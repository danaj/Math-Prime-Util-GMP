
/* GMP version of Racing SQUFOF for up to 126-bit inputs.
 * Oct 2017-2026, Dana Jacobsen
 *
 * Version 1.1.0, October 2026
 *
 * Based heavily on Ben Buhrow's racing SQUFOF implementation.
 * The SQUFOF recurrences use native unsigned 64-bit arithmetic.
 *
 * Uses lazy multiplier initialization and interleaved racing, with four
 * lanes by default, arranged to expose compiler vectorization opportunities.
 *
 * The 126-bit bound is a limit for this implementation, not a performance
 * recommendation.  Crossovers with rho, ECM, and SIQS depend on input shape,
 * work budget, their implementations, and platform.
 */

#include <gmp.h>
#include <limits.h>
#include <math.h>
#include "squfof126.h"

#define TEST_FOR_2357(n, f) \
  { \
    if (mpz_divisible_ui_p(n, 2)) { mpz_set_ui(f, 2); return 1; } \
    if (mpz_divisible_ui_p(n, 3)) { mpz_set_ui(f, 3); return 1; } \
    if (mpz_divisible_ui_p(n, 5)) { mpz_set_ui(f, 5); return 1; } \
    if (mpz_divisible_ui_p(n, 7)) { mpz_set_ui(f, 7); return 1; } \
    if (mpz_cmp_ui(n, 121) < 0) { return 0; } \
  }

/* Low 64 magnitude bits, independent of unsigned long and GMP limb width.
 * GMP_NUMB_BITS also accounts for any unused nail bits in a limb. */
static inline uint64_t mpz_get64(const mpz_t n) {
  uint64_t v = (uint64_t)mpz_getlimbn(n, 0);
#if GMP_NUMB_BITS < 64
  unsigned shift;
  mp_size_t i;
  for (i = 1, shift = GMP_NUMB_BITS; shift < 64; i++, shift += GMP_NUMB_BITS)
    v |= (uint64_t)mpz_getlimbn(n, i) << shift;
#endif
  return v;
}

/* GMP's _ui arguments are unsigned long, which may have only 32 bits even
 * when uint64_t and GMP limbs have 64. Keep this conversion local. */
static inline void mpz_set64(mpz_t n, uint64_t v) {
#if ULONG_MAX >= UINT64_MAX
  mpz_set_ui(n, (unsigned long)v);
#else
  if (v <= ULONG_MAX) {
    mpz_set_ui(n, (unsigned long)v);
  } else {
    mpz_set_ui(n, (unsigned long)(v >> 32));
    mpz_mul_2exp(n, n, 32);
    mpz_add_ui(n, n, (unsigned long)(uint32_t)v);
  }
#endif
}

typedef struct
{
  int valid;
  uint64_t P;
  uint64_t bn;
  uint64_t Qn;
  uint64_t Q0;
  uint64_t b0;
  uint64_t it;
  uint64_t imax;
  uint64_t mult;
} mult_t;

/* Bounded groups for autovectorization, with scalar cleanup for short groups. */
#ifndef SQUFOF_RACE_LANES
#define SQUFOF_RACE_LANES 4
#endif
#if SQUFOF_RACE_LANES != 2 && SQUFOF_RACE_LANES != 4 && SQUFOF_RACE_LANES != 8
#error SQUFOF_RACE_LANES must be 2 or 4 or 8
#endif

static inline unsigned squfof_ctz64(uint64_t n)
{
#if defined(__GNUC__) || defined(__clang__)
  return (unsigned)__builtin_ctzll((unsigned long long)n);
#else
  unsigned count = 0;
  while ((n & 1U) == 0) { n >>= 1; count++; }
  return count;
#endif
}

static uint64_t squfof_gcd64(uint64_t a, uint64_t b)
{
  unsigned shift;
  uint64_t tmp;
  if (a == 0) return b;
  if (b == 0) return a;
  /* Reduce disparate widths once before the division-free binary loop. */
  if (a < b) { tmp = a; a = b; b = tmp; }
  a %= b;
  if (a == 0) return b;
  shift = squfof_ctz64(a | b);
  a >>= squfof_ctz64(a);
  do {
    b >>= squfof_ctz64(b);
    if (a > b) { tmp = a; a = b; b = tmp; }
    b -= a;
  } while (b != 0);
  return a << shift;
}

/* Scalar reverse search and extraction, shared by the lane race and cleanup.
 * imax is the absolute forward limit for this visit, as in the scalar code. */
static uint64_t squfof_extract(const mpz_t n, uint64_t n64,
                                 mult_t* mult_save, uint64_t P,
                                 uint64_t Q0, uint64_t S,
                                 uint64_t imax, mpz_t t)
{
  uint64_t j,bbn,Ro,So,t1,t2,k,g;
  const uint64_t b0 = mult_save->b0;

  /* Reduce to G0 */
  k = (b0 - P)/S;
  Ro = P + S*k;
  /* D=P*P+Q0*S*S, so (D-Ro*Ro)/S=S*Q0-k*(P+Ro).
   * The true result is in (0,2*b0], below 2^64.  Unsigned wrapping
   * of the intermediate products therefore gives that exact result. */
  So = S*Q0 - k*(P+Ro);
  bbn = (b0+Ro)/So;

#define SYMMETRY_POINT_ITERATION \
      t1 = Ro; \
      Ro = bbn*So - Ro; \
      if (Ro == t1) break; \
      t2 = So; \
      So = S + bbn*(t1-Ro); \
      S = t2; \
      bbn = (b0+Ro)/So;

  /* Search for symmetry point, occurs at approximately i/2 */
  j = 0;
  while (1) {
    SYMMETRY_POINT_ITERATION;
    SYMMETRY_POINT_ITERATION;
    SYMMETRY_POINT_ITERATION;
    SYMMETRY_POINT_ITERATION;
    if (j++ > imax) {
      mult_save->valid = 0;
      return 0;
    }
  }
#undef SYMMETRY_POINT_ITERATION

  /* gcd(Ro,m*n)/gcd(gcd(Ro,m*n),m) = gcd(Ro/g,n), where
   * g=gcd(Ro,m).  Keep multiplier-only outcomes distinct from gcd 1. */
  g = squfof_gcd64(Ro, mult_save->mult);
  /* Racing can expose a full-input GCD before a useful multiplier wins.
   * Treat it as a non-split and keep searching this live state. */
  if (n64 != 0) {
    t1 = squfof_gcd64(Ro/g, n64);
    if (t1 == n64)
      return 0;
  } else {
    mpz_set64(t, Ro/g);
    mpz_gcd(t, t, n);
    if (mpz_cmp(t, n) == 0)
      return 0;
    t1 = mpz_get64(t);
  }
  if (t1 > 1)
    return t1;
  if (g > 1)
    mult_save->valid = 0;
  return 0;
}

/* Return 0 or a factor with the multiplier already removed.
 * n is the original input; n64 is zero when n needs the GMP GCD path.
 * An odd starting iteration has already had its square test/extraction. */
static uint64_t squfof_unit(const mpz_t n, uint64_t n64,
                              mult_t* mult_save, uint64_t imax, mpz_t t)
{
  uint64_t i,Q0,Qn,bn,b0,P,t1,t2,f64;

  P = mult_save->P;
  bn = mult_save->bn;
  Qn = mult_save->Qn;
  Q0 = mult_save->Q0;
  b0 = mult_save->b0;
  i  = mult_save->it;

#define SQUARE_SEARCH_ITERATION \
      t1 = P; \
      P = bn*Qn - P; \
      t2 = Qn; \
      Qn = Q0 + bn*(t1-P); \
      Q0 = t2; \
      bn = (b0 + P) / Qn; \
      i++;

  while (1) {
    if (i & 0x1) {
      SQUARE_SEARCH_ITERATION;
    }
    /* i is now even */
    while (1) {
      /* We need to know P, bn, Qn, Q0, iteration count, i  from prev */
      if (i >= imax) {
        /* save state and try another multiplier. */
        mult_save->P = P;
        mult_save->bn = bn;
        mult_save->Qn = Qn;
        mult_save->Q0 = Q0;
        mult_save->it = i;
        return 0;
      }

      SQUARE_SEARCH_ITERATION;  /* Even iteration */

      /* Check if Qn is a perfect square */
      /* This residue filter is not a win at this location. Our real-cofactor
       * samples fail the mod-64 test only about 22% of the time; removing it
       * reduced splitter CPU time by about 5% on M1 Pro, despite extra roots.
       * Keep it commented out for comparison on other machines.
      if (!((UINT64_C(1) << (Qn & 63)) & UINT64_C(0xfdfdfdedfdfcfdec))) {
      */
      /* The FP root may round to 2^32; keep it wide for the square test. */
      t1 = (uint64_t) sqrt((double) Qn);
      if (Qn == t1*t1)
        break;
      /* } */

      SQUARE_SEARCH_ITERATION;  /* Odd iteration */
    }
    mult_save->it = i;
    f64 = squfof_extract(n, n64, mult_save, P, Q0, t1, imax, t);
    if (f64 > 1 || !mult_save->valid)
      return f64;
  }
}
#undef SQUARE_SEARCH_ITERATION

/* Full groups race until a lane exhausts its batch or retires.  The SoA
 * scratch and fixed-length loops expose independent exact recurrences;
 * reverse searches and all GMP operations remain outside those loops.
 * Short groups and uneven tails use the same scalar recurrence as before. */
static uint64_t squfof_race(const mpz_t n, uint64_t n64,
                              mult_t *states[SQUFOF_RACE_LANES],
                              unsigned count, mpz_t t)
{
  uint64_t P[SQUFOF_RACE_LANES], bn[SQUFOF_RACE_LANES];
  uint64_t Qn[SQUFOF_RACE_LANES], Q0[SQUFOF_RACE_LANES];
  uint64_t b0[SQUFOF_RACE_LANES], root[SQUFOF_RACE_LANES];
  uint64_t start[SQUFOF_RACE_LANES], limit[SQUFOF_RACE_LANES];
  uint64_t used = 0, common = states[0]->imax, f64;
  unsigned lane, squares, retired;

  for (lane = 0; lane < count; lane++) {
    start[lane] = states[lane]->it;
    limit[lane] = start[lane] + states[lane]->imax;
    if (states[lane]->imax < common)
      common = states[lane]->imax;
  }
  if (count == SQUFOF_RACE_LANES) {
    for (lane = 0; lane < SQUFOF_RACE_LANES; lane++) {
      P[lane] = states[lane]->P;
      bn[lane] = states[lane]->bn;
      Qn[lane] = states[lane]->Qn;
      Q0[lane] = states[lane]->Q0;
      b0[lane] = states[lane]->b0;
    }

#define RACE_ITERATION \
      uint64_t oldP = P[lane], oldQn = Qn[lane]; \
      P[lane] = bn[lane]*Qn[lane] - P[lane]; \
      Qn[lane] = Q0[lane] + bn[lane]*(oldP-P[lane]); \
      Q0[lane] = oldQn; \
      bn[lane] = (b0[lane]+P[lane]) / Qn[lane];

    /* Every live state starts/resumes even.  Each pass tests the odd state
     * after one step and then completes its pair, just like squfof_unit. */
    while (used < common) {
      for (lane = 0; lane < SQUFOF_RACE_LANES; lane++) {
        RACE_ITERATION;
      }
      /* Keep the independent conversion/root stage separate from division. */
#if defined(__clang__)
      /* Otherwise this short loop unrolls before Clang considers vectors. */
#pragma clang loop unroll(disable) vectorize(enable) vectorize_width(2) interleave_count(2)
#endif
      for (lane = 0; lane < SQUFOF_RACE_LANES; lane++)
        root[lane] = (uint64_t) sqrt((double) Qn[lane]);
      used++;
      squares = 0;
      for (lane = 0; lane < SQUFOF_RACE_LANES; lane++)
        squares |= (Qn[lane] == root[lane]*root[lane]) << lane;
      if (squares) {
        retired = 0;
        for (lane = 0; lane < SQUFOF_RACE_LANES; lane++) {
          if (!(squares & (1U << lane)))
            continue;
          states[lane]->it = start[lane] + used;
          f64 = squfof_extract(n, n64, states[lane], P[lane], Q0[lane],
                               root[lane], limit[lane], t);
          if (f64 > 1)
            return f64;
          retired |= !states[lane]->valid;
        }
        if (retired)
          break;
      }
      for (lane = 0; lane < SQUFOF_RACE_LANES; lane++) {
        RACE_ITERATION;
      }
      used++;
    }
#undef RACE_ITERATION

    for (lane = 0; lane < SQUFOF_RACE_LANES; lane++) {
      states[lane]->P = P[lane];
      states[lane]->bn = bn[lane];
      states[lane]->Qn = Qn[lane];
      states[lane]->Q0 = Q0[lane];
      states[lane]->it = start[lane] + used;
    }
  }
  for (lane = 0; lane < count; lane++) {
    if (!states[lane]->valid)
      continue;
    /* Odd states have already been square-tested.  Even an exhausted odd
     * tail must finish its pair before the scalar loop checks its limit. */
    f64 = squfof_unit(n, n64, states[lane], limit[lane], t);
    if (f64 > 1)
      return f64;
  }
  return 0;
}

/* Gower and Wagstaff 2008:
 *    http://www.ams.org/journals/mcom/2008-77-261/S0025-5718-07-02010-8/
 * Section 5.3.  I've added some with 13,17,19.  Sorted by F(). */
/* These fit GMP's unsigned-long multiplier argument on 32-bit hosts too. */
static const uint32_t squfof_multipliers[] =
  { 33*1680, 11*1680, 66*1680,  3*1680,  2*1680,  6*1680, 22*1680, 78*1680,
     1*1680, 26*1680, 39*1680, 13*1680,102*1680, 30*1680, 34*1680, 10*1680,
    15*1680, 51*1680,  5*1680, 57*1680, 17*1680, 19*1680,
    3*5*7*11, 3*5*7,  3*5*7*11*13, 3*5*7*13, 3*5*7*11*17, 3*5*11,
    3*5*7*17, 3*5,    3*5*7*11*19, 3*5*11*13,3*5*7*19,    3*5*7*13*17,
    3*5*13,   3*7*11, 3*7,         5*7*11,   3*7*13,      5*7,
    3*5*17,   5*7*13, 3*5*19,      3*11,     3*7*17,      3,
    3*11*13,  5*11,   3*7*19,      3*13,     5,           5*11*13,
    5*7*19,   5*13,   7*11,        7,        3*17,        7*13,
    11,       1 };
#define NSQUFOF_MULT (sizeof(squfof_multipliers)/sizeof(squfof_multipliers[0]))

int squfof126(const mpz_t n, mpz_t f, uint64_t rounds)
{
  mpz_t t, nn64;
  mult_t mult_save[NSQUFOF_MULT];
  mult_t *states[SQUFOF_RACE_LANES];
  size_t i;
  uint32_t mult;
  uint64_t f64, sqrtnn64, n64, rounds_done = 0;
  uint64_t reserved, available;
  unsigned count, lane;
  int mults_racing = NSQUFOF_MULT;
  int first_pass = 1;
  const uint32_t max_bits = 126;
  const size_t nbits = mpz_sizeinbase(n, 2);
  const double batch_scale = nbits <= 60 ? 0.3 : 0.5;

  if (nbits > max_bits) {
    mpz_set(f, n);
    return 0;
  }
  TEST_FOR_2357(n, f);
  n64 = nbits <= 64 ? mpz_get64(n) : 0;

  mpz_init(t);  mpz_init(nn64);

  /* Short batches help through 60 bits; larger inputs favor fewer long tails. */
  while (mults_racing > 0 && rounds_done < rounds) {
    for (i = 0; i < NSQUFOF_MULT && rounds_done < rounds;) {
      count = 0;
      reserved = 0;
      available = rounds - rounds_done;
      /* Reserve normal batches in multiplier order.  Their sum, rather
       * than the number of parallel steps, consumes the shared budget. */
      for (; i < NSQUFOF_MULT && count < SQUFOF_RACE_LANES
             && reserved < available; i++) {
        if (!first_pass && mult_save[i].valid == 0)  continue;
        if (first_pass) {
          mult = squfof_multipliers[i];
          mpz_mul_ui(nn64, n, (unsigned long)mult);
          if (mpz_sizeinbase(nn64,2) > max_bits) {
            mult_save[i].valid = 0; /* This multiplier would overflow */
            mults_racing--;
            continue;
          }
          mpz_sqrt(t,nn64);
          sqrtnn64 = mpz_get64(t);
          mult_save[i].valid = 1;
          mult_save[i].Q0    = 1;
          mult_save[i].b0    = sqrtnn64;
          mult_save[i].P     = sqrtnn64;
          /* 0 <= D-b0*b0 <= 2*b0 < 2^64: the low-word subtraction
           * is exact even when b0*b0 wraps. */
          mult_save[i].Qn    = mpz_get64(nn64) - sqrtnn64*sqrtnn64;
          if (mult_save[i].Qn == 0) {
            mpz_clear(t); mpz_clear(nn64);
            mpz_set64(f, sqrtnn64);
            return 1;  /* nn64 is a perfect square */
          }
          mult_save[i].bn    = (2 * sqrtnn64) / mult_save[i].Qn;
          mult_save[i].it    = 0;
          mult_save[i].mult  = mult;
          /* The fifth root only sizes a work batch; no algebra needs it exact.
           * Preserve the integer estimate before applying the scale. */
          mult_save[i].imax  = (uint64_t)
            (batch_scale * (uint64_t) pow(mpz_get_d(nn64), 0.2));
          if (mult_save[i].imax < 20)
            mult_save[i].imax = 20;
        }
        if (mults_racing == 1 || mult_save[i].imax > available-reserved)
          mult_save[i].imax = available-reserved;
        states[count++] = &mult_save[i];
        reserved += mult_save[i].imax;
      }
      if (count == 0)
        continue;
      f64 = squfof_race(n, n64, states, count, t);
      if (f64 > 1) {
        mpz_clear(t); mpz_clear(nn64);
        mpz_set64(f, f64);
        return 1;
      }
      for (lane = 0; lane < count; lane++)
        if (!states[lane]->valid)
          mults_racing--;
      rounds_done += reserved;   /* Charge each complete reserved batch. */
    }
    /* A later sweep can only start after every record was visited:
     * exhausting the budget also ends the outer loop. */
    first_pass = 0;
  }

  /* No factors found */
  mpz_clear(t);
  mpz_clear(nn64);
  mpz_set(f, n);
  return 0;
}
