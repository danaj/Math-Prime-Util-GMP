
/* GMP version of Racing SQUFOF for up to 126-bit inputs.
 * Oct 2017 - Dana Jacobsen
 *
 * Based heavily on Ben Buhrow's racing SQUFOF implementation.
 * All factoring operations use 64-bit unsigned longs so it's quite fast.
 *
 * Realistically it has decent performance up to about 80 bits.
 * From 54 to 80 it is faster to first try p-1 to B1 = 1k-16k.
 *
 * As of 2017, the fastest method for 64-bit is x86-64 Pollard-Rho.
 * After that, a tinyqs such as Jason P's cofactorize-siqs seems fastest.
 */

#include <gmp.h>
#include <math.h>
#include "ptypes.h"
#include "squfof126.h"

#define TEST_FOR_2357(n, f) \
  { \
    if (mpz_divisible_ui_p(n, 2)) { mpz_set_ui(f, 2); return 1; } \
    if (mpz_divisible_ui_p(n, 3)) { mpz_set_ui(f, 3); return 1; } \
    if (mpz_divisible_ui_p(n, 5)) { mpz_set_ui(f, 5); return 1; } \
    if (mpz_divisible_ui_p(n, 7)) { mpz_set_ui(f, 7); return 1; } \
    if (mpz_cmp_ui(n, 121) < 0) { return 0; } \
  }

/* Pick type for 64-bit core, plus methods to get/set from GMP */

#if HAVE_STD_U64
#define SQUFOF_TYPE uint64_t
#elif BITS_PER_WORD == 64
#define SQUFOF_TYPE UV
#else
#define SQUFOF_TYPE unsigned long long
#endif

static INLINE SQUFOF_TYPE mpz_get64(mpz_t n) {
  SQUFOF_TYPE v = mpz_getlimbn(n,0);
  if (GMP_LIMB_BITS < 64 || sizeof(mp_limb_t) < sizeof(SQUFOF_TYPE))
    v |= ((SQUFOF_TYPE)mpz_getlimbn(n,1)) << 32;
  return v;
}
static INLINE void mpz_set64(mpz_t n, SQUFOF_TYPE v) {
  if (v <= 0xFFFFFFFFUL || sizeof(unsigned long int) >= sizeof(SQUFOF_TYPE)) {
    mpz_set_ui(n, v);
  } else {
    uint32_t upper = (v >> 32), lower = v & 0xFFFFFFFFUL;
    mpz_set_ui(n, upper);
    mpz_mul_2exp(n, n, 32);
    mpz_add_ui(n, n, lower);
  }
}

typedef struct
{
  int valid;
  SQUFOF_TYPE P;
  SQUFOF_TYPE bn;
  SQUFOF_TYPE Qn;
  SQUFOF_TYPE Q0;
  SQUFOF_TYPE b0;
  SQUFOF_TYPE it;
  SQUFOF_TYPE imax;
  SQUFOF_TYPE mult;
} mult_t;

/* Return 0 or factor */
static SQUFOF_TYPE squfof_unit(mpz_t n, mult_t* mult_save, mpz_t t)
{
  SQUFOF_TYPE imax,i,j,Q0,Qn,bn,b0,P,bbn,Ro,S,So,t1,t2;

  P = mult_save->P;
  bn = mult_save->bn;
  Qn = mult_save->Qn;
  Q0 = mult_save->Q0;
  b0 = mult_save->b0;
  i  = mult_save->it;
  imax = i + mult_save->imax;

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
      t2 = Qn & 127;
      if (!((t2*0x8bc40d7d) & (t2*0xa1e2f5d1) & 0x14020a)) {
        t1 = (uint32_t) sqrt(Qn);
        if (Qn == t1*t1)
          break;
      }

      SQUARE_SEARCH_ITERATION;  /* Odd iteration */
    }
    S = t1;
    mult_save->it = i;

    /* Reduce to G0 */
    Ro = P + S*((b0 - P)/S);
    /* So = (n - (UV)Ro*(UV)Ro)/(UV)S; */
    mpz_set64(t, Ro);
    mpz_mul(t, t, t);
    mpz_sub(t, n, t);
    mpz_div_ui(t, t, S);   /* S is 32-bit so this is ok */
    So = mpz_get64(t);
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

    mpz_set64(t, Ro);
    mpz_gcd(t, t, n);
    t1 = mpz_get64(t);
    if (t1 > 1)
      return t1;
  }
}

/* Gower and Wagstaff 2008:
 *    http://www.ams.org/journals/mcom/2008-77-261/S0025-5718-07-02010-8/
 * Section 5.3.  I've added some with 13,17,19.  Sorted by F(). */
static const SQUFOF_TYPE squfof_multipliers[] =
#if 0
  { 3*5*7*11, 3*5*7, 3*5*11, 3*5, 3*7*11, 3*7, 5*7*11, 5*7,
    3*11,     3,     5*11,   5,   7*11,   7,   11,     1   };
#endif
#if 0
  { 3*5*7*11, 3*5*7,  3*5*7*11*13, 3*5*7*13, 3*5*7*11*17, 3*5*11,
    3*5*7*17, 3*5,    3*5*7*11*19, 3*5*11*13,3*5*7*19,    3*5*7*13*17,
    3*5*13,   3*7*11, 3*7,         5*7*11,   3*7*13,      5*7,
    3*5*17,   5*7*13, 3*5*19,      3*11,     3*7*17,      3,
    3*11*13,  5*11,   3*7*19,      3*13,     5,           5*11*13,
    5*7*19,   5*13,   7*11,        7,        3*17,        7*13,
    11,       1 };
#endif
#if 1
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
#endif
#define NSQUFOF_MULT (sizeof(squfof_multipliers)/sizeof(squfof_multipliers[0]))

int squfof126(mpz_t n, mpz_t f, UV rounds)
{
  mpz_t t, nn64;
  mult_t mult_save[NSQUFOF_MULT];
  SQUFOF_TYPE i, mult, f64, f64red, sqrtnn64, rounds_done = 0;
  int mults_racing = NSQUFOF_MULT;
  const uint32_t max_bits = 2 * sizeof(SQUFOF_TYPE)*8 - 2;

  if (sizeof(SQUFOF_TYPE) <  8 || mpz_sizeinbase(n,2) > max_bits) {
    mpz_set(f, n);
    return 0;
  }
  TEST_FOR_2357(n, f);

  for (i = 0; i < NSQUFOF_MULT; i++) {
    mult_save[i].valid = -1;
    mult_save[i].it = 0;
  }
  mpz_init(t);  mpz_init(nn64);

  /* Process the multipliers a little at a time: 0.5 * (n*mult)^1/5 */
  while (mults_racing > 0 && rounds_done < rounds) {
    for (i = 0; i < NSQUFOF_MULT && rounds_done < rounds; i++) {
      if (mult_save[i].valid == 0)  continue;
      mult = squfof_multipliers[i];
      mpz_mul_ui(nn64, n, mult);
      if (mult_save[i].valid == -1) {
        if (mpz_sizeinbase(nn64,2) > max_bits) {
          mult_save[i].valid = 0; /* This multiplier would overflow */
          mults_racing--;
          continue;
        }
        mpz_sqrt(t,nn64);
        sqrtnn64 = mpz_get64(t);
        mpz_mul(t,t,t);
        mpz_sub(t, nn64, t);
        mult_save[i].valid = 1;
        mult_save[i].Q0    = 1;
        mult_save[i].b0    = sqrtnn64;
        mult_save[i].P     = sqrtnn64;
        mult_save[i].Qn    = mpz_get64(t);  /* nn64 - isqrt(nn64)^2 */
        if (mult_save[i].Qn == 0) {
          mpz_clear(t); mpz_clear(nn64);
          mpz_set64(f, sqrtnn64);
          return 1;  /* nn64 is a perfect square */
        }
        mpz_root(t, nn64, 5);
        mult_save[i].bn    = (2 * sqrtnn64) / mult_save[i].Qn; /* n < 127-bit */
        mult_save[i].it    = 0;
        mult_save[i].mult  = mult;
        mult_save[i].imax  = (SQUFOF_TYPE) (0.5 * mpz_get64(t));
        if (mult_save[i].imax < 20)
          mult_save[i].imax = 20;
      }
      if (mults_racing == 1 || mult_save[i].imax > (rounds-rounds_done))
        mult_save[i].imax = (rounds - rounds_done);
      f64 = squfof_unit(nn64, &mult_save[i], t);
      if (f64 > 1) {
        mpz_set64(t, f64);
        f64red = f64 / mpz_gcd_ui(t, t, mult);
        if (f64red > 1) {
          mpz_clear(t); mpz_clear(nn64);
          mpz_set64(f, f64red);
          return 1;
        }
        /* Found trivial factor.  Quit working with this multiplier. */
        mult_save[i].valid = 0;
      }
      if (mult_save[i].valid == 0)
        mults_racing--;
      rounds_done += mult_save[i].imax;   /* Assume we did all rounds */
    }
  }

  /* No factors found */
  mpz_clear(t);
  mpz_clear(nn64);
  mpz_set(f, n);
  return 0;
}
