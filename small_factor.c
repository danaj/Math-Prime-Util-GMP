
/* Racing SQUFOF, based in part on Ben Buhrow's code. */

#include <math.h>

#include "ptypes.h"
#include "small_factor.h"

static UV isqrt(UV n) {
  UV root;
#if BITS_PER_WORD == 32
  if (n >= UVCONST(4294836225)) return UVCONST(65535);
#else
  if (n >= UVCONST(18446744065119617025)) return UVCONST(4294967295);
#endif
  root = (UV) sqrt((double)n);
  while (root*root > n)  root--;
  while ((root+1)*(root+1) <= n)  root++;
  return root;
}

/* Return 0 if n is not a perfect square.  Set sqrtn to int(sqrt(n)) if so.
 * See Math::Prime::Util factor.c for more info on this.
 */
static int is_perfect_square(UV n, UV* sqrtn)
{
  UV m;
  m = n & 127;
  if ((m*0x8bc40d7d) & (m*0xa1e2f5d1) & 0x14020a)  return 0;
  /* This cuts out another 80%: */
  m = n % 63; if ((m*0x3d491df7) & (m*0xc824a9f9) & 0x10f14008) return 0;
  /* m = n % 25; if ((m*0x1929fc1b) & (m*0x4c9ea3b2) & 0x51001005) return 0; */
  m = isqrt(n);
  if (n != m*m) return 0;
  if (sqrtn != 0) *sqrtn = m;
  return 1;
}

static UV gcd_ui(UV x, UV y) {
  UV t;
  if (y < x) { t = x; x = y; y = t; }
  while (y > 0) {
    t = y;  y = x % y;  x = t;  /* y1 <- x0 % y0 ; x1 <- y0 */
  }
  return x;
}

typedef struct
{
  int valid;
  UV P;
  UV bn;
  UV Qn;
  UV Q0;
  UV b0;
  UV it;
  UV imax;
} mult_t;

/* N < 2^63 (or 2^31).  Returns 0 or a factor */
static UV squfof_unit(UV n, mult_t* mult_save)
{
  UV imax,i,Q0,b0,Qn,bn,P,bbn,Ro,S,So,t1,t2;

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
    int j = 0;
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

      SQUARE_SEARCH_ITERATION;

      /* Even iteration.  Check for square: Qn = S*S */
      if (is_perfect_square(Qn, &S))
        break;

      /* Odd iteration. */
      SQUARE_SEARCH_ITERATION;
    }
    /* printf("found square %lu after %lu iterations with mult %d\n", Qn, i, mult_save->mult); */

    /* Reduce to G0 */
    Ro = P + S*((b0 - P)/S);
    t1 = Ro;
    So = (n - t1*t1)/S;
    bbn = (b0+Ro)/So;

    /* Search for symmetry point */
#define SYMMETRY_POINT_ITERATION \
      t1 = Ro; \
      Ro = bbn*So - Ro; \
      t2 = So; \
      So = S + bbn*(t1-Ro); \
      S = t2; \
      bbn = (b0+Ro)/So; \
      if (Ro == t1) break;

    j = 0;
    while (1) {
      SYMMETRY_POINT_ITERATION;
      SYMMETRY_POINT_ITERATION;
      SYMMETRY_POINT_ITERATION;
      SYMMETRY_POINT_ITERATION;
      if (j++ > 2000000) {
         mult_save->valid = 0;
         return 0;
      }
    }

    t1 = gcd_ui(Ro, n);
    if (t1 > 1)
      return t1;
  }
}

/* Gower and Wagstaff 2008:
 *    http://www.ams.org/journals/mcom/2008-77-261/S0025-5718-07-02010-8/
 * Section 5.3.  I've added some with 13,17,19.  Sorted by F(). */
static const UV squfof_multipliers[] =
  /* { 3*5*7*11, 3*5*7, 3*5*11, 3*5, 3*7*11, 3*7, 5*7*11, 5*7,
       3*11,     3,     5*11,   5,   7*11,   7,   11,     1   }; */
  { 3*5*7*11, 3*5*7,  3*5*7*11*13, 3*5*7*13, 3*5*7*11*17, 3*5*11,
    3*5*7*17, 3*5,    3*5*7*11*19, 3*5*11*13,3*5*7*19,    3*5*7*13*17,
    3*5*13,   3*7*11, 3*7,         5*7*11,   3*7*13,      5*7,
    3*5*17,   5*7*13, 3*5*19,      3*11,     3*7*17,      3,
    3*11*13,  5*11,   3*7*19,      3*13,     5,           5*11*13,
    5*7*19,   5*13,   7*11,        7,        3*17,        7*13,
    11,       1 };
#define NSQUFOF_MULT (sizeof(squfof_multipliers)/sizeof(squfof_multipliers[0]))

int racing_squfof_factor(UV n, UV *factors, UV rounds)
{
  const UV big2 = UV_MAX;
  mult_t mult_save[NSQUFOF_MULT];
  int still_racing;
  UV i, nn64, mult, f64;
  UV rounds_done = 0;

  /* Caller should have handled these trivial cases */
  MPUassert( (n >= 3) && ((n%2) != 0) , "bad n in squfof_factor");

  /* Too big */
  if (n > big2) {
    factors[0] = n;  return 1;
  }

  for (i = 0; i < NSQUFOF_MULT; i++)
    mult_save[i].valid = -1;

  /* Process the multipliers a little at a time: 0.33*(n*mult)^1/4: 20-20k */
  do {
    still_racing = 0;
    for (i = 0; i < NSQUFOF_MULT; i++) {
      if (mult_save[i].valid == 0)  continue;
      mult = squfof_multipliers[i];
      nn64 = n * mult;
      if (mult_save[i].valid == -1) {
        if ((big2 / mult) < n) {
          mult_save[i].valid = 0; /* This multiplier would overflow 64-bit */
          continue;
        }
        mult_save[i].valid = 1;
        mult_save[i].b0 = isqrt(nn64);
        mult_save[i].imax = (UV) (sqrt(mult_save[i].b0) / 16);
        if (mult_save[i].imax < 20)     mult_save[i].imax = 20;
        if (mult_save[i].imax > rounds) mult_save[i].imax = rounds;
        mult_save[i].Q0 = 1;
        mult_save[i].P  = mult_save[i].b0;
        mult_save[i].Qn = nn64 - (mult_save[i].b0 * mult_save[i].b0);
        if (mult_save[i].Qn == 0) {
          factors[0] = mult_save[i].b0;
          factors[1] = n / mult_save[i].b0;
          MPUassert( factors[0] * factors[1] == n , "incorrect factoring");
          return 2;
        }
        mult_save[i].bn = (mult_save[i].b0 + mult_save[i].P) / mult_save[i].Qn;
        mult_save[i].it = 0;
      }
      f64 = squfof_unit(nn64, &mult_save[i]);
      if (f64 > 1) {
        if (f64 != mult) {
          f64 /= gcd_ui(f64, mult);
          if (f64 != 1) {
            factors[0] = f64;
            factors[1] = n / f64;
            MPUassert( factors[0] * factors[1] == n , "incorrect factoring");
            return 2;
          }
        }
        /* Found trivial factor.  Quit working with this multiplier. */
        mult_save[i].valid = 0;
      }
      if (mult_save[i].valid == 1)
        still_racing = 1;
      rounds_done += mult_save[i].imax;
      if (rounds_done >= rounds)
        break;
    }
  } while (still_racing && rounds_done < rounds);

  /* No factors found */
  factors[0] = n;
  return 1;
}
