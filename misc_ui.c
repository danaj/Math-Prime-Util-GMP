/* Native functions, either ui or UV.
 *
 * Utility as well as some more specialized.
 */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <gmp.h>
#include <math.h>

#include "ptypes.h"
#include "prime_iterator.h"
#include "gmp_main.h"

#define FUNC_ipow 1
#define FUNC_isqrt 1
#define FUNC_icbrt 1
#define FUNC_is_perfect_square
#include "misc_ui.h"


#if BITS_PER_WORD == 64
  #define MPU_MAX_POW3 40
  static const uint32_t root_max[1+MPU_MAX_POW3] = {0,0,4294967295U,2642245,65535,7131,1625,565,255,138,84,56,40,30,23,19,15,13,11,10,9,8,7,6,6,5,5,5,4,4,4,4,3,3,3,3,3,3,3,3,3};
  static const unsigned char _debruijn64[64] = {
    63, 0,58, 1,59,47,53, 2, 60,39,48,27,54,33,42, 3, 61,51,37,40,49,18,28,20,
    55,30,34,11,43,14,22, 4, 62,57,46,52,38,26,32,41, 50,36,17,19,29,10,13,21,
    56,45,25,31,35,16, 9,12, 44,24,15, 8,23, 7, 6, 5 };
  unsigned int log2_ui(UV n) {
    if (n == 0) return 0;
    n |= n >> 1;   n |= n >> 2;   n |= n >> 4;
    n |= n >> 8;   n |= n >> 16;  n |= n >> 32;
    return _debruijn64[((n-(n>>1))*UVCONST(0x07EDD5E59A4E28C2)) >> 58];
  }
#else
  #define MPU_MAX_POW3 20
  static const uint32_t root_max[1+MPU_MAX_POW3] = {0,0,65535,1625,255,84,40,23,15,11,9,7,6,5,4,4,3,3,3,3,3};
  static const unsigned char _lead_debruijn32[32] = {
    0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30,
    8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31 };
  unsigned int log2_ui(UV n) {
    if (n == 0) return 0;
    n |= n >> 1;   n |= n >> 2;   n |= n >> 4;   n |= n >> 8;   n |= n >> 16;
    return _lead_debruijn32[(n * UVCONST(0x07C4ACDD)) >> 27];
  }
#endif

UV rootint_ui(UV n, UV k)
{
  UV lo, hi, max;
  if (k == 0) return 0;
  if (k == 1) return n;
  if (k == 2) return isqrt(n);
  if (k == 3) return icbrt(n);

  /* Bracket between powers of 2, but never exceed max power so ipow works */
  max = 1 + ((k > MPU_MAX_POW3) ? 2 : root_max[k]);
  lo = UVCONST(1) << (log2_ui(n)/k);
  hi = ((lo*2) < max) ? lo*2 : max;

  /* Binary search */
  while (lo < hi) {
    UV mid = lo + (hi-lo)/2;
    if (ipow(mid,k) <= n) lo = mid+1;
    else                  hi = mid;
  }
  return lo-1;
}
UV logint_ui(UV n, UV b)
{
  /* UV e;  for (e=0; n; n /= b) e++;  return e-1; */
  UV v, e = 0;
  if (b == 2)
    return log2_ui(n);
  if (b > n)
    return 0;
  if (n > UV_MAX/b) {
    n /= b;
    e = 1;
  }
  for (v = b; v <= n; v *= b)
    e++;
  return e;
}

/* Modular inversion: invert a mod p.
 * This implementation from William Hart, using extended gcd.
 */
unsigned long modinverse(unsigned long a, unsigned long p)
{
  long u1 = 1;
  long u3 = a;
  long v1 = 0;
  long v3 = p;
  long t1 = 0;
  long t3 = 0;
  long quot;
  while (v3) {
    quot = u3 - v3;
    if (u3 < (v3<<2)) {
      if (quot < v3) {
        if (quot < 0) {
          t1 = u1; u1 = v1; v1 = t1;
          t3 = u3; u3 = v3; v3 = t3;
        } else {
          t1 = u1 - v1; u1 = v1; v1 = t1;
          t3 = u3 - v3; u3 = v3; v3 = t3;
        }
      } else if (quot < (v3<<1)) {
        t1 = u1 - (v1<<1); u1 = v1; v1 = t1;
        t3 = u3 - (v3<<1); u3 = v3; v3 = t3;
      } else {
        t1 = u1 - v1*3; u1 = v1; v1 = t1;
        t3 = u3 - v3*3; u3 = v3; v3 = t3;
      }
    } else {
      quot = u3 / v3;
      t1 = u1 - v1*quot; u1 = v1; v1 = t1;
      t3 = u3 - v3*quot; u3 = v3; v3 = t3;
    }
 }
 if (u1 < 0) u1 += p;
 return u1;
}

UV next_prime_ui(UV n)
{
  mpz_t t;
  if (n < 13)
    return (n<3) ? 3 : (n<5) ? 5 : (n<7) ? 7 : (n<11) ? 11 : 13;
  mpz_init_set_ui(t, n);
  _GMP_next_prime(t);
  n = mpz_get_ui(t);
  mpz_clear(t);
  return n;
}

/******************************************************************************/

#define P_GT_LO(f,p,lo)  ( ((f)>=(lo)) ? (f) : (lo)+(((p)-((lo)%(p)))%(p)) )

/* Return a char array with lo-hi+1 elements. mu[k-lo] = µ(k) for k = lo .. hi.
 * It is the callers responsibility to call Safefree on the result. */
signed char* range_moebius(UV lo, UV hi)
{
  signed char* mu;
  UV p, i, sqrtn = isqrt(hi), count = hi-lo+1;

  /* Kuznetsov indicates that the Deléglise & Rivat (1996) method can be
   * modified to work on logs, which allows us to operate with no
   * intermediate memory at all.  Same time as the D&R method, less memory. */
  unsigned char logp;
  UV nextlog, nextlogi;

  if (hi < lo) croak("range_mobius error hi %"UVuf" < lo %"UVuf"\n", hi, lo);

  Newz(0, mu, count, signed char);
  if (sqrtn*sqrtn != hi && sqrtn < (UVCONST(1)<<(BITS_PER_WORD/2))-1) sqrtn++;

  {
    PRIME_ITERATOR(iter);

    logp = 1; nextlog = 3; /* 2+1 */
    for (p = 2; p <= sqrtn; p = prime_iterator_next(&iter)) {
      UV p2 = p*p;
      if (p > nextlog) {
        logp += 2;   /* logp is 1 | ceil(log(p)/log(2)) */
        nextlog = ((nextlog-1)*4)+1;
      }
      for (i = P_GT_LO(p, p, lo); i >= lo && i <= hi; i += p)
        mu[i-lo] += logp;
      for (i = P_GT_LO(p2, p2, lo); i >= lo && i <= hi; i += p2)
        mu[i-lo] = 0x80;
    }
    prime_iterator_destroy(&iter);
  }

  logp = log2_ui(lo);
  nextlogi = (UVCONST(2) << logp) - lo;
  for (i = 0; i < count; i++) {
    unsigned char a = mu[i];
    if (i >= nextlogi) nextlogi = (UVCONST(2) << ++logp) - lo;
    if (a & 0x80)       { a = 0; }
    else if (a >= logp) { a =  1 - 2*(a&1); }
    else                { a = -1 + 2*(a&1); }
    mu[i] = a;
  }
  if (lo == 0)  mu[0] = 0;

  return mu;
}

/******************************************************************************/

static short* mertens_array(UV hi)
{
  signed char* mu;
  short* M;
  UV i;

  /* We could blend this with range_moebius but it seems not worth it. */
  mu = range_moebius(0, hi);
  New(0, M, hi+1, short);
  M[0] = 0;
  for (i = 1; i <= hi; i++)
    M[i] = M[i-1] + mu[i];
  Safefree(mu);

  return M;
}

typedef struct {
  UV n;
  IV sum;
} mertens_value_t;
static void _insert_mert_hash(mertens_value_t *H, UV hsize, UV n, IV sum) {
  UV idx = n % hsize;
  H[idx].n = n;
  H[idx].sum = sum;
}
static int _get_mert_hash(mertens_value_t *H, UV hsize, UV n, IV *sum) {
  UV idx = n % hsize;
  if (H[idx].n == n) {
    *sum = H[idx].sum;
    return 1;
  }
  return 0;
}

typedef struct {
  UV               msize;
  const short     *M;       /* 16 bits is enough range 32-bit M => 64-bit n */
  UV               hsize;
  mertens_value_t *H;       /* cache of calculated values */
} mertens_t;

void* hmertens_create(UV n)
{
  mertens_t *HM;
  UV j = 2 * icbrt(n);   /* Biased toward a lot of memory */
  UV maxmu = 1 * j * j;
  UV hsize = 100 + 8*j;

  /* At large sizes, start clamping memory use. */
  if (maxmu > 100000000UL) {
    /* Exponential decay, reduce by factor of 1 to 8 */
    float rfactor = 1.0 + 7.0 * (1.0 - exp(-(float)maxmu/8000000000.0));
    maxmu /= rfactor;
    hsize = hsize * 16;  /* Increase the result cache size */
  }

  hsize = next_prime_ui(hsize-1);  /* Make hash size a prime */

#if BITS_PER_WORD == 64
  /* A 16-bit signed short will overflow at maxmu > 7613644883 */
  if (maxmu > UVCONST(7613644883))  maxmu = UVCONST(7613644883);
#endif

  New(0, HM, 1, mertens_t);
  HM->msize = maxmu;
  HM->M     = mertens_array(maxmu);
  HM->hsize = hsize;
  Newz(0, HM->H, hsize, mertens_value_t);

  return HM;
}

UV hmertens_value(void* ctx, UV n)
{
  mertens_t *HM = ctx;

  UV msize = HM->msize,  hsize = HM->hsize;
  const short* M = HM->M;
  mertens_value_t *H = HM->H;

  UV s, k, ns, nk, nk1, mk, mnk;
  IV sum;

  if (n <= msize)
    return M[n];

  if (_get_mert_hash(H, hsize, n, &sum))
    return sum;

  s = isqrt(n);
  ns = n / (s+1);
  sum = 1;

#if 0
  for (k = 2; k <= ns; k++)
    sum -= _rmertens(n/k, msize, M, H, hsize);
  for (k = 1; k <= s; k++)
    sum -= M[k] * (n/k - n/(k+1));
#else
  /* Take the above: merge the loops and iterate the divides. */
  if (s != ns && s != ns+1) croak("mertens  s / ns");
  nk  = n;
  nk1 = n/2;
  sum -= (nk - nk1);
  for (k = 2; k <= ns; k++) {
    nk = nk1;
    nk1 = n/(k+1);
    mnk = (nk <= msize)  ?  M[nk]  :  hmertens_value(ctx, nk);
    mk  = (k  <= msize)  ?  M[k]   :  hmertens_value(ctx, k);
    sum -= mnk + mk * (nk-nk1);
  }
  if (s > ns)
    sum -= hmertens_value(ctx, s) * (n/s - n/(s+1));
#endif

  _insert_mert_hash(H, hsize, n, sum);
  return sum;
}

void  hmertens_destroy(void* ctx)
{
  mertens_t *HM = ctx;
  Safefree(HM->H);
  Safefree(HM->M);
  Safefree(ctx);
}

IV mertens_ui(UV n) {
  void* mctx;
  IV sum;

  if (n <= 512) {
    signed char* mu = range_moebius(0,n);
    UV j;
    for (j = 0, sum = 0; j <= n; j++)
      sum += mu[j];
    Safefree(mu);
    return sum;
  }

  mctx = hmertens_create(n);
  sum = hmertens_value(mctx, n);
  hmertens_destroy(mctx);

  return sum;
}

/******************************************************************************/

static const signed char _small_liouville[16] = {-1,1,-1,-1,1,-1,1,-1,-1,1,1,-1,-1,-1,1,1};
IV sumliouville_ui(UV n) {
  void* mctx;
  UV j, k, nk, sqrtn;
  IV sum;

  if (n <= 15) {
    for (sum = 0, j = 1; j <= n; j++)
      sum += _small_liouville[j];
    return sum;
  }
  mctx = hmertens_create(n);

  sqrtn = isqrt(n);
  sum = hmertens_value(mctx, n);
  for (k = 2; k <= sqrtn; k++) {
    nk = n / (k*k);
    if (nk == 1) break;
    sum += hmertens_value(mctx, nk);
  }
  sum += (sqrtn + 1 - k);  /* all k where n/(k*k) == 1 */
  /* TODO: find method to get exact number of n/(k*k)==1 .. 4.  Halves k */
  /*       Ends up with method like Lehmer's g. */
  hmertens_destroy(mctx);

  return sum;
}

/******************************************************************************/

/* This method that doesn't get the divisors of n is just as good under 1e6.
 * Above that it is much slower.
 */

IV hclassno_ui(UV n) {
  UV h, b, b2, i, lim;
  int square;

  if (n == 0) return -1;
  if ((n % 4) == 1 || (n % 4) == 2) return 0;

  h = 0;
  square = 0;
  b = n & 1;
  b2 = (n+1) >> 2;

  if (b == 0) {
    lim = isqrt(b2);
    if (lim*lim == b2) {
      square = 1;
      lim--;
    }
    for (i = 1; i <= lim; i++)
      if ((b2 % i) == 0)
        h++;
    b = 2;
    b2 = (n+4) >> 2;
  }
  while (b2*3 < n) {
    if ((b2 % b) == 0)  h++;
    lim = isqrt(b2);
    if (lim*lim == b2) {
      h++;
      lim--;
    }
    for (i = b+1; i <= lim; i++)
      if ((b2 % i) == 0)
        h += 2;
    b += 2;
    b2 = (n+b*b) >> 2;
  }
  return 2 * ((b2*3 == n)  ?  2*(3*h+1)  :  (square  ?  3*(2*h+1)  :  6*h));
}
