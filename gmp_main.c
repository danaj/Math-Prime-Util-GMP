#include <string.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
#include <gmp.h>
#include "ptypes.h"

#include "gmp_main.h"
#include "primality.h"
#include "prime_iterator.h"
#include "ecpp.h"
#include "factor.h"
#include "real.h"
#include "random_prime.h"

#define FUNC_gcd_ui 1
#define FUNC_mpz_logn 1
#define FUNC_isqrt 1
#include "utility.h"

static mpz_t _bgcd;
static mpz_t _bgcd2;
static mpz_t _bgcd3;
#define BGCD_PRIMES       168
#define BGCD_LASTPRIME    997
#define BGCD_NEXTPRIME   1009
#define BGCD2_PRIMES     1229
#define BGCD2_NEXTPRIME 10007
#define BGCD3_PRIMES     4203
#define BGCD3_NEXTPRIME 40009

#define NSMALLPRIMES 168
static const unsigned short sprimes[NSMALLPRIMES] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997};

#define TSTAVAL(arr, val)   (arr[(val) >> 6] & (1U << (((val)>>1) & 0x1F)))
#define SETAVAL(arr, val)   arr[(val) >> 6] |= 1U << (((val)>>1) & 0x1F)


void _GMP_init(void)
{
  /* For real use of randomness we need to be seeded properly.
   * This gives us a start until someone calls seed_csprng().
   * We could try to improve this duct-tape in various ways. */
  unsigned long seed = time(NULL);
  init_randstate(seed);
  prime_iterator_global_startup();
  mpz_init(_bgcd);
  _GMP_pn_primorial(_bgcd, BGCD_PRIMES);   /* mpz_primorial_ui(_bgcd, 1000) */
  mpz_init_set_ui(_bgcd2, 0);
  mpz_init_set_ui(_bgcd3, 0);
  _init_factor();
}

void _GMP_destroy(void)
{
  _GMP_memfree();
  prime_iterator_global_shutdown();
  clear_randstate();
  mpz_clear(_bgcd);
  mpz_clear(_bgcd2);
  mpz_clear(_bgcd3);
}

void _GMP_memfree(void)
{
  free_float_constants();
  destroy_ecpp_gcds();
  free_borwein_zeta();
  free_bernoulli();
}

static const unsigned char next_wheel[30] =
  {1,7,7,7,7,7,7,11,11,11,11,13,13,17,17,17,17,19,19,23,23,23,23,29,29,29,29,29,29,1};
static const unsigned char prev_wheel[30] =
  {29,29,1,1,1,1,1,1,7,7,7,7,11,11,13,13,13,13,17,17,19,19,19,19,23,23,23,23,23,23};
static const unsigned char wheel_advance[30] =
  {1,6,5,4,3,2,1,4,3,2,1,2,1,4,3,2,1,2,1,4,3,2,1,6,5,4,3,2,1,2};
static const unsigned char wheel_retreat[30] =
  {1,2,1,2,3,4,5,6,1,2,3,4,1,2,1,2,3,4,1,2,1,2,3,4,1,2,3,4,5,6};



/*****************************************************************************/

static int is_tiny_prime(uint32_t n) {
  uint32_t f, limit;
  if (n < 11) {
    if (n == 2 || n == 3 || n == 5 || n == 7)      return 2;
    else                                           return 0;
  }
  if (!(n%2) || !(n%3) || !(n%5) || !(n%7))        return 0;
  if (n <  121) /* 11*11 */                        return 2;
  if (!(n%11) || !(n%13) || !(n%17) || !(n%19) ||
       !(n%23) || !(n%29) || !(n%31) || !(n%37) ||
       !(n%41) || !(n%43) || !(n%47) || !(n%53))   return 0;
  if (n < 3481) /* 59*59 */                        return 2;

  f = 59;
  limit = (uint32_t) (sqrt((double)n));
  while (f <= limit) {
    if ( !(n%f) || !(n%(f+2)) || !(n%(f+8)) || !(n%(f+12)) ) return 0;
    f += 14;
    if ( !(n%f) || !(n%(f+4)) || !(n%(f+6)) || !(n%(f+10)) ) return 0;
    f += 16;
  }
  return 2;
}

int primality_pretest(mpz_t n)
{
  if (mpz_cmp_ui(n, 100000) < 0)
    return is_tiny_prime((uint32_t)mpz_get_ui(n));

  /* Check for tiny divisors */
  if (mpz_even_p(n)) return 0;
  if (sizeof(unsigned long) < 8) {
    if (mpz_gcd_ui(NULL, n, 3234846615UL) != 1) return 0;           /*  3-29 */
  } else {
    if (mpz_gcd_ui(NULL, n, 4127218095UL*3948078067UL)!=1) return 0;/*  3-53 */
    if (mpz_gcd_ui(NULL, n, 4269855901UL*1673450759UL)!=1) return 0;/* 59-101 */
  }

  {
    UV log2n = mpz_sizeinbase(n,2);
    mpz_t t;
    mpz_init(t);

    /* Do a GCD with all primes < 1009 */
    mpz_gcd(t, n, _bgcd);
    if (mpz_cmp_ui(t, 1))
      { mpz_clear(t); return 0; }

    /* No divisors under 1009 */
    if (mpz_cmp_ui(n, BGCD_NEXTPRIME*BGCD_NEXTPRIME) < 0)
      { mpz_clear(t); return 2; }

    /* If we're reasonably large, do a gcd with more primes */
    if (log2n > 700) {
      if (mpz_sgn(_bgcd3) == 0) {
        _GMP_pn_primorial(_bgcd3, BGCD3_PRIMES);
        mpz_divexact(_bgcd3, _bgcd3, _bgcd);
      }
      mpz_gcd(t, n, _bgcd3);
      if (mpz_cmp_ui(t, 1))
        { mpz_clear(t); return 0; }
    } else if (log2n > 300) {
      if (mpz_sgn(_bgcd2) == 0) {
        _GMP_pn_primorial(_bgcd2, BGCD2_PRIMES);
        mpz_divexact(_bgcd2, _bgcd2, _bgcd);
      }
      mpz_gcd(t, n, _bgcd2);
      if (mpz_cmp_ui(t, 1))
        { mpz_clear(t); return 0; }
    }
    mpz_clear(t);
    /* Do more trial division if we think we should.
     * According to Menezes (section 4.45) as well as Park (ISPEC 2005),
     * we want to select a trial limit B such that B = E/D where E is the
     * time for our primality test (one M-R test) and D is the time for
     * one trial division.  Example times on my machine came out to
     *   log2n = 840375, E= 6514005000 uS, D=1.45 uS, E/D = 0.006 * log2n
     *   log2n = 465618, E= 1815000000 uS, D=1.05 uS, E/D = 0.008 * log2n
     *   log2n = 199353, E=  287282000 uS, D=0.70 uS, E/D = 0.01  * log2n
     *   log2n =  99678, E=   56956000 uS, D=0.55 uS, E/D = 0.01  * log2n
     *   log2n =  33412, E=    4289000 uS, D=0.30 uS, E/D = 0.013 * log2n
     *   log2n =  13484, E=     470000 uS, D=0.21 uS, E/D = 0.012 * log2n
     * Our trial division could also be further improved for large inputs.
     */
    if (log2n > 16000) {
      double dB = (double)log2n * (double)log2n * 0.005;
      if (BITS_PER_WORD == 32 && dB > 4200000000.0) dB = 4200000000.0;
      if (_GMP_trial_factor(n, BGCD3_NEXTPRIME, (UV)dB))  return 0;
    } else if (log2n > 4000) {
      if (_GMP_trial_factor(n, BGCD3_NEXTPRIME, 80*log2n))  return 0;
    } else if (log2n > 1600) {
      if (_GMP_trial_factor(n, BGCD3_NEXTPRIME, 30*log2n))  return 0;
    }
  }
  return 1;
}

/* Primality using purely trial division.
 * We basically never want to use this for anything practical,
 * but it is good for testing and comparison.
 */
int is_trial_prime(mpz_t n)
{
  mpz_t sqrtn;
  uint32_t res, p, lim;
  int _verbose = get_verbose_level();

  if (_GMP_trial_factor(n, 2, 256))
    return 0;
  if (mpz_cmp_ui(n, 257*257) < 0)
    return 1;

  mpz_init(sqrtn);
  mpz_sqrt(sqrtn, n);
  if (_verbose >= 2) gmp_printf("    trial division to %Zd\n", sqrtn);
  /* Next prime is 257 */
  {
    PRIME_ITERATOR(iter);

    prime_iterator_setprime(&iter, 256);
    p = prime_iterator_next(&iter);
    lim = (mpz_cmp_ui(sqrtn, 4294967279) >= 0) ? 4294967279 : mpz_get_ui(sqrtn);
    if (_verbose >= 3) gmp_printf("    ... using mpz_ui and prime iterator to %lu\n", lim);
    while (p <= lim) {
      if (mpz_divisible_ui_p(n, p))
        break;
      p = prime_iterator_next(&iter);
    }
    res = (p > lim);
    prime_iterator_destroy(&iter);

    /* After benchmarking, sieve_primes is faster than using a 29-rough
     * iterator (basically next_prime without the final primality test).
     * A remainder tree would be faster, or even vecmul + gcd. */

    if (res == 1 && mpz_cmp_ui(sqrtn, p) >= 0) {   /* We need to keep going */
      mpz_t zp, nlo, nhi;
      UV i, nprimes, *list;
      mpz_init_set_ui(nlo, p);
      mpz_init(nhi);
      mpz_init(zp);
      if (_verbose >= 3) gmp_printf("    ... using mpz    and sieve_primes from %Zd to %Zd\n", nlo, sqrtn);
      while (res == 1 && mpz_cmp(nlo, sqrtn) <= 0) {
        mpz_add_ui(nhi, nlo, 10000000);
        if (mpz_cmp(nhi, sqrtn) > 0) mpz_set(nhi, sqrtn);
        list = sieve_primes(nlo, nhi, 0, &nprimes);
        if (list != 0) {
          for (i = 0; i < nprimes; i++) {
            mpz_add_ui(zp, nlo, list[i]);
            if (mpz_divisible_p(n, zp))
              break;
          }
          res = (i >= nprimes);
          Safefree(list);
        }
        mpz_add_ui(nlo, nhi, mpz_odd_p(nhi) ? 2 : 1);
      }
      mpz_clear(zp);
      mpz_clear(nhi);
      mpz_clear(nlo);
    }
  }
  mpz_clear(sqrtn);
  return res;
}


/*****************************************************************************/

/* Controls how many numbers to sieve.  Little time impact. */
#define NPS_MERIT  30.0
/* Controls how many primes to use.  Big time impact. */
static UV _nps_depth(UV log2n, UV log2log2n) {
  double d;
  if (log2n < 100) return 1000;
  d = 0.75 * (double)log2n * (double)(log2n >> 5) * (double)log2log2n;
  /* Make sure we don't try to sieve too far. */
  if (d >= (double)(UV_MAX>>1)) return (UV_MAX>>1);
  return (UV) d;
}

static void next_prime_with_sieve(mpz_t n) {
  UV i, log2n, log2log2n, width, depth;
  uint32_t* comp;
  mpz_t t, base;
  log2n = mpz_sizeinbase(n, 2);
  for (log2log2n = 1, i = log2n; i >>= 1; ) log2log2n++;
  width = (UV) (NPS_MERIT/1.4427 * (double)log2n + 0.5);
  depth = _nps_depth(log2n, log2log2n);

  if (width & 1) width++;                     /* Make width even */
  mpz_add_ui(n, n, mpz_even_p(n) ? 1 : 2);    /* Set n to next odd */
  mpz_init(t);  mpz_init(base);
  while (1) {
    mpz_set(base, n);
    comp = partial_sieve(base, width, depth); /* sieve range to depth */
    for (i = 1; i <= width; i += 2) {
      if (!TSTAVAL(comp, i)) {
        mpz_add_ui(t, base, i);               /* We found a candidate */
        if (_GMP_BPSW(t)) {
          mpz_set(n, t);
          mpz_clear(t);  mpz_clear(base);
          Safefree(comp);
          return;
        }
      }
    }
    Safefree(comp);       /* A huge gap found, so sieve another range */
    mpz_add_ui(n, n, width);
  }
}

static void prev_prime_with_sieve(mpz_t n) {
  UV i, j, log2n, log2log2n, width, depth;
  uint32_t* comp;
  mpz_t t, base;
  log2n = mpz_sizeinbase(n, 2);
  for (log2log2n = 1, i = log2n; i >>= 1; ) log2log2n++;
  width = (UV) (NPS_MERIT/1.4427 * (double)log2n + 0.5);
  depth = _nps_depth(log2n, log2log2n);

  mpz_sub_ui(n, n, mpz_even_p(n) ? 1 : 2);       /* Set n to prev odd */
  width = 64 * ((width+63)/64);                /* Round up to next 64 */
  mpz_init(t);  mpz_init(base);
  while (1) {
    mpz_sub_ui(base, n, width-2);
    /* gmp_printf("sieve from %Zd to %Zd width %lu\n", base, n, width); */
    comp = partial_sieve(base, width, depth); /* sieve range to depth */
    /* if (mpz_odd_p(base)) croak("base off after partial");
       if (width & 1) croak("width is odd after partial"); */
    for (j = 1; j < width; j += 2) {
      i = width - j;
      if (!TSTAVAL(comp, i)) {
        mpz_add_ui(t, base, i);               /* We found a candidate */
        if (_GMP_BPSW(t)) {
          mpz_set(n, t);
          mpz_clear(t);  mpz_clear(base);
          Safefree(comp);
          return;
        }
      }
    }
    Safefree(comp);       /* A huge gap found, so sieve another range */
    mpz_sub_ui(n, n, width);
  }
}

/* Modifies argument */
void _GMP_next_prime(mpz_t n)
{
  if (mpz_cmp_ui(n, 29) < 0) { /* small inputs */

    UV m = mpz_get_ui(n);
    m = (m < 2) ? 2 : (m < 3) ? 3 : (m < 5) ? 5 : next_wheel[m];
    mpz_set_ui(n, m);

  } else if (mpz_sizeinbase(n,2) > 120) {

    next_prime_with_sieve(n);

  } else {

    UV m23 = mpz_fdiv_ui(n, 223092870UL);  /* 2*3*5*7*11*13*17*19*23 */
    UV m = m23 % 30;
    do {
      uint32_t skip = wheel_advance[m];
      mpz_add_ui(n, n, skip);
      m23 += skip;
      m = next_wheel[m];
    } while ( !(m23% 7) || !(m23%11) || !(m23%13) || !(m23%17) ||
              !(m23%19) || !(m23%23) || !_GMP_is_prob_prime(n) );

  }
}

/* Modifies argument */
void _GMP_prev_prime(mpz_t n)
{
  UV m, m23;
  if (mpz_cmp_ui(n, 29) <= 0) { /* small inputs */

    m = mpz_get_ui(n);
    m = (m < 3) ? 0 : (m < 4) ? 2 : (m < 6) ? 3 : (m < 8) ? 5 : prev_wheel[m];
    mpz_set_ui(n, m);

  } else if (mpz_sizeinbase(n,2) > 200) {

    prev_prime_with_sieve(n);

  } else {

    m23 = mpz_fdiv_ui(n, 223092870UL);  /* 2*3*5*7*11*13*17*19*23 */
    m = m23 % 30;
    m23 += 223092870UL;  /* No need to re-mod inside the loop */
    do {
      uint32_t skip = wheel_retreat[m];
      mpz_sub_ui(n, n, skip);
      m23 -= skip;
      m = prev_wheel[m];
    } while ( !(m23% 7) || !(m23%11) || !(m23%13) || !(m23%17) ||
              !(m23%19) || !(m23%23) || !_GMP_is_prob_prime(n) );

  }
}

void surround_primes(mpz_t n, UV* prev, UV* next, UV skip_width) {
  UV i, j, log2n, log2log2n, width, depth, fprev, fnext, search_merits;
  uint32_t* comp;
  mpz_t t, base;
  int neven, found;

  log2n = mpz_sizeinbase(n, 2);
  for (log2log2n = 1, i = log2n; i >>= 1; ) log2log2n++;

  if (log2n < 8*sizeof(unsigned long)) {
    mpz_init(t);
    mpz_set(t, n);
    _GMP_prev_prime(t);
    mpz_sub(t, n, t);
    *prev = mpz_get_ui(t);
    mpz_set(t, n);
    _GMP_next_prime(t);
    mpz_sub(t, t, n);
    *next = mpz_get_ui(t);
    mpz_clear(t);
    return;
  }

  mpz_init(t);
  mpz_init(base);
  fprev = fnext = 0;
  neven = mpz_even_p(n);
  j = 1 + !neven;         /* distance from n we're looking. */

  for (found = 0, search_merits = 20; !found; search_merits *= 2) {
    double logn = mpz_logn(n);

    if (BITS_PER_WORD == 32 && log2n >  16600)
      depth = UVCONST(   2500000000);
    else if (BITS_PER_WORD == 64 && log2n > 203600)
      depth = UVCONST(6000000000000);
    else if (log2n > 900)
      depth = (UV) ((-.05L+(log2n/8000.0L)) * logn * logn * log(logn));
    else
      depth = _nps_depth(log2n, log2log2n);

    width = (UV) (search_merits * logn + 0.5);
    width = 64 * ((width+63)/64);    /* Round up to next 64 */
    if (neven) width++;              /* base will always be odd */
    mpz_sub_ui(base, n, width);
    /* printf("merits %lu  width %lu  depth %lu  skip_width %lu\n", search_merits, width, depth, skip_width); */

    /* gmp_printf("partial sieve width %lu  depth %lu\n", 2*width+1, depth); */
    comp = partial_sieve(base, 2*width+1, depth);

    for (; j < width; j += 2) {
      if (!fprev) {
        if (!TSTAVAL(comp, width+1-j)) {
          mpz_sub_ui(t, n, j);
          if ( (skip_width == 0) ? _GMP_BPSW(t) : miller_rabin_ui(t,2) ) {
            fprev = j;
            if (fnext || (skip_width != 0 && j <= skip_width))
              break;
          }
        }
      }
      if (!fnext) {
        if (!TSTAVAL(comp, width+1+j)) {
          mpz_add_ui(t, n, j);
          if ( (skip_width == 0) ? _GMP_BPSW(t) : miller_rabin_ui(t,2) ) {
            fnext = j;
            if (fprev || (skip_width != 0 && j <= skip_width))
              break;
          }
        }
      }
    }

    Safefree(comp);
    if ( (fprev && fnext) ||
         (skip_width != 0 && j <= skip_width && (fprev || fnext)) )
      found = 1;
  }

  mpz_clear(base);
  mpz_clear(t);

  *prev = fprev;
  *next = fnext;
}

/*****************************************************************************/


#define LAST_TRIPLE_PROD \
  ((ULONG_MAX <= 4294967295UL) ? UVCONST(1619) : UVCONST(2642231))
#define LAST_DOUBLE_PROD \
  ((ULONG_MAX <= 4294967295UL) ? UVCONST(65521) : UVCONST(4294967291))
void _GMP_pn_primorial(mpz_t prim, UV n)
{
  UV i = 0, al = 0, p = 2;
  mpz_t* A;

  if (n <= 4) {                 /* tiny input */

    p = (n == 0) ? 1 : (n == 1) ? 2 : (n == 2) ? 6 : (n == 3) ? 30 : 210;
    mpz_set_ui(prim, p);

  } else if (n < 200) {         /* simple linear multiply */

    PRIME_ITERATOR(iter);
    mpz_set_ui(prim, 1);
    while (n-- > 0) {
      if (n > 0) { p *= prime_iterator_next(&iter); n--; }
      mpz_mul_ui(prim, prim, p);
      p = prime_iterator_next(&iter);
    }
    prime_iterator_destroy(&iter);

  } else {                      /* tree mult array of products of 8 UVs */

    PRIME_ITERATOR(iter);
    New(0, A, n, mpz_t);
    while (n-- > 0) {
      if (p <= LAST_TRIPLE_PROD && n > 0)
        { p *= prime_iterator_next(&iter); n--; }
      if (p <= LAST_DOUBLE_PROD && n > 0)
        { p *= prime_iterator_next(&iter); n--; }
      /* each array entry holds the product of 8 UVs */
      if ((i & 7) == 0) mpz_init_set_ui( A[al++], p );
      else              mpz_mul_ui(A[al-1],A[al-1], p );
      i++;
      p = prime_iterator_next(&iter);
    }
    mpz_product(A, 0, al-1);
    mpz_set(prim, A[0]);
    for (i = 0; i < al; i++)  mpz_clear(A[i]);
    Safefree(A);
    prime_iterator_destroy(&iter);

  }
}
void _GMP_primorial(mpz_t prim, UV n)
{
#if (__GNU_MP_VERSION > 5) || (__GNU_MP_VERSION == 5 && __GNU_MP_VERSION_MINOR >= 1)
  mpz_primorial_ui(prim, n);
#else
  if (n <= 4) {

    UV p = (n == 0) ? 1 : (n == 1) ? 1 : (n == 2) ? 2 : (n == 3) ? 6 : 6;
    mpz_set_ui(prim, p);

  } else {

    mpz_t *A;
    UV nprimes, i, al;
    UV *primes = sieve_to_n(n, &nprimes);

    /* Multiply native pairs until we overflow the native type */
    while (nprimes > 1 && ULONG_MAX/primes[0] > primes[nprimes-1]) {
      i = 0;
      while (nprimes > i+1 && ULONG_MAX/primes[i] > primes[nprimes-1])
        primes[i++] *= primes[--nprimes];
    }

    if (nprimes <= 8) {
      /* Just multiply if there are only a few native values left */
      mpz_set_ui(prim, primes[0]);
      for (i = 1; i < nprimes; i++)
        mpz_mul_ui(prim, prim, primes[i]);
    } else {
      /* Create n/4 4-way products, then use product tree */
      New(0, A, nprimes/4 + 1, mpz_t);
      for (i = 0, al = 0; i < nprimes; al++) {
        mpz_init_set_ui(A[al], primes[i++]);
        if (i < nprimes) mpz_mul_ui(A[al], A[al], primes[i++]);
        if (i < nprimes) mpz_mul_ui(A[al], A[al], primes[i++]);
        if (i < nprimes) mpz_mul_ui(A[al], A[al], primes[i++]);
      }
      mpz_product(A, 0, al-1);
      mpz_set(prim, A[0]);
      for (i = 0; i < al; i++)  mpz_clear(A[i]);
      Safefree(A);
    }
    Safefree(primes);

  }
#endif
}

/*****************************************************************************/

void stirling(mpz_t r, unsigned long n, unsigned long m, UV type)
{
  mpz_t t, t2;
  unsigned long j;
  if (type < 1 || type > 3) croak("stirling type must be 1, 2, or 3");
  if (n == m) {
    mpz_set_ui(r, 1);
  } else if (n == 0 || m == 0 || m > n) {
    mpz_set_ui(r,0);
  } else if (m == 1) {
    switch (type) {
      case 1:  mpz_fac_ui(r, n-1);  if (!(n&1)) mpz_neg(r, r); break;
      case 2:  mpz_set_ui(r, 1); break;
      case 3:
      default: mpz_fac_ui(r, n); break;
    }
  } else {
    mpz_init(t);  mpz_init(t2);
    mpz_set_ui(r,0);
    if (type == 3) { /* Lah: binomial(n k) * binomial(n-1 k-1) * (n-k)!*/
      mpz_bin_uiui(t, n, m);
      mpz_bin_uiui(t2, n-1, m-1);
      mpz_mul(r, t, t2);
      mpz_fac_ui(t2, n-m);
      mpz_mul(r, r, t2);
    } else if (type == 2) {
      mpz_t binom;
      mpz_init_set_ui(binom, m);
      mpz_ui_pow_ui(r, m, n);
      /* Use symmetry to halve the number of loops */
      for (j = 1; j <= ((m-1)>>1); j++) {
        mpz_ui_pow_ui(t, j, n);
        mpz_ui_pow_ui(t2, m-j, n);
        if (m&1) mpz_sub(t, t2, t);
        else     mpz_add(t, t2, t);
        mpz_mul(t, t, binom);
        if (j&1) mpz_sub(r, r, t);
        else     mpz_add(r, r, t);
        mpz_mul_ui(binom, binom, m-j);
        mpz_divexact_ui(binom, binom, j+1);
      }
      if (!(m&1)) {
        mpz_ui_pow_ui(t, j, n);
        mpz_mul(t, t, binom);
        if (j&1) mpz_sub(r, r, t);
        else     mpz_add(r, r, t);
      }
      mpz_clear(binom);
      mpz_fac_ui(t, m);
      mpz_divexact(r, r, t);
    } else {
      mpz_bin_uiui(t,  n-1+1, n-m+1);
      mpz_bin_uiui(t2, n-m+n, n-m-1);
      mpz_mul(t2, t2, t);
      for (j = 1; j <= n-m; j++) {
        stirling(t, n-m+j, j, 2);
        mpz_mul(t, t, t2);
        if (j & 1)      mpz_sub(r, r, t);
        else            mpz_add(r, r, t);
        mpz_mul_ui(t2, t2, n+j);
        mpz_divexact_ui(t2, t2, n-m+j+1);
        mpz_mul_ui(t2, t2, n-m-j);
        mpz_divexact_ui(t2, t2, n+j+1);
      }
    }
    mpz_clear(t2);  mpz_clear(t);
  }
}

/* Goetgheluck method.  Also thanks to Peter Luschny. */
void binomial(mpz_t r, UV n, UV k)
{
  UV fi, nk, sqrtn, piN, prime, i, j;
  UV* primes;
  mpz_t* mprimes;

  if (k > n)            { mpz_set_ui(r, 0); return; }
  if (k == 0 || k == n) { mpz_set_ui(r, 1); return; }

  if (k > n/2)  k = n-k;

  sqrtn = (UV) (sqrt((double)n));
  fi = 0;
  nk = n-k;
  primes = sieve_to_n(n, &piN);

#define PUSHP(p) \
 do { \
   if ((j++ % 8) == 0) mpz_init_set_ui(mprimes[fi++], p); \
   else                mpz_mul_ui(mprimes[fi-1], mprimes[fi-1], p); \
 } while (0)

  New(0, mprimes, (piN+7)/8, mpz_t);

  for (i = 0, j = 0; i < piN; i++) {
    prime = primes[i];
    if (prime > nk) {
      PUSHP(prime);
    } else if (prime > n/2) {
      /* nothing */
    } else if (prime > sqrtn) {
      if (n % prime < k % prime)
        PUSHP(prime);
    } else {
      UV N = n, K = k, p = 1, s = 0;
      while (N > 0) {
        s = (N % prime) < (K % prime + s) ? 1 : 0;
        if (s == 1)  p *= prime;
        N /= prime;
        K /= prime;
      }
      if (p > 1)
        PUSHP(p);
    }
  }
  Safefree(primes);
  mpz_product(mprimes, 0, fi-1);
  mpz_set(r, mprimes[0]);
  for (i = 0; i < fi; i++)
    mpz_clear(mprimes[i]);
  Safefree(mprimes);
}

void multifactorial(mpz_t r, unsigned long n, unsigned long k)
{
  if (k == 0) {  mpz_set_ui(r, 1); return;  }
  if (k == 1) {  mpz_fac_ui(r, n); return;  }
#if (__GNU_MP_VERSION > 5) || (__GNU_MP_VERSION == 5 && __GNU_MP_VERSION_MINOR >= 1)
  mpz_mfac_uiui(r, n, k);
#else
  /* Naive code.  Slow with large n and small k.
   * See Luschny or mpz_mfac_uiui for better. */
  mpz_set_ui(r, (n > 1) ? n : 1);
  while (n > k) {
    n -= k;
    mpz_mul_ui(r, r, n);
  }
#endif
}

void factorial_sum(mpz_t r, unsigned long n)
{
  mpz_t t;
  unsigned long k;
  if (n == 0) { mpz_set_ui(r,0); return; }

  mpz_set_ui(r,1);
  mpz_init_set_ui(t,1);
  for (k = 1; k < n; k++) {
    mpz_mul_ui(t, t, k);
    mpz_add(r, r, t);
  }
  mpz_clear(t);
}

void subfactorial(mpz_t r, unsigned long n)
{
  unsigned long k;
  if (n == 0) { mpz_set_ui(r,1); return; }
  if (n == 1) { mpz_set_ui(r,0); return; }
  /* We could loop using Pochhammer, but that's much slower. */
  mpz_set_ui(r,0);
  for (k = 2; k <= n; k++) {
    mpz_mul_ui(r, r, k);
    if (k & 1) mpz_sub_ui(r, r, 1);
    else       mpz_add_ui(r, r, 1);
  }
}

void falling_factorial(mpz_t r, unsigned long x, unsigned long n)
{
  mpz_t t;
  if (n == 0) { mpz_set_ui(r,1); return; }
  mpz_init(t);
  mpz_bin_uiui(t, x, n);
  mpz_fac_ui(r, n);
  mpz_mul(r, r, t);
  mpz_clear(t);
}
void rising_factorial(mpz_t r, unsigned long x, unsigned long n) {
  falling_factorial(r, x+n-1, n);
}


void factorialmod(mpz_t r, UV N, mpz_t m)
{
  int m_is_prime;
  mpz_t t, t2;
  UV D = N, i, p;

  if (mpz_cmp_ui(m,N) <= 0 || mpz_cmp_ui(m,1) <= 0) {
    mpz_set_ui(r,0);
    return;
  }

  m_is_prime = _GMP_is_prime(m);
  mpz_init(t);
  mpz_tdiv_q_2exp(t, m, 1);
  if (mpz_cmp_ui(t, N) < 0 && m_is_prime)
    D = mpz_get_ui(m) - N - 1;

  if (D < 2 && N > D) {
    if (D == 0) mpz_sub_ui(r, m, 1);
    else        mpz_set_ui(r, 1);
    mpz_clear(t);
    return;
  }

  if (D > 500 && !m_is_prime) {
    mpz_t *factors;
    int j, nfactors, *exponents, reszero;
    nfactors = factor(m, &factors, &exponents);
    /* Find max factor */
    mpz_set_ui(t, 0);
    for (j = 0; j < nfactors; j++) {
      if (exponents[j] > 1)
        mpz_mul_ui(factors[j], factors[j], exponents[j]);
      if (mpz_cmp(factors[j], t) > 0)
        mpz_set(t, factors[j]);
    }
    /* for m=p^k * p^k ..., t is max(p*k,p*k,...).  This is >= S(m), where
     * S(m) is the smallest value where m divides S(m)!.  Hence, every
     * n! mod m will be zero at that value or higher.  We could calculate
     * the exact value of S(m), then we would know there are no zero results
     * for the larger case. */
    reszero = (mpz_cmp_ui(t, N) <= 0);
    clear_factors(nfactors, &factors, &exponents);
    if (reszero) { mpz_clear(t); mpz_set_ui(r,0); return; }
  }

  /* Accumulate into t, then mod into r at the end. */
  mpz_set_ui(t,1);

  /* For small D, naive method. */
  if (D <= 1000) {
    for (i = 2; i <= D && mpz_sgn(t); i++) {
      mpz_mul_ui(t, t, i);
      if ((i & 15) == 0) mpz_mod(t, t, m);
    }
  } else {
    UV j, sd = isqrt(D);
    PRIME_ITERATOR(iter);

    mpz_init(t2);
    mpz_set_ui(t,1);
    /* Group into powers of primes */
    for (p = 2, i = 0; p <= D/sd; p = prime_iterator_next(&iter)) {
      UV td = D/p,  e = td;
      do { td /= p; e += td; } while (td > 0);
      mpz_set_ui(t2, p);
      mpz_powm_ui(t2, t2, e, m);
      mpz_mul(t, t, t2);
      if ((i++ & 15) == 0) {
        mpz_mod(t, t, m);
        if (!mpz_sgn(t)) break;
      }
    }
    /* Further group by primes with the same power. */
    for (j = sd-1; j >= 1 && mpz_sgn(t); j--) {
      UV lo = D / (j+1)+1,  hi = D / j;
      MPUassert(p >= lo, "factorialmod prime loop p should be in range");
      /* while (p < lo) p = prime_iterator_next(&iter); */
      for (mpz_set_ui(t2,1), i=0;  p <= hi;  p = prime_iterator_next(&iter)) {
        mpz_mul_ui(t2, t2, p);
        if ((i++ & 15) == 0) mpz_mod(t2, t2, m);
      }
      mpz_powm_ui(t2, t2, j, m);
      mpz_mul(t, t, t2);
      if ((j & 15) == 0) mpz_mod(t, t, m);
    }
    mpz_clear(t2);
    prime_iterator_destroy(&iter);
  }
  mpz_mod(r, t, m);
  mpz_clear(t);

  /* If we used Wilson's theorem, turn the result for D! into N! */
  if (D != N && mpz_sgn(r)) {
    if (!(D&1)) mpz_sub(r, m, r);
    mpz_invert(r, r, m);
  }
}

void partitions(mpz_t npart, UV n)
{
  mpz_t psum, *part;
  UV *pent, i, j, k, d = (UV) sqrt(n+1);

  if (n <= 3) {
    mpz_set_ui(npart, (n == 0) ? 1 : n);
    return;
  }

  New(0, pent, 2*d+2, UV);
  pent[0] = 0;
  pent[1] = 1;
  for (i = 1; i <= d; i++) {
    pent[2*i  ] = ( i   *(3*i+1)) / 2;
    pent[2*i+1] = ((i+1)*(3*i+2)) / 2;
  }
  New(0, part, n+1, mpz_t);
  mpz_init_set_ui(part[0], 1);
  mpz_init(psum);
  for (j = 1; j <= n; j++) {
    mpz_set_ui(psum, 0);
    for (k = 1; pent[k] <= j; k++) {
      if ((k+1) & 2) mpz_add(psum, psum, part[ j - pent[k] ]);
      else           mpz_sub(psum, psum, part[ j - pent[k] ]);
    }
    mpz_init_set(part[j], psum);
  }

  mpz_set(npart, part[n]);

  mpz_clear(psum);
  for (i = 0; i <= n; i++)
    mpz_clear(part[i]);
  Safefree(part);
  Safefree(pent);
}

void faulhaber_sum(mpz_t sum, mpz_t zn, unsigned long p) /*Sum_1^n(k^p)*/
{
  const mpz_t *N, *D;
  mpz_t t, nj, num, den;
  unsigned long j;

  if (mpz_cmp_ui(zn, 1) <= 0) {
    mpz_set_ui(sum, (mpz_sgn(zn) > 0) ? 1 : 0);
    return;
  }

  mpz_init(t);

  /* Use the polynomials directly for tiny powers */
  if (p <= 3) {
    mpz_add_ui(t, zn, 1);
    switch (p) {
      case 0: mpz_set(sum, zn);
              break;
      case 1: mpz_mul(sum, zn, t);
              mpz_divexact_ui(sum, sum, 2);
              break;
      case 2: mpz_mul(sum, zn, t);
              mpz_mul_ui(t, t, 2); mpz_sub_ui(t, t, 1);
              mpz_mul(sum, sum, t);
              mpz_divexact_ui(sum, sum, 6);
              break;
      case 3: mpz_mul(sum, zn, t);
              mpz_divexact_ui(sum, sum, 2);
              mpz_mul(sum, sum, sum);
              break;
    }
    mpz_clear(t);
    return;
  }
  /* If n is small, doing directly can be much faster.  Cheating.... */
  if (mpz_cmp_ui(zn, p) <= 0) {
    unsigned long n = mpz_get_ui(zn);
    mpz_set_ui(sum, 1);
    for (j = 1; j < n; j++) {
      mpz_ui_pow_ui(t, j+1, p);
      mpz_add(sum, sum, t);
    }
    mpz_clear(t);
    return;
  }

  bernvec(&N, &D, p >> 1);
  mpz_init_set_ui(num,0);
  mpz_init_set_ui(den,1);
  mpz_init_set(nj, zn);

  /* Loop backwards so we build up n^(p+1-j) as we go */
  for (j = p; j >= 2; j--) {
    if (!(j&1)) {  /* j is even so B_j is non-zero */
      mpz_bin_uiui(t, p+1, j);
      mpz_mul(t, t, nj);
      mpz_mul(t, t, N[j>>1]);
      mpz_mul(num, num, D[j>>1]);
      mpz_addmul(num, t, den);
      mpz_mul(den, den, D[j>>1]);
      /* Reduce num/den */
      mpz_gcd(t, num, den);
      mpz_divexact(num, num, t);
      mpz_divexact(den, den, t);
    }
    /* Update n^j */
    mpz_mul(nj, nj, zn);
  }
  /* j = 1 */
  if (p >= 1) {
    mpz_mul_ui(t, nj, p+1);
    mpz_mul_ui(num, num, 2);
    mpz_addmul(num, t, den);
    mpz_mul_ui(den, den, 2);
    /* Skip reduction */
    mpz_mul(nj, nj, zn);
  }
  /* j = 0 */
  mpz_addmul(num, nj, den);
  /* Finalize */
  mpz_mul_ui(den, den, p+1);
  mpz_divexact(sum, num, den);

  mpz_clear(nj);
  mpz_clear(den);
  mpz_clear(num);
  mpz_clear(t);
}


static void _powerful_count_recurse(mpz_t sum, mpz_t n, unsigned long k,
                                   mpz_t m, unsigned long r, mpz_t t)
{
  mpz_t i, lim, newm;

  mpz_fdiv_q(t, n, m);
  mpz_root(t, t, r);
  if (r <= k) {
    mpz_add(sum, sum, t);
    return;
  }
  mpz_init_set(lim, t);
  mpz_init(newm);
  mpz_init(i);

  for (mpz_set_ui(i,1);  mpz_cmp(i,lim) <= 0;  mpz_add_ui(i,i,1)) {
    mpz_gcd(t, m, i);
    if (mpz_cmp_ui(t,1) == 0 && moebius(i) != 0) {
      mpz_pow_ui(t, i, r);
      mpz_mul(newm, m, t);
      if (r-1 == k) {
        mpz_fdiv_q(t, n, newm);
        mpz_root(t, t, k);
        mpz_add(sum, sum, t);
      } else {
        _powerful_count_recurse(sum, n, k, newm, r-1, t);
      }
    }
  }
  mpz_clear(i);
  mpz_clear(lim);
  mpz_clear(newm);
}

void powerful_count(mpz_t r, mpz_t n, unsigned long k)
{
  mpz_t m, t;

  mpz_set_ui(r, 0);
  if (k == 0)
    return;
  if (k == 1 || mpz_cmp_ui(n,1) <= 0) {
    mpz_set(r, n);
    return;
  }
  mpz_init_set_ui(m, 1);
  mpz_init(t);
  _powerful_count_recurse(r, n, k, m, 2*k-1, t);
  mpz_clear(t);
  mpz_clear(m);
}

void perfect_power_count(mpz_t r, mpz_t n)
{
  unsigned long k, log2n;
  mpz_t t;

  if (mpz_cmp_ui(n,1) <= 0) {
    mpz_set(r, n);
    return;
  }

  log2n = mpz_sizeinbase(n,2);
  mpz_init(t);
  mpz_set_ui(r, 1);
  for (k = 2; k <= log2n; k++) {
    int m;
    mpz_set_ui(t,k);
    m = moebius(t);
    if (m != 0) {
      mpz_root(t, n, k);
      mpz_sub_ui(t, t, 1);
      if (m < 0) mpz_add(r, r, t);
      else       mpz_sub(r, r, t);
    }
  }
  mpz_clear(t);
}
void perfect_power_count_range(mpz_t r, mpz_t lo, mpz_t hi) {
  if (mpz_cmp(lo, hi) > 0 || mpz_cmp_ui(hi,1) < 0) {
    mpz_set_ui(r, 0);
    return;
  }

  perfect_power_count(r, hi);

  if (mpz_cmp_ui(lo, 1) > 0) {
    mpz_t locount, lom1;
    mpz_init(locount);
    mpz_init(lom1);
    mpz_sub_ui(lom1, lo, 1);
    perfect_power_count(locount, lom1);
    mpz_sub(r, r, locount);
    mpz_clear(lom1);
    mpz_clear(locount);
  }
}

void prime_power_count(mpz_t r, mpz_t n)
{
  unsigned long k, log2n;
  mpz_t t1, t2;

  if (mpz_cmp_ui(n,5) <= 0) {
    if (mpz_cmp_ui(n,1) <= 0) mpz_set_ui(r, 0);
    else                      mpz_sub_ui(r, n, 1);
    return;
  }

  log2n = mpz_sizeinbase(n,2);
  mpz_init(t1);
  mpz_init(t2);
  prime_count(r, n);

  for (k = 2; k <= log2n; k++) {
    mpz_root(t1, n, k);
    prime_count(t2, t1);
    mpz_add(r, r, t2);
  }
  mpz_clear(t2);  mpz_clear(t1);
}
void prime_power_count_range(mpz_t r, mpz_t lo, mpz_t hi) {
  if (mpz_cmp(lo, hi) > 0 || mpz_cmp_ui(hi,2) < 0) {
    mpz_set_ui(r, 0);
    return;
  }

  prime_power_count(r, hi);

  if (mpz_cmp_ui(lo, 2) > 0) {
    mpz_t locount, lom1;
    mpz_init(locount);
    mpz_init(lom1);
    mpz_sub_ui(lom1, lo, 1);
    prime_power_count(locount, lom1);
    mpz_sub(r, r, locount);
    mpz_clear(lom1);
    mpz_clear(locount);
  }
}

void consecutive_integer_lcm(mpz_t m, unsigned long B)
{
  unsigned long i, p, p_power, pmin;
  mpz_t t[8];
  PRIME_ITERATOR(iter);

  /* Create eight sub-products to combine at the end. */
  for (i = 0; i < 8; i++)  mpz_init_set_ui(t[i], 1);
  i = 0;
  /* For each prime, multiply m by p^floor(log B / log p), which means
   * raise p to the largest power e such that p^e <= B.
   */
  if (B >= 2) {
    p_power = 2;
    while (p_power <= B/2)
      p_power *= 2;
    mpz_mul_ui(t[i&7], t[i&7], p_power);  i++;
  }
  p = prime_iterator_next(&iter);
  while (p <= B) {
    pmin = B/p;
    if (p > pmin)
      break;
    p_power = p*p;
    while (p_power <= pmin)
      p_power *= p;
    mpz_mul_ui(t[i&7], t[i&7], p_power);  i++;
    p = prime_iterator_next(&iter);
  }
  while (p <= B) {
    mpz_mul_ui(t[i&7], t[i&7], p);  i++;
    p = prime_iterator_next(&iter);
  }
  /* combine the eight sub-products */
  for (i = 0; i < 4; i++)  mpz_mul(t[i], t[2*i], t[2*i+1]);
  for (i = 0; i < 2; i++)  mpz_mul(t[i], t[2*i], t[2*i+1]);
  mpz_mul(m, t[0], t[1]);
  for (i = 0; i < 8; i++)  mpz_clear(t[i]);
  prime_iterator_destroy(&iter);
}


void exp_mangoldt(mpz_t res, mpz_t n)
{
  if (prime_power(res, n) < 1)
    mpz_set_ui(res, 1);
}

int is_carmichael(mpz_t n)
{
  mpz_t nm1, base, t, *factors;
  int i, res, nfactors, *exponents;

  /* small or even */
  if (mpz_cmp_ui(n,1105) < 0 || mpz_even_p(n))
    return mpz_cmp_ui(n,561)==0;

  mpz_init(nm1);
  mpz_sub_ui(nm1, n, 1);

  for (i = 1; i < 30; i++) {
    if (mpz_divisible_ui_p(n, sprimes[i])) {
      if (    mpz_divisible_ui_p(n, sprimes[i]*sprimes[i])
           || !mpz_divisible_ui_p(nm1, sprimes[i]-1)       )
        { mpz_clear(nm1); return 0; }
    }
  }

  mpz_init(t);
  mpz_init(base);

  /* Check 2^(n-1) mod n = 1 */
  res = (mpz_set_ui(t,2), mpz_powm(t,t,nm1,n), mpz_cmp_ui(t,1) == 0);

  /* Check with a few random small primes */
  for (i = 0; res && i < 5; i++) {
    mpz_random_nbit_prime(base, 32);
    if (mpz_divisible_p(n, base)) {
      if ( (mpz_mul(t, base, base), mpz_divisible_p(n, t)) ||
           (mpz_sub_ui(t, base, 1), !mpz_divisible_p(nm1, t)) )
        res = 0;
    } else {
      mpz_powm(t, t, nm1, n),
      res = (mpz_cmp_ui(t, 1) == 0);
    }
  }
  if (!res)
    { mpz_clear(base); mpz_clear(t);  mpz_clear(nm1);  return 0; }

  /* Decent chance that at this point it's prime or a Carmichael number. */

  if (mpz_sizeinbase(n,10) > 50) {      /* Probabilistic test */

    res = !_GMP_is_prime(n);  /* It must be a composite */
    for (i = 0; res && i < 128; i++) {
      mpz_sub_ui(t, n, 4);
      mpz_isaac_urandomm(base, t);
      mpz_add_ui(base, base, 3);    /* random base between 3 and n-2 */
      mpz_powm(t, base, nm1, n);
      res = (mpz_cmp_ui(t, 1) == 0);  /* if base^(n-1) mod n != 1, fail */
    }

  } else {                              /* Deterministic test (factor n) */

    nfactors = factor(n, &factors, &exponents);
    res = (nfactors > 2);                       /* must have 3+ factors */
    for (i = 0; res && i < nfactors; i++)       /* must be square free  */
      if (exponents[i] > 1)
        res = 0;
    for (i = 0; res && i < nfactors; i++)       /* p-1 | n-1 for all p */
      if (mpz_sub_ui(t, factors[i], 1), !mpz_divisible_p(nm1, t))
        res = 0;
    clear_factors(nfactors, &factors, &exponents);

  }
  mpz_clear(base);  mpz_clear(t);  mpz_clear(nm1);
  return res;
}

int is_fundamental(mpz_t n)
{
  int r, neg = (mpz_sgn(n) < 0), ret = 0;
  if (neg) mpz_neg(n,n);

  r = mpz_fdiv_ui(n,16);
  if (r != 0) {
    int r4 = r & 3;
    if (!neg && r4 == 1 && moebius(n) != 0) ret = 1;
    if ( neg && r4 == 3 && moebius(n) != 0) ret = 1;
    if (r4 == 0 && ((!neg && r != 4) || (neg && r != 12))) {
      mpz_t t;
      mpz_init(t);
      mpz_tdiv_q_2exp(t, n, 2);
      if (moebius(t) != 0)
        ret = 1;
      mpz_clear(t);
    }
  }
  if (neg) mpz_neg(n,n);
  return ret;
}

int is_practical(mpz_t n)
{
  mpz_t pke, fmult, prod, *factors;
  int i, j, nfactors, *exponents;

  if (mpz_cmp_ui(n,1) == 0) return 1;
  if (mpz_sgn(n) <= 0 || mpz_odd_p(n))  return 0;
  if (mpz_popcount(n) == 1) return 1;   /* n > 0 is a power of 2 */
  if (    !mpz_divisible_ui_p(n, 6)
       && !mpz_divisible_ui_p(n, 20)
       && !mpz_divisible_ui_p(n, 28)
       && !mpz_divisible_ui_p(n, 88)
       && !mpz_divisible_ui_p(n, 104)
       && !mpz_divisible_ui_p(n, 16) )
    return 0;

  nfactors = factor(n, &factors, &exponents);
  mpz_init(pke);
  mpz_init(fmult);
  mpz_init_set_ui(prod, 1);

  for (i = 1; i < nfactors; i++) {
    /* prod *= ipow(fac[i-1],exp[i-1]);  sum = 1 + divisor_sum(prod,1); */
    mpz_set(pke, factors[i-1]);
    mpz_add_ui(fmult, factors[i-1], 1);
    for (j = 1; j < exponents[i-1]; j++) {
      mpz_mul(pke, pke, factors[i-1]);
      mpz_add(fmult, fmult, pke);
    }
    mpz_mul(prod, prod, fmult);
    mpz_add_ui(fmult, prod, 1);
    if (mpz_cmp(factors[i], fmult) > 0)
      break;
  }
  clear_factors(nfactors, &factors, &exponents);
  mpz_clear(prod);
  mpz_clear(fmult);
  mpz_clear(pke);
  return (i >= nfactors);
}

int _totpred(mpz_t n, mpz_t maxd)
{
    int i, res, ndivisors;
    mpz_t N, r, d, p, *D;

    if (mpz_odd_p(n)) return 0;
    if (mpz_cmp_ui(n,2) == 0) return 1;
    if (mpz_popcount(n) == 1) return 1;   /* n > 0 is a power of 2 */

    mpz_init(N);  mpz_init(p);
    mpz_tdiv_q_2exp(N, n, 1);
    res = 0;
    mpz_add_ui(p, n, 1);
    if (mpz_cmp(N, maxd) < 0 && _GMP_is_prime(p)) {
      res = 1;
    } else {
      mpz_init(d);  mpz_init(r);
      D = divisor_list(&ndivisors, N);
      for (i = 0; res == 0 && i < ndivisors && mpz_cmp(D[i],maxd) < 0; i++) {
        mpz_set(d, D[i]);
        mpz_mul_2exp(p,d,1);  mpz_add_ui(p,p,1);
        if (!_GMP_is_prime(p))
          continue;
        mpz_divexact(r, N, d);
        while (1) {
          if (mpz_cmp(r, p) == 0 || _totpred(r, d)) {
            res = 1;
            break;
          }
          if (!mpz_divisible_p(r, p))
            break;
          mpz_divexact(r, r, p);
        }
      }
      mpz_clear(r);  mpz_clear(d);
      for (i = 0; i < ndivisors; i++)
        mpz_clear(D[i]);
      Safefree(D);
    }
    mpz_clear(p);  mpz_clear(N);
    return res;
}

int is_totient(mpz_t n)
{
  if (mpz_sgn(n) == 0 || mpz_odd_p(n))
    return !mpz_cmp_ui(n,1) ? 1 : 0;
  return _totpred(n, n);
}

void polygonal_nth(mpz_t r, mpz_t n, mpz_t k)
{
  mpz_t D, t;
  uint32_t ksmall = mpz_fits_uint_p(k) ? mpz_get_ui(k) : 0xFFFFFFFFU;

  if (mpz_sgn(n) < 0 || mpz_sgn(k) < 0) { mpz_set_ui(r,0); return; }
  if (ksmall < 3)                       { mpz_set_ui(r,0); return; }
  if (mpz_cmp_ui(n,1) <= 0)             { mpz_set_ui(r,1); return; }

  if (ksmall == 4) {
    if (mpz_perfect_square_p(n)) mpz_sqrt(r, n);
    else                         mpz_set_ui(r, 0);
    return;
  }

  mpz_init(D);
  mpz_init(t);

  if (ksmall == 3) {
    mpz_mul_2exp(D, n, 3);
    mpz_add_ui(D, D, 1);
  } else if (ksmall == 5) {
    mpz_mul_ui(D, n, 24); /* 8*k-16 = 24 if k=5 */
    mpz_add_ui(D, D, 1);
  } else {
    mpz_mul_ui(D, k, 8);
    mpz_sub_ui(D, D, 16);
    mpz_mul(D, D, n);
    mpz_sub_ui(t, k, 4);
    mpz_mul(t, t, t);
    mpz_add(D, D, t);
  }
  if (mpz_perfect_square_p(D)) {
    mpz_sqrt(D, D);
    if (ksmall == 3)  { mpz_sub_ui(D, D, 1); }
    else              { mpz_sub_ui(t, k, 4); mpz_add(D, D, t); }
    mpz_mul_ui(t, k, 2);
    mpz_sub_ui(t, t, 4);
    if (mpz_divisible_p(D, t)) {
      mpz_divexact(r, D, t);
      mpz_clear(t);
      mpz_clear(D);
      return;
    }
  }
  mpz_clear(t);
  mpz_clear(D);
  mpz_set_ui(r, 0);
}


static void word_tile(uint32_t* source, uint32_t from, uint32_t to) {
  while (from < to) {
    uint32_t words = (2*from > to) ? to-from : from;
    memcpy(source+from, source, sizeof(uint32_t)*words);
    from += words;
  }
}
static void sievep_ui(uint32_t* comp, UV pos, UV p, UV len, int verbose) {
  if (!(pos & 1)) pos += p;
  if (verbose > 3) {
    for ( ; pos < len; pos += 2*p ) {
      if (!TSTAVAL(comp, pos)) {
        printf("factor: %"UVuf" at %"UVuf"\n", p, pos);
        SETAVAL(comp, pos);
      }
    }
  } else {
    for ( ; pos < len; pos += 2*p )
      SETAVAL(comp, pos);
  }
}
/* Find first multiple of p after start */
#define sievep(comp, start, p, len, verbose) \
  sievep_ui(comp, (p) - mpz_fdiv_ui((start),(p)), p, len, verbose)

uint32_t* partial_sieve(mpz_t start, UV length, UV maxprime)
{
  uint32_t* comp;
  UV p, wlen, pwlen;
  int _verbose = get_verbose_level();
  PRIME_ITERATOR(iter);

  /* mpz_init(t);
     mpz_add_ui(t, start, (length & 1) ? length-1 : length-2);
     gmp_printf("partial sieve start %Zd  length %lu mark %Zd to %Zd\n", start, length, start, t); */
  MPUassert(mpz_odd_p(start), "partial sieve given even start");
  MPUassert(length > 0, "partial sieve given zero length");
  mpz_sub_ui(start, start, 1);
  if (length & 1) length++;

  /* Possibly reduce maxprime */
  if (mpz_cmp_ui(start, maxprime) <= 0) {
    mpz_t t;
    mpz_init(t);
    mpz_add_ui(t, start, length+1);
    mpz_sqrt(t, t);
    if (maxprime > mpz_get_ui(t))
      maxprime = mpz_get_ui(t);
    mpz_clear(t);
  }

  /* Allocate odds-only array in uint32_t units */
  wlen = (length+63)/64;
  New(0, comp, wlen, uint32_t);
  p = prime_iterator_next(&iter);

  /* Mark 3, 5, ... by tiling as long as we can. */
  pwlen = (wlen < 3) ? wlen : 3;
  memset(comp, 0x00, pwlen*sizeof(uint32_t));
  while (p <= maxprime) {
    sievep(comp, start, p, pwlen*64, _verbose);
    p = prime_iterator_next(&iter);
    if (pwlen*p >= wlen) break;
    word_tile(comp, pwlen, pwlen*p);
    pwlen *= p;
  }
  word_tile(comp, pwlen, wlen);

  /* Sieve up to their limit.
   *
   * Simple code for this:
   *    for ( ; p <= maxprime; p = prime_iterator_next(&iter))
   *      sievep(comp, start, p, length);
   * We'll save some time for large start values by doubling up primes.
   */
  {
    UV p1, p2;
    UV doublelim = (1UL << (sizeof(unsigned long) * 4)) - 1;
    UV ulim = (maxprime > ULONG_MAX) ? ULONG_MAX : maxprime;
    if (doublelim > maxprime) doublelim = maxprime;
    /* Do 2 primes at a time.  Fewer mpz remainders. */
    for ( p1 = p, p2 = prime_iterator_next(&iter);
          p2 <= doublelim;
          p1 = prime_iterator_next(&iter), p2 = prime_iterator_next(&iter) ) {
      UV p1p2 = p1 * p2;
      UV ddiv = mpz_fdiv_ui(start, p1p2);
      sievep_ui(comp, p1 - (ddiv % p1), p1, length, _verbose);
      sievep_ui(comp, p2 - (ddiv % p2), p2, length, _verbose);
    }
    if (p1 <= maxprime) sievep(comp, start, p1, length, _verbose);
    for (p = p2; p <= ulim; p = prime_iterator_next(&iter))
      sievep(comp, start, p, length, _verbose);
    if (p < maxprime) {
      /* UV is 64-bit, GMP's ui functions are 32-bit.  Sigh. */
      UV lastp, pos;
      mpz_t mp, rem;
      mpz_init(rem);
      mpz_init_set_ui(mp, (p >> 16) >> 16);
      mpz_mul_2exp(mp, mp, 32);
      mpz_add_ui(mp, mp, p & 0xFFFFFFFFUL);
      for (lastp = p;  p <= maxprime;  lastp=p, p=prime_iterator_next(&iter)) {
        mpz_add_ui(mp, mp, p-lastp);                 /* Calc mp = p */
        mpz_fdiv_r(rem, start, mp);                  /* Calc start % mp */
        if (mpz_cmp_ui(rem, ULONG_MAX) <= 0) {       /* pos = p - (start % p) */
          pos = p - mpz_get_ui(rem);
        } else {
          p1 = mpz_fdiv_q_ui(rem, rem, 2147483648UL);
          p1 += ((UV)mpz_get_ui(rem)) << 31;
          pos = p - p1;
        }
        sievep_ui(comp, pos, p, length, _verbose);
      }
      mpz_clear(mp);
      mpz_clear(rem);
    }
  }

  prime_iterator_destroy(&iter);
  return comp;
}


/*****************************************************************************/

static unsigned short small_prime_count(unsigned short n)
{
  unsigned long i;
  for (i = 0; i < NSMALLPRIMES; i++)
    if (n < sprimes[i])
      break;
  return i;
}

void prime_count_lower(mpz_t pc, mpz_t n)
{
  mpf_t x, logx, logx2, t, s;
  unsigned long bits = 7 + DIGS2BITS(mpz_sizeinbase(n,10));
  unsigned long N = mpz_get_ui(n);

  if (mpz_cmp_ui(n, 1009) < 0)
    { mpz_set_ui(pc, small_prime_count(N)); return; }

  mpf_init2(x,     bits);
  mpf_init2(logx,  bits);
  mpf_init2(logx2, bits);
  mpf_init2(t,     bits);
  mpf_init2(s,     bits);

  mpf_set_z(x, n);
  mpf_log(logx, x);
  mpf_mul(logx2, logx, logx);

  if (mpz_cmp_ui(n, 300000) < 0) {
    double a = (N <  33000) ? 1190
             : (N <  70200) ? 947
             : (N < 176000) ? 904
                            : 829;
    mpf_set(s, logx);
    mpf_sub_ui(s, s, 1);
    mpf_ui_div(t, 1, logx);
    mpf_sub(s, s, t);
    mpf_set_d(t, 2.85);
    mpf_div(t, t, logx2);
    mpf_sub(s, s, t);
    mpf_set_d(t, 13.15);
    mpf_mul(logx2, logx2, logx);
    mpf_div(t, t, logx2);
    mpf_sub(s, s, t);
    mpf_set_d(t, a);
    mpf_mul(logx2, logx2, logx);
    mpf_div(t, t, logx2);
    mpf_add(s, s, t);
    mpf_div(x, x, s);
  } else if (mpf_cmp_d(x, 1e19) < 0) { /* Büthe 2015 1.9    1511.02032v2.pdf */
    double a = 2.50;
    double b = (N <     88783) ?   4.0
             : (N <    300000) ?  -3.0
             : (N <    303000) ?   5.0
             : (N <   1100000) ?  -7.0
             : (N <   4500000) ? -37.0
             : (N <  10200000) ? -70.0
             : (N <  36900000) ? -53.0
             : (N <  38100000) ? -29.0
             :                   -84.0;
    mpf_set_str(s, "1.95", 10);
    if (N >= 4000000000UL) {
      mpf_set_str(t, "3.9", 10);
      mpf_div(t, t, logx);
      mpf_add(s, s, t);
      mpf_set_str(t, "19.5", 10);
      mpf_div(t, t, logx2);
      mpf_add(s, s, t);
    } else {
      mpf_set_d(t, a);
      mpf_div(t, t, logx);
      mpf_add(s, s, t);
      mpf_set_d(t, b);
      mpf_div(t, t, logx2);
      mpf_add(s, s, t);
    }
    mpf_sqrt(t, x);
    mpf_div(t, t, logx);
    mpf_mul(s, s, t);

    li(t, x, 20);
    mpf_sub(x, t, s);

  } else if (mpf_cmp_d(x, 5.5e25) < 0) { /* Büthe 2014 v3 7.2 1410.7015v3.pdf */
    mpf_sqrt(t, x);                      /* Axler 2017 2.2    1703.08032.pdf */
    mpf_mul(s, logx, t);
    const_pi(t, 30);
    mpf_mul_2exp(t, t, 3);
    mpf_div(s, s, t);
    li(t, x, 30);
    mpf_sub(x, t, s);
  } else { /* Axler 2017 1.3  https://arxiv.org/pdf/1703.08032.pdf */
    mpf_set(s, logx);
    mpf_sub_ui(s, s, 1);
    mpf_ui_div(t, 1, logx);
    mpf_sub(s, s, t);
    mpf_set_str(t, "2.85", 10);
    mpf_div(t, t, logx2);
    mpf_sub(s, s, t);
    mpf_set_str(t, "13.15", 10);
    mpf_mul(logx2, logx2, logx);
    mpf_div(t, t, logx2);
    mpf_sub(s, s, t);
    mpf_set_str(t, "70.7", 10);
    mpf_mul(logx2, logx2, logx);
    mpf_div(t, t, logx2);
    mpf_sub(s, s, t);
    mpf_set_str(t, "458.7275", 10);
    mpf_mul(logx2, logx2, logx);
    mpf_div(t, t, logx2);
    mpf_sub(s, s, t);
    mpf_set_str(t, "3428.7225", 10);
    mpf_mul(logx2, logx2, logx);
    mpf_div(t, t, logx2);
    mpf_sub(s, s, t);
    mpf_div(x, x, s);
  }
  if (!mpf_integer_p(x)) mpf_add_ui(x, x, 1);  /* ceil */
  mpz_set_f(pc,x);
  mpf_clear(logx2); mpf_clear(logx); mpf_clear(x); mpf_clear(t); mpf_clear(s);
}
void prime_count_upper(mpz_t pc, mpz_t n)
{
  mpf_t x, logx, logx2, t, s;
  unsigned long bits = 7 + DIGS2BITS(mpz_sizeinbase(n,10));
  unsigned long N = mpz_get_ui(n);

  if (mpz_cmp_ui(n, 1009) < 0)
    { mpz_set_ui(pc, small_prime_count(N)); return; }

  if (mpz_cmp_ui(n, 15900) < 0) {
    if (N < 7) {
      mpz_set_ui(pc, 0 + (N >= 2) + (N >= 3) + (N >= 5));
    } else {
      double a = (N < 1621) ? 1.048 : (N < 5000) ? 1.071 : 1.098;
      double X = 1.0 + ((double)N) / (log((double)N) - a);
      mpz_set_d(pc, X);
    }
    return;
  }

  mpf_init2(x,     bits);
  mpf_init2(logx,  bits);
  mpf_init2(logx2, bits);
  mpf_init2(t,     bits);
  mpf_init2(s,     bits);

  mpf_set_z(x, n);
  mpf_log(logx, x);
  mpf_mul(logx2, logx, logx);

  if (mpz_cmp_ui(n, 821800000UL) < 0) {
    double a = (N <    356000) ? 2.54
             : (N <  48000000) ? 2.51
             : (N < 400000000) ? 2.47
                               : 2.37;
    mpf_set_ui(s, 1);
    mpf_ui_div(t, 1, logx);
    mpf_add(s, s, t);
    mpf_set_d(t, a);
    mpf_div(t, t, logx2);
    mpf_add(s, s, t);

    mpf_div(t, x, logx);
    mpf_mul(x, t, s);
  } else if (mpf_cmp_d(x, 1e19) < 0) { /* Büthe 2015 1.10    1511.02032v2.pdf */
    double a = (mpf_cmp_d(x,   1100000000.0) < 0) ? 0.032
             : (mpf_cmp_d(x,  10010000000.0) < 0) ? 0.027
             : (mpf_cmp_d(x, 101260000000.0) < 0) ? 0.021
                                                  : 0.0;
    if (a > 0) {
      mpf_sqrt(t, x);
      mpf_mul(s, logx, t);
      mpf_set_d(t, a);
      mpf_mul(s, s, t);
      const_pi(t, 25);
      mpf_mul_2exp(t, t, 3);
      mpf_div(s, s, t);
      li(t, x, 25);
      mpf_sub(x, t, s);
    } else {
      li(x, x, 25);
    }
  } else if (mpf_cmp_d(x, 5.5e25) < 0) { /* Büthe 2014 v3 7.2 1410.7015v3.pdf */
    mpf_sqrt(t, x);                      /* Axler 2017 2.2    1703.08032.pdf */
    mpf_mul(s, logx, t);
    const_pi(t, 30);
    mpf_mul_2exp(t, t, 3);
    mpf_div(s, s, t);
    li(t, x, 30);
    mpf_add(x, t, s);
  } else { /* Axler 2017 1.2  https://arxiv.org/pdf/1703.08032.pdf */
    mpf_set(s, logx);
    mpf_sub_ui(s, s, 1);
    mpf_ui_div(t, 1, logx);
    mpf_sub(s, s, t);
    mpf_set_str(t, "3.15", 10);
    mpf_div(t, t, logx2);
    mpf_sub(s, s, t);
    mpf_set_str(t, "12.85", 10);
    mpf_mul(logx2, logx2, logx);
    mpf_div(t, t, logx2);
    mpf_sub(s, s, t);
    mpf_set_str(t, "71.3", 10);
    mpf_mul(logx2, logx2, logx);
    mpf_div(t, t, logx2);
    mpf_sub(s, s, t);
    mpf_set_str(t, "463.2275", 10);
    mpf_mul(logx2, logx2, logx);
    mpf_div(t, t, logx2);
    mpf_sub(s, s, t);
    mpf_set_str(t, "4585", 10);
    mpf_mul(logx2, logx2, logx);
    mpf_div(t, t, logx2);
    mpf_sub(s, s, t);
    mpf_div(x, x, s);
  }
  /* floor */
  mpz_set_f(pc,x);
  mpf_clear(logx2); mpf_clear(logx); mpf_clear(x); mpf_clear(t); mpf_clear(s);
}
/*****************************************************************************/

void prime_count_range(mpz_t count, mpz_t ilo, mpz_t ihi)
{
  uint32_t* comp;
  UV i, cnt, hbits, depth, width;
  mpz_t lo, hi, shi, t;

  mpz_set_ui(count, 0);

  if (mpz_cmp(ilo, ihi) > 0 || mpz_cmp_ui(ihi,2) < 0) return;

  if (mpz_cmp_ui(ihi, 1009) < 0) {
     uint32_t ulo = mpz_get_ui(ilo), uhi = mpz_get_ui(ihi);
     uint32_t cnt = small_prime_count(uhi)
                  - ((ulo <= 2) ? 0 : small_prime_count(ulo-1));
     mpz_set_ui(count, cnt);
     return;
  }

  mpz_init_set(lo, ilo);
  mpz_init_set(hi, ihi);

  if (mpz_cmp_ui(lo,2) <= 0) {
    if (mpz_cmp_ui(hi,2) >= 0) mpz_add_ui(count,count,1);
    mpz_set_ui(lo,3);
  }
  if (mpz_cmp(lo,hi) > 0) { mpz_clear(lo); mpz_clear(hi); return; }

  mpz_init(t);
  if (mpz_add_ui(t, lo, 100000), mpz_cmp(t,hi) > 0) {
    mpz_sub_ui(lo,lo,1);
    for (cnt = 0; mpz_cmp(lo, hi) <= 0; _GMP_next_prime(lo))
      cnt++;
    mpz_add_ui(count,count,cnt-1);
    mpz_clear(t); mpz_clear(lo); mpz_clear(hi);
    return;
  }

  /* Determine a reasonable depth for sieve. */
  hbits = mpz_sizeinbase(hi,2);
  depth = (hbits < 100) ? 50000000UL : ((UV)hbits)*UVCONST(500000);
  if (hbits < 64) {
    mpz_sqrt(t, hi);
    if (mpz_cmp_ui(t,depth) < 0)
      depth = mpz_get_ui(t);
  }

  /* Count small inputs directly.  We should do this faster. */
  if (mpz_cmp_ui(lo, depth) <= 0) {
    mpz_sub_ui(lo,lo,1);
    for (cnt = 0; mpz_cmp_ui(lo, depth) <= 0; _GMP_next_prime(lo))
      cnt++;
    mpz_add_ui(count,count,cnt-1);
  }

  /* Ensure endpoints are odd */
  if (mpz_even_p(lo)) mpz_add_ui(lo,lo,1);
  if (mpz_even_p(hi)) mpz_sub_ui(hi,hi,1);

  mpz_init(shi);  /* segment hi */
  while (mpz_cmp(lo,hi) <= 0) {
    mpz_add_ui(shi, lo, 1000000000UL);
    if (mpz_cmp(shi, hi) > 0)
      mpz_set(shi, hi);
    mpz_sub(t, shi, lo);
    width = mpz_get_ui(t) + 1;

    comp = partial_sieve(lo, width, depth);
    for (i = 1, cnt = 0; i <= width; i += 2) {
      if (!TSTAVAL(comp, i) && (mpz_add_ui(t, lo, i), _GMP_BPSW(t)))
        cnt++;
    }
    Safefree(comp);

    mpz_add_ui(lo, shi, 2);
    mpz_add_ui(count, count, cnt);
  }
  mpz_clear(shi);
  mpz_clear(t); mpz_clear(lo); mpz_clear(hi);
}

void prime_count(mpz_t count, mpz_t hi)
{
  mpz_t lo;
  mpz_init_set_ui(lo, 0);
  prime_count_range(count, lo, hi);
  mpz_clear(lo);
}


typedef struct {
  uint32_t nmax;
  uint32_t nsize;
  UV* list;
} vlist;
#define INIT_VLIST(v) \
  v.nsize = 0; \
  v.nmax = 1024; \
  New(0, v.list, v.nmax, UV);
#define RESIZE_VLIST(v, size) \
  do { if (v.nmax < size) Renew(v.list, v.nmax = size, UV); } while (0)
#define PUSH_VLIST(v, n) \
  do { \
    if (v.nsize >= v.nmax) \
      Renew(v.list, v.nmax += 1024, UV); \
    v.list[v.nsize++] = n; \
  } while (0)

#define ADDVAL32(v, n, max, val) \
  do { if (n >= max) Renew(v, max += 1024, UV);  v[n++] = val; } while (0)
#define SWAPL32(l1, n1, m1,  l2, n2, m2) \
  { UV t_, *u_ = l1;  l1 = l2;  l2 = u_; \
            t_ = n1;  n1 = n2;  n2 = t_; \
            t_ = m1;  m1 = m2;  m2 = t_; }

UV* sieve_primes(mpz_t inlow, mpz_t high, UV k, UV *rn) {
  mpz_t t, low;
  int test_primality = 0, k_primality = 0, force_full = 0, width2, hbits;
  uint32_t* comp;
  vlist retlist;

  if (mpz_cmp_ui(inlow, 2) < 0) mpz_set_ui(inlow, 2);
  if (mpz_cmp(inlow, high) > 0) { *rn = 0; return 0; }

  mpz_init(t);
  mpz_sub(t, high, inlow);
  width2 = mpz_sizeinbase(t,2);
  hbits = mpz_sizeinbase(high,2);

  mpz_sqrt(t, high);           /* No need for k to be > sqrt(high) */
  /* If auto-setting k or k >= sqrt(n), pick a good depth and test primality */
  if (k == 0 || mpz_cmp_ui(t, k) <= 0) {
    test_primality = 1;
    k = (hbits < 100) ? 50000000UL : ((UV)hbits)*UVCONST(500000);
    /* For small widths we're much better off with smaller sieves. */
    if (width2 <= 23)  k /= 2;
    if (width2 <= 21)  k /= 2;
    if (width2 <= 19)  k /= 2;
    if (width2 <= 17)  k /= 2;
    if (width2 <= 15)  k /= 2;
    if (width2 <= 13)  k /= 2;
    if (width2 <= 11)  k /= 2;
  }
  /* For smallist n and large ranges, force full sieve */
  if (test_primality && hbits >= 40 && hbits <  56 && (width2 * 2.8) >= hbits)
    force_full = 1;
  if (test_primality && hbits >= 56 && hbits <  64 && (width2 * 2.6) >= hbits)
    force_full = 1;
  if (test_primality && hbits >= 64 && hbits <= 82 && (width2 * 2.5) >= hbits && sizeof(unsigned long) > 4)
    force_full = 1;
  /* If k >= sqrtn, sieving is enough.  Use k=sqrtn, turn off post-sieve test */
  if (force_full || mpz_cmp_ui(t, k) <= 0) {
    k = mpz_get_ui(t);
    k_primality = 1;           /* Our sieve is complete */
    test_primality = 0;        /* Don't run BPSW */
  }

  INIT_VLIST(retlist);

  /* If we want small primes, do it quickly */
  if ( (k_primality || test_primality) && mpz_cmp_ui(high,2000000000U) <= 0 ) {
    UV ulow = mpz_get_ui(inlow), uhigh = mpz_get_ui(high);
    if (uhigh < 1000000U || uhigh/ulow >= 4) {
      UV n, Pi, *primes;
      primes = sieve_to_n(mpz_get_ui(high), &Pi);
      RESIZE_VLIST(retlist, Pi);
      for (n = 0; n < Pi; n++) {
        if (primes[n] >= ulow)
          PUSH_VLIST(retlist, primes[n]-ulow);
      }
      Safefree(primes);
      mpz_clear(t);
      *rn = retlist.nsize;
      return retlist.list;
    }
  }

  mpz_init_set(low, inlow);
  if (k < 2) k = 2;   /* Should have been handled by quick return */

  /* Include all primes up to k, since they will get filtered */
  if (mpz_cmp_ui(low, k) <= 0) {
    UV n, Pi, *primes, ulow = mpz_get_ui(low);
    primes = sieve_to_n(k, &Pi);
    RESIZE_VLIST(retlist, Pi);
    for (n = 0; n < Pi; n++) {
      if (primes[n] >= ulow)
        PUSH_VLIST(retlist, primes[n]-ulow);
    }
    Safefree(primes);
  }

  if (mpz_even_p(low))           mpz_add_ui(low, low, 1);
  if (mpz_even_p(high))          mpz_sub_ui(high, high, 1);

  if (mpz_cmp(low, high) <= 0) {
    UV i, length, offset;
    mpz_sub(t, high, low); length = mpz_get_ui(t) + 1;
    /* Get bit array of odds marked with composites(k) marked with 1 */
    comp = partial_sieve(low, length, k);
    mpz_sub(t, low, inlow); offset = mpz_get_ui(t);
    for (i = 1; i <= length; i += 2) {
      if (!TSTAVAL(comp, i)) {
        if (!test_primality || (mpz_add_ui(t,low,i),_GMP_BPSW(t)))
          PUSH_VLIST(retlist, i - offset);
      }
    }
    Safefree(comp);
  }

  mpz_clear(low);
  mpz_clear(t);
  *rn = retlist.nsize;
  return retlist.list;
}

void next_twin_prime(mpz_t res, mpz_t n) {
  mpz_t low, t;

  mpz_init(t);
  if (mpz_cmp_ui(n, 1000000) < 0) {
    uint32_t p, ulow = mpz_get_ui(n), last = 0;
    PRIME_ITERATOR(iter);
    prime_iterator_setprime(&iter, ulow);
    while (1) {
      p = prime_iterator_next(&iter);
      if (p-last == 2 && last > 0) {
        mpz_set_ui(res, last);
        break;
      }
      last = p;
    }
    prime_iterator_destroy(&iter);
  } else {
    UV i, starti, l2, length, depth, found = 0;
    uint32_t* comp;

    mpz_init(low);
    mpz_add_ui(low, n, 1);
    if (mpz_even_p(low))  mpz_add_ui(low, low, 1);

    l2 = mpz_sizeinbase(low,2);
    if (l2 <= 6000) {
      depth = (l2/160.0) * l2 * l2;
      length = 3 * 0.634 * l2 * l2;  /* we will resieve sometimes */
      if (length % 2) length++;  /* low is odd, length must be even */
    } else {
      depth  = UVCONST(1350000000);
      length = UVCONST(  91296000);
    }

    while (!found) {
      starti = (6 - mpz_fdiv_ui(low,6)) % 6;
      comp = partial_sieve(low, length + 2, depth);
      for (i = starti; i <= length && !found; i += 6) {
        if (!TSTAVAL(comp, i) && !TSTAVAL(comp, i+2)) {
          if ( (mpz_add_ui(t,low,i),  miller_rabin_ui(t,2)) &&
               (mpz_add_ui(t,t,2),    miller_rabin_ui(t,2)) &&
               (mpz_add_ui(t,low,i),  _GMP_is_lucas_pseudoprime(t,2)) &&
               (mpz_add_ui(t,t,2),    _GMP_is_lucas_pseudoprime(t,2)) ) {
            mpz_add_ui(res, low, i);
            found = 1;
          }
        }
      }
      Safefree(comp);
      mpz_add_ui(low, low, length+1);
    }
    mpz_clear(low);
  }
  mpz_clear(t);
}

UV* sieve_twin_primes(mpz_t low, mpz_t high, UV twin, UV *rn) {
  mpz_t t;
  UV i, length, k, starti = 1, skipi = 2;
  uint32_t* comp;
  vlist retlist;

  /* Twin offset will always be even, so we will never return 2 */
  MPUassert( !(twin & 1), "twin prime offset is even" );
  if (mpz_cmp_ui(low, 3) <= 0) mpz_set_ui(low, 3);

  if (mpz_even_p(low))  mpz_add_ui(low, low, 1);
  if (mpz_even_p(high)) mpz_sub_ui(high, high, 1);

  i = twin % 6;
  if (i == 2 || i == 4) { skipi = 6; starti = ((twin%6)==2) ? 5 : 1; }

  /* If no way to return any more results, leave now */
  if (mpz_cmp(low, high) > 0 || (i == 1 || i == 3 || i == 5))
    { *rn = 0; return 0; }

  INIT_VLIST(retlist);
  mpz_init(t);

  /* Use a much higher k value than we do for primes */
  k = 80000 * mpz_sizeinbase(high,2);
  /* No need for k to be > sqrt(high) */
  mpz_sqrt(t, high);
  if (mpz_cmp_ui(t, k) < 0)
    k = mpz_get_ui(t);

  /* Handle small primes that will get sieved out */
  if (mpz_cmp_ui(low, k) <= 0) {
    uint32_t p, ulow = mpz_get_ui(low);
    PRIME_ITERATOR(iter);
    for (p = 2; p <= k; p = prime_iterator_next(&iter)) {
      if (p < ulow) continue;
      if (mpz_set_ui(t, p+twin), _GMP_BPSW(t))
        PUSH_VLIST(retlist, p-ulow+1);
    }
    prime_iterator_destroy(&iter);
  }

  mpz_sub(t, high, low);
  length = mpz_get_uv(t) + 1;
  starti = ((starti+skipi) - mpz_fdiv_ui(low,skipi) + 1) % skipi;

  /* Get bit array of odds marked with composites(k) marked with 1 */
  comp = partial_sieve(low, length + twin, k);
  for (i = starti; i <= length; i += skipi) {
    if (!TSTAVAL(comp, i) && !TSTAVAL(comp, i+twin)) {
      /* Add to list if both t,t+2 pass MR and if both pass ES Lucas */
      if ( (mpz_add_ui(t,low,i),  miller_rabin_ui(t,2)) &&
           (mpz_add_ui(t,t,twin), miller_rabin_ui(t,2)) &&
           (mpz_add_ui(t,low,i),  _GMP_is_lucas_pseudoprime(t,2)) &&
           (mpz_add_ui(t,t,twin), _GMP_is_lucas_pseudoprime(t,2)) )
        PUSH_VLIST(retlist, i);
    }
  }
  Safefree(comp);
  mpz_clear(t);
  *rn = retlist.nsize;
  return retlist.list;
}


#define addmodded(r,a,b,n)  do { r = a + b; if (r >= n) r -= n; } while(0)

UV* sieve_cluster(mpz_t low, mpz_t high, uint32_t* cl, UV nc, UV *rn) {
  mpz_t t, savelow;
  vlist retlist;
  UV i, ppr, nres, allocres;
  uint32_t const targres = 4000000;
  uint32_t const maxpi = 168;
  UV *residues, *cres;
  uint32_t pp_0, pp_1, pp_2, *resmod_0, *resmod_1, *resmod_2;
  uint32_t rem_0, rem_1, rem_2, remadd_0, remadd_1, remadd_2;
  uint32_t pi, startpi = 1;
  uint32_t lastspr = sprimes[maxpi-1];
  uint32_t c, smallnc;
  UV ibase = 0, num_mr = 0, num_lucas = 0;
  char crem_0[53*59], crem_1[61*67], crem_2[71*73], *VPrem;
  int run_pretests = 0;
  int _verbose = get_verbose_level();

  if (nc == 1) return sieve_primes(low, high, 0, rn);
  if (nc == 2) return sieve_twin_primes(low, high, cl[1], rn);

  if (mpz_even_p(low))           mpz_add_ui(low, low, 1);
  if (mpz_even_p(high))          mpz_sub_ui(high, high, 1);

  if (mpz_cmp(low, high) > 0) { *rn = 0; return 0; }

  INIT_VLIST(retlist);
  mpz_init(t);

  /* Handle small values that would get sieved away */
  if (mpz_cmp_ui(low, lastspr) <= 0) {
    UV ui_low = mpz_get_ui(low);
    UV ui_high = (mpz_cmp_ui(high,lastspr) > 0) ? lastspr : mpz_get_ui(high);
    for (pi = 0; pi < maxpi; pi++) {
      UV p = sprimes[pi];
      if (p > ui_high) break;
      if (p < ui_low) continue;
      for (c = 1; c < nc; c++)
        if (!(mpz_set_ui(t, p+cl[c]), _GMP_is_prob_prime(t))) break;
      if (c == nc)
        PUSH_VLIST(retlist, p-ui_low+1);
    }
  }
  if (mpz_odd_p(low)) mpz_sub_ui(low, low, 1);
  if (mpz_cmp_ui(high, lastspr) <= 0) {
    mpz_clear(t);
    *rn = retlist.nsize;
    return retlist.list;
  }

  /* Determine the primorial size and acceptable residues */
  New(0, residues, allocres = 1024, UV);
  {
    UV remr, *res2, allocres2, nres2, maxppr;
    /* Calculate residues for a small primorial */
    for (pi = 2, ppr = 1, i = 0;  i <= pi;  i++) ppr *= sprimes[i];
    remr = mpz_fdiv_ui(low, ppr);
    nres = 0;
    for (i = 1; i <= ppr; i += 2) {
      UV remi = remr + i;
      for (c = 0; c < nc; c++) {
        if (gcd_ui(remi + cl[c], ppr) != 1) break;
      }
      if (c == nc)
        ADDVAL32(residues, nres, allocres, i);
    }
    /* Raise primorial size until we have plenty of residues */
    New(0, res2, allocres2 = 1024, UV);
    mpz_sub(t, high, low);
    maxppr = (mpz_sizeinbase(t,2) >= BITS_PER_WORD) ? UV_MAX : (UVCONST(1) << mpz_sizeinbase(t,2));
#if BITS_PER_WORD == 64
    while (pi++ < 14) {
#else
    while (pi++ < 8) {
#endif
      uint32_t j, p = sprimes[pi];
      UV r, newppr = ppr * p;
      if (nres == 0 || nres > targres/(p/2) || newppr > maxppr) break;
      if (_verbose > 1) printf("cluster sieve found %"UVuf" residues mod %"UVuf"\n", nres, ppr);
      remr = mpz_fdiv_ui(low, newppr);
      nres2 = 0;
      for (i = 0; i < p; i++) {
        for (j = 0; j < nres; j++) {
          r = i*ppr + residues[j];
          for (c = 0; c < nc; c++) {
            UV v = remr + r + cl[c];
            if ((v % p) == 0) break;
          }
          if (c == nc)
            ADDVAL32(res2, nres2, allocres2, r);
        }
      }
      ppr = newppr;
      SWAPL32(residues, nres, allocres,  res2, nres2, allocres2);
    }
    startpi = pi;
    Safefree(res2);
  }
  if (_verbose) printf("cluster sieve using %"UVuf" residues mod %"UVuf"\n", nres, ppr);

  /* Return if not admissible, maybe with a single small value */
  if (nres == 0) {
    Safefree(residues);
    mpz_clear(t);
    *rn = retlist.nsize;
    return retlist.list;
  }

  mpz_init_set(savelow, low);
  if (mpz_sizeinbase(low, 2) > 310) run_pretests = 1;
  if (run_pretests && mpz_sgn(_bgcd2) == 0) {
    _GMP_pn_primorial(_bgcd2, BGCD2_PRIMES);
    mpz_divexact(_bgcd2, _bgcd2, _bgcd);
  }

  /* Pre-mod the residues with first two primes for fewer modulos every chunk */
  {
    uint32_t p1 = sprimes[startpi+0], p2 = sprimes[startpi+1];
    uint32_t p3 = sprimes[startpi+2], p4 = sprimes[startpi+3];
    uint32_t p5 = sprimes[startpi+4], p6 = sprimes[startpi+5];
    pp_0 = p1*p2; pp_1 = p3*p4; pp_2 = p5*p6;
    memset(crem_0, 1, pp_0);
    memset(crem_1, 1, pp_1);
    memset(crem_2, 1, pp_2);
    /* Mark remainders that indicate a composite for this residue. */
    for (i = 0; i < p1; i++) { crem_0[i*p1]=0; crem_0[i*p2]=0; }
    for (     ; i < p2; i++) { crem_0[i*p1]=0;                }
    for (i = 0; i < p3; i++) { crem_1[i*p3]=0; crem_1[i*p4]=0; }
    for (     ; i < p4; i++) { crem_1[i*p3]=0;                }
    for (i = 0; i < p5; i++) { crem_2[i*p5]=0; crem_2[i*p6]=0; }
    for (     ; i < p6; i++) { crem_2[i*p5]=0;                }
    for (c = 1; c < nc; c++) {
      uint32_t c1=cl[c], c2=cl[c], c3=cl[c], c4=cl[c], c5=cl[c], c6=cl[c];
      if (c1 >= p1) c1 %= p1;
      if (c2 >= p2) c2 %= p2;
      for (i = 1; i <= p1; i++) { crem_0[i*p1-c1]=0; crem_0[i*p2-c2]=0; }
      for (     ; i <= p2; i++) { crem_0[i*p1-c1]=0;                   }
      if (c3 >= p3) c3 %= p3;
      if (c4 >= p4) c4 %= p4;
      for (i = 1; i <= p3; i++) { crem_1[i*p3-c3]=0; crem_1[i*p4-c4]=0; }
      for (     ; i <= p4; i++) { crem_1[i*p3-c3]=0;                   }
      if (c5 >= p5) c5 %= p5;
      if (c6 >= p6) c6 %= p6;
      for (i = 1; i <= p5; i++) { crem_2[i*p5-c5]=0; crem_2[i*p6-c6]=0; }
      for (     ; i <= p6; i++) { crem_2[i*p5-c5]=0;                   }
    }
    New(0, resmod_0, nres, uint32_t);
    New(0, resmod_1, nres, uint32_t);
    New(0, resmod_2, nres, uint32_t);
    for (i = 0; i < nres; i++) {
      resmod_0[i] = residues[i] % pp_0;
      resmod_1[i] = residues[i] % pp_1;
      resmod_2[i] = residues[i] % pp_2;
    }
  }

  /* Precalculate acceptable residues for more primes */
  MPUassert( lastspr <= 1024, "cluster sieve internal" );
  New(0, VPrem, maxpi * 1024, char);
  memset(VPrem, 1, maxpi * 1024);
  for (pi = startpi+6; pi < maxpi; pi++)
    VPrem[pi*1024] = 0;
  for (pi = startpi+6, smallnc = 0; pi < maxpi; pi++) {
    uint32_t p = sprimes[pi];
    char* prem = VPrem + pi*1024;
    while (smallnc < nc && cl[smallnc] < p)   smallnc++;
    for (c = 1; c < smallnc; c++) prem[p-cl[c]] = 0;
    for (     ; c <      nc; c++) prem[p-(cl[c]%p)] = 0;
  }

  New(0, cres, nres, UV);

  rem_0 = mpz_fdiv_ui(low,pp_0);  remadd_0 = ppr % pp_0;
  rem_1 = mpz_fdiv_ui(low,pp_1);  remadd_1 = ppr % pp_1;
  rem_2 = mpz_fdiv_ui(low,pp_2);  remadd_2 = ppr % pp_2;

  /* Loop over their range in chunks of size 'ppr' */
  while (mpz_cmp(low, high) <= 0) {

    uint32_t r, nr, remr, ncres;
    unsigned long ui_low = (mpz_sizeinbase(low,2) > 8*sizeof(unsigned long)) ? 0 : mpz_get_ui(low);

    /* Reduce the allowed residues for this chunk using more primes */

    { /* Start making a list of this chunk's residues using three pairs */
      for (r = 0, ncres = 0; r < nres; r++) {
        addmodded(remr, rem_0, resmod_0[r], pp_0);
        if (crem_0[remr]) {
          addmodded(remr, rem_1, resmod_1[r], pp_1);
          if (crem_1[remr]) {
            addmodded(remr, rem_2, resmod_2[r], pp_2);
            if (crem_2[remr]) {
              cres[ncres++] = residues[r];
            }
          }
        }
      }
      addmodded(rem_0, rem_0, remadd_0, pp_0);
      addmodded(rem_1, rem_1, remadd_1, pp_1);
      addmodded(rem_2, rem_2, remadd_2, pp_2);
    }

    /* Sieve through more primes one at a time, removing residues. */
    for (pi = startpi+6; pi < maxpi && ncres > 0; pi++) {
      uint32_t p = sprimes[pi];
      uint32_t rem = (ui_low) ? (ui_low % p) : mpz_fdiv_ui(low,p);
      char* prem = VPrem + pi*1024;
      /* Check divisibility of each remaining residue with this p */
      if (startpi <= 9 || cres[ncres-1] < 4294967295U) {   /* Residues are 32-bit */
        for (r = 0, nr = 0; r < ncres; r++) {
          if (prem[ (rem+(uint32_t)cres[r]) % p ])
            cres[nr++] = cres[r];
        }
      } else {              /* Residues are 64-bit */
        for (r = 0, nr = 0; r < ncres; r++) {
          if (prem[ (rem+cres[r]) % p ])
            cres[nr++] = cres[r];
        }
      }
      ncres = nr;
    }
    if (_verbose > 2) printf("cluster sieve range has %u residues left\n", ncres);

    /* Now check each of the remaining residues for inclusion */
    for (r = 0; r < ncres; r++) {
      i = cres[r];
      mpz_add_ui(t, low, i);
      if (mpz_cmp(t, high) > 0) break;
      /* Pretest each element if the input is large enough */
      if (run_pretests) {
        for (c = 0; c < nc; c++)
          if (mpz_add_ui(t, low, i+cl[c]), mpz_gcd(t,t,_bgcd2), mpz_cmp_ui(t,1)) break;
        if (c != nc) continue;
      }
      /* PRP test.  Split BPSW in two for faster rejection. */
      for (c = 0; c < nc; c++)
        if (! (mpz_add_ui(t, low, i+cl[c]), num_mr++, miller_rabin_ui(t,2)) ) break;
      if (c != nc) continue;
      for (c = 0; c < nc; c++)
        if (! (mpz_add_ui(t, low, i+cl[c]), num_lucas++, _GMP_is_lucas_pseudoprime(t,2)) ) break;
      if (c != nc) continue;
      PUSH_VLIST(retlist, ibase + i);
    }
    ibase += ppr;
    mpz_add_ui(low, low, ppr);
  }

  if (_verbose) printf("cluster sieve ran %"UVuf" MR and %"UVuf" Lucas tests (pretests %s)\n", num_mr, num_lucas, run_pretests ? "on" : "off");
  mpz_set(low, savelow);
  Safefree(cres);
  Safefree(VPrem);
  Safefree(resmod_0);
  Safefree(resmod_1);
  Safefree(resmod_2);
  Safefree(residues);
  mpz_clear(savelow);
  mpz_clear(t);
  *rn = retlist.nsize;
  return retlist.list;
}

static uint32_t* _todigits32(uint32_t *ndigits, uint32_t n, uint32_t base) {
  uint32_t bits[32];
  uint32_t *digits, i, d;

  for (d = 0; n; n /= base)
    bits[d++] = n % base;

  New(0, digits, d, uint32_t);
  for (i = 0; i < d; i++)
    digits[i] = bits[d-i-1];
  *ndigits = d;
  return digits;
}
static uint32_t* _todigits_base2(uint32_t *ndigits, mpz_t n) {
  uint32_t *digits, d, bits = mpz_sizeinbase(n, 2);
  New(0, digits, bits, uint32_t);
  for (d = 0; d < bits; d++)
    digits[d] = mpz_tstbit(n,bits-d-1);
  *ndigits = d;
  return digits;
}
/* Algorithm 1.26 FastIntegerOutput from MCA v0.5.9 */
uint32_t* todigits(uint32_t *ndigits, mpz_t n, uint32_t base) {
  uint32_t* digits, *rdigits, *qdigits;
  uint32_t  k, nlen, rlen, qlen, zlen, i, j;
  mpz_t Q, R;

  if (base == 2)
    return _todigits_base2(ndigits, n);

  if (mpz_cmp_ui(n, 0xFFFFFFFFUL) <= 0)
    return _todigits32(ndigits, mpz_get_ui(n), base);

  mpz_init(Q); mpz_init(R);
  nlen = logint(n, base) + 1;

  /* Find k s.t.  B^(2k-2) <= n <= B^(2k) */
  k = ((nlen-1) >> 1) + 1;

  /* Compute Q,R = DivRem(n, base^k) */
  mpz_ui_pow_ui(Q, base, k);
  mpz_tdiv_qr(Q, R, n, Q);

  /* In theory we could do this all in place, avoiding the allocations
   * and copying. */

  qdigits = todigits(&qlen, Q, base);
  rdigits = todigits(&rlen, R, base);

  zlen = k - rlen;
  if (qlen + rlen + zlen != nlen) croak("todigits: internal sizes wrong!");
  mpz_clear(Q);  mpz_clear(R);

  /* Allocate exactly as much space as needed. */
  New(0, digits, nlen, uint32_t);
  j = 0;

  for (i = 0; i < qlen; i++)
    digits[j++] = qdigits[i];
  for (i = 0; i < zlen; i++)
    digits[j++] = 0;
  for (i = 0; i < rlen; i++)
    digits[j++] = rdigits[i];

  Safefree(rdigits);
  Safefree(qdigits);

  *ndigits = nlen;
  return digits;
}

/* Algorithm 1.25 FastIntegerInput from MCA v0.5.9 */
void fromdigits(mpz_t n, uint32_t *d, uint32_t len, uint32_t base)
{
  if (len == 0) {
    mpz_set_ui(n, 0);
  } else if (len == 1) {
    mpz_set_ui(n, d[0]);
  } else {
    mpz_t *l, b;
    uint32_t i, k = len;
    mpz_init_set_ui(b, base);
    New(0, l, len, mpz_t);
    /* We expect input in lowest -> highest order, i.e. right to left */
    for (i = 0; i < len; i++)
      mpz_init_set_ui(l[i], d[i]);
    while (k > 1) {
      for (i = 0; i < k-1; i += 2) {
        /* l[i>>1] = l[i] + b * l[i+1] */
        mpz_mul(l[i+1], l[i+1], b);
        mpz_add(l[i>>1], l[i], l[i+1]);
      }
      if (k&1)
        mpz_set(l[k>>1], l[k-1]);
      k = (k+1)>>1;
      mpz_mul(b, b, b);
    }
    mpz_clear(b);
    mpz_set(n, l[0]);
    for (i = 0; i < len; i++)
      mpz_clear(l[i]);
    Safefree(l);
  }
}

void fromdigits_str(mpz_t n, const char* s, uint32_t base)
{
  uint32_t i, len, *dig;

  if (s[0] == '-' || s[0] == '+') s++;
  while (s[0] == '0') s++;

  len = strlen(s);
  New(0, dig, len, uint32_t);
  for (i = 0; i < len; i++) {
    const unsigned char c = s[i];
    uint32_t d = !isalnum(c) ? 255 : (c <= '9') ? c-'0' : (c <= 'Z') ? c-'A'+10 : c-'a'+10;
    if (d >= base) croak("Invalid digit for base %u", base);
    dig[len-1-i] = d;
  }
  fromdigits(n, dig, len, base);
  Safefree(dig);
}
