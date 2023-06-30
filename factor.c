#include <gmp.h>
#include "ptypes.h"

#include "factor.h"
#include "primality.h"
#include "prime_iterator.h"
#include "utility.h"
#include "pbrent63.h"
#include "squfof126.h"
#include "ecm.h"
#include "tinyqs.h"
#include "simpqs.h"
#include "lucas_seq.h"

#define _GMP_ECM_FACTOR(n, f, b1, ncurves) \
   _GMP_ecm_factor_projective(n, f, b1, 0, ncurves)

#define NPRIMES_SMALL 3450
static unsigned short primes_small[NPRIMES_SMALL];
static mpz_t _gcd_1k;
static mpz_t _gcd_4k;
static mpz_t _gcd_16k;
static mpz_t _gcd_32k;
void _init_factor(void) {
  uint32_t pn;
  PRIME_ITERATOR(iter);
  primes_small[0] = 0;
  primes_small[1] = 2;
  mpz_init_set_ui(_gcd_1k, 1);
  mpz_init_set_ui(_gcd_4k, 1);
  mpz_init_set_ui(_gcd_16k, 1);
  mpz_init_set_ui(_gcd_32k, 1);
  for (pn = 2; pn < NPRIMES_SMALL; pn++) {
    unsigned long p = prime_iterator_next(&iter);
    primes_small[pn] = p;
    if (p >     2 && p <=  1000)  mpz_mul_ui(_gcd_1k, _gcd_1k, p);
    if (p >  1000 && p <=  4000)  mpz_mul_ui(_gcd_4k, _gcd_4k, p);
    if (p >  4000 && p <= 16000)  mpz_mul_ui(_gcd_16k, _gcd_16k, p);
    if (p > 16000 && p <= 32000)  mpz_mul_ui(_gcd_32k, _gcd_32k, p);
  }
  prime_iterator_destroy(&iter);
}

/* Max number of factors on the unfactored stack, not the max total factors.
 * This is used when we split n into two or more composites.  Since we work
 * on the smaller of the composites first, this rarely goes above 10 even
 * with thousands of non-trivial factors. */
#define MAX_FACTORS 128

static int add_factor(int nfactors, mpz_t f, int e, mpz_t** pfactors, int** pexponents)
{
  int i, j, cmp = 0;
  MPUassert(e >= 1, "Adding factor with 0 exponent");
  if (nfactors == 0) {                      /* First factor */
    mpz_t *factors;
    int* exponents;
    New(0, factors, 10, mpz_t);
    New(0, exponents, 10, int);
    mpz_init_set(factors[0], f);
    exponents[0] = e;
    *pfactors = factors;
    *pexponents = exponents;
    return 1;
  } else if (mpz_cmp((*pfactors)[nfactors-1],f) < 0) {  /* New biggest factor */
    if (!(nfactors % 10)) {
      Renew(*pfactors, nfactors+10, mpz_t);
      Renew(*pexponents, nfactors+10, int);
    }
    mpz_init_set((*pfactors)[nfactors], f);
    (*pexponents)[nfactors] = e;
    return nfactors+1;
  }
  /* Insert in sorted order.  Find out where we will put it. */
  for (i = 0; i < nfactors; i++)
    if ((cmp = mpz_cmp((*pfactors)[i], f)) >= 0)
      break;
  if (cmp == 0) {                           /* Duplicate factor */
    (*pexponents)[i] += e;
    return nfactors;
  }
  /* factor[i] > f.  Move everything from i to nfactors up. */
  if (!(nfactors % 10)) {
    Renew(*pfactors, nfactors+10, mpz_t);
    Renew(*pexponents, nfactors+10, int);
  }
  mpz_init((*pfactors)[nfactors]);
  for (j = nfactors; j > i; j--) {
    mpz_set( (*pfactors)[j], (*pfactors)[j-1] );
    (*pexponents)[j] = (*pexponents)[j-1];
  }
  mpz_set((*pfactors)[i], f);
  (*pexponents)[i] = e;
  return nfactors+1;
}

#define ADD_FACTOR_UI(f, t) \
  do { \
    mpz_set_uv(f, t); \
    nfactors = add_factor(nfactors, f, 1, &factors, &exponents); \
  } while (0)

#define ADD_FACTOR(f) \
  do { nfactors = add_factor(nfactors, f, 1, &factors, &exponents); } while (0)

#define ADD_FACTORS(f, e) \
  do { nfactors = add_factor(nfactors, f, e, &factors, &exponents); } while (0)

#define TRIAL_DIVIDE_SMALL(n, pn_lo, pn_hi) \
  { unsigned long sp, p; \
    for (sp = pn_lo, p = primes_small[sp]; \
         sp <= pn_hi && mpz_cmp_ui(n,p*p) >= 0; \
         p = primes_small[++sp]) { \
      if (mpz_divisible_ui_p(n, p)) { \
        mpz_set_ui(f, p); \
        ADD_FACTORS(f, mpz_remove(n,n,f)); \
      } \
    } \
  }

int factor(mpz_t input_n, mpz_t* pfactors[], int* pexponents[])
{
  mpz_t tofac_stack[MAX_FACTORS];
  int ntofac = 0;
  mpz_t* factors = 0;
  int* exponents = 0;
  int nfactors = 0;
  mpz_t f, n;
  UV tf, tlim;

  mpz_init_set(n, input_n);
  mpz_init(f);
  if (mpz_cmp_ui(n, 4) < 0) {
    if (mpz_cmp_ui(n, 1) != 0)    /* 1 should return no results */
      ADD_FACTOR(n);
    goto DONE;
  }

  /* Trial factor to small limit */
  if (mpz_even_p(n)) {
    mpz_set_ui(f,2);
    ADD_FACTORS(f, mpz_remove(n,n,f));
  }
  /* Using gcd to detect any factors in the range, then remove them if found */
  if (mpz_gcd(f,n,_gcd_1k), mpz_cmp_ui(f,1)>0) TRIAL_DIVIDE_SMALL(n,   2,  168);
  if (mpz_gcd(f,n,_gcd_4k), mpz_cmp_ui(f,1)>0) TRIAL_DIVIDE_SMALL(n, 169,  550);
  if (mpz_gcd(f,n,_gcd_16k),mpz_cmp_ui(f,1)>0) TRIAL_DIVIDE_SMALL(n, 551, 1862);
  if (mpz_gcd(f,n,_gcd_32k),mpz_cmp_ui(f,1)>0) TRIAL_DIVIDE_SMALL(n,1863, 3432);

  tlim = 32003;
  if (mpz_cmp_ui(n,tlim*tlim) < 0) {
    if (mpz_cmp_ui(n,1) > 0)
      ADD_FACTOR(n);
    goto DONE;
  }

  /* Power factor */
  tf = power_factor(n, f);
  if (tf) {
    mpz_t* pow_factors;
    int* pow_exponents;
    int pow_nfactors, i, j;

    pow_nfactors = factor(f, &pow_factors, &pow_exponents);
    for (i = 0; i < pow_nfactors; i++)
      pow_exponents[i] *= tf;
    for (i = 0; i < pow_nfactors; i++)
      for (j = 0; j < pow_exponents[i]; j++)
        ADD_FACTOR(pow_factors[i]);
    clear_factors(pow_nfactors, &pow_factors, &pow_exponents);
    goto DONE;
  }

  do { /* loop over each remaining factor */
    while ( mpz_cmp_ui(n, tlim*tlim) > 0 && !_GMP_is_prob_prime(n) ) {
      int success = 0;
      int o = get_verbose_level();
      UV B1 = 5000;
      UV nbits = mpz_sizeinbase(n, 2);

      /*
       * This set of operations is meant to provide good performance for
       * "random" numbers as input.  Hence we stack lots of effort up front
       * looking for small factors: prho and pbrent are ~ O(f^1/2) where
       * f is the smallest factor.  SQUFOF is O(N^1/4), so arguable not
       * any better.  p-1 and ECM are quite useful for pulling out small
       * factors (6-20 digits).
       *
       * Factoring a 778-digit number consisting of 101 8-digit factors
       * should complete in under 3 seconds.  Factoring numbers consisting
       * of many 12-digit or 14-digit primes should take under 10 seconds.
       */

      /* Handle small inputs here */
      if (nbits <= 63) {
        if (!success) success = pbrent63(n, f, 400000);
        if (success&&o) {gmp_printf("UV Rho-Brent found factor %Zd\n", f);o=0;}
      }
      if (nbits >= 65 && nbits <= 126) {
        if (!success) success = _GMP_pminus1_factor(n, f, 5000, 5000);
        if (success&&o) {gmp_printf("p-1 (%dk) found factor %Zd\n",5,f);o=0;}
        if (!success) success = tinyqs(n, f);
        if (success&&o) {gmp_printf("tinyqs found factor %Zd\n", f);o=0;}
      }

      /* It's possible the previous calls failed or weren't available */
      if (nbits <= 53) {
        if (!success)  success = squfof126(n, f, 400000);
        if (success&&o) {gmp_printf("UV SQUFOF126 found factor %Zd\n", f);o=0;}
      } else if (nbits <= 77) {
        int sb1 = (nbits < 58) ?  1
                : (nbits < 63) ?  2
                : (nbits < 72) ?  4
                               : 10;
        if (!success)  success = _GMP_pminus1_factor(n, f, sb1*1000, sb1*10000);
        if (success&&o) {gmp_printf("p-1 (%dk) found factor %Zd\n",sb1,f);o=0;}

        if (!success)  success = squfof126(n, f, 1000000);
        if (success&&o) {gmp_printf("SQUFOF126 found factor %Zd\n", f);o=0;}
      }

      /* Make sure it isn't a perfect power */
      if (!success)  success = (int)power_factor(n, f);
      if (success&&o) {gmp_printf("perfect power found factor %Zd\n", f);o=0;}

      if (!success)  success = _GMP_pminus1_factor(n, f, 20000, 200000);
      if (success&&o) {gmp_printf("p-1 (20k) found factor %Zd\n", f);o=0;}

      /* Small ECM to find small factors */
      if (!success)  success = _GMP_ECM_FACTOR(n, f, 200, 4);
      if (success&&o) {gmp_printf("tiny ecm (200) found factor %Zd\n", f);o=0;}
      if (!success)  success = _GMP_ECM_FACTOR(n, f, 600, 20);
      if (success&&o) {gmp_printf("tiny ecm (600) found factor %Zd\n", f);o=0;}
      if (!success)  success = _GMP_ECM_FACTOR(n, f, 2000, 10);
      if (success&&o) {gmp_printf("tiny ecm (2000) found factor %Zd\n", f);o=0;}

      /* Small p-1 */
      if (!success) {
        if (nbits < 100 || nbits >= 160) {
          success = _GMP_pminus1_factor(n, f, 200000, 3000000);
          if (success&&o) {gmp_printf("p-1 (200k) found factor %Zd\n", f);o=0;}
        }
      }

      /* Set ECM parameters that have a good chance of success */
      if (!success) {
        UV curves;
        if      (nbits < 100){ B1 =   5000; curves =  20; }
        else if (nbits < 128){ B1 =  10000; curves =   2; } /* go to QS */
        else if (nbits < 160){ B1 =  20000; curves =   2; } /* go to QS */
        else if (nbits < 192){ B1 =  30000; curves =  20; }
        else if (nbits < 224){ B1 =  40000; curves =  40; }
        else if (nbits < 256){ B1 =  80000; curves =  40; }
        else if (nbits < 512){ B1 = 160000; curves =  80; }
        else                 { B1 = 320000; curves = 160; }
        if (curves > 0) {
          success = _GMP_ECM_FACTOR(n, f, B1, curves);
          if (success&&o) {gmp_printf("small ecm (%luk,%lu) found factor %Zd\n", B1/1000, curves, f);o=0;}
        }
      }

      /* QS (30+ digits).  Fantastic if it is a semiprime, but can be
       * slow and a memory hog if not (compared to ECM).  Restrict to
       * reasonable size numbers (< 91 digits).  Because of the way it
       * works, it will generate (possibly) multiple factors for the same
       * amount of work.  Go to some trouble to use them. */
      if (!success && mpz_sizeinbase(n,10) >= 30 && nbits < 300) {
        mpz_t farray[66];
        int i, qs_nfactors;
        for (i = 0; i < 66; i++)
          mpz_init(farray[i]);
        qs_nfactors = _GMP_simpqs(n, farray);
        mpz_set(f, farray[0]);
        if (qs_nfactors > 2) {
          /* We found multiple factors */
          for (i = 2; i < qs_nfactors; i++) {
            if (o){gmp_printf("SIMPQS found extra factor %Zd\n",farray[i]);}
            if (ntofac >= MAX_FACTORS-1) croak("Too many factors\n");
            mpz_init_set(tofac_stack[ntofac], farray[i]);
            ntofac++;
            mpz_divexact(n, n, farray[i]);
          }
          /* f = farray[0], n = farray[1], farray[2..] pushed */
        }
        for (i = 0; i < 66; i++)
          mpz_clear(farray[i]);
        success = qs_nfactors > 1;
        if (success&&o) {gmp_printf("SIMPQS found factor %Zd\n", f);o=0;}
      }

      if (!success)  success = _GMP_ECM_FACTOR(n, f, 2*B1, 20);
      if (success&&o) {gmp_printf("ecm (%luk,20) found factor %Zd\n",2*B1/1000,f);o=0;}

      if (!success)  success = _GMP_pbrent_factor(n, f, 1, 1*1024*1024);
      if (success&&o) {gmp_printf("pbrent (1,1M) found factor %Zd\n", f);o=0;}

      if (!success)  success = _GMP_ECM_FACTOR(n, f, 4*B1, 20);
      if (success&&o) {gmp_printf("ecm (%luk,20) ecm found factor %Zd\n", 4*B1,f);o=0;}

      if (!success)  success = _GMP_ECM_FACTOR(n, f, 8*B1, 20);
      if (success&&o) {gmp_printf("ecm (%luk,20) ecm found factor %Zd\n", 8*B1,f);o=0;}

      /* HOLF in case it's a near-ratio-of-perfect-square */
      if (!success)  success = _GMP_holf_factor(n, f, 1*1024*1024);
      if (success&&o) {gmp_printf("holf found factor %Zd\n", f);o=0;}

      /* Chebyshev poly in case p+1 is smooth */
      if (!success) success = _GMP_cheb_factor(n, f, 50000, 0);
      if (success&&o) {gmp_printf("cheb (%dk) found factor %Zd\n",50,f);o=0;}

      /* Large p-1 with stage 2: B2 = 20*B1 */
      if (!success)  success = _GMP_pminus1_factor(n,f,5000000,5000000*20);
      if (success&&o) {gmp_printf("p-1 (5M) found factor %Zd\n", f);o=0;}

      if (!success)  success = _GMP_ECM_FACTOR(n, f, 32*B1, 40);
      if (success&&o) {gmp_printf("ecm (%luk,40) ecm found factor %Zd\n", 32*B1,f);o=0;}

      /*
      if (!success)  success = _GMP_pbrent_factor(n, f, 2, 512*1024*1024);
      if (success&&o) {gmp_printf("pbrent (2,512M) found factor %Zd\n", f);o=0;}
      */

      /* Our method of last resort: ECM with high bmax and many curves*/
      if (!success) {
        int i;
        if (get_verbose_level()) gmp_printf("starting large ECM on %Zd\n",n);
        B1 *= 8;
        for (i = 0; i < 10; i++) {
          success = _GMP_ECM_FACTOR(n, f, B1, 100);
          if (success) break;
          B1 *= 2;
        }
        if (success&&o) {gmp_printf("ecm (%luk,100) ecm found factor %Zd\n", B1,f);o=0;}
      }

      if (success) {
        if (!mpz_divisible_p(n, f) || !mpz_cmp_ui(f, 1) || !mpz_cmp(f, n)) {
          gmp_printf("n = %Zd  f = %Zd\n", n, f);
          croak("Incorrect factoring");
        }
      }

      if (!success) {
        /* TODO: What to do with composites we can't factor?
         *       Push them as "C#####" ?
         *       For now, just push them as if we factored.
         */
        if (get_verbose_level()) gmp_printf("gave up on %Zd\n", n);
        ADD_FACTOR(n);
        mpz_set_ui(n, 1);
      } else {
        int ndiv = mpz_remove(n, n, f);
        if (_GMP_is_prob_prime(f)) { /* prime factor */
          while (ndiv-- > 0)
            ADD_FACTOR(f);
        } else if (ndiv > 1) {       /* Repeated non-trivial composite factor */
          mpz_t* pow_factors;
          int* pow_exponents;
          int pow_nfactors, i, j;
          pow_nfactors = factor(f, &pow_factors, &pow_exponents);
          for (i = 0; i < pow_nfactors; i++)
            pow_exponents[i] *= ndiv;
          for (i = 0; i < pow_nfactors; i++)
            for (j = 0; j < pow_exponents[i]; j++)
              ADD_FACTOR(pow_factors[i]);
          clear_factors(pow_nfactors, &pow_factors, &pow_exponents);
        } else {                     /* One non-trivial composite factor */
          if (ntofac >= MAX_FACTORS-1) croak("Too many factors\n");
          /* If f < n and both are composites, put n on stack and work on f */
          if (mpz_cmp(f, n) < 0 && !_GMP_is_prob_prime(n)) {
            mpz_init_set(tofac_stack[ntofac++], n);
            mpz_set(n, f);
          } else {
            mpz_init_set(tofac_stack[ntofac++], f);
          }
        }
      }
    }
    /* n is now prime or 1 */
    if (mpz_cmp_ui(n, 1) > 0) {
      ADD_FACTOR(n);
      mpz_set_ui(n, 1);
    }
    if (ntofac-- > 0) {
      mpz_set(n, tofac_stack[ntofac]);
      mpz_clear(tofac_stack[ntofac]);
    }
  } while (mpz_cmp_ui(n, 1) > 0);

DONE:
  mpz_clear(f);
  mpz_clear(n);
  *pfactors = factors;
  *pexponents = exponents;
  return nfactors;
}

void clear_factors(int nfactors, mpz_t* pfactors[], int* pexponents[])
{
  while (nfactors > 0)
    mpz_clear((*pfactors)[--nfactors]);
  Safefree(*pfactors);
  Safefree(*pexponents);
}


/*****************************************************************************/


/* TODO: Optimize to factor down to the last semiprime. */
static void _omega(mpz_t n, uint32_t* omega, uint32_t* bigomega)
{
  mpz_t* factors;
  int i, nfactors, result, *exponents;

  mpz_abs(n,n);
  nfactors = factor(n, &factors, &exponents);
  if (bigomega != 0) {
    for (i = 0, result = 0; i < nfactors; i++)
      result += exponents[i];
    *bigomega = result;
  }
  if (omega != 0)
    *omega = nfactors;
  clear_factors(nfactors, &factors, &exponents);
}

uint32_t omega(mpz_t n)
{
  uint32_t o;
  _omega(n, &o, 0);
  return o;
}
uint32_t bigomega(mpz_t n)
{
  uint32_t bo;
  _omega(n, 0, &bo);
  return bo;
}

void sigma(mpz_t res, mpz_t n, unsigned long k)
{
  mpz_t* factors;
  mpz_t pk, pke, fmult;
  int* exponents;
  int i, j, nfactors;

  if (mpz_cmp_ui(n, 1) <= 0) {
    mpz_set_ui(res, (k == 0 && mpz_cmp_ui(n,1) < 0) ? 2 : 1);
    return;
  }

  if (_GMP_is_prob_prime(n)) {
    mpz_pow_ui(res, n, k);
    mpz_add_ui(res, res, 1);
    return;
  }

  nfactors = factor(n, &factors, &exponents);
  if (k == 0) {
    for (i = 0; i < nfactors; i++) {
      mpz_set_ui(factors[i], exponents[i]+1);
    }
  } else if (k == 1) {
    mpz_init(pke);
    mpz_init(fmult);
    for (i = 0; i < nfactors; i++) {
      mpz_set(pke, factors[i]);
      mpz_add_ui(fmult, factors[i], 1);
      for (j = 1; j < exponents[i]; j++) {
        mpz_mul(pke, pke, factors[i]);
        mpz_add(fmult, fmult, pke);
      }
      mpz_set(factors[i], fmult);
    }
    mpz_clear(fmult);
    mpz_clear(pke);
  } else {
    mpz_init(pk);
    mpz_init(pke);
    mpz_init(fmult);
    for (i = 0; i < nfactors; i++) {
      mpz_pow_ui(pk, factors[i], k);
      mpz_add_ui(fmult, pk, 1);
      mpz_set(pke, pk);
      for (j = 1; j < exponents[i]; j++) {
        mpz_mul(pke, pke, pk);
        mpz_add(fmult, fmult, pke);
      }
      mpz_set(factors[i], fmult);
    }
    mpz_clear(fmult);
    mpz_clear(pke);
    mpz_clear(pk);
  }
  mpz_product(factors, 0, nfactors-1);
  mpz_set(res, factors[0]);
  clear_factors(nfactors, &factors, &exponents);
}


static const unsigned long smalldiv[] = {4, 9, 25, 49, 121, 169, 289};
int moebius(mpz_t n)
{
  uint32_t i, o, bo;

  if (mpz_sgn(n) < 0) {
    mpz_neg(n,n);
    o = moebius(n);
    mpz_neg(n,n);
    return o;
  }
  if (mpz_sgn(n) == 0) return 0;
  if (mpz_cmp_ui(n, 1) == 0) return 1;

  for (i = 0; i < 7; i++)
    if (mpz_divisible_ui_p(n, smalldiv[i]))
      return 0;

  _omega(n, &o, &bo);
  return (o != bo) ?   0
      :  (o % 2)   ?  -1
                   :   1;
}

int liouville(mpz_t n)
{
  uint32_t result = bigomega(n);
  return (result & 1)  ?  -1  : 1;
}

void totient(mpz_t tot, mpz_t n_input)
{
  mpz_t t, n;
  mpz_t* factors;
  int* exponents;
  int i, j, nfactors;

  if (mpz_cmp_ui(n_input, 1) <= 0) {
    mpz_set(tot, n_input);
    return;
  }
  mpz_init_set(n, n_input);
  mpz_set_ui(tot, 1);
  /* Fast reduction of multiples of 2 */
  i = mpz_scan1(n, 0);
  if (i > 0) {
    if (i > 1)  mpz_mul_2exp(tot, tot, i-1);
    mpz_tdiv_q_2exp(n, n, i);
  }
  /* Now factor and calculate totient */
  nfactors = factor(n, &factors, &exponents);
  mpz_init(t);
  for (i = 0; i < nfactors; i++) {
    mpz_sub_ui(t, factors[i], 1);
    for (j = 1; j < exponents[i]; j++)
      mpz_mul(t, t, factors[i]);
    mpz_mul(tot, tot, t);
  }
  mpz_clear(t);
  clear_factors(nfactors, &factors, &exponents);
  mpz_clear(n);
}

void jordan_totient(mpz_t tot, mpz_t n, unsigned long k)
{
  if (k == 0) {
    mpz_set_ui(tot, (mpz_cmp_ui(n, 1) == 0) ? 1 : 0);
  } else if (k == 1) {
    totient(tot, n);
  } else if (mpz_cmp_ui(n, 1) <= 0) {
    mpz_set_ui(tot, (mpz_cmp_ui(n, 1) == 0) ? 1 : 0);
  } else {
    mpz_t t;
    mpz_t* factors;
    int* exponents;
    int i, j, nfactors;

    nfactors = factor(n, &factors, &exponents);
    mpz_init(t);
    for (i = 0; i < nfactors; i++) {
      mpz_pow_ui(t, factors[i], k);
      mpz_sub_ui(factors[i], t, 1);
      for (j = 1; j < exponents[i]; j++)
        mpz_mul(factors[i], factors[i], t);
    }
    mpz_product(factors, 0, nfactors-1);
    mpz_set(tot, factors[0]);
    mpz_clear(t);
    clear_factors(nfactors, &factors, &exponents);
  }
}

void carmichael_lambda(mpz_t lambda, mpz_t n)
{
  if (mpz_cmp_ui(n, 8) < 0) {
    totient(lambda, n);
  } else if (mpz_scan1(n, 0) == mpz_sizeinbase(n, 2)-1) {
    mpz_tdiv_q_2exp(lambda, n, 2);
  } else {
    mpz_t t;
    mpz_t* factors;
    int* exponents;
    int i, j, nfactors;

    nfactors = factor(n, &factors, &exponents);
    mpz_init(t);
    mpz_set_ui(lambda, 1);
    if (exponents[0] > 2 && mpz_cmp_ui(factors[0], 2) == 0) exponents[0]--;
    for (i = 0; i < nfactors; i++) {
      mpz_sub_ui(t, factors[i], 1);
      for (j = 1; j < exponents[i]; j++)
        mpz_mul(t, t, factors[i]);
      mpz_lcm(lambda, lambda, t);
    }
    mpz_clear(t);
    clear_factors(nfactors, &factors, &exponents);
  }
}

static const unsigned char _zntiny[7][7] = {
  {2},
  {0,2},
  {4,4,2},
  {0,0,0,2},
  {3,6,3,6,2},
  {0,2,0,2,0,2},
  {6,0,3,6,0,3,2} };
static const unsigned char _zn11[ 9] = {10,5,5,5,10,10,10,5,2};
static const unsigned char _zn13[11] = {12,3,6,4,12,12,4,3,6,12,2};
static const unsigned char _zn16[14] = {0,4,0,4,0,2,0,2,0,4,0,4,0,2};
static const unsigned char _zn17[15] = {8,16,4,16,16,16,8,8,16,16,16,4,16,8,2};

/* Do not alias any of the arguments */
static void _znorder1(mpz_t order, mpz_t a, mpz_t p, int e, mpz_t t, mpz_t n)
{
  mpz_t phi, *factors;
  int* exponents;
  int i, j, nfactors;

  mpz_set_ui(order, 1);
  mpz_pow_ui(n, p, e);

  /* Remove some simple cases */
  if (mpz_cmp_ui(n, 1) <= 0)  { mpz_set(order, n); return; }
  mpz_mod(t, a, n);
  if (mpz_cmp_ui(t, 1) <= 0)  { mpz_set(order, t); return; }
  /* These are purely for performance. */
  if (mpz_cmp_ui(n,17) <= 0) {
    unsigned long int nn = mpz_get_ui(n), aa = mpz_get_ui(t);
    if (nn <=  9)  { mpz_set_ui(order, _zntiny[nn-3][aa-2]); return; }
    if (nn == 11)  { mpz_set_ui(order, _zn11[aa-2]); return; }
    if (nn == 13)  { mpz_set_ui(order, _zn13[aa-2]); return; }
    if (nn == 16)  { mpz_set_ui(order, _zn16[aa-2]); return; }
    if (nn == 17)  { mpz_set_ui(order, _zn17[aa-2]); return; }
  }

  /* Abhijit Das, algorithm 1.7 */
  /* This could be further simplified / optimized */
  mpz_init(phi);
  mpz_sub_ui(phi, p, 1);
  nfactors = factor(phi, &factors, &exponents);
  if (e > 1) {
    mpz_pow_ui(t, p, e-1);
    mpz_mul(phi, phi, t);
    ADD_FACTORS(p, e-1);
  }
  for (i = 0; i < nfactors; i++) {
    mpz_divexact(t, phi, factors[i]);
    for (j = 1; j < exponents[i]; j++)
      mpz_divexact(t, t, factors[i]);
    mpz_powm(t, a, t, n);
    for (j = 0;  mpz_cmp_ui(t, 1) != 0;  mpz_powm(t, t, factors[i], n)) {
      if (j++ >= exponents[i]) {
        mpz_set_ui(order, 0);
        break;
      }
      mpz_mul(order, order, factors[i]);
    }
    if (j > exponents[i]) break;
  }
  mpz_clear(phi);
  clear_factors(nfactors, &factors, &exponents);
}

void znorder(mpz_t res, mpz_t a, mpz_t n)
{
  mpz_t t, order;

  /* TODO: Usually we don't want to modify their inputs */
  mpz_abs(n,n);
  if (mpz_cmp_ui(n, 1) <= 0) { mpz_set(res, n); return; }
  mpz_mod(a, a, n);
  if (mpz_cmp_ui(a, 1) <= 0) { mpz_set(res, a); return; }

  mpz_init(t);
  mpz_gcd(t, a, n);
  if (mpz_cmp_ui(t, 1) != 0) {
    mpz_set_ui(res, 0);
    mpz_clear(t);
    return;
  }
  mpz_init_set_ui(order, 1);

  { /* Factor n, then lcm all the znorder(p^e). */
    mpz_t order1, t2, *factors;
    int* exponents;
    int i, nfactors;

    mpz_init(order1);
    mpz_init(t2);
    nfactors = factor(n, &factors, &exponents);
    for (i = 0; i < nfactors; i++) {
      _znorder1(order1, a, factors[i], exponents[i], t, t2);
      mpz_lcm(order, order, order1);
    }
    clear_factors(nfactors, &factors, &exponents);
    mpz_clear(t2);
    mpz_clear(order1);
  }
  mpz_set(res, order);
  mpz_clear(order);
  mpz_clear(t);
}

void znprimroot(mpz_t root, mpz_t n)
{
  mpz_t t, phi, a, on, r, *factors;
  int i, nfactors, *exponents, oddprime, k;

  mpz_set_ui(root, 0);
  if (mpz_cmp_ui(n, 4) <= 0) {
    if (mpz_sgn(n) > 0)
      mpz_sub_ui(root, n, 1);
    return;
  }
  if (mpz_divisible_ui_p(n, 4))
    return;

  mpz_init(r);
  mpz_init_set(on, n);
  if (mpz_even_p(on))
    mpz_tdiv_q_2exp(on, on, 1);
  if (!power_factor(on, r))
    mpz_set(r, on);
  if (!_GMP_is_prob_prime(r)) {
    mpz_clear(on);
    mpz_clear(r);
    return;
  }

  mpz_init(phi);
  mpz_sub_ui(phi, r, 1);
  mpz_divexact(on, on, r);
  mpz_mul(phi, phi, on);

  mpz_sub_ui(r,n,1);
  oddprime = (mpz_cmp(r,phi) == 0);

  mpz_clear(on);
  mpz_clear(r);

  mpz_init(t);
  mpz_init(a);
  nfactors = factor(phi, &factors, &exponents);

  /* Replace each factor with phi/factor */
  for (i = 0; i < nfactors; i++)
    mpz_divexact(factors[i], phi, factors[i]);

  for (mpz_set_ui(a,2);  mpz_cmp(a,n) < 0;  mpz_add_ui(a,a,1)) {
    if (!mpz_cmp_ui(a,4) || !mpz_cmp_ui(a,8) || !mpz_cmp_ui(a,9)) continue;
    k = mpz_kronecker(a, n);
    if ( (oddprime && k != -1) || (!oddprime && k == 0) ) continue;
    for (i = 0; i < nfactors; i++) {
      mpz_powm(t, a, factors[i], n);
      if (mpz_cmp_ui(t, 1) == 0)
        break;
    }
    if (i == nfactors) { mpz_set(root, a); break; }
  }
  clear_factors(nfactors, &factors, &exponents);
  mpz_clear(a);
  mpz_clear(t);
  mpz_clear(phi);
}

static const int32_t tau_table[] = {
  0,1,-24,252,-1472,4830,-6048,-16744,84480,-113643,-115920,534612,-370944,-577738,401856,1217160,987136,-6905934,2727432,10661420,-7109760,-4219488,-12830688,18643272,21288960,-25499225,13865712,-73279080,24647168,128406630,-29211840,-52843168,-196706304,134722224,165742416,-80873520,167282496,-182213314,-255874080,-145589976,408038400,308120442,101267712,-17125708,-786948864,-548895690,-447438528
};
#define NTAU (sizeof(tau_table)/sizeof(tau_table[0]))
void ramanujan_tau(mpz_t res, mpz_t n)
{
  mpz_t t, t1, t2, t3, t4, *factors;
  int i, nfactors, *exponents;
  unsigned long j, p2;

  if (mpz_cmp_ui(n, NTAU) < 0) {
    if (mpz_sgn(n) <= 0) mpz_set_si(res, 0);
    else                 mpz_set_si(res, tau_table[mpz_get_ui(n)]);
    return;
  }

  /* We are doing far too much work here for sigma5.  We could do it just
   * for primes then use the multiplicative property.  However that works
   * for prime *powers*, so it isn't quite so simple.  This solution also
   * gets to be high memory use. */

  /* Pari/GP does this using Hurwitz class numbers.  That is a more
   * complicated but far more efficient solution. */

  mpz_init(t); mpz_init(t1); mpz_init(t2); mpz_init(t3); mpz_init(t4);
  nfactors = factor(n, &factors, &exponents);
  for (i = 0; i < nfactors; i++) {
    /* t = tau(p) */
    if (mpz_cmp_ui(factors[i], NTAU) < 0) {
      mpz_set_si(t, tau_table[mpz_get_ui(factors[i])]);
    } else {
      mpz_pow_ui(t, factors[i], 11);   mpz_add_ui(t, t, 1); /* sigma(t,f,11) */
      mpz_mul_ui(t1, t, 65);
      mpz_pow_ui(t, factors[i],  5);   mpz_add_ui(t, t, 1); /* sigma(t,f, 5) */
      mpz_mul_ui(t2, t, 691);
      mpz_add(t1, t1, t2);

      /* t1 in use. t2 accumulate.  t3, t4, t free. */
      mpz_sub_ui(t, factors[i], 1);
      mpz_tdiv_q_2exp(t, t, 1);
      p2 = mpz_get_ui(t);
      mpz_set_ui(t2, 0);
      for (j = 1; j <= p2; j++) {
        mpz_set_uv(t, j);
        sigma(t3, t, 5);
        mpz_sub_ui(t, factors[i], j);
        sigma(t, t, 5);
        mpz_mul(t4, t3, t);
        mpz_add(t2, t2, t4);
      }
      mpz_mul_ui(t2, t2, 2*691*252);
      mpz_sub(t, t1, t2);
      mpz_tdiv_q_ui(t, t, 756);
    }
    /* t holds tau(p), all other temps are free. */

    if (exponents[i] > 1) {
      mpz_pow_ui(t1, t, exponents[i]);
      if (exponents[i] == 2) {
        mpz_pow_ui(t2, factors[i], 11);
      } else if (exponents[i] == 3) {
        mpz_pow_ui(t2, factors[i], 11);
        mpz_mul(t2, t2, t);
        mpz_mul_ui(t2, t2, 2);
      } else {
        /* t1 = t^e  t2 = neg sum,  t3 = prod,  t4 = temp */
        mpz_set_ui(t2, 0);
        for (j = 1; j <= ((unsigned long)exponents[i])>>1; j++) {
          mpz_set_si(t3, (j&1) ? -1 : 1);
          mpz_pow_ui(t4, factors[i], 11*j);
          mpz_mul(t3, t3, t4);
          mpz_bin_uiui(t4, exponents[i]-j, exponents[i]-2*j);
          mpz_mul(t3, t3, t4);
          mpz_pow_ui(t4, t, exponents[i]-2*j);
          mpz_mul(t3, t3, t4);
          mpz_sub(t2, t2, t3);
        }
      }
      mpz_sub(t, t1, t2);
    }
    mpz_set(factors[i], t);
  }
  mpz_product(factors, 0, nfactors-1);
  mpz_set(res, factors[0]);
  clear_factors(nfactors, &factors, &exponents);
  mpz_clear(t4); mpz_clear(t3); mpz_clear(t2); mpz_clear(t1); mpz_clear(t);
}

int is_semiprime(mpz_t n)
{
  int ret;
  unsigned long k, div, lim = 6000;
  mpz_t t;

  if (mpz_cmp_ui(n,6) < 0)
    return (mpz_cmp_ui(n,4) == 0);

  mpz_init(t);

  div = _GMP_trial_factor(n, 2, lim);
  if (div > 0) {
    mpz_divexact_ui(t, n, div);
    ret = _GMP_is_prime(t);
    mpz_clear(t);
    return !!ret;
  }
  /* No small divisors */
  if (_GMP_BPSW(n))
    { mpz_clear(t); return 0; }
  if (mpz_ui_pow_ui(t,lim,3), mpz_cmp(n,t) < 0)
    { mpz_clear(t); return 1; }

  /* Quick check for power */
  k = power_factor(n, t);
  if (k >= 2) {
    ret = (k == 2 && _GMP_is_prime(t));
    mpz_clear(t);
    return ret;
  }

  /* Number is composite, isn't tiny, and has no small divisors */
  if (    0
       || _GMP_pbrent_factor(n, t, 1, 15000)
       || _GMP_pminus1_factor(n, t, 50000, 500000)
       || _GMP_ECM_FACTOR(n, t,     800, 10)
       || _GMP_ECM_FACTOR(n, t,    8000, 20)
       || _GMP_ECM_FACTOR(n, t,   80000, 40)
       || _GMP_ECM_FACTOR(n, t,  320000, 40)
       || _GMP_ECM_FACTOR(n, t, 1000000, 80)
     ) {
    ret = _GMP_BPSW(t);
    if (ret) {
      mpz_divexact(t, n, t);
      ret = _GMP_BPSW(t);
    }
    mpz_clear(t);
    return !!ret;
  }
  /* No luck finding a small factor.  Do it the hard way. */
  {
    mpz_t* factors;
    int* exponents;
    int nfactors, i, j;

    nfactors = factor(n, &factors, &exponents);
    for (i = 0, j = 0; i < nfactors; i++)
      j += exponents[i];
    clear_factors(nfactors, &factors, &exponents);
    mpz_clear(t);
    return (j == 2);
  }
}


/*****************************************************************************/

UV _GMP_trial_factor(mpz_t n, UV from_n, UV to_n)
{
  size_t log2n = mpz_sizeinbase(n, 2);
  UV p = 0;
  PRIME_ITERATOR(iter);

  if (mpz_cmp_ui(n, 6) < 0) {
    unsigned long un = mpz_get_ui(n);
    if (un == 1) p = 1;
    else if (un == 4 && from_n <= 2 && to_n >= 2) p = 2;
    prime_iterator_destroy(&iter);
    return p;
  }
  if      (from_n <= 2 && to_n >= 2 && mpz_even_p(n))             p = 2;
  else if (from_n <= 3 && to_n >= 3 && mpz_divisible_ui_p(n, 3))  p = 3;
  else if (from_n <= 5 && to_n >= 5 && mpz_divisible_ui_p(n, 5))  p = 5;
  if (p != 0) {
    prime_iterator_destroy(&iter);
    return p;
  }

  if (from_n < 7)
    from_n = 7;
  if (from_n > to_n)
    { prime_iterator_destroy(&iter);  return 0; }
  /* p will be the next prime >= from_n */
  prime_iterator_setprime(&iter, from_n-1);
  p = prime_iterator_next(&iter);

  /* All native math if n fits in an unsigned long */
  if (log2n <= sizeof(unsigned long)*8) {
    unsigned long un = mpz_get_ui(n);
    unsigned long sqrtn = (unsigned long) sqrt((double)un);
    /* Be extra careful here, as we are using unsigned long, which may not
     * match a UV.  But GMP's ui is 'unsigned long' so that's what we have
     * to deal with.  We want to make sure we get the correct integer sqrt,
     * but also watch out for overflow. */
    while (sqrtn*sqrtn > un) sqrtn--;
    while ( (sqrtn+1)*(sqrtn+1) <= un
            && sqrtn < (1UL << 4*sizeof(unsigned long)) )
      sqrtn++;
    if (to_n > sqrtn)
      to_n = sqrtn;
    while (p <= to_n) {
      if ((un % p) == 0)
        break;
      p = prime_iterator_next(&iter);
    }
    prime_iterator_destroy(&iter);
    return (p <= to_n) ? p : 0;
  }

  /* For "small" numbers, this simple method is best. */
  {
    UV small_to = (log2n < 3000)  ?  to_n  :  30000;
    while (p <= small_to) {
      if (mpz_divisible_ui_p(n, p))
        break;
      p = prime_iterator_next(&iter);
    }
    if (p <= small_to || p > to_n) {
      prime_iterator_destroy(&iter);
      return (p <= small_to) ? p : 0;
    }
  }

  /* Simple treesieve.
   * This is much faster than simple divisibility for really big numbers.
   * Credit to Jens K Andersen for writing up the generic algorithm.
   *
   * This will search until the first group element is > to_n, which means
   * we will search a bit farther than to_n.
   */
  {
    UV found = 0;
    unsigned long* xn;        /* leaves */
    mpz_t* xtree[16+1];       /* the tree (maxdepth = 16) */
    mpz_t* xtemp;
    unsigned int i, j, d, depth, leafsize, nleaves;

    /* Decide on the tree depth (3-16) and number of leaves (10-31) */
    {
      unsigned int dp = log2n >> 10;
      depth = 0;
      while (dp >>= 1) depth++;
      if (depth < 3) depth = 3;
      if (depth > 16) depth = 16;
    }
    leafsize = log2n / (1U << depth) / 68;
    nleaves = 1 << depth;
    /* printf("log2n %lu  depth %u  leafsize %u  nleaves %u\n",log2n,depth,leafsize,nleaves); */

    New(0, xn, nleaves * leafsize, unsigned long);
    for (d = 0; d <= depth; d++) {
      unsigned int nodes = 1 << (depth - d);
      New(0, xtree[d], nodes, mpz_t);
      for (j = 0; j < nodes; j++)
        mpz_init(xtree[d][j]);
    }
    xtemp = xtree[1];   /* implies mindepth = 3 */

    while (!found && p <= to_n) {
      /* Create nleaves x[0] values, each the product of leafsize primes */
      for (i = 0; i < nleaves; i++) {
        for (j = 0; j < 4; j++)                  /* Create 4 sub-products */
          mpz_set_ui(xtemp[j], 1);
        for (j = 0; j < leafsize; j++) {
          xn[i*leafsize+j] = p;
          mpz_mul_ui(xtemp[j&3], xtemp[j&3], p);
          p = prime_iterator_next(&iter);
        }
        mpz_mul(xtemp[0], xtemp[0], xtemp[1]);   /* Combine for final product*/
        mpz_mul(xtemp[2], xtemp[2], xtemp[3]);
        mpz_mul(xtree[0][i], xtemp[0], xtemp[2]);
      }
      /* Multiply product tree, xtree[depth][0] has nleaves*leafsize product */
      for (d = 1; d <= depth; d++)
        for (i = 0; i < (1U << (depth-d)); i++)
          mpz_mul(xtree[d][i], xtree[d-1][2*i], xtree[d-1][2*i+1]);
      /* Go backwards replacing the products with remainders */
      mpz_tdiv_r(xtree[depth][0], n, xtree[depth][0]);
      for (d = 1; d <= depth; d++)
        for (i = 0; i < (1U << d); i++)
          mpz_tdiv_r(xtree[depth-d][i], xtree[depth-d+1][i>>1], xtree[depth-d][i]);
      /* Search each leaf for divisors */
      for (i = 0; !found && i < nleaves; i++)
        for (j = 0; j < leafsize; j++)
          if (mpz_divisible_ui_p(xtree[0][i], xn[i*leafsize+j]))
            { found = xn[i*leafsize+j]; break; }
    }
    p = found;
    for (d = 0; d <= depth; d++) {
      unsigned int nodes = 1U << (depth - d);
      for (j = 0; j < nodes; j++)
        mpz_clear(xtree[d][j]);
      Safefree(xtree[d]);
    }
    Safefree(xn);
    if (p > 0 && !mpz_divisible_ui_p(n, p))
      croak("incorrect trial factor\n");
  }
  prime_iterator_destroy(&iter);
  return p;
}


#define TEST_FOR_2357(n, f) \
  { \
    if (mpz_divisible_ui_p(n, 2)) { mpz_set_ui(f, 2); return 1; } \
    if (mpz_divisible_ui_p(n, 3)) { mpz_set_ui(f, 3); return 1; } \
    if (mpz_divisible_ui_p(n, 5)) { mpz_set_ui(f, 5); return 1; } \
    if (mpz_divisible_ui_p(n, 7)) { mpz_set_ui(f, 7); return 1; } \
    if (mpz_cmp_ui(n, 121) < 0) { return 0; } \
  }

int _GMP_prho_factor(mpz_t n, mpz_t f, UV a, UV rounds)
{
  mpz_t U, V, oldU, oldV, m;
  int i;
  const UV inner = 256;

  TEST_FOR_2357(n, f);
  rounds = (rounds + inner - 1) / inner;
  mpz_init_set_ui(U, 7);
  mpz_init_set_ui(V, 7);
  mpz_init(m);
  mpz_init(oldU);
  mpz_init(oldV);
  while (rounds-- > 0) {
    mpz_set_ui(m, 1); mpz_set(oldU, U);  mpz_set(oldV, V);
    for (i = 0; i < (int)inner; i++) {
      mpz_mul(U, U, U);  mpz_add_ui(U, U, a);  mpz_tdiv_r(U, U, n);
      mpz_mul(V, V, V);  mpz_add_ui(V, V, a);  mpz_tdiv_r(V, V, n);
      mpz_mul(V, V, V);  mpz_add_ui(V, V, a);  mpz_tdiv_r(V, V, n);
      mpz_sub(f, U, V);
      mpz_mul(m, m, f);
      if ((i%4) == ((inner-1)%4)) mpz_tdiv_r(m, m, n);
    }
    mpz_gcd(f, m, n);
    if (!mpz_cmp_ui(f, 1))
      continue;
    if (!mpz_cmp(f, n)) {
      /* f == n, so we have to back up to see what factor got found */
      mpz_set(U, oldU); mpz_set(V, oldV);
      i = inner;
      do {
        mpz_mul(U, U, U);  mpz_add_ui(U, U, a);  mpz_tdiv_r(U, U, n);
        mpz_mul(V, V, V);  mpz_add_ui(V, V, a);  mpz_tdiv_r(V, V, n);
        mpz_mul(V, V, V);  mpz_add_ui(V, V, a);  mpz_tdiv_r(V, V, n);
        if (mpz_cmp(U, V) >= 0)  mpz_sub(f, U, V);
        else                     mpz_sub(f, V, U);
        mpz_gcd(f, f, n);
      } while (!mpz_cmp_ui(f, 1) && i-- != 0);
      if ( (!mpz_cmp_ui(f, 1)) || (!mpz_cmp(f, n)) )  break;
    }
    mpz_clear(U); mpz_clear(V); mpz_clear(m); mpz_clear(oldU); mpz_clear(oldV);
    return 1;
  }
  mpz_clear(U); mpz_clear(V); mpz_clear(m); mpz_clear(oldU); mpz_clear(oldV);
  mpz_set(f, n);
  return 0;
}

int _GMP_pbrent_factor(mpz_t n, mpz_t f, UV a, UV rounds)
{
  mpz_t Xi, Xm, saveXi, m, t;
  UV i, r;
  const UV inner = 256;

  TEST_FOR_2357(n, f);
  mpz_init_set_ui(Xi, 2);
  mpz_init_set_ui(Xm, 2);
  mpz_init(m);
  mpz_init(t);
  mpz_init(saveXi);

  r = 1;
  mpz_set_ui(f,1);
  while (rounds > 0) {
    UV rleft = (r > rounds) ? rounds : r;
    while (rleft > 0) {   /* Do rleft rounds, inner at a time */
      UV dorounds = (rleft > inner) ? inner : rleft;
      mpz_set_ui(m, 1);
      mpz_set(saveXi, Xi);
      for (i = 0; i < dorounds; i++) {
        mpz_mul(t, Xi, Xi);  mpz_add_ui(t, t, a);  mpz_tdiv_r(Xi, t, n);
        mpz_sub(f, Xm, Xi);
        mpz_mul(m, m, f);
        if ((i%4) == ((dorounds-1)%4)) mpz_tdiv_r(m, m, n);
      }
      rleft -= dorounds;
      rounds -= dorounds;
      mpz_gcd(f, m, n);
      if (mpz_cmp_ui(f, 1) != 0)
        break;
    }
    if (!mpz_cmp_ui(f, 1)) {
      r *= 2;
      mpz_set(Xm, Xi);
      continue;
    }
    if (!mpz_cmp(f, n)) {
      /* f == n, so we have to back up to see what factor got found */
      mpz_set(Xi, saveXi);
      do {
        mpz_mul(t, Xi, Xi);  mpz_add_ui(t, t, a);  mpz_tdiv_r(Xi, t, n);
        mpz_sub(f, Xm, Xi); if (mpz_sgn(f) < 0) mpz_add(f,f,n);
        mpz_gcd(f, f, n);
      } while (!mpz_cmp_ui(f, 1) && r-- != 0);
    }
    break;
  }
  mpz_clear(Xi); mpz_clear(Xm); mpz_clear(m); mpz_clear(saveXi); mpz_clear(t);
  if (!mpz_cmp_ui(f, 1) || !mpz_cmp(f, n)) {
    mpz_set(f, n);
    return 0;
  }
  return 1;
}

/* References for P-1:
 *  Montgomery 1987:  https://cr.yp.to/bib/1987/montgomery.pdf
 *  Brent 1990:       http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.127.4316
 *
 * The main advantage of this over ECM is that it is *much* lower overhead,
 * so very cheap to run with relatively small B1,B2 values.  A disadvantage
 * is no continuation method, so subsequent calls with larger B1,B2 will
 * repeat all the previous work.  ECM is much better for harder factorisations,
 * so we typically want to try a little p-1 then move to ECM (or QS).
 */
int _GMP_pminus1_factor(mpz_t n, mpz_t f, UV B1, UV B2)
{
  mpz_t a, savea, t;
  UV q, saveq, j, sqrtB1;
  int _verbose = get_verbose_level();
  PRIME_ITERATOR(iter);

  TEST_FOR_2357(n, f);
  if (B1 < 7) return 0;

  mpz_init(a);
  mpz_init(savea);
  mpz_init(t);

  if (_verbose>2) gmp_printf("# p-1 trying %Zd (B1=%"UVuf" B2=%"UVuf")\n", n, (unsigned long)B1, (unsigned long)B2);

  /* STAGE 1
   * See Montgomery 1987 p249-250 or Brent 1990 p5.  We can take E to be the
   * lcm of integers to B1, then gcd(a^E-1,n) may be a factor of n.  While
   * we could actually calculate the LCM, it is quite inefficient to do so.
   * There are various ways to speed it up, but generally we prefer to do it
   * the way Brent indicates, which is one powmod for each prime p below B1,
   * a = a^(p^e) mod n, where e is the largest e such that p^e <= B1.
   *
   * One thing that can speed things up quite a bit is not running the GCD
   * on every step.  However with small factors this means we can easily end
   * up with multiple factors between GCDs, so we allow backtracking.  This
   * could also be added to stage 2, but it's far less likely to happen there.
   */
  j = 15;
  mpz_set_ui(a, 2);
  mpz_set_ui(savea, 2);
  saveq = 2;
  /* We could wrap this in a loop trying a few different a values, in case
   * the current one ended up going to 0. */
  q = 2;
  mpz_set_ui(t, 1);
  sqrtB1 = (UV) sqrt(B1);
  while (q <= B1) {
    UV k = q;
    if (q <= sqrtB1) {
      UV kmin = B1/q;
      while (k <= kmin)
        k *= q;
    }
    mpz_mul_ui(t, t, k);        /* Accumulate powers for a */
    if ( (j++ % 32) == 0) {
      mpz_powm(a, a, t, n);     /* a=a^(k1*k2*k3*...) mod n */
      if (mpz_sgn(a))  mpz_sub_ui(t, a, 1);
      else             mpz_sub_ui(t, n, 1);
      mpz_gcd(f, t, n);         /* f = gcd(a-1, n) */
      mpz_set_ui(t, 1);
      if (mpz_cmp(f, n) == 0)
        break;
      if (mpz_cmp_ui(f, 1) != 0)
        goto end_success;
      saveq = q;
      mpz_set(savea, a);
    }
    q = prime_iterator_next(&iter);
  }
  mpz_powm(a, a, t, n);
  if (mpz_sgn(a))  mpz_sub_ui(t, a, 1);
  else             mpz_sub_ui(t, n, 1);
  mpz_gcd(f, t, n);
  if (mpz_cmp(f, n) == 0) {
    /* We found multiple factors.  Loop one at a time. */
    prime_iterator_setprime(&iter, saveq);
    mpz_set(a, savea);
    for (q = saveq; q <= B1; q = prime_iterator_next(&iter)) {
      UV k = q;
      if (q <= sqrtB1) {
        UV kmin = B1/q;
        while (k <= kmin)
          k *= q;
      }
      mpz_powm_ui(a, a, k, n );
      mpz_sub_ui(t, a, 1);
      mpz_gcd(f, t, n);
      if (mpz_cmp(f, n) == 0)
        goto end_fail;
      if (mpz_cmp_ui(f, 1) != 0)
        goto end_success;
    }
  }
  if ( (mpz_cmp_ui(f, 1) != 0) && (mpz_cmp(f, n) != 0) )
    goto end_success;

  /* STAGE 2
   * This is the standard continuation which replaces the powmods in stage 1
   * with two mulmods, with a GCD every 64 primes (no backtracking).
   * This is quite a bit faster than stage 1.
   * See Montgomery 1987, p250-253 for possible optimizations.
   * We quickly precalculate a few of the prime gaps, and lazily cache others
   * up to a gap of 222.  That's enough for a B2 value of 189 million.  We
   * still work above that, we just won't cache the value for big gaps.
   */
  if (B2 > B1) {
    mpz_t b, bm, bmdiff;
    mpz_t precomp_bm[111];
    int   is_precomp[111] = {0};
    UV* primes = 0;
    UV sp = 1;

    mpz_init(bmdiff);
    mpz_init_set(bm, a);
    mpz_init_set_ui(b, 1);

    /* Set the first 20 differences */
    mpz_powm_ui(bmdiff, bm, 2, n);
    mpz_init_set(precomp_bm[0], bmdiff);
    is_precomp[0] = 1;
    for (j = 1; j < 22; j++) {
      mpz_mul(bmdiff, bmdiff, bm);
      mpz_mul(bmdiff, bmdiff, bm);
      mpz_tdiv_r(bmdiff, bmdiff, n);
      mpz_init_set(precomp_bm[j], bmdiff);
      is_precomp[j] = 1;
    }

    mpz_powm_ui(a, a, q, n );
    if (B2 < 10000000) {
      /* grab all the primes at once.  Hack around non-perfect iterator. */
      primes = sieve_to_n(B2+300, 0);
      for (sp = B1>>4; primes[sp] <= q; sp++)  ;
      /* q is primes <= B1, primes[sp] is the next prime */
    }

    j = 31;
    while (q <= B2) {
      UV lastq, qdiff;

      lastq = q;
      q = primes ? primes[sp++] : prime_iterator_next(&iter);
      qdiff = (q - lastq) / 2 - 1;

      if (qdiff < 111 && is_precomp[qdiff]) {
        mpz_mul(t, a, precomp_bm[qdiff]);
      } else if (qdiff < 111) {
        mpz_powm_ui(bmdiff, bm, q-lastq, n);
        mpz_init_set(precomp_bm[qdiff], bmdiff);
        is_precomp[qdiff] = 1;
        mpz_mul(t, a, bmdiff);
      } else {
        mpz_powm_ui(bmdiff, bm, q-lastq, n);  /* Big gap */
        mpz_mul(t, a, bmdiff);
      }
      mpz_tdiv_r(a, t, n);
      if (mpz_sgn(a))  mpz_sub_ui(t, a, 1);
      else             mpz_sub_ui(t, n, 1);
      mpz_mul(b, b, t);
      if ((j % 2) == 0)           /* put off mods a little */
        mpz_tdiv_r(b, b, n);
      if ( (j++ % 64) == 0) {     /* GCD every so often */
        mpz_gcd(f, b, n);
        if ( (mpz_cmp_ui(f, 1) != 0) && (mpz_cmp(f, n) != 0) )
          break;
      }
    }
    mpz_gcd(f, b, n);
    mpz_clear(b);
    mpz_clear(bm);
    mpz_clear(bmdiff);
    for (j = 0; j < 111; j++) {
      if (is_precomp[j])
        mpz_clear(precomp_bm[j]);
    }
    if (primes != 0) Safefree(primes);
    if ( (mpz_cmp_ui(f, 1) != 0) && (mpz_cmp(f, n) != 0) )
      goto end_success;
  }

  end_fail:
    mpz_set(f,n);
  end_success:
    prime_iterator_destroy(&iter);
    mpz_clear(a);
    mpz_clear(savea);
    mpz_clear(t);
    if ( (mpz_cmp_ui(f, 1) != 0) && (mpz_cmp(f, n) != 0) ) {
      if (_verbose>2) gmp_printf("# p-1: %Zd\n", f);
      return 1;
    }
    if (_verbose>2) gmp_printf("# p-1: no factor\n");
    mpz_set(f, n);
    return 0;
}

static void pp1_pow(mpz_t X, mpz_t Y, unsigned long exp, mpz_t n)
{
  mpz_t x0;
  unsigned long bit;
  {
    unsigned long v = exp;
    unsigned long b = 1;
    while (v >>= 1) b++;
    bit = 1UL << (b-2);
  }
  mpz_init_set(x0, X);
  mpz_mul(Y, X, X);
  mpz_sub_ui(Y, Y, 2);
  mpz_tdiv_r(Y, Y, n);
  while (bit) {
    if ( exp & bit ) {
      mpz_mul(X, X, Y);
      mpz_sub(X, X, x0);
      mpz_mul(Y, Y, Y);
      mpz_sub_ui(Y, Y, 2);
    } else {
      mpz_mul(Y, X, Y);
      mpz_sub(Y, Y, x0);
      mpz_mul(X, X, X);
      mpz_sub_ui(X, X, 2);
    }
    mpz_mod(X, X, n);
    mpz_mod(Y, Y, n);
    bit >>= 1;
  }
  mpz_clear(x0);
}

int _GMP_pplus1_factor(mpz_t n, mpz_t f, UV P0, UV B1, UV B2)
{
  UV j, q, saveq, sqrtB1;
  mpz_t X, Y, saveX;
  PRIME_ITERATOR(iter);

  TEST_FOR_2357(n, f);
  if (B1 < 7) return 0;

  mpz_init_set_ui(X, P0);
  mpz_init(Y);
  mpz_init(saveX);

  /* Montgomery 1987 */
  if (P0 == 0) {
    mpz_set_ui(X, 7);
    if (mpz_invert(X, X, n)) {
      mpz_mul_ui(X, X, 2);
      mpz_mod(X, X, n);
    } else
      P0 = 1;
  }
  if (P0 == 1) {
    mpz_set_ui(X, 5);
    if (mpz_invert(X, X, n)) {
      mpz_mul_ui(X, X, 6);
      mpz_mod(X, X, n);
    } else
      P0 = 2;
  }
  if (P0 == 2) {
    mpz_set_ui(X, 11);
    if (mpz_invert(X, X, n)) {
      mpz_mul_ui(X, X, 23);
      mpz_mod(X, X, n);
    }
  }

  sqrtB1 = (UV) sqrt(B1);
  j = 8;
  q = 2;
  saveq = q;
  mpz_set(saveX, X);
  while (q <= B1) {
    UV k = q;
    if (q <= sqrtB1) {
      UV kmin = B1/q;
      while (k <= kmin)
        k *= q;
    }
    pp1_pow(X, Y, k, n);
    if ( (j++ % 16) == 0) {
      mpz_sub_ui(f, X, 2);
      if (mpz_sgn(f) == 0)        break;
      mpz_gcd(f, f, n);
      if (mpz_cmp(f, n) == 0)     break;
      if (mpz_cmp_ui(f, 1) > 0)   goto end_success;
      saveq = q;
      mpz_set(saveX, X);
    }
    q = prime_iterator_next(&iter);
  }
  mpz_sub_ui(f, X, 2);
  mpz_gcd(f, f, n);
  if (mpz_cmp_ui(X, 2) == 0 || mpz_cmp(f, n) == 0) {
    /* Backtrack */
    prime_iterator_setprime(&iter, saveq);
    mpz_set(X, saveX);
    for (q = saveq; q <= B1; q = prime_iterator_next(&iter)) {
      UV k = q;
      if (q <= sqrtB1) {
        UV kmin = B1/q;
        while (k <= kmin)
          k *= q;
      }
      pp1_pow(X, Y, k, n);
      mpz_sub_ui(f, X, 2);
      if (mpz_sgn(f) == 0)        goto end_fail;
      mpz_gcd(f, f, n);
      if (mpz_cmp(f, n) == 0)     break;
      if (mpz_cmp_ui(f, 1) > 0)   goto end_success;
    }
  }
  if ( (mpz_cmp_ui(f, 1) > 0) && (mpz_cmp(f, n) != 0) )
    goto end_success;
  /* TODO: stage 2 */
  end_fail:
    mpz_set(f,n);
  end_success:
    prime_iterator_destroy(&iter);
    mpz_clear(X);  mpz_clear(Y);  mpz_clear(saveX);
    return (mpz_cmp_ui(f, 1) != 0) && (mpz_cmp(f, n) != 0);
}

int _GMP_cheb_factor(mpz_t n, mpz_t f, UV B, UV initx)
{
  unsigned long p;
  double logB;
  mpz_t x, inv, t, k, P, Q;
  PRIME_ITERATOR(iter);

  if (B == 0) { B = mpz_sizeinbase(n,2);  B = B*B*B; }

  TEST_FOR_2357(n, f);
  if (B < 7) return 0;

  logB = logl(B);
  mpz_init_set_ui(inv, 2);
  mpz_invert(inv, inv, n);   /* multiplying by this will divide by two */
  mpz_init_set_ui(x, (initx == 0) ? 72 : initx);
  mpz_init(t);
  mpz_init(k);
  mpz_init(P);
  mpz_init_set_ui(Q, 1);

  mpz_set_ui(f, 1);
  for (p = 2; p <= B && mpz_cmp_ui(f,1) <= 0; p = prime_iterator_next(&iter)) {
    unsigned long lgbp = (unsigned long) (logB / logl(p));   /* Alternately logint(mpzB,p) */
    if (lgbp > 1)  mpz_ui_pow_ui(k, p, lgbp);
    else           mpz_set_uv(k, p);
    mpz_mul_2exp(P, x, 1);
    lucasvmod(x, P, Q, k, n, t);
    mpz_mul(x, x, inv);
    mpz_mod(x, x, n);
    mpz_sub_ui(t, x, 1);
    mpz_gcd(f, t, n);
  }
  if (mpz_cmp_ui(f,1) <= 0)
    mpz_set(f,n);
  prime_iterator_destroy(&iter);
  mpz_clear(Q);  mpz_clear(P);  mpz_clear(k);  mpz_clear(t);
  mpz_clear(x);  mpz_clear(inv);
  return (mpz_cmp_ui(f, 1) > 0) && (mpz_cmp(f, n) < 0);
}

int _GMP_holf_factor(mpz_t n, mpz_t f, UV rounds)
{
  mpz_t s, m;
  UV i;

#define PREMULT 480   /* 1  2  6  12  480  151200 */

  TEST_FOR_2357(n, f);
  if (mpz_perfect_square_p(n)) {
    mpz_sqrt(f, n);
    return 1;
  }

  mpz_mul_ui(n, n, PREMULT);
  mpz_init(s);
  mpz_init(m);
  for (i = 1; i <= rounds; i++) {
    mpz_mul_ui(f, n, i);    /* f = n*i */
    if (mpz_perfect_square_p(f)) {
      /* s^2 = n*i, so m = s^2 mod n = 0.  Hence f = GCD(n, s) = GCD(n, n*i) */
      mpz_divexact_ui(n, n, PREMULT);
      mpz_gcd(f, f, n);
      mpz_clear(s); mpz_clear(m);
      if (mpz_cmp(f, n) == 0)  return 0;
      return 1;
    }
    mpz_sqrt(s, f);
    mpz_add_ui(s, s, 1);    /* s = ceil(sqrt(n*i)) */
    mpz_mul(m, s, s);
    mpz_sub(m, m, f);       /* m = s^2 mod n = s^2 - n*i */
    if (mpz_perfect_square_p(m)) {
      mpz_divexact_ui(n, n, PREMULT);
      mpz_sqrt(f, m);
      mpz_sub(s, s, f);
      mpz_gcd(f, s, n);
      mpz_clear(s); mpz_clear(m);
      return (mpz_cmp_ui(f, 1) > 0);
    }
  }
  mpz_divexact_ui(n, n, PREMULT);
  mpz_set(f, n);
  mpz_clear(s); mpz_clear(m);
  return 0;
}

/* See if n is a perfect power */
unsigned long power_factor(mpz_t n, mpz_t f)
{
  unsigned long k = 1, b = 2;
  if (mpz_cmp_ui(n, 1) > 0 && mpz_perfect_power_p(n)) {
    mpz_t nf, tf;
    PRIME_ITERATOR(iter);

    mpz_init_set(nf, n);
    mpz_init(tf);
    while (1) {
      unsigned long ok = k;
      while (mpz_root(tf, nf, b)) {
        mpz_set(f, tf);
        mpz_set(nf, tf);
        k *= b;
      }
      if (ok != k && !mpz_perfect_power_p(nf)) break;
      if (mpz_cmp_ui(tf, 1) <= 0) break; /* Exit if we can't find the power */
      b = prime_iterator_next(&iter);
    }
    mpz_clear(tf);  mpz_clear(nf);
    prime_iterator_destroy(&iter);
  }
  return (k == 1) ? 0 : k;
}

static int numcmp(const void *av, const void *bv)
  { return mpz_cmp(*(const mpz_t*)av, *(const mpz_t*)bv); }

mpz_t * divisor_list(int *num_divisors, mpz_t n)
{
  mpz_t *factors, *divs, mult;
  int nfactors, ndivisors, i, j, k, count, *exponents;

  nfactors = factor(n, &factors, &exponents);
  ndivisors = 1 + ((nfactors > 0) ? exponents[0] : 0);
  for (i = 1; i < nfactors; i++)
    ndivisors *= (exponents[i] + 1);

  mpz_init(mult);
  New(0, divs, ndivisors, mpz_t);
  mpz_init_set_ui(divs[0], 1);
  for (count = 1, k = 0; k < nfactors; k++) {
    int scount = count;
    mpz_set_ui(mult, 1);
    for (j = 0; j < exponents[k]; j++) {
      mpz_mul(mult, mult, factors[k]);
      for (i = 0; i < scount; i++) {
        mpz_init(divs[count]);
        mpz_mul(divs[count], divs[i], mult);
        count++;
      }
    }
  }
  mpz_clear(mult);
  clear_factors(nfactors, &factors, &exponents);

  qsort(divs, ndivisors, sizeof(mpz_t), numcmp);

  *num_divisors = ndivisors;
  return divs;
}

int is_smooth(mpz_t n, mpz_t k) {
  mpz_t *factors;
  mpz_t N;
  int i, nfactors, *exponents;
  uint32_t klo, khi, div;

  if (mpz_cmp_ui(n,1) <= 0) return 1;
  if (mpz_cmp_ui(k,1) <= 0) return 0;
  if (mpz_cmp(n, k) <= 0) return 1;

  mpz_init_set(N, n);
  klo = 2;
  khi = (mpz_cmp_ui(k, 10000000) >= 0) ? 10000000 : mpz_get_ui(k);

  while (klo <= khi && (div = _GMP_trial_factor(N, klo, khi)) > 0) {
    do {
      mpz_divexact_ui(N, N, div);
    } while (mpz_divisible_ui_p(N, div));
    if (mpz_cmp(N, k) <= 0) {
      mpz_clear(N);
      return 1;
    }
    klo = div+1;
  }
  /* N still has at least one factor, and no factors in N <= khi */
  if (mpz_cmp_ui(k, khi) <= 0) {
    mpz_clear(N);
    return 0;
  }

  nfactors = factor(N, &factors, &exponents);
  for (i = 0; i < nfactors; i++) {
    if (mpz_cmp(factors[i], k) > 0)
      break;
  }
  clear_factors(nfactors, &factors, &exponents);
  mpz_clear(N);
  return i >= nfactors;
}

int is_rough(mpz_t n, mpz_t k) {
  mpz_t *factors;
  mpz_t N, f;
  int i, nfactors, *exponents;
  uint32_t khi, stage;

  if (mpz_sgn(n) == 0) return 0 + (mpz_sgn(k) == 0);
  if (mpz_cmp_ui(n,1) == 0 || mpz_cmp_ui(k,1) <= 0) return 1;
  if (mpz_cmp_ui(k,2) == 0) return 0 + (mpz_cmp_ui(n,1) >= 0);
  if (mpz_cmp_ui(k,3) == 0) return (mpz_sgn(n) > 0 && mpz_odd_p(n));
  if (mpz_cmp(n, k) < 0) return 0;  /* 1 < n < k all are zero */
  /* k > 3;  n >= k;  n >= 2 */

  if (mpz_even_p(n) || mpz_divisible_ui_p(n,3)) return 0;
  if (mpz_cmp_ui(k,5) <= 0) return 1;
  if (mpz_divisible_ui_p(n,5)) return 0;
  /* k > 5;  n >= k;  n not divisible by 2, 3, or 5 */

  /* 1. Trial division up to a limit */
  khi = (mpz_cmp_ui(k, 1000000) >= 0) ? 1000000 : mpz_get_ui(k);

  if (_GMP_trial_factor(n, 2, khi-1) > 0)
    return 0;
  if (mpz_cmp_ui(k, khi) <= 0)
    return 1;

  /* 2. Try to efficiently pull out small factors */
  mpz_init_set(N, n);
  mpz_init(f);
  for (stage = 1; stage <= 7; stage++) {
    int success = 0;
    if (stage == 5 && _GMP_BPSW(N)) { mpz_clear(f); mpz_clear(N); return 1; }
    switch (stage) {
      case 1: success = _GMP_pminus1_factor(N, f,   13,   130); break;
      case 2: success = _GMP_pminus1_factor(N, f,  150,  1500); break;
      case 3: success = _GMP_pminus1_factor(N, f, 1500, 15000); break;
      case 4: success = _GMP_ECM_FACTOR(N, f,  100, 1); break;
      case 5: success = _GMP_ECM_FACTOR(N, f,  400, 1); break;
      case 6: success = _GMP_ECM_FACTOR(N, f,  800, 1); break;
      case 7: success = _GMP_ECM_FACTOR(N, f, 1600, 1); break;
      default: break;
    }
    if (!success) continue;
    nfactors = factor(f, &factors, &exponents);
    for (i = 0; i < nfactors; i++) {
      if (mpz_cmp(factors[i], k) < 0)
        break;
    }
    clear_factors(nfactors, &factors, &exponents);
    if (i < nfactors) break;       /* Return 0: Found a small factor. */
    mpz_divexact(N, N, f);         /* Divide out the factors all checked */
    if (mpz_cmp(N, k) < 0) break;  /* Return 0: Reduced N to lower than k */
    if (stage >= 5 && _GMP_BPSW(N)) { mpz_clear(f); mpz_clear(N); return 1; }
  }
  mpz_clear(f);
  if (stage <= 7) {
    mpz_clear(N);
    return 0;
  }

  /* 3. Fully factor */
  nfactors = factor(N, &factors, &exponents);
  for (i = 0; i < nfactors; i++) {
    if (mpz_cmp(factors[i], k) < 0)
      break;
  }
  clear_factors(nfactors, &factors, &exponents);
  mpz_clear(N);
  return i >= nfactors;
}

int is_powerful(mpz_t n, uint32_t k) {
  mpz_t *factors;
  mpz_t N, f;
  int i, nfactors, *exponents;
  uint32_t e, klo, khi, div;

  if (k == 0) k = 2;   /* API */

  if (k <= 1 || mpz_cmp_ui(n,1) <= 0) return 1;

  mpz_init_set(N, n);
  e = mpz_scan1(N, 0);
  if (e > 0) {
    if (e < k) { mpz_clear(N); return 0; }
    mpz_tdiv_q_2exp(N, N, e);
    if (mpz_cmp_ui(N,1) == 0) { mpz_clear(N); return 1; }
  }

  for (i = 2; i <= 8; i++) {
    static const unsigned char pr[9] = {0,2,3,5,7,11,13,17,19};
    uint32_t p = pr[i];
    uint32_t pk = (k == 2) ? p*p : (k == 3) ? p*p*p : p*p*p*p;
    if (mpz_divisible_ui_p(N, p) && !mpz_divisible_ui_p(N, pk))
      { mpz_clear(N); return 0; }
  }

  mpz_init(f);
  mpz_root(f, N, k);
  klo = 3;
  khi = (mpz_cmp_ui(f, 1000000) >= 0) ? 1000000 : mpz_get_ui(f);

  while (klo <= khi && (div = _GMP_trial_factor(N, klo, khi)) > 0) {
    mpz_set_ui(f, div);
    if (mpz_remove(N, N, f) < k)
      { mpz_clear(N); mpz_clear(f); return 0; }
    if (mpz_cmp_ui(N,1) == 0 || power_factor(N, f) >= k)
      { mpz_clear(N); mpz_clear(f); return 1; }
    mpz_ui_pow_ui(f, div, 2*k);
    if (mpz_cmp(N,f) < 0)
      { mpz_clear(N); mpz_clear(f); return 0; }
    klo = div+1;
  }
  if (mpz_cmp_ui(N,1) == 0 || power_factor(N, f) >= k)
    { mpz_clear(N); mpz_clear(f); return 1; }
  mpz_ui_pow_ui(f, khi, 2*k);
  if (mpz_cmp(N,f) < 0)
    { mpz_clear(N); mpz_clear(f); return 0; }
  mpz_clear(f);

  /* 3. Fully factor */
  nfactors = factor(N, &factors, &exponents);
  for (i = 0; i < nfactors; i++) {
    if ((uint32_t)exponents[i] < k)
      break;
  }
  clear_factors(nfactors, &factors, &exponents);
  mpz_clear(N);
  return i >= nfactors;
}

int is_almost_prime(uint32_t k, mpz_t n)
{
#if 0
  if (k == 0)  return (mpz_cmp_ui(n,1) == 0) ? 1 : 0;
  if (k == 1)  return _GMP_is_prime(n) ? 1 : 0;
  if (k == 2)  return is_semiprime(n);
#endif
  while (k > 0 && mpz_even_p(n))           { k--; mpz_divexact_ui(n,n,2); }
  while (k > 0 && mpz_divisible_ui_p(n,3)) { k--; mpz_divexact_ui(n,n,3); }
  while (k > 0 && mpz_divisible_ui_p(n,5)) { k--; mpz_divexact_ui(n,n,5); }
  while (k > 0 && mpz_divisible_ui_p(n,7)) { k--; mpz_divexact_ui(n,n,7); }
  if (k == 0)  return (mpz_cmp_ui(n,1) == 0) ? 1 : 0;
  if (k == 1)  return _GMP_is_prime(n) ? 1 : 0;
  if (k == 2)  return is_semiprime(n);

  /* Optimally, we should be factoring one at a time.  This will let us
   * exit early if we find too many factors, and stop when we get a final
   * semiprime.  This is all rather tedious. */

  return bigomega(n) == k;
}
