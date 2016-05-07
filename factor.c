#include <gmp.h>
#include "ptypes.h"

#include "factor.h"
#include "primality.h"
#include "prime_iterator.h"
#include "utility.h"
#include "small_factor.h"
#include "ecm.h"
#include "simpqs.h"

#define _GMP_ECM_FACTOR(n, f, b1, ncurves) \
   _GMP_ecm_factor_projective(n, f, b1, 0, ncurves)

#define NPRIMES_SMALL 2000
static unsigned short primes_small[NPRIMES_SMALL];
void _init_factor(void) {
  UV pn;
  PRIME_ITERATOR(iter);
  primes_small[0] = 0;
  primes_small[1] = 2;
  for (pn = 2; pn < NPRIMES_SMALL; pn++) {
    primes_small[pn] = prime_iterator_next(&iter);
  }
  prime_iterator_destroy(&iter);
}

/* Max number of factors on the unfactored stack, not the max total factors.
 * This is used when we split n into two or more composites.  Since we work
 * on the smaller of the composites first, this rarely goes above 10 even
 * with thousands of non-trivial factors. */
#define MAX_FACTORS 128

static int add_factor(int nfactors, mpz_t f, mpz_t** pfactors, int** pexponents)
{
  int i, j, cmp = 0;
  if (nfactors == 0) {                      /* First factor */
    mpz_t *factors;
    int* exponents;
    New(0, factors, 10, mpz_t);
    New(0, exponents, 10, int);
    mpz_init_set(factors[0], f);
    exponents[0] = 1;
    *pfactors = factors;
    *pexponents = exponents;
    return 1;
  } else if (mpz_cmp((*pfactors)[nfactors-1],f) < 0) {  /* New biggest factor */
    if (!(nfactors % 10)) {
      Renew(*pfactors, nfactors+10, mpz_t);
      Renew(*pexponents, nfactors+10, int);
    }
    mpz_init_set((*pfactors)[nfactors], f);
    (*pexponents)[nfactors] = 1;
    return nfactors+1;
  }
  /* Insert in sorted order.  Find out where we will put it. */
  for (i = 0; i < nfactors; i++)
    if ((cmp = mpz_cmp((*pfactors)[i], f)) >= 0)
      break;
  if (cmp == 0) {                           /* Duplicate factor */
    (*pexponents)[i]++;
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
  (*pexponents)[i] = 1;
  return nfactors+1;
}

#define ADD_FACTOR_UI(f, t) \
  do { \
    mpz_set_ui(f, t); \
    nfactors = add_factor(nfactors, f, &factors, &exponents); \
  } while (0)

#define ADD_FACTOR(f) \
  do { nfactors = add_factor(nfactors, f, &factors, &exponents); } while (0)

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
  while (mpz_even_p(n)) {
    ADD_FACTOR_UI(f, 2);
    mpz_divexact_ui(n, n, 2);
  }
  tlim = (mpz_sizeinbase(n,2) > 80)  ?  4001  :  16001;
  {
    UV sp, p, un;
    un = (mpz_cmp_ui(n,2*tlim*tlim) >= 0) ? 2*tlim*tlim : mpz_get_ui(n);

    for (sp = 2, p = primes_small[sp];
         p < tlim && p*p <= un;
         p = primes_small[++sp]) {
      while (mpz_divisible_ui_p(n, p)) {
        ADD_FACTOR_UI(f, p);
        mpz_divexact_ui(n, n, p);
        un = (mpz_cmp_ui(n,2*tlim*tlim) > 0) ? 2*tlim*tlim : mpz_get_ui(n);
      }
    }

    if (un < p*p) {
      if (un > 1)
        ADD_FACTOR(n);
      goto DONE;
    }
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
      UV nbits, B1 = 5000;

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

      if (mpz_cmp_ui(n, (unsigned long)(UV_MAX>>4)) < 0) {
        UV ui_n = mpz_get_ui(n);
        UV ui_factors[2];
        if (!mpz_cmp_ui(n, ui_n)) {
          success = racing_squfof_factor(ui_n, ui_factors, 200000)-1;
          if (success) {
            mpz_set_ui(f, ui_factors[0]);
          } else {
            if (o > 2) {gmp_printf("UV SQUFOF failed %Zd\n", n);}
          }
        }
        if (success&&o) {gmp_printf("UV SQUFOF found factor %Zd\n", f);o=0;}
      }

      /* Make sure it isn't a perfect power */
      if (!success)  success = (int)power_factor(n, f);
      if (success&&o) {gmp_printf("perfect power found factor %Zd\n", f);o=0;}

      if (!success)  success = _GMP_pminus1_factor(n, f, 15000, 150000);
      if (success&&o) {gmp_printf("p-1 (15k) found factor %Zd\n", f);o=0;}

      /* Small ECM to find small factors */
      if (!success)  success = _GMP_ECM_FACTOR(n, f, 200, 4);
      if (success&&o) {gmp_printf("tiny ecm (200) found factor %Zd\n", f);o=0;}
      if (!success)  success = _GMP_ECM_FACTOR(n, f, 600, 20);
      if (success&&o) {gmp_printf("tiny ecm (600) found factor %Zd\n", f);o=0;}
      if (!success)  success = _GMP_ECM_FACTOR(n, f, 2000, 10);
      if (success&&o) {gmp_printf("tiny ecm (2000) found factor %Zd\n", f);o=0;}

      /* Small p-1 */
      if (!success) {
        nbits = mpz_sizeinbase(n, 2);
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

      /* Large p-1 with stage 2: B2 = 20*B1 */
      if (!success)  success = _GMP_pminus1_factor(n,f,5000000,5000000*20);
      if (success&&o) {gmp_printf("p-1 (5M) found factor %Zd\n", f);o=0;}

      if (!success)  success = _GMP_ECM_FACTOR(n, f, 32*B1, 40);
      if (success&&o) {gmp_printf("ecm (%luk,40) ecm found factor %Zd\n", 32*B1,f);o=0;}

      /*
      if (!success)  success = _GMP_pbrent_factor(n, f, 2, 512*1024*1024);
      if (success&&o) {gmp_printf("pbrent (2,512M) found factor %Zd\n", f);o=0;}

      if (!success)  success = _GMP_squfof_factor(n, f, 256*1024*1024);
      if (success&&o) {gmp_printf("squfof found factor %Zd\n", f);o=0;}
      */

      /* Our method of last resort: ECM with high bmax and many curves*/
      if (!success) {
        UV i;
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


void sigma(mpz_t res, mpz_t n, UV k)
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
  mpz_t* factors;
  int* exponents;
  int i, nfactors, result;

  if (mpz_sgn(n) <= 0)       return 0;
  if (mpz_cmp_ui(n, 1) == 0) return 1;

  for (i = 0; i < 7; i++)
    if (mpz_divisible_ui_p(n, smalldiv[i]))
      return 0;

  nfactors = factor(n, &factors, &exponents);
  result = (nfactors % 2) ? -1 : 1;
  for (i = 0; i < nfactors; i++)
    if (exponents[i] > 1)
      { result = 0; break; }
  clear_factors(nfactors, &factors, &exponents);
  return result;
}

int liouville(mpz_t n)
{
  mpz_t* factors;
  int* exponents;
  int i, nfactors, result;

  nfactors = factor(n, &factors, &exponents);
  for (i = 0, result = 0; i < nfactors; i++)
    result += exponents[i];
  result = (result & 1) ? -1 : 1;
  clear_factors(nfactors, &factors, &exponents);
  return result;
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
    mpz_set_ui(tot, 1);
    for (i = 0; i < nfactors; i++) {
      mpz_pow_ui(t, factors[i], k);
      mpz_sub_ui(t, t, 1);
      mpz_mul(tot, tot, t);
      mpz_add_ui(t, t, 1);
      for (j = 1; j < exponents[i]; j++)
        mpz_mul(tot, tot, t);
    }
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

void znorder(mpz_t res, mpz_t a, mpz_t n)
{
  mpz_t t;

  mpz_init(t);
  mpz_gcd(t, a, n);

  if (mpz_cmp_ui(n, 1) <= 0) {
    mpz_set(res, n);
  } else if (mpz_cmp_ui(a, 1) <= 0) {
    mpz_set(res, a);
  } else if (mpz_cmp_ui(t, 1) != 0) {
    mpz_set_ui(res, 0);
  } else {
    mpz_t order, phi;
    mpz_t* factors;
    int* exponents;
    int i, j, nfactors;

    mpz_init_set_ui(order, 1);
    mpz_init(phi);
    /* Abhijit Das, algorithm 1.7, applied to Carmichael Lambda */
    carmichael_lambda(phi, n);
    nfactors = factor(phi, &factors, &exponents);
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
    mpz_set(res, order);
    mpz_clear(phi);
    mpz_clear(order);
    clear_factors(nfactors, &factors, &exponents);
  }
  mpz_clear(t);
}

void znprimroot(mpz_t root, mpz_t n)
{
  if (mpz_cmp_ui(n, 4) <= 0) {
    if (mpz_sgn(n) <= 0)  mpz_set_ui(root, 0);
    else                  mpz_sub_ui(root, n, 1);
  } else if (mpz_divisible_ui_p(n, 4)) {
    mpz_set_ui(root, 0);
  } else {
    mpz_t* factors;
    int* exponents;
    int i, nfactors;
    mpz_t t, phi, a;

    mpz_init(phi);  mpz_init(t);
    nfactors = 1;
    mpz_sub_ui(phi, n, 1);
    if (!_GMP_is_prob_prime(n)) {
      if (mpz_even_p(n)) mpz_tdiv_q_2exp(t, n, 1);
      else               mpz_set(t, n);
      nfactors = factor(t, &factors, &exponents);
      mpz_sub_ui(phi, factors[0], 1);
      for (i = 1; i < exponents[0]; i++)
        mpz_mul(phi, phi, factors[0]);
      clear_factors(nfactors, &factors, &exponents);
    }
    if (nfactors != 1) {
      mpz_set_ui(root, 0);
    } else {
      nfactors = factor(phi, &factors, &exponents);
      i = 0;
      for (mpz_init_set_ui(a,2);  mpz_cmp(a, n) < 0;  mpz_add_ui(a, a, 1)) {
        if (mpz_kronecker(a, n) == 0) continue;
        for (i = 0; i < nfactors; i++) {
          mpz_divexact(t, phi, factors[i]);
          mpz_powm(t, a, t, n);
          if (mpz_cmp_ui(t, 1) == 0)
            break;
        }
        if (i == nfactors) break;
      }
      if (i == nfactors)  mpz_set(root, a);
      else                mpz_set_ui(root, 0);
      mpz_clear(a);
      clear_factors(nfactors, &factors, &exponents);
    }
    mpz_clear(t);  mpz_clear(phi);
  }
}

static const int32_t tau_table[] = {
  0,1,-24,252,-1472,4830,-6048,-16744,84480,-113643,-115920,534612,-370944,-577738,401856,1217160,987136,-6905934,2727432,10661420,-7109760,-4219488,-12830688,18643272,21288960,-25499225,13865712,-73279080,24647168,128406630,-29211840,-52843168,-196706304,134722224,165742416,-80873520,167282496,-182213314,-255874080,-145589976,408038400,308120442,101267712,-17125708,-786948864,-548895690,-447438528
};
#define NTAU (sizeof(tau_table)/sizeof(tau_table[0]))
void ramanujan_tau(mpz_t res, mpz_t n)
{
  mpz_t t, t1, t2, t3, t4, *factors;
  int i, nfactors, *exponents;
  UV j, p2;

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
        mpz_set_ui(t, j);
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
        for (j = 1; j <= (UV) (exponents[i]>>1); j++) {
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
  mpz_clear(t1); mpz_clear(t2); mpz_clear(t3); mpz_init(t4); mpz_clear(t);
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
      if (mpz_cmp(U, V) >= 0)  mpz_sub(f, U, V);
      else                     mpz_sub(f, V, U);
      mpz_mul(m, m, f);
      mpz_tdiv_r(m, m, n);
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
  while (rounds > 0) {
    UV rleft = (r > rounds) ? rounds : r;
    while (rleft > 0) {   /* Do rleft rounds, inner at a time */
      UV dorounds = (rleft > inner) ? inner : rleft;
      mpz_set_ui(m, 1);
      mpz_set(saveXi, Xi);
      for (i = 0; i < dorounds; i++) {
        mpz_mul(t, Xi, Xi);  mpz_add_ui(t, t, a);  mpz_tdiv_r(Xi, t, n);
        if (mpz_cmp(Xi, Xm) >= 0)  mpz_sub(f, Xi, Xm);
        else                       mpz_sub(f, Xm, Xi);
        mpz_mul(t, m, f);
        mpz_tdiv_r(m, t, n);
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
        if (mpz_cmp(Xi, Xm) >= 0)  mpz_sub(f, Xi, Xm);
        else                       mpz_sub(f, Xm, Xi);
        mpz_gcd(f, f, n);
      } while (!mpz_cmp_ui(f, 1) && r-- != 0);
      if ( (!mpz_cmp_ui(f, 1)) || (!mpz_cmp(f, n)) )  break;
    }
    mpz_clear(Xi); mpz_clear(Xm); mpz_clear(m); mpz_clear(saveXi); mpz_clear(t);
    return 1;
  }
  mpz_clear(Xi); mpz_clear(Xm); mpz_clear(m); mpz_clear(saveXi); mpz_clear(t);
  mpz_set(f, n);
  return 0;
}

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
   * Montgomery 1987 p249-250 and Brent 1990 p5 both indicate we can calculate
   * a^m mod n where m is the lcm of the integers to B1.  This can be done
   * using either
   *    m = calc_lcm(B), b = a^m mod n
   * or
   *    calculate_b_lcm(b, B1, a, n);
   *
   * The first means raising a to a huge power then doing the mod, which is
   * inefficient and can be _very_ slow on some machines.  The latter does
   * one powmod for each prime power, which works pretty well.  Yet another
   * way to handle this is to loop over each prime p below B1, calculating
   * a = a^(p^e) mod n, where e is the largest e such that p^e <= B1.
   * My experience with GMP is that this last method is faster with large B1,
   * sometimes a lot faster.
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


/*----------------------------------------------------------------------
 * GMP version of Ben Buhrow's public domain 9/24/09 implementation.
 * It uses ideas and code from Jason Papadopoulos, Scott Contini, and
 * Tom St. Denis.  Also see the papers of Stephen McMath, Daniel Shanks,
 * and Jason Gower.  Gower and Wagstaff is particularly useful:
 *    http://homes.cerias.purdue.edu/~ssw/squfof.pdf
 *--------------------------------------------------------------------*/

static int shanks_mult(mpz_t n, mpz_t f)
{
   /*
    * use shanks SQUFOF to factor N.
    *
    * return 0 if no factor found, 1 if found with factor in f1.
    *
    * Input should have gone through trial division to 5.
    */

   int result = 0;
   unsigned long j=0;
   mpz_t b0, bn, imax, tmp, Q0, Qn, P, i, t1, t2, S, Ro, So, bbn;

   if (mpz_cmp_ui(n, 3) <= 0)
     return 0;

   if (mpz_perfect_square_p(n)) {
     mpz_sqrt(f, n);
     return 1;
   }

   mpz_init(b0);
   mpz_init(bn);
   mpz_init(imax);
   mpz_init(tmp);
   mpz_init(Q0);
   mpz_init(Qn);
   mpz_init(P);
   mpz_init(i);
   mpz_init(t1);
   mpz_init(t2);
   mpz_init(S);
   mpz_init(Ro);
   mpz_init(So);
   mpz_init(bbn);

   mpz_sqrt(b0, n);
   mpz_sqrt(tmp, b0);
   mpz_mul_ui(imax, tmp, 3);

   /* set up recurrence */
   mpz_set_ui(Q0, 1);
   mpz_set(P, b0);
   mpz_mul(tmp, b0, b0);
   mpz_sub(Qn, n, tmp);

   mpz_add(tmp, b0, P);
   mpz_tdiv_q(bn, tmp, Qn);

   mpz_set_ui(i, 0);
   while (1) {
      j=0;
      while (1) {
         mpz_set(t1, P);   /* hold Pn for this iteration */
         mpz_mul(tmp, bn, Qn);
         mpz_sub(P, tmp, P);
         mpz_set(t2, Qn);  /* hold Qn for this iteration */
         mpz_sub(tmp, t1, P);
         mpz_mul(tmp, tmp, bn);
         mpz_add(Qn, Q0, tmp);
         mpz_set(Q0, t2);  /* remember last Q */
         mpz_add(tmp, b0, P);
         mpz_tdiv_q(bn, tmp, Qn);

         if (mpz_even_p(i)) {
           if (mpz_perfect_square_p(Qn)) {
             mpz_add_ui(i, i, 1);
             break;
           }
         }
         mpz_add_ui(i, i, 1);

         if (mpz_cmp(i, imax) >= 0) {
           result = 0;
           goto end;
         }
      }

      /* reduce to G0 */
      mpz_sqrt(S, Qn);
      mpz_sub(tmp, b0, P);
      mpz_tdiv_q(tmp, tmp, S);
      mpz_mul(tmp, S, tmp);
      mpz_add(Ro, P, tmp);
      mpz_mul(tmp, Ro, Ro);
      mpz_sub(tmp, n, tmp);
      mpz_tdiv_q(So, tmp, S);
      mpz_add(tmp, b0, Ro);
      mpz_tdiv_q(bbn, tmp, So);

      /* search for symmetry point */
      while (1) {
         mpz_set(t1, Ro);  /* hold Ro for this iteration */
         mpz_mul(tmp, bbn, So);
         mpz_sub(Ro, tmp, Ro);
         mpz_set(t2, So);  /* hold So for this iteration */
         mpz_sub(tmp, t1, Ro);
         mpz_mul(tmp, bbn, tmp);
         mpz_add(So, S, tmp);
         mpz_set(S, t2);   /* remember last S */
         mpz_add(tmp, b0, Ro);
         mpz_tdiv_q(bbn, tmp, So);

         /* check for symmetry point */
         if (mpz_cmp(Ro, t1) == 0)
            break;

         /* this gets stuck very rarely, but it does happen. */
         if (++j > 1000000000)
         {
            result = -1;
            goto end;
         }
      }

      mpz_gcd(t1, Ro, n);
      if (mpz_cmp_ui(t1, 1) > 0) {
         mpz_set(f, t1);
         /* gmp_printf("GMP SQUFOF found factor after %Zd/%lu rounds: %Zd\n", i, j, f); */
         result = 1;
         goto end;
      }
   }

   end:
   mpz_clear(b0);
   mpz_clear(bn);
   mpz_clear(imax);
   mpz_clear(tmp);
   mpz_clear(Q0);
   mpz_clear(Qn);
   mpz_clear(P);
   mpz_clear(i);
   mpz_clear(t1);
   mpz_clear(t2);
   mpz_clear(S);
   mpz_clear(Ro);
   mpz_clear(So);
   mpz_clear(bbn);
   return result;
}

int _GMP_squfof_factor(mpz_t n, mpz_t f, UV rounds)
{
   const UV multipliers[] = {
      3*5*7*11, 3*5*7,  3*5*7*11*13, 3*5*7*13, 3*5*7*11*17, 3*5*11,
      3*5*7*17, 3*5,    3*5*7*11*19, 3*5*11*13,3*5*7*19,    3*5*7*13*17,
      3*5*13,   3*7*11, 3*7,         5*7*11,   3*7*13,      5*7,
      3*5*17,   5*7*13, 3*5*19,      3*11,     3*7*17,      3,
      3*11*13,  5*11,   3*7*19,      3*13,     5,           5*11*13,
      5*7*19,   5*13,   7*11,        7,        3*17,        7*13,
      11,       1 };
   const size_t sz_mul = sizeof(multipliers)/sizeof(multipliers[0]);
   size_t i;
   int result;
   mpz_t t;

   TEST_FOR_2357(n, f);
   mpz_init(t);
   mpz_set_ui(f, 1);

   for (i = 0; i < sz_mul; i++) {
      UV mult = multipliers[i];
      /* Only run when 64*m^3 < n */
      mpz_set_ui(t, mult);
      mpz_pow_ui(t, t, 3);
      mpz_mul_ui(t, t, 64);
      if (mpz_cmp(t, n) >= 0)
        continue;
      /* Run with this multiplier */
      mpz_mul_ui(t, n, mult);
      result = shanks_mult(t, f);
      if (result == -1)
        continue;
      if ( (result == 1) && (mpz_cmp_ui(f, mult) != 0) ) {
        unsigned long gcdf = mpz_gcd_ui(NULL, f, mult);
        mpz_divexact_ui(f, f, gcdf);
        if (mpz_cmp_ui(f, 1) > 0)
          break;
      }
   }
   mpz_clear(t);
   return (mpz_cmp_ui(f, 1) > 0);
}

/* See if n is a perfect power */
UV power_factor(mpz_t n, mpz_t f)
{
  UV k = 1, b = 2;
  if (mpz_cmp_ui(n, 1) > 0 && mpz_perfect_power_p(n)) {
    mpz_t nf, tf;
    PRIME_ITERATOR(iter);

    mpz_init_set(nf, n);
    mpz_init(tf);
    while (1) {
      UV ok = k;
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
