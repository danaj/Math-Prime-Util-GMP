//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <math.h>
//#include <time.h>

#include "factor.h"
#include "gmp_main.h"
#include "prime_iterator.h"

#define TRIAL_LIM 2000
#define MAX_FACTORS 256
#define _GMP_ECM_FACTOR(n, f, b1, ncurves) \
   _GMP_ecm_factor_projective(n, f, b1, 0, ncurves)

static int add_factor(int nfactors, mpz_t f, mpz_t** pfactors, int** pexponents)
{
  int i, j, cmp;
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

int factor(mpz_t n, mpz_t* pfactors[], int* pexponents[])
{
  mpz_t tofac_stack[MAX_FACTORS];
  int ntofac = 0;
  mpz_t* factors;
  int* exponents;
  int nfactors = 0;
  mpz_t f;
  UV tf;

  mpz_init(f);
  if (mpz_cmp_ui(n, 4) < 0) {
    ADD_FACTOR(n);
    goto DONE;
  }

  /* Trial factor to small limit */
  while (mpz_even_p(n)) {
    ADD_FACTOR_UI(f, 2);
    mpz_divexact_ui(n, n, 2);
  }
  {
    UV p;
    PRIME_ITERATOR(iter);

    for (p = prime_iterator_next(&iter);
         p < TRIAL_LIM && mpz_cmp_ui(n, p*p) >= 0;
         p = prime_iterator_next(&iter)) {
      while (mpz_divisible_ui_p(n, p)) {
        ADD_FACTOR_UI(f, p);
        mpz_divexact_ui(n, n, p);
      }
    }
    if (mpz_cmp_ui(n, p*p) < 0) {
      if (mpz_cmp_ui(n, 1) > 0)
        ADD_FACTOR(n);
      goto DONE;
    }
  }

  /* Power factor */
  tf = power_factor(n, f);
  if (tf) {
    mpz_t* pow_factors;
    int* pow_exponents;
    int pow_nfactors;
    int i, j;

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
    while ( mpz_cmp_ui(n, TRIAL_LIM*TRIAL_LIM) > 0 && !_GMP_is_prob_prime(n) ) {
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

      if (mpz_cmp_ui(n, (unsigned long)(UV_MAX>>2)) < 0) {
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

      if (!success)  success = _GMP_pminus1_factor(n, f, 10000, 150000);
      if (success&&o) {gmp_printf("p-1 (10k) found factor %Zd\n", f);o=0;}

      /* Really small ECM to find small factors */
      if (!success)  success = _GMP_ECM_FACTOR(n, f, 150, 30);
      if (success&&o) {gmp_printf("tiny ecm (150) found factor %Zd\n", f);o=0;}
      if (!success)  success = _GMP_ECM_FACTOR(n, f, 500, 25);
      if (success&&o) {gmp_printf("tiny ecm (500) found factor %Zd\n", f);o=0;}
      if (!success)  success = _GMP_ECM_FACTOR(n, f, 2000, 10);
      if (success&&o) {gmp_printf("tiny ecm (2000) found factor %Zd\n", f);o=0;}

      /* Small p-1 */
      if (!success)  success = _GMP_pminus1_factor(n, f, 200000, 3000000);
      if (success&&o) {gmp_printf("p-1 (200k) found factor %Zd\n", f);o=0;}

      /* Set ECM parameters that have a good chance of success */
      if (!success) {
        UV curves;
        nbits = mpz_sizeinbase(n, 2);
        if      (nbits < 100){ B1 =   5000; curves =  20; }
        else if (nbits < 128){ B1 =  10000; curves =   5; } /* go to QS */
        else if (nbits < 160){ B1 =  20000; curves =   5; } /* go to QS */
        else if (nbits < 192){ B1 =  30000; curves =  20; }
        else if (nbits < 224){ B1 =  40000; curves =  40; }
        else if (nbits < 256){ B1 =  80000; curves =  40; }
        else if (nbits < 512){ B1 = 160000; curves =  80; }
        else                 { B1 = 320000; curves = 160; }
        success = _GMP_ECM_FACTOR(n, f, B1, curves);
        if (success&&o) {gmp_printf("small ecm (%luk,%lu) found factor %Zd\n", B1/1000, curves, f);o=0;}
      }

      /* QS (30+ digits).  Fantastic if it is a semiprime, but can be
       * slow and a memory hog if not (compared to ECM).  Restrict to
       * reasonable size numbers (< 91 digits).  Because of the way it
       * works, it will generate (possibly) multiple factors for the same
       * amount of work.  Go to some trouble to use them. */
      if (!success && mpz_sizeinbase(n,10) >= 30 && nbits < 300) {
        mpz_t farray[66];
        int i, nfactors;
        for (i = 0; i < 66; i++)
          mpz_init(farray[i]);
        nfactors = _GMP_simpqs(n, farray);
        mpz_set(f, farray[0]);
        if (nfactors > 2) {
          /* We found multiple factors */
          for (i = 2; i < nfactors; i++) {
            if (o){gmp_printf("SIMPQS found extra factor %Zd\n",farray[i]);}
            if (ntofac == MAX_FACTORS-1) croak("Too many factors\n");
            mpz_init_set(tofac_stack[ntofac], farray[i]);
            ntofac++;
            mpz_divexact(n, n, farray[i]);
          }
          /* f = farray[0], n = farray[1], farray[2..] pushed */
        }
        for (i = 0; i < 66; i++)
          mpz_clear(farray[i]);
        success = nfactors > 1;
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
      } else if (_GMP_is_prob_prime(f)) {
        int ndiv = mpz_remove(n, n, f);
        while (ndiv-- > 0)
          ADD_FACTOR(f);
      } else {
        int ndiv = mpz_remove(n, n, f);
        if (ntofac == MAX_FACTORS-ndiv)
          croak("Too many factors\n");
        while (ndiv-- > 0)
          mpz_init_set(tofac_stack[ntofac++], f);
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
