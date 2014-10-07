#include "factor.h"
#include "gmp_main.h"
#include "prime_iterator.h"
#include "utility.h"
#include "small_factor.h"
#include "ecm.h"
#include "simpqs.h"

#define _GMP_ECM_FACTOR(n, f, b1, ncurves) \
   _GMP_ecm_factor_projective(n, f, b1, 0, ncurves)

static const unsigned short primes_small[] =
  {0,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
   101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,
   193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,
   293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,
   409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,
   521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,
   641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,
   757,761,769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,
   881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997,1009,
   1013,1019,1021,1031,1033,1039,1049,1051,1061,1063,1069,1087,1091,1093,
   1097,1103,1109,1117,1123,1129,1151,1153,1163,1171,1181,1187,1193,1201,
   1213,1217,1223,1229,1231,1237,1249,1259,1277,1279,1283,1289,1291,1297,
   1301,1303,1307,1319,1321,1327,1361,1367,1373,1381,1399,1409,1423,1427,
   1429,1433,1439,1447,1451,1453,1459,1471,1481,1483,1487,1489,1493,1499,
   1511,1523,1531,1543,1549,1553,1559,1567,1571,1579,1583,1597,1601,1607,
   1609,1613,1619,1621,1627,1637,1657,1663,1667,1669,1693,1697,1699,1709,
   1721,1723,1733,1741,1747,1753,1759,1777,1783,1787,1789,1801,1811,1823,
   1831,1847,1861,1867,1871,1873,1877,1879,1889,1901,1907,1913,1931,1933,
   1949,1951,1973,1979,1987,1993,1997,1999,2003,2011};
#define NPRIMES_SMALL (sizeof(primes_small)/sizeof(primes_small[0]))

#define TRIAL_LIM 2000
#define MAX_FACTORS 256

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
  UV tf;

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
  {
    UV sp, p, un;
    un = (mpz_cmp_ui(n,2*TRIAL_LIM*TRIAL_LIM) >= 0) ? 2*TRIAL_LIM*TRIAL_LIM
                                                    : mpz_get_ui(n);

    for (sp = 2, p = primes_small[sp];
         p < TRIAL_LIM && p*p <= un;
         p = primes_small[++sp]) {
      while (mpz_divisible_ui_p(n, p)) {
        ADD_FACTOR_UI(f, p);
        mpz_divexact_ui(n, n, p);
        un = (mpz_cmp_ui(n,2*TRIAL_LIM*TRIAL_LIM) > 0) ? 2*TRIAL_LIM*TRIAL_LIM
                                                       : mpz_get_ui(n);
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

      if (!success)  success = _GMP_pminus1_factor(n, f, 15000, 150000);
      if (success&&o) {gmp_printf("p-1 (10k) found factor %Zd\n", f);o=0;}

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
