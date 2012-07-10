
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include "gmp_main.h"
#include <math.h>

static const unsigned short primes_small[] =
  {0,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
   101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,
   193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,
   293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,
   409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,
   521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,
   641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,
   757,761,769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,
   881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997,1009
  };
#define NPRIMES_SMALL (sizeof(primes_small)/sizeof(primes_small[0]))

static const unsigned char next_wheel[30] =
  {1,7,7,7,7,7,7,11,11,11,11,13,13,17,17,17,17,19,19,23,23,23,23,29,29,29,29,29,29,1};
static const unsigned char prev_wheel[30] =
  {29,29,1,1,1,1,1,1,7,7,7,7,11,11,13,13,13,13,17,17,19,19,19,19,23,23,23,23,23,23};


static int _is_small_prime7(UV n)
{
  UV limit = sqrt(n);
  UV i = 7;
  while (1) {   /* trial division, skipping multiples of 2/3/5 */
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 4;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 2;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 4;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 2;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 4;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 6;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 2;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 6;
  }
  return 2;
}

static UV next_small_prime(UV n)
{
  if (n < 7)
    return (n < 2) ? 2 : (n < 3) ? 3 : (n < 5) ? 5 : 7;

  UV d = n/30;
  UV m = n - d*30;
    /* Move forward one, knowing we may not be on the wheel */
  if (m == 29) { d++; m = 1; } else  { m = next_wheel[m]; }
  while (!_is_small_prime7(d*30+m)) {
    m = next_wheel[m];  if (m == 1) d++;
  }
  return(d*30+m);
}

static int _GMP_miller_rabin_ui(mpz_t n, UV base)
{
  int rval;
  mpz_t a;
  mpz_init_set_ui(a, base);
  rval = _GMP_miller_rabin(n, a);
  mpz_clear(a);
  return rval;
}

int _GMP_miller_rabin(mpz_t n, mpz_t a)
{
  mpz_t nminus1, d, x;
  UV s, r;
  int rval;

  /* gmp_printf("Computing MR base %Zd on %Zd\n", a, n); */
  {
    int cmpr = mpz_cmp_ui(n, 2);
    if (cmpr == 0)     return 1;  /* 2 is prime */
    if (cmpr < 0)      return 0;  /* below 2 is composite */
    if (mpz_even_p(n)) return 0;  /* multiple of 2 is composite */
  }
  mpz_init_set(nminus1, n);
  mpz_sub_ui(nminus1, nminus1, 1);
  mpz_init_set(d, nminus1);
  s = 0;
  while (mpz_even_p(d)) {
    s++;
    mpz_divexact_ui(d, d, 2);
  }
  /* faster way, verify s is identical
   *    s = mpz_scan1(d, 0);
   *    mpz_tdiv_q_2exp(d, d, s);
   */
  mpz_init(x);
  mpz_powm(x, a, d, n);
  mpz_clear(d); /* done with a and d */
  rval = 0;
  if (!mpz_cmp_ui(x, 1) || !mpz_cmp(x, nminus1)) {
    rval = 1;
  } else {
    for (r = 0; r < s; r++) {
      mpz_powm_ui(x, x, 2, n);
      if (!mpz_cmp_ui(x, 1)) {
        break;
      }
      if (!mpz_cmp(x, nminus1)) {
        rval = 1;
        break;
      }
    }
  }
  mpz_clear(nminus1); mpz_clear(x);
  return rval;
}

int _GMP_is_strong_lucas_pseudoprime(mpz_t n)
{
  mpz_t d, U, V, Qkd;
  IV D;
  UV P = 1;
  IV Q;
  UV s;
  int rval;

  {
    int cmpr = mpz_cmp_ui(n, 2);
    if (cmpr == 0)     return 1;  /* 2 is prime */
    if (cmpr < 0)      return 0;  /* below 2 is composite */
    if (mpz_even_p(n)) return 0;  /* multiple of 2 is composite */
    if (mpz_perfect_square_p(n)) return 0;  /* n perfect square is composite */
  }
  /* Determine Selfridge D, P, Q parameters */
  {
    mpz_t t;
    UV D_ui = 5;
    IV sign = 1;
    mpz_init(t);
    while (1) {
      UV gcd, j;
      gcd = mpz_gcd_ui(NULL, n, D_ui);
      if ((gcd > 1) && mpz_cmp_ui(n, gcd) != 0) {
        D_ui = 0;
        break;
      }
      mpz_set_si(t, (IV)D_ui * sign);
      j = mpz_jacobi(t, n);
      if (j == -1)  break;
      D_ui += 2;
      sign = -sign;
    }
    mpz_clear(t);
    if (D_ui == 0)  return 0;
    D = (IV)D_ui * sign;
  }
  Q = (1 - D) / 4;
  //gmp_printf("N: %Zd  D: %ld  P: %lu  Q: %ld\n", n, D, P, Q);
  if (D != P*P - 4*Q)  croak("incorrect DPQ\n");
  /* Now start on the lucas sequence */
  mpz_init_set(d, n);
  mpz_add_ui(d, d, 1);
  s = 0;
  while (mpz_even_p(d)) {
    s++;
    mpz_divexact_ui(d, d, 2);
  }
  mpz_init_set_ui(U, 1);
  mpz_init_set_ui(V, P);
  {
    mpz_t U2m, V2m, Qm, T1, T2;
    mpz_init_set(U2m, U);
    mpz_init_set(V2m, V);
    mpz_init_set_si(Qm, Q);
    mpz_init_set(Qkd, Qm);
    mpz_tdiv_q_ui(d, d, 2);
    mpz_init(T1);
    mpz_init(T2);
    while (mpz_sgn(d) > 0) {
      //gmp_printf("U=%Zd  V=%Zd  Qm=%Zd\n", U, V, Qm);
      mpz_mul(U2m, U2m, V2m);
      mpz_mod(U2m, U2m, n);
      mpz_powm_ui(V2m, V2m, 2, n);
      mpz_submul_ui(V2m, Qm, 2);
      mpz_mod(V2m, V2m, n);
      //gmp_printf("  l  U2m=%Zd  V2m=%Zd\n", U2m, V2m);
      mpz_powm_ui(Qm, Qm, 2, n);
      if (mpz_odd_p(d)) {
        mpz_mul(T1, U2m, V);
        mpz_mul(T2, U2m, U);
        mpz_mul_si(T2, T2, D);
        //gmp_printf("      T1 %Zd  T2 %Zd\n", T1, T2);
        /* U */
        mpz_mul(U, U, V2m);
        mpz_add(U, U, T1);
        if (mpz_odd_p(U)) mpz_add(U, U, n);
        mpz_fdiv_q_ui(U, U, 2);
        mpz_mod(U, U, n);
        /* V */
        mpz_mul(V, V, V2m);
        mpz_add(V, V, T2);
        if (mpz_odd_p(V)) mpz_add(V, V, n);
        mpz_fdiv_q_ui(V, V, 2);
        mpz_mod(V, V, n);
        /* Qkd */
        mpz_mul(Qkd, Qkd, Qm);
        mpz_mod(Qkd, Qkd, n);
      }
      mpz_tdiv_q_ui(d, d, 2);
    }
    mpz_clear(U2m); mpz_clear(V2m); mpz_clear(Qm); mpz_clear(T1); mpz_clear(T2);
  }
  rval = 0;
  //gmp_printf("l0 U=%Zd  V=%Zd\n", U, V);
  if ( (mpz_sgn(U) == 0) || (mpz_sgn(V) == 0) ) {
    rval = 1;
    s = 0;
  }
  /* Powers of V */
  while (s--) {
    mpz_mul(V, V, V);
    mpz_submul_ui(V, Qkd, 2);
    mpz_mod(V, V, n);
    if (mpz_sgn(V) == 0) {
      rval = 1;
      break;
    }
    if (s)
      mpz_powm_ui(Qkd, Qkd, 2, n);
  }
  mpz_clear(d); mpz_clear(U); mpz_clear(V); mpz_clear(Qkd);
  return rval;
}

int _GMP_trial_div(mpz_t n, UV to_n)
{
  int small_n = 0;
  int primei = 2;
  UV f;
  {
    int cmpr = mpz_cmp_ui(n, 2);
    if (cmpr == 0)     return 2;  /* 2 is prime */
    if (cmpr < 0)      return 0;  /* below 2 is composite */
    if (mpz_even_p(n)) return 0;  /* multiple of 2 is composite */
  }

  /* This is a really simple function which should be replaced with a
   * proper version sometime.  It's not too bad up to the small primes
   * limit, but is shockingly stupid after that.  There are also much
   * faster ways to test the small numbers (e.g. gcd within a limb) */
  if (mpz_cmp_ui(n, to_n*to_n) < 0)
    small_n = 1;
  f = primes_small[primei];
  while (f <= to_n) {
    if (small_n && mpz_cmp_ui(n, f*f) < 0) return 2;
    //gmp_printf("  is %Zd divisible by %lu?\n", n, primes_small[primei]);
    if (mpz_divisible_ui_p(n, primes_small[primei])) return 0;
    if (++primei < NPRIMES_SMALL) {
      f = primes_small[primei];
    } else {
      f = next_small_prime(f);
    }
  }
  return 1;
}

int _GMP_is_prob_prime(mpz_t n)
{
  int rval;

  /* Trial divide some small primes. Returns:
   *     0  Small factor found, number is composite.
   *     1  No small factors found.  No answer.
   *     2  All primes to sqrt(n) tested, number is definitely prime.
   */
  rval = _GMP_trial_div(n, 400);
  if (rval != 1)  return rval;

  /* Miller Rabin with base 2 */
  if (_GMP_miller_rabin_ui(n, 2) == 0)
    return 0;

  /* Strong Lucas-Selfridge */
  if (_GMP_is_strong_lucas_pseudoprime(n) == 0)
    return 0;

  /* BPSW is deterministic below 2^64 */
  if (mpz_sizeinbase(n, 2) <= 64)
    return 2;

  return 1;
}

int _GMP_is_prime(mpz_t n)
{
  return _GMP_is_prob_prime(n);
}

/* Modifies argument */
void _GMP_next_prime(mpz_t n)
{
  mpz_t d;
  UV m;

  /* small inputs */
  if (mpz_cmp_ui(n, 1) <= 0) { mpz_set_ui(n, 2); return; }
  if (mpz_cmp_ui(n, 2) <= 0) { mpz_set_ui(n, 3); return; }
  if (mpz_cmp_ui(n, 4) <= 0) { mpz_set_ui(n, 5); return; }

  mpz_init(d);
  m = mpz_fdiv_q_ui(d, n, 30);
  
  if (m == 29) {
    mpz_add_ui(d, d, 1);
    m = 1;
  } else {
    m = next_wheel[m];
  }
  while (1) {
    mpz_mul_ui(n, d, 30);
    mpz_add_ui(n, n, m);
    if (_GMP_is_prob_prime(n))
      break;
    m = next_wheel[m];
    if (m == 1)
      mpz_add_ui(d, d, 1);
  }
  mpz_clear(d);
}

/* Modifies argument */
void _GMP_prev_prime(mpz_t n)
{
  mpz_t d;
  UV m;

  /* small inputs */
  if (mpz_cmp_ui(n, 2) <= 0) { mpz_set_ui(n, 0); return; }
  if (mpz_cmp_ui(n, 3) <= 0) { mpz_set_ui(n, 2); return; }
  if (mpz_cmp_ui(n, 5) <= 0) { mpz_set_ui(n, 3); return; }
  if (mpz_cmp_ui(n, 7) <= 0) { mpz_set_ui(n, 5); return; }

  mpz_init(d);
  m = mpz_fdiv_q_ui(d, n, 30);
  
  while (1) {
    m = prev_wheel[m];
    if (m == 29)
      mpz_sub_ui(d, d, 1);
    mpz_mul_ui(n, d, 30);
    mpz_add_ui(n, n, m);
    if (_GMP_is_prob_prime(n))
      break;
  }
  mpz_clear(d);
}

static UV basic_factor(mpz_t n)
{
  if (mpz_cmp_ui(n, 4) < 0)      return mpz_get_ui(n);
  if (mpz_even_p(n))             return 2;
  if (mpz_divisible_ui_p(n, 3))  return 3;
  if (mpz_divisible_ui_p(n, 5))  return 5;
  return 0;
}

int _GMP_prho_factor(mpz_t n, mpz_t f, UV a, UV rounds)
{
  mpz_t U, V;
  int inloop = 0;

  mpz_init(U);
  mpz_init(V);
  while (rounds-- > 0) {
    mpz_mul(U, U, U);  mpz_add_ui(U, U, a);  mpz_mod(U, U, n);
    mpz_mul(V, V, V);  mpz_add_ui(V, V, a);  mpz_mod(V, V, n);
    mpz_mul(V, V, V);  mpz_add_ui(V, V, a);  mpz_mod(V, V, n);
    mpz_sub(f, U, V);
    mpz_gcd(f, f, n);
    if (!mpz_cmp(f, n)) {
      if (inloop++) break;
    } else if (mpz_cmp_ui(f, 1) != 0) {
      mpz_clear(U); mpz_clear(V);
      return 1;
    }
  }
  mpz_clear(U); mpz_clear(V);
  mpz_set(f, n);
  return 0;
}

static void lcm_to_B(UV B, mpz_t m)
{
  mpz_t t;
  double logB = log(B);
  UV p = 1;

  mpz_init(t);
  mpz_set_ui(m, 1);
  while ( (p = next_small_prime(p)) <= B) {
    /* mpz_mul_ui(m, m, (UV)pow(p, (UV)( logB / log(p) )) ); */
    mpz_set_ui(t, p);
    mpz_pow_ui(t, t, (UV)( logB / log(p) ));
    mpz_mul(m, m, t);
  }
  mpz_clear(t);
}

int _GMP_pminus1_factor(mpz_t n, mpz_t f, UV smoothness_bound)
{
  mpz_t a, m, x;
  UV i, p;
  UV B;

  mpz_init(a);
  mpz_init(m);
  mpz_init(x);

  B = 5;
  while (B <= smoothness_bound) {
    /* gmp_printf("   calculating new m...\n"); */
    lcm_to_B(B, m);
    /* gmp_printf("trying %Zd  with B=%lu", n, B); if (B<100) gmp_printf(" m=%Zd", m); gmp_printf("\n"); */
    /* Use primes for a's to try.  We rarely make it past the first couple. */
    p = 1;
    while ( (p = next_small_prime(p)) < 200000 ) {
      mpz_set_ui(a, p);
      mpz_powm(x, a, m, n);
      if (mpz_sgn(x) == 0)
        mpz_set(x, n);
      mpz_sub_ui(x, x, 1);
      mpz_gcd(f, x, n);
      if (mpz_cmp_ui(f, 1) == 0)
        break;
      if (mpz_cmp(f, n) != 0) {
        mpz_clear(a); mpz_clear(m); mpz_clear(x);
        return 1;
      }
    }
    if (B == smoothness_bound) break;
    B *= 3; if (B > smoothness_bound) B = smoothness_bound;
  }
  mpz_clear(a); mpz_clear(m); mpz_clear(x);
  mpz_set(f, n);
  return 0;
}

int _GMP_holf_factor(mpz_t n, mpz_t f, UV rounds)
{
  mpz_t s, m;
  UV i;

#define PREMULT 480   /* 1  2  6  12  480  151200 */

  mpz_mul_ui(n, n, PREMULT);
  mpz_init(s);
  mpz_init(m);
  for (i = 1; i <= rounds; i++) {
    mpz_mul_ui(s, n, i);
    mpz_sqrtrem(s, m, s);    /* s = sqrt(n*i), m = remainder */
    if (mpz_sgn(m) != 0)
      mpz_add_ui(s, s, 1);   /* s++ if s*s != n*i */
    mpz_powm_ui(m, s, 2, n); /* m = s^2 % n */
    if (mpz_perfect_square_p(m)) {
      mpz_sqrt(f, m);
      mpz_sub(s, s, f);
      mpz_divexact_ui(n, n, PREMULT);
      mpz_gcd(f, s, n);
      mpz_clear(s); mpz_clear(m); return 1;
    }
  }
  mpz_clear(s); mpz_clear(m);
  mpz_divexact_ui(n, n, PREMULT);
  mpz_set(f, n);
  return 0;
}


UV _GMP_trial_factor(mpz_t n, UV from_n, UV to_n)
{
  int small_n = 0;
  int primei = 2;
  UV f;

  f = basic_factor(n);
  if (f > 0)
    return f;

  if (from_n > to_n)
    croak("GMP_trial_factor from > to");
  if (to_n > primes_small[NPRIMES_SMALL-1])
    croak("GMP_trial_factor too large");

  if (mpz_cmp_ui(n, to_n*to_n) < 0)
    small_n = 1;

  f = primes_small[primei];
  while (f < from_n) {
    f = primes_small[++primei];
  }

  while (f <= to_n) {
    if (small_n && mpz_cmp_ui(n, f*f) < 0) return 0;
    if (mpz_divisible_ui_p(n, f)) return f;
    if (primei == NPRIMES_SMALL-1)  return 0;
    f = primes_small[++primei];
  }
  return 0;
}
