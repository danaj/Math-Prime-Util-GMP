
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
  UV i, limit;

  if (( 7 *  7) > n)  return 2;  if (!(n %  7)) return 0;
  if ((11 * 11) > n)  return 2;  if (!(n % 11)) return 0;
  if ((13 * 13) > n)  return 2;  if (!(n % 13)) return 0;
  if ((17 * 17) > n)  return 2;  if (!(n % 17)) return 0;
  if ((19 * 19) > n)  return 2;  if (!(n % 18)) return 0;
  if ((23 * 23) > n)  return 2;  if (!(n % 23)) return 0;
  if ((29 * 29) > n)  return 2;  if (!(n % 29)) return 0;

  limit = sqrt(n);
  i = 31;
  while (1) {   /* trial division, skipping multiples of 2/3/5 */
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 6;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 4;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 2;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 4;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 2;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 4;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 6;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 2;
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

/* Really simple little odd-only sieve. */
static unsigned char* sieve_erat(UV end)
{
  unsigned char* mem;
  UV n, s;
  UV last = (end+1)/2;

  Newz(0, mem, (last+7)/8, unsigned char);
  if (mem == 0) return 0;

  n = 3;
  while ( (n*n) <= end ) {
    for (s = n*n; s <= end; s += 2*n)
      mem[s/16] |= (1 << ((s/2) % 8));
    do { n += 2; } while (mem[n/16] & (1 << ((n/2) % 8)));
  }
  mem[0] |= 1;  /* 1 is composite */
  return mem;
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
      mpz_mul(T1, U2m, V2m);
      mpz_mod(U2m, T1, n);
      mpz_mul(T1, V2m, V2m);
      mpz_submul_ui(T1, Qm, 2);
      mpz_mod(V2m, T1, n);
      //gmp_printf("  l  U2m=%Zd  V2m=%Zd\n", U2m, V2m);
      mpz_mul(T1, Qm, Qm); mpz_mod(Qm, T1, n);
      if (mpz_odd_p(d)) {
        /* Save T1 and T2 for later operations in this block */
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
        mpz_fdiv_q_ui(T1, V, 2);
        mpz_mod(V, T1, n);
        /* Qkd */
        mpz_mul(T1, Qkd, Qm);
        mpz_mod(Qkd, T1, n);
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

#if 0
int _GMP_prho_factor(mpz_t n, mpz_t f, UV a, UV rounds)
{
  mpz_t U, V;
  int inloop = 0;

  mpz_init_set_ui(U, 7);
  mpz_init_set_ui(V, 7);
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
#else
int _GMP_prho_factor(mpz_t n, mpz_t f, UV a, UV rounds)
{
  mpz_t U, V, oldU, oldV, m;
  int i;
  const UV inner = 256;

  rounds = (rounds + inner - 1) / inner;
  mpz_init_set_ui(U, 7);
  mpz_init_set_ui(V, 7);
  mpz_init(m);
  mpz_init(oldU);
  mpz_init(oldV);
  while (rounds-- > 0) {
    mpz_set_ui(m, 1); mpz_set(oldU, U);  mpz_set(oldV, V);
    for (i = 0; i < inner; i++) {
      mpz_mul(U, U, U);  mpz_add_ui(U, U, a);  mpz_tdiv_r(U, U, n);
      mpz_mul(V, V, V);  mpz_add_ui(V, V, a);  mpz_tdiv_r(V, V, n);
      mpz_mul(V, V, V);  mpz_add_ui(V, V, a);  mpz_tdiv_r(V, V, n);
      mpz_sub(f, U, V);
      if (mpz_sgn(f) < 0)  mpz_add(f, f, n);
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
        mpz_sub(f, U, V);
        if (mpz_sgn(f) < 0)  mpz_add(f, f, n);
        mpz_gcd(f, f, n);
      } while (!mpz_cmp_ui(f, 1) && i-- != 0);
      if ( (!mpz_cmp_ui(f, 1)) || (!mpz_cmp(f, n)) )
        break;
    }
    mpz_clear(U); mpz_clear(V); mpz_clear(m); mpz_clear(oldU); mpz_clear(oldV);
    return 1;
  }
  mpz_clear(U); mpz_clear(V); mpz_clear(m); mpz_clear(oldU); mpz_clear(oldV);
  mpz_set(f, n);
  return 0;
}
#endif

#if 0
int _GMP_pbrent_factor(mpz_t n, mpz_t f, UV a, UV rounds)
{
  mpz_t Xi, Xm, saveXi, m;
  UV i, j, r;

  mpz_init_set_ui(Xi, 2);
  mpz_init_set_ui(Xm, 2);
  mpz_init(m);
  mpz_init(saveXi);

  r = 1;
  while (rounds > 0) {
    mpz_set_ui(m, 1); mpz_set(saveXi, Xi);
    for (i = 0; i < r; i++) {
      mpz_mul(Xi, Xi, Xi);  mpz_add_ui(Xi, Xi, a);  mpz_tdiv_r(Xi, Xi, n);
      mpz_sub(f, Xi, Xm);
      if (mpz_sgn(f) < 0)  mpz_add(f, f, n);
      mpz_mul(m, m, f);
      mpz_tdiv_r(m, m, n);
      if (i > rounds) break;
    }
    rounds = (rounds <= r) ? 0 : rounds-r;
    mpz_gcd(f, m, n);
    if (!mpz_cmp_ui(f, 1)) {
      r *= 2;
      mpz_set(Xm, Xi);
      continue;
    }
    if (!mpz_cmp(f, n)) {
      /* f == n, so we have to back up to see what factor got found */
      mpz_set(Xi, saveXi);
      do {
        mpz_mul(Xi, Xi, Xi);  mpz_add_ui(Xi, Xi, a);  mpz_tdiv_r(Xi, Xi, n);
        mpz_sub(f, Xi, Xm);
        if (mpz_sgn(f) < 0)  mpz_add(f, f, n);
        mpz_gcd(f, f, n);
      } while (!mpz_cmp_ui(f, 1) && r-- != 0);
    }
    mpz_clear(Xi); mpz_clear(Xm); mpz_clear(m); mpz_clear(saveXi);
    return (mpz_cmp_ui(f, 1) != 0);
  }
  mpz_clear(Xi); mpz_clear(Xm); mpz_clear(m); mpz_clear(saveXi);
  mpz_set(f, n);
  return 0;
}
#else
int _GMP_pbrent_factor(mpz_t n, mpz_t f, UV a, UV rounds)
{
  mpz_t Xi, Xm, saveXi, m;
  UV i, r;
  const UV inner = 256;

  mpz_init_set_ui(Xi, 2);
  mpz_init_set_ui(Xm, 2);
  mpz_init(m);
  mpz_init(saveXi);

  r = 1;
  while (rounds > 0) {
    UV rleft = (r > rounds) ? rounds : r;
    while (rleft > 0) {   /* Do rleft rounds, inner at a time */
      UV dorounds = (rleft > inner) ? inner : rleft;
      mpz_set_ui(m, 1);
      mpz_set(saveXi, Xi);
      for (i = 0; i < dorounds; i++) {
        mpz_mul(Xi, Xi, Xi);  mpz_add_ui(Xi, Xi, a);  mpz_tdiv_r(Xi, Xi, n);
        mpz_sub(f, Xi, Xm);
        if (mpz_sgn(f) < 0)  mpz_add(f, f, n);
        mpz_mul(m, m, f);
        mpz_tdiv_r(m, m, n);
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
        mpz_mul(Xi, Xi, Xi);  mpz_add_ui(Xi, Xi, a);  mpz_tdiv_r(Xi, Xi, n);
        mpz_sub(f, Xi, Xm);
        if (mpz_sgn(f) < 0)  mpz_add(f, f, n);
        mpz_gcd(f, f, n);
      } while (!mpz_cmp_ui(f, 1) && r-- != 0);
      if (!mpz_cmp_ui(f, 1)) break;
    }
    mpz_clear(Xi); mpz_clear(Xm); mpz_clear(m); mpz_clear(saveXi);
    return 1;
  }
  mpz_clear(Xi); mpz_clear(Xm); mpz_clear(m); mpz_clear(saveXi);
  mpz_set(f, n);
  return 0;
}
#endif


#if 0
static void lcm_to_B(UV B, mpz_t m)
{
  double logB = log(B);
  UV p = 1;
  UV exponent;

  mpz_set_ui(m, 1);
  while ( (p = next_small_prime(p)) <= B) {
    exponent = logB / log(p);
    if (exponent == 1) {
      mpz_mul_ui(m, m, p);
      break;
    }
    /* pow(p, exponent) < B, so fits in a UV */
    mpz_mul_ui(m, m, (UV)pow(p, exponent) );
  }
  /* All the exponent = 1 portion.  */
  while ( (p = next_small_prime(p)) <= B) {
    mpz_mul_ui(m, m, p);
  }
}
#else
static void lcm_to_B(UV B, mpz_t m)
{
  double logB = log(B);
  UV p, exponent;

  /* Simple sieve to B */
  unsigned char* s = sieve_erat(B);

  mpz_set_ui(m, 1);
  exponent = logB / log(2);
  if (B >= 2)  mpz_mul_ui(m, m, (UV)pow(2, exponent));
  p = 3;
  while ( (p <= B) && (exponent != 1) ) {
    exponent = logB / log(p);
    mpz_mul_ui(m, m, (UV)pow(p, exponent) );
    do { p += 2; } while (s[p/16] & (1UL << ((p/2) % 8)));
  }
  /* All the exponent = 1 portion.  */
  while ( p <= B ) {
    mpz_mul_ui(m, m, p);
    do { p += 2; } while (s[p/16] & (1UL << ((p/2) % 8)));
  }
  Safefree(s);
}
#endif

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
      /* gmp_printf("trying with a=%lu\n", p); */
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
    B *= 10; if (B > smoothness_bound) B = smoothness_bound;
  }
  mpz_clear(a); mpz_clear(m); mpz_clear(x);
  mpz_set(f, n);
  return 0;
}

/* Alternate way.  Much less memory for really big B.  On some machines this
 * is slower.  On others the lcm process above is agonizingly slow, so large
 * smoothness factors work far better with this.
 *
 * The above method does the textbook p-1, with a smoothness bound B generating
 * a stupendously large M, then we pick an a to examine GCD(a^M - 1 mod n, n).
 * In contrast, this function examines GCD(b^M - 1 mod n, n) where b is
 * semi-random between 1 and n, and M is just the loop counter.  This doesn't
 * seem like the same thing at all.
 */
int _GMP_pminus1_factor2(mpz_t n, mpz_t f, UV rounds)
{
  mpz_t b;
  UV loops;

  mpz_init_set_ui(b, 13);

  for (loops = 1; loops <= rounds; loops++) {
    mpz_add_ui(b, b, 1);
    mpz_powm_ui(b, b, loops, n);
    if (mpz_sgn(b) == 0) mpz_set(b, n);
    mpz_sub_ui(b, b, 1);
    mpz_gcd(f, b, n);
    if (mpz_cmp(f, n) == 0)
      break;
    if (mpz_cmp_ui(f, 1) > 0) {
      mpz_clear(b);
      return 1;
    }
  }
  mpz_clear(b);
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
    mpz_mul_ui(f, n, i);    /* f = n*i */
    if (mpz_perfect_square_p(f)) {
      /* s^2 = n*i, so m = s^2 mod n = 0.  Hence f = GCD(n, s) = GCD(n, n*i) */
      mpz_divexact_ui(n, n, PREMULT);
      mpz_gcd(f, f, n);
      mpz_clear(s); mpz_clear(m); return 1;
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
      mpz_clear(s); mpz_clear(m); return 1;
    }
  }
  mpz_divexact_ui(n, n, PREMULT);
  mpz_set(f, n);
  mpz_clear(s); mpz_clear(m);
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
    croak("GMP_trial_factor from > to: %"UVuf" - %"UVuf, from_n, to_n);
  if (to_n > primes_small[NPRIMES_SMALL-1])
    croak("GMP_trial_factor too large: %"UVuf" - %"UVuf, from_n, to_n);

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


/*----------------------------------------------------------------------
 * GMP version of Ben Buhrow's public domain 9/24/09 implementation.
 * It uses ideas and code from Jason Papadopoulos, Scott Contini, and
 * Tom St. Denis.  Also see the papers of Stephen McMath, Daniel Shanks,
 * and Jason Gower.
 *--------------------------------------------------------------------*/

static int shanks_mult(mpz_t n, mpz_t f)
{
   // use shanks SQUFOF to factor N.
   //
   // return 0 if no factor found, 1 if found with factor in f1.
   //
   // Input should have gone through trial division to 5.

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

   //set up recurrence
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
         mpz_set(t1, P);  //hold Pn for this iteration
         mpz_mul(tmp, bn, Qn);
         mpz_sub(P, tmp, P);
         mpz_set(t2, Qn);  //hold Qn for this iteration
         mpz_sub(tmp, t1, P);
         mpz_mul(tmp, tmp, bn);
         mpz_add(Qn, Q0, tmp);
         mpz_set(Q0, t2);  //remember last Q
         mpz_add(tmp, b0, P);
         mpz_tdiv_q(bn, tmp, Qn);

         if (mpz_even_p(i)) {   // i even
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

      //reduce to G0
      mpz_sqrt(S, Qn);
      mpz_sub(tmp, b0, P);
      mpz_tdiv_q(tmp, tmp, S);
      mpz_mul(tmp, S, tmp);
      mpz_add(Ro, P, tmp);
      //mpz_set(t1, Ro);
      //mpz_mul(tmp, t1, t1);
      mpz_mul(tmp, Ro, Ro);
      mpz_sub(tmp, n, tmp);
      mpz_tdiv_q(So, tmp, S);
      mpz_add(tmp, b0, Ro);
      mpz_tdiv_q(bbn, tmp, So);

      //search for symmetry point
      while (1) {
         mpz_set(t1, Ro);  //hold Ro for this iteration
         mpz_mul(tmp, bbn, So);
         mpz_sub(Ro, tmp, Ro);
         mpz_set(t2, So);  //hold So for this iteration
         mpz_sub(tmp, t1, Ro);
         mpz_mul(tmp, bbn, tmp);
         mpz_add(So, S, tmp);
         mpz_set(S, t2);   //remember last S
         mpz_add(tmp, b0, Ro);
         mpz_tdiv_q(bbn, tmp, So);

         //check for symmetry point
         if (mpz_cmp(Ro, t1) == 0)
            break;

         //this gets stuck very rarely, but it does happen.
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
   //call shanks with multiple small multipliers
   const int multipliers[] = {3, 5, 7, 
            11, 3*5, 3*7, 3*11, 
            5*7, 5*11, 7*11, 
            3*5*7, 3*5*11, 3*7*11, 
            5*7*11, 3*5*7*11, 1};

   const int sz_mul = 16;
   int i;
   int result;
   mpz_t nm;

   if (mpz_cmp_ui(n, 1) <= 0)
     return 0;

   mpz_init(nm);
   mpz_set_ui(f, 1);

   for (i=sz_mul-1;i>=0;i--) {
      mpz_mul_ui(nm, n, multipliers[i]);
      result = shanks_mult(nm, f);
      if (result == -1)
        break;
      if ( (result == 1) && (mpz_cmp_ui(f, multipliers[i]) != 0) ) {
        unsigned long gcdf = mpz_gcd_ui(NULL, f, multipliers[i]);
        mpz_divexact_ui(f, f, gcdf);
        if (mpz_cmp_ui(f, 1) > 0)
          break;
      }
   }
   mpz_clear(nm);
   return (mpz_cmp_ui(f, 1) > 0);
}
