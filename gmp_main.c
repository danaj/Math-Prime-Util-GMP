
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include "gmp_main.h"
#include "prime_iterator.h"
#include "small_factor.h"
#include "simpqs.h"
#include "ecm.h"
#define _GMP_ECM_FACTOR _GMP_ecm_factor_projective
#include "utility.h"

static int _verbose = 0;
void _GMP_set_verbose(int v) { _verbose = v; }
int _GMP_get_verbose(void) { return _verbose; }

static gmp_randstate_t _randstate;
gmp_randstate_t* _GMP_get_randstate(void) { return &_randstate; }

void _GMP_init(void)
{
  /* We should  not use this random number system for crypto, so
   * using this lousy seed is ok.  We just would like something a
   * bit different every run. */
  unsigned long seed = time(NULL);
  gmp_randinit_mt(_randstate);
  gmp_randseed_ui(_randstate, seed);
  prime_iterator_global_startup();
}

void _GMP_destroy(void)
{
  prime_iterator_global_shutdown();
  gmp_randclear(_randstate);
}


static const unsigned char next_wheel[30] =
  {1,7,7,7,7,7,7,11,11,11,11,13,13,17,17,17,17,19,19,23,23,23,23,29,29,29,29,29,29,1};
static const unsigned char prev_wheel[30] =
  {29,29,1,1,1,1,1,1,7,7,7,7,11,11,13,13,13,13,17,17,19,19,19,19,23,23,23,23,23,23};
static const unsigned char wheel_advance[30] =
  {0,6,0,0,0,0,0,4,0,0,0,2,0,4,0,0,0,2,0,4,0,0,0,6,0,0,0,0,0,2};


static INLINE int _GMP_miller_rabin_ui(mpz_t n, UV base)
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

  {
    int cmpr = mpz_cmp_ui(n, 2);
    if (cmpr == 0)     return 1;  /* 2 is prime */
    if (cmpr < 0)      return 0;  /* below 2 is composite */
    if (mpz_even_p(n)) return 0;  /* multiple of 2 is composite */
  }
  mpz_init_set(nminus1, n);
  mpz_sub_ui(nminus1, nminus1, 1);
  mpz_init_set(x, a);

  /* Handle large and small bases.  Use x so we don't modify their input a. */
  if (mpz_cmp(x, n) >= 0)
    mpz_mod(x, x, n);
  if ( (mpz_cmp_ui(x, 1) <= 0) || (mpz_cmp(x, nminus1) >= 0) ) {
    mpz_clear(nminus1);
    mpz_clear(x);
    return 1;
  }

  mpz_init_set(d, nminus1);
  s = mpz_scan1(d, 0);
  mpz_tdiv_q_2exp(d, d, s);

  mpz_powm(x, x, d, n);
  mpz_clear(d); /* done with a and d */
  rval = 0;
  if (!mpz_cmp_ui(x, 1) || !mpz_cmp(x, nminus1)) {
    rval = 1;
  } else {
    for (r = 1; r < s; r++) {
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
      UV gcd = mpz_gcd_ui(NULL, n, D_ui);
      if ((gcd > 1) && mpz_cmp_ui(n, gcd) != 0) {
        D_ui = 0;
        break;
      }
      mpz_set_si(t, (IV)D_ui * sign);
      if (mpz_jacobi(t, n) == -1)  break;
      D_ui += 2;
      sign = -sign;
    }
    mpz_clear(t);
    if (D_ui == 0)  return 0;
    D = (IV)D_ui * sign;
  }
  Q = (1 - D) / 4;
  if (_verbose>3) gmp_printf("N: %Zd  D: %ld  P: %lu  Q: %ld\n", n, D, P, Q);
  if (D != ((IV)(P*P)) - 4*Q)  croak("incorrect DPQ\n");
  /* Now start on the lucas sequence */
  mpz_init_set(d, n);
  mpz_add_ui(d, d, 1);

  s = mpz_scan1(d, 0);
  mpz_tdiv_q_2exp(d, d, s);

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
      if (_verbose>3) gmp_printf("U=%Zd  V=%Zd  Qm=%Zd\n", U, V, Qm);
      mpz_mul(T1, U2m, V2m);
      mpz_mod(U2m, T1, n);
      mpz_mul(T1, V2m, V2m);
      mpz_submul_ui(T1, Qm, 2);
      mpz_mod(V2m, T1, n);
      if (_verbose>3) gmp_printf("  l  U2m=%Zd  V2m=%Zd\n", U2m, V2m);
      mpz_mul(T1, Qm, Qm); mpz_mod(Qm, T1, n);
      if (mpz_odd_p(d)) {
        /* Save T1 and T2 for later operations in this block */
        mpz_mul(T1, U2m, V);
        mpz_mul(T2, U2m, U);
        mpz_mul_si(T2, T2, D);
        if (_verbose>3) gmp_printf("      T1 %Zd  T2 %Zd\n", T1, T2);
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
  if (_verbose>3) gmp_printf("l0 U=%Zd  V=%Zd\n", U, V);
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


UV _GMP_trial_factor(mpz_t n, UV from_n, UV to_n)
{
  int small_n = 0;
  UV f;
  PRIME_ITERATOR(iter);

  if (mpz_cmp_ui(n, 4) < 0) {
    return (mpz_cmp_ui(n, 1) <= 0) ? 1 : 0;   /* 0,1 => 1   2,3 => 0 */
  }
  if ( (from_n <= 2) && mpz_even_p(n) )   return 2;

  if (from_n > to_n)
    croak("GMP_trial_factor from > to: %"UVuf" - %"UVuf, from_n, to_n);

  if (mpz_cmp_ui(n, to_n*to_n) < 0)
    small_n = 1;

  for (f = 2; f <= to_n; f = prime_iterator_next(&iter)) {
    if (small_n && mpz_cmp_ui(n, f*f) < 0) break;
    if (mpz_divisible_ui_p(n, f)) { prime_iterator_destroy(&iter); return f; }
  }
  prime_iterator_destroy(&iter);
  return 0;
}


int _GMP_is_prob_prime(mpz_t n)
{
  if (_GMP_trial_factor(n, 2, 400))
    return 0;
  if (mpz_cmp_ui(n, 400*400) <= 0)
    return 2;

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
  int prob_prime;

  prob_prime = _GMP_is_prob_prime(n);
  if (prob_prime != 1)
    return prob_prime;
  /*
   * BPSW indicated this is a prime, so we're pretty darn sure it is.  For
   * smaller numbers, try a quick proof using the Pocklington or BLS test.
   * We only do this for reasonably small numbers, and we're rather
   * half-hearted about the test, as it requires partially factoring n-1,
   * which may be easy or may be very hard.  They might not care at all
   * about this, so don't spend too much time.
   *
   * The idea is that they can call:
   *    is_prob_prime      BPSW
   *    is_prime           is_prob_prime + a little extra
   *    is_provable_prime  really prove it, which could take hours or days
   *
   * If we had better large-number factoring, we could extend the reasonable
   * range for this test.  In the long run, ECPP is the better method.  What
   * Pocklington/BLS gives us is a quite fast and very easy (in code) way to
   * prove primality for a lot of numbers in the 2^64 - 2^200 bit range (up
   * to about 50 digits).
   */
  if (mpz_sizeinbase(n, 2) <= 200) {
    prob_prime = _GMP_primality_bls(n, 1, 0);
    if (prob_prime < 0)
      return 0;
    if (prob_prime > 0)
      return 2;
  }
  return 1;
}

int _GMP_is_provable_prime(mpz_t n, char* prooftext)
{
  int prob_prime;

  prob_prime = _GMP_is_prob_prime(n);
  if (prob_prime != 1)
    return prob_prime;

  /* Miller-Rabin tests are fast, and the primality proving methods are
   * really, really, really slow for composites.  So run a few more MR tests
   * just in case.  These bases aren't special.
   */
  if (    (_GMP_miller_rabin_ui(n, 325) == 0)
       || (_GMP_miller_rabin_ui(n, 131) == 0)
       || (_GMP_miller_rabin_ui(n, 223) == 0)
       || (_GMP_miller_rabin_ui(n, 887) == 0)
     )
    return 0;

  prob_prime = _GMP_primality_bls(n, 0, &prooftext);
  if (prob_prime < 0)
    return 0;
  if (prob_prime > 0)
    return 2;
  return 1;
}

/*****************************************************************************/
/*          AKS.    This implementation is quite slow, but useful to have.   */

static int test_anr(UV a, mpz_t n, UV r, mpz_t* px, mpz_t* py, mpz_t* ptmp)
{
  int retval = 1;
  UV i, n_mod_r;
  mpz_t t;

  for (i = 0; i < r; i++)
    mpz_set_ui(px[i], 0);

  a %= r;
  mpz_set_ui(px[0], a);
  mpz_set_ui(px[1], 1);

  poly_mod_pow(py, px, ptmp, n, r, n);

  mpz_init(t);
  n_mod_r = mpz_mod_ui(t, n, r);
  if (n_mod_r >= r)  croak("n % r >= r ?!");
  mpz_sub_ui(t, py[n_mod_r], 1);
  mpz_mod(py[n_mod_r], t, n);
  mpz_sub_ui(t, py[0], a);
  mpz_mod(py[0], t, n);
  mpz_clear(t);

  for (i = 0; i < r; i++)
    if (mpz_sgn(py[i]))
      retval = 0;
  return retval;
}

int _GMP_is_aks_prime(mpz_t n)
{
  mpz_t sqrtn;
  mpz_t *px, *py, *pt;
  int retval;
  UV i, log2n, limit, rlimit, r, a;
  PRIME_ITERATOR(iter);

  if (mpz_cmp_ui(n, 4) < 0) {
    return (mpz_cmp_ui(n, 1) <= 0) ? 0 : 1;
  }

  if (mpz_perfect_power_p(n)) {
    return 0;
  }

  mpz_init(sqrtn);
  mpz_sqrt(sqrtn, n);
  /* limit should be floor(log2(n) * log2(n)).  This overcalculates it. */
  log2n = mpz_sizeinbase(n, 2);  /* ceil log2(n) */
  limit = log2n * log2n;

  if (_verbose>1) gmp_printf("# AKS checking order_r(%Zd) to %lu\n", n, (unsigned long) limit);

  /* Using a native r limits us to ~2000 digits in the worst case (r ~ log^5n)
   * but would typically work for 100,000+ digits (r ~ log^3n).  This code is
   * far too slow to matter either way. */

  for (r = 2; mpz_cmp_ui(n, r) >= 0; r = prime_iterator_next(&iter)) {
    if (mpz_divisible_ui_p(n, r)) {   /* r divides n.  composite. */
      prime_iterator_destroy(&iter);
      mpz_clear(sqrtn);
      return 0;
    }
    if (mpz_cmp_ui(sqrtn, r) < 0) {  /* no r <= sqrtn divides n.  prime. */
      prime_iterator_destroy(&iter);
      mpz_clear(sqrtn);
      return 1;
    }
    if (mpz_order_ui(r, n, limit) > limit)
      break;
  }
  prime_iterator_destroy(&iter);
  mpz_clear(sqrtn);

  if (mpz_cmp_ui(n, r) <= 0) {
    return 1;
  }

  rlimit = (UV) floor( sqrt(r-1) * log2n );

  if (_verbose) gmp_printf("# AKS %Zd.  r = %lu rlimit = %lu\n", n, (unsigned long) r, (unsigned long) rlimit);

  /* Create the three polynomials we will use */
  New(0, px, r, mpz_t);
  New(0, py, r, mpz_t);
  New(0, pt, r, mpz_t);
  if ( !px || !py || !pt )
    croak("allocation failure\n");
  for (i = 0; i < r; i++) {
    mpz_init(px[i]);
    mpz_init(py[i]);
    mpz_init(pt[i]);
  }

  retval = 1;
  for (a = 1; a <= rlimit; a++) {
    if (! test_anr(a, n, r, px, py, pt) ) {
      retval = 0;
      break;
    }
    if (_verbose>1) { printf("."); fflush(stdout); }
  }
  if (_verbose>1) { printf("\n"); fflush(stdout); };

  /* Free the polynomials */
  for (i = 0; i < r; i++) {
    mpz_clear(px[i]);
    mpz_clear(py[i]);
    mpz_clear(pt[i]);
  }
  Safefree(px);
  Safefree(py);
  Safefree(pt);

  return retval;
}

/*****************************************************************************/


/*
 * Lucas (1876): Given a completely factored n-1, if there exists an a s.t.
 *     a^(n-1) % n = 1
 *     a^((n-1/f) % n != 1 for ALL factors f of n-1
 * then n is prime.
 *
 * PPBLS:, given n-1 = A*B, A > sqrt(n), if we can find an a s.t.
 *     a^A % n = 1
 *     gcd(a^(A/f)-1,n) = 1 for ALL factors f of A
 * then n is prime.
 *
 * Generalized Pocklington: given n-1 = A*B, gcd(A,B)=1, A > sqrt(n), then if
 *     for each each factor f of A, there exists an a (1 < a < n-1) s.t.
 *         a^(n-1) % n = 1
 *         gcd(a^((n-1)/f)-1,n) = 1
 * then n is prime.
 *
 * BLS: given n-1 = A*B, factored A, s=B/2A r=B mod (2A), and an a, then if:
 *   - A is even, B is odd, and AB=n-1 (all implied by n = odd and the above),
 *   - n < (A+1) * (2*A*A + (r-1) * A + 1)
 *   - for each each factor f of A, there exists an a (1 < a < n-1) s.t.
 *     - a^(n-1) % n = 1
 *     - gcd(a^((n-1)/f)-1,n) = 1  for ALL factors f of A
 * then:
 *     if s = 0 or r*r - 8*s is not a perfect square
 *         n is prime
 *     else
 *         n is composite
 *
 * The generalized Pocklington test is also sometimes known as the
 * Pocklington-Lehmer test.  It's definitely an improvement over Lucas
 * since we only have to find factors up to sqrt(n), _and_ we can choose
 * a different 'a' value for each factor.  This is corollary 1 from BLS75.
 *
 * BLS is the Brillhart-Lehmer-Selfridge 1975 theorem 5 (see link below).
 * We can factor even less of n, and the test lets us kick out some
 * composites early, without having to test n-3 different 'a' values.
 *
 * Once we've found the factors of n-1 (or enough of them), verification
 * usually happens really fast.  a=2 works for most, and few seem to require
 * more than ~ log2(n).  However all but BLS75 require testing all integers
 * 1 < a < n-1 before answering in the negative, which is impractical.
 *
 *
 * AKS is not too hard to implement, but it's impractically slow.
 *
 * ECPP is very fast and definitely the best method for most numbers.
 *
 * BLS75:  http://www.ams.org/journals/mcom/1975-29-130/S0025-5718-1975-0384673-1/S0025-5718-1975-0384673-1.pdf
 *
 */

/* For these primality functions:
 *   -1   n is definitely composite
 *    1   n is definitely prime
 *    0   we can't decide
 *
 * You really should run is_prob_prime on n first, so we only have to run
 * these tests on numbers that are very probably prime.
 */

/* FIXME:
 *   (1) too much repetitious overhead code in these
 *   (2) way too much copy/paste between Pocklington and BLS
 */
#define PRIM_STACK_SIZE 128

#define primality_handle_factor(f, primality_func, factor_prob) \
  { \
    int f_prob_prime = _GMP_is_prob_prime(f); \
    if ( (f_prob_prime == 1) && (primality_func(f, do_quick, prooftextptr)) ) \
      f_prob_prime = 2; \
    if (f_prob_prime == 2) { \
      if (fsp >= PRIM_STACK_SIZE) { success = 0; } \
      else                        { mpz_init_set(fstack[fsp++], f); } \
      while (mpz_divisible_p(B, f)) { \
        mpz_mul(A, A, f); \
        mpz_divexact(B, B, f); \
      } \
    } else if ( (f_prob_prime == 0) || (factor_prob) ) { \
      if (msp >= PRIM_STACK_SIZE) { success = 0; } \
      else                        { mpz_init_set(mstack[msp++], f); } \
    } \
  }

#if 0
int _GMP_primality_pocklington(mpz_t n, int do_quick)
{
  mpz_t nm1, A, B, sqrtn, t, m, f;
  mpz_t mstack[PRIM_STACK_SIZE];
  mpz_t fstack[PRIM_STACK_SIZE];
  int msp = 0;
  int fsp = 0;
  int success = 1;

  mpz_init(nm1);
  mpz_sub_ui(nm1, n, 1);
  mpz_init_set_ui(A, 1);
  mpz_init_set(B, nm1);
  mpz_init(sqrtn);
  mpz_sqrt(sqrtn, n);
  mpz_init(m);
  mpz_init(f);
  mpz_init(t);

  { /* Pull small factors out */
    UV tf = 2;
    while ( (tf = _GMP_trial_factor(B, tf, 1000)) != 0 ) {
      if (fsp >= PRIM_STACK_SIZE) { success = 0; break; }
      mpz_init_set_ui(fstack[fsp++], tf);
      while (mpz_divisible_ui_p(B, tf)) {
        mpz_mul_ui(A, A, tf);
        mpz_divexact_ui(B, B, tf);
      }
      tf++;
    }
  }

  if (success) {
    mpz_set(f, B);
    primality_handle_factor(f, _GMP_primality_pocklington, 1);
  }

  while (success) {
    mpz_gcd(t, A, B);
    if ( (mpz_cmp(A, sqrtn) > 0) && (mpz_cmp_ui(t, 1) == 0) )
      break;
    success = 0;
    /* If the stack is empty, we have failed. */
    if (msp == 0)
      break;
    /* pop a component off the stack */
    mpz_set(m, mstack[--msp]); mpz_clear(mstack[msp]);

    /* Try to factor it without trying too hard */
    if (!success)  success = _GMP_power_factor(m, f);
    if (do_quick) {
      UV log2m = mpz_sizeinbase(m, 2);
      UV rounds = (log2m <= 64) ? 300000 : 300000 / (log2m-63);
      if (!success)  success = _GMP_pbrent_factor(m, f, 3, rounds);
    } else {
      if (!success)  success = _GMP_pbrent_factor(m, f, 3, 64*1024);
      if (!success)  success = _GMP_pbrent_factor(m, f, 5, 64*1024);
      if (!success)  success = _GMP_pbrent_factor(m, f, 7, 64*1024);
      if (!success)  success = _GMP_pbrent_factor(m, f,11, 64*1024);
      if (!success)  success = _GMP_pbrent_factor(m, f,13, 64*1024);
      if (!success)  success = _GMP_pbrent_factor(m, f, 1, 16*1024*1024);
      if (!success)  success = _GMP_ECM_FACTOR   (m, f, 12500, 4);
      if (!success)  success = _GMP_ECM_FACTOR   (m, f, 3125000, 10);
      if (!success)  success = _GMP_pbrent_factor(m, f, 3, 256*1024*1024);
      if (!success)  success = _GMP_ECM_FACTOR   (m, f, 400000000, 200);
    }
    /* If we couldn't factor m and the stack is empty, we've failed. */
    if ( (!success) && (msp == 0) )
      break;
    /* Put the two factors f and m/f into the stacks */
    primality_handle_factor(f, _GMP_primality_pocklington, 0);
    mpz_divexact(f, m, f);
    primality_handle_factor(f, _GMP_primality_pocklington, 0);
  }
  if (success) {
    int pcount, a;
    int const alimit = do_quick ? 200 : 10000;
    mpz_t p, ap;

    mpz_init(p);
    mpz_init(ap);

    for (pcount = 0; success && pcount < fsp; pcount++) {
      mpz_set(p, fstack[pcount]);
      success = 0;
      for (a = 2; !success && a <= alimit; a = next_small_prime(a)) {
        mpz_set_ui(ap, a);
        /* Does a^(n-1) % n = 1 ? */
        mpz_powm(t, ap, nm1, n);
        if (mpz_cmp_ui(t, 1) != 0)
          continue;
        /* Does gcd(a^((n-1)/f)-1,n) = 1 ? */
        mpz_divexact(B, nm1, p);
        mpz_powm(t, ap, B, n);
        mpz_sub_ui(t, t, 1);
        mpz_gcd(t, t, n);
        if (mpz_cmp_ui(t, 1) != 0)
          continue;
        success = 1;   /* We found an a for this p */
      }
    }
    mpz_clear(p);
    mpz_clear(ap);
  }
  while (msp-- > 0) {
    mpz_clear(mstack[msp]);
  }
  while (fsp-- > 0) {
    mpz_clear(fstack[fsp]);
  }
  mpz_clear(nm1);
  mpz_clear(A);
  mpz_clear(B);
  mpz_clear(sqrtn);
  mpz_clear(m);
  mpz_clear(f);
  mpz_clear(t);
  return success;
}
#endif

int _GMP_primality_bls(mpz_t n, int do_quick, char** prooftextptr)
{
  mpz_t nm1, A, B, t, m, f, r, s;
  mpz_t mstack[PRIM_STACK_SIZE];
  mpz_t fstack[PRIM_STACK_SIZE];
  int msp = 0;
  int fsp = 0;
  int success = 1;

  /* We need to do this for BLS */
  if (mpz_even_p(n)) return -1;

  mpz_init(nm1);
  mpz_sub_ui(nm1, n, 1);
  mpz_init_set_ui(A, 1);
  mpz_init_set(B, nm1);
  mpz_init(m);
  mpz_init(f);
  mpz_init(t);
  mpz_init(r);
  mpz_init(s);

#define TRIAL_B1 2000
  { /* Pull small factors out */
    PRIME_ITERATOR(iter);
    UV tf;
    for (tf = 2; tf < TRIAL_B1; tf = prime_iterator_next(&iter)) {
      if (mpz_cmp_ui(B, tf*tf) < 0) break;
      if (mpz_divisible_ui_p(B, tf)) {
        if (fsp >= PRIM_STACK_SIZE) { success = 0; break; }
        mpz_init_set_ui(fstack[fsp++], tf);
        do {
          mpz_mul_ui(A, A, tf);
          mpz_divexact_ui(B, B, tf);
        } while (mpz_divisible_ui_p(B, tf));
      }
    }
    prime_iterator_destroy(&iter);
  }

  if (success) {
    mpz_set(f, B);
    primality_handle_factor(f, _GMP_primality_bls, 1);
  }

  while (success) {
    mpz_mul_ui(t, A, 2);
    mpz_tdiv_q(s, B, t);
    mpz_mod(r, B, t);

    mpz_mul(m, t, A);    /* m = 2*A*A              */
    mpz_sub_ui(t, r, 1); /* t = r-1                */
    mpz_mul(t, t, A);    /* t = A*(r-1)            */
    mpz_add(m, m, t);    /* m = 2A^2 + A(r-1)      */
    mpz_add_ui(m, m, 1); /* m = 2A^2 + A(r-1) + 1  */
    mpz_add_ui(t, A, 1);
    mpz_mul(m, m, t);

    if (mpz_cmp(m, n) > 0)
      break;
#if 0
    /* Theorem 7 would let us do this, though we also have to do one extra
     * test at the end finding an a for B (treat the unfactored part B as if
     * it was the last factor). */
    mpz_mul(m, t, A);           /* m = 2*A*A               */
    mpz_sub_ui(t, r, TRIAL_B1); /* t = r-B1                */
    mpz_mul(t, t, A);           /* t = A*(r-B1)            */
    mpz_add(m, m, t);           /* m = 2A^2 + A(r-B1)      */
    mpz_add_ui(m, m, 1);        /* m = 2A^2 + A(r-B1) + 1  */
    mpz_mul_ui(t, A, TRIAL_B1);
    mpz_add_ui(t, t, 1);        /* t = B1*A+1              */
    mpz_mul(m, m, t);           /* m = (B1*A+1)*(2A^2+(r-B1)A+1) */
    if (mpz_cmp(m, n) > 0)
      break;
#endif
    success = 0;
    /* If the stack is empty, we have failed. */
    if (msp == 0)
      break;
    /* pop a component off the stack */
    mpz_set(m, mstack[--msp]); mpz_clear(mstack[msp]);

    /* Try to factor it without trying too hard */
    if (!success)  success = _GMP_power_factor(m, f);

    if (!success && mpz_cmp_ui(m, (unsigned long)(UV_MAX>>2)) < 0) {
      UV ui_m = mpz_get_ui(m);
      UV ui_factors[2];
      if (!mpz_cmp_ui(m, ui_m)) {
        success = racing_squfof_factor(ui_m, ui_factors, 200000)-1;
        if (success)
          mpz_set_ui(f, ui_factors[0]);
      }
    }
    if (do_quick) {
      /* These keep the time in the realm of 2ms per input number, which
       * means they adjust for size.  Success varies by size:
       *   99% of  70-bit inputs
       *   90% of  90-bit inputs
       *   58% of 120-bit inputs
       *   25% of 160-bit inputs
       *   10% of 190-bit inputs
       */
      UV log2m = mpz_sizeinbase(m, 2);
      UV brent_rounds = (log2m <= 64) ? 100000 : 100000 / (log2m-63);
      int final_B2 = 1000 * (150-(int)log2m);
      if (log2m < 70) brent_rounds *= 3;
      if (!success)  success = _GMP_pminus1_factor(m, f, 100, 1000);
      if (!success)  success = _GMP_pminus1_factor(m, f, 1000, 10000);
      if (!success && log2m < 80)  success = _GMP_ECM_FACTOR(m, f, 150, 5);
      if (!success)  success = _GMP_pbrent_factor(m, f, 3, brent_rounds);
      if (!success && final_B2 > 10000)  success = _GMP_pminus1_factor(m, f, 10000, final_B2);
    } else {
      if (!success)  success = _GMP_pminus1_factor(m, f, 100, 1000);
      if (!success)  success = _GMP_pminus1_factor(m, f, 1000, 10000);
      if (!success)  success = _GMP_pminus1_factor(m, f, 10000, 200000);
      if (!success)  success = _GMP_ECM_FACTOR(m, f, 500, 30);
      if (!success)  success = _GMP_ECM_FACTOR(m, f, 2000, 20);
      if (!success)  success = _GMP_pminus1_factor(m, f, 200000, 4000000);
      if (!success)  success = _GMP_ECM_FACTOR(m, f, 10000, 10);
      /* QS.  Uses lots of memory, but finds multiple factors quickly */
      if (!success && mpz_sizeinbase(m,10)>=30 && mpz_sizeinbase(m,2)<300) {
            mpz_t farray[66];
            int i, nfactors;
            for (i = 0; i < 66; i++)
              mpz_init(farray[i]);
            nfactors = _GMP_simpqs(m, farray);
            /* Insert all found factors */
            if (nfactors > 1) {
              success = 1;
              for (i = 0; i < nfactors; i++) {
                primality_handle_factor(farray[i], _GMP_primality_bls, 0);
              }
            }
            for (i = 0; i < 66; i++)
              mpz_clear(farray[i]);
            if (success)
              continue;
      }
      if (!success) {
        UV i;
        UV B1 = 10000;
        UV curves = 10;
        for (i = 1; i < 18; i++) {
          B1 *= 2;
          success = _GMP_ECM_FACTOR(m, f, B1, curves);
          if (success) break;
        }
      }
    }
    /* If we couldn't factor m and the stack is empty, we've failed. */
    if ( (!success) && (msp == 0) )
      break;
    /* Put the two factors f and m/f into the stacks, smallest first */
    mpz_divexact(m, m, f);
    if (mpz_cmp(m, f) < 0) { mpz_set(t, m); mpz_set(m, f); mpz_set(f, t); }
    primality_handle_factor(f, _GMP_primality_bls, 0);
    primality_handle_factor(m, _GMP_primality_bls, 0);
  }

  /* clear mstack since we don't care about it.  Use to hold a values. */
  while (msp-- > 0)
    mpz_clear(mstack[msp]);
  msp = 0;

  if (success > 0) {
    mpz_mul(t, r, r);
    mpz_submul_ui(t, s, 8);   /* t = r^2 - 8s */
    /* N is prime if and only if s=0 OR t not a perfect square */
    /* So, N is a composite if s != 0 and t is perfect square */
    if (mpz_sgn(s) != 0 && mpz_perfect_square_p(t))
      success = -1;
  }

  if (success > 0) {
    int pcount, a;
    int const alimit = do_quick ? 200 : 10000;
    mpz_t p, ap;

    mpz_init(p);
    mpz_init(ap);

    for (pcount = 0; success && pcount < fsp; pcount++) {
      PRIME_ITERATOR(iter);
      mpz_set(p, fstack[pcount]);
      success = 0;
      for (a = 2; !success && a <= alimit; a = prime_iterator_next(&iter)) {
        mpz_set_ui(ap, a);
        /* Does a^(n-1) % n = 1 ? */
        mpz_powm(t, ap, nm1, n);
        if (mpz_cmp_ui(t, 1) != 0)
          continue;
        /* Does gcd(a^((n-1)/f)-1,n) = 1 ? */
        mpz_divexact(B, nm1, p);
        mpz_powm(t, ap, B, n);
        mpz_sub_ui(t, t, 1);
        mpz_gcd(t, t, n);
        if (mpz_cmp_ui(t, 1) != 0)
          continue;
        success = 1;   /* We found an a for this p */
        mpz_init_set(mstack[msp++], ap);
      }
      prime_iterator_destroy(&iter);
    }
    mpz_clear(p);
    mpz_clear(ap);
  }
  if (success > 0 && prooftextptr != 0 && *prooftextptr != 0) {
    int i;
    if (fsp != msp) croak("Different f and a counts\n");
    *prooftextptr += gmp_sprintf(*prooftextptr, "%Zd :", n);
    for (i = 0; i < fsp; i++)
      *prooftextptr += gmp_sprintf(*prooftextptr, " %Zd", fstack[i]);
    *prooftextptr += gmp_sprintf(*prooftextptr, " :");
    for (i = 0; i < msp; i++)
      *prooftextptr += gmp_sprintf(*prooftextptr, " %Zd", mstack[i]);
    *prooftextptr += gmp_sprintf(*prooftextptr, "\n");
  }
  while (fsp-- > 0)
    mpz_clear(fstack[fsp]);
  while (msp-- > 0)
    mpz_clear(mstack[msp]);
  mpz_clear(nm1);
  mpz_clear(A);
  mpz_clear(B);
  mpz_clear(m);
  mpz_clear(f);
  mpz_clear(t);
  mpz_clear(r);
  mpz_clear(s);
  return success;
}


/* Modifies argument */
void _GMP_next_prime(mpz_t n)
{
  mpz_t d;
  UV m;

  /* small inputs */
  if (mpz_cmp_ui(n, 7) < 0) {
    if      (mpz_cmp_ui(n, 2) < 0) { mpz_set_ui(n, 2); }
    else if (mpz_cmp_ui(n, 3) < 0) { mpz_set_ui(n, 3); }
    else if (mpz_cmp_ui(n, 5) < 0) { mpz_set_ui(n, 5); }
    else                           { mpz_set_ui(n, 7); }
    return;
  }

  mpz_init(d);
  m = mpz_fdiv_q_ui(d, n, 30);

  if (m == 29) {
    mpz_add_ui(d, d, 1);
    m = 1;
  } else {
    m = next_wheel[m];
  }
  mpz_mul_ui(n, d, 30);
  mpz_add_ui(n, n, m);
  while (1) {
    if (_GMP_is_prob_prime(n))
      break;
    mpz_add_ui(n, n, wheel_advance[m]);
    m = next_wheel[m];
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

void _GMP_pn_primorial(mpz_t prim, UV n)
{
  UV p = 2;
  PRIME_ITERATOR(iter);

  mpz_set_ui(prim, 1);
  while (n--) {
    mpz_mul_ui(prim, prim, p);
    p = prime_iterator_next(&iter);
  }
  prime_iterator_destroy(&iter);
}
void _GMP_primorial(mpz_t prim, mpz_t n)
{
  UV p = 2;
  PRIME_ITERATOR(iter);

  mpz_set_ui(prim, 1);
  while (mpz_cmp_ui(n, p) >= 0) {
    mpz_mul_ui(prim, prim, p);
    p = prime_iterator_next(&iter);
  }
  prime_iterator_destroy(&iter);
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


void _GMP_lcm_of_consecutive_integers(UV B, mpz_t m)
{
  UV p, p_power, pmin;
  PRIME_ITERATOR(iter);

  /* For each prime, multiply m by p^floor(log B / log p), which means
   * raise p to the largest power e such that p^e <= B.
   */
  mpz_set_ui(m, 1);
  if (B >= 2) {
    p_power = 2;
    while (p_power <= B/2)
      p_power *= 2;
    mpz_mul_ui(m, m, p_power);
  }
  p = prime_iterator_next(&iter);
  while (p <= B) {
    pmin = B/p;
    if (p > pmin)
      break;
    p_power = p*p;
    while (p_power <= pmin)
      p_power *= p;
    mpz_mul_ui(m, m, p_power);
    p = prime_iterator_next(&iter);
  }
  while (p <= B) {
    mpz_mul_ui(m, m, p);
    p = prime_iterator_next(&iter);
  }
  prime_iterator_destroy(&iter);
}


int _GMP_pminus1_factor(mpz_t n, mpz_t f, UV B1, UV B2)
{
  mpz_t a, savea, t;
  UV q, saveq, j;
  PRIME_ITERATOR(iter);

  TEST_FOR_2357(n, f);
  if (B1 < 7) return 0;

  mpz_init(a);
  mpz_init(savea);
  mpz_init(t);

  if (_verbose>2) gmp_printf("# p-1 trying %Zd (B1=%lu B2=%lu)\n", n, (unsigned long)B1, (unsigned long)B2);

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
  j = 1;
  mpz_set_ui(a, 2);
  mpz_set_ui(savea, 2);
  saveq = 2;
  /* We could wrap this in a loop trying a few different a values, in case
   * the current one ended up going to 0. */
  q = 2;
  mpz_set_ui(t, 1);
  while (q <= B1) {
    UV k, kmin;
    k = q;
    kmin = B1/q;
    while (k <= kmin)
      k *= q;
    mpz_mul_ui(t, t, k);        /* Accumulate powers for a */
    if ( (j++ % 32) == 0) {
      mpz_powm(a, a, t, n);     /* a=a^(k1*k2*k3*...) mod n */
      if ( !mpz_sgn(a) )
        goto end_fail;
      mpz_sub_ui(t, a, 1);
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
  if ( !mpz_sgn(a) )
    goto end_fail;
  mpz_sub_ui(t, a, 1);
  mpz_gcd(f, t, n);
  if (mpz_cmp(f, n) == 0) {
    /* We found multiple factors.  Loop one at a time. */
    prime_iterator_setprime(&iter, saveq);
    mpz_set(a, savea);
    for (q = saveq; q <= B1; q = prime_iterator_next(&iter)) {
      UV k = q;
      UV kmin = B1/q;
      while (k <= kmin)
        k *= q;
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
   * See Montgomery 1987, p249-250 for what one _should_ do.
   * This is the standard continuation which replaces the powmods in stage 1
   * with two mulmods, with a GCD every 64 primess (no backtracking
   * implemented).  This is quite a bit faster than stage 1.
   * We quickly precalculate a few of the prime gaps, and lazily cache others
   * up to a gap of 222.  That's enough for a B2 value of 189 million.  We
   * still work above that, we just won't cache the value.
   */
  if (B2 > B1) {
    mpz_t b, bm, bmdiff;
    mpz_t precomp_bm[111];
    int   is_precomp[111] = {0};

    mpz_init(bmdiff);
    mpz_init_set(bm, a);
    mpz_init_set_ui(b, 1);

    /* Set the first 20 differences */
    mpz_powm_ui(bmdiff, bm, 2, n);
    mpz_init_set(precomp_bm[0], bmdiff);
    is_precomp[0] = 1;
    for (j = 1; j < 20; j++) {
      mpz_mul(bmdiff, bmdiff, bm);
      mpz_mul(bmdiff, bmdiff, bm);
      mpz_tdiv_r(bmdiff, bmdiff, n);
      mpz_init_set(precomp_bm[j], bmdiff);
      is_precomp[j] = 1;
    }

    mpz_powm_ui(a, a, q, n );

    j = 1;
    while (q <= B2) {
      UV lastq = q;
      UV qdiff;

      q = prime_iterator_next(&iter);
      qdiff = (q - lastq) / 2 - 1;
      if (qdiff >= 111) {
        mpz_powm_ui(bmdiff, bm, q-lastq, n);  /* Big gap */
      } else if (is_precomp[qdiff]) {
        mpz_set(bmdiff, precomp_bm[qdiff]);
      } else {
        mpz_powm_ui(bmdiff, bm, q-lastq, n);
        mpz_init_set(precomp_bm[qdiff], bmdiff);
        is_precomp[qdiff] = 1;
      }
      mpz_mul(t, a, bmdiff);
      mpz_tdiv_r(a, t, n);
      if ( !mpz_sgn(a) )
        break;
      mpz_sub_ui(t, a, 1);
      mpz_mul(t, b, t);
      mpz_tdiv_r(b, t, n);
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
      mpz_clear(s); mpz_clear(m); return 1;
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

   mpz_mod_ui(t1, n, 4);
   if (mpz_cmp_ui(t1, 3))
     croak("Incorrect call to shanks\n");

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
      3*5*7*11, 3*5*7, 3*5*11, 3*5, 3*7*11, 3*7, 5*7*11, 5*7,
      3*11,     3,     5*11,   5,   7*11,   7,   11,     1   };
   const size_t sz_mul = sizeof(multipliers)/sizeof(multipliers[0]);
   size_t i;
   int result;
   UV nmod4;
   mpz_t nm, t;

   TEST_FOR_2357(n, f);
   mpz_init(nm);
   mpz_init(t);
   mpz_set_ui(f, 1);
   nmod4 = mpz_tdiv_r_ui(nm, n, 4);

   for (i = 0; i < sz_mul; i++) {
      UV mult = multipliers[i];
      /* Only use multipliers where n*m = 3 mod 4 */
      if (nmod4 == (mult % 4))
        continue;
      /* Only run when 64*m^3 < n */
      mpz_set_ui(t, mult);
      mpz_pow_ui(t, t, 3);
      mpz_mul_ui(t, t, 64);
      if (mpz_cmp(t, n) >= 0)
        continue;
      /* Run with this multiplier */
      mpz_mul_ui(nm, n, mult);
      result = shanks_mult(nm, f);
      if (result == -1)
        break;
      if ( (result == 1) && (mpz_cmp_ui(f, mult) != 0) ) {
        unsigned long gcdf = mpz_gcd_ui(NULL, f, mult);
        mpz_divexact_ui(f, f, gcdf);
        if (mpz_cmp_ui(f, 1) > 0)
          break;
      }
   }
   mpz_clear(t);
   mpz_clear(nm);
   return (mpz_cmp_ui(f, 1) > 0);
}

/* See if n is a perfect power */
int _GMP_power_factor(mpz_t n, mpz_t f)
{
  if (mpz_perfect_power_p(n)) {
    unsigned long k;

    mpz_set_ui(f, 1);
    for (k = 2; mpz_sgn(f); k++) {
      if (mpz_root(f, n, k))
        return 1;
    }
    /* GMP says it's a perfect power, but we couldn't find an integer root? */
  }
  return 0;
}
