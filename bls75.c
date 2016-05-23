#include <string.h>
#include <gmp.h>
#include "ptypes.h"

#include "bls75.h"
#include "primality.h"
#include "prime_iterator.h"
#include "small_factor.h"
#include "factor.h"
#include "simpqs.h"
#include "ecm.h"
#define _GMP_ECM_FACTOR(n, f, b1, ncurves) \
   _GMP_ecm_factor_projective(n, f, b1, 0, ncurves)
#include "utility.h"

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
 * BLS T5: given n-1 = A*B, factored A, s=B/2A r=B mod (2A), and an a, then if:
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
 * BLS75 theorem 7 is the final n-1 theorem and takes into account any
 * knowledge that the remaining factor is not below a threshold B.  Since
 * we do initial trial division this helps.  It is usually of only small
 * benefit.
 *
 *
 * AKS is not too hard to implement, but it's impractically slow.
 *
 * ECPP is very fast and definitely the best method for most numbers.
 *
 * BLS75:  http://www.ams.org/journals/mcom/1975-29-130/S0025-5718-1975-0384673-1/S0025-5718-1975-0384673-1.pdf
 *
 */

/* Like all the primality functions:
 *   2 = definitely prime, 1 = maybe prime, 0 = definitely composite
 *
 * You really should run is_prob_prime on n first, so we only have to run
 * these tests on numbers that are very probably prime.
 */


static int try_factor(mpz_t f, mpz_t n, int effort)
{
  int success = 0;
  UV log2n = mpz_sizeinbase(n, 2);

  if (!success && mpz_cmp_ui(n, (unsigned long)(UV_MAX>>5)) < 0) {
    UV ui_n = mpz_get_ui(n);
    UV ui_factors[2];
    if (!mpz_cmp_ui(n, ui_n)) {
      success = racing_squfof_factor(ui_n, ui_factors, 200000)-1;
      if (success)
        mpz_set_ui(f, ui_factors[0]);
    }
  }

  if (effort >= 1) {
    if (!success)  success = _GMP_pminus1_factor(n, f, 1000, 10000);
    //if (!success)  success = _GMP_pminus1_factor(n, f, 4000, 4000);
    //if (!success)  success = _GMP_pminus1_factor(n, f, 2000, 2000);
    //if (!success)  success = _GMP_pminus1_factor(n, f, 1000, 1000);
    //if (!success)  success = _GMP_pminus1_factor(n, f, 400, 400);
  }

  if (!success)  success = (int)power_factor(n, f);

  if (!success && effort == 2) {
    UV brent_rounds = (log2n <= 64) ? 100000 : 100000 / (log2n-63);
    int final_B2 = 1000 * (150-(int)log2n);
    if (log2n < 70) brent_rounds *= 3;
    if (!success && log2n < 80)  success = _GMP_ECM_FACTOR(n, f, 150, 5);
    if (!success)  success = _GMP_pbrent_factor(n, f, 3, brent_rounds);
    if (!success && final_B2 > 10000)  success = _GMP_pminus1_factor(n, f, 10000, final_B2);
  }

  if (!success && effort >= 3) {
    if (!success)  success = _GMP_pminus1_factor(n, f, 10000, 200000);
    if (!success)  success = _GMP_ECM_FACTOR(n, f, 500, 30);
    if (!success)  success = _GMP_ECM_FACTOR(n, f, 2000, 20);
  }
  if (!success && effort >= 5 && log2n > 170) {
    UV B1 = (log2n > 2500) ? 10000000 : 4000 * log2n;
    if (!success)  success = _GMP_pminus1_factor(n, f, B1, 20*B1);
    /* To head off expensive QS, do these early */
    if (!success && log2n > 210) success = _GMP_ECM_FACTOR(n, f, 20000, 10);
    if (!success && log2n > 240) success = _GMP_ECM_FACTOR(n, f, 40000, 10);
    if (!success && log2n > 240) success = _GMP_ECM_FACTOR(n, f, 80000,  5);
    if (!success && log2n > 270) success = _GMP_ECM_FACTOR(n, f,160000, 20);
  }

  return success;
}
static int try_factor2(mpz_t f, mpz_t n, int effort)
{
  int success = 0;

  if (!success && effort >= 4) {
    if (!success)  success = _GMP_pminus1_factor(n, f, 200000, 4000000);
    if (!success)  success = _GMP_ECM_FACTOR(n, f, 10000, 10);
  }
  if (!success && effort >= 5) {
    UV i;
    UV ecm_B1 = 10000;
    UV curves = 10;
    if (_GMP_is_prob_prime(n)) croak("Internal error in BLS75\n");
    for (i = 1; i < 18 && !success; i++) {
      if ((4+i) > (UV)effort) break;
      ecm_B1 *= 2;
      success = _GMP_ECM_FACTOR(n, f, ecm_B1, curves);
    }
  }
  return success;
}

/* F*R = n, F is factored part, R is remainder */
static void small_factor(mpz_t F, mpz_t R, UV B1)
{
  PRIME_ITERATOR(iter);
  UV tf;
  for (tf = 2; tf < B1; tf = prime_iterator_next(&iter)) {
    if (mpz_cmp_ui(R, tf*tf) < 0) break;
    if (mpz_divisible_ui_p(R, tf)) {
      do {
        mpz_mul_ui(F, F, tf);
        mpz_divexact_ui(R, R, tf);
      } while (mpz_divisible_ui_p(R, tf));
    }
  }
  prime_iterator_destroy(&iter);
}


/* FIXME:
 *   (1) too much repetitious overhead code in these
 *   (2) way too much copy/paste between Pocklington and BLS
 *   (3) Pocklington has code rotted, so fix before using
 */
#define ADD_TO_STACK(val, stack, cur, max) \
  if (cur == max) \
    Renew(stack, max += 10, mpz_t); \
  mpz_init_set( stack[cur++], val );
#define ADD_TO_STACK_UI(val, stack, cur, max) \
  if (cur == max) \
    Renew(stack, max += 10, mpz_t); \
  mpz_init_set_ui( stack[cur++], val );


#define test_ui_factor(f, R, F, fstack, fsp, fsmax) \
  { \
    if (mpz_divisible_ui_p(R, f)) { \
      ADD_TO_STACK_UI( f, fstack, fsp, fsmax ); \
      do { \
        mpz_mul_ui(F, F, f); \
        mpz_divexact_ui(R, R, f); \
      } while (mpz_divisible_ui_p(R, f)); \
    } \
  }

#define primality_handle_factor(f, R, F, fstack, fsp, fsmax, primality_func, factor_prob) \
  { \
    int f_prob_prime = _GMP_is_prob_prime(f); \
    if ( (f_prob_prime == 1) && (primality_func(f, effort, prooftextptr) == 2) ) \
      f_prob_prime = 2; \
    if (f_prob_prime == 2) { \
      ADD_TO_STACK( f, fstack, fsp, fsmax ); \
      while (mpz_divisible_p(R, f)) { \
        mpz_mul(F, F, f); \
        mpz_divexact(R, R, f); \
      } \
    } else if ( (f_prob_prime == 0) || (factor_prob) ) { \
      ADD_TO_STACK( f, mstack, msp, msmax ); \
    } \
  }

#define INNER_QS_FACTOR(qn, primality_func) \
  { \
    mpz_t farray[66]; \
    int i, nfactors; \
    for (i = 0; i < 66; i++)  mpz_init(farray[i]); \
    nfactors = _GMP_simpqs(qn, farray); \
    /* Insert all found factors */ \
    if (nfactors > 1) { \
      success = 1; \
      for (i = 0; i < nfactors; i++) \
        primality_handle_factor(farray[i], B, A, fstack, fsp, fsmax, primality_func, 0); \
    } \
    for (i = 0; i < 66; i++)  mpz_clear(farray[i]); \
    if (success) \
      continue; \
  }

/* Sort factors found from largest to smallest, but 2 must be at start. */
static void sort_and_trim_factors(int* fsp, mpz_t* fstack)
{
  int i, j;
  for (i = 2; i < *fsp; i++)
    for (j = i; j > 1 && mpz_cmp(fstack[j-1], fstack[j]) < 0; j--)
      mpz_swap( fstack[j-1], fstack[j] );
  for (i = 2; i < *fsp; i++) { /* Remove any duplicate factors */
    if (mpz_cmp(fstack[i], fstack[i-1]) == 0) {
      for (j = i+1; j < *fsp; j++)
        mpz_set(fstack[j-1], fstack[j]);
      *fsp -= 1;
    }
  }
}



/******************************************************************************/



static int bls_theorem5_limit(mpz_t n, mpz_t A, mpz_t B,
                              mpz_t t, mpz_t y, mpz_t r, mpz_t s)
{
  mpz_mul(t, A, B);
  mpz_add_ui(t, t, 1);
  if (mpz_cmp(t, n) != 0) croak("BLS75 internal error: A*B != n-1\n");

  mpz_mul_ui(t, A, 2);
  mpz_tdiv_qr(s, r, B, t);

  mpz_mul(y, t, A);     /* y = 2*A*A              */
  mpz_sub_ui(t, r, 1);  /* t = r-1                */
  mpz_mul(t, t, A);     /* t = A*(r-1)            */
  mpz_add(y, y, t);     /* y = 2A^2 + A(r-1)      */
  mpz_add_ui(y, y, 1);  /* y = 2A^2 + A(r-1) + 1  */
  mpz_add_ui(t, A, 1);  /* t = A+1                */
  mpz_mul(y, y, t);     /* y = (A+1)*(2A^2+(r-1)A+1) */

  return (mpz_cmp(n, y) < 0) ? 1 : 0;
}
static int bls_theorem17_limit(mpz_t n, mpz_t A, mpz_t B,
                               mpz_t t, mpz_t y, mpz_t r, mpz_t s)
{
  mpz_mul(t, A, B);
  mpz_sub_ui(t, t, 1);
  if (mpz_cmp(t, n) != 0) croak("BLS75 internal error: A*B != n+1\n");

  mpz_mul_ui(t, A, 2);
  mpz_tdiv_qr(s, r, B, t);

  mpz_mul(y, t, A);     /* y = 2*A*A              */
  mpz_add(y, y, A);     /* y = 2A^2 + A           */
  mpz_mul(t, A, r);     /* t = rA                 */
  mpz_sub(y, y, t);     /* y = 2A^2 + A(1-r)      */
  mpz_add_ui(y, y, 1);  /* y = 2A^2 + A(1-r) + 1  */
  mpz_sub_ui(t, A, 1);  /* t = A-1                */
  mpz_mul(y, y, t);     /* y = (A-1)*(2A^2+(1-r)A+1) */

  return (mpz_cmp(n, y) < 0) ? 1 : 0;
}

static int bls_corollary11_limit(mpz_t n, mpz_t F1, mpz_t F2, UV B,
                                 mpz_t m, mpz_t t)
{
  mpz_mul(m, F1, F2);
  mpz_mul_ui(m, m, B);
  mpz_mul_ui(m, m, B);
  mpz_mul_ui(m, m, B);
  mpz_tdiv_q_2exp(m, m, 1);
  mpz_mul(t, m, F1);
  if (mpz_cmp(n, t) >= 0) {
    mpz_mul(t, m, F2);
    if (mpz_cmp(n, t) >= 0) {
      return 0;  /* We need more factoring or a higher B */
    }
  }
  return 1;
}
static int bls_theorem20_limit(mpz_t n, mpz_t R1, mpz_t F1, mpz_t F2,
                               UV B, UV m,
                               mpz_t t, mpz_t g, mpz_t r, mpz_t s)
{
  mpz_tdiv_q_2exp(t, F2, 1);
  mpz_tdiv_qr(s, r, R1, t);

  mpz_mul_ui(t, F1, B);
  mpz_add_ui(g, t, 1);

  mpz_mul_ui(t, F2, B);
  mpz_sub_ui(t, t, 1);

  if (mpz_cmp(t, g) > 0)  mpz_set(g, t);

  mpz_mul(t, F1, F2);
  mpz_tdiv_q_2exp(t, t, 1);
  mpz_mul_ui(t, t, B);
  mpz_mul_ui(t, t, B);
  mpz_add_ui(s, t, 1);    /* s = B1*B2*F1*F2/2+1 */

  mpz_mul(g, g, s);
  if (mpz_cmp(n, g) < 0) {
    mpz_set_ui(s, 0);     /* Use s to signal whether we must test for m. */
    return 1;
  }

  mpz_mul(t, F1, F2);
  mpz_mul_ui(t, t, m);
  mpz_tdiv_q_2exp(t, t, 1);
  mpz_mul(r, r, F1);           /* r *= F1;  t += r,  r /= F1; */
  mpz_add(t, t, r);
  mpz_divexact(r, r, F1);
  mpz_add_ui(t, t, 1);         /* t = m * F1 * F2/2 + r * F1 + 1 */
  mpz_mul(g, s, t);
  mpz_set_ui(s, 1);

  return (mpz_cmp(n, g) < 0) ? 1 : 0;
}


/******************************************************************************/


/* (I) For each prime p_i dividing F1 [N-1 = F1R1] there exists an a_i
 *     such that N is a psp base a_i and gcd(A_i^{(N-1)/p_i}-1,N) = 1.
 */
static int _verify_cond_I_p(mpz_t n, mpz_t pi, mpz_t ap, mpz_t t, int alimit, char* pspcache)
{
  int a, success = 0;
  PRIME_ITERATOR(iter);

  for (a = 2; !success && a <= alimit; a = prime_iterator_next(&iter)) {
     int psp = -1;
     mpz_set_ui(ap, a);

     if (pspcache)  psp = pspcache[a];
     if (psp == -1) {
       mpz_sub_ui(t, n, 1);
       mpz_powm(t, ap, t, n);
       psp = (mpz_cmp_ui(t, 1) == 0);
     }
     if (pspcache)  pspcache[a] = psp;
     if (!psp)
       continue;

     mpz_sub_ui(t, n, 1);
     mpz_divexact(t, t, pi);
     mpz_powm(t, ap, t, n);
     mpz_sub_ui(t, t, 1);
     mpz_gcd(t, t, n);
     if (mpz_cmp_ui(t, 1) != 0)
       continue;

     success = 1;   /* We found an a for this p */
  }
  prime_iterator_destroy(&iter);
  return success;
}

/* (III) For each prime q_i dividing F2 [N+1 = F2R2] there exists a Lucas
 *       sequence U_k with discriminant D for which D/N = -1, N divides
 *       U_{N+1} and gcd(U_{(N+1)/q_i},N} = 1.
 */
static int _verify_cond_III_q(mpz_t n, mpz_t qi, IV* lp, IV* lq)
{
  int rval = 0;
  IV d, p, q;
  mpz_t U, V, k, t1, t2;

  mpz_init(U);  mpz_init(V);  mpz_init(k);  mpz_init(t1);  mpz_init(t2);

  for (q = 2; q < 1000; q++) {
    p = (q % 2) ? 2 : 1;
    d = p*p - 4*q;
    mpz_set_si(t1, d);
    if (mpz_jacobi(t1, n) != -1)
      continue;

    /* we have a d/p/q where d = -1.  Check the Lucas sequence. */

    mpz_add_ui(k, n, 1);
    lucas_seq(U, V, n, p, q, k,    t1, t2);
    if (mpz_sgn(U) == 0) {
      mpz_divexact(k, k, qi);
      lucas_seq(U, V, n, p, q, k,  t1, t2);
      mpz_gcd(t1, V, n);
      if (mpz_cmp_ui(t1, 1) == 0) {
        rval = 2;
        if (lp) *lp = p;
        if (lq) *lq = q;
        break;
      }
    }
  }
  mpz_clear(U);  mpz_clear(V);  mpz_clear(k);  mpz_clear(t1);  mpz_clear(t2);
  return rval;
}



/******************************************************************************/



int _GMP_primality_bls_nm1(mpz_t n, int effort, char** prooftextptr)
{
  mpz_t nm1, A, B, t, m, f, r, s;
  mpz_t* fstack;
  mpz_t* mstack;
  int fsp = 0, fsmax = 10;
  int msp = 0, msmax = 10;
  int success = 1;
  UV B1 = (mpz_sizeinbase(n,10) > 1000) ? 100000 : 2000;

  /* We need to do this for BLS */
  if (mpz_even_p(n)) return 0;

  mpz_init(nm1);
  mpz_sub_ui(nm1, n, 1);
  mpz_init_set_ui(A, 1);
  mpz_init_set(B, nm1);
  mpz_init(m);
  mpz_init(f);
  mpz_init(t);
  mpz_init(r);
  mpz_init(s);

  New(0, fstack, fsmax, mpz_t);
  New(0, mstack, msmax, mpz_t);

  { /* Pull small factors out */
    PRIME_ITERATOR(iter);
    UV tf;
    for (tf = 2; tf < B1; tf = prime_iterator_next(&iter)) {
      if (mpz_cmp_ui(B, tf*tf) < 0) break;
      test_ui_factor(tf, B, A, fstack, fsp, fsmax);
    }
    prime_iterator_destroy(&iter);
  }

  if (success) {
    mpz_set(f, B);
    primality_handle_factor(f, B, A, fstack, fsp, fsmax, _GMP_primality_bls_nm1, 1);
  }

  while (success) {

    if (bls_theorem5_limit(n, A, B, t, m, r, s))
      break;

    success = 0;
    /* If the stack is empty, we have failed. */
    if (msp == 0)
      break;
    /* pop a component off the stack */
    mpz_set(m, mstack[--msp]); mpz_clear(mstack[msp]);

    success = try_factor(f, m, effort);

    /* QS.  Uses lots of memory, but finds multiple factors quickly */
    if (!success && effort >= 5 &&
        mpz_sizeinbase(m,10) >= 30 && mpz_sizeinbase(m,10) <= 90) {
      if (effort > 5 || (effort == 5 && mpz_sizeinbase(m,10) < 55) ) {
        INNER_QS_FACTOR(m, _GMP_primality_bls_nm1);
      }
    }

    if (!success)
      success = try_factor2(f, m, effort);

    /* If we couldn't factor m and the stack is empty, we've failed. */
    if ( (!success) && (msp == 0) )
      break;
    /* Put the two factors f and m/f into the stacks, smallest first */
    mpz_divexact(m, m, f);
    if (mpz_cmp(m, f) < 0)
      mpz_swap(m, f);
    primality_handle_factor(f, B, A, fstack, fsp, fsmax, _GMP_primality_bls_nm1, 0);
    primality_handle_factor(m, B, A, fstack, fsp, fsmax, _GMP_primality_bls_nm1, 0);
  }

  /* clear mstack since we don't care about it.  Use to hold a values. */
  while (msp-- > 0)
    mpz_clear(mstack[msp]);
  msp = 0;

  sort_and_trim_factors(&fsp, fstack);

  /* Shrink to smallest set and verify conditions. */
  if (success > 0) {
    int i;
    mpz_set_ui(A, 1);
    mpz_set(B, nm1);
    for (i = 0; i < fsp; i++) {
      if (bls_theorem5_limit(n, A, B, t, m, r, s))
        break;
      do {
        mpz_mul(A, A, fstack[i]);
        mpz_divexact(B, B, fstack[i]);
      } while (mpz_divisible_p(B, fstack[i]));
    }
    /* Delete any extra factors */
    while (i < fsp)
      mpz_clear(fstack[--fsp]);
    /* Verify Q[0] = 2 */
    if (mpz_cmp_ui(fstack[0], 2) != 0)
      croak("BLS75 internal error: 2 not at start of fstack");
    /* Verify conditions */
    success = 0;
    if (bls_theorem5_limit(n, A, B, t, m, r, s)) {
      mpz_mul(t, r, r);
      mpz_submul_ui(t, s, 8);   /* t = r^2 - 8s */
      /* N is prime if and only if s=0 OR t not a perfect square */
      success = (mpz_sgn(s) == 0 || !mpz_perfect_square_p(t))  ?  1  :  -1;
    }
  }

  if (success > 0) {
    int pcount, a;
    int const alimit = (effort <= 2) ? 200 : 10000;
    char afermat[10000+1];
    mpz_t p, ap;

    mpz_init(p);
    mpz_init(ap);

    /* Cache result that doesn't depend on factor */
    for (a = 0; a <= alimit; a++)  afermat[a] = -1;

    for (pcount = 0; success && pcount < fsp; pcount++) {
      success = _verify_cond_I_p(n, fstack[pcount], ap, t, alimit, afermat);
      if (success) {
        ADD_TO_STACK( ap, mstack, msp, msmax );
      }
    }
    /* If we could not find 'a' values, then we should return 1 (maybe prime)
     * since we did not perform an exhaustive search.  It would be quite
     * unusual to find a prime that didn't have an 'a' in the first 10,000
     * primes, but it could happen.  It's a "dubiously prime" :) */
    if (!success && get_verbose_level() > 0)
      printf("N-1 factored but failed to prove.  Perhaps composite.\n");
    mpz_clear(p);
    mpz_clear(ap);
  }
  if (success > 0 && prooftextptr != 0) {
    int i;
    char *proofstr, *proofptr;
    int curprooflen = (*prooftextptr == 0) ? 0 : strlen(*prooftextptr);
    int myprooflen = (5 + mpz_sizeinbase(n, 10)) * (2 + fsp + msp) + 200;

    if (fsp != msp) croak("Different f and a counts\n");
    New(0, proofstr, myprooflen + curprooflen + 1, char);
    proofptr = proofstr;
    proofptr += gmp_sprintf(proofptr, "Type BLS5\nN  %Zd\n", n);
    /* Q[0] is always 2 */
    for (i = 1; i < fsp; i++)
      proofptr += gmp_sprintf(proofptr, "Q[%d]  %Zd\n", i, fstack[i]);
    /* A[i] only printed if not 2 */
    for (i = 0; i < msp; i++)
      if (mpz_cmp_ui(mstack[i], 2) != 0)
        proofptr += gmp_sprintf(proofptr, "A[%d]  %Zd\n", i, mstack[i]);
    proofptr += gmp_sprintf(proofptr, "----\n");
    /* Set or prepend */
    if (*prooftextptr) {
      proofptr += gmp_sprintf(proofptr, "\n");
      strcat(proofptr, *prooftextptr);
      Safefree(*prooftextptr);
    }
    *prooftextptr = proofstr;
  }
  while (fsp-- > 0)
    mpz_clear(fstack[fsp]);
  Safefree(fstack);
  while (msp-- > 0)
    mpz_clear(mstack[msp]);
  Safefree(mstack);
  mpz_clear(nm1);
  mpz_clear(A);
  mpz_clear(B);
  mpz_clear(m);
  mpz_clear(f);
  mpz_clear(t);
  mpz_clear(r);
  mpz_clear(s);
  if (success < 0) return 0;
  if (success > 0) return 2;
  return 1;
}


int _GMP_primality_bls_np1(mpz_t n, int effort, char** prooftextptr)
{
  mpz_t np1, A, B, t, m, f, r, s;
  mpz_t* fstack;
  mpz_t* mstack;
  int fsp = 0, fsmax = 10;
  int msp = 0, msmax = 10;
  int success = 1;
  UV B1 = (mpz_sizeinbase(n,10) > 1000) ? 100000 : 2000;

  /* We need to do this for BLS */
  if (mpz_even_p(n)) return 0;

  mpz_init(np1);
  mpz_add_ui(np1, n, 1);
  mpz_init_set_ui(A, 1);
  mpz_init_set(B, np1);
  mpz_init(m);
  mpz_init(f);
  mpz_init(t);
  mpz_init(r);
  mpz_init(s);

  New(0, fstack, fsmax, mpz_t);
  New(0, mstack, msmax, mpz_t);

  { /* Pull small factors out */
    PRIME_ITERATOR(iter);
    UV tf;
    for (tf = 2; tf < B1; tf = prime_iterator_next(&iter)) {
      if (mpz_cmp_ui(B, tf*tf) < 0) break;
      test_ui_factor(tf, B, A, fstack, fsp, fsmax);
    }
    prime_iterator_destroy(&iter);
  }

  if (success) {
    mpz_set(f, B);
    primality_handle_factor(f, B, A, fstack, fsp, fsmax, _GMP_primality_bls_np1, 1);
  }

  while (success) {

    if (bls_theorem17_limit(n, A, B, t, m, r, s))
      break;

    success = 0;
    /* If the stack is empty, we have failed. */
    if (msp == 0)
      break;
    /* pop a component off the stack */
    mpz_set(m, mstack[--msp]); mpz_clear(mstack[msp]);

    success = try_factor(f, m, effort);

    /* QS.  Uses lots of memory, but finds multiple factors quickly */
    if (!success && effort >= 5 &&
        mpz_sizeinbase(m,10) >= 30 && mpz_sizeinbase(m,10) <= 90) {
      if (effort > 5 || (effort == 5 && mpz_sizeinbase(m,10) < 55) ) {
        INNER_QS_FACTOR(m, _GMP_primality_bls_np1);
      }
    }

    if (!success)
      success = try_factor2(f, m, effort);

    /* If we couldn't factor m and the stack is empty, we've failed. */
    if ( (!success) && (msp == 0) )
      break;
    /* Put the two factors f and m/f into the stacks, smallest first */
    mpz_divexact(m, m, f);
    if (mpz_cmp(m, f) < 0)
      mpz_swap(m, f);
    primality_handle_factor(f, B, A, fstack, fsp, fsmax, _GMP_primality_bls_np1, 0);
    primality_handle_factor(m, B, A, fstack, fsp, fsmax, _GMP_primality_bls_np1, 0);
  }

  /* clear mstack since we don't care about it.  Use to hold a values. */
  while (msp-- > 0)
    mpz_clear(mstack[msp]);
  msp = 0;

  sort_and_trim_factors(&fsp, fstack);

  /* Shrink to smallest set and verify conditions. */
  if (success > 0) {
    int i;
    mpz_set_ui(A, 1);
    mpz_set(B, np1);
    for (i = 0; i < fsp; i++) {
      if (bls_theorem17_limit(n, A, B, t, m, r, s))
        break;
      do {
        mpz_mul(A, A, fstack[i]);
        mpz_divexact(B, B, fstack[i]);
      } while (mpz_divisible_p(B, fstack[i]));
    }
    /* Delete any extra factors */
    while (i < fsp)
      mpz_clear(fstack[--fsp]);
    /* Verify Q[0] = 2 */
    if (mpz_cmp_ui(fstack[0], 2) != 0)
      croak("BLS75 internal error: 2 not at start of fstack");
    /* Verify conditions */
    success = 0;
    if (bls_theorem17_limit(n, A, B, t, m, r, s)) {
      mpz_mul(t, r, r);
      mpz_submul_ui(t, s, 8);   /* t = r^2 - 8s */
      /* N is prime if and only if s=0 OR t not a perfect square */
      success = (mpz_sgn(s) == 0 || !mpz_perfect_square_p(t))  ?  1  :  -1;
    }
  }

  if (success > 0) {
    int pcount;
    for (pcount = 0; success && pcount < fsp; pcount++) {
      success = _verify_cond_III_q(n, fstack[pcount], 0, 0);
    }
    if (!success && get_verbose_level() > 0)
      printf("N+1 factored but failed to prove.  Perhaps composite.\n");
  }

  /* TODO: Proof text */

  while (fsp-- > 0)
    mpz_clear(fstack[fsp]);
  Safefree(fstack);
  while (msp-- > 0)
    mpz_clear(mstack[msp]);
  Safefree(mstack);
  mpz_clear(np1);
  mpz_clear(A);
  mpz_clear(B);
  mpz_clear(m);
  mpz_clear(f);
  mpz_clear(t);
  mpz_clear(r);
  mpz_clear(s);
  if (success < 0) return 0;
  if (success > 0) return 2;
  return 1;
}

/******************************************************************************/

#define PRINT_PCT 0

/* Use Theorem 20 */
int bls75_hybrid(mpz_t n, int effort, char** prooftextptr)
{
  mpz_t nm1, np1, F1, F2, R1, R2;
  mpz_t r, s, t, u;
  mpz_t *f1stack, *f2stack, *mstack;
  int f1sp = 0, f1smax = 10;
  int f2sp = 0, f2smax = 10;
  int msp = 0, msmax = 10;
  int pcount;
  int success = 1;
  UV B1 = (mpz_sizeinbase(n,10) > 1000) ? 100000 : 10000;
  UV m = 8000;   /* m should be less than B1 */
#if PRINT_PCT
  double trial_pct, prime_pct, nm1_pct, np1_pct;
#endif

  /* We need to do this for BLS */
  if (mpz_even_p(n)) return 0;

  mpz_init(nm1); mpz_sub_ui(nm1, n, 1);
  mpz_init(np1); mpz_add_ui(np1, n, 1);

  mpz_init_set_ui(F1, 1); mpz_init_set(R1, nm1);
  mpz_init_set_ui(F2, 1); mpz_init_set(R2, np1);

  mpz_init(r);
  mpz_init(s);
  mpz_init(u);
  mpz_init(t);

  New(0, f1stack, f1smax, mpz_t);
  New(0, f2stack, f2smax, mpz_t);
  New(0, mstack, msmax, mpz_t);

  { /* Pull small factors out */
    PRIME_ITERATOR(iter);
    UV tf;
    for (tf = 2; tf < B1; tf = prime_iterator_next(&iter)) {
      /* Divisibility on R1 and R2 is faster than the page 635 suggestion */
      if (mpz_cmp_ui(R1, tf*tf) >= 0)
        test_ui_factor(tf, R1, F1, f1stack, f1sp, f1smax);
      if (mpz_cmp_ui(R2, tf*tf) >= 0)
        test_ui_factor(tf, R2, F2, f2stack, f2sp, f2smax);
    }
    prime_iterator_destroy(&iter);
  }

#if PRINT_PCT
  trial_pct = (100.0 * (mpz_sizeinbase(F1,2) + mpz_sizeinbase(F2,2))) / (mpz_sizeinbase(nm1,2) + mpz_sizeinbase(np1,2));
#endif

  if (0) { /* Check primality of R2 immediately */
    /* This leads to less time factoring n-1, but also takes recurse time */
    mpz_set(t, R2);
    primality_handle_factor(t, R2, F2, f2stack, f2sp, f2smax, bls75_hybrid, 1);
    /* But remove from work stack */
    while (msp-- > 0) mpz_clear(mstack[msp]);
    msp = 0;
  }

#if PRINT_PCT
  prime_pct = (100.0 * (mpz_sizeinbase(F1,2) + mpz_sizeinbase(F2,2))) / (mpz_sizeinbase(nm1,2) + mpz_sizeinbase(np1,2));
#endif

  /* TODO: We'd like to cut down the time on low-effort runs.  Ideas:
   *
   *  1) put off recursive primality checks until needed, run from smallest.
   *
   *  2) make more progressive effort levels, and run both n-1 and n+1 in a
   *     loop over them.  So we do both at effort 1, then bump up to effort 2
   *     for both, etc.
   *
   */
  if (1) {
  //if (mpz_cmp_ui(R1,1) == 0 || mpz_cmp_ui(R2,1) == 0) {
    /* More factoring for N-1 */
    mpz_set(t, R1);
    primality_handle_factor(t, R1, F1, f1stack, f1sp, f1smax, bls75_hybrid, 1);
    while (success) {
      if (bls_theorem20_limit(n, R1, F1, F2, B1, m, t, u, r, s))
        break;
      if (mpz_cmp_ui(R1,1) == 0 || mpz_cmp_ui(R2,1) == 0) croak("bad");
      success = 0;
      if (msp == 0) break; /* If the stack is empty, we have failed. */
      /* pop a component off the stack */
      mpz_set(u, mstack[--msp]); mpz_clear(mstack[msp]);

      success = try_factor(t, u, effort) || try_factor2(t, u, effort);

      /* If we couldn't factor u and the stack is empty, we've failed. */
      if (!success && msp == 0) break;
      /* Put the two factors f and u/f into the stacks, smallest first */
      mpz_divexact(u, u, t);
      if (mpz_cmp(u, t) < 0)
        mpz_swap(u, t);
      primality_handle_factor(t, R1, F1, f1stack, f1sp, f1smax, bls75_hybrid, 0);
      primality_handle_factor(u, R1, F1, f1stack, f1sp, f1smax, bls75_hybrid, 0);
    }
    success = 1;
    while (msp-- > 0) mpz_clear(mstack[msp]);
    msp = 0;
  }
#if PRINT_PCT
  nm1_pct = (100.0 * (mpz_sizeinbase(F1,2) + mpz_sizeinbase(F2,2))) / (mpz_sizeinbase(nm1,2) + mpz_sizeinbase(np1,2));
#endif
  if (1) {
  //if (mpz_cmp_ui(R1,1) == 0 || mpz_cmp_ui(R2,1) == 0) {
    /* More factoring for N+1 */
    mpz_set(t, R2);
    primality_handle_factor(t, R2, F2, f2stack, f2sp, f2smax, bls75_hybrid, 1);
    while (success) {
      if (bls_theorem20_limit(n, R1, F1, F2, B1, m, t, u, r, s))
        break;
      if (mpz_cmp_ui(R1,1) == 0 || mpz_cmp_ui(R2,1) == 0) croak("bad");
      success = 0;
      if (msp == 0) break; /* If the stack is empty, we have failed. */
      /* pop a component off the stack */
      mpz_set(u, mstack[--msp]); mpz_clear(mstack[msp]);

      success = try_factor(t, u, effort) || try_factor2(t, u, effort);

      /* If we couldn't factor u and the stack is empty, we've failed. */
      if (!success && msp == 0) break;
      /* Put the two factors f and u/f into the stacks, smallest first */
      mpz_divexact(u, u, t);
      if (mpz_cmp(u, t) < 0)
        mpz_swap(u, t);
      primality_handle_factor(t, R2, F2, f2stack, f2sp, f2smax, bls75_hybrid, 0);
      primality_handle_factor(u, R2, F2, f2stack, f2sp, f2smax, bls75_hybrid, 0);
    }
    success = 1;
    while (msp-- > 0) mpz_clear(mstack[msp]);
    msp = 0;
  }
#if PRINT_PCT
  np1_pct = (100.0 * (mpz_sizeinbase(F1,2) + mpz_sizeinbase(F2,2))) / (mpz_sizeinbase(nm1,2) + mpz_sizeinbase(np1,2));
  printf("%6.2f .. %6.2f .. %6.2f .. %6.2f\n", trial_pct, prime_pct, nm1_pct, np1_pct);  fflush(stdout);
#endif

  /* If we've completely factored either one, quickly handle it. */
  if (!mpz_cmp_ui(R1,1)) { /* Theorem 1 */
    for (pcount = 0; success && pcount < f1sp; pcount++)
      success = _verify_cond_I_p(n, f1stack[pcount], u, t, 10000, 0);
    goto end_hybrid;
  }
  if (!mpz_cmp_ui(R2,1)) { /* Theorem 13 */
    for (pcount = 0; success && pcount < f2sp; pcount++)
      success = _verify_cond_III_q(n, f2stack[pcount], 0, 0);
    goto end_hybrid;
  }

  mpz_mul(t, F1, R1); if (mpz_cmp(nm1, t) != 0) croak("Bad n-1 factor");
  mpz_mul(t, F2, R2); if (mpz_cmp(np1, t) != 0) croak("Bad n+1 factor");

  sort_and_trim_factors(&f1sp, f1stack);
  sort_and_trim_factors(&f2sp, f2stack);

  /* Check N < B^3 F1*F2*F2/2  or  N < B^3 F1*F1*F2/2 */
  if (success > 0) 
    success = bls_theorem20_limit(n, R1, F1, F2, B1, m, t, u, r, s);

  /* Check lambda divisibility if needed */
  if (success > 0 && mpz_sgn(s)) {
    UV lambda;
    mpz_mul(t, F1, F2);
    mpz_tdiv_q_2exp(t, t, 1);
    mpz_mul(u, r, F1);
    mpz_add_ui(u, u, 1);
    for (lambda = 0; success > 0 && lambda < m; lambda++) {
      if (mpz_cmp(n, u) >= 0)
        break;
      if (lambda > 0 || mpz_sgn(r))
        if (mpz_divisible_p(n, u))
          success = -1;
      mpz_add(u, u, t);
    }
  }

  /* Verify (I)   page 623 */
  for (pcount = 0; success > 0 && pcount < f1sp; pcount++) {
    success = _verify_cond_I_p(n, f1stack[pcount], u, t, 10000, 0);
  }

  /* Verify (II)  page 625 (II applied to R1) */
  if (success > 0)
    success = _verify_cond_I_p(n, R1, u, t, 100, 0);

  /* Verify (III) page 631 */
  for (pcount = 0; success > 0 && pcount < f2sp; pcount++) {
    success = _verify_cond_III_q(n, f2stack[pcount], 0, 0);
  }
  /* Verify (IV)  page 633 (III applied to R2) */
  if (success > 0)
    success = _verify_cond_III_q(n, R2, 0, 0);

end_hybrid:
  while (f1sp-- > 0)
    mpz_clear(f1stack[f1sp]);
  Safefree(f1stack);
  while (f2sp-- > 0)
    mpz_clear(f2stack[f2sp]);
  Safefree(f2stack);
  while (msp-- > 0)
    mpz_clear(mstack[msp]);
  Safefree(mstack);
  mpz_clear(nm1); mpz_clear(np1);
  mpz_clear(F1);  mpz_clear(F2);
  mpz_clear(R1);  mpz_clear(R2);
  mpz_clear(r);
  mpz_clear(s);
  mpz_clear(u);
  mpz_clear(t);
  if (success < 0) return 0;
  if (success > 0) return 2;
  return 1;
}


/* Given an n where we're factored n-1 down to p, check BLS theorem 3 */
int _GMP_primality_bls_3(mpz_t n, mpz_t p, UV* reta)
{
  mpz_t nm1, m, t, t2;
  int rval = 0;

  if (reta) *reta = 0;
  if (mpz_cmp_ui(n, 2) <= 0 || mpz_even_p(n) || mpz_even_p(p))
    return 0;                 /* n is <= 2, n is even, or p is even */
  if (!_GMP_is_prob_prime(p))
    return 0;                 /* p is not a probable prime */

  mpz_init(nm1);  mpz_init(m);  mpz_init(t);  mpz_init(t2);
  mpz_sub_ui(nm1, n, 1);
  mpz_divexact(m, nm1, p);
  mpz_mul(t, m, p);
  if (mpz_cmp(nm1, t) != 0)
    goto end_bls3;           /* m*p != n+1 */

  mpz_mul_ui(t, p, 2);
  mpz_add_ui(t, t, 1);
  mpz_sqrt(t2, n);
  if (mpz_cmp(t, t2) <= 0)
    goto end_bls3;           /* 2p+1 <= sqrt(n) */

  {
    /* N-1 = mp, p is an odd probable prime, and 2p+1 > sqrt(n).
     * Now find an 'a' where a^(n-1)/2 = -1 mod n, a^(m/2) != -1 mod n. */
    PRIME_ITERATOR(iter);
    UV const alimit = 1000;
    UV a;
    for (a = 2; a <= alimit; a = prime_iterator_next(&iter)) {
      mpz_set_ui(t2, a);
      mpz_divexact_ui(t, m, 2);
      mpz_powm(t, t2, t, n);       /* a^(m/2) mod n */
      if (mpz_cmp(t, nm1) == 0)
        continue;
      mpz_divexact_ui(t, nm1, 2);
      mpz_powm(t, t2, t, n);       /* a^((n-1)/2) mod n */
      if (mpz_cmp(t, nm1) != 0)
        continue;
      rval = 2;
      if (reta) *reta = a;
      break;
    }
    prime_iterator_destroy(&iter);
  }

end_bls3:
  mpz_clear(nm1);  mpz_clear(m);  mpz_clear(t);  mpz_clear(t2);
  return rval;
}

/* Given an n where we're factored n+1 down to f, check BLS theorem 15 */
int _GMP_primality_bls_15(mpz_t n, mpz_t f, IV* lp, IV* lq)
{
  mpz_t np1, m, t, t2;
  int rval = 0;

  if (lp) *lp = 0;
  if (lq) *lq = 0;
  if (mpz_cmp_ui(n, 2) <= 0 || mpz_even_p(n) || mpz_even_p(f))
    return 0;                 /* n is <= 2, n is even, or f is even */
  if (!_GMP_is_prob_prime(f))
    return 0;                 /* f is not a probable prime */

  mpz_init(np1);  mpz_init(m);  mpz_init(t);  mpz_init(t2);
  mpz_add_ui(np1, n, 1);
  mpz_divexact(m, np1, f);
  mpz_mul(t, m, f);
  if (mpz_cmp(np1, t) != 0)
    goto end_bls15;           /* m*f != n+1 */

  mpz_mul_ui(t, f, 2);
  mpz_sub_ui(t, t, 1);
  mpz_sqrt(t2, n);
  if (mpz_cmp(t, t2) <= 0)
    goto end_bls15;           /* 2f-1 <= sqrt(n) */

  {
    /* N+1 = mf, f is an odd probable prime, and 2f-1 > sqrt(n).
     * Now find a Lucas sequence V_k with discriminant D s.t. D/N = -1
     * where N divides V_(N+1)/2 and N does not divide V_m/2. */
    IV d, p, q;
    mpz_t U, V, k;

    mpz_init(U);  mpz_init(V);  mpz_init(k);

    /* Primo gave me the idea of this p/q selection method */
    for (q = 2; q < 1000; q++) {
      p = (q % 2) ? 2 : 1;
      d = p*p - 4*q;
      mpz_set_si(t, d);
      if (mpz_jacobi(t, n) != -1)
        continue;
      /* we have a d/p/q where d = -1.  Check the Lucas sequences. */
      mpz_divexact_ui(k, m, 2);
      lucas_seq(U, V, n, p, q, k,    t, t2);
      if (mpz_sgn(V) != 0) {
        mpz_divexact_ui(k, np1, 2);
        lucas_seq(U, V, n, p, q, k,    t, t2);
        if (mpz_sgn(V) == 0) {
          rval = 2;
          if (lp) *lp = p;
          if (lq) *lq = q;
          break;
        }
      }
    }
    mpz_clear(U);  mpz_clear(V);  mpz_clear(k);
  }

end_bls15:
  /* Somehow there is a tester getting 0 for LQ */
  if (rval && lq && *lq < 2) croak("Internal error in BLS15\n");
  mpz_clear(np1);  mpz_clear(m);  mpz_clear(t);  mpz_clear(t2);
  return rval;
}

/* Given an n, try using BLS75 theorem 15, N+1 = mq.
 * Note: this does _not_ prove n is prime!  If it returns 1, then we have
 * found a q/D that satisfy theorem 15, but we leave proving q for the caller.
 */
int _GMP_primality_bls_np1_split(mpz_t n, int effort, mpz_t q, IV* lp, IV* lq)
{
  mpz_t np1, m, f, sqrtn, t;
  int success = 1;
  UV B1 = 2000;

  /* We need to do this for BLS */
  if (mpz_even_p(n)) return 0;

  mpz_init(np1);  mpz_init(m);  mpz_init(f);  mpz_init(sqrtn);  mpz_init(t);
  mpz_add_ui(np1, n, 1);
  mpz_set_ui(m, 1);
  mpz_set(q, np1);
  mpz_sqrt(sqrtn, n);

  small_factor(m, q, B1);

  while (success) {
    success = 0;
    mpz_mul_ui(t, q, 2);
    mpz_sub_ui(t, t, 1);
    if (mpz_cmp(t, sqrtn) <= 0)
      break;
    if (_GMP_is_prob_prime(q)) {
      success = 1;
      break;
    }
    success = try_factor(f, q, effort);
    if (!success)
      success = try_factor2(f, q, effort);
    if (success) {
      mpz_divexact(q, q, f);
      if (mpz_cmp(q, f) < 0)
        mpz_swap(q, f);
      mpz_mul(m, m, f);
    }
  }

  if (success)
    success = _GMP_primality_bls_15(n, q, lp, lq);

  mpz_clear(np1);
  mpz_clear(m);
  mpz_clear(f);
  mpz_clear(sqrtn);
  mpz_clear(t);
  return success;
}

/* Given an n, try using BLS75 theorem 3, N-1 = mp. */
int _GMP_primality_bls_nm1_split(mpz_t n, int effort, mpz_t p, UV *reta)
{
  mpz_t nm1, m, f, sqrtn, t;
  int success = 1;
  UV B1 = 2000;

  /* We need to do this for BLS */
  if (mpz_even_p(n)) return 0;

  mpz_init(nm1);  mpz_init(m);  mpz_init(f);  mpz_init(sqrtn);  mpz_init(t);
  mpz_sub_ui(nm1, n, 1);
  mpz_set_ui(m, 1);
  mpz_set(p, nm1);
  mpz_sqrt(sqrtn, n);

  small_factor(m, p, B1);

  while (success) {
    success = 0;
    mpz_mul_ui(t, p, 2);
    mpz_add_ui(t, t, 1);
    if (mpz_cmp(t, sqrtn) <= 0)
      break;
    if (_GMP_is_prob_prime(p)) {
      success = 1;
      break;
    }
    success = try_factor(f, p, effort);
    if (!success)
      success = try_factor2(f, p, effort);
    if (success) {
      mpz_divexact(p, p, f);
      if (mpz_cmp(p, f) < 0)
        mpz_swap(p, f);
      mpz_mul(m, m, f);
    }
  }

  if (success)
    success = _GMP_primality_bls_3(n, p, reta);

  mpz_clear(nm1);
  mpz_clear(m);
  mpz_clear(f);
  mpz_clear(sqrtn);
  mpz_clear(t);
  return success;
}
