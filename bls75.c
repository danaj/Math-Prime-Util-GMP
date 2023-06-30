#include <string.h>
#include <gmp.h>
#include "ptypes.h"
static const int ev = 0;

#include "bls75.h"
#include "primality.h"
#include "lucas_seq.h"
#include "prime_iterator.h"
#include "pbrent63.h"
#include "squfof126.h"
#include "factor.h"
#include "simpqs.h"
#include "ecm.h"
#define _GMP_ECM_FACTOR(n, f, b1, ncurves) \
   _GMP_ecm_factor_projective(n, f, b1, 0, ncurves)
#include "utility.h"

static const int maxuibits = (sizeof(unsigned long) > BITS_PER_WORD)
                           ? BITS_PER_WORD : (8*sizeof(unsigned long));

/*
 *
 * Theorems where we are factoring n-1.
 *
 *
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
 * BLS T5: given n-1 = A*B, factored A, s=B/2A r=B mod (2A), then if:
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
 * BLS T7: given n-1 = A*B, factored A, a B1 where all factors of B are >= B1,
 *         s=B/2A r=B mod (2A), then if:
 *   - A is even, B is odd, and AB=n-1 (all implied by n = odd and the above),
 *   - n < (B1*A+1) * (2*A*A + (r-B1)*A + 1)
 *   - for each f in {B, factors of A}, there exists an a (1 < a < n-1) s.t.
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
 * BLS T5 is the Brillhart-Lehmer-Selfridge 1975 theorem 5 (see link below).
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
 * benefit, but it costs nothing.
 *
 *
 * AKS is not too hard to implement, but it's impractically slow.
 *
 * ECPP is very fast and definitely the best method for most numbers.
 *
 * APR-CL is very practical for numbers of a few hundred digits.
 *
 * BLS75:  http://www.ams.org/journals/mcom/1975-29-130/S0025-5718-1975-0384673-1/S0025-5718-1975-0384673-1.pdf
 *
 */

/*
 * Theorems where we are factoring n+1.
 *
 * Corollary 8 is analagous to Corollary 1 (generalized Pocklington).

 * Theorem 15 is analagous to Theorem 3.  These are simple proofs for the case
 * of n+1 (theorem 15) / n-1 (theorem 3) having a large prime factor.  These
 * are used by ECPP.  T3 is slightly better than the original Proth theorem.

 * Theorem 17 is analagous to Theorem 5.
 * Theorem 19 is analagous to Theorem 7.
 *
 * Theorem 20 is a hybrid, combining theorems 7 and 19.
 */

/* Like all the primality functions:
 *   2 = definitely prime, 1 = maybe prime, 0 = definitely composite
 *
 * You really should run is_prob_prime on n first, so we only have to run
 * these tests on numbers that are very probably prime.
 */




static int tfe(mpz_t f, mpz_t n, int effort)
{
  int success = 0;
  UV log2n = mpz_sizeinbase(n, 2);

  if (mpz_cmp_ui(n,3) <= 0) {
    mpz_set(f,n);
    return 1;
  }

  if (effort == 0) {
    if (!success && log2n <= 63) success = pbrent63(n, f, 1000000);
    if (!success) success = (int)power_factor(n, f);
    if (success) return success;
  }
  /* TODO: Use tinyqs effectively here, e.g. stage 3 for 50-90 bit */

  switch (effort) {
    case 0: success = _GMP_pminus1_factor(n, f, 400, 400);
            break;
    case 1: success = _GMP_pminus1_factor(n, f, 1000, 11000);
            break;
    case 2: { UV brent_rounds = (log2n <= 64) ? 100000 : 100000 / (log2n-63);
              int final_B2 = 1000 * (150-(int)log2n);
              if (log2n < 70) brent_rounds *= 3;
              if (log2n < 80)
                success = _GMP_ECM_FACTOR(n, f, 150, 5);
              if (!success)
                success = _GMP_pbrent_factor(n, f, 3, brent_rounds);
              if (!success && final_B2 > 11000)
                success = _GMP_pminus1_factor(n, f, 10000, final_B2);
            } break;
    case 3: success = _GMP_pminus1_factor(n, f, 20000, 200000);
            if (!success && log2n <= 80) success = squfof126(n, f, 2000000);
            break;
    case 4: success = _GMP_ECM_FACTOR(n, f, 500, 30);
            break;
    case 5: success = _GMP_ECM_FACTOR(n, f, 2000, 20);
            break;
    case 6: if (log2n > 170) {
              UV B1 = (log2n > 2500) ? 10000000 : 4000 * log2n;
              success = _GMP_pminus1_factor(n, f, B1, 20*B1);
            } break;
    case 7: if (log2n > 210) { success = _GMP_ECM_FACTOR(n, f, 20000, 10); }
            break;
    case 8: success = _GMP_cheb_factor(n, f, 50000, 0);
            break;
    case 9: if (log2n > 240) { success = _GMP_ECM_FACTOR(n, f, 40000, 10); }
            break;
    case 10:if (log2n > 240) { success = _GMP_ECM_FACTOR(n, f, 80000,  5); }
            break;
    case 11:if (log2n > 270) { success = _GMP_ECM_FACTOR(n, f,160000, 20); }
            break;

    /* QS for sizes 30-90 digits */
    case 20:
    case 21:{ UV log10n = mpz_sizeinbase(n, 10);
              if (log10n >= 30 && log10n <= ((effort == 20) ? 54 : 90)) {
                mpz_t farray[66];
                int i, nfactors;
                for (i = 0; i < 66; i++)  mpz_init(farray[i]);
                nfactors = _GMP_simpqs(n, farray);
                /* TODO: Return all factors */
                if (nfactors > 1) {
                 success = 1;
                  mpz_set(f, farray[nfactors-1]);   /* Return largest */
                }
                for (i = 0; i < 66; i++)  mpz_clear(farray[i]);
              }
            } break;

    case 30:success = _GMP_pminus1_factor(n, f, 200000, 4000000);
            break;

    case 40:
    case 41:
    case 42:
    case 43:
    case 44:
    case 45:
    case 46:
    case 47:
    case 48:
    case 49:
    case 50:
    case 51:
    case 52:
    case 53:
    case 54:
    case 55:
    case 56:
    case 57:{ UV B1 = UVCONST(1) << (13 + (effort-40));
              success = _GMP_ECM_FACTOR(n, f, B1, 10);
            } break;

    default: break;
  }

  /* if (success) printf("   bls75 factored %lu-bit at effort %d\n", log2n, effort); */
  return success;
}


typedef struct {
  int    cur;
  int    max;
  mpz_t* stack;
} fstack_t;

#define FACTOR_STACK(name)  fstack_t name = {0, 0, 0}

static int nstack(fstack_t* s) { return s->cur; }
static void push_fstack(fstack_t* s, mpz_t v) {
  if (s->stack == 0)    New(0,s->stack, s->max = 10, mpz_t);
  if (s->cur == s->max) Renew(s->stack, s->max += 10, mpz_t);
  mpz_init_set(s->stack[(s->cur)++], v);
}
static void push_fstack_ui(fstack_t* s, unsigned long v) {
  if (s->stack == 0)    New(0,s->stack, s->max = 10, mpz_t);
  if (s->cur == s->max) Renew(s->stack, s->max += 10, mpz_t);
  mpz_init_set_ui(s->stack[(s->cur)++], v);
}
static void pop_fstack(mpz_t rv, fstack_t* s) {
  mpz_set(rv, s->stack[--(s->cur)]);
  mpz_clear(s->stack[s->cur]);
}
static void clear_fstack(fstack_t* s) {
  while (s->cur > 0)
    mpz_clear(s->stack[--(s->cur)]);
}
static void destroy_fstack(fstack_t* s) {
  clear_fstack(s);
  Safefree(s->stack);
  s->stack = 0;
}

static void factor_out(mpz_t R, mpz_t F, mpz_t v) {
  int ndiv = mpz_remove(R, R, v);
  while (ndiv-- > 0)
    mpz_mul(F, F, v);
}
static void factor_out_ui(mpz_t R, mpz_t F, unsigned long v) {
  while (mpz_divisible_ui_p(R, v)) {
    mpz_mul_ui(F, F, v);
    mpz_divexact_ui(R, R, v);
  }
}

static int factor_test_ui(unsigned long f, mpz_t R, mpz_t F, fstack_t* s) {
  if (mpz_divisible_ui_p(R, f)) {
    push_fstack_ui(s, f);
    factor_out_ui(R, F, f);
    return 1;
  }
  return 0;
}

typedef int (*bls_func_t)(mpz_t, int, char**);
typedef int (*limit_func_t)(mpz_t, mpz_t, mpz_t, UV, mpz_t,mpz_t,mpz_t,mpz_t);

static void handle_factor(mpz_t f, mpz_t R, mpz_t F,
                          fstack_t* sf, fstack_t* sm,
                          int effort, char** prtext,
                          int push_if_probable,
                          bls_func_t func) {
  int pr = _GMP_BPSW(f);
  if (pr == 1) { /* Try to prove */
    if (effort > 1 || mpz_sizeinbase(f,2) < 200) {
      pr = (*func)(f, effort, prtext);
    }
  }
  if (pr == 2) {
    push_fstack(sf, f);
    factor_out(R, F, f);
  } else if (pr == 0 || push_if_probable) {
    push_fstack(sm, f);
  }
}
static void handle_factor2(mpz_t f, mpz_t R, mpz_t F,
                           fstack_t* sf, fstack_t* sp, fstack_t* sm,
                           int effort, char** prtext,
                           bls_func_t func) {
  int pr = _GMP_BPSW(f);
  if (pr == 1) { /* Try to prove */
    pr = (*func)(f, effort, prtext);
  }
  if (pr == 0) {
    push_fstack(sm, f);
  } else if (pr == 2) {
    push_fstack(sf, f);
    factor_out(R, F, f);
  } else {
    /* Save actual proof for later, but for now assume we can do it */
    push_fstack(sp, f);
    factor_out(R, F, f);
  }
}

static void trim_factors(mpz_t F, mpz_t R, mpz_t n, mpz_t n_pm_one, UV B, fstack_t* fs, limit_func_t func, mpz_t t, mpz_t m, mpz_t r, mpz_t s) {
  int i;
  if (ev > 1) gmp_printf("Starting trim with F %Zd R %Zd\n", F, R);
  if (fs->cur > 1) {
    if (mpz_odd_p(n_pm_one)) croak("n-1 / n+1 isn't even in trim_factors");
    mpz_set_ui(F, 1);
    mpz_set(R, n_pm_one);
    for (i = 0; i < fs->cur; i++) {
      /* 2 always goes into F and all factors < B go into F */
      if (i > 0 && mpz_cmp_ui(fs->stack[i],B) >= 0
                && func(n, F, R, B, t, m, r, s))
        break;
      factor_out(R, F, fs->stack[i]);
      if (ev > 1) gmp_printf("   %Zd -> F %Zd R %Zd\n", fs->stack[i], F, R);
    }
    /* Remove excess factors */
    while (i < fs->cur)
      pop_fstack(t, fs);
  }
  /* Verify Q[0] = 2 */
  if (mpz_cmp_ui(fs->stack[0], 2) != 0)
    croak("BLS75 internal error: 2 not at start of fstack");
  /* r and s have been set by func */
}

static int numcmp(const void *av, const void *bv)
  { return mpz_cmp(*(const mpz_t*)av, *(const mpz_t*)bv); }

/* Remove duplicates, sort smallest factors < B, reverse sort rest. */
static void _sort_and_remove_dup_factors(int* fsp, mpz_t* fstack, UV B)
{
  int i, j;
  if (B < 3) B = 3;
  /* 1.  Sort everything smallest to largest */
  qsort(fstack, *fsp, sizeof(mpz_t), numcmp);
  /* 2.  Remove any duplicate factors */
  for (i = 1; i < *fsp; i++) {
    if (mpz_cmp(fstack[i], fstack[i-1]) == 0) {
      for (j = i+1; j < *fsp; j++)
        mpz_set(fstack[j-1], fstack[j]);
      *fsp -= 1;
      mpz_clear(fstack[*fsp]);
    }
  }
  /* 3.  Find first factor >= B */
  for (i = 1; i < *fsp; i++)
    if (mpz_cmp_ui(fstack[i], B) >= 0)
      break;
  /* 4.  Reverse all remaining factors */
  if ((*fsp-i) > 1) {
    int ihead, itail;
    for (ihead = i, itail = *fsp-1;  ihead < itail;  ihead++, itail--)
      mpz_swap(fstack[ihead], fstack[itail]);
  }
}
static void fstack_tidy(fstack_t* s, UV B) {
  _sort_and_remove_dup_factors(&(s->cur), s->stack, B);
}



/******************************************************************************/


static int bls_theorem5_limit(mpz_t n, mpz_t A, mpz_t B, UV dummy,
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
static int bls_theorem7_limit(mpz_t n, mpz_t F1, mpz_t R1, UV B1,
                              mpz_t t, mpz_t y, mpz_t r, mpz_t s)
{
  mpz_mul(t, F1, R1);
  mpz_add_ui(t, t, 1);
  if (mpz_cmp(t, n) != 0) croak("BLS75 internal error: F1*R1 != n-1\n");

  mpz_mul_ui(t, F1, 2);
  mpz_tdiv_qr(s, r, R1, t);

  mpz_add(y, t, r);
  mpz_sub_ui(y, y, B1);
  mpz_mul(y, y, F1);
  mpz_add_ui(y, y, 1);   /* y = 2F1^2 + (r - B1)F1 + 1 */
  mpz_mul_ui(t, F1, B1);
  mpz_add_ui(t, t, 1);
  mpz_mul(y, y, t);      /* times (B1F1+1) */

  return (mpz_cmp(n, y) < 0) ? 1 : 0;
}

/*
17:  N < (m F2 - 1)  ( 2 F2 F2 + m F2 - |r| F2 + 1 )
     (III) test (l*F2+1) doesn't divide N for 1 .. m

19:  N < (B2 F2 - 1) ( 2 F2 F2 + B2 F2 - |r| F2 + 1 )
     (III) (IV) R2 factors > B2
*/

#if 0
static int bls_theorem17_limit(mpz_t n, mpz_t F2, mpz_t R2, UV dummy,
                               mpz_t t, mpz_t y, mpz_t r, mpz_t s)
{
  mpz_mul(t, F2, R2);
  mpz_sub_ui(t, t, 1);
  if (mpz_cmp(t, n) != 0) croak("BLS75 internal error: F2*R2 != n+1\n");

  mpz_mul_ui(t, F2, 2);
  mpz_tdiv_qr(s, r, R2, t);
  if (mpz_cmp(r, F2) >= 0) {
    mpz_add_ui(s, s, 1);
    mpz_sub(r, r, t);
  }
  /* Let m = 1 */
  mpz_add_ui(y, t, 1);
  mpz_abs(t, r);
  mpz_sub(y, y, t);
  mpz_mul(y, y, F2);
  mpz_add_ui(y, y, 1);   /* y = 2F2^2 + (m-r)F2 + 1 */
  mpz_sub_ui(t, F2, 1);
  mpz_mul(y, y, t);      /* times (mF2-1) */

  return (mpz_cmp(n, y) < 0) ? 1 : 0;
}
#endif
static int bls_theorem19_limit(mpz_t n, mpz_t F2, mpz_t R2, UV B2,
                               mpz_t t, mpz_t y, mpz_t r, mpz_t s)
{
  mpz_mul(t, F2, R2);
  mpz_sub_ui(t, t, 1);
  if (mpz_cmp(t, n) != 0) croak("BLS75 internal error: F2*R2 != n+1\n");

  mpz_mul_ui(t, F2, 2);
  mpz_tdiv_qr(s, r, R2, t);
  if (mpz_cmp(r, F2) >= 0) {
    mpz_add_ui(s, s, 1);
    mpz_sub(r, r, t);
  }

  mpz_add_ui(y, t, B2);
  mpz_abs(t, r);
  mpz_sub(y, y, t);
  mpz_mul(y, y, F2);
  mpz_add_ui(y, y, 1);   /* y = 2F2^2 + (B2 - r)F2 + 1 */
  mpz_mul_ui(t, F2, B2);
  mpz_sub_ui(t, t, 1);
  mpz_mul(y, y, t);      /* times (B2F2-1) */

  return (mpz_cmp(n, y) < 0) ? 1 : 0;
}

static int bls_corollary11_limit(mpz_t n, mpz_t R1, mpz_t F1, mpz_t F2, UV B,
                                 mpz_t t, mpz_t g, mpz_t r, mpz_t s)
{
  if (mpz_cmp(F1,F2) >= 0) {
    mpz_tdiv_q_2exp(t, F2, 1);
    mpz_mul(t, t, F1);
    mpz_mul(t, t, F1);
  } else {
    mpz_tdiv_q_2exp(t, F1, 1);
    mpz_mul(t, t, F2);
    mpz_mul(t, t, F2);
  }
  mpz_ui_pow_ui(g, B, 3);
  mpz_mul(t, t, g);
  return (mpz_cmp(n, t) < 0);
}

static int bls_theorem20_limit(mpz_t n, mpz_t R1, mpz_t F1, mpz_t F2,
                               UV B, UV m,
                               mpz_t t, mpz_t g, mpz_t r, mpz_t s)
{
  int m_used = 0;

  if (bls_corollary11_limit(n,R1,F1,F2,B,t,g,r,s)) {
    mpz_set_ui(s,0);   /* No test for m needed */
    return 1;
  }

  mpz_mul_ui(t, F1, B);
  mpz_add_ui(g, t, 1);

  mpz_mul_ui(t, F2, B);
  mpz_sub_ui(t, t, 1);

  if (mpz_cmp(t, g) > 0)  mpz_set(g, t);

  mpz_tdiv_q_2exp(t, F2, 1);
  mpz_tdiv_qr(s, r, R1, t);
  mpz_mul(t, F1, F2);
  mpz_mul_ui(t, t, m);
  mpz_tdiv_q_2exp(t, t, 1);
  mpz_mul(r, r, F1);           /* r *= F1;  t += r */
  mpz_add(t, t, r);
  mpz_add_ui(t, t, 1);         /* t = m * F1 * F2/2 + r * F1 + 1 */

  if (mpz_cmp(t, g) > 0) {
    m_used = 1;
    mpz_set(g, t);
  }

  mpz_mul(t, F1, F2);
  mpz_tdiv_q_2exp(t, t, 1);
  mpz_mul_ui(t, t, B);
  mpz_mul_ui(t, t, B);
  mpz_add_ui(s, t, 1);    /* s = B1*B2*F1*F2/2+1 */
  mpz_mul(g, g, s);

  mpz_set_ui(s, m_used);   /* Use s to signal whether we must test for m. */
  return (mpz_cmp(n, g) < 0);
}


/******************************************************************************/


/* (I) For each prime p_i dividing F1 [N-1 = F1R1] there exists an a_i
 *     such that N is a psp base a_i and gcd(a_i^{(N-1)/p_i}-1,N) = 1.
 */
static int _verify_cond_I_p(mpz_t n, mpz_t pi, mpz_t ap, mpz_t t, int alimit, char* pspcache)
{
  int a, success = 0;
  PRIME_ITERATOR(iter);

  for (a = 2; !success && a <= alimit; a = prime_iterator_next(&iter)) {
     int psp = pspcache  ?  pspcache[a]  :  -1;
     mpz_set_ui(ap, a);

     if (psp == -1) {
       mpz_sub_ui(t, n, 1);
       mpz_powm(t, ap, t, n);
       psp = (mpz_cmp_ui(t, 1) == 0);
     }
     if (!psp)
       return -1;  /* We failed a Fermat test, n is composite. */
     if (pspcache)  pspcache[a] = psp;

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

  /* No result found in first <alimit> primes.  Check random values. */
  if (!success) {
    int starta = a;
    for (a = 0; !success && a < 100; a++) {
      if (a < 10) {
        mpz_set_ui(ap, irand64(32));    /* ap: 0 to 2^32-1 */
        mpz_add_ui(ap, ap, starta);     /* ap: starta to 2^32-1+starta */
      } else {
        mpz_sub_ui(t, n, starta+2);
        mpz_isaac_urandomm(ap, t);      /* ap: 0 to n-starta-2 */
        mpz_add_ui(ap, ap, 2);          /* ap: starta to n-2 */
      }

      mpz_sub_ui(t, n, 1);
      mpz_powm(t, ap, t, n);
      if (mpz_cmp_ui(t,1) != 0)
       return -1;  /* We failed a Fermat test, n is composite. */

      mpz_sub_ui(t, n, 1);
      mpz_divexact(t, t, pi);
      mpz_powm(t, ap, t, n);
      mpz_sub_ui(t, t, 1);
      mpz_gcd(t, t, n);
      if (mpz_cmp_ui(t, 1) != 0)
        continue;
      success = 1;  /* we found an a for this p */
    }
  }
  if (ev>1) gmp_printf("  %s pi %Zd  a=%Zd\n", success ? "PASS" : "FAIL", pi, ap);
  return success;
}

/* (III) For each prime q_i dividing F2 [N+1 = F2R2] there exists a Lucas
 *       sequence U_k with discriminant D for which D/N = -1, N divides
 *       U_{N+1} and gcd(U_{(N+1)/q_i},N} = 1.
 *       [Note: All sequences used must have the same D.]
 */

static int _test_III_D(mpz_t n, mpz_t np1, IV D, IV P, IV Q, fstack_t* qi, mpz_t R2, mpz_t t, mpz_t U, mpz_t t1, mpz_t t2)
{
  mpz_t t3;
  char *qdone;
  int i, ipq, nqi, nqileft;

  nqi = nstack(qi)+1;
  nqileft = nqi;
  Newz(0, qdone, nqi, char);

  if (mpz_cmp_ui(R2,1) == 0) { qdone[nqi-1] = 1; nqileft--; }

  for (ipq = 0; nqileft > 0 && ipq < 12; ipq++) {
    if (ipq > 0) {  Q += P+1;  P += 2;  }
    if (ev > 1) gmp_printf("   n %Zd   III start  D %ld  P %ld Q %ld\n", n, D, P, Q);
    if ((P*P-4*Q) != D) croak("test_III_D bad D construction");
    mpz_set_iv(t1,P);
    mpz_set_iv(t2,Q);
    mpz_gcd(t, n, t2);
    if (mpz_cmp_ui(t,1) > 0 && mpz_cmp(t,n) < 0) {  /* Found a factor of n */
      Safefree(qdone);
      return -1;
    }
    if (mpz_cmp_ui(t,1) != 0)
      continue;
    lucasumod(U, t1, t2, np1, n, t);
    if (mpz_sgn(U) != 0)
      continue;

    if (ev > 1) gmp_printf("   n %Zd   III facts  D %ld  P %ld Q %ld\n", n, D, P, Q);
    /* This P/Q sequence with the given D is acceptable.  Now check each Qi. */
    mpz_init(t3);
    for (i = 0; i < nqi; i++) {
      if (qdone[i])
        continue;
      mpz_divexact(t3, np1, (i == nqi-1) ? R2 : qi->stack[i]);
      lucasumod(U, t1, t2, t3, n, t);
      mpz_gcd(t, n, U);
      if (mpz_cmp_ui(t,1) > 0 && mpz_cmp(t,n) < 0) {  /* Found a factor of n */
        mpz_clear(t3);
        Safefree(qdone);
        return -1;   /* We found a factor of n */
      }
      if (mpz_cmp_ui(t,1) == 0) {
        qdone[i] = 1;
        nqileft--;
        continue;
      }
    }
    mpz_clear(t3);
  }
  Safefree(qdone);
  return (nqileft == 0)  ?  2 /* Prime */  :  1 /* Not sure */;
}

static int _test_III_IV(mpz_t n, mpz_t np1, mpz_t R2, fstack_t* fstack,
                        mpz_t t, mpz_t m, mpz_t r, mpz_t s)
{
  PRIME_ITERATOR(iter);
  IV P, Q, D;
  int i, knd, success = 1;

  for (Q = 2; success == 1 && Q < 72; Q = prime_iterator_next(&iter)) {
    P = 83;
    D = P*P - 4*Q;
    /* if (D == 0) continue; */   /* None with this method. */
    knd = mpz_si_kronecker(D,n);
    if (knd == 0) {
      mpz_set_iv(t,D);
      mpz_gcd(t,t,n);
      if (mpz_cmp_ui(t,1) > 0 && mpz_cmp(t,n) < 0)
        success = -1;  /* We found a factor */
    }
    if (knd != -1) continue;
    /* -1=comp,2=prime,1=dunno */
    success = _test_III_D(n, np1, D, P, Q, fstack, R2, /*temps*/ t, m, r, s);
  }
  for (i = 0; success == 1 && i < 12; i++) {
    P = 100 + irand64((maxuibits/2)-4);
    Q = 100 + irand64((maxuibits/2)-4);
    D = P*P - 4*Q;
    if (D == 0) continue;
    knd = mpz_si_kronecker(D,n);
    if (knd == 0) {
      mpz_set_iv(t,D);
      mpz_gcd(t,t,n);
      if (mpz_cmp_ui(t,1) > 0 && mpz_cmp(t,n) < 0)
        success = -1;  /* We found a factor */
    }
    if (knd != -1) continue;
    /* -1=comp,2=prime,1=dunno */
    success = _test_III_D(n, np1, D, P, Q, fstack, R2, /*temps*/ t, m, r, s);
  }
  prime_iterator_destroy(&iter);
  if (success == 1)
    success = 0;
  return success;  /* -1 = composite, 2 = prime, 0 = no proof found */
}

  /*
   *   -1   definitely composite.
   *    0   no proof.  We can't say anything about n.
   *    2   definitely prime.
   */
static int _prove_T19(mpz_t n, mpz_t np1, mpz_t F2, mpz_t R2, UV B2,
                      fstack_t* fstack, mpz_t t, mpz_t m, mpz_t r, mpz_t s)
{
  trim_factors(F2, R2, n, np1, B2, fstack, &bls_theorem19_limit, t, m, r, s);
  if (!bls_theorem19_limit(n, F2, R2, B2, t,m,r,s))
    return 0;
   mpz_mul(t, r, r);
   mpz_addmul_ui(t, s, 8);   /* t = r^2 + 8s */
   /* N is prime if and only if s=0 OR t not a perfect square */
   if (mpz_sgn(s) != 0 && mpz_perfect_square_p(t))
     return -1;
  if (ev) gmp_printf("N %Zd  F2 %Zd  R2 %Zd  B2 %lu\n", n, F2, R2, B2);

  /* Now verify (III) and (IV). */
  return _test_III_IV(n, np1, R2, fstack, t, m, r, s);
}


/******************************************************************************/



int BLS_primality_nm1(mpz_t n, int effort, char** prooftextptr)
{
  mpz_t nm1, F1, R1, t, m, f, r, s;
  FACTOR_STACK(fstack);
  FACTOR_STACK(mstack);
  int e, success = 1;
  UV B1 = (effort < 2 && mpz_sizeinbase(n,2) < 160) ?  6000 :
          (mpz_sizeinbase(n,2) < 1000)              ? 20000 : 200000;
  limit_func_t limitfunc = (prooftextptr) ? bls_theorem5_limit : bls_theorem7_limit;

  /* We need to do this for BLS */
  if (mpz_even_p(n)) return (mpz_cmp_ui(n,2) == 0);

  mpz_init(nm1);
  mpz_sub_ui(nm1, n, 1);
  mpz_init_set_ui(F1, 1);
  mpz_init_set(R1, nm1);
  mpz_init(m);
  mpz_init(f);
  mpz_init(t);
  mpz_init(r);
  mpz_init(s);

  { /* Pull small factors out */
    PRIME_ITERATOR(iter);
    UV tf;
    for (tf = 2; tf < B1; tf = prime_iterator_next(&iter)) {
      if (mpz_cmp_ui(R1, tf*tf) < 0)
        break;
      if (factor_test_ui(tf, R1, F1, &fstack)) {
        if (limitfunc(n, F1, R1, B1, t, m, r, s))
          break;
      }
    }
    B1 = (tf >= B1) ?  B1  :  (tf == 2)  ?  3  :  tf+2;
    prime_iterator_destroy(&iter);
  }

  if (prooftextptr) B1 = 3;  /* BLS Theorem 5 doesn't care about B1 */

  if (ev) gmp_printf("n-1 after trial  %d n %Zd  nm1 %Zd  F1 %Zd R1 %Zd B1 %lu\n",success,n,nm1,F1,R1,B1);

  /* If we got enough from trial factoring then no more factoring needed */
  if (limitfunc(n, F1, R1, B1, t, m, r, s))
    goto start_nm1_proof;

  if (success && mpz_cmp_ui(R1,1) > 0) {
    mpz_set(f, R1);
    handle_factor(f, R1, F1, &fstack, &mstack, effort, prooftextptr, 1, &BLS_primality_nm1);
  }

  if (ev) gmp_printf("n-1 after R1     %d n %Zd  nm1 %Zd  F1 %Zd R1 %Zd B1 %lu\n",success,n,nm1,F1,R1,B1);

  while (success) {

    if (limitfunc(n, F1, R1, B1, t, m, r, s))
      break;

    success = 0;
    if (nstack(&mstack) == 0) /* If the stack is empty, we have failed. */
      break;
    pop_fstack(m, &mstack);   /* pop a component off the stack */

    for (e = 0; !success && e <= effort; e++)
      success = tfe(f, m, e);

    if (ev) gmp_printf("n-1 factored    %d n %Zd  nm1 %Zd  F1 %Zd R1 %Zd B1 %lu\n",success,n,nm1,F1,R1,B1);
    /* If we couldn't factor m and the stack is empty, we've failed. */
    if ( (!success) && (nstack(&mstack) == 0) )
      break;
    /* Put the two factors f and m/f into the stacks, smallest first */
    mpz_divexact(m, m, f);
    if (mpz_cmp(m, f) < 0)
      mpz_swap(m, f);
    handle_factor(f, R1, F1, &fstack, &mstack, effort, prooftextptr, 0, &BLS_primality_nm1);
    handle_factor(m, R1, F1, &fstack, &mstack, effort, prooftextptr, 0, &BLS_primality_nm1);
  }

start_nm1_proof:

  if (ev) gmp_printf("n-1 start proof  %d n %Zd  nm1 %Zd  F1 %Zd R1 %Zd B1 %lu\n",success,n,nm1,F1,R1,B1);

  /* clear mstack since we don't care about it.  Use to hold a values. */
  clear_fstack(&mstack);

  fstack_tidy(&fstack, B1);

  /* Shrink to smallest set and verify conditions. */
  if (success > 0) {
    trim_factors(F1, R1, n, nm1, B1, &fstack, limitfunc, t, m, r, s);
    /* Verify conditions */
    success = 0;
    if (limitfunc(n, F1, R1, B1, t, m, r, s)) {
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

    for (pcount = 0; success > 0 && pcount < fstack.cur; pcount++) {
      success = _verify_cond_I_p(n, fstack.stack[pcount], ap, t, alimit, afermat);
      if (success > 0)
        push_fstack(&mstack, ap);
    }
    /* Using T7 proof, need to test R1 */
    if (success > 0 && !prooftextptr && mpz_cmp_ui(R1,1) > 0) {
      success = _verify_cond_I_p(n, R1, ap, t, alimit, afermat);
      if (success > 0)
        push_fstack(&mstack, ap);
    }

    /* If we could not find 'a' values, then we should return 1 (maybe prime)
     * since we did not perform an exhaustive search.  It would be quite
     * unusual to find a prime that didn't have an 'a' in the first 10,000
     * primes, but it could happen.  It's a "dubiously prime" :) */
    if (success == 0 && get_verbose_level() > 0)
      printf("N-1 factored but failed to prove.  Perhaps composite.\n");
    mpz_clear(p);
    mpz_clear(ap);
  }

  if (success > 0 && prooftextptr != 0) {
    int i;
    char *proofstr, *proofptr;
    int curprooflen = (*prooftextptr == 0) ? 0 : strlen(*prooftextptr);
    int fsp = nstack(&fstack);
    int msp = nstack(&mstack);
    int myprooflen = (5 + mpz_sizeinbase(n, 10)) * (2 + fsp + msp) + 200;

    if (fsp != msp) croak("Different f and a counts\n");
    New(0, proofstr, myprooflen + curprooflen + 1, char);
    proofptr = proofstr;
    proofptr += gmp_sprintf(proofptr, "Type BLS5\nN  %Zd\n", n);
    /* Q[0] is always 2 */
    for (i = 1; i < fsp; i++)
      proofptr += gmp_sprintf(proofptr, "Q[%d]  %Zd\n", i, fstack.stack[i]);
    /* A[i] only printed if not 2 */
    for (i = 0; i < msp; i++)
      if (mpz_cmp_ui(mstack.stack[i], 2) != 0)
        proofptr += gmp_sprintf(proofptr, "A[%d]  %Zd\n", i, mstack.stack[i]);
    proofptr += gmp_sprintf(proofptr, "----\n");
    /* Set or prepend */
    if (*prooftextptr) {
      proofptr += gmp_sprintf(proofptr, "\n");
      strcat(proofptr, *prooftextptr);
      Safefree(*prooftextptr);
    }
    *prooftextptr = proofstr;
  }

  destroy_fstack(&fstack);
  destroy_fstack(&mstack);
  mpz_clear(nm1);
  mpz_clear(F1);
  mpz_clear(R1);
  mpz_clear(m);
  mpz_clear(f);
  mpz_clear(t);
  mpz_clear(r);
  mpz_clear(s);
  if (success < 0) return 0;
  if (success > 0) return 2;
  return 1;
}


int BLS_primality_np1(mpz_t n, int effort, char** prooftextptr)
{
  mpz_t np1, F2, R2, t, m, f, r, s;
  FACTOR_STACK(fstack);
  FACTOR_STACK(mstack);
  int e, success = 1;
  UV B2 = (effort < 2 && mpz_sizeinbase(n,2) < 160) ?  6000 :
          (mpz_sizeinbase(n,2) < 1000)              ? 20000 : 200000;
  limit_func_t limitfunc = bls_theorem19_limit;

  /* We need to do this for BLS */
  if (mpz_even_p(n)) return (mpz_cmp_ui(n,2) == 0) ? 2 : 0;

  mpz_init(np1);
  mpz_add_ui(np1, n, 1);
  mpz_init_set_ui(F2, 1);
  mpz_init_set(R2, np1);
  mpz_init(m);
  mpz_init(f);
  mpz_init(t);
  mpz_init(r);
  mpz_init(s);

  { /* Pull small factors out */
    PRIME_ITERATOR(iter);
    UV tf;
    for (tf = 2; tf < B2; tf = prime_iterator_next(&iter)) {
      if (mpz_cmp_ui(R2, tf*tf) < 0)
        break;
      if (factor_test_ui(tf, R2, F2, &fstack)) {
        if (limitfunc(n, F2, R2, B2, t, m, r, s))
          break;
      }
    }
    B2 = (tf >= B2) ?  B2  :  (tf == 2)  ?  3  :  tf+2;
    prime_iterator_destroy(&iter);
  }
  if (ev) gmp_printf("N %Zd  F2 %Zd  R2 %Zd  B2 %lu\n", n, F2, R2, B2);

  /* printf("trial  np1: %lu bits, %lu bits factored\n", mpz_sizeinbase(np1,2), mpz_sizeinbase(F2,2)); */

  /* If we got enough from trial factoring then no more factoring needed */
  if (limitfunc(n, F2, R2, B2, t, m, r, s))
    goto start_np1_proof;

  if (success && mpz_cmp_ui(R2,1) > 0) {
    mpz_set(f, R2);
    handle_factor(f, R2, F2, &fstack, &mstack, effort, prooftextptr, 1, &BLS_primality_np1);
  }

  while (success) {

    if (limitfunc(n, F2, R2, B2, t, m, r, s))
      break;

    success = 0;

    if (nstack(&mstack) == 0) /* If the stack is empty, we have failed. */
      break;
    pop_fstack(m, &mstack);   /* pop a component off the stack */

    for (e = 0; !success && e <= effort; e++)
      success = tfe(f, m, e);

    /* If we couldn't factor m and the stack is empty, we've failed. */
    if ( (!success) && (nstack(&mstack) == 0) )
      break;
    /* Put the two factors f and m/f into the stacks, smallest first */
    mpz_divexact(m, m, f);
    if (mpz_cmp(m, f) < 0)
      mpz_swap(m, f);
    handle_factor(f, R2, F2, &fstack, &mstack, effort, prooftextptr, 0, &BLS_primality_np1);
    handle_factor(m, R2, F2, &fstack, &mstack, effort, prooftextptr, 0, &BLS_primality_np1);
  }
  /* TODO: success = 0 here just means we couldn't factor */

start_np1_proof:

  /* clear mstack since we don't care about it.  Use to hold a values. */
  clear_fstack(&mstack);

  /* printf("factor  np1: %lu bits, %lu bits factored\n", mpz_sizeinbase(np1,2), mpz_sizeinbase(F2,2)); */

  fstack_tidy(&fstack, B2);
  if (ev) gmp_printf("end trim:  N %Zd  F2 %Zd  R2 %Zd  B2 %lu\n", n, F2, R2, B2);

  if (success > 0)
    success = _prove_T19(n, np1, F2, R2, B2, &fstack, t, m, r, s);
  /* -1 = composite, 2 = prime, 0 = no proof found */

  /* TODO: Proof text */

  destroy_fstack(&fstack);
  destroy_fstack(&mstack);
  mpz_clear(np1);
  mpz_clear(F2);
  mpz_clear(R2);
  mpz_clear(m);
  mpz_clear(f);
  mpz_clear(t);
  mpz_clear(r);
  mpz_clear(s);
  if (success < 0) return 0;   /* Composite */
  if (success > 0) return 2;   /* Prime */
  return 1;                    /* Not sure */
}

/******************************************************************************/

#define PRINT_PCT 0

/* Will use one of:
 *    N-1   Corollary 1
 *    N-1   Theorem 5
 *    N-1   Theorem 7
 *    N+1   Corollary 8      not used now
 *    N+1   Theorem 17       not used now
 *    N+1   Theorem 19       not used now
 *    Comb  Theorem 20
 *    Comb  Corollary 11
 */
int BLS_primality(mpz_t n, int effort, char** prooftextptr)
{
  mpz_t nm1, np1, F1, F2, R1, R2;
  mpz_t r, s, t, u, f, c1, c2;
  /* fstack:  definite prime factors
   * pstack:  probable prime factors   product of fstack and pstack = F
   * mstack:  composite remainders     product of mstack = R
   */
  FACTOR_STACK(f1stack);
  FACTOR_STACK(f2stack);
  FACTOR_STACK(p1stack);
  FACTOR_STACK(p2stack);
  FACTOR_STACK(m1stack);
  FACTOR_STACK(m2stack);
  int pcount, e, success = 1;
  int low_effort = (effort < 1) ? 0 : 1;
  UV m, B1 = (effort < 2 && mpz_sizeinbase(n,2) < 160) ?  6000 :
             (mpz_sizeinbase(n,2) < 1024)              ? 20000 : 200000;
#if PRINT_PCT
  double trial_pct, prime_pct, fac_pct, fin_pct;
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
  mpz_init(f); mpz_init(c1); mpz_init(c2);

  { /* Pull small factors out */
    PRIME_ITERATOR(iter);
    UV tf;
    for (tf = 2; tf < B1; tf = prime_iterator_next(&iter)) {
      /* Page 635 of BLS75 describes an optimization for divisibility
       * testing.  It seems slower than just doing two UI div tests. */
      int testr1 = (mpz_cmp_ui(R1, tf*tf) >= 0);
      int testr2 = (mpz_cmp_ui(R2, tf*tf) >= 0);
      if (!testr1 || !testr2) { B1 = tf; break; }
      if (ev) printf("   factoring out %lu\n", tf);
      if (testr1) factor_test_ui(tf, R1, F1, &f1stack);
      if (testr2) factor_test_ui(tf, R2, F2, &f2stack);
    }
    prime_iterator_destroy(&iter);
  }
  m = B1-1;   /* m should be less than B1 */
  if (ev) gmp_printf("N %Zd  F1 %Zd  R1 %Zd  B1 %lu\n", n, F1, R1, B1);
  if (ev) gmp_printf("N %Zd  F2 %Zd  R2 %Zd  B1 %lu\n", n, F2, R2, B1);

#if PRINT_PCT
  trial_pct = (100.0 * (mpz_sizeinbase(F1,2) + mpz_sizeinbase(F2,2))) / (mpz_sizeinbase(nm1,2) + mpz_sizeinbase(np1,2));
  printf("\n%6.2f .. ", trial_pct);  fflush(stdout);
#endif

  if ( bls_theorem7_limit(n, F1, R1, B1, t, u, r, s) ||
       bls_theorem20_limit(n, R1, F1, F2, B1, m, t, u, r, s) )
    goto start_hybrid_proof;

  if (mpz_cmp_ui(R1,1) > 0) {
    mpz_set(f, R1);
    handle_factor2(f, R1, F1, &f1stack, &p1stack, &m1stack, low_effort, prooftextptr, &BLS_primality);
  }
  if (mpz_cmp_ui(R2,1) > 0) {
    mpz_set(f, R2);
    handle_factor2(f, R2, F2, &f2stack, &p2stack, &m2stack, low_effort, prooftextptr, &BLS_primality);
  }

#if PRINT_PCT
  prime_pct = (100.0 * (mpz_sizeinbase(F1,2) + mpz_sizeinbase(F2,2))) / (mpz_sizeinbase(nm1,2) + mpz_sizeinbase(np1,2));
  printf("\n%6.2f .. ", prime_pct);  fflush(stdout);
#endif

  while (1) {
    int d1, d2;

    success = 1;
    if ( bls_theorem7_limit(n, F1, R1, B1, t, u, r, s) ||
         bls_theorem20_limit(n, R1, F1, F2, B1, m, t, u, r, s) )
      break;

    success = 0;
    mpz_set_ui(c1, 0);
    mpz_set_ui(c2, 0);
    d1 = nstack(&m1stack) > 0;
    d2 = nstack(&m2stack) > 0;
    if (!d1 && !d2) break;

    if (d1) pop_fstack(c1, &m1stack);
    if (d2) pop_fstack(c2, &m2stack);
    for (e = 0; !success && e <= effort; e++) {
      if (d1 && tfe(f, c1, e)) {
        if (d2) push_fstack(&m2stack, c2);
        mpz_set(u, c1);
        success = 1;
      } else if (d2 && tfe(f, c2, e)) {
        if (d1) push_fstack(&m1stack, c1);
        mpz_set(u, c2);
        success = 2;
      }
    }

    /* No success for this set of composites.  Move on. */
    if (!success) continue;

    if (success == 1) {
      mpz_divexact(u, u, f);
      if (mpz_cmp(u, f) < 0)
        mpz_swap(u, f);
      handle_factor2(f, R1, F1, &f1stack, &p1stack, &m1stack, low_effort, prooftextptr, &BLS_primality);
      handle_factor2(u, R1, F1, &f1stack, &p1stack, &m1stack, low_effort, prooftextptr, &BLS_primality);
    } else if (success == 2) {
      mpz_divexact(u, u, f);
      if (mpz_cmp(u, f) < 0)
        mpz_swap(u, f);
      handle_factor2(f, R2, F2, &f2stack, &p2stack, &m2stack, low_effort, prooftextptr, &BLS_primality);
      handle_factor2(u, R2, F2, &f2stack, &p2stack, &m2stack, low_effort, prooftextptr, &BLS_primality);
    }
#if PRINT_PCT
    fac_pct = (100.0 * (mpz_sizeinbase(F1,2) + mpz_sizeinbase(F2,2))) / (mpz_sizeinbase(nm1,2) + mpz_sizeinbase(np1,2));
    printf("%6.2f .. ", fac_pct);  fflush(stdout);
#endif
  }

start_hybrid_proof:
  mpz_mul(t, F1, R1); if (mpz_cmp(nm1, t) != 0) croak("Bad n-1 factor");
  mpz_mul(t, F2, R2); if (mpz_cmp(np1, t) != 0) croak("Bad n+1 factor");

  /* We've done all the factoring we need or can. */
  if (ev) gmp_printf("start hybrid:  N %Zd  F2 %Zd  R2 %Zd  B1 %lu\n", n, F2, R2, B1);

  /* Finish proofs for p{1,2}stack as needed. */
  /* TODO: optimize for cases of both n-1 and n+1 working */
  if (nstack(&p1stack) > 0) {
    while (nstack(&p1stack) > 0) {
      int pr = 1;
      pop_fstack(f, &p1stack);
      if (effort > low_effort)
        pr = BLS_primality(f, effort, prooftextptr);
      if      (pr == 0) croak("probable prime factor proved composite");
      else if (pr == 2) push_fstack(&f1stack, f); /* Proved, put on F stack */
      else              factor_out(F1, R1, f);    /* No proof.  Move to R */
    }
  }
  if (nstack(&p2stack) > 0) {
    while (nstack(&p2stack) > 0) {
      int pr = 1;
      pop_fstack(f, &p2stack);
      if (effort > low_effort)
        pr = BLS_primality(f, effort, prooftextptr);
      if      (pr == 0) croak("probable prime factor proved composite");
      else if (pr == 2) push_fstack(&f2stack, f); /* Proved, put on F stack */
      else              factor_out(F2, R2, f);    /* No proof.  Move to R */
    }
  }

  if (ev) gmp_printf("start tidy:  N %Zd  F2 %Zd  R2 %Zd  B1 %lu\n", n, F2, R2, B1);
  fstack_tidy(&f1stack, B1);
  fstack_tidy(&f2stack, B1);
  if (ev) gmp_printf("end tidy:  N %Zd  F2 %Zd  R2 %Zd  B1 %lu\n", n, F2, R2, B1);

#if PRINT_PCT
  fin_pct = (100.0 * (mpz_sizeinbase(F1,2) + mpz_sizeinbase(F2,2))) / (mpz_sizeinbase(nm1,2) + mpz_sizeinbase(np1,2));
  printf("%6.2f .. ", fin_pct);  fflush(stdout);
  printf("\n");  fflush(stdout);
#endif

  if (ev) gmp_printf("check:  N %Zd  F2 %Zd  R2 %Zd  B1 %lu\n", n, F2, R2, B1);
  /* Check the theorems we have available */

  /* If we can do a standard n-1 proof, do that. */
  if (bls_theorem7_limit(n, F1, R1, B1, t, u, r, s)) {
    if (get_verbose_level() > 0) printf("BLS75 proof using N-1\n");
    if (ev) gmp_printf("N %Zd  F1 %Zd  R1 %Zd  B1 %lu\n", n, F1, R1, B1);
    trim_factors(F1, R1, n, nm1, B1, &f1stack, &bls_theorem7_limit, t, u, r, s);
    if (ev) gmp_printf("N %Zd  F1 %Zd  R1 %Zd  B1 %lu\n", n, F1, R1, B1);
    for (pcount = 0; success > 0 && pcount < f1stack.cur; pcount++)
      success = _verify_cond_I_p(n, f1stack.stack[pcount], u, t, 1000, 0);
    if (success > 0 && (mpz_mul(t, F1, F1), mpz_cmp(t,n) > 0))
      goto end_hybrid;  /* Corollary 1, n-1 factored more than sqrt(n) */
    if (success > 0 && !bls_theorem5_limit(n, F1, R1, B1, t, u, r, s))
      success = _verify_cond_I_p(n, R1, u, t, 1000, 0);
    if (success > 0) {
      mpz_mul(t, r, r);
      mpz_submul_ui(t, s, 8);   /* t = r^2 - 8s */
      /* N is prime if and only if s=0 OR t not a perfect square */
      success = (mpz_sgn(s) == 0 || !mpz_perfect_square_p(t))  ?  1  :  -1;
    }
    goto end_hybrid;    /* Theorem 5 or 7 */
  }

  /* Rather than looking at theorem 19 for an N+1 proof, we'll just go to
   * theorem 20 (or corollary 11).  T20 is faster. */

  /* Check N < B^3 F1*F2*F2/2  or  N < B^3 F1*F1*F2/2 */
  success = bls_theorem20_limit(n, R1, F1, F2, B1, m, t, u, r, s);

  if (get_verbose_level() > 0) printf("BLS75 proof using N-1 / N+1 (T20)\n");

  /* Trim some factors from f2stack if possible */
  if (nstack(&f2stack) > 1) {
    int i;
    mpz_set_ui(F2, 1);
    mpz_set(R2, np1);
    for (i = 0; i < f2stack.cur; i++) {
      if (i > 0 && bls_theorem20_limit(n, R1, F1, F2, B1, m, t, u, r, s))
        break;
      factor_out(R2, F2, f2stack.stack[i]);
    }
    /* Remove excess factors */
    while (i < f2stack.cur)
      pop_fstack(t, &f2stack);
    /* Verify Q[0] = 2 */
    if (mpz_cmp_ui(f2stack.stack[0], 2) != 0)
      croak("BLS75 internal error: 2 not at start of fstack");
  }

  /* Check lambda divisibility if needed */
  if (success > 0 && mpz_sgn(s)) {
    UV lambda;
    mpz_mul(t, F1, F2);
    mpz_tdiv_q_2exp(t, t, 1);
    mpz_mul(u, r, F1);
    mpz_add_ui(u, u, 1);
    for (lambda = 0; success > 0 && lambda < m; lambda++, mpz_add(u,u,t)) {
      if (lambda > 0 || mpz_sgn(r))
        if (mpz_divisible_p(n, u)) {
          /* Check that we found a non-trivial divisor */
          mpz_gcd(t, u, n);
          success = (mpz_cmp_ui(t,1) > 0 && mpz_cmp(t,n) < 0) ? -1 : 0;
          break;
        }
    }
  }

  /* Verify (I)   page 623 and (II) page 625 */
  if (success > 0) {
    for (pcount = 0; success > 0 && pcount < f1stack.cur; pcount++)
      success = _verify_cond_I_p(n, f1stack.stack[pcount], u, t, 1000, 0);
    if (success > 0)
      success = _verify_cond_I_p(n, R1, u, t, 1000, 0);
  }

  /* Verify (III) page 631 and (IV) page 633 */
  if (success > 0) {
    success = _test_III_IV(n, np1, R2, &f2stack, t, u, r, s);
  }

#if 0
  { double p1 = (100.0 * mpz_sizeinbase(F1,2) / mpz_sizeinbase(nm1,2));
    double p2 = (100.0 * mpz_sizeinbase(F2,2) / mpz_sizeinbase(np1,2));
    printf("%6.2f  %6.2f\n", p1, p2);  fflush(stdout); }
  //{ double pct = (100.0 * (mpz_sizeinbase(R1,2) + mpz_sizeinbase(R2,2))) / (mpz_sizeinbase(nm1,2) + mpz_sizeinbase(np1,2)); printf("%6.2f\n", 100.0-pct);  fflush(stdout); }
#endif

end_hybrid:
  destroy_fstack(&f1stack);
  destroy_fstack(&f2stack);
  destroy_fstack(&p1stack);
  destroy_fstack(&p2stack);
  destroy_fstack(&m1stack);
  destroy_fstack(&m2stack);
  mpz_clear(nm1); mpz_clear(np1);
  mpz_clear(F1);  mpz_clear(F2);
  mpz_clear(R1);  mpz_clear(R2);
  mpz_clear(r);
  mpz_clear(s);
  mpz_clear(u);
  mpz_clear(t);
  mpz_clear(f); mpz_clear(c1); mpz_clear(c2);
  if (success < 0) return 0;
  if (success > 0) return 2;
  return 1;
}



/******************************************************************************/

/* Given an n where we're factored n-1 down to p, check BLS theorem 3 */
int BLS_check_T3(mpz_t n, mpz_t p, UV* reta)
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
      /* should check kronecker(a,n) == -1 here */
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
int BLS_check_T15(mpz_t n, mpz_t f, IV* lp, IV* lq)
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
      mpz_set_iv(t, d);
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
