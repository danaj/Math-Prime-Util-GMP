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


static int tfe(mpz_t f, mpz_t n, int effort)
{
  int success = 0;
  UV log2n = mpz_sizeinbase(n, 2);

  if (!success && mpz_cmp_ui(n, (unsigned long)(UV_MAX>>4)) < 0) {
    UV ui_n = mpz_get_ui(n);
    UV ui_factors[2];
    if (!mpz_cmp_ui(n, ui_n)) {
      success = racing_squfof_factor(ui_n, ui_factors, 200000)-1;
      if (success)
        mpz_set_ui(f, ui_factors[0]);
    }
  }

  if (!success && effort == 0)  success = (int)power_factor(n, f);

  if (success) return success;

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
    case 8: if (log2n > 240) { success = _GMP_ECM_FACTOR(n, f, 40000, 10); }
            break;
    case 9: if (log2n > 240) { success = _GMP_ECM_FACTOR(n, f, 80000,  5); }
            break;
    case 10:if (log2n > 270) { success = _GMP_ECM_FACTOR(n, f,160000, 20); }
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

static void factor_test_ui(unsigned long f, mpz_t R, mpz_t F, fstack_t* s) {
  if (mpz_divisible_ui_p(R, f)) {
    push_fstack_ui(s, f);
    factor_out_ui(R, F, f);
  }
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

static void trim_factors(mpz_t F, mpz_t R, mpz_t n, mpz_t none, UV flim, fstack_t* fs, limit_func_t func, mpz_t t, mpz_t m, mpz_t r, mpz_t s) {
  if (fs->cur > 1) {
    int i;
    mpz_set_ui(F, 1);
    mpz_set(R, none);
    for (i = 0; i < fs->cur; i++) {
      if (i > 0 && func(n, F, R, flim, t, m, r, s))
        break;
      factor_out(R, F, fs->stack[i]);
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
static void fstack_sort_trim(fstack_t* s) {
  sort_and_trim_factors(&(s->cur), s->stack);
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

/*
17:  N < (m F2 - 1)  ( 2 F2 F2 + m F2 - |r| F2 + 1 )
     (III) test (l*F2+1) doesn't divide N for 1 .. m

19:  N < (B2 F2 - 1) ( 2 F2 F2 + B2 F2 - |r| F2 + 1 )
     (III) (IV) R2 factors > B2
*/

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
static int _verify_cond_III_q2(mpz_t n, mpz_t qi, IV p, IV q, mpz_t U, mpz_t V, mpz_t k, mpz_t t1, mpz_t t2)
{
  mpz_add_ui(k, n, 1);
  mpz_divexact(k, k, qi);
  lucas_seq(U, V, n, p, q, k,  t1, t2);
  mpz_gcd(k, U, n);
  return (mpz_cmp_ui(k, 1) == 0) ? 1 : 0;
}

#define MAXQV 50

static int _verify_cond_III_q(mpz_t n, mpz_t qi, IV* qv, int* pnumqv, IV* lp, IV* lq)
{
  int i, numqv = *pnumqv, rval = 0;
  IV d, p, q, startq;
  mpz_t U, V, k, t1, t2;

  mpz_init(U);  mpz_init(V);  mpz_init(k);  mpz_init(t1);  mpz_init(t2);

  /* Try previous q values with (D|n)=-1 and U_(n+1) = 0 mod n */
  for (i = 0; i < numqv; i++) {
    q = qv[i];
    p = (q % 2) ? 2 : 1;
    if (_verify_cond_III_q2(n, qi, p, q, U, V, k, t1, t2)) {
      rval = 2;
      break;
    }
  }

  if (rval == 0) {
    /* Search for a q value */
    startq = (numqv > 0) ? qv[numqv-1]+1 : 2;
    for (q = startq; q < startq+1000; q++) {
      if (mpz_cmp_ui(n, (unsigned long)q) <= 0) break;
      p = (q % 2) ? 2 : 1;
      d = p*p - 4*q;
      mpz_set_si(t1, d);
      if (mpz_jacobi(t1, n) != -1)
        continue;
      /* we have a d/p/q where d = -1.  Check the first Lucas sequence. */
      mpz_add_ui(k, n, 1);
      lucas_seq(U, V, n, p, q, k,  t1, t2);
      if (mpz_sgn(U) != 0)
        continue;
      /* Passed first test, add to qv list */
      if (numqv < MAXQV) {
        qv[numqv] = q;
        *pnumqv = ++numqv;
      }
      /* Verify second Lucas sequence */
      if (_verify_cond_III_q2(n, qi, p, q, U, V, k, t1, t2)) {
        rval = 2;
        break;
      }
    }
  }
  mpz_clear(U);  mpz_clear(V);  mpz_clear(k);  mpz_clear(t1);  mpz_clear(t2);
  if (lp) *lp = p;
  if (lq) *lq = q;
  return rval;
}
#if 0
/* N divides V_{(N+1)/2} and for q > 2 does not divide V{(N+1)/(2q)} */
static int _verify_theorem14_q(mpz_t n, mpz_t qi, IV* lastq, IV* lp, IV* lq)
{
  int rval = 0;
  IV d, p, q;
  mpz_t U, V, k, t1, t2;

  mpz_init(U);  mpz_init(V);  mpz_init(k);  mpz_init(t1);  mpz_init(t2);

  if (lastq && *lastq > 0 && mpz_cmp_ui(qi,2) > 0) {
    q = *lastq;
    p = (q % 2) ? 2 : 1;
    d = p*p - 4*q;
    mpz_set_si(t1, d);
    if (mpz_jacobi(t1, n) == -1) {
      /* q passed the first tst.  Do the second that depends on the factor. */
      mpz_add_ui(k, n, 1);
      mpz_tdiv_q_2exp(k, k, 1);
      mpz_divexact(k, k, qi);
      lucas_seq(U, V, n, p, q, k,  t1, t2);
      if (mpz_sgn(V) != 0) {
        rval = 2;
      }
    }
  }

  if (!rval) {
    for (q = (lastq && *lastq > 0) ? *lastq + 1 : 2; q < 1000; q++) {
      p = (q % 2) ? 2 : 1;
      d = p*p - 4*q;
      mpz_set_si(t1, d);
      if (mpz_jacobi(t1, n) != -1)
        continue;

      /* we have a d/p/q where d = -1.  Check the Lucas sequence. */

      mpz_add_ui(k, n, 1);
      mpz_tdiv_q_2exp(k, k, 1);
      lucas_seq(U, V, n, p, q, k,    t1, t2);
      if (mpz_sgn(V) == 0) {
        if (mpz_cmp_ui(qi, 2) <= 0) {
          rval = 2; break;
        } else {
          mpz_divexact(k, k, qi);
          lucas_seq(U, V, n, p, q, k,  t1, t2);
          if (mpz_sgn(V) != 0) {
            rval = 2; break;
          }
        }
      }
    }
  }

  mpz_clear(U);  mpz_clear(V);  mpz_clear(k);  mpz_clear(t1);  mpz_clear(t2);
  if (lp) *lp = p;
  if (lq) *lq = q;
  if (lastq)  *lastq = q;
  return rval;
}
#endif



/******************************************************************************/



int _GMP_primality_bls_nm1(mpz_t n, int effort, char** prooftextptr)
{
  mpz_t nm1, A, B, t, m, f, r, s;
  FACTOR_STACK(fstack);
  FACTOR_STACK(mstack);
  int e, success = 1;
  UV B1 = (mpz_sizeinbase(n,10) > 1000) ? 100000 : 10000;

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

  { /* Pull small factors out */
    PRIME_ITERATOR(iter);
    UV tf;
    for (tf = 2; tf < B1; tf = prime_iterator_next(&iter)) {
      if (mpz_cmp_ui(B, tf*tf) < 0) break;
      factor_test_ui(tf, B, A, &fstack);
    }
    prime_iterator_destroy(&iter);
  }

  if (success && mpz_cmp_ui(B,1) > 0) {
    mpz_set(f, B);
    handle_factor(f, B, A, &fstack, &mstack, effort, prooftextptr, 1, &_GMP_primality_bls_nm1);
  }

  while (success) {

    if ( (prooftextptr && bls_theorem5_limit(n, A, B, B1, t, m, r, s)) || (!prooftextptr && bls_theorem7_limit(n, A, B, B1, t, m, r, s)) )
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
    handle_factor(f, B, A, &fstack, &mstack, effort, prooftextptr, 0, &_GMP_primality_bls_nm1);
    handle_factor(m, B, A, &fstack, &mstack, effort, prooftextptr, 0, &_GMP_primality_bls_nm1);
  }

  /* clear mstack since we don't care about it.  Use to hold a values. */
  clear_fstack(&mstack);

  fstack_sort_trim(&fstack);

  /* Shrink to smallest set and verify conditions. */
  if (success > 0) {
    if (prooftextptr)
      trim_factors(A, B, n, nm1, B1, &fstack, &bls_theorem5_limit, t, m, r, s);
    else
      trim_factors(A, B, n, nm1, B1, &fstack, &bls_theorem7_limit, t, m, r, s);
    /* Verify conditions */
    success = 0;
    if ( (prooftextptr && bls_theorem5_limit(n, A, B, B1, t, m, r, s)) || (!prooftextptr && bls_theorem7_limit(n, A, B, B1, t, m, r, s)) ) {
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

    for (pcount = 0; success && pcount < fstack.cur; pcount++) {
      success = _verify_cond_I_p(n, fstack.stack[pcount], ap, t, alimit, afermat);
      if (success)
        push_fstack(&mstack, ap);
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
  mpz_t np1, F2, R2, t, m, f, r, s;
  FACTOR_STACK(fstack);
  FACTOR_STACK(mstack);
  int e, success = 1;
  UV B2 = (mpz_sizeinbase(n,10) > 1000) ? 100000 : 10000;

  /* TODO: T19 doesn't seem right, and we're still
   * passing some composites like 14299. */

  /* We need to do this for BLS */
  if (mpz_even_p(n)) return 0;

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
      if (mpz_cmp_ui(R2, tf*tf) < 0) break;
      factor_test_ui(tf, R2, F2, &fstack);
    }
    prime_iterator_destroy(&iter);
  }

  /* printf("trial  np1: %lu bits, %lu bits factored\n", mpz_sizeinbase(np1,2), mpz_sizeinbase(F2,2)); */

  if (success && mpz_cmp_ui(R2,1) > 0) {
    mpz_set(f, R2);
    handle_factor(f, R2, F2, &fstack, &mstack, effort, prooftextptr, 1, &_GMP_primality_bls_np1);
  }

  while (success) {

    if (bls_theorem17_limit(n, F2, R2, B2, t, m, r, s))
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
    handle_factor(f, R2, F2, &fstack, &mstack, effort, prooftextptr, 0, &_GMP_primality_bls_np1);
    handle_factor(m, R2, F2, &fstack, &mstack, effort, prooftextptr, 0, &_GMP_primality_bls_np1);
  }

  /* clear mstack since we don't care about it.  Use to hold a values. */
  clear_fstack(&mstack);

  /* printf("factor  np1: %lu bits, %lu bits factored\n", mpz_sizeinbase(np1,2), mpz_sizeinbase(F2,2)); */

  fstack_sort_trim(&fstack);

  /* Shrink to smallest set and verify conditions. */
  if (success > 0) {
    trim_factors(F2, R2, n, np1, B2, &fstack, &bls_theorem17_limit, t, m, r, s);
    /* Verify conditions */
    success = 0;
    if (bls_theorem17_limit(n, F2, R2, B2, t, m, r, s)) {
      mpz_mul(t, r, r);
      mpz_addmul_ui(t, s, 8);   /* t = r^2 + 8s */
      /* N is prime if and only if s=0 OR t not a perfect square */
      success = (mpz_sgn(s) == 0 || !mpz_perfect_square_p(t))  ?  1  :  -1;
    }
  }

  if (success > 0) {
    IV qv[MAXQV];
    int pcount, numqv = 0;
    for (pcount = 0; success && pcount < fstack.cur; pcount++) {
      success = _verify_cond_III_q(n, fstack.stack[pcount], qv, &numqv, 0, 0);
    }
    /* If we meet theorem 17 limits, then no need to test R2 */
    if (success && !bls_theorem17_limit(n, F2, R2, B2, t, m, r, s)) {
      success = _verify_cond_III_q(n, R2, qv, &numqv, 0, 0);
    }
    if (!success && get_verbose_level() > 0)
      printf("N+1 factored but failed to prove.  Perhaps composite.\n");
  }

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
  if (success < 0) return 0;
  if (success > 0) return 2;
  return 1;
}

/******************************************************************************/

#define PRINT_PCT 0

/* Will use one of:
 *    N-1   Corollary 1
 *    N-1   Theorem 5
 *    N-1   Theorem 7
 *    N+1   Corollary 8
 *    N+1   Theorem 17
 *    N+1   Theorem 19
 *    Comb  Theorem 20
 */
int bls75_hybrid(mpz_t n, int effort, char** prooftextptr)
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
  UV B1 = (effort < 2 && mpz_sizeinbase(n,2) < 160) ?  6000 :
          (mpz_sizeinbase(n,2) < 1024)              ? 20000 : 200000;
  UV m = B1-1;   /* m should be less than B1 */
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
      if (mpz_cmp_ui(R1, tf*tf) >= 0)
        factor_test_ui(tf, R1, F1, &f1stack);
      if (mpz_cmp_ui(R2, tf*tf) >= 0)
        factor_test_ui(tf, R2, F2, &f2stack);
    }
    prime_iterator_destroy(&iter);
  }

#if PRINT_PCT
  trial_pct = (100.0 * (mpz_sizeinbase(F1,2) + mpz_sizeinbase(F2,2))) / (mpz_sizeinbase(nm1,2) + mpz_sizeinbase(np1,2));
  printf("\n%6.2f .. ", trial_pct);  fflush(stdout);
#endif

  if ( bls_theorem7_limit(n, F1, R1, B1, t, u, r, s) ||
       bls_theorem19_limit(n, F2, R2, B1, t, u, r, s) ||
       bls_theorem20_limit(n, R1, F1, F2, B1, m, t, u, r, s) )
    goto start_hybrid_proof;

  if (mpz_cmp_ui(R1,1) > 0) {
    mpz_set(f, R1);
    handle_factor2(f, R1, F1, &f1stack, &p1stack, &m1stack, low_effort, prooftextptr, &bls75_hybrid);
  }
  if (mpz_cmp_ui(R2,1) > 0) {
    mpz_set(f, R2);
    handle_factor2(f, R2, F2, &f2stack, &p2stack, &m2stack, low_effort, prooftextptr, &bls75_hybrid);
  }

#if PRINT_PCT
  prime_pct = (100.0 * (mpz_sizeinbase(F1,2) + mpz_sizeinbase(F2,2))) / (mpz_sizeinbase(nm1,2) + mpz_sizeinbase(np1,2));
  printf("\n%6.2f .. ", prime_pct);  fflush(stdout);
#endif

  while (1) {
    int d1, d2;

    success = 1;
    if ( bls_theorem7_limit(n, F1, R1, B1, t, u, r, s) ||
         bls_theorem19_limit(n, F2, R2, B1, t, u, r, s) ||
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
      handle_factor2(f, R1, F1, &f1stack, &p1stack, &m1stack, low_effort, prooftextptr, &bls75_hybrid);
      handle_factor2(u, R1, F1, &f1stack, &p1stack, &m1stack, low_effort, prooftextptr, &bls75_hybrid);
    } else if (success == 2) {
      mpz_divexact(u, u, f);
      if (mpz_cmp(u, f) < 0)
        mpz_swap(u, f);
      handle_factor2(f, R2, F2, &f2stack, &p2stack, &m2stack, low_effort, prooftextptr, &bls75_hybrid);
      handle_factor2(u, R2, F2, &f2stack, &p2stack, &m2stack, low_effort, prooftextptr, &bls75_hybrid);
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

  /* Finish proofs for p{1,2}stack as needed. */
  /* TODO: optimize for cases of both n-1 and n+1 working */
  if (nstack(&p1stack) > 0) {
    while (nstack(&p1stack) > 0) {
      int pr = 1;
      pop_fstack(f, &p1stack);
      if (effort > low_effort)
        pr = bls75_hybrid(f, effort, prooftextptr);
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
        pr = bls75_hybrid(f, effort, prooftextptr);
      if      (pr == 0) croak("probable prime factor proved composite");
      else if (pr == 2) push_fstack(&f2stack, f); /* Proved, put on F stack */
      else              factor_out(F2, R2, f);    /* No proof.  Move to R */
    }
  }

  fstack_sort_trim(&f1stack);
  fstack_sort_trim(&f2stack);

#if PRINT_PCT
  fin_pct = (100.0 * (mpz_sizeinbase(F1,2) + mpz_sizeinbase(F2,2))) / (mpz_sizeinbase(nm1,2) + mpz_sizeinbase(np1,2));
  printf("%6.2f .. ", fin_pct);  fflush(stdout);
  printf("\n");  fflush(stdout);
#endif

  /* Check the theorems we have available */

  if (bls_theorem7_limit(n, F1, R1, B1, t, u, r, s)) {
    if (get_verbose_level() > 0) printf("BLS75 proof using N-1\n");
    trim_factors(F1, R1, n, nm1, B1, &f1stack, &bls_theorem7_limit, t, u, r, s);
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


  if (bls_theorem19_limit(n, F2, R2, B1, t, u, r, s)) {
    IV qv[MAXQV];
    int numqv = 0;
    if (get_verbose_level() > 0) printf("BLS75 proof using N+1\n");
    trim_factors(F2, R2, n, np1, B1, &f2stack, &bls_theorem19_limit, t, u, r, s);
    for (pcount = 0; success > 0 && pcount < f2stack.cur; pcount++)
      success = _verify_cond_III_q(n, f2stack.stack[pcount], qv, &numqv, 0, 0);
    if (success > 0 && (mpz_mul(t,F2,F2), mpz_add_ui(t,t,1), mpz_cmp(t,n) > 0))
      goto end_hybrid;  /* Corollary 8, n+1 factored more than sqrt(n)+1 */
    if (success > 0 && !bls_theorem17_limit(n, F2, R2, B1, t, u, r, s))
      success = _verify_cond_III_q(n, R2, qv, &numqv, 0, 0);
    if (success > 0) {
      mpz_mul(t, r, r);
      mpz_submul_ui(t, s, 8);   /* t = r^2 - 8s */
      /* N is prime if and only if s=0 OR t not a perfect square */
      success = (mpz_sgn(s) == 0 || !mpz_perfect_square_p(t))  ?  1  :  -1;
    }
    goto end_hybrid;   /* Theorem 17 or 19 */
  }

  if (get_verbose_level() > 0) printf("BLS75 proof using N-1 / N+1 (T20)\n");

  /* Check N < B^3 F1*F2*F2/2  or  N < B^3 F1*F1*F2/2 */
  success = bls_theorem20_limit(n, R1, F1, F2, B1, m, t, u, r, s);

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
          success = (mpz_cmp_ui(t,1) > 0 || mpz_cmp(t,n) < 0) ? -1 : 0;
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
    IV qv[MAXQV];
    int numqv = 0;
    for (pcount = 0; success > 0 && pcount < f2stack.cur; pcount++)
      success = _verify_cond_III_q(n, f2stack.stack[pcount], qv, &numqv, 0, 0);
    if (success > 0)
      success = _verify_cond_III_q(n, R2, qv, &numqv, 0, 0);
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
  int e, success = 1;
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
    for (e = 0; !success && e <= effort; e++)
      success = tfe(f, q, e);
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
  int e, success = 1;
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
    for (e = 0; !success && e <= effort; e++)
      success = tfe(f, p, e);
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
