
/* If set, then tries to apply theorem 7 in addition to theorem 5.
 * Normally I would just have this on, but then we'd produce certificates
 * that Math::Prime::Util 0.26 couldn't understand. :(
 */
#define BLS_THEOREM7 0

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

/* FIXME:
 *   (1) too much repetitious overhead code in these
 *   (2) way too much copy/paste between Pocklington and BLS
 *   (3) Pocklington has code rotted, so fix before using
 */
#define PRIM_STACK_SIZE 128

#define primality_handle_factor(f, primality_func, factor_prob) \
  { \
    int f_prob_prime = _GMP_is_prob_prime(f); \
    if ( (f_prob_prime == 1) && (primality_func(f, effort, prooftextptr) == 2) ) \
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

static int bls_theorem5_limit(mpz_t n, mpz_t A, mpz_t B,
                              mpz_t t, mpz_t y, mpz_t r, mpz_t s)
{
  mpz_mul(t, A, B);
  mpz_add_ui(t, t, 1);
  if (mpz_cmp(t, n) != 0) croak("BLS75:  A*B != n\n");

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

#if BLS_THEOREM7
static int bls_theorem7_limit(mpz_t n, mpz_t A, mpz_t B, UV B1,
                              mpz_t t, mpz_t y, mpz_t r, mpz_t s)
{
  mpz_mul(t, A, B);
  mpz_add_ui(t, t, 1);
  if (mpz_cmp(t, n) != 0) croak("BLS75:  A*B != n-1\n");

  mpz_mul_ui(t, A, 2);
  mpz_tdiv_qr(s, r, B, t);

  mpz_mul(y, t, A);     /* y = 2*A*A              */
  mpz_sub_ui(t, r, B1); /* t = r-B1               */
  mpz_mul(t, t, A);     /* t = A*(r-B1)           */
  mpz_add(y, y, t);     /* y = 2A^2 + A(r-B1)     */
  mpz_add_ui(y, y, 1);  /* y = 2A^2 + A(r-B1) + 1 */
  mpz_mul_ui(t, A, B1); /* t = A*B1               */
  mpz_add_ui(t, t, 1);  /* t = A*B1+1             */
  mpz_mul(y, y, t);     /* y = (A*B1+1)*(2A^2+(r-B1)A+1) */

  return (mpz_cmp(n, y) < 0) ? 1 : 0;
}
#endif

int _GMP_primality_bls_nm1(mpz_t n, int effort, char** prooftextptr)
{
  mpz_t nm1, A, B, t, m, f, r, s;
  mpz_t mstack[PRIM_STACK_SIZE];
  mpz_t fstack[PRIM_STACK_SIZE];
  int msp = 0;
  int fsp = 0;
  int success = 1;
  int theorem7 = 0;
  UV B1 = 2000;

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
    primality_handle_factor(f, _GMP_primality_bls_nm1, 1);
  }

  while (success) {

    if (bls_theorem5_limit(n, A, B, t, m, r, s))
      break;

#if BLS_THEOREM7
    /* With theorem 7, we can exit even earlier if we know the lower limit
     * on the size of the factors of B.  We need to do an additional test
     * at the end (treat B as if it was the last factor). */
    if (bls_theorem7_limit(n, A, B, B1, t, m, r, s)) {
      theorem7 = 1;
      break;
    }
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
    if (effort >= 1) {
      if (!success)  success = _GMP_pminus1_factor(m, f, 100, 1000);
      if (!success)  success = _GMP_pminus1_factor(m, f, 1000, 10000);
    }
    if (!success && effort == 2) {
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
      if (!success && log2m < 80)  success = _GMP_ECM_FACTOR(m, f, 150, 5);
      if (!success)  success = _GMP_pbrent_factor(m, f, 3, brent_rounds);
      if (!success && final_B2 > 10000)  success = _GMP_pminus1_factor(m, f, 10000, final_B2);
    }
    if (!success && effort >= 3) {
      if (!success)  success = _GMP_pminus1_factor(m, f, 10000, 200000);
      if (!success)  success = _GMP_ECM_FACTOR(m, f, 500, 30);
      if (!success)  success = _GMP_ECM_FACTOR(m, f, 2000, 20);
    }
    if (!success && effort >= 4) {
      if (!success)  success = _GMP_pminus1_factor(m, f, 200000, 4000000);
      if (!success)  success = _GMP_ECM_FACTOR(m, f, 10000, 10);
    }
#if BLS_THEOREM7
    /* Could we exit now if B1 was larger? */
    if (!success && effort > 2) {
      UV newB1 = B1;
      if      (bls_theorem7_limit(n,A,B,    10000,t,f,r,s)) newB1 = 10000;
      else if (bls_theorem7_limit(n,A,B,  1000000,t,f,r,s)) newB1 = 1000000;
      else if (bls_theorem7_limit(n,A,B,100000000,t,f,r,s)) newB1 = 100000000;
      if (newB1 > B1) {
        UV tf = _GMP_trial_factor(m, B1, newB1);
        if (tf == 0) {
          B1 = newB1;
          continue; /* Go back to the top and we'll retest with this B1 */
        } else {
          /* Holy trial factorization, Batman! */
          B1 = tf;
          mpz_set_ui(f, tf);
          success = 1;
        }
      }
    }
#endif
    if (!success && effort > 5) {  /* do here only if effort > 5 */
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
                primality_handle_factor(farray[i], _GMP_primality_bls_nm1, 0);
              }
            }
            for (i = 0; i < 66; i++)
              mpz_clear(farray[i]);
            if (success)
              continue;
      }
    }
    if (!success && effort >= 5) {
      UV i;
      UV B1 = 10000;
      UV curves = 10;
      for (i = 1; i < 18; i++) {
        if ((4+i) > (UV)effort) break;
        B1 *= 2;
        success = _GMP_ECM_FACTOR(m, f, B1, curves);
        if (success) break;
      }
    }
    /* If we couldn't factor m and the stack is empty, we've failed. */
    if ( (!success) && (msp == 0) )
      break;
    /* Put the two factors f and m/f into the stacks, smallest first */
    mpz_divexact(m, m, f);
    if (mpz_cmp(m, f) < 0) { mpz_set(t, m); mpz_set(m, f); mpz_set(f, t); }
    primality_handle_factor(f, _GMP_primality_bls_nm1, 0);
    primality_handle_factor(m, _GMP_primality_bls_nm1, 0);
  }

  /* clear mstack since we don't care about it.  Use to hold a values. */
  while (msp-- > 0)
    mpz_clear(mstack[msp]);
  msp = 0;

  if (success > 0) {
    /* Calculate r and s from final A/B values */
    mpz_mul_ui(t, A, 2);
    mpz_tdiv_qr(s, r, B, t);

    mpz_mul(t, r, r);
    mpz_submul_ui(t, s, 8);   /* t = r^2 - 8s */
    /* N is prime if and only if s=0 OR t not a perfect square */
    /* So, N is a composite if s != 0 and t is perfect square */
    if (mpz_sgn(s) != 0 && mpz_perfect_square_p(t))
      success = -1;
  }

  if (success > 0) {
    int pcount, a;
    int const alimit = (effort <= 2) ? 200 : 10000;
    mpz_t p, ap;

    mpz_init(p);
    mpz_init(ap);

    if (theorem7) {
      /* Theorem 7 used, so add B to the factor list so it gets checked */
      if (fsp >= PRIM_STACK_SIZE) { croak("BLS75 stack overflow\n"); }
      mpz_init_set(fstack[fsp++], B);
    }

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
  if (success > 0 && prooftextptr != 0) {
    int i, prooflen;
    char *proofstr, *proofptr;
    if (fsp != msp) croak("Different f and a counts\n");
    prooflen = mpz_sizeinbase(n, 10) * (1 + fsp + msp) + 200;
    New(0, proofstr, prooflen, char);
    proofptr = proofstr;
    *proofptr = 0;
    proofptr += gmp_sprintf(proofptr, "%Zd :", n);
    if (theorem7) {
      proofptr += gmp_sprintf(proofptr, " N-1 T7 :");
      proofptr += gmp_sprintf(proofptr, " %lu %Zd %Zd :",
                       B1, fstack[--fsp], mstack[--msp]);
       mpz_clear(fstack[fsp]);  mpz_clear(mstack[msp]);
    } else {
      proofptr += gmp_sprintf(proofptr, " N-1 T5 :");
    }
    for (i = 0; i < fsp; i++)
      proofptr += gmp_sprintf(proofptr, " %Zd", fstack[i]);
    proofptr += gmp_sprintf(proofptr, " :");
    for (i = 0; i < msp; i++)
      proofptr += gmp_sprintf(proofptr, " %Zd", mstack[i]);
    proofptr += gmp_sprintf(proofptr, "\n");
    /* Set or append */
    if (*prooftextptr == 0) {
      *prooftextptr = proofstr;
    } else {
      Renew(*prooftextptr, strlen(*prooftextptr) + strlen(proofstr) + 1, char);
      (void) strcat(*prooftextptr, proofstr);
      Safefree(proofstr);
    }
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
  if (success < 0) return 0;
  if (success > 0) return 2;
  return 1;
}
