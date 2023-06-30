#include <gmp.h>
#include "ptypes.h"

#include "lucas_seq.h"
#include "utility.h"


void lucasuv(mpz_t Uh, mpz_t Vl, mpz_t P, mpz_t Q, mpz_t k)
{
  mpz_t Vh, Ql, Qh, t;
  int j, s, n;

  if (mpz_sgn(k) <= 0) {
    mpz_set_ui(Uh, 0);
    mpz_set_ui(Vl, 2);
    return;
  }

  mpz_set_ui(Uh, 1);
  mpz_set_ui(Vl, 2);
  mpz_init_set(Vh, P);
  mpz_init(t);

  s = mpz_scan1(k, 0);     /* number of zero bits at the end */
  n = mpz_sizeinbase(k,2);

  /* It is tempting to try to pull out the various Q operations when Q=1 or
   * Q=-1.  This doesn't lead to any immediate savings.  Don't bother unless
   * there is a way to reduce the actual operations involving U and V. */
  mpz_init_set_ui(Ql,1);
  mpz_init_set_ui(Qh,1);

  for (j = n-1; j > s; j--) {
    mpz_mul(Ql, Ql, Qh);
    if (mpz_tstbit(k, j)) {
      mpz_mul(Qh, Ql, Q);
      mpz_mul(Uh, Uh, Vh);
      mpz_mul(t, Ql, P);  mpz_mul(Vl, Vl, Vh); mpz_sub(Vl, Vl, t);
      mpz_mul(Vh, Vh, Vh); mpz_sub(Vh, Vh, Qh); mpz_sub(Vh, Vh, Qh);
    } else {
      mpz_set(Qh, Ql);
      mpz_mul(Uh, Uh, Vl);  mpz_sub(Uh, Uh, Ql);
      mpz_mul(t, Ql, P);  mpz_mul(Vh, Vh, Vl); mpz_sub(Vh, Vh, t);
      mpz_mul(Vl, Vl, Vl);  mpz_sub(Vl, Vl, Ql);  mpz_sub(Vl, Vl, Ql);
    }
  }
  mpz_mul(Ql, Ql, Qh);
  mpz_mul(Qh, Ql, Q);
  mpz_mul(Uh, Uh, Vl);  mpz_sub(Uh, Uh, Ql);
  mpz_mul(t, Ql, P);  mpz_mul(Vl, Vl, Vh);  mpz_sub(Vl, Vl, t);
  mpz_mul(Ql, Ql, Qh);
  mpz_clear(Qh);  mpz_clear(t);  mpz_clear(Vh);
  for (j = 0; j < s; j++) {
    mpz_mul(Uh, Uh, Vl);
    mpz_mul(Vl, Vl, Vl);  mpz_sub(Vl, Vl, Ql);  mpz_sub(Vl, Vl, Ql);
    mpz_mul(Ql, Ql, Ql);
  }
  mpz_clear(Ql);
}

/* No argument checking is done in this internal function. */
void internal_lucas_vmod_q1(mpz_t V, mpz_t W, mpz_t P, mpz_t k, mpz_t n)
{
  int b = mpz_sizeinbase(k, 2);

  mpz_mod(P, P, n);
  mpz_set(V, P);
  mpz_mul(W, P, P);
  mpz_sub_ui(W, W, 2);  /*   V_{k} = V = P,  V_{k+1} = W = P*P-2   */
  mpz_mod(W, W, n);
  while (b-- > 1) {
    if (mpz_tstbit(k, b-1)) {
      mpz_mul(V, V, W);  mpz_sub(V, V, P);
      mpz_mul(W, W, W);  mpz_sub_ui(W, W, 2);
    } else {
      mpz_mul(W, V, W);  mpz_sub(W, W, P);
      mpz_mul(V, V, V);  mpz_sub_ui(V, V, 2);
    }
    mpz_mod(V, V, n);  mpz_mod(W, W, n);
  }
  /* V = V_{k} and W = V_{k+1} */
  /* Therefore U = D^-1 * (2W - V*P) where D = P*P-4Q */
}


void lucasuvmod(mpz_t U, mpz_t V, mpz_t P, mpz_t Q, mpz_t k, mpz_t n, mpz_t t)
{
#if 0
  { lucasuv(U, V, P, Q, k);   /* This will be horribly slow with large k */
    mpz_mod(U, U, n);
    mpz_mod(V, V, n);
    return;
  }
#endif
  mpz_t Pmod, Qmod, Dmod;
  unsigned long b;
  int qsignorig, qsign;

  if (mpz_cmp_ui(n, 1) < 0) croak("Lucas sequence modulus n must be > 0");
  if (mpz_cmp_ui(n, 1) == 0) { mpz_set_ui(V,0); return; }
  if (mpz_cmp_ui(k, 0) <= 0) { mpz_set_ui(U,0); mpz_set_ui(V,2); mpz_mod(V,V,n); return; }

  mpz_init(Pmod);
  mpz_init(Qmod);
  mpz_init(Dmod);
  mpz_mod(Pmod, P, n);
  mpz_mod(Qmod, Q, n);

#if 0
  /* A proof of concept of the idea from C&P page 147.
     They also how we can turn any Q w/gcd(Q,n)=1 into a square.
     We could run all gcd(Q,n)=1 cases through the fast V inversion code. */
  if (mpz_cmp_ui(Qmod,3) > 0 && mpz_perfect_square_p(Qmod)) {
    mpz_t root, invroot, Q1;
    mpz_init(root); mpz_init(invroot); mpz_init(Q1);
    mpz_sqrt(root, Qmod);
    mpz_invert(invroot, root, n);
    mpz_set_ui(Q1, 1);
    mpz_mul(invroot, invroot, Pmod);

    lucasuvmod(U, V, invroot, Q1, k, n, t);

    mpz_sub_ui(t,k,1);
    mpz_powm(t, root, t, n);
    mpz_mul(U, U, t);
    mpz_mod(U, U, n);
    mpz_mul(V, V, t);
    mpz_mul(V, V, root);
    mpz_mod(V, V, n);

    mpz_clear(root);  mpz_clear(invroot);  mpz_clear(Q1);
    return;
  }
#endif

  mpz_sub_ui(t, n, 1);
  qsignorig = qsign = !mpz_cmp_si(Qmod,1) ? 1 : !mpz_cmp(Qmod,t) ? -1 : 0;

  mpz_mul(Dmod, Pmod, Pmod);
  mpz_submul_ui(Dmod, Qmod, 4);
  mpz_mod(Dmod, Dmod, n);    /* D = (P*P - 4*Q) mod n */

  mpz_set_ui(t, 2);
  if (mpz_sgn(Dmod) == 0 && mpz_invert(t, t, n)) {
    mpz_mul(t, t, Pmod);
    mpz_mod(U, t, n);     /* S is P/2 mod n */
    mpz_set(V, U);
    mpz_sub_ui(t, k, 1);
    mpz_powm(U, U, t, n);
    mpz_mul(U, U, k);
    mpz_mod(U, U, n);     /* U = (k * S^(k-1)) mod n */

    mpz_powm(V, V, k, n);
    mpz_mul_ui(V, V, 2);
    mpz_mod(V, V, n);     /* V = (2 * S^k) mod n */

    mpz_clear(Dmod); mpz_clear(Qmod); mpz_clear(Pmod);
    return;
  }

  /* Mulmods used:
   *
   * Q=1 invertible      2 even 2 odd
   * odd N, Q=1 or Q=-1  2 even 5 odd
   * odd N  any P/Q      3 even 7 odd
   * generic             5 even 6 odd + 6 + 3 per trailing zero
   */

  b = mpz_sizeinbase(k, 2);

  if (qsign == 1 && mpz_invert(t, Dmod, n)) {
    /* Compute V_k and V_{k+1}, then compute U_k from them. */
    mpz_set(V, Pmod);
    mpz_mul(U, V, V);
    mpz_sub_ui(U, U, 2);  /*   V = P; U = P*P-2   */
    while (b-- > 1) {
      if (mpz_tstbit(k, b-1)) {
        mpz_mul(V, V, U);  mpz_sub(V, V, Pmod);  mpz_mod(V, V, n);
        mpz_mul(U, U, U);  mpz_sub_ui(U, U, 2);  mpz_mod(U, U, n);
      } else {
        mpz_mul(U, V, U);  mpz_sub(U, U, Pmod);  mpz_mod(U, U, n);
        mpz_mul(V, V, V);  mpz_sub_ui(V, V, 2);  mpz_mod(V, V, n);
      }
    }
    mpz_mul_ui(U, U, 2);
    mpz_submul(U, V, Pmod);
    mpz_mul(U, U, t);
  } else if (mpz_odd_p(n) && qsign) {   /* Odd N with Q=1 or Q=-1 */

    mpz_set_ui(U, 1);
    mpz_set   (V, Pmod);
    while (b-- > 1) {
      mpz_mulmod(U, U, V, n, t);     /* U2k = Uk * Vk */

      mpz_mul(V, V, V);
      if (qsign == 1)  mpz_sub_ui(V, V, 2);  else  mpz_add_ui(V, V, 2);
      mpz_mod(V, V, n);              /* V2k = Vk^2 - 2 Q^k */

      qsign = 1;
      if (mpz_tstbit(k, b-1)) {
        mpz_mul(t, U, Dmod);
                                     /* U:  U2k+1 = (P*U2k + V2k)/2 */
        mpz_mul(U, U, Pmod);
        mpz_add(U, U, V);
        if (mpz_odd_p(U)) mpz_add(U, U, n);
        mpz_fdiv_q_2exp(U, U, 1);
                                     /* V:  V2k+1 = (D*U2k + P*V2k)/2 */
        mpz_mul(V, V, Pmod);
        mpz_add(V, V, t);
        if (mpz_odd_p(V)) mpz_add(V, V, n);
        mpz_fdiv_q_2exp(V, V, 1);

        qsign = qsignorig;
      }
    }

  } else if (mpz_odd_p(n)) {   /* Odd N with arbitrary Q */

    mpz_t Qk;
    mpz_init_set(Qk, Qmod);
    mpz_set_ui(U, 1);
    mpz_set   (V, Pmod);

    while (b-- > 1) {
      mpz_mulmod(U, U, V, n, t);     /* U2k = Uk * Vk */

      mpz_mul(V, V, V);
      mpz_submul_ui(V, Qk, 2);
      mpz_mod(V, V, n);              /* V2k = Vk^2 - 2 Q^k */

      mpz_mul(Qk, Qk, Qk);           /* Q2k = Qk^2 */
      if (mpz_tstbit(k, b-1)) {
        mpz_mul(t, U, Dmod);
                                     /* U:  U2k+1 = (P*U2k + V2k)/2 */
        mpz_mul(U, U, Pmod);
        mpz_add(U, U, V);
        if (mpz_odd_p(U)) mpz_add(U, U, n);
        mpz_fdiv_q_2exp(U, U, 1);
                                     /* V:  V2k+1 = (D*U2k + P*V2k)/2 */
        mpz_mul(V, V, Pmod);
        mpz_add(V, V, t);
        if (mpz_odd_p(V)) mpz_add(V, V, n);
        mpz_fdiv_q_2exp(V, V, 1);

        mpz_mul(Qk, Qk, Qmod);
      }
      mpz_mod(Qk, Qk, n);
    }
    mpz_clear(Qk);

  } else {    /* Even N.  This is generic code for any input. */
    mpz_t Uh, Vl, Vh, Ql, Qh;
    int j, s = mpz_scan1(k, 0);

    mpz_init_set_ui(Uh, 1);
    mpz_init_set_ui(Vl, 2);
    mpz_init_set   (Vh, Pmod);
    mpz_init_set_ui(Ql, 1);
    mpz_init_set_ui(Qh, 1);

    for (j = b-1; j > s; j--) {
      mpz_mul(Ql, Ql, Qh);
      mpz_mod(Ql, Ql, n);
      if (mpz_tstbit(k, j)) {
        mpz_mul(Qh, Ql, Qmod);
        mpz_mul(Uh, Uh, Vh);
        mpz_mul(t, Ql, Pmod); mpz_mul(Vl, Vl, Vh);  mpz_sub(Vl, Vl, t);
        mpz_mul(Vh, Vh, Vh);  mpz_sub(Vh, Vh, Qh);  mpz_sub(Vh, Vh, Qh);
      } else {
        mpz_set(Qh, Ql);
        mpz_mul(Uh, Uh, Vl);  mpz_sub(Uh, Uh, Ql);
        mpz_mul(t, Ql, Pmod); mpz_mul(Vh, Vh, Vl);  mpz_sub(Vh, Vh, t);
        mpz_mul(Vl, Vl, Vl);  mpz_sub(Vl, Vl, Ql);  mpz_sub(Vl, Vl, Ql);
      }
      mpz_mod(Qh, Qh, n);
      mpz_mod(Uh, Uh, n);
      mpz_mod(Vh, Vh, n);
      mpz_mod(Vl, Vl, n);
    }
    mpz_mul(Ql, Ql, Qh);
    mpz_mul(Qh, Ql, Qmod);
    mpz_mul(Uh, Uh, Vl);  mpz_sub(Uh, Uh, Ql);
    mpz_mul(t, Ql, Pmod);  mpz_mul(Vl, Vl, Vh);  mpz_sub(Vl, Vl, t);
    mpz_mul(Ql, Ql, Qh);
    mpz_clear(Qh);  mpz_clear(Vh);
    mpz_mod(Ql, Ql, n);
    mpz_mod(Uh, Uh, n);
    mpz_mod(Vl, Vl, n);
    for (j = 0; j < s; j++) {
      mpz_mul(Uh, Uh, Vl);
      mpz_mul(Vl, Vl, Vl);  mpz_sub(Vl, Vl, Ql);  mpz_sub(Vl, Vl, Ql);
      mpz_mul(Ql, Ql, Ql);
      mpz_mod(Ql, Ql, n);
      mpz_mod(Uh, Uh, n);
      mpz_mod(Vl, Vl, n);
    }
    mpz_set(U, Uh);
    mpz_set(V, Vl);
    mpz_clear(Ql); mpz_clear(Vl); mpz_clear(Uh);
  }

  mpz_mod(U, U, n);
  mpz_mod(V, V, n);
  mpz_clear(Dmod); mpz_clear(Qmod); mpz_clear(Pmod);
}


/* Returns Lucas sequence  U_k mod n and V_k mod n  defined by P,Q */
void lucas_seq(mpz_t U, mpz_t V, mpz_t n, IV P, IV Q, mpz_t k,
               mpz_t Qk, mpz_t t)
{
  mpz_t mP, mQ;
  mpz_init(mP);  mpz_init(mQ);
  mpz_set_iv(mP, P);
  mpz_set_iv(mQ, Q);
  lucasuvmod(U, V, mP, mQ, k, n, t);
  mpz_powm(Qk, mQ, k, n);
  mpz_clear(mQ); mpz_clear(mP);
}

#if 1
void lucasvmod(mpz_t V, mpz_t P, mpz_t Q, mpz_t k, mpz_t n, mpz_t t)
{
  mpz_t U, Pmod, Qmod, Dmod;

  /* This doesn't save much vs. calling lucasuvmod, but it's a tiny help
   * with the Q=1 case. */

  if (mpz_cmp_ui(n, 1) < 0) croak("Lucas sequence modulus n must be > 0");
  if (mpz_cmp_ui(n, 1) == 0) { mpz_set_ui(V,0); return; }
  if (mpz_cmp_ui(k, 0) <= 0) { mpz_set_ui(V,2); mpz_mod(V, V, n); return; }

  mpz_init(Pmod);
  mpz_init(Qmod);
  mpz_init(Dmod);
  mpz_mod(Pmod, P, n);
  mpz_mod(Qmod, Q, n);

  mpz_mul(Dmod, Pmod, Pmod);
  mpz_submul_ui(Dmod, Qmod, 4);
  mpz_mod(Dmod, Dmod, n);    /* D = (P*P - 4*Q) mod n */

  mpz_set_ui(t, 2);
  if (mpz_sgn(Dmod) == 0 && mpz_invert(t, t, n)) {
    mpz_mul(t, t, Pmod);
    mpz_mod(t, t, n);     /* S is P/2 mod n */

    mpz_powm(V, t, k, n);
    mpz_mul_ui(V, V, 2);
    mpz_mod(V, V, n);     /* V = (2 * S^k) mod n */

    mpz_clear(Dmod); mpz_clear(Qmod); mpz_clear(Pmod);
    return;
  }

  if (!mpz_cmp_ui(Q,1)) {
    internal_lucas_vmod_q1(V, t, Pmod, k, n);
  } else {
    mpz_init(U);
    lucasuvmod(U, V, P, Q, k, n, t);
    mpz_clear(U);
  }
}
#else
void lucasvmod(mpz_t V, mpz_t P, mpz_t Q, mpz_t k, mpz_t n, mpz_t t)
{
  mpz_t U;
  mpz_init(U);
  lucasuvmod(U, V, P, Q, k, n, t);
  mpz_clear(U);
}
#endif
void lucasumod(mpz_t U, mpz_t P, mpz_t Q, mpz_t k, mpz_t n, mpz_t t)
{
  mpz_t V;
  mpz_init(V);
  lucasuvmod(U, V, P, Q, k, n, t);
  mpz_clear(V);
}
