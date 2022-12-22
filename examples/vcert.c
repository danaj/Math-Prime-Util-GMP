/*
 * Verify Cert
 * version 0.96
 *
 * Copyright (c) 2013-2016 Dana Jacobsen (dana@acm.org).
 * This is free software; you can redistribute it and/or modify it under
 * the same terms as the Perl 5 programming language system itself.
 *
 * Verifies Primo v3, Primo v4, and MPU certificates.
 *
 * Return values:
 *   0  all numbers are verified prime.
 *   1  at least one number was verified composite.
 *   2  the certificate does not provide a complete proof.
 *   3  there is an error in the certificate.
 *
 * TODO: Allow multiple proofs per input file
 * TODO: Projective EC for ~4x faster operation
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <gmp.h>

/*****************************************************************************/
/* Preliminary definitions                                                   */
/*****************************************************************************/

/* Projective doesn't work yet */
#define USE_AFFINE_EC 1

#define RET_PRIME 0
#define RET_COMPOSITE 1
#define RET_INVALID 2
#define RET_ERROR 3

#define CERT_UNKNOWN 0
#define CERT_PRIMO   1
#define CERT_PRIMO42 2
#define CERT_MPU     3

#define MAX_LINE_LEN 60000
#define MAX_STEPS    20000
#define MAX_QARRAY   100
#define BAD_LINES_ALLOWED  5    /* Similar to WraithX's verifier */

typedef unsigned long UV;
typedef   signed long IV;
#define croak(fmt,...)  { gmp_printf(fmt,##__VA_ARGS__); exit(RET_ERROR); }
#define MPUassert(c,text) if (!(c)) { croak("Internal error: " text); }

#define BGCD_PRIMES      168
#define BGCD_LASTPRIME   997
#define BGCD_NEXTPRIME  1009

void GMP_pn_primorial(mpz_t prim, UV n);
UV trial_factor(mpz_t n, UV from_n, UV to_n);
int miller_rabin_ui(mpz_t n, UV base);
int miller_rabin(mpz_t n, mpz_t a);
void lucas_seq(mpz_t U, mpz_t V, mpz_t n, IV P, IV Q, mpz_t k, mpz_t Qk, mpz_t t);
#define mpz_mulmod(r, a, b, n, t)  \
  do { mpz_mul(t, a, b); mpz_mod(r, t, n); } while (0)

/*****************************************************************************/
/* Some global variables and functions we'll use                             */
/*****************************************************************************/

int _verbose = 0;
int _quiet = 0;
int _testcount = 0;
int _base = 10;
int _step = 0;
int _format = CERT_UNKNOWN;
char _line[MAX_LINE_LEN+1];
char _vstr[MAX_LINE_LEN+1];
const char* _filename;
FILE* _fh;
mpz_t PROOFN, N, A, B, M, Q, X, Y, LQ, LP, S, R, T, J, W, T1, T2;
mpz_t QARRAY[MAX_QARRAY];
mpz_t AARRAY[MAX_QARRAY];
mpz_t _bgcd;

int   _num_chains = 0;
mpz_t _chain_n[MAX_STEPS];
mpz_t _chain_q[MAX_STEPS];

static void var_init(void) {
  int i;
  mpz_init(PROOFN);
  mpz_init(N);
  mpz_init(A);
  mpz_init(B);
  /* MPU: */
  mpz_init(M);
  mpz_init(Q);
  mpz_init(X);
  mpz_init(Y);
  mpz_init(LQ);
  mpz_init(LP);
  /* Primo: */
  mpz_init(S);
  mpz_init(R);
  mpz_init(T);
  mpz_init(J);
  mpz_init(W);

  mpz_init(_bgcd);
  GMP_pn_primorial(_bgcd, BGCD_PRIMES);

  mpz_init(T1);
  mpz_init(T2);

  for (i = 0; i < MAX_QARRAY; i++) {
    mpz_init(QARRAY[i]);
    mpz_init(AARRAY[i]);
  }
}
static void var_free(void) {
  int i;
  mpz_clear(PROOFN);
  mpz_clear(N);
  mpz_clear(A);
  mpz_clear(B);
  mpz_clear(M);
  mpz_clear(Q);
  mpz_clear(X);
  mpz_clear(Y);
  mpz_clear(LQ);
  mpz_clear(LP);
  mpz_clear(S);
  mpz_clear(R);
  mpz_clear(T);
  mpz_clear(J);
  mpz_clear(W);

  mpz_clear(_bgcd);

  mpz_clear(T1);
  mpz_clear(T2);

  for (i = 0; i < MAX_QARRAY; i++) {
    mpz_clear(QARRAY[i]);
    mpz_clear(AARRAY[i]);
  }
}

static void quit_prime(void) {
  if (!_quiet) printf("                                             \r");
  if (!_quiet) printf("PRIME\n");
  var_free();
  exit(RET_PRIME);
}
static void quit_composite(void) {
  if (!_quiet) printf("                                             \r");
  if (!_quiet) printf("COMPOSITE\n");
  var_free();
  exit(RET_COMPOSITE);
}
static void quit_invalid(const char* type, const char* msg) {
  if (!_quiet) printf("\n");
  if (!_quiet) gmp_printf("%s: step %d, %Zd failed condition %s\n", type, _step, N, msg);
  var_free();
  exit(RET_INVALID);
}
static void quit_error(const char* msg1, const char* msg2) {
  if (!_quiet) printf("\n");
  if (!_quiet) gmp_printf("ERROR: step %d, %s%s\n", _step, msg1, msg2);
  var_free();
  exit(RET_ERROR);
}

#if USE_AFFINE_EC
/*****************************************************************************/
/* EC: affine with point (x,y,1)                                             */
/*****************************************************************************/

struct ec_affine_point  { mpz_t x, y; };
/* P3 = P1 + P2 */
static void _ec_add_AB(mpz_t n,
                    struct ec_affine_point P1,
                    struct ec_affine_point P2,
                    struct ec_affine_point *P3,
                    mpz_t m,
                    mpz_t t1,
                    mpz_t t2)
{
  if (!mpz_cmp(P1.x, P2.x)) {
    mpz_add(t2, P1.y, P2.y);
    mpz_mod(t1, t2, n);
    if (!mpz_cmp_ui(t1, 0) ) {
      mpz_set_ui(P3->x, 0);
      mpz_set_ui(P3->y, 1);
      return;
    }
  }

  mpz_sub(t1, P2.x, P1.x);
  mpz_mod(t2, t1, n);

  /* t1 = 1/deltay mod n */
  if (!mpz_invert(t1, t2, n)) {
    /* We've found a factor!  In multiply, gcd(mult,n) will be a factor. */
    mpz_set_ui(P3->x, 0);
    mpz_set_ui(P3->y, 1);
    return;
  }

  mpz_sub(m, P2.y, P1.y);
  mpz_mod(t2, m, n);        /* t2 = deltay   mod n */
  mpz_mul(m, t1, t2);
  mpz_mod(m, m, n);         /* m = deltay / deltax   mod n */

  /* x3 = m^2 - x1 - x2 mod n */
  mpz_mul(t1, m, m);
  mpz_sub(t2, t1, P1.x);
  mpz_sub(t1, t2, P2.x);
  mpz_mod(P3->x, t1, n);
  /* y3 = m(x1 - x3) - y1 mod n */
  mpz_sub(t1, P1.x, P3->x);
  mpz_mul(t2, m, t1);
  mpz_sub(t1, t2, P1.y);
  mpz_mod(P3->y, t1, n);
}

/* P3 = 2*P1 */
static void _ec_add_2A(mpz_t a,
                    mpz_t n,
                    struct ec_affine_point P1,
                    struct ec_affine_point *P3,
                    mpz_t m,
                    mpz_t t1,
                    mpz_t t2)
{
  /* m = (3x1^2 + a) * (2y1)^-1 mod n */
  mpz_mul_ui(t1, P1.y, 2);
  if (!mpz_invert(m, t1, n)) {
    mpz_set_ui(P3->x, 0);
    mpz_set_ui(P3->y, 1);
    return;
  }
  mpz_mul_ui(t1, P1.x, 3);
  mpz_mul(t2, t1, P1.x);
  mpz_add(t1, t2, a);
  mpz_mul(t2, m, t1);
  mpz_tdiv_r(m, t2, n);

  /* x3 = m^2 - 2x1 mod n */
  mpz_mul(t1, m, m);
  mpz_mul_ui(t2, P1.x, 2);
  mpz_sub(t1, t1, t2);
  mpz_tdiv_r(P3->x, t1, n);

  /* y3 = m(x1 - x3) - y1 mod n */
  mpz_sub(t1, P1.x, P3->x);
  mpz_mul(t2, t1, m);
  mpz_sub(t1, t2, P1.y);
  mpz_tdiv_r(P3->y, t1, n);
}

static int ec_affine_multiply(mpz_t a, mpz_t k, mpz_t n, struct ec_affine_point P, struct ec_affine_point *R, mpz_t d)
{
  int found = 0;
  struct ec_affine_point A, B, C;
  mpz_t t, t2, t3, mult;

  mpz_init(A.x); mpz_init(A.y);
  mpz_init(B.x); mpz_init(B.y);
  mpz_init(C.x); mpz_init(C.y);
  mpz_init(t);   mpz_init(t2);   mpz_init(t3);
  mpz_init_set_ui(mult, 1);  /* holds intermediates, gcd at end */

  mpz_set(A.x, P.x);  mpz_set(A.y, P.y);
  mpz_set_ui(B.x, 0); mpz_set_ui(B.y, 1);

  /* Binary ladder multiply. */
  while (mpz_cmp_ui(k, 0) > 0) {
    if (mpz_odd_p(k)) {
      mpz_sub(t, B.x, A.x);
      mpz_mul(t2, mult, t);
      mpz_mod(mult, t2, n);
      if ( !mpz_cmp_ui(A.x, 0) && !mpz_cmp_ui(A.y, 1) ) {
        /* nothing */
      } else if ( !mpz_cmp_ui(B.x, 0) && !mpz_cmp_ui(B.y, 1) ) {
        mpz_set(B.x, A.x);  mpz_set(B.y, A.y);
      } else {
        _ec_add_AB(n, A, B, &C, t, t2, t3);
        /* If the add failed to invert, then we have a factor. */
        mpz_set(B.x, C.x);  mpz_set(B.y, C.y);
      }
      mpz_sub_ui(k, k, 1);
    } else {
      mpz_mul_ui(t, A.y, 2);
      mpz_mul(t2, mult, t);
      mpz_mod(mult, t2, n);

      _ec_add_2A(a, n, A, &C, t, t2, t3);
      mpz_set(A.x, C.x);  mpz_set(A.y, C.y);
      mpz_tdiv_q_2exp(k, k, 1);
    }
  }
  mpz_gcd(d, mult, n);
  found = (mpz_cmp_ui(d, 1) && mpz_cmp(d, n));

  mpz_tdiv_r(R->x, B.x, n);
  mpz_tdiv_r(R->y, B.y, n);

  mpz_clear(mult);
  mpz_clear(t);   mpz_clear(t2);   mpz_clear(t3);
  mpz_clear(A.x); mpz_clear(A.y);
  mpz_clear(B.x); mpz_clear(B.y);
  mpz_clear(C.x); mpz_clear(C.y);

  return found;
}

#else

/*****************************************************************************/
/* EC: projective with point (X,1,Z) (Montgomery)                            */
/*****************************************************************************/

/* (xout:zout) = (x1:z1) + (x2:z2) */
static void pec_add3(mpz_t xout, mpz_t zout,
                     mpz_t x1, mpz_t z1,
                     mpz_t x2, mpz_t z2,
                     mpz_t xin, mpz_t zin,
                     mpz_t n, mpz_t u, mpz_t v, mpz_t w)
{
  mpz_sub(u, x2, z2);
  mpz_add(v, x1, z1);
  mpz_mulmod(u, u, v, n, w);     /* u = (x2 - z2) * (x1 + z1) % n */

  mpz_add(v, x2, z2);
  mpz_sub(w, x1, z1);
  mpz_mulmod(v, v, w, n, v);     /* v = (x2 + z2) * (x1 - z1) % n */

  mpz_add(w, u, v);              /* w = u+v */
  mpz_sub(v, u, v);              /* v = u-v */

  mpz_mulmod(w, w, w, n, u);     /* w = (u+v)^2 % n */
  mpz_mulmod(v, v, v, n, u);     /* v = (u-v)^2 % n */

  mpz_set(u, xin);
  mpz_mulmod(xout, w, zin, n, w);
  mpz_mulmod(zout, v, u,   n, w);
  /* 6 mulmods, 6 adds */
}

/* (x2:z2) = 2(x1:z1) */
static void pec_double(mpz_t x2, mpz_t z2, mpz_t x1, mpz_t z1,
                       mpz_t b, mpz_t n, mpz_t u, mpz_t v, mpz_t w)
{
  mpz_add(u, x1, z1);
  mpz_mulmod(u, u, u, n, w);     /* u = (x1+z1)^2 % n */

  mpz_sub(v, x1, z1);
  mpz_mulmod(v, v, v, n, w);     /* v = (x1-z1)^2 % n */

  mpz_mulmod(x2, u, v, n, w);    /* x2 = uv % n */

  mpz_sub(w, u, v);              /* w = u-v = 4(x1 * z1) */
  mpz_mulmod(u, b, w, n, z2);
  mpz_add(u, u, v);              /* u = (v+b*w) mod n */
  mpz_mulmod(z2, w, u, n, v);    /* z2 = (w*u) mod n */
  /* 5 mulmods, 4 adds */
}

#define NORMALIZE(f, u, v, x, z, n) \
    mpz_gcdext(f, u, NULL, z, n); \
    mpz_mulmod(x, x, u, n, v); \
    mpz_set_ui(z, 1);

static void pec_mult(mpz_t a, mpz_t b, mpz_t k, mpz_t n, mpz_t x, mpz_t z)
{
  mpz_t u, v, w, x1, x2, z1, z2, r;
  int l = -1;

  mpz_init(u); mpz_init(v); mpz_init(w);
  mpz_init(x1);  mpz_init(x2);  mpz_init(z1);  mpz_init(z2);

  mpz_sub_ui(k, k, 1);
  mpz_init_set(r, k);
  while (mpz_cmp_ui(r, 1) > 0) {
    mpz_tdiv_q_2exp(r, r, 1);
    l++;
  }
  mpz_clear(r);

  if (mpz_tstbit(k, l)) {
    pec_double(x2, z2, x, z, b, n, u, v, w);
    pec_add3(x1, z1, x2, z2, x, z, x, z, n, u, v, w);
    pec_double(x2, z2, x2, z2, b, n, u, v, w);
  } else {
    pec_double(x1, z1, x, z, b, n, u, v, w);
    pec_add3(x2, z2, x, z, x1, z1, x, z, n, u, v, w);
  }
  l--;
  while (l >= 1) {
    if (mpz_tstbit(k, l)) {
      pec_add3(x1, z1, x1, z1, x2, z2, x, z, n, u, v, w);
      pec_double(x2, z2, x2, z2, b, n, u, v, w);
    } else {
      pec_add3(x2, z2, x2, z2, x1, z1, x, z, n, u, v, w);
      pec_double(x1, z1, x1, z1, b, n, u, v, w);
    }
    l--;
  }
  if (mpz_tstbit(k, 0)) {
    pec_double(x, z, x2, z2, b, n, u, v, w);
  } else {
    pec_add3(x, z, x2, z2, x1, z1, x, z, n, u, v, w);
  }
  mpz_clear(u);  mpz_clear(v);  mpz_clear(w);
  mpz_clear(x1);  mpz_clear(x2);  mpz_clear(z1);  mpz_clear(z2);
}

#endif

/*****************************************************************************/
/* M-R, Lucas, BPSW                                                          */
/*****************************************************************************/
int miller_rabin_ui(mpz_t n, UV base)
{
  int rval;
  mpz_t a;
  mpz_init_set_ui(a, base);
  rval = miller_rabin(n, a);
  mpz_clear(a);
  return rval;
}

int miller_rabin(mpz_t n, mpz_t a)
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
  if (mpz_cmp_ui(a, 1) <= 0)
    croak("Base %ld is invalid", mpz_get_si(a));
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

/* Returns Lucas sequence  U_k mod n and V_k mod n  defined by P,Q */
void lucas_seq(mpz_t U, mpz_t V, mpz_t n, IV P, IV Q, mpz_t k,
               mpz_t Qk, mpz_t t)
{
  UV b = mpz_sizeinbase(k, 2);
  IV D = P*P - 4*Q;

  MPUassert( mpz_cmp_ui(n, 2) >= 0, "lucas_seq: n is less than 2" );
  MPUassert( mpz_cmp_ui(k, 0) >= 0, "lucas_seq: k is negative" );
  MPUassert( P >= 0 && mpz_cmp_si(n, P) >= 0, "lucas_seq: P is out of range" );
  MPUassert( mpz_cmp_si(n, Q) >= 0, "lucas_seq: Q is out of range" );
  MPUassert( D != 0, "lucas_seq: D is zero" );

  if (mpz_cmp_ui(k, 0) <= 0) {
    mpz_set_ui(U, 0);
    mpz_set_ui(V, 2);
    return;
  }

  MPUassert( mpz_odd_p(n), "lucas_seq: implementation is for odd n" );

  mpz_set_ui(U, 1);
  mpz_set_si(V, P);
  mpz_set_si(Qk, Q);

  if (Q == 1) {
    /* Use the fast V method if possible.  Much faster with small n. */
    mpz_set_si(t, P*P-4);
    if (P > 2 && mpz_invert(t, t, n)) {
      /* Compute V_k and V_{k+1}, then computer U_k from them. */
      mpz_set_si(V, P);
      mpz_init_set_si(U, P*P-2);
      while (b > 1) {
        b--;
        if (mpz_tstbit(k, b-1)) {
          mpz_mul(V, V, U);  mpz_sub_ui(V, V, P);  mpz_mod(V, V, n);
          mpz_mul(U, U, U);  mpz_sub_ui(U, U, 2);  mpz_mod(U, U, n);
        } else {
          mpz_mul(U, V, U);  mpz_sub_ui(U, U, P);  mpz_mod(U, U, n);
          mpz_mul(V, V, V);  mpz_sub_ui(V, V, 2);  mpz_mod(V, V, n);
        }
      }
      mpz_mul_ui(U, U, 2);
      mpz_submul_ui(U, V, P);
      mpz_mul(U, U, t);
    } else {
      /* Fast computation of U_k and V_k, specific to Q = 1 */
      while (b > 1) {
        mpz_mulmod(U, U, V, n, t);     /* U2k = Uk * Vk */
        mpz_mul(V, V, V);
        mpz_sub_ui(V, V, 2);
        mpz_mod(V, V, n);               /* V2k = Vk^2 - 2 Q^k */
        b--;
        if (mpz_tstbit(k, b-1)) {
          mpz_mul_si(t, U, D);
                                      /* U:  U2k+1 = (P*U2k + V2k)/2 */
          mpz_mul_si(U, U, P);
          mpz_add(U, U, V);
          if (mpz_odd_p(U)) mpz_add(U, U, n);
          mpz_fdiv_q_2exp(U, U, 1);
                                      /* V:  V2k+1 = (D*U2k + P*V2k)/2 */
          mpz_mul_si(V, V, P);
          mpz_add(V, V, t);
          if (mpz_odd_p(V)) mpz_add(V, V, n);
          mpz_fdiv_q_2exp(V, V, 1);
        }
      }
    }
  } else {
    while (b > 1) {
      mpz_mulmod(U, U, V, n, t);     /* U2k = Uk * Vk */
      mpz_mul(V, V, V);
      mpz_submul_ui(V, Qk, 2);
      mpz_mod(V, V, n);               /* V2k = Vk^2 - 2 Q^k */
      mpz_mul(Qk, Qk, Qk);            /* Q2k = Qk^2 */
      b--;
      if (mpz_tstbit(k, b-1)) {
        mpz_mul_si(t, U, D);
                                    /* U:  U2k+1 = (P*U2k + V2k)/2 */
        mpz_mul_si(U, U, P);
        mpz_add(U, U, V);
        if (mpz_odd_p(U)) mpz_add(U, U, n);
        mpz_fdiv_q_2exp(U, U, 1);
                                    /* V:  V2k+1 = (D*U2k + P*V2k)/2 */
        mpz_mul_si(V, V, P);
        mpz_add(V, V, t);
        if (mpz_odd_p(V)) mpz_add(V, V, n);
        mpz_fdiv_q_2exp(V, V, 1);

        mpz_mul_si(Qk, Qk, Q);
      }
      mpz_mod(Qk, Qk, n);
    }
  }
  mpz_mod(U, U, n);
  mpz_mod(V, V, n);
}

static int lucas_selfridge_params(IV* P, IV* Q, mpz_t n, mpz_t t)
{
  IV D = 5;
  UV Dui = (UV) D;
  while (1) {
    UV gcd = mpz_gcd_ui(NULL, n, Dui);
    if ((gcd > 1) && mpz_cmp_ui(n, gcd) != 0)
      return 0;
    mpz_set_si(t, D);
    if (mpz_jacobi(t, n) == -1)
      break;
    if (Dui == 21 && mpz_perfect_square_p(n))
      return 0;
    Dui += 2;
    D = (D > 0)  ?  -Dui  :  Dui;
    if (Dui > 1000000)
      croak("lucas_selfridge_params: D exceeded 1e6");
  }
  if (P) *P = 1;
  if (Q) *Q = (1 - D) / 4;
  return 1;
}

static int lucas_extrastrong_params(IV* P, IV* Q, mpz_t n, mpz_t t, UV inc)
{
  UV tP = 3;
  if (inc < 1 || inc > 256)
    croak("Invalid lucas parameter increment: %lu\n", (unsigned long)inc);
  while (1) {
    UV D = tP*tP - 4;
    UV gcd = mpz_gcd_ui(NULL, n, D);
    if (gcd > 1 && mpz_cmp_ui(n, gcd) != 0)
      return 0;
    mpz_set_ui(t, D);
    if (mpz_jacobi(t, n) == -1)
      break;
    if (tP == (3+20*inc) && mpz_perfect_square_p(n))
      return 0;
    tP += inc;
    if (tP > 65535)
      croak("lucas_extrastrong_params: P exceeded 65535");
  }
  if (P) *P = (IV)tP;
  if (Q) *Q = 1;
  return 1;
}
int is_lucas_pseudoprime(mpz_t n, int strength)
{
  mpz_t d, U, V, Qk, t;
  IV P = 0, Q = 0;
  UV s = 0;
  int rval;

  {
    int cmpr = mpz_cmp_ui(n, 2);
    if (cmpr == 0)     return 1;  /* 2 is prime */
    if (cmpr < 0)      return 0;  /* below 2 is composite */
    if (mpz_even_p(n)) return 0;  /* multiple of 2 is composite */
  }

  mpz_init(t);
  rval = (strength < 2) ? lucas_selfridge_params(&P, &Q, n, t)
                        : lucas_extrastrong_params(&P, &Q, n, t, 1);
  if (!rval) {
    mpz_clear(t);
    return 0;
  }

  mpz_init(U);  mpz_init(V);  mpz_init(Qk);
  mpz_init_set(d, n);
  mpz_add_ui(d, d, 1);

  if (strength > 0) {
    s = mpz_scan1(d, 0);
    mpz_tdiv_q_2exp(d, d, s);
  }

  lucas_seq(U, V, n, P, Q, d, Qk, t);
  mpz_clear(d);

  rval = 0;
  if (strength == 0) {
    /* Standard checks U_{n+1} = 0 mod n. */
    rval = (mpz_sgn(U) == 0);
  } else if (strength == 1) {
    if (mpz_sgn(U) == 0) {
      rval = 1;
    } else {
      while (s--) {
        if (mpz_sgn(V) == 0) {
          rval = 1;
          break;
        }
        if (s) {
          mpz_mul(V, V, V);
          mpz_submul_ui(V, Qk, 2);
          mpz_mod(V, V, n);
          mpz_mulmod(Qk, Qk, Qk, n, t);
        }
      }
    }
  } else {
    mpz_sub_ui(t, n, 2);
    if ( mpz_sgn(U) == 0 && (mpz_cmp_ui(V, 2) == 0 || mpz_cmp(V, t) == 0) ) {
      rval = 1;
    } else {
      s--;  /* The extra strong test tests r < s-1 instead of r < s */
      while (s--) {
        if (mpz_sgn(V) == 0) {
          rval = 1;
          break;
        }
        if (s) {
          mpz_mul(V, V, V);
          mpz_sub_ui(V, V, 2);
          mpz_mod(V, V, n);
        }
      }
    }
  }
  mpz_clear(Qk); mpz_clear(V); mpz_clear(U); mpz_clear(t);
  return rval;
}

int is_prob_prime(mpz_t n)
{
  /*  Step 1: Look for small divisors.  This is done purely for performance.
   *          It is *not* a requirement for the BPSW test. */

  /* If less than 1009, make trial factor handle it. */
  if (mpz_cmp_ui(n, BGCD_NEXTPRIME) < 0)
    return trial_factor(n, 2, BGCD_LASTPRIME) ? 0 : 2;

  /* Check for tiny divisors (GMP can do these really fast) */
  if ( mpz_even_p(n)
    || mpz_divisible_ui_p(n, 3)
    || mpz_divisible_ui_p(n, 5) ) return 0;
  /* Do a big GCD with all primes < 1009 */
  {
    mpz_t t;
    mpz_init(t);
    mpz_gcd(t, n, _bgcd);
    if (mpz_cmp_ui(t, 1) != 0) { mpz_clear(t); return 0; }
    mpz_clear(t);
  }
  /* No divisors under 1009 */
  if (mpz_cmp_ui(n, BGCD_NEXTPRIME*BGCD_NEXTPRIME) < 0)
    return 2;

  /*  Step 2: The BPSW test.  psp base 2 and slpsp. */

  /* Miller Rabin with base 2 */
  if (miller_rabin_ui(n, 2) == 0)
    return 0;

  /* Extra-Strong Lucas test */
  if (is_lucas_pseudoprime(n, 2 /*extra strong*/) == 0)
    return 0;
  /* BPSW is deterministic below 2^64 */
  if (mpz_sizeinbase(n, 2) <= 64)
    return 2;

  return 1;
}

/* These primorial and trial factor functions are really slow for numerous
 * reasons, but most of all because mpz_nextprime is dog slow.  We don't
 * really use them, so don't worry about it too much. */
void GMP_pn_primorial(mpz_t prim, UV n)
{
  mpz_t p;
  mpz_init_set_ui(p, 2);
  mpz_set_ui(prim, 1);
  while (n--) {
    mpz_mul(prim, prim, p);
    mpz_nextprime(p, p);
  }
  mpz_clear(p);
}
UV trial_factor(mpz_t n, UV from_n, UV to_n)
{
  mpz_t p;
  UV f = 0;

  if (mpz_cmp_ui(n, 4) < 0)
    return (mpz_cmp_ui(n, 1) <= 0) ? 1 : 0;   /* 0,1 => 1   2,3 => 0 */
  
  if      (from_n <= 2 && to_n >= 2 && mpz_even_p(n)           ) return 2;
  else if (from_n <= 3 && to_n >= 3 && mpz_divisible_ui_p(n, 3)) return 2;
  if (from_n < 5)
    from_n = 5;
  if (from_n > to_n)
    return 0;

  mpz_init(p);
  mpz_sqrt(p, n);
  if (mpz_cmp_ui(p, to_n) < 0)
    to_n = mpz_get_ui(p);      /* limit to_n to sqrtn */
  mpz_set_ui(p, from_n-1);
  mpz_nextprime(p, p);         /* Set p to the first prime >= from_n */

  while (mpz_cmp_ui(p, to_n) <= 0) {
    if (mpz_divisible_p(n, p)) {
      f = mpz_get_ui(p);
      break;
    }
    mpz_nextprime(p, p);
  }
  mpz_clear(p);
  return f;
}

/*****************************************************************************/
/* Proof verification                                                        */
/*****************************************************************************/

/* What each of these does is verify:
 *      Assume Q is prime.
 *      Then N is prime based on the proof given.
 * We verify any necessary conditions on Q (e.g. it must be odd, or > 0, etc.
 * but do not verify Q prime.  That is done in another proof step.
 */

/* ECPP using N, A, B, M, Q, X, Y
 *
 * A.O.L. Atkin and F. Morain, "Elliptic Curves and primality proving"
 * Mathematics of Computation, v61, 1993, pages 29-68.
 * http://www.ams.org/journals/mcom/1993-61-203/S0025-5718-1993-1199989-X/
 *
 * Page 10, Theorem 5.2:
 *  "Let N be an integer prime to 6, E an elliptic curve over Z/NZ, together
 *   with a point P on E and m and s two integers with s | m.  For each prime
 *   divisor q of s, we put (m/q)P = (x_q : y_q : z_q).  We assume that
 *   mP = O_E and gcd(z_q,N) = 1 for all q.  Then, if p is a prime divisor
 *   of N, one has #E(Z/pZ) = 0 mod s."
 * Page 10, Corollary 5.1:
 *  "With the same conditions, if s > (N^(1/4) + 1)^2, then N is prime."
 *
 * Basically this same result is repeated in Crandall and Pomerance 2005,
 * Theorem 7.6.1 "Goldwasser-Kilian ECPP theorem".
 *
 * Wikipedia, "Elliptic curve primality testing":
 *  "Let N be a positive integer, and E be the set which is defined by the
 *   equation y^2 = x^3 + ax + b (mod N).  Consider E over Z/NZ, use the
 *   usual addition law on E, and write O for the neutral element on E.
 *   Let m be an integer.  If there is a prime q which divides m, and is
 *   greater than (N^(1/4) + 1)^2 and there exists a point P on E such that
 *   (1) mP = O, (2) (m/q)P is defined and not equal to O, then N is prime."
 *
 * We use the restricted form as stated by Wikipedia and used in the
 * Atkin/Morain ECPP algorithm, where s is a prime (hence the "for each prime
 * divisor q of s" of the general theorem is just s).
 */
void verify_ecpp(void) {
  mpz_mod(A, A, N);
  mpz_mod(B, B, N);

  if (mpz_cmp_ui(N, 0) <= 0)   quit_invalid("ECPP", "N > 0");
  if (mpz_gcd_ui(NULL, N, 6) != 1)  quit_invalid("ECPP", "gcd(N, 6) = 1");
  mpz_mul(T1, A, A);
  mpz_mul(T1, T1, A);
  mpz_mul_ui(T1, T1, 4);
  mpz_mul(T2, B, B);
  mpz_mul_ui(T2, T2, 27);
  mpz_add(T1, T1, T2);
  mpz_gcd(T1, T1, N);
  if (mpz_cmp_ui(T1, 1) != 0)  quit_invalid("ECPP", "gcd(4*a^3 + 27*b^2, N) = 1");
  mpz_mul(T1, X, X);
  mpz_add(T1, T1, A);
  mpz_mul(T1, T1, X);
  mpz_add(T1, T1, B);
  mpz_mod(T1, T1, N);
  mpz_mul(T2, Y, Y);
  mpz_mod(T2, T2, N);
  if (mpz_cmp(T1, T2) != 0)   quit_invalid("ECPP", "Y^2 = X^3 + A*X + B mod N");
  mpz_mul_ui(T2, N, 4);
  mpz_sqrt(T2, T2);
  mpz_add_ui(T1, N, 1);
  mpz_sub(T1, T1, T2);
  if (mpz_cmp(M, T1) < 0)     quit_invalid("ECPP", "M >= N + 1 - 2*sqrt(N)");
  mpz_add_ui(T1, N, 1);
  mpz_add(T1, T1, T2);
  if (mpz_cmp(M, T1) > 0)     quit_invalid("ECPP", "M <= N + 1 + 2*sqrt(N)");
  mpz_root(T1, N, 4);
  mpz_add_ui(T1, T1, 1);
  mpz_mul(T1, T1, T1);
  if (mpz_cmp(Q, T1) <= 0)    quit_invalid("ECPP", "Q > (N^(1/4)+1)^2");
  if (mpz_cmp(Q, N) >= 0)     quit_invalid("ECPP", "Q < N");
  /* While M = Q is odd to compute in a proof, it is allowed.
   * In Primo terms, this means S=1 is allowed.
   *   if (mpz_cmp(M, Q) == 0)     quit_invalid("ECPP", "M != Q");
   */
  if (!mpz_divisible_p(M, Q)) quit_invalid("ECPP", "Q divides M");

  {
#if USE_AFFINE_EC
    struct ec_affine_point P0, P1, P2;
    mpz_init_set(P0.x, X);  mpz_init_set(P0.y, Y);
    mpz_init(P1.x); mpz_init(P1.y);
    mpz_init(P2.x); mpz_init(P2.y);
    mpz_divexact(T1, M, Q);
    if (ec_affine_multiply(A, T1, N, P0, &P2, T2))
      quit_invalid("ECPP", "Factor found for N");
    /* Check that P2 is not (0,1) */
    if (mpz_cmp_ui(P2.x, 0) == 0 && mpz_cmp_ui(P2.y, 1) == 0)
      quit_invalid("ECPP", "(M/Q) * EC(A,B,N,X,Y) is not identity");
    mpz_set(T1, Q);
    if (ec_affine_multiply(A, T1, N, P2, &P1, T2))
      quit_invalid("ECPP", "Factor found for N");
    /* Check that P1 is (0, 1) */
    if (! (mpz_cmp_ui(P1.x, 0) == 0 && mpz_cmp_ui(P1.y, 1) == 0) )
      quit_invalid("ECPP", "M * EC(A,B,N,X,Y) is identity");
    mpz_clear(P0.x); mpz_clear(P0.y);
    mpz_clear(P1.x); mpz_clear(P1.y);
    mpz_clear(P2.x); mpz_clear(P2.y);
#else
    mpz_t PX, PY, PA, PB;
    mpz_init(PX);  mpz_init(PY);  mpz_init(PA);  mpz_init(PB);

    /* We have A,B,X,Y in affine coordinates, for the curve:
     *   Y^2 = X^3 + AX + B
     * and want to turn this into points on a Montgomery curve:
     *   by^2 = x^3 + ax^2 + x
     * so we can use the much faster (~4x) multiplication routines.
     * The inverse of this operation is:
     *    X = (3x+a)/3b
     *    Y = y/b
     *    A = (3-a^2)/(3b^2)
     *    B = (2a^3-9a)/27b^3
     * In our case we need to do the harder job of going the other direction.
     */

    /* Make Montgomery variables from affine (TODO: make this work) */
    mpz_add(PB, X, A);
    mpz_mul(PB, PB, X);
    mpz_add_ui(PB, PB, 1);
    mpz_mul(PB, PB, X);
    mpz_mod(PB, PB, N);

    mpz_mul_ui(T2, PB, 3);
    mpz_mul(T2, T2, PB);
    mpz_mod(T2, T2, N);
    mpz_gcdext(T2, T1, NULL, T2, N); /* T1 = 1/3g^2 */
    if (mpz_cmp_ui(T2,1) != 0) quit_invalid("ECPP", "Factor found during gcd");

    mpz_mul_ui(PX, X, 3);
    mpz_add(PX, PX, A);
    mpz_mul(PX, PX, PB);
    mpz_mul(PX, PX, T1);
    mpz_mod(PX, PX, N);

    mpz_set(PY, Y);
    mpz_mul_ui(PY, PY, 3);
    mpz_mul(PY, PY, PB);
    mpz_mul(PY, PY, T1);
    mpz_mod(PY, PY, N);  /* y = (3gY)/(3g^2) = Y/g */

    mpz_mul(PA, A, A);
    mpz_sub_ui(PA, PA, 3);
    mpz_neg(PA, PA);
    mpz_mul(PA, PA, T1);
    mpz_mod(PA, PA, N);

    mpz_divexact(T1, M, Q);
    pec_mult(PA, PB, T1, N, PX, PY);
    /* Check that point is not (0,0) */
    if (mpz_cmp_ui(PX, 0) == 0 && mpz_cmp_ui(PY, 0) == 0)
      quit_invalid("ECPP", "(M/Q) * EC(A,B,N,X,Y) is not identity");
    mpz_set(T1, Q);
    pec_mult(PA, PB, T1, N, PX, PY);
    /* Check that point is (0, 0) */
    if (! (mpz_cmp_ui(PX, 0) == 0 && mpz_cmp_ui(PY, 0) == 0) )
      quit_invalid("ECPP", "M * EC(A,B,N,X,Y) is identity");
    mpz_clear(PX);  mpz_clear(PY);  mpz_clear(PA);  mpz_clear(PB);
#endif
  }
}

/* Basic N+1 using N, Q, LP, LQ
 *
 * John Brillhart, D.H. Lehmer, J.L. Selfridge,
 *   "New Primality Criteria and Factorizations of 2^m +/- 1"
 * Mathematics of Computation, v29, n130, April 1975, pp 620-647.
 * http://www.ams.org/journals/mcom/1975-29-130/S0025-5718-1975-0384673-1/S0025-5718-1975-0384673-1.pdf
 *
 * Page 631, Theorem 15:
 *  "Let N+1 = mq, where q is an odd prime such that 2q-1 > sqrt(N).
 *   If there exists a Lucas sequence {V_k} of discriminant D with
 *   (D|N) = -1 for which N|V_{(N+1)/2}, but N∤V_{m/2}, then N is prime."
 */
void verify_bls15(void) {
  if (mpz_even_p(Q))            quit_invalid("BLS15", "Q odd");
  if (mpz_cmp_ui(Q, 2) <= 0)    quit_invalid("BLS15", "Q > 2");
  mpz_add_ui(T2, N, 1);
  if (!mpz_divisible_p(T2, Q))  quit_invalid("BLS15", "Q divides N+1");
  mpz_divexact(M, T2, Q);
  mpz_mul(T1, M, Q);
  mpz_sub_ui(T1, T1, 1);
  if (mpz_cmp(T1, N) != 0)      quit_invalid("BLS15", "MQ-1 = N");
  if (mpz_cmp_ui(M, 0) <= 0)    quit_invalid("BLS15", "M > 0");
  mpz_mul_ui(T1, Q, 2);
  mpz_sub_ui(T1, T1, 1);
  mpz_sqrt(T2, N);
  if (mpz_cmp(T1, T2) <= 0)     quit_invalid("BLS15", "2Q-1 > sqrt(N)");
  mpz_mul(T1, LP, LP);
  mpz_mul_ui(T2, LQ, 4);
  mpz_sub(T1, T1, T2);
  if (mpz_sgn(T1) == 0)         quit_invalid("BLS15", "D != 0");
  if (mpz_jacobi(T1, N) != -1)  quit_invalid("BLS15", "jacobi(D,N) = -1");
  {
    mpz_t U, V, k;
    IV iLP, iLQ;
    mpz_init(U);  mpz_init(V);  mpz_init(k);

    iLP = mpz_get_si(LP);
    iLQ = mpz_get_si(LQ);
    if (mpz_cmp_si(LP, iLP) != 0)  quit_error("BLS15 LP out of range", "");
    if (mpz_cmp_si(LQ, iLQ) != 0)  quit_error("BLS15 LQ out of range", "");

    mpz_tdiv_q_2exp(k, M, 1);
    lucas_seq(U, V, N, iLP, iLQ, k, T1, T2);
    if (mpz_sgn(V) == 0)        quit_invalid("BLS15", "V_{m/2} mod N != 0");
    mpz_add_ui(k, N, 1);
    mpz_tdiv_q_2exp(k, k, 1);
    lucas_seq(U, V, N, iLP, iLQ, k, T1, T2);
    if (mpz_sgn(V) != 0)        quit_invalid("BLS15", "V_{(N+1)/2} mod N == 0");

    mpz_clear(U);  mpz_clear(V);  mpz_clear(k);
  }
}

/* Simplistic N-1 using N, Q, A
 *
 * Hans Riesel, "Prime Numbers and Computer Methods for Factorization"
 * page 103-104, Theorem 4.6:
 *  "Suppose N-1 = R*F = R prod(j=1,n,q_j^{B_j}), with all q_j's distinct
 *   primes, with GCD(R,F) = 1 and R < F.  If an integer a can be found, s.t.
 *     GCD(A^((N-1)/q_j)-1,N) = 1 for all j=1..n
 *   and satisfying
 *     a^(N-1) = 1 mod N
 *   then N is a prime."
 *
 * Now make the severe restriction that F must be a single prime q.  This then
 * reduces to the Wikipedia "Pocklington primality test":
 *  "Let N > 1 be an integer, and suppose there exist numbers a and q such that
 *   (1) q is prime, q|N-1 and q > sqrt(N)-1
 *   (2) a^(N-1) = 1 mod N
 *   (3) gcd(a^((N-1)/q)-1,N) = 1
 *   Then N is prime."
 *
 * Note that BLS75 theorem 3 is similar, but also allows a smaller q.  Also,
 * BLS75 theorem 5 is a much better method than generalized Pocklington.
 */
void verify_pocklington(void)
{
  mpz_sub_ui(T2, N, 1);
  if (!mpz_divisible_p(T2, Q))  quit_invalid("Pocklington", "Q divides N-1");
  mpz_divexact(M, T2, Q);
  if (mpz_odd_p(M))             quit_invalid("Pocklington", "M is even");
  if (mpz_cmp_ui(M, 0) <= 0)    quit_invalid("Pocklington", "M > 0");
  if (mpz_cmp(M, Q) >= 0)       quit_invalid("Pocklington", "M < Q");
  mpz_mul(T1, M, Q);
  mpz_add_ui(T1, T1, 1);
  if (mpz_cmp(T1, N) != 0)      quit_invalid("Pocklington", "MQ+1 = N");
  if (mpz_cmp_ui(A, 1) <= 0)    quit_invalid("Pocklington", "A > 1");
  if (mpz_cmp(A, N) >= 0)       quit_invalid("Pocklington", "A < N");
  mpz_powm(T1, A, T2, N);
  if (mpz_cmp_ui(T1, 1) != 0)   quit_invalid("Pocklington", "A^(N-1) mod N = 1");
  mpz_powm(T1, A, M, N);
  if (mpz_sgn(T1)) mpz_sub_ui(T1, T1, 1);
  else             mpz_set(T1, T2);
  mpz_gcd(T1, T1, N);
  if (mpz_cmp_ui(T1, 1) != 0)   quit_invalid("Pocklington", "gcd(A^M - 1, N) = 1");
}

/* Basic N-1 using N, Q, A
 *
 * John Brillhart, D.H. Lehmer, J.L. Selfridge,
 *   "New Primality Criteria and Factorizations of 2^m +/- 1"
 * Mathematics of Computation, v29, n130, April 1975, pp 620-647.
 * http://www.ams.org/journals/mcom/1975-29-130/S0025-5718-1975-0384673-1/S0025-5718-1975-0384673-1.pdf
 *
 * Page 622-623, Theorem 3:
 *  "Let N-1 = mp, where p is an odd prime such that 2p+1 > sqrt(N).
 *   If there exists an a for which a^((N-1)/2) = -1 mod N,
 *   but a^(m/2) != -1 mod N, then N is prime."
 */
void verify_bls3(void) {
  if (mpz_even_p(Q))            quit_invalid("BLS3", "Q odd");
  if (mpz_cmp_ui(Q, 2) <= 0)    quit_invalid("BLS3", "Q > 2");
  mpz_sub_ui(T2, N, 1);
  if (!mpz_divisible_p(T2, Q))  quit_invalid("BLS3", "Q divides N-1");
  mpz_divexact(M, T2, Q);
  mpz_mul(T1, M, Q);
  mpz_add_ui(T1, T1, 1);
  if (mpz_cmp(T1, N) != 0)      quit_invalid("BLS3", "MQ+1 = N");
  if (mpz_cmp_ui(M, 0) <= 0)    quit_invalid("BLS3", "M > 0");
  mpz_mul_ui(T1, Q, 2);
  mpz_add_ui(T1, T1, 1);
  mpz_sqrt(T2, N);
  if (mpz_cmp(T1, T2) <= 0)     quit_invalid("BLS3", "2Q+1 > sqrt(N)");
  mpz_sub_ui(T2, N, 1);
  mpz_divexact_ui(T1, T2, 2);
  mpz_powm(T1, A, T1, N);
  if (mpz_cmp(T1, T2) != 0)     quit_invalid("BLS3", "A^((N-1)/2) = N-1 mod N");
  mpz_divexact_ui(T1, M, 2);
  mpz_powm(T1, A, T1, N);
  if (mpz_cmp(T1, T2) == 0)     quit_invalid("BLS3", "A^(M/2) != N-1 mod N");
}

/* Sophisticated N-1 using N, QARRAY, AARRAY
 *
 * John Brillhart, D.H. Lehmer, J.L. Selfridge,
 *   "New Primality Criteria and Factorizations of 2^m +/- 1"
 * Mathematics of Computation, v29, n130, April 1975, pp 620-647.
 * http://www.ams.org/journals/mcom/1975-29-130/S0025-5718-1975-0384673-1/S0025-5718-1975-0384673-1.pdf
 *
 * Page 621-622: "The expression 'N is a psp base a' will be used for a
 * number N which satisfies the congruence a^(N-1) = 1 mod N, 1 < a < N-1,
 * i.e., N is a "pseudoprime" base a."
 * Page 623: "Throughout the rest of this paper the notation N-1 = F_1 R_1
 * will be used, where F_1 is the even factored portion of N-1, R_1 is > 1,
 * and (F_1,R_1) = 1."
 * Page 623: "(I) For each prime p_i dividing F_1 there exists an a_i such
 * that N is a psp base a_i and (a_i^((N-1)/p_i),N) = 1."
 * Page 624, Theorem 5.
 *  "Assume (I) and let m be >= 1.
 *   When m > 1, assume further that λF_1 + 1 ∤ N for 1 <= λ < m.
 *   If N < (mF_1 + 1) [2(F_1)^2 + (r-m)F_1 + 1],
 *   where r and s are defined by R_1 = (N-1)/F_1 = 2(F_1)s + r, 1 <= r < 2F_1,
 *   then N is prime if and only if s = 0 or r^2 - 8s is not a perfect square.
 *   r != 0 since R_1 is odd."
 *
 * Note that we are using m=1, which is simple for the verifier.
 *
 */
void verify_bls5(int num_qs) {
  int i;
  mpz_t F, R, s, r;

  if (mpz_cmp_ui(N, 2) <= 0)   quit_invalid("BLS5", "N > 2");
  if (mpz_even_p(N))           quit_invalid("BLS5", "N odd");
  mpz_sub_ui(T2, N, 1);
  mpz_init_set_ui(F, 1);
  mpz_init_set(R, T2);
  mpz_init(s);  mpz_init(r);
  for (i = 0; i < num_qs; i++) {
    if (mpz_cmp_ui(QARRAY[i], 1 ) <= 0)  quit_invalid("BLS5", "Q > 1");
    if (mpz_cmp(   QARRAY[i], T2) >= 0)  quit_invalid("BLS5", "Q < N-1");
    if (mpz_cmp_ui(AARRAY[i], 1 ) <= 0)  quit_invalid("BLS5", "A > 1");
    if (mpz_cmp(   AARRAY[i], T2) >= 0)  quit_invalid("BLS5", "A < N-1");
    if (!mpz_divisible_p(T2, QARRAY[i])) quit_invalid("BLS5", "Q divides N-1");
    while (mpz_divisible_p(R, QARRAY[i])) {
      mpz_mul(F, F, QARRAY[i]);
      mpz_divexact(R, R, QARRAY[i]);
    }
  }
  mpz_mul(T1, R, F);
  if (mpz_cmp(T1, T2) != 0)    quit_invalid("BLS5", "R == (N-1)/F");
  if (mpz_odd_p(F))            quit_invalid("BLS5", "F is even");
  mpz_gcd(T1, F, R);
  if (mpz_cmp_ui(T1, 1) != 0)  quit_invalid("BLS5", "gcd(F, R) = 1");
  mpz_mul_ui(T1, F, 2);
  mpz_tdiv_qr(s, r, R, T1);

  mpz_mul_ui(T1, F, 2);   /* T1 = 2*F */
  mpz_sub_ui(T2, r, 1);
  mpz_add(T1, T1, T2);    /* T1 = 2*F + (r-1) */
  mpz_mul(T1, T1, F);     /* T1 = 2*F*F + (r-1)*F */
  mpz_add_ui(T1, T1, 1);  /* T1 = 2*F*F + (r-1)*F + 1 */
  mpz_add_ui(T2, F, 1);
  mpz_mul(T1, T1, T2);    /* T1 = (F+1) * (2*F*F + (r-1)*F + 1) */
  if (mpz_cmp(N, T1) >= 0) quit_invalid("BLS5", "N < (F+1)(2F^2+(r-1)F+1)");
  if (mpz_sgn(s) != 0) {
    mpz_mul(T2, r, r);
    mpz_submul_ui(T2, s, 8);  /* T2 = r*r - 8*s */
    if (mpz_perfect_square_p(T2))  quit_invalid("BLS5", "S=0 OR R^2-8S not a perfect square");
  }
  mpz_clear(F); mpz_clear(R); mpz_clear(s); mpz_clear(r);
  mpz_sub_ui(T2, N, 1);
  for (i = 0; i < num_qs; i++) {
    mpz_powm(T1, AARRAY[i], T2, N);
    if (mpz_cmp_ui(T1, 1) != 0)   quit_invalid("BLS5", "A[i]^(N-1) mod N = 1");
    mpz_divexact(T1, T2, QARRAY[i]);
    mpz_powm(T1, AARRAY[i], T1, N);
    if (mpz_sgn(T1)) mpz_sub_ui(T1, T1, 1);
    else             mpz_set(T1, T2);
    mpz_gcd(T1, T1, N);
    if (mpz_cmp_ui(T1, 1) != 0)   quit_invalid("BLS5", "gcd(A[i]^((N-1)/Q[i]) - 1, N) = 1");
  }
}

/* Most basic N-1 using N, QARRAY, A
 *
 * D.H. Lehmer, "Tests for Primality by the Converse of Fermat's Theorem"
 * Bull. AMS, v33, n3, 1927, pp 327-340.
 * http://projecteuclid.org/DPubS?service=UI&version=1.0&verb=Display&handle=euclid.bams/1183492108
 *
 * Page 330, Theorem 2:
 *  "If a^x = 1 mod N for x = N-1, but not for x a quotient of N-1 on
 *   division by any of its prime factors, then N is a prime."
 */
void verify_lucas(int num_qs) {
  int i;
  mpz_sub_ui(T2, N, 1);
  mpz_set(R, T2);
  if (mpz_cmp_ui(A, 1) <= 0)  quit_invalid("Lucas", "A > 1");
  if (mpz_cmp(   A, N) >= 0)  quit_invalid("Lucas", "A < N");

  mpz_powm(T1, A, T2, N);
  if (mpz_cmp_ui(T1, 1) != 0)   quit_invalid("Lucas", "A^(N-1) mod N = 1");

  for (i = 1; i < num_qs; i++) {
    if (mpz_cmp_ui(QARRAY[i], 1 ) <= 0)  quit_invalid("Lucas", "Q > 1");
    if (mpz_cmp(   QARRAY[i], T2) >= 0)  quit_invalid("Lucas", "Q < N-1");
    if (!mpz_divisible_p(T2, QARRAY[i])) quit_invalid("Lucas", "Q divides N-1");

    mpz_divexact(T1, T2, QARRAY[i]);
    mpz_powm(T1, A, T1, N);
    if (mpz_cmp_ui(T1, 1) == 0)   quit_invalid("Lucas", "A^((N-1)/Q[i]) mod N != 1");
    while (mpz_divisible_p(R, QARRAY[i]))
      mpz_divexact(R, R, QARRAY[i]);
  }
  if (mpz_cmp_ui(R, 1) != 0)
    quit_invalid("Lucas", "N-1 has only factors Q[i]");
}

void verify_primo_ecpp(void) {
  /* Calculate a,b,x,y from A, B, T, then verify */
  mpz_mul(T1, T, T);
  mpz_add(T1, T1, A);
  mpz_mul(T1, T1, T);
  mpz_add(T1, T1, B);
  mpz_mod(T1, T1, N);   /* L = T1 = T^3 + AT + B mod N */
  if (mpz_sgn(T1) <= 0) quit_invalid("Primo ECPP", "L > 0");

  mpz_mul(T2, T1, T1);
  mpz_mul(A, A, T2);
  mpz_mod(A, A, N);     /* a = AL^2 mod N */

  mpz_mul(T2, T2, T1);
  mpz_mul(B, B, T2);
  mpz_mod(B, B, N);     /* b = BL^3 mod N */

  mpz_mul(X, T, T1);
  mpz_mod(X, X, N);     /* x = TL mod N */
  mpz_mul(Y, T1, T1);
  mpz_mod(Y, Y, N);     /* y = L^2 mod N */

  mpz_mul(M, R, S);     /* M = R*S */
  mpz_set(Q, R);        /* Q = R */

  verify_ecpp();  /* N, A, B, M, Q, X, Y */
}

void verify_ecpp4(void) {
  mpz_mul_ui(T1, J, 2);
  if (mpz_cmpabs(T1, N) > 0) quit_invalid("Primo Type 4", "|J| <= N/2");
  if (mpz_cmp_ui(T, 0) < 0)  quit_invalid("Primo Type 4", "T >= 0");
  if (mpz_cmp(T, N) >= 0)    quit_invalid("Primo Type 4", "T < N");

  mpz_set_ui(T2, 1728);
  mpz_sub(T2, T2, J);
  mpz_mul(A, T2, J);
  mpz_mul_ui(A, A, 3);    /* A = 3 J (1728-J)   */

  mpz_mul(T2, T2, T2);
  mpz_mul(B, T2, J);
  mpz_mul_ui(B, B, 2);    /* B = 2 J (1728-J)^2 */

  verify_primo_ecpp();  /* (A,B,T)->(a,b,x,y)  (R,S)->(M,Q)  Verify */
}

void verify_ecpp3(void) {
  mpz_mul_ui(T1, A, 2);
  mpz_mul_ui(T2, B, 2);
  if (mpz_cmpabs(T1, N) > 0) quit_invalid("Primo Type 3", "|A| <= N/2");
  if (mpz_cmpabs(T2, N) > 0) quit_invalid("Primo Type 3", "|B| <= N/2");
  if (mpz_cmp_ui(T, 0) < 0)  quit_invalid("Primo Type 3", "T >= 0");
  if (mpz_cmp(T, N) >= 0)    quit_invalid("Primo Type 3", "T < N");

  verify_primo_ecpp();  /* (A,B,T)->(a,b,x,y)  (R,S)->(M,Q)  Verify */
}

void verify_primo2(void) {
  mpz_set(LQ, Q);
  mpz_set_ui(LP, mpz_odd_p(LQ) ? 2 : 1);
  mpz_set(Q, R);
  verify_bls15();  /* N, Q, LP, LQ */
}

void verify_primo1(void) {
  mpz_set(Q, R);
  mpz_set(A, B);
  mpz_mul(T1, S, R);
  mpz_add_ui(T1, T1, 1);
  if (mpz_cmp(T1, N) != 0)  quit_invalid("Primo Type 1", "SR+1 = N");
  verify_pocklington();  /* N, Q, A */
}

void verify_ecpp4_42(void) {
  if (mpz_cmp_ui(S, 0) <= 0)  quit_invalid("Primo Type 4", "S > 0");
  mpz_mul(T1, W, W);
  mpz_mul_ui(T2, N, 4);
  if (mpz_cmp(T1, T2) >= 0)   quit_invalid("Primo Type 4", "W^2 < 4*N");
  mpz_add_ui(R, N, 1);
  mpz_sub(R, R, W);
  if (!mpz_divisible_p(R, S)) quit_invalid("Primo Type 4", "(N+1-W) mod S = 0");
  mpz_divexact(R, R, S);

  verify_ecpp4();
}
void verify_ecpp3_42(void) {
  if (mpz_cmp_ui(S, 0) <= 0)  quit_invalid("Primo Type 4", "S > 0");
  mpz_mul(T1, W, W);
  mpz_mul_ui(T2, N, 4);
  if (mpz_cmp(T1, T2) >= 0)   quit_invalid("Primo Type 4", "W^2 < 4*N");
  mpz_add_ui(R, N, 1);
  mpz_sub(R, R, W);
  if (!mpz_divisible_p(R, S)) quit_invalid("Primo Type 4", "(N+1-W) mod S = 0");
  mpz_divexact(R, R, S);

  verify_ecpp3();
}
void verify_primo2_42(void) {
  mpz_add_ui(R, N, 1);
  if (!mpz_even_p(S))          quit_invalid("Primo N+1", "S is even");
  if (mpz_cmp_ui(S,1) <= 0)    quit_invalid("Primo N+1", "S > 1");
  if (!mpz_divisible_p(R, S))  quit_invalid("Primo N+1", "(N+1) mod S = 0");
  if (mpz_cmp_ui(Q,0) <= 0)    quit_invalid("Primo N+1", "Q > 0");
  if (mpz_cmp(Q,N) >= 0)       quit_invalid("Primo N+1", "Q < N");
  mpz_divexact(R, R, S);
  if (mpz_jacobi(Q, N) != -1)  quit_invalid("Primo N+1", "jacobi(Q,N) = -1");

  mpz_set(LQ, Q);
  mpz_set_ui(LP, mpz_odd_p(LQ) ? 2 : 1);
  mpz_set(Q, R);
  verify_bls15();  /* N, Q, LP, LQ */
}

void verify_primo1_42(void) {
  mpz_sub_ui(R, N, 1);
  if (!mpz_even_p(S))          quit_invalid("Primo N-1", "S is even");
  if (mpz_cmp_ui(S,1) <= 0)    quit_invalid("Primo N-1", "S > 1");
  if (!mpz_divisible_p(R, S))  quit_invalid("Primo N-1", "(N-1) mod S = 0");
  if (mpz_cmp_ui(B,1) <= 0)    quit_invalid("Primo N-1", "B > 1");
  if (mpz_cmp(B,N) >= 0)       quit_invalid("Primo N-1", "B < N");
  mpz_divexact(R, R, S);

  mpz_set(Q, R);
  mpz_set(A, B);
  verify_pocklington();  /* N, Q, A */
}



void verify_small(void) {
  if (mpz_sizeinbase(N, 2) > 64)  quit_invalid("Small", "N <= 2^64");
  if (is_prob_prime(N) != 2)      quit_invalid("Small", "N does not pass BPSW");
}

void add_chain(mpz_t n, mpz_t q) {
  mpz_init_set(_chain_n[_num_chains], n);
  mpz_init_set(_chain_q[_num_chains], q);
  _num_chains++;
}
void free_chains(void) {
  while (_num_chains > 0) {
    _num_chains--;
    mpz_clear(_chain_n[_num_chains]);
    mpz_clear(_chain_q[_num_chains]);
  }
}
void verify_chain(mpz_t n) {
  int i, found;
  if (mpz_sizeinbase(n, 2) <= 64) {
    mpz_set(N, n);
    verify_small();
    return;
  }
  found = 0;
  for (i = 0; i < _num_chains; i++) {
    if (mpz_cmp(n, _chain_n[i]) == 0) {
      found = 1;
      verify_chain(_chain_q[i]);
    }
  }
  mpz_set(N, n);
  if (!found)  quit_invalid("Final", "q value has no proof");
}


void verify_final(void) {
  if (_format == CERT_PRIMO) {
    if (mpz_cmp_ui(N, 18) <= 0)  quit_invalid("Primo Type 0", "N > 18");
    if (mpz_cmp_ui(N, 340000000000000UL) >= 0)  quit_invalid("Primo Type 0", "N < 34 * 10^13");
    if (!miller_rabin_ui(N,  2) || !miller_rabin_ui(N,  3) ||
        !miller_rabin_ui(N,  5) || !miller_rabin_ui(N,  7) ||
        !miller_rabin_ui(N, 11) || !miller_rabin_ui(N, 13) ||
        !miller_rabin_ui(N, 17))
      quit_invalid("Primo Type 0", "N is SPSP(2,3,5,7,11,13,17)");
  } else if (_format == CERT_PRIMO42) {
    verify_small();
  } else {
    verify_chain(PROOFN);
    free_chains();
  }
}

/*****************************************************************************/
/* File parsing                                                              */
/*****************************************************************************/

static int get_line(int signaleof) {
  size_t i;
  /* Read in a line */
  if (fgets(_line, MAX_LINE_LEN, _fh) != _line) {
    if (signaleof) return 1;
    else           quit_error("Error reading from file: ", _filename);
  }
  _line[MAX_LINE_LEN] = '\0';
  /* Remove trailing newlines and spaces */
  i = strlen(_line);
  while ( i > 0 && isspace(_line[i-1]) )
    _line[--i] = '\0';
  return 0;
}

#define PROCESS_VAR(v) \
  do { \
    mpz_set_str(v, _vstr, _base); \
    for (i = 0; i < nargs; i++) { \
      if (vlist[i] != 0 && strcmp(vlist[i], #v) == 0) { \
        vfound++; \
        vlist[i] = 0; \
        break; \
      } \
    } \
    if (i >= nargs) \
      quit_error("Unknown variable: ", #v); \
  } while (0)

void read_vars(const char* vars) {
  char* varstring = strdup(vars);
  char* vlist[10];
  char  varname;
  int i;
  int nargs = 0;
  int vfound = 0;
  int bad_lines = 0;

  vlist[0] = strtok(varstring, " ");
  while (vlist[nargs] != 0)
    vlist[++nargs] = strtok(NULL, " ");
  while (vfound < nargs) {
    get_line(0);
    if (strlen(_line) == 0)    /* Skip extrenuous blank lines */
      continue;
    if (_format == CERT_PRIMO) {
      int varnum = 0;  /* Did we read a variable properly this line? */
      if (sscanf(_line, "%c$=%s", &varname, _vstr) == 2) {
        for (i = 0; i < nargs && varnum == 0; i++) {
          if (vlist[i] != 0 && varname == vlist[i][0]) {
            switch (varname) {
              case 'S':  mpz_set_str(S, _vstr, 16);  break;
              case 'R':  mpz_set_str(R, _vstr, 16);  break;
              case 'A':  mpz_set_str(A, _vstr, 16);  break;
              case 'B':  mpz_set_str(B, _vstr, 16);  break;
              case 'Q':  mpz_set_str(Q, _vstr, 16);  break;
              case 'T':  mpz_set_str(T, _vstr, 16);  break;
              case 'J':  mpz_set_str(J, _vstr, 16);  break;
              default: quit_error("Internal error: bad Primo variable type","");
                       break;
            }
            varnum = i+1;
          }
        }
      } else if (sscanf(_line, "[%d]", &i) == 1) {
        quit_error("Variables missing from proof step", "");
      }
      if (varnum != 0) {
        vfound++;             /* We found a variable on the list */
        vlist[varnum-1] = 0;  /* It should only appear once */
        bad_lines = 0;
      } else {
        if (_verbose) { printf("%60s\r", ""); printf("skipping bad line: %s\n", _line); }
        if (bad_lines++ >= BAD_LINES_ALLOWED)
          quit_error("Too many bad lines reading variables", "");
      }
    } else {
      if      (sscanf(_line, "N %s", _vstr) == 1) PROCESS_VAR(N);
      else if (sscanf(_line, "A %s", _vstr) == 1) PROCESS_VAR(A);
      else if (sscanf(_line, "B %s", _vstr) == 1) PROCESS_VAR(B);
      else if (sscanf(_line, "M %s", _vstr) == 1) PROCESS_VAR(M);
      else if (sscanf(_line, "Q %s", _vstr) == 1) PROCESS_VAR(Q);
      else if (sscanf(_line, "X %s", _vstr) == 1) PROCESS_VAR(X);
      else if (sscanf(_line, "Y %s", _vstr) == 1) PROCESS_VAR(Y);
      else if (sscanf(_line, "LQ %s", _vstr) == 1) PROCESS_VAR(LQ);
      else if (sscanf(_line, "LP %s", _vstr) == 1) PROCESS_VAR(LP);
      /* ECPP3 and ECPP4 */
      else if (sscanf(_line, "T %s", _vstr) == 1) PROCESS_VAR(T);
      else if (sscanf(_line, "J %s", _vstr) == 1) PROCESS_VAR(J);
      else if (sscanf(_line, "S %s", _vstr) == 1) PROCESS_VAR(S);
      else if (sscanf(_line, "R %s", _vstr) == 1) PROCESS_VAR(R);
      else
        quit_error("Internal error: bad MPU variable type", "");
    }
  }
  free(varstring);
}

/* Primo version 4.2 doesn't give us a step type, but infers it from which
 * variables appear.  We have to use a different processing method. */
int read42_vars(void) {
  char  varname;
  int i;
  int bad_lines = 0;
  int varsign, fS = 0, fW = 0, fT = 0, fJ = 0, fA = 0, fB = 0, fQ = 0;

  MPUassert( _format == CERT_PRIMO42, "wrong parse routine for cert type");
  while (1) {
    varsign = 0;
    get_line(0);
    if (strlen(_line) == 0)    /* Skip extrenuous blank lines */
      continue;

    if (sscanf(_line, "[%d]", &i) == 1)
      quit_error("Variables missing from proof step", "");

    /* TODO: Hacky even for this file.  More robust method would be nice. */
    if      (sscanf(_line, "%c=-$%s",  &varname, _vstr) == 2)  varsign =-1;
    else if (sscanf(_line, "%c=$%s",   &varname, _vstr) == 2)  varsign = 1;
    else if (sscanf(_line, "%c=-0x%s", &varname, _vstr) == 2)  varsign =-1;
    else if (sscanf(_line, "%c=0x%s",  &varname, _vstr) == 2)  varsign = 1;
    else if (sscanf(_line, "%c=%s",    &varname, _vstr) == 2)  varsign = 2;

    if (varsign != 0) { /* Process variable found */
      mpz_set_str(T1, _vstr, (varsign == 2) ? 10 : 16);
      if (varsign < 0) mpz_neg(T1, T1);
      switch (varname) {
        case 'S': mpz_set(S,T1);  fS=1;  break;
        case 'W': mpz_set(W,T1);  fW=1;  break;
        case 'T': mpz_set(T,T1);  fT=1;  break;
        case 'J': mpz_set(J,T1);  fJ=1;  break;
        case 'A': mpz_set(A,T1);  fA=1;  break;
        case 'B': mpz_set(B,T1);  fB=1;  break;
        case 'Q': mpz_set(Q,T1);  fQ=1;  break;
        default: quit_error("Internal error: bad Primo variable type","");
                 break;
      }
      bad_lines = 0;
    } else {
      if (_verbose) { printf("%60s\r", ""); printf("skipping bad line: %s\n", _line); }
      if (bad_lines++ >= BAD_LINES_ALLOWED)
        quit_error("Too many bad lines reading variables", "");
    }
    if (fS && fW && fT && fJ) break;
    if (fS && fW && fT && fA && fB) break;
    if (fS && !fW && fQ) break;
    if (fS && !fW && fB) break;
  }
  if (!fS) quit_error("Invalid file: no S", "");
  if (fW) {
    if (!fT) quit_error("Invalid file: W no T", "");
    if (fJ)
      return 4;     /* type 4: ECPP */
    if (!fA) quit_error("Invalid file: W no A", "");
    if (!fB) quit_error("Invalid file: W no B", "");
      return 3;     /* type 3: ECPP */
  } else if (fQ) {
    return 2;       /* type 2: n+1 */
  } else if (fB) {
    return 1;       /* type 2: n-1 */
  }
  quit_error("Invalid file: no {W,Q,B}", "");
  return 0;
}

/* TODO:
 * rearrange so we (1) open and read everything up to proof for / candidate.
 * then (2) func that does the proof
 * this should let us call #2 if we hit another proof, so we can do multiple
 * proof in a file.
 */

void parse_top(void)
{
  do {
    if (get_line(1))
      quit_error("Count not find primality certificate indicator", "");
  } while (strstr(_line, "Primality Certificate") == 0);
  if (strcmp(_line, "[PRIMO - Primality Certificate]") == 0)
    _format = CERT_PRIMO;
  else if (strcmp(_line, "[MPU - Primality Certificate]") == 0)
    _format = CERT_MPU;
  else
    quit_error("First line in file is not primality certificate indicator","");

  if (_format == CERT_PRIMO) {
    int items_found = 0;
    int format_ver;
    while (items_found < 3) {
      get_line(0);
      if (sscanf(_line, "Format=%d", &format_ver) == 1) {
        switch (format_ver) {
          case 3: break;
          case 4: _format = CERT_PRIMO42; break;
          default: quit_error("Unknown version number", "");
        }
      }
      if (sscanf(_line, "TestCount=%d", &_testcount) == 1) items_found++;
      if (strcmp(_line, "[Candidate]") == 0)  items_found++;
      if (_format == CERT_PRIMO42) {
        if (sscanf(_line, "N=$%s", _vstr) == 1) items_found++;
        else if (sscanf(_line, "N=0x%s", _vstr) == 1) items_found++;
      } else {
        if (sscanf(_line, "N$=%s", _vstr) == 1) items_found++;
      }
    }
    mpz_set_str(PROOFN, _vstr, 16);
  } else {
    while (1) {
      get_line(0);
      if (_line[0] == '#')  continue;
      if (sscanf(_line, "Base %d", &_base) == 1) continue;
      if (strcmp(_line, "Proof for:") == 0) {
        read_vars("N");
        mpz_set(PROOFN, N);
        break;
      }
    }
  }
  if (!_quiet) {
    printf("%60s\r", "");
    printf("N: ");
    if (_verbose) gmp_printf("%Zd ", PROOFN);
    printf("(%d digits)\n", (int)mpz_sizeinbase(PROOFN, 10));
    printf("Verifying probable prime status.\r");
    fflush(stdout);
  }
  if (is_prob_prime(PROOFN) == 0)
    quit_composite();
  mpz_set(N, PROOFN);
}

void process_file(const char* filename)
{
  _step = 0;
  if (strcmp(filename, "-") == 0)
    _fh = stdin;
  else if ((_fh = fopen(filename, "r")) == NULL)
    quit_error("Unable to open file: ", filename);
  _filename = filename;

  parse_top();

  if (_format == CERT_PRIMO) {
    int type;
    while (1) {
      int rstep;
      get_line(0);
      if (sscanf(_line, "[%d]", &rstep) == 1) {
        if (rstep != _step+1)
          quit_error("Wrong step number found", "");
        _step++;
      }
      if (sscanf(_line, "Type=%d", &type) == 1) {
        if (!_quiet) { printf("%60s\r", ""); printf("Step %3d/%-3d %5d digits  Type %d\r", _step, _testcount, (int)mpz_sizeinbase(N,10), type); fflush(stdout); }
        switch (type) {
          case 4:  read_vars("S R J T");   verify_ecpp4();  break;
          case 3:  read_vars("S R A B T"); verify_ecpp3();  break;
          case 2:  read_vars("S R Q");     verify_primo2(); break;
          case 1:  read_vars("S R B");     verify_primo1(); break;
          case 0:  /* verify_small */ break;
          default: quit_error("Parsing", "Unknown type");   break;
        }
        if (type == 0) break;
        mpz_set(N, R);
      }
    }
  } else if (_format == CERT_PRIMO42) {
    int type, rstep, rvars = 0;
    while (_step < _testcount) {
      get_line(0);
      if (sscanf(_line, "[%d]", &rstep) == 1) {
        if (rstep != _step+1)
          quit_error("Wrong step number found", "");
        _step++;
        rvars = 1;
      }
      if (!rvars) continue;  /* skip until we find a section */
      type = read42_vars();
      if (!_quiet) { printf("%60s\r", ""); printf("Step %3d/%-3d %5d digits  Type %d\r", _step, _testcount, (int)mpz_sizeinbase(N,10), type); fflush(stdout); }
      switch (type) {
        case 4:  verify_ecpp4_42();  break;
        case 3:  verify_ecpp3_42();  break;
        case 2:  verify_primo2_42(); break;
        case 1:  verify_primo1_42(); break;
        case 0:  /* verify_small */ break;
        default: quit_error("Parsing", "Unknown type");   break;
      }
      mpz_set(N, R);
      rvars = 0;
      if (type == 0) break;
    }
  } else {
    char type[MAX_LINE_LEN+1];
    while (1) {
      if (get_line(1)) break;
      if (sscanf(_line, "Type %s", type) == 1) {
        { /* Convert type to upper case */
          char* s = type;
          while (*s != '\0') {
            if (islower(*s))  *s = toupper(*s);
            s++;
          }
        }
        _step++;
        /* TODO: Quick parse of the file to count steps? */
        if (!_quiet) { printf("%60s\r", ""); printf("Step %-4d %5d digits  Type %s\r", _step, (int)mpz_sizeinbase(N,10), type); fflush(stdout); }
        if        (strcmp(type, "ECPP" ) == 0) { read_vars("N A B M Q X Y");
                                                 verify_ecpp();
                                                 add_chain(N, Q);
        } else if (strcmp(type, "ECPP3") == 0) { read_vars("N S R A B T");
                                                 verify_ecpp3();
                                                 add_chain(N, Q);
        } else if (strcmp(type, "ECPP4") == 0) { read_vars("N S R J T");
                                                 verify_ecpp4();
                                                 add_chain(N, Q);
        } else if (strcmp(type, "BLS15") == 0) { read_vars("N Q LP LQ");
                                                 verify_bls15();
                                                 add_chain(N, Q);
        } else if (strcmp(type, "BLS3" ) == 0) { read_vars("N Q A");
                                                 verify_bls3();
                                                 add_chain(N, Q);
        } else if (strcmp(type, "POCKLINGTON") == 0) { read_vars("N Q A");
                                                 verify_pocklington();
                                                 add_chain(N, Q);
        } else if (strcmp(type, "SMALL") == 0) { read_vars("N");
                                                 verify_small();
        } else if (strcmp(type, "BLS5") == 0) {
          int i, index;
          for (index = 0; index < MAX_QARRAY; index++) {
            mpz_set_ui(QARRAY[index], 0);
            mpz_set_ui(AARRAY[index], 2);
          }
          mpz_set_ui(QARRAY[0], 2);
          index = 1;
          while (1) {
            get_line(0);
            if (_line[0] == '-') {
              break;
            } else if (sscanf(_line, "N %s", _vstr) == 1) {
              mpz_set_str(N, _vstr, _base);
            } else if (sscanf(_line, "Q[%d] %s", &i, _vstr) == 2) {
              if (i != index) quit_error("BLS5", "Invalid Q index");
              mpz_set_str(QARRAY[i], _vstr, _base);
              index++;
            } else if (sscanf(_line, "A[%d] %s", &i, _vstr) == 2) {
              if (i < 0 || i > index) quit_error("BLS5", "Invalid A index");
              mpz_set_str(AARRAY[i], _vstr, _base);
            }
          }
          verify_bls5(index);
          for (i = 0; i < index; i++)
            add_chain(N, QARRAY[i]);
        } else if (strcmp(type, "LUCAS") == 0) {
          int i, index;
          for (index = 0; index < MAX_QARRAY; index++)
            mpz_set_ui(QARRAY[index], 0);
          index = 1;
          while (1) {
            get_line(0);
            if        (sscanf(_line, "N %s", _vstr) == 1) {
              mpz_set_str(N, _vstr, _base);
            } else if (sscanf(_line, "Q[%d] %s", &i, _vstr) == 2) {
              if (i != index) quit_error("Lucas", "Invalid Q index");
              mpz_set_str(QARRAY[i], _vstr, _base);
              index++;
            } else if (sscanf(_line, "A %s", _vstr) == 1) {
              mpz_set_str(A, _vstr, _base);
              break;
            }
          }
          verify_lucas(index);
          for (i = 1; i < index; i++)
            add_chain(N, QARRAY[i]);
        } else {
          quit_error("Parsing", "Unknown type");
        }
      }
    }
  }
  verify_final();
}

static void dieusage(const char* prog) {
  printf("Verify Cert version 0.96.  Dana Jacobsen\n\n");
  printf("Usage: %s [options] <file>\n\n", prog);
  printf("Options:\n");
  printf("   -v     set verbose\n");
  printf("   -q     set quiet (no output, only exit code)\n");
  printf("   -help  this message\n");
  var_free();
  exit(RET_INVALID);
}


int main(int argc, char *argv[])
{
  int i;
  int optdone = 0;
  if (argc < 2) dieusage(argv[0]);

  var_init();
  mpz_set_ui(N, 0);

  for (i = 1; i < argc; i++) {
    if (!optdone && argv[i][0] == '-') {
      if (strcmp(argv[i], "--") == 0) {
        optdone = 1;
      } else if (argv[i][1] == '\0') {
        process_file("-");
      } else if (strcmp(argv[i], "-v") == 0) {
        _verbose++;
      } else if (strcmp(argv[i], "-q") == 0) {
        _quiet++;
      } else if (strcmp(argv[i], "-help") == 0 || strcmp(argv[i], "--help") == 0) {
        dieusage(argv[0]);
      } else {
        printf("Unknown option: %s\n\n", argv[i]);
        dieusage(argv[0]);
      }
      continue;
    }
    /* process_file will exit if not verified prime */
    process_file(argv[i]);
  }
  if (mpz_sgn(N) == 0)
    dieusage(argv[0]);
  quit_prime();
  exit(RET_PRIME);
}
