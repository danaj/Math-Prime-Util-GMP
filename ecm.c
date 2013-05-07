
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include "ecm.h"
#include "prime_iterator.h"
#include "gmp_main.h"

#define USE_PRAC

#define TEST_FOR_2357(n, f) \
  { \
    if (mpz_divisible_ui_p(n, 2)) { mpz_set_ui(f, 2); return 1; } \
    if (mpz_divisible_ui_p(n, 3)) { mpz_set_ui(f, 3); return 1; } \
    if (mpz_divisible_ui_p(n, 5)) { mpz_set_ui(f, 5); return 1; } \
    if (mpz_divisible_ui_p(n, 7)) { mpz_set_ui(f, 7); return 1; } \
    if (mpz_cmp_ui(n, 121) < 0) { return 0; } \
  }

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

int ec_affine_multiply(mpz_t a, mpz_t k, mpz_t n, struct ec_affine_point P, struct ec_affine_point *R, mpz_t d)
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

int _GMP_ecm_factor_affine(mpz_t n, mpz_t f, UV B1, UV ncurves)
{
  mpz_t a, mk;
  struct ec_affine_point X, Y;
  UV B, curve, q;
  gmp_randstate_t* p_randstate = _GMP_get_randstate();

  TEST_FOR_2357(n, f);

  mpz_init(a);   mpz_init(mk);
  mpz_init(X.x); mpz_init(X.y);
  mpz_init(Y.x); mpz_init(Y.y);

  for (B = 100; B < B1*5; B *= 5) {
    if (B*5 > 2*B1) B = B1;
    for (curve = 0; curve < ncurves; curve++) {
      PRIME_ITERATOR(iter);
      mpz_urandomm(a, *p_randstate, n);
      mpz_set_ui(X.x, 0); mpz_set_ui(X.y, 1);
      for (q = 2; q < B; q = prime_iterator_next(&iter)) {
        UV k = q;
        UV kmin = B / q;

        while (k <= kmin)
          k *= q;

        mpz_set_ui(mk, k);
        if (ec_affine_multiply(a, mk, n, X, &Y, f)) {
          prime_iterator_destroy(&iter);
          mpz_clear(a);
          mpz_clear(X.x); mpz_clear(X.y);
          mpz_clear(Y.x); mpz_clear(Y.y);
          return 1;
        }
        mpz_set(X.x, Y.x);  mpz_set(X.y, Y.y);
        /* Check that we're not starting over */
        if ( !mpz_cmp_ui(X.x, 0) && !mpz_cmp_ui(X.y, 1) )
          break;
      }
      prime_iterator_destroy(&iter);
    }
  }

  mpz_clear(a);   mpz_clear(mk);
  mpz_clear(X.x); mpz_clear(X.y);
  mpz_clear(Y.x); mpz_clear(Y.y);

  return 0;
}


/*******************************************************************/

/* A better ECM, with a stage 2.
 * Heavily inspired by GMP-ECM, written by Paul Zimmermann (1998),
 * especially the stage 2 method.  Also see "The elliptic curve
 * integer factorization method" by Bosma and Lenstra as well as many
 * other articles.
 */

static mpz_t b, ecn;          /* used throughout ec mult */
static mpz_t u, v, w;         /* temporaries */
static mpz_t x1, z1, x2, z2;  /* used by ec_mult and stage2 */

#define mpz_mulmod(r, a, b, n, t)  \
  do { mpz_mul(t, a, b); mpz_mod(r, t, n); } while (0)

/* (x2:z2) = (x1:z1) + (x2:z2) */
static void ec_add(mpz_t x2, mpz_t z2, mpz_t x1, mpz_t z1, mpz_t xinit)
{
  mpz_sub(u, x2, z2);
  mpz_add(v, x1, z1);
  mpz_mulmod(u, u, v, ecn, w);   /* u = (x2 - z2) * (x1 + z1) % n */

  mpz_add(v, x2, z2);
  mpz_sub(w, x1, z1);
  mpz_mulmod(v, v, w, ecn, x2);  /* v = (x2 + z2) * (x1 - z1) % n */

  mpz_add(w, u, v);
  mpz_mulmod(x2, w, w, ecn, z2); /* x2 = (u+v)^2 % n */

  mpz_sub(w, u, v);
  mpz_mulmod(z2, w, w, ecn, v);  /* z2 = (u-v)^2 % n */

  mpz_mulmod(z2, xinit, z2, ecn, v); /* z2 *= X1. */
  /* Per Montgomery 1987, we set Z1 to 1, so no need for x2 *= Z1 */
  /* 5 mulmods, 6 adds */
}

/* This version assumes no normalization, so uses an extra mulmod. */
/* (xout:zout) = (x1:z1) + (x2:z2) */
static void ec_add3(mpz_t xout, mpz_t zout,
                    mpz_t x1, mpz_t z1,
                    mpz_t x2, mpz_t z2,
                    mpz_t xin, mpz_t zin)
{
  mpz_sub(u, x2, z2);
  mpz_add(v, x1, z1);
  mpz_mulmod(u, u, v, ecn, w);   /* u = (x2 - z2) * (x1 + z1) % n */

  mpz_add(v, x2, z2);
  mpz_sub(w, x1, z1);
  mpz_mulmod(v, v, w, ecn, v);   /* v = (x2 + z2) * (x1 - z1) % n */

  mpz_add(w, u, v);              /* w = u+v */
  mpz_sub(v, u, v);              /* v = u-v */

  mpz_mulmod(w, w, w, ecn, u);   /* w = (u+v)^2 % n */
  mpz_mulmod(v, v, v, ecn, u);   /* v = (u-v)^2 % n */

  mpz_set(u, xin);
  mpz_mulmod(xout, w, zin, ecn, w);
  mpz_mulmod(zout, v, u,   ecn, w);
  /* 6 mulmods, 6 adds */
}

/* (x2:z2) = 2(x1:z1) */
static void ec_double(mpz_t x2, mpz_t z2, mpz_t x1, mpz_t z1)
{
  mpz_add(u, x1, z1);
  mpz_mulmod(u, u, u, ecn, w);   /* u = (x1+z1)^2 % n */

  mpz_sub(v, x1, z1);
  mpz_mulmod(v, v, v, ecn, w);   /* v = (x1-z1)^2 % n */

  mpz_mulmod(x2, u, v, ecn, w);  /* x2 = uv % n */

  mpz_sub(w, u, v);              /* w = u-v = 4(x1 * z1) */
  mpz_mulmod(u, b, w, ecn, z2);
  mpz_add(u, u, v);              /* u = (v+b*w) mod n */
  mpz_mulmod(z2, w, u, ecn, v);  /* z2 = (w*u) mod n */
  /* 5 mulmods, 4 adds */
}

/* See http://alexandria.tue.nl/extra2/200311829.pdf for lots of discussion
 * and algorithms for various addition chains.
 */

#ifndef USE_PRAC

static void ec_mult(UV k, mpz_t x, mpz_t z)
{
  int l, r;

  r = --k; l = -1; while (r != 1) { r >>= 1; l++; }
  if (k & ( UVCONST(1)<<l)) {
    ec_double(x2, z2, x, z);
    ec_add3(x1, z1, x2, z2, x, z, x, z);
    ec_double(x2, z2, x2, z2);
  } else {
    ec_double(x1, z1, x, z);
    ec_add3(x2, z2, x, z, x1, z1, x, z);
  }
  l--;
  while (l >= 1) {
    if (k & ( UVCONST(1)<<l)) {
      ec_add3(x1, z1, x1, z1, x2, z2, x, z);
      ec_double(x2, z2, x2, z2);
    } else {
      ec_add3(x2, z2, x2, z2, x1, z1, x, z);
      ec_double(x1, z1, x1, z1);
    }
    l--;
  }
  if (k & 1) {
    ec_double(x, z, x2, z2);
  } else {
    ec_add3(x, z, x2, z2, x1, z1, x, z);
  }
}

#else

/* PRAC
 *
 * "Evaluating recurrences of form X_{m+n} = f(X_m, X_n, X_{m-n}) via
 *  Lucas chains, Peter L. Montgomery, Dec 1983 (revised Jan 1992).
 *  http://research.microsoft.com/en-us/um/people/petmon/lucas.pdf
 *
 * "20 years of ECM" by Paul Zimmerman, 2006.
 *
 * Code derived from GMP-ECM (Zimmermann & Kruppa)
 */

/* PRAC, details from GMP-ECM, algorithm from Montgomery */
/* See "20 years of ECM" by Paul Zimmermann for more info */
static mpz_t x3, z3, x4, z4;  /* used by prac */
#define ADD 6 /* number of multiplications in an addition */
#define DUP 5 /* number of multiplications in a double */

/* Returns the number of mulmods */
static UV lucas_cost(UV n, double v)
{
  UV c, d, e, r;

  d = n;
  r = (UV) ( ((double)d / v) + 0.5 );
  if (r >=n )
    return(ADD*n);
  d = n - r;
  e = 2 * r - n;
  c = DUP + ADD;  /* initial double and final add */
  while (d != e) {
    if (d < e) { r = d;  d = e;  e = r; }
    if (4 * d <= 5 * e && ((d + e) % 3) == 0) {
 /*C1*/ d = (2 * d - e) / 3;  e = (e - d) / 2;  c += 3 * ADD; /* 3 adds */
    } else if (4 * d <= 5 * e && (d - e) % 6 == 0) {
 /*C2*/ d = (d - e) / 2;  c += ADD + DUP; /* one add, one double */
    } else if (d <= 4 * e) {
 /*C3*/ d -= e;  c += ADD; /* one add */
    } else if ((d + e) % 2 == 0) {
 /*C4*/ d = (d - e) / 2;  c += ADD + DUP; /* one add, one double */
    } else if (d % 2 == 0) {
 /*C5*/ d /= 2;  c += ADD + DUP; /* one add, one double */
    } else if (d % 3 == 0) {
 /*C6*/ d = d / 3 - e;  c += 3 * ADD + DUP; /* three adds, one double */
    } else if ((d + e) % 3 == 0) {
 /*C7*/ d = (d - 2 * e) / 3;  c += 3 * ADD + DUP; /* three adds, one double */
    } else if ((d - e) % 3 == 0) {
 /*C8*/ d = (d - e) / 3;  c += 3 * ADD + DUP; /* three adds, one double */
    } else {
 /*C9*/ e /= 2;  c += ADD + DUP; /* one add, one double */
    }
  }
  return(c);
}

#define SWAP(a, b) \
  t = x##a; x##a = x##b; x##b = t;  t = z##a; z##a = z##b; z##b = t;

/* PRAC: computes kP from P=(x:z) and puts the result in (x:z). Assumes k>2.*/
static void ec_mult(UV k, mpz_t x, mpz_t z)
{
   unsigned int  d, e, r, i;
   __mpz_struct *xA, *zA, *xB, *zB, *xC, *zC, *xT, *zT, *xT2, *zT2, *t;

   static double const v[] =
     {1.61803398875, 1.72360679775, 1.618347119656, 1.617914406529,
      1.58017872826};

   /* chooses the best value of v */
   r = ADD * k;
   i = 0;
   for (d = 0; d < 5; d++) {
     e = lucas_cost(k, v[d]);
     if (e < r) { r = e;  i = d; }
   }
   r = (unsigned int)((double)k / v[i] + 0.5);
   /* A=(x:z) B=(x1:z1) C=(x2:z2) T=T1=(x3:z3) T2=(x4:z4) */
   xA=x; zA=z; xB=x1; zB=z1; xC=x2; zC=z2; xT=x3; zT=z3; xT2=x4; zT2=z4;
   /* first iteration always begins by Condition 3, then a swap */
   d = k - r;
   e = 2 * r - k;
   mpz_set(xB,xA); mpz_set(zB,zA); /* B=A */
   mpz_set(xC,xA); mpz_set(zC,zA); /* C=A */
   ec_double(xA,zA,xA,zA);         /* A=2*A */
   while (d != e) {
      if (d < e) {
         r = d;  d = e;  e = r;
         SWAP(A,B);
      }
      /* do the first line of Table 4 whose condition qualifies */
      if (4 * d <= 5 * e && ((d + e) % 3) == 0) { /* condition 1 */
         d = (2 * d - e) / 3;
         e = (e - d) / 2;
         ec_add3(xT,zT,xA,zA,xB,zB,xC,zC);   /* T = f(A,B,C) */
         ec_add3(xT2,zT2,xT,zT,xA,zA,xB,zB); /* T2= f(T,A,B) */
         ec_add3(xB,zB,xB,zB,xT,zT,xA,zA);   /* B = f(B,T,A) */
         SWAP(A,T2);
      } else if (4 * d <= 5 * e && (d - e) % 6 == 0) { /* condition 2 */
         d = (d - e) / 2;
         ec_add3(xB,zB,xA,zA,xB,zB,xC,zC);   /* B = f(A,B,C) */
         ec_double(xA,zA,xA,zA);             /* A = 2*A */
      } else if (d <= (4 * e)) { /* condition 3 */
         d -= e;
         ec_add3(xC,zC,xB,zB,xA,zA,xC,zC);   /* C = f(B,A,C) */
         SWAP(B,C);
      } else if ((d + e) % 2 == 0) { /* condition 4 */
         d = (d - e) / 2;
         ec_add3(xB,zB,xB,zB,xA,zA,xC,zC);   /* B = f(B,A,C) */
         ec_double(xA,zA,xA,zA);             /* A = 2*A */
      } else if (d % 2 == 0) { /* condition 5 */
         d /= 2;
         ec_add3(xC,zC,xC,zC,xA,zA,xB,zB);   /* C = f(C,A,B) */
         ec_double(xA,zA,xA,zA);             /* A = 2*A */
      } else if (d % 3 == 0) { /* condition 6 */
         d = d / 3 - e;
         ec_double(xT,zT,xA,zA);             /* T = 2*A */
         ec_add3(xT2,zT2,xA,zA,xB,zB,xC,zC); /* T2= f(A,B,C) */
         ec_add3(xA,zA,xT,zT,xA,zA,xA,zA);   /* A = f(T,A,A) */
         ec_add3(xC,zC,xT,zT,xT2,zT2,xC,zC); /* C = f(T,T2,C) */
         SWAP(B,C);
      } else if ((d + e) % 3 == 0) { /* condition 7 */
         d = (d - 2 * e) / 3;
         ec_add3(xT,zT,xA,zA,xB,zB,xC,zC);   /* T = f(A,B,C) */
         ec_add3(xB,zB,xT,zT,xA,zA,xB,zB);   /* B = f(T1,A,B) */
         ec_double(xT,zT,xA,zA);
         ec_add3(xA,zA,xA,zA,xT,zT,xA,zA);   /* A = 3*A */
      } else if ((d - e) % 3 == 0) { /* condition 8 */
         d = (d - e) / 3;
         ec_add3(xT,zT,xA,zA,xB,zB,xC,zC);   /* T = f(A,B,C) */
         ec_add3(xC,zC,xC,zC,xA,zA,xB,zB);   /* C = f(A,C,B) */
         SWAP(B,T);
         ec_double(xT,zT,xA,zA);
         ec_add3(xA,zA,xA,zA,xT,zT,xA,zA);   /* A = 3*A */
      } else { /* condition 9 */
         e /= 2;
         ec_add3(xC,zC,xC,zC,xB,zB,xA,zA);   /* C = f(C,B,A) */
         ec_double(xB,zB,xB,zB);             /* B = 2*B */
      }
   }
   ec_add3(xA,zA,xA,zA,xB,zB,xC,zC);
   if (x!=xA) { mpz_set(x,xA); mpz_set(z,zA); }
}

#endif /* PRAC */

#define NORMALIZE(f, u, v, x, z, n) \
    mpz_gcdext(f, u, NULL, z, n); \
    found = mpz_cmp_ui(f, 1); \
    if (found) break; \
    mpz_mulmod(x, x, u, n, v); \
    mpz_set_ui(z, 1);

static int ec_stage2(UV B1, UV B2, mpz_t x, mpz_t z, mpz_t f)
{
  UV D, i, m;
  mpz_t* nqx = 0;
  mpz_t g, one;
  int found;
  PRIME_ITERATOR(iter);

  do {
    NORMALIZE(f, u, v, x, z, ecn);

    D = sqrt( (double)B2 / 2.0 );
    if (D%2) D++;

    /* We really only need half of these. Only even values used. */
    Newz(0, nqx, 2*D+1, mpz_t);
    mpz_init_set(nqx[1], x);
    mpz_init_set_ui(g, 1);
    mpz_init_set_ui(one, 1);

    for (i = 2; i <= 2*D; i++) {
      if (i % 2) {
        mpz_set(x2, nqx[(i+1)/2]);  mpz_set_ui(z2, 1);
        ec_add(x2, z2, nqx[(i-1)/2], one, x);
      } else {
        ec_double(x2, z2, nqx[i/2], one);
      }
      mpz_init_set(nqx[i], x2);
      NORMALIZE(f, u, v, nqx[i], z2, ecn);
    }
    if (found) break;

    mpz_set(x1, x);
    mpz_set(z1, z);
    mpz_set(x, nqx[2*D-1]);
    mpz_set_ui(z, 1);

    /* See Zimmermann, "20 Years of ECM" slides, 2006, page 11-12 */
    for (m = 1; m < B2+D; m += 2*D) {
      if (m != 1) {
        mpz_set(x2, x1);
        mpz_set(z2, z1);
        ec_add(x1, z1, nqx[2*D], one, x);
        NORMALIZE(f, u, v, x1, z1, ecn);
        mpz_set(x, x2);  mpz_set(z, z2);
      }
      if (m+D > B1) {
        prime_iterator_setprime(&iter, m-D-1);
        for (i = prime_iterator_next(&iter); i < m; i = prime_iterator_next(&iter)) {
          /* if (m+D-i<1 || m+D-i>2*D) croak("index %lu range\n",i-(m-D)); */
          mpz_sub(w, x1, nqx[m+D-i]);
          mpz_mulmod(g, g, w, ecn, u);
        }
        for ( ; i <= m+D; i = prime_iterator_next(&iter)) {
          if (i > m && !prime_iterator_isprime(&iter, m+m-i)) {
            /* if (i-m<1 || i-m>2*D) croak("index %lu range\n",i-(m-D)); */
            mpz_sub(w, x1, nqx[i-m]);
            mpz_mulmod(g, g, w, ecn, u);
          }
        }
        mpz_gcd(f, g, ecn);
        found = mpz_cmp_ui(f, 1);
        if (found) break;
      }
    }
  } while (0);
  prime_iterator_destroy(&iter);

  if (nqx != 0) {
    for (i = 1; i <= 2*D; i++) {
      if (nqx[i] != 0)
        mpz_clear(nqx[i]);
    }
    Safefree(nqx);
    mpz_clear(g);
    mpz_clear(one);
  }
  if (found && !mpz_cmp(f, ecn)) found = 0;
  return (found) ? 2 : 0;
}

int _GMP_ecm_factor_projective(mpz_t n, mpz_t f, UV B1, UV ncurves)
{
  mpz_t sigma, a, x, z;
  UV i, curve, q, k;
  UV B2 = 100*B1;  /* time(S1) == time(S2) ~ 125 */
  int found = 0;
  gmp_randstate_t* p_randstate = _GMP_get_randstate();
  int _verbose = _GMP_get_verbose();

  TEST_FOR_2357(n, f);

  mpz_init_set(ecn, n);
  mpz_init(x);   mpz_init(z);   mpz_init(a);   mpz_init(b);
  mpz_init(u);   mpz_init(v);   mpz_init(w);   mpz_init(sigma);
  mpz_init(x1);  mpz_init(z1);  mpz_init(x2);  mpz_init(z2);
#ifdef USE_PRAC
  mpz_init(x3);  mpz_init(z3);  mpz_init(x4);  mpz_init(z4);
#endif

  if (_verbose>2) gmp_printf("# ecm trying %Zd (B1=%lu B2=%lu ncurves=%lu)\n", n, (unsigned long)B1, (unsigned long)B2, (unsigned long)ncurves);

  for (curve = 0; curve < ncurves; curve++) {
    PRIME_ITERATOR(iter);
    do {
      mpz_urandomm(sigma, *p_randstate, n);
    } while (mpz_cmp_ui(sigma, 5) <= 0);
    mpz_mul_ui(w, sigma, 4);
    mpz_mod(v, w, n);             /* v = 4σ */

    mpz_mul(x, sigma, sigma);
    mpz_sub_ui(w, x, 5);
    mpz_mod(u, w, n);             /* u = σ^2-5 */

    mpz_mul(x, u, u);
    mpz_mulmod(x, x, u, n, w);    /* x = u^3 */

    mpz_mul(z, v, v);
    mpz_mulmod(z, z, v, n, w);    /* z = v^3 */

    mpz_mul(b, x, v);
    mpz_mul_ui(w, b, 4);
    mpz_mod(b, w, n);             /* b = 4 u^3 v */

    mpz_sub(a, v, u);
    mpz_mul(w, a, a);
    mpz_mulmod(w, w, a, n, w);

    mpz_mul_ui(a, u, 3);
    mpz_add(a, a, v);
    mpz_mul(w, w, a);
    mpz_mod(a, w, n);             /* a = ((v-u)^3 * (3*u + v)) % n */

    mpz_gcdext(f, u, NULL, b, n);
    found = mpz_cmp_ui(f, 1);
    if (found) { if (!mpz_cmp(f, n)) { found = 0; continue; } break; }
    mpz_mul(a, a, u);

    mpz_sub_ui(a, a, 2);
    mpz_mod(a, a, n);

    mpz_add_ui(b, a, 2);
    if (mpz_mod_ui(w, b, 2)) mpz_add(b, b, n);
    mpz_tdiv_q_2exp(b, b, 1);
    if (mpz_mod_ui(w, b, 2)) mpz_add(b, b, n);
    mpz_tdiv_q_2exp(b, b, 1);

    /* Use sigma to collect possible factors */
    mpz_set_ui(sigma, 1);

    /* Stage 1 */
    for (q = 2; q < B1; q *= 2)
      ec_double(x, z, x, z);
    mpz_mulmod(sigma, sigma, x, ecn, w);
    i = 15;
    for (q = prime_iterator_next(&iter); q < B1; q = prime_iterator_next(&iter)) {
      /* PRAC is a little faster with:
       *     for (k = 1; k <= B1/q; k *= q)
       *       ec_mult(q, x, z);
       * but binary multiplication is much slower that way. */
      for (k = q; k <= B1/q; k *= q) ;
      ec_mult(k, x, z);
      mpz_mulmod(sigma, sigma, x, ecn, w);
      if (i++ % 32 == 0) {
        mpz_gcd(f, sigma, ecn);
        if (mpz_cmp_ui(f, 1))  break;
      }
    }
    prime_iterator_destroy(&iter);

    /* Find factor in S1 */
    do { NORMALIZE(f, u, v, x, z, n); } while (0);
    if (!found) {
      mpz_gcd(f, sigma, ecn);
      found = mpz_cmp_ui(f, 1);
    }
    if (found) { if (!mpz_cmp(f, n)) { found = 0; continue; } break; }

    /* Stage 2 */
    if (!found && B2 > B1)
      found = ec_stage2(B1, B2, x, z, f);

    if (found) { if (!mpz_cmp(f, n)) { found = 0; continue; } break; }
  }
  if (_verbose>2) {
    if (found) gmp_printf("# ecm: %Zd in stage %d\n", f, found);
    else       gmp_printf("# ecm: no factor\n");
  }

  mpz_clear(ecn);
  mpz_clear(x);   mpz_clear(z);   mpz_clear(a);   mpz_clear(b);
  mpz_clear(u);   mpz_clear(v);   mpz_clear(w);   mpz_clear(sigma);
  mpz_clear(x1);  mpz_clear(z1);  mpz_clear(x2);  mpz_clear(z2);
#ifdef USE_PRAC
  mpz_clear(x3);  mpz_clear(z3);  mpz_clear(x4);  mpz_clear(z4);
#endif

  return found;
}
