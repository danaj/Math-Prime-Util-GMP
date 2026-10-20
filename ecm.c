
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <gmp.h>

#include "ptypes.h"
#include "ecm.h"
#include "utility.h"
#include "prime_iterator.h"

#ifndef SIZE_MAX
#define SIZE_MAX ((size_t)-1)
#endif

#define TEST_FOR_2357(n, f) \
  { \
    if (mpz_divisible_ui_p(n, 2)) { mpz_set_ui(f, 2); return 1; } \
    if (mpz_divisible_ui_p(n, 3)) { mpz_set_ui(f, 3); return 1; } \
    if (mpz_divisible_ui_p(n, 5)) { mpz_set_ui(f, 5); return 1; } \
    if (mpz_divisible_ui_p(n, 7)) { mpz_set_ui(f, 7); return 1; } \
    if (mpz_cmp_ui(n, 121) < 0) { return 0; } \
  }

/* D is below 2^(UV bits/2).  Reserve three such widths so the stage-2
 * endpoints B2+D, m+D and the final m+=2*D cannot wrap.  Saturate both
 * explicit B2 values and the default without overflowing its multiplication. */
static UV ecm_resolve_b2(UV B1, UV B2)
{
  const UV max_B2 = UV_MAX - 3 * ((UV)1 << (sizeof(UV)*CHAR_BIT/2));
  if (B2 < B1)
    B2 = B1 > max_B2 / 100 ? max_B2 : 100 * B1;
  return B2 > max_B2 ? max_B2 : B2;
}

static void ec_affine_set(struct ec_affine_point *R,
                          const struct ec_affine_point P)
{
  mpz_set(R->x, P.x);  mpz_set(R->y, P.y);
  R->infinity = P.infinity;
}

static void ec_affine_set_infinity(struct ec_affine_point *R)
{
  mpz_set_ui(R->x, 0);  mpz_set_ui(R->y, 1);
  R->infinity = 1;
}

/* P3 = 2*P1.  P3 must be distinct from P1. */
static int _ec_add_2A(const mpz_t a,
                    const mpz_t n,
                    const struct ec_affine_point P1,
                    struct ec_affine_point *P3,
                    mpz_t d,
                    mpz_t m,
                    mpz_t t1,
                    mpz_t t2)
{
  if (P1.infinity || mpz_sgn(P1.y) == 0) {
    ec_affine_set_infinity(P3);
    return 0;
  }
  /* m = (3x1^2 + a) * (2y1)^-1 mod n */
  mpz_mul_ui(t1, P1.y, 2);
  if (!mpz_invert(m, t1, n)) {
    mpz_gcd(d, t1, n);
    return 1;
  }
  mpz_mul_ui(t1, P1.x, 3);
  mpz_mul(t2, t1, P1.x);
  mpz_add(t1, t2, a);
  mpz_mul(t2, m, t1);
  mpz_mod(m, t2, n);

  /* x3 = m^2 - 2x1 mod n */
  mpz_mul(t1, m, m);
  mpz_mul_ui(t2, P1.x, 2);
  mpz_sub(t1, t1, t2);
  mpz_mod(P3->x, t1, n);

  /* y3 = m(x1 - x3) - y1 mod n */
  mpz_sub(t1, P1.x, P3->x);
  mpz_mul(t2, t1, m);
  mpz_sub(t1, t2, P1.y);
  mpz_mod(P3->y, t1, n);
  P3->infinity = 0;
  return 0;
}

/* P3 = P1 + P2.  P3 must be distinct from both inputs. */
static int _ec_add_AB(const mpz_t a, const mpz_t n,
                     const struct ec_affine_point P1,
                     const struct ec_affine_point P2,
                     struct ec_affine_point *P3,
                     mpz_t d, mpz_t m, mpz_t t1, mpz_t t2)
{
  if (P1.infinity) { ec_affine_set(P3, P2); return 0; }
  if (P2.infinity) { ec_affine_set(P3, P1); return 0; }
  if (mpz_cmp(P1.x, P2.x) == 0) {
    if (mpz_cmp(P1.y, P2.y) == 0)
      return _ec_add_2A(a, n, P1, P3, d, m, t1, t2);
    mpz_add(t1, P1.y, P2.y);
    mpz_mod(t1, t1, n);
    if (mpz_sgn(t1) == 0) {
      ec_affine_set_infinity(P3);
      return 0;
    }
    /* Equal x but neither equal nor opposite y: the points agree in
     * some CRT components and are opposites in the others. */
    mpz_sub(t1, P1.y, P2.y);
    mpz_gcd(d, t1, n);
    return 1;
  }
  mpz_sub(t2, P2.x, P1.x);
  if (!mpz_invert(t1, t2, n)) {
    mpz_gcd(d, t2, n);
    return 1;
  }
  mpz_sub(t2, P2.y, P1.y);
  mpz_mul(m, t1, t2);
  mpz_mod(m, m, n);

  mpz_mul(t1, m, m);
  mpz_sub(t1, t1, P1.x);
  mpz_sub(t1, t1, P2.x);
  mpz_mod(P3->x, t1, n);
  mpz_sub(t1, P1.x, P3->x);
  mpz_mul(t2, m, t1);
  mpz_sub(t1, t2, P1.y);
  mpz_mod(P3->y, t1, n);
  P3->infinity = 0;
  return 0;
}

/* k >= 0 is consumed to zero on success, or partially consumed if an
 * inversion fails.  Copying P first permits R=P; see ecm.h for aliasing. */
int ec_affine_multiply(const mpz_t a, mpz_t k, const mpz_t n, const struct ec_affine_point P, struct ec_affine_point *R, mpz_t d)
{
  int found = 0;
  struct ec_affine_point A, B, C;
  mpz_t t, t2, t3;

  mpz_init(A.x); mpz_init(A.y);
  mpz_init(B.x); mpz_init(B.y);
  mpz_init(C.x); mpz_init(C.y);
  mpz_init(t);   mpz_init(t2);   mpz_init(t3);
  mpz_set_ui(d, 1);
  mpz_mod(A.x, P.x, n);  mpz_mod(A.y, P.y, n);
  A.infinity = P.infinity;
  ec_affine_set_infinity(&B);

  /* Binary ladder multiply. */
  while (mpz_cmp_ui(k, 0) > 0) {
    if (mpz_odd_p(k)) {
      found = _ec_add_AB(a, n, A, B, &C, d, t, t2, t3);
      if (found) break;
      ec_affine_set(&B, C);
      mpz_sub_ui(k, k, 1);
    } else {
      found = _ec_add_2A(a, n, A, &C, d, t, t2, t3);
      if (found) break;
      ec_affine_set(&A, C);
      mpz_tdiv_q_2exp(k, k, 1);
    }
  }
  /* A failed inversion is not a valid point at infinity. */
  if (!found) ec_affine_set(R, B);
  mpz_clear(t);   mpz_clear(t2);   mpz_clear(t3);
  mpz_clear(A.x); mpz_clear(A.y);
  mpz_clear(B.x); mpz_clear(B.y);
  mpz_clear(C.x); mpz_clear(C.y);

  return found;
}

int _GMP_ecm_factor_affine(const mpz_t n, mpz_t f, UV B1, UV ncurves)
{
  mpz_t a, mk;
  struct ec_affine_point X, Y;
  UV B, curve, q;

  TEST_FOR_2357(n, f);

  mpz_init(a);   mpz_init(mk);
  mpz_init(X.x); mpz_init(X.y);
  mpz_init(Y.x); mpz_init(Y.y);

  for (B = 100; B < B1*5; B *= 5) {
    if (B*5 > 2*B1) B = B1;
    for (curve = 0; curve < ncurves; curve++) {
      PRIME_ITERATOR(iter);
      mpz_isaac_urandomm(a, n);
      mpz_set_ui(X.x, 0); mpz_set_ui(X.y, 1);
      X.infinity = 0;
      for (q = 2; q < B; q = prime_iterator_next(&iter)) {
        UV k = q;
        UV kmin = B / q;

        while (k <= kmin)
          k *= q;

        mpz_set_uv(mk, k);
        if (ec_affine_multiply(a, mk, n, X, &Y, f)) {
          prime_iterator_destroy(&iter);
          mpz_clear(a);   mpz_clear(mk);
          mpz_clear(X.x); mpz_clear(X.y);
          mpz_clear(Y.x); mpz_clear(Y.y);
          return 1;
        }
        ec_affine_set(&X, Y);
        /* Check that we're not starting over */
        if (X.infinity)
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

/* Projective arithmetic shares these scratch variables, including their
 * initialization/clear in the driver.  Serialize entire projective calls:
 * this implementation is neither reentrant nor thread-safe. */
static mpz_t b, ecn;          /* used throughout ec mult */
static mpz_t u, v, w;         /* temporaries */
static mpz_t x1, z1, x2, z2;  /* used by ec_mult and stage2 */

#define mpz_mulmod(r, a, b, n, t)  \
  do { mpz_mul(t, a, b); mpz_mod(r, t, n); } while (0)

/* (x2:z2) = (x1:z1) + (x2:z2).  xinit must not alias an output coordinate. */
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
/* (xout:zout) = (x1:z1) + (x2:z2).  The output point may alias either
 * addend or the difference (xin:zin); its x and z must be distinct. */
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

/* (x2:z2) = 2(x1:z1), permitting an in-place point double. */
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

/* Full-width binary ladder, also used outside PRAC's arithmetic range.
 * The input point remains the difference of the two working points.
 * x and z are overwritten, must be distinct, and must not be scratch globals. */
static void ec_mult_binary(UV k, mpz_t x, mpz_t z)
{
  UV bit = 1;
  if (k == 0) {
    mpz_set_ui(x, 1); mpz_set_ui(z, 0);
    return;
  }
  if (k == 1 || mpz_sgn(z) == 0) return;
  /* (0:1) is a point of order two, not a valid differential input. */
  if (mpz_sgn(x) == 0) {
    if (!(k & 1)) { mpz_set_ui(x, 1); mpz_set_ui(z, 0); }
    return;
  }
  while (bit <= k/2) bit *= 2;
  mpz_set(x1, x); mpz_set(z1, z);
  ec_double(x2, z2, x, z);
  for (bit >>= 1; bit != 0; bit >>= 1) {
    if (k & bit) {
      ec_add3(x1, z1, x1, z1, x2, z2, x, z);
      ec_double(x2, z2, x2, z2);
    } else {
      ec_add3(x2, z2, x2, z2, x1, z1, x, z);
      ec_double(x1, z1, x1, z1);
    }
  }
  mpz_set(x, x1); mpz_set(z, z1);
}

/* PRAC
 *
 * "Evaluating recurrences of form X_{m+n} = f(X_m, X_n, X_{m-n}) via
 *  Lucas chains, Peter L. Montgomery, Dec 1983 (revised Jan 1992).
 *  https://cr.yp.to/bib/1992/montgomery-lucas.pdf
 *  http://research.microsoft.com/en-us/um/people/petmon/lucas.pdf
 *
 * "20 years of ECM" by Paul Zimmermann, 2006.
 *
 * Code derived from GMP-ECM (Zimmermann & Kruppa)
 */

/* PRAC, details from GMP-ECM, algorithm from Montgomery */
/* See "20 years of ECM" by Paul Zimmermann for more info */
static mpz_t x3, z3, x4, z4;  /* used by prac */
#define ADD 6 /* number of multiplications in an addition */
#define DUP 5 /* number of multiplications in a double */

/* Return the mulmod cost, or UV_MAX for an invalid chain.  Only a terminal
 * d=e=1 computes the full nP; some reductions absorb a common factor.
 * Caller limits n to UV_MAX/6, keeping the small multiples below safe. */
static UV lucas_cost(UV n, UV r)
{
  UV c, d, e, t;

  if (r <= n/2 || r >= n)
    return UV_MAX;
  d = n - r;
  e = 2 * r - n;
  c = DUP + ADD;  /* initial double and final add */
  while (d != e) {
    if (d < e) { t = d;  d = e;  e = t; }
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
  return d == 1 ? c : UV_MAX;
}

#define SWAP(a, b) \
  t = x##a; x##a = x##b; x##b = t;  t = z##a; z##a = z##b; z##b = t;

/* PRAC for odd scalars in its safe arithmetic range, binary otherwise.
 * x and z have the same destructive/nonaliasing contract as the binary ladder. */
static void ec_mult(UV k, mpz_t x, mpz_t z)
{
   UV d, e, r, cost, best_cost = UV_MAX, best_r = 0;
   unsigned int i;
   __mpz_struct *xA, *zA, *xB, *zB, *xC, *zC, *xT, *zT, *xT2, *zT2, *t;

   static double const val[] =
     {1.61803398875, 1.72360679775, 1.618347119656, 1.617914406529,
      1.58017872826};

   /* Bound all 2*r, 4*d, 5*e and cost arithmetic before choosing a chain.
    * Keep full UV scalars working without depending on floating rounding. */
   if (k <= 2 || !(k & 1) || k > UV_MAX/6 ||
       mpz_sgn(x) == 0 || mpz_sgn(z) == 0) {
     ec_mult_binary(k, x, z);
     return;
   }
   for (i = 0; i < sizeof(val)/sizeof(val[0]); i++) {
     r = (UV)((double)k / val[i] + 0.5);
     cost = lucas_cost(k, r);
     if (cost < best_cost) { best_cost = cost; best_r = r; }
   }
   /* All five can be invalid (e.g. k=7^4).  For odd k this reserve r
    * always satisfies gcd(k,r)=1, since 2*r-k=1. */
   r = best_r != 0 ? best_r : k/2 + 1;
   /* A=(x:z) B=(x1:z1) C=(x2:z2) T=T1=(x3:z3) T2=(x4:z4) */
   xA=x; zA=z; xB=x1; zB=z1; xC=x2; zC=z2; xT=x3; zT=z3; xT2=x4; zT2=z4;
   /* Initially kP = (k-r)*(2P) + (2*r-k)*P, with C=A-B=P. */
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
   /* Costing guarantees d=e=1, so the final addition gives kP. */
   ec_add3(xA,zA,xA,zA,xB,zB,xC,zC);
   if (x!=xA) { mpz_set(x,xA); mpz_set(z,zA); }
}

/* Normalize in place.  Return 1 on success (f=1), or 0 with gcd(z,n)
 * in f.  The latter can be n, meaning a lost curve rather than a factor.
 * x, z and f must be distinct and must not alias scratch u, v or modulus ecn. */
static int ec_normalize(mpz_t x, mpz_t z, mpz_t f)
{
  mpz_gcdext(f, u, NULL, z, ecn);
  if (mpz_cmp_ui(f, 1) != 0) return 0;
  mpz_mulmod(x, x, u, ecn, v);
  mpz_set_ui(z, 1);
  return 1;
}

static int ec_stage2(UV B1, UV B2, mpz_t x, mpz_t z, mpz_t f)
{
  UV D, i, m, count, initialized = 0;
  mpz_t* nqx;
  mpz_t g, one;
  PRIME_ITERATOR(iter);

  if (B2 <= B1 || B2 < 2) return 0;
  D = sqrt((double)B2 / 2.0);
  if (D < 2) D = 2;
  else if (D%2) D++;
  count = 2*D+1;  /* sqrt(UV_MAX/2), rounded up even, keeps this in UV. */
  /* Compare before Newz can narrow the count to size_t (or Perl MEM_SIZE).
   * In particular, UV can be 64 bits while allocation sizes are 32 bits. */
  if (count > SIZE_MAX / sizeof(mpz_t))
    croak("ECM stage-2 allocation size overflow");

  if (!ec_normalize(x, z, f))
    return mpz_cmp(f, ecn) != 0 ? 2 : 0;

  Newz(0, nqx, (size_t)count, mpz_t);
  mpz_init_set(nqx[1], x);
  initialized = 1;
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
    initialized = i;
    if (!ec_normalize(nqx[i], z2, f)) goto cleanup;
  }

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
      if (!ec_normalize(x1, z1, f)) goto cleanup;
      mpz_set(x, x2);  mpz_set(z, z2);
    }
    if (m+D > B1 && m >= D) {
      prime_iterator_setprime(&iter, m-D-1);
      for (i = prime_iterator_next(&iter); i < m; i = prime_iterator_next(&iter)) {
        /* Giant step mQ, baby step (m-i)Q; 1 <= m-i <= D. */
        mpz_sub(w, x1, nqx[m-i]);
        mpz_mulmod(g, g, w, ecn, u);
      }
      for ( ; i <= m+D; i = prime_iterator_next(&iter)) {
        if (i > m && !prime_iterator_isprime(&iter, m-(i-m))) {
          /* Right-half baby step (i-m)Q; 1 <= i-m <= D. */
          mpz_sub(w, x1, nqx[i-m]);
          mpz_mulmod(g, g, w, ecn, u);
        }
      }
      mpz_gcd(f, g, ecn);
      if (mpz_cmp_ui(f, 1) != 0) goto cleanup;
    }
  }

cleanup:
  prime_iterator_destroy(&iter);
  for (i = 1; i <= initialized; i++) mpz_clear(nqx[i]);
  Safefree(nqx);
  mpz_clear(g);
  mpz_clear(one);
  return mpz_cmp_ui(f, 1) != 0 && mpz_cmp(f, ecn) != 0 ? 2 : 0;
}

int _GMP_ecm_factor_projective(const mpz_t n, mpz_t f, UV B1, UV B2, UV ncurves)
{
  mpz_t sigma, a, x, z;
  UV i, curve, q, k;
  int found = 0;
  int _verbose = get_verbose_level();

  TEST_FOR_2357(n, f);

  B2 = ecm_resolve_b2(B1, B2);  /* Default B2/B1 ratio: 100. */

  mpz_init_set(ecn, n);
  mpz_init(x);   mpz_init(z);   mpz_init(a);   mpz_init(b);
  mpz_init(u);   mpz_init(v);   mpz_init(w);   mpz_init(sigma);
  mpz_init(x1);  mpz_init(z1);  mpz_init(x2);  mpz_init(z2);
  mpz_init(x3);  mpz_init(z3);  mpz_init(x4);  mpz_init(z4);

  if (_verbose>2) gmp_printf("# ecm trying %Zd (B1=%"UVuf" B2=%"UVuf" ncurves=%"UVuf")\n", n, B1, B2, ncurves);

  for (curve = 0; curve < ncurves; curve++) {
    PRIME_ITERATOR(iter);
    do {
      mpz_isaac_urandomm(sigma, n);
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
    found = mpz_cmp_ui(f, 1) != 0;
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
    for (q = 2; q < B1; q *= 2) {
      ec_double(x, z, x, z);
      if (q > (B1-1)/2) break; /* The next power reaches B1 or would wrap. */
    }
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
    found = !ec_normalize(x, z, f);
    if (!found) {
      mpz_gcd(f, sigma, ecn);
      found = mpz_cmp_ui(f, 1) != 0;
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
  mpz_clear(x3);  mpz_clear(z3);  mpz_clear(x4);  mpz_clear(z4);

  return found;
}
