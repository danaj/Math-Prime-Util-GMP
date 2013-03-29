
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

#define TEST_FOR_2357(n, f) \
  { \
    if (mpz_divisible_ui_p(n, 2)) { mpz_set_ui(f, 2); return 1; } \
    if (mpz_divisible_ui_p(n, 3)) { mpz_set_ui(f, 3); return 1; } \
    if (mpz_divisible_ui_p(n, 5)) { mpz_set_ui(f, 5); return 1; } \
    if (mpz_divisible_ui_p(n, 7)) { mpz_set_ui(f, 7); return 1; } \
    if (mpz_cmp_ui(n, 121) < 0) { return 0; } \
  }

struct _ec_point  { mpz_t x, y; };

/* P3 = P1 + P2 */
static void _ec_add_AB(mpz_t n,
                    struct _ec_point P1,
                    struct _ec_point P2,
                    struct _ec_point *P3,
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

  /* m = (y2 - y1) * (x2 - x1)^-1 mod n */
  if (!mpz_invert(t1, t2, n)) {
    /* We've found a factor!  In multiply, gcd(mult,n) will be a factor. */
    mpz_set_ui(P3->x, 0);
    mpz_set_ui(P3->y, 1);
    return;
  }

  mpz_sub(m, P2.y, P1.y);
  mpz_mod(t2, m, n);        /* t2 = deltay */
  mpz_mul(m, t1, t2);
  mpz_mod(m, m, n);         /* m = deltay / deltax */

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
                    struct _ec_point P1,
                    struct _ec_point *P3,
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

static int _ec_multiply(mpz_t a, UV k, mpz_t n, struct _ec_point P, struct _ec_point *R, mpz_t d)
{
  int found = 0;
  struct _ec_point A, B, C;
  mpz_t t, t2, t3, mult;

  mpz_init(A.x); mpz_init(A.y);
  mpz_init(B.x); mpz_init(B.y);
  mpz_init(C.x); mpz_init(C.y);
  mpz_init(t);   mpz_init(t2);   mpz_init(t3);
  mpz_init_set_ui(mult, 1);  /* holds intermediates, gcd at end */

  mpz_set(A.x, P.x);  mpz_set(A.y, P.y);
  mpz_set_ui(B.x, 0); mpz_set_ui(B.y, 1);

  /* Binary ladder multiply.  Should investigate Lucas chains. */
  while (k > 0) {
    if ( k & 1 ) {
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
      k--;
    } else {
      mpz_mul_ui(t, A.y, 2);
      mpz_mul(t2, mult, t);
      mpz_mod(mult, t2, n);

      _ec_add_2A(a, n, A, &C, t, t2, t3);
      mpz_set(A.x, C.x);  mpz_set(A.y, C.y);
      k >>= 1;
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
  mpz_t a;
  struct _ec_point X, Y;
  UV B, curve, q;
  gmp_randstate_t* p_randstate = _GMP_get_randstate();

  TEST_FOR_2357(n, f);

  mpz_init(a);
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

        if (_ec_multiply(a, k, n, X, &Y, f)) {
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

  mpz_clear(a);
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
static mpz_t x1pz1, x1mz1;    /* intermediates */

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
}

/* #define ec_adddouble(x1, z1, x2, z2, x) ec_add(x1, z1, x2, z2, x); ec_double(x2, z2, x2, z2); */
static void ec_adddouble(mpz_t x2, mpz_t z2, mpz_t x1, mpz_t z1, mpz_t xinit)
{
  mpz_add(x1pz1, x1, z1);
  mpz_sub(x1mz1, x1, z1);

  /* ADD:  x2:z2 += x1:z1 */
  mpz_sub(u, x2, z2);
  mpz_mulmod(u, u, x1pz1, ecn, w); /* u = (x2 - z2) * (x1 + z1) % n */

  mpz_add(v, x2, z2);
  mpz_mulmod(v, v, x1mz1, ecn, w); /* v = (x2 + z2) * (x1 - z1) % n */

  mpz_add(w, u, v);
  mpz_mulmod(x2, w, w, ecn, z2); /* x2 = (u+v)^2 % n */

  mpz_sub(w, u, v);
  mpz_mulmod(z2, w, w, ecn, v);  /* z2 = (u-v)^2 % n */
  mpz_mulmod(z2, xinit, z2, ecn, v); /* z2 *= X1. */

  /* DOUBLE: x1:z1 = 2*(x1:z1) */
  mpz_mulmod(u, x1pz1, x1pz1, ecn, w); /* u = (x1+z1)^2 % n */
  mpz_mulmod(v, x1mz1, x1mz1, ecn, w); /* v = (x1-z1)^2 % n */

  mpz_mulmod(x1, u, v, ecn, w);  /* x1 = uv % n */
  mpz_sub(w, u, v);              /* w = u-v = 4(x1 * z1) */
  mpz_mulmod(u, b, w, ecn, z1);
  mpz_add(u, u, v);              /* u = (v+b*w) mod n */
  mpz_mulmod(z1, w, u, ecn, v);  /* z1 = (w*u) mod n */
}

static void ec_mult(UV k, mpz_t x, mpz_t z)
{
  int l, r;

  r = --k; l = -1; while (r != 1) { r >>= 1; l++; }
  if (k & ( UVCONST(1)<<l)) {
    mpz_set(x1, x);
    mpz_set(z1, z);
    ec_double(x2, z2, x1, z1);
    /* (x1:z1) = (x1:z1) + (x2:z2);  (x2:z2) = 2(x2:z2) */
    ec_adddouble(x1, z1, x2, z2, x);
  } else {
    mpz_set(x2, x);
    mpz_set(z2, z);
    ec_double(x1, z1, x2, z2);
    ec_add(x2, z2, x1, z1, x);
  }
  l--;
  while (l >= 1) {
    if (k & ( UVCONST(1)<<l)) {
      /* (x1:z1) = (x1:z1) + (x2:z2);  (x2:z2) = 2(x2:z2) */
      ec_adddouble(x1, z1, x2, z2, x);
    } else {
      /* (x2:z2) = (x1:z1) + (x2:z2);  (x1:z1) = 2(x1:z1) */
      ec_adddouble(x2, z2, x1, z1, x);
    }
    l--;
  }
  if (k & 1) {
    ec_double(x, z, x2, z2);
  } else {
    ec_add(x2, z2, x1, z1, x);
    mpz_set(x, x2);
    mpz_set(z, z2);
  }
}

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
  UV curve, q;
  UV B2 = 100*B1;
  int found = 0;
  gmp_randstate_t* p_randstate = _GMP_get_randstate();

  TEST_FOR_2357(n, f);

  mpz_init_set(ecn, n);
  mpz_init(x);  mpz_init(z);  mpz_init(sigma);  mpz_init(a);

  mpz_init(u);  mpz_init(v);  mpz_init(w);
  mpz_init(x1);  mpz_init(z1);
  mpz_init(x2);  mpz_init(z2);
  mpz_init(x1pz1);  mpz_init(x1mz1);
  mpz_init(b);

  for (curve = 0; curve < ncurves; curve++) {
    PRIME_ITERATOR(iter);
    do {
      mpz_urandomm(sigma, *p_randstate, n);
    } while (mpz_cmp_ui(sigma, 5) <= 0);
    mpz_mul_ui(w, sigma, 4);
    mpz_mod(v, w, n);

    mpz_mul(x, sigma, sigma);
    mpz_sub_ui(w, x, 5);
    mpz_mod(u, w, n);

    mpz_mul(x, u, u);
    mpz_mulmod(x, x, u, n, w);

    mpz_mul(z, v, v);
    mpz_mulmod(z, z, v, n, w);

    mpz_mul(b, x, v);
    mpz_mul_ui(w, b, 4);
    mpz_mod(b, w, n);

    mpz_sub(a, v, u);
    mpz_mul(w, a, a);
    mpz_mulmod(w, w, a, n, w);

    mpz_mul_ui(a, u, 3);
    mpz_add(a, a, v);
    mpz_mul(w, w, a);
    mpz_mod(a, w, n);     /* a = ((v-u)^3 * (3*u + v)) % n */

    mpz_gcdext(f, u, NULL, b, n);
    found = mpz_cmp_ui(f, 1);
    if (found) {
      if (!mpz_cmp(f, n)) { found = 0; continue; }
      break;
    }

    mpz_mul(a, a, u);
    mpz_sub_ui(a, a, 2);
    mpz_mod(a, a, n);

    mpz_add_ui(b, a, 2);
    if (mpz_mod_ui(w, b, 2)) mpz_add(b, b, n);
    mpz_tdiv_q_2exp(b, b, 1);
    if (mpz_mod_ui(w, b, 2)) mpz_add(b, b, n);
    mpz_tdiv_q_2exp(b, b, 1);

    do { NORMALIZE(f, u, v, x, z, n); } while (0);
    if (found) {
      if (!mpz_cmp(f, n)) { found = 0; continue; }
      break;
    }

    /* Stage 1 */
    for (q = 2; q < B1; q = prime_iterator_next(&iter)) {
      UV k = q;
      UV kmin = B1 / q;
      while (k <= kmin)  k *= q;

      NORMALIZE(f, u, v, x, z, n);
      ec_mult(k, x, z);
    }
    prime_iterator_destroy(&iter);

    if (!found)
      do { NORMALIZE(f, u, v, x, z, n); } while (0);

    /* Stage 2 */
    if (!found && B2 > B1)
      found = ec_stage2(B1, B2, x, z, f);

    if (found) {
      if (!mpz_cmp(f, n)) { found = 0; continue; }
      /* gmp_printf("S%d\n", found); */
      break;
    }
  }

  mpz_clear(a);  mpz_clear(b);  mpz_clear(ecn);
  mpz_clear(u);  mpz_clear(v);  mpz_clear(w);
  mpz_clear(x);  mpz_clear(z);
  mpz_clear(x1);  mpz_clear(z1);
  mpz_clear(x2);  mpz_clear(z2);
  mpz_clear(x1pz1);  mpz_clear(x1mz1);
  mpz_clear(sigma);

  return found;
}
