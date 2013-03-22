
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
