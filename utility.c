/*
 * Utility functions, such as sqrt mod p, polynomial manipulation, etc.
 */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <gmp.h>

#ifdef STANDALONE
  typedef unsigned long UV;
  typedef   signed long IV;
  #define INLINE
  #define UV_MAX ULONG_MAX
  #define UVCONST(x) ((unsigned long)x##UL)
  #define croak(fmt,...)            { printf(fmt,##__VA_ARGS__); exit(1); }
  #define New(id, mem, size, type)  mem = (type*) malloc((size)*sizeof(type))
  #define Newz(id, mem, size, type) mem = (type*) calloc(size, sizeof(type))
  #define Safefree(mem)             free((void*)mem)
  #define PRIME_ITERATOR(i) mpz_t i; mpz_init_set_ui(i, 2)
  /*
  static UV prime_iterator_next(mpz_t *iter) { mpz_nextprime(*iter, *iter); return mpz_get_ui(*iter); }
  static void prime_iterator_destroy(mpz_t *iter) { mpz_clear(*iter); }
  static void prime_iterator_setprime(mpz_t *iter, UV n) {mpz_set_ui(*iter, n);}
  static int prime_iterator_isprime(mpz_t *iter, UV n) {int isp; mpz_t t; mpz_init_set_ui(t, n); isp = mpz_probab_prime_p(t, 10); mpz_clear(t); return isp;}
  */
#else
  #include "EXTERN.h"
  #include "perl.h"
  #include "XSUB.h"
#endif

#include "utility.h"


/* set x to sqrt(a) mod p.  Returns 0 if a is not a square root mod p */
/* See Cohen section 1.5.
 * See http://www.math.vt.edu/people/brown/doc/sqrts.pdf
 */
int sqrtmod(mpz_t x, mpz_t a, mpz_t p,
            mpz_t t, mpz_t q, mpz_t b, mpz_t z) /* 4 temp variables */
{
  int r, e, m;

  if (mpz_kronecker(a, p) != 1) { /* No solution exists */
    mpz_set_ui(x, 0);
    return 0;
  }

  /* Easy cases from page 31 */
  if (mpz_congruent_ui_p(p, 3, 4)) {
    mpz_add_ui(t, p, 1);
    mpz_divexact_ui(t, t, 4);
    mpz_powm(x, a, t, p);
    return 1;
  }

  if (mpz_congruent_ui_p(p, 5, 8)) {
    mpz_sub_ui(t, p, 1);
    mpz_divexact_ui(t, t, 4);
    mpz_powm(q, a, t, p);
    if (mpz_cmp_si(q, 1) == 0) { /* s = a^((p+3)/8) mod p */
      mpz_add_ui(t, p, 3);
      mpz_divexact_ui(t, t, 8);
      mpz_powm(x, a, t, p);
    } else {                      /* s = 2a * (4a)^((p-5)/8) mod p */
      mpz_sub_ui(t, p, 5);
      mpz_divexact_ui(t, t, 8);
      mpz_mul_ui(q, a, 4);
      mpz_powm(x, q, t, p);
      mpz_mul_ui(x, x, 2);
      mpz_mulmod(x, x, a, p, x);
    }
    return 1;
  }

  mpz_sub_ui(q, p, 1);
  e = mpz_scan1(q, 0);              /* Remove 2^e from q */
  mpz_tdiv_q_2exp(q, q, e);
  mpz_set_ui(t, 2);
  while (mpz_legendre(t, p) != -1)  /* choose t "at random" */
    mpz_add_ui(t, t, 1);
  mpz_powm(z, t, q, p);                     /* Step 1 complete */
  r = e;

  mpz_powm(b, a, q, p);
  mpz_add_ui(q, q, 1);
  mpz_divexact_ui(q, q, 2);
  mpz_powm(x, a, q, p);   /* Done with q, will use it for y now */

  while (mpz_cmp_ui(b, 1)) {
    /* calculate how many times b^2 mod p == 1 */
    mpz_set(t, b);
    m = 0;
    do {
      mpz_powm_ui(t, t, 2, p);
      m++;
    } while (m < r && mpz_cmp_ui(t, 1));
    if (m == r) return 0;
    mpz_ui_pow_ui(t, 2, r-m-1);
    mpz_powm(t, z, t, p);
    mpz_mulmod(x, x, t, p, x);
    mpz_powm_ui(z, t, 2, p);
    mpz_mulmod(b, b, z, p, b);
    r = m;
  }
  return 1;
}

/* Smith-Cornacchia: Solve x,y for x^2 + |D|y^2 = p given prime p */
/* See Cohen 1.5.2 */
int cornacchia(mpz_t x, mpz_t y, mpz_t D, mpz_t p)
{
  int result = 0;
  mpz_t a, b, c, d;

  if (mpz_jacobi(D, p) < 0)     /* No solution */
    return 0;

  mpz_init(a); mpz_init(b); mpz_init(c); mpz_init(d);

  sqrtmod(x, D, p, a, b, c, d);
  mpz_set(a, p);
  mpz_set(b, x);
  mpz_sqrt(c, p);

  while (mpz_cmp(b,c) > 0) {
    mpz_set(d, a);
    mpz_set(a, b);
    mpz_mod(b, d, b);
  }

  mpz_mul(a, b, b);
  mpz_sub(a, p, a);   /* a = p - b^2 */
  mpz_abs(d, D);      /* d = |D| */

  if (mpz_divisible_p(a, d)) {
    mpz_divexact(c, a, d);
    if (mpz_perfect_square_p(c)) {
      mpz_set(x, b);
      mpz_sqrt(y, c);
      result = 1;
    }
  }

  mpz_clear(a); mpz_clear(b); mpz_clear(c); mpz_clear(d);

  return result;
}

/* Modified Cornacchia, Solve x,y for x^2 + |D|y^2 = 4p given prime p */
/* See Cohen 1.5.3 */
int modified_cornacchia(mpz_t x, mpz_t y, mpz_t D, mpz_t p)
{
  int result = 0;
  mpz_t a, b, c, d;

  if (mpz_cmp_ui(p, 2) == 0) {
    mpz_add_ui(x, D, 8);
    if (mpz_perfect_square_p(x)) {
      mpz_sqrt(x, x);
      mpz_set_ui(y, 1);
      result = 1;
    }
    return result;
  }
  if (mpz_jacobi(D, p) == -1)     /* No solution */
    return 0;

  mpz_init(a); mpz_init(b); mpz_init(c); mpz_init(d);

  sqrtmod(x, D, p, a, b, c, d);
  if ( (mpz_even_p(D) && mpz_odd_p(x)) || (mpz_odd_p(D) && mpz_even_p(x)) )
    mpz_sub(x, p, x);

  mpz_mul_ui(a, p, 2);
  mpz_set(b, x);
  mpz_sqrt(c, p);
  mpz_mul_ui(c, c, 2);

  /* Euclidean algorithm */
  while (mpz_cmp(b, c) > 0) {
    mpz_set(d, a);
    mpz_set(a, b);
    mpz_mod(b, d, b);
  }

  mpz_mul_ui(c, p, 4);
  mpz_mul(a, b, b);
  mpz_sub(a, c, a);   /* a = 4p - b^2 */
  mpz_abs(d, D);      /* d = |D| */

  if (mpz_divisible_p(a, d)) {
    mpz_divexact(c, a, d);
    if (mpz_perfect_square_p(c)) {
      mpz_set(x, b);
      mpz_sqrt(y, c);
      result = 1;
    }
  }

  mpz_clear(a); mpz_clear(b); mpz_clear(c); mpz_clear(d);

  return result;
}


/* Modular inversion: invert a mod p.
 * This implementation from William Hart, using extended gcd.
 */
unsigned long modinverse(unsigned long a, unsigned long p)
{
  long u1 = 1;
  long u3 = a;
  long v1 = 0;
  long v3 = p;
  long t1 = 0;
  long t3 = 0;
  long quot;
  while (v3) {
    quot = u3 - v3;
    if (u3 < (v3<<2)) {
      if (quot < v3) {
        if (quot < 0) {
          t1 = u1; u1 = v1; v1 = t1;
          t3 = u3; u3 = v3; v3 = t3;
        } else {
          t1 = u1 - v1; u1 = v1; v1 = t1;
          t3 = u3 - v3; u3 = v3; v3 = t3;
        }
      } else if (quot < (v3<<1)) {
        t1 = u1 - (v1<<1); u1 = v1; v1 = t1;
        t3 = u3 - (v3<<1); u3 = v3; v3 = t3;
      } else {
        t1 = u1 - v1*3; u1 = v1; v1 = t1;
        t3 = u3 - v3*3; u3 = v3; v3 = t3;
      }
    } else {
      quot = u3 / v3;
      t1 = u1 - v1*quot; u1 = v1; v1 = t1;
      t3 = u3 - v3*quot; u3 = v3; v3 = t3;
    }
 }
 if (u1 < 0) u1 += p;
 return u1;
}



UV mpz_order_ui(UV r, mpz_t n, UV limit) {
  UV j;
  mpz_t t;

  mpz_init_set_ui(t, 1);
  for (j = 1; j <= limit; j++) {
    mpz_mul(t, t, n);
    mpz_mod_ui(t, t, r);
    if (!mpz_cmp_ui(t, 1))
      break;
  }
  mpz_clear(t);
  return j;
}




void poly_mod_mul(mpz_t* px, mpz_t* py, mpz_t* ptmp, UV r, mpz_t mod)
{
  UV i, j, prindex;

  for (i = 0; i < r; i++)
    mpz_set_ui(ptmp[i], 0);
  for (i = 0; i < r; i++) {
    if (!mpz_sgn(px[i])) continue;
    for (j = 0; j < r; j++) {
      if (!mpz_sgn(py[j])) continue;
      prindex = (i+j) % r;
      mpz_addmul( ptmp[prindex], px[i], py[j] );
    }
  }
  /* Put ptmp into px and mod n */
  for (i = 0; i < r; i++)
    mpz_mod(px[i], ptmp[i], mod);
}

void poly_mod_sqr(mpz_t* px, mpz_t* ptmp, UV r, mpz_t mod)
{
  UV i, d, s;
  UV degree = r-1;

  for (i = 0; i < r; i++)
    mpz_set_ui(ptmp[i], 0);
  for (d = 0; d <= 2*degree; d++) {
    UV prindex = d % r;
    for (s = (d <= degree) ? 0 : d-degree; s <= (d/2); s++) {
      if (s*2 == d) {
        mpz_addmul( ptmp[prindex], px[s], px[s] );
      } else {
        mpz_addmul( ptmp[prindex], px[s], px[d-s] );
        mpz_addmul( ptmp[prindex], px[s], px[d-s] );
      }
    }
  }
  /* Put ptmp into px and mod n */
  for (i = 0; i < r; i++)
    mpz_mod(px[i], ptmp[i], mod);
}

void poly_mod_pow(mpz_t *pres, mpz_t *pn, mpz_t *ptmp, mpz_t power, UV r, mpz_t mod)
{
  UV i;
  mpz_t mpow;

  for (i = 0; i < r; i++)
    mpz_set_ui(pres[i], 0);
  mpz_set_ui(pres[0], 1);

  mpz_init_set(mpow, power);

  while (mpz_cmp_ui(mpow, 0) > 0) {
    if (mpz_odd_p(mpow))            poly_mod_mul(pres, pn, ptmp, r, mod);
    mpz_tdiv_q_2exp(mpow, mpow, 1);
    if (mpz_cmp_ui(mpow, 0) > 0)    poly_mod_sqr(pn, ptmp, r, mod);
  }
  mpz_clear(mpow);
}
