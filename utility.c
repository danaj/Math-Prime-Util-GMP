/*
 * Utility functions, such as sqrt mod p, polynomial manipulation, etc.
 */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <gmp.h>

#include "ptypes.h"

/* includes mpz_mulmod(r, a, b, n, temp) */
#include "utility.h"

static int _verbose = 0;
int get_verbose_level(void) { return _verbose; }
void set_verbose_level(int level) { _verbose = level; }

static gmp_randstate_t _randstate;
gmp_randstate_t* get_randstate(void) { return &_randstate; }
void init_randstate(unsigned long seed) {
  gmp_randinit_mt(_randstate);
  gmp_randseed_ui(_randstate, seed);
}
void clear_randstate(void) {  gmp_randclear(_randstate);  }


int mpz_divmod(mpz_t r, mpz_t a, mpz_t b, mpz_t n, mpz_t t)
{
  int invertible;
  invertible = mpz_invert(t, b, n);
  if (!invertible)
    return 0;
  mpz_mulmod(r, t, a, n, t); /* mpz_mul(t,t,a); mpz_mod(r,t,n); */
  return 1;
}

/* Returns 1 if x^2 = a mod p, otherwise set x to 0 and return 0. */
static int verify_sqrt(mpz_t x, mpz_t a, mpz_t p, mpz_t t, mpz_t t2) {
  mpz_mulmod(t, x, x, p, t2);
  mpz_mod(t2, a, p);
  if (mpz_cmp(t, t2) == 0) return 1;
  mpz_set_ui(x, 0);
  return 0;
}

/* set x to sqrt(a) mod p.  Returns 0 if a is not a square root mod p */
/* See Cohen section 1.5.
 * See http://www.math.vt.edu/people/brown/doc/sqrts.pdf
 */
int sqrtmod(mpz_t x, mpz_t a, mpz_t p,
            mpz_t t, mpz_t q, mpz_t b, mpz_t z) /* 4 temp variables */
{
  int r, e, m;

  /* Easy cases from page 31 (or Menezes 3.36, 3.37) */
  if (mpz_congruent_ui_p(p, 3, 4)) {
    mpz_add_ui(t, p, 1);
    mpz_tdiv_q_2exp(t, t, 2);
    mpz_powm(x, a, t, p);
    return verify_sqrt(x, a, p, t, q);
  }

  if (mpz_congruent_ui_p(p, 5, 8)) {
    mpz_sub_ui(t, p, 1);
    mpz_tdiv_q_2exp(t, t, 2);
    mpz_powm(q, a, t, p);
    if (mpz_cmp_si(q, 1) == 0) { /* s = a^((p+3)/8) mod p */
      mpz_add_ui(t, p, 3);
      mpz_tdiv_q_2exp(t, t, 3);
      mpz_powm(x, a, t, p);
    } else {                      /* s = 2a * (4a)^((p-5)/8) mod p */
      mpz_sub_ui(t, p, 5);
      mpz_tdiv_q_2exp(t, t, 3);
      mpz_mul_ui(q, a, 4);
      mpz_powm(x, q, t, p);
      mpz_mul_ui(x, x, 2);
      mpz_mulmod(x, x, a, p, x);
    }
    return verify_sqrt(x, a, p, t, q);
  }

  if (mpz_kronecker(a, p) != 1) {
    /* Possible no solution exists.  Check Euler criterion. */
    mpz_sub_ui(t, p, 1);
    mpz_tdiv_q_2exp(t, t, 1);
    mpz_powm(x, a, t, p);
    if (mpz_cmp_si(x, 1) != 0) {
      mpz_set_ui(x, 0);
      return 0;
    }
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
    if (m == r) break;
    mpz_ui_pow_ui(t, 2, r-m-1);
    mpz_powm(t, z, t, p);
    mpz_mulmod(x, x, t, p, x);
    mpz_powm_ui(z, t, 2, p);
    mpz_mulmod(b, b, z, p, b);
    r = m;
  }
  return verify_sqrt(x, a, p, t, q);
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

  /* If n < limit, set limit to n */
  if (mpz_cmp_ui(n, limit) < 0)
    limit = mpz_get_ui(n);
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

void mpz_arctan(mpz_t r, unsigned long base, mpz_t pow, mpz_t t1, mpz_t t2)
{
  unsigned long i = 1;
  mpz_tdiv_q_ui(r, pow, base);
  mpz_set(t1, r);
  do {
    mpz_ui_pow_ui(t2, base, 2);
    mpz_tdiv_q(t1, t1, t2);
    mpz_tdiv_q_ui(t2, t1, 2*i+1);
    if (i++ & 1) mpz_sub(r, r, t2); else mpz_add(r, r, t2);
  } while (mpz_sgn(t2));
}



#if 0
/* Simple polynomial multiplication */
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
#endif

#if 0
/* Binary segmentation, using simple shift+add method for processing p.
 * Faster than twiddling bits, but not nearly as fast as import/export.
 */
void poly_mod_mul(mpz_t* px, mpz_t* py, UV r, mpz_t mod, mpz_t p, mpz_t p2, mpz_t t)
{
  UV i, d, bits;
  UV degree = r-1;

  mpz_mul(t, mod, mod);
  mpz_mul_ui(t, t, r);
  bits = mpz_sizeinbase(t, 2);

  mpz_set_ui(p, 0);
  for (i = 0; i < r; i++) {
    mpz_mul_2exp(p, p, bits);
    mpz_add(p, p, px[r-i-1]);
  }

  if (px == py) {
    mpz_mul(p, p, p);
  } else {
    mpz_set_ui(p2, 0);
    for (i = 0; i < r; i++) {
      mpz_mul_2exp(p2, p2, bits);
      mpz_add(p2, p2, py[r-i-1]);
    }
    mpz_mul(p, p, p2);
  }

  for (d = 0; d <= 2*degree; d++) {
    mpz_tdiv_r_2exp(t, p, bits);
    mpz_tdiv_q_2exp(p, p, bits);
    if (d < r)
      mpz_set(px[d], t);
    else
      mpz_add(px[d-r], px[d-r], t);
  }
  for (i = 0; i < r; i++)
    mpz_mod(px[i], px[i], mod);
}
#endif

/* Binary segmentation, using import/export method for processing p.
 * Thanks to Dan Bernstein's 2007 Quartic paper.
 */
void poly_mod_mul(mpz_t* px, mpz_t* py, UV r, mpz_t mod, mpz_t p, mpz_t p2, mpz_t t)
{
  UV i, bytes;
  char* s;

  mpz_mul(t, mod, mod);
  mpz_mul_ui(t, t, r);
  bytes = mpz_sizeinbase(t, 256);
  mpz_set_ui(p, 0);
  mpz_set_ui(p2, 0);

  /* 1. Create big integers from px and py with padding. */
  {
    Newz(0, s, r*bytes, char);
    for (i = 0; i < r; i++)
      mpz_export(s + i*bytes, NULL, -1, 1, 0, 0, px[i]);
    mpz_import(p, r*bytes, -1, 1, 0, 0, s);
    Safefree(s);
  }
  if (px != py) {
    Newz(0, s, r*bytes, char);
    for (i = 0; i < r; i++)
      mpz_export(s + i*bytes, NULL, -1, 1, 0, 0, py[i]);
    mpz_import(p2, r*bytes, -1, 1, 0, 0, s);
    Safefree(s);
  }

  /* 2. Multiply using the awesomeness that is GMP. */
  mpz_mul( p, p, (px == py) ? p : p2 );

  /* 3. Pull out the parts of the result, add+mod, and put in px. */
  {
    Newz(0, s, 2*r*bytes, char);
    /* fill s with data from p */
    mpz_export(s, NULL, -1, 1, 0, 0, p);
    for (i = 0; i < r; i++) {
      /* Set px[i] to the top part, t to the bottom. */
      mpz_import(px[i], bytes, -1, 1, 0, 0, s + (i+r)*bytes);
      mpz_import(t,     bytes, -1, 1, 0, 0, s +     i*bytes);
      /* Add and mod */
      mpz_add(px[i], px[i], t);
      mpz_mod(px[i], px[i], mod);
    }
    Safefree(s);
  }
}

void poly_mod_pow(mpz_t *pres, mpz_t *pn, mpz_t power, UV r, mpz_t mod)
{
  UV i;
  mpz_t mpow, t1, t2, t3;

  for (i = 0; i < r; i++)
    mpz_set_ui(pres[i], 0);
  mpz_set_ui(pres[0], 1);

  mpz_init_set(mpow, power);
  mpz_init(t1);  mpz_init(t2);  mpz_init(t3);

  while (mpz_cmp_ui(mpow, 0) > 0) {
    if (mpz_odd_p(mpow))            poly_mod_mul(pres, pn, r, mod, t1, t2, t3);
    mpz_tdiv_q_2exp(mpow, mpow, 1);
    if (mpz_cmp_ui(mpow, 0) > 0)    poly_mod_mul(pn, pn, r, mod, t1, t2, t3);
  }
  mpz_clear(t1);  mpz_clear(t2);  mpz_clear(t3);
  mpz_clear(mpow);
}

void poly_mod(mpz_t *pres, mpz_t *pn, UV *dn, mpz_t mod)
{
  UV i;
  for (i = 0; i < *dn; i++) {
    mpz_mod(pres[i], pn[i], mod);
  }
  while (*dn > 0 && mpz_sgn(pres[*dn-1]) == 0)
    *dn -= 1;
}
void polyz_mod(mpz_t *pres, mpz_t *pn, long *dn, mpz_t mod)
{
  long i;
  for (i = 0; i <= *dn; i++) {
    mpz_mod(pres[i], pn[i], mod);
  }
  while (*dn > 0 && mpz_sgn(pres[*dn]) == 0)
    *dn -= 1;
}

void polyz_set(mpz_t* pr, long* dr, mpz_t* ps, long ds)
{
  *dr = ds;
  do {
    mpz_set(pr[ds], ps[ds]);
  } while (ds-- > 0);
}

void polyz_print(const char* header, mpz_t* pn, long dn)
{
  gmp_printf("%s", header);
  do { gmp_printf("%Zd ", pn[dn]); } while (dn-- > 0);
  gmp_printf("\n");
}

/* Multiply polys px and py to create new poly pr, all modulo 'mod' */
#if 0
void polyz_mulmod(mpz_t* pr, mpz_t* px, mpz_t *py, long *dr, long dx, long dy, mpz_t mod)
{
  long i, j;
  *dr = dx + dy;
  for (i = 0; i <= *dr; i++)
    mpz_set_ui(pr[i], 0);
  for (i = 0; i <= dx; i++) {
    if (!mpz_sgn(px[i])) continue;
    for (j = 0; j <= dy; j++) {
      if (!mpz_sgn(py[j])) continue;
      mpz_addmul( pr[i+j], px[i], py[j] );
    }
  }
  for (i = 0; i <= *dr; i++)
    mpz_mod(pr[i], pr[i], mod);
  while (*dr > 0 && mpz_sgn(pr[*dr]) == 0)  dr[0]--;
}
#endif
#if 1
void polyz_mulmod(mpz_t* pr, mpz_t* px, mpz_t *py, long *dr, long dx, long dy, mpz_t mod)
{
  UV i, bits, r;
  mpz_t p, p2, t;

  mpz_init(p); mpz_init(t);
  *dr = dx+dy;
  r = *dr+1;
  mpz_mul(t, mod, mod);
  mpz_mul_ui(t, t, r);
  bits = mpz_sizeinbase(t, 2);
  mpz_set_ui(p, 0);

  /* Create big integers p and p2 from px and py, with padding */
  {
    for (i = 0; i <= (UV)dx; i++) {
      mpz_mul_2exp(p, p, bits);
      mpz_add(p, p, px[dx-i]);
    }
  }
  if (px == py) {
    mpz_pow_ui(p, p, 2);
  } else {
    mpz_init_set_ui(p2, 0);
    for (i = 0; i <= (UV)dy; i++) {
      mpz_mul_2exp(p2, p2, bits);
      mpz_add(p2, p2, py[dy-i]);
    }
    mpz_mul(p, p, p2);
    mpz_clear(p2);
  }

  /* Pull out parts of result p to pr */
  for (i = 0; i < r; i++) {
    mpz_tdiv_r_2exp(t, p, bits);
    mpz_tdiv_q_2exp(p, p, bits);
    mpz_mod(pr[i], t, mod);
  }

  mpz_clear(p); mpz_clear(t);
}
#endif
#if 0
void polyz_mulmod(mpz_t* pr, mpz_t* px, mpz_t *py, long *dr, long dx, long dy, mpz_t mod)
{
  UV i, bytes, r;
  char* s;
  mpz_t p, p2, t;

  mpz_init(p); mpz_init(p2); mpz_init(t);
  *dr = dx+dy;
  r = *dr+1;
  mpz_mul(t, mod, mod);
  mpz_mul_ui(t, t, r);
  bytes = mpz_sizeinbase(t, 256);
  mpz_set_ui(p, 0);
  mpz_set_ui(p2, 0);

  /* Create big integers p and p2 from px and py, with padding */
  {
    Newz(0, s, (dx+1)*bytes, char);
    for (i = 0; i <= dx; i++)
      mpz_export(s + i*bytes, NULL, -1, 1, 0, 0, px[i]);
    mpz_import(p, (dx+1)*bytes, -1, 1, 0, 0, s);
    Safefree(s);
  }
  if (px != py) {
    Newz(0, s, (dy+1)*bytes, char);
    for (i = 0; i <= dy; i++)
      mpz_export(s + i*bytes, NULL, -1, 1, 0, 0, py[i]);
    mpz_import(p2, (dy+1)*bytes, -1, 1, 0, 0, s);
    Safefree(s);
  }

  /* Multiply! */
  mpz_mul( p ,p, (px == py) ? p : p2 );

  /* Pull out parts of result p to pr */
  {
    Newz(0, s, r*bytes, char);
    /* fill s with data from p */
    mpz_export(s, NULL, -1, 1, 0, 0, p);
    for (i = 0; i < r; i++) {
      mpz_import(t,     bytes, -1, 1, 0, 0, s +     i*bytes);
      mpz_mod(pr[i], t, mod);
    }
    Safefree(s);
  }

  mpz_clear(p); mpz_clear(p2); mpz_clear(t);
}
#endif

/* Polynomial division modulo N.
 * This is Cohen algorithm 3.1.2 "pseudo-division". */
void polyz_div(mpz_t *pq, mpz_t *pr, mpz_t *pn, mpz_t *pd,
               long *dq, long *dr, long dn, long dd, mpz_t NMOD)
{
  long i, j;
  mpz_t t;

  /* Ensure n and d are reduced */
  while (dn > 0 && mpz_sgn(pn[dn]) == 0)   dn--;
  while (dd > 0 && mpz_sgn(pd[dd]) == 0)   dd--;
  if (dd == 0 && mpz_sgn(pd[0]) == 0)
    croak("polyz_divmod: divide by zero\n");

  /* Q = 0 */
  *dq = 0;
  mpz_set_ui(pq[0], 0);

  /* R = N */
  *dr = dn;
  for (i = 0; i <= dn; i++)
    mpz_set(pr[i], pn[i]);  /* pn should already be mod */

  if (*dr < dd)
    return;
  if (dd == 0) {
    *dq = 0; *dr = 0;
    mpz_tdiv_qr( pq[0], pr[0], pn[0], pd[0] );
    return;
  }

  *dq = dn - dd;
  *dr = dd-1;

  if (mpz_cmp_ui(pd[dd], 1) == 0) {
    for (i = *dq; i >= 0; i--) {
      long di = dd + i;
      mpz_set(pq[i], pr[di]);
      for (j = di-1; j >= i; j--) {
        mpz_submul(pr[j], pr[di], pd[j-i]);
        mpz_mod(pr[j], pr[j], NMOD);
      }
    }
  } else {
    mpz_init(t);
    for (i = *dq; i >= 0; i--) {
      long di = dd + i;
      mpz_powm_ui(t, pd[dd], i, NMOD);
      mpz_mul(t, t, pr[di]);
      mpz_mod(pq[i], t, NMOD);
      for (j = di-1; j >= 0; j--) {
        mpz_mul(pr[j], pr[j], pd[dd]);   /* j != di so this is safe */
        if (j >= i)
          mpz_submul(pr[j], pr[di], pd[j-i]);
        mpz_mod(pr[j], pr[j], NMOD);
      }
    }
    mpz_clear(t);
  }
  /* Reduce R and Q. */
  while (*dr > 0 && mpz_sgn(pr[*dr]) == 0)  dr[0]--;
  while (*dq > 0 && mpz_sgn(pq[*dq]) == 0)  dq[0]--;
}

/* Raise poly pn to the power, modulo poly pmod and coefficient NMOD. */
void polyz_pow_polymod(mpz_t* pres,  mpz_t* pn,  mpz_t* pmod,
                              long *dres,   long   dn,  long   dmod,
                              mpz_t power, mpz_t NMOD)
{
  mpz_t mpow;
  long dProd, dQ, dX, maxd, i;
  mpz_t *pProd, *pQ, *pX;

  /* Product = res*x.  With a prediv this would be dmod+dmod, but without it
   * is max(dmod,dn)+dmod. */
  maxd = (dn > dmod) ? dn+dmod : dmod+dmod;
  New(0, pProd, maxd+1, mpz_t);
  New(0, pQ, maxd+1, mpz_t);
  New(0, pX, maxd+1, mpz_t);
  for (i = 0; i <= maxd; i++) {
    mpz_init(pProd[i]);
    mpz_init(pQ[i]);
    mpz_init(pX[i]);
  }

  *dres = 0;
  mpz_set_ui(pres[0], 1);

  dX = dn;
  for (i = 0; i <= dX; i++)
    mpz_set(pX[i], pn[i]);

  mpz_init_set(mpow, power);
  while (mpz_cmp_ui(mpow, 0) > 0) {
    if (mpz_odd_p(mpow)) {
      polyz_mulmod(pProd, pres, pX, &dProd, *dres, dX, NMOD);
      polyz_div(pQ, pres,  pProd, pmod,  &dQ, dres, dProd, dmod, NMOD);
    }
    mpz_tdiv_q_2exp(mpow, mpow, 1);
    if (mpz_cmp_ui(mpow, 0) > 0) {
      polyz_mulmod(pProd, pX, pX, &dProd, dX, dX, NMOD);
      polyz_div(pQ, pX,  pProd, pmod,  &dQ, &dX, dProd, dmod, NMOD);
    }
  }
  mpz_clear(mpow);
  for (i = 0; i <= maxd; i++) {
    mpz_clear(pProd[i]);
    mpz_clear(pQ[i]);
    mpz_clear(pX[i]);
  }
  Safefree(pProd);
  Safefree(pQ);
  Safefree(pX);
}

void polyz_gcd(mpz_t* pres, mpz_t* pa, mpz_t* pb, long* dres, long da, long db, mpz_t MODN)
{
  long i;
  long dr1, dq, dr, maxd;
  mpz_t *pr1, *pq, *pr;

  /* Early reduction so we're not wasteful with messy input. */
  while (da > 0 && mpz_sgn(pa[da]) == 0)   da--;
  while (db > 0 && mpz_sgn(pb[db]) == 0)   db--;

  /* Swap a/b so degree a >= degree b */
  if (db > da) {
    mpz_t* mtmp;
    long ltmp;
    mtmp = pa; pa = pb; pb = mtmp;
    ltmp = da; da = db; db = ltmp;
  }
  /* TODO: Should set c=pa[da], then after Euclid loop, res = c^-1 * res */

  /* Allocate temporary polys */
  maxd = da;
  New(0, pr1, maxd+1, mpz_t);
  New(0, pq, maxd+1, mpz_t);
  New(0, pr, maxd+1, mpz_t);
  for (i = 0; i <= maxd; i++) {
    mpz_init(pr1[i]);
    mpz_init(pq[i]);
    mpz_init(pr[i]);
  }

  /* R0 = a (we use pres) */
  *dres = da;
  for (i = 0; i <= da; i++)
    mpz_mod(pres[i], pa[i], MODN);
  while (*dres > 0 && mpz_sgn(pres[*dres]) == 0)   dres[0]--;

  /* R1 = b */
  dr1 = db;
  for (i = 0; i <= db; i++)
    mpz_mod(pr1[i], pb[i], MODN);
  while (dr1 > 0 && mpz_sgn(pr1[dr1]) == 0)   dr1--;

  while (dr1 > 0 || mpz_sgn(pr1[dr1]) != 0) {
    polyz_div(pq, pr,  pres, pr1,  &dq, &dr,  *dres, dr1, MODN);
    if (dr < 0 || dq < 0 || dr > maxd || dq > maxd)
      croak("division error: dq %ld dr %ld maxd %ld\n", dq, dr, maxd);
    /* pr0 = pr1.  pr1 = pr */
    *dres = dr1;
    for (i = 0; i <= dr1; i++)
      mpz_set(pres[i], pr1[i]);  /* pr1 is mod MODN already */
    dr1 = dr;
    for (i = 0; i <= dr; i++)
      mpz_set(pr1[i], pr[i]);
  }
  /* return pr0 */
  while (*dres > 0 && mpz_sgn(pres[*dres]) == 0)   dres[0]--;

  for (i = 0; i <= maxd; i++) {
    mpz_clear(pr1[i]);
    mpz_clear(pq[i]);
    mpz_clear(pr[i]);
  }
  Safefree(pr1);
  Safefree(pq);
  Safefree(pr);
}



void polyz_root_deg1(mpz_t root, mpz_t* pn, mpz_t NMOD)
{
  mpz_invert(root, pn[1], NMOD);
  mpz_mul(root, root, pn[0]);
  mpz_neg(root, root);
  mpz_mod(root, root, NMOD);
}
void polyz_root_deg2(mpz_t root1, mpz_t root2, mpz_t* pn, mpz_t NMOD)
{
  mpz_t e, d, t, t2, t3, t4;

  mpz_init(e); mpz_init(d);
  mpz_init(t); mpz_init(t2); mpz_init(t3); mpz_init(t4);

  mpz_mul(t, pn[0], pn[2]);
  mpz_mul_ui(t, t, 4);
  mpz_mul(d, pn[1], pn[1]);
  mpz_sub(d, d, t);
  sqrtmod(e, d, NMOD, t, t2, t3, t4);

  mpz_neg(t4, pn[1]);                    /* t4 = -a_1      */
  mpz_mul_ui(t3, pn[2], 2);              /* t3 = 2a_2      */
  mpz_add(t, t4, e);                     /* t = -a_1 + e   */
  mpz_divmod( root1, t, t3, NMOD, t2);   /* (-a_1+e)/2a_2  */
  mpz_sub(t, t4, e);                     /* t = -a_1 - e   */
  mpz_divmod( root2, t, t3, NMOD, t2);   /* (-a_1-e)/2a_2  */
  mpz_clear(e); mpz_clear(d);
  mpz_clear(t); mpz_clear(t2); mpz_clear(t3); mpz_clear(t4);
}

/* Algorithm 2.3.10 procedure "roots" from Crandall & Pomerance.
 * Step 3/4 of Cohen Algorithm 1.6.1.
 * Uses some hints from Pate Williams (1997-1998) for the poly math */
static void polyz_roots(mpz_t* roots, long *nroots,
                        long maxroots, mpz_t* pg, long dg, mpz_t NMOD,
                        gmp_randstate_t* p_randstate)
{
  long i, ntries, maxtries, maxd, dxa, dt, dh, dq, dup;
  mpz_t t, power;
  mpz_t pxa[2];
  mpz_t *pt, *ph, *pq;

  if (*nroots >= maxroots || dg <= 0)  return;

  mpz_init(t);
  mpz_init(pxa[0]);  mpz_init(pxa[1]);

  if (dg <= 2) {
    if (dg == 1) polyz_root_deg1( pxa[0], pg, NMOD );
    else         polyz_root_deg2( pxa[0], pxa[1], pg, NMOD );
    for (dt = 0; dt < dg; dt++) {
      mpz_set(t, pxa[dt]);
      dup = 0; /* don't add duplicate roots */
      for (i = 0; i < *nroots; i++)
        if (mpz_cmp(t, roots[i]) == 0)
          { dup = 1; break; }
      if (!dup)
        mpz_set(roots[ (*nroots)++ ], t);
    }
    mpz_clear(t);
    mpz_clear(pxa[0]);  mpz_clear(pxa[1]);
    return;
  }

  /* If not a monic polynomial, divide by leading coefficient */
  if (mpz_cmp_ui(pg[dg], 1) != 0) {
    if (!mpz_invert(t, pg[dg], NMOD)) {
      mpz_clear(t);
      return;
    }
    for (i = 0; i <= dg; i++)
      mpz_mulmod(pg[i], pg[i], t, NMOD, pg[i]);
  }

  /* Try hard to find a single root, work less after we got one or two.
   * In a generic routine we would want to try hard all the time. */
  ntries = 0;
  maxtries = (*nroots == 0) ? 200 : (*nroots == 1) ? 50 : 10;

  mpz_init(power);
  mpz_set_ui(pxa[1], 1);
  dxa = 1;

  maxd = 2 * dg;
  New(0, pt, maxd+1, mpz_t);
  New(0, ph, maxd+1, mpz_t);
  New(0, pq, maxd+1, mpz_t);
  for (i = 0; i <= maxd; i++) {
    mpz_init(pt[i]);
    mpz_init(ph[i]);
    mpz_init(pq[i]);
  }

  mpz_sub_ui(t, NMOD, 1);
  mpz_tdiv_q_2exp(power, t, 1);
  /* We'll pick random "a" values from 1 to 1000M */
  mpz_set_ui(t, 1000000000UL);
  if (mpz_cmp(t, NMOD) > 0) mpz_set(t, NMOD);

  while (ntries++ < maxtries) {
    /* pxa = X+a for randomly selected a */
    if (ntries <= 2)  mpz_set_ui(pxa[0], ntries);  /* Simple small values */
    else              mpz_urandomm(pxa[0], *p_randstate, t);

    /* Raise pxa to (NMOD-1)/2, all modulo NMOD and g(x) */
    polyz_pow_polymod(pt, pxa, pg, &dt, dxa, dg, power, NMOD);

    /* Subtract 1 and gcd */
    mpz_sub_ui(pt[0], pt[0], 1);
    polyz_gcd(ph, pt, pg, &dh, dt, dg, NMOD);

    if (dh >= 1 && dh < dg)
      break;
  }

  if (dh >= 1 && dh < dg) {
    /* Pick the smaller of the two splits to process first */
    if (dh <= 2 || dh <= (dg-dh)) {
      polyz_roots(roots, nroots, maxroots, ph, dh, NMOD, p_randstate);
      if (*nroots < maxroots) {
        /* q = g/h, and recurse */
        polyz_div(pq, pt,  pg, ph,  &dq, &dt, dg, dh, NMOD);
        polyz_roots(roots, nroots, maxroots, pq, dq, NMOD, p_randstate);
      }
    } else {
      polyz_div(pq, pt,  pg, ph,  &dq, &dt, dg, dh, NMOD);
      polyz_roots(roots, nroots, maxroots, pq, dq, NMOD, p_randstate);
      if (*nroots < maxroots) {
        polyz_roots(roots, nroots, maxroots, ph, dh, NMOD, p_randstate);
      }
    }
  }

  mpz_clear(t);
  mpz_clear(power);
  mpz_clear(pxa[0]);  mpz_clear(pxa[1]);

  for (i = 0; i <= maxd; i++) {
    mpz_clear(pt[i]);
    mpz_clear(ph[i]);
    mpz_clear(pq[i]);
  }
  Safefree(pt);
  Safefree(ph);
  Safefree(pq);
}


/* Algorithm 1.6.1 from Cohen, minus step 1. */
void polyz_roots_modp(mpz_t** roots, long *nroots, long maxroots,
                      mpz_t *pP, long dP, mpz_t NMOD,
                      gmp_randstate_t* p_randstate)
{
  long i;

  *nroots = 0;
  *roots = 0;

  if (dP == 0)
    return;

  /* Allocate space for the maximum number of roots (plus 1 for safety) */
  New(0, *roots, dP+1, mpz_t);
  for (i = 0; i <= dP; i++)
    mpz_init((*roots)[i]);

  if (maxroots > dP || maxroots == 0)
    maxroots = dP;

  /* Do degree 1 or 2 explicitly */
  if (dP == 1) {
    polyz_root_deg1( (*roots)[0], pP, NMOD );
    *nroots = 1;
    return;
  }
  if (dP == 2) {
    polyz_root_deg2( (*roots)[0], (*roots)[1], pP, NMOD );
    *nroots = 2;
    return;
  }

  polyz_roots(*roots, nroots, maxroots, pP, dP, NMOD, p_randstate);
  /* This could be just really bad luck.  Let the caller handle it. */
  /* if (*nroots == 0) croak("failed to find roots\n"); */

  /* Clear out space for roots we didn't find */
  for (i = *nroots; i <= dP; i++)
    mpz_clear((*roots)[i]);
}


#include "class_poly_data.h"

int* poly_class_nums(void)
{
  int* dlist;
  UV i;
  int degree_offset[256] = {0};

  for (i = 1; i < NUM_CLASS_POLYS; i++)
    if (_class_poly_data[i].D < _class_poly_data[i-1].D)
      croak("Problem with data file, out of order at D=%d\n", (int)_class_poly_data[i].D);

  Newz(0, dlist, NUM_CLASS_POLYS + 1, int);
  /* init degree_offset to total number of this degree */
  for (i = 0; i < NUM_CLASS_POLYS; i++)
    degree_offset[_class_poly_data[i].degree]++;
  /* set degree_offset to sum of this and all previous degrees. */
  for (i = 1; i < 256; i++)
    degree_offset[i] += degree_offset[i-1];
  /* Fill in dlist, sorted */
  for (i = 0; i < NUM_CLASS_POLYS; i++) {
    int position = degree_offset[_class_poly_data[i].degree-1]++;
    dlist[position] = i+1;
  }
  /* Null terminate */
  dlist[NUM_CLASS_POLYS] = 0;
  return dlist;
}

UV poly_class_poly_num(int i, int *D, mpz_t**T, int* type)
{
  UV degree, j;
  int ctype;
  mpz_t t;
  const char* s;

  if (i < 1 || i > (int)NUM_CLASS_POLYS) { /* Invalid number */
     if (D != 0) *D = 0;
     if (T != 0) *T = 0;
     return 0;
  }
  i--; /* i now is the index into our table */

  degree = _class_poly_data[i].degree;
  ctype  = _class_poly_data[i].type;
  s = _class_poly_data[i].coefs;

  if (D != 0)  *D = -_class_poly_data[i].D;
  if (type != 0)  *type = ctype;
  if (T == 0) return degree;

  New(0, *T, degree+1, mpz_t);
  mpz_init(t);
  for (j = 0; j < degree; j++) {
    unsigned char signcount = (*s++) & 0xFF;
    unsigned char sign = signcount >> 7;
    unsigned long count = signcount & 0x7F;
    if (count == 127) {
      do {
        signcount = (*s++) & 0xFF;
        count += signcount;
      } while (signcount == 127);
    }
    mpz_set_ui(t, 0);
    while (count-- > 0) {
      mpz_mul_2exp(t, t, 8);
      mpz_add_ui(t, t, (unsigned long) (*s++) & 0xFF);
    }
    /* Cube the last coefficient of Hilbert polys */
    if (j == 0 && ctype == 1) mpz_pow_ui(t, t, 3);
    if (sign) mpz_neg(t, t);
    mpz_init_set( (*T)[j], t );
  }
  mpz_clear(t);
  mpz_init_set_ui( (*T)[degree], 1 );
  return degree;
}
