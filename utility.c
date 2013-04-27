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

/* includes mpz_mulmod(r, a, b, n, temp) */
#include "utility.h"

int mpz_divmod(mpz_t r, mpz_t a, mpz_t b, mpz_t n, mpz_t t)
{
  int invertible;
  invertible = mpz_invert(t, b, n);
  if (!invertible)
    return 0;
  mpz_mulmod(r, t, a, n, t); /* mpz_mul(t,t,a); mpz_mod(r,t,n); */
  return 1;
}

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

/* Ex:
 *   -4    x - 1728
 *   -7    x + 3375
 *   -23   x^3 + 3491750*x^2 - 5151296875*x + 12771880859375
 *   -148  x^2 - 39660183801072000*x - 7898242515936467904000000
 *   -163  x + 262537412640768000
 *
 * Some good ways to get these:
 *   1. Cohen 7.6.2 (page 415)
 *   2. NZMATH:  http://hilbert-class-polynomial.appspot.com
 *   3. SAGE
 * Missing values are usually because of duplication from earlier values
 */
struct _hilbert_poly {
  unsigned short  D;
  unsigned short  degree;
  const char *coefs;
};

/* Storage,
 *   1. The first term (x^degree) is not stored since the coefficient is 1.
 *   2. Numbers are in base 58 for a little more efficient storage.
 *   3. Coefficients are separated by nulls 
 *   4. The final term (x^0) should be cubed.
 *
 * Split into two arrays so we can disable the big one to save space.
 */
static const struct _hilbert_poly _hilbert_data1[] =
{
  { 3, 1, "0"},
  { 4, 1, "-C"},
  { 7, 1, "F"},
  { 8, 1, "-K"},
  { 11, 1, "W"},
  { 15, 2, "-8V\x00" "ujV"},
  { 19, 1, "1c"},
  { 20, 2, "-FA\x00" "-6Rh6"},
  { 23, 3, "6t1\x00" "-7nBh5L\x00" "HpuY"},
  { 24, 2, "gC\x00" "-OjF6"},
  { 31, 3, "YUB\x00" "-1VNXnn0\x00" "3SNMZ"},
  { 35, 2, "-1UG\x00" "AOYlo"},
  { 39, 4, "E74T\x00" "Ejtqjjo5Dj\x00" "-BGss2JC\x00" "THAk0"}, /* SAGE */
  { 40, 2, "6C0\x00" "-bZjb2"},
  { 43, 1, "GW"},
  { 47, 5, "MGb3T\x00" "-AJZrA0OkAMQNq\x00" "BocU01cGqJt\x00" "-4ULG9Kcf\x00" "3PTvOX"},
  { 51, 2, "5Rk\x00" "8PbYPs"},
  { 52, 2, "-OZY\x00" "-ATQI88"},
  { 55, 4, "-2Jh3T\x00" "6qYJaDWKBIoj\x00" "-9SGADmpf\x00" "K0mscT"},
  { 56, 4, "1qc24\x00" "5Du6pXTKXDQ\x00" "G4lPRoiLk\x00" "-OfJfBs"},
  { 59, 3, "3eHE\x00" "-15ip6mG2C\x00" "k0Qva4"},
  { 67, 1, "1X2"},
  { 68, 4, "-BHCLs\x00" "-ChJCUYM0Oa8W\x00" "-ACDhAfc5WC\x00" "-4dTq9Ns"},
  { 71, 7, "NgbHFId\x00" "-7f7oZ1NsODjsjUM3KR5\x00" "5NJarvPWDgsdbrvhne\x00" "-2qrWEVPTmTbNGfm9\x00" "19poArmToUCnJa\x00" "-O8LJYUSN4\x00" "8DnkhuN"},
  { 79, 5, "Qln89T\x00" "-KfpJMZZ5qoLp7UPq\x00" "LJU08d2sFJ1Gmt\x00" "-nfSq3f0KC\x00" "ZAJcEU4"},
  { 83, 3, "fvBM\x00" "-5XurMpAFMe\x00" "1CfGc048"},
  { 84, 4, "-FE73s\x00" "13GhPRbcHIoQ3Q\x00" "-D8TeTGLeVqC\x00" "-1PuUKQbs"},
  { 87, 6, "LU86us5\x00" "7mKFR07ZhmVDDCQLruT\x00" "UJf9U40p39vXM8QXp\x00" "5ksD96BtkpuvB8q\x00" "BUHelecuva\x00" "2Nk1nlVf"},
  { 88, 2, "Cnq4\x00" "-2nKXdrnA"},
  { 91, 2, "-kXE\x00" "4e6dksHk"},
  { 95, 8, "3fVLhjbaT\x00" "10Ga3dZ43XFmeh3UfEvgQuJ5B\x00" "-2HV1BjMGNhYTivh40CBUV3n\x00" "1q3D78YJklMMHASPjo2mh\x00" "-Di8dRLIUktoVVMV8fr\x00" "1N1Qim3vW9oRcvq9\x00" "-1YbdjZiFrik\x00" "903vWH3A"},
  { 103, 5, "838kCLD\x00" "E2KUL9h0ammE7VNUfd\x00" "HRStVfdoeSFG3Oam\x00" "BTPqi8ZZkR\x00" "VmQaaiYv"},
  { 104, 6, "AXtEie0\x00" "-R3BIJDCSIVHsrTqn0C\x00" "1Q2X8iikb9ehTu59KC\x00" "6Ku5shcjY8CdFHo\x00" "TYcMnDSqe9dc\x00" "-b8h3v6iu"},
  { 107, 3, "68p9g\x00" "-FefkejoHLDo\x00" "10jAlvoG0"},
  { 111, 8, "NXQCHjXVl\x00" "-UmjCXuQTTqH5BotolS83sHO8n\x00" "2QRe7J2iBNp1UQcRnIoVIE0S\x00" "-2MJnc2U1MT4bJYnqOsclUY\x00" "s8jE0RnLXXFQTIUXWnI\x00" "3OTIrXRGOYC1hIL07\x00" "SQGQR70LO02\x00" "1nHOAGMdT"},
  { 115, 2, "2YdI\x00" "3JjG6VK9k"},
  { 116, 6, "-26CrM7A0\x00" "iklJdQupkaaqvshAm000\x00" "-1BsUhl2JNORU0i1hXmq\x00" "-2tHAPj8nK0281TQW\x00" "-7aTYeGUoCMQ6C\x00" "-3oG7pNQXs"},
  { 120, 4, "PsOOH6\x00" "-98sjoVDOsDZje4RY\x00" "I9goh3tpAonbg\x00" "-6pseSdX6e"},
  { 123, 2, "R9FU\x00" "AXH8BE7jo"},
  { 127, 5, "35a5EtM7\x00" "-K62lFUJN30ZNf4PAk5p6\x00" "5s3PepgPfV4dLr2eqa\x00" "-1D3bZGBM2LLP\x00" "IVmGdqC78"},
  { 131, 5, "1LtNFSa\x00" "-3dPtKcJ9r0BLh7vJ6\x00" "LvOJUv4DD9RYsQa\x00" "-QntoLOVb6bts\x00" "WEenhFL68"},
  { 132, 4, "I1YW3o\x00" "3KQOvanr0dDCHhXWK\x00" "-3oNODN7vQt3OTc\x00" "-avJLWgRIm"},
  { 136, 4, "KQnFQu\x00" "-6SIEYD9K5kLeE0nM\x00" "8hm7rkIve90naa\x00" "-15bgGYCTQ8"},
  { 139, 3, "3YY3A\x00" "-2776mZ1lrLoW\x00" "1b7iuZaSC8"},
  { 148, 2, "-HYdp6\x00" "-5Je9j2Xgho"},
  { 151, 7, "1vpjRYh9E7\x00" "2D8f61QRkqN4iLKM391p0YNbup\x00" "LULAIWnbZ7Gfga8vrmDNP0pf\x00" "6R2TZvD6oOhsNOIcKK9gPW\x00" "13aL64SnFK8J405aPs74\x00" "iHbveCfE0eLJ\x00" "7nIKcX6bbG"},
  { 152, 6, "ZFl7JGpg\x00" "-2CELUvC4aSEWaofkFgBteiW\x00" "Evh2IUXurN1dD32r5lrGa\x00" "7BtQMh6aDLrjhYJoEK\x00" "3SV0n8D2kZEjHRs\x00" "-8rH349ZARo"},
  { 155, 4, "otll60\x00" "1E7dPGqFP3RUopWb6\x00" "-Ue3lPS270HsdA\x00" "D2eeGn7VOC"},
  { 163, 1, "3GK0"},
  { 164, 8, "-CiQEj5nA00\x00" "1hMAHYBpoOHjV7PVYCnbm4000000\x00" "S1rFHRkWi9L58WOCYMeebvlIO\x00" "-2GlWc350NgLhfBVgKrjMjKAa\x00" "-6e2WllAEvU0gtPFqRk0W5U\x00" "-1tTX0Af8Mv97l1OuDDQ\x00" "-XCVtsD2335gPoiO\x00" "-du1ZfckCm0"},
  { 168, 4, "KkVJ5gW\x00" "-4kFWDuF3h9cRds5piFY\x00" "1B1R48gXEDJ1DWRg\x00" "-174uqMmNYUK"},
  { 179, 5, "1ns1DPhc\x00" "-7IHVPHmSnvr13nXgT9IW\x00" "2ktcYPFh3PjWVcl5TY\x00" "-QADuR4FptntZhE\x00" "49e2mjsiArs"},
  { 183, 8, "1ESCgpX4EZQn\x00" "-5kg8Q7iWDW9rg1a5jq6718CQle7GS7pf\x00" "DDM6jJVHspS7OpnavRKVi7RSus8qjf\x00" "1Ym2O7uUGWUL2nn2LN9NkvdoSdlh\x00" "5huBNPW0S9pDWu80DO9jNQmEU\x00" "iB0caVWnmobaNhBiNrIie\x00" "oAt25HonUc0rmP\x00" "6bWJiCdRakr"},
  { 184, 4, "Chv7iA4\x00" "ee3h5JqplhpKaiOmlQ\x00" "KMroR71hhP8Va4Dk\x00" "-7Qtiia8sGns"},
  { 187, 2, "-1MHI8\x00" "AVst1CXO5s0"},
  { 195, 4, "-1ZlNOi0\x00" "6MV3ZsjJqJPfC53cG\x00" "5BXTa59uoTqSFtk\x00" "QBE2FFmQr60"},
  { 203, 4, "8JFLMcq\x00" "4VSOGETrtr94R5hUtvQ\x00" "-H24gbf53k9PmENY\x00" "15meh7TeBe7o"},
  { 211, 3, "2cA2K0\x00" "3HNW7SKtX8cK64\x00" "2aqdq88Q2iPs"},
  { 212, 6, "-5Rnrj1uqUq\x00" "5bh9m6vFQbfpWPeru8AOO8FSk9k\x00" "-5E0BZu2Tuj00iASuupusW63se\x00" "-E9tmpkvkFgUttLdf3mfl6\x00" "-cKPfVfkRANFIocllk\x00" "-2sKE5kBLaSuG"},
  { 219, 4, "Ca9R3no\x00" "-GkL3kW4Ulct6JlDH3M\x00" "2sQmjBmX0QeMikuO\x00" "6CGUuU1ufm9k"},
  { 223, 7, "VskAFuDSFu3\x00" "-bhFYY00NWYSuAFNRVH74p0V5Ro1KcgO\x00" "c4OuPV92igN2jcBc2pl5tm6WCnI7e\x00" "-1P3jvPQIhUTs3vMTbnL47pPXDZr\x00" "1Vneq5F7M4OC9hKn5bcEnSi5\x00" "-8mPAIbPLFndgu81\x00" "9RkJeuqfJhMA"},
  { 227, 5, "7jg8hkWG\x00" "-BLhZTYku5n2OqnOT94Fvg\x00" "5eMGnioA5vXBQrST62EC\x00" "-93XqbGS0c8dII1K4\x00" "ENmbro5OKQFo"},
  { 228, 4, "MJghYEHk\x00" "5QuLi4T8DiJXAC4Pgpa6lo\x00" "-8Mq7rfNsESXkBWqLiK\x00" "-FvXa0mQt1bHQ"},
  { 232, 2, "bR2aEq\x00" "-OBftZLKGhClo"},
  { 235, 2, "KBbcu\x00" "WsjthVrmHvD2"},
  { 244, 6, "-GfeVU16aO\x00" "-BaHgiJUI14PuI0JAY9Saler60\x00" "-o74CrA61RJABNRRlOaSY78a\x00" "-ET3N47YArXT10kONnau96\x00" "-1dh4OWNKLGOOIMIFE1M\x00" "-1O7oUoRAS8vJk"},
  { 247, 6, "-HC0EXjk25JT\x00" "50OJOemJ58mDr0fu1apiqoGITtvZsgH\x00" "62ILCt0avrWc8gSbtGNcsCr5KOdf\x00" "3cF9h0tSmItscF1NMQetsoffo\x00" "38FvjU0kEeAPOf59\x00" "1qtOv3ctBLqe2"},
  { 248, 8, "4MnQP4rANclg\x00" "-43jjYSPOWdie1bFs3qFi66pEn1XCNqsbs\x00" "3NHsaiIUMAhuFoYOs7riPdFIrkbMtSi\x00" "-19fkqT8GYPk5l0l5PhKKGHUsLXu40\x00" "53Xnm0nSpCvc0k0QgbFMb9Yi3GK\x00" "3sdaOThqDC6NXW9tnARQACW\x00" "34CMDusgo6MXtqNNE5U\x00" "-26Z9BoF900hT2"},
  { 251, 7, "GjSBabcem0\x00" "-4La7u8KqUlg9FrJJbq5ZVnrlJ5bY\x00" "3O7iVpcWoJavD0TSYYP5Tj7tFM\x00" "Ci51XDbftsHrYv6nTCnZfSS\x00" "5fpfLI8BQvjqAeNbkK4OW\x00" "-427keG5IYg36qAEGO\x00" "2nD6cJY48nOi0"},
  { 259, 4, "4HKolu4\x00" "5itO38hCf0fYtvYTTQ\x00" "-1HUJT8j7aAI1Eh4C\x00" "6EtIplr7N4ko4"},
  { 260, 8, "-5te2OU3K11MO\x00" "JXjQRIOhPje9WZjFLX8MLcieFCoZR90ANo\x00" "-IVLRltlo3N21joE8XWniMB7lab4Lt90K\x00" "AjQT1Hfhog4TRGnI2ICsYPO09aVrCS\x00" "-17WlFPQSYe06sF19oRgTZIhisMue\x00" "-YPvKdjldfClKfQKqGntPnii\x00" "-I7Af9XtVDG6Hs5ra8CO\x00" "-6q7R21Gej3DNk"},
  { 264, 8, "2i3gJhCrVvI0\x00" "-b264i4Wnia65WfIc8OMKqmbiJ7OrBTTI\x00" "2LcRE7NgfNPLm771Ud4H0SaakmbI4Pc\x00" "-kbKvTRBEGTv3u9rZ2lUeAmhFvHCi\x00" "-CshICnKgrD6W4PHH5eXkec4JsO\x00" "-1CAFRJP4cUEHhd4WOf2ifosK\x00" "WPeK57HuH3hSAjOIUUC\x00" "-A9qaY5m6ZK4um"},
  { 267, 2, "CJiIdI\x00" "DXhAWVqTa4jrg"},
  { 276, 8, "9395iR1uDM00\x00" "-OEhIM1O2BsAv8NLjtMaRgv6OQsLpmC000\x00" "9YfRgdgsdc96QMTgAaCZfS6ma04jQiO\x00" "H81F2C4TfFfNJDtMUmU2eAfmYEsLM\x00" "-14K64OkZ7NQJA0KcLmOLZtBB97U\x00" "9s6MkqqTD7h0PSP1qlfChfAW\x00" "-37FBr6CZ6IYLPL0uMU0O\x00" "-W1oDlZrIvg7Q0"},
  { 280, 4, "sma1OCM4\x00" "-17lkqPPlRWOWFREZI9QAXkW\x00" "5T3O9tSTjGPDcZ2ORu1k\x00" "-ke9YeIBdIKFJs"},
  { 283, 3, "8rt4oq\x00" "IafDePEn7j86DJE\x00" "13mLN29jn0qr0C"},
  { 291, 4, "CbZIWDmC\x00" "1VApsBTVdE3bukZsffA4u\x00" "BJfJLu3Mt27umfZ1jY\x00" "2DmEEJZdKiNZ4K"},
  { 292, 4, "-Wl53gO00\x00" "4DG4TOnStsgQgYalc1M000\x00" "-TGFUT0fT3auilPYTUDts\x00" "-2QJrHb3OdSg0U4"},
  { 295, 8, "16bIHTZtKg6pk9\x00" "-6KfbcSLGDl14rOnRti7Mcl5kWtOjdioCa3p4nV\x00" "YUUgpD1NX9TLROpSs8KCMUajKtugZmAQedsO\x00" "2ds0oFJNL7T2g6u5j5C1vZAbtbLKUf3dLP\x00" "v4hiRqNNlGDJerkD8aTLKpqBcE518T\x00" "APXG3OsnlXj9WgSpeV5a8HttAUb\x00" "HbAMHWCSDK4jLk0sp\x00" "3DNq868R6USfGs"},
  { 299, 8, "-1389Ghrl5ojk\x00" "9okJj16sLTPnYUI8F5LZqmkUucCJ82NI\x00" "-Ds8vaGdKkLnQfJTM8ICXh1fS7WLVfo\x00" "FcFTlpoipfRDfqa6XgcO9jEHQZFE\x00" "-1UlajfUCiH3GIUcKpR8ab2rXqC\x00" "3KV5YD6t9dJBdMNKANKTAZ2\x00" "-U629reWUJnb74TjbaC\x00" "4bnlESoAY1hTjo"},
  { 307, 3, "VdElse\x00" "-Huht4Y350B17psG0\x00" "9XRpHcljImVkP6"},
  { 308, 8, "-8rhuhJDNjtk00\x00" "1A9srbBkoaAHEresIPYOo1Guf7YkVGK000000\x00" "-1H3vUNWQjOr8hcM5VZT7OA6BC5FZmuokWpY\x00" "PQ0vGHV4ZDuh5SCYfY72odH9naBZVD3M\x00" "-5TF00dFqCQKujNb9L0aKKgecUHcA1U\x00" "-beUNmvnq4FP0BsbrGBvdSn5PQ\x00" "-4QCCiO1U6jJOMrkoOW52O\x00" "-ARV4YRElBQPNc0"},
  { 312, 4, "1ZA88JL168\x00" "-1EUvQFWdITnpVPG65pQPO6Wojc\x00" "7Z0SeVfPo9TphHhYNfe8u\x00" "-Eu48DeLQ7fQVOq"},
  { 323, 4, "-3oZIk7aUi\x00" "22CMQsiLI5oVtDDQshnB7Mca\x00" "-8uDF3jSHki8YME8HaK0\x00" "dRJXCXvE465JpU"},
  { 328, 4, "3SAAIr6vM\x00" "1Tj4Z0HQbW0DVUKdgfdfsslY\x00" "14WTTuH3qpc5HHb0EaYQbc\x00" "-13067FZpD7iT1O4"},
  { 331, 3, "A3TFLg0\x00" "MS1q2iSL8DbiGogiG\x00" "1L4oIO79o4s9ee8"},
  { 339, 6, "A4Nb8uFWhg\x00" "G8GPsHlSFs5VvH1XmcAYN6iRYpo\x00" "FsuYrIklrapj68LdJrit5Ddaa\x00" "-n2W3C3jSJ1XAFUms56fs8m\x00" "uVo0uiTHXOZSRN2blNc\x00" "2f7hD8RFi9apcRo"},
  { 340, 4, "2gfbACNm4\x00" "2l6RsMcuLPJEUMYpgVTMKmMpQ\x00" "-53rCUKlG0gM5HXS3MIOqAC\x00" "-2t6jKmn3JcMlS5s"},
  { 347, 5, "iSDiHTqNY\x00" "-aYZUeWhCDXkqRAeCJbivYmHhU\x00" "158T052DV7XEs6Zd4KlFHoG7g\x00" "-2Nn6HjHOqJmVBHhLqo8O\x00" "5Jf1SJmK7PWtJNc"},
  { 355, 4, "2SkFR6Ii\x00" "-7TFU95gEeEAqhFUk6j3k\x00" "7AMokgX8COdpO3SJom\x00" "APeSpnX4IT528tg"},
  { 371, 8, "1GZteMcNh1SJA\x00" "-AOt3mVZG4dklPGIuEIgNSjvQDCXIO4HTXEO\x00" "TXRatnOVCfbgEMHDrA7aJ97hGQMiSMGtQ\x00" "PpjkIackd3JNtsQLcQ0ihrBvavrSqCe\x00" "HT8FL96HSvhVWif6HDYfWNO5Fh6m\x00" "U0MEdU4CpCV1fI6ZhnQ5uYdbo\x00" "-YDTBc0Wm2Gt7TivoqkMi\x00" "d3Y06RRovEF4EkC"},
  { 372, 4, "kX9u02jSJQ\x00" "4GRHbArCCoBbL9EIo9O6oG9jutmq\x00" "-4hoKE62kAVJvYXebsQWuSmm\x00" "-gLtRvR7Y7hvqcEi"},
  { 376, 8, "HHrO2FHeI6Wu\x00" "-p5ZUTH2ZevYnZgMcVglI596dU15v7TAs40\x00" "NaHqd7mZeisbUHObtChaIPSeTp4iP3Gpg\x00" "2EUV0Jbq4YsNh0kRScTSi0TXnmVv2R2\x00" "7j21S93Bplh0M6aHdI4XvV12qZRbQ\x00" "-4l7eeribKKEISXtUIkCDdJqGcq\x00" "7h7Qh0hZAMHiHMeuMjpnMbE\x00" "-10bV7LCP11EB8sUi"},
  { 379, 3, "17O6B8oq\x00" "-2BkYo0sbfO8atVvdiMC\x00" "1GhBTI0M447AoWu0"},
  { 388, 4, "-1N7hMMi3Rs\x00" "1f8EorHcJdKU1e17aLbrQCdffQ\x00" "-WqKKIC5cTNJVp7WmAQ5tsYC\x00" "-2bnQr23K5tBTskDs"},
  { 395, 8, "2rNX5DXnOB38i\x00" "-1EQE6Qqm3bkLqGMVIkjNiJ9e2JNRH9gdAeoC\x00" "YY0635HKH8h07fA7CHt0jlKdC8UBl9vjI\x00" "LNRqTnmqUCdr210GNNbtcBvvd9eZovs\x00" "HbhrOLGEcEB6PMNbEpOEEFgEhVtqq\x00" "CmC5tGSCiXeL8qarrYeCFJPlLg\x00" "-7eeoDkoGFqTVcMjLjleCW\x00" "4aBKIoWlrkkgvPj2"},
  { 403, 2, "-7Fqu32\x00" "8d70FAdo0cksb2J2"},
  { 408, 4, "5TY4Wg6aAo4\x00" "-E5HYk57B9nlUuPi3Zk4QCY8bkdIWW\x00" "61GAJTP8NLZjTi9YiWl0cmLk\x00" "-Cl713PHOvtuBQWfs"},
  { 411, 6, "9ehOK3r7aBg\x00" "-3MpmHDc5mXRsjXY3ZoRQK6OeSbc1rQ\x00" "27BV3USHXXHOrEccKXioMLNeoqKu\x00" "gc3EnFv9O1EK371n6AXPFU3I\x00" "eqHCpdZ4IRvmuPsgbS7rI\x00" "G9upnrdoiU5ZcPC4"},
  { 420, 8, "-2EWpcYYqqe4llU\x00" "IlPPuiFNF9VsIRTLN3JOjS1Q3GYXqqdoVZNLdE\x00" "1sfuG2lF0A4X7O6c4UKolY7trTXKTTk3fFJRLg\x00" "-SC4U3KOf3MdGphLChZN6v9OWeNEBd4DtDma\x00" "-2P5h4AZQj3oR8jAZB0VLfjqsC5d3S35dU\x00" "1J1ZR96Y3Nvqt5r50GLMcb9bLF7oVA\x00" "-O9gi788KVsFrUBiJ3vD7hbPQ\x00" "-WKgXNpsBIBDJAM5c"},
  { 424, 6, "2VNusf7T7440\x00" "-1i0rdiNiWO0jPFDVs1SIZDvI4DlmIFaEC\x00" "2tjhdlbDa8A84qnisTprgRhIUkZUJaC\x00" "2F9EsGNLYOgtNNgU2Gj2XCYXacVo\x00" "cDljpEma3mUC965pBGLNiKnc\x00" "-hsCe1VqGLsXoJNou"},
  { 427, 2, "1NnBaoq\x00" "tCCq1M2062uATRg8"},
  { 435, 4, "-I7OoKoUi0\x00" "KNjHFlSbDPNDCejspt6o1bg00\x00" "87G44RUT4VQf9vGFTqIM4G\x00" "1h7ov2ZFI1rEG1Rg8"},
  { 436, 6, "-5trOSonvfeeu\x00" "FXa46AVSOmCkGmWAJuM8uQmvb7T2EvMrU\x00" "-GrGereJBLNa3XfYKka2g6UeVCiM1mpg\x00" "47dYis9EsnhEQsGZihbG2KOfYuaO\x00" "-2XTgivr0fTn5T6DGlhCJfPGqW\x00" "-1p2Q4qsLBkOhYqgbM"},
  { 443, 5, "2elXQfT9hia\x00" "-R9bgRDurvYs1oKQnAg8RRZuQT5er6\x00" "1inIlqrsipINVQH00gSrP98Ae5Ka\x00" "-I4KQDsv2SRTutrYQS9CWHs\x00" "3ADmDbAnaVUqhNWu8"},
  { 451, 6, "8R4slKsQ7o\x00" "FtOXSD5JR7gmOtsFnFqMYge3TJia\x00" "4kGX7RqB1eZLq1Qkqj0ttvU9S4\x00" "11fMXYY0HHN3hLCslJo0TadI\x00" "BjrkfVV1AtvAOJmebJEW\x00" "5hlaqlImkJmPAvsTk"},
  { 452, 8, "-kAnWc5ueVZGXUdg\x00" "13b7F6JOA009rKemMqro4nRv6MNEfHbFaMltrApCMkbM\x00" "-1ffdT7KZNfTir6PlFPJXvqI0juf73Loj1iShW0m1km\x00" "4tMMdNbhePdaW9Lb6mld6n8oHIVA9784s46Su64\x00" "-DPopU2JppDmG55LXJCh2fctDah70KaGm7cFg\x00" "-1pdqj7opVFEpQfkCh9BtiBXG8AZK0Ai\x00" "-FOjmToBlRa6FTvQL9SeqEk7k8\x00" "-6BPAGMbZWSa3Oba7o"},
  { 456, 8, "2C6qpRuHg1Ogq00\x00" "oEAoIsojuITrIaXccnFmjqTrFYkIrqN9iEC000000\x00" "EErmCK7hX7gV06ncepYtQG8LNdKkbSsWvj0L0uiO\x00" "-f33dS7jf1ki7N84P6sF5MRtMHPHV7uZ1L0K6a\x00" "87qNCioBYnOGlp3enPmE89m68aMVEbNRuS\x00" "-3EKs9banircgmlHNhuq3D48sDhZQYhQ\x00" "O0XpSsXVaaZWoqlgcUlCL7gJY\x00" "-8IfANXBKORY76YYa0"},
  { 463, 7, "2ZM2ApcCJ5QJsJnv\x00" "QdFSMo5Z1bZRlATaE8PIWZZUhZu7aoeJJ2goPNRaMZKM\x00" "CO66g7slLrln7VfC0kNQQa2c4XrTmjhiN9jiZTcPsR\x00" "-95VGGeTtHo6OVTQTDaUBXtOL80UNnlASbeinYI\x00" "3JD1qn8B6KUdjmShRl0vmWfV7TVpbDeECM\x00" "-btFMjujI5L0sfK41YBWaF\x00" "DqDM0bmO8LeLnW7er"},
  { 467, 7, "24EmRu7iAbc00\x00" "-Mdhs4mXikI7SvurMFbgvkTDcCu4SO000000\x00" "GcavXnEEYV2m9a3DgraUsD3pa8vTOmRlk\x00" "-2OGlKuUrceDFtTO7QrPO2946qPPAsXk\x00" "XgsO3IoVhFj3RnavWsNKEgbvkNjo\x00" "-3GkKAIXHhB2UlgFpQa7TdiK\x00" "IYtFY7KrUfQZqJJ12"},
  { 472, 6, "JdEcm7kWqJY\x00" "1dJhhpenjqKrvSUsOnE6An12ljYAVrI0\x00" "172Thn66YVoNNalYM3tWM1Tq0MFuOnM\x00" "-G9F3iI9Wu6VohYGL0o4dhUOMWv76\x00" "2MBGBpWaQcvjb5mM08njGvW2Aa\x00" "-QgIbXuea7kp9qOvZk"},
  { 483, 4, "-q39DXCgjM\x00" "1IQhPt76jMMG22ZFD7VVrS4XgO\x00" "FHoqG7RUTG9QhWsVFJb0mPc\x00" "10sSoIWMKuknj7rH2G"},
  { 487, 7, "DbJtP2dOicPmE2Lr\x00" "EJZXanKQF97d8gA1Ibk7UdZV0k83kcU2uTT1STWZXH3gX\x00" "GAOLucmF9KU52Zo0FAe3KAhVFojHhGEb1blUNnX21ES\x00" "7duIrHP0skGrFNeiDG18IKeCK33kqooJ3GmoF20\x00" "1luSmKb1p75UqHQshDLWNsJAA2UJlA0VdmQ\x00" "5drHWsPTEVD2eLRm9ho3iS\x00" "1KNAEXB0hsA8JgWXQu"},
  { 499, 3, "1Hf5S40h2\x00" "-WcvR2fG9N1SbIigqXgB3c\x00" "39DeQQCh7siagGlDFk"},
};
#define NUM_HILBERT1_POLYS (sizeof(_hilbert_data1)/sizeof(_hilbert_data1[0]))

static const struct _hilbert_poly _hilbert_data2[] =
{
  /* Adds 39 more values to fill out through 500, at the cost of 21k of ROM */
#if 0
  { 119, 10, "-5FLqclQ9lhT\x00" "4JK6gMJSZruI6cJ7Pdg4SOpiYc3aET0\x00" "-CH6esl4n86hvsmmVIaKAb6I1HL6tp\x00" "4A74F8Nd4a9Z9MpZbMTeGj4j7gp\x00" "5scJDNoirHssuoLOKD8i7FLpd\x00" "7dPoGd6ojlk5D5SdM2PLKc4\x00" "-GU2pERaCPtq2GlbT600D\x00" "Zd0WtCC9MlMSLRFah\x00" "-2l2gGGqS6Dpb\x00" "5uNpqY8Ub"},
  { 143, 10, "-26RbV6Y9cbhT\x00" "bGujjtraDL51hR8gv7seS09VOE9uAohT\x00" "-1fnLCn6Q1edoYuWu0pTVJ2dvGeghUeaI\x00" "4M72vgQdZ8F7QACSR5EnhCVe5QBXVH\x00" "-1aiooQv6bh0hHNJiRQhSHa7lUuFe\x00" "GVm9YL8cVMI4hFigBUSLj6K8P\x00" "-Bi9p0Tccl7TfsN5VXXv1AZ\x00" "7hatPkT2Mv8p5sHHZHD\x00" "-4EDft4GmDsfc3\x00" "2jSMZpBX3V"},
  { 159, 10, "EcbTogjLL7mT\x00" "-829HfjDIQAAm3s6nr1pvtf3YeQ5qcUFtK1\x00" "SHBjhJamQleYbqm9rsCf3LkkWeY3NVv8\x00" "-c49PKi7NqoL0uFCXqUcarPlDb5isjs\x00" "MHISbMZFgiVUrIdv9J0DesYYrULt\x00" "-469F1JrGN3I3jBFv8WM7LbAOcJ\x00" "Kej0Pld3DBpH4dMJWdhXNVh\x00" "802TXg2IWi6A9Y82XdnB\x00" "tXTp6Sj7IUosr\x00" "LVO9KqDkSf"},
  { 167, 11, "LUBo3QOB7dUml\x00" "Df7GMLrC1CT0QPFYCnM1rF3q34k6fIMCndkmR\x00" "2Kil9T72C6pVl3Z8uCu67EKJVNfhB71IqNfB\x00" "-CHeDAveVE3iljmpmj9LcqTWS7lSVBZoeSv\x00" "BnETDIEe5bBhRZ9gPnM24FUiOuZMtJFI\x00" "2JmJiYtLSLRmEKb78aCoUWNjAn6DM2\x00" "E8cEUQZAnKbBm9E1C3piv5X9qoX\x00" "-3cRPAWZWaSZuRkrImFCTotQv\x00" "vHA5eXI9tNnE9ACpDQSG\x00" "-3eCMf7pKD79Bb8\x00" "vbUI3EE4fA"},
  { 191, 13, "1LTtEREasoMYFBhT\x00" "6H9pAmd0TS1QtbUacQMMsdBrZ09XcTcee0lRUtUl7ET0\x00" "72DFp8aGqbLPruEq6C9aMcpFfbFIF8Mh2Ndnu4InrP\x00" "-1P8FWc2XreFNKvuRvX2CPjFcKQNIl0DE5s0ZBNqBp\x00" "1SNUSHUbuLRL2LLQ4kHoAUAKrYIGs7TqBB77gdC\x00" "-FoPpCUhdnVI80T12PNbafKsUATU9WA7vqgrK\x00" "YgOHPcjoCSXY7QWZso17ncfIGb0jm9U1ZR\x00" "-3AfP5QcJPJMNDjt9mlX0WjJQ6qcUIfdG\x00" "8MpMDTLa0I24r6tcnF58Jb63kLJ8u\x00" "-nD4oSZpTpfkcUr4mmTPa37bdo\x00" "4jdhiSRIgmIoKl0P4ih2dL\x00" "-2tQfJmmLhJFY2Hn\x00" "GcTJ9ZZFBn8"},
  { 199, 9, "7HU1GQu8fS4F\x00" "-1dbjpUu0GkfObiGJbRZ5F1Rn4IgBSJWujP\x00" "HUuktVs98G3Mkc5TgNKcKAWLLVRTNffG\x00" "-B8JWNHRah9i949oXMc5v2NpLSeFORi\x00" "2GGvCZL3TojkI4KUekPmIbfeXPOZ\x00" "BNrpYsQQLHq79ZJZn6siV5b2h\x00" "Stf87HqPnaPXsU5K00KHrA\x00" "FmaY88WPVMU7tO\x00" "ev52J1QBcQm"},
  { 215, 14, "4csTYQlgNUUXvsT0T\x00" "-83FL7iiOuWsaZbp5PCSfAWDveHoqiXO2D8OkGFFMbNjWpvsZ\x00" "2ipo2YK0AUC8780CKSjeN6PdbcobgGivL9OZhNecWe8eIp\x00" "OgED7Acb2mjqfqAefOEa5g4pjvaQb4Jq8iiKqtE81INC\x00" "tGD7aSTgL9u0tmj0ccbKI6LkBEbji78GUm1pLm5r9F\x00" "-2XCPnnoNeM1rHLdGI8Qt6HM8cepJtUnuWN4W8535Z\x00" "n9T0hBZskf0sWAoFM8joCTGgCBsPYDviQc71cR\x00" "ikoAH4kWtlLIHJDOXfoQ3v0ir88gILsvOtAf\x00" "22HSsvfgCdKhqQMPB5q2H8eBXBJMc9dEFc\x00" "2eqoiQZmPqOpjU0Y9OiB4pJpmSSiTsL\x00" "-6WigZgXHsLZSnbqMA3gcuSHEDMS\x00" "GPae2EjPCRoVFPcsmXO4B8l\x00" "-1ZsFapMhDljVodkR\x00" "43AX0DADeDrG"},
  { 231, 12, "-iYNKXLmYCe12lZr\x00" "mUSZY95JM5kefgOL6bj8ETthQcNosli8AgeQIsnCMAY\x00" "hv6HLRl4jutW2s271lDF48C1UUHPNEWFNBR6I7mq7B\x00" "-IvWsAhrpgXcW9E48GWh9je4vuGnahoBsfhRhe8Zd\x00" "11aOSgYHvX02jHDcTVqq6BFmeN9sZH0mPOq5cRC\x00" "84aLPJ71u96oQAlmRraN48r07Tc29AYivBS9\x00" "3mbHd6dDna6PeEn54A2oqtZePVpkYGf4Io\x00" "6jo1fW6a0XAosloIh1lRr0CY6qMTshM\x00" "65EdAfA7rud03OKbr2EretkQrV7I\x00" "8CRBm3Cp3QkbaCmAJk6JpKYF\x00" "GAR94K7TVugWrn17\x00" "Lm1u1II0MvH8"},
  { 239, 15, "9sA39TadkD0lNT7psl\x00" "-LktKkBdDIR2fbVqCocWIbSTaAC7crURn4Bj39uGQei7B9m61bFa\x00" "E4blHcGf6oQh08btaDNcMQpdHgho6bQrvlb7s80no8ZnvViW0KZ\x00" "1V6V9hpdQSKFIjtmLOW7gtvZlBLNJQFLGe71vgD1WMBaN3bos\x00" "-D4Xl4VWWCv12tm0cYe1Naa2CEEDeOoi71BXf4378oUb4AZS\x00" "-9ivvg4qNkFgK2m96VEk7v6rK29PQVlilYmb4fYkmQDl8T\x00" "KE12ErkYMOJP18eNu7o0hL3io41rquRPlL8Wo6F0Xb7\x00" "-6B4Co7QCNoNO03Hvv811mAlvpibi5UGoGBEoTY0Fq\x00" "iEpYmu8vNXXoY8YUNv9Yi8svoGmGvFSrUr2lWc\x00" "-111COcvKb1ILq10hQhQatmA04cW317hWNqFj\x00" "dAlp63rKuLWKdhadLarMSNDefZDArkCL\x00" "-fGXOBtmhHrJWlbbqHBpkY6BTh6sK\x00" "gJ8QnPDm77aCsWagPIDeHgus\x00" "-mJJvtuhCut3OmoVI\x00" "nVoFg7CjQOvM"},
  { 255, 12, "-B2G36rN3eEBuAnTP\x00" "LEnJOqhsHAn6OivOHGWgnI5mGNb8PcXDfI7tWcb0ioFfg\x00" "-4WKfa7E0LYCb2m4tatBe61Tf8gNuUXFSQi788BvQjMAq\x00" "BfpgO515mYUOXQTRVBFcM0oiqWiPVMJPG7Qk7KK1j9\x00" "-3r9UYBgmd5M9pMvs9orYEbpNQ9tcB6aJieHpgnuG\x00" "5t8flTu9FKao5IE4oX36Geor8utXeA2J8HOstX\x00" "rqruI2BTa9hNZ8RhnY9s7oKbkbffDGTaR9g\x00" "-tPEWUQj1BSnCQYGCBkTsimWpjDnlb0jQ\x00" "SIRC03cUeb5Rv7EaHoljCL3itcNRb\x00" "Hp8m8BvDe5RjqXSc2795JN9fP\x00" "6eLpPvZVU4lO6Xg76\x00" "4DFKmm640b6gp"},
  { 263, 13, "18K2ItVBfA7sM0DRhT\x00" "-HVFCLtEWm1GRQYKRhuAcUZhRZtmekcpo2sP2L2F9EYbi0qWaET\x00" "1N1qAeoEmTms1STQWH7H781m3kKd9Z9KIgmc3LX0Z5r7OYJDb\x00" "QN0hmNk331GZEGqb20qdCvcvhI9g2F63RrZWKdL53NcW8o\x00" "-WMr70g4qZDWuVJAdicXZ2F1S1tLBLpZBt1JctE3PQ1lk\x00" "-LpFqcPBTHd81cJQ1JGriQEKhgRA7FnQrssgBTbu0il\x00" "JH30d9TqsTG7iRgtdhUt9H9N6M3N2q54cBc0Y83Y\x00" "FZg5f7c9QvMjXb1qhsGp91G9bNCrjNSEjQYZo\x00" "6Jg4nrP9qkm4nGSmYVnFToGXuI1v8dPFB8\x00" "-31EmJqfWZMBOvPtTQk3ccmdisl2Sbe\x00" "1RE6bKeZ15LlO2NqKvZatGsrro\x00" "-Iv1V7Ba6ol9SqUvSl\x00" "9DSqrGHIYXsRh"},
  { 271, 11, "5BSP1K1sn3ZJ7hT\x00" "coYnBTp9ql2KqlhA3WBdSiib9VvGBfSGWRNsbDnET0\x00" "BeHvaaCp4oquU6stTcCoMktE4TGblTFbXSMeVi7cF\x00" "-QQ1GJlgAcBHpLL4UKQYQ7ADdoKhLXo2iDrIv0gb\x00" "h4PgkYUEujoLaBTbCZHKmK2NUrmH50KlYdE22\x00" "1ti74sJK49Gq2suiKovG069ETKZcXvZD7DV\x00" "T89ugk0ilDHLIuYuffc5NclmDP8K7Xb1\x00" "-JsvhXiLrr7RA7gFZ6C7QdCG4g3QnC\x00" "6mpQRJ069KNo4WEFA4RmFjYeB1\x00" "-10DD6As6j32e1UWUr\x00" "JrReaN6BKrt3Z"},
  { 287, 14, "9pOCMatr1jOUQiOt8ET\x00" "-LXhdsUANXHiBWkm9Kq2s09nkkou7c6SDSqkUU3phrXPKIbLoCvR000\x00" "-1iqoZltAklReKArnPPB86kPBk3pdoUceipPFtIu6hFduVYNoE4d4G\x00" "2a7pAOIAagSErIr2FhYVttdAVQ7UWpTdQO2FRiSs1TboiNlB93O\x00" "aivpGbdHFeeOdvcQEf7ilBJBYcX1eP4a2DDKpia3iVQtEukQ\x00" "-UTIgp80k1MapHlbI0KPoStY1nFB0317fE5p9nTSEqMJe2d\x00" "-BWrihh5vBGDtLP94cAm3GTT2tgePbrSXdlh2i9smr1ev\x00" "7sMW7sPh127UXXmI2uP7gB8tCl7X9aV1bWrZpv4NEo\x00" "-3c4ncJsvVISjNvYuAeaJQBBWIXkOp0ame7YRuYr\x00" "miSY5GhJJc0VtuFOkWtSPa7vg76qEIhkVc1\x00" "-AphSqO7vPcuOkTZlcc27BUs8imHbFhr\x00" "2MijpGarM1Yg2cGoV78SYVmEtNS\x00" "-72P5MCVNBaVullUA0H\x00" "1VfRPbhZ2dlPuf"},
  { 296, 10, "7IkKs8I5ALHLKC\x00" "-7FYOkD7Y4B7nBHLaKiZE6Osi0kv2tGvGQfaKClg\x00" "Acv8HgIsKeetSJo6k77fTTQqovaFmEnUGchuHsm\x00" "3RnlKqs5KUeJ6ZBMI8IrUbcOJNXefIckejujc\x00" "1nbD5A0e79mvLe25cpKv5WtmcvfM9C9WMT6\x00" "-3tOiuAOt4NDW08DjPdDYBcYBh90PnAgS\x00" "UP8UU9McFUHCSKhf2MSWcDSMSW6Mi\x00" "598smd5fvM6uorDv2E44boDou\x00" "ogAstF4WqbVg1QkNqJQW\x00" "-3VJtpdgqZ4WMA4"},
  { 303, 10, "2HffhKpI6fJmCWSP\x00" "2OU81CCZOlmWIDQMib6EVPbB7H3ZTOCpaj8rXWhgsrjTH\x00" "Dkud9Dt1CI3JQE5pReUaDS6l0oWmmu0aUMC2HYXOZI0\x00" "-XuAEbR93R45PjuBuC1cZOs4Q7WlWasOFQ7NQaDFI\x00" "1qnqlgpXMZMeOLAR1pARonqgvMLHZUvMBkmWUE8\x00" "-Z2PEONNMchY0r3WcIAEuSb2Nu95sDKRlhlL\x00" "5UmPgj2eLVT31HCb3cANDqZn36CnOm2b\x00" "icd0Z0jfLS8vGEgMiZWiR7MEZpY\x00" "m0F5iVs93sYYkEOg0E\x00" "6daduG4mhEAq12"},
  { 311, 19, "8mCX4ZSFUOhosQMZhNaQnhT\x00" "-5iQVMCiaXGJQsVk735vQMbIEgMe3lZJdKJ1kngL4PVZ8ZNMDnG6Y1N7DH0btBuBe7ET\x00" "23ZbnRs3CrfLsos5S0oF2LvCC1F3vno7djAma7egalFGRAcYsFo3skYj2ZlennYAZp\x00" "-56f7CpTt3VeVN29j7k5qd3oPoIWHVdQ6gvtihCQLvlEXGW4Hb5jpGiOMSo1PXWFr\x00" "4O0Qp6pQi4Y6J4BAIRWg4C4VcWZD4AtccUckRXKMu2RMfTVAf5hn5VMjMp79lK\x00" "-2VVE9lhp9vZmFBFpuamS5CQCR1JJA1MFYG32bJciJITafdWhkH6FCJXfOHu2\x00" "L6eRDskIc9gToc2O7U6JCZnFpGHZQXb7OBtnn2dPT1DTb1H3ptIM6OeFQ2\x00" "22IhTO0XK2AsMZ2RKUbSL32TqsLdDGcHL6blWZ6OvTqS2vYBYuqS1uaf\x00" "-25o4Rq9cforcS8A1pZv0svOJTopuskL4I8MOBgJnXbplsC6PIGUoMA\x00" "-Ji4ZcXj6rBfGCD19Mde3O45gtgKj5k90OvrrvX3nQ3s8YVvec9P\x00" "Do9Fl8OPYq9aMqK9F2vcLrmpEm9sR7ovMCvU9nOrRQoL0m0Mv\x00" "Rpd0S9dilmjbUQK1HnO8RuB9ZrorqnBLrJEVORXCFReCLt\x00" "2PPZ3BCcduvgL8lVZ0ZZ2YOrHXavTMGlpceBP3NZ0fok\x00" "aPWBiqWfGLcSAM9LqkX9hWmHOQaevjO7NnfoqoNd\x00" "4fFNFrevI6MD3WV5eCJ4oteeRg9rpndnmkALi\x00" "-TTJi7NTC5UppObjG7KjR4Ac2W0gZtgqM\x00" "3DRg1LHZWmWl5bui8AJTFLj4aIs5\x00" "-28vGqPjpIWpAkCjpFQ4\x00" "De8iNIZY37Zn4K"},
  { 319, 10, "-NMX6e07BhuAsT3T\x00" "RmI6FIsbjeXYFErIUTEToXVZ9ml4WvqTCuLq94ITULET\x00" "-1akH0fuWURS8NAHj4YlQaNf2KFRNvkhSkHNjEP3dvRp\x00" "1T1AIP04k8EglWE4mAbOAuOka0ltOpAc9dGcRSSIp\x00" "-6I4poh5bBtO8nAqf2JVKKjg786oTIiplsI1Fcm\x00" "AYmSd1oCHUSXbbXAaX5B0sbJQK9fO6qrrd4\x00" "-2C5ecOPCMFPDGUjWJKKWQBlThPI1apeg\x00" "DIYfdfAJ6BJNCQavN6MdFmvZCsHE\x00" "-4ZTUKELVuj5p9THfi8\x00" "Rk8YNS0d9AsF28"},
  { 327, 12, "2ebPM4cTtsWYBi2cET\x00" "-hh3iI8rOHaUEXfgd5JJW6dFGSD4TaXIpZOl28rt7nR0lMeKLhT\x00" "4VhC2jRCYfKX6nkF4pfeFCTPfgGQWBaRRTgbBSmjG0D0FnQTV\x00" "-3d3vni9I7hqVIThBMdO0mqmqEDXtS6LCisBFCkge6FtVhLk\x00" "3n98fiIAcfGGgIQSuqUItbmdE91fUSMADHXlAsQ33BRW0\x00" "-ItjP89vOnGdHdtpdQ8JrBGI6nMDDn2FSKBEpQgVmvm\x00" "DPleOlvBLY32lKj2HgcPcZtUkdGULJ8fgnp1f7N6\x00" "2VLktVVuvaBUpgMAOtJf7scFG4PnMQ0t1SHXe\x00" "DsH8VhvmW2dKalXJ5lrC6d8F0r7YLUYlo\x00" "rsAuKrVuTp1WnX0lBWJNs4KSbEFo\x00" "EQBaThZkSUX7bmq0lqT\x00" "trtESKJKU9nnCA"},
  { 335, 18, "9VdNeHBGGai9oLt4SU6cilB\x00" "-14WnHcJnI2GjcNqFPHfgIm1Y4VTHL6VVUJLhkAqNeWqpPSYo975PfUUI7QUUEJsrn2X\x00" "-3gqgGRMaUorippXL7bRCRoFfVCjcolNZ7XDArNUdZYPZiJckVofNTYYYuF3hXItIG\x00" "DjVK5Ue3KjXDkYZO8OsET5cLMdpqr1eROMYoSm6but4XKSmUPodLorSa5LUHCNU\x00" "2UYl6L5Pu5O5XC4AcvCBdllpdoUGdLTeu6YfQ3GitQeIL0InGLGTYsds3U3qOa\x00" "-3ekTTonQcUeC4WIcYXga42R3V0mlbCU0uXSfaRU4d7VNoSYRvXnOk9rDE8i3\x00" "1EQ6tg45cf4cWrIGYroBjRLM8UunNiBKrbP26kHK3NYsWPBJbpp0nqKain\x00" "GA96DOSeaAk9WeEP6gbsaDZhNk08Fc1bb0foVf1em15FfrvWnS3U5Jr\x00" "-3SclEllL87IjCvJbRuf1ZLqamh0g9PsegQ56evRfUhhQffpWohkSr\x00" "-102ckEUkaYJPSXXk9V8NN3dt2P39pFerjKPB9TbLF78AfWJsjcc\x00" "Do8j3uv1LScM3T7a1VKdPieeDjvWkal1pDh7G8UlAFOtqSoH\x00" "SXP6QNFVPWTpClqnQGFIOK6Fhv0RT4QX1cJigTqAoWEbv\x00" "-4eurI6SXgb4TgveIBJrEpAVZFtNLIrWKJrNSsvnNY0\x00" "MOCJ9WGhXMWS4e9O4ZqYmuVsddsnviQRbmSB9I\x00" "-1BiLXX5ruAFE03r6nK01FEj5CS03U6EQYv\x00" "3eeDRc9PDEZsZPnhr2V9p0n7SCUTf\x00" "-aF9V6bH0jvVHte18IIY\x00" "1rYBIJ1Lr0Ylu6W"},
  { 344, 10, "t7qrJId4B1cfFTg\x00" "-5dEZnhebQD8VEA6e2jeVZHbMXPjCQjsZTUbgNqYq8GFQ\x00" "A7sApOsW6kXX7ZsONfUl6nnobGje7935vSdiOV9sK8\x00" "2R6Sj0fTJGQkc73bICiugQg53GpLIgXrGHSYUgK0\x00" "YuePK5iS7ODg8Bbsp43WRtToZh1KA8g9ISosS\x00" "1RVYZZTMCb14cH5402RJTOJgALI5VWUmHmS\x00" "1QHmbYcOnrl7ELgIZracRA1BbkNpoDLA\x00" "3T9LlCTotjS9eEtpl4VLfm6Vj7U\x00" "8PEUEA3Z3uefPYlhTL4InE\x00" "-48KNv6H1ZVbL0nk"},
  { 356, 12, "-GD3vkKsZe2YG9aq00\x00" "3siPLd5S3OvXWIj6CgNG2RdKr74P7fduf26W5N01WY8shA000\x00" "4A6qFIAdMSgRmG42172O9084n3oL6RFleYYjhenbPDT61s\x00" "-524RvAWbLob7fJe4XhZgq9VEbbLlYnneeWUpSd5qclkSu\x00" "-348n9hVBbgvr9kZ6ZuGLnETSBYKAIjU8Opdo6b6CPXY\x00" "aL43j8OVBWYXveruH8nJUf73j5SvCcKmo9DMn2Ro\x00" "-7JsHRZlfSKXuVC3O8NNKG6vjOYDFNbaKgDV7vo\x00" "6tYcceas3TWN0E2RsAXEjYL7ac6LbYSCtc8\x00" "-CpM98YJCaMqulXKITVPrdGI8DDWTVpk0\x00" "-MBLRDH332GdaLSBqFXtVaTORZ7c\x00" "-cE0sM2Mk9n7q5QH08n9CEi\x00" "-BKIZT7rF87cpJDg"},
  { 359, 19, "GUvrAqUEM53j3ReKrCroSLhT\x00" "-7HSvDm8YODDrmuRmCHSF27sHDlmEZvRRHSvNM6jKCD9TbTsF835P5jgPhnu5TJ0Hb93aET\x00" "3ibrbOgS2SPETHaUvEVpYLc933elZHtU3cZFYP7uJhOWJ0lDaDCfFUaPNSKXtnrr6aLLhT\x00" "-FBWrtrYAgMmjFhS1HCq6CeqmvtEBIQQdBmkM3eU1d2GSU02lU0hJE8td4S8EUsGcYjKc\x00" "a4FHbrmD2atNPMEWq8QEELhaXcISWEgRr7UkdlQlRVDkgqN6TT2cSaEnPf45HrLhgI\x00" "-1KPXpp3siODqoUIjKALCaXKn6iGheaMq8tJ90W90FZi7BkAM8tkgTZKhBi4P3YLWh\x00" "1H7HVdTXOTMGL0RY6RQ1j7cYij4DhlihvoSnV4bNKNfAPhaVGIcVGglgJCciQIW\x00" "Ne3WgFOpsV9g9sXOpYqWruQmTn3PUkk6TuV463E37fEphg2357M39rtBKQ7\x00" "-20gHY1TIhJg2BF3IC48KvEiQ4G9uTLMMBjJSlGmgdDrTCPel4jqeWSloH6\x00" "-2BZNr0P5gEFIrhg1KojvibmJDsIv4LW9qEmqslCB0NLvt3V8DQh9ouU\x00" "3MSRBkT1SFMRdV8d4rQ0sj9Hha7MNqmo83qs5uadQDiP9QODfYu2q\x00" "-6fjo3lC7nLOHS6hKfgiK5KEjHnlotY36nGh1a0l1fpaP8SIWsZ\x00" "4qEJeWBu2PnX0i64TQPlLX8UdYRIJIoYtVPh3nYpZ0dNvQ8\x00" "SaQDUMAiKa3LvZknJg23QZJ8JrnERuAb9biQ3u1a6YO\x00" "1PaF9j6P9aOKROnpJZVXbJJfZu5cENq26SGXNAhB\x00" "-2GVokGIMql7oJYfXuTagWCdh4ct6fcfXalv\x00" "3c4arLTjpWRTUY493c2QEs45qi0Zqj\x00" "-95qc4iVG1eIS6ipvA1Fd\x00" "EWbl9RqkMCBaZpG"},
  { 367, 9, "g2ptfSc1sGXPhhT\x00" "202judZCopXDdSW2rIjbWU3uFXp98utegE9CScENGhT0\x00" "49AHCN13l5CGbFb8BIL86fQ3b9bIEGLk6e8FUmAEs9\x00" "2i6IsY1F9uFG3SWBpkf9bosXCqUQuoWArSV1U3dP\x00" "1EM279pd4fPl0JlEo9SeJZEfjsjOWshDhXiTBl\x00" "-5mr1Jksr6neQ15XQNl3DUsM7eLhFYPYcNl\x00" "DdDBPEjrZM3Uc8WGMqpSXNrAlZRKWq\x00" "-C292MSg8U9HIhhIJV8r\x00" "S9Ue0ZuOHHOM9j4"},
  { 383, 17, "1cdZCJ8U83AqfSaHEPCQ7d2t\x00" "3n3ofKv2Z03mPcSCZ8c3b7mjbXqZV7jnc2ivLcOHRF1SWLpBN1Xr39NFvc3CpZc8MTqU\x00" "38R5EMomR1UlOMqSFs9eU6BOskeJQhI4c3UqHZb2DBk6HMUeteIoYB2MkpOt7TjZNhn\x00" "-8bfGoNVEf9vp7MCvEdefToYQuf7pDrcXgqr9rBdV08oGF1vJZ8XW19dvkmgNN35A4\x00" "7tfFDU2Nhfc9BD5DOE6V2C1PbrT5ljkPYhPqVNqttl2XONvFQJ3JPPaNObKKKWc\x00" "-3IfM2clfTWpKN6ZXHgpcBmb2uT9lhSLpWeYVAuaJtkfAYF1XbBAHeDumAMZ8o\x00" "13XsEYhmCDodZgKQldFV995NE1OZIvM4pU2nKFYpUhaJSfTPJJcQN05jMWS\x00" "33L8FIdcoCbvJuq6RBn2ecoUZ8TFJjuVoeN1ZXc7p31lDnZdDTPFlLqF\x00" "-nTYYRRWaqFuFNTurR0NhW1jajQhtXGMlYNfIjvEee33u2NKnrACDB\x00" "-32cF0R5rnXvLqa45IQd5AZLjXOd5SRZY0eMv9GpklQsI4INREWu\x00" "lI7hCq7ZcldlPJZACeVCsrDDMNW9uTMLqoHfpZ0ZcGufBKl6\x00" "-2iFulrSipnFpNG5KRRhqVZolUqHa5QSsFCKCoiZnJfI15\x00" "4dNvE1gAr0BXnpaHvjYd0mvY4MMbv90Hqk50NUTTU\x00" "-3o1V7DH1pZRYtTefRLTjQgR0jdfuVY0fSm6g\x00" "39LV6Nd0eUVWakMKuSXfXsBHtjbkLev\x00" "-29vaR5aBDMr4d4FfiiOA1\x00" "1j7R46YJIP5Wnm2X"},
  { 391, 14, "-13LgDcHiZ7YBHVCFJFhT\x00" "RdI62Vl3fFEqDKdPrfPeJGXRdvttkNsiDfHrJMVmCWgEEPefWnli1hT0\x00" "-4pTI0lPVkS3c2vDshUSK4EqSBcqVcJ0maipOMrTLYgnidY9vp6Sonjj\x00" "CtBaEX6dmC5qjC5UGDCKOS2YYeADWXkWA7B4bfe0mdRniROQKVijr\x00" "-H0vBA5WFrWcKMbDclXiLWv6k1ObI0VGXRirZ5VnaaGgZmDboU2W\x00" "Wesvb3J76d59O03AKOmR5veiTepsm5dnJ5JV6fI2OJlE3OLLm\x00" "-KKBe3au0Wf5gucacBQtWbDLOkrlgRWZbAnTfO6Yrhh6b5uS\x00" "4EXTqq44Eu6ONp3dvIQNZs18Ekt7X02THEn5OD8BhR6Fo\x00" "27DjcFgVJILPihSBcWfTkE0hF6l9ce9BtGl2Mqh01Y\x00" "2mJH8tXvRAnbqArd3vCs7dKcNt7sVR98hUWbf9U\x00" "80ZF6WDgFllP1hqGJgWd55hmeo3JPYRUH85\x00" "BJru4Z6g5b4kh7ilfpAY3DpcQacEbXS\x00" "2LrNZZG7s8flBRh0Fqpf\x00" "3LK38Mj4VrGqpE6f"},
  { 399, 16, "-IU2r1jl7TXYMfDv3RnJHGJ\x00" "-17jXZReNJKnDAPBJSIRNC18o9unHv6bUkHD8OJuaBA8htuUWeeO8vanFOq3RjunPc\x00" "-YRPMYlH96UGKrga6mYPOOPPNjRShDGaMMcBC7XT9fJFqU5a6WqrYA655Kvrh69O\x00" "-jCaafqAsLSk6j3btJrdFXuvGR3BKiJ9B4nV0kMJdVrQmaTGTv36YgTErNtH75\x00" "1b5UTnMRE7VTMK9UJ3o0a9QSa03des7XeBJU2XFSFS59e9ERYOoUWunp4SOv\x00" "4IaQYsftK1ufi6dpl4LnY7muipYhBG59pOAi8X1UTdZHKRCkREpR9vep6g\x00" "1WPe0N3PuW2MWuNghv9uTjA8lbm0cb1d0Hrtd39AuK1Fdiee7lvgF4St\x00" "5YDMPvu2KXe9nPYJIQh9hDFTVN0ohm3soTG1iUV5rHL2QNDbTOTs\x00" "oEi7W7IhEc7OnifGQlqfT4v2dWiE1praoBKNdJG7mGQ6MRtceq\x00" "-10rM6ANMr9ZZIJ9bKvB0sD5t0dMN8Y4JIOqmD0JLLTLVp47u\x00" "R1ViiN942qhDZaBVKGOd7SEIOpf6gNt5PQiCjKDeLtO8\x00" "-15A6d3lBP1To5cfkWVfoBjuLmiA6IIMVhQUIS1qIl\x00" "1FhvrCkCoSZ9DXOOqiV5ec8Q6KefqcnR8qIfh\x00" "eA7hZZ4bYYZ3Qpg8WCt9RUD8AjoTD9U\x00" "BasMMUHdlOuvogtJutBKF\x00" "6JaCfnuk14ppVVKH"},
  { 404, 14, "-9m37ksYmj8avDMD5g00\x00" "jD4ic0ChH6DIQHDP5fn1cnhjKMDqksHfPo2XNET745jWBdLIIEC000\x00" "9R6qiK2AU9J5q4iTGrMKrrjEXmsIZX6I9C3nbFiZqtgdoKJcqSNXg\x00" "-DSqAmUb1m5VFg6l9KA5HJs641qC8RFX5ACcEv3Ulm0AUMWjrlPY\x00" "-8lYV5nqBJ8uECVU3dcQ36rusjFII7u3ZAvWA10sIEoTUDJsWGW\x00" "-M5ANLgTnP1Qo921vS1YNY3vB1C196EtY44fuVcQTp6rjH1A\x00" "-1SqBJ7TC6PWTlQkDtqQ2Z3ik0hFEd0SJbK8koIO3TSRN8S\x00" "650lrpPksUl4B0KltrI6PL3dHWmmEqoR2Pf8mLW4O3o\x00" "-XaLAqXVAlp2DKhCTclFMhAWrfHKchrWNc9LjQLYO\x00" "-dD0Nvtug2IJlb2NLMuT9QAZZVPXWPf5uonvSS\x00" "-GjOXHKbMLoH6AjhU6KQUmPLGHvfqULN8GO\x00" "-7tXTUqlcgq18ZYMJKX29fKGfKatlU\x00" "-3ik0FbfPj7oApd8OkQnD6kXk\x00" "-9M2TMrZA1gs6bkDk"},
  { 407, 16, "9NdpgjVublBhWnEt0IUTW9V\x00" "2Wbjano5GPJOghXU9paRc3qY6ifiA8Aao9Jv1ndjDQnJK8eEYaJ3duqd7VAbcJoE98A1\x00" "-AbNfPa7gZOWlRkYPiunjWksMgPA7q4IsP7PUdqf4oqMcJGWm02ttWKaGG3k0v5fA2v\x00" "bH63M3AKTNlgCveDlH10nRjLHc4fJCuAbbON0k462hrR9hdSNBS6eUPalCUisEXM\x00" "-eaK5RCU3cDDu4PM777So1LJF3QheiEN9HEi6bUXvm43L3grJPQThtYeYVHC21n\x00" "T8PjKaO8HksS5EY7V67C4voXYfGuLqUbYmvkv9F3F8MnMSXPu4i0R8634SIJ\x00" "OVCWYo1YNCJ4BCkB3FYhNXSnLCLS3rmrHbkYqhUgPvdmn5DhpCfEdfYNg\x00" "-CeRl8Xl4tv9NCrMnZumVhOsABpsUB2vJiMBSrPZdD3Ma7rkUaU8CFhr\x00" "-DKWRmUPmR2FuKnhOO65DAmG0tIslCfELuM4SnaODu43GYslNmc1d\x00" "5ddBdacEoRbT8nEVOuOrBWH3kB4JUp5gM0fUojFpKN03kcfrHm\x00" "BrPO9eUD1lUG1SHT77Fr9kQJTNblMqEbYHHoLNXDTP7kBl\x00" "CTAV49CLQPCuA31oIOkhUCv6Sre7T0ZRNsrNBT2WNo\x00" "-5S2cdPTj0avTZSgkR59C4PMkeVQEA6C38aTMZ\x00" "2OSuTOOdhWHf1A7KobP1KKRkJBRITg9k\x00" "-QmEl4P0mfLGMIgrXhIY7V\x00" "BnSJS2jH528D3HoY"},
  { 415, 10, "CNUPgritrDJdgOqlR\x00" "-5hnKToTLRQ5JsbSXfWqgcR1OAoNc0AEZs9PnjVaSeNVUdCFbX\x00" "KMvKrcdn547dZScDrKHsI2afiTSAISVTTbenOqjQgY9QqCf\x00" "-NsIi3pPo52ufMNXkhfLlIbipE3N9OVZCWBM8koBRNnFZ6\x00" "9JqSb44tb4K0tDL46JVdVKO7f8EMMPH3ebnHsRC2R0S\x00" "4oL3sQ6UgaJX9kgZ3Lc7KSNFDpQAlrft9cXeHKoN\x00" "-91l8RnXIHWHW2rkdWAITnNLjnA7D9eoeJ1t7\x00" "8LLnqLShDrd6HcW0IOoLselLpBdQNu7q\x00" "-Nj9kOD2DiJZLTuMfdUuv\x00" "M1lBKUj3UJSAcSVa"},
  { 419, 9, "8VtpElNPR10gSG\x00" "Vu4FuQNkjATVUcvi46GCWH2AkDC9MdU3SCrZvMu\x00" "2oCvKfAoP6vhRZrCU4akLFJFKLdOZpOAPfr9YG\x00" "-1p7kDqpfb1stPuDmn66WCSajTjRZJg2WUfsG\x00" "L0L3ZTFGM0T6CVMEnYLOv7Ek6EXbrYkBA\x00" "-GP7qG7sN9b4UqsWoOLCSQSugZOQtJ6\x00" "4vDXuRqlAGGfWDqW3EjlWMvi0AC\x00" "-1Z6iNGWmjUoD94islSiugS\x00" "TuAZvNACkiHMo240"},
  { 431, 21, "54043IoSbW7RXqg6BagrPAuXWUqF\x00" "-SlGSJE7SDOBEePF8kUaDTcdd1X6RFDXMpG4FJiUZp0oX4imO0VesprVao1cacMg64jF8lKggab26RYnLu\x00" "4fBqC2McTCoJ6QsuLu8tIgRQbJhXmDBbp56ABZ98fFV4beHRfStFF5Opn9ceJgUoEa26Tm0RZPYE0FjH\x00" "-IEQclsbQM72vIUtFPfL3XgDWjbRsSiIrlfDELEus7QZNgdjbu55QhDZYW6qrBuQpAb0XiZmt5PNZ9s\x00" "Jpq4pWUZgGQgZAEXCcPob1TRL0r3frcb93akDcTOpq73OXXcnXs6XJMVr5N8s0HimOG97llPoTcW\x00" "Q9cE5FD6d5NAK7IaJpNgTDtLuSnRK5Df30QL3VEpGtYQ6qriUekRnr98Wu6u9RJIKACdeM4Zsh\x00" "-fOZmbXOYN7sHiQK6j1n4eoMhtDpKahL8rEHXPMRHNhbSPonJMZEEV1OWUN2EC3rBuFS8Dd7k\x00" "-dXetP0FTh6spXfMlGUdd9BPZSZRgMo1nLvqE4l3sS9bMb6huAnoNtDJJgLVlApFMhH7HsX\x00" "vQmukC7lqGIeepmnCE83Ti542esZARDKqKAbMPrjMnK1UV9WJLPjWR6EAXBeNAqvkbDW\x00" "SrmX6JvU5LosYtpgs0gQf5DHShB4fdFsXEs2Z6mg51dG99EeQ5qjoFHDlUkoiqnpI\x00" "-DO54O5YcfkRD0GH8UIYAShhsBcILumluTNDVOYV2Qgf71UPds5NhOvQNgR0mSqT\x00" "-8FSjRQBpd6Cu3UoHVmie948eu8Rlf61ALB92cZm5OWo4Ln96fM6M2o91qYVP\x00" "3APLgVUltUDZukfPuX0BeNjmZg74K6sNjgrHUNJqtat7Oa2kYL0a4pIB2Y\x00" "hdS8nsrFLATc88Zk1vctWnevt6TP8hlIqKLP1XiSfsqs93ZsHHDKjM\x00" "WQ8Tl7JmUZZ912uUFGiKJXdM2gcskcBUb4WSZR6aJXZabJAdpUH\x00" "-hoLdJGMHRsUP97KYbjWoUmk3ODasv7uYjgmCnXWFEJqHNb0\x00" "TdX6jK0LsWNmDj2uL4gg8AAuGJ5XmT6Oekc5tk1Y6Sf\x00" "-72aM6Nu8asTNgI7KmqiSo8XmYYhTbW6NC3518C\x00" "1cNaaR4aAHRRBj56obP6u2PsuAd9ImRHk\x00" "-5Qsqmb7VhjjW8Hg5H3eZdU\x00" "1Gj4VTWWXgpCh5cbk"},
  { 439, 15, "HDKqjg1bf34jFXbIoAiQB\x00" "O7fftICNPUm4efYek6M0iU8t6JLNa6sUlRHpSihP6O34B7gD40Yta4mDNscZ7\x00" "2rjUmkiK5QMINPKVOYMbCrbGqIk9H3EOPeIeE04bsWXDmJCfRbIBTPsu10os\x00" "-M44cYKfPMWMSnEqCvI64XR5hq2tdNLLKv21BDUbjMDPgqBEqHTBgcbjLmg\x00" "hRnLMRJGnA8m2Ctfhcfot56GTioWuUbgsctjZH2QmcHcrAqG5NMgRpL3\x00" "-dSRspMH4mhnYlJDrH9HTW3qhqVrr5bLHPl95WCshtW5IAp81RNhhKZ\x00" "MDIFJ1fYuIZmUKGLt5KfP5EiTuSt1gdntl3X345WT8sfp1HULS9s\x00" "-4fjIZ2iTXp7f203Q6furpKOQrJI85Bpvv9JQESMnDpQY7iIiuS\x00" "15cBsPjkZhVSFlEUJ3v3co2RZFNeIhJEqFt7clnfkMWIEHnP\x00" "-aEtZslqJTOivOq5JaRkUIVrgdTeCGRiQWODkWbVZuntl\x00" "7jVknEHC1nmpfFIFvhAnhqPVXG3gFdfKoDjXrUYjv\x00" "9H6K08fUOfMpSaba4aa8jIQTEeIcuY98tP8Zp\x00" "5VeY8520NeXmXgtYcSfjUB8mkhY9FHooc\x00" "3spi9GGHveSunNLoClEl9\x00" "2KYYPQA4mP8QYTVeX"},
  { 440, 12, "ZNg9OYE7nHeNSUZmOu\x00" "16aNFMLEqZGHSmisOpXoLTJLCJe6EUsZFS5fcrSfUqJ6lZT3aL31o\x00" "rA31YGtHjILCdRNGG7itWVbR7faq8HspL5D7EH1RbHMCev2GUK\x00" "-5KI8rJrkuaXLkPPOC5rLNpL34Rj3RtZlHQ8nXZKcr0RsGtqFo\x00" "-38o4P2ppQ1E1ZsMrpdgctC9b5W1Q4OjKdW187Z5O1LLeskC\x00" "-anKcQ3ZIPIQXGSFZGWD377hnLSEoVIX8AAFJeHBkggXE\x00" "-1oDeRR2JWdAdWq08GOQ8Z3rNP74EOZMcmsZfu1q0iO\x00" "-4TRihUcDo2h1PKTWFEZ9NbMRMm0VvBHsTmW5GGO\x00" "1sh6i8Bh67rneLMXuvOB1Z80uoXcGMW793Ie\x00" "LKLcEup7DUqIViYMJrAhUO2sB4Lu9k\x00" "42WHvqhnGY47CshU6GaZttpdM\x00" "-2VD12l9paTSZlAjhE"},
  { 447, 14, "mlEEGR1q25KICW4kKlaET\x00" "-QPIcjEYJgrtshPa8GpQKGi6Nh3AEISvNvGMGVBlBeZFZ8sPXcQBmiW6DtGI7ET\x00" "G1ACNZ3U0mMIlVW930Qdpn4DbiFAO5feI1eBGPluphRmao8Vs65Sm4CrCeET0\x00" "-RZbkCaUiYaO4s3ZSpbC99kltHKesTtFlneIH3WvHvSQi4nCktYSVtteqNcX\x00" "HdHmNQ9LZsMvei15QorKHNgBbS4XRPPQ0UhQbfUuRNaPDX5t5UEvBbfea\x00" "-4bunEmm1bLOflagX828L3GaIiOBGYZtVe0vYEDLp6pJAkXbgovZtN4k\x00" "Sv4LK83fFhueDAWqaTOJHH7C8ZZR2Srjdq3i7Q6QQrL7npIYqamc\x00" "-iIc1hW9ud6GUG2AoHD7XKT6EtelXf7oIQg34MHoUvDL33YlQV\x00" "46g8T45vuJ4vrFlKpIrblUHLoG7HgZnnKNHligr1aG9IFeD\x00" "-41fdNTf5X3YrFTiUNbb1YdN7pDg0JFbrhla6c8VgGb2\x00" "1uNCIIo3hLe0blqlLFAUrrYWDYuY8Y2pmvTLdT8\x00" "IHkWnsjQDr3YVQvDAXmJHrBoTWnE8WY8o\x00" "QgiUB3qefIPgn716KEXFm9\x00" "4G9OEAIb4SGajc6N8"},
  { 455, 20, "-Kq0CmOUDbT7LW0DuZLbKsfOt87H\x00" "XY0AifeLmd1OCBApZipW7P06Q1llVRSGFR1QG9VIAmgGJ6gHhJnu4KKn7RgukeATg5Zh6TT830NAlbLo\x00" "-1Ea1EFDhqJEbkOF7LarFfb266P1iegcgjp4llmj5IQSh6NVHKJ1AHC2Y2HKYo7BXODZmsulakMc6lT2\x00" "1p9XHiZ5a27IVa12In8O8hI1QLk0vbaZ0B0uksgZk91P4vOZ3UeaTr3hX0B7DSn3VPRpmIlDKrCGm\x00" "-1di3lk5UsEYmPQ9F28viBbonGnKES9v2pnF7nlg10nIEVd2dKR6YHuMe5PNgMRE4hnS6PZ02lcm\x00" "11PoJ6LYgNDQDHLqYOQShgG5PAcBXPFqCcumpf39vt1ebiFXI58Nhpp4Ne9PhSpuD9iBfjEUHT\x00" "a7krgWbGgZRXEDp3aGFCXDoZGSQJEmR2dUVSDXPJi41ZTcW1eZ5bajlJtRhTPatReJqcZpg\x00" "jPNbUZhbiDf1W6TJQA0APe3un8Km8YKc27rVtgQuoQ3jF0W3JAO0bJfrNH3teCuH6WqU8\x00" "A6TJCgahEUdmDofWgoB9pvg6gAKsPPWCjUnnQ0ANRRBFXpRB770Z83QhGjA1j0ucSRE\x00" "-2GBdmjWicGpnWSHttGqhvsSiijH5EHFIq9HRHYDd7944LsqIA5VblfJKBEoVu0Ce\x00" "-1IfH7oEXEndICfVaZVePIZNN1hZAYPfOh39YQesPvsfEVumvXO3D0p8jqHrBMe\x00" "CrcIVnWkTjhlOnSeDC41VYLEKfqYqGBCsSmPIjXdFEjoN0Bl4btHrugGJM\x00" "9vdUtSNteeEqSc6QPoP5aN8uD6Ml51HfLojt4ShAV2XmbR2pFG0raJYc\x00" "32LTfquLjaFbFOFpSGf06FKbpfpU73jRXf9Gbhe0HvX1rs72bYhSZ\x00" "2TZ8ABkv8RgHroFfFQbtFTYvL9o4Fo9S4uipLDW4XjhmrAfp0\x00" "123dbWoISBGLmtUsI89QqhDcvVY8cchuWstvVTJD0rFDn\x00" "-7n6poaOaj27edolZ4sIfLN8dLXv94duUdpncUKE\x00" "11jT79B9IXmNe7FTBfXIeg5kW3nfhpP44M\x00" "-10nthj0hmeTM09m5mlYFvcF\x00" "7gQYhvH9pAFtLTT1R"},
  { 471, 16, "CQr0i7BpYN9Va8NbKPWJBEBp\x00" "-3XNdMokp6v4HLpVCOLDsm5hbf6VuWPQhIWk31HTgCWk815RVtfBuCl8vKfpLJMGINOLRn5\x00" "QunrXLHd4L55PbLo0suLvlM4rGsUCoJHsuY3j615aA3VYWgJEuRo2DhbIvZUFl6daJtf\x00" "-1eCgnMvNOp5B5FvSH2E11Ka3BMS7OAcGFAHTiBdn8s6Bnb25gbjONAeH4ZglnjAckVr\x00" "2cB3s1IBpPvJNarvMIDpQiqRQUToDvR0LiOvtT3VcdjHitUSkjmVmt97hkh8Gq65e\x00" "-1XAjMgFlT41MkAIkSNUIsBTm9OFmGcVdfDDETXW7iM8cWqhCWd5KH3O0hb1L5pS\x00" "b2GFaXLnlAB4Lb3sG0h4Q3gZqUhg5h28lWYXnXLkDt72bElr6akYrJPY2PDf\x00" "-3n6PJFMgT00Lrnc8iYK6WG5WiB8Gqlga2Js3WkLMCXaji3UT3mjiTIASl\x00" "3GRZkjPWeNOg1XkKbKdSJkfisGDUol0nIAjFKBoRfKbSrZrLpHvc8kj\x00" "SpphVXGZJCZqgvlTTmfsLU0aL75gjBBYH44JODJjgdPQiigKd8P\x00" "Ac9IOfmjl6vjrfoZKLSDhq0qGBkImb4XYkFHGS7f6RvrKBA4\x00" "6cZWbjIoKA3o9LntsO19lOC3DLusm4hMHjMoZaJ2g3AI\x00" "250tWTcYFKNruqq5XtGCJsgdrbmKsE2mjrQTIhEU\x00" "AcDbIteq4OOiPEtKQCIXDqjL0Cff06nqNA\x00" "4oGlth624OOa29ImSXl1CT2\x00" "Oo7aPUkRoYUeYknFU"},
  { 479, 25, "4p16uBP5KYZvG3LaEKpPIDjWm5TEsLhT\x00" "-FWIWUn1fYVtBiphNhCt7h9qDWA2X1fegLH7uoiDWrRAOETt76ISr9mHJvtQs9kVmY2GNcG90TM2enTlcC3Zbb93aET000\x00" "-13ccZ4veMrM6SRFVYNrJOl6TvDWPQDnOTXQ5doZObe8Cdo7KFZRVGnaAeqGLde8JGOMZaLqGjG8RnKjQUOS4aAAohT00\x00" "Ae6Z5blos3MUH1JD2lJimpq95pqjWcrhHmKQGESUm8R95na3RV7PTuRP63WPZGrUHetaqLmj3S3hU00tuHM51PSaET\x00" "2Ier04SnqHUS6iTRoYDiD0dL2hunME7CrChScPc8ZXvtIfPped38oa999DZpP8IYMnfDnSYMMM0QrZNJa5CjDT04o\x00" "-DGm0MWKoRgCJ49TceL4kQctajC4MjrNrhWVf5b13ZiTi8UYFgP9MGPsknCBG1fHa6RSfkOLZrAJaIQCWQnFrhn5\x00" "WSWcISoUjbH6madGldt1V7flTCbGHA9OlXophpedi991T9KDEYYUbq2XMYbZZkLnqK8laQPVmoU26Drevc0no\x00" "-sYJu0n2Im83UjVS2f4TN2ipvE1c65NEa0S2mKl3NgssWjWh5BUKXsVYpX8vE0q3qXuFhbagLg6PYjD9H1Cq\x00" "umc0AR2e9HT7otAn6MeFdgaMtOoWI4Nf0hDRs8hum7pb7YIObjum5SLESgpCJ918Ioh3h4qi4ijiSFYQY\x00" "-OcvllK1hTjaLR7Ufc2IvFCItfKlg7PO2RY6igjXsHF5ibT9Hh0VrH0DZA2MCZoUPNQcledLIrHtbhDs\x00" "1vIlfjavqob0vRffm2g6l97diEGEJ3So33aXTXptXT0LD4H0LZEgPF3Xp7hsepHVgHRU97faAE1mk\x00" "RoXXJqt0DNGAj50W0EImSJf9rOYQLQmsXBXmaanTnmhhcqsdNF6Qsql5CoJN40CgGVlf4nthrI\x00" "FnhpFm9LlhVLDTLllQ1UccADGIEoD8kO4sE6FDvXGWsjWK18IIVdPWiOt2uiB3um5KOvT1P0\x00" "-7585VWGQo2kUle90JKnjm8EdOnCNn7iia91tsC5OMBsMRl3L8uAWF7qEdo8dl8ClgpCKI\x00" "-18Zi9L0ZJlcFNAuqJTZv3NlO6psrWfPAOPrhr2MjKVcjN846GnstsCeXfU4gnv5JHOW\x00" "REYV8J2aXREPh3CGrOuJoeBjTbM9bYf8uJ6fkGeKWrmrEeE6XRVRb6tCLX5LS4s\x00" "4rDGAg3DY0lWA6UXSciQe3ge9ETbIjGQ20LZAbRQ9p0KGQjAKSdPlZAnmmI8B\x00" "-1g1InXReuKWNGVaWFtQALgfJmL2EcEmP8jEdD1ai6tVDt1qrHeeYBCIaXg\x00" "DQPu5bB6O8gOqMlTYqSeBDrQSPskVc6k7baeJJakRSgODHRKsJ1otO\x00" "-77kHq0kkQScSknXPaVM6SBDH7O8DeXQHGNgafe2ZeprA75m46b\x00" "1pdQNtCCUe56mkgDCQ9qLjoD6uOOXKYcoXr6VVR0lAa7AB\x00" "-80YnWoPiQkKkrfbH3fW88emRsa95G6DnalH9pYIS\x00" "XifvNIupF0uh2tjRFkY6I0Tc7EPQ67LZCj\x00" "-ASoCIpElkSj51k9XB0dsGCR\x00" "iEs89hoU6TcINYsEI"},
  { 488, 10, "1lIu0R6UUHlH0nGUSW\x00" "-7OgSNJVFUudWY6ATH3ZebqLYr1P1Vs5TY9oIUfOppBBeIRDA0W\x00" "BBADCKnT53cvZ3tcVtjXnDSCFk3lncNt80mXSMS24tA6Yonk\x00" "oKSnKq7VuriE0hDZeCqMfX27ruqSk0FW38mIvIm7RF1su\x00" "KlCbZGob0cRXNeB8HUB4TX1cPO8PVYvLv36RbvtNGAK\x00" "b0DkF8cUDWtrdgdv6fIAYSbHfdVdb5Io8u8p1LTY\x00" "184UUjajia9cnfWXqR0NFpG7jO1HY0dMcjRQFg\x00" "3pqhsqWvSDnmHsXfA4bs9hkTRvchNkPM\x00" "DIK5koAkcZ8itlQbpUJcbWch8K\x00" "-1QAQkQXFJVqaq9j7Ns"},
  { 491, 9, "a3614GQ3HG2q00\x00" "CJAZNiNdvg0GSihMvAmKHpvG7M1Cjf9luS000000\x00" "1ZRBAClZIVZCW6oCqKFkSPWkB2EROliUslIcejQe\x00" "3c434UvBOKgkRBmtPs012QUS8JRuXEe7aghU9c\x00" "5iTV7AsHRRugS6U47jnV8q2N3AcN5TqVFWvE\x00" "aDuTONCrHj4NuG6nMoonA0ed4hHApsJY\x00" "AFdvUQbpImNvlk8UQ6arOBkAPaHVc\x00" "-Wf06eP6ICOFoSkK7evG9VZk\x00" "1k9PYKGMEQ0DiVcHDA"},
#endif
};
#define NUM_HILBERT2_POLYS (sizeof(_hilbert_data2)/sizeof(_hilbert_data2[0]))

UV poly_hilbert(IV D, mpz_t**T)
{
  UV i, j;
  if (T == 0) return 0;
  *T = 0;
  for (i = 0; i < NUM_HILBERT1_POLYS; i++) {
    if ((IV)_hilbert_data1[i].D == -D) {
      UV degree = _hilbert_data1[i].degree;
      const char* s = _hilbert_data1[i].coefs;
      New(0, *T, degree+1, mpz_t);
      for (j = 0; j < degree; j++) {
        mpz_init_set_str( (*T)[j], s, 58 );
        if (j == 0) mpz_pow_ui( (*T)[j], (*T)[j], 3);
        s += strlen(s)+1;
      }
      mpz_init_set_ui( (*T)[degree], 1 );
      return degree;
    }
  }
  for (i = 0; i < NUM_HILBERT2_POLYS; i++) {
    if ((IV)_hilbert_data2[i].D == -D) {
      UV degree = _hilbert_data2[i].degree;
      const char* s = _hilbert_data2[i].coefs;
      New(0, *T, degree+1, mpz_t);
      for (j = 0; j < degree; j++) {
        mpz_init_set_str( (*T)[j], s, 58 );
        if (j == 0) mpz_pow_ui( (*T)[j], (*T)[j], 3);
        s += strlen(s)+1;
      }
      mpz_init_set_ui( (*T)[degree], 1 );
      return degree;
    }
  }
  return 0;
}
