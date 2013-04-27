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
 * See NZMATH's great site: http://hilbert-class-polynomial.appspot.com
 * Nice because these are *really* slow to calculate.
 */
struct _hilbert_poly_small {
  unsigned short  D;
  unsigned short  degree;
  const char *coef[4];
};

struct _hilbert_poly_medium {
  unsigned short  D;
  unsigned short  degree;
  const char *coef[8];
};

struct _hilbert_poly_large {
  unsigned short  D;
  unsigned short  degree;
  const char *coef[25];
};

/* Storage,
 *   1. The first term (x^degree) is not stored since the coefficient is 1.
 *   2. Numbers are in base 58 for a little more efficient storage.
 *   3. The final term (x^0) should be cubed.
 *
 * Storing as something like a stream of Delta codes would be more efficient,
 * but starts getting really wonky to read (simple with Perl+Data::BitStream).
 * This 3-struct system works, but it really wastes space and just is a PITA.
 * I think the answer will be to merge all coefficients into one string.  So:
 *    { 23, 3, "6t1\x00-7nBh5L\x00HpuY" },
 *
 * Split into three arrays so we can disable the third one to save space.
 */
static const struct _hilbert_poly_small _hilbert_data1[] =
{
  { 3, 1, {"0"}},
  { 4, 1, {"-C"}},
  { 7, 1, {"F"}},
  { 8, 1, {"-K"}},
  { 11, 1, {"W"}},
  { 12, 2, {"0","-G32"}},
  { 15, 2, {"-8V","ujV"}},
  { 16, 2, {"Dc","-1Rua"}},
  { 19, 1, {"1c"}},
  { 20, 2, {"-FA","-6Rh6"}},
  { 23, 3, {"6t1","-7nBh5L","HpuY"}},
  { 24, 2, {"gC","-OjF6"}},
  { 27, 2, {"0","14uk4"}},
  { 28, 2, {"-17t","-1Qu3Y"}},
  { 31, 3, {"YUB","-1VNXnn0","3SNMZ"}},
  { 32, 3, {"-Dd6","BHMMIHk","-4ZmS0"}},
  { 35, 2, {"-1UG","AOYlo"}},
  { 36, 3, {"hIW","-e4NHjVE","-DWtBs"}},
  { 39, 4, {"E74T","Ejtqjjo5Dj","-BGss2JC","THAk0"}},
  { 40, 2, {"6C0","-bZjb2"}},
  { 43, 1, {"GW"}},
  { 44, 4, {"-EDM0","165FbrvvfQ","-GVDphf20","-1fBigq"}},
  { 48, 4, {"0","-lZfXocuCeW","1EIZGfoW0","-4IYVWq"}},
  { 51, 2, {"5Rk","8PbYPs"}},
  { 52, 2, {"-OZY","-ATQI88"}},
  { 55, 4, {"-2Jh3T","6qYJaDWKBIoj","-9SGADmpf","K0mscT"}},
  { 56, 4, {"1qc24","5Du6pXTKXDQ","G4lPRoiLk","-OfJfBs"}},
  { 59, 3, {"3eHE","-15ip6mG2C","k0Qva4"}},
  { 60, 4, {"-2Jh3T","1KJk5leif8pQ","-s1FTcL4nT","-uN8Jeq"}},
  { 64, 4, {"-1Kv4K","-1Voj5megCI0u","3BbSdf5OCe","-29G4n4C"}}, /* NA */
  { 67, 1, {"1X2"}},
  { 68, 4, {"-BHCLs","-ChJCUYM0Oa8W","-ACDhAfc5WC","-4dTq9Ns"}},
  { 72, 3, {"-151A0","VeANIiMCZo","-9rNlsqu"}}, /* NA */
  { 75, 3, {"0","edGEDKNbo","HB1OdT2"}}, /* NA */
  { 76, 4, {"-7tp7s","ds1gRfCEUL5Q","-1ZE2GoMvE2C","-KYXpqrs"}}, /* NA */
  { 83, 3, {"fvBM","-5XurMpAFMe","1CfGc048"}},
  { 84, 4, {"-FE73s","13GhPRbcHIoQ3Q","-D8TeTGLeVqC","-1PuUKQbs"}},
  { 88, 2, {"Cnq4","-2nKXdrnA"}},
  { 91, 2, {"-kXE","4e6dksHk"}},
  { 99, 3, {"-14l36","2gMK3PaZjFk","H26MSuG0"}}, /* NA */
  { 100, 3, {"725YO","-Bdtag3q7fW6C","-Jsaeko8u"}}, /* NA */
  { 107, 3, {"68p9g","-FefkejoHLDo","10jAlvoG0"}},
  { 112, 4, {"-6OOIF3","-1MjTfulSIiR1Pku","38QJjC5NUXgc4","-28TZTuSZI"}}, /* NA */
  { 115, 2, {"2YdI","3JjG6VK9k"}},
  { 120, 4, {"PsOOH6","-98sjoVDOsDZje4RY","I9goh3tpAonbg","-6pseSdX6e"}},
  { 123, 2, {"R9FU","AXH8BE7jo"}},
  { 132, 4, {"I1YW3o","3KQOvanr0dDCHhXWK","-3oNODN7vQt3OTc","-avJLWgRIm"}},
  { 136, 4, {"KQnFQu","-6SIEYD9K5kLeE0nM","8hm7rkIve90naa","-15bgGYCTQ8"}},
  { 139, 3, {"3YY3A","-2776mZ1lrLoW","1b7iuZaSC8"}},
  { 148, 2, {"-HYdp6","-5Je9j2Xgho"}},
  { 155, 4, {"otll60","1E7dPGqFP3RUopWb6","-Ue3lPS270HsdA","D2eeGn7VOC"}},
  { 163, 1, {"3GK0"}},
  { 168, 4, {"KkVJ5gW","-4kFWDuF3h9cRds5piFY","1B1R48gXEDJ1DWRg","-174uqMmNYUK"}},
  { 184, 4, {"Chv7iA4","ee3h5JqplhpKaiOmlQ","KMroR71hhP8Va4Dk","-7Qtiia8sGns"}},
  { 187, 2, {"-1MHI8","AVst1CXO5s0"}},
  { 195, 4, {"-1ZlNOi0","6MV3ZsjJqJPfC53cG","5BXTa59uoTqSFtk","QBE2FFmQr60"}},
  { 203, 4, {"8JFLMcq","4VSOGETrtr94R5hUtvQ","-H24gbf53k9PmENY","15meh7TeBe7o"}},
  { 211, 3, {"2cA2K0","3HNW7SKtX8cK64","2aqdq88Q2iPs"}},
  { 219, 4, {"Ca9R3no","-GkL3kW4Ulct6JlDH3M","2sQmjBmX0QeMikuO","6CGUuU1ufm9k"}},
  { 228, 4, {"MJghYEHk","5QuLi4T8DiJXAC4Pgpa6lo","-8Mq7rfNsESXkBWqLiK","-FvXa0mQt1bHQ"}},
  { 232, 2, {"bR2aEq","-OBftZLKGhClo"}},
  { 235, 2, {"KBbcu","WsjthVrmHvD2"}},
  { 259, 4, {"4HKolu4","5itO38hCf0fYtvYTTQ","-1HUJT8j7aAI1Eh4C","6EtIplr7N4ko4"}},
  { 267, 2, {"CJiIdI","DXhAWVqTa4jrg"}},
  { 280, 4, {"sma1OCM4","-17lkqPPlRWOWFREZI9QAXkW","5T3O9tSTjGPDcZ2ORu1k","-ke9YeIBdIKFJs"}},
  { 283, 3, {"8rt4oq","IafDePEn7j86DJE","13mLN29jn0qr0C"}},
  { 291, 4, {"CbZIWDmC","1VApsBTVdE3bukZsffA4u","BJfJLu3Mt27umfZ1jY","2DmEEJZdKiNZ4K"}},
  { 292, 4, {"-Wl53gO00","4DG4TOnStsgQgYalc1M000","-TGFUT0fT3auilPYTUDts","-2QJrHb3OdSg0U4"}},
  { 307, 3, {"VdElse","-Huht4Y350B17psG0","9XRpHcljImVkP6"}},
  { 312, 4, {"1ZA88JL168","-1EUvQFWdITnpVPG65pQPO6Wojc","7Z0SeVfPo9TphHhYNfe8u","-Eu48DeLQ7fQVOq"}},
  { 323, 4, {"-3oZIk7aUi","22CMQsiLI5oVtDDQshnB7Mca","-8uDF3jSHki8YME8HaK0","dRJXCXvE465JpU"}},
  { 328, 4, {"3SAAIr6vM","1Tj4Z0HQbW0DVUKdgfdfsslY","14WTTuH3qpc5HHb0EaYQbc","-13067FZpD7iT1O4"}},
  { 331, 3, {"A3TFLg0","MS1q2iSL8DbiGogiG","1L4oIO79o4s9ee8"}},
  { 340, 4, {"2gfbACNm4","2l6RsMcuLPJEUMYpgVTMKmMpQ","-53rCUKlG0gM5HXS3MIOqAC","-2t6jKmn3JcMlS5s"}},
  { 355, 4, {"2SkFR6Ii","-7TFU95gEeEAqhFUk6j3k","7AMokgX8COdpO3SJom","APeSpnX4IT528tg"}},
  { 372, 4, {"kX9u02jSJQ","4GRHbArCCoBbL9EIo9O6oG9jutmq","-4hoKE62kAVJvYXebsQWuSmm","-gLtRvR7Y7hvqcEi"}},
  { 379, 3, {"17O6B8oq","-2BkYo0sbfO8atVvdiMC","1GhBTI0M447AoWu0"}},
  { 388, 4, {"-1N7hMMi3Rs","1f8EorHcJdKU1e17aLbrQCdffQ","-WqKKIC5cTNJVp7WmAQ5tsYC","-2bnQr23K5tBTskDs"}},
  { 403, 2, {"-7Fqu32","8d70FAdo0cksb2J2"}},
  { 408, 4, {"5TY4Wg6aAo4","-E5HYk57B9nlUuPi3Zk4QCY8bkdIWW","61GAJTP8NLZjTi9YiWl0cmLk","-Cl713PHOvtuBQWfs"}},
  { 427, 2, {"1NnBaoq","tCCq1M2062uATRg8"}},
  { 435, 4, {"-I7OoKoUi0","KNjHFlSbDPNDCejspt6o1bg00","87G44RUT4VQf9vGFTqIM4G","1h7ov2ZFI1rEG1Rg8"}},
  { 483, 4, {"-q39DXCgjM","1IQhPt76jMMG22ZFD7VVrS4XgO","FHoqG7RUTG9QhWsVFJb0mPc","10sSoIWMKuknj7rH2G"}},
  { 499, 3, {"1Hf5S40h2","-WcvR2fG9N1SbIigqXgB3c","39DeQQCh7siagGlDFk"}},
};
#define NUM_HILBERT1_POLYS (sizeof(_hilbert_data1)/sizeof(_hilbert_data1[0]))

static const struct _hilbert_poly_medium _hilbert_data2[] =
{
  { 47, 5, {"MGb3T","-AJZrA0OkAMQNq","BocU01cGqJt","-4ULG9Kcf","3PTvOX"}},
  { 71, 7, {"NgbHFId","-7f7oZ1NsODjsjUM3KR5","5NJarvPWDgsdbrvhne","-2qrWEVPTmTbNGfm9","19poArmToUCnJa","-O8LJYUSN4","8DnkhuN"}},
  { 79, 5, {"Qln89T","-KfpJMZZ5qoLp7UPq","LJU08d2sFJ1Gmt","-nfSq3f0KC","ZAJcEU4"}},
  { 80, 6, {"-1gY5HRM","-POt9PRcqEK5KfYXPI","kI4WtSjFTrB06Fo","BvjLCeA1UH4uO","4c2ZtkHrF0C","-ftNAu9M"}}, /* NA */
  { 87, 6, {"LU86us5","7mKFR07ZhmVDDCQLruT","UJf9U40p39vXM8QXp","5ksD96BtkpuvB8q","BUHelecuva","2Nk1nlVf"}},
  { 92, 6, {"-BIZcWQj","UNkDCRkvBlDiUJLhRE","-1F92JUPM9kvK9AGNL","-AEdVI0FBUZtsfM","-1fWIhEbsqNog","-5UdLGIKm"}}, /* NA */
  { 95, 8, {"3fVLhjbaT","10Ga3dZ43XFmeh3UfEvgQuJ5B","-2HV1BjMGNhYTivh40CBUV3n","1q3D78YJklMMHASPjo2mh","-Di8dRLIUktoVVMV8fr","1N1Qim3vW9oRcvq9","-1YbdjZiFrik","903vWH3A"}},
  { 96, 6, {"-6MvJPGK","Bqf8fKXsX0VZlhWee0","R74cSutnanG2rtpE","-gaaiAWS2qTvVYu","4VTE4cPZtB0C","-AX6938FI"}}, /* NA */
  { 103, 5, {"838kCLD","E2KUL9h0ammE7VNUfd","HRStVfdoeSFG3Oam","BTPqi8ZZkR","VmQaaiYv"}},
  { 104, 6, {"AXtEie0","-R3BIJDCSIVHsrTqn0C","1Q2X8iikb9ehTu59KC","6Ku5shcjY8CdFHo","TYcMnDSqe9dc","-b8h3v6iu"}},
  { 108, 6, {"0","MZQk20gqO6i5Td8RNhk","-AROHlhWPiYUa5Sr3k","L6h5VfvI66oa83o","-1FrtgvCVcIh8K","-1AMnslTcu"}}, /* NA */
  { 111, 8, {"NXQCHjXVl","-UmjCXuQTTqH5BotolS83sHO8n","2QRe7J2iBNp1UQcRnIoVIE0S","-2MJnc2U1MT4bJYnqOsclUY","s8jE0RnLXXFQTIUXWnI","3OTIrXRGOYC1hIL07","SQGQR70LO02","1nHOAGMdT"}},
  { 116, 6, {"-26CrM7A0","iklJdQupkaaqvshAm000","-1BsUhl2JNORU0i1hXmq","-2tHAPj8nK0281TQW","-7aTYeGUoCMQ6C","-3oG7pNQXs"}},
  { 124, 6, {"-PfoJfKT","-OA7onDHGMGO043g6QMC","-2VEe537dBheLP8QGSk","BeDJ40cTX4Q32Jo","-gTBSLC5kLWf8n","-CANfGFEKu"}}, /* NA */
  { 127, 5, {"35a5EtM7","-K62lFUJN30ZNf4PAk5p6","5s3PepgPfV4dLr2eqa","-1D3bZGBM2LLP","IVmGdqC78"}},
  { 128, 7, {"2U9j9rI00","-3tHamsJB6LtpRBosoq000000","-VAeAPOBGYaF6oIh3qQ0ZY","nSDG2ArLrqYHlAAnsEq","1AeXKhEP7GU0p65Ka","1e9fpsndYuifK4","-LINP4RQho"}}, /* NA */
  { 131, 5, {"1LtNFSa","-3dPtKcJ9r0BLh7vJ6","LvOJUv4DD9RYsQa","-QntoLOVb6bts","WEenhFL68"}},
  { 140, 8, {"-4HOp3sROsi","5ppD14JNd0CPOqvPrMj9JApAC0G","-5a67iSG4htsqLtDQitLCqPjW4","f3mWKesahOsRUcKJjQtYn6","-T11nhMjOlrSC4tPkaFQW","-NQSnv5ql32cjhU74S","-JX0WFLNVdekpom","-1ocbJvaT5Q"}}, /* NA */
  { 151, 7, {"1vpjRYh9E7","2D8f61QRkqN4iLKM391p0YNbup","LULAIWnbZ7Gfga8vrmDNP0pf","6R2TZvD6oOhsNOIcKK9gPW","13aL64SnFK8J405aPs74","iHbveCfE0eLJ","7nIKcX6bbG"}},
  { 152, 6, {"ZFl7JGpg","-2CELUvC4aSEWaofkFgBteiW","Evh2IUXurN1dD32r5lrGa","7BtQMh6aDLrjhYJoEK","3SV0n8D2kZEjHRs","-8rH349ZARo"}},
  { 164, 8, {"-CiQEj5nA00","1hMAHYBpoOHjV7PVYCnbm4000000","S1rFHRkWi9L58WOCYMeebvlIO","-2GlWc350NgLhfBVgKrjMjKAa","-6e2WllAEvU0gtPFqRk0W5U","-1tTX0Af8Mv97l1OuDDQ","-XCVtsD2335gPoiO","-du1ZfckCm0"}},
  { 179, 5, {"1ns1DPhc","-7IHVPHmSnvr13nXgT9IW","2ktcYPFh3PjWVcl5TY","-QADuR4FptntZhE","49e2mjsiArs"}},
  { 183, 8, {"1ESCgpX4EZQn","-5kg8Q7iWDW9rg1a5jq6718CQle7GS7pf","DDM6jJVHspS7OpnavRKVi7RSus8qjf","1Ym2O7uUGWUL2nn2LN9NkvdoSdlh","5huBNPW0S9pDWu80DO9jNQmEU","iB0caVWnmobaNhBiNrIie","oAt25HonUc0rmP","6bWJiCdRakr"}},
  { 212, 6, {"-5Rnrj1uqUq","5bh9m6vFQbfpWPeru8AOO8FSk9k","-5E0BZu2Tuj00iASuupusW63se","-E9tmpkvkFgUttLdf3mfl6","-cKPfVfkRANFIocllk","-2sKE5kBLaSuG"}},
  { 223, 7, {"VskAFuDSFu3","-bhFYY00NWYSuAFNRVH74p0V5Ro1KcgO","c4OuPV92igN2jcBc2pl5tm6WCnI7e","-1P3jvPQIhUTs3vMTbnL47pPXDZr","1Vneq5F7M4OC9hKn5bcEnSi5","-8mPAIbPLFndgu81","9RkJeuqfJhMA"}},
  { 227, 5, {"7jg8hkWG","-BLhZTYku5n2OqnOT94Fvg","5eMGnioA5vXBQrST62EC","-93XqbGS0c8dII1K4","ENmbro5OKQFo"}},
  { 244, 6, {"-GfeVU16aO","-BaHgiJUI14PuI0JAY9Saler60","-o74CrA61RJABNRRlOaSY78a","-ET3N47YArXT10kONnau96","-1dh4OWNKLGOOIMIFE1M","-1O7oUoRAS8vJk"}},
  { 247, 6, {"-HC0EXjk25JT","50OJOemJ58mDr0fu1apiqoGITtvZsgH","62ILCt0avrWc8gSbtGNcsCr5KOdf","3cF9h0tSmItscF1NMQetsoffo","38FvjU0kEeAPOf59","1qtOv3ctBLqe2"}},
  { 248, 8, {"4MnQP4rANclg","-43jjYSPOWdie1bFs3qFi66pEn1XCNqsbs","3NHsaiIUMAhuFoYOs7riPdFIrkbMtSi","-19fkqT8GYPk5l0l5PhKKGHUsLXu40","53Xnm0nSpCvc0k0QgbFMb9Yi3GK","3sdaOThqDC6NXW9tnARQACW","34CMDusgo6MXtqNNE5U","-26Z9BoF900hT2"}},
  { 251, 7, {"GjSBabcem0","-4La7u8KqUlg9FrJJbq5ZVnrlJ5bY","3O7iVpcWoJavD0TSYYP5Tj7tFM","Ci51XDbftsHrYv6nTCnZfSS","5fpfLI8BQvjqAeNbkK4OW","-427keG5IYg36qAEGO","2nD6cJY48nOi0"}},
  { 260, 8, {"-5te2OU3K11MO","JXjQRIOhPje9WZjFLX8MLcieFCoZR90ANo","-IVLRltlo3N21joE8XWniMB7lab4Lt90K","AjQT1Hfhog4TRGnI2ICsYPO09aVrCS","-17WlFPQSYe06sF19oRgTZIhisMue","-YPvKdjldfClKfQKqGntPnii","-I7Af9XtVDG6Hs5ra8CO","-6q7R21Gej3DNk"}},
  { 264, 8, {"2i3gJhCrVvI0","-b264i4Wnia65WfIc8OMKqmbiJ7OrBTTI","2LcRE7NgfNPLm771Ud4H0SaakmbI4Pc","-kbKvTRBEGTv3u9rZ2lUeAmhFvHCi","-CshICnKgrD6W4PHH5eXkec4JsO","-1CAFRJP4cUEHhd4WOf2ifosK","WPeK57HuH3hSAjOIUUC","-A9qaY5m6ZK4um"}},
  { 276, 8, {"9395iR1uDM00","-OEhIM1O2BsAv8NLjtMaRgv6OQsLpmC000","9YfRgdgsdc96QMTgAaCZfS6ma04jQiO","H81F2C4TfFfNJDtMUmU2eAfmYEsLM","-14K64OkZ7NQJA0KcLmOLZtBB97U","9s6MkqqTD7h0PSP1qlfChfAW","-37FBr6CZ6IYLPL0uMU0O","-W1oDlZrIvg7Q0"}},
  { 347, 5, {"iSDiHTqNY","-aYZUeWhCDXkqRAeCJbivYmHhU","158T052DV7XEs6Zd4KlFHoG7g","-2Nn6HjHOqJmVBHhLqo8O","5Jf1SJmK7PWtJNc"}},
};
#define NUM_HILBERT2_POLYS (sizeof(_hilbert_data2)/sizeof(_hilbert_data2[0]))

static const struct _hilbert_poly_large _hilbert_data3[] =
{
#if 0
  { 119, 10, {"-5FLqclQ9lhT","4JK6gMJSZruI6cJ7Pdg4SOpiYc3aET0","-CH6esl4n86hvsmmVIaKAb6I1HL6tp","4A74F8Nd4a9Z9MpZbMTeGj4j7gp","5scJDNoirHssuoLOKD8i7FLpd","7dPoGd6ojlk5D5SdM2PLKc4","-GU2pERaCPtq2GlbT600D","Zd0WtCC9MlMSLRFah","-2l2gGGqS6Dpb","5uNpqY8Ub"}},
  { 143, 10, {"-26RbV6Y9cbhT","bGujjtraDL51hR8gv7seS09VOE9uAohT","-1fnLCn6Q1edoYuWu0pTVJ2dvGeghUeaI","4M72vgQdZ8F7QACSR5EnhCVe5QBXVH","-1aiooQv6bh0hHNJiRQhSHa7lUuFe","GVm9YL8cVMI4hFigBUSLj6K8P","-Bi9p0Tccl7TfsN5VXXv1AZ","7hatPkT2Mv8p5sHHZHD","-4EDft4GmDsfc3","2jSMZpBX3V"}},
  { 159, 10, {"EcbTogjLL7mT","-829HfjDIQAAm3s6nr1pvtf3YeQ5qcUFtK1","SHBjhJamQleYbqm9rsCf3LkkWeY3NVv8","-c49PKi7NqoL0uFCXqUcarPlDb5isjs","MHISbMZFgiVUrIdv9J0DesYYrULt","-469F1JrGN3I3jBFv8WM7LbAOcJ","Kej0Pld3DBpH4dMJWdhXNVh","802TXg2IWi6A9Y82XdnB","tXTp6Sj7IUosr","LVO9KqDkSf"}},
  { 167, 11, {"LUBo3QOB7dUml","Df7GMLrC1CT0QPFYCnM1rF3q34k6fIMCndkmR","2Kil9T72C6pVl3Z8uCu67EKJVNfhB71IqNfB","-CHeDAveVE3iljmpmj9LcqTWS7lSVBZoeSv","BnETDIEe5bBhRZ9gPnM24FUiOuZMtJFI","2JmJiYtLSLRmEKb78aCoUWNjAn6DM2","E8cEUQZAnKbBm9E1C3piv5X9qoX","-3cRPAWZWaSZuRkrImFCTotQv","vHA5eXI9tNnE9ACpDQSG","-3eCMf7pKD79Bb8","vbUI3EE4fA"}},
  { 191, 13, {"1LTtEREasoMYFBhT","6H9pAmd0TS1QtbUacQMMsdBrZ09XcTcee0lRUtUl7ET0","72DFp8aGqbLPruEq6C9aMcpFfbFIF8Mh2Ndnu4InrP","-1P8FWc2XreFNKvuRvX2CPjFcKQNIl0DE5s0ZBNqBp","1SNUSHUbuLRL2LLQ4kHoAUAKrYIGs7TqBB77gdC","-FoPpCUhdnVI80T12PNbafKsUATU9WA7vqgrK","YgOHPcjoCSXY7QWZso17ncfIGb0jm9U1ZR","-3AfP5QcJPJMNDjt9mlX0WjJQ6qcUIfdG","8MpMDTLa0I24r6tcnF58Jb63kLJ8u","-nD4oSZpTpfkcUr4mmTPa37bdo","4jdhiSRIgmIoKl0P4ih2dL","-2tQfJmmLhJFY2Hn","GcTJ9ZZFBn8"}},
  { 199, 9, {"7HU1GQu8fS4F","-1dbjpUu0GkfObiGJbRZ5F1Rn4IgBSJWujP","HUuktVs98G3Mkc5TgNKcKAWLLVRTNffG","-B8JWNHRah9i949oXMc5v2NpLSeFORi","2GGvCZL3TojkI4KUekPmIbfeXPOZ","BNrpYsQQLHq79ZJZn6siV5b2h","Stf87HqPnaPXsU5K00KHrA","FmaY88WPVMU7tO","ev52J1QBcQm"}},
  { 215, 14, {"4csTYQlgNUUXvsT0T","-83FL7iiOuWsaZbp5PCSfAWDveHoqiXO2D8OkGFFMbNjWpvsZ","2ipo2YK0AUC8780CKSjeN6PdbcobgGivL9OZhNecWe8eIp","OgED7Acb2mjqfqAefOEa5g4pjvaQb4Jq8iiKqtE81INC","tGD7aSTgL9u0tmj0ccbKI6LkBEbji78GUm1pLm5r9F","-2XCPnnoNeM1rHLdGI8Qt6HM8cepJtUnuWN4W8535Z","n9T0hBZskf0sWAoFM8joCTGgCBsPYDviQc71cR","ikoAH4kWtlLIHJDOXfoQ3v0ir88gILsvOtAf","22HSsvfgCdKhqQMPB5q2H8eBXBJMc9dEFc","2eqoiQZmPqOpjU0Y9OiB4pJpmSSiTsL","-6WigZgXHsLZSnbqMA3gcuSHEDMS","GPae2EjPCRoVFPcsmXO4B8l","-1ZsFapMhDljVodkR","43AX0DADeDrG"}},
  { 231, 12, {"-iYNKXLmYCe12lZr","mUSZY95JM5kefgOL6bj8ETthQcNosli8AgeQIsnCMAY","hv6HLRl4jutW2s271lDF48C1UUHPNEWFNBR6I7mq7B","-IvWsAhrpgXcW9E48GWh9je4vuGnahoBsfhRhe8Zd","11aOSgYHvX02jHDcTVqq6BFmeN9sZH0mPOq5cRC","84aLPJ71u96oQAlmRraN48r07Tc29AYivBS9","3mbHd6dDna6PeEn54A2oqtZePVpkYGf4Io","6jo1fW6a0XAosloIh1lRr0CY6qMTshM","65EdAfA7rud03OKbr2EretkQrV7I","8CRBm3Cp3QkbaCmAJk6JpKYF","GAR94K7TVugWrn17","Lm1u1II0MvH8"}},
  { 239, 15, {"9sA39TadkD0lNT7psl","-LktKkBdDIR2fbVqCocWIbSTaAC7crURn4Bj39uGQei7B9m61bFa","E4blHcGf6oQh08btaDNcMQpdHgho6bQrvlb7s80no8ZnvViW0KZ","1V6V9hpdQSKFIjtmLOW7gtvZlBLNJQFLGe71vgD1WMBaN3bos","-D4Xl4VWWCv12tm0cYe1Naa2CEEDeOoi71BXf4378oUb4AZS","-9ivvg4qNkFgK2m96VEk7v6rK29PQVlilYmb4fYkmQDl8T","KE12ErkYMOJP18eNu7o0hL3io41rquRPlL8Wo6F0Xb7","-6B4Co7QCNoNO03Hvv811mAlvpibi5UGoGBEoTY0Fq","iEpYmu8vNXXoY8YUNv9Yi8svoGmGvFSrUr2lWc","-111COcvKb1ILq10hQhQatmA04cW317hWNqFj","dAlp63rKuLWKdhadLarMSNDefZDArkCL","-fGXOBtmhHrJWlbbqHBpkY6BTh6sK","gJ8QnPDm77aCsWagPIDeHgus","-mJJvtuhCut3OmoVI","nVoFg7CjQOvM"}},
  { 255, 12, {"-B2G36rN3eEBuAnTP","LEnJOqhsHAn6OivOHGWgnI5mGNb8PcXDfI7tWcb0ioFfg","-4WKfa7E0LYCb2m4tatBe61Tf8gNuUXFSQi788BvQjMAq","BfpgO515mYUOXQTRVBFcM0oiqWiPVMJPG7Qk7KK1j9","-3r9UYBgmd5M9pMvs9orYEbpNQ9tcB6aJieHpgnuG","5t8flTu9FKao5IE4oX36Geor8utXeA2J8HOstX","rqruI2BTa9hNZ8RhnY9s7oKbkbffDGTaR9g","-tPEWUQj1BSnCQYGCBkTsimWpjDnlb0jQ","SIRC03cUeb5Rv7EaHoljCL3itcNRb","Hp8m8BvDe5RjqXSc2795JN9fP","6eLpPvZVU4lO6Xg76","4DFKmm640b6gp"}},
  { 263, 13, {"18K2ItVBfA7sM0DRhT","-HVFCLtEWm1GRQYKRhuAcUZhRZtmekcpo2sP2L2F9EYbi0qWaET","1N1qAeoEmTms1STQWH7H781m3kKd9Z9KIgmc3LX0Z5r7OYJDb","QN0hmNk331GZEGqb20qdCvcvhI9g2F63RrZWKdL53NcW8o","-WMr70g4qZDWuVJAdicXZ2F1S1tLBLpZBt1JctE3PQ1lk","-LpFqcPBTHd81cJQ1JGriQEKhgRA7FnQrssgBTbu0il","JH30d9TqsTG7iRgtdhUt9H9N6M3N2q54cBc0Y83Y","FZg5f7c9QvMjXb1qhsGp91G9bNCrjNSEjQYZo","6Jg4nrP9qkm4nGSmYVnFToGXuI1v8dPFB8","-31EmJqfWZMBOvPtTQk3ccmdisl2Sbe","1RE6bKeZ15LlO2NqKvZatGsrro","-Iv1V7Ba6ol9SqUvSl","9DSqrGHIYXsRh"}},
  { 271, 11, {"5BSP1K1sn3ZJ7hT","coYnBTp9ql2KqlhA3WBdSiib9VvGBfSGWRNsbDnET0","BeHvaaCp4oquU6stTcCoMktE4TGblTFbXSMeVi7cF","-QQ1GJlgAcBHpLL4UKQYQ7ADdoKhLXo2iDrIv0gb","h4PgkYUEujoLaBTbCZHKmK2NUrmH50KlYdE22","1ti74sJK49Gq2suiKovG069ETKZcXvZD7DV","T89ugk0ilDHLIuYuffc5NclmDP8K7Xb1","-JsvhXiLrr7RA7gFZ6C7QdCG4g3QnC","6mpQRJ069KNo4WEFA4RmFjYeB1","-10DD6As6j32e1UWUr","JrReaN6BKrt3Z"}},
  { 287, 14, {"9pOCMatr1jOUQiOt8ET","-LXhdsUANXHiBWkm9Kq2s09nkkou7c6SDSqkUU3phrXPKIbLoCvR000","-1iqoZltAklReKArnPPB86kPBk3pdoUceipPFtIu6hFduVYNoE4d4G","2a7pAOIAagSErIr2FhYVttdAVQ7UWpTdQO2FRiSs1TboiNlB93O","aivpGbdHFeeOdvcQEf7ilBJBYcX1eP4a2DDKpia3iVQtEukQ","-UTIgp80k1MapHlbI0KPoStY1nFB0317fE5p9nTSEqMJe2d","-BWrihh5vBGDtLP94cAm3GTT2tgePbrSXdlh2i9smr1ev","7sMW7sPh127UXXmI2uP7gB8tCl7X9aV1bWrZpv4NEo","-3c4ncJsvVISjNvYuAeaJQBBWIXkOp0ame7YRuYr","miSY5GhJJc0VtuFOkWtSPa7vg76qEIhkVc1","-AphSqO7vPcuOkTZlcc27BUs8imHbFhr","2MijpGarM1Yg2cGoV78SYVmEtNS","-72P5MCVNBaVullUA0H","1VfRPbhZ2dlPuf"}},
  { 295, 8, {"16bIHTZtKg6pk9","-6KfbcSLGDl14rOnRti7Mcl5kWtOjdioCa3p4nV","YUUgpD1NX9TLROpSs8KCMUajKtugZmAQedsO","2ds0oFJNL7T2g6u5j5C1vZAbtbLKUf3dLP","v4hiRqNNlGDJerkD8aTLKpqBcE518T","APXG3OsnlXj9WgSpeV5a8HttAUb","HbAMHWCSDK4jLk0sp","3DNq868R6USfGs"}},
  { 296, 10, {"7IkKs8I5ALHLKC","-7FYOkD7Y4B7nBHLaKiZE6Osi0kv2tGvGQfaKClg","Acv8HgIsKeetSJo6k77fTTQqovaFmEnUGchuHsm","3RnlKqs5KUeJ6ZBMI8IrUbcOJNXefIckejujc","1nbD5A0e79mvLe25cpKv5WtmcvfM9C9WMT6","-3tOiuAOt4NDW08DjPdDYBcYBh90PnAgS","UP8UU9McFUHCSKhf2MSWcDSMSW6Mi","598smd5fvM6uorDv2E44boDou","ogAstF4WqbVg1QkNqJQW","-3VJtpdgqZ4WMA4"}},
  { 299, 8, {"-1389Ghrl5ojk","9okJj16sLTPnYUI8F5LZqmkUucCJ82NI","-Ds8vaGdKkLnQfJTM8ICXh1fS7WLVfo","FcFTlpoipfRDfqa6XgcO9jEHQZFE","-1UlajfUCiH3GIUcKpR8ab2rXqC","3KV5YD6t9dJBdMNKANKTAZ2","-U629reWUJnb74TjbaC","4bnlESoAY1hTjo"}},
  { 303, 10, {"2HffhKpI6fJmCWSP","2OU81CCZOlmWIDQMib6EVPbB7H3ZTOCpaj8rXWhgsrjTH","Dkud9Dt1CI3JQE5pReUaDS6l0oWmmu0aUMC2HYXOZI0","-XuAEbR93R45PjuBuC1cZOs4Q7WlWasOFQ7NQaDFI","1qnqlgpXMZMeOLAR1pARonqgvMLHZUvMBkmWUE8","-Z2PEONNMchY0r3WcIAEuSb2Nu95sDKRlhlL","5UmPgj2eLVT31HCb3cANDqZn36CnOm2b","icd0Z0jfLS8vGEgMiZWiR7MEZpY","m0F5iVs93sYYkEOg0E","6daduG4mhEAq12"}},
  { 308, 8, {"-8rhuhJDNjtk00","1A9srbBkoaAHEresIPYOo1Guf7YkVGK000000","-1H3vUNWQjOr8hcM5VZT7OA6BC5FZmuokWpY","PQ0vGHV4ZDuh5SCYfY72odH9naBZVD3M","-5TF00dFqCQKujNb9L0aKKgecUHcA1U","-beUNmvnq4FP0BsbrGBvdSn5PQ","-4QCCiO1U6jJOMrkoOW52O","-ARV4YRElBQPNc0"}},
  { 311, 19, {"8mCX4ZSFUOhosQMZhNaQnhT","-5iQVMCiaXGJQsVk735vQMbIEgMe3lZJdKJ1kngL4PVZ8ZNMDnG6Y1N7DH0btBuBe7ET","23ZbnRs3CrfLsos5S0oF2LvCC1F3vno7djAma7egalFGRAcYsFo3skYj2ZlennYAZp","-56f7CpTt3VeVN29j7k5qd3oPoIWHVdQ6gvtihCQLvlEXGW4Hb5jpGiOMSo1PXWFr","4O0Qp6pQi4Y6J4BAIRWg4C4VcWZD4AtccUckRXKMu2RMfTVAf5hn5VMjMp79lK","-2VVE9lhp9vZmFBFpuamS5CQCR1JJA1MFYG32bJciJITafdWhkH6FCJXfOHu2","L6eRDskIc9gToc2O7U6JCZnFpGHZQXb7OBtnn2dPT1DTb1H3ptIM6OeFQ2","22IhTO0XK2AsMZ2RKUbSL32TqsLdDGcHL6blWZ6OvTqS2vYBYuqS1uaf","-25o4Rq9cforcS8A1pZv0svOJTopuskL4I8MOBgJnXbplsC6PIGUoMA","-Ji4ZcXj6rBfGCD19Mde3O45gtgKj5k90OvrrvX3nQ3s8YVvec9P","Do9Fl8OPYq9aMqK9F2vcLrmpEm9sR7ovMCvU9nOrRQoL0m0Mv","Rpd0S9dilmjbUQK1HnO8RuB9ZrorqnBLrJEVORXCFReCLt","2PPZ3BCcduvgL8lVZ0ZZ2YOrHXavTMGlpceBP3NZ0fok","aPWBiqWfGLcSAM9LqkX9hWmHOQaevjO7NnfoqoNd","4fFNFrevI6MD3WV5eCJ4oteeRg9rpndnmkALi","-TTJi7NTC5UppObjG7KjR4Ac2W0gZtgqM","3DRg1LHZWmWl5bui8AJTFLj4aIs5","-28vGqPjpIWpAkCjpFQ4","De8iNIZY37Zn4K"}},
  { 319, 10, {"-NMX6e07BhuAsT3T","RmI6FIsbjeXYFErIUTEToXVZ9ml4WvqTCuLq94ITULET","-1akH0fuWURS8NAHj4YlQaNf2KFRNvkhSkHNjEP3dvRp","1T1AIP04k8EglWE4mAbOAuOka0ltOpAc9dGcRSSIp","-6I4poh5bBtO8nAqf2JVKKjg786oTIiplsI1Fcm","AYmSd1oCHUSXbbXAaX5B0sbJQK9fO6qrrd4","-2C5ecOPCMFPDGUjWJKKWQBlThPI1apeg","DIYfdfAJ6BJNCQavN6MdFmvZCsHE","-4ZTUKELVuj5p9THfi8","Rk8YNS0d9AsF28"}},
  { 327, 12, {"2ebPM4cTtsWYBi2cET","-hh3iI8rOHaUEXfgd5JJW6dFGSD4TaXIpZOl28rt7nR0lMeKLhT","4VhC2jRCYfKX6nkF4pfeFCTPfgGQWBaRRTgbBSmjG0D0FnQTV","-3d3vni9I7hqVIThBMdO0mqmqEDXtS6LCisBFCkge6FtVhLk","3n98fiIAcfGGgIQSuqUItbmdE91fUSMADHXlAsQ33BRW0","-ItjP89vOnGdHdtpdQ8JrBGI6nMDDn2FSKBEpQgVmvm","DPleOlvBLY32lKj2HgcPcZtUkdGULJ8fgnp1f7N6","2VLktVVuvaBUpgMAOtJf7scFG4PnMQ0t1SHXe","DsH8VhvmW2dKalXJ5lrC6d8F0r7YLUYlo","rsAuKrVuTp1WnX0lBWJNs4KSbEFo","EQBaThZkSUX7bmq0lqT","trtESKJKU9nnCA"}},
  { 335, 18, {"9VdNeHBGGai9oLt4SU6cilB","-14WnHcJnI2GjcNqFPHfgIm1Y4VTHL6VVUJLhkAqNeWqpPSYo975PfUUI7QUUEJsrn2X","-3gqgGRMaUorippXL7bRCRoFfVCjcolNZ7XDArNUdZYPZiJckVofNTYYYuF3hXItIG","DjVK5Ue3KjXDkYZO8OsET5cLMdpqr1eROMYoSm6but4XKSmUPodLorSa5LUHCNU","2UYl6L5Pu5O5XC4AcvCBdllpdoUGdLTeu6YfQ3GitQeIL0InGLGTYsds3U3qOa","-3ekTTonQcUeC4WIcYXga42R3V0mlbCU0uXSfaRU4d7VNoSYRvXnOk9rDE8i3","1EQ6tg45cf4cWrIGYroBjRLM8UunNiBKrbP26kHK3NYsWPBJbpp0nqKain","GA96DOSeaAk9WeEP6gbsaDZhNk08Fc1bb0foVf1em15FfrvWnS3U5Jr","-3SclEllL87IjCvJbRuf1ZLqamh0g9PsegQ56evRfUhhQffpWohkSr","-102ckEUkaYJPSXXk9V8NN3dt2P39pFerjKPB9TbLF78AfWJsjcc","Do8j3uv1LScM3T7a1VKdPieeDjvWkal1pDh7G8UlAFOtqSoH","SXP6QNFVPWTpClqnQGFIOK6Fhv0RT4QX1cJigTqAoWEbv","-4eurI6SXgb4TgveIBJrEpAVZFtNLIrWKJrNSsvnNY0","MOCJ9WGhXMWS4e9O4ZqYmuVsddsnviQRbmSB9I","-1BiLXX5ruAFE03r6nK01FEj5CS03U6EQYv","3eeDRc9PDEZsZPnhr2V9p0n7SCUTf","-aF9V6bH0jvVHte18IIY","1rYBIJ1Lr0Ylu6W"}},
  { 339, 6, {"A4Nb8uFWhg","G8GPsHlSFs5VvH1XmcAYN6iRYpo","FsuYrIklrapj68LdJrit5Ddaa","-n2W3C3jSJ1XAFUms56fs8m","uVo0uiTHXOZSRN2blNc","2f7hD8RFi9apcRo"}},
  { 344, 10, {"t7qrJId4B1cfFTg","-5dEZnhebQD8VEA6e2jeVZHbMXPjCQjsZTUbgNqYq8GFQ","A7sApOsW6kXX7ZsONfUl6nnobGje7935vSdiOV9sK8","2R6Sj0fTJGQkc73bICiugQg53GpLIgXrGHSYUgK0","YuePK5iS7ODg8Bbsp43WRtToZh1KA8g9ISosS","1RVYZZTMCb14cH5402RJTOJgALI5VWUmHmS","1QHmbYcOnrl7ELgIZracRA1BbkNpoDLA","3T9LlCTotjS9eEtpl4VLfm6Vj7U","8PEUEA3Z3uefPYlhTL4InE","-48KNv6H1ZVbL0nk"}},
  { 347, 5, {"iSDiHTqNY","-aYZUeWhCDXkqRAeCJbivYmHhU","158T052DV7XEs6Zd4KlFHoG7g","-2Nn6HjHOqJmVBHhLqo8O","5Jf1SJmK7PWtJNc"}},
  { 356, 12, {"-GD3vkKsZe2YG9aq00","3siPLd5S3OvXWIj6CgNG2RdKr74P7fduf26W5N01WY8shA000","4A6qFIAdMSgRmG42172O9084n3oL6RFleYYjhenbPDT61s","-524RvAWbLob7fJe4XhZgq9VEbbLlYnneeWUpSd5qclkSu","-348n9hVBbgvr9kZ6ZuGLnETSBYKAIjU8Opdo6b6CPXY","aL43j8OVBWYXveruH8nJUf73j5SvCcKmo9DMn2Ro","-7JsHRZlfSKXuVC3O8NNKG6vjOYDFNbaKgDV7vo","6tYcceas3TWN0E2RsAXEjYL7ac6LbYSCtc8","-CpM98YJCaMqulXKITVPrdGI8DDWTVpk0","-MBLRDH332GdaLSBqFXtVaTORZ7c","-cE0sM2Mk9n7q5QH08n9CEi","-BKIZT7rF87cpJDg"}},
  { 359, 19, {"GUvrAqUEM53j3ReKrCroSLhT","-7HSvDm8YODDrmuRmCHSF27sHDlmEZvRRHSvNM6jKCD9TbTsF835P5jgPhnu5TJ0Hb93aET","3ibrbOgS2SPETHaUvEVpYLc933elZHtU3cZFYP7uJhOWJ0lDaDCfFUaPNSKXtnrr6aLLhT","-FBWrtrYAgMmjFhS1HCq6CeqmvtEBIQQdBmkM3eU1d2GSU02lU0hJE8td4S8EUsGcYjKc","a4FHbrmD2atNPMEWq8QEELhaXcISWEgRr7UkdlQlRVDkgqN6TT2cSaEnPf45HrLhgI","-1KPXpp3siODqoUIjKALCaXKn6iGheaMq8tJ90W90FZi7BkAM8tkgTZKhBi4P3YLWh","1H7HVdTXOTMGL0RY6RQ1j7cYij4DhlihvoSnV4bNKNfAPhaVGIcVGglgJCciQIW","Ne3WgFOpsV9g9sXOpYqWruQmTn3PUkk6TuV463E37fEphg2357M39rtBKQ7","-20gHY1TIhJg2BF3IC48KvEiQ4G9uTLMMBjJSlGmgdDrTCPel4jqeWSloH6","-2BZNr0P5gEFIrhg1KojvibmJDsIv4LW9qEmqslCB0NLvt3V8DQh9ouU","3MSRBkT1SFMRdV8d4rQ0sj9Hha7MNqmo83qs5uadQDiP9QODfYu2q","-6fjo3lC7nLOHS6hKfgiK5KEjHnlotY36nGh1a0l1fpaP8SIWsZ","4qEJeWBu2PnX0i64TQPlLX8UdYRIJIoYtVPh3nYpZ0dNvQ8","SaQDUMAiKa3LvZknJg23QZJ8JrnERuAb9biQ3u1a6YO","1PaF9j6P9aOKROnpJZVXbJJfZu5cENq26SGXNAhB","-2GVokGIMql7oJYfXuTagWCdh4ct6fcfXalv","3c4arLTjpWRTUY493c2QEs45qi0Zqj","-95qc4iVG1eIS6ipvA1Fd","EWbl9RqkMCBaZpG"}},
  { 367, 9, {"g2ptfSc1sGXPhhT","202judZCopXDdSW2rIjbWU3uFXp98utegE9CScENGhT0","49AHCN13l5CGbFb8BIL86fQ3b9bIEGLk6e8FUmAEs9","2i6IsY1F9uFG3SWBpkf9bosXCqUQuoWArSV1U3dP","1EM279pd4fPl0JlEo9SeJZEfjsjOWshDhXiTBl","-5mr1Jksr6neQ15XQNl3DUsM7eLhFYPYcNl","DdDBPEjrZM3Uc8WGMqpSXNrAlZRKWq","-C292MSg8U9HIhhIJV8r","S9Ue0ZuOHHOM9j4"}},
  { 371, 8, {"1GZteMcNh1SJA","-AOt3mVZG4dklPGIuEIgNSjvQDCXIO4HTXEO","TXRatnOVCfbgEMHDrA7aJ97hGQMiSMGtQ","PpjkIackd3JNtsQLcQ0ihrBvavrSqCe","HT8FL96HSvhVWif6HDYfWNO5Fh6m","U0MEdU4CpCV1fI6ZhnQ5uYdbo","-YDTBc0Wm2Gt7TivoqkMi","d3Y06RRovEF4EkC"}},
  { 376, 8, {"HHrO2FHeI6Wu","-p5ZUTH2ZevYnZgMcVglI596dU15v7TAs40","NaHqd7mZeisbUHObtChaIPSeTp4iP3Gpg","2EUV0Jbq4YsNh0kRScTSi0TXnmVv2R2","7j21S93Bplh0M6aHdI4XvV12qZRbQ","-4l7eeribKKEISXtUIkCDdJqGcq","7h7Qh0hZAMHiHMeuMjpnMbE","-10bV7LCP11EB8sUi"}},
  { 383, 17, {"1cdZCJ8U83AqfSaHEPCQ7d2t","3n3ofKv2Z03mPcSCZ8c3b7mjbXqZV7jnc2ivLcOHRF1SWLpBN1Xr39NFvc3CpZc8MTqU","38R5EMomR1UlOMqSFs9eU6BOskeJQhI4c3UqHZb2DBk6HMUeteIoYB2MkpOt7TjZNhn","-8bfGoNVEf9vp7MCvEdefToYQuf7pDrcXgqr9rBdV08oGF1vJZ8XW19dvkmgNN35A4","7tfFDU2Nhfc9BD5DOE6V2C1PbrT5ljkPYhPqVNqttl2XONvFQJ3JPPaNObKKKWc","-3IfM2clfTWpKN6ZXHgpcBmb2uT9lhSLpWeYVAuaJtkfAYF1XbBAHeDumAMZ8o","13XsEYhmCDodZgKQldFV995NE1OZIvM4pU2nKFYpUhaJSfTPJJcQN05jMWS","33L8FIdcoCbvJuq6RBn2ecoUZ8TFJjuVoeN1ZXc7p31lDnZdDTPFlLqF","-nTYYRRWaqFuFNTurR0NhW1jajQhtXGMlYNfIjvEee33u2NKnrACDB","-32cF0R5rnXvLqa45IQd5AZLjXOd5SRZY0eMv9GpklQsI4INREWu","lI7hCq7ZcldlPJZACeVCsrDDMNW9uTMLqoHfpZ0ZcGufBKl6","-2iFulrSipnFpNG5KRRhqVZolUqHa5QSsFCKCoiZnJfI15","4dNvE1gAr0BXnpaHvjYd0mvY4MMbv90Hqk50NUTTU","-3o1V7DH1pZRYtTefRLTjQgR0jdfuVY0fSm6g","39LV6Nd0eUVWakMKuSXfXsBHtjbkLev","-29vaR5aBDMr4d4FfiiOA1","1j7R46YJIP5Wnm2X"}},
  { 391, 14, {"-13LgDcHiZ7YBHVCFJFhT","RdI62Vl3fFEqDKdPrfPeJGXRdvttkNsiDfHrJMVmCWgEEPefWnli1hT0","-4pTI0lPVkS3c2vDshUSK4EqSBcqVcJ0maipOMrTLYgnidY9vp6Sonjj","CtBaEX6dmC5qjC5UGDCKOS2YYeADWXkWA7B4bfe0mdRniROQKVijr","-H0vBA5WFrWcKMbDclXiLWv6k1ObI0VGXRirZ5VnaaGgZmDboU2W","Wesvb3J76d59O03AKOmR5veiTepsm5dnJ5JV6fI2OJlE3OLLm","-KKBe3au0Wf5gucacBQtWbDLOkrlgRWZbAnTfO6Yrhh6b5uS","4EXTqq44Eu6ONp3dvIQNZs18Ekt7X02THEn5OD8BhR6Fo","27DjcFgVJILPihSBcWfTkE0hF6l9ce9BtGl2Mqh01Y","2mJH8tXvRAnbqArd3vCs7dKcNt7sVR98hUWbf9U","80ZF6WDgFllP1hqGJgWd55hmeo3JPYRUH85","BJru4Z6g5b4kh7ilfpAY3DpcQacEbXS","2LrNZZG7s8flBRh0Fqpf","3LK38Mj4VrGqpE6f"}},
  { 395, 8, {"2rNX5DXnOB38i","-1EQE6Qqm3bkLqGMVIkjNiJ9e2JNRH9gdAeoC","YY0635HKH8h07fA7CHt0jlKdC8UBl9vjI","LNRqTnmqUCdr210GNNbtcBvvd9eZovs","HbhrOLGEcEB6PMNbEpOEEFgEhVtqq","CmC5tGSCiXeL8qarrYeCFJPlLg","-7eeoDkoGFqTVcMjLjleCW","4aBKIoWlrkkgvPj2"}},
  { 399, 16, {"-IU2r1jl7TXYMfDv3RnJHGJ","-17jXZReNJKnDAPBJSIRNC18o9unHv6bUkHD8OJuaBA8htuUWeeO8vanFOq3RjunPc","-YRPMYlH96UGKrga6mYPOOPPNjRShDGaMMcBC7XT9fJFqU5a6WqrYA655Kvrh69O","-jCaafqAsLSk6j3btJrdFXuvGR3BKiJ9B4nV0kMJdVrQmaTGTv36YgTErNtH75","1b5UTnMRE7VTMK9UJ3o0a9QSa03des7XeBJU2XFSFS59e9ERYOoUWunp4SOv","4IaQYsftK1ufi6dpl4LnY7muipYhBG59pOAi8X1UTdZHKRCkREpR9vep6g","1WPe0N3PuW2MWuNghv9uTjA8lbm0cb1d0Hrtd39AuK1Fdiee7lvgF4St","5YDMPvu2KXe9nPYJIQh9hDFTVN0ohm3soTG1iUV5rHL2QNDbTOTs","oEi7W7IhEc7OnifGQlqfT4v2dWiE1praoBKNdJG7mGQ6MRtceq","-10rM6ANMr9ZZIJ9bKvB0sD5t0dMN8Y4JIOqmD0JLLTLVp47u","R1ViiN942qhDZaBVKGOd7SEIOpf6gNt5PQiCjKDeLtO8","-15A6d3lBP1To5cfkWVfoBjuLmiA6IIMVhQUIS1qIl","1FhvrCkCoSZ9DXOOqiV5ec8Q6KefqcnR8qIfh","eA7hZZ4bYYZ3Qpg8WCt9RUD8AjoTD9U","BasMMUHdlOuvogtJutBKF","6JaCfnuk14ppVVKH"}},
  { 404, 14, {"-9m37ksYmj8avDMD5g00","jD4ic0ChH6DIQHDP5fn1cnhjKMDqksHfPo2XNET745jWBdLIIEC000","9R6qiK2AU9J5q4iTGrMKrrjEXmsIZX6I9C3nbFiZqtgdoKJcqSNXg","-DSqAmUb1m5VFg6l9KA5HJs641qC8RFX5ACcEv3Ulm0AUMWjrlPY","-8lYV5nqBJ8uECVU3dcQ36rusjFII7u3ZAvWA10sIEoTUDJsWGW","-M5ANLgTnP1Qo921vS1YNY3vB1C196EtY44fuVcQTp6rjH1A","-1SqBJ7TC6PWTlQkDtqQ2Z3ik0hFEd0SJbK8koIO3TSRN8S","650lrpPksUl4B0KltrI6PL3dHWmmEqoR2Pf8mLW4O3o","-XaLAqXVAlp2DKhCTclFMhAWrfHKchrWNc9LjQLYO","-dD0Nvtug2IJlb2NLMuT9QAZZVPXWPf5uonvSS","-GjOXHKbMLoH6AjhU6KQUmPLGHvfqULN8GO","-7tXTUqlcgq18ZYMJKX29fKGfKatlU","-3ik0FbfPj7oApd8OkQnD6kXk","-9M2TMrZA1gs6bkDk"}},
  { 407, 16, {"9NdpgjVublBhWnEt0IUTW9V","2Wbjano5GPJOghXU9paRc3qY6ifiA8Aao9Jv1ndjDQnJK8eEYaJ3duqd7VAbcJoE98A1","-AbNfPa7gZOWlRkYPiunjWksMgPA7q4IsP7PUdqf4oqMcJGWm02ttWKaGG3k0v5fA2v","bH63M3AKTNlgCveDlH10nRjLHc4fJCuAbbON0k462hrR9hdSNBS6eUPalCUisEXM","-eaK5RCU3cDDu4PM777So1LJF3QheiEN9HEi6bUXvm43L3grJPQThtYeYVHC21n","T8PjKaO8HksS5EY7V67C4voXYfGuLqUbYmvkv9F3F8MnMSXPu4i0R8634SIJ","OVCWYo1YNCJ4BCkB3FYhNXSnLCLS3rmrHbkYqhUgPvdmn5DhpCfEdfYNg","-CeRl8Xl4tv9NCrMnZumVhOsABpsUB2vJiMBSrPZdD3Ma7rkUaU8CFhr","-DKWRmUPmR2FuKnhOO65DAmG0tIslCfELuM4SnaODu43GYslNmc1d","5ddBdacEoRbT8nEVOuOrBWH3kB4JUp5gM0fUojFpKN03kcfrHm","BrPO9eUD1lUG1SHT77Fr9kQJTNblMqEbYHHoLNXDTP7kBl","CTAV49CLQPCuA31oIOkhUCv6Sre7T0ZRNsrNBT2WNo","-5S2cdPTj0avTZSgkR59C4PMkeVQEA6C38aTMZ","2OSuTOOdhWHf1A7KobP1KKRkJBRITg9k","-QmEl4P0mfLGMIgrXhIY7V","BnSJS2jH528D3HoY"}},
  { 411, 6, {"9ehOK3r7aBg","-3MpmHDc5mXRsjXY3ZoRQK6OeSbc1rQ","27BV3USHXXHOrEccKXioMLNeoqKu","gc3EnFv9O1EK371n6AXPFU3I","eqHCpdZ4IRvmuPsgbS7rI","G9upnrdoiU5ZcPC4"}},
  { 415, 10, {"CNUPgritrDJdgOqlR","-5hnKToTLRQ5JsbSXfWqgcR1OAoNc0AEZs9PnjVaSeNVUdCFbX","KMvKrcdn547dZScDrKHsI2afiTSAISVTTbenOqjQgY9QqCf","-NsIi3pPo52ufMNXkhfLlIbipE3N9OVZCWBM8koBRNnFZ6","9JqSb44tb4K0tDL46JVdVKO7f8EMMPH3ebnHsRC2R0S","4oL3sQ6UgaJX9kgZ3Lc7KSNFDpQAlrft9cXeHKoN","-91l8RnXIHWHW2rkdWAITnNLjnA7D9eoeJ1t7","8LLnqLShDrd6HcW0IOoLselLpBdQNu7q","-Nj9kOD2DiJZLTuMfdUuv","M1lBKUj3UJSAcSVa"}},
  { 419, 9, {"8VtpElNPR10gSG","Vu4FuQNkjATVUcvi46GCWH2AkDC9MdU3SCrZvMu","2oCvKfAoP6vhRZrCU4akLFJFKLdOZpOAPfr9YG","-1p7kDqpfb1stPuDmn66WCSajTjRZJg2WUfsG","L0L3ZTFGM0T6CVMEnYLOv7Ek6EXbrYkBA","-GP7qG7sN9b4UqsWoOLCSQSugZOQtJ6","4vDXuRqlAGGfWDqW3EjlWMvi0AC","-1Z6iNGWmjUoD94islSiugS","TuAZvNACkiHMo240"}},
  { 420, 8, {"-2EWpcYYqqe4llU","IlPPuiFNF9VsIRTLN3JOjS1Q3GYXqqdoVZNLdE","1sfuG2lF0A4X7O6c4UKolY7trTXKTTk3fFJRLg","-SC4U3KOf3MdGphLChZN6v9OWeNEBd4DtDma","-2P5h4AZQj3oR8jAZB0VLfjqsC5d3S35dU","1J1ZR96Y3Nvqt5r50GLMcb9bLF7oVA","-O9gi788KVsFrUBiJ3vD7hbPQ","-WKgXNpsBIBDJAM5c"}},
  { 424, 6, {"2VNusf7T7440","-1i0rdiNiWO0jPFDVs1SIZDvI4DlmIFaEC","2tjhdlbDa8A84qnisTprgRhIUkZUJaC","2F9EsGNLYOgtNNgU2Gj2XCYXacVo","cDljpEma3mUC965pBGLNiKnc","-hsCe1VqGLsXoJNou"}},
  { 431, 21, {"54043IoSbW7RXqg6BagrPAuXWUqF","-SlGSJE7SDOBEePF8kUaDTcdd1X6RFDXMpG4FJiUZp0oX4imO0VesprVao1cacMg64jF8lKggab26RYnLu","4fBqC2McTCoJ6QsuLu8tIgRQbJhXmDBbp56ABZ98fFV4beHRfStFF5Opn9ceJgUoEa26Tm0RZPYE0FjH","-IEQclsbQM72vIUtFPfL3XgDWjbRsSiIrlfDELEus7QZNgdjbu55QhDZYW6qrBuQpAb0XiZmt5PNZ9s","Jpq4pWUZgGQgZAEXCcPob1TRL0r3frcb93akDcTOpq73OXXcnXs6XJMVr5N8s0HimOG97llPoTcW","Q9cE5FD6d5NAK7IaJpNgTDtLuSnRK5Df30QL3VEpGtYQ6qriUekRnr98Wu6u9RJIKACdeM4Zsh","-fOZmbXOYN7sHiQK6j1n4eoMhtDpKahL8rEHXPMRHNhbSPonJMZEEV1OWUN2EC3rBuFS8Dd7k","-dXetP0FTh6spXfMlGUdd9BPZSZRgMo1nLvqE4l3sS9bMb6huAnoNtDJJgLVlApFMhH7HsX","vQmukC7lqGIeepmnCE83Ti542esZARDKqKAbMPrjMnK1UV9WJLPjWR6EAXBeNAqvkbDW","SrmX6JvU5LosYtpgs0gQf5DHShB4fdFsXEs2Z6mg51dG99EeQ5qjoFHDlUkoiqnpI","-DO54O5YcfkRD0GH8UIYAShhsBcILumluTNDVOYV2Qgf71UPds5NhOvQNgR0mSqT","-8FSjRQBpd6Cu3UoHVmie948eu8Rlf61ALB92cZm5OWo4Ln96fM6M2o91qYVP","3APLgVUltUDZukfPuX0BeNjmZg74K6sNjgrHUNJqtat7Oa2kYL0a4pIB2Y","hdS8nsrFLATc88Zk1vctWnevt6TP8hlIqKLP1XiSfsqs93ZsHHDKjM","WQ8Tl7JmUZZ912uUFGiKJXdM2gcskcBUb4WSZR6aJXZabJAdpUH","-hoLdJGMHRsUP97KYbjWoUmk3ODasv7uYjgmCnXWFEJqHNb0","TdX6jK0LsWNmDj2uL4gg8AAuGJ5XmT6Oekc5tk1Y6Sf","-72aM6Nu8asTNgI7KmqiSo8XmYYhTbW6NC3518C","1cNaaR4aAHRRBj56obP6u2PsuAd9ImRHk","-5Qsqmb7VhjjW8Hg5H3eZdU","1Gj4VTWWXgpCh5cbk"}},
  { 436, 6, {"-5trOSonvfeeu","FXa46AVSOmCkGmWAJuM8uQmvb7T2EvMrU","-GrGereJBLNa3XfYKka2g6UeVCiM1mpg","47dYis9EsnhEQsGZihbG2KOfYuaO","-2XTgivr0fTn5T6DGlhCJfPGqW","-1p2Q4qsLBkOhYqgbM"}},
  { 439, 15, {"HDKqjg1bf34jFXbIoAiQB","O7fftICNPUm4efYek6M0iU8t6JLNa6sUlRHpSihP6O34B7gD40Yta4mDNscZ7","2rjUmkiK5QMINPKVOYMbCrbGqIk9H3EOPeIeE04bsWXDmJCfRbIBTPsu10os","-M44cYKfPMWMSnEqCvI64XR5hq2tdNLLKv21BDUbjMDPgqBEqHTBgcbjLmg","hRnLMRJGnA8m2Ctfhcfot56GTioWuUbgsctjZH2QmcHcrAqG5NMgRpL3","-dSRspMH4mhnYlJDrH9HTW3qhqVrr5bLHPl95WCshtW5IAp81RNhhKZ","MDIFJ1fYuIZmUKGLt5KfP5EiTuSt1gdntl3X345WT8sfp1HULS9s","-4fjIZ2iTXp7f203Q6furpKOQrJI85Bpvv9JQESMnDpQY7iIiuS","15cBsPjkZhVSFlEUJ3v3co2RZFNeIhJEqFt7clnfkMWIEHnP","-aEtZslqJTOivOq5JaRkUIVrgdTeCGRiQWODkWbVZuntl","7jVknEHC1nmpfFIFvhAnhqPVXG3gFdfKoDjXrUYjv","9H6K08fUOfMpSaba4aa8jIQTEeIcuY98tP8Zp","5VeY8520NeXmXgtYcSfjUB8mkhY9FHooc","3spi9GGHveSunNLoClEl9","2KYYPQA4mP8QYTVeX"}},
  { 440, 12, {"ZNg9OYE7nHeNSUZmOu","16aNFMLEqZGHSmisOpXoLTJLCJe6EUsZFS5fcrSfUqJ6lZT3aL31o","rA31YGtHjILCdRNGG7itWVbR7faq8HspL5D7EH1RbHMCev2GUK","-5KI8rJrkuaXLkPPOC5rLNpL34Rj3RtZlHQ8nXZKcr0RsGtqFo","-38o4P2ppQ1E1ZsMrpdgctC9b5W1Q4OjKdW187Z5O1LLeskC","-anKcQ3ZIPIQXGSFZGWD377hnLSEoVIX8AAFJeHBkggXE","-1oDeRR2JWdAdWq08GOQ8Z3rNP74EOZMcmsZfu1q0iO","-4TRihUcDo2h1PKTWFEZ9NbMRMm0VvBHsTmW5GGO","1sh6i8Bh67rneLMXuvOB1Z80uoXcGMW793Ie","LKLcEup7DUqIViYMJrAhUO2sB4Lu9k","42WHvqhnGY47CshU6GaZttpdM","-2VD12l9paTSZlAjhE"}},
  { 443, 5, {"2elXQfT9hia","-R9bgRDurvYs1oKQnAg8RRZuQT5er6","1inIlqrsipINVQH00gSrP98Ae5Ka","-I4KQDsv2SRTutrYQS9CWHs","3ADmDbAnaVUqhNWu8"}},
  { 447, 14, {"mlEEGR1q25KICW4kKlaET","-QPIcjEYJgrtshPa8GpQKGi6Nh3AEISvNvGMGVBlBeZFZ8sPXcQBmiW6DtGI7ET","G1ACNZ3U0mMIlVW930Qdpn4DbiFAO5feI1eBGPluphRmao8Vs65Sm4CrCeET0","-RZbkCaUiYaO4s3ZSpbC99kltHKesTtFlneIH3WvHvSQi4nCktYSVtteqNcX","HdHmNQ9LZsMvei15QorKHNgBbS4XRPPQ0UhQbfUuRNaPDX5t5UEvBbfea","-4bunEmm1bLOflagX828L3GaIiOBGYZtVe0vYEDLp6pJAkXbgovZtN4k","Sv4LK83fFhueDAWqaTOJHH7C8ZZR2Srjdq3i7Q6QQrL7npIYqamc","-iIc1hW9ud6GUG2AoHD7XKT6EtelXf7oIQg34MHoUvDL33YlQV","46g8T45vuJ4vrFlKpIrblUHLoG7HgZnnKNHligr1aG9IFeD","-41fdNTf5X3YrFTiUNbb1YdN7pDg0JFbrhla6c8VgGb2","1uNCIIo3hLe0blqlLFAUrrYWDYuY8Y2pmvTLdT8","IHkWnsjQDr3YVQvDAXmJHrBoTWnE8WY8o","QgiUB3qefIPgn716KEXFm9","4G9OEAIb4SGajc6N8"}},
  { 451, 6, {"8R4slKsQ7o","FtOXSD5JR7gmOtsFnFqMYge3TJia","4kGX7RqB1eZLq1Qkqj0ttvU9S4","11fMXYY0HHN3hLCslJo0TadI","BjrkfVV1AtvAOJmebJEW","5hlaqlImkJmPAvsTk"}},
  { 452, 8, {"-kAnWc5ueVZGXUdg","13b7F6JOA009rKemMqro4nRv6MNEfHbFaMltrApCMkbM","-1ffdT7KZNfTir6PlFPJXvqI0juf73Loj1iShW0m1km","4tMMdNbhePdaW9Lb6mld6n8oHIVA9784s46Su64","-DPopU2JppDmG55LXJCh2fctDah70KaGm7cFg","-1pdqj7opVFEpQfkCh9BtiBXG8AZK0Ai","-FOjmToBlRa6FTvQL9SeqEk7k8","-6BPAGMbZWSa3Oba7o"}},
  { 455, 20, {"-Kq0CmOUDbT7LW0DuZLbKsfOt87H","XY0AifeLmd1OCBApZipW7P06Q1llVRSGFR1QG9VIAmgGJ6gHhJnu4KKn7RgukeATg5Zh6TT830NAlbLo","-1Ea1EFDhqJEbkOF7LarFfb266P1iegcgjp4llmj5IQSh6NVHKJ1AHC2Y2HKYo7BXODZmsulakMc6lT2","1p9XHiZ5a27IVa12In8O8hI1QLk0vbaZ0B0uksgZk91P4vOZ3UeaTr3hX0B7DSn3VPRpmIlDKrCGm","-1di3lk5UsEYmPQ9F28viBbonGnKES9v2pnF7nlg10nIEVd2dKR6YHuMe5PNgMRE4hnS6PZ02lcm","11PoJ6LYgNDQDHLqYOQShgG5PAcBXPFqCcumpf39vt1ebiFXI58Nhpp4Ne9PhSpuD9iBfjEUHT","a7krgWbGgZRXEDp3aGFCXDoZGSQJEmR2dUVSDXPJi41ZTcW1eZ5bajlJtRhTPatReJqcZpg","jPNbUZhbiDf1W6TJQA0APe3un8Km8YKc27rVtgQuoQ3jF0W3JAO0bJfrNH3teCuH6WqU8","A6TJCgahEUdmDofWgoB9pvg6gAKsPPWCjUnnQ0ANRRBFXpRB770Z83QhGjA1j0ucSRE","-2GBdmjWicGpnWSHttGqhvsSiijH5EHFIq9HRHYDd7944LsqIA5VblfJKBEoVu0Ce","-1IfH7oEXEndICfVaZVePIZNN1hZAYPfOh39YQesPvsfEVumvXO3D0p8jqHrBMe","CrcIVnWkTjhlOnSeDC41VYLEKfqYqGBCsSmPIjXdFEjoN0Bl4btHrugGJM","9vdUtSNteeEqSc6QPoP5aN8uD6Ml51HfLojt4ShAV2XmbR2pFG0raJYc","32LTfquLjaFbFOFpSGf06FKbpfpU73jRXf9Gbhe0HvX1rs72bYhSZ","2TZ8ABkv8RgHroFfFQbtFTYvL9o4Fo9S4uipLDW4XjhmrAfp0","123dbWoISBGLmtUsI89QqhDcvVY8cchuWstvVTJD0rFDn","-7n6poaOaj27edolZ4sIfLN8dLXv94duUdpncUKE","11jT79B9IXmNe7FTBfXIeg5kW3nfhpP44M","-10nthj0hmeTM09m5mlYFvcF","7gQYhvH9pAFtLTT1R"}},
  { 456, 8, {"2C6qpRuHg1Ogq00","oEAoIsojuITrIaXccnFmjqTrFYkIrqN9iEC000000","EErmCK7hX7gV06ncepYtQG8LNdKkbSsWvj0L0uiO","-f33dS7jf1ki7N84P6sF5MRtMHPHV7uZ1L0K6a","87qNCioBYnOGlp3enPmE89m68aMVEbNRuS","-3EKs9banircgmlHNhuq3D48sDhZQYhQ","O0XpSsXVaaZWoqlgcUlCL7gJY","-8IfANXBKORY76YYa0"}},
  { 463, 7, {"2ZM2ApcCJ5QJsJnv","QdFSMo5Z1bZRlATaE8PIWZZUhZu7aoeJJ2goPNRaMZKM","CO66g7slLrln7VfC0kNQQa2c4XrTmjhiN9jiZTcPsR","-95VGGeTtHo6OVTQTDaUBXtOL80UNnlASbeinYI","3JD1qn8B6KUdjmShRl0vmWfV7TVpbDeECM","-btFMjujI5L0sfK41YBWaF","DqDM0bmO8LeLnW7er"}},
  { 467, 7, {"24EmRu7iAbc00","-Mdhs4mXikI7SvurMFbgvkTDcCu4SO000000","GcavXnEEYV2m9a3DgraUsD3pa8vTOmRlk","-2OGlKuUrceDFtTO7QrPO2946qPPAsXk","XgsO3IoVhFj3RnavWsNKEgbvkNjo","-3GkKAIXHhB2UlgFpQa7TdiK","IYtFY7KrUfQZqJJ12"}},
  { 471, 16, {"CQr0i7BpYN9Va8NbKPWJBEBp","-3XNdMokp6v4HLpVCOLDsm5hbf6VuWPQhIWk31HTgCWk815RVtfBuCl8vKfpLJMGINOLRn5","QunrXLHd4L55PbLo0suLvlM4rGsUCoJHsuY3j615aA3VYWgJEuRo2DhbIvZUFl6daJtf","-1eCgnMvNOp5B5FvSH2E11Ka3BMS7OAcGFAHTiBdn8s6Bnb25gbjONAeH4ZglnjAckVr","2cB3s1IBpPvJNarvMIDpQiqRQUToDvR0LiOvtT3VcdjHitUSkjmVmt97hkh8Gq65e","-1XAjMgFlT41MkAIkSNUIsBTm9OFmGcVdfDDETXW7iM8cWqhCWd5KH3O0hb1L5pS","b2GFaXLnlAB4Lb3sG0h4Q3gZqUhg5h28lWYXnXLkDt72bElr6akYrJPY2PDf","-3n6PJFMgT00Lrnc8iYK6WG5WiB8Gqlga2Js3WkLMCXaji3UT3mjiTIASl","3GRZkjPWeNOg1XkKbKdSJkfisGDUol0nIAjFKBoRfKbSrZrLpHvc8kj","SpphVXGZJCZqgvlTTmfsLU0aL75gjBBYH44JODJjgdPQiigKd8P","Ac9IOfmjl6vjrfoZKLSDhq0qGBkImb4XYkFHGS7f6RvrKBA4","6cZWbjIoKA3o9LntsO19lOC3DLusm4hMHjMoZaJ2g3AI","250tWTcYFKNruqq5XtGCJsgdrbmKsE2mjrQTIhEU","AcDbIteq4OOiPEtKQCIXDqjL0Cff06nqNA","4oGlth624OOa29ImSXl1CT2","Oo7aPUkRoYUeYknFU"}},
  { 472, 6, {"JdEcm7kWqJY","1dJhhpenjqKrvSUsOnE6An12ljYAVrI0","172Thn66YVoNNalYM3tWM1Tq0MFuOnM","-G9F3iI9Wu6VohYGL0o4dhUOMWv76","2MBGBpWaQcvjb5mM08njGvW2Aa","-QgIbXuea7kp9qOvZk"}},
  { 479, 25, {"4p16uBP5KYZvG3LaEKpPIDjWm5TEsLhT","-FWIWUn1fYVtBiphNhCt7h9qDWA2X1fegLH7uoiDWrRAOETt76ISr9mHJvtQs9kVmY2GNcG90TM2enTlcC3Zbb93aET000","-13ccZ4veMrM6SRFVYNrJOl6TvDWPQDnOTXQ5doZObe8Cdo7KFZRVGnaAeqGLde8JGOMZaLqGjG8RnKjQUOS4aAAohT00","Ae6Z5blos3MUH1JD2lJimpq95pqjWcrhHmKQGESUm8R95na3RV7PTuRP63WPZGrUHetaqLmj3S3hU00tuHM51PSaET","2Ier04SnqHUS6iTRoYDiD0dL2hunME7CrChScPc8ZXvtIfPped38oa999DZpP8IYMnfDnSYMMM0QrZNJa5CjDT04o","-DGm0MWKoRgCJ49TceL4kQctajC4MjrNrhWVf5b13ZiTi8UYFgP9MGPsknCBG1fHa6RSfkOLZrAJaIQCWQnFrhn5","WSWcISoUjbH6madGldt1V7flTCbGHA9OlXophpedi991T9KDEYYUbq2XMYbZZkLnqK8laQPVmoU26Drevc0no","-sYJu0n2Im83UjVS2f4TN2ipvE1c65NEa0S2mKl3NgssWjWh5BUKXsVYpX8vE0q3qXuFhbagLg6PYjD9H1Cq","umc0AR2e9HT7otAn6MeFdgaMtOoWI4Nf0hDRs8hum7pb7YIObjum5SLESgpCJ918Ioh3h4qi4ijiSFYQY","-OcvllK1hTjaLR7Ufc2IvFCItfKlg7PO2RY6igjXsHF5ibT9Hh0VrH0DZA2MCZoUPNQcledLIrHtbhDs","1vIlfjavqob0vRffm2g6l97diEGEJ3So33aXTXptXT0LD4H0LZEgPF3Xp7hsepHVgHRU97faAE1mk","RoXXJqt0DNGAj50W0EImSJf9rOYQLQmsXBXmaanTnmhhcqsdNF6Qsql5CoJN40CgGVlf4nthrI","FnhpFm9LlhVLDTLllQ1UccADGIEoD8kO4sE6FDvXGWsjWK18IIVdPWiOt2uiB3um5KOvT1P0","-7585VWGQo2kUle90JKnjm8EdOnCNn7iia91tsC5OMBsMRl3L8uAWF7qEdo8dl8ClgpCKI","-18Zi9L0ZJlcFNAuqJTZv3NlO6psrWfPAOPrhr2MjKVcjN846GnstsCeXfU4gnv5JHOW","REYV8J2aXREPh3CGrOuJoeBjTbM9bYf8uJ6fkGeKWrmrEeE6XRVRb6tCLX5LS4s","4rDGAg3DY0lWA6UXSciQe3ge9ETbIjGQ20LZAbRQ9p0KGQjAKSdPlZAnmmI8B","-1g1InXReuKWNGVaWFtQALgfJmL2EcEmP8jEdD1ai6tVDt1qrHeeYBCIaXg","DQPu5bB6O8gOqMlTYqSeBDrQSPskVc6k7baeJJakRSgODHRKsJ1otO","-77kHq0kkQScSknXPaVM6SBDH7O8DeXQHGNgafe2ZeprA75m46b","1pdQNtCCUe56mkgDCQ9qLjoD6uOOXKYcoXr6VVR0lAa7AB","-80YnWoPiQkKkrfbH3fW88emRsa95G6DnalH9pYIS","XifvNIupF0uh2tjRFkY6I0Tc7EPQ67LZCj","-ASoCIpElkSj51k9XB0dsGCR","iEs89hoU6TcINYsEI"}},
  { 487, 7, {"DbJtP2dOicPmE2Lr","EJZXanKQF97d8gA1Ibk7UdZV0k83kcU2uTT1STWZXH3gX","GAOLucmF9KU52Zo0FAe3KAhVFojHhGEb1blUNnX21ES","7duIrHP0skGrFNeiDG18IKeCK33kqooJ3GmoF20","1luSmKb1p75UqHQshDLWNsJAA2UJlA0VdmQ","5drHWsPTEVD2eLRm9ho3iS","1KNAEXB0hsA8JgWXQu"}},
  { 488, 10, {"1lIu0R6UUHlH0nGUSW","-7OgSNJVFUudWY6ATH3ZebqLYr1P1Vs5TY9oIUfOppBBeIRDA0W","BBADCKnT53cvZ3tcVtjXnDSCFk3lncNt80mXSMS24tA6Yonk","oKSnKq7VuriE0hDZeCqMfX27ruqSk0FW38mIvIm7RF1su","KlCbZGob0cRXNeB8HUB4TX1cPO8PVYvLv36RbvtNGAK","b0DkF8cUDWtrdgdv6fIAYSbHfdVdb5Io8u8p1LTY","184UUjajia9cnfWXqR0NFpG7jO1HY0dMcjRQFg","3pqhsqWvSDnmHsXfA4bs9hkTRvchNkPM","DIK5koAkcZ8itlQbpUJcbWch8K","-1QAQkQXFJVqaq9j7Ns"}},
  { 491, 9, {"a3614GQ3HG2q00","CJAZNiNdvg0GSihMvAmKHpvG7M1Cjf9luS000000","1ZRBAClZIVZCW6oCqKFkSPWkB2EROliUslIcejQe","3c434UvBOKgkRBmtPs012QUS8JRuXEe7aghU9c","5iTV7AsHRRugS6U47jnV8q2N3AcN5TqVFWvE","aDuTONCrHj4NuG6nMoonA0ed4hHApsJY","AFdvUQbpImNvlk8UQ6arOBkAPaHVc","-Wf06eP6ICOFoSkK7evG9VZk","1k9PYKGMEQ0DiVcHDA"}},
#endif
};
#define NUM_HILBERT3_POLYS (sizeof(_hilbert_data3)/sizeof(_hilbert_data3[0]))

UV poly_hilbert(IV D, mpz_t**T)
{
  UV i, j;
  if (T == 0) return 0;
  *T = 0;
  for (i = 0; i < NUM_HILBERT1_POLYS; i++) {
    if ((IV)_hilbert_data1[i].D == -D) {
      UV degree = _hilbert_data1[i].degree;
      New(0, *T, degree+1, mpz_t);
      for (j = 0; j < degree; j++) {
        mpz_init_set_str( (*T)[j], _hilbert_data1[i].coef[j], 58 );
        if (j == 0) mpz_pow_ui( (*T)[j], (*T)[j], 3);
        gmp_printf("j %d T %Zd\n", j, (*T)[j]);
      }
      mpz_init_set_ui( (*T)[degree], 1 );
      return degree;
    }
  }
  for (i = 0; i < NUM_HILBERT2_POLYS; i++) {
    if ((IV)_hilbert_data2[i].D == -D) {
      UV degree = _hilbert_data2[i].degree;
      New(0, *T, degree+1, mpz_t);
      for (j = 0; j < degree; j++) {
        mpz_init_set_str( (*T)[j], _hilbert_data2[i].coef[j], 58 );
        if (j == 0) mpz_pow_ui( (*T)[j], (*T)[j], 3);
        gmp_printf("j %d T %Zd\n", j, (*T)[j]);
      }
      mpz_init_set_ui( (*T)[degree], 1 );
      return degree;
    }
  }
  for (i = 0; i < NUM_HILBERT3_POLYS; i++) {
    if ((IV)_hilbert_data3[i].D == -D) {
      UV degree = _hilbert_data3[i].degree;
      New(0, *T, degree+1, mpz_t);
      for (j = 0; j < degree; j++) {
        mpz_init_set_str( (*T)[j], _hilbert_data3[i].coef[j], 58 );
        if (j == 0) mpz_pow_ui( (*T)[j], (*T)[j], 3);
        gmp_printf("j %d T %Zd\n", j, (*T)[j]);
      }
      mpz_init_set_ui( (*T)[degree], 1 );
      return degree;
    }
  }
  return 0;
}
