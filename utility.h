#ifndef MPU_UTILITY_H
#define MPU_UTILITY_H

#include <math.h>
#include <gmp.h>
#include "ptypes.h"

extern int get_verbose_level(void);
extern void set_verbose_level(int level);

extern gmp_randstate_t* get_randstate(void);
extern void init_randstate(unsigned long seed);
extern void clear_randstate(void);
extern void mpz_isaac_urandomb(mpz_t rop, int nbits);
extern void mpz_isaac_urandomm(mpz_t rop, mpz_t n);
extern UV irand64(int nbits);
extern NV drand64(void);

extern void mpz_set_uv(mpz_t n, UV v);
extern void mpz_set_iv(mpz_t n, IV v);
extern UV   mpz_get_uv(mpz_t n);
/* TODO: cmp_uv, cmp_iv, mul_iv, mul_uv, etc. */

extern UV   is_power(mpz_t n, UV a);
extern UV   prime_power(mpz_t prime, mpz_t n);
extern int  is_primitive_root(mpz_t a, mpz_t b, int nprime);

/* tdiv_r is faster, but we'd need to guarantee the input is positive */
#define mpz_mulmod(r, a, b, n, t)  \
  do { mpz_mul(t, a, b); mpz_mod(r, t, n); } while (0)

#undef mpz_divmod
extern int mpz_divmod(mpz_t r, mpz_t a, mpz_t b, mpz_t n, mpz_t t);

extern unsigned long modinverse(unsigned long a, unsigned long p);

#if __GNU_MP_VERSION < 5
/* Older versions left out a normalization step */
extern void gcdext(mpz_t g, mpz_t s, mpz_t t, const mpz_t a, const mpz_t b);
#else
#define gcdext(g,s,t,a,b) mpz_gcdext(g,s,t,a,b)
#endif

extern int chinese(mpz_t ret, mpz_t lcm, mpz_t *a, mpz_t *m, int items);

extern UV mpz_order_ui(unsigned long r, mpz_t n, unsigned long limit);

extern void mpz_arctan(mpz_t r, unsigned long base, mpz_t pow, mpz_t t1, mpz_t t2);
extern void mpz_arctanh(mpz_t r, unsigned long base, mpz_t pow, mpz_t t1, mpz_t t2);
extern void mpz_product(mpz_t* A, UV a, UV b);
extern void mpz_product_ui(mpz_t prod, unsigned long *v, unsigned long n);
extern void mpz_veclcm(mpz_t* A, UV a, UV b);

extern void poly_mod_mul(mpz_t* px, mpz_t* py, UV r, mpz_t mod, mpz_t t1, mpz_t t2, mpz_t t3);
extern void poly_mod_pow(mpz_t *pres, mpz_t *pn, mpz_t power, UV r, mpz_t mod);

extern void poly_mod(mpz_t *pres, mpz_t *pn, UV *dn, mpz_t mod);
extern void polyz_mod(mpz_t *pres, mpz_t *pn, long *dn, mpz_t mod);

extern void polyz_set(mpz_t* pr, long* dr, mpz_t* ps, long ds);
extern void polyz_print(const char* header, mpz_t* pn, long dn);
extern void polyz_mulmod(mpz_t* pr, mpz_t* px, mpz_t *py, long *dr, long dx, long dy, mpz_t mod);
extern void polyz_div(mpz_t *pq, mpz_t *pr, mpz_t *pn, mpz_t *pd,
                      long *dq, long *dr, long dn, long dd, mpz_t NMOD);
extern void polyz_pow_polymod(mpz_t* pres,  mpz_t* pn,  mpz_t* pmod,
                              long *dres,   long   dn,  long   dmod,
                              mpz_t power, mpz_t NMOD);
extern void polyz_gcd(mpz_t* pres, mpz_t* pa, mpz_t* pb, long* dres, long da, long db, mpz_t MODN);

extern void polyz_root_deg1(mpz_t root, mpz_t* pn, mpz_t NMOD);
extern void polyz_root_deg2(mpz_t root1, mpz_t root2, mpz_t* pn, mpz_t NMOD);

/* Find roots of a polynomial mod a prime, slightly modified. */
/* We will stop if we've found at least maxroots unique roots. */
extern void polyz_roots_modp(mpz_t** roots, long *nroots, long maxroots,
                             mpz_t *pP, long dP, mpz_t NMOD);

/* Solve x^2 + |D|y^2 = p */
extern int cornacchia(mpz_t x, mpz_t y, mpz_t D, mpz_t p);
/* Solve x^2 + |D|y^2 = 4p */
extern int modified_cornacchia(mpz_t x, mpz_t y, mpz_t D, mpz_t p);

/* return a class poly (Hilbert [type 1] or Weber [type 2]) */
extern UV poly_class_poly(IV D, mpz_t**T, int* type);

extern const char* poly_class_type_name(int type);

/* return a 0 terminated list of all D's sorted by degree */
extern IV* poly_class_degrees(int insert_1s);

/* List of class polynomial indices in order */
extern int* poly_class_nums(void);
/* Given a class poly index, return the degree and fill in (if not null):
 *   D     the discriminant number
 *   T     the polynomial coefficients
 *   type  the poly type:  1 Hilber, 2 Weber
 */
extern UV poly_class_poly_num(int i, int *D, mpz_t**T, int* type);

#define BITS2DIGS(bits) ceil(bits/3.3219281)
#define DIGS2BITS(digs) ceil(digs*3.3219281)

extern void mpf_log(mpf_t logx, mpf_t x);
extern void mpf_exp(mpf_t expx, mpf_t x);
extern void mpf_pow(mpf_t powx, mpf_t b, mpf_t x);
extern void mpf_root(mpf_t rootx, mpf_t x, mpf_t n);
extern void mpf_agm(mpf_t r, mpf_t a, mpf_t b);

extern UV logint(mpz_t n, UV base);

#if defined(FUNC_isqrt) || defined(FUNC_is_perfect_square)
static UV isqrt(UV n) {
  UV root;
#if BITS_PER_WORD == 32
  if (n >= UVCONST(4294836225)) return UVCONST(65535);
#else
  if (n >= UVCONST(18446744065119617025)) return UVCONST(4294967295);
#endif
  root = (UV) sqrt((double)n);
  while (root*root > n)  root--;
  while ((root+1)*(root+1) <= n)  root++;
  return root;
}
#endif

#if defined(FUNC_gcd_ui) || defined(FUNC_lcm_ui)
static UV gcd_ui(UV x, UV y) {
  UV t;
  if (y < x) { t = x; x = y; y = t; }
  while (y > 0) {
    t = y;  y = x % y;  x = t;  /* y1 <- x0 % y0 ; x1 <- y0 */
  }
  return x;
}
#endif

#ifdef FUNC_lcm_ui
static UV lcm_ui(UV x, UV y) {
  /* Can overflow if lcm(x,y) > 2^64 (e.g. two primes each > 2^32) */
  return x * (y / gcd_ui(x,y));
}
#endif

#ifdef FUNC_is_perfect_square
/* Return 0 if n is not a perfect square.  Set sqrtn to int(sqrt(n)) if so.
 * See:  http://mersenneforum.org/showpost.php?p=110896
 */
static int is_perfect_square(UV n, UV* sqrtn)
{
  UV m;
  m = n & 127;
  if ((m*0x8bc40d7d) & (m*0xa1e2f5d1) & 0x14020a)  return 0;
  /* This cuts out another 80%: */
  m = n % 63; if ((m*0x3d491df7) & (m*0xc824a9f9) & 0x10f14008) return 0;
  /* m = n % 25; if ((m*0x1929fc1b) & (m*0x4c9ea3b2) & 0x51001005) return 0; */
  m = isqrt(n);
  if (n != m*m) return 0;
  if (sqrtn != 0) *sqrtn = m;
  return 1;
}
#endif

#ifdef FUNC_mpz_logn
static double mpz_logn(mpz_t n)
{
  long exp;
  double logn = mpz_get_d_2exp(&exp, n);
  logn = log(logn) + (log(2) * exp);
  return logn;
}
#endif

#ifdef FUNC_mpz_log2
static double mpz_log2(mpz_t n)
{
  long exp;
  double logn = mpz_get_d_2exp(&exp, n);
  logn = exp + log(logn)/log(2);
  return logn;
}
#endif

#endif
