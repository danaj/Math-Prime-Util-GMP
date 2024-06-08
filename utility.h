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
extern void mpz_isaac_urandomm(mpz_t rop, const mpz_t n);
extern UV irand64(int nbits);
extern NV drand64(void);

extern int  mpz_fits_uv_p(const mpz_t n);
extern void mpz_set_uv(mpz_t n, UV v);
extern void mpz_set_iv(mpz_t n, IV v);
extern UV   mpz_get_uv(const mpz_t n);
/* TODO: cmp_uv, cmp_iv, mul_iv, mul_uv, etc. */

extern UV   is_power(mpz_t n, UV a);
extern UV   prime_power(mpz_t prime, mpz_t n);
extern int  is_primitive_root(mpz_t a, mpz_t b, int nprime);
extern int  is_qr(mpz_t a, mpz_t n);

/* tdiv_r is faster, but we'd need to guarantee the input is positive */
#define mpz_mulmod(r, a, b, n, t)  \
  do { mpz_mul(t, a, b); mpz_mod(r, t, n); } while (0)

#undef mpz_divmod
extern int mpz_divmod(mpz_t r, const mpz_t a, const mpz_t b, const mpz_t n, mpz_t t);

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

/* Solve x^2 + |D|y^2 = p */
extern int cornacchia(mpz_t x, mpz_t y, const mpz_t D, const mpz_t p);
/* Solve x^2 + |D|y^2 = 4p */
extern int modified_cornacchia(mpz_t x, mpz_t y, const mpz_t D, const mpz_t p);

#define BITS2DIGS(bits) ceil(bits/3.3219281)
#define DIGS2BITS(digs) ceil(digs*3.3219281)

extern void mpf_log(mpf_t logx, mpf_t x);
extern void mpf_exp(mpf_t expx, mpf_t x);
extern void mpf_pow(mpf_t powx, mpf_t b, mpf_t x);
extern void mpf_root(mpf_t rootx, mpf_t x, mpf_t n);
extern void mpf_agm(mpf_t r, mpf_t a, mpf_t b);

extern UV logint(const mpz_t n, UV base);

#ifdef FUNC_mpz_logn
static double mpz_logn(const mpz_t n)
{
  long exp;
  double logn = mpz_get_d_2exp(&exp, n);
  logn = log(logn) + (log(2) * exp);
  return logn;
}
#endif

#ifdef FUNC_mpz_log2
static double mpz_log2(const mpz_t n)
{
  long exp;
  double logn = mpz_get_d_2exp(&exp, n);
  logn = exp + log(logn)/log(2);
  return logn;
}
#endif

#endif
