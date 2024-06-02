/*
 * Utility functions, such as sqrt mod p, polynomial manipulation, etc.
 */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <gmp.h>
#include <math.h>

#include "ptypes.h"

/* includes mpz_mulmod(r, a, b, n, temp) */
#include "utility.h"
#include "rootmod.h"
#include "factor.h"
#include "primality.h"
#include "isaac.h"

static int _verbose = 0;
int get_verbose_level(void) { return _verbose; }
void set_verbose_level(int level) { _verbose = level; }

static gmp_randstate_t _randstate;
gmp_randstate_t* get_randstate(void) { return &_randstate; }

#if __LITTLE_ENDIAN__ || (defined(BYTEORDER) && (BYTEORDER == 0x1234 || BYTEORDER == 0x12345678))
#define LESWAP(mem, val)
#else
#if !defined(__x86_64__)
#undef U8TO32_LE
#undef U32TO8_LE
#endif
#ifndef U32TO8_LE
#define U32TO8_LE(p, v) \
  do { \
    uint32_t _v = v; \
    (p)[0] = (((_v)      ) & 0xFFU); \
    (p)[1] = (((_v) >>  8) & 0xFFU); \
    (p)[2] = (((_v) >> 16) & 0xFFU); \
    (p)[3] = (((_v) >> 24) & 0xFFU); \
  } while (0)
#endif
#define LESWAP(mem, val) U32TO8_LE(mem,val)
#endif

void init_randstate(unsigned long seed) {
  unsigned char seedstr[8] = {0};
#if (__GNU_MP_VERSION > 4) || (__GNU_MP_VERSION == 4 && __GNU_MP_VERSION_MINOR >= 2)
  /* MT was added in GMP 4.2 released in 2006. */
  gmp_randinit_mt(_randstate);
#else
  gmp_randinit_default(_randstate);
#endif
  gmp_randseed_ui(_randstate, seed);

#if BITS_PER_WORD == 64
  if (seed > UVCONST(4294967295)) {
    LESWAP(seedstr, seed);
    LESWAP(seedstr + 4, (seed >> 16) >> 16);
    isaac_init(8, seedstr);
  } else
#endif
  {
    LESWAP(seedstr, seed);
    isaac_init(4, seedstr);
  }
}
void clear_randstate(void) {  gmp_randclear(_randstate);  }

UV irand64(int nbits)
{
  if (nbits ==  0) return 0;
  if (nbits <= 32) return( isaac_rand32() >> (32-nbits) );
#if BITS_PER_WORD == 64
  if (nbits <= 64) return( ((UV)isaac_rand32() | ((UV)isaac_rand32() << 32)) >> (64-nbits) );
#endif
  croak("irand64 too many bits for UV");
}
static NV _tonv_32 = -1.0;
static NV _tonv_64;
NV drand64(void)
{
  if (_tonv_32 < 0) {
    int i;
    NV t64, t32 = 1.0;
    for (i = 0; i < 32; i++)
      t32 /= 2.0;
    t64 = t32;
    for (i = 0; i < 32; i++)
      t64 /= 2.0;
    _tonv_64 = t64;
    _tonv_32 = t32;
  }
  return isaac_rand32() * _tonv_32 + isaac_rand32() * _tonv_64;
}

void mpz_isaac_urandomb(mpz_t rop, int nbits)
{
  if (nbits <= 32) {
    mpz_set_ui(rop, irand64(nbits));
  } else {
    unsigned char* d;
    int nbytes = (nbits+7)/8;
    New(0, d, nbytes, unsigned char);
    isaac_rand_bytes(nbytes, d);
    mpz_import(rop, nbytes, 1, sizeof(unsigned char), 0, 0, d);
    Safefree(d);
    if (nbits != nbytes*8)
      mpz_tdiv_r_2exp(rop, rop, nbits);
  }
}

void mpz_isaac_urandomm(mpz_t rop, mpz_t n)
{
  int count = 80;
  unsigned long nbits = mpz_sizeinbase(n,2);

  if (mpz_sgn(n) <= 0) {
    mpz_set_ui(rop,0);
    return;
  } else if (nbits <= 32) {
    mpz_set_ui(rop, isaac_rand(mpz_get_ui(n)));
  } else if (nbits < 3000) {
    /* Just like GMP, try until we're in range or we're tired. */
    while (count-- > 0) {
      mpz_isaac_urandomb(rop, nbits);
      if (mpz_cmp(rop,n) < 0)
        return;
    }
    mpz_mod(rop, rop, n);
  } else {
    /* Reduce tries needed by selecting from a range that is a multiple of n
     * (so no bias) and uses the max space inside the power-of-2 range.
     * Downside is that we do an alloc and two mods.  For large values
     * it can be much faster however. */
    mpz_t rmax;
    mpz_init(rmax);
    mpz_setbit(rmax, nbits+8);
    mpz_sub_ui(rmax,rmax,1);
    mpz_tdiv_q(rmax, rmax, n);
    mpz_mul(rmax, rmax, n);
    do {
      mpz_isaac_urandomb(rop, nbits+8);
    } while (mpz_cmp(rop, rmax) >= 0 && count-- > 0);
    mpz_clear(rmax);
    mpz_mod(rop, rop, n);
  }
}

int mpz_fits_uv_p(const mpz_t n)
{
  if (sizeof(UV) == sizeof(unsigned long int))
    return mpz_fits_ulong_p(n);
  else
    return (mpz_sgn(n) >= 0 && mpz_sizeinbase(n,2) <= BITS_PER_WORD);
}

void mpz_set_uv(mpz_t n, UV v)
{
#if BITS_PER_WORD == 32
  mpz_set_ui(n, v);
#else
  if (v <= 0xFFFFFFFFUL || sizeof(unsigned long int) >= sizeof(UV)) {
    mpz_set_ui(n, v);
  } else {
    mpz_set_ui(n, (v >> 32));
    mpz_mul_2exp(n, n, 32);
    mpz_add_ui(n, n, v & 0xFFFFFFFFUL);
  }
#endif
}
void mpz_set_iv(mpz_t n, IV v)
{
#if BITS_PER_WORD == 32
    mpz_set_si(n, v);
#else
  if ((v <= 0x7FFFFFFFL && v >= -0x80000000L) || sizeof(unsigned long int) >= sizeof(UV)) {
    mpz_set_si(n, v);
  } else if (v >= 0) {
    mpz_set_uv(n, v);
  } else {
    mpz_set_uv(n, -v);
    mpz_neg(n, n);
  }
#endif
}
UV mpz_get_uv(const mpz_t n)
{
#if BITS_PER_WORD == 32
  return mpz_get_ui(n);
#else
  UV v = mpz_getlimbn(n,0);
  if (GMP_LIMB_BITS < 64 || sizeof(mp_limb_t) < sizeof(UV))
    v |= ((UV)mpz_getlimbn(n,1)) << 32;
  return v;
#endif
}


/* a=0, return power.  a>1, return bool if an a-th power */
UV is_power(mpz_t n, UV a)
{
  if (mpz_cmp_ui(n,3) <= 0 && a == 0)
    return 0;
  else if (a == 1 || mpz_cmp_ui(n,1) == 0)
    return 1;
  else if (a == 2)
    return mpz_perfect_square_p(n);
  else if (a >= 2 && mpz_sizeinbase(n,2) < a)
    return 0;
  else {
    UV result;
    mpz_t t;
    mpz_init(t);
    result = (a == 0)  ?  power_factor(n, t)  :  (UV)mpz_root(t, n, a);
    mpz_clear(t);
    return result;
  }
}

UV prime_power(mpz_t prime, mpz_t n)
{
  UV k;
  if (mpz_even_p(n)) {
    k = mpz_scan1(n, 0);
    if (k+1 == mpz_sizeinbase(n, 2)) {
      mpz_set_ui(prime, 2);
      return k;
    }
    return 0;
  }
  if (_GMP_is_prob_prime(n)) {
    mpz_set(prime, n);
    return 1;
  }
  k = power_factor(n, prime);
  if (k && !_GMP_is_prob_prime(prime))
    k = 0;
  return k;
}

int is_primitive_root(mpz_t ina, mpz_t n, int nprime)
{
  mpz_t a, s, r, sreduced, t, *factors;
  int ret, i, nfactors, *exponents;

  if (mpz_sgn(n) == 0)
    return 0;
  if (mpz_sgn(n) < 0)
    mpz_neg(n,n);
  if (mpz_cmp_ui(n,1) == 0)
    return 1;
  if (mpz_cmp_ui(n,4) <= 0)
    return (mpz_cmp_ui(ina,4) <= 0) && (mpz_get_ui(ina) == mpz_get_ui(n)-1);
  if (mpz_divisible_2exp_p(n,2))
    return 0;

  mpz_init(a);
  mpz_mod(a,ina,n);
  mpz_init(s);
  mpz_gcd(s, a, n);

  if (mpz_cmp_ui(s,1) != 0) {
    mpz_clear(s);
    mpz_clear(a);
    return 0;
  }

  mpz_init(t);
  if (nprime) {
    mpz_sub_ui(s, n, 1);
  } else { /* totient(s, n); */   /* Fine, but slow. */
    UV k;
    mpz_init(r);
    if (mpz_odd_p(n)) mpz_set(t, n); else mpz_fdiv_q_2exp(t, n, 1);
    k = prime_power(r, t);
    if (!k) {  /* Not of form p^a or 2p^a */
      mpz_clear(r); mpz_clear(t); mpz_clear(s); mpz_clear(a); return 0;
    }
    mpz_divexact(t, t, r);
    mpz_mul(s, t, r);  mpz_sub(s, s, t);
    mpz_clear(r);
  }
  mpz_init_set(sreduced, s);

  ret = 0;
  mpz_sub_ui(t, n, 1);
  if (mpz_cmp(s,t) == 0 && mpz_kronecker(a,n) != -1)
    goto DONE_IPR;
  /* Unclear if this is worth doing.
  i = is_power(a, 0);
  if (i > 1 && mpz_gcd_ui(NULL, s, i) != 1)
    goto DONE_IPR;
  */

#define IPR_TEST_UI(s, p, a, n, t, ret) \
  mpz_divexact_ui(t, s, p); \
  mpz_powm(t, a, t, n); \
  if (mpz_cmp_ui(t, 1) == 0) { ret = 0; }

#define IPR_TEST(s, p, a, n, t, ret) \
  mpz_divexact(t, s, p); \
  mpz_powm(t, a, t, n); \
  if (mpz_cmp_ui(t, 1) == 0) { ret = 0; }

  ret = 1;
  /* Pull out small factors and test */
  {
    void *iter = trial_factor_iterator_create(sreduced, 64000);
    unsigned long f;
    uint32_t e;

    while (ret && mpz_cmp_ui(sreduced,1) > 0 && trial_factor_iterator_next(&f, &e, iter)) {
      IPR_TEST_UI(s, f, a, n, t, ret);
      while (e-- > 0)
        mpz_divexact_ui(sreduced,sreduced,f);
    }
    trial_factor_iterator_destroy(iter);
  }
  if (ret == 0 || mpz_cmp_ui(sreduced,1) == 0)
    goto DONE_IPR;
  if (_GMP_BPSW(sreduced)) {
    IPR_TEST(s, sreduced, a, n, t, ret);
    goto DONE_IPR;
  }

  /* We have a composite and so far it could be a primitive root.  Factor. */
  nfactors = factor(sreduced, &factors, &exponents);
  for (i = 0; ret == 1 && i < nfactors; i++) {
    IPR_TEST(s, factors[i], a, n, t, ret);
  }
  clear_factors(nfactors, &factors, &exponents);

DONE_IPR:
  mpz_clear(sreduced);  mpz_clear(t);
  mpz_clear(s);  mpz_clear(a);
  return ret;
}

int  is_qr(mpz_t a, mpz_t n)
{
  if (mpz_sgn(n) == 0) return -1;
  if (mpz_sgn(n) < 0) mpz_abs(n,n);
  if (mpz_cmp_ui(n,2) <= 0) return 1;
  if (mpz_cmp(a,n) > 0) mpz_mod(a,a,n);
  if (mpz_cmp_ui(a,1) <= 0) return 1;

  if (_GMP_is_prime(n))
    return (mpz_kronecker(a,n) == 1);

#if 0
  return sqrtmod(a,a,n);
#else
  {
    mpz_t r, *fac;
    int ret, i, nfactors, *exp;

    mpz_init(r);
    nfactors = factor(n, &fac, &exp);
    for (i = 0, ret = 1; ret && i < nfactors; i++) {
      int gcd_is_1, fac_is_2;
      mpz_gcd(r, a,fac[i]);
      gcd_is_1 = mpz_cmp_ui(r,1) == 0;
      fac_is_2 = mpz_cmp_ui(fac[i],2) == 0;
      if      (exp[i] == 1 && (fac_is_2 || !gcd_is_1))
        ret = 1;
      else if (exp[i] == 1 || (!fac_is_2 && gcd_is_1))
        ret = mpz_kronecker(a,fac[i]) == 1;
      else {
        mpz_pow_ui(fac[i], fac[i], exp[i]);
        ret = sqrtmod(r, a, fac[i]);
      }
    }
    clear_factors(nfactors, &fac, &exp);
    mpz_clear(r);
    return ret;
  }
#endif
}


int mpz_divmod(mpz_t r, mpz_t a, mpz_t b, mpz_t n, mpz_t t)
{
  if (mpz_invert(t, b, n)) {
    mpz_mulmod(r, t, a, n, t); /* mpz_mul(t,t,a); mpz_mod(r,t,n); */
    return 1;
  }
  return 0;  /* Cannot invert */
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

  sqrtmodp_t(x, D, p, a, b, c, d);
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

  sqrtmodp_t(x, D, p, a, b, c, d);

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

#if __GNU_MP_VERSION < 5
extern void gcdext(mpz_t g, mpz_t s, mpz_t t, const mpz_t ia, const mpz_t ib)
{
  mpz_t a, b;
  mpz_init_set(a, ia);
  mpz_init_set(b, ib);

  if (mpz_sgn(a) == 0 || mpz_cmp(a,b) == 0) {
    mpz_set_si(s, 0);
    mpz_set_si(t, mpz_sgn(b));
    mpz_abs(g, b);
  } else if (mpz_sgn(b) == 0) {
    mpz_set_si(s, mpz_sgn(a));
    mpz_set_si(t, 0);
    mpz_abs(g, a);
  } else {
    mpz_t os, ot, or, r, q, tmp;
    mpz_init(os);  mpz_init(ot);  mpz_init(or);
    mpz_init(r);  mpz_init(q);  mpz_init(tmp);
    mpz_set_ui(s,0);  mpz_set_ui(os,1);
    mpz_set_ui(t,1);  mpz_set_ui(ot,0);
    mpz_set(r,b);  mpz_set(or,a);
    while (mpz_sgn(r)) {
      mpz_tdiv_q(q, or, r);
      mpz_set(tmp, r); mpz_mul(r, r, q); mpz_sub(r, or, r); mpz_set(or, tmp);
      mpz_set(tmp, s); mpz_mul(s, s, q); mpz_sub(s, os, s); mpz_set(os, tmp);
      mpz_set(tmp, t); mpz_mul(t, t, q); mpz_sub(t, ot, t); mpz_set(ot, tmp);
    }
    mpz_set(s, os);
    mpz_set(t, ot);
    mpz_set(g, or);
    mpz_clear(r);  mpz_clear(q);  mpz_clear(tmp);
    mpz_clear(os);  mpz_clear(ot);  mpz_clear(or);

    if (mpz_sgn(g) < 0) {
      mpz_neg(s, s);
      mpz_neg(t, t);
      mpz_neg(g, g);
    }
  }
  mpz_clear(a); mpz_clear(b);
}
#endif

int chinese(mpz_t ret, mpz_t lcm, mpz_t *a, mpz_t *m, int items)
{
  mpz_t sum, gcd, u, v, s, t, temp1, temp2;
  int i, rval = 1;

#if 0
  if (items >= 128) {
    int first = items/2;
    mpz_t ca[2], cm[2];

    for (i = 0; i < 2; i++)
      { mpz_init(ca[i]); mpz_init(cm[i]); }
    rval = chinese(ca[0], cm[0], a, m, first);
    if (rval == 1)
      rval = chinese(ca[1], cm[1], a+first, m+first, items-first);
    if (rval == 1)
      rval = chinese(ret, lcm, ca, cm, 2);
    for (i = 0; i < 2; i++)
      { mpz_clear(ca[i]); mpz_clear(cm[i]); }
    return rval;
  }
#else
#define CRTN 8
  if (items >= 64) {
    int step = items/CRTN;
    mpz_t ca[CRTN], cm[CRTN];

    for (i = 0; i < CRTN; i++)
      { mpz_init(ca[i]); mpz_init(cm[i]); }
    for (i = 0; rval && i < CRTN; i++) {
      int citems = (i==CRTN-1) ? items-(CRTN-1)*step : step;
      rval = chinese(ca[i], cm[i], a+i*step, m+i*step, citems);
    }
    if (rval) rval = chinese(ret, lcm, ca, cm, CRTN);
    for (i = 0; i < CRTN; i++)
      { mpz_clear(ca[i]); mpz_clear(cm[i]); }
    return rval;
  }
#endif

  /* Avoid dividing by zero and use absolute value of modulus */
  for (i = 0; i < items; i++) {
    mpz_abs(m[i],m[i]);
    if (mpz_sgn(m[i]) == 0)
      return 0;
  }

  mpz_init(temp1); mpz_init(temp2);
  mpz_init(sum); mpz_init(gcd);
  mpz_init(s); mpz_init(t);
  mpz_init(u); mpz_init(v);

  mpz_set(lcm, m[0]);
  mpz_mod(sum, a[0], m[0]);
  for (i = 1; i < items; i++) {
    mpz_gcdext(gcd, u, v, lcm, m[i]);
    mpz_divexact(s, m[i], gcd);
    mpz_divexact(t, lcm, gcd);
    if (mpz_cmp_ui(gcd,1) != 0) {
      mpz_mod(temp1, sum, gcd);
      mpz_mod(temp2, a[i], gcd);
      if (mpz_cmp(temp1, temp2) != 0) {
        rval = 0;
        break;
      }
    }
    if (mpz_sgn(s) < 0) mpz_neg(s,s);
    if (mpz_sgn(t) < 0) mpz_neg(t,t);
    mpz_mul(lcm, lcm, s);
    if (mpz_sgn(u) < 0) mpz_add(u, u, lcm);
    if (mpz_sgn(v) < 0) mpz_add(v, v, lcm);
    mpz_mul(temp1, v, s);
    mpz_mul(v, temp1, sum);
    mpz_mul(temp1, u, t);
    mpz_mul(u, temp1, a[i]);
    mpz_add(temp1, v, u);
    mpz_mod(sum, temp1, lcm);
  }
  mpz_set(ret, sum);
  mpz_clear(sum); mpz_clear(gcd);
  mpz_clear(s); mpz_clear(t);
  mpz_clear(u); mpz_clear(v);
  mpz_clear(temp1); mpz_clear(temp2);
  return rval;
}


UV mpz_order_ui(unsigned long r, mpz_t n, unsigned long limit) {
  unsigned long j;
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
    if (base > 65535) { mpz_ui_pow_ui(t2, base, 2); mpz_tdiv_q(t1, t1, t2); }
    else               mpz_tdiv_q_ui(t1, t1, base*base);
    mpz_tdiv_q_ui(t2, t1, 2*i+1);
    if (i++ & 1) mpz_sub(r, r, t2); else mpz_add(r, r, t2);
  } while (mpz_sgn(t2));
}
void mpz_arctanh(mpz_t r, unsigned long base, mpz_t pow, mpz_t t1, mpz_t t2)
{
  unsigned long i = 1;
  mpz_tdiv_q_ui(r, pow, base);
  mpz_set(t1, r);
  do {
    if (base > 65535) { mpz_ui_pow_ui(t2, base, 2); mpz_tdiv_q(t1, t1, t2); }
    else               mpz_tdiv_q_ui(t1, t1, base*base);
    mpz_tdiv_q_ui(t2, t1, 1 + (i++ << 1));
    mpz_add(r, r, t2);
  } while (mpz_sgn(t2));
}

void mpz_product(mpz_t* A, UV a, UV b) {
  if (b <= a) {
    /* nothing */
  } else if (b == a+1) {
    mpz_mul(A[a], A[a], A[b]);
  } else if (b == a+2) {
    mpz_mul(A[a+1], A[a+1], A[a+2]);
    mpz_mul(A[a], A[a], A[a+1]);
  } else {
    UV c = a + (b-a+1)/2;
    mpz_product(A, a, c-1);
    mpz_product(A, c, b);
    mpz_mul(A[a], A[a], A[c]);
  }
}

void mpz_product_ui(mpz_t prod, unsigned long *v, unsigned long n) {
  mpz_set_ui(prod, 1);
  while (n > 0) {
    unsigned long p = v[--n];
    while (n > 0 && v[n-1] < ULONG_MAX/p)
      p *= v[--n];
    mpz_mul_ui(prod, prod, p);
  }
}

void mpz_veclcm(mpz_t* A, UV a, UV b) {
  if (b <= a) {
    /* nothing */
  } else if (b == a+1) {
    mpz_lcm(A[a], A[a], A[b]);
  } else if (b == a+2) {
    mpz_lcm(A[a+1], A[a+1], A[a+2]);
    mpz_lcm(A[a], A[a], A[a+1]);
  } else {
    UV c = a + (b-a+1)/2;
    mpz_veclcm(A, a, c-1);
    mpz_veclcm(A, c, b);
    mpz_lcm(A[a], A[a], A[c]);
  }
}

/* TODO: possible UV / unsigned long mismatch */
UV logint(mpz_t n, UV base) {
  mpz_t nt;
  double logn, logbn, coreps;
  UV res, nbits;

  if (mpz_cmp_ui(n,0) <= 0 || base <= 1)
    croak("mpz_logint: bad input\n");

  /* If base is a small power of 2, then this is exact */
  if (base <= 62 && (base & (base-1)) == 0)
    return mpz_sizeinbase(n, base)-1;

  if (mpz_cmp_ui(n,base) < 0)
    return 0;

#if 0  /* example using mpf_log for high precision.  Slow. */
  {
    mpf_t fr, fn;
    mpf_init(fr); mpf_init(fn);
    mpf_set_z(fn, n);
    mpf_log(fr, fn);
    logn = mpf_get_d(fr);
    mpf_clear(fn); mpf_clear(fr);
    coreps = 1e-8;
  }
#endif

  /* A typical way to do this is to start with base, then compare
   * base^2, base^4, base^8, ... until larger than n.  Then either work
   * back down or do a binary search
   * It uses O(log2(log2(n)) integer squares+multiplies plus some space.
   *
   * However, libc gives us the very fast log() function for doubles.  While
   * reducing the argument as needed to make sure we fit inside a double,
   * we can use this to give us a result extremely close to the right
   * answer, then adjust if we're not sure of the result.
   *
   * My benchmarks show it as about 2-10x faster than the all-integer method.
   */

  nbits = mpz_sizeinbase(n,2);
  mpz_init(nt);

  /* Step 1, get an approximation of log(n) */
  if (nbits < 512) {
    logn = log(mpz_get_d(n));
    coreps = 1e-8;
  } else {
    /* Reduce bits using log(x * 2^k) = log(x) + k * log(2) */
    uint32_t redbits = nbits - 256;
    mpz_tdiv_q_2exp(nt, n, redbits);
    logn = log(mpz_get_d(nt)) + redbits * 0.69314718055994530941723212145818L;
    coreps = 1e-7;
  }

  /* Step 2, approximate log_base(n) */
  logbn = logn / log(base);
  res = (UV) logbn;

  /* Step 3, correct if logbn might be rounded wrong */
  if (res != (UV)(logbn+coreps) || res != (UV)(logbn-coreps)) {
    mpz_ui_pow_ui(nt, base, res);
    if (mpz_cmp(nt, n) > 0) {
      res--;
    } else if (mpz_cmp(nt, n) < 0) {
      mpz_mul_ui(nt, nt, base);
      if (mpz_cmp(nt, n) <= 0)
        res++;
    }
  }
  mpz_clear(nt);
  /* res is largest res such that base^res <= n */
  return res;
}

/******************************************************************************/
/*
 * Floating point routines.
 * These are not rigorously accurate.  Use MPFR if possible.
 *
 * See also: http://fredrikj.net/math/elefun.pdf
 * for how to really look at this correctly.
 *
 * Many ideas from Brent's presentation:
 * https://pdfs.semanticscholar.org/8aec/ea97b8f2f23d4f09ec8f69025598f742ae9e.pdf
 */


extern void const_pi(mpf_t pi, unsigned long prec);
extern void const_log2(mpf_t logn, unsigned long prec);

/* Log using Brent's second AGM algorithm (Sasaki and Kanada theta) */
void mpf_log(mpf_t logn, mpf_t n)
{
  mpf_t N, t, q, theta2, theta3, logdn;
  unsigned long k, bits = mpf_get_prec(logn);
  int neg = (mpf_sgn(n) < 0);

  if (mpf_sgn(n) == 0)
    croak("mpf_log(0)");
  if (mpf_cmp_ui(n,2) == 0)
    { const_log2(logn,BITS2DIGS(bits)); return; }
  if ((neg && !mpf_cmp_si(n,-1)) || (!neg && !mpf_cmp_si(n,1)))
    { mpf_set_ui(logn,0); return; }

  mpf_init2(N, bits);
  mpf_set(N, n);
  if (neg) mpf_neg(N, N);

  mpf_init2(t, 64 + bits);
  mpf_init2(q,      64 + bits);
  mpf_init2(theta2, 64 + bits);
  mpf_init2(theta3, 64 + bits);
  mpf_init2(logdn,  64 + bits);
  mpf_set_ui(logn, 0);

  /* ensure N >> 1 */
  mpf_set_ui(t, 1);  mpf_mul_2exp(t, t, 1+(35+bits)/36);
  if (mpf_cmp(N, t) <= 0) {
    /* log(n) = log(n*2^k) - k*log(2) */
    for (k = 0; mpf_cmp(N, t) <= 0; k += 16)
      mpf_mul_2exp(N, N, 16);
    if (k > 0) {
      const_log2(t, BITS2DIGS(bits));
      mpf_mul_ui(logn, t, k);
      mpf_neg(logn, logn);
    }
  }

  mpf_ui_div(q, 1, N);
  mpf_pow_ui(t,q,9); mpf_add(theta2, q, t);
  mpf_pow_ui(t,q,25); mpf_add(theta2, theta2, t);
  mpf_mul_2exp(theta2, theta2, 1);
  mpf_pow_ui(theta3,q,4);
  mpf_pow_ui(t,q,16); mpf_add(theta3, theta3, t);
  mpf_mul_2exp(theta3, theta3, 1);
  mpf_add_ui(theta3, theta3, 1);

  /* Normally we would do:
   *   mpf_mul(theta2, theta2, theta2); mpf_mul(theta3, theta3, theta3);
   *   mpf_agm(t, theta2, theta3); mpf_mul_2exp(t, t, 2);
   * but Brent points out the one term acceleration:
   *   AGM(t2^2,t3^2)  =>  AGM(2*t2*t3,t2^2+t3^2)/2
   */
  mpf_mul(t, theta2, theta3);
  mpf_mul_2exp(q, t, 1);  /* q = 2*t2*t3 */
  mpf_add(t, theta2, theta3);
  mpf_mul(t, t, t);       /* t = (t2+t3)^2 = t2^2 + t3^2 + 2*t2*t3 */
  mpf_sub(theta3, t, q);
  mpf_set(theta2, q);
  mpf_agm(t, theta2, theta3);
  mpf_mul_2exp(t, t, 1);

  const_pi(logdn, BITS2DIGS(bits));
  mpf_div(logdn, logdn, t);

  mpf_add(logn, logn, logdn);
  mpf_clear(logdn); mpf_clear(theta3); mpf_clear(theta2); mpf_clear(q);
  mpf_clear(t); mpf_clear(N);
  if (neg) mpf_neg(logn, logn);
}

/* x should be 0 < x < 1 */
static void _exp_sinh(mpf_t expx, mpf_t x, unsigned long bits)
{
  unsigned long k;
  mpf_t t, s, N, D, X;

  mpf_init2(t, 10 + bits);
  mpf_init2(s, 10 + bits);
  mpf_init2(N, 10 + bits);
  mpf_init2(D, 10 + bits);
  mpf_init2(X, 10 + bits);

  /* 1. Compute s =~ sinh(x). */
  mpf_set(s,x);
  mpf_set(N,x);
  mpf_mul(X,x,x);
  mpf_set_ui(D,1);
  for (k = 1; k < bits; k++) {
    mpf_mul(N, N, X);
    mpf_mul_ui(D, D, 2*k);
    mpf_mul_ui(D, D, 2*k+1);
    mpf_div(t, N, D);
    mpf_add(s, s, t);

    mpf_abs(t, t);
    mpf_mul_2exp(t, t, bits);
    if (mpf_cmp_d(t, .5) < 0)
      break;
  }
  mpf_clear(X); mpf_clear(D); mpf_clear(N);

  /* 2. Compute s =~ e(x) from sinh(x). */
  mpf_mul(t, s, s);
  mpf_add_ui(t, t, 1);
  mpf_sqrt(t, t);
  mpf_add(s, s, t);

  mpf_set(expx, s);
  mpf_clear(s); mpf_clear(t);
}

static void _exp_lift(mpf_t expx, mpf_t x, unsigned long bits)
{
  mpf_t s, t1, t2;
  unsigned long k;

  mpf_init2(s,  10 + bits);
  mpf_init2(t1, 10 + bits);
  mpf_init2(t2, 10 + bits);

  mpf_set(s, expx);
  mpf_log(t1, s);
  mpf_sub(t2, x, t1);     /* t2 = s-ln(x) */
  mpf_mul(t1, s, t2);     /* t1 = s(s-ln(x) */
  mpf_add(s, s, t1);
  /* third and higher orders */
  for (k = 3; k <= 8; k++) {
    mpf_mul(t1, t1, t2);
    mpf_div_ui(t1, t1, k-1);
    mpf_add(s, s, t1);
  }
  mpf_set(expx, s);
  mpf_clear(t2); mpf_clear(t1); mpf_clear(s);
}

void mpf_exp(mpf_t expn, mpf_t x)
{
  mpf_t t;
  unsigned long k, r, rbits, bits = mpf_get_prec(expn);

  if (mpf_sgn(x) == 0) { mpf_set_ui(expn, 1); return; }

  mpf_init2(t, 10 + bits);

  if (mpf_sgn(x) < 0) { /* As per Brent, exp(x) = 1/exp(-x) */
    mpf_neg(t, x);
    mpf_exp(t, t);
    if (mpf_sgn(t) != 0) mpf_ui_div(expn, 1, t);
    else                 mpf_set_ui(expn, 0);
    mpf_clear(t);
    return;
  }

  /* Doubling rule, to make -.25 < x < .25.  Speeds convergence. */
  mpf_set(t, x);
  for (k = 0; mpf_cmp_d(t, 1.0L/8192.0L) > 0; k++)
    mpf_div_2exp(t, t, 1);

  /* exp with sinh method, then Newton lift to final bits */
  for (rbits = bits, r = 0;  rbits > 4000;  rbits = (rbits+7)/8)
    r++;
  _exp_sinh(expn, t, rbits);
  while (r-- > 0) {
    rbits *= 8;
    _exp_lift(expn, t, rbits);
  }
  if (rbits < bits)
    _exp_lift(expn, t, bits);

  if (k > 0) {
    const unsigned long maxred = 8*sizeof(unsigned long)-1;
    for ( ; k > maxred; k -= maxred)
      mpf_pow_ui(expn, expn, 1UL << maxred);
    mpf_pow_ui(expn, expn, 1UL << k);
  }
  mpf_clear(t);
}


/* Negative b with non-int x usually gives a complex result.
 * We try to at least give a consistent result. */

void mpf_pow(mpf_t powx, mpf_t b, mpf_t x)
{
  mpf_t t;
  int neg = (mpf_sgn(b) < 0);

  if (mpf_sgn(x) == 0) { mpf_set_ui(powx, 1); return; }
  if (mpf_sgn(b) == 0) { mpf_set_ui(powx, 0); return; }
  if (mpf_cmp_ui(b,1) == 0) { mpf_set_ui(powx, 1); return; }
  if (mpf_cmp_ui(x,1) == 0) { mpf_set(powx, b); return; }
  if (mpf_cmp_si(x,-1) == 0) { mpf_ui_div(powx, 1, b); return; }

  if (mpf_integer_p(x) && mpf_fits_ulong_p(x)) {
    mpf_pow_ui(powx, b, mpf_get_ui(x));
    return;
  }

  if (neg) mpf_neg(b,b);
  mpf_init2(t, mpf_get_prec(powx));
  mpf_log(t, b);
  mpf_mul(t, t, x);
  mpf_exp(powx, t);
  if (neg) mpf_neg(powx,powx);
  mpf_clear(t);
}

void mpf_root(mpf_t rootx, mpf_t x, mpf_t n)
{
  if (mpf_sgn(n) == 0) {
    mpf_set_ui(rootx, 0);
  } else if (mpf_cmp_ui(n, 2) == 0 && mpf_sgn(x) >= 0) {
    mpf_sqrt(rootx, x);
  } else {
    mpf_t t;
    mpf_init2(t, mpf_get_prec(rootx));
    mpf_ui_div(t, 1, n);
    mpf_pow(rootx, x, t);
    mpf_clear(t);
  }
}

void mpf_agm(mpf_t r, mpf_t a, mpf_t b)
{
  mpf_t t;
  unsigned long bits = mpf_get_prec(r);

  if (mpf_cmp(a,b) > 0) mpf_swap(a,b);

  mpf_init2(t, 6+bits);
  while (1) {
    mpf_sub(t, b, a);
    mpf_abs(t, t);
    mpf_mul_2exp(t, t, bits);
    mpf_sub_ui(t,t,1);
    if (mpf_sgn(t) < 0)
      break;
    mpf_sub_ui(t,t,1);
    mpf_set(t, a);
    mpf_add(a, a, b);
    mpf_div_2exp(a, a, 1);
    mpf_mul(b, b, t);
    mpf_sqrt(b, b);
  }
  mpf_set(r, b);
  mpf_clear(t);
}
