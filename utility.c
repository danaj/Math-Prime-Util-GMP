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
UV mpz_get_uv(mpz_t n)
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
  { /* Pull out small factors and test */
    UV p, fromp = 0;
    while (ret == 1 && (p = _GMP_trial_factor(sreduced, fromp, 60))) {
      if (mpz_cmp_ui(sreduced,p) == 0) break;
      IPR_TEST_UI(s, p, a, n, t, ret);
      mpz_set_ui(t, p);
      (void) mpz_remove(sreduced, sreduced, t);
      fromp = p+1;
    }
  }
  if (ret == 0 || mpz_cmp_ui(sreduced,1) == 0)
    goto DONE_IPR;
  if (_GMP_BPSW(sreduced)) {
    IPR_TEST(s, sreduced, a, n, t, ret);
    goto DONE_IPR;
  }

  /* Pull out more factors, noting they can be composites. */
  while (_GMP_pminus1_factor(sreduced, t, 100000, 100000)) {
    mpz_divexact(sreduced, sreduced, t);
    if (_GMP_BPSW(t)) {
      IPR_TEST(s, t, a, n, t, ret);
    } else {
      nfactors = factor(t, &factors, &exponents);
      for (i = 0; ret == 1 && i < nfactors; i++) {
        IPR_TEST(s, factors[i], a, n, t, ret);
      }
      clear_factors(nfactors, &factors, &exponents);
    }
    if (ret == 0 || mpz_cmp_ui(sreduced,1) == 0)
      goto DONE_IPR;
    if (_GMP_BPSW(sreduced)) {
      IPR_TEST(s, sreduced, a, n, t, ret);
      goto DONE_IPR;
    }
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


int mpz_divmod(mpz_t r, mpz_t a, mpz_t b, mpz_t n, mpz_t t)
{
  int invertible;
  invertible = mpz_invert(t, b, n);
  if (!invertible)
    return 0;
  mpz_mulmod(r, t, a, n, t); /* mpz_mul(t,t,a); mpz_mod(r,t,n); */
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

/******************************************************************************/


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

#if (__GNU_MP_VERSION < 4) || (__GNU_MP_VERSION == 4 && __GNU_MP_VERSION_MINOR < 1)
/* Binary segmentation, using simple shift+add method for processing p.
 * Faster than twiddling bits, but not nearly as fast as import/export.
 * mpz_import and mpz_export were added in GMP 4.1 released in 2002.
 */
void poly_mod_mul(mpz_t* px, mpz_t* py, UV r, mpz_t mod, mpz_t p, mpz_t p2, mpz_t t)
{
  UV i, d, bits;
  UV degree = r-1;

  mpz_mul(t, mod, mod);
  mpz_mul_ui(t, t, r);    /* TODO: possible UV / ui mismatch */
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

#else

/* Binary segmentation, using import/export method for processing p.
 * Thanks to Dan Bernstein's 2007 Quartic paper.
 */
void poly_mod_mul(mpz_t* px, mpz_t* py, UV r, mpz_t mod, mpz_t p, mpz_t p2, mpz_t t)
{
  UV i, bytes;
  char* s;

  mpz_mul(t, mod, mod);
  mpz_mul_ui(t, t, r);    /* TODO: possible UV / ui mismatch */
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
#endif

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
  r = *dr+1;    /* TODO: long to UV */
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
  sqrtmodp_t(e, d, NMOD, t, t2, t3, t4);

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
                        long maxroots, mpz_t* pg, long dg, mpz_t NMOD)
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
    else              mpz_isaac_urandomm(pxa[0], t);

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
      polyz_roots(roots, nroots, maxroots, ph, dh, NMOD);
      if (*nroots < maxroots) {
        /* q = g/h, and recurse */
        polyz_div(pq, pt,  pg, ph,  &dq, &dt, dg, dh, NMOD);
        polyz_roots(roots, nroots, maxroots, pq, dq, NMOD);
      }
    } else {
      polyz_div(pq, pt,  pg, ph,  &dq, &dt, dg, dh, NMOD);
      polyz_roots(roots, nroots, maxroots, pq, dq, NMOD);
      if (*nroots < maxroots) {
        polyz_roots(roots, nroots, maxroots, ph, dh, NMOD);
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
                      mpz_t *pP, long dP, mpz_t NMOD)
{
  long i;

  *nroots = 0;
  *roots = 0;

  if (dP == 0)
    return;

  /* Do degree 1 or 2 explicitly */
  if (dP == 1) {
    New(0, *roots, 1, mpz_t);
    mpz_init((*roots)[0]);
    polyz_root_deg1( (*roots)[0], pP, NMOD );
    *nroots = 1;
    return;
  }
  if (dP == 2) {
    New(0, *roots, 2, mpz_t);
    mpz_init((*roots)[0]);
    mpz_init((*roots)[1]);
    polyz_root_deg2( (*roots)[0], (*roots)[1], pP, NMOD );
    *nroots = 2;
    return;
  }

  /* Allocate space for the maximum number of roots (plus 1 for safety) */
  New(0, *roots, dP+1, mpz_t);
  for (i = 0; i <= dP; i++)
    mpz_init((*roots)[i]);

  if (maxroots > dP || maxroots == 0)
    maxroots = dP;

  polyz_roots(*roots, nroots, maxroots, pP, dP, NMOD);
  /* This could be just really bad luck.  Let the caller handle it. */
  /* if (*nroots == 0) croak("failed to find roots\n"); */

  /* Clear out space for roots we didn't find */
  for (i = *nroots; i <= dP; i++)
    mpz_clear((*roots)[i]);
}


#include "class_poly_data.h"

const char* poly_class_type_name(int type)
{
  switch (type) {
    case 1: return "Hilbert";
    case 2: return "Weber";
    case 3: return "Ramanujan";
    default: return "Unknown";
  }
}

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
