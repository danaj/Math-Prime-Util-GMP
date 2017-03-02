#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <gmp.h>
#include "ptypes.h"
#include "random_prime.h"
#include "utility.h"
#include "primality.h"
#include "gmp_main.h"
#include "isaac.h"

static char pr[31] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127};

void mpz_random_nbit_prime(mpz_t p, UV n)
{
  switch (n) {
    case 0:
    case 1:   mpz_set_ui(p, 0);  return;
    case 2:   mpz_set_ui(p, pr[ 0+isaac_rand(2)]);  return;
    case 3:   mpz_set_ui(p, pr[ 2+isaac_rand(2)]);  return;
    case 4:   mpz_set_ui(p, pr[ 4+isaac_rand(2)]);  return;
    case 5:   mpz_set_ui(p, pr[ 6+isaac_rand(5)]);  return;
    case 6:   mpz_set_ui(p, pr[11+isaac_rand(7)]);  return;
    case 7:   mpz_set_ui(p, pr[18+isaac_rand(13)]); return;
    default:  break;
  }
  /* For 32-bit inputs, use fast trivial method */
  if (n <= 32) {
    uint32_t mask = (0xFFFFFFFFU >> (34-n)) << 1,  base = mask+3;
    do {
      mpz_set_ui(p, base | (isaac_rand32() & mask));
    } while (!_GMP_is_prob_prime(p));
    return;
  }

#if 0
  do {                         /* Trivial method. */
    mpz_isaac_urandomb(p, n);
    mpz_setbit(p, n-1);
    mpz_setbit(p, 0);
  } while (!_GMP_is_prob_prime(p));
#else
  {                            /* Fouque+Tibouchi Alg 1, without modulo checks */
    mpz_t base;
    mpz_init(base);
    if (n > 33) { mpz_isaac_urandomb(base, n-33); mpz_mul_2exp(base,base,1); }
    mpz_setbit(base, n-1);
    mpz_setbit(base, 0);
    do {
      mpz_set_ui(p, isaac_rand32());
      mpz_mul_2exp(p, p, n-32);
      mpz_ior(p, p, base);
    } while (!_GMP_is_prob_prime(p));
    mpz_clear(base);
  }
#endif
}

/* PRIMEINC: pick random value, select next prime. */
/* Fast but bad distribution. */
static int _random_prime_primeinc(mpz_t p, mpz_t lo, mpz_t hi)
{
  mpz_t r, t;
  mpz_init(t);
  mpz_init(r);
  mpz_sub(r, hi, lo);
  mpz_isaac_urandomm(t, r);
  mpz_clear(r);
  mpz_add(t, t, lo);
  mpz_sub_ui(t, t, 1);
  _GMP_next_prime(t);
  if (mpz_cmp(t,hi) > 0) {
     mpz_sub_ui(t, lo, 1);
    _GMP_next_prime(t);
    if (mpz_cmp(t,hi) > 0) {
      mpz_clear(t);
      return 0;
    }
  }
  mpz_set(p, t);
  mpz_clear(t);
  return 1;
}

/* TRIVIAL: pick random values until one is prime */
/* Perfect distribution. */
static int _random_prime_trivial(mpz_t p, mpz_t lo_in, mpz_t hi_in)
{
  mpz_t r, t, lo, hi;
  int res = 0, tries = 10000;

  if (mpz_cmp_ui(hi_in,2) < 0 || mpz_cmp(lo_in,hi_in) > 0)
    return 0;

  mpz_init_set(lo, lo_in);
  mpz_init_set(hi, hi_in);
  if (mpz_cmp_ui(lo,2) <= 0) {
    mpz_set_ui(lo,1);
  } else if (mpz_even_p(lo)) {
    mpz_add_ui(lo,lo,1);
  }
  if (mpz_cmp_ui(hi,2) <= 0) {
    mpz_set_ui(hi,1);
  } else if (mpz_even_p(hi)) {
    mpz_sub_ui(hi,hi,1);
  }
  /* lo and hi are now odd */
  if (mpz_cmp(lo,hi) >= 0) {
    if (mpz_cmp(lo,hi) > 0) {
      /* null range */
    } else if (mpz_cmp_ui(lo,1) == 0) {
      mpz_set_ui(p,2);
      res = 1;
    } else if (_GMP_is_prob_prime(lo)) {
      mpz_set(p,lo);
      res = 1;
    }
    return res;
  }
  /* lo and hi are now odd and at least one odd between them */

  mpz_init(t);
  mpz_init(r);
  mpz_sub(r, hi, lo);
  mpz_tdiv_q_2exp(r, r, 1);
  mpz_add_ui(r,r,1);
  do {
    mpz_isaac_urandomm(t, r);
    mpz_mul_2exp(t, t, 1);
    mpz_add(t, t, lo);
    if (mpz_cmp_ui(t,1) == 0) mpz_set_ui(t,2);  /* map 1 back to 2 */
  } while (!_GMP_is_prob_prime(t) && --tries > 0);

  if (tries > 0) {
    mpz_set(p, t);
    res = 1;
  } else {
    /* We couldn't find anything.  Perhaps no primes in range. */
    res = _random_prime_primeinc(p, lo, hi);
  }
  mpz_clear(r);
  mpz_clear(t);
  mpz_clear(lo);
  mpz_clear(hi);
  return res;
}

/* Set p to a random prime between lo and hi inclusive */
int mpz_random_prime(mpz_t p, mpz_t lo, mpz_t hi)
{
  return _random_prime_trivial(p,lo,hi);
}

void mpz_random_ndigit_prime(mpz_t p, UV n)
{
  mpz_t lo, hi;
  switch (n) {
    case 0:   mpz_set_ui(p,0); return;
    case 1:   mpz_set_ui(p, pr[isaac_rand(4)]);  return;
    case 2:   mpz_set_ui(p, pr[4+isaac_rand(21)]);  return;
    default:  break;
  }
  mpz_init_set_ui(lo,10);
  mpz_pow_ui(lo, lo, n-1);
  mpz_init(hi);
  mpz_mul_ui(hi, lo, 10);

  if (!mpz_random_prime(p, lo, hi))
    croak("Failed to find %lu digit prime\n", n);

  mpz_clear(lo);
  mpz_clear(hi);
}

/* Random number rop such that 2*mult*rop+1 has nbits bits. */
static void _rand_in_bit_interval(mpz_t rop, UV nbits, mpz_t mult)
{
  mpz_t t, lo, hi;
  mpz_init(t); mpz_init(lo); mpz_init(hi);

  mpz_mul_ui(t, mult, 2);

  mpz_setbit(lo, nbits-1);
  mpz_sub_ui(lo, lo, 1);
  mpz_cdiv_q(lo, lo, t);   /* lo = ceil(2^(nbits-1)-1 / (2*mult)) */

  mpz_setbit(hi, nbits);
  mpz_sub_ui(hi, hi, 2);
  mpz_fdiv_q(hi, hi, t);   /* hi = floor(2^nbits-2 / (2*mult)) */

  mpz_sub(t, hi, lo);
  mpz_isaac_urandomm(rop, t);
  mpz_add(rop, rop, lo);

  mpz_clear(t); mpz_clear(lo); mpz_clear(hi);
}

/* Gordon's algorithm */
void mpz_random_strong_prime(mpz_t p, UV nbits)
{
  mpz_t S, T, R, P0, t, i, j;
  UV rbits, sbits, tbits;

  if (nbits < 128)  croak("random_strong_prime, bits must be >= 128");

  if (nbits < 256) {
    rbits = ((nbits+1) >> 1) - 2;
    sbits = (nbits >> 1) - 20;
    tbits = rbits - 20;
  } else {
    UV N1, N2;
    { /* Calculate FIPS 186-4 C.10 recommended parameter */
      UV t_, l2_;
      for (l2_ = 1, t_ = nbits; t_ >>= 1; ) l2_++;
      N1 =  (nbits/2)-l2_-7;
      N2 = N1/2;
    }
    if (N1 > 200) N1 = 201;
    if (N2 > 100) N2 = 101;
    if (N2 < 100) N2 += N1/4;
    rbits = sbits = N1;
    tbits = N2;
  }

  mpz_init(S);  mpz_init(T);  mpz_init(R);  mpz_init(P0);
  mpz_init(t);  mpz_init(i);  mpz_init(j);

  while (1) {
    mpz_random_nbit_prime(S, sbits);
    mpz_random_nbit_prime(T, tbits);

    _rand_in_bit_interval(i, rbits, T);
    while (1) {
      mpz_mul(t, i, T);
      mpz_mul_ui(t, t, 2);
      mpz_add_ui(R, t, 1);                 /* R = 2*i*T+1 */
      if (_GMP_is_prob_prime(R)) break;
      mpz_add_ui(i,i,1);
    }

    mpz_sub_ui(t, R, 2);
    mpz_powm(P0, S, t, R);
    mpz_mul_ui(P0, P0, 2);
    mpz_mul(P0, P0, S);
    mpz_sub_ui(P0, P0, 1);

    mpz_mul(i, R, S);
    mpz_mul_ui(t, i, 2);
    _rand_in_bit_interval(j, nbits, i);
    while (1) {
      mpz_mul(p, j, t);
      mpz_add(p, p, P0);                 /* p = 2*j*R*S+p0 */
      if (mpz_sizeinbase(p,2) > nbits) break;
      if (_GMP_is_prob_prime(p)) {
        mpz_clear(t);  mpz_clear(i);  mpz_clear(j);
        mpz_clear(S);  mpz_clear(T);  mpz_clear(R);  mpz_clear(P0);
        /* p-1 has factor R.  p+1 has factor S.  r-1 has factor T. */
        return;
      }
      mpz_add_ui(j,j,1);
    }
  }
}
