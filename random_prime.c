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
  /* Trivial method.  We should use Fouque/Tibouchi Algorithm 1. */
  do {
    mpz_isaac_urandomb(p, n);
    mpz_setbit(p, n-1);
    if (n > 2) mpz_setbit(p, 0);
  } while (!_GMP_is_prob_prime(p));
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
