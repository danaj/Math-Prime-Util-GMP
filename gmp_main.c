
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>

#define FUNC_gcd_ui 1
#define FUNC_is_perfect_square 1
#include "ptypes.h"
#include "gmp_main.h"
#include "prime_iterator.h"
#include "bls75.h"
#include "ecpp.h"
#include "utility.h"
#include "factor.h"

#define AKS_VARIANT_V6        1    /* The V6 paper with Lenstra impr */
#define AKS_VARIANT_BERN21    2    /* AKS-Bernstein-Morain theorem 2.1 */
#define AKS_VARIANT_BERN22    3    /* AKS-Bernstein-Morain theorem 2.2 */
#define AKS_VARIANT_BERN23    4    /* AKS-Bernstein-Morain theorem 2.3 */
#define AKS_VARIANT_BORNEMANN 5    /* Based on Folkmar Bornemann's impl */
#define AKS_VARIANT_BERN41    6    /* Bernstein 2003, theorem 4.1 */

#define AKS_VARIANT  AKS_VARIANT_BERN41

static mpz_t _bgcd;
static mpz_t _bgcd2;
static mpz_t _bgcd3;
#define BGCD_PRIMES       168
#define BGCD_LASTPRIME    997
#define BGCD_NEXTPRIME   1009
#define BGCD2_PRIMES     1229
#define BGCD2_NEXTPRIME 10007
#define BGCD3_PRIMES     4203
#define BGCD3_NEXTPRIME 40009

#define NSMALLPRIMES 168
static const unsigned short sprimes[NSMALLPRIMES] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997};

#define TSTAVAL(arr, val)   (arr[(val) >> 6] & (1U << (((val)>>1) & 0x1F)))
#define SETAVAL(arr, val)   arr[(val) >> 6] |= 1U << (((val)>>1) & 0x1F)


void _GMP_init(void)
{
  /* We should  not use this random number system for crypto, so
   * using this lousy seed is ok.  We just would like something a
   * bit different every run.  Using Perl_seed(aTHX) would be better. */
  unsigned long seed = time(NULL);
  init_randstate(seed);
  prime_iterator_global_startup();
  mpz_init(_bgcd);
  _GMP_pn_primorial(_bgcd, BGCD_PRIMES);   /* mpz_primorial_ui(_bgcd, 1000) */
  mpz_init_set_ui(_bgcd2, 0);
  mpz_init_set_ui(_bgcd3, 0);
  _init_factor();
}

void _GMP_destroy(void)
{
  prime_iterator_global_shutdown();
  clear_randstate();
  mpz_clear(_bgcd);
  mpz_clear(_bgcd2);
  mpz_clear(_bgcd3);
  destroy_ecpp_gcds();
}


static const unsigned char next_wheel[30] =
  {1,7,7,7,7,7,7,11,11,11,11,13,13,17,17,17,17,19,19,23,23,23,23,29,29,29,29,29,29,1};
static const unsigned char prev_wheel[30] =
  {29,29,1,1,1,1,1,1,7,7,7,7,11,11,13,13,13,13,17,17,19,19,19,19,23,23,23,23,23,23};
static const unsigned char wheel_advance[30] =
  {1,6,5,4,3,2,1,4,3,2,1,2,1,4,3,2,1,2,1,4,3,2,1,6,5,4,3,2,1,2};
static const unsigned char wheel_retreat[30] =
  {1,2,1,2,3,4,5,6,1,2,3,4,1,2,1,2,3,4,1,2,1,2,3,4,1,2,3,4,5,6};


static INLINE int _GMP_miller_rabin_ui(mpz_t n, UV base)
{
  int rval;
  mpz_t a;
  mpz_init_set_ui(a, base);
  rval = _GMP_miller_rabin(n, a);
  mpz_clear(a);
  return rval;
}

int _GMP_miller_rabin_random(mpz_t n, UV numbases, char* seedstr)
{
  gmp_randstate_t* p_randstate = get_randstate();
  mpz_t t, base;
  UV i;

  if (numbases == 0)  return 1;
  if (mpz_cmp_ui(n, 100) < 0)     /* tiny n */
    return (_GMP_is_prob_prime(n) > 0);

  mpz_init(base);  mpz_init(t);

  if (seedstr != 0) { /* Set the RNG seed if they gave us a seed */
    mpz_set_str(t, seedstr, 0);
    gmp_randseed(*p_randstate, t);
  }

  mpz_sub_ui(t, n, 3);
  for (i = 0; i < numbases; i++) {
    mpz_urandomm(base, *p_randstate, t);  /* base = 0 .. (n-3)-1 */
    mpz_add_ui(base, base, 2);            /* base = 2 .. n-2     */
    if (_GMP_miller_rabin(n, base) == 0)
      break;
  }
  mpz_clear(base);  mpz_clear(t);
  return (i >= numbases);
}

int _GMP_miller_rabin(mpz_t n, mpz_t a)
{
  mpz_t nminus1, d, x;
  UV s, r;
  int rval;

  {
    int cmpr = mpz_cmp_ui(n, 2);
    if (cmpr == 0)     return 1;  /* 2 is prime */
    if (cmpr < 0)      return 0;  /* below 2 is composite */
    if (mpz_even_p(n)) return 0;  /* multiple of 2 is composite */
  }
  if (mpz_cmp_ui(a, 1) <= 0)
    croak("Base %ld is invalid", mpz_get_si(a));
  mpz_init_set(nminus1, n);
  mpz_sub_ui(nminus1, nminus1, 1);
  mpz_init_set(x, a);

  /* Handle large and small bases.  Use x so we don't modify their input a. */
  if (mpz_cmp(x, n) >= 0)
    mpz_mod(x, x, n);
  if ( (mpz_cmp_ui(x, 1) <= 0) || (mpz_cmp(x, nminus1) >= 0) ) {
    mpz_clear(nminus1);
    mpz_clear(x);
    return 1;
  }

  mpz_init_set(d, nminus1);
  s = mpz_scan1(d, 0);
  mpz_tdiv_q_2exp(d, d, s);

  mpz_powm(x, x, d, n);
  mpz_clear(d); /* done with a and d */
  rval = 0;
  if (!mpz_cmp_ui(x, 1) || !mpz_cmp(x, nminus1)) {
    rval = 1;
  } else {
    for (r = 1; r < s; r++) {
      mpz_powm_ui(x, x, 2, n);
      if (!mpz_cmp_ui(x, 1)) {
        break;
      }
      if (!mpz_cmp(x, nminus1)) {
        rval = 1;
        break;
      }
    }
  }
  mpz_clear(nminus1); mpz_clear(x);
  return rval;
}

int is_miller_prime(mpz_t n, int assume_grh)
{
  mpz_t nminus1, d, x;
  UV s, r, logn, maxa, a;
  int rval;

  {
    int cmpr = mpz_cmp_ui(n, 2);
    if (cmpr == 0)     return 1;  /* 2 is prime */
    if (cmpr < 0)      return 0;  /* below 2 is composite */
    if (mpz_even_p(n)) return 0;  /* multiple of 2 is composite */
  }
  mpz_init_set(nminus1, n);
  mpz_sub_ui(nminus1, nminus1, 1);

  mpz_init_set(d, nminus1);
  s = mpz_scan1(d, 0);
  mpz_tdiv_q_2exp(d, d, s);
  logn = mpz_sizeinbase(n, 2);

  if (mpz_cmp_ui(n, 1373653) < 0) {
    maxa = 3;
  } else if (assume_grh) {
    /* Rough computation of ln(n), slightly high */
    logn = 1 + (0.69314718056 * ((double)logn));
    if (logn > ((BITS_PER_WORD == 32) ? 46340UL : 3037000499UL))
      croak("is_miller_prime: n is too large for GRH DMR");
    maxa = 2*logn*logn;
  } else { /* Bober and Goldmakher 2015 (http://arxiv.org/abs/1311.7556) */
    logn = 1 + ((1/6.59) * ((double)logn));
    if (logn >= BITS_PER_WORD)
      croak("is_miller_prime: n is too large for unconditional DMR");
    maxa = UVCONST(1) << logn;
  }
  if (mpz_cmp_ui(n, maxa) <= 0)
    maxa = mpz_get_ui(n) - 1;
  if (get_verbose_level() > 1)
    printf("Deterministic Miller-Rabin testing bases from 2 to %"UVuf"\n", maxa);

  mpz_init(x);
  rval = 1;
  for (a = 2; rval && a <= maxa; a++) {
    rval = 0;
    mpz_set_ui(x, a);
    mpz_powm(x, x, d, n);
    if (!mpz_cmp_ui(x, 1) || !mpz_cmp(x, nminus1)) {
      rval = 1;
    } else {
      for (r = 1; r < s; r++) {
        mpz_powm_ui(x, x, 2, n);
        if (!mpz_cmp_ui(x, 1)) {
          break;
        }
        if (!mpz_cmp(x, nminus1)) {
          rval = 1;
          break;
        }
      }
    }
  }
  mpz_clear(nminus1); mpz_clear(x); mpz_clear(d);
  return rval;
}

/* Returns Lucas sequence  U_k mod n and V_k mod n  defined by P,Q */
void lucas_seq(mpz_t U, mpz_t V, mpz_t n, IV P, IV Q, mpz_t k,
               mpz_t Qk, mpz_t t)
{
  UV b = mpz_sizeinbase(k, 2);
  IV D = P*P - 4*Q;

  if (mpz_cmp_ui(n, 2) < 0) croak("Lucas sequence modulus n must be > 1");
  MPUassert( mpz_cmp_ui(k, 0) >= 0, "lucas_seq: k is negative" );
  MPUassert( mpz_cmp_si(n,(P>=0) ? P : -P) > 0, "lucas_seq: P is out of range");
  MPUassert( mpz_cmp_si(n,(Q>=0) ? Q : -Q) > 0, "lucas_seq: Q is out of range");
  MPUassert( D != 0, "lucas_seq: D is zero" );

  if (mpz_cmp_ui(k, 0) <= 0) {
    mpz_set_ui(U, 0);
    mpz_set_ui(V, 2);
    return;
  }
  if (mpz_even_p(n)) {
    /* If n is even, then we can't divide by 2.  Do it differently. */
    alt_lucas_seq(U, V, n, P, Q, k, Qk, t);
    return;
  }

  mpz_set_ui(U, 1);
  mpz_set_si(V, P);
  mpz_set_si(Qk, Q);

  if (Q == 1) {
    /* Use the fast V method if possible.  Much faster with small n. */
    mpz_set_si(t, P*P-4);
    if (P > 2 && mpz_invert(t, t, n)) {
      /* Compute V_k and V_{k+1}, then computer U_k from them. */
      mpz_set_si(V, P);
      mpz_set_si(U, P*P-2);
      while (b > 1) {
        b--;
        if (mpz_tstbit(k, b-1)) {
          mpz_mul(V, V, U);  mpz_sub_ui(V, V, P);  mpz_mod(V, V, n);
          mpz_mul(U, U, U);  mpz_sub_ui(U, U, 2);  mpz_mod(U, U, n);
        } else {
          mpz_mul(U, V, U);  mpz_sub_ui(U, U, P);  mpz_mod(U, U, n);
          mpz_mul(V, V, V);  mpz_sub_ui(V, V, 2);  mpz_mod(V, V, n);
        }
      }
      mpz_mul_ui(U, U, 2);
      mpz_submul_ui(U, V, P);
      mpz_mul(U, U, t);
    } else {
      /* Fast computation of U_k and V_k, specific to Q = 1 */
      while (b > 1) {
        mpz_mulmod(U, U, V, n, t);     /* U2k = Uk * Vk */
        mpz_mul(V, V, V);
        mpz_sub_ui(V, V, 2);
        mpz_mod(V, V, n);               /* V2k = Vk^2 - 2 Q^k */
        b--;
        if (mpz_tstbit(k, b-1)) {
          mpz_mul_si(t, U, D);
                                      /* U:  U2k+1 = (P*U2k + V2k)/2 */
          mpz_mul_si(U, U, P);
          mpz_add(U, U, V);
          if (mpz_odd_p(U)) mpz_add(U, U, n);
          mpz_fdiv_q_2exp(U, U, 1);
                                      /* V:  V2k+1 = (D*U2k + P*V2k)/2 */
          mpz_mul_si(V, V, P);
          mpz_add(V, V, t);
          if (mpz_odd_p(V)) mpz_add(V, V, n);
          mpz_fdiv_q_2exp(V, V, 1);
        }
      }
    }
  } else {
    while (b > 1) {
      mpz_mulmod(U, U, V, n, t);     /* U2k = Uk * Vk */
      mpz_mul(V, V, V);
      mpz_submul_ui(V, Qk, 2);
      mpz_mod(V, V, n);               /* V2k = Vk^2 - 2 Q^k */
      mpz_mul(Qk, Qk, Qk);            /* Q2k = Qk^2 */
      b--;
      if (mpz_tstbit(k, b-1)) {
        mpz_mul_si(t, U, D);
                                    /* U:  U2k+1 = (P*U2k + V2k)/2 */
        mpz_mul_si(U, U, P);
        mpz_add(U, U, V);
        if (mpz_odd_p(U)) mpz_add(U, U, n);
        mpz_fdiv_q_2exp(U, U, 1);
                                    /* V:  V2k+1 = (D*U2k + P*V2k)/2 */
        mpz_mul_si(V, V, P);
        mpz_add(V, V, t);
        if (mpz_odd_p(V)) mpz_add(V, V, n);
        mpz_fdiv_q_2exp(V, V, 1);

        mpz_mul_si(Qk, Qk, Q);
      }
      mpz_mod(Qk, Qk, n);
    }
  }
  mpz_mod(U, U, n);
  mpz_mod(V, V, n);
}

void alt_lucas_seq(mpz_t Uh, mpz_t Vl, mpz_t n, IV P, IV Q, mpz_t k,
                   mpz_t Ql, mpz_t t)
{
  mpz_t Vh, Qh;
  int j, s = mpz_scan1(k,0), b = mpz_sizeinbase(k,2);

  if (mpz_sgn(k) <= 0) {
    mpz_set_ui(Uh, 0);
    mpz_set_ui(Vl, 2);
    return;
  }

  mpz_set_ui(Uh, 1);
  mpz_set_ui(Vl, 2);
  mpz_set_ui(Ql, 1);
  mpz_init_set_si(Vh,P);
  mpz_init_set_ui(Qh,1);

  for (j = b; j > s; j--) {
    mpz_mul(Ql, Ql, Qh);
    if (mpz_tstbit(k, j)) {
      mpz_mul_si(Qh, Ql, Q);
      mpz_mul(Uh, Uh, Vh);
      mpz_mul_si(t, Ql, P);  mpz_mul(Vl, Vl, Vh); mpz_sub(Vl, Vl, t);
      mpz_mul(Vh, Vh, Vh); mpz_sub(Vh, Vh, Qh); mpz_sub(Vh, Vh, Qh);
    } else {
      mpz_set(Qh, Ql);
      mpz_mul(Uh, Uh, Vl);  mpz_sub(Uh, Uh, Ql);
      mpz_mul_si(t, Ql, P);  mpz_mul(Vh, Vh, Vl); mpz_sub(Vh, Vh, t);
      mpz_mul(Vl, Vl, Vl);  mpz_sub(Vl, Vl, Ql);  mpz_sub(Vl, Vl, Ql);
    }
    mpz_mod(Qh, Qh, n);
    mpz_mod(Uh, Uh, n);
    mpz_mod(Vh, Vh, n);
    mpz_mod(Vl, Vl, n);
  }
  mpz_mul(Ql, Ql, Qh);
  mpz_mul_si(Qh, Ql, Q);
  mpz_mul(Uh, Uh, Vl);  mpz_sub(Uh, Uh, Ql);
  mpz_mul_si(t, Ql, P);  mpz_mul(Vl, Vl, Vh);  mpz_sub(Vl, Vl, t);
  mpz_mul(Ql, Ql, Qh);
  mpz_clear(Qh);  mpz_clear(Vh);
  mpz_mod(Ql, Ql, n);
  mpz_mod(Uh, Uh, n);
  mpz_mod(Vl, Vl, n);
  for (j = 0; j < s; j++) {
    mpz_mul(Uh, Uh, Vl);
    mpz_mul(Vl, Vl, Vl);  mpz_sub(Vl, Vl, Ql);  mpz_sub(Vl, Vl, Ql);
    mpz_mul(Ql, Ql, Ql);
    mpz_mod(Ql, Ql, n);
    mpz_mod(Uh, Uh, n);
    mpz_mod(Vl, Vl, n);
  }
}

void lucasuv(mpz_t Uh, mpz_t Vl, IV P, IV Q, mpz_t k)
{
  mpz_t Vh, Ql, Qh, t;
  int j, s, n;

  if (mpz_sgn(k) <= 0) {
    mpz_set_ui(Uh, 0);
    mpz_set_ui(Vl, 2);
    return;
  }

  mpz_set_ui(Uh, 1);
  mpz_set_ui(Vl, 2);
  mpz_init_set_si(Vh,P);
  mpz_init(t);

  s = mpz_scan1(k, 0);     /* number of zero bits at the end */
  n = mpz_sizeinbase(k,2);

  /* It is tempting to try to pull out the various Q operations when Q=1 or
   * Q=-1.  This doesn't lead to any immediate savings.  Don't bother unless
   * there is a way to reduce the actual operations involving U and V. */
  mpz_init_set_ui(Ql,1);
  mpz_init_set_ui(Qh,1);

  for (j = n; j > s; j--) {
    mpz_mul(Ql, Ql, Qh);
    if (mpz_tstbit(k, j)) {
      mpz_mul_si(Qh, Ql, Q);
      mpz_mul(Uh, Uh, Vh);
      mpz_mul_si(t, Ql, P);  mpz_mul(Vl, Vl, Vh); mpz_sub(Vl, Vl, t);
      mpz_mul(Vh, Vh, Vh); mpz_sub(Vh, Vh, Qh); mpz_sub(Vh, Vh, Qh);
    } else {
      mpz_set(Qh, Ql);
      mpz_mul(Uh, Uh, Vl);  mpz_sub(Uh, Uh, Ql);
      mpz_mul_si(t, Ql, P);  mpz_mul(Vh, Vh, Vl); mpz_sub(Vh, Vh, t);
      mpz_mul(Vl, Vl, Vl);  mpz_sub(Vl, Vl, Ql);  mpz_sub(Vl, Vl, Ql);
    }
  }
  mpz_mul(Ql, Ql, Qh);
  mpz_mul_si(Qh, Ql, Q);
  mpz_mul(Uh, Uh, Vl);  mpz_sub(Uh, Uh, Ql);
  mpz_mul_si(t, Ql, P);  mpz_mul(Vl, Vl, Vh);  mpz_sub(Vl, Vl, t);
  mpz_mul(Ql, Ql, Qh);
  mpz_clear(Qh);  mpz_clear(t);  mpz_clear(Vh);
  for (j = 0; j < s; j++) {
    mpz_mul(Uh, Uh, Vl);
    mpz_mul(Vl, Vl, Vl);  mpz_sub(Vl, Vl, Ql);  mpz_sub(Vl, Vl, Ql);
    mpz_mul(Ql, Ql, Ql);
  }
  mpz_clear(Ql);
}

int lucas_lehmer(UV p)
{
  UV k, tlim;
  int res, pbits;
  mpz_t V, mp, t;

  if (p == 2) return 1;
  if (!(p&1)) return 0;

  mpz_init_set_ui(t, p);
  if (!_GMP_is_prob_prime(t))    /* p must be prime */
    { mpz_clear(t); return 0; }
  if (p < 23)
    { mpz_clear(t); return (p != 11); }

  pbits = mpz_sizeinbase(t,2);
  mpz_init(mp);
  mpz_setbit(mp, p);
  mpz_sub_ui(mp, mp, 1);

  /* If p=3 mod 4 and p,2p+1 both prime, then 2p+1 | 2^p-1.  Cheap test. */
  if (p > 3 && p % 4 == 3) {
    mpz_mul_ui(t, t, 2);
    mpz_add_ui(t, t, 1);
    if (_GMP_is_prob_prime(t) && mpz_divisible_p(mp, t))
      { mpz_clear(mp); mpz_clear(t); return 0; }
  }

  /* Do a little trial division first.  Saves quite a bit of time. */
  tlim = (p < 1500) ? p/2 : (p < 5000) ? p : 2*p;
  if (tlim > UV_MAX/(2*p)) tlim = UV_MAX/(2*p);
  for (k = 1; k < tlim; k++) {
    UV q = 2*p*k+1;
    if ( (q%8==1 || q%8==7) &&                 /* factor must be 1 or 7 mod 8 */
         q % 3 && q % 5 && q % 7 && q % 11 && q % 13) {  /* and must be prime */
      if (1 && q < (UVCONST(1) << (BITS_PER_WORD/2)) ) {
        UV b = 1, j = pbits;
        while (j--) {
          b = (b*b) % q;
          if (p & (UVCONST(1) << j)) { b *= 2; if (b >= q) b -= q; }
        }
        if (b == 1)
          { mpz_clear(mp); mpz_clear(t); return 0; }
      } else {
        if( mpz_divisible_ui_p(mp, q) )
          { mpz_clear(mp); mpz_clear(t); return 0; }
      }
    }
  }
  /* We could do some specialized p+1 factoring here. */

  mpz_init_set_ui(V, 4);

  for (k = 3; k <= p; k++) {
    mpz_mul(V, V, V);
    mpz_sub_ui(V, V, 2);
    /* mpz_mod(V, V, mp) but more efficiently done given mod 2^p-1 */
    if (mpz_sgn(V) < 0) mpz_add(V, V, mp);
    /* while (n > mp) { n = (n >> p) + (n & mp) } if (n==mp) n=0 */
    /* but in this case we can have at most one loop plus a carry */
    mpz_tdiv_r_2exp(t, V, p);
    mpz_tdiv_q_2exp(V, V, p);
    mpz_add(V, V, t);
    while (mpz_cmp(V, mp) >= 0) mpz_sub(V, V, mp);
  }
  res = !mpz_sgn(V);
  mpz_clear(t); mpz_clear(mp); mpz_clear(V);
  return res;
}

/* Returns:  -1 unknown, 0 composite, 2 definitely prime */
int llr(mpz_t N)
{
  mpz_t v, k, V, U, Qk, t;
  UV i, n, P;
  int res = -1;

  if (mpz_cmp_ui(N,100) <= 0) return (_GMP_is_prob_prime(N) ? 2 : 0);
  if (mpz_even_p(N) || mpz_divisible_ui_p(N, 3)) return 0;
  mpz_init(v); mpz_init(k);
  mpz_add_ui(v, N, 1);
  n = mpz_scan1(v, 0);
  mpz_tdiv_q_2exp(k, v, n);
  /* N = k * 2^n - 1 */
  if (mpz_cmp_ui(k,1) == 0) {
    res = lucas_lehmer(n) ? 2 : 0;
    goto DONE_LLR;
  }
  if (mpz_sizeinbase(k,2) > n)
    goto DONE_LLR;

  mpz_init(V);
  mpz_init(U); mpz_init(Qk); mpz_init(t);
  if (!mpz_divisible_ui_p(k, 3)) { /* Select V for 3 not divis k */
    lucas_seq(U, V, N, 4, 1, k, Qk, t);
  } else if ((n % 4 == 0 || n % 4 == 3) && mpz_cmp_ui(k,3)==0) {
    mpz_set_ui(V, 5778);
  } else {
    /* Öystein J. Rödseth: http://www.uib.no/People/nmaoy/papers/luc.pdf */
    for (P=5; P < 1000; P++) {
      mpz_set_ui(t, P-2);
      if (mpz_jacobi(t, N) == 1) {
        mpz_set_ui(t, P+2);
        if (mpz_jacobi(t, N) == -1) {
          break;
        }
      }
    }
    if (P >= 1000) {
      mpz_clear(t);  mpz_clear(Qk); mpz_clear(U);
      mpz_clear(V);
      goto DONE_LLR;
    }
    lucas_seq(U, V, N, P, 1, k, Qk, t);
  }
  mpz_clear(t);  mpz_clear(Qk); mpz_clear(U);

  for (i = 3; i <= n; i++) {
    mpz_mul(V, V, V);
    mpz_sub_ui(V, V, 2);
    mpz_mod(V, V, N);
  }
  res = mpz_sgn(V) ? 0 : 2;
  mpz_clear(V);

DONE_LLR:
  if (res != -1 && get_verbose_level() > 1) printf("N shown %s with LLR\n", res ? "prime" : "composite");
  mpz_clear(k); mpz_clear(v);
  return res;
}
/* Returns:  -1 unknown, 0 composite, 2 definitely prime */
int proth(mpz_t N)
{
  mpz_t v, k, a;
  UV n;
  int i, res = -1;
  if (mpz_cmp_ui(N,100) <= 0) return (_GMP_is_prob_prime(N) ? 2 : 0);
  if (mpz_even_p(N) || mpz_divisible_ui_p(N, 3)) return 0;
  mpz_init(v); mpz_init(k);
  mpz_sub_ui(v, N, 1);
  n = mpz_scan1(v, 0);
  mpz_tdiv_q_2exp(k, v, n);
  /* N = k * 2^n + 1 */
  if (mpz_sizeinbase(k,2) <= n) {
    mpz_init(a);
    for (i = 0; i < 25 && res == -1; i++) {
      mpz_set_ui(a, sprimes[i]);
      if (mpz_jacobi(a, N) == -1)
        res = 0;
    }
    /* (a,N) = -1 if res=0.  Do Proth primality test. */
    if (res == 0) {
      mpz_tdiv_q_2exp(k, v, 1);
      mpz_powm(a, a, k, N);
      if (mpz_cmp(a, v) == 0)
        res = 2;
    }
    mpz_clear(a);
  }
  mpz_clear(k); mpz_clear(v);
  if (res != -1 && get_verbose_level() > 1) printf("N shown %s with Proth\n", res ? "prime" : "composite");
  fflush(stdout);
  return res;
}
static int is_proth_form(mpz_t N)
{
  mpz_t v, k;
  UV n;
  int res = 0;
  if (mpz_cmp_ui(N,100) <= 0) return (_GMP_is_prob_prime(N) ? 2 : 0);
  if (mpz_even_p(N) || mpz_divisible_ui_p(N, 3)) return 0;
  mpz_init(v); mpz_init(k);
  mpz_sub_ui(v, N, 1);
  n = mpz_scan1(v, 0);
  mpz_tdiv_q_2exp(k, v, n);
  /* N = k * 2^n + 1 */
  if (mpz_sizeinbase(k,2) <= n)  res = 1;
  mpz_clear(k); mpz_clear(v);
  return res;
}

static int lucas_selfridge_params(IV* P, IV* Q, mpz_t n, mpz_t t)
{
  IV D = 5;
  UV Dui = (UV) D;
  while (1) {
    UV gcd = mpz_gcd_ui(NULL, n, Dui);
    if ((gcd > 1) && mpz_cmp_ui(n, gcd) != 0)
      return 0;
    mpz_set_si(t, D);
    if (mpz_jacobi(t, n) == -1)
      break;
    if (Dui == 21 && mpz_perfect_square_p(n))
      return 0;
    Dui += 2;
    D = (D > 0)  ?  -Dui  :  Dui;
    if (Dui > 1000000)
      croak("lucas_selfridge_params: D exceeded 1e6");
  }
  if (P) *P = 1;
  if (Q) *Q = (1 - D) / 4;
  return 1;
}

static int lucas_extrastrong_params(IV* P, IV* Q, mpz_t n, mpz_t t, UV inc)
{
  UV tP = 3;
  if (inc < 1 || inc > 256)
    croak("Invalid lucas parameter increment: %"UVuf"\n", inc);
  while (1) {
    UV D = tP*tP - 4;
    UV gcd = mpz_gcd_ui(NULL, n, D);
    if (gcd > 1 && mpz_cmp_ui(n, gcd) != 0)
      return 0;
    mpz_set_ui(t, D);
    if (mpz_jacobi(t, n) == -1)
      break;
    if (tP == (3+20*inc) && mpz_perfect_square_p(n))
      return 0;
    tP += inc;
    if (tP > 65535)
      croak("lucas_extrastrong_params: P exceeded 65535");
  }
  if (P) *P = (IV)tP;
  if (Q) *Q = 1;
  return 1;
}



/* This code was verified against Feitsma's psps-below-2-to-64.txt file.
 * is_strong_pseudoprime reduced it from 118,968,378 to 31,894,014.
 * all three variations of the Lucas test reduce it to 0.
 * The test suite should check that they generate the correct pseudoprimes.
 *
 * The standard and strong versions use the method A (Selfridge) parameters,
 * while the extra strong version uses Baillie's parameters from OEIS A217719.
 *
 * Using the strong version, we can implement the strong BPSW test as
 * specified by Baillie and Wagstaff, 1980, page 1401.
 *
 * Testing on my x86_64 machine, the strong Lucas code is over 35% faster than
 * T.R. Nicely's implementation, and over 40% faster than David Cleaver's.
 */
int _GMP_is_lucas_pseudoprime(mpz_t n, int strength)
{
  mpz_t d, U, V, Qk, t;
  IV P, Q;
  UV s = 0;
  int rval;
  int _verbose = get_verbose_level();

  {
    int cmpr = mpz_cmp_ui(n, 2);
    if (cmpr == 0)     return 1;  /* 2 is prime */
    if (cmpr < 0)      return 0;  /* below 2 is composite */
    if (mpz_even_p(n)) return 0;  /* multiple of 2 is composite */
  }

  mpz_init(t);
  rval = (strength < 2) ? lucas_selfridge_params(&P, &Q, n, t)
                        : lucas_extrastrong_params(&P, &Q, n, t, 1);
  if (!rval) {
    mpz_clear(t);
    return 0;
  }
  if (_verbose>3) gmp_printf("N: %Zd  D: %"IVdf"  P: %"UVuf"  Q: %"IVdf"\n", n, P*P-4*Q, P, Q);

  mpz_init(U);  mpz_init(V);  mpz_init(Qk);
  mpz_init_set(d, n);
  mpz_add_ui(d, d, 1);

  if (strength > 0) {
    s = mpz_scan1(d, 0);
    mpz_tdiv_q_2exp(d, d, s);
  }

  lucas_seq(U, V, n, P, Q, d, Qk, t);
  mpz_clear(d);

  rval = 0;
  if (strength == 0) {
    /* Standard checks U_{n+1} = 0 mod n. */
    rval = (mpz_sgn(U) == 0);
  } else if (strength == 1) {
    if (mpz_sgn(U) == 0) {
      rval = 1;
    } else {
      while (s--) {
        if (mpz_sgn(V) == 0) {
          rval = 1;
          break;
        }
        if (s) {
          mpz_mul(V, V, V);
          mpz_submul_ui(V, Qk, 2);
          mpz_mod(V, V, n);
          mpz_mulmod(Qk, Qk, Qk, n, t);
        }
      }
    }
  } else {
    mpz_sub_ui(t, n, 2);
    if ( mpz_sgn(U) == 0 && (mpz_cmp_ui(V, 2) == 0 || mpz_cmp(V, t) == 0) ) {
      rval = 1;
    } else {
      s--;  /* The extra strong test tests r < s-1 instead of r < s */
      while (s--) {
        if (mpz_sgn(V) == 0) {
          rval = 1;
          break;
        }
        if (s) {
          mpz_mul(V, V, V);
          mpz_sub_ui(V, V, 2);
          mpz_mod(V, V, n);
        }
      }
    }
  }
  mpz_clear(Qk); mpz_clear(V); mpz_clear(U); mpz_clear(t);
  return rval;
}

/* Pari's clever method.  It's an extra-strong Lucas test, but without
 * computing U_d.  This makes it faster, but yields more pseudoprimes.
 *
 * increment:  1 for Baillie OEIS, 2 for Pari.
 *
 * With increment = 1, these results will be a subset of the extra-strong
 * Lucas pseudoprimes.  With increment = 2, we produce Pari's results (we've
 * added the necessary GCD with D so we produce somewhat fewer).
 */
int _GMP_is_almost_extra_strong_lucas_pseudoprime(mpz_t n, UV increment)
{
  mpz_t d, V, W, t;
  UV P, s;
  int rval;

  {
    int cmpr = mpz_cmp_ui(n, 2);
    if (cmpr == 0)     return 1;  /* 2 is prime */
    if (cmpr < 0)      return 0;  /* below 2 is composite */
    if (mpz_even_p(n)) return 0;  /* multiple of 2 is composite */
  }

  mpz_init(t);
  {
    IV PP;
    if (! lucas_extrastrong_params(&PP, 0, n, t, increment) ) {
      mpz_clear(t);
      return 0;
    }
    P = (UV) PP;
  }

  mpz_init(d);
  mpz_add_ui(d, n, 1);

  s = mpz_scan1(d, 0);
  mpz_tdiv_q_2exp(d, d, s);

  /* Calculate V_d */
  {
    UV b = mpz_sizeinbase(d, 2);
    mpz_init_set_ui(V, P);
    mpz_init_set_ui(W, P*P-2);   /* V = V_{k}, W = V_{k+1} */

    while (b > 1) {
      b--;
      if (mpz_tstbit(d, b-1)) {
        mpz_mul(V, V, W);
        mpz_sub_ui(V, V, P);

        mpz_mul(W, W, W);
        mpz_sub_ui(W, W, 2);
      } else {
        mpz_mul(W, V, W);
        mpz_sub_ui(W, W, P);

        mpz_mul(V, V, V);
        mpz_sub_ui(V, V, 2);
      }
      mpz_mod(V, V, n);
      mpz_mod(W, W, n);
    }
    mpz_clear(W);
  }
  mpz_clear(d);

  rval = 0;
  mpz_sub_ui(t, n, 2);
  if ( mpz_cmp_ui(V, 2) == 0 || mpz_cmp(V, t) == 0 ) {
    rval = 1;
  } else {
    s--;  /* The extra strong test tests r < s-1 instead of r < s */
    while (s--) {
      if (mpz_sgn(V) == 0) {
        rval = 1;
        break;
      }
      if (s) {
        mpz_mul(V, V, V);
        mpz_sub_ui(V, V, 2);
        mpz_mod(V, V, n);
      }
    }
  }
  mpz_clear(V); mpz_clear(t);
  return rval;
}

static void mat_mulmod_3x3(mpz_t* a, mpz_t* b, mpz_t n, mpz_t* t, mpz_t t2) {
  int i, row, col;
  for (row = 0; row < 3; row++) {
    for (col = 0; col < 3; col++) {
      mpz_mul(t[3*row+col], a[3*row+0], b[0+col]);
      mpz_mul(t2, a[3*row+1], b[3+col]);
      mpz_add(t[3*row+col], t[3*row+col], t2);
      mpz_mul(t2, a[3*row+2], b[6+col]);
      mpz_add(t[3*row+col], t[3*row+col], t2);
    }
  }
  for (i = 0; i < 9; i++) mpz_mod(a[i], t[i], n);
}
static void mat_powmod_3x3(mpz_t* m, mpz_t kin, mpz_t n) {
  mpz_t k, t2, t[9], res[9];
  int i;
  mpz_init_set(k, kin);
  mpz_init(t2);
  for (i = 0; i < 9; i++) { mpz_init(t[i]); mpz_init(res[i]); }
  mpz_set_ui(res[0],1);  mpz_set_ui(res[4],1);  mpz_set_ui(res[8],1);
  while (mpz_sgn(k)) {
    if (mpz_odd_p(k))  mat_mulmod_3x3(res, m, n, t, t2);
    mpz_fdiv_q_2exp(k, k, 1);
    if (mpz_sgn(k))    mat_mulmod_3x3(m, m, n, t, t2);
  }
  for (i = 0; i < 9; i++)
    { mpz_set(m[i],res[i]);  mpz_clear(res[i]);  mpz_clear(t[i]); }
  mpz_clear(t2);
  mpz_clear(k);
}
int is_perrin_pseudoprime(mpz_t n)
{
  int P[9] = {0,1,0, 0,0,1, 1,1,0};
  mpz_t m[9];
  int i, rval;
  {
    int cmpr = mpz_cmp_ui(n, 2);
    if (cmpr == 0)     return 1;  /* 2 is prime */
    if (cmpr < 0)      return 0;  /* below 2 is composite */
  }
  for (i = 0; i < 9; i++) mpz_init_set_ui(m[i], P[i]);
  mat_powmod_3x3(m, n, n);
  mpz_add(m[1], m[0], m[4]);
  mpz_add(m[2], m[1], m[8]);
  mpz_mod(m[0], m[2], n);
  rval = mpz_sgn(m[0]) ? 0 : 1;
  for (i = 0; i < 9; i++) mpz_clear(m[i]);
  return rval;
}

int is_frobenius_pseudoprime(mpz_t n, IV P, IV Q)
{
  mpz_t t, Vcomp, d, U, V, Qk;
  IV D;
  int k = 0;
  int rval;

  {
    int cmpr = mpz_cmp_ui(n, 2);
    if (cmpr == 0)     return 1;  /* 2 is prime */
    if (cmpr < 0)      return 0;  /* below 2 is composite */
    if (mpz_even_p(n)) return 0;  /* multiple of 2 is composite */
  }
  mpz_init(t);
  if (P == 0 && Q == 0) {
    P = 1;  Q = 2;
    do {
      P += 2;
      if (P == 3) P = 5;  /* P=3,Q=2 -> D=9-8=1 => k=1, so skip */
      if (P == 21 && mpz_perfect_square_p(n))
        { mpz_clear(t); return 0; }
      D = P*P-4*Q;
      if (mpz_cmp_ui(n, P >= 0 ? P : -P) <= 0) break;
      if (mpz_cmp_ui(n, D >= 0 ? D : -D) <= 0) break;
      mpz_set_si(t, D);
      k = mpz_jacobi(t, n);
    } while (k == 1);
  } else {
    D = P*P-4*Q;
    if (is_perfect_square( D >= 0 ? D : -D, 0 ))
      croak("Frobenius invalid P,Q: (%"IVdf",%"IVdf")", P, Q);
    mpz_set_si(t, D);
    k = mpz_jacobi(t, n);
  }

  /* Check initial conditions */
  {
    UV Pu = P >= 0 ? P : -P;
    UV Qu = Q >= 0 ? Q : -Q;
    UV Du = D >= 0 ? D : -D;

    /* If abs(P) or abs(Q) or abs(D) >= n, exit early. */
    if (mpz_cmp_ui(n, Pu) <= 0 || mpz_cmp_ui(n, Qu) <= 0 || mpz_cmp_ui(n, Du) <= 0) {
      mpz_clear(t);
      return _GMP_trial_factor(n, 2, Du+Pu+Qu) ? 0 : 1;
    }
    /* If k = 0, then D divides n */
    if (k == 0) {
      mpz_clear(t);
      return 0;
    }
    /* If n is not coprime to P*Q*D then we found a factor */
    if (mpz_gcd_ui(NULL, n, Du*Pu*Qu) > 1) {
      mpz_clear(t);
      return 0;
    }
  }

  mpz_init(Vcomp);
  if (k == 1) {
    mpz_set_si(Vcomp, 2);
  } else {
    mpz_set_si(Vcomp, Q);
    mpz_mul_ui(Vcomp, Vcomp, 2);
    mpz_mod(Vcomp, Vcomp, n);
  }

  mpz_init(U);  mpz_init(V);  mpz_init(Qk);  mpz_init(d);
  if (k == 1) mpz_sub_ui(d, n, 1);
  else        mpz_add_ui(d, n, 1);

  lucas_seq(U, V, n, P, Q, d, Qk, t);
  rval = ( mpz_sgn(U) == 0 && mpz_cmp(V, Vcomp) == 0 );

  mpz_clear(d); mpz_clear(Qk); mpz_clear(V); mpz_clear(U);
  mpz_clear(Vcomp); mpz_clear(t);

  return rval;
}

/* Use Crandall/Pomerance, steps from Loebenberger 2008 */
int is_frobenius_cp_pseudoprime(mpz_t n, UV ntests)
{
  gmp_randstate_t* p_randstate = get_randstate();
  mpz_t t, a, b, d, w1, wm, wm1, m;
  UV tn;
  int j;
  int result = 1;

  if (mpz_cmp_ui(n, 100) < 0)     /* tiny n */
    return (_GMP_is_prob_prime(n) > 0);
  if (mpz_even_p(n))
    return 0;

  mpz_init(t);  mpz_init(a);  mpz_init(b);  mpz_init(d);
  mpz_init(w1);  mpz_init(wm);  mpz_init(wm1);  mpz_init(m);
  for (tn = 0; tn < ntests; tn++) {
    /* Step 1: choose a and b in 1..n-1 and d=a^2-4b not square and coprime */
    do {
      mpz_sub_ui(t, n, 1);
      mpz_urandomm(a, *p_randstate, t);
      mpz_add_ui(a, a, 1);
      mpz_urandomm(b, *p_randstate, t);
      mpz_add_ui(b, b, 1);
      /* Check d and gcd */
      mpz_mul(d, a, a);
      mpz_mul_ui(t, b, 4);
      mpz_sub(d, d, t);
    } while (mpz_perfect_square_p(d));
    mpz_mul(t, a, b);
    mpz_mul(t, t, d);
    mpz_gcd(t, t, n);
    if (mpz_cmp_ui(t, 1) != 0 && mpz_cmp(t, n) != 0)
      { result = 0; break; }
    /* Step 2: W1 = a^2b^-1 - 2 mod n */
    if (!mpz_invert(t, b, n))
      { result = 0; break; }
    mpz_mul(w1, a, a);
    mpz_mul(w1, w1, t);
    mpz_sub_ui(w1, w1, 2);
    mpz_mod(w1, w1, n);
    /* Step 3: m = (n-(d|n))/2 */
    j = mpz_jacobi(d, n);
    if (j == -1)     mpz_add_ui(m, n, 1);
    else if (j == 0) mpz_set(m, n);
    else if (j == 1) mpz_sub_ui(m, n, 1);
    mpz_tdiv_q_2exp(m, m, 1);
    /* Step 8 here:  B = b^(n-1)/2 mod n  (stored in d) */
    mpz_sub_ui(t, n, 1);
    mpz_tdiv_q_2exp(t, t, 1);
    mpz_powm(d, b, t, n);
    /* Quick Euler test */
    mpz_sub_ui(t, n, 1);
    if (mpz_cmp_ui(d, 1) != 0 && mpz_cmp(d, t) != 0)
      { result = 0; break; }
    /* Step 4: calculate Wm,Wm+1 */
    mpz_set_ui(wm, 2);
    mpz_set(wm1, w1);
    {
      UV bit = mpz_sizeinbase(m, 2);
      while (bit-- > 0) {
        if (mpz_tstbit(m, bit)) {
          mpz_mul(t, wm, wm1);
          mpz_sub(wm, t, w1);
          mpz_mul(t, wm1, wm1);
          mpz_sub_ui(wm1, t, 2);
        } else {
          mpz_mul(t, wm, wm1);
          mpz_sub(wm1, t, w1);
          mpz_mul(t, wm, wm);
          mpz_sub_ui(wm, t, 2);
        }
        mpz_mod(wm, wm, n);
        mpz_mod(wm1, wm1, n);
      }
    }
    /* Step 5-7: compare w1/wm */
    mpz_mul(t, w1, wm);
    mpz_mod(t, t, n);
    mpz_mul_ui(wm1, wm1, 2);
    mpz_mod(wm1, wm1, n);
    if (mpz_cmp(t, wm1) != 0)
      { result = 0; break; }
    /* Step 8 was earlier */
    /* Step 9: See if Bwm = 2 mod n */
    mpz_mul(wm, wm, d);
    mpz_mod(wm, wm, n);
    if (mpz_cmp_ui(wm, 2) != 0)
      { result = 0; break; }
  }
  mpz_clear(w1);  mpz_clear(wm);  mpz_clear(wm1);  mpz_clear(m);
  mpz_clear(t);  mpz_clear(a);  mpz_clear(b);  mpz_clear(d);
  return result;
}

/* New code based on draft paper */
int _GMP_is_frobenius_underwood_pseudoprime(mpz_t n)
{
  mpz_t temp1, temp2, n_plus_1, s, t;
  unsigned long a, ap2, len;
  int bit, j, rval = 0;
  int _verbose = get_verbose_level();

  {
    int cmpr = mpz_cmp_ui(n, 2);
    if (cmpr == 0)     return 1;  /* 2 is prime */
    if (cmpr < 0)      return 0;  /* below 2 is composite */
    if (mpz_even_p(n)) return 0;  /* multiple of 2 is composite */
  }

  mpz_init(temp1);

  for (a = 0; a < 1000000; a++) {
    if (a==2 || a==4 || a==7 || a==8 || a==10 || a==14 || a==16 || a==18)
      continue;
    mpz_set_si(temp1, (long)(a*a) - 4);
    j = mpz_jacobi(temp1, n);
    if (j == -1) break;
    if (j == 0 || (a == 20 && mpz_perfect_square_p(n)))
      { mpz_clear(temp1); return 0; }
  }
  if (a >= 1000000)
    { mpz_clear(temp1); croak("FU test failure, unable to find suitable a"); }
  if (mpz_gcd_ui(NULL, n, (a+4)*(2*a+5)) != 1)
    { mpz_clear(temp1); return 0; }

  mpz_init(temp2); mpz_init(n_plus_1),

  ap2 = a+2;
  mpz_add_ui(n_plus_1, n, 1);
  len = mpz_sizeinbase(n_plus_1, 2);
  mpz_init_set_ui(s, 1);
  mpz_init_set_ui(t, 2);

  for (bit = len-2; bit >= 0; bit--) {
    mpz_add(temp2, t, t);
    if (a != 0) {
      mpz_mul_ui(temp1, s, a);
      mpz_add(temp2, temp1, temp2);
    }
    mpz_mul(temp1, temp2, s);
    mpz_sub(temp2, t, s);
    mpz_add(s, s, t);
    mpz_mul(t, s, temp2);
    mpz_mod(t, t, n);
    mpz_mod(s, temp1, n);
    if ( mpz_tstbit(n_plus_1, bit) ) {
      if (a == 0)   mpz_add(temp1, s, s);
      else          mpz_mul_ui(temp1, s, ap2);
      mpz_add(temp1, temp1, t);
      mpz_add(temp2, t, t);
      mpz_sub(t, temp2, s);
      mpz_set(s, temp1);
    }
  }
  /* n+1 always has an even last bit, so s and t always modded */
  mpz_set_ui(temp1, 2*a+5);
  mpz_mod(temp1, temp1, n);
  if (mpz_cmp_ui(s, 0) == 0 && mpz_cmp(t, temp1) == 0)
    rval = 1;
  if (_verbose>1) gmp_printf("%Zd is %s with a = %"UVuf"\n", n, (rval) ? "probably prime" : "composite", a);

  mpz_clear(temp1); mpz_clear(temp2); mpz_clear(n_plus_1);
  mpz_clear(s); mpz_clear(t);
  return rval;
}

int _GMP_is_frobenius_khashin_pseudoprime(mpz_t n)
{
  mpz_t t, ta, tb, ra, rb, a, b, n_minus_1;
  unsigned long c = 1;
  int bit, len, k, rval = 0;

  {
    int cmpr = mpz_cmp_ui(n, 2);
    if (cmpr == 0)     return 1;  /* 2 is prime */
    if (cmpr < 0)      return 0;  /* below 2 is composite */
    if (mpz_even_p(n)) return 0;  /* multiple of 2 is composite */
  }
  if (mpz_perfect_square_p(n)) return 0;

  mpz_init(t);
  do {
    c += 2;
    mpz_set_ui(t, c);
    k = mpz_jacobi(t, n);
  } while (k == 1);
  if (k == 0) {
    mpz_clear(t);
    return 0;
  }

  mpz_init_set_ui(ra, 1);   mpz_init_set_ui(rb, 1);
  mpz_init_set_ui(a, 1);    mpz_init_set_ui(b, 1);
  mpz_init(ta);   mpz_init(tb);
  mpz_init(n_minus_1);
  mpz_sub_ui(n_minus_1, n, 1);

  len = mpz_sizeinbase(n_minus_1, 2);
  for (bit = 0; bit < len; bit++) {
    if ( mpz_tstbit(n_minus_1, bit) ) {
      mpz_mul(ta, ra, a);
      mpz_mul(tb, rb, b);
      mpz_add(t, ra, rb);
      mpz_add(rb, a, b);
      mpz_mul(rb, rb, t);
      mpz_sub(rb, rb, ta);
      mpz_sub(rb, rb, tb);
      mpz_mod(rb, rb, n);
      mpz_mul_ui(tb, tb, c);
      mpz_add(ra, ta, tb);
      mpz_mod(ra, ra, n);
    }
    if (bit < len-1) {
      mpz_mul(t, b, b);
      mpz_mul_ui(t, t, c);
      mpz_mul(b, b, a);
      mpz_add(b, b, b);
      mpz_mod(b, b, n);
      mpz_mul(a, a, a);
      mpz_add(a, a, t);
      mpz_mod(a, a, n);
    }
  }
  if ( (mpz_cmp_ui(ra,1) == 0) && (mpz_cmp(rb, n_minus_1) == 0) )
    rval = 1;

  mpz_clear(n_minus_1);
  mpz_clear(ta); mpz_clear(tb);
  mpz_clear(a);  mpz_clear(b);
  mpz_clear(ra); mpz_clear(rb);
  mpz_clear(t);
  return rval;
}


UV _GMP_trial_factor(mpz_t n, UV from_n, UV to_n)
{
  size_t log2n = mpz_sizeinbase(n, 2);
  UV p = 0;
  PRIME_ITERATOR(iter);

  if (mpz_cmp_ui(n, 6) < 0) {
    unsigned long un = mpz_get_ui(n);
    if (un == 1) p = 1;
    else if (un == 4 && from_n <= 2 && to_n >= 2) p = 2;
    prime_iterator_destroy(&iter);
    return p;
  }
  if      (from_n <= 2 && to_n >= 2 && mpz_even_p(n))             p = 2;
  else if (from_n <= 3 && to_n >= 3 && mpz_divisible_ui_p(n, 3))  p = 3;
  else if (from_n <= 5 && to_n >= 5 && mpz_divisible_ui_p(n, 5))  p = 5;
  if (p != 0) {
    prime_iterator_destroy(&iter);
    return p;
  }

  if (from_n < 7)
    from_n = 7;
  if (from_n > to_n)
    { prime_iterator_destroy(&iter);  return 0; }
  /* p will be the next prime >= from_n */
  prime_iterator_setprime(&iter, from_n-1);
  p = prime_iterator_next(&iter);

  /* All native math if n fits in an unsigned long */
  if (log2n <= sizeof(unsigned long)*8) {
    unsigned long un = mpz_get_ui(n);
    unsigned long sqrtn = (unsigned long) sqrt((double)un);
    /* Be extra careful here, as we are using unsigned long, which may not
     * match a UV.  But GMP's ui is 'unsigned long' so that's what we have
     * to deal with.  We want to make sure we get the correct integer sqrt,
     * but also watch out for overflow. */
    while (sqrtn*sqrtn > un) sqrtn--;
    while ( (sqrtn+1)*(sqrtn+1) <= un
            && sqrtn < (1UL << 4*sizeof(unsigned long)) )
      sqrtn++;
    if (to_n > sqrtn)
      to_n = sqrtn;
    while (p <= to_n) {
      if ((un % p) == 0)
        break;
      p = prime_iterator_next(&iter);
    }
    prime_iterator_destroy(&iter);
    return (p <= to_n) ? p : 0;
  }

  /* For "small" numbers, this simple method is best. */
  {
    UV small_to = (log2n < 3000)  ?  to_n  :  30000;
    while (p <= small_to) {
      if (mpz_divisible_ui_p(n, p))
        break;
      p = prime_iterator_next(&iter);
    }
    if (p <= small_to || p > to_n) {
      prime_iterator_destroy(&iter);
      return (p <= small_to) ? p : 0;
    }
  }

  /* Simple treesieve.
   * This is much faster than simple divisibility for really big numbers.
   * Credit to Jens K Andersen for writing up the generic algorithm.
   *
   * This will search until the first group element is > to_n, which means
   * we will search a bit farther than to_n.
   */
  {
    UV found = 0;
    unsigned long* xn;        /* leaves */
    mpz_t* xtree[16+1];       /* the tree (maxdepth = 16) */
    mpz_t* xtemp;
    unsigned int i, j, d, depth, leafsize, nleaves;

    /* Decide on the tree depth (3-16) and number of leaves (10-31) */
    {
      unsigned int dp = log2n >> 10;
      depth = 0;
      while (dp >>= 1) depth++;
      if (depth < 3) depth = 3;
      if (depth > 16) depth = 16;
    }
    leafsize = log2n / (1U << depth) / 68;
    nleaves = 1 << depth;
    /* printf("log2n %lu  depth %u  leafsize %u  nleaves %u\n",log2n,depth,leafsize,nleaves); */

    New(0, xn, nleaves * leafsize, unsigned long);
    for (d = 0; d <= depth; d++) {
      unsigned int nodes = 1 << (depth - d);
      New(0, xtree[d], nodes, mpz_t);
      for (j = 0; j < nodes; j++)
        mpz_init(xtree[d][j]);
    }
    xtemp = xtree[1];   /* implies mindepth = 3 */

    while (!found && p <= to_n) {
      /* Create nleaves x[0] values, each the product of leafsize primes */
      for (i = 0; i < nleaves; i++) {
        for (j = 0; j < 4; j++)                  /* Create 4 sub-products */
          mpz_set_ui(xtemp[j], 1);
        for (j = 0; j < leafsize; j++) {
          xn[i*leafsize+j] = p;
          mpz_mul_ui(xtemp[j&3], xtemp[j&3], p);
          p = prime_iterator_next(&iter);
        }
        mpz_mul(xtemp[0], xtemp[0], xtemp[1]);   /* Combine for final product*/
        mpz_mul(xtemp[2], xtemp[2], xtemp[3]);
        mpz_mul(xtree[0][i], xtemp[0], xtemp[2]);
      }
      /* Multiply product tree, xtree[depth][0] has nleaves*leafsize product */
      for (d = 1; d <= depth; d++)
        for (i = 0; i < (1U << (depth-d)); i++)
          mpz_mul(xtree[d][i], xtree[d-1][2*i], xtree[d-1][2*i+1]);
      /* Go backwards replacing the products with remainders */
      mpz_tdiv_r(xtree[depth][0], n, xtree[depth][0]);
      for (d = 1; d <= depth; d++)
        for (i = 0; i < (1U << d); i++)
          mpz_tdiv_r(xtree[depth-d][i], xtree[depth-d+1][i>>1], xtree[depth-d][i]);
      /* Search each leaf for divisors */
      for (i = 0; !found && i < nleaves; i++)
        for (j = 0; j < leafsize; j++)
          if (mpz_divisible_ui_p(xtree[0][i], xn[i*leafsize+j]))
            { found = xn[i*leafsize+j]; break; }
    }
    p = found;
    for (d = 0; d <= depth; d++) {
      unsigned int nodes = 1U << (depth - d);
      for (j = 0; j < nodes; j++)
        mpz_clear(xtree[d][j]);
      Safefree(xtree[d]);
    }
    Safefree(xn);
    if (p > 0 && !mpz_divisible_ui_p(n, p))
      croak("incorrect trial factor\n");
  }
  prime_iterator_destroy(&iter);
  return p;
}

/*
 * is_prob_prime      BPSW -- fast, no known counterexamples
 * is_prime           is_prob_prime + a little extra
 * is_provable_prime  really prove it, which could take a very long time
 *
 * They're all identical for numbers <= 2^64.
 *
 * The extra M-R tests in is_prime start actually costing something after
 * 1000 bits or so.  Primality proving will get quite a bit slower as the
 * number of bits increases.
 *
 * What are these extra M-R tests getting us?  The primary reference is
 * Damgård, Landrock, and Pomerance, 1993.  From Rabin-Monier, with random
 * bases we have p <= 4^-t.  So one extra test gives us p = 0.25, and four
 * tests gives us p = 0.00390625.  But for larger k (bits in n) this is very
 * conservative.  Since the value has passed BPSW, we are only interested in
 * k > 64.  See the calculate-mr-probs.pl script in the xt/ directory.
 * For a 256-bit input, we get:
 *   1 test:  p < 0.000244140625
 *   2 tests: p < 0.00000000441533
 *   3 tests: p < 0.0000000000062550875
 *   4 tests: p < 0.000000000000028421709
 *
 * It's even more extreme as k goes higher.  Also recall that this is the
 * probability once we've somehow found a BPSW pseudoprime.
 */

int primality_pretest(mpz_t n)
{
  /* If less than 1009, make trial factor handle it. */
  if (mpz_cmp_ui(n, BGCD_NEXTPRIME) < 0)
    return _GMP_trial_factor(n, 2, BGCD_LASTPRIME) ? 0 : 2;

  /* Check for tiny divisors */
  if (mpz_even_p(n)) return 0;
  if (sizeof(unsigned long) < 8) {
    if (mpz_gcd_ui(NULL, n, 3234846615UL) != 1) return 0;           /*  3-29 */
  } else {
    if (mpz_gcd_ui(NULL, n, 4127218095UL*3948078067UL)!=1) return 0;/*  3-53 */
    if (mpz_gcd_ui(NULL, n, 4269855901UL*1673450759UL)!=1) return 0;/* 59-101 */
  }

  {
    UV log2n = mpz_sizeinbase(n,2);
    mpz_t t;
    mpz_init(t);

    /* Do a GCD with all primes < 1009 */
    mpz_gcd(t, n, _bgcd);
    if (mpz_cmp_ui(t, 1))
      { mpz_clear(t); return 0; }

    /* No divisors under 1009 */
    if (mpz_cmp_ui(n, BGCD_NEXTPRIME*BGCD_NEXTPRIME) < 0)
      { mpz_clear(t); return 2; }

    /* If we're reasonably large, do a gcd with more primes */
    if (log2n > 700) {
      if (mpz_sgn(_bgcd3) == 0) {
        _GMP_pn_primorial(_bgcd3, BGCD3_PRIMES);
        mpz_divexact(_bgcd3, _bgcd3, _bgcd);
      }
      mpz_gcd(t, n, _bgcd3);
      if (mpz_cmp_ui(t, 1))
        { mpz_clear(t); return 0; }
    } else if (log2n > 300) {
      if (mpz_sgn(_bgcd2) == 0) {
        _GMP_pn_primorial(_bgcd2, BGCD2_PRIMES);
        mpz_divexact(_bgcd2, _bgcd2, _bgcd);
      }
      mpz_gcd(t, n, _bgcd2);
      if (mpz_cmp_ui(t, 1))
        { mpz_clear(t); return 0; }
    }
    mpz_clear(t);
    /* Do more trial division if we think we should.
     * According to Menezes (section 4.45) as well as Park (ISPEC 2005),
     * we want to select a trial limit B such that B = E/D where E is the
     * time for our primality test (one M-R test) and D is the time for
     * one trial division.  Example times on my machine came out to
     *   log2n = 840375, E= 6514005000 uS, D=1.45 uS, E/D = 0.006 * log2n
     *   log2n = 465618, E= 1815000000 uS, D=1.05 uS, E/D = 0.008 * log2n
     *   log2n = 199353, E=  287282000 uS, D=0.70 uS, E/D = 0.01  * log2n
     *   log2n =  99678, E=   56956000 uS, D=0.55 uS, E/D = 0.01  * log2n
     *   log2n =  33412, E=    4289000 uS, D=0.30 uS, E/D = 0.013 * log2n
     *   log2n =  13484, E=     470000 uS, D=0.21 uS, E/D = 0.012 * log2n
     * Our trial division could also be further improved for large inputs.
     */
    if (log2n > 16000) {
      double dB = (double)log2n * (double)log2n * 0.005;
      if (BITS_PER_WORD == 32 && dB > 4200000000.0) dB = 4200000000.0;
      if (_GMP_trial_factor(n, BGCD3_NEXTPRIME, (UV)dB))  return 0;
    } else if (log2n > 4000) {
      if (_GMP_trial_factor(n, BGCD3_NEXTPRIME, 80*log2n))  return 0;
    } else if (log2n > 1600) {
      if (_GMP_trial_factor(n, BGCD3_NEXTPRIME, 30*log2n))  return 0;
    }
  }
  return 1;
}

int _GMP_BPSW(mpz_t n)
{
  if (mpz_cmp_ui(n, 4) < 0)
    return (mpz_cmp_ui(n, 1) <= 0) ? 0 : 1;

  if (_GMP_miller_rabin_ui(n, 2) == 0)   /* Miller Rabin with base 2 */
    return 0;

  if (_GMP_is_lucas_pseudoprime(n, 2 /*extra strong*/) == 0)
    return 0;

  if (mpz_sizeinbase(n, 2) <= 64)        /* BPSW is deterministic below 2^64 */
    return 2;

  return 1;
}

/* Assume n is a BPSW PRP, return 1 (no result), 0 (composite), 2 (prime) */
static int is_deterministic_miller_rabin_prime(mpz_t n)
{
  mpz_t t;
  int i, res = 1, maxp = 0;

  if (mpz_sizeinbase(n, 2) <= 82) {
    mpz_init(t);
    /* n < 3825123056546413051  =>  maxp=9, but BPSW should have handled */
    if      (mpz_set_str(t, "318665857834031151167461",10), mpz_cmp(n,t) < 0)
      maxp = 12;
    else if (mpz_set_str(t,"3317044064679887385961981",10), mpz_cmp(n,t) < 0)
      maxp = 13;
    if (maxp > 0) {
      for (i = 1; i < maxp && res; i++) {
        mpz_set_ui(t, sprimes[i]);
        res = _GMP_miller_rabin(n, t);
      }
      if (res == 1) res = 2;
    }
    mpz_clear(t);
  }
  return res;
}


int _GMP_is_prob_prime(mpz_t n)
{
  /*  Step 1: Look for small divisors.  This is done purely for performance.
   *          It is *not* a requirement for the BPSW test. */
  int res = primality_pretest(n);
  if (res != 1)  return res;

  /* We'd like to do the LLR test here, but it screws with certificates. */

  /*  Step 2: The BPSW test.  spsp base 2 and slpsp. */
  return _GMP_BPSW(n);
}

int _GMP_is_prime(mpz_t n)
{
  UV nbits;
  /* Similar to is_prob_prime, but put LLR before BPSW, then do more testing */

  /* First, simple pretesting */
  int prob_prime = primality_pretest(n);
  if (prob_prime != 1)  return prob_prime;

  /* If the number is of form N=k*2^n-1 and we have a fast proof, do it. */
  prob_prime = llr(n);
  if (prob_prime == 0 || prob_prime == 2) return prob_prime;

  /* If the number is of form N=k*2^n+1 and we have a fast proof, do it. */
  prob_prime = proth(n);
  if (prob_prime == 0 || prob_prime == 2) return prob_prime;

  /* Start with BPSW */
  prob_prime = _GMP_BPSW(n);
  nbits = mpz_sizeinbase(n, 2);

  /* Use Sorenson/Webster 2015 deterministic M-R if possible */
  if (prob_prime == 1) {
    prob_prime = is_deterministic_miller_rabin_prime(n);
    if (prob_prime == 0) gmp_printf("\n\n**** BPSW counter-example found?  ****\n**** N = %Zd ****\n\n", n);
  }

  /* n has passed the ES BPSW test, making it quite unlikely it is a
   * composite (and it cannot be if n < 2^64). */

  /* For small numbers, try a quick BLS75 n-1 proof. */
  if (prob_prime == 1) {
    if (is_proth_form(n))
      prob_prime = _GMP_primality_bls_nm1(n, 2 /* effort */, 0 /* cert */);
    else if (nbits <= 200)
      prob_prime = _GMP_primality_bls_nm1(n, 1 /* effort */, 0 /* cert */);
  }

  /* If prob_prime is still 1, let's run some extra tests.
   * Argument against:  Nobody has yet found a BPSW counterexample, so
   *                    this is wasted time.  They can make a call to one of
   *                    the other tests themselves if they care.
   * Argument for:      is_prime() should be as correct as reasonably possible.
   *                    They can call is_prob_prime() to skip extra tests.
   *
   * Choices include:
   *
   *   - A number of fixed-base Miller-Rabin tests.
   *   - A number of random-base Miller-Rabin tests.
   *   - A Frobenius-Underwood test.
   *   - A Frobenius-Khashin test.
   *
   * The Miller-Rabin tests have better performance.
   *
   * We will use random bases to make it more difficult to get a false
   * result even if someone has a way to generate BPSW pseudoprimes.
   * We use enough bases to keep the probability below 1/595,000.
   */

  if (prob_prime == 1) {
    UV ntests;
    if      (nbits <  80) ntests = 5;  /* p < .00000168 */
    else if (nbits < 105) ntests = 4;  /* p < .00000156 */
    else if (nbits < 160) ntests = 3;  /* p < .00000164 */
    else if (nbits < 413) ntests = 2;  /* p < .00000156 */
    else                  ntests = 1;  /* p < .00000159 */
    prob_prime = _GMP_miller_rabin_random(n, ntests, 0);
    /* prob_prime = _GMP_is_frobenius_underwood_pseudoprime(n); */
    /* prob_prime = _GMP_is_frobenius_khashin_pseudoprime(n); */
  }

  /* Using Damgård, Landrock, and Pomerance, we get upper bounds:
   * k <=  64      p = 0
   * k <   80      p < 7.53e-08
   * k <  105      p < 3.97e-08
   * k <  160      p < 1.68e-08
   * k >= 160      p < 2.57e-09
   * This (1) pretends the SPSP-2 is included in the random tests, (2) doesn't
   * take into account the trial division, (3) ignores the observed
   * SPSP-2 / Strong Lucas anti-correlation that makes the BPSW test so
   * useful.  Hence these numbers are still extremely conservative.
   *
   * For an idea of how conservative we are being, we have now exceeded the
   * tests used in Mathematica, Maple, Pari, and SAGE.  Note: Pari pre-2.3
   * used just M-R tests.  Pari 2.3+ uses BPSW with no extra M-R checks for
   * is_pseudoprime, but isprime uses APRCL which, being a proof, could only
   * output a pseudoprime through a coding error.
   */

  return prob_prime;
}


int _GMP_is_provable_prime(mpz_t n, char** prooftext)
{
  int prob_prime = primality_pretest(n);
  if (prob_prime != 1)  return prob_prime;

  /* Try LLR and Proth if they don't need a proof certificate. */
  if (prooftext == 0) {
    prob_prime = llr(n);
    if (prob_prime == 0 || prob_prime == 2) return prob_prime;
    prob_prime = proth(n);
    if (prob_prime == 0 || prob_prime == 2) return prob_prime;
  }

  /* Start with BPSW */
  prob_prime = _GMP_BPSW(n);
  if (prob_prime != 1)  return prob_prime;

  /* Use Sorenson/Webster 2015 deterministic M-R if possible */
  if (prooftext == 0) {
    prob_prime = is_deterministic_miller_rabin_prime(n);
    if (prob_prime != 1)  return prob_prime;
  }

  /* Run one more M-R test, just in case. */
  prob_prime = _GMP_miller_rabin_random(n, 1, 0);
  if (prob_prime != 1)  return prob_prime;

  /* We can choose a primality proving algorithm:
   *   AKS    _GMP_is_aks_prime       really slow, don't bother
   *   N-1    _GMP_primality_bls_nm1  small or special numbers
   *   ECPP   _GMP_ecpp               fastest in general
   */

  /* Give n-1 a small go */
  prob_prime = _GMP_primality_bls_nm1(n, is_proth_form(n) ? 3 : 1, prooftext);
  if (prob_prime != 1)  return prob_prime;

  /* ECPP */
  prob_prime = _GMP_ecpp(n, prooftext);

  return prob_prime;
}


/*****************************************************************************/
/*          AKS.    This implementation is quite slow, but useful to have.   */

static int test_anr(UV a, mpz_t n, UV r, mpz_t* px, mpz_t* py)
{
  int retval = 1;
  UV i, n_mod_r;
  mpz_t t;

  for (i = 0; i < r; i++)
    mpz_set_ui(px[i], 0);

  a %= r;
  mpz_set_ui(px[0], a);
  mpz_set_ui(px[1], 1);

  poly_mod_pow(py, px, n, r, n);

  mpz_init(t);
  n_mod_r = mpz_mod_ui(t, n, r);
  if (n_mod_r >= r)  croak("n mod r >= r ?!");
  mpz_sub_ui(t, py[n_mod_r], 1);
  mpz_mod(py[n_mod_r], t, n);
  mpz_sub_ui(t, py[0], a);
  mpz_mod(py[0], t, n);
  mpz_clear(t);

  for (i = 0; i < r; i++)
    if (mpz_sgn(py[i]))
      retval = 0;
  return retval;
}

/* TODO: is_primitive_root, mpz_logn, mpz_log2 should move to utilities */

#if AKS_VARIANT != AKS_VARIANT_V6
static int is_primitive_root(mpz_t n, UV r)
{
  mpz_t m, modr;
  UV p, rm1 = r-1;
  UV lim = (UV) (sqrt(rm1) + 0.00001);
  UV ret = 0;

  mpz_init(m);
  mpz_init_set_ui(modr, r);
  if ((rm1 % 2) == 0) {
    mpz_powm_ui(m, n, (rm1/2), modr);
    if (!mpz_cmp_ui(m,1))
      goto END_PRIMROOT;
  }
  for (p = 3; p <= lim; p += 2) {
    if ((rm1 % p) == 0) {
      mpz_powm_ui(m, n, (rm1/p), modr);
      if (!mpz_cmp_ui(m,1))
        goto END_PRIMROOT;
    }
  }
  ret = 1;
END_PRIMROOT:
  mpz_clear(m);
  mpz_clear(modr);
  return ret;
}
#endif

static double mpz_logn(mpz_t n)
{
  long exp;
  double logn = mpz_get_d_2exp(&exp, n);
  logn = log(logn) + (log(2) * exp);
  return logn;
}

static double mpz_log2(mpz_t n) {
  long exp;
  double logn = mpz_get_d_2exp(&exp, n);
  logn = exp + log(logn)/log(2);
  return logn;
}

#if AKS_VARIANT == AKS_VARIANT_BERN21
static UV largest_factor(UV n) {
  UV p = 2;
  PRIME_ITERATOR(iter);
  while (n >= p*p && !prime_iterator_isprime(&iter, n)) {
    while ( (n % p) == 0  &&  n >= p*p ) { n /= p; }
    p = prime_iterator_next(&iter);
  }
  prime_iterator_destroy(&iter);
  return n;
}
#endif


int _GMP_is_aks_prime(mpz_t n)
{
  mpz_t *px, *py;
  int retval;
  UV i, s, r, a;
  UV starta = 1;
  int _verbose = get_verbose_level();

  if (mpz_cmp_ui(n, 4) < 0)
    return (mpz_cmp_ui(n, 1) <= 0) ? 0 : 1;

  /* Just for performance: check small divisors: 2*3*5*7*11*13*17*19*23 */
  if (mpz_gcd_ui(0, n, 223092870UL) != 1 && mpz_cmp_ui(n, 23) > 0)
    return 0;

  if (mpz_perfect_power_p(n))
    return 0;

#if AKS_VARIANT == AKS_VARIANT_V6    /* From the V6 AKS paper */
  {
    mpz_t sqrtn, t;
    double log2n;
    UV limit, startr;
    PRIME_ITERATOR(iter);

    mpz_init(sqrtn);
    mpz_sqrt(sqrtn, n);

    log2n = mpz_log2(n);
    limit = (UV) floor( log2n * log2n );

    if (_verbose>1) gmp_printf("# AKS checking order_r(%Zd) to %"UVuf"\n", n, (unsigned long) limit);

    /* Using a native r limits us to ~2000 digits in the worst case (r ~ log^5n)
     * but would typically work for 100,000+ digits (r ~ log^3n).  This code is
     * far too slow to matter either way.  Composite r is ok here, but it will
     * always end up prime, so save time and just check primes. */
    retval = 0;
    /* Start order search at a good spot.  Idea from Nemana and Venkaiah. */
    startr = (mpz_sizeinbase(n,2)-1) * (mpz_sizeinbase(n,2)-1);
    startr = (startr < 1002) ? 2 : startr - 100;
    for (r = 2; /* */; r = prime_iterator_next(&iter)) {
      if (mpz_divisible_ui_p(n, r) ) /* r divides n.  composite. */
        { retval = 0; break; }
      if (mpz_cmp_ui(sqrtn, r) <= 0) /* no r <= sqrtn divides n.  prime. */
        { retval = 1; break; }
      if (r < startr) continue;
      if (mpz_order_ui(r, n, limit) > limit)
        { retval = 2; break; }
    }
    prime_iterator_destroy(&iter);
    mpz_clear(sqrtn);
    if (retval != 2) return retval;

    /* Since r is prime, phi(r) = r-1. */
    s = (UV) floor( sqrt(r-1) * log2n );
  }
#elif AKS_VARIANT == AKS_VARIANT_BORNEMANN /* Bernstein + Voloch */
  {
    UV slim;
    double c2, x;
    /* small t = few iters of big poly.  big t = many iters of small poly */
    double const t = (mpz_sizeinbase(n, 2) <= 64) ? 32 : 40;
    double const t1 = (1.0/((t+1)*log(t+1)-t*log(t)));
    double const dlogn = mpz_logn(n);
    mpz_t tmp;
    PRIME_ITERATOR(iter);

    mpz_init(tmp);
    prime_iterator_setprime(&iter, (UV) (t1*t1 * dlogn*dlogn) );
    r = prime_iterator_next(&iter);
    while (!is_primitive_root(n,r))
      r = prime_iterator_next(&iter);
    prime_iterator_destroy(&iter);

    slim = (UV) (2*t*(r-1));
    c2 = dlogn * floor(sqrt(r));
    { /* Binary search for first s in [1,slim] where x >= 0 */
      UV bi = 1;
      UV bj = slim;
      while (bi < bj) {
        s = bi + (bj-bi)/2;
        mpz_bin_uiui(tmp, r+s-1, s);
        x = mpz_logn(tmp) / c2 - 1.0;
        if (x < 0)  bi = s+1;
        else        bj = s;
      }
      s = bi-1;
    }
    s = (s+3) >> 1;
    /* Bornemann checks factors up to (s-1)^2, we check to max(r,s) */
    /* slim = (s-1)*(s-1); */
    slim = (r > s) ? r : s;
    if (_verbose > 1) printf("# aks trial to %"UVuf"\n", slim);
    if (_GMP_trial_factor(n, 2, slim) > 1)
      { mpz_clear(tmp); return 0; }
    mpz_sqrt(tmp, n);
    if (mpz_cmp_ui(tmp, slim) <= 0)
      { mpz_clear(tmp); return 1; }
    mpz_clear(tmp);
  }
#elif AKS_VARIANT == AKS_VARIANT_BERN21
  { /* Bernstein 2003, theorem 2.1 (simplified) */
    UV q;
    double slim, scmp, x;
    mpz_t t, t2;
    PRIME_ITERATOR(iter);
    mpz_init(t);  mpz_init(t2);
    r = s = 0;
    while (1) {
      /* todo: Check r|n and r >= sqrt(n) here instead of waiting */
      if (mpz_cmp_ui(n, r) <= 0) break;
      r = prime_iterator_next(&iter);
      q = largest_factor(r-1);
      mpz_set_ui(t, r);
      mpz_powm_ui(t, n, (r-1)/q, t);
      if (mpz_cmp_ui(t, 1) <= 0) continue;
      scmp = 2 * floor(sqrt(r)) * mpz_log2(n);

      slim = 20 * (r-1);

      /* Check viability */
      mpz_bin_uiui(t, q+slim-1, slim); if (mpz_log2(t) < scmp) continue;

      for (s = 2; s < slim; s++) {
        mpz_bin_uiui(t, q+s-1, s);
        if (mpz_log2(t) > scmp) break;
      }
      if (s < slim) break;
    }
    mpz_clear(t);  mpz_clear(t2);
    prime_iterator_destroy(&iter);
    if (_GMP_trial_factor(n, 2, s) > 1)
      return 0;
  }
#elif AKS_VARIANT == AKS_VARIANT_BERN22
  { /* Bernstein 2003, theorem 2.2 (simplified) */
    UV q;
    double slim, scmp, x;
    mpz_t t, t2;
    PRIME_ITERATOR(iter);
    mpz_init(t);  mpz_init(t2);
    r = s = 0;
    while (1) {
      /* todo: Check r|n and r >= sqrt(n) here instead of waiting */
      if (mpz_cmp_ui(n, r) <= 0) break;
      r = prime_iterator_next(&iter);
      if (!is_primitive_root(n,r)) continue;
      q = r-1;   /* Since r is prime, phi(r) = r-1 */
      scmp = 2 * floor(sqrt(r-1)) * mpz_log2(n);

      slim = 20 * (r-1);

      /* Check viability */
      mpz_bin_uiui(t, q+slim-1, slim); if (mpz_log2(t) < scmp) continue;

      for (s = 2; s < slim; s++) {
        mpz_bin_uiui(t, q+s-1, s);
        if (mpz_log2(t) > scmp) break;
      }
      if (s < slim) break;
    }
    mpz_clear(t);  mpz_clear(t2);
    prime_iterator_destroy(&iter);
    if (_GMP_trial_factor(n, 2, s) > 1)
      return 0;
  }
#elif AKS_VARIANT == AKS_VARIANT_BERN23
  { /* Bernstein 2003, theorem 2.3 (simplified) */
    UV q, d, limit;
    double slim, scmp, sbin, x, log2n;
    mpz_t t, t2;
    PRIME_ITERATOR(iter);
    mpz_init(t);  mpz_init(t2);
    log2n = mpz_log2(n);
    limit = (UV) floor( log2n * log2n );
    r = 2;
    s = 0;
    while (1) {
      /* todo: Check r|n and r >= sqrt(n) here instead of waiting */
      if (mpz_cmp_ui(n, r) <= 0) break;
      r++;
      UV gcd = mpz_gcd_ui(NULL, n, r);
      if (gcd != 1) { mpz_clear(t); mpz_clear(t2); return 0; }
      UV v = mpz_order_ui(r, n, limit);
      if (v >= limit) continue;

      mpz_set_ui(t2, r);
      totient(t, t2);
      q = mpz_get_ui(t);
      UV phiv = q/v;
      /* printf("phi(%lu)/v = %lu/%lu = %lu\n", r, q, v, phiv); */

      /* This is extremely inefficient. */

      /* Choose an s value we'd be happy with */
      slim = 20 * (r-1);

      /* Quick check to see if it could work with s=slim, d=1 */
      mpz_bin_uiui(t, q+slim-1, slim);
      sbin = mpz_log2(t);
      if (sbin < 2*floor(sqrt(q))*log2n)
        continue;

      for (s = 2; s < slim; s++) {
        mpz_bin_uiui(t, q+s-1, s);
        sbin = mpz_log2(t);
        if (sbin < 2*floor(sqrt(q))*log2n) continue;   /* d=1 */
        /* Check each number dividing phi(r)/v */
        for (d = 2; d < phiv; d++) {
          if ((phiv % d) != 0) continue;
          scmp = 2 * d * floor(sqrt(q/d)) * log2n;
          if (sbin < scmp) break;
        }
        /* if we did not exit early, this s worked for each d.  This s wins. */
        if (d >= phiv) break;
      }
      if (s < slim) break;
    }
    mpz_clear(t);  mpz_clear(t2);
    prime_iterator_destroy(&iter);
    if (_GMP_trial_factor(n, 2, s) > 1)
      return 0;
  }
#elif AKS_VARIANT == AKS_VARIANT_BERN41
  {
    double const log2n = mpz_log2(n);
    double const r0 = 0.010 * log2n * log2n;
    double scmp;
    UV d, j, slim;
    mpz_t tmp, tmp2;
    PRIME_ITERATOR(iter);

    mpz_init(tmp);  mpz_init(tmp2);
    prime_iterator_setprime(&iter, (UV) r0);
    r = prime_iterator_next(&iter);
    while (!is_primitive_root(n,r) || r < 3)
      r = prime_iterator_next(&iter);
    prime_iterator_destroy(&iter);

    d = (UV) (0.5 * (r-1) + 0.5);
    i = (UV) (0.475 * (r-1) + 0.5);
    j = i;

    /*gmp_printf("for n=%Zd chose r=%lu d=%lu i=%lu j=%lu\n", n, r, d, i, j);*/
    /* Ensure conditions are correct */
    if (d > r-2)     d = r-2;
    if (i > d)       i = d;
    if (j > (r-2-d)) j = r-2-d;

    /* Find a suitable s */

    slim = (UV) (8*(r-1));   /* This could be tightened up.  Depends on r0! */
    scmp = ceil(sqrt( (r-1)/3.0 )) * mpz_log2(n);
    { /* Binary search for first s in [1,slim] where x >= 0 */
      UV bi = 1;
      UV bj = slim;
      while (bi < bj) {
        s = bi + (bj-bi)/2;
        mpz_bin_uiui(tmp2, 2*s, i);
        mpz_bin_uiui(tmp, d, i);       mpz_mul(tmp2, tmp2, tmp);
        mpz_bin_uiui(tmp, 2*s-i, j);   mpz_mul(tmp2, tmp2, tmp);
        mpz_bin_uiui(tmp, r-2-d, j);   mpz_mul(tmp2, tmp2, tmp);
        double left = mpz_log2(tmp2);
        if (left < scmp)  bi = s+1;
        else              bj = s;
      }
      /* Our S goes from 2 to s+1. */
      starta = 2;
      s = s+1;
    }
    /* printf("chose r=%lu s=%lu d = %lu i = %lu j = %lu\n", r, s, d, i, j); */

    /* Check divisibility to (s-1)^2 to cover both gcd conditions */
    slim = (s-1)*(s-1);
    if (_verbose > 1) printf("# aks trial to %"UVuf"\n", slim);
    if (_GMP_trial_factor(n, 2, slim) > 1)
      { mpz_clear(tmp); mpz_clear(tmp2); return 0; }
    mpz_sqrt(tmp, n);
    if (mpz_cmp_ui(tmp, slim) <= 0)
      { mpz_clear(tmp); mpz_clear(tmp2); return 1; }

    /* Check b^(n-1) = 1 mod n for b in [2..s] */
    mpz_sub_ui(tmp2, n, 1);
    for (i = 2; i <= s; i++) {
      mpz_set_ui(tmp, i);
      mpz_powm(tmp, tmp, tmp2, n);
      if (mpz_cmp_ui(tmp, 1) != 0)
        { mpz_clear(tmp); mpz_clear(tmp2); return 0; }
    }

    mpz_clear(tmp);  mpz_clear(tmp2);
  }

#endif

  if (_verbose) gmp_printf("# AKS %Zd.  r = %"UVuf" s = %"UVuf"\n", n, (unsigned long) r, (unsigned long) s);

  /* Create the three polynomials we will use */
  New(0, px, r, mpz_t);
  New(0, py, r, mpz_t);
  if ( !px || !py )
    croak("allocation failure\n");
  for (i = 0; i < r; i++) {
    mpz_init(px[i]);
    mpz_init(py[i]);
  }

  retval = 1;
  for (a = starta; a <= s; a++) {
    if (! test_anr(a, n, r, px, py) ) {
      retval = 0;
      break;
    }
    if (_verbose>1) { printf("."); fflush(stdout); }
  }
  if (_verbose>1) { printf("\n"); fflush(stdout); };

  /* Free the polynomials */
  for (i = 0; i < r; i++) {
    mpz_clear(px[i]);
    mpz_clear(py[i]);
  }
  Safefree(px);
  Safefree(py);

  return retval;
}

/*****************************************************************************/

/* Controls how many numbers to sieve.  Little time impact. */
#define NPS_MERIT  30.0
/* Controls how many primes to use.  Big time impact. */
#define NPS_DEPTH(log2n, log2log2n) \
  (log2n < 100) ? 1000 : \
  (BITS_PER_WORD == 32 && log2n > 9000U) ? UVCONST(2500000000) : \
  (BITS_PER_WORD == 64 && log2n > 4294967294U) ? UVCONST(9300000000000000000) :\
  ((log2n * (log2n >> 5) * log2log2n) >> 1)

static void next_prime_with_sieve(mpz_t n) {
  UV i, log2n, log2log2n, width, depth;
  uint32_t* comp;
  mpz_t t, base;
  log2n = mpz_sizeinbase(n, 2);
  for (log2log2n = 1, i = log2n; i >>= 1; ) log2log2n++;
  width = (UV) (NPS_MERIT/1.4427 * (double)log2n + 0.5);
  depth = NPS_DEPTH(log2n, log2log2n);

  if (width & 1) width++;                     /* Make width even */
  mpz_add_ui(n, n, mpz_even_p(n) ? 1 : 2);    /* Set n to next odd */
  mpz_init(t);  mpz_init(base);
  while (1) {
    mpz_set(base, n);
    comp = partial_sieve(base, width, depth); /* sieve range to depth */
    for (i = 1; i <= width; i += 2) {
      if (!TSTAVAL(comp, i)) {
        mpz_add_ui(t, base, i);               /* We found a candidate */
        if (_GMP_BPSW(t)) {
          mpz_set(n, t);
          mpz_clear(t);  mpz_clear(base);
          Safefree(comp);
          return;
        }
      }
    }
    Safefree(comp);       /* A huge gap found, so sieve another range */
    mpz_add_ui(n, n, width);
  }
}

static void prev_prime_with_sieve(mpz_t n) {
  UV i, j, log2n, log2log2n, width, depth;
  uint32_t* comp;
  mpz_t t, base;
  log2n = mpz_sizeinbase(n, 2);
  for (log2log2n = 1, i = log2n; i >>= 1; ) log2log2n++;
  width = (UV) (NPS_MERIT/1.4427 * (double)log2n + 0.5);
  depth = NPS_DEPTH(log2n, log2log2n);

  mpz_sub_ui(n, n, mpz_even_p(n) ? 1 : 2);       /* Set n to prev odd */
  width = 64 * ((width+63)/64);                /* Round up to next 64 */
  mpz_init(t);  mpz_init(base);
  while (1) {
    mpz_sub_ui(base, n, width-2);
    /* gmp_printf("sieve from %Zd to %Zd width %lu\n", base, n, width); */
    comp = partial_sieve(base, width, depth); /* sieve range to depth */
    /* if (mpz_odd_p(base)) croak("base off after partial");
       if (width & 1) croak("width is odd after partial"); */
    for (j = 1; j < width; j += 2) {
      i = width - j;
      if (!TSTAVAL(comp, i)) {
        mpz_add_ui(t, base, i);               /* We found a candidate */
        if (_GMP_BPSW(t)) {
          mpz_set(n, t);
          mpz_clear(t);  mpz_clear(base);
          Safefree(comp);
          return;
        }
      }
    }
    Safefree(comp);       /* A huge gap found, so sieve another range */
    mpz_sub_ui(n, n, width);
  }
}

/* Modifies argument */
void _GMP_next_prime(mpz_t n)
{
  if (mpz_cmp_ui(n, 29) < 0) { /* small inputs */

    UV m = mpz_get_ui(n);
    m = (m < 2) ? 2 : (m < 3) ? 3 : (m < 5) ? 5 : next_wheel[m];
    mpz_set_ui(n, m);

  } else if (mpz_sizeinbase(n,2) > 120) {

    next_prime_with_sieve(n);

  } else {

    UV m23 = mpz_fdiv_ui(n, 223092870UL);  /* 2*3*5*7*11*13*17*19*23 */
    UV m = m23 % 30;
    do {
      UV skip = wheel_advance[m];
      mpz_add_ui(n, n, skip);
      m23 += skip;
      m = next_wheel[m];
    } while ( !(m23% 7) || !(m23%11) || !(m23%13) || !(m23%17) ||
              !(m23%19) || !(m23%23) || !_GMP_is_prob_prime(n) );

  }
}

/* Modifies argument */
void _GMP_prev_prime(mpz_t n)
{
  UV m, m23;
  if (mpz_cmp_ui(n, 29) <= 0) { /* small inputs */

    m = mpz_get_ui(n);
    m = (m < 3) ? 0 : (m < 4) ? 2 : (m < 6) ? 3 : (m < 8) ? 5 : prev_wheel[m];
    mpz_set_ui(n, m);

  } else if (mpz_sizeinbase(n,2) > 200) {

    prev_prime_with_sieve(n);

  } else {

    m23 = mpz_fdiv_ui(n, 223092870UL);  /* 2*3*5*7*11*13*17*19*23 */
    m = m23 % 30;
    m23 += 223092870UL;  /* No need to re-mod inside the loop */
    do {
      UV skip = wheel_retreat[m];
      mpz_sub_ui(n, n, skip);
      m23 -= skip;
      m = prev_wheel[m];
    } while ( !(m23% 7) || !(m23%11) || !(m23%13) || !(m23%17) ||
              !(m23%19) || !(m23%23) || !_GMP_is_prob_prime(n) );

  }
}

#define LAST_DOUBLE_PROD \
  ((BITS_PER_WORD == 32) ? UVCONST(65521) : UVCONST(4294967291))
void _GMP_pn_primorial(mpz_t prim, UV n)
{
  UV p = 2;
  PRIME_ITERATOR(iter);

  if (n < 800) {  /* Don't go above 6500 to prevent overflow below */
    /* Simple linear multiplication, two at a time */
    mpz_set_ui(prim, 1);
    while (n-- > 0) {
      if (n > 0) { p *= prime_iterator_next(&iter); n--; }
      mpz_mul_ui(prim, prim, p);
      p = prime_iterator_next(&iter);
    }
  } else {
    /* Shallow product tree, ~10x faster for large values */
    mpz_t t[16];
    UV i;
    for (i = 0; i < 16; i++)  mpz_init_set_ui(t[i], 1);
    i = 0;
    while (n-- > 0) {
      if (p <= LAST_DOUBLE_PROD && n > 0)
        { p *= prime_iterator_next(&iter); n--; }
      mpz_mul_ui(t[i&15], t[i&15], p);
      p = prime_iterator_next(&iter);
      i++;
    }
    mpz_product(t, 0, 16-1);
    mpz_set(prim, t[0]);
    for (i = 0; i < 16; i++)  mpz_clear(t[i]);
  }
  prime_iterator_destroy(&iter);
}
void _GMP_primorial(mpz_t prim, UV n)
{
#if (__GNU_MP_VERSION > 5) || (__GNU_MP_VERSION == 5 && __GNU_MP_VERSION_MINOR >= 1)
  mpz_primorial_ui(prim, n);
#else
  UV p = 2;
  PRIME_ITERATOR(iter);

  if (n < 1000) {
    /* Simple linear multiplication, one at a time */
    mpz_set_ui(prim, 1);
    while (p <= n) {
      mpz_mul_ui(prim, prim, p);
      p = prime_iterator_next(&iter);
    }
  } else {
    /* Shallow product tree, ~10x faster for large values */
    mpz_t t[16];
    UV i;
    for (i = 0; i < 16; i++)  mpz_init_set_ui(t[i], 1);
    i = 0;
    while (p <= n) {
      mpz_mul_ui(t[i&15], t[i&15], p);
      p = prime_iterator_next(&iter);
      i++;
    }
    mpz_product(t, 0, 16-1);
    mpz_set(prim, t[0]);
    for (i = 0; i < 16; i++)  mpz_clear(t[i]);
  }
  prime_iterator_destroy(&iter);
#endif
}

/* Luschny's version of the "Brent-Harvey" method */
void bernfrac(mpz_t num, mpz_t den, mpz_t zn)
{
  UV k, j, n = mpz_get_ui(zn);
  mpz_t* T;
  mpz_t t;

  if      (n == 0) { mpz_set_ui(num, 1);  mpz_set_ui(den, 1); return; }
  else if (n == 1) { mpz_set_ui(num, 1);  mpz_set_ui(den, 2); return; }
  else if (n & 1)  { mpz_set_ui(num, 0);  mpz_set_ui(den, 1); return; }
  n >>= 1;
  New(0, T, n+1, mpz_t);
  for (k = 1; k <= n; k++)  mpz_init(T[k]);
  mpz_set_ui(T[1], 1);

  mpz_init(t);

  for (k = 2; k <= n; k++)
    mpz_mul_ui(T[k], T[k-1], k-1);

  for (k = 2; k <= n; k++) {
    for (j = k; j <= n; j++) {
      mpz_mul_ui(t, T[j], j-k+2);
      mpz_mul_ui(T[j], T[j-1], j-k);
      mpz_add(T[j], T[j], t);
    }
  }

  mpz_mul_ui(num, T[n], n);
  mpz_mul_si(num, num, (n & 1) ? 2 : -2);
  mpz_set_ui(t, 1);
  mpz_mul_2exp(den, t, 2*n);  /* den = U = 1 << 2n  */
  mpz_sub_ui(t, den, 1);      /* t = U-1            */
  mpz_mul(den, den, t);       /* den = U*(U-1)      */
  mpz_gcd(t, num, den);
  mpz_divexact(num, num, t);
  mpz_divexact(den, den, t);
  mpz_clear(t);
  for (k = 1; k <= n; k++)  mpz_clear(T[k]);
  Safefree(T);
}
static void _harmonic(mpz_t a, mpz_t b, mpz_t t) {
  mpz_sub(t, b, a);
  if (mpz_cmp_ui(t, 1) == 0) {
    mpz_set(b, a);
    mpz_set_ui(a, 1);
  } else {
    mpz_t q, r;
    mpz_add(t, a, b);
    mpz_tdiv_q_2exp(t, t, 1);
    mpz_init_set(q, t); mpz_init_set(r, t);
    _harmonic(a, q, t);
    _harmonic(r, b, t);
    mpz_mul(a, a, b);
    mpz_mul(t, q, r);
    mpz_add(a, a, t);
    mpz_mul(b, b, q);
    mpz_clear(q); mpz_clear(r);
  }
}

void harmfrac(mpz_t num, mpz_t den, mpz_t zn)
{
  mpz_t t;
  mpz_init(t);
  mpz_add_ui(den, zn, 1);
  mpz_set_ui(num, 1);
  _harmonic(num, den, t);
  mpz_gcd(t, num, den);
  mpz_divexact(num, num, t);
  mpz_divexact(den, den, t);
  mpz_clear(t);
}

char* harmreal(mpz_t zn, unsigned long prec) {
  char* out;
  mpz_t num, den;
  mpf_t fnum, fden, res;
  unsigned long numsize, densize;

  mpz_init(num); mpz_init(den);
  harmfrac(num, den, zn);

  numsize = mpz_sizeinbase(num, 10);
  densize = mpz_sizeinbase(den, 10);

  mpf_init2(fnum, 1 + mpz_sizeinbase(num,2));
  mpf_init2(fden, 1 + mpz_sizeinbase(den,2));
  mpf_set_z(fnum, num);  mpf_set_z(fden, den);
  mpz_clear(den); mpz_clear(num);

  mpf_init2(res, (unsigned long) (8+prec*3.4) );
  mpf_div(res, fnum, fden);
  mpf_clear(fnum);  mpf_clear(fden);

  New(0, out, (10+numsize-densize)+prec, char);
  gmp_sprintf(out, "%.*Ff", (int)(prec), res);
  mpf_clear(res);

  return out;
}

void stirling(mpz_t r, unsigned long n, unsigned long m, UV type)
{
  mpz_t t, t2;
  unsigned long j;
  if (type < 1 || type > 3) croak("stirling type must be 1, 2, or 3");
  if (n == m) {
    mpz_set_ui(r, 1);
  } else if (n == 0 || m == 0 || m > n) {
    mpz_set_ui(r,0);
  } else if (m == 1) {
    switch (type) {
      case 1:  mpz_fac_ui(r, n-1);  if (!(n&1)) mpz_neg(r, r); break;
      case 2:  mpz_set_ui(r, 1); break;
      case 3:
      default: mpz_fac_ui(r, n); break;
    }
  } else {
    mpz_init(t);  mpz_init(t2);
    mpz_set_ui(r,0);
    if (type == 3) {
      mpz_bin_uiui(t, n-1, m-1);
      mpz_fac_ui(t2, n);
      mpz_mul(r, t, t2);
      mpz_fac_ui(t2, m);
      mpz_divexact(r, r, t2);
    } else if (type == 2) {
      for (j = 1; j <= m; j++) {
        mpz_bin_uiui(t, m, j);
        mpz_ui_pow_ui(t2, j, n);
        mpz_mul(t, t, t2);
        if ((m-j) & 1)  mpz_sub(r, r, t);
        else            mpz_add(r, r, t);
      }
      mpz_fac_ui(t, m);
      mpz_divexact(r, r, t);
    } else {
      for (j = 1; j <= n-m; j++) {
        mpz_bin_uiui(t,  n+j-1, n+j-m);
        mpz_bin_uiui(t2, n+n-m, n-j-m);
        mpz_mul(t, t, t2);
        stirling(t2, n+j-m, j, 2);
        mpz_mul(t, t, t2);
        if (j & 1)      mpz_sub(r, r, t);
        else            mpz_add(r, r, t);
      }
    }
    mpz_clear(t2);  mpz_clear(t);
  }
}

/* Goetgheluck method.  Also thanks to Peter Luschny */
void binomial(mpz_t r, UV n, UV k)
{
  UV fi, nk, sqrtn, piN, prime, i, j;
  UV* primes;
  mpz_t* mprimes;

  if (k > n)            { mpz_set_ui(r, 0); return; }
  if (k == 0 || k == n) { mpz_set_ui(r, 1); return; }

  if (k > n/2)  k = n-k;

  sqrtn = (UV) (sqrt((double)n));
  fi = 0;
  nk = n-k;
  primes = sieve_to_n(n, &piN);

#define PUSHP(p) \
 do { \
   if ((j++ % 8) == 0) mpz_init_set_ui(mprimes[fi++], p); \
   else                mpz_mul_ui(mprimes[fi-1], mprimes[fi-1], p); \
 } while (0)

  New(0, mprimes, (piN+7)/8, mpz_t);

  for (i = 0, j = 0; i < piN; i++) {
    prime = primes[i];
    if (prime > nk) {
      PUSHP(prime);
    } else if (prime > n/2) {
      /* nothing */
    } else if (prime > sqrtn) {
      if (n % prime < k % prime)
        PUSHP(prime);
    } else {
      UV N = n, K = k, p = 1, s = 0;
      while (N > 0) {
        s = (N % prime) < (K % prime + s) ? 1 : 0;
        if (s == 1)  p *= prime;
        N /= prime;
        K /= prime;
      }
      if (p > 1)
        PUSHP(p);
    }
  }
  Safefree(primes);
  mpz_product(mprimes, 0, fi-1);
  mpz_set(r, mprimes[0]);
  for (i = 0; i < fi; i++)
    mpz_clear(mprimes[i]);
  Safefree(mprimes);
}



#define TEST_FOR_2357(n, f) \
  { \
    if (mpz_divisible_ui_p(n, 2)) { mpz_set_ui(f, 2); return 1; } \
    if (mpz_divisible_ui_p(n, 3)) { mpz_set_ui(f, 3); return 1; } \
    if (mpz_divisible_ui_p(n, 5)) { mpz_set_ui(f, 5); return 1; } \
    if (mpz_divisible_ui_p(n, 7)) { mpz_set_ui(f, 7); return 1; } \
    if (mpz_cmp_ui(n, 121) < 0) { return 0; } \
  }



int _GMP_prho_factor(mpz_t n, mpz_t f, UV a, UV rounds)
{
  mpz_t U, V, oldU, oldV, m;
  int i;
  const UV inner = 256;

  TEST_FOR_2357(n, f);
  rounds = (rounds + inner - 1) / inner;
  mpz_init_set_ui(U, 7);
  mpz_init_set_ui(V, 7);
  mpz_init(m);
  mpz_init(oldU);
  mpz_init(oldV);
  while (rounds-- > 0) {
    mpz_set_ui(m, 1); mpz_set(oldU, U);  mpz_set(oldV, V);
    for (i = 0; i < (int)inner; i++) {
      mpz_mul(U, U, U);  mpz_add_ui(U, U, a);  mpz_tdiv_r(U, U, n);
      mpz_mul(V, V, V);  mpz_add_ui(V, V, a);  mpz_tdiv_r(V, V, n);
      mpz_mul(V, V, V);  mpz_add_ui(V, V, a);  mpz_tdiv_r(V, V, n);
      if (mpz_cmp(U, V) >= 0)  mpz_sub(f, U, V);
      else                     mpz_sub(f, V, U);
      mpz_mul(m, m, f);
      mpz_tdiv_r(m, m, n);
    }
    mpz_gcd(f, m, n);
    if (!mpz_cmp_ui(f, 1))
      continue;
    if (!mpz_cmp(f, n)) {
      /* f == n, so we have to back up to see what factor got found */
      mpz_set(U, oldU); mpz_set(V, oldV);
      i = inner;
      do {
        mpz_mul(U, U, U);  mpz_add_ui(U, U, a);  mpz_tdiv_r(U, U, n);
        mpz_mul(V, V, V);  mpz_add_ui(V, V, a);  mpz_tdiv_r(V, V, n);
        mpz_mul(V, V, V);  mpz_add_ui(V, V, a);  mpz_tdiv_r(V, V, n);
        if (mpz_cmp(U, V) >= 0)  mpz_sub(f, U, V);
        else                     mpz_sub(f, V, U);
        mpz_gcd(f, f, n);
      } while (!mpz_cmp_ui(f, 1) && i-- != 0);
      if ( (!mpz_cmp_ui(f, 1)) || (!mpz_cmp(f, n)) )  break;
    }
    mpz_clear(U); mpz_clear(V); mpz_clear(m); mpz_clear(oldU); mpz_clear(oldV);
    return 1;
  }
  mpz_clear(U); mpz_clear(V); mpz_clear(m); mpz_clear(oldU); mpz_clear(oldV);
  mpz_set(f, n);
  return 0;
}

int _GMP_pbrent_factor(mpz_t n, mpz_t f, UV a, UV rounds)
{
  mpz_t Xi, Xm, saveXi, m, t;
  UV i, r;
  const UV inner = 256;

  TEST_FOR_2357(n, f);
  mpz_init_set_ui(Xi, 2);
  mpz_init_set_ui(Xm, 2);
  mpz_init(m);
  mpz_init(t);
  mpz_init(saveXi);

  r = 1;
  while (rounds > 0) {
    UV rleft = (r > rounds) ? rounds : r;
    while (rleft > 0) {   /* Do rleft rounds, inner at a time */
      UV dorounds = (rleft > inner) ? inner : rleft;
      mpz_set_ui(m, 1);
      mpz_set(saveXi, Xi);
      for (i = 0; i < dorounds; i++) {
        mpz_mul(t, Xi, Xi);  mpz_add_ui(t, t, a);  mpz_tdiv_r(Xi, t, n);
        if (mpz_cmp(Xi, Xm) >= 0)  mpz_sub(f, Xi, Xm);
        else                       mpz_sub(f, Xm, Xi);
        mpz_mul(t, m, f);
        mpz_tdiv_r(m, t, n);
      }
      rleft -= dorounds;
      rounds -= dorounds;
      mpz_gcd(f, m, n);
      if (mpz_cmp_ui(f, 1) != 0)
        break;
    }
    if (!mpz_cmp_ui(f, 1)) {
      r *= 2;
      mpz_set(Xm, Xi);
      continue;
    }
    if (!mpz_cmp(f, n)) {
      /* f == n, so we have to back up to see what factor got found */
      mpz_set(Xi, saveXi);
      do {
        mpz_mul(t, Xi, Xi);  mpz_add_ui(t, t, a);  mpz_tdiv_r(Xi, t, n);
        if (mpz_cmp(Xi, Xm) >= 0)  mpz_sub(f, Xi, Xm);
        else                       mpz_sub(f, Xm, Xi);
        mpz_gcd(f, f, n);
      } while (!mpz_cmp_ui(f, 1) && r-- != 0);
      if ( (!mpz_cmp_ui(f, 1)) || (!mpz_cmp(f, n)) )  break;
    }
    mpz_clear(Xi); mpz_clear(Xm); mpz_clear(m); mpz_clear(saveXi); mpz_clear(t);
    return 1;
  }
  mpz_clear(Xi); mpz_clear(Xm); mpz_clear(m); mpz_clear(saveXi); mpz_clear(t);
  mpz_set(f, n);
  return 0;
}


void _GMP_lcm_of_consecutive_integers(UV B, mpz_t m)
{
  UV i, p, p_power, pmin;
  mpz_t t[8];
  PRIME_ITERATOR(iter);

  /* Create eight sub-products to combine at the end. */
  for (i = 0; i < 8; i++)  mpz_init_set_ui(t[i], 1);
  i = 0;
  /* For each prime, multiply m by p^floor(log B / log p), which means
   * raise p to the largest power e such that p^e <= B.
   */
  if (B >= 2) {
    p_power = 2;
    while (p_power <= B/2)
      p_power *= 2;
    mpz_mul_ui(t[i&7], t[i&7], p_power);  i++;
  }
  p = prime_iterator_next(&iter);
  while (p <= B) {
    pmin = B/p;
    if (p > pmin)
      break;
    p_power = p*p;
    while (p_power <= pmin)
      p_power *= p;
    mpz_mul_ui(t[i&7], t[i&7], p_power);  i++;
    p = prime_iterator_next(&iter);
  }
  while (p <= B) {
    mpz_mul_ui(t[i&7], t[i&7], p);  i++;
    p = prime_iterator_next(&iter);
  }
  /* combine the eight sub-products */
  for (i = 0; i < 4; i++)  mpz_mul(t[i], t[2*i], t[2*i+1]);
  for (i = 0; i < 2; i++)  mpz_mul(t[i], t[2*i], t[2*i+1]);
  mpz_mul(m, t[0], t[1]);
  for (i = 0; i < 8; i++)  mpz_clear(t[i]);
  prime_iterator_destroy(&iter);
}


int _GMP_pminus1_factor(mpz_t n, mpz_t f, UV B1, UV B2)
{
  mpz_t a, savea, t;
  UV q, saveq, j, sqrtB1;
  int _verbose = get_verbose_level();
  PRIME_ITERATOR(iter);

  TEST_FOR_2357(n, f);
  if (B1 < 7) return 0;

  mpz_init(a);
  mpz_init(savea);
  mpz_init(t);

  if (_verbose>2) gmp_printf("# p-1 trying %Zd (B1=%"UVuf" B2=%"UVuf")\n", n, (unsigned long)B1, (unsigned long)B2);

  /* STAGE 1
   * Montgomery 1987 p249-250 and Brent 1990 p5 both indicate we can calculate
   * a^m mod n where m is the lcm of the integers to B1.  This can be done
   * using either
   *    m = calc_lcm(B), b = a^m mod n
   * or
   *    calculate_b_lcm(b, B1, a, n);
   *
   * The first means raising a to a huge power then doing the mod, which is
   * inefficient and can be _very_ slow on some machines.  The latter does
   * one powmod for each prime power, which works pretty well.  Yet another
   * way to handle this is to loop over each prime p below B1, calculating
   * a = a^(p^e) mod n, where e is the largest e such that p^e <= B1.
   * My experience with GMP is that this last method is faster with large B1,
   * sometimes a lot faster.
   *
   * One thing that can speed things up quite a bit is not running the GCD
   * on every step.  However with small factors this means we can easily end
   * up with multiple factors between GCDs, so we allow backtracking.  This
   * could also be added to stage 2, but it's far less likely to happen there.
   */
  j = 15;
  mpz_set_ui(a, 2);
  mpz_set_ui(savea, 2);
  saveq = 2;
  /* We could wrap this in a loop trying a few different a values, in case
   * the current one ended up going to 0. */
  q = 2;
  mpz_set_ui(t, 1);
  sqrtB1 = (UV) sqrt(B1);
  while (q <= B1) {
    UV k = q;
    if (q <= sqrtB1) {
      UV kmin = B1/q;
      while (k <= kmin)
        k *= q;
    }
    mpz_mul_ui(t, t, k);        /* Accumulate powers for a */
    if ( (j++ % 32) == 0) {
      mpz_powm(a, a, t, n);     /* a=a^(k1*k2*k3*...) mod n */
      if (mpz_sgn(a))  mpz_sub_ui(t, a, 1);
      else             mpz_sub_ui(t, n, 1);
      mpz_gcd(f, t, n);         /* f = gcd(a-1, n) */
      mpz_set_ui(t, 1);
      if (mpz_cmp(f, n) == 0)
        break;
      if (mpz_cmp_ui(f, 1) != 0)
        goto end_success;
      saveq = q;
      mpz_set(savea, a);
    }
    q = prime_iterator_next(&iter);
  }
  mpz_powm(a, a, t, n);
  if (mpz_sgn(a))  mpz_sub_ui(t, a, 1);
  else             mpz_sub_ui(t, n, 1);
  mpz_gcd(f, t, n);
  if (mpz_cmp(f, n) == 0) {
    /* We found multiple factors.  Loop one at a time. */
    prime_iterator_setprime(&iter, saveq);
    mpz_set(a, savea);
    for (q = saveq; q <= B1; q = prime_iterator_next(&iter)) {
      UV k = q;
      if (q <= sqrtB1) {
        UV kmin = B1/q;
        while (k <= kmin)
          k *= q;
      }
      mpz_powm_ui(a, a, k, n );
      mpz_sub_ui(t, a, 1);
      mpz_gcd(f, t, n);
      if (mpz_cmp(f, n) == 0)
        goto end_fail;
      if (mpz_cmp_ui(f, 1) != 0)
        goto end_success;
    }
  }
  if ( (mpz_cmp_ui(f, 1) != 0) && (mpz_cmp(f, n) != 0) )
    goto end_success;

  /* STAGE 2
   * This is the standard continuation which replaces the powmods in stage 1
   * with two mulmods, with a GCD every 64 primes (no backtracking).
   * This is quite a bit faster than stage 1.
   * See Montgomery 1987, p250-253 for possible optimizations.
   * We quickly precalculate a few of the prime gaps, and lazily cache others
   * up to a gap of 222.  That's enough for a B2 value of 189 million.  We
   * still work above that, we just won't cache the value for big gaps.
   */
  if (B2 > B1) {
    mpz_t b, bm, bmdiff;
    mpz_t precomp_bm[111];
    int   is_precomp[111] = {0};
    UV* primes = 0;
    UV sp = 1;

    mpz_init(bmdiff);
    mpz_init_set(bm, a);
    mpz_init_set_ui(b, 1);

    /* Set the first 20 differences */
    mpz_powm_ui(bmdiff, bm, 2, n);
    mpz_init_set(precomp_bm[0], bmdiff);
    is_precomp[0] = 1;
    for (j = 1; j < 22; j++) {
      mpz_mul(bmdiff, bmdiff, bm);
      mpz_mul(bmdiff, bmdiff, bm);
      mpz_tdiv_r(bmdiff, bmdiff, n);
      mpz_init_set(precomp_bm[j], bmdiff);
      is_precomp[j] = 1;
    }

    mpz_powm_ui(a, a, q, n );
    if (B2 < 10000000) {
      /* grab all the primes at once.  Hack around non-perfect iterator. */
      primes = sieve_to_n(B2+300, 0);
      for (sp = B1>>4; primes[sp] <= q; sp++)  ;
      /* q is primes <= B1, primes[sp] is the next prime */
    }

    j = 31;
    while (q <= B2) {
      UV lastq, qdiff;

      lastq = q;
      q = primes ? primes[sp++] : prime_iterator_next(&iter);
      qdiff = (q - lastq) / 2 - 1;

      if (qdiff < 111 && is_precomp[qdiff]) {
        mpz_mul(t, a, precomp_bm[qdiff]);
      } else if (qdiff < 111) {
        mpz_powm_ui(bmdiff, bm, q-lastq, n);
        mpz_init_set(precomp_bm[qdiff], bmdiff);
        is_precomp[qdiff] = 1;
        mpz_mul(t, a, bmdiff);
      } else {
        mpz_powm_ui(bmdiff, bm, q-lastq, n);  /* Big gap */
        mpz_mul(t, a, bmdiff);
      }
      mpz_tdiv_r(a, t, n);
      if (mpz_sgn(a))  mpz_sub_ui(t, a, 1);
      else             mpz_sub_ui(t, n, 1);
      mpz_mul(b, b, t);
      if ((j % 2) == 0)           /* put off mods a little */
        mpz_tdiv_r(b, b, n);
      if ( (j++ % 64) == 0) {     /* GCD every so often */
        mpz_gcd(f, b, n);
        if ( (mpz_cmp_ui(f, 1) != 0) && (mpz_cmp(f, n) != 0) )
          break;
      }
    }
    mpz_gcd(f, b, n);
    mpz_clear(b);
    mpz_clear(bm);
    mpz_clear(bmdiff);
    for (j = 0; j < 111; j++) {
      if (is_precomp[j])
        mpz_clear(precomp_bm[j]);
    }
    if (primes != 0) Safefree(primes);
    if ( (mpz_cmp_ui(f, 1) != 0) && (mpz_cmp(f, n) != 0) )
      goto end_success;
  }

  end_fail:
    mpz_set(f,n);
  end_success:
    prime_iterator_destroy(&iter);
    mpz_clear(a);
    mpz_clear(savea);
    mpz_clear(t);
    if ( (mpz_cmp_ui(f, 1) != 0) && (mpz_cmp(f, n) != 0) ) {
      if (_verbose>2) gmp_printf("# p-1: %Zd\n", f);
      return 1;
    }
    if (_verbose>2) gmp_printf("# p-1: no factor\n");
    mpz_set(f, n);
    return 0;
}

static void pp1_pow(mpz_t X, mpz_t Y, unsigned long exp, mpz_t n)
{
  mpz_t x0;
  unsigned long bit;
  {
    unsigned long v = exp;
    unsigned long b = 1;
    while (v >>= 1) b++;
    bit = 1UL << (b-2);
  }
  mpz_init_set(x0, X);
  mpz_mul(Y, X, X);
  mpz_sub_ui(Y, Y, 2);
  mpz_tdiv_r(Y, Y, n);
  while (bit) {
    if ( exp & bit ) {
      mpz_mul(X, X, Y);
      mpz_sub(X, X, x0);
      mpz_mul(Y, Y, Y);
      mpz_sub_ui(Y, Y, 2);
    } else {
      mpz_mul(Y, X, Y);
      mpz_sub(Y, Y, x0);
      mpz_mul(X, X, X);
      mpz_sub_ui(X, X, 2);
    }
    mpz_mod(X, X, n);
    mpz_mod(Y, Y, n);
    bit >>= 1;
  }
  mpz_clear(x0);
}

int _GMP_pplus1_factor(mpz_t n, mpz_t f, UV P0, UV B1, UV B2)
{
  UV j, q, saveq, sqrtB1;
  mpz_t X, Y, saveX;
  PRIME_ITERATOR(iter);

  TEST_FOR_2357(n, f);
  if (B1 < 7) return 0;

  mpz_init_set_ui(X, P0);
  mpz_init(Y);
  mpz_init(saveX);

  /* Montgomery 1987 */
  if (P0 == 0) {
    mpz_set_ui(X, 7);
    if (mpz_invert(X, X, n)) {
      mpz_mul_ui(X, X, 2);
      mpz_mod(X, X, n);
    } else
      P0 = 1;
  }
  if (P0 == 1) {
    mpz_set_ui(X, 5);
    if (mpz_invert(X, X, n)) {
      mpz_mul_ui(X, X, 6);
      mpz_mod(X, X, n);
    } else
      P0 = 2;
  }
  if (P0 == 2) {
    mpz_set_ui(X, 11);
    if (mpz_invert(X, X, n)) {
      mpz_mul_ui(X, X, 23);
      mpz_mod(X, X, n);
    }
  }

  sqrtB1 = (UV) sqrt(B1);
  j = 8;
  q = 2;
  saveq = q;
  mpz_set(saveX, X);
  while (q <= B1) {
    UV k = q;
    if (q <= sqrtB1) {
      UV kmin = B1/q;
      while (k <= kmin)
        k *= q;
    }
    pp1_pow(X, Y, k, n);
    if ( (j++ % 16) == 0) {
      mpz_sub_ui(f, X, 2);
      if (mpz_sgn(f) == 0)        break;
      mpz_gcd(f, f, n);
      if (mpz_cmp(f, n) == 0)     break;
      if (mpz_cmp_ui(f, 1) > 0)   goto end_success;
      saveq = q;
      mpz_set(saveX, X);
    }
    q = prime_iterator_next(&iter);
  }
  mpz_sub_ui(f, X, 2);
  mpz_gcd(f, f, n);
  if (mpz_cmp_ui(X, 2) == 0 || mpz_cmp(f, n) == 0) {
    /* Backtrack */
    prime_iterator_setprime(&iter, saveq);
    mpz_set(X, saveX);
    for (q = saveq; q <= B1; q = prime_iterator_next(&iter)) {
      UV k = q;
      if (q <= sqrtB1) {
        UV kmin = B1/q;
        while (k <= kmin)
          k *= q;
      }
      pp1_pow(X, Y, k, n);
      mpz_sub_ui(f, X, 2);
      if (mpz_sgn(f) == 0)        goto end_fail;
      mpz_gcd(f, f, n);
      if (mpz_cmp(f, n) == 0)     break;
      if (mpz_cmp_ui(f, 1) > 0)   goto end_success;
    }
  }
  if ( (mpz_cmp_ui(f, 1) > 0) && (mpz_cmp(f, n) != 0) )
    goto end_success;
  /* TODO: stage 2 */
  end_fail:
    mpz_set(f,n);
  end_success:
    prime_iterator_destroy(&iter);
    mpz_clear(X);  mpz_clear(Y);  mpz_clear(saveX);
    return (mpz_cmp_ui(f, 1) != 0) && (mpz_cmp(f, n) != 0);
}

int _GMP_holf_factor(mpz_t n, mpz_t f, UV rounds)
{
  mpz_t s, m;
  UV i;

#define PREMULT 480   /* 1  2  6  12  480  151200 */

  TEST_FOR_2357(n, f);
  if (mpz_perfect_square_p(n)) {
    mpz_sqrt(f, n);
    return 1;
  }

  mpz_mul_ui(n, n, PREMULT);
  mpz_init(s);
  mpz_init(m);
  for (i = 1; i <= rounds; i++) {
    mpz_mul_ui(f, n, i);    /* f = n*i */
    if (mpz_perfect_square_p(f)) {
      /* s^2 = n*i, so m = s^2 mod n = 0.  Hence f = GCD(n, s) = GCD(n, n*i) */
      mpz_divexact_ui(n, n, PREMULT);
      mpz_gcd(f, f, n);
      mpz_clear(s); mpz_clear(m);
      if (mpz_cmp(f, n) == 0)  return 0;
      return 1;
    }
    mpz_sqrt(s, f);
    mpz_add_ui(s, s, 1);    /* s = ceil(sqrt(n*i)) */
    mpz_mul(m, s, s);
    mpz_sub(m, m, f);       /* m = s^2 mod n = s^2 - n*i */
    if (mpz_perfect_square_p(m)) {
      mpz_divexact_ui(n, n, PREMULT);
      mpz_sqrt(f, m);
      mpz_sub(s, s, f);
      mpz_gcd(f, s, n);
      mpz_clear(s); mpz_clear(m);
      return (mpz_cmp_ui(f, 1) > 0);
    }
  }
  mpz_divexact_ui(n, n, PREMULT);
  mpz_set(f, n);
  mpz_clear(s); mpz_clear(m);
  return 0;
}


/*----------------------------------------------------------------------
 * GMP version of Ben Buhrow's public domain 9/24/09 implementation.
 * It uses ideas and code from Jason Papadopoulos, Scott Contini, and
 * Tom St. Denis.  Also see the papers of Stephen McMath, Daniel Shanks,
 * and Jason Gower.  Gower and Wagstaff is particularly useful:
 *    http://homes.cerias.purdue.edu/~ssw/squfof.pdf
 *--------------------------------------------------------------------*/

static int shanks_mult(mpz_t n, mpz_t f)
{
   /*
    * use shanks SQUFOF to factor N.
    *
    * return 0 if no factor found, 1 if found with factor in f1.
    *
    * Input should have gone through trial division to 5.
    */

   int result = 0;
   unsigned long j=0;
   mpz_t b0, bn, imax, tmp, Q0, Qn, P, i, t1, t2, S, Ro, So, bbn;

   if (mpz_cmp_ui(n, 3) <= 0)
     return 0;

   if (mpz_perfect_square_p(n)) {
     mpz_sqrt(f, n);
     return 1;
   }

   mpz_init(b0);
   mpz_init(bn);
   mpz_init(imax);
   mpz_init(tmp);
   mpz_init(Q0);
   mpz_init(Qn);
   mpz_init(P);
   mpz_init(i);
   mpz_init(t1);
   mpz_init(t2);
   mpz_init(S);
   mpz_init(Ro);
   mpz_init(So);
   mpz_init(bbn);

   mpz_mod_ui(t1, n, 4);
   if (mpz_cmp_ui(t1, 3))
     croak("Incorrect call to shanks\n");

   mpz_sqrt(b0, n);
   mpz_sqrt(tmp, b0);
   mpz_mul_ui(imax, tmp, 3);

   /* set up recurrence */
   mpz_set_ui(Q0, 1);
   mpz_set(P, b0);
   mpz_mul(tmp, b0, b0);
   mpz_sub(Qn, n, tmp);

   mpz_add(tmp, b0, P);
   mpz_tdiv_q(bn, tmp, Qn);

   mpz_set_ui(i, 0);
   while (1) {
      j=0;
      while (1) {
         mpz_set(t1, P);   /* hold Pn for this iteration */
         mpz_mul(tmp, bn, Qn);
         mpz_sub(P, tmp, P);
         mpz_set(t2, Qn);  /* hold Qn for this iteration */
         mpz_sub(tmp, t1, P);
         mpz_mul(tmp, tmp, bn);
         mpz_add(Qn, Q0, tmp);
         mpz_set(Q0, t2);  /* remember last Q */
         mpz_add(tmp, b0, P);
         mpz_tdiv_q(bn, tmp, Qn);

         if (mpz_even_p(i)) {
           if (mpz_perfect_square_p(Qn)) {
             mpz_add_ui(i, i, 1);
             break;
           }
         }
         mpz_add_ui(i, i, 1);

         if (mpz_cmp(i, imax) >= 0) {
           result = 0;
           goto end;
         }
      }

      /* reduce to G0 */
      mpz_sqrt(S, Qn);
      mpz_sub(tmp, b0, P);
      mpz_tdiv_q(tmp, tmp, S);
      mpz_mul(tmp, S, tmp);
      mpz_add(Ro, P, tmp);
      mpz_mul(tmp, Ro, Ro);
      mpz_sub(tmp, n, tmp);
      mpz_tdiv_q(So, tmp, S);
      mpz_add(tmp, b0, Ro);
      mpz_tdiv_q(bbn, tmp, So);

      /* search for symmetry point */
      while (1) {
         mpz_set(t1, Ro);  /* hold Ro for this iteration */
         mpz_mul(tmp, bbn, So);
         mpz_sub(Ro, tmp, Ro);
         mpz_set(t2, So);  /* hold So for this iteration */
         mpz_sub(tmp, t1, Ro);
         mpz_mul(tmp, bbn, tmp);
         mpz_add(So, S, tmp);
         mpz_set(S, t2);   /* remember last S */
         mpz_add(tmp, b0, Ro);
         mpz_tdiv_q(bbn, tmp, So);

         /* check for symmetry point */
         if (mpz_cmp(Ro, t1) == 0)
            break;

         /* this gets stuck very rarely, but it does happen. */
         if (++j > 1000000000)
         {
            result = -1;
            goto end;
         }
      }

      mpz_gcd(t1, Ro, n);
      if (mpz_cmp_ui(t1, 1) > 0) {
         mpz_set(f, t1);
         /* gmp_printf("GMP SQUFOF found factor after %Zd/%lu rounds: %Zd\n", i, j, f); */
         result = 1;
         goto end;
      }
   }

   end:
   mpz_clear(b0);
   mpz_clear(bn);
   mpz_clear(imax);
   mpz_clear(tmp);
   mpz_clear(Q0);
   mpz_clear(Qn);
   mpz_clear(P);
   mpz_clear(i);
   mpz_clear(t1);
   mpz_clear(t2);
   mpz_clear(S);
   mpz_clear(Ro);
   mpz_clear(So);
   mpz_clear(bbn);
   return result;
}

int _GMP_squfof_factor(mpz_t n, mpz_t f, UV rounds)
{
   const UV multipliers[] = {
      3*5*7*11, 3*5*7, 3*5*11, 3*5, 3*7*11, 3*7, 5*7*11, 5*7,
      3*11,     3,     5*11,   5,   7*11,   7,   11,     1   };
   const size_t sz_mul = sizeof(multipliers)/sizeof(multipliers[0]);
   size_t i;
   int result;
   UV nmod4;
   mpz_t nm, t;

   TEST_FOR_2357(n, f);
   mpz_init(nm);
   mpz_init(t);
   mpz_set_ui(f, 1);
   nmod4 = mpz_tdiv_r_ui(nm, n, 4);

   for (i = 0; i < sz_mul; i++) {
      UV mult = multipliers[i];
      /* Only use multipliers where n*m = 3 mod 4 */
      if (nmod4 == (mult % 4))
        continue;
      /* Only run when 64*m^3 < n */
      mpz_set_ui(t, mult);
      mpz_pow_ui(t, t, 3);
      mpz_mul_ui(t, t, 64);
      if (mpz_cmp(t, n) >= 0)
        continue;
      /* Run with this multiplier */
      mpz_mul_ui(nm, n, mult);
      result = shanks_mult(nm, f);
      if (result == -1)
        break;
      if ( (result == 1) && (mpz_cmp_ui(f, mult) != 0) ) {
        unsigned long gcdf = mpz_gcd_ui(NULL, f, mult);
        mpz_divexact_ui(f, f, gcdf);
        if (mpz_cmp_ui(f, 1) > 0)
          break;
      }
   }
   mpz_clear(t);
   mpz_clear(nm);
   return (mpz_cmp_ui(f, 1) > 0);
}

/* See if n is a perfect power */
UV power_factor(mpz_t n, mpz_t f)
{
  UV k = 1, b = 2;
  if (mpz_cmp_ui(n, 1) > 0 && mpz_perfect_power_p(n)) {
    mpz_t nf, tf;
    PRIME_ITERATOR(iter);

    mpz_init_set(nf, n);
    mpz_init(tf);
    while (1) {
      UV ok = k;
      while (mpz_root(tf, nf, b)) {
        mpz_set(f, tf);
        mpz_set(nf, tf);
        k *= b;
      }
      if (ok != k && !mpz_perfect_power_p(nf)) break;
      if (mpz_cmp_ui(tf, 1) <= 0) break; /* Exit if we can't find the power */
      b = prime_iterator_next(&iter);
    }
    mpz_clear(tf);  mpz_clear(nf);
    prime_iterator_destroy(&iter);
  }
  return (k == 1) ? 0 : k;
}

/* a=0, return power.  a>1, return bool if an a-th power */
UV is_power(mpz_t n, UV a)
{
  if (mpz_cmp_ui(n,3) <= 0)
    return 0;
  else if (a == 1)
    return 1;
  else if (a == 2)
    return mpz_perfect_square_p(n);
  else {
    UV result;
    mpz_t t;
    mpz_init(t);
    result = (a == 0)  ?  power_factor(n, t)  :  (UV)mpz_root(t, n, a);
    mpz_clear(t);
    return result;
  }
}

void exp_mangoldt(mpz_t res, mpz_t n)
{
  UV k;
  mpz_set_ui(res, 1);
  if (mpz_cmp_ui(n, 1) <= 0)
    return;
  k = mpz_scan1(n, 0);
  if (k > 0) {
    if (k+1 == mpz_sizeinbase(n, 2))
      mpz_set_ui(res, 2);
    return;
  }
  if (_GMP_is_prob_prime(n)) {
    mpz_set(res, n);
    return;
  }
  k = power_factor(n, res);
  if (k > 1 && _GMP_is_prob_prime(res))
    return;
  mpz_set_ui(res, 1);
}


static void word_tile(uint32_t* source, uint32_t from, uint32_t to) {
  while (from < to) {
    uint32_t words = (2*from > to) ? to-from : from;
    memcpy(source+from, source, sizeof(uint32_t)*words);
    from += words;
  }
}
static void sievep(uint32_t* comp, mpz_t start, UV p, UV len) {
  UV pos = p - mpz_fdiv_ui(start,p);  /* First multiple of p after start  */
  if (!(pos & 1)) pos += p;           /* Make sure it is odd.             */
  for ( ; pos < len; pos += 2*p )
    SETAVAL(comp, pos);
}
static void sievep_ui(uint32_t* comp, UV pos, UV p, UV len) {
  if (!(pos & 1)) pos += p;
  for ( ; pos < len; pos += 2*p )
    SETAVAL(comp, pos);
}
uint32_t* partial_sieve(mpz_t start, UV length, UV maxprime)
{
  uint32_t* comp;
  UV p, wlen, pwlen;
  PRIME_ITERATOR(iter);

  /* mpz_init(t);
     mpz_add_ui(t, start, (length & 1) ? length-1 : length-2);
     gmp_printf("partial sieve start %Zd  length %lu mark %Zd to %Zd\n", start, length, start, t); */
  MPUassert(mpz_odd_p(start), "partial sieve given even start");
  MPUassert(length > 0, "partial sieve given zero length");
  mpz_sub_ui(start, start, 1);
  if (length & 1) length++;

  /* Allocate odds-only array in uint32_t units */
  wlen = (length+63)/64;
  New(0, comp, wlen, uint32_t);
  p = prime_iterator_next(&iter);

  /* Mark 3, 5, ... by tiling as long as we can. */
  pwlen = (wlen < 3) ? wlen : 3;
  memset(comp, 0x00, pwlen*sizeof(uint32_t));
  while (p <= maxprime) {
    sievep(comp, start, p, pwlen*64);
    p = prime_iterator_next(&iter);
    if (pwlen*p >= wlen) break;
    word_tile(comp, pwlen, pwlen*p);
    pwlen *= p;
  }
  word_tile(comp, pwlen, wlen);

  /* Sieve up to their limit.
   *
   * Simple code for this:
   *    for ( ; p <= maxprime; p = prime_iterator_next(&iter))
   *      sievep(comp, start, p, length);
   * We'll save some time for large start values by doubling up primes.
   */
  {
    UV p1, p2, doublelim = (1UL << (sizeof(unsigned long) * 4)) - 1;
    if (doublelim > maxprime) doublelim = maxprime;
    /* Do 2 primes at a time.  Fewer mpz remainders. */
    for ( p1 = p, p2 = prime_iterator_next(&iter);
          p2 <= doublelim;
          p1 = prime_iterator_next(&iter), p2 = prime_iterator_next(&iter) ) {
      UV p1p2 = p1 * p2;
      UV ddiv = mpz_fdiv_ui(start, p1p2);
      sievep_ui(comp, p1 - (ddiv % p1), p1, length);
      sievep_ui(comp, p2 - (ddiv % p2), p2, length);
    }
    if (p1 <= maxprime) sievep(comp, start, p1, length);
    for (p = p2; p <= maxprime; p = prime_iterator_next(&iter))
      sievep(comp, start, p, length);
  }

  prime_iterator_destroy(&iter);
  return comp;
}


char* pidigits(UV n) {
  char* out;

  New(0, out, n+4, char);
  out[0] = '3';  out[1] = '\0';
  if (n <= 1)
    return out;

#if 0
  /* Spigot method.
   * ~40x slower than the Machin formulas, 2x slower than spigot in plain C */
  {
    mpz_t t1, t2, acc, den, num;
    UV i, k, d, d2;

    mpz_init(t1);  mpz_init(t2);
    mpz_init_set_ui(acc, 0);
    mpz_init_set_ui(den, 1);
    mpz_init_set_ui(num, 1);
    n++;   /* rounding */
    for (i = k = 0; i < n; ) {
      {
        UV k2 = ++k * 2 + 1;
        mpz_mul_2exp(t1, num, 1);
        mpz_add(acc, acc, t1);
        mpz_mul_ui(acc, acc, k2);
        mpz_mul_ui(den, den, k2);
        mpz_mul_ui(num, num, k);
      }
      if (mpz_cmp(num, acc) > 0)  continue;
      {
        mpz_mul_ui(t1, num, 3);
        mpz_add(t2, t1, acc);
        mpz_tdiv_q(t1, t2, den);
        d = mpz_get_ui(t1);
      }
      {
        mpz_mul_ui(t1, num, 4);
        mpz_add(t2, t1, acc);
        mpz_tdiv_q(t1, t2, den);
        d2 = mpz_get_ui(t1);
      }
      if (d != d2)  continue;
      out[++i] = '0' + d;
      {
        mpz_submul_ui(acc, den, d);
        mpz_mul_ui(acc, acc, 10);
        mpz_mul_ui(num, num, 10);
      }
    }
    mpz_clear(num); mpz_clear(den); mpz_clear(acc); mpz_clear(t2);mpz_clear(t1);
    if (out[n] >= '5') out[n-1]++;  /* Round */
    for (i = n-1; out[i] == '9'+1; i--)    /* Keep rounding */
      { out[i] = '0';  out[i-1]++; }
    n--;  /* Undo the extra digit we used for rounding */
    out[1] = '.';
    out[n+1] = '\0';
  }
#elif 0
  /* https://en.wikipedia.org/wiki/Machin-like_formula
   * Thanks to Ledrug from RosettaCode for the simple code for base 10.
   * Pretty fast, but growth is a lot slower than AGM. */
  {
    mpz_t t1, t2, term1, term2, pows;
    UV i, k;

    mpz_init(t1); mpz_init(t2); mpz_init(term1); mpz_init(term2); mpz_init(pows);
    n++;   /* rounding */
    mpz_ui_pow_ui(pows, 10, n+20);

#if 0
    /* Machin 1706 */
    mpz_arctan(term1,       5, pows, t1, t2);  mpz_mul_ui(term1, term1, 4);
    mpz_arctan(term2,     239, pows, t1, t2);
    mpz_sub(term1, term1, term2);
#elif 0
    /* Störmer 1896 */
    mpz_arctan(term1,      57, pows, t1, t2);  mpz_mul_ui(term1, term1, 44);
    mpz_arctan(term2,     239, pows, t1, t2);  mpz_mul_ui(term2, term2, 7);
    mpz_add(term1, term1, term2);
    mpz_arctan(term2,     682, pows, t1, t2);  mpz_mul_ui(term2, term2, 12);
    mpz_sub(term1, term1, term2);
    mpz_arctan(term2,   12943, pows, t1, t2);  mpz_mul_ui(term2, term2, 24);
    mpz_add(term1, term1, term2);
#else
    /* Chien-Lih 1997 */
    mpz_arctan(term1,     239, pows, t1, t2);  mpz_mul_ui(term1, term1, 183);
    mpz_arctan(term2,    1023, pows, t1, t2);  mpz_mul_ui(term2, term2,  32);
    mpz_add(term1, term1, term2);
    mpz_arctan(term2,    5832, pows, t1, t2);  mpz_mul_ui(term2, term2,  68);
    mpz_sub(term1, term1, term2);
    mpz_arctan(term2,  110443, pows, t1, t2);  mpz_mul_ui(term2, term2,  12);
    mpz_add(term1, term1, term2);
    mpz_arctan(term2, 4841182, pows, t1, t2);  mpz_mul_ui(term2, term2,  12);
    mpz_sub(term1, term1, term2);
    mpz_arctan(term2, 6826318, pows, t1, t2);  mpz_mul_ui(term2, term2, 100);
    mpz_sub(term1, term1, term2);
#endif
    mpz_mul_ui(term1, term1, 4);

    mpz_ui_pow_ui(pows, 10, 20);
    mpz_tdiv_q(term1, term1, pows);

    mpz_clear(t1); mpz_clear(t2); mpz_clear(term2); mpz_clear(pows);

    k = mpz_sizeinbase(term1, 10);           /* Copy result to out string */
    while (k > n+1) {                        /* making sure we don't overflow */
      mpz_tdiv_q_ui(term1, term1, 10);
      k = mpz_sizeinbase(term1, 10);
    }
    (void) mpz_get_str(out+1, 10, term1);
    mpz_clear(term1);

    if (out[n] >= '5') out[n-1]++;  /* Round */
    for (i = n-1; out[i] == '9'+1; i--)    /* Keep rounding */
      { out[i] = '0';  out[i-1]++; }
    n--;  /* Undo the extra digit we used for rounding */
    out[1] = '.';
    out[n+1] = '\0';
  }
#else
  /* AGM using GMP's floating point.  Fast and very good growth. */
  {
    mpf_t t, an, bn, tn, prev_an;
    UV k = 0;
    unsigned long oldprec = mpf_get_default_prec();
    mpf_set_default_prec(10 + n * 3.322);

    mpf_init(t);  mpf_init(prev_an);
    mpf_init_set_d(an, 1);  mpf_init_set_d(bn, 0.5);  mpf_init_set_d(tn, 0.25);
    mpf_sqrt(bn, bn);

    while ((n >> k) > 0) {
      mpf_set(prev_an, an);
      mpf_add(t, an, bn);
      mpf_div_ui(an, t, 2);
      mpf_mul(t, bn, prev_an);
      mpf_sqrt(bn, t);
      mpf_sub(prev_an, prev_an, an);
      mpf_mul(t, prev_an, prev_an);
      mpf_mul_2exp(t, t, k);
      mpf_sub(tn, tn, t);
      k++;
    }
    mpf_add(t, an, bn);
    mpf_mul(an, t, t);
    mpf_mul_2exp(t, tn, 2);
    mpf_div(bn, an, t);
    gmp_sprintf(out, "%.*Ff", (int)(n-1), bn);
    mpf_clear(tn); mpf_clear(bn); mpf_clear(an);
    mpf_clear(prev_an); mpf_clear(t);
    mpf_set_default_prec(oldprec);
  }
#endif
  return out;
}

typedef struct {
  uint32_t nmax;
  uint32_t nsize;
  UV* list;
} vlist;
#define INIT_VLIST(v) \
  v.nsize = 0; \
  v.nmax = 1024; \
  New(0, v.list, v.nmax, UV);
#define RESIZE_VLIST(v, size) \
  do { if (v.nmax < size) Renew(v.list, v.nmax = size, UV); } while (0)
#define PUSH_VLIST(v, n) \
  do { \
    if (v.nsize >= v.nmax) \
      Renew(v.list, v.nmax += 1024, UV); \
    v.list[v.nsize++] = n; \
  } while (0)

#define ADDVAL32(v, n, max, val) \
  do { if (n >= max) Renew(v, max += 1024, UV);  v[n++] = val; } while (0)
#define SWAPL32(l1, n1, m1,  l2, n2, m2) \
  { UV t_, *u_ = l1;  l1 = l2;  l2 = u_; \
            t_ = n1;  n1 = n2;  n2 = t_; \
            t_ = m1;  m1 = m2;  m2 = t_; }

UV* sieve_primes(mpz_t inlow, mpz_t high, UV k, UV *rn) {
  mpz_t t, low;
  int test_primality = 0;
  int k_primality = 0;
  uint32_t* comp;
  vlist retlist;

  if (k < 2) {
    test_primality = 1;
    k = 5000 * mpz_sizeinbase(high,2);
  }

  if (mpz_cmp_ui(inlow, 2) < 0) mpz_set_ui(inlow, 2);
  if (mpz_cmp(inlow, high) > 0) { *rn = 0; return 0; }

  mpz_init(t);
  mpz_sqrt(t, high);           /* No need for k to be > sqrt(high) */
  if (mpz_cmp_ui(t, k) < 0)
    k = mpz_get_ui(t);
  if (mpz_cmp_ui(t, k) <= 0) {
    k_primality = 1;           /* Our sieve is complete */
    test_primality = 0;        /* Don't run BPSW */
  }

  INIT_VLIST(retlist);

  /* If we want small primes, do it quickly */
  if ( (k_primality || test_primality) && mpz_cmp_ui(high,2000000000U) <= 0 ) {
    UV ulow = mpz_get_ui(inlow), uhigh = mpz_get_ui(high);
    if (uhigh < 1000000U || uhigh/ulow >= 4) {
      UV n, Pi, *primes;
      primes = sieve_to_n(mpz_get_ui(high), &Pi);
      RESIZE_VLIST(retlist, Pi);
      for (n = 0; n < Pi; n++) {
        if (primes[n] >= ulow)
          PUSH_VLIST(retlist, primes[n]-ulow);
      }
      Safefree(primes);
      mpz_clear(t);
      *rn = retlist.nsize;
      return retlist.list;
    }
  }

  mpz_init_set(low, inlow);
  if (k < 2) k = 2;   /* Should have been handled by quick return */

  /* Include all primes up to k, since they will get filtered */
  if (mpz_cmp_ui(low, k) <= 0) {
    UV n, Pi, *primes, ulow = mpz_get_ui(low);
    primes = sieve_to_n(k, &Pi);
    RESIZE_VLIST(retlist, Pi);
    for (n = 0; n < Pi; n++) {
      if (primes[n] >= ulow)
        PUSH_VLIST(retlist, primes[n]-ulow);
    }
    Safefree(primes);
  }

  if (mpz_even_p(low))           mpz_add_ui(low, low, 1);
  if (mpz_even_p(high))          mpz_sub_ui(high, high, 1);

  if (mpz_cmp(low, high) <= 0) {
    UV i, length, offset;
    mpz_sub(t, high, low); length = mpz_get_ui(t) + 1;
    /* Get bit array of odds marked with composites(k) marked with 1 */
    comp = partial_sieve(low, length, k);
    mpz_sub(t, low, inlow); offset = mpz_get_ui(t);
    for (i = 1; i <= length; i += 2) {
      if (!TSTAVAL(comp, i)) {
        if (!test_primality || (mpz_add_ui(t,low,i),_GMP_BPSW(t)))
          PUSH_VLIST(retlist, i - offset);
      }
    }
    Safefree(comp);
  }

  mpz_clear(low);
  mpz_clear(t);
  *rn = retlist.nsize;
  return retlist.list;
}

UV* sieve_twin_primes(mpz_t low, mpz_t high, UV twin, UV *rn) {
  mpz_t t;
  UV i, length, k, starti = 1, skipi = 2;
  uint32_t* comp;
  vlist retlist;

  /* Twin offset will always be even, so we will never return 2 */
  MPUassert( !(twin & 1), "twin prime offset is even" );
  if (mpz_cmp_ui(low, 3) <= 0) mpz_set_ui(low, 3);

  if (mpz_even_p(low))  mpz_add_ui(low, low, 1);
  if (mpz_even_p(high)) mpz_sub_ui(high, high, 1);

  i = twin % 6;
  if (i == 2 || i == 4) { skipi = 6; starti = ((twin%6)==2) ? 5 : 1; }

  /* If no way to return any more results, leave now */
  if (mpz_cmp(low, high) > 0 || (i == 1 || i == 3 || i == 5))
    { *rn = 0; return 0; }

  INIT_VLIST(retlist);
  mpz_init(t);

  /* Use a much higher k value than we do for primes */
  k = 80000 * mpz_sizeinbase(high,2);
  /* No need for k to be > sqrt(high) */
  mpz_sqrt(t, high);
  if (mpz_cmp_ui(t, k) < 0)
    k = mpz_get_ui(t);

  /* Handle small primes that will get sieved out */
  if (mpz_cmp_ui(low, k) <= 0) {
    UV p, ulow = mpz_get_ui(low);
    PRIME_ITERATOR(iter);
    for (p = 2; p <= k; p = prime_iterator_next(&iter)) {
      if (p < ulow) continue;
      if (mpz_set_ui(t, p+twin), _GMP_BPSW(t))
        PUSH_VLIST(retlist, p-ulow+1);
    }
    prime_iterator_destroy(&iter);
  }

  mpz_sub(t, high, low);
  length = mpz_get_ui(t) + 1;
  starti = ((starti+skipi) - mpz_fdiv_ui(low,skipi) + 1) % skipi;

  /* Get bit array of odds marked with composites(k) marked with 1 */
  comp = partial_sieve(low, length + twin, k);
  for (i = starti; i <= length; i += skipi) {
    if (!TSTAVAL(comp, i) && !TSTAVAL(comp, i+twin)) {
      mpz_add_ui(t, low, i);
      if (_GMP_BPSW(t)) {
        mpz_add_ui(t, t, twin);
        if (_GMP_BPSW(t)) {
          PUSH_VLIST(retlist, i);
        }
      }
    }
  }
  Safefree(comp);
  mpz_clear(t);
  *rn = retlist.nsize;
  return retlist.list;
}


#define addmodded(r,a,b,n)  do { r = a + b; if (r >= n) r -= n; } while(0)

UV* sieve_cluster(mpz_t low, mpz_t high, uint32_t* cl, UV nc, UV *rn) {
  mpz_t t, savelow;
  vlist retlist;
  UV i, ppr, nres, allocres;
  uint32_t const targres = 4000000;
  uint32_t const maxpi = 168;
  UV *residues, *cres;
  uint32_t pp_0, pp_1, pp_2, *resmod_0, *resmod_1, *resmod_2;
  uint32_t rem_0, rem_1, rem_2, remadd_0, remadd_1, remadd_2;
  uint32_t pi, startpi = 1;
  uint32_t lastspr = sprimes[maxpi-1];
  uint32_t c, smallnc;
  UV ibase = 0, nprps = 0;
  char crem_0[53*59], crem_1[61*67], crem_2[71*73], *VPrem;
  int run_pretests = 0;
  int _verbose = get_verbose_level();

  if (nc == 1) return sieve_primes(low, high, 0, rn);
  if (nc == 2) return sieve_twin_primes(low, high, cl[1], rn);

  if (mpz_even_p(low))           mpz_add_ui(low, low, 1);
  if (mpz_even_p(high))          mpz_sub_ui(high, high, 1);

  if (mpz_cmp(low, high) > 0) { *rn = 0; return 0; }

  INIT_VLIST(retlist);
  mpz_init(t);

  /* Handle small values that would get sieved away */
  if (mpz_cmp_ui(low, lastspr) <= 0) {
    UV ui_low = mpz_get_ui(low);
    UV ui_high = (mpz_cmp_ui(high,lastspr) > 0) ? lastspr : mpz_get_ui(high);
    for (pi = 0; pi < maxpi; pi++) {
      UV p = sprimes[pi];
      if (p > ui_high) break;
      if (p < ui_low) continue;
      for (c = 1; c < nc; c++)
        if (!(mpz_set_ui(t, p+cl[c]), _GMP_is_prob_prime(t))) break;
      if (c == nc)
        PUSH_VLIST(retlist, p-ui_low+1);
    }
  }
  if (mpz_odd_p(low)) mpz_sub_ui(low, low, 1);
  if (mpz_cmp_ui(high, lastspr) <= 0) {
    mpz_clear(t);
    *rn = retlist.nsize;
    return retlist.list;
  }

  /* Determine the primorial size and acceptable residues */
  New(0, residues, allocres = 1024, UV);
  {
    UV remr, *res2, allocres2, nres2, maxppr;
    /* Calculate residues for a small primorial */
    for (pi = 2, ppr = 1, i = 0;  i <= pi;  i++) ppr *= sprimes[i];
    remr = mpz_fdiv_ui(low, ppr);
    nres = 0;
    for (i = 1; i <= ppr; i += 2) {
      UV remi = remr + i;
      for (c = 0; c < nc; c++) {
        if (gcd_ui(remi + cl[c], ppr) != 1) break;
      }
      if (c == nc)
        ADDVAL32(residues, nres, allocres, i);
    }
    /* Raise primorial size until we have plenty of residues */
    New(0, res2, allocres2 = 1024, UV);
    mpz_sub(t, high, low);
    maxppr = (mpz_sizeinbase(t,2) >= BITS_PER_WORD) ? UV_MAX : (UVCONST(1) << mpz_sizeinbase(t,2));
#if BITS_PER_WORD == 64
    while (pi++ < 14) {
#else
    while (pi++ < 8) {
#endif
      uint32_t j, p = sprimes[pi];
      UV r, newppr = ppr * p;
      if (nres == 0 || nres > targres/(p/2) || newppr > maxppr) break;
      if (_verbose > 1) printf("cluster sieve found %"UVuf" residues mod %"UVuf"\n", nres, ppr);
      remr = mpz_fdiv_ui(low, newppr);
      nres2 = 0;
      for (i = 0; i < p; i++) {
        for (j = 0; j < nres; j++) {
          r = i*ppr + residues[j];
          for (c = 0; c < nc; c++) {
            UV v = remr + r + cl[c];
            if ((v % p) == 0) break;
          }
          if (c == nc)
            ADDVAL32(res2, nres2, allocres2, r);
        }
      }
      ppr = newppr;
      SWAPL32(residues, nres, allocres,  res2, nres2, allocres2);
    }
    startpi = pi;
    Safefree(res2);
  }
  if (_verbose) printf("cluster sieve using %"UVuf" residues mod %"UVuf"\n", nres, ppr);

  /* Return if not admissible, maybe with a single small value */
  if (nres == 0) {
    Safefree(residues);
    mpz_clear(t);
    *rn = retlist.nsize;
    return retlist.list;
  }

  mpz_init_set(savelow, low);
  if (mpz_sizeinbase(low, 2) > 260) run_pretests = 1;
  if (run_pretests && mpz_sgn(_bgcd2) == 0) {
    _GMP_pn_primorial(_bgcd2, BGCD2_PRIMES);
    mpz_divexact(_bgcd2, _bgcd2, _bgcd);
  }

  /* Pre-mod the residues with first two primes for fewer modulos every chunk */
  {
    uint32_t p1 = sprimes[startpi+0], p2 = sprimes[startpi+1];
    uint32_t p3 = sprimes[startpi+2], p4 = sprimes[startpi+3];
    uint32_t p5 = sprimes[startpi+4], p6 = sprimes[startpi+5];
    pp_0 = p1*p2; pp_1 = p3*p4; pp_2 = p5*p6;
    memset(crem_0, 1, pp_0);
    memset(crem_1, 1, pp_1);
    memset(crem_2, 1, pp_2);
    /* Mark remainders that indicate a composite for this residue. */
    for (i = 0; i < p1; i++) { crem_0[i*p1]=0; crem_0[i*p2]=0; }
    for (     ; i < p2; i++) { crem_0[i*p1]=0;                }
    for (i = 0; i < p3; i++) { crem_1[i*p3]=0; crem_1[i*p4]=0; }
    for (     ; i < p4; i++) { crem_1[i*p3]=0;                }
    for (i = 0; i < p5; i++) { crem_2[i*p5]=0; crem_2[i*p6]=0; }
    for (     ; i < p6; i++) { crem_2[i*p5]=0;                }
    for (c = 1; c < nc; c++) {
      uint32_t c1=cl[c], c2=cl[c], c3=cl[c], c4=cl[c], c5=cl[c], c6=cl[c];
      if (c1 >= p1) c1 %= p1;
      if (c2 >= p2) c2 %= p2;
      for (i = 1; i <= p1; i++) { crem_0[i*p1-c1]=0; crem_0[i*p2-c2]=0; }
      for (     ; i <= p2; i++) { crem_0[i*p1-c1]=0;                   }
      if (c3 >= p3) c3 %= p3;
      if (c4 >= p4) c4 %= p4;
      for (i = 1; i <= p3; i++) { crem_1[i*p3-c3]=0; crem_1[i*p4-c4]=0; }
      for (     ; i <= p4; i++) { crem_1[i*p3-c3]=0;                   }
      if (c5 >= p5) c5 %= p5;
      if (c6 >= p6) c6 %= p6;
      for (i = 1; i <= p5; i++) { crem_2[i*p5-c5]=0; crem_2[i*p6-c6]=0; }
      for (     ; i <= p6; i++) { crem_2[i*p5-c5]=0;                   }
    }
    New(0, resmod_0, nres, uint32_t);
    New(0, resmod_1, nres, uint32_t);
    New(0, resmod_2, nres, uint32_t);
    for (i = 0; i < nres; i++) {
      resmod_0[i] = residues[i] % pp_0;
      resmod_1[i] = residues[i] % pp_1;
      resmod_2[i] = residues[i] % pp_2;
    }
  }

  /* Precalculate acceptable residues for more primes */
  MPUassert( lastspr <= 1024, "cluster sieve internal" );
  New(0, VPrem, maxpi * 1024, char);
  memset(VPrem, 1, maxpi * 1024);
  for (pi = startpi+6; pi < maxpi; pi++)
    VPrem[pi*1024] = 0;
  for (pi = startpi+6, smallnc = 0; pi < maxpi; pi++) {
    uint32_t p = sprimes[pi];
    char* prem = VPrem + pi*1024;
    while (smallnc < nc && cl[smallnc] < p)   smallnc++;
    for (c = 1; c < smallnc; c++) prem[p-cl[c]] = 0;
    for (     ; c <      nc; c++) prem[p-(cl[c]%p)] = 0;
  }

  New(0, cres, nres, UV);

  rem_0 = mpz_fdiv_ui(low,pp_0);  remadd_0 = ppr % pp_0;
  rem_1 = mpz_fdiv_ui(low,pp_1);  remadd_1 = ppr % pp_1;
  rem_2 = mpz_fdiv_ui(low,pp_2);  remadd_2 = ppr % pp_2;

  /* Loop over their range in chunks of size 'ppr' */
  while (mpz_cmp(low, high) <= 0) {

    uint32_t r, nr, remr, ncres;
    unsigned long ui_low = (mpz_sizeinbase(low,2) > 8*sizeof(unsigned long)) ? 0 : mpz_get_ui(low);

    /* Reduce the allowed residues for this chunk using more primes */

    { /* Start making a list of this chunk's residues using three pairs */
      for (r = 0, ncres = 0; r < nres; r++) {
        addmodded(remr, rem_0, resmod_0[r], pp_0);
        if (crem_0[remr]) {
          addmodded(remr, rem_1, resmod_1[r], pp_1);
          if (crem_1[remr]) {
            addmodded(remr, rem_2, resmod_2[r], pp_2);
            if (crem_2[remr]) {
              cres[ncres++] = residues[r];
            }
          }
        }
      }
      addmodded(rem_0, rem_0, remadd_0, pp_0);
      addmodded(rem_1, rem_1, remadd_1, pp_1);
      addmodded(rem_2, rem_2, remadd_2, pp_2);
    }

    /* Sieve through more primes one at a time, removing residues. */
    for (pi = startpi+6; pi < maxpi && ncres > 0; pi++) {
      uint32_t p = sprimes[pi];
      uint32_t rem = (ui_low) ? (ui_low % p) : mpz_fdiv_ui(low,p);
      char* prem = VPrem + pi*1024;
      /* Check divisibility of each remaining residue with this p */
      if (startpi <= 9 || cres[ncres-1] < 4294967295U) {   /* Residues are 32-bit */
        for (r = 0, nr = 0; r < ncres; r++) {
          if (prem[ (rem+(uint32_t)cres[r]) % p ])
            cres[nr++] = cres[r];
        }
      } else {              /* Residues are 64-bit */
        for (r = 0, nr = 0; r < ncres; r++) {
          if (prem[ (rem+cres[r]) % p ])
            cres[nr++] = cres[r];
        }
      }
      ncres = nr;
    }
    if (_verbose > 2) printf("cluster sieve range has %u residues left\n", ncres);

    /* Now check each of the remaining residues for inclusion */
    for (r = 0; r < ncres; r++) {
      i = cres[r];
      mpz_add_ui(t, low, i);
      if (mpz_cmp(t, high) > 0) break;
      /* Pretest each element if the input is large enough */
      if (run_pretests) {
        for (c = 0; c < nc; c++)
          if (mpz_add_ui(t, low, i+cl[c]), mpz_gcd(t,t,_bgcd2), mpz_cmp_ui(t,1)) break;
        if (c != nc) continue;
      }
      /* PRP test */
      for (c = 0; c < nc; c++)
        if (! (mpz_add_ui(t, low, i+cl[c]), nprps++, _GMP_BPSW(t)) ) break;
      if (c != nc) continue;
      PUSH_VLIST(retlist, ibase + i);
    }
    ibase += ppr;
    mpz_add_ui(low, low, ppr);
  }

  if (_verbose) printf("cluster sieve ran %"UVuf" BPSW tests (pretests %s)\n", nprps, run_pretests ? "on" : "off");
  mpz_set(low, savelow);
  Safefree(cres);
  Safefree(VPrem);
  Safefree(resmod_0);
  Safefree(resmod_1);
  Safefree(resmod_2);
  Safefree(residues);
  mpz_clear(savelow);
  mpz_clear(t);
  *rn = retlist.nsize;
  return retlist.list;
}
