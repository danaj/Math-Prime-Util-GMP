#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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
  } else {
#if 0
    do {                        /* Trivial method. */
      mpz_isaac_urandomb(p, n);
      mpz_setbit(p, n-1);
      mpz_setbit(p, 0);
    } while (!_GMP_is_prob_prime(p));
#else
    mpz_t base;                 /* Fouque+Tibouchi Alg 1, without modulo checks */
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
#endif
  }
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
    croak("Failed to find %"UVuf" digit prime\n", n);

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

/*===========================================================================*/
/* Proven primes (Maurer and Shawe-Taylor */
/*===========================================================================*/

#define MAKE_PROOF_START(proofptr, n, nums) \
  if (proofptr) { \
    char* thisproof, *thisptr; \
    int prevlen = (*proofptr == 0) ? 0 : strlen(*proofptr); \
    int thislen = (5 + mpz_sizeinbase(n,10)) * nums + 200; \
    New(0, thisproof, thislen + prevlen + 1, char); \
    thisptr = thisproof; \
    thisptr += gmp_sprintf(thisptr,
#define MAKE_PROOF_END(proofptr) \
               ); \
    if (*proofptr) { \
      thisptr += gmp_sprintf(thisptr,"\n"); \
      strcat(thisptr, *proofptr); \
      Safefree(*proofptr); \
    } \
    *proofptr = thisproof; \
  }

#define USE_THEOREM5 0
void mpz_random_maurer_prime(mpz_t n, UV k, char** proofptr)
{
  mpz_t t, a, q, I, R;
  double m, r, minr = USE_THEOREM5 ? 0.334 : 0.5;
  int i, verbose = get_verbose_level();

  /* We could use safely use k <= 64. */
  if (k <= 32)
    return mpz_random_nbit_prime(n, k);

  r = minr;        /* size of q relative to size of n */
  m = 20;          /* always use at least this many bits of randomness */
  if (k > 2*m) {
    do {
      double s = ((double)isaac_rand32()) / ((double)4294967295.0);  /* [0,1] */
      r = pow(2,s-1);  /* exp2 is C99 */
#if USE_THEOREM5
      r = 0.334 + 1.332 * (r-0.5);  /* Stretch r to cover 0.334 - 1 */
#endif
    } while ((k-r*k) <= m);
  }
#if 0  /* Improve efficiency for less than ideal distribution */
  r -= 0.25;  if (r < minr) r = minr;
#endif
  mpz_init(t); mpz_init(a); mpz_init(q); mpz_init(I); mpz_init(R);
  mpz_random_maurer_prime(q, (UV)(r*k)+1, proofptr);
  mpz_setbit(I, k-1);
  mpz_mul_ui(t, q, 2);
  mpz_fdiv_q(I, I, t);   /* I = floor(2^(k-1) / 2q) */
  if (verbose && verbose != 3)
    { gmp_printf("r = %lf  k = %lu  q = %Zd  I = %Zd\n",r,k,q,I); fflush(stdout); }

  while (1) {
    if (verbose > 2) { printf("."); fflush(stdout); }
    mpz_isaac_urandomm(R, I);  /* [0, I-1]  */
    mpz_add(R, R, I);          /* [I, 2I-I] */
    mpz_add_ui(R, R, 1);       /* [I+1,2I]  */
#if USE_THEOREM5
    mpz_setbit(R, 0);          /* We need R to be odd */
#endif
    mpz_mul(n, R, q);
    mpz_mul_ui(n,n,2);
    mpz_add_ui(n,n,1);         /* n = 2Rq+1 */

    if (!primality_pretest(n)) continue;
    if (verbose > 2) { printf("+"); fflush(stdout); }
    /* if (!is_euler_plumb_pseudoprime(n)) continue; */
    if (!miller_rabin_ui(n,2)) continue;
    /* n is a base-2 psp and probably prime */

    if (verbose > 2) { printf("*"); fflush(stdout); }

    /* See if we can use BLS75 theorem 3 */
    mpz_mul_ui(t, q, 2);
    mpz_add_ui(t, t, 1);
    mpz_mul(t, t, t);
    if (mpz_cmp(t, n) > 0) {
      for (i = 0; i < 20; i++) {
        mpz_set_ui(a, pr[i]);
        /* Check A^R mod N != N-1 */
        mpz_powm(a, a, R, n);
        mpz_add_ui(t,a,1);
        if (mpz_cmp(t, n) == 0) continue;
        /* Check A^{Rq} mod N == N-1 */
        mpz_powm(a, a, q, n);
        mpz_add_ui(t,a,1);
        if (mpz_cmp(t, n) != 0) continue;
        if (verbose > 2) { printf("(%"UVuf")",k); fflush(stdout); }
        /* Ensure all results passed BPSW.  ~20% speed penalty. */
        if (!_GMP_is_lucas_pseudoprime(n,2)) croak("Maurer internal failure");
        MAKE_PROOF_START(proofptr, n, 3)
           "Type BLS3\nN %Zd\nQ %Zd\nA %u\n", n, q, pr[i]
        MAKE_PROOF_END(proofptr)
        mpz_clear(t); mpz_clear(a); mpz_clear(q); mpz_clear(I); mpz_clear(R);
        return;
      }
      /* Blast, we couldn't find the right 'a' value fast enough.  Try a new n. */
      continue;
    }

    /* Our q is smaller than sqrt(n)/2-1, so use BLS75 theorem 5. */
#if !USE_THEOREM5
    croak("random_maurer_prime: internal bit size error");
#else
    /* Check for obvious generation problems. */
    if (mpz_even_p(R)) continue;
    if (mpz_cmp_ui(R, 1) <= 0) continue;
    mpz_gcd(t, q, R);  if (mpz_cmp_ui(t, 1) != 0) continue;

    /* Theorem 5 with m = 2, assuming (I) which we'll check after this. */
    {
      mpz_t ts, tr, F;
      mpz_init(ts); mpz_init(tr); mpz_init(F);
      mpz_mul_ui(F, q, 2);
      /* Calculate r,s from page 624 of BLS75 */
      mpz_mul_ui(t, F, 2);
      mpz_tdiv_qr(ts, tr, R, t);
      /* Verify the r,s condition */
      mpz_mul(t, tr, tr);
      mpz_submul_ui(t, ts, 8);   /* t = r^2-8s */
      if (mpz_sgn(ts) != 0 && mpz_perfect_square_p(t)) {
        /* printf("fail r/s check\n"); */
        mpz_clear(ts);  mpz_clear(tr);  mpz_clear(F);
        continue;
      }
      /* Verify size of N with m=2. a,t are temps.  Should not fail. */
      mpz_mul(t, F, tr);
      mpz_add_ui(a, t, 1);         /* a = rF + 1 */
      mpz_sub_ui(tr, tr, 1);
      mpz_mul(t, F, F);
      mpz_mul_ui(t, t, 2);
      mpz_mul(t, t, tr);
      mpz_add(a, a, t);           /* a = (r-1)2FF + rF + 1 */
      mpz_mul(t, F, F);
      mpz_mul(t, t, F);
      mpz_mul_ui(t, t, 4);
      mpz_add(a, a, t);           /* a = 4FFF + (r-1)2FF + rF + 1 */
      mpz_add_ui(t, F, 1);
      mpz_clear(tr);  mpz_clear(ts);  mpz_clear(F);
      if (mpz_cmp(n,a) >= 0) {
        /* printf("fail N size check\n"); */
        continue;
      }
      /* Check divisibility required to use m=2 */
      if (mpz_divisible_p(n,t)) {
        /* printf("fail N divisiblity check\n"); */
        continue;
      }
    }

#define SET_A_CHECK_PSP(i) \
   mpz_set_ui(a, pr[i]); \
   if (apsp[i] == -1) \
     { mpz_sub_ui(t,n,1); mpz_powm(t,a,t,n); apsp[i] = (mpz_cmp_ui(t,1) == 0); } \
   if (apsp[i] == 0) continue;
#define CHECK_GCD(t) \
    mpz_powm(t, a, t, n); mpz_sub_ui(t, t, 1); mpz_gcd(t, t, n); \
    if (mpz_cmp_ui(t,1) != 0) continue;

    {
      int j, apsp[20];  /* apsp caches psp check.  Init all to -1. */
      for (i = 0; i < 20; i++) apsp[i] = -1;
      apsp[0] = 1;  /* We passed a base 2 psp test to get here */
      /* Find an a that works for p=2 */
      for (i = 0; i < 20; i++) {
        SET_A_CHECK_PSP(i);
        mpz_mul(t, q, R);
        CHECK_GCD(t);
        /* We are good for p=2.  Find an a for p=q. */
        for (j = 0; j < 20; j++) {
          SET_A_CHECK_PSP(j);
          mpz_mul_ui(t, R, 2);
          CHECK_GCD(t);
          /* Success */
          if (verbose > 2) { printf("(%lu)",k); fflush(stdout); }
          if (i == 0 && j == 0) {
            MAKE_PROOF_START(proofptr, n, 2)
              "Type BLS5\nN %Zd\nQ[1] %Zd\n----\n", n, q
            MAKE_PROOF_END(proofptr)
          } else {
            MAKE_PROOF_START(proofptr, n, 2)
              "Type BLS5\nN %Zd\nQ[1] %Zd\nA[0]  %lu\nA[1]  %lu\n----\n", n, q, pr[i], pr[j]
            MAKE_PROOF_END(proofptr)
          }
          /* Ensure all results passed BPSW.  ~20% speed penalty. */
          if (!_GMP_is_lucas_pseudoprime(n,2)) croak("Maurer internal failure");
          mpz_clear(t); mpz_clear(a); mpz_clear(q); mpz_clear(I); mpz_clear(R);
          return;
        }
        break; /* Failed for p=q */
      }
    }
    /* Blast, we couldn't find the right 'a' value fast enough.  Try a new n. */
#endif
  }
}

/* FIPS 186-4 algorithm but using our CSPRNG (ISAAC) instead of SHA-256 */
void mpz_random_shawe_taylor_prime(mpz_t c, UV k, char** proofptr)
{
  mpz_t c0, t, u, a, z;

  if (k <= 32)
    return mpz_random_nbit_prime(c, k);

  mpz_init(c0); mpz_init(t); mpz_init(u); mpz_init(a); mpz_init(z);

  mpz_random_shawe_taylor_prime(c0, 1 + (k+1)/2, proofptr);

  mpz_isaac_urandomb(t, k-1);
  mpz_setbit(t,k-1);       /* Steps 18-21: t a random k-bit integer */
  mpz_mul_ui(u, c0, 2);
  mpz_cdiv_q(t, t, u);     /* Step 22: set t based on random integer */

  while (1) {
    /* Steps 23-24 */
    mpz_mul_ui(u, c0, 2);  /* u = 2c0 */
    mpz_mul(c, u, t);
    mpz_add_ui(c, c, 1);   /* c = 2tc0+1 */
    if (mpz_sizeinbase(c,2) > k) {
      mpz_set_ui(t,0);
      mpz_setbit(t,k-1);
      mpz_cdiv_q(t, t, u);
      mpz_mul(c, u, t);
      mpz_add_ui(c, c, 1);
    }
    /* Don't bother with Steps 26-31 for obvious composites */
    if (primality_pretest(c) && miller_rabin_ui(c,2)) {
      /* Steps 26-29 */
      mpz_sub_ui(u, c, 3);
      mpz_isaac_urandomm(a, u);
      mpz_add_ui(a, a, 2);

      /* Step 30 */
      mpz_mul_ui(u, t, 2);
      mpz_powm(z, a, u, c);

      /* Step 31 */
      mpz_sub_ui(u, z, 1);
      mpz_gcd(u, u, c);
      if (mpz_cmp_ui(u, 1) == 0) {
        mpz_powm(u, z, c0, c);
        if (mpz_cmp_ui(u, 1) == 0) {
          /* Ensure all results passed BPSW.  ~20% speed penalty. */
          if (!_GMP_is_lucas_pseudoprime(c,2)) croak("ST internal failure");
          MAKE_PROOF_START(proofptr, c, 3)
             "Type Pocklington\nN %Zd\nQ %Zd\nA %Zd\n", c, c0, a
          MAKE_PROOF_END(proofptr)
          mpz_clear(c0); mpz_clear(t); mpz_clear(u); mpz_clear(a); mpz_clear(z);
          return;
        }
      }
    }
    mpz_add_ui(t,t,1);
  }
}

/*===========================================================================*/
