#include <gmp.h>
#include "ptypes.h"

#include "primality.h"
#include "lucas_seq.h"
#include "gmp_main.h"  /* primality_pretest */
#include "bls75.h"
#include "ecpp.h"
#include "factor.h"

#define FUNC_is_perfect_square 1
#define FUNC_mpz_logn
#include "utility.h"

#define NSMALLPRIMES 54
static const unsigned char sprimes[NSMALLPRIMES] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251};


/* 0 fail, 1 pass, -1 nothing.  Modifies base a */
static int _preprocess_base(mpz_t n, mpz_t a)
{
  if (mpz_cmp_ui(a, 1) <= 0)
    croak("Base %ld is invalid", mpz_get_si(a));
  if (mpz_cmp_ui(n, 2) <= 0)
    return (mpz_cmp_ui(n, 2) >= 0);

  if (mpz_cmp_ui(a, 2) > 0) {
    if (mpz_cmp(a, n) >= 0) {
      mpz_mod(a, a, n);
      if (mpz_cmp_ui(a, 1) <= 0)
        return mpz_sgn(a);
    }
  }
  return -1;
}

int is_pseudoprime(mpz_t n, mpz_t a)
{
  mpz_t nm1;
  int res;

  if ((res = _preprocess_base(n, a)) >= 0)
    return res;

  mpz_init(nm1);
  mpz_sub_ui(nm1, n, 1);
  mpz_powm(nm1, a, nm1, n);
  res = (mpz_cmp_ui(nm1, 1) == 0);
  mpz_clear(nm1);
  return res;
}

int is_euler_pseudoprime(mpz_t n, mpz_t a)
{
  mpz_t nm1, ap;
  int res;

  if (mpz_even_p(n))
    return (mpz_cmp_ui(n,2) == 0);
  if ((res = _preprocess_base(n, a)) >= 0)
    return res;

  mpz_init(ap);
  if (mpz_gcd(ap, a, n), mpz_cmp_ui(ap, 1) != 0) {
    mpz_clear(ap);
    return 0;
  }
  mpz_init(nm1);
  mpz_sub_ui(nm1, n, 1);
  mpz_tdiv_q_2exp(ap, nm1, 1);
  mpz_powm(ap, a, ap, n);

  if (mpz_cmp_ui(ap, 1) && mpz_cmp(ap, nm1))
    res = 0;
  else if (mpz_kronecker(a, n) >= 0)
    res = (mpz_cmp_ui(ap, 1) == 0);
  else
    res = (mpz_cmp(ap, nm1) == 0);

  mpz_clear(nm1);
  mpz_clear(ap);
  return res;
}

static int mrx(/*destroyed*/mpz_t x, /*destroyed*/ mpz_t d, mpz_t n, UV s)
{
  UV r;
  mpz_powm(x, x, d, n);
  mpz_sub_ui(d, n, 1);
  if (!mpz_cmp_ui(x, 1) || !mpz_cmp(x, d))
    return 1;
  for (r = 1; r < s; r++) {
    mpz_powm_ui(x, x, 2, n);
    if (!mpz_cmp_ui(x, 1))
      break;
    if (!mpz_cmp(x, d))
      return 1;
  }
  return 0;
}


int miller_rabin(mpz_t n, mpz_t a)
{
  mpz_t d, x;
  int cmpr, rval = 1;

  cmpr = mpz_cmp_ui(n, 2);
  if (cmpr == 0)     return 1;  /* 2 is prime */
  if (cmpr < 0)      return 0;  /* below 2 is not prime */
  if (mpz_even_p(n)) return 0;  /* multiple of 2 is composite */
  if (mpz_cmp_ui(a, 1) <= 0) croak("Base %ld is invalid", mpz_get_si(a));

  mpz_init_set(x, a);
  mpz_init_set(d, n);
  mpz_sub_ui(d, d, 1);

  /* Handle large and small bases.  Use x so we don't modify their input a. */
  if (mpz_cmp(x, n) >= 0)
    mpz_mod(x, x, n);
  if ( (mpz_cmp_ui(x, 1) > 0) && (mpz_cmp(x, d) < 0) ) {
    UV s = mpz_scan1(d, 0);
    mpz_tdiv_q_2exp(d, d, s);
    rval = mrx(x, d, n, s);
  }
  mpz_clear(d);
  mpz_clear(x);
  return rval;
}
int miller_rabin_ui(mpz_t n, unsigned long a)
{
  mpz_t d, x;
  int cmpr, rval = 1;

  cmpr = mpz_cmp_ui(n, 2);
  if (cmpr == 0)     return 1;  /* 2 is prime */
  if (cmpr < 0)      return 0;  /* below 2 is not prime */
  if (mpz_even_p(n)) return 0;  /* multiple of 2 is composite */
  if (a <= 1) croak("Base %lu is invalid", a);

  mpz_init_set_ui(x, a);
  mpz_init_set(d, n);
  mpz_sub_ui(d, d, 1);

  if (mpz_cmp(x, n) >= 0)
    mpz_mod(x, x, n);
  if ( (mpz_cmp_ui(x, 1) > 0) && (mpz_cmp(x, d) < 0) ) {
    UV s = mpz_scan1(d, 0);
    mpz_tdiv_q_2exp(d, d, s);
    rval = mrx(x, d, n, s);
  }
  mpz_clear(d);
  mpz_clear(x);
  return rval;
}

int is_miller_prime(mpz_t n, int assume_grh)
{
  mpz_t d, x, D;
  UV s, maxa, a;
  int rval;

  {
    int cmpr = mpz_cmp_ui(n, 2);
    if (cmpr == 0)     return 1;  /* 2 is prime */
    if (cmpr < 0)      return 0;  /* below 2 is not prime */
    if (mpz_even_p(n)) return 0;  /* multiple of 2 is composite */
  }

  if (mpz_cmp_ui(n, 1373653) < 0) {
    maxa = 3;
  } else if (assume_grh) {
    double logn = mpz_logn(n);
    double dmaxa = 2 * logn * logn;  /* Bach (1990) */
    /* Wedeniwski (2001) claims the following, but it might be wrong
     * double dmaxa = 1.5L * logn * logn - 44.0L/5.0L * logn + 13; */
    if (dmaxa >= (double)ULONG_MAX)
      croak("is_miller_prime: n is too large for GRH DMR");
    maxa = ceil(dmaxa);
  } else { /* Bober and Goldmakher 2015 (http://arxiv.org/abs/1311.7556) */
    /* n_p < p^(1/(4*sqrt(e))+epsilon).  Do it with logs */
    double dmaxa = exp( (1.0L/6.5948850828L) * mpz_logn(n) );
    if (dmaxa >= (double)ULONG_MAX)
      croak("is_miller_prime: n is too large for unconditional DMR");
    maxa = ceil(dmaxa);
  }

  if (mpz_cmp_ui(n, maxa) <= 0)
    maxa = mpz_get_ui(n) - 1;
  if (get_verbose_level() > 1)
    printf("Deterministic Miller-Rabin testing bases from 2 to %"UVuf"\n", maxa);

  mpz_init_set(d, n);
  mpz_sub_ui(d, d, 1);
  s = mpz_scan1(d, 0);
  mpz_tdiv_q_2exp(d, d, s);
  mpz_init(D);
  mpz_init(x);

  for (a = 2, rval = 1; rval && a <= maxa; a++) {
    mpz_set_ui(x, a);
    mpz_set(D, d);
    rval = mrx(x, D, n, s);
  }
  mpz_clear(x);
  mpz_clear(D);
  mpz_clear(d);
  return rval;
}

int miller_rabin_random(mpz_t n, UV numbases, char* seedstr)
{
  gmp_randstate_t randstate;
  mpz_t t, base;
  UV i;

  if (numbases == 0)  return 1;
  if (mpz_cmp_ui(n, 100) < 0)     /* tiny n */
    return (_GMP_is_prob_prime(n) > 0);

  /* See if they've asked for a ludicrous number of bases */
  mpz_init(t);
  mpz_mul_ui(t, n, 3);
  mpz_div_ui(t, t, 4);
  if (mpz_cmp_ui(t, numbases) <= 0) {
    int res = is_bpsw_dmr_prime(n);
    if (res != 1) {
      mpz_clear(t);
      return !!res;
    }
    numbases = mpz_get_ui(t);
  }

  mpz_init(base);
  mpz_sub_ui(t, n, 3);

  if (seedstr == 0) { /* Normal case with no seed string.  Use our CSPRNG. */
    for (i = 0; i < numbases; i++) {
      mpz_isaac_urandomm(base, t);    /* base 0 .. (n-3)-1 */
      mpz_add_ui(base, base, 2);      /* base 2 .. n-2     */
      if (miller_rabin(n, base) == 0)
        break;
    }
  } else {
    gmp_randinit_default(randstate);
    mpz_set_str(base, seedstr, 0);
    gmp_randseed(randstate, base);
    for (i = 0; i < numbases; i++) {
      mpz_urandomm(base, randstate, t); /* base 0 .. (n-3)-1 */
      mpz_add_ui(base, base, 2);        /* base 2 .. n-2     */
      if (miller_rabin(n, base) == 0)
        break;
    }
    gmp_randclear(randstate);
  }
  mpz_clear(base);  mpz_clear(t);
  return (i >= numbases);
}


int is_euler_plumb_pseudoprime(mpz_t n)
{
  unsigned int nmod8;
  mpz_t x, two;
  int result = 0;
  if (mpz_cmp_ui(n,5) < 0)
    return (mpz_cmp_ui(n,2) == 0 || mpz_cmp_ui(n,3) == 0);
  if (mpz_even_p(n))
    return 0;
  nmod8 = mpz_fdiv_ui(n, 8);
  mpz_init(x);  mpz_init_set_ui(two,2);
  mpz_sub_ui(x, n, 1);
  mpz_fdiv_q_2exp(x, x, 1 + (nmod8 == 1));
  mpz_powm(x, two, x, n);
  if (mpz_cmp_ui(x,1) == 0) {
    result = (nmod8 == 1 || nmod8 == 7);
  } else {
    mpz_add_ui(x,x,1);
    result = (mpz_cmp(x,n) == 0 && (nmod8 == 1 || nmod8 == 3 || nmod8 == 5));
  }
  mpz_clear(two); mpz_clear(x);
  return result;
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
  /* TODO: Should have a flag to skip pretests */
  if (mpz_cmp_ui(N,100) <= 0) return (_GMP_is_prob_prime(N) ? 2 : 0);
  if (mpz_even_p(N) || mpz_divisible_ui_p(N, 3)) return 0;
  mpz_init(v); mpz_init(k);
  mpz_sub_ui(v, N, 1);
  n = mpz_scan1(v, 0);
  mpz_tdiv_q_2exp(k, v, n);
  /* N = k * 2^n + 1 */
  if (mpz_sizeinbase(k,2) <= n) {
    mpz_init(a);
    /* Sze (2018) form without Jacobi tests */
    mpz_tdiv_q_2exp(k, v, 1);
    for (i = 0; i < 30 && res == -1; i++) {
      mpz_set_ui(a, sprimes[i]);
      mpz_powm(a, a, k, N);
      if (mpz_cmp_ui(a,1) != 0)
        res = (mpz_cmp(a, v) == 0)  ?  2  :  0;
    }
    mpz_clear(a);
  }
  /* TODO: look into Rao (2018): k*2^n+1 for n>1, k prime */
  /* if (n > 1 && res == -1 && _GMP_BPSW(k)) { ... } */
  mpz_clear(k); mpz_clear(v);
  if (res != -1 && get_verbose_level() > 1) {
    printf("N shown %s with Proth\n", res ? "prime" : "composite");
    fflush(stdout);
  }
  return res;
}
int is_proth_form(mpz_t N)
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
    mpz_set_uv(t, D);
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
    if (cmpr < 0)      return 0;  /* below 2 is not prime */
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
    if (cmpr < 0)      return 0;  /* below 2 is not prime */
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
  }
  mpz_clear(d);

  rval = 0;
  mpz_sub_ui(t, n, 2);

#if 0
  /* According to Stan Wagon's 1999 "Mathematica in Action", this method,
   * equivalent to our Extra Strong test, is used, though with a different
   * (unspecified) parameter selection method.  I'm not convinced Mathematica
   * actually uses this method, without any confirmation from Wolfram.
   */
  {
    /* V must be all 2 or x,x,-2,2,2,2.  First V=+/-2 must have a W=P. */
    int must_have_2 = 0, bad = 0;

    if (mpz_cmp_ui(V,2) == 0) {
      must_have_2 = 1;
      if (mpz_cmp_ui(W, P) != 0)
        bad = 1;
    }
    while (s-- && !bad) {
      if (must_have_2) {
        if (mpz_cmp_ui(V,2) != 0)
          bad = 1;
      } else {
        if (mpz_cmp(V,t) == 0) {  /* Found -2. Check W=-P.  Look for 2's. */
          mpz_sub(W, n, W);
          if (mpz_cmp_ui(W, P) != 0)
            bad = 1;
          must_have_2 = 1;
        } else {                  /* Keep looking for a -2. */
          mpz_mul(W, V, W);
          mpz_sub_ui(W, W, P);
          mpz_mod(W, W, n);
        }
      }
      mpz_mul(V, V, V);
      mpz_sub_ui(V, V, 2);
      mpz_mod(V, V, n);
    }
    /* Fail if bad val found or didn't have a 2 at the end */
    rval = !bad && must_have_2 && !mpz_cmp_ui(V,2);
  }
#else
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
#endif
  mpz_clear(W); mpz_clear(V); mpz_clear(t);
  return rval;
}

int is_perrin_pseudoprime(mpz_t n, int restricted)
{
  mpz_t S[6], T[6], T01, T34, T45, t;
  int cmpr, i, j, rval;

  cmpr = mpz_cmp_ui(n, 2);
  if (cmpr == 0)     return 1;  /* 2 is prime */
  if (cmpr < 0)      return 0;  /* below 2 is not prime */
  if (restricted > 2 && mpz_even_p(n)) return 0;

  { /* Simple filter for composites */
    uint32_t n32 = mpz_fdiv_ui(n, 2762760);
    if (!(n32& 1) && !((   22 >> (n32% 7)) & 1)) return 0;
    if (!(n32% 3) && !((  523 >> (n32%13)) & 1)) return 0;
    if (!(n32% 5) && !((65890 >> (n32%24)) & 1)) return 0;
    if (!(n32% 4) && !((  514 >> (n32%14)) & 1)) return 0;
    if (!(n32%23) && !((    2 >> (n32%22)) & 1)) return 0;
  }

  /* Calculate signature using Adams/Shanks doubling rule. */
  mpz_init(t);
  mpz_init_set_ui(S[0], 1);
  mpz_init(S[1]); mpz_sub_ui(S[1], n, 1);
  mpz_init_set_ui(S[2], 3);
  mpz_init_set_ui(S[3], 3);
  mpz_init_set_ui(S[4], 0);
  mpz_init_set_ui(S[5], 2);

  for (i=0; i < 6; i++)
    mpz_init(T[i]);
  mpz_init(T01); mpz_init(T34); mpz_init(T45);

  for (i = mpz_sizeinbase(n,2)-2; i >= 0; i--) {
    mpz_mul(t,S[0],S[0]); mpz_sub(t,t,S[5]); mpz_sub(t,t,S[5]); mpz_mod(T[0],t,n);
    mpz_mul(t,S[1],S[1]); mpz_sub(t,t,S[4]); mpz_sub(t,t,S[4]); mpz_mod(T[1],t,n);
    mpz_mul(t,S[2],S[2]); mpz_sub(t,t,S[3]); mpz_sub(t,t,S[3]); mpz_mod(T[2],t,n);
    mpz_mul(t,S[3],S[3]); mpz_sub(t,t,S[2]); mpz_sub(t,t,S[2]); mpz_mod(T[3],t,n);
    mpz_mul(t,S[4],S[4]); mpz_sub(t,t,S[1]); mpz_sub(t,t,S[1]); mpz_mod(T[4],t,n);
    mpz_mul(t,S[5],S[5]); mpz_sub(t,t,S[0]); mpz_sub(t,t,S[0]); mpz_mod(T[5],t,n);
    mpz_sub(t,T[2],T[1]); mpz_mod(T01,t,n);
    mpz_sub(t,T[5],T[4]); mpz_mod(T34,t,n);
    mpz_add(t,T34, T[3]); mpz_mod(T45,t,n);
    if (mpz_tstbit(n, i)) {
      mpz_set(S[0],T[0]); mpz_set(S[1],T01); mpz_set(S[2],T[1]);
      mpz_set(S[3],T[4]); mpz_set(S[4],T45); mpz_set(S[5],T[5]);
    } else {
      mpz_add(t,T01,T[0]);
      mpz_set(S[0],T01); mpz_set(S[1],T[1]); mpz_mod(S[2],t,n);
      mpz_set(S[3],T34); mpz_set(S[4],T[4]); mpz_set(S[5],T45);
    }
  }

  for (i=0; i < 6; i++)
    mpz_clear(T[i]);
  mpz_clear(T01); mpz_clear(T34); mpz_clear(T45);

  rval = !mpz_sgn(S[4]);
  if (rval == 0 || restricted == 0)
    goto DONE_PERRIN;

  mpz_sub_ui(t,n,1);
  rval = !mpz_cmp(S[1],t);
  if (rval == 0 || restricted == 1)
    goto DONE_PERRIN;

  /* Adams/Shanks or Arno,Grantham full signature test */
  rval = 0;
  j = mpz_si_kronecker(-23, n);

  if (j == -1) {
    mpz_t A, B, C;
    mpz_init_set(B, S[2]); mpz_init(A); mpz_init(C);

    mpz_mul(t,B,B); mpz_mod(t,t,n);
    mpz_mul_ui(A,B,3); mpz_add_ui(A,A,1); mpz_sub(A,A,t); mpz_mod(A,A,n);
    mpz_mul_ui(C,t,3); mpz_sub_ui(C,C,2); mpz_mod(C,C,n);

    mpz_mul(t,t,B); mpz_sub(t,t,B); mpz_mod(t,t,n);
    rval = !mpz_cmp(S[0],A) && !mpz_cmp(S[2],B) && !mpz_cmp(S[3],B) && !mpz_cmp(S[5],C) && mpz_cmp_ui(B,3) && !mpz_cmp_ui(t,1);
  } else if (restricted > 2 && j == 0 && mpz_cmp_ui(n,23)) {
    /* Jacobi symbol says 23|n.  n is composite if != 23 */
    rval = 0;
  } else {
    if (!mpz_cmp(S[2],S[3])) {
      rval = !mpz_cmp_ui(S[0],1) && !mpz_cmp_ui(S[2],3) && !mpz_cmp_ui(S[3],3) && !mpz_cmp_ui(S[5],2);
    } else {
      mpz_sub_ui(t, n, 1);
      rval = !mpz_cmp_ui(S[0],0) && !mpz_cmp(S[5],t) &&
             (mpz_add(t,S[2],S[3]),mpz_add_ui(t,t,3),mpz_mod(t,t,n),!mpz_sgn(t)) &&
             (mpz_sub(t,S[2],S[3]),mpz_mul(t,t,t),mpz_add_ui(t,t,23),mpz_mod(t,t,n),!mpz_sgn(t));
    }
  }

DONE_PERRIN:
  for (i=0; i < 6; i++)
    mpz_clear(S[i]);
  mpz_clear(t);
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
    if (cmpr < 0)      return 0;  /* below 2 is not prime */
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
    mpz_set_iv(Vcomp, Q);
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
      mpz_isaac_urandomm(a, t);
      mpz_add_ui(a, a, 1);
      mpz_isaac_urandomm(b, t);
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
    if (cmpr < 0)      return 0;  /* below 2 is not prime */
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
  if (_verbose>1) { gmp_printf("%Zd is %s with a = %"UVuf"\n", n, (rval) ? "probably prime" : "composite", a); fflush(stdout); }

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
    if (cmpr < 0)      return 0;  /* below 2 is not prime */
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




int _GMP_BPSW(mpz_t n)
{
  if (mpz_cmp_ui(n, 4) < 0)
    return (mpz_cmp_ui(n, 1) <= 0) ? 0 : 2;

  if (miller_rabin_ui(n, 2) == 0)   /* Miller Rabin with base 2 */
    return 0;

  if (_GMP_is_lucas_pseudoprime(n, 2 /*extra strong*/) == 0)
    return 0;

  if (mpz_sizeinbase(n, 2) <= 64)        /* BPSW is deterministic below 2^64 */
    return 2;

  return 1;
}

/* Assume n is a BPSW PRP, return 1 (no result), 0 (composite), 2 (prime) */
int is_deterministic_miller_rabin_prime(mpz_t n)
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
        res = miller_rabin_ui(n, sprimes[i]);
      }
      if (res == 1) res = 2;
    }
    mpz_clear(t);
  }
  return res;
}


/*
 * is_prob_prime      BPSW -- fast, no known counterexamples
 * is_prime           is_prob_prime + a little extra
 * is_provable_prime  really prove it, which could take a very long time
 *
 * They're all identical for numbers <= 2^64.
 *
 * The extra M-R test in is_prime start actually costing something after
 * 1000 bits or so.  Primality proving will get quite a bit slower as the
 * number of bits increases.
 *
 * Both is_prime and is_provable prime start with some trial division
 * followed by an extra strong BPSW test, just like is_prob_prime.  For
 * 65+ bit inputs, they do more.  A single extra random-base M-R test is
 * next done.  In the worst case this has a 1/4 chance of rejecting the
 * composite, but an average change of < 1e-12.  is_prime will return at
 * this point.  is_provable_prime tries various special proofs, followed
 * by ECPP.  ECPP returns either:
 *    2  "definitely prime",
 *    0  "definitely composite, or
 *    1  "probably prime"
 * where the latter doesn't really tell us anything.  It's unusual, but
 * possible.  If this even happs, we do a couple Frobenius-type tests, so
 * even an answer of 1 (not proven) will have less chance of false results.
 * A composite would have to have no small factors, be the first known
 * counterexample to each of three different and distinct tests, pass a
 * random-base M-R test, and manage to have the ECPP implementation not
 * find it a composite.
 */

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

int is_bpsw_dmr_prime(mpz_t n)
{
  int prob_prime = _GMP_BPSW(n);
  if (prob_prime == 1) {
    prob_prime = is_deterministic_miller_rabin_prime(n);
    if (prob_prime == 0) gmp_printf("\n\n**** BPSW counter-example found?  ****\n**** N = %Zd ****\n\n", n);
  }
  return prob_prime;
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

  /* For small numbers, try a quick BLS75 proof. */
  if (prob_prime == 1) {
    if (is_proth_form(n))
      prob_prime = BLS_primality_nm1(n, 2 /* effort */, 0 /* cert */);
    else if (nbits <= 150)
      prob_prime = BLS_primality_nm1(n, 0 /* effort */, 0 /* cert */);
    /* We could do far better by calling bls75_hybrid, especially with a
     * larger effort.  But that takes more time.  I've decided it isn't worth
     * the extra time here.  We'll still try a very quick N-1 proof. */
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
   *
   * We can estimate probabilities of random odd 'k'-bit composites passing
   * 't' random-base Miller-Rabin tests using papers such as Kim and
   * Pomerance; Damg??rd, Landrock, and Pomerance; Lichtman and Pomerance.
   * We can also perform random sampling to get empirical data, which shows
   * a probability of less than 1e-12 at 65 bits, and lowering as the number
   * of bits increases.  One extra test is all that is necessary.
   */

  if (prob_prime == 1) {
    UV ntests = 1;
    prob_prime = miller_rabin_random(n, ntests, 0);
    /* prob_prime = _GMP_is_frobenius_underwood_pseudoprime(n); */
    /* prob_prime = _GMP_is_frobenius_khashin_pseudoprime(n); */
  }

  /* We have now done trial division, an SPSP-2 test, an extra-strong Lucas
   * test, and a random-base Miller-Rabin test.  We have exceeded the testing
   * done in most math packages.  Any composites will have no small factors,
   * be the first known BPSW counterexample, and gotten past the random-base
   * Miller-Rabin test.
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
  prob_prime = miller_rabin_random(n, 1, 0);
  if (prob_prime != 1)  return prob_prime;

  /* We can choose a primality proving algorithm:
   *   AKS      _GMP_is_aks_prime       really slow, don't bother
   *   N-1      BLS_primality_nm1       small or special numbers
   *   N+1      BLS_primality_np1       small or special numbers
   *   N-1/N+1  BLS_primality           small or special numbers
   *   ECPP     _GMP_ecpp               fastest in general
   */

  if (prooftext) {
    /* Give n-1 a small go */
    prob_prime = BLS_primality_nm1(n, is_proth_form(n) ? 3 : 1, prooftext);
    if (prob_prime != 1)  return prob_prime;
  } else {
    /* See if there's an easy n-1 / n+1 hybrid proof */
    prob_prime = BLS_primality(n, is_proth_form(n) ? 3 : 1, prooftext);
    if (prob_prime != 1)  return prob_prime;
  }

  /* ECPP */
  prob_prime = _GMP_ecpp(n, prooftext);

  /* While extremely unusual, it is not impossible for our ECPP implementation
   * to give up.  If this happens, we won't return 2 with a proof, but let's
   * at least run some more tests. */
  if (prob_prime == 1)
    prob_prime = _GMP_is_frobenius_underwood_pseudoprime(n);
  if (prob_prime == 1)
    prob_prime = _GMP_is_frobenius_khashin_pseudoprime(n);

  return prob_prime;
}
