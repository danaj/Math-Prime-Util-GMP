#include <string.h>
#include <math.h>
#include <float.h>
#include <gmp.h>
#include "ptypes.h"

#include "real.h"
#include "prime_iterator.h"
#include "primality.h"
#include "factor.h"
#define FUNC_mpz_logn 1
#include "utility.h"

static unsigned long precbits(mpf_t x, unsigned long prec, unsigned long add) {
  unsigned long bits1 = mpf_get_prec(x), bits2 = DIGS2BITS(prec);
  return add + ((bits1 > bits2) ? bits1 : bits2);
}

/*****************************************************************************/

/* Put result into char with correct number of digits */
static char* _str_real(mpf_t f, unsigned long prec) {
  char* out;
  unsigned long k;
  int neg = (mpf_sgn(f) < 0);

  if (neg)
    mpf_neg(f, f);

  for (k = 0;  mpf_cmp_ui(f, 1000000000U) >= 0;  k += 9)
    mpf_div_ui(f, f, 1000000000U);
  for (;  mpf_cmp_ui(f, 1) >= 0;  k++)
    mpf_div_ui(f, f, 10);

  New(0, out, 10+((k>prec) ? k : prec), char);
  gmp_sprintf(out, "%.*Ff", prec, f);
  if (out[0] == '0') {
    memmove(out, out+2, prec);
  } else { /* We rounded up.  Treat like 0.1 with one larger k */
    memmove(out+1, out+2, prec);
    k++;
  }

  if (k >= prec) { /* No decimal */
    if (k-prec < 10) {
      memset(out+prec, '0', k-prec);
      prec=k-1;
    } else {
      out[prec++] = 'E';
      prec += sprintf(out+prec, "%lu", k-prec+1);
    }
  } else {        /* insert decimal in correct place */
    memmove(out+k+1, out+k, prec-k);
    out[k] = '.';
  }
  out[prec+1]='\0';
  if (neg) {
    memmove(out+1, out, prec+2);
    out[0] = '-';
  }
  return out;
}

static char* _frac_real(mpz_t num, mpz_t den, unsigned long prec) {
#if 0
  char* out;
  mpf_t fnum, fden, res;
  unsigned long numbits = mpz_sizeinbase(num,  2);
  unsigned long denbits = mpz_sizeinbase(den,  2);
  unsigned long numdigs = mpz_sizeinbase(num, 10);
  unsigned long dendigs = mpz_sizeinbase(den, 10);

  mpf_init2(fnum, 1 + numbits);  mpf_set_z(fnum, num);
  mpf_init2(fden, 1 + denbits);  mpf_set_z(fden, den);
  mpf_init2(res, (unsigned long) (8 + (numbits-denbits+1) + prec*3.4) );
  mpf_div(res, fnum, fden);
  mpf_clear(fnum);  mpf_clear(fden);

  New(0, out, (10+numdigs-dendigs)+prec, char);
  gmp_sprintf(out, "%.*Ff", (int)(prec), res);
  mpf_clear(res);

  return out;
#else
  char* out;
  mpf_t fnum, fden;
  unsigned long bits = 32+(unsigned long)(prec*3.32193);
  mpf_init2(fnum, bits);   mpf_set_z(fnum, num);
  mpf_init2(fden, bits);   mpf_set_z(fden, den);
  mpf_div(fnum, fnum, fden);
  out = _str_real(fnum, prec);
  mpf_clear(fden);
  mpf_clear(fnum);
  return out;
#endif
}


/*********************     Riemann Zeta and Riemann R     *********************/

static void _bern_real_zeta(mpf_t bn, unsigned long s, unsigned long prec);
static unsigned long zeta_n = 0;
static mpz_t* zeta_d = 0;

void free_borwein_zeta(void) {
  unsigned long i;
  if (zeta_n > 0) {
    for (i = 0; i <= zeta_n; i++)
      mpz_clear(zeta_d[i]);
    Safefree(zeta_d);
    zeta_n = 0;
  }
}

static void _borwein_d(unsigned long D) {
  mpz_t t1, t2, t3, sum;
  unsigned long i, n = 3 + (1.31 * D);

  if (zeta_n >= n)
    return;

  free_borwein_zeta();

  n += 10;   /* Add some in case we want a few more digits later */
  zeta_n = n;
  New(0, zeta_d, n+1, mpz_t);
  mpz_init(t1); mpz_init(t2); mpz_init(t3);

  mpz_init_set_ui(sum, 1);
  mpz_init_set(zeta_d[0], sum);

  mpz_fac_ui(t1, n);
  mpz_fac_ui(t2, n);
  for (i = 1; i <= n; i++) {
    mpz_mul_ui(t1, t1, 2*(n+i-1));    /* We've pulled out a 2 from t1 and t2 */
    mpz_divexact_ui(t2, t2, n-i+1);
    mpz_mul_ui(t2, t2, (2*i-1) * i);
    mpz_divexact(t3, t1, t2);
    mpz_add(sum, sum, t3);
    mpz_init_set(zeta_d[i], sum);
  }
  mpz_clear(sum); mpz_clear(t3); mpz_clear(t2); mpz_clear(t1);
}

/* MPFR does some shortcuts, then does an in-place version of Borwein 1991.
 * It's quite clever, and has the advantage of not using the statics.  Our
 * code can be a little faster in some cases, slower in others.  They
 * certainly have done more rigorous error bounding, allowing fewer guard
 * bits and earlier loop exits.
 *
 * For real values, we use our home-grown mpf_pow function, which is slower
 * at high precisions compared to MPFR or Pari.
 */

static void _zeta(mpf_t z, mpf_t f, unsigned long prec)
{
  unsigned long k, S, p;
  mpf_t s, tf, term;
  mpz_t t1;

  if (mpf_cmp_ui(f,1) == 0) {
   mpf_set_ui(z, 0);
   return;
 }

  /* Shortcut if we know all prec terms are zeros. */
  if (mpf_cmp_ui(f, 1+3.3219281*prec) >= 0 || mpf_cmp_ui(f, mpf_get_prec(z)) > 0) {
    mpf_set_ui(z,1);
    return;
  }

  S = (mpf_integer_p(f) && mpf_fits_ulong_p(f))  ?  mpf_get_ui(f)  :  0;

  /* Negative integers using Bernoulli */
  if (S == 0 && mpf_integer_p(f) && mpf_fits_slong_p(f) && mpf_sgn(f) != 0) {
    S = -mpf_get_si(f);
    if (!(S & 1)) { /* negative even integers are zero */
      mpf_set_ui(z,0);
    } else {        /* negative odd integers are -B_(n+1)/(n+1) */
      _bern_real_zeta(z, S+1, prec);
      mpf_div_ui(z, z, S+1);
      mpf_neg(z,z);
    }
    return;
  }

  mpf_init2(s,    96+mpf_get_prec(z));   mpf_set(s, f);
  mpf_init2(tf,   96+mpf_get_prec(z));
  mpf_init2(term, 96+mpf_get_prec(z));
  mpz_init(t1);

  if (S && S <= 14 && !(S & 1)) {         /* Small even S can be done with Pi */
    unsigned long div[]={0,6,90,945,9450,93555,638512875,18243225};
    const_pi(z, prec);
    mpf_pow_ui(z, z, S);
    if (S == 12) mpf_mul_ui(z, z, 691);
    if (S == 14) mpf_mul_ui(z, z, 2);
    mpf_div_ui(z, z, div[S/2]);
  } else if (mpf_cmp_ui(f, 3+prec*2.15) > 0) {  /* Only one term (3^s < prec) */
    if (S) {
      mpf_set_ui(term, 1);
      mpf_mul_2exp(term, term, S);
    } else {
      mpf_set_ui(term, 2);
      mpf_pow(term, term, s);
    }
    mpf_sub_ui(tf, term, 1);
    mpf_div(z, term, tf);
  } else if ( (mpf_cmp_ui(f,20) > 0 && mpf_cmp_ui(f, prec/3.5) > 0) ||
              (prec > 500 && (mpz_ui_pow_ui(t1, 8*prec, S), mpz_sizeinbase(t1,2) > (20+3.3219281*prec))) ) {
    /* Basic formula, for speed (also note only valid for > 1) */
    PRIME_ITERATOR(iter);
    mpf_set_ui(z, 1);
    for (p = 2; p <= 1000000000; p = prime_iterator_next(&iter)) {
      if (S) {
        mpz_ui_pow_ui(t1, p, S);
        mpf_set_z(term, t1);
      } else {
        mpf_set_ui(tf, p);
        mpf_pow(term, tf, s);
        mpz_set_f(t1, term);
      }
      if (mpz_sizeinbase(t1,2) > (20+3.3219281*prec)) break;
      mpf_sub_ui(tf, term, 1);
      mpf_div(term, term, tf);
      mpf_mul(z, z, term);
    }
    prime_iterator_destroy(&iter);
  } else {
    /* TODO: negative non-integer inputs past -20 or so are very wrong. */
    _borwein_d( (mpf_cmp_d(f,-3.0) >= 0)  ?  prec  :  80+2*prec );

    mpf_set_ui(z, 0);
    for (k = 0; k <= zeta_n-1; k++) {
      if (S) {
        mpz_ui_pow_ui(t1, k+1, S);
        mpf_set_z(term, t1);
      } else {
        mpf_set_ui(tf, k+1);
        mpf_pow(term, tf, s);
      }

      mpz_sub(t1, zeta_d[k], zeta_d[zeta_n]);
      mpf_set_z(tf, t1);

      mpf_div(term, tf, term);

      if (k&1) mpf_sub(z, z, term);
      else     mpf_add(z, z, term);
    }

    mpf_set_z(tf, zeta_d[zeta_n]);
    mpf_div(z, z, tf);

    if (S) {
      mpf_set_ui(tf, 1);
      mpf_div_2exp(tf, tf, S-1);
    } else {
      mpf_set_ui(term, 2);
      mpf_ui_sub(tf, 1, s);
      mpf_pow(tf, term, tf);
    }

    mpf_ui_sub(tf, 1, tf);
    mpf_div(z, z, tf);

    mpf_neg(z, z);
  }
  mpz_clear(t1);
  mpf_clear(term); mpf_clear(tf); mpf_clear(s);
}

static void _zetaint(mpf_t z, unsigned long s, unsigned long prec)
{
  mpf_t f;

  if (s <= 1) {
    mpf_set_ui(z, 0);
  } else if (s >= (1+3.3219281*prec) || s > mpf_get_prec(z)) {
    /* Shortcut if we know all prec terms are zeros. */
    mpf_set_ui(z,1);
  } else {
    mpf_init2(f, mpf_get_prec(z));
    mpf_set_ui(f, s);
    _zeta(z, f, prec);
    mpf_clear(f);
  }
}

static void _riemann_r(mpf_t r, mpf_t n, unsigned long prec)
{
  mpf_t logn, sum, term, part_term, tol, tf;
  unsigned long k, bits = precbits(r, prec, 7);

  mpf_init2(logn,      bits);
  mpf_init2(sum,       bits);
  mpf_init2(term,      bits);
  mpf_init2(part_term, bits);
  mpf_init2(tol,       bits);
  mpf_init2(tf,        bits);

  mpf_log(logn, n);
  mpf_set_ui(tol, 10);  mpf_pow_ui(tol, tol, prec);  mpf_ui_div(tol,1,tol);

#if 1 /* Standard Gram Series */
  mpf_set_ui(part_term, 1);
  mpf_set_ui(sum, 1);
  for (k = 1; k < 10000; k++) {
    mpf_mul(part_term, part_term, logn);
    mpf_div_ui(part_term, part_term, k);

    _zetaint(tf, k+1, prec+1);
    mpf_mul_ui(tf, tf, k);
    mpf_div(term, part_term, tf);
    mpf_add(sum, sum, term);

    mpf_abs(term, term);
    mpf_mul(tf, sum, tol);
    if (mpf_cmp(term, tf) <= 0) break;
  }
#else  /* Accelerated (about half the number of terms needed) */
  /* See:   http://mathworld.wolfram.com/GramSeries.html (5)
   * Ramanujan's G (3) and its restatements (5) and (6) are equal,
   * but G is only asymptotically equal to R (restated in (4) as per Gram).
   * To avoid confusion we won't use this.  Too bad, as it can be 2x faster.
   */
  mpf_set(part_term, logn);
  mpf_set_ui(sum, 0);
  _zetaint(tf, 2, prec);
  mpf_div(term, part_term, tf);
  mpf_add(sum, sum, term);
  for (k = 2; k < 1000000; k++) {
    if (mpf_cmp(part_term, tol) <= 0) break;

    mpf_mul(tf, logn, logn);
    if (k < 32768) {
      mpf_div_ui(tf, tf, (2*k-2) * (2*k-1));
    } else {
      mpf_div_ui(tf, tf, 2*k-2);  mpf_div_ui(tf, tf, 2*k-1);
    }
    mpf_mul(part_term, part_term, tf);

    _zetaint(tf, 2*k, prec);
    mpf_mul_ui(tf, tf, 2*k-1);
    mpf_div(term, part_term, tf);
    mpf_add(sum, sum, term);
  }
  mpf_mul_ui(sum, sum, 2);
  mpf_add_ui(sum, sum, 1);
#endif

  mpf_set(r, sum);

  mpf_clear(tf); mpf_clear(tol); mpf_clear(part_term);
  mpf_clear(term); mpf_clear(sum); mpf_clear(logn);
}

/***********************     Constants: Euler, Pi      ***********************/

/* See:
 *   http://numbers.computation.free.fr/Constants/Gamma/gamma.pdf
 *   https://www.ginac.de/CLN/binsplit.pdf (3.1)
 *   https://www-fourier.ujf-grenoble.fr/~demailly/manuscripts/gamma_gazmath_eng.pdf
 *   Pari/GP trans1.c
 *
 * Mortici and Chen (2013) have a O(n^-12) method, but it still too slow.
 * https://link.springer.com/content/pdf/10.1186/1029-242X-2013-222.pdf
 *
 * The Stieltjes zeta method isn't terrible but too slow for large n.
 *
 * We use the same series method as Pari, running about 3x faster.
 * We should use binary splitting as it still isn't really fast.
 * Each doubling of digits increases time by 4x.
 */
static void _const_euler(mpf_t gamma, unsigned long prec)
{
  const double log10 = 2.3025850929940456840179914546843642076L;
  const unsigned long maxsqr = (1UL << (4*sizeof(unsigned long))) - 1;
  unsigned long bits = 40 + DIGS2BITS(prec);
  unsigned long x = floor(2 + prec * log10/4);
  unsigned long N = ceil(1 + 3.591121477*x - 0.195547*log(x));
  unsigned long xx = x*x;
  unsigned long k;
  mpf_t u, v, a, b, fxx;

  if (prec <= 100) {
    mpf_set_str(gamma, "0.5772156649015328606065120900824024310421593359399235988057672348848677267776646709369470632917467495", 10);
    return;
  }

  mpf_init2(u,    bits);
  mpf_init2(v,    bits);
  mpf_init2(a,    bits);
  mpf_init2(b,    bits);

  /*
   * Brent and McMillan (1980) algorithm B1.
   * http://www.ams.org/journals/mcom/1980-34-149/S0025-5718-1980-0551307-4/
   */
  mpf_set_ui(u, x);
  mpf_log(u, u);
  mpf_neg(u, u);
  mpf_set(a, u);
  mpf_set_ui(b, 1);
  mpf_set_ui(v, 1);

  if (x <= maxsqr && N <= maxsqr) {
    for (k = 1; k <= N; k++) {
      mpf_mul_ui(b, b, xx);
      mpf_div_ui(b, b, k*k);   /* B_k = B_{k-1} * x^2 / k^2 */
      mpf_mul_ui(a, a, xx);
      mpf_div_ui(a, a, k);
      mpf_add(a, a, b);
      mpf_div_ui(a, a, k);     /* A_k = (A_{k-1} * x^2 / k + B_k) / k */
      mpf_add(u, u, a);
      mpf_add(v, v, b);
    }
  } else {
    mpf_init2(fxx,bits);
    mpf_set_ui(fxx, x);
    mpf_mul(fxx, fxx, fxx);
    for (k = 1; k <= N; k++) {
      mpf_mul(b,b,fxx);
      if (k <= maxsqr) { mpf_div_ui(b,b,k*k); }
      else             { mpf_div_ui(b,b,k);  mpf_div_ui(b,b,k); }
      mpf_mul(a,a,fxx);
      mpf_div_ui(a, a, k);
      mpf_add(a, a, b);
      mpf_div_ui(a, a, k);
      mpf_add(u, u, a);
      mpf_add(v, v, b);
    }
    mpf_clear(fxx);
  }
  mpf_div(gamma, u, v);
  mpf_clear(u); mpf_clear(v); mpf_clear(a); mpf_clear(b);
}

/* There are a plethora of interesting ways to calculate Pi.
 *
 * - Spigot.  Far too slow, though nice in plain C for a few hundred digits.
 *
 * - Machin-like, using Machin, Störmer, Chien-lih, Arndt, etc.
 *   See http://www.jjj.de/arctan/arctanpage.html for best arctan series.
 *
 * - AGM.  Quite good, and this seems to be best for relatively small sizes.
 *   See: https://arxiv.org/abs/1802.07558
 *
 * - Ramanujan / Chudnovsky with binary splitting.
 *   About 2-4x faster than AGM for large enough sizes.  This version is based
 *   on Alexander Yee's example.  I have tested vs a port of Pari/GP's abpq_sum
 *   and it came out about the same speed but uses a lot more memory.
 *   There are many more optimizations that can be done for this.  Xue's code
 *   on the GMP page uses quite a bit of code to do running reduction of P,Q
 *   which makes it about 1.5x faster.
 */

static void _set_pqr(mpz_t P, mpz_t Q, mpz_t R, unsigned long b)
{
  if (sizeof(unsigned long) < 8 || b > 630000) {
    mpz_set_ui(P, b);
    mpz_mul(Q, P, P);
    mpz_mul_ui(R, P, 26726400UL);
    mpz_mul_ui(R, R, 409297880UL);
    mpz_mul(Q, Q, R);
    mpz_set_ui(R,    2*b-1);
    mpz_mul_ui(R, R, 6*b-5);
    mpz_mul_ui(R, R, 6*b-1);
    mpz_mul_ui(P, P, 545140134UL);
    mpz_add_ui(P, P, 13591409UL);
  } else {
    mpz_set_ui(Q, b*b*b);
    mpz_mul_ui(Q, Q, 26726400UL*409297880UL); /* 10939058860032000UL */
    mpz_set_ui(R, (2*b-1) * (6*b-5) * (6*b-1));
    mpz_set_ui(P, b*545140134UL);
    mpz_add_ui(P, P, 13591409UL);
  }
  mpz_mul(P, P, R);
  if (b % 2 == 1) mpz_neg(P, P);
}

static void _sum_pqr(mpz_t P, mpz_t Q, mpz_t R, mpz_t u, unsigned long a, unsigned long b)
{
  if (b-a == 1) {
    _set_pqr(P, Q, R, b);
  } else {
    mpz_t P1, Q1, R1;
    unsigned long m;

    mpz_init(P1); mpz_init(Q1); mpz_init(R1);

    if (b-a == 2) {
      _set_pqr(P, Q, R, b-1);
      _set_pqr(P1, Q1, R1, b);
    } else if (b-a == 3) {
      m = a+2;
      _sum_pqr(P, Q, R, u, a, m);
      _set_pqr(P1, Q1, R1, b);
    } else {
       m = a + (b-a)*0.54;   /* Biased splitting */
      _sum_pqr(P, Q, R, u, a, m);
      _sum_pqr(P1, Q1, R1, u, m, b);
    }

    /* P = P0*Q1+P1*R0  Q = Q0*Q1  R = R0*R1 */
    mpz_mul(u, P1, R);
    mpz_mul(P, P, Q1);  mpz_add(P, P, u);
    mpz_mul(Q, Q, Q1);
    mpz_mul(R, R, R1);

    mpz_clear(P1); mpz_clear(Q1); mpz_clear(R1);
  }
}

static void _ramanujan_pi(mpf_t pi, unsigned long prec)
{
  unsigned long terms = (1 + DIGS2BITS(prec)/47.11041314);
  mpz_t P, Q, R, u;
  mpf_t t;

  mpz_init(P); mpz_init(Q); mpz_init(R);
  mpz_init(u);
  _sum_pqr(P, Q, R, u, 0, terms);
  mpz_clear(u);

  mpz_mul_ui(R, Q, 13591409UL);
  mpz_add(P, P, R);
  mpz_mul_ui(Q, Q, 4270934400UL);

  /* pi = Q / (P * sqrt(10005)) */
  mpf_init2(t, mpf_get_prec(pi));
  mpf_set_ui(t, 10005);
  mpf_sqrt(t, t);
  mpf_set_z(pi, P);
  mpf_mul(t, t, pi);
  mpf_set_z(pi, Q);
  mpf_div(pi, pi, t);
  mpf_clear(t);
  mpz_clear(R); mpz_clear(Q); mpz_clear(P);
}

static void _agm_pi(mpf_t pi, unsigned long prec)
{
  mpf_t t, an, bn, tn, prev_an;
  unsigned long k, bits = ceil(prec * 3.322);

  mpf_init2(t,       10+bits);
  mpf_init2(an,      10+bits);
  mpf_init2(bn,      10+bits);
  mpf_init2(tn,      10+bits);
  mpf_init2(prev_an, 10+bits);

  mpf_set_ui(an, 1);
  mpf_div_2exp(bn, an, 1);
  mpf_div_2exp(tn, an, 2);
  mpf_sqrt(bn, bn);
                                    /* Comments from Brent 1976 */
  for (k = 0; (prec >> (k+1)) > 0; k++) {
    mpf_set(prev_an, an);           /* Y <- A */
    mpf_add(t, an, bn);
    mpf_div_2exp(an, t, 1);         /* A <- (A+B)/2 */
    mpf_mul(t, bn, prev_an);
    mpf_sqrt(bn, t);                /* B <- (BY)^(1/2) */
    mpf_sub(prev_an, prev_an, an);
    mpf_mul(t, prev_an, prev_an);
    mpf_mul_2exp(t, t, k);
    mpf_sub(tn, tn, t);             /* T <- T - 2^k (A-Y)^2 */
#if 0 /* Instead of doing the comparison, we assume doubling per iteration */
    mpf_sub(t, an, bn);
    mpf_mul_2exp(t, t, bits);
    if (mpf_cmp_ui(t,1) <= 0) break;
#endif
  }
  mpf_add(t, an, bn);
  mpf_mul(an, t, t);
  mpf_mul_2exp(t, tn, 2);
  mpf_div(pi, an, t);               /* return (A+B)^2 / 4T */
  mpf_clear(tn); mpf_clear(bn); mpf_clear(an);
  mpf_clear(prev_an); mpf_clear(t);
}

static void _const_pi(mpf_t pi, unsigned long prec)
{
  if (prec <= 100) {
    mpf_set_str(pi, "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798215", 10);
  } else if (prec <= 3000) {
    _agm_pi(pi, prec);
  } else {
    _ramanujan_pi(pi, prec);
  }
}

/* http://numbers.computation.free.fr/Constants/Log2/log2.ps
 * Machin-like formula 25. */
static void _const_log2(mpf_t logn, unsigned long prec)
{
  mpf_t t;
  mpz_t t1, t2, term1, term2, pows;
  unsigned long bits = precbits(logn, prec, 64);

  mpz_init(t1); mpz_init(t2); mpz_init(term1); mpz_init(term2); mpz_init(pows);
  mpf_init2(t, bits);
  mpz_ui_pow_ui(pows, 10, 20+prec);
  mpz_arctanh(term1,   26, pows, t1, t2);  mpz_mul_ui(term1, term1,  18);
  mpz_arctanh(term2, 4801, pows, t1, t2);  mpz_mul_ui(term2, term2,   2);
  mpz_sub(term1, term1, term2);
  mpz_arctanh(term2, 8749, pows, t1, t2);  mpz_mul_ui(term2, term2,   8);
  mpz_add(term1, term1, term2);
  /* term1 = 69313147... */
  mpf_set_z(logn, term1);
  mpf_set_z(t, pows);
  mpf_div(logn, logn, t);
  /* logn = .69313147... */
  mpf_clear(t);
  mpz_clear(t1); mpz_clear(t2); mpz_clear(term1); mpz_clear(term2); mpz_clear(pows);
}

/* Cache constants.  We should thread lock these. */
static mpf_t _fconst_euler, _fconst_pi, _fconst_log2;
static unsigned long _prec_euler = 0, _prec_pi = 0, _prec_log2 = 0;

#define CONST_FUNC(name) \
  void const_##name(mpf_t c, unsigned long prec) { \
    if (prec > _prec_##name) { \
      prec += 10; \
      if (_prec_##name == 0) mpf_init2(_fconst_##name, 7+DIGS2BITS(prec)); \
      else                   mpf_set_prec(_fconst_##name, 7+DIGS2BITS(prec)); \
      _const_##name(_fconst_##name, prec); \
      _prec_##name = prec; \
    } \
    mpf_set(c, _fconst_##name); \
  }

CONST_FUNC(euler);
CONST_FUNC(pi);
CONST_FUNC(log2);

void free_float_constants(void) {
  if (_prec_euler) { _prec_euler = 0;  mpf_clear(_fconst_euler); }
  if (_prec_pi)    { _prec_pi    = 0;  mpf_clear(_fconst_pi);    }
  if (_prec_log2)  { _prec_log2  = 0;  mpf_clear(_fconst_log2);  }
}

/*****************     Exponential / Logarithmic Integral     *****************/

void li(mpf_t r, mpf_t n, unsigned long prec)
{
  mpz_t factorial;
  mpf_t logn, sum, inner_sum, term, p, q, tol;
  unsigned long j, k, bits = precbits(r, prec, 10);

  mpf_init2(logn,      bits);
  mpf_log(logn, n);

  mpf_init2(sum,       bits);
  mpf_init2(inner_sum, bits);
  mpf_init2(term,      bits);
  mpf_init2(p,         bits);
  mpf_init2(q,         bits);
  mpf_init2(tol,       bits);

  mpf_set_ui(tol, 10);  mpf_pow_ui(tol, tol, prec);  mpf_ui_div(tol,1,tol);

  mpz_init_set_ui(factorial, 1);

  mpf_set_si(p, -1);
  for (j = 1, k = 0; j < 1000000; j++) {
    mpz_mul_ui(factorial, factorial, j);
    mpf_mul(p, p, logn);
    mpf_neg(p, p);
    for (; k <= (j - 1) / 2; k++) {
      mpf_set_ui(q, 1);
      mpf_div_ui(q, q, 2*k+1);
      mpf_add(inner_sum, inner_sum, q);
    }
    mpf_set_z(q, factorial);
    mpf_mul_2exp(q, q, j-1);
    mpf_mul(term, p, inner_sum);
    mpf_div(term, term, q);
    mpf_add(sum, sum, term);

    mpf_abs(term, term);
    mpf_mul(q, sum, tol);
    mpf_abs(q, q);
    if (mpf_cmp(term, q) <= 0) break;
  }
  mpf_sqrt(q, n);
  mpf_mul(r, sum, q);

  mpf_abs(logn,logn);
  mpf_log(q, logn);
  mpf_add(r, r, q);

  /* Find out roughly how many digits of C we need, then get it and add */
  mpf_set(q, r);
  for (k = prec; k > 100 && mpf_cmp_ui(q, 1024*1024) >= 0; k -= 6)
    mpf_div_2exp(q, q, 20);
  const_euler(q, k);
  mpf_add(r, r, q);

  mpz_clear(factorial);
  mpf_clear(tol);
  mpf_clear(q);
  mpf_clear(p);
  mpf_clear(term);
  mpf_clear(inner_sum);
  mpf_clear(sum);
  mpf_clear(logn);
}

void ei(mpf_t r, mpf_t x, unsigned long prec)
{
#if 0
  #include <mpfr.h>
  mpfr_t C, Cin;
  mpfr_init2(C,    10+DIGS2BITS(prec));
  mpfr_init2(Cin,  mpf_get_prec(x));
  mpfr_set_f(Cin, x, MPFR_RNDN);
  mpfr_eint(C, Cin, MPFR_RNDN);
  mpfr_get_f(r, C, MPFR_RNDN);
  mpfr_clear(Cin);
  mpfr_clear(C);
  return;
#endif
  if (mpf_sgn(x) > 0 && mpf_cmp_ui(x, 100) < 0) {
    /* x > 0 only */
    mpf_t factn, invn, term, sum, t, tol;
    unsigned long n, bits = precbits(r, prec, 14);

    mpf_init2(factn, bits);
    mpf_init2(invn,  bits);
    mpf_init2(term,  bits);
    mpf_init2(sum,   bits);
    mpf_init2(t,     bits);
    mpf_init2(tol,   bits);

    mpf_set_ui(tol, 10);  mpf_pow_ui(tol, tol, prec+4);  mpf_ui_div(tol,1,tol);
    mpf_set(factn, x);

    for (n = 2; n <= 1000000; n++) {
      mpf_set_ui(t, n);
      mpf_ui_div(invn, 1, t);
      mpf_mul(t, x, invn);
      mpf_mul(factn, factn, t);
      mpf_mul(term, factn, invn);
      mpf_add(sum, sum, term);

      mpf_abs(term, term);
      mpf_mul(t, sum, tol);
      mpf_abs(t, t);
      if (mpf_cmp(term, t) <= 0) break;
    }
    const_euler(t, prec+4);  mpf_add(sum, sum, t);
    mpf_log(t, x);           mpf_add(sum, sum, t);
    mpf_add(sum, sum, x);
    mpf_set(r, sum);

    mpf_clear(tol); mpf_clear(t); mpf_clear(sum);
    mpf_clear(term); mpf_clear(invn); mpf_clear(factn);
  } else {
    mpf_exp(r, x);
    li(r, r, prec+3);
  }
}


/***************************        Harmonic        ***************************/

static void _harmonic(mpz_t a, mpz_t b, mpz_t t) {
  mpz_sub(t, b, a);
  if (mpz_cmp_ui(t, 1) <= 0) {
    mpz_set(b, a);
    mpz_set(a, t);
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

/**************************        Bernoulli        **************************/

static void _bern_real_zeta(mpf_t bn, unsigned long s, unsigned long prec)
{
  mpf_t tf;

  if (s & 1) {
    mpf_set_d(bn, (s == 1) ? 0.5 : 0.0);
    return;
  }

  mpf_init2(tf, mpf_get_prec(bn));

  /* For large values with low precision, we should look at approximations.
   *   http://www.ebyte.it/library/downloads/2008_MTH_Nemes_GammaApproximationUpdate.pdf
   *   http://www.luschny.de/math/primes/bernincl.html
   *   http://arxiv.org/pdf/math/0702300.pdf
   */

  _zetaint(bn, s, prec);

  /* We should be using an approximation here, e.g. Pari's mpfactr.  For
   * large values this is the majority of time taken for this function. */
  { mpz_t t; mpz_init(t); mpz_fac_ui(t, s); mpf_set_z(tf, t); mpz_clear(t);}
  mpf_mul(bn, bn, tf);
  /* bn = s! * zeta(s) */

  const_pi(tf, prec);
  mpf_mul_ui(tf, tf, 2);
  mpf_pow_ui(tf, tf, s);
  mpf_div(bn, bn, tf);
  /* bn = s! * zeta(s) / (2Pi)^s */

  mpf_mul_2exp(bn, bn, 1);
  if ((s & 3) == 0) mpf_neg(bn, bn);
  /* bn = (-1)^(n-1) * 2 * s! * zeta(s) / (2Pi)^s */
  mpf_clear(tf);
}

/* Compute a vector of n+1 even Bernoulli numbers:  B[0],B[2],...,B[2n] */
static void _bernoulli_vector(mpz_t** pN, mpz_t **pD, unsigned long n) {
  mpz_t *T, *N, *D, g, den, p;
  unsigned long i, j, k, h;

  New(0, T, n+1, mpz_t);
  New(0, N, n+1, mpz_t);
  New(0, D, n+1, mpz_t);

  for (i = 0; i <= n; i++) {
    mpz_init(T[i]);
    mpz_init(N[i]);
    mpz_init(D[i]);
  }
  mpz_set_ui(N[0], 1);
  mpz_set_ui(D[0], 1);
  if (n >= 1)
    mpz_set_ui(T[1],1);

  /* Use Luschny's Seidel method */

  mpz_init_set_ui(den, 1);
  mpz_init_set_ui(p, 1);
  mpz_init(g);

  for (i = 1, j = 1, h = 0; i <= 2*n; i++) {
    if (i & 1) {
      mpz_mul_ui(p, p, 4);
      mpz_sub_ui(den, p, 1);
      mpz_mul_2exp(den, den, 1);
      for (k = h++; k > 0; k--)
        mpz_add(T[k], T[k], T[k+1]);
    } else {
      for (k = 1; k <= h; k++)
        mpz_add(T[k], T[k], T[k-1]);
      mpz_gcd(g, T[h], den);
      mpz_divexact(N[j], T[h], g);
      mpz_divexact(D[j], den, g);
      if (!(j&1)) mpz_neg(N[j],N[j]);
      j++;
    }
  }
  mpz_clear(g); mpz_clear(p); mpz_clear(den);

  for (i = 0; i <= n; i++)
    mpz_clear(T[i]);
  Safefree(T);
  *pN = N;
  *pD = D;
}

static int _bern_cache_init = 0;
static unsigned long _bern_cache_n = 0;
static mpz_t *_bern_cache_NUM = 0;
static mpz_t *_bern_cache_DEN = 0;

void free_bernoulli(void) {
  if (_bern_cache_init) {
    unsigned long i, n = _bern_cache_n;
    mpz_t *N = _bern_cache_NUM, *D = _bern_cache_DEN;
    _bern_cache_NUM = _bern_cache_DEN = 0;
    _bern_cache_n = 0;
    _bern_cache_init = 0;
    for (i = 0; i <= n; i++) {
      mpz_clear(N[i]);
      mpz_clear(D[i]);
    }
    Safefree(N);
    Safefree(D);
  }
}

static void _fill_bern_cache(unsigned long n) {
  if (n < 100)
     n = 100; /* Make it at least this large. */
  if (_bern_cache_init && _bern_cache_n >= n)
    return;
  free_bernoulli();
  _bernoulli_vector( &_bern_cache_NUM, &_bern_cache_DEN, n);
  _bern_cache_n = n;
  _bern_cache_init = 1;
}

static int _get_bern_cache(mpz_t num, mpz_t den, unsigned long n) {
  unsigned long k = n >> 1;
  if (n <= 1 || (n & 1)) {
    mpz_set_ui(num, (n<=1) ? 1 : 0);
    mpz_set_ui(den, (n==1) ? 2 : 1);
    return 1;
  } else if (_bern_cache_init && k <= _bern_cache_n) {
    mpz_set(num, _bern_cache_NUM[k]);
    mpz_set(den, _bern_cache_DEN[k]);
    return 1;
  }
  return 0;
}

/* Return first n even Bernoulli numbers: B[0], B[2], ... B[2n] as READ ONLY */
void bernvec(const mpz_t **N, const mpz_t **D, unsigned long n) {
  _fill_bern_cache(n);
  *N = _bern_cache_NUM;
  *D = _bern_cache_DEN;
}


static void _bernfrac_comb(mpz_t num, mpz_t den, unsigned long n, mpz_t t)
{
  unsigned long k, j;
  mpz_t* T;

  if (n <= 1 || (n & 1)) {
    mpz_set_ui(num, (n<=1) ? 1 : 0);
    mpz_set_ui(den, (n==1) ? 2 : 1);
    return;
  }

  /* Denominator */
  mpz_set_ui(t, 1);
  mpz_mul_2exp(den, t, n);    /* den = U = 1 << n  */
  mpz_sub_ui(t, den, 1);      /* t = U-1            */
  mpz_mul(den, den, t);       /* den = U*(U-1)      */

  n >>= 1;

  /* Luschny's version of the "Brent-Harvey" method */
  /* Algorithm TangentNumbers from https://arxiv.org/pdf/1108.0286.pdf */
  New(0, T, n+1, mpz_t);
  for (k = 1; k <= n; k++)  mpz_init(T[k]);
  mpz_set_ui(T[1], 1);

  for (k = 2; k <= n; k++)
    mpz_mul_ui(T[k], T[k-1], k-1);

  for (k = 2; k <= n; k++) {
    for (j = k; j <= n; j++) {
      mpz_mul_ui(t, T[j], j-k+2);
      mpz_mul_ui(T[j], T[j-1], j-k);
      mpz_add(T[j], T[j], t);
    }
  }

  /* (14), also last line of Algorithm FastTangentNumbers from paper */
  mpz_mul_ui(num, T[n], n);
  mpz_mul_si(num, num, (n & 1) ? 2 : -2);

  for (k = 1; k <= n; k++)  mpz_clear(T[k]);
  Safefree(T);
}

static void _bernfrac_zeta(mpz_t num, mpz_t den, unsigned long n, mpz_t t)
{
  unsigned long prec;
  double nbits;
  mpf_t bn, tf;
  /* Compute integer numerator by getting the real bn first. */

  if (n <= 1 || (n & 1)) {
    mpz_set_ui(num, (n<=1) ? 1 : 0);
    mpz_set_ui(den, (n==1) ? 2 : 1);
    return;
  }
  if (n == 2) { mpz_set_ui(num, 1); mpz_set_ui(den, 6); return; }

  /* Calculate denominator */
  {
    int i, ndivisors;
    mpz_t *D;

    mpz_set_ui(t, n >> 1);
    D = divisor_list(&ndivisors, t);
    mpz_set_ui(den, 6);
    for (i = 1; i < ndivisors; i++) {
      mpz_mul_2exp(t,D[i],1);  mpz_add_ui(t,t,1);
      if (_GMP_is_prime(t))
        mpz_mul(den, den, t);
    }
    for (i = 0; i < ndivisors; i++)
      mpz_clear(D[i]);
    Safefree(D);
  }

  /* Estimate number of bits, from Pari, also see Stein 2006 */
  nbits = mpz_logn(den) + (n+0.5) * log((double)n) - n*2.8378770664093454835606594728L + 1.712086L;
  nbits /= log(2);
  nbits += 32;
  prec = (unsigned long)(nbits/3.32193 + 1);

  mpf_init2(bn, nbits);
  mpf_init2(tf, nbits);
  _bern_real_zeta(bn, n, prec);
  mpf_set_z(tf, den);
  mpf_mul(bn, bn, tf);

  mpf_set_d(tf, (mpf_sgn(bn) < 0) ? -0.5 : 0.5);
  mpf_add(bn, bn, tf);

  mpz_set_f(num, bn);

  mpf_clear(tf);
  mpf_clear(bn);
}


void bernfrac(mpz_t num, mpz_t den, mpz_t zn)
{
  unsigned long n = mpz_get_ui(zn);
  mpz_t t;

  if (n < 100)
    _fill_bern_cache(n);
  if (_get_bern_cache(num, den, n))
    return;

  mpz_init(t);
  if (n < 46) {
    _bernfrac_comb(num, den, n, t);
  } else {
    _bernfrac_zeta(num, den, n, t);
  }
  mpz_gcd(t, num, den);
  mpz_divexact(num, num, t);
  mpz_divexact(den, den, t);
  mpz_clear(t);
}

/***************************       Lambert W       ***************************/

static double _lambertw_approx(double x) {
  double w, k1, k2, k3;

  if (x < -0.312) {
    /* Near the branch point.  See Fukushima (2013) section 2.5. */
    k2 = 2.0 * (1.0 + 2.71828182845904523536 * x);
    if (k2 <= 0) return -1.0 + 1*DBL_EPSILON;
    k1 = sqrt(k2);
    /* Use Puiseux series, e.g. Verberic 2009, Boost, Johannson (2020). */
    w = -1.0 + (1.0 + (-1.0/3.0 + (11.0/72.0 + (-43.0/540.0 + (769.0/17280.0 + (-221.0/8505.0 + (680863.0/43545600.0 + (-1963.0/204120.0 + 226287557.0/37623398400.0
    * k1) * k1) * k1) * k1) * k1) * k1) * k1) * k1) * k1;

  } else if (x > -0.14 && x < 0.085) {
    /* Around zero.  See Fukushima (2013) section 2.6. */
    w = (1.0 + (-1.0 + (3.0/2.0 + (-8.0/3.0 + (125.0/24.0 + (-54.0/5.0 + (16807.0/720.0 + (-16384.0/315.0 + 531441.0/4480.0
        * x) * x) * x) * x) * x) * x) * x) * x) * x;

  } else if (x < 1) {
    /* This and the rest from Vazquez-Leal et al. (2019). */
    k1 = sqrt(1.0 + 2.71828182845904523536 * x);
    k2 = 0.33333333333333333333 + 0.7071067811865475244 / k1 - 0.058925565098880 * k1 +
         (x + 0.36787944117144) * (0.050248489761611 + (0.11138904851051 + 0.040744556245195 * x) * x)
         /
         (1.0 + (2.7090878606183 + (1.5510922597820 + 0.095477712183841 * x) * x) * x);
    w = -(k2-1)/k2;

  } else if (x < 40) {
    k1 = 1.0 + (5.950065500550155 + (13.96586471370701 + (10.52192021050505 + (3.065294254265870 + 0.1204576876518760 * x) * x) * x) * x) * x;
    w = 0.1600049638651493 * log(k1);
  } else if (x < 20000) {
    k1 = 1.0 + (-3.16866642511229e11 + (3.420439800038598e10 +
         (-1.501433652432257e9 + (3.44887729947585e7 + (-4.453783741137856e5 +
         (3257.926478908996 + (-10.82545259305382 + (0.6898058947898353e-1 +
         0.4703653406071575e-4 * x) * x) * x) * x) * x) * x) * x) * x) * x;
    w = 0.9898045358731312e-1 * log(k1);

  } else {
    k1 = 1.0 / (1.0 + log(1.0 + x));
    k2 = 1.0 / k1;
    k3 = log(k2);
    w = k2-1-k3+(1+k3+(-1/2+(1/2)*k3*k3 +(-1/6+(-1+(-1/2+
        (1/3) * k3) * k3) * k3) * k1) * k1) * k1;
  }

  /* Improve the FP estimate using two simple Halley iterations. */
  if (x >= -0.36728) {
    if (w != 0) w = (w/(1.0+w)) * (1.0+log(x/w));
    if (w != 0) w = (w/(1.0+w)) * (1.0+log(x/w));
    if (isnan(w)) w = DBL_EPSILON;
  }

  /* Result is over 15 digits, 16 when x > 0 */
  return w;
}


static void _lambertw(mpf_t r, mpf_t x, unsigned long prec)
{
  int i;
  unsigned long bits = 96+mpf_get_prec(r);  /* More bits for intermediate */
  mpf_t w, t, tol, w1, zn, qn, en;

  if (mpf_cmp_d(x, -0.36787944117145) < 0)
    croak("Invalid input to LambertW:  x must be >= -1/e");
  if (mpf_sgn(x) == 0)
    { mpf_set(r, x); return; }

  /* Use Fritsch rather than Halley. */
  mpf_init2(w,   bits);
  mpf_init2(t,   bits);
  mpf_init2(tol, bits);
  mpf_init2(w1,  bits);
  mpf_init2(zn,  bits);
  mpf_init2(qn,  bits);
  mpf_init2(en,  bits);

  /* Initial estimate done in FP instead of mpf. */
  mpf_set_d(w, _lambertw_approx(mpf_get_d(x)));

  /* Divide prec by 2 since t should be have 4x number of zeros each round */
  mpf_set_ui(tol, 10);
  mpf_pow_ui(tol, tol, (mpf_cmp_d(x, -.36) < 0) ? prec : prec/2);
  mpf_ui_div(tol,1,tol);

  for (i = 0; i < 500 && mpz_sgn(w) != 0; i++) {
    mpf_add_ui(w1, w, 1);

    mpf_div(t, x, w);
    mpf_log(zn, t);
    mpf_sub(zn, zn, w);

    mpf_mul_ui(t, zn, 2);
    mpf_div_ui(t, t, 3);
    mpf_add(t, t, w1);
    mpf_mul(t, t, w1);
    mpf_mul_ui(qn, t, 2);

    mpf_sub(en, qn, zn);
    mpf_mul_ui(t, zn, 2);
    mpf_sub(t, qn, t);
    mpf_div(en, en, t);
    mpf_div(t, zn, w1);
    mpf_mul(en, en, t);

    mpf_mul(t, w, en);
    mpf_add(w, w, t);

    mpf_abs(t, t);
    if (mpf_cmp(t, tol) <= 0) break;
    if (mpf_cmp_d(w,-1) <= 0) break;
  }

  mpf_clear(en); mpf_clear(qn); mpf_clear(zn);
  mpf_clear(w1); mpf_clear(tol); mpf_clear(t);

  if (mpf_cmp_d(w, -1) <= 0)
    mpf_set_si(r, -1);
  else
    mpf_set(r, w);
  mpf_clear(w);
}

/*****************************************************************************/

/*****************************************************************************/

char* zetareal(mpf_t z, unsigned long prec)
{
  size_t est_digits = 10+prec;
  char* out;
  if (mpf_cmp_ui(z,1) == 0) return 0;
  if (mpz_sgn(z) < 0) est_digits += -mpf_get_si(z);
  _zeta(z, z, prec);
  New(0, out, est_digits, char);
  gmp_sprintf(out, "%.*Ff", (int)(prec), z);
  return out;
}
char* riemannrreal(mpf_t r, unsigned long prec)
{
  if (mpf_cmp_ui(r,0) <= 0) return 0;
  _riemann_r(r, r, prec);
  return _str_real(r, prec);
}
char* lambertwreal(mpf_t x, unsigned long prec) {
  _lambertw(x, x, prec);
  return _str_real(x, prec);
}
char* lireal(mpf_t r, unsigned long prec)
{
  if (mpf_cmp_ui(r,0) < 0) return 0;
  if (mpf_cmp_ui(r,1) == 0) return 0;
  li(r, r, prec);
  return _str_real(r, prec);
}
char* eireal(mpf_t r, unsigned long prec)
{
  if (mpf_cmp_ui(r,0) == 0) return 0;
  ei(r, r, prec);
  return _str_real(r, prec);
}
#define DEFINE_REAL_1ARG(func, mfunc) \
  char* func(mpf_t r, unsigned long prec) { \
    mfunc(r, r); \
    return _str_real(r, prec); \
  }
#define DEFINE_REAL_2ARG(func, mfunc) \
  char* func(mpf_t r, mpf_t x, unsigned long prec) { \
    mfunc(r, r, x); \
    return _str_real(r, prec); \
  }
DEFINE_REAL_1ARG(logreal, mpf_log)
DEFINE_REAL_1ARG(expreal, mpf_exp)
DEFINE_REAL_2ARG(powreal, mpf_pow)
DEFINE_REAL_2ARG(rootreal, mpf_root)
DEFINE_REAL_2ARG(addreal, mpf_add)
DEFINE_REAL_2ARG(subreal, mpf_sub)
DEFINE_REAL_2ARG(mulreal, mpf_mul)
DEFINE_REAL_2ARG(divreal, mpf_div)

char* agmreal(mpf_t a, mpf_t b, unsigned long prec)
{
  if (mpz_sgn(a) == 0 || mpz_sgn(b) == 0) {
     mpf_set_ui(a,0);
  } else if (mpz_sgn(a) < 0 || mpz_sgn(b) < 0) {
     return 0;  /* NaN */
  }
  mpf_agm(a, a, b);
  return _str_real(a, prec);
}

char* eulerconst(unsigned long prec) {
  char* out;
  mpf_t gamma;
  unsigned long bits = 7 + DIGS2BITS(prec);

  mpf_init2(gamma, bits);
  const_euler(gamma, prec);
  New(0, out, prec+4, char);
  gmp_sprintf(out, "%.*Ff", (int)(prec), gamma);
  mpf_clear(gamma);
  return out;
}
char* piconst(unsigned long prec) {
  char* out;
  mpf_t pi;
  unsigned long bits = 7 + DIGS2BITS(prec);

  mpf_init2(pi, bits);
  const_pi(pi, prec);
  New(0, out, prec+4, char);
  gmp_sprintf(out, "%.*Ff", (int)(prec-1), pi);
  mpf_clear(pi);
  return out;
}

char* harmreal(mpz_t zn, unsigned long prec) {
  char* out;
  mpz_t num, den;

  mpz_init(num); mpz_init(den);
  harmfrac(num, den, zn);
  out = _frac_real(num, den, prec);
  mpz_clear(den); mpz_clear(num);

  return out;
}

char* bernreal(mpz_t zn, unsigned long prec) {
  mpz_t num, den;
  unsigned long n = mpz_get_ui(zn);
  char* out;

  if (n < 100)
    _fill_bern_cache(n);

  mpz_init(num); mpz_init(den);
  if (_get_bern_cache(num, den, n)) {
    out = _frac_real(num, den, prec);
  } else {
    mpf_t z;
    unsigned long bits = 32 + DIGS2BITS(prec);
    mpf_init2(z, bits);
    _bern_real_zeta(z, n, prec);
    out = _str_real(z, prec);
    mpf_clear(z);
  }
  mpz_clear(den); mpz_clear(num);
  return out;
}
