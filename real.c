#include <string.h>
#include <math.h>
#include <gmp.h>
#include "ptypes.h"

#include "real.h"
#include "prime_iterator.h"
#include "primality.h"
#include "factor.h"
#define FUNC_mpz_logn 1
#include "utility.h"


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

static void _bern_real_zeta(mpf_t bn, mpz_t zn, unsigned long prec);
static unsigned long zeta_n = 0;
static mpz_t* zeta_d = 0;

static void _borwein_d(unsigned long D) {
  mpz_t t1, t2, t3, sum;
  unsigned long i, n = 3 + (1.31 * D);

  if (zeta_n >= n)
    return;

  if (zeta_n > 0) {
    for (i = 0; i <= zeta_n; i++)
      mpz_clear(zeta_d[i]);
    Safefree(zeta_d);
  }

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
      mpz_t n;
      mpz_init_set_ui(n, S+1);
      _bern_real_zeta(z, n, prec);
      mpf_div_ui(z, z, S+1);
      mpf_neg(z,z);
      mpz_clear(n);
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
  unsigned long k, bits = mpf_get_prec(n);

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
 * We'll use the series method as Pari does.  We should use binary splitting,
 * as it still isn't really fast.  For high precision it's about 2x slower
 * than Pari due to mpf_log / mpf_exp.
 */
static void _const_euler(mpf_t gamma, unsigned long prec)
{
  const double log2 = 0.693147180559945309417232121458176568L;
  const unsigned long maxsqr = (1UL << (4*sizeof(unsigned long))) - 1;
  unsigned long bits = ceil(40 + prec * 3.322);
  unsigned long x = ceil((2 + bits) * log2/4);
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

  mpf_set_ui(u, x);
  mpf_log(u, u);
  mpf_neg(u, u);
  mpf_set(a, u);
  mpf_set_ui(b, 1);
  mpf_set_ui(v, 1);

  if (x <= maxsqr && N <= maxsqr) {
    /*  v_k = x^2 / k^2  |  u_k = x^2 / k + v_k) / k  */
    for (k = 1; k <= N; k++) {
      mpf_mul_ui(b, b, xx);
      mpf_div_ui(b, b, k*k);
      mpf_mul_ui(a, a, xx);
      mpf_div_ui(a, a, k);
      mpf_add(a, a, b);
      mpf_div_ui(a, a, k);
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
 *
 * - Ramanujan / Chudnovsky with binary splitting.
 *   About 2-4x faster than AGM for large enough sizes.  This version is
 *   based on Alexander Yee's example.  I have tested with a port of Pari/GP's
 *   abpq_sum and it came out about the same speed (albeit is more generic).
 *   There are many more optimizations that can be done for this.
 */

static void _sum_pqr(mpz_t P, mpz_t Q, mpz_t R, unsigned long a, unsigned long b)
{
  if (b-a == 1) {
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
    mpz_mul(P, P, R);
    if (b % 2 == 1) mpz_neg(P, P);
  } else {
    mpz_t P0, Q0, R0, P1, Q1, R1;
    unsigned long m = a + (b-a)/2;

    mpz_init(P0); mpz_init(Q0); mpz_init(R0);
    mpz_init(P1); mpz_init(Q1); mpz_init(R1);

    _sum_pqr(P0, Q0, R0, a, m);
    _sum_pqr(P1, Q1, R1, m, b);

    mpz_mul(Q, P1, R0);
    mpz_mul(P, P0, Q1);  mpz_add(P, P, Q);
    mpz_mul(Q, Q0, Q1);
    mpz_mul(R, R0, R1);

    mpz_clear(P0); mpz_clear(Q0); mpz_clear(R0);
    mpz_clear(P1); mpz_clear(Q1); mpz_clear(R1);
  }
}

static void _ramanujan_pi(mpf_t pi, unsigned long prec)
{
  unsigned long terms = (1 + DIGS2BITS(prec)/47.11041314);
  mpz_t P, Q, R;
  mpf_t t;

  mpz_init(P); mpz_init(Q); mpz_init(R);
  _sum_pqr(P, Q, R, 0, terms);

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
  mpz_clear(R); mpz_clear(Q); mpz_clear(P);
  mpf_clear(t);
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

  mpf_set_d(an, 1);
  mpf_set_d(bn, 0.5);
  mpf_set_d(tn, 0.25);
  mpf_sqrt(bn, bn);
                                    /* Comments from Brent 1976 */
  for (k = 0; (prec >> (k+1)) > 0; k++) {
    mpf_set(prev_an, an);           /* Y <- A */
    mpf_add(t, an, bn);
    mpf_div_ui(an, t, 2);           /* A <- (A+B)/2 */
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
  } else if (prec <= 6000) {
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
  unsigned long bits = mpf_get_prec(logn);

  mpz_init(t1); mpz_init(t2); mpz_init(term1); mpz_init(term2); mpz_init(pows);
  mpf_init2(t, 64+bits);
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
  mpz_clear(t1); mpz_clear(t2); mpz_clear(term2); mpz_clear(pows);
}

/* Cache constants.  We should thread lock these. */
static mpf_t _fconst_euler, _fconst_pi, _fconst_log2;
static unsigned long _prec_euler = 0, _prec_pi = 0, _prec_log2 = 0;

#define CONST_FUNC(name) \
  void const_##name(mpf_t c, unsigned long prec) { \
    if (prec > _prec_##name) { \
      prec += 10; \
      if (_prec_##name == 0) mpf_init2(_fconst_##name, DIGS2BITS(prec)); \
      else                   mpf_set_prec(_fconst_##name, DIGS2BITS(prec)); \
      _const_##name(_fconst_##name, prec); \
      _prec_##name = prec; \
    } \
    mpf_set(c, _fconst_##name); \
  }

CONST_FUNC(euler);
CONST_FUNC(pi);
CONST_FUNC(log2);

void free_constants(void) {
  _prec_euler = 0;  mpf_clear(_fconst_euler);
  _prec_pi    = 0;  mpf_clear(_fconst_pi);
  _prec_log2  = 0;  mpf_clear(_fconst_log2);
}

/*****************     Exponential / Logarithmic Integral     *****************/

static void _li_r(mpf_t r, mpf_t n, unsigned long prec)
{
  mpz_t factorial;
  mpf_t logn, sum, inner_sum, term, p, q, tol;
  unsigned long j, k, bits = 10 + mpf_get_prec(n);

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
  for (k = prec; mpf_cmp_ui(q, 1024*1024) >= 0; k -= 6)
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

static void _ei_r(mpf_t r, mpf_t n, unsigned long prec)
{
  mpf_exp(r, n);
  _li_r(r, r, prec+3);
}


/***************************        Harmonic        ***************************/

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

/**************************        Bernoulli        **************************/

static void _bern_real_zeta(mpf_t bn, mpz_t zn, unsigned long prec)
{
  unsigned long s = mpz_get_ui(zn);
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


static void _bernfrac_comb(mpz_t num, mpz_t den, mpz_t zn, mpz_t t)
{
  unsigned long k, j, n = mpz_get_ui(zn);
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


static void _bernfrac_zeta(mpz_t num, mpz_t den, mpz_t zn, mpz_t t)
{
  unsigned long prec, n = mpz_get_ui(zn);
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
  _bern_real_zeta(bn, zn, prec);
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
  mpz_t t;
  mpz_init(t);

  if (mpz_cmp_ui(zn,46) < 0) {
    _bernfrac_comb(num, den, zn, t);
  } else {
    _bernfrac_zeta(num, den, zn, t);
  }

  mpz_gcd(t, num, den);
  mpz_divexact(num, num, t);
  mpz_divexact(den, den, t);
  mpz_clear(t);
}

/***************************       Lambert W       ***************************/

static void _lambertw(mpf_t w, mpf_t x, unsigned long prec)
{
  int i;
  unsigned long bits = 96+mpf_get_prec(x);  /* More bits for intermediate */
  mpf_t t, w1, zn, qn, en, tol;

  if (mpf_cmp_d(x, -0.36787944117145) < 0)
    croak("Invalid input to LambertW:  x must be >= -1/e");
  if (mpf_sgn(x) == 0)
    { mpf_set(w, x); return; }

  /* Use Fritsch rather than Halley. */
  mpf_init2(t,   bits);
  mpf_init2(w1,  bits);
  mpf_init2(zn,  bits);
  mpf_init2(qn,  bits);
  mpf_init2(en,  bits);
  mpf_init2(tol, bits);

  /* Initial estimate */
  if (mpf_cmp_d(x, -0.06) < 0) {  /* Pade(3,2) */
    mpf_set_d(t, 5.4365636569180904707205749);
    mpf_mul(t, t, x);
    mpf_add_ui(t, t, 2);
    if (mpf_sgn(t) <= 0) { mpf_set_ui(t, 0); } else { mpf_sqrt(t, t); }
    mpf_mul(zn, t, t);
    mpf_mul(qn, zn, t);

    mpf_set_d(w, -1);
    mpf_set_d(w1, 1.0L/6.0L);      mpf_mul(w1, w1,  t);  mpf_add(w, w, w1);
    mpf_set_d(w1, 257.0L/720.0L);  mpf_mul(w1, w1, zn);  mpf_add(w, w, w1);
    mpf_set_d(w1, 13.0L/720.0L);   mpf_mul(w1, w1, zn);  mpf_add(w, w, w1);
    mpf_set(en, w);  /* numerator */

    mpf_set_d(w, 1);
    mpf_set_d(w1, 5.0L/6.0L);      mpf_mul(w1, w1,  t);  mpf_add(w, w, w1);
    mpf_set_d(w1, 103.0L/720.0L);  mpf_mul(w1, w1, zn);  mpf_add(w, w, w1);

    mpf_div(w, en, w);
  } else if (mpf_cmp_d(x, 1.363) < 0) {  /* Winitzki 2003 */
    mpf_add_ui(t, x, 1);
    mpf_log(w1, t);
    mpf_add_ui(zn, w1, 1);
    mpf_log(zn, zn);
    mpf_add_ui(qn, w1, 2);
    mpf_div(t, zn, qn);
    mpf_ui_sub(t, 1, t);
    mpf_mul(w, w1, t);
  } else if (mpf_cmp_d(x, 3.7) < 0) {  /* Vargas 2013 modified */
    mpf_log(w, x);
    mpf_log(w1, w);
    mpf_div(t, w1, w);
    mpf_ui_sub(t, 1, t);
    mpf_log(t, t);
    mpf_div_ui(t, t, 2);
    mpf_sub(w, w, w1);
    mpf_sub(w, w, t);
  } else {  /* Corless et al. 1993 */
    mpf_t l1, l2, d1, d2, d3;
    mpf_init2(l1,bits); mpf_init2(l2,bits);
    mpf_init2(d1,bits); mpf_init2(d2,bits); mpf_init2(d3,bits);

    mpf_log(l1, x);
    mpf_log(l2, l1);
    mpf_mul(d1, l1, l1);  mpf_mul_ui(d1, d1, 2);
    mpf_mul(d2, l1, d1);  mpf_mul_ui(d2, d2, 3);
    mpf_mul(d3, l1, d2);  mpf_mul_ui(d3, d3, 2);

    mpf_sub(w, l1, l2);

    mpf_div(t, l2, l1);
    mpf_add(w, w, t);

    mpf_sub_ui(t, l2, 2);
    mpf_mul(t, t, l2);
    mpf_div(t, t, d1);
    mpf_add(w, w, t);

    mpf_mul_ui(t, l2, 2);
    mpf_sub_ui(t, t, 9);
    mpf_mul(t, t, l2);
    mpf_add_ui(t, t, 6);
    mpf_mul(t, t, l2);
    mpf_div(t, t, d2);
    mpf_add(w, w, t);

    mpf_mul_ui(t, l2, 3);
    mpf_sub_ui(t, t, 22);
    mpf_mul(t, t, l2);
    mpf_add_ui(t, t, 36);
    mpf_mul(t, t, l2);
    mpf_sub_ui(t, t, 12);
    mpf_mul(t, t, l2);
    mpf_div(t, t, d3);
    mpf_add(w, w, t);

    mpf_clear(l1); mpf_clear(l2);
    mpf_clear(d1); mpf_clear(d2); mpf_clear(d3);
  }

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

  if (mpf_cmp_d(w, -1) <= 0)
    mpf_set_si(w, -1);

  mpf_clear(en); mpf_clear(qn); mpf_clear(zn);
  mpf_clear(w1); mpf_clear(t); mpf_clear(tol);
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
char* lireal(mpf_t r, unsigned long prec)
{
  if (mpf_cmp_ui(r,0) < 0) return 0;
  if (mpf_cmp_ui(r,1) == 0) return 0;
  _li_r(r, r, prec);
  return _str_real(r, prec);
}
char* eireal(mpf_t r, unsigned long prec)
{
  if (mpf_cmp_ui(r,0) == 0) return 0;
  _ei_r(r, r, prec);
  return _str_real(r, prec);
}
char* logreal(mpf_t r, unsigned long prec)
{
  mpf_log(r, r);
  return _str_real(r, prec);
}
char* expreal(mpf_t r, unsigned long prec)
{
  mpf_exp(r, r);
  return _str_real(r, prec);
}
char* powreal(mpf_t r, mpf_t x, unsigned long prec)
{
  mpf_pow(r, r, x);
  return _str_real(r, prec);
}
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
  unsigned long bits = 7 + (unsigned long)(prec*3.32193);

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
  unsigned long bits = 7 + (unsigned long)(prec*3.32193);

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
  char* out;

  if (mpz_cmp_ui(zn,40) < 0) {
    mpz_t num, den, t;
    mpz_init(num); mpz_init(den); mpz_init(t);
    _bernfrac_comb(num, den, zn, t);
    out = _frac_real(num, den, prec);
    mpz_clear(t); mpz_clear(den); mpz_clear(num);
  } else {
    mpf_t z;
    unsigned long bits = 32+(unsigned long)(prec*3.32193);
    mpf_init2(z, bits);
    _bern_real_zeta(z, zn, prec);
    out = _str_real(z, prec);
    mpf_clear(z);
  }
  return out;
}

char* lambertwreal(mpf_t x, unsigned long prec) {
  char* out;
  mpf_t w;
  unsigned long bits = 64+(unsigned long)(prec*3.32193);

  mpf_init2(w, bits);
  _lambertw(w, x, 10+prec);
  out = _str_real(w, prec);
  mpf_clear(w);
  return out;
}
