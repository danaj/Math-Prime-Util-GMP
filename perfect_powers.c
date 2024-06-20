#include <gmp.h>
#include "ptypes.h"

#include "perfect_powers.h"
#include "misc_ui.h"
#include "utility.h"

int is_perfect_power(const mpz_t n)
{
  const unsigned char smallres[10] = {0,1,0,0,1,0,0,0,1,1};
  UV res;
  if (mpz_sgn(n) < 0) {
    mpz_t N;
    if (mpz_cmp_si(n,-1) == 0) return 1;
    mpz_init(N);
    mpz_neg(N,n);
    res = is_power(N,0);
    mpz_clear(N);
    return (res > 2 && ((res & (res-1)) != 0));
  }
  if (mpz_cmp_ui(n,9) <= 0)  return smallres[mpz_get_ui(n)];
  return (is_power(n,0) > 1);
}

void next_perfect_power(mpz_t next, const mpz_t n)
{
  UV power, log2n, k;
  mpz_t N, best, r;

  if (mpz_sgn(n) < 0) {
    if (mpz_cmp_si(n,-1) == 0) { mpz_set_ui(next,1);  return; }
    mpz_init(N);
    mpz_neg(N, n);
    do {
      prev_perfect_power(N, N);
      power = is_power(N, 0);
    } while (mpz_cmp_ui(N,1) > 0 && (power <= 2 || (power & (power-1)) == 0));
    mpz_neg(next, N);
    mpz_clear(N);
    return;
  }

  if (mpz_sgn(n) == 0) { mpz_set_ui(next, 1); return; }
  if (mpz_cmp_ui(n,1) == 0) { mpz_set_ui(next,4); return; }

  mpz_init(r);
  mpz_init(best);

  mpz_sqrt(r, n);
  mpz_add_ui(r, r, 1);
  mpz_pow_ui(best, r, 2);
  log2n = mpz_sizeinbase(n,2);
  for (k = 3; k <= 1+log2n; k++) {
    mpz_root(r, n, k);
    mpz_add_ui(r, r, 1);
    mpz_pow_ui(r, r, k);
    if (mpz_cmp(r, best) < 0 && mpz_cmp(r, n) > 0)
      mpz_set(best, r);
  }
  mpz_set(next, best);

  mpz_clear(best);
  mpz_clear(r);
}
void prev_perfect_power(mpz_t prev, const mpz_t n)
{
  UV power, log2n, k;
  mpz_t N, best, r, c;

  if (mpz_sgn(n) < 0) {
    mpz_init(N);
    mpz_neg(N, n);
    do {
      next_perfect_power(N, N);
      power = is_power(N, 0);
    } while (mpz_cmp_ui(N,1) > 0 && (power <= 2 || (power & (power-1)) == 0));
    mpz_neg(prev, N);
    mpz_clear(N);
    return;
  }

  if (mpz_cmp_ui(n,4) <= 0) {
    mpz_set_si(prev, (mpz_cmp_ui(n,1) > 0) ? 1 : -1);
    return;
  }

  mpz_init(r);
  mpz_init(c);
  mpz_init(best);

  mpz_set_ui(best, 4);
  log2n = mpz_sizeinbase(n,2);
  for (k = 2; k <= log2n; k++) {
    mpz_root(r, n, k);
    if (mpz_cmp_ui(r,1) > 0) {
      mpz_pow_ui(c, r, k);
      if (mpz_cmp(c,n) >= 0) {
        mpz_sub_ui(r, r, 1);
        mpz_pow_ui(c, r, k);
      }
      if (mpz_cmp(c, best) > 0 && mpz_cmp(c, n) < 0)
        mpz_set(best, c);
    }
  }
  mpz_set(prev, best);

  mpz_clear(best);
  mpz_clear(c);
  mpz_clear(r);
}



void perfect_power_count(mpz_t r, const mpz_t n)
{
  signed char* mu;
  unsigned long k, log2n;
  mpz_t t, count;

  if (mpz_cmp_ui(n,1) <= 0) {
    mpz_set(r, n);
    return;
  }

  log2n = mpz_sizeinbase(n,2);
  mpz_init(t);
  mpz_init_set_ui(count, 1);
  mu = range_moebius(0, log2n);
  for (k = 2; k <= log2n; k++) {
    int m = mu[k];
    if (m != 0) {
      mpz_root(t, n, k);
      mpz_sub_ui(t, t, 1);
      if (m < 0) mpz_add(count, count, t);
      else       mpz_sub(count, count, t);
    }
  }
  Safefree(mu);
  mpz_set(r, count);
  mpz_clear(count);
  mpz_clear(t);
}
void perfect_power_count_range(mpz_t r, const mpz_t lo, const mpz_t hi) {
  if (mpz_cmp(lo, hi) > 0 || mpz_cmp_ui(hi,1) < 0) {
    mpz_set_ui(r, 0);
    return;
  }

  perfect_power_count(r, hi);

  if (mpz_cmp_ui(lo, 1) > 0) {
    mpz_t locount, lom1;
    mpz_init(locount);
    mpz_init(lom1);
    mpz_sub_ui(lom1, lo, 1);
    perfect_power_count(locount, lom1);
    mpz_sub(r, r, locount);
    mpz_clear(lom1);
    mpz_clear(locount);
  }
}



/* r and t must not alias */
static void mpf_mul_d_pow_d(mpf_t r, double m, mpf_t b, double e,  mpf_t t)
{
  mpf_set_d(t, e);
  mpf_pow(r, b, t);    /* This is our function from utility.c */
  mpf_set_d(t, m);
  mpf_mul(r, r, t);
}

void nth_perfect_power_approx(mpz_t nth, const mpz_t n)
{
  mpf_t nf, pp, t, u;

  unsigned long bits = 2 * mpz_sizeinbase(n,2) + 8;

  if (mpz_cmp_ui(n,1) <= 0) {
    mpz_set(nth, n);
    return;
  }

  mpf_init2(nf, bits);
  mpf_init2(pp, bits);
  mpf_init2(t, bits);
  mpf_init2(u, bits);

  mpf_set_z(nf, n);
  mpf_mul(pp, nf, nf);

  mpf_mul_d_pow_d(t, 13./ 3., nf,  8./ 6., u);  mpf_add(pp, pp, t);
  mpf_mul_d_pow_d(t, 32./15., nf, 32./30., u);  mpf_add(pp, pp, t);

  mpf_mul_d_pow_d(t, -2., nf,  5./ 3., u);  mpf_add(pp, pp, t);
  mpf_mul_d_pow_d(t, -2., nf,  7./ 5., u);  mpf_add(pp, pp, t);
  mpf_mul_d_pow_d(t, -2., nf,  9./ 7., u);  mpf_add(pp, pp, t);
  mpf_mul_d_pow_d(t,  2., nf, 12./10., u);  mpf_add(pp, pp, t);
  mpf_mul_d_pow_d(t, -2., nf, 13./11., u);  mpf_add(pp, pp, t);
  mpf_mul_d_pow_d(t, -2., nf, 15./13., u);  mpf_add(pp, pp, t);
  mpf_mul_d_pow_d(t,  2., nf, 16./14., u);  mpf_add(pp, pp, t);
  mpf_mul_d_pow_d(t,  2., nf, 17./15., u);  mpf_add(pp, pp, t);
#if 1
  mpf_mul_d_pow_d(t,-.48, nf, 19./17., u);  mpf_add(pp, pp, t);
#else
  mpf_mul_d_pow_d(t, -2., nf, 19./17., u);  mpf_add(pp, pp, t);
  mpf_mul_d_pow_d(t, -2., nf, 21./19., u);  mpf_add(pp, pp, t);
  mpf_mul_d_pow_d(t,  2., nf, 23./21., u);  mpf_add(pp, pp, t);
  mpf_mul_d_pow_d(t,  2., nf, 24./22., u);  mpf_add(pp, pp, t);
  mpf_mul_d_pow_d(t, -2., nf, 25./23., u);  mpf_add(pp, pp, t);
  mpf_mul_d_pow_d(t,  2., nf, 28./26., u);  mpf_add(pp, pp, t);
#endif
  mpf_set_d(t, 1.5);
  mpf_sub(pp, pp, t);

  mpz_set_f(nth, pp);
  mpf_clear(u);  mpf_clear(t);  mpf_clear(pp);  mpf_clear(nf);
}

void nth_perfect_power(mpz_t nth, const mpz_t n)
{
  const unsigned char smallres[10] = {0, 1, 4, 8, 9, 16, 25, 27, 32, 36};
  mpz_t g, c, apn, diff, t;
  int gn;

  if (mpz_cmp_ui(n,10) < 0) {
    if (mpz_sgn(n) <= 0) mpz_set_ui(nth, 0);
    else                 mpz_set_ui(nth, smallres[mpz_get_ui(n)]);
    return;
  }

  mpz_init(g);
  mpz_init(c);
  mpz_init(apn);
  mpz_init(t);
  mpz_init(diff);
  nth_perfect_power_approx(apn, n);
  mpz_set(g, apn);
  perfect_power_count(c, g);
  for (gn = 0; gn < 10; gn++) {
    mpz_sub(diff, c, n);
    mpz_abs(diff, diff);
    if (mpz_cmp_ui(diff,1000) <= 0) break;
    nth_perfect_power_approx(t, c);
    mpz_sub(t, apn, t);
    mpz_add(g, g, t);
    perfect_power_count(c, g);
  }
  if (mpz_cmp(c, n) >= 0) {
    mpz_add_ui(g,g,1);
    for (prev_perfect_power(g,g);  mpz_cmp(c,n) > 0;  mpz_sub_ui(c,c,1))
      prev_perfect_power(g,g);
  } else {
    for ( ; mpz_cmp(c,n) < 0; mpz_add_ui(c,c,1))
      next_perfect_power(g,g);
  }
  mpz_set(nth, g);
  mpz_clear(diff);
  mpz_clear(t);
  mpz_clear(apn);
  mpz_clear(c);
  mpz_clear(g);
}
