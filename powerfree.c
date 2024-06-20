#include <gmp.h>
#include "ptypes.h"

#include "powerfree.h"
#include "factor.h"
#include "real.h"
#include "utility.h"
#define FUNC_isqrt 1
#define FUNC_ipow 1
#include "misc_ui.h"

int is_powerfree(const mpz_t n, uint32_t k)
{
  int ret = 1;

  if (k < 2 || mpz_cmp_ui(n, 1) <= 0)
    return (mpz_cmp_ui(n,1) == 0);

  if (mpz_sizeinbase(n,2) < k) return 1;  /* Too small */
  if (mpz_scan1(n,0) >= k) return 0;      /* n = 2^k * N */

  if (k == 2)  return is_square_free(n);

#if 0
  if (k == 3) {
    if (mpz_divisible_ui_p(n,  27) || mpz_divisible_ui_p(n, 125) ||
        mpz_divisible_ui_p(n, 343) || mpz_divisible_ui_p(n,1331) ||
        mpz_divisible_ui_p(n,2197))
      return 0;
    if (mpz_cmp_ui(n,4913) < 0)
      return 1;
  }
#endif

  {
    unsigned long f;
    uint32_t e;
    void *iter = trial_factor_iterator_create(n, 1000);
    while (ret && trial_factor_iterator_next(&f, &e, iter))
      if (e >= k)
        ret = 0;
    trial_factor_iterator_destroy(iter);
    if (ret == 0) return 0;
    if (mpz_cmp_ui(n, 1027243729UL) < 0) return 1;  /* next_prime(1000)^3 */
  }

  {
    mpz_t *factors;
    int i, nfactors, *exponents;

    nfactors = factor(n, &factors, &exponents);
    for (i = 0; ret == 1 && i < nfactors; i++)
      if ((uint32_t)exponents[i] >= k)
        ret = 0;
    clear_factors(nfactors, &factors, &exponents);
    return ret;
  }
}

void next_powerfree(mpz_t next, const mpz_t n, uint32_t k)
{
  mpz_t N;

  if (mpz_sgn(n) <= 0) {
    mpz_set_ui(next, 1);
    return;
  }

  mpz_init_set(N, n);
  do {  mpz_add_ui(N, N, 1);  } while (!is_powerfree(N, k));
  mpz_set(next, N);
  mpz_clear(N);
}
void prev_powerfree(mpz_t prev, const mpz_t n, uint32_t k)
{
  mpz_t N;

  if (mpz_cmp_ui(n, 1) <= 0) {
    mpz_set_ui(prev, 0);      /* Nothing before 1 */
    return;
  }

  mpz_init_set(N, n);
  do {  mpz_sub_ui(N, N, 1);  } while (!is_powerfree(N, k));
  mpz_set(prev, N);
  mpz_clear(N);
}



static UV squarefree_count_ui(UV n)
{
  signed char* mu;
  IV *M, *Mx, Mxisum, mert;
  UV I, D, i, j, S1 = 0, S2 = 0;

  if (n < 4) return n;

  I = rootint_ui(n, 5);   /* times loglogn ^ (4/5) */
  D = isqrt(n / I);
  mu = range_moebius(0, D);

  S1 += n;
  New(0, M, D+1, IV);
  M[0] = 0;
  M[1] = 1;
  mert = 1;
  for (i = 2; i <= D; i++) {
    if (mu[i] != 0) {
      S1 += mu[i] * (n/(i*i));
      mert += mu[i];
    }
    M[i] = mert;
  }
  Safefree(mu);

  Newz(0, Mx, I+1, IV);
  Mxisum = 0;
  for (i = I-1; i > 0; i--) {
    IV Mxi = 1;
    UV xi = isqrt(n/i);
    UV L = isqrt(xi);
    for (j = 1; j <= xi/(L+1); j++)
      Mxi -= M[j] * (xi/j - xi/(j+1));
    for (j = 2; j <= L; j++)
      Mxi -=  (xi/j <= D)  ?  M[xi/j]  :  Mx[j*j*i];
    Mx[i] = Mxi;
    Mxisum += Mxi;
  }
  S2 = Mxisum - (I - 1) * M[D];
  Safefree(Mx);
  Safefree(M);

  return S1 + S2;
}

static UV powerfree_count_ui(UV n, uint32_t k)
{
  UV i, nk, count;

  if (k < 2) return (n >= 1);
  if (n < 4) return n;
  if (k == 2) return squarefree_count_ui(n);

  count = n;
  nk = rootint_ui(n, k);

  {
    signed char* mu = range_moebius(0, nk);
    for (i = 2; i <= nk; i++)
      if (mu[i] != 0)
        count += mu[i] * (n/ipow(i,k));
    Safefree(mu);
  }
  return count;
}

void powerfree_count(mpz_t count, const mpz_t n, uint32_t k)
{
  mpz_t t, i, c, c1, nk;
  static const int PERF_MM = 4;  /* Larger = more moebius, less mertens */

  if (k < 2 || mpz_sgn(n) <= 0) {
    mpz_set_ui(count, (mpz_cmp_ui(n,1) >= 0));
    return;
  }
  if (mpz_cmp_ui(n, 4) < 0) {
    mpz_set(count, n);
    return;
  }
  if (mpz_fits_uv_p(n)) {
    UV cnt = powerfree_count_ui(mpz_get_uv(n), k);
    mpz_set_uv(count, cnt);
    return;
  }

  mpz_init(t);
  mpz_init(i);
  mpz_init_set_ui(c,  0);
  mpz_init_set_ui(c1, 0);
  mpz_init(nk);
  mpz_root(nk, n, k);

  if (mpz_sizeinbase(nk,2) <= 50) {
    unsigned long int j, LA, wbeg, wend, nku = mpz_get_ui(nk);
    unsigned long int A = isqrt(nku) / PERF_MM;
    signed char* mu;

    if (A <= 1) { A = 1; LA = nku; }
    else        { mpz_fdiv_q_ui(t,n,A); mpz_root(t,t,k); LA = mpz_get_ui(t); }
    /* gmp_printf("   A %lu  LA %lu\n",A,LA); */

    /* Dense region.  Use a windowed ranged moebius and sum each each value. */
    for (wbeg = 2;  wbeg <= LA;  wbeg = wend+1) {
      wend = ((LA-wbeg) > 10000000UL)  ?  wbeg + 10000000UL - 1  :  LA;
      mu = range_moebius(wbeg, wend);
      for (j = wbeg; j <= wend; j++) {
        if (mu[j-wbeg] != 0) {
          mpz_ui_pow_ui(t, j, k);
          mpz_fdiv_q(t, n, t);
          if (mu[j-wbeg] > 0) mpz_add(c, c, t);
          else                mpz_sub(c, c, t);
        }
      }
      Safefree(mu);
    }
    /* printf("   MOEBIUS DONE\n"); fflush(stdout); */

    /* Sparse region.  Many of the j from LA to nk will result in the same
     * value for n/(j^k).  We walk the possible values instead, multiplying
     * by the moebius values.  Furthermore, we use a fast caching Mertens
     * function instead of summing the moebius values directly.
     */
    if (A > 1) {
      void *mctx = hmertens_create(nku);
      signed long int *M;

      New(0, M, A+1, signed long int);
      M[0] = 0;
      M[1] = hmertens_value(mctx, nku);
      for (j = 2; j <= A; j++) {
        mpz_fdiv_q_ui(t, n, j);
        mpz_root(t, t, k);
        M[j] = hmertens_value(mctx, mpz_get_ui(t));
      }
      hmertens_destroy(mctx);

      for (j = 2; j <= A; j++) {
        mpz_set_si(t, M[j-1] - M[j]);
        mpz_mul_ui(t, t, j-1);
        mpz_add(c1, c1, t);
      }
      Safefree(M);
    }

  } else {
    mpz_t L1;
    signed long int c1i;

    mpz_init(L1);
    mpz_fdiv_q_2exp(t, n, 1);
    mpz_root(L1, t, k);
    for (mpz_set_ui(i,2);  mpz_cmp(i,L1) <= 0;  mpz_add_ui(i,i,1)) {
      int m = moebius(i);
      if (m != 0) {
        mpz_pow_ui(t, i, k);
        mpz_fdiv_q(t, n, t);
        if (m > 0)  mpz_add(c, c, t);
        else        mpz_sub(c, c, t);
      }
    }
    mpz_clear(L1);
    for (c1i = 0;  mpz_cmp(i,nk) <= 0;  mpz_add_ui(i,i,1))
      c1i += moebius(i);
    mpz_set_si(c1, c1i);
  }
  mpz_set(count, n);
  mpz_add(count, count, c);
  mpz_add(count, count, c1);
  mpz_clear(nk);  mpz_clear(c);  mpz_clear(c1);  mpz_clear(i);  mpz_clear(t);
}



void nth_powerfree(mpz_t nth, const mpz_t n, uint32_t k)
{
  uint32_t tol = (k>=6) ? 7 : (k==5) ? 10 : (k==4) ? 20 : (k==3) ? 400 : 2000;
  unsigned long prec;
  int cmp;
  mpf_t fzm, fqk, t;
  mpz_t v, count, diff;

  if (k < 2 || mpz_sgn(n) <= 0) {
    mpz_set_ui(nth,0);
    return;
  }
  if (mpz_cmp_ui(n,4) < 0) {
    mpz_set(nth,n);
    return;
  }

  prec = mpz_sizeinbase(n,2) + 10;
  mpf_init2(fzm, prec);
  mpf_init2(fqk, prec);
  mpf_init2(t, prec);
  mpf_set_ui(fzm, k);
  (void) zetareal(fzm, prec);
  mpf_set_z(fqk, n);
  mpf_mul(fqk, fqk, fzm);

  mpz_init(v);
  mpz_init(count);
  mpz_init(diff);
  mpz_set_f(v, fqk);

  while (1) {
    powerfree_count(count, v, k);
    if (mpz_sizeinbase(n,2) <= 43) break;
    mpz_sub(diff, n, count);
    /* gmp_printf("qk %Zd  count %Zd  diff %Zd\n", v, count, diff); fflush(stdout); */
    if (mpz_cmpabs_ui(diff, tol) <= 0) break;

    /* Adjust difference for expected number of k-powerfree values */
    mpf_set_z(fqk,diff);
    mpf_mul(fqk,fqk,fzm);
    mpf_set_d(t, (mpz_sgn(diff) > 0) ? 0.5 : -0.5);
    mpf_add(fqk,fqk,t);
    mpz_set_f(diff,fqk);

    mpz_add(v, v, diff);
  }
  mpf_clear(t);  mpf_clear(fzm);  mpf_clear(fqk);

  while (!is_powerfree(v, k))
    mpz_sub_ui(v,v,1);
  cmp = mpz_cmp(n,count);
  while (cmp != 0 && mpz_cmp(n,count) != 0) {
    do {
      if (cmp > 0) mpz_add_ui(v,v,1);
      else         mpz_sub_ui(v,v,1);
    } while (!is_powerfree(v,k));
    if (cmp > 0) mpz_add_ui(count,count,1);
    else         mpz_sub_ui(count,count,1);
  }
  mpz_set(nth,v);
  mpz_clear(diff);  mpz_clear(count);  mpz_clear(v);
}
