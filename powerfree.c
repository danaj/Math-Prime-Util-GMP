#include <gmp.h>
#include "ptypes.h"

#include "powerfree.h"
#include "factor.h"
#include "real.h"
#include "prime_iterator.h"
#include "utility.h"
#define FUNC_isqrt 1
#define FUNC_ipow 1
#include "misc_ui.h"

int is_powerfree(mpz_t n, uint32_t k)
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
      if (exponents[i] >= k)
        ret = 0;
    clear_factors(nfactors, &factors, &exponents);
    return ret;
  }
}

void next_powerfree(mpz_t next, mpz_t n, uint32_t k)
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
void prev_powerfree(mpz_t prev, mpz_t n, uint32_t k)
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

#if 1
/*
 * Powerfree count really needs a range moebius for performance.
 */

void powerfree_count(mpz_t count, mpz_t n, uint32_t k)
{
#if 0
  mpz_t t, i, nk, c;

  if (k < 2 || mpz_sgn(n) <= 0) {
    mpz_set_ui(count, (mpz_cmp_ui(n,1) >= 0));
    return;
  }
  if (mpz_cmp_ui(n, 4) < 0) {
    mpz_set(count, n);
    return;
  }

  mpz_init(t);
  mpz_init(i);
  mpz_init(nk);
  mpz_root(nk, n, k);
  mpz_init_set(c, n);

  for (mpz_set_ui(i,2);  mpz_cmp(i,nk) <= 0;  mpz_add_ui(i,i,1)) {
    int m = moebius(i);
    if (m == 0) continue;
    mpz_pow_ui(t, i, k);
    mpz_fdiv_q(t, n, t);
    if (m > 0)  mpz_add(c, c, t);
    else        mpz_sub(c, c, t);
  }
  mpz_set(count, c);
  mpz_clear(c);  mpz_clear(nk);  mpz_clear(i);  mpz_clear(t);
#else
  mpz_t t, i, L1, nk, c;
  signed long c1;

  if (k < 2 || mpz_sgn(n) <= 0) {
    mpz_set_ui(count, (mpz_cmp_ui(n,1) >= 0));
    return;
  }
  if (mpz_cmp_ui(n, 4) < 0) {
    mpz_set(count, n);
    return;
  }

  mpz_init(t);
  mpz_init(i);
  mpz_init(L1);
  mpz_init(nk);
  mpz_init(c);

  mpz_fdiv_q_2exp(t, n, 1);
  mpz_root(L1, t, k);
  mpz_root(nk, n, k);
  mpz_init_set_ui(c, 0);

  for (mpz_set_ui(i,2);  mpz_cmp(i,L1) <= 0;  mpz_add_ui(i,i,1)) {
    int m = moebius(i);
    if (m != 0) {
      mpz_pow_ui(t, i, k);
      mpz_fdiv_q(t, n, t);
      if (m > 0)  mpz_add(c, c, t);
      else        mpz_sub(c, c, t);
    }
  }
  for (c1 = 0;  mpz_cmp(i,nk) <= 0;  mpz_add_ui(i,i,1))
    c1 += moebius(i);
  mpz_set(count, n);
  mpz_add(count, count, c);
  if (c1 >= 0) mpz_add_ui(count, count, (unsigned long) c1);
  else         mpz_sub_ui(count, count, (unsigned long) (-c1));
  mpz_clear(c);  mpz_clear(nk);  mpz_clear(L1);  mpz_clear(i);  mpz_clear(t);
#endif
}
#else

#define P_GT_LO(f,p,lo)  ( ((f)>=(lo)) ? (f) : (lo)+(((p)-((lo)%(p)))%(p)) )

/* Return a char array with lo-hi+1 elements. mu[k-lo] = µ(k) for k = lo .. hi.
 * It is the callers responsibility to call Safefree on the result. */
signed char* range_moebius(UV lo, UV hi)
{
  signed char* mu;
  UV p, i, sqrtn = isqrt(hi), count = hi-lo+1;

  /* Kuznetsov indicates that the Deléglise & Rivat (1996) method can be
   * modified to work on logs, which allows us to operate with no
   * intermediate memory at all.  Same time as the D&R method, less memory. */
  unsigned char logp;
  UV nextlog, nextlogi;

  if (hi < lo) croak("range_mobius error hi %"UVuf" < lo %"UVuf"\n", hi, lo);

  Newz(0, mu, count, signed char);
  if (sqrtn*sqrtn != hi && sqrtn < (UVCONST(1)<<(BITS_PER_WORD/2))-1) sqrtn++;

  {
    PRIME_ITERATOR(iter);

    logp = 1; nextlog = 3; /* 2+1 */
    for (p = 2; p <= sqrtn; p = prime_iterator_next(&iter)) {
      UV p2 = p*p;
      if (p > nextlog) {
        logp += 2;   /* logp is 1 | ceil(log(p)/log(2)) */
        nextlog = ((nextlog-1)*4)+1;
      }
      for (i = P_GT_LO(p, p, lo); i >= lo && i <= hi; i += p)
        mu[i-lo] += logp;
      for (i = P_GT_LO(p2, p2, lo); i >= lo && i <= hi; i += p2)
        mu[i-lo] = 0x80;
    }
    prime_iterator_destroy(&iter);
  }

  logp = log2_ui(lo);
  nextlogi = (UVCONST(2) << logp) - lo;
  for (i = 0; i < count; i++) {
    unsigned char a = mu[i];
    if (i >= nextlogi) nextlogi = (UVCONST(2) << ++logp) - lo;
    if (a & 0x80)       { a = 0; }
    else if (a >= logp) { a =  1 - 2*(a&1); }
    else                { a = -1 + 2*(a&1); }
    mu[i] = a;
  }
  if (lo == 0)  mu[0] = 0;

  return mu;
}

static UV squarefree_count_ui(UV n)
{
  signed char* mu;
  IV *M, *Mx, Mxisum, mert;
  UV sqrtn, I, D, i, j, S1 = 0, S2 = 0;

  if (n < 4) return n;

  sqrtn = isqrt(n);
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


void powerfree_count(mpz_t count, mpz_t n, uint32_t k)
{
  mpz_t t, i, L1, nk, c;
  signed long int c1;

  if (k < 2 || mpz_sgn(n) <= 0) {
    mpz_set_ui(count, (mpz_cmp_ui(n,1) >= 0));
    return;
  }
  if (mpz_cmp_ui(n, 4) < 0) {
    mpz_set(count, n);
    return;
  }
  if (mpz_fits_uv_p(n)) {
    UV c = powerfree_count_ui(mpz_get_uv(n), k);
    mpz_set_uv(count, c);
    return;
  }

  mpz_init(t);
  mpz_init(i);
  mpz_init(L1);
  mpz_init(nk);
  mpz_init(c);

  mpz_fdiv_q_2exp(t, n, 1);
  mpz_root(L1, t, k);
  mpz_root(nk, n, k);
  mpz_init_set_ui(c, 0);

  gmp_printf("   nk %Zd   L1 %Zd\n", nk, L1);
  if (0 && mpz_cmp_ui(nk, 1e8) < 0) {
    unsigned long int j, nku = mpz_get_ui(nk), L1u = mpz_get_ui(L1);
    signed char* mu = range_moebius(0, nku);
    for (j = 2; j <= L1u; j++) {
      if (mu[j] != 0) {
        mpz_ui_pow_ui(t, j, k);
        mpz_fdiv_q(t, n, t);
        if (mu[j] > 0) mpz_add(c, c, t);
        else           mpz_sub(c, c, t);
      }
    }
    for (c1 = 0; j <= nku; j++)
      c1 += mu[j];
    Safefree(mu);
  } else {
    for (mpz_set_ui(i,2);  mpz_cmp(i,L1) <= 0;  mpz_add_ui(i,i,1)) {
      int m = moebius(i);
      if (m != 0) {
        mpz_pow_ui(t, i, k);
        mpz_fdiv_q(t, n, t);
        if (m > 0)  mpz_add(c, c, t);
        else        mpz_sub(c, c, t);
      }
    }
    for (c1 = 0;  mpz_cmp(i,nk) <= 0;  mpz_add_ui(i,i,1))
      c1 += moebius(i);
  }
  mpz_set(count, n);
  mpz_add(count, count, c);
  if (c1 >= 0) mpz_add_ui(count, count, (unsigned long) c1);
  else         mpz_sub_ui(count, count, (unsigned long) (-c1));
  mpz_clear(c);  mpz_clear(nk);  mpz_clear(L1);  mpz_clear(i);  mpz_clear(t);
}
#endif

void nth_powerfree(mpz_t nth, mpz_t n, uint32_t k)
{
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
    if (mpz_sizeinbase(n,10) <= 16) break;
    mpz_sub(diff, n, count);
    /* gmp_printf("qk %Zd  count %Zd  diff %Zd\n", v, count, diff); fflush(stdout); */
    if (mpz_cmpabs_ui(diff, 20) <= 0) break;

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
