#include <gmp.h>
#include "ptypes.h"

#include "powerfree.h"
#include "utility.h"
#include "factor.h"
#include "real.h"

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
