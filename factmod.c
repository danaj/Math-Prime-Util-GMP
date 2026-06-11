#include <gmp.h>
#include "ptypes.h"

#include "factmod.h"
#include "primality.h"
#include "prime_iterator.h"
#include "factor.h"
#include "utility.h"
#define FUNC_isqrt 1
#include "misc_ui.h"

#if BITS_PER_WORD == 32
#define FACTORIALMOD_MAX_D UV_MAX
#else
#define FACTORIALMOD_MAX_D UVCONST(500000000000)
#endif

void factorialmod(mpz_t r, const mpz_t n, const mpz_t m)
{
  mpz_t D, t, t2;
  UV du, i, p;
  int m_is_prime = -1, use_wilson = 0, check_for_zero = 0;

  if (mpz_sgn(n) < 0)
    croak("factorialmod: n must be >= 0");

  if (mpz_cmp_ui(m,1) <= 0 || mpz_cmp(n,m) >= 0) {
    mpz_set_ui(r,0);
    return;
  }

  mpz_init_set(D, n);
  mpz_init(t);
  mpz_tdiv_q_2exp(t, m, 1);
  if (mpz_cmp(t, n) < 0) {
    m_is_prime = _GMP_is_prime(m);
    if (m_is_prime) {
      mpz_sub(D, m, n);
      mpz_sub_ui(D, D, 1);
      use_wilson = 1;
    }
  }

  if (mpz_cmp_ui(D, 2) < 0) {
    if (use_wilson && mpz_sgn(D) == 0) mpz_sub_ui(r, m, 1);
    else                               mpz_set_ui(r, 1);
    mpz_clear(t); mpz_clear(D);
    return;
  }

  /* For large n, check for zero if m is composite and not too big.  This
   * must happen before the D -> UV conversion so huge but trivially zero
   * cases can return 0 instead of croaking. */
  if (mpz_cmp_ui(n, 500) > 0 && mpz_sizeinbase(m, 2) <= 150) {
    if (m_is_prime < 0) m_is_prime = _GMP_is_prime(m);
    if (!m_is_prime) check_for_zero = 1;
  }
  if (check_for_zero) {
    mpz_t *factors;
    int j, nfactors, *exponents, reszero;

    nfactors = factor(m, &factors, &exponents);
    /* Find max factor */
    mpz_set_ui(t, 0);
    for (j = 0; j < nfactors; j++) {
      if (exponents[j] > 1)
        mpz_mul_ui(factors[j], factors[j], exponents[j]);
      if (mpz_cmp(factors[j], t) > 0)
        mpz_set(t, factors[j]);
    }
    /* for m=p^k * p^k ..., t is max(p*k,p*k,...).  This is >= S(m), where
     * S(m) is the smallest value where m divides S(m)!.  Hence, every
     * n! mod m will be zero at that value or higher.  We could calculate
     * the exact value of S(m), then we would know there are no zero results
     * for the larger case. */
    reszero = (mpz_cmp(t, n) <= 0);
    clear_factors(nfactors, &factors, &exponents);
    if (reszero) { mpz_clear(t); mpz_clear(D); mpz_set_ui(r,0); return; }
  }

  /* Set du = D = reduced n if small enough, otherwise error. */
  if (!mpz_fits_uv_p(D) || (du = mpz_get_uv(D)) >= FACTORIALMOD_MAX_D) {
    mpz_clear(t); mpz_clear(D);
    croak("factorialmod: reduced n too large for processing");
  }

  /* Accumulate into t, then mod into r at the end. */
  mpz_set_ui(t,1);

  /* For small D, naive method. */
  if (du <= 1000) {
    for (i = 2; i <= du && mpz_sgn(t); i++) {
      mpz_mul_ui(t, t, i);
      if ((i & 15) == 0) mpz_mod(t, t, m);
    }
  } else {
    UV j, sd = isqrt(du);
    PRIME_ITERATOR(iter);

    mpz_init(t2);
    mpz_set_ui(t,1);
    /* Group into powers of primes */
    for (p = 2, i = 0; p <= du/sd; p = prime_iterator_next(&iter)) {
      UV td = du/p,  e = td;
      do { td /= p; e += td; } while (td > 0);
      mpz_set_ui(t2, p);
      mpz_powm_ui(t2, t2, e, m);
      mpz_mul(t, t, t2);
      if ((i++ & 15) == 0) {
        mpz_mod(t, t, m);
        if (!mpz_sgn(t)) break;
      }
    }
    /* Further group by primes with the same power. */
    for (j = sd-1; j >= 1 && mpz_sgn(t); j--) {
      UV lo = du / (j+1)+1,  hi = du / j;
      MPUassert(p >= lo, "factorialmod prime loop p should be in range");
      /* while (p < lo) p = prime_iterator_next(&iter); */
      for (mpz_set_ui(t2,1), i=0;  p <= hi;  p = prime_iterator_next(&iter)) {
        mpz_mul_ui(t2, t2, p);
        if ((i++ & 15) == 0) mpz_mod(t2, t2, m);
      }
      mpz_powm_ui(t2, t2, j, m);
      mpz_mul(t, t, t2);
      if ((j & 15) == 0) mpz_mod(t, t, m);
    }
    mpz_clear(t2);
    prime_iterator_destroy(&iter);
  }
  mpz_mod(r, t, m);
  mpz_clear(t);

  /* If we used Wilson's theorem, turn the result for D! into N! */
  if (use_wilson && mpz_sgn(r)) {
    if (!(du&1)) mpz_sub(r, m, r);
    mpz_invert(r, r, m);
  }
  mpz_clear(D);
}

/* Compute C(n,k) modulo p^e, where k fits in an unsigned long.
 * This avoids creating a potentially enormous intermediate binomial. */
#define BINOMIALMOD_SMALLK_MAX  UVCONST(100000)
#define BINOMIALMOD_BIN_UIUI_MAX  UVCONST(1000000)
#define BINOMIALMOD_FULL_UIUI_MAX  UVCONST(100000)
static void _binomialmod_native_digit_prime(mpz_t r, unsigned long n, unsigned long k, const mpz_t p)
{
  unsigned long i;
  mpz_t num, den, inv;

  if (k > n) {
    mpz_set_ui(r, 0);
    return;
  }
  if (k == 0 || k == n) {
    mpz_set_ui(r, 1);
    return;
  }
  if (k > n/2)
    k = n-k;
  if (k == 1) {
    mpz_set_ui(r, n);
    return;
  }

  mpz_init_set_ui(num, 1);
  mpz_init_set_ui(den, 1);
  mpz_init(inv);

  for (i = 1; i <= k; i++) {
    mpz_mul_ui(num, num, n-k+i);
    mpz_mod(num, num, p);
    mpz_mul_ui(den, den, i);
    mpz_mod(den, den, p);
  }
  if (!mpz_invert(inv, den, p))
    croak("binomialmod: internal digit inverse failure");
  mpz_mul(r, num, inv);
  mpz_mod(r, r, p);

  mpz_clear(inv);
  mpz_clear(den);
  mpz_clear(num);
}

static void _binomialmod_smallk_prime_power(mpz_t r, const mpz_t n, unsigned long k, const mpz_t p, unsigned int e, const mpz_t pk)
{
  unsigned long i;
  long v = 0;
  mpz_t nmk, t, num, den, inv;

  mpz_init(nmk);
  mpz_init(t);
  mpz_init_set_ui(num, 1);
  mpz_init_set_ui(den, 1);
  mpz_init(inv);

  mpz_sub_ui(nmk, n, k);
  for (i = 1; i <= k; i++) {
    mpz_add_ui(t, nmk, i);
    v += (long) mpz_remove(t, t, p);
    mpz_mod(t, t, pk);
    mpz_mul(num, num, t);
    mpz_mod(num, num, pk);

    mpz_set_ui(t, i);
    v -= (long) mpz_remove(t, t, p);
    mpz_mod(t, t, pk);
    mpz_mul(den, den, t);
    mpz_mod(den, den, pk);
  }

  if (v >= (long)e) {
    mpz_set_ui(r, 0);
  } else {
    if (!mpz_invert(inv, den, pk))
      croak("binomialmod: internal inverse failure");
    mpz_mul(r, num, inv);
    mpz_mod(r, r, pk);
    if (v > 0) {
      mpz_powm_ui(t, p, (unsigned long)v, pk);
      mpz_mul(r, r, t);
      mpz_mod(r, r, pk);
    }
  }

  mpz_clear(inv);
  mpz_clear(den);
  mpz_clear(num);
  mpz_clear(t);
  mpz_clear(nmk);
}

static void _binomialmod_lucas_prime(mpz_t r, const mpz_t n, unsigned long k, const mpz_t p)
{
  mpz_t N, np, term;
  unsigned long pu = 0;

  mpz_init_set(N, n);
  mpz_init(np);
  mpz_init(term);
  mpz_set_ui(r, 1);

  if (mpz_fits_ulong_p(p))
    pu = mpz_get_ui(p);

  while (k > 0 && mpz_sgn(r) != 0) {
    unsigned long kp;

    if (pu) {
      kp = k % pu;
      k /= pu;
    } else {
      kp = k;
      k = 0;
    }

    mpz_mod(np, N, p);
    if (mpz_cmp_ui(np, kp) < 0) {
      mpz_set_ui(r, 0);
      break;
    }

    if (kp == 1) {
      mpz_set(term, np);
    } else if (mpz_fits_ulong_p(np) && mpz_get_ui(np) <= BINOMIALMOD_BIN_UIUI_MAX) {
      mpz_bin_uiui(term, mpz_get_ui(np), kp);
      mpz_mod(term, term, p);
    } else if (kp <= BINOMIALMOD_SMALLK_MAX && mpz_fits_ulong_p(np)) {
      _binomialmod_native_digit_prime(term, mpz_get_ui(np), kp, p);
    } else if (kp > 1) {
      if (kp <= BINOMIALMOD_SMALLK_MAX) {
        _binomialmod_smallk_prime_power(term, np, kp, p, 1, p);
      } else {
        mpz_bin_ui(term, np, kp);
        mpz_mod(term, term, p);
      }
    } else {
      mpz_set_ui(term, 1);
    }

    mpz_mul(r, r, term);
    mpz_mod(r, r, p);

    if (pu)
      mpz_tdiv_q_ui(N, N, pu);
  }

  mpz_clear(term);
  mpz_clear(np);
  mpz_clear(N);
}

static void _factorial_without_prime_square(mpz_t r, const mpz_t n, UV p, const mpz_t pk)
{
  mpz_t a, t, pm1_fact, b_fact, hb, term;
  UV b, j;

  if (mpz_cmp_ui(n, 1) <= 0) {
    mpz_set_ui(r, 1);
    return;
  }
  if (mpz_cmp_ui(n, p) < 0) {
    factorialmod(r, n, pk);
    return;
  }

  mpz_init(a);
  mpz_init(t);
  mpz_init(pm1_fact);
  mpz_init(b_fact);
  mpz_init_set_ui(hb, 0);
  mpz_init(term);

  b = mpz_tdiv_q_ui(a, n, p);
  mpz_set_ui(t, p-1);
  factorialmod(pm1_fact, t, pk);
  mpz_set_ui(t, b);
  factorialmod(b_fact, t, pk);

  for (j = 1; j <= b; j++) {
    mpz_set_ui(t, j);
    if (!mpz_invert(term, t, pk))
      croak("binomialmod: internal harmonic inverse failure");
    mpz_add(hb, hb, term);
    mpz_mod_ui(hb, hb, p);
  }

  mpz_powm(r, pm1_fact, a, pk);
  mpz_mul(r, r, b_fact);
  mpz_mod(r, r, pk);

  if (mpz_sgn(a) > 0 && mpz_sgn(hb) > 0) {
    mpz_mul_ui(term, a, p);
    mpz_mul(term, term, hb);
    mpz_add_ui(term, term, 1);
    mpz_mod(term, term, pk);
    mpz_mul(r, r, term);
    mpz_mod(r, r, pk);
  }

  mpz_clear(term);
  mpz_clear(hb);
  mpz_clear(b_fact);
  mpz_clear(pm1_fact);
  mpz_clear(t);
  mpz_clear(a);
}

static void _binomialmod_lucas_prime_square(mpz_t r, const mpz_t n, unsigned long k, const mpz_t pz, const mpz_t pk)
{
  mpz_t Ntmp, Ktmp, Rtmp, NminusK, prq, x, y, z, den, inv, pterm;
  UV p, i, d, rq;
  UV *carry;

  p = mpz_get_uv(pz);
  mpz_init_set(Ntmp, n);
  mpz_init_set_ui(Ktmp, k);
  mpz_init(Rtmp);
  mpz_init(NminusK);
  mpz_sub(NminusK, n, Ktmp);
  mpz_init(x);

  for (d = 0, mpz_set(x, n); mpz_sgn(x) > 0; d++)
    mpz_tdiv_q_ui(x, x, p);
  Newz(0, carry, d+1, UV);

  for (i = 0; i < d; i++) {
    UV ni = mpz_tdiv_q_ui(Ntmp, Ntmp, p);
    UV ki = mpz_tdiv_q_ui(Ktmp, Ktmp, p);
    carry[i] = (ni < ki + (i > 0 ? carry[i-1] : 0));
  }
  for (i = d; i-- > 1; )
    carry[i-1] += carry[i];

  if (carry[0] >= 2) {
    mpz_set_ui(r, 0);
    Safefree(carry);
    mpz_clear(x); mpz_clear(NminusK); mpz_clear(Rtmp); mpz_clear(Ktmp); mpz_clear(Ntmp);
    return;
  }

  rq = 2 - carry[0];
  mpz_init(prq);
  mpz_pow_ui(prq, pz, rq);
  mpz_init(y);
  mpz_init(z);
  mpz_init(den);
  mpz_init(inv);
  mpz_init(pterm);

  if ((rq < d) && (carry[rq-1] & 1))
    mpz_sub_ui(r, pk, 1);
  else
    mpz_set_ui(r, 1);
  if (carry[0]) {
    mpz_powm_ui(pterm, pz, carry[0], pk);
    mpz_mul(r, r, pterm);
    mpz_mod(r, r, pk);
  }

  mpz_set(Ntmp, n);
  mpz_set_ui(Ktmp, k);
  mpz_set(Rtmp, NminusK);
  for (i = 0; i < d && mpz_sgn(r) != 0; i++) {
    mpz_mod(x, Ntmp, prq);
    mpz_mod(y, Ktmp, prq);
    mpz_mod(z, Rtmp, prq);

    if (rq == 2) {
      _factorial_without_prime_square(x, x, p, pk);
      _factorial_without_prime_square(y, y, p, pk);
      _factorial_without_prime_square(z, z, p, pk);
    } else {
      factorialmod(x, x, pk);
      factorialmod(y, y, pk);
      factorialmod(z, z, pk);
    }

    mpz_mul(den, y, z);
    mpz_mod(den, den, pk);
    if (!mpz_invert(inv, den, pk))
      croak("binomialmod: internal prime-square inverse failure");
    mpz_mul(x, x, inv);
    mpz_mod(x, x, pk);
    mpz_mul(r, r, x);
    mpz_mod(r, r, pk);

    mpz_tdiv_q_ui(Ntmp, Ntmp, p);
    mpz_tdiv_q_ui(Ktmp, Ktmp, p);
    mpz_tdiv_q_ui(Rtmp, Rtmp, p);
  }

  mpz_clear(pterm);
  mpz_clear(inv);
  mpz_clear(den);
  mpz_clear(z);
  mpz_clear(y);
  mpz_clear(prq);
  Safefree(carry);
  mpz_clear(x);
  mpz_clear(NminusK);
  mpz_clear(Rtmp);
  mpz_clear(Ktmp);
  mpz_clear(Ntmp);
}

static void _poly_init(mpz_t **pa, unsigned int e)
{
  unsigned int i;
  mpz_t *a;
  New(0, a, e, mpz_t);
  for (i = 0; i < e; i++)
    mpz_init(a[i]);
  *pa = a;
}

static void _poly_clear(mpz_t **pa, unsigned int e)
{
  unsigned int i;
  for (i = 0; i < e; i++)
    mpz_clear((*pa)[i]);
  Safefree(*pa);
}

static void _poly_mul_mod(mpz_t *C, mpz_t *A, mpz_t *B, unsigned int e, const mpz_t mod)
{
  unsigned int i, j;
  mpz_t t;

  mpz_init(t);
  for (i = 0; i < e; i++)
    mpz_set_ui(C[i], 0);

  for (i = 0; i < e; i++) {
    if (mpz_sgn(A[i]) == 0) continue;
    for (j = 0; j < e-i; j++) {
      if (mpz_sgn(B[j]) == 0) continue;
      mpz_mul(t, A[i], B[j]);
      mpz_add(C[i+j], C[i+j], t);
      mpz_mod(C[i+j], C[i+j], mod);
    }
  }
  mpz_clear(t);
}

static void _poly_shift_mod(mpz_t *B, mpz_t *A, const mpz_t h, unsigned int e, const mpz_t mod)
{
  int j, m;
  mpz_t hpow, bin, term;

  mpz_init(hpow);
  mpz_init(bin);
  mpz_init(term);
  for (j = 0; j < (int)e; j++)
    mpz_set_ui(B[j], 0);

  for (j = 0; j < (int)e; j++) {
    if (mpz_sgn(A[j]) == 0) continue;
    mpz_set_ui(hpow, 1);
    for (m = j; m >= 0; m--) {
      mpz_bin_uiui(bin, j, m);
      mpz_mul(term, bin, hpow);
      mpz_mod(term, term, mod);
      mpz_mul(term, term, A[j]);
      mpz_add(B[m], B[m], term);
      mpz_mod(B[m], B[m], mod);
      if (m > 0) {
        mpz_mul(hpow, hpow, h);
        mpz_mod(hpow, hpow, mod);
      }
    }
  }

  mpz_clear(term);
  mpz_clear(bin);
  mpz_clear(hpow);
}

static void _poly_build_mod(mpz_t *R, const mpz_t q, mpz_t *poly, unsigned int e, const mpz_t mod)
{
  unsigned int i;
  mpz_t h, two_h;
  mpz_t *P_h, *P_h_shifted, *P_2h, *poly_shifted;

  if (mpz_sgn(q) == 0) {
    mpz_set_ui(R[0], 1);
    for (i = 1; i < e; i++)
      mpz_set_ui(R[i], 0);
    return;
  }
  if (mpz_cmp_ui(q, 1) == 0) {
    for (i = 0; i < e; i++)
      mpz_set(R[i], poly[i]);
    return;
  }

  mpz_init(h);
  mpz_init(two_h);
  mpz_fdiv_q_2exp(h, q, 1);
  _poly_init(&P_h, e);
  _poly_init(&P_h_shifted, e);
  _poly_init(&P_2h, e);

  _poly_build_mod(P_h, h, poly, e, mod);
  _poly_shift_mod(P_h_shifted, P_h, h, e, mod);
  _poly_mul_mod(P_2h, P_h, P_h_shifted, e, mod);

  if (mpz_odd_p(q)) {
    _poly_init(&poly_shifted, e);
    mpz_mul_2exp(two_h, h, 1);
    _poly_shift_mod(poly_shifted, poly, two_h, e, mod);
    _poly_mul_mod(R, P_2h, poly_shifted, e, mod);
    _poly_clear(&poly_shifted, e);
  } else {
    for (i = 0; i < e; i++)
      mpz_set(R[i], P_2h[i]);
  }

  _poly_clear(&P_2h, e);
  _poly_clear(&P_h_shifted, e);
  _poly_clear(&P_h, e);
  mpz_clear(two_h);
  mpz_clear(h);
}

static void _factorial_without_prime_pe(mpz_t r, const mpz_t n, UV p, unsigned int e, const mpz_t pk)
{
  unsigned int idx;
  UV j, rem;
  mpz_t fact, inv, t, ppow, q, pq;
  mpz_t *c, *poly, *Pq;

  if (mpz_cmp_ui(n, 1) <= 0) {
    mpz_set_ui(r, 1);
    return;
  }
  if (mpz_cmp_ui(n, p) < 0) {
    factorialmod(r, n, pk);
    return;
  }

  _poly_init(&c, e);
  _poly_init(&poly, e);
  _poly_init(&Pq, e);
  mpz_init_set_ui(fact, 1);
  mpz_init(inv);
  mpz_init(t);
  mpz_init_set_ui(ppow, 1);
  mpz_init(q);
  mpz_init(pq);
  mpz_set_ui(c[0], 1);

  for (j = 1; j < p; j++) {
    mpz_mul_ui(fact, fact, j);
    mpz_mod(fact, fact, pk);
    mpz_set_ui(t, j);
    if (!mpz_invert(inv, t, pk))
      croak("binomialmod: internal polynomial inverse failure");
    for (idx = e; idx-- > 1; ) {
      mpz_mul(t, c[idx-1], inv);
      mpz_add(c[idx], c[idx], t);
      mpz_mod(c[idx], c[idx], pk);
    }
  }

  for (idx = 0; idx < e; idx++) {
    mpz_mul(poly[idx], c[idx], fact);
    mpz_mod(poly[idx], poly[idx], pk);
    mpz_mul(poly[idx], poly[idx], ppow);
    mpz_mod(poly[idx], poly[idx], pk);
    mpz_mul_ui(ppow, ppow, p);
    mpz_mod(ppow, ppow, pk);
  }

  rem = mpz_tdiv_q_ui(q, n, p);
  _poly_build_mod(Pq, q, poly, e, pk);
  mpz_set(r, Pq[0]);

  mpz_mul_ui(pq, q, p);
  for (j = 1; j <= rem; j++) {
    mpz_add_ui(t, pq, j);
    mpz_mul(r, r, t);
    mpz_mod(r, r, pk);
  }

  mpz_clear(pq);
  mpz_clear(q);
  mpz_clear(ppow);
  mpz_clear(t);
  mpz_clear(inv);
  mpz_clear(fact);
  _poly_clear(&Pq, e);
  _poly_clear(&poly, e);
  _poly_clear(&c, e);
}

static void _binomialmod_lucas_prime_power(mpz_t r, const mpz_t n, unsigned long k, const mpz_t pz, unsigned int q, const mpz_t pk)
{
  mpz_t Ntmp, Ktmp, Rtmp, NminusK, prq, x, y, z, den, inv, pterm;
  UV p, i, d, rq;
  UV *carry;

  p = mpz_get_uv(pz);
  mpz_init_set(Ntmp, n);
  mpz_init_set_ui(Ktmp, k);
  mpz_init(Rtmp);
  mpz_init(NminusK);
  mpz_sub(NminusK, n, Ktmp);
  mpz_init(x);

  for (d = 0, mpz_set(x, n); mpz_sgn(x) > 0; d++)
    mpz_tdiv_q_ui(x, x, p);
  Newz(0, carry, d+1, UV);

  for (i = 0; i < d; i++) {
    UV ni = mpz_tdiv_q_ui(Ntmp, Ntmp, p);
    UV ki = mpz_tdiv_q_ui(Ktmp, Ktmp, p);
    carry[i] = (ni < ki + (i > 0 ? carry[i-1] : 0));
  }
  for (i = d; i-- > 1; )
    carry[i-1] += carry[i];

  if (carry[0] >= q) {
    mpz_set_ui(r, 0);
    Safefree(carry);
    mpz_clear(x); mpz_clear(NminusK); mpz_clear(Rtmp); mpz_clear(Ktmp); mpz_clear(Ntmp);
    return;
  }

  rq = q - carry[0];
  mpz_init(prq);
  mpz_pow_ui(prq, pz, rq);
  mpz_init(y);
  mpz_init(z);
  mpz_init(den);
  mpz_init(inv);
  mpz_init(pterm);

  if ((p > 2 || rq < 3) && rq <= d && (carry[rq-1] & 1))
    mpz_sub_ui(r, pk, 1);
  else
    mpz_set_ui(r, 1);
  if (carry[0]) {
    mpz_powm_ui(pterm, pz, carry[0], pk);
    mpz_mul(r, r, pterm);
    mpz_mod(r, r, pk);
  }

  mpz_set(Ntmp, n);
  mpz_set_ui(Ktmp, k);
  mpz_set(Rtmp, NminusK);
  for (i = 0; i < d && mpz_sgn(r) != 0; i++) {
    mpz_mod(x, Ntmp, prq);
    mpz_mod(y, Ktmp, prq);
    mpz_mod(z, Rtmp, prq);

    _factorial_without_prime_pe(x, x, p, q, pk);
    _factorial_without_prime_pe(y, y, p, q, pk);
    _factorial_without_prime_pe(z, z, p, q, pk);

    mpz_mul(den, y, z);
    mpz_mod(den, den, pk);
    if (!mpz_invert(inv, den, pk))
      croak("binomialmod: internal prime-power inverse failure");
    mpz_mul(x, x, inv);
    mpz_mod(x, x, pk);
    mpz_mul(r, r, x);
    mpz_mod(r, r, pk);

    mpz_tdiv_q_ui(Ntmp, Ntmp, p);
    mpz_tdiv_q_ui(Ktmp, Ktmp, p);
    mpz_tdiv_q_ui(Rtmp, Rtmp, p);
  }

  mpz_clear(pterm);
  mpz_clear(inv);
  mpz_clear(den);
  mpz_clear(z);
  mpz_clear(y);
  mpz_clear(prq);
  Safefree(carry);
  mpz_clear(x);
  mpz_clear(NminusK);
  mpz_clear(Rtmp);
  mpz_clear(Ktmp);
  mpz_clear(Ntmp);
}

static void _binomialmod_positive(mpz_t r, const mpz_t n, unsigned long k, const mpz_t m)
{
  mpz_t *factors, *mods, *residues, lcm, fullbinom;
  int i, nfactors, *exponents, use_full = 0, ninit = 0;

  if (k == 0) {
    mpz_set_ui(r, 1);
    return;
  }
  if (k == 1) {
    mpz_mod(r, n, m);
    return;
  }
  if (mpz_fits_ulong_p(n) && mpz_get_ui(n) <= BINOMIALMOD_FULL_UIUI_MAX) {
    mpz_bin_uiui(r, mpz_get_ui(n), k);
    mpz_mod(r, r, m);
    return;
  }
  if (_GMP_is_prime(m)) {
    _binomialmod_lucas_prime(r, n, k, m);
    return;
  }

  nfactors = factor(m, &factors, &exponents);
  New(0, mods, nfactors, mpz_t);
  New(0, residues, nfactors, mpz_t);

  for (i = 0; i < nfactors; i++) {
    mpz_init(mods[i]);
    mpz_init(residues[i]);
    ninit++;
    mpz_pow_ui(mods[i], factors[i], exponents[i]);

    if (exponents[i] == 1) {
      _binomialmod_lucas_prime(residues[i], n, k, factors[i]);
    } else if (exponents[i] == 2 && mpz_fits_uv_p(factors[i]) && mpz_get_uv(factors[i]) > 2 && mpz_get_uv(factors[i]) <= UVCONST(5000000)) {
      _binomialmod_lucas_prime_square(residues[i], n, k, factors[i], mods[i]);
    } else if (exponents[i] <= 4 && mpz_fits_uv_p(factors[i]) && mpz_get_uv(factors[i]) > 2 && mpz_get_uv(factors[i]) <= UVCONST(5000000)) {
      _binomialmod_lucas_prime_power(residues[i], n, k, factors[i], exponents[i], mods[i]);
    } else if (k <= BINOMIALMOD_SMALLK_MAX) {
      _binomialmod_smallk_prime_power(residues[i], n, k, factors[i], exponents[i], mods[i]);
    } else {
      use_full = 1;
      break;
    }
  }

  if (use_full) {
    mpz_init(fullbinom);
    if (mpz_fits_ulong_p(n))
      mpz_bin_uiui(fullbinom, mpz_get_ui(n), k);
    else
      mpz_bin_ui(fullbinom, n, k);
    mpz_mod(r, fullbinom, m);
    mpz_clear(fullbinom);
  } else {
    mpz_init(lcm);
    if (!chinese(r, lcm, residues, mods, nfactors))
      croak("binomialmod: internal CRT failure");
    mpz_mod(r, r, m);
    mpz_clear(lcm);
  }

  for (i = 0; i < ninit; i++) {
    mpz_clear(residues[i]);
    mpz_clear(mods[i]);
  }
  Safefree(residues);
  Safefree(mods);
  clear_factors(nfactors, &factors, &exponents);
}

void binomialmod(mpz_t r, const mpz_t n, const mpz_t k, const mpz_t m)
{
  mpz_t N, K, T;
  unsigned long ku;
  int neg = 0;

  if (mpz_cmp_ui(m, 1) <= 0) {
    mpz_set_ui(r, 0);
    return;
  }

  mpz_init(N);
  mpz_init(K);
  mpz_init(T);

  if (mpz_sgn(n) >= 0) {
    if (mpz_sgn(k) < 0 || mpz_cmp(k,n) > 0) {
      mpz_set_ui(r, 0);
      mpz_clear(T); mpz_clear(K); mpz_clear(N);
      return;
    }
    mpz_set(N, n);
    mpz_set(K, k);
  } else if (mpz_sgn(k) >= 0) {
    mpz_neg(N, n);
    mpz_add(N, N, k);
    mpz_sub_ui(N, N, 1);
    mpz_set(K, k);
    neg = mpz_odd_p(k);
  } else {
    if (mpz_cmp(k,n) > 0) {
      mpz_set_ui(r, 0);
      mpz_clear(T); mpz_clear(K); mpz_clear(N);
      return;
    }
    mpz_neg(N, k);
    mpz_sub_ui(N, N, 1);
    mpz_sub(K, n, k);
    neg = mpz_odd_p(K);
  }

  mpz_sub(T, N, K);
  if (mpz_cmp(T, K) < 0)
    mpz_swap(T, K);

  if (!mpz_fits_ulong_p(K)) {
    mpz_clear(T); mpz_clear(K); mpz_clear(N);
    croak("binomialmod: k too large");
  }
  ku = mpz_get_ui(K);

  _binomialmod_positive(r, N, ku, m);
  if (neg && mpz_sgn(r) != 0)
    mpz_sub(r, m, r);

  mpz_clear(T);
  mpz_clear(K);
  mpz_clear(N);
}
