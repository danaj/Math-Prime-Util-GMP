#include <math.h>
#include <gmp.h>
#include "ptypes.h"

/*****************************************************************************
 *
 *            The AKS polynomial-time deterministic primality test.
 *
 * Multiple variants are implemented.  The original algorithm as shown in
 * articles such as Rotella (2005) is not here as it is much slower than
 * the essentially identical algorithm in the updated (V6) AKS paper.  The
 * updated paper uses improved theorems based on Lenstra et al. which allows
 * loweing some limits used in the algorithm.
 *
 * All versions have a relatively similar O(log^{6.x}(n)) asymptotic growth,
 * which is what we expect from AKS.
 *
 * A version with improvements from Voloch and Bornemann is included, and
 * is similar to Bornemann's 2002 Pari/GP implementation.  It is *much* faster
 * than the V6 algorithm, and to the best of my knowledge is the fastest
 * publicly available AKS implementation in early 2016.
 *
 * The Bernstein 4.1 algorithm implements theorem 4.1 from Bernstein's 2003
 * paper, which has the Voloch improvements as well as many more.  It is
 * another 10-20x faster than the Bornemann version, hence is substantially
 * faster than any other known implementation in mid 2016.
 *
 * Copyright (2012-2016) Dana Jacobsen.
 *
 *****************************************************************************/


/* In approximate order of performance. */
#define AKS_VARIANT_V6        1    /* The V6 paper with Lenstra impr */
#define AKS_VARIANT_BERN21    2    /* AKS-Bernstein-Morain theorem 2.1 */
#define AKS_VARIANT_BERN22    3    /* AKS-Bernstein-Morain theorem 2.2 */
#define AKS_VARIANT_BERN23    4    /* AKS-Bernstein-Morain theorem 2.3 */
#define AKS_VARIANT_BORNEMANN 5    /* Based on Folkmar Bornemann's impl */
#define AKS_VARIANT_BERN41    6    /* Bernstein 2003, theorem 4.1 */

#define AKS_VARIANT  AKS_VARIANT_BERN41


#include "aks.h"
#include "prime_iterator.h"
#include "factor.h"

#if AKS_VARIANT == AKS_VARIANT_BORNEMANN
#define FUNC_mpz_logn 1
#endif
#define FUNC_mpz_log2 1
#include "utility.h"


static int test_anr(UV a, mpz_t n, UV r, mpz_t* px, mpz_t* py)
{
  int retval = 1;
  UV i, n_mod_r;
  mpz_t t;

  for (i = 0; i < r; i++)
    mpz_set_ui(px[i], 0);

  mpz_set_uv(px[0], a);
  mpz_set_ui(px[1], 1);

  poly_mod_pow(py, px, n, r, n);

  mpz_init(t);
  n_mod_r = mpz_fdiv_ui(n, r);
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

#if AKS_VARIANT != AKS_VARIANT_V6
static int is_primitive_root_uiprime(mpz_t n, UV r)
{
  int res;
  mpz_t zr;
  mpz_init(zr);
  mpz_set_uv(zr, r);
  res = is_primitive_root(n, zr, 1);
  mpz_clear(zr);
  return res;
}
#endif
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
#if AKS_VARIANT == AKS_VARIANT_BERN41
int bern41_acceptable(mpz_t n, UV r, UV s, mpz_t t1, mpz_t t2)
{
  double scmp = ceil(sqrt( (r-1)/3.0 )) * mpz_log2(n);
  UV d = (UV) (0.5 * (r-1));
  UV i = (UV) (0.475 * (r-1));
  UV j = i;
  /* Ensure conditions are correct */
  if (d > r-2)     d = r-2;
  if (i > d)       i = d;
  if (j > (r-2-d)) j = r-2-d;

  mpz_bin_uiui(t2, 2*s, i);
  mpz_bin_uiui(t1, d, i);       mpz_mul(t2, t2, t1);
  mpz_bin_uiui(t1, 2*s-i, j);   mpz_mul(t2, t2, t1);
  mpz_bin_uiui(t1, r-2-d, j);   mpz_mul(t2, t2, t1);
  return (mpz_log2(t2) >= scmp);
}
#endif



int is_aks_prime(mpz_t n)
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
    unsigned long limit, startr;
    PRIME_ITERATOR(iter);

    mpz_init(sqrtn);
    mpz_sqrt(sqrtn, n);

    log2n = mpz_log2(n);
    limit = (unsigned long) floor( log2n * log2n );

    if (_verbose>1) gmp_printf("# AKS checking order_r(%Zd) to %lu\n", n, limit);

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
    s = (unsigned long) floor( sqrt(r-1) * log2n );
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
    while (!is_primitive_root_uiprime(n,r))
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
      mpz_set_uv(t, r);
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
      if (!is_primitive_root_uiprime(n,r)) continue;
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
    unsigned long q, d, limit;
    unsigned long limit;
    double slim, scmp, sbin, x, log2n;
    mpz_t t, t2;
    PRIME_ITERATOR(iter);
    mpz_init(t);  mpz_init(t2);
    log2n = mpz_log2(n);
    limit = (unsigned long) floor( log2n * log2n );
    r = 2;
    s = 0;
    while (1) {
      /* todo: Check r|n and r >= sqrt(n) here instead of waiting */
      if (mpz_cmp_ui(n, r) <= 0) break;
      r++;
      unsigned long gcd = mpz_gcd_ui(NULL, n, r);
      if (gcd != 1) { mpz_clear(t); mpz_clear(t2); return 0; }
      unsigned long v = mpz_order_ui(r, n, limit);
      if (v >= limit) continue;

      mpz_set_ui(t2, r);
      totient(t, t2);
      q = mpz_get_ui(t);
      unsigned long phiv = q/v;
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
    /* Tuning: Initial 'r' selection */
    double const r0 = 0.008 * log2n * log2n;
    /* Tuning: Try a larger 'r' if 's' looks very large */
    UV const rmult = 8;
    UV slim;
    mpz_t tmp, tmp2;
    PRIME_ITERATOR(iter);

    mpz_init(tmp);  mpz_init(tmp2);
    /* r has to be at least 3. */
    prime_iterator_setprime(&iter, (r0 < 2) ? 2 : (UV) r0);
    r = prime_iterator_next(&iter);

    /* r must be a primitive root.  For performance, skip if s looks too big. */
    while ( !is_primitive_root_uiprime(n, r) ||
            !bern41_acceptable(n, r, rmult*(r-1), tmp, tmp2) )
      r = prime_iterator_next(&iter);
    prime_iterator_destroy(&iter);

    { /* Binary search for first s in [1,lim] where conditions met */
      UV bi = 1;
      UV bj = rmult * (r-1);
      while (bi < bj) {
        s = bi + (bj-bi)/2;
        if (!bern41_acceptable(n,r,s,tmp,tmp2))  bi = s+1;
        else                                     bj = s;
      }
      s = bj;
      /* Our S goes from 2 to s+1. */
      starta = 2;
      s = s+1;
    }
    /* printf("chose r=%lu s=%lu d = %lu i = %lu j = %lu\n", r, s, d, i, j); */

    /* Check divisibility to s(s-1) to cover both gcd conditions */
    slim = s * (s-1);
    if (_verbose > 1) printf("# aks trial to %"UVuf"\n", slim);
    if (_GMP_trial_factor(n, 2, slim) > 1)
      { mpz_clear(tmp); mpz_clear(tmp2); return 0; }
    /* If we checked divisibility to sqrt(n), then it is prime. */
    mpz_sqrt(tmp, n);
    /* TODO: slim UV / ui */
    if (mpz_cmp_ui(tmp, slim) <= 0)
      { mpz_clear(tmp); mpz_clear(tmp2); return 1; }

    /* Check b^(n-1) = 1 mod n for b in [2..s] */
    if (_verbose > 1) printf("# aks checking fermat to %"UVuf"\n", s);
    mpz_sub_ui(tmp2, n, 1);
    for (i = 2; i <= s; i++) {
      mpz_set_uv(tmp, i);
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
    retval = test_anr(a, n, r, px, py);
    if (!retval) break;
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
