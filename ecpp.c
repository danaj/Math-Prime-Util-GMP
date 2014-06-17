/*****************************************************************************
 *
 * ECPP - Elliptic Curve Primality Proving
 *
 * Copyright (c) 2013-2014 Dana Jacobsen (dana@acm.org).
 * This is free software; you can redistribute it and/or modify it under
 * the same terms as the Perl 5 programming language system itself.
 *
 * This file is part of the Math::Prime::Util::GMP Perl module.  A script
 * is included to build this as a standalone program (see the README file).
 *
 * This is pretty good for numbers less than 800 digits.  Over that, it needs
 * larger discriminant sets.  Comparing to other contemporary software:
 *
 *   - Primo is much faster for inputs over 300 digits.  Not open source.
 *   - mpz_aprcl 1.1 (APR-CL).  Nearly the same speed to ~600 digits, with
 *     very little speed variation.  Faster over 800 digits.  No certificate.
 *   - GMP-ECPP is much slower at all sizes, and nearly useless > 300 digits.
 *   - AKS is stupendously slow, even with Bernstein and Voloch improvements.
 *   - François Morain's 10-20 year old work describes optimizations not
 *     present here, but his (very old!) binaries run slower than this code at
 *     all sizes.  Not open source.
 *
 * A set of fixed discriminants are used, rather than calculating them as
 * needed.  Having a way to calculate values as needed would be a big help.
 * In the interests of space for the MPU package, I've chosen ~600 values which
 * compile into about 35k of data.  This is about 1/5 of the entire code size
 * for the MPU package.  The github repository includes an expanded set of 5271
 * discriminants that compile to 2MB.  This is recommended if proving 300+
 * digit numbers is a regular occurance.  There is a set available for download
 * with almost 15k polys, taking 15.5MB.
 *
 * This version uses the FAS "factor all strategy", meaning it first constructs
 * the entire factor chain, with backtracking if necessary, then will do the
 * elliptic curve proof as it recurses back.
 *
 * If your goal is primality proofs for very large numbers, use Primo.  It's
 * free, it is very fast, it is widely used, it can process batch results,
 * and it makes independently verifiable certificates (including the verifier
 * included in this package).  MPU's ECPP (this software) is an open source
 * alternative with many of the same features for "small" numbers of <1000
 * digits.  Improvements are possible since it is open source.
 *
 * Another open source alternative if one does not need certificates is the
 * mpz_aprcl code from David Cleaver.  To about 600 digits the speeds are
 * very similar, but past that this ECPP code starts slowing down.
 *
 * Thanks to H. Cohen, R. Crandall & C. Pomerance, and H. Riesel for their
 * text books.  Thanks to the authors of open source software who allow me
 * to compare and contrast (GMP-ECM, GMP-ECPP).  Thanks to the authors of GMP.
 * Thanks to Schoof, Goldwasser, Kilian, Atkin, Morain, Lenstra, etc. for all
 * the math and publications.  Thanks to Gauss, Euler, et al.
 *
 * The ECM code in ecm.c was heavily influenced by early GMP-ECM work by Phil
 * Zimmerman, as well as all the articles of Montgomery, Bosma, Lentra,
 * Cohen, and others.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

#include "ptypes.h"
#include "ecpp.h"
#include "gmp_main.h"  /* is_prob_prime, pminus1_factor, miller_rabin_random */
#include "ecm.h"
#include "utility.h"
#include "prime_iterator.h"
#include "bls75.h"

#define MAX_SFACS 1000

#ifdef USE_LIBECM
 #include <ecm.h>
#endif

#ifdef USE_APRCL
 #include "mpz_aprcl.h"
 #include "mpz_aprcl.c"
#endif

/*********** big primorials and lcm for divisibility tests  **********/
static int _gcdinit = 0;
static mpz_t _gcd_small;
static mpz_t _gcd_large;

void init_ecpp_gcds(UV nsize) {
  if (_gcdinit == 0) {
    mpz_init(_gcd_small);
    mpz_init(_gcd_large);
    _GMP_pn_primorial(_gcd_small,  3000);
    /* This is never re-adjusted -- first number proved sets the size */
    nsize *= 20;
    if      (nsize < 20000) nsize = 20000;
    else if (nsize > 500000) nsize = 500000;
    _GMP_pn_primorial(_gcd_large, nsize);
    mpz_divexact(_gcd_large, _gcd_large, _gcd_small);
    mpz_divexact_ui(_gcd_small, _gcd_small, 2*3*5);
    _gcdinit = 1;
  }
}

void destroy_ecpp_gcds(void) {
  if (!_gcdinit) return;
  mpz_clear(_gcd_small);
  mpz_clear(_gcd_large);
  _gcdinit = 0;
}

/* We could use a function with a prefilter here, but my tests are showing
 * that adding a Fermat test (ala GMP's is_probab_prime) is slower than going
 * straight to the base-2 Miller-Rabin test we use in BPSW. */
#define is_bpsw_prime(n) _GMP_BPSW(n)

static int check_for_factor(mpz_t f, mpz_t inputn, mpz_t fmin, mpz_t n, int stage, mpz_t* sfacs, int* nsfacs, int degree)
{
  int success, sfaci;
  UV B1;

  /* Use this so we don't modify their input value */
  mpz_set(n, inputn);

  if (mpz_cmp(n, fmin) <= 0) return 0;

#if 0
  /* Use this to really encourage n-1 / n+1 proof types */
  if (degree <= 0) {
    if (stage == 1) return -1;
    if (stage == 0) stage = 1;
  }
#endif

  /* Utilize GMP's fast gcd algorithms.  Trial to 224737+ with two gcds. */
  mpz_tdiv_q_2exp(n, n, mpz_scan1(n, 0));
  while (mpz_divisible_ui_p(n, 3))  mpz_divexact_ui(n, n, 3);
  while (mpz_divisible_ui_p(n, 5))  mpz_divexact_ui(n, n, 5);
  if (mpz_cmp(n, fmin) <= 0) return 0;
  mpz_gcd(f, n, _gcd_small);
  while (mpz_cmp_ui(f, 1) > 0) {
    mpz_divexact(n, n, f);
    mpz_gcd(f, f, n);
  }
  if (mpz_cmp(n, fmin) <= 0) return 0;
  mpz_gcd(f, n, _gcd_large);
  while (mpz_cmp_ui(f, 1) > 0) {
    mpz_divexact(n, n, f);
    mpz_gcd(f, f, n);
  }

  sfaci = 0;
  success = 1;
  while (success) {
    UV nsize = mpz_sizeinbase(n, 2);
    const int do_pm1 = 1;
    const int do_pp1 = 1;
    const int do_pbr = 0;
    const int do_ecm = 0;

    if (mpz_cmp(n, fmin) <= 0) return 0;
    if (is_bpsw_prime(n)) { mpz_set(f, n); return (mpz_cmp(f, fmin) > 0); }

    success = 0;
    B1 = 300 + 3 * nsize;
    if (degree <= 2) B1 += nsize;             /* D1 & D2 are cheap to prove. */
    if (degree <= 0) B1 += 2*nsize;         /* N-1 and N+1 are really cheap. */
    if (degree > 20 && stage <= 1) B1 -= nsize;   /* Less time on big polys. */
    if (degree > 40) B1 -= nsize/2;               /* Less time on big polys. */
    if (stage == 0) {
      /* A relatively small performance hit, makes slightly smaller proofs. */
      if (nsize < 900 && degree <= 2) B1 *= 1.8;
      /* We need to try a bit harder for the large sizes :( */
      if (nsize > 1400)  B1 *= 2;
      if (nsize > 2000)  B1 *= 2;
      if (!success)
        success = _GMP_pminus1_factor(n, f, 100+B1/8, 100+B1);
    } else if (stage >= 1) {
      /* P-1 */
      if ((!success && do_pm1))
        success = _GMP_pminus1_factor(n, f, B1, 6*B1);
      /* Pollard's Rho */
      if ((!success && do_pbr && nsize < 500))
        success = _GMP_pbrent_factor(n, f, nsize % 53, 1000-nsize);
      /* P+1 */
      if ((!success && do_pp1)) {
        UV ppB = (nsize < 2000) ? B1/4 : B1/16;
        success = _GMP_pplus1_factor(n, f, 0, ppB, ppB);
      }
      if ((!success && do_ecm))
        success = _GMP_ecm_factor_projective(n, f, 400, 2000, 1);
#ifdef USE_LIBECM
      /* TODO: LIBECM in other stages */
      /* Note: this will be substantially slower than our code for small sizes
       *       and the small B1/B2 values we're using. */
      if (!success && degree <= 2 && nsize > 600) {
        ecm_params params;
        ecm_init(params);
        params->method = ECM_ECM;
        mpz_set_ui(params->B2, 10*B1);
        mpz_set_ui(params->sigma, 0);
        success = ecm_factor(f, n, B1/4, params);
        ecm_clear(params);
        if (mpz_cmp(f, n) == 0)  success = 0;
        if (success) { printf("ECM FOUND FACTOR\n"); }
      }
#endif
    }
    /* Try any factors found in previous stage 2+ calls */
    while (!success && sfaci < *nsfacs) {
      if (mpz_divisible_p(n, sfacs[sfaci])) {
        mpz_set(f, sfacs[sfaci]);
        success = 1;
      }
      sfaci++;
    }
    if (stage > 1 && !success) {
      if (stage == 2) {
        /* if (!success) success = _GMP_pbrent_factor(n, f, nsize-1, 8192); */
        if (!success) success = _GMP_pminus1_factor(n, f, 6*B1, 60*B1);
        /* p+1 with different initial point and searching farther */
        if (!success) success = _GMP_pplus1_factor(n, f, 1, B1/2, B1/2);
        if (!success) success = _GMP_ecm_factor_projective(n, f, 250, 2500, 8);
      } else if (stage == 3) {
        if (!success) success = _GMP_pbrent_factor(n, f, nsize+1, 16384);
        if (!success) success = _GMP_pminus1_factor(n, f, 60*B1, 600*B1);
        /* p+1 with a third initial point and searching farther */
        if (!success) success = _GMP_pplus1_factor(n, f, 2, 1*B1, 1*B1);
        if (!success) success = _GMP_ecm_factor_projective(n, f, B1/4, B1*4, 5);
      } else if (stage == 4) {
        if (!success) success = _GMP_pminus1_factor(n, f, 300*B1, 300*20*B1);
        if (!success) success = _GMP_ecm_factor_projective(n, f, B1/2, B1*8, 4);
      } else if (stage >= 5) {
        UV B = B1 * (stage-4) * (stage-4) * (stage-4);
        if (!success) success = _GMP_ecm_factor_projective(n, f, B, 10*B, 8+stage);
      }
    }
    if (success) {
      if (mpz_cmp_ui(f, 1) == 0 || mpz_cmp(f, n) == 0) {
        gmp_printf("factoring %Zd resulted in factor %Zd\n", n, f);
        croak("internal error in ECPP factoring");
      }
      /* Add the factor to the saved factors list */
      if (stage > 1 && *nsfacs < MAX_SFACS) {
        /* gmp_printf(" ***** adding factor %Zd ****\n", f); */
        mpz_init_set(sfacs[*nsfacs], f);
        nsfacs[0]++;
      }
      /* Is the factor f what we want? */
      if ( mpz_cmp(f, fmin) > 0 && is_bpsw_prime(f) )  return 1;
      /* Divide out f */
      mpz_divexact(n, n, f);
    }
  }
  /* n is larger than fmin and not prime */
  mpz_set(f, n);
  return -1;
}

/* See:
 *   (1) Kaltofen, Valente, Yui 1989
 *   (2) Valente 1992 (Thesis)
 *   (3) Konstantinou, Stamatiou, and Zaroliagis (CHES 2002)
 * This code is performing table 1 of reference 3.
 */
static void weber_root_to_hilbert_root(mpz_t r, mpz_t N, long D)
{
  mpz_t A, t;

  if (D < 0) D = -D;
  D = ((D % 4) == 0)  ?  D/4  :  D;
  if ( (D % 8) == 0 )
    return;

  mpz_init(A);  mpz_init(t);

  switch (D % 8) {
    case 1:  if ((D % 3) != 0)  mpz_powm_ui(t, r, 12, N);
             else               mpz_powm_ui(t, r,  4, N);
             mpz_mul_ui(A, t, 64);
             mpz_sub_ui(t, A, 16);
             break;
    case 2:
    case 6:  if ((D % 3) != 0)  mpz_powm_ui(t, r, 12, N);
             else               mpz_powm_ui(t, r,  4, N);
             mpz_mul_ui(A, t, 64);
             mpz_add_ui(t, A, 16);
             break;
    case 5:  if ((D % 3) != 0)  mpz_powm_ui(t, r, 6, N);
             else               mpz_powm_ui(t, r, 2, N);
             mpz_mul_ui(A, t, 64);
             mpz_sub_ui(t, A, 16);
             break;
    case 7:  if (!mpz_invert(t, r, N)) mpz_set_ui(t, 0);
             if ((D % 3) != 0)  mpz_powm_ui(A, t, 24, N);
             else               mpz_powm_ui(A, t,  8, N);
             mpz_sub_ui(t, A, 16);
             break;
    /* Results in degree 3x Hilbert, so typically not used */
    case 3:  if (!mpz_invert(t, r, N)) mpz_set_ui(t, 0);
             if ((D % 3) != 0) {
               mpz_powm_ui(t, t, 24, N);
               mpz_mul_2exp(A, t, 12);
             } else {
               mpz_powm_ui(t, t, 8, N);
               mpz_mul_2exp(A, t, 4);
             }
             mpz_sub_ui(t, A, 16);
             break;
    default: break;
  }
  /* r = t^3 / A */
  mpz_powm_ui(t, t, 3, N);
  if ( ! mpz_divmod(r, t, A, N, r) )
    mpz_set_ui(r, 0);
  mpz_clear(A);  mpz_clear(t);
}


static int find_roots(long D, int poly_index, mpz_t N, mpz_t** roots, int maxroots)
{
  mpz_t* T;
  UV degree;
  long dT, i, nroots;
  int poly_type;
  gmp_randstate_t* p_randstate = get_randstate();

  if (D == -3 || D == -4) {
    *roots = 0;
    return 1;
  }

  degree = poly_class_poly_num(poly_index, NULL, &T, &poly_type);
  if (degree == 0 || (poly_type != 1 && poly_type != 2))
    return 0;

  dT = degree;
  polyz_mod(T, T, &dT, N);

  polyz_roots_modp(roots, &nroots, maxroots, T, dT, N, p_randstate);
  if (nroots == 0) {
    gmp_printf("N = %Zd\n", N);
    croak("Failed to find roots for D = %ld\n", D);
  }
  for (i = 0; i <= dT; i++)
    mpz_clear(T[i]);
  Safefree(T);
#if 0
  if (nroots != dT && get_verbose_level())
    printf("  found %ld roots of the %ld degree poly\n", nroots, dT);
#endif

  /* Convert Weber roots to Hilbert roots */
  if (poly_type == 2)
    for (i = 0; i < nroots; i++)
      weber_root_to_hilbert_root((*roots)[i], N, D);

  return nroots;
}

static void select_curve_params(mpz_t a, mpz_t b, mpz_t g,
                                long D, mpz_t *roots, long i, mpz_t N, mpz_t t)
{
  int N_is_not_1_congruent_3;

  mpz_set_ui(a, 0);
  mpz_set_ui(b, 0);
  if      (D == -3) { mpz_set_si(b, -1); }
  else if (D == -4) { mpz_set_si(a, -1); }
  else {
    mpz_sub_ui(t, roots[i], 1728);
    mpz_mod(t, t, N);
    /* c = (j * inverse(j-1728)) mod n */
    if (mpz_divmod(b, roots[i], t, N, b)) {
      mpz_mul_si(a, b, -3);   /* r = -3c */
      mpz_mul_si(b, b, 2);    /* s =  2c */
    }
  }
  mpz_mod(a, a, N);
  mpz_mod(b, b, N);

  /* g:  1 < g < Ni && (g/Ni) != -1 && (g%3!=1 || cubic non-residue) */
  N_is_not_1_congruent_3 = ! mpz_congruent_ui_p(N, 1, 3);
  for ( mpz_set_ui(g, 2);  mpz_cmp(g, N) < 0;  mpz_add_ui(g, g, 1) ) {
    if (mpz_jacobi(g, N) != -1)
      continue;
    if (N_is_not_1_congruent_3)
      break;
    mpz_sub_ui(t, N, 1);
    mpz_tdiv_q_ui(t, t, 3);
    mpz_powm(t, g, t, N);   /* t = g^((Ni-1)/3) mod Ni */
    if (mpz_cmp_ui(t, 1) == 0)
      continue;
    if (D == -3) {
      mpz_powm_ui(t, t, 3, N);
      if (mpz_cmp_ui(t, 1) != 0)   /* Additional check when D == -3 */
        continue;
    }
    break;
  }
  if (mpz_cmp(g, N) >= 0)    /* No g can be found: N is composite */
    mpz_set_ui(g, 0);
}

static void select_point(mpz_t x, mpz_t y, mpz_t a, mpz_t b, mpz_t N,
                         mpz_t t, mpz_t t2)
{
  mpz_t Q, t3, t4;
  gmp_randstate_t* p_randstate = get_randstate();

  mpz_init(Q); mpz_init(t3); mpz_init(t4);
  mpz_set_ui(y, 0);

  while (mpz_sgn(y) == 0) {
    /* select a Q s.t. (Q,N) != -1 */
    do {
      do {
        /* mpz_urandomm(x, *p_randstate, N); */
        mpz_urandomb(x, *p_randstate, 32);   /* May as well make x small */
        mpz_mod(x, x, N);
      } while (mpz_sgn(x) == 0);
      mpz_mul(t, x, x);
      mpz_add(t, t, a);
      mpz_mul(t, t, x);
      mpz_add(t, t, b);
      mpz_mod(Q, t, N);
    } while (mpz_jacobi(Q, N) == -1);
    /* Select Y */
    sqrtmod(y, Q, N, t, t2, t3, t4);
    /* TODO: if y^2 mod Ni != t, return composite */
    if (mpz_sgn(y) == 0) croak("y == 0 in point selection\n");
  }
  mpz_clear(Q); mpz_clear(t3); mpz_clear(t4);
}

/* Returns 0 (composite), 1 (didn't find a point), 2 (found point) */
int ecpp_check_point(mpz_t x, mpz_t y, mpz_t m, mpz_t q, mpz_t a, mpz_t N,
                     mpz_t t, mpz_t t2)
{
  struct ec_affine_point P, P1, P2;
  int result = 1;

  mpz_init_set(P.x, x);  mpz_init_set(P.y, y);
  mpz_init(P1.x); mpz_init(P1.y);
  mpz_init(P2.x); mpz_init(P2.y);

  mpz_tdiv_q(t, m, q);
  if (!ec_affine_multiply(a, t, N, P, &P2, t2)) {
    mpz_tdiv_q(t, m, q);
    /* P2 should not be (0,1) */
    if (!(mpz_cmp_ui(P2.x, 0) == 0 && mpz_cmp_ui(P2.y, 1) == 0)) {
      mpz_set(t, q);
      if (!ec_affine_multiply(a, t, N, P2, &P1, t2)) {
        /* P1 should be (0,1) */
        if (mpz_cmp_ui(P1.x, 0) == 0 && mpz_cmp_ui(P1.y, 1) == 0) {
          result = 2;
        }
      } else result = 0;
    }
  } else result = 0;

  mpz_clear(P.x);  mpz_clear(P.y);
  mpz_clear(P1.x); mpz_clear(P1.y);
  mpz_clear(P2.x); mpz_clear(P2.y);
  return result;
}

static void update_ab(mpz_t a, mpz_t b, long D, mpz_t g, mpz_t N)
{
  if      (D == -3) { mpz_mul(b, b, g); }
  else if (D == -4) { mpz_mul(a, a, g); }
  else {
    mpz_mul(a, a, g);
    mpz_mul(a, a, g);
    mpz_mul(b, b, g);
    mpz_mul(b, b, g);
    mpz_mul(b, b, g);
  }
  mpz_mod(a, a, N);
  mpz_mod(b, b, N);
}

/* Once we have found a D and q, this will find a curve and point.
 * Returns: 0 (composite), 1 (didn't work), 2 (success)
 * It's debatable what to do with a 1 return.
 */
static int find_curve(mpz_t a, mpz_t b, mpz_t x, mpz_t y,
                      long D, int poly_index, mpz_t m, mpz_t q, mpz_t N, int maxroots)
{
  long nroots, npoints, i, rooti, unity, result;
  mpz_t g, t, t2;
  mpz_t* roots = 0;

  /* TODO: A better way to do this, I believe, would be to have the root
   *       finder set up as an iterator.  That way we'd get the first root,
   *       try to find a curve, and probably we'd be done.  Only if we tried
   *       10+ points on that root would we get another root.  This would
   *       probably be set up as a stack (array) of polynomials plus one
   *       saved root (for when we solve a degree 2 poly).
   */
  /* Step 1: Get the roots of the Hilbert class polynomial. */
  nroots = find_roots(D, poly_index, N, &roots, maxroots);
  if (nroots == 0)
    return 1;

  /* Step 2: Loop selecting curves and trying points.
   *         On average it takes about 3 points, but we'll try 100+. */

  mpz_init(g);  mpz_init(t);  mpz_init(t2);
  npoints = 0;
  result = 1;
  for (rooti = 0; result == 1 && rooti < 50*nroots; rooti++) {
    /* Given this D and root, select curve a,b */
    select_curve_params(a, b, g,  D, roots, rooti % nroots, N, t);
    if (mpz_sgn(g) == 0) { result = 0; break; }

    /* See Cohen 5.3.1, page 231 */
    unity = (D == -3) ? 6 : (D == -4) ? 4 : 2;
    for (i = 0; result == 1 && i < unity; i++) {
      if (i > 0)
        update_ab(a, b, D, g, N);
      npoints++;
      select_point(x, y,  a, b, N, t, t2);
      result = ecpp_check_point(x, y, m, q, a, N, t, t2);
    }
  }
  if (npoints > 10 && get_verbose_level() > 0)
    printf("  # point finding took %ld points\n", npoints);

  if (roots != 0) {
    for (rooti = 0; rooti < nroots; rooti++)
      mpz_clear(roots[rooti]);
    Safefree(roots);
  }
  mpz_clear(g);  mpz_clear(t);  mpz_clear(t2);

  return result;
}

/* Select the 2, 4, or 6 numbers we will try to factor. */
static void choose_m(mpz_t* mlist, long D, mpz_t u, mpz_t v, mpz_t N,
                     mpz_t t, mpz_t Nplus1)
{
  int i, j;
  mpz_add_ui(Nplus1, N, 1);

  mpz_sub(mlist[0], Nplus1, u);     /* N+1-u */
  mpz_add(mlist[1], Nplus1, u);     /* N+1+u */
  for (i = 2; i < 6; i++)
    mpz_set_ui(mlist[i], 0);

  if (D == -3) {
    /* If reading Cohen, be sure to see the errata for page 474. */
    mpz_mul_si(t, v, 3);
    mpz_add(t, t, u);
    mpz_tdiv_q_2exp(t, t, 1);
    mpz_sub(mlist[2], Nplus1, t);   /* N+1-(u+3v)/2 */
    mpz_add(mlist[3], Nplus1, t);   /* N+1+(u+3v)/2 */
    mpz_mul_si(t, v, -3);
    mpz_add(t, t, u);
    mpz_tdiv_q_2exp(t, t, 1);
    mpz_sub(mlist[4], Nplus1, t);   /* N+1-(u-3v)/2 */
    mpz_add(mlist[5], Nplus1, t);   /* N+1+(u-3v)/2 */
  } else if (D == -4) {
    mpz_mul_ui(t, v, 2);
    mpz_sub(mlist[2], Nplus1, t);   /* N+1-2v */
    mpz_add(mlist[3], Nplus1, t);   /* N+1+2v */
  }
  /* m must not be prime */
  for (i = 0; i < 6; i++)
    if (mpz_sgn(mlist[i]) && _GMP_is_prob_prime(mlist[i]))
      mpz_set_ui(mlist[i], 0);
  /* Sort the m values so we test the smallest first */
  for (i = 0; i < 5; i++)
    if (mpz_sgn(mlist[i]))
      for (j = i+1; j < 6; j++)
        if (mpz_sgn(mlist[j]) && mpz_cmp(mlist[i],mlist[j]) > 0)
          mpz_swap( mlist[i], mlist[j] );
}





/* This is the "factor all strategy" FAS version, which ends up being a lot
 * simpler than the FPS code.
 *
 * It should have a little more smarts for not repeating work when repeating
 * steps.  This could be complicated trying to save all state, but I think we
 * could get most of the benefit by keeping a simple list of all factors
 * found after stage 1, and we just try each of them.
 */

#define VERBOSE_PRINT_N(step, ndigits, maxH, factorstage) \
  if (verbose) { \
    printf("%*sN[%d] (%d dig)", i, "", step, ndigits); \
    if (factorstage > 1) printf(" [FS %d]", factorstage); \
    fflush(stdout); \
  }

/* Recursive routine to prove via ECPP */
static int ecpp_down(int i, mpz_t Ni, int facstage, int *pmaxH, int* dilist, mpz_t* sfacs, int* nsfacs, char** prooftextptr)
{
  mpz_t a, b, u, v, m, q, minfactor, sqrtn, mD, t, t2;
  mpz_t mlist[6];
  mpz_t qlist[6];
  UV nm1a;
  IV np1lp, np1lq;
  struct ec_affine_point P;
  int k, dindex, pindex, nidigits, facresult, curveresult, downresult, stage, D;
  int verbose = get_verbose_level();

  nidigits = mpz_sizeinbase(Ni, 10);

  downresult = _GMP_is_prob_prime(Ni);
  if (downresult == 0)  return 0;
  if (downresult == 2) {
    /* No need to put anything in the proof */
    if (verbose) printf("%*sN[%d] (%d dig)  PRIME\n", i, "", i, nidigits);
    return 2;
  }
  if (i == 0 && facstage == 2 && _GMP_miller_rabin_random(Ni, 2, 0) == 0) {
    gmp_printf("\n\n**** BPSW counter-example found?  ****\n**** N = %Zd ****\n\n", Ni);
    return 0;
  }

  VERBOSE_PRINT_N(i, nidigits, *pmaxH, facstage);

  mpz_init(a);  mpz_init(b);
  mpz_init(u);  mpz_init(v);
  mpz_init(m);  mpz_init(q);
  mpz_init(mD); mpz_init(minfactor);  mpz_init(sqrtn);
  mpz_init(t);  mpz_init(t2);
  mpz_init(P.x);mpz_init(P.y);
  for (k = 0; k < 6; k++) {
    mpz_init(mlist[k]);
    mpz_init(qlist[k]);
  }

  /* Any factors q found must be strictly > minfactor.
   * See Atkin and Morain, 1992, section 6.4 */
  mpz_root(minfactor, Ni, 4);
  mpz_add_ui(minfactor, minfactor, 1);
  mpz_mul(minfactor, minfactor, minfactor);
  mpz_sqrt(sqrtn, Ni);

  stage = 0;
  if (nidigits > 700) stage = 1;  /* Too rare to find them */
  if (i == 0 && facstage > 1)  stage = facstage;
  for ( ; stage <= facstage; stage++) {
    int next_stage = (stage > 1) ? stage : 1;
    for (dindex = -1; dindex < 0 || dilist[dindex] != 0; dindex++) {
      int poly_type;  /* just for debugging/verbose */
      int poly_degree;
      int allq = (nidigits < 400);  /* Do all q values together, or not */

      if (dindex == -1) {   /* n-1 and n+1 tests */
        int nm1_success = 0;
        int np1_success = 0;
        const char* ptype = "";
        mpz_sub_ui(m, Ni, 1);
        mpz_sub_ui(t2, sqrtn, 1);
        mpz_tdiv_q_2exp(t2, t2, 1);    /* t2 = minfactor */
        nm1_success = check_for_factor(u, m, t2, t, stage, sfacs, nsfacs, 0);
        mpz_add_ui(m, Ni, 1);
        mpz_add_ui(t2, sqrtn, 1);
        mpz_tdiv_q_2exp(t2, t2, 1);    /* t2 = minfactor */
        np1_success = check_for_factor(v, m, t2, t, stage, sfacs, nsfacs, 0);
        /* If both successful, pick smallest */
        if (nm1_success > 0 && np1_success > 0) {
          if (mpz_cmp(u, v) <= 0) np1_success = 0;
          else                    nm1_success = 0;
        }
        if      (nm1_success > 0) {  ptype = "n-1";  mpz_set(q, u);  D =  1; }
        else if (np1_success > 0) {  ptype = "n+1";  mpz_set(q, v);  D = -1; }
        else                      continue;
        if (verbose) { printf(" %s\n", ptype); fflush(stdout); }
        downresult = ecpp_down(i+1, q, next_stage, pmaxH, dilist, sfacs, nsfacs, prooftextptr);
        if (downresult == 0) goto end_down;   /* composite */
        if (downresult == 1) {   /* nothing found at this stage */
          VERBOSE_PRINT_N(i, nidigits, *pmaxH, facstage);
          continue;
        }
        if (verbose)
          { printf("%*sN[%d] (%d dig) %s", i, "", i, nidigits, ptype); fflush(stdout); }
        curveresult = (nm1_success > 0)
                    ? _GMP_primality_bls_3(Ni, q, &nm1a)
                    : _GMP_primality_bls_15(Ni, q, &np1lp, &np1lq);
        if (verbose) { printf("  %d\n", curveresult); fflush(stdout); }
        if ( ! curveresult ) { /* This ought not happen */
          if (verbose)
            gmp_printf("\n  Could not prove %s with N = %Zd\n", ptype, Ni);
          downresult = 1;
          continue;
        }
        goto end_down;
      }

      pindex = dilist[dindex];
      if (pindex < 0) continue;  /* We marked this for skip */
      /* Get the values for D, degree, and poly type */
      poly_degree = poly_class_poly_num(pindex, &D, NULL, &poly_type);
      if (poly_degree == 0)
        croak("Unknown value in dilist[%d]: %d\n", dindex, pindex);

      if ( (-D % 4) != 3 && (-D % 16) != 4 && (-D % 16) != 8 )
        croak("Invalid discriminant '%d' in list\n", D);
      /* D must also be squarefree in odd divisors, but assume it. */
      /* Make sure we can get a class polynomial for this D. */
      if (poly_degree > 16 && stage == 0) {
        if (verbose) printf(" [1]");
        break;
      }
      /* Make the continue-search vs. backtrack decision */
      if (*pmaxH > 0 && poly_degree > *pmaxH)  break;
      mpz_set_si(mD, D);
      /* (D/N) must be 1, and we have to have a u,v solution */
      if (mpz_jacobi(mD, Ni) != 1)
        continue;
      if ( ! modified_cornacchia(u, v, mD, Ni) )
        continue;

      if (verbose > 1)
        { printf(" %d", D); fflush(stdout); }

      /* We're going to factor all the values for this discriminant then pick
       * the smallest.  This adds a little time, but it means we go down
       * faster.  This makes smaller proofs, and might even save time. */

      choose_m(mlist, D, u, v, Ni, t, t2);
      if (allq) {
        int i, j;
        /* We have 0 to 6 m values.  Try to factor them, put in qlist. */
        for (k = 0; k < 6; k++) {
          mpz_set_ui(qlist[k], 0);
          if (mpz_sgn(mlist[k])) {
            facresult = check_for_factor(qlist[k], mlist[k], minfactor, t, stage, sfacs, nsfacs, poly_degree);
            /* -1 = couldn't find, 0 = no big factors, 1 = found */
            if (facresult <= 0)
              mpz_set_ui(qlist[k], 0);
          }
        }
        /* Sort any q values by size, so we work on the smallest first */
        for (i = 0; i < 5; i++)
          if (mpz_sgn(qlist[i]))
            for (j = i+1; j < 6; j++)
              if (mpz_sgn(qlist[j]) && mpz_cmp(qlist[i],qlist[j]) > 0) {
                mpz_swap( qlist[i], qlist[j] );
                mpz_swap( mlist[i], mlist[j] );
              }
      }
      /* Try to make a proof with the first (smallest) q value.
       * Repeat for others if we have to. */
      for (k = 0; k < 6; k++) {
        int maxH = *pmaxH;
        int minH = (nidigits <= 240) ? 7 : (nidigits+39)/40;

        if (allq) {
          if (mpz_sgn(qlist[k]) == 0) continue;
          mpz_set(m, mlist[k]);
          mpz_set(q, qlist[k]);
        } else {
          if (mpz_sgn(mlist[k]) == 0) continue;
          mpz_set(m, mlist[k]);
          facresult = check_for_factor(q, m, minfactor, t, stage, sfacs, nsfacs, poly_degree);
          if (facresult <= 0) continue;
        }

        if (verbose)
          { printf(" %d (%s %d)\n", D, (poly_type == 1) ? "Hilbert" : "Weber", poly_degree); fflush(stdout); }
        if (maxH == 0) {
          maxH = minH-1 + poly_degree;
          if (facstage > 1)              /* We worked hard to get here, */
            maxH = 2*maxH + 10;          /* try hard to make use of it. */
        } else if (maxH > minH && maxH > (poly_degree+2)) {
          maxH--;
        }
        /* Great, now go down. */
        downresult = ecpp_down(i+1, q, next_stage, &maxH, dilist, sfacs, nsfacs, prooftextptr);
        /* Nothing found, look at more polys in the future */
        if (downresult == 1 && *pmaxH > 0)  *pmaxH = maxH;

        if (downresult == 0) goto end_down;   /* composite */
        if (downresult == 1) {   /* nothing found at this stage */
          VERBOSE_PRINT_N(i, nidigits, *pmaxH, facstage);
          continue;
        }

        /* Awesome, we found the q chain and are in STAGE 2 */
        if (verbose)
          { printf("%*sN[%d] (%d dig) %d (%s %d)", i, "", i, nidigits, D, (poly_type == 1) ? "Hilbert" : "Weber", poly_degree); fflush(stdout); }

        /* Try with only one or two roots, then 8 if that didn't work. */
        /* TODO: This should be done using a root iterator in find_curve() */
        curveresult = find_curve(a, b, P.x, P.y, D, pindex, m, q, Ni, 1);
        if (curveresult == 1) {
          if (verbose) { printf(" [redo roots]"); fflush(stdout); }
          curveresult = find_curve(a, b, P.x, P.y, D, pindex, m, q, Ni, 8);
        }
        if (verbose) { printf("  %d\n", curveresult); fflush(stdout); }
        if (curveresult == 1) {
          /* Something is wrong.  Very likely the class poly coefficients are
             incorrect.  We've wasted lots of time, and need to try again. */
          dilist[dindex] = -2; /* skip this D value from now on */
          if (verbose) gmp_printf("\n  Invalidated D = %d with N = %Zd\n", D, Ni);
          downresult = 1;
          continue;
        }
        /* We found it was composite or proved it */
        goto end_down;
      } /* k loop for D */
    } /* D */
  } /* fac stage */
  /* Nothing at this level */
  if (downresult != 1) croak("ECPP internal error: downresult is %d at end\n", downresult);
  if (verbose) {
    if (*pmaxH > 0) printf(" (max %d)", *pmaxH);
    printf(" ---\n");
    fflush(stdout);
  }
  if (*pmaxH > 0) *pmaxH = *pmaxH + 2;

end_down:

  if (downresult == 2) {
    if (0 && verbose > 1) {
      gmp_printf("\n");
      if (D == 1) {
        gmp_printf("Type BLS3\nN  %Zd\nQ  %Zd\nA  %"UVuf"\n", Ni, q, nm1a);
      } else if (D == -1) {
        gmp_printf("Type BLS15\nN  %Zd\nQ  %Zd\nLP %"IVdf"\nLQ %"IVdf"\n", Ni, q, np1lp, np1lq);
      } else {
        gmp_printf("Type ECPP\nN  %Zd\nA  %Zd\nB  %Zd\nM  %Zd\nQ  %Zd\nX  %Zd\nY  %Zd\n", Ni, a, b, m, q, P.x, P.y);
      }
      gmp_printf("\n");
      fflush(stdout);
    }
    /* Prepend our proof to anything that exists. */
    if (prooftextptr != 0) {
      char *proofstr, *proofptr;
      int curprooflen = (*prooftextptr == 0) ? 0 : strlen(*prooftextptr);

      if (D == 1) {
        int myprooflen = 20 + 2*(4 + mpz_sizeinbase(Ni, 10)) + 1*21;
        New(0, proofstr, myprooflen + curprooflen + 1, char);
        proofptr = proofstr;
        proofptr += gmp_sprintf(proofptr, "Type BLS3\nN  %Zd\nQ  %Zd\nA  %"UVuf"\n", Ni, q, nm1a);
      } else if (D == -1) {
        int myprooflen = 20 + 2*(4 + mpz_sizeinbase(Ni, 10)) + 2*21;
        New(0, proofstr, myprooflen + curprooflen + 1, char);
        proofptr = proofstr;
        proofptr += gmp_sprintf(proofptr, "Type BLS15\nN  %Zd\nQ  %Zd\nLP %"IVdf"\nLQ %"IVdf"\n", Ni, q, np1lp, np1lq);
      } else {
        int myprooflen = 20 + 7*(4 + mpz_sizeinbase(Ni, 10)) + 0;
        New(0, proofstr, myprooflen + curprooflen + 1, char);
        proofptr = proofstr;
        mpz_sub_ui(t, Ni, 1);
        if (mpz_cmp(a, t) == 0)  mpz_set_si(a, -1);
        if (mpz_cmp(b, t) == 0)  mpz_set_si(b, -1);
        proofptr += gmp_sprintf(proofptr, "Type ECPP\nN  %Zd\nA  %Zd\nB  %Zd\nM  %Zd\nQ  %Zd\nX  %Zd\nY  %Zd\n", Ni, a, b, m, q, P.x, P.y);
      }
      if (*prooftextptr) {
        proofptr += gmp_sprintf(proofptr, "\n");
        strcat(proofptr, *prooftextptr);
        Safefree(*prooftextptr);
      }
      *prooftextptr = proofstr;
    }
  }

  /* Ni passed BPSW, so it's highly unlikely to be composite */
  if (downresult == 0) {
    if (mpz_probab_prime_p(Ni, 2) == 0) {
      gmp_printf("\n\n**** BPSW counter-example found?  ****\n**** N = %Zd ****\n\n", Ni);
    } else {
      /* Q was composite, but we don't seem to be. */
      downresult = 1;
    }
  }

  mpz_clear(a);  mpz_clear(b);
  mpz_clear(u);  mpz_clear(v);
  mpz_clear(m);  mpz_clear(q);
  mpz_clear(mD); mpz_clear(minfactor);  mpz_clear(sqrtn);
  mpz_clear(t);  mpz_clear(t2);
  mpz_clear(P.x);mpz_clear(P.y);
  for (k = 0; k < 6; k++) {
    mpz_clear(mlist[k]);
    mpz_clear(qlist[k]);
  }

  return downresult;
}

/* returns 2 if N is proven prime, 1 if probably prime, 0 if composite */
int _GMP_ecpp(mpz_t N, char** prooftextptr)
{
  int* dilist;
  mpz_t* sfacs;
  int i, fstage, result, nsfacs;
  UV nsize = mpz_sizeinbase(N,2);

  /* We must check gcd(N,6), let's check 2*3*5*7*11*13*17*19*23. */
  if (nsize <= 64 || mpz_gcd_ui(NULL, N, 223092870UL) != 1) {
    result = _GMP_is_prob_prime(N);
    if (result != 1) return result;
  }

  init_ecpp_gcds( nsize );

  if (prooftextptr)
    *prooftextptr = 0;

  New(0, sfacs, MAX_SFACS, mpz_t);
  dilist = poly_class_nums();
  nsfacs = 0;
  result = 1;
  for (fstage = 1; fstage < 20; fstage++) {
    int maxH = 0;
    if (fstage == 3 && get_verbose_level())
      gmp_printf("Working hard on: %Zd\n", N);
    result = ecpp_down(0, N, fstage, &maxH, dilist, sfacs, &nsfacs, prooftextptr);
    if (result != 1)
      break;
  }
  Safefree(dilist);
  for (i = 0; i < nsfacs; i++)
    mpz_clear(sfacs[i]);
  Safefree(sfacs);

  return result;
}


#ifdef STANDALONE_ECPP
static void dieusage(char* prog) {
  printf("ECPP-DJ version 1.03.  Dana Jacobsen\n\n");
  printf("Usage: %s [options] <number>\n\n", prog);
  printf("Options:\n");
  printf("   -v     set verbose\n");
  printf("   -V     set extra verbose\n");
  printf("   -c     print certificate\n");
  printf("   -bpsw  use the extra strong BPSW test (probable prime test)\n");
  printf("   -nm1   use n-1 proof only (BLS75 theorem 5)\n");
  printf("   -aks   use AKS for proof\n");
#ifdef USE_APRCL
  printf("   -aprcl use APR-CL for proof\n");
#endif
  printf("   -help  this message\n");
  printf("\n");
  printf("Return codes: 0 prime, 1 composite, 2 prp, 3 error\n");
  exit(3);
}

int main(int argc, char **argv)
{
  mpz_t n;
  int isprime, i, do_printcert;
  int do_nminus1 = 0;
  int do_aks = 0;
  int do_aprcl = 0;
  int do_bpsw = 0;
  int retcode = 3;
  char* cert = 0;

  if (argc < 2) dieusage(argv[0]);
  _GMP_init();
  mpz_init(n);
  set_verbose_level(0);
  do_printcert = 0;

  /* Braindead hacky option parsing */
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      if (strcmp(argv[i], "-v") == 0) {
        set_verbose_level(1);
      } else if (strcmp(argv[i], "-V") == 0) {
        set_verbose_level(2);
      } else if (strcmp(argv[i], "-c") == 0) {
        do_printcert = 1;
      } else if (strcmp(argv[i], "-nm1") == 0) {
        do_nminus1 = 1;
      } else if (strcmp(argv[i], "-aks") == 0) {
        do_aks = 1;
      } else if (strcmp(argv[i], "-aprcl") == 0) {
        do_aprcl = 1;
      } else if (strcmp(argv[i], "-bpsw") == 0) {
        do_bpsw = 1;
      } else if (strcmp(argv[i], "-help") == 0 || strcmp(argv[i], "--help") == 0) {
        dieusage(argv[0]);
      } else {
        printf("Unknown option: %s\n\n", argv[i]);
        dieusage(argv[0]);
      }
      continue;
    }
    mpz_set_str(n, argv[i], 10);
    /* gmp_printf("%Zd\n", n); */

    isprime = _GMP_is_prob_prime(n);
    /* If isprime = 2 here, that means it's so small it fits in the
     * deterministic M-R or BPSW range. */
    if (isprime == 2) {
      Newz(0, cert, 20 + mpz_sizeinbase(n, 10), char);
      gmp_sprintf(cert, "Type Small\nN  %Zd\n", n);
    } else if (isprime == 1) {
      if (do_bpsw) {
        /* Done */
      } else if (do_nminus1) {
        isprime = _GMP_primality_bls_nm1(n, 100, &cert);
      } else if (do_aks) {
        isprime = 2 * _GMP_is_aks_prime(n);
        do_printcert = 0;
      } else if (do_aprcl) {
#ifdef USE_APRCL
        /* int i; for (i = 0; i < 10000; i++) */
        isprime = mpz_aprtcle(n, get_verbose_level());
        do_printcert = 0;
#else
        croak("Compiled without USE_APRCL.  Sorry.");
#endif
      } else {
        /* Quick n-1 test */
        isprime = _GMP_primality_bls_nm1(n, 1, &cert);
        if (isprime == 1)
          isprime = _GMP_ecpp(n, &cert);
      }
    }

    /* printf("(%d digit) ", (int)mpz_sizeinbase(n, 10)); */
    if (isprime == 0) {
      printf("COMPOSITE\n");
      retcode = 1;
    } else if (isprime == 1) {
      /* This would normally only be from BPSW */
      printf("PROBABLY PRIME\n");
      retcode = 2;
    } else if (isprime == 2) {
      if (do_printcert) {
        gmp_printf("[MPU - Primality Certificate]\n");
        gmp_printf("Version 1.0\n");
        gmp_printf("\n");
        gmp_printf("Proof for:\n");
        gmp_printf("N %Zd\n", n);
        gmp_printf("\n");
        printf("%s", cert);
      } else {
        printf("PRIME\n");
      }
      retcode = 0;
    } else {
      /* E.g. APRCL returns -1 for error */
      croak("Unknown return code, probable error.\n");
    }
    if (cert != 0) {
      Safefree(cert);
      cert = 0;
    }
  }
  mpz_clear(n);
  _GMP_destroy();
  return retcode;
}
#endif
