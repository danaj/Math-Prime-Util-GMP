/*****************************************************************************
 *
 * ECPP - Elliptic Curve Primality Proving
 *
 * Copyright (c) 2013 Dana Jacobsen (dana@acm.org).
 * This is free software; you can redistribute it and/or modify it under
 * the same terms as the Perl 5 programming language system itself.
 *
 * This file is part of the Math::Prime::Util Perl module.  It is not too hard
 * to build this as a standalone program (see the README file).
 *
 * This is pretty good for numbers less than 800 digits.  It is many orders
 * of magnitude faster than the contemporary GMP-ECPP.  It is far, far slower
 * than PRIMO for numbers over 300 digits (it's a bit faster below that).  It
 * is missing optimizations that were in François Morain's 10-20 year old work.
 *
 * A set of fixed discriminants are used, rather than calculating them as
 * needed.  Having a way to calculate values as needed would be a big help.
 * In the interests of space for the MPU package, I've chosen ~500 values which
 * compile into about 35k of data.  This is about 1/5 of the entire code size
 * for the MPU package.  The github repository includes a alternate set of 3182
 * discriminants that compile to 2.3MB.  This is recommended if proving 300+
 * digit numbers is a regular occurance.
 *
 * This version uses the FAS "factor all strategy", meaning it first constructs
 * the entire factor chain, with backtracking if necessary, then will do the
 * elliptic curve proof as it recurses back.
 *
 * If your goal is primality proofs for very large numbers, use Primo.  It's
 * free, it is really fast, it is widely used, it can process batch results,
 * and it makes independently verifiable certificates.  This software also
 * makes certificates in a described format, works well under 800 or so digits,
 * and is open source.
 *
 * There are other proof methods, such as the open source APR-CL code from
 * WraithX.  For numbers over 1000 digits, the APR-CL implementation will
 * likely be faster than this, though does not produce a certificate.
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

#undef USE_NM1
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
static mpz_t _lcm_small;

void init_ecpp_gcds(void) {
  if (_gcdinit == 0) {
    mpz_init(_gcd_small);
    mpz_init(_gcd_large);
    _GMP_pn_primorial(_gcd_small,  3000);
    _GMP_pn_primorial(_gcd_large, 20000);
    mpz_divexact(_gcd_large, _gcd_large, _gcd_small);
    mpz_divexact_ui(_gcd_small, _gcd_small, 2*3*5);
    mpz_init(_lcm_small);
    _GMP_lcm_of_consecutive_integers(300, _lcm_small);
    _gcdinit = 1;
  }
}

void destroy_ecpp_gcds(void) {
  if (!_gcdinit) return;
  mpz_clear(_gcd_small);
  mpz_clear(_gcd_large);
  mpz_clear(_lcm_small);
  _gcdinit = 0;
}



static int check_for_factor2(mpz_t f, mpz_t inputn, mpz_t fmin, mpz_t n, int stage, mpz_t* sfacs, int* nsfacs, int degree)
{
  int success, sfaci;
  UV B1;

  /* Use this so we don't modify their input value */
  mpz_set(n, inputn);

  if (mpz_cmp(n, fmin) <= 0) return 0;

#if 0
  {
    /* Straightforward trial division up to 3000. */
    PRIME_ITERATOR(iter);
    UV tf;
    UV const trial_limit = 3000;
    for (tf = 2; tf < trial_limit; tf = prime_iterator_next(&iter)) {
      if (mpz_cmp_ui(n, tf*tf) < 0) break;
      while (mpz_divisible_ui_p(n, tf))
        mpz_divexact_ui(n, n, tf);
    }
    prime_iterator_destroy(&iter);
  }
#else
  /* Utilize GMP's fast gcd algorithms.  Trial to 224737 with two gcds. */
  mpz_tdiv_q_2exp(n, n, mpz_scan1(n, 0));
  while (mpz_divisible_ui_p(n, 3))  mpz_divexact_ui(n, n, 3);
  while (mpz_divisible_ui_p(n, 5))  mpz_divexact_ui(n, n, 5);
  if (mpz_cmp(n, fmin) <= 0) return 0;
  mpz_gcd(f, n, _gcd_small);
  while (mpz_cmp_ui(f, 1) > 0) {
    mpz_divexact(n, n, f);
    mpz_gcd(f, n, _gcd_small);
  }
  if (mpz_cmp(n, fmin) <= 0) return 0;
  mpz_gcd(f, n, _gcd_large);
  while (mpz_cmp_ui(f, 1) > 0) {
    mpz_divexact(n, n, f);
    mpz_gcd(f, n, _gcd_large);
  }
  /* Quick stage 1 n-1 using a single big powm + gcd. */
  if (stage == 0) {
    if (mpz_cmp(n, fmin) <= 0) return 0;
    mpz_set_ui(f, 2);
    mpz_powm(f, f, _lcm_small, n);
    mpz_sub_ui(f, f, 1);
    mpz_gcd(f, f, n);
    if (mpz_cmp_ui(f, 1) != 0 && mpz_cmp(f, n) != 0) {
      mpz_divexact(n, n, f);
      if (mpz_cmp(f, n) > 0)
        mpz_set(n, f);
    }
  }
#endif

  sfaci = 0;
  success = 1;
  while (success) {
    UV nsize = mpz_sizeinbase(n, 2);

    if (mpz_cmp(n, fmin) <= 0) return 0;
    if (_GMP_is_prob_prime(n)) { mpz_set(f, n); return (mpz_cmp(f, fmin) > 0); }

    success = 0;
    B1 = 300 + 3 * nsize;
    if (degree <= 2) B1 += nsize; /* D1 & D2 are cheap to prove.  Encourage. */
    if (degree <= 0) B1 += nsize; /* N-1 and N+1 are really cheap. */
    if (degree > 20 && stage <= 1) B1 -= nsize;   /* Less time on big polys. */
    if (degree > 40) B1 -= nsize/2;               /* Less time on big polys. */
    if (stage >= 1) {
#ifdef USE_LIBECM
      /* TODO: Tune stage 1 (PM1?) */
      /* TODO: LIBECM in other stages */
      if (!success) {
        ecm_params params;
        ecm_init(params);
        params->method = ECM_ECM;
        mpz_set_ui(params->B2, 10*B1);
        mpz_set_ui(params->sigma, 0);
        success = ecm_factor(f, n, B1/4, params);
        ecm_clear(params);
        if (mpz_cmp(f, n) == 0)  success = 0;
      }
#else
      if (!success) success = _GMP_pminus1_factor(n, f, B1, 6*B1);
      if (!success) success = _GMP_pplus1_factor(n, f, 0, B1/8, B1/8);
      if (!success && nsize < 500) success = _GMP_pbrent_factor(n, f, nsize, 1024-nsize);
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
        if (!success) success = _GMP_pbrent_factor(n, f, nsize-1, 8192);
        if (!success) success = _GMP_pminus1_factor(n, f, 6*B1, 60*B1);
        /* p+1 with different initial point and searching farther */
        if (!success) success = _GMP_pplus1_factor(n, f, 1, B1/2, B1/2);
        if (!success) success = _GMP_ecm_factor_projective(n, f, 250, 3);
      } else if (stage == 3) {
        if (!success) success = _GMP_pbrent_factor(n, f, nsize+1, 16384);
        if (!success) success = _GMP_pminus1_factor(n, f, 60*B1, 600*B1);
        /* p+1 with a third initial point and searching farther */
        if (!success) success = _GMP_pplus1_factor(n, f, 2, 1*B1, 1*B1);
        if (!success) success = _GMP_ecm_factor_projective(n, f, B1/4, 4);
      } else if (stage == 4) {
        if (!success) success = _GMP_pminus1_factor(n, f, 300*B1, 300*20*B1);
        if (!success) success = _GMP_ecm_factor_projective(n, f, B1/2, 4);
      } else if (stage >= 5) {
        UV B = B1 * (stage-4) * (stage-4) * (stage-4);
        if (!success) success = _GMP_ecm_factor_projective(n, f, B, 8+stage);
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
      if ( mpz_cmp(f, fmin) > 0 && _GMP_is_prob_prime(f) )  return 1;
      /* Divide out f */
      mpz_divexact(n, n, f);
    }
  }
  /* n is larger than fmin and not prime */
  mpz_set(f, n);
  return -1;
}

/* return false or true with f (prime) > fmin */
/* Set up for clever FPS, where the reduced n is remembered each time */
static int check_for_factor(mpz_t f, mpz_t inputn, mpz_t fmin, mpz_t n, long stage)
{
  int success;

  /* Use this so we don't modify their input value */
  mpz_set(n, inputn);

  if (mpz_cmp(n, fmin) <= 0) return 0;

  if (stage == 1) {
    /* simple trial division */
    while (mpz_divisible_ui_p(n,  2)) mpz_divexact_ui(n, n,  2);
    while (mpz_divisible_ui_p(n,  3)) mpz_divexact_ui(n, n,  3);
    if (mpz_gcd_ui(NULL, n, 2850092245UL) != 1) {
      while (mpz_divisible_ui_p(n,  5)) mpz_divexact_ui(n, n,  5);
      while (mpz_divisible_ui_p(n,  7)) mpz_divexact_ui(n, n,  7);
      while (mpz_divisible_ui_p(n, 11)) mpz_divexact_ui(n, n, 11);
      while (mpz_divisible_ui_p(n, 13)) mpz_divexact_ui(n, n, 13);
      while (mpz_divisible_ui_p(n, 17)) mpz_divexact_ui(n, n, 17);
      while (mpz_divisible_ui_p(n, 19)) mpz_divexact_ui(n, n, 19);
      while (mpz_divisible_ui_p(n, 41)) mpz_divexact_ui(n, n, 41);
      while (mpz_divisible_ui_p(n, 43)) mpz_divexact_ui(n, n, 43);
    }
    if (mpz_gcd_ui(NULL, n, 2392308223UL) != 1) {
      while (mpz_divisible_ui_p(n, 23)) mpz_divexact_ui(n, n, 23);
      while (mpz_divisible_ui_p(n, 29)) mpz_divexact_ui(n, n, 29);
      while (mpz_divisible_ui_p(n, 31)) mpz_divexact_ui(n, n, 31);
      while (mpz_divisible_ui_p(n, 37)) mpz_divexact_ui(n, n, 37);
      while (mpz_divisible_ui_p(n, 53)) mpz_divexact_ui(n, n, 53);
      while (mpz_divisible_ui_p(n, 59)) mpz_divexact_ui(n, n, 59);
    }
  }

  success = 1;
  while (success) {

    if (mpz_cmp(n, fmin) <= 0) return 0;
    if (_GMP_is_prob_prime(n)) { mpz_set(f, n); return (mpz_cmp(f, fmin) > 0); }

    success = 0;
    if (stage == 1) {
      /* We should probably do something like B1 = 10 * ndigits */
      UV B1 = (mpz_sizeinbase(n, 2) > 1200) ? 6000 : 3000;
      if (!success) success = _GMP_pminus1_factor(n, f, B1, 10*B1);
    } else if (stage == 2) {
      if (!success) success = _GMP_pminus1_factor(n, f, 10000, 200000);
      if (!success) success = _GMP_ecm_factor_projective(n, f, 250, 4);
    } else if (stage == 3) {
      if (!success) success = _GMP_pminus1_factor(n, f, 40000, 800000);
      if (!success) success = _GMP_ecm_factor_projective(n, f, 500, 4);
    } else if (stage == 4) {
      if (!success) success = _GMP_pminus1_factor(n, f, 100000, 5000000);
      if (!success) success = _GMP_ecm_factor_projective(n, f, 1000, 10);
    } else if (stage >= 5) {
      UV B = 4000 * (stage-4) * (stage-4) * (stage-4);
      if (!success) success = _GMP_ecm_factor_projective(n, f, B, 20);
    }
    if (success) {
      if (mpz_cmp_ui(f, 1) == 0 || mpz_cmp(f, n) == 0) {
        gmp_printf("factoring %Zd resulted in factor %Zd\n", n, f);
        croak("internal error in ECPP factoring");
      }
      /* Is the factor f what we want? */
      if ( mpz_cmp(f, fmin) > 0 && _GMP_is_prob_prime(f) )  return 1;
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


static int find_roots(long D, mpz_t N, mpz_t** roots)
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

  degree = poly_class_poly(D, &T, &poly_type);
  if (degree == 0 || (poly_type != 1 && poly_type != 2))
    return 0;

  dT = degree;
  polyz_mod(T, T, &dT, N);

  /* Don't bother getting more than 4 roots */
  polyz_roots_modp(roots, &nroots,  4, T, dT, N, p_randstate);
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
                      long D, mpz_t m, mpz_t q, mpz_t N)
{
  long nroots, npoints, i, rooti, unity, result;
  mpz_t g, t, t2;
  mpz_t* roots = 0;
  int verbose = get_verbose_level();

  /* Step 1: Get the roots of the Hilbert class polynomial. */
  nroots = find_roots(D, N, &roots);
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
  if (verbose && npoints > 10)
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
                     mpz_t t, mpz_t Nminus1)
{
  int i;
  mpz_add_ui(Nminus1, N, 1);

  mpz_add(mlist[0], Nminus1, u);
  mpz_sub(mlist[1], Nminus1, u);
  for (i = 2; i < 6; i++)
    mpz_set_ui(mlist[i], 0);

  if (D == -3) {
    /* If reading Cohen, be sure to see the errata for page 474. */
    mpz_mul_si(t, v, 3);
    mpz_add(t, t, u);
    mpz_tdiv_q_2exp(t, t, 1);
    mpz_add(mlist[2], Nminus1, t);
    mpz_sub(mlist[3], Nminus1, t);
    mpz_mul_si(t, v, -3);
    mpz_add(t, t, u);
    mpz_tdiv_q_2exp(t, t, 1);
    mpz_add(mlist[4], Nminus1, t);
    mpz_sub(mlist[5], Nminus1, t);
  } else if (D == -4) {
    mpz_mul_ui(t, v, 2);
    mpz_add(mlist[2], Nminus1, t);
    mpz_sub(mlist[3], Nminus1, t);
  }
  /* m must not be prime */
  for (i = 0; i < 6; i++)
    if (mpz_sgn(mlist[i]) && _GMP_is_prob_prime(mlist[i]))
      mpz_set_ui(mlist[i], 0);
}





/* This is the "factor all strategy" FAS version, which ends up being a lot
 * simpler than the FPS code.
 *
 * It should have a little more smarts for not repeating work when repeating
 * steps.  This could be complicated trying to save all state, but I think we
 * could get most of the benefit by keeping a simple list of all factors
 * found after stage 1, and we just try each of them.
 */

/* Recursive routine to prove via ECPP */
static int ecpp_down(int i, mpz_t Ni, int facstage, IV* dlist, mpz_t* sfacs, int* nsfacs, char** prooftextptr)
{
  mpz_t a, b, u, v, m, q, minfactor, sqrtn, mD, t, t2;
  mpz_t mlist[6];
  UV nm1a;
  IV np1lp, np1lq;
  struct ec_affine_point P;
  int k, dnum, nidigits, facresult, curveresult, downresult, stage, D;
  int verbose = get_verbose_level();

  nidigits = mpz_sizeinbase(Ni, 10);

  k = _GMP_is_prob_prime(Ni);
  if (k == 0)  return 0;
  if (k == 2) {
    /* No need to put anything in the proof */
    if (verbose) printf("%*sN[%d] (%d dig)  PRIME\n", i, "", i, nidigits);
    return 2;
  }
  downresult = 0;

  if (verbose) {
    printf("%*sN[%d] (%d dig)", i, "", i, nidigits);
    if (facstage > 1) printf(" FS %d", facstage);
    fflush(stdout);
  }

  mpz_init(a);  mpz_init(b);
  mpz_init(u);  mpz_init(v);
  mpz_init(m);  mpz_init(q);
  mpz_init(mD); mpz_init(minfactor);  mpz_init(sqrtn);
  mpz_init(t);  mpz_init(t2);
  mpz_init(P.x);mpz_init(P.y);
  for (k = 0; k < 6; k++)
    mpz_init(mlist[k]);

  /* Any factors q found must be strictly > minfactor.
   * See Atkin and Morain, 1992, section 6.4 */
  mpz_root(minfactor, Ni, 4);
  mpz_add_ui(minfactor, minfactor, 1);
  mpz_mul(minfactor, minfactor, minfactor);
  mpz_sqrt(sqrtn, Ni);

  for (stage = (i == 0) ? facstage : 1; stage <= facstage; stage++) {
    int next_stage = (stage > 1) ? stage : 1;
    for (dnum = 0; dlist[dnum] != 0; dnum++) {
      int poly_type;  /* just for debugging/verbose */
      int poly_degree;
      D = -dlist[dnum];
      if (D > 1) continue;  /* Marked for skip */

      if (D == 1) {   /* n-1 test */
        mpz_sub_ui(m, Ni, 1);          /* m = N-1 */
        mpz_sub_ui(t2, sqrtn, 1);
        mpz_tdiv_q_2exp(t2, t2, 1);    /* t2 = minfactor */

        facresult = check_for_factor2(q, m, t2, t, stage, sfacs, nsfacs, 0);
        if (facresult <= 0)
          continue;
        if (verbose)
          { printf(" n-1\n"); fflush(stdout); }
        downresult = ecpp_down(i+1, q, next_stage, dlist, sfacs, nsfacs, prooftextptr);
        if (downresult == 0)     /* composite */
          goto end_down;
        if (downresult == 1) {   /* nothing found at this stage */
          if (verbose) {
            printf("%*sN[%d] (%d dig)", i, "", i, nidigits);
            if (facstage > 1) printf(" FS %d", facstage);
            fflush(stdout);
          }
          continue;
        }
        if (verbose)
          { printf("%*sN[%d] (%d dig) n-1", i, "", i, nidigits); fflush(stdout); }
        curveresult = _GMP_primality_bls_3(Ni, q, &nm1a);
        if (verbose) { printf("  %d\n", curveresult); fflush(stdout); }
        if ( ! curveresult ) {
          /* This ought not happen */
          dlist[dnum] = -2; /* skip this D value from now on */
          if (verbose) gmp_printf("\n  Could not prove n-1 with N = %Zd\n", D, Ni);
          downresult = 0;
          continue;
        }
        goto end_down;
      }
      if (D == -1) {  /* n+1 test */
        mpz_add_ui(m, Ni, 1);          /* m = N+1 */
        mpz_add_ui(t2, sqrtn, 1);
        mpz_tdiv_q_2exp(t2, t2, 1);    /* t2 = minfactor */

        facresult = check_for_factor2(q, m, t2, t, stage, sfacs, nsfacs, 0);
        if (facresult <= 0)
          continue;
        if (verbose)
          { printf(" n+1\n"); fflush(stdout); }
        downresult = ecpp_down(i+1, q, next_stage, dlist, sfacs, nsfacs, prooftextptr);
        if (downresult == 0)     /* composite */
          goto end_down;
        if (downresult == 1) {   /* nothing found at this stage */
          if (verbose) {
            printf("%*sN[%d] (%d dig)", i, "", i, nidigits);
            if (facstage > 1) printf(" FS %d", facstage);
            fflush(stdout);
          }
          continue;
        }
        if (verbose)
          { printf("%*sN[%d] (%d dig) n-1", i, "", i, nidigits); fflush(stdout); }
        curveresult = _GMP_primality_bls_15(Ni, q, &np1lp, &np1lq);
        if (verbose) { printf("  %d\n", curveresult); fflush(stdout); }
        if ( ! curveresult ) {
          /* This ought not happen */
          dlist[dnum] = -2; /* skip this D value from now on */
          if (verbose) gmp_printf("\n  Could not prove n-1 with N = %Zd\n", D, Ni);
          downresult = 0;
          continue;
        }
        goto end_down;
      }

      if ( (-D % 4) != 3 && (-D % 16) != 4 && (-D % 16) != 8 )
        croak("Invalid discriminant '%d' in list\n", D);
      /* D must also be squarefree in odd divisors, but assume it. */
      /* Make sure we can get a class polynomial for this D. */
      poly_degree = poly_class_poly(D, NULL, &poly_type);
      if (poly_degree == 0)  continue;
      /* We'll save time in the long run by not looking at big polys once
       * we've found a good path from the start.  TODO: Needs more tuning. */
      if (stage == 0 && i >= 0 && poly_degree > 2) break;
      if (facstage == 1) {
        if (i >  2 && nidigits < 1100 && poly_degree > 24)  break;
        if (i >  3 && nidigits <  950 && poly_degree > 15)  break;
        if (i >  4 && nidigits <  800 && poly_degree > 11)  break;
        if (i >  8 && nidigits <  700 && poly_degree >  9)  break;
        if (i > 16 && nidigits <  600 && poly_degree >  8)  break;
      }
      mpz_set_si(mD, D);
      /* (D/N) must be 1, and we have to have a u,v solution */
      if (mpz_jacobi(mD, Ni) != 1)
        continue;
      if ( ! modified_cornacchia(u, v, mD, Ni) )
        continue;

      if (verbose > 1)
        { printf(" %d", D); fflush(stdout); }

      choose_m(mlist, D, u, v, Ni, t, t2);
      for (k = 0; k < 6; k++) {
        facresult = check_for_factor2(q, mlist[k], minfactor, t, stage, sfacs, nsfacs, poly_degree);
        /* -1 = couldn't find, 0 = no big factors, 1 = found */
        if (facresult <= 0)
          continue;
        mpz_set(m, mlist[k]);
        if (verbose)
          { printf(" %d (%s %d)\n", D, (poly_type == 1) ? "Hilbert" : "Weber", poly_degree); fflush(stdout); }
        /* Great, now go down. */
        downresult = ecpp_down(i+1, q, next_stage, dlist, sfacs, nsfacs, prooftextptr);
        if (downresult == 0)     /* composite */
          goto end_down;
        if (downresult == 1) {   /* nothing found at this stage */
          if (verbose) {
            printf("%*sN[%d] (%d dig)", i, "", i, nidigits);
            if (facstage > 1) printf(" FS %d", facstage);
            fflush(stdout);
          }
          continue;
        }

        /* Awesome, we found the q chain and are in STAGE 2 */
        if (verbose)
          { printf("%*sN[%d] (%d dig) %d (%s %d)", i, "", i, nidigits, D, (poly_type == 1) ? "Hilbert" : "Weber", poly_degree); fflush(stdout); }

        curveresult = find_curve(a, b, P.x, P.y, D, m, q, Ni);
        if (verbose) { printf("  %d\n", curveresult); fflush(stdout); }
        if (curveresult == 1) {
          /* Oh no!  We can't find a point on the curve.  Something is right
           * messed up, and we've wasted a lot of time.  Sigh. */
          dlist[dnum] = -2; /* skip this D value from now on */
          if (verbose) gmp_printf("\n  Invalidated D = %d with N = %Zd\n", D, Ni);
          downresult = 0;
          continue;
        }
        /* We found it was composite or proved it */
        goto end_down;
      } /* k loop for D */
    } /* D */
  } /* fac stage */
  /* Nothing at this level */
  downresult = 1;
  if (verbose) { printf(" ---\n"); fflush(stdout); }

end_down:

  if (downresult == 2) {
    if (0 && verbose > 1) {
      if (D == 1) {
        gmp_printf("\n");
        gmp_printf("Type BLS3\n");
        gmp_printf("N  %Zd\n", Ni);
        gmp_printf("Q  %Zd\n", q);
        gmp_printf("A  %lu\n", (unsigned long) nm1a);
        gmp_printf("\n");
        fflush(stdout);
      } else if (D == -1) {
        gmp_printf("\n");
        gmp_printf("Type BLS15\n");
        gmp_printf("N  %Zd\n", Ni);
        gmp_printf("Q  %Zd\n", q);
        gmp_printf("LP %ld\n", (signed long) np1lp);
        gmp_printf("LQ %ld\n", (signed long) np1lq);
        gmp_printf("\n");
        fflush(stdout);
      } else {
        gmp_printf("\n");
        gmp_printf("Type ECPP\n");
        gmp_printf("N  %Zd\n", Ni);
        gmp_printf("A  %Zd\n", a);
        gmp_printf("B  %Zd\n", b);
        gmp_printf("M  %Zd\n", m);
        gmp_printf("Q  %Zd\n", q);
        gmp_printf("X  %Zd\n", P.x);
        gmp_printf("Y  %Zd\n", P.y);
        gmp_printf("\n");
        fflush(stdout);
      }
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

  mpz_clear(a);  mpz_clear(b);
  mpz_clear(u);  mpz_clear(v);
  mpz_clear(m);  mpz_clear(q);
  mpz_clear(mD); mpz_clear(minfactor);  mpz_clear(sqrtn);
  mpz_clear(t);  mpz_clear(t2);
  mpz_clear(P.x);mpz_clear(P.y);
  for (k = 0; k < 6; k++)
    mpz_clear(mlist[k]);

  /* TODO: proof text */
  return downresult;
}

/* returns 2 if N is proven prime, 1 if probably prime, 0 if composite */
int _GMP_ecpp(mpz_t N, char** prooftextptr)
{
  IV* dlist;
  mpz_t* sfacs;
  int i, fstage, result, nsfacs;
  int verbose = get_verbose_level();

  /* We must check gcd(N,6), let's check 2*3*5*7*11*13*17*19*23. */
  if (mpz_gcd_ui(NULL, N, 223092870UL) != 1)
    return _GMP_is_prob_prime(N);

  init_ecpp_gcds();

  if (prooftextptr)
    *prooftextptr = 0;

  New(0, sfacs, MAX_SFACS, mpz_t);
  dlist = poly_class_degrees(1);
  nsfacs = 0;
  result = 1;
  for (fstage = 1; fstage < 20; fstage++) {
    if (verbose && fstage == 3) gmp_printf("Working hard on: %Zd\n", N);
    result = ecpp_down(0, N, fstage, dlist, sfacs, &nsfacs, prooftextptr);
    if (result != 1)
      break;
  }
  Safefree(dlist);
  for (i = 0; i < nsfacs; i++)
    mpz_clear(sfacs[i]);
  Safefree(sfacs);

  return result;
}




/* This is the "factor and prove strategy" version, with no backtracking.
 * It tries to be efficient with factoring, but with a limited set of
 * discriminants and no backtracking it tends to get stuck factoring once
 * into the 250+ digit range.  More discriminants helps, but it just pushes
 * the problem out some.  The FAS version tends to be much faster.
 */
typedef struct {
  long D;
  mpz_t m;
  mpz_t q;
} dmqlist_t;

/* returns 2 if N is proven prime, 1 if probably prime, 0 if composite */
int _GMP_ecpp_fps(mpz_t N, char** prooftextptr)
{
  mpz_t Ni, a, b, u, v, m, q, mD, t, t2, minfactor;
  mpz_t mlist[6];
  IV* dlist;
  struct ec_affine_point P;
  UV i, dnum;
  long D, k, numQs, numDs, stage;
  int result = 0;
  int facresult, curveresult, niresult;
  dmqlist_t *dmqlist;
  char *proofstr, *proofptr;
  size_t proofstr_size;
  int verbose = get_verbose_level();

  /* We must check gcd(N,6), let's check 2*3*5*7*11*13*17*19*23. */
  if (mpz_gcd_ui(NULL, N, 223092870UL) != 1)
    return _GMP_is_prob_prime(N);

  init_ecpp_gcds();

  if (prooftextptr) {
    proofstr_size = 4096;
    New(0, proofstr, proofstr_size, char);
    *proofstr = 0;
  } else {
    proofstr_size = 0;
    proofstr = 0;
    proofptr = 0;
  }
  proofptr = proofstr;
  mpz_init_set(Ni, N);
  mpz_init(a); mpz_init(b); mpz_init(u); mpz_init(v);
  mpz_init(m); mpz_init(q);
  mpz_init(mD);
  mpz_init(t); mpz_init(t2);
  mpz_init(minfactor);
  mpz_init(P.x); mpz_init(P.y);
  dlist = poly_class_degrees(0);
  for (k = 0; k < 6; k++)
    mpz_init(mlist[k]);
  /* TODO: make this 100 and realloc if needed */
  New(0, dmqlist, 1000, dmqlist_t);

  numQs = 0;
  for (i = 0; i < 10000000; i++) {
    if (verbose) {
      int nidigits = mpz_sizeinbase(Ni, 10);
      printf("  N[%d] (%d digits)\n", (int)i, nidigits);
    }
    /* Stop now if Ni is small enough to be known prime. */
    k = _GMP_is_prob_prime(Ni);
    if (k != 1) {
      if (k == 2) result = 2;
      goto end_ecpp;
    }
#ifdef USE_NM1
    /* Try a simple n-1 test.  Note that Atkin and Morain discuss using
     * n-1 and n+1 tests at each level, so this is fairly common.  It probably
     * made even more sense in 1992 when the usual cutoff was a sieve.
     *
     * For not-big numbers, Effort 1 is nearly free, 2 is very fast, 3 is ok.
     * The likelihood of this succeeding goes down rapidly with size, and the
     * cost goes up.
     * For effort 2, it is about 50% at 128 bits, 25% at 160, 10% at 190. */
    {
      UV log2ni = mpz_sizeinbase(Ni, 2);
      int effort = (log2ni <= 128) ? 2
                 : (log2ni <= 256) ? 1
                 :                   0;
      if (effort > 0) {
        char* nm1_proof = 0;
        /* -1 = composite, 0 = dunno, 1 = proved prime */
        k = _GMP_primality_bls_nm1(Ni, effort, &nm1_proof);
        if (k > 0) {
          result = 2;
          if (prooftextptr && nm1_proof != 0) {
            size_t new_size = strlen(proofstr) + strlen(nm1_proof) + 1;
            if ( new_size > proofstr_size ) {
              proofstr_size = new_size;
              Renew(proofstr, proofstr_size, char);
            }
            (void) strcat(proofstr, nm1_proof);
          }
        }
        if (nm1_proof != 0)
          Safefree(nm1_proof);
        if (k != 0)
          goto end_ecpp;
      }
    }
#endif

    /* Any factors q found must be strictly > minfactor.
     * See Atkin and Morain, 1992, section 6.4 */
    mpz_root(minfactor, Ni, 4);
    mpz_add_ui(minfactor, minfactor, 1);
    mpz_mul(minfactor, minfactor, minfactor);

    niresult = 0;
    numDs = 0;
    if (numQs != 0) croak("assert numQs == 0\n");
    for (stage = 1; stage <= 20; stage++) {
      if (stage == 1) {
        /* Walk the degree-sorted discriminants */
        for (dnum = 0; dlist[dnum] > 0; dnum++) {
          int poly_type;  /* just for debugging/verbose */
          int poly_degree;
          D = -dlist[dnum];
          if ( (-D % 4) != 3 && (-D % 16) != 4 && (-D % 16) != 8 )
            croak("Invalid discriminant '%ld' in list\n", D);
          /* D must also be squarefree in odd divisors, but assume it. */
          /* Make sure we can get a class polynomial for this D. */
          poly_degree = poly_class_poly(D, NULL, &poly_type);
          if (poly_degree == 0)  continue;
          /* For debugging, one can put things here like:
           * if poly_degree < 4) continue;
           */
          mpz_set_si(mD, D);
          /* (D/N) must be 1, and we have to have a u,v solution */
          if (mpz_jacobi(mD, Ni) != 1)
            continue;
          if ( ! modified_cornacchia(u, v, mD, Ni) )
            continue;
          numDs++;

          if (verbose > 1)
            { printf("  Candidate D %ld (%s poly)\n",
                     D, (poly_type == 1) ? "Hilbert" : "Weber");
              fflush(stdout); }

          choose_m(mlist, D, u, v, Ni, t, t2);
          for (k = 0; k < 6; k++) {
            facresult = check_for_factor(q, mlist[k], minfactor, t, stage);
            if (facresult == 0)
              continue;
            if (facresult > 0) {
              mpz_set(m, mlist[k]);
              if (verbose)
                printf("  Found factor with D = %ld (class %d %s poly)\n",
                       D, poly_degree, (poly_type == 1) ? "Hilbert" : "Weber");
              curveresult = find_curve(a, b, P.x, P.y, D, m, q, Ni);
              if (curveresult == 0)
                goto end_ecpp;
              /* If curveresult=1, then we're stumped.  Skip it? */
              if (curveresult == 2) { niresult = 1; break; }
            }
            /* Add to list for future stages */
            dmqlist[numQs].D = D;
            mpz_init_set(dmqlist[numQs].m, mlist[k]);
            mpz_init_set(dmqlist[numQs].q, q);
            numQs++;
          }
          if (niresult) break;
        }
      } else {
        /* We may want to consider sorting by size instead.  Certainly after
         * stage 2 we're far more concerned with factoring than root finding. */
        if (verbose)
          { printf("  At stage %ld with %ld D values and %ld q values\n", stage, numDs, numQs); fflush(stdout); }
        for (k = 0; k < numQs; k++) {
          mpz_set(m, dmqlist[k].m);
          facresult = check_for_factor(q, dmqlist[k].q, minfactor, t, stage);
          if (facresult == 0) {
            if (verbose > 1)
              printf("  q has no large factor\n");
            mpz_set(dmqlist[k].q, q);
          } else if (facresult == -1) {
            if (verbose > 1) {
              unsigned long qdigits = mpz_sizeinbase(q, 10);
              if (mpz_cmp(dmqlist[k].q,q) != 0)
                printf("  q reduced (%lu digits)\n", qdigits);
              else
                printf("  q no luck (%lu digits)\n", qdigits);
            }
            mpz_set(dmqlist[k].q, q);
          } else {
            D = dmqlist[k].D;
            curveresult = find_curve(a, b, P.x, P.y, D, m, q, Ni);
            if (curveresult == 0)
              goto end_ecpp;
            /* If k=1, then we're stumped.  Skip it? */
            if (curveresult == 2) { niresult = 1; break; }
          }
        }
      }

      if (niresult)  break;

      if (stage == 1 && numDs == 0) {
        /* If no Ds were found, make sure we don't have a composite. */
        if (_GMP_miller_rabin_random(Ni, 10) == 0) {
          if (verbose) gmp_printf("composite %Zd found in ecpp!\n", Ni);
          goto end_ecpp;
        }
      }

    } /* end stage loop */
    for ( ; numQs > 0; numQs--) {
      mpz_clear(dmqlist[numQs-1].m);
      mpz_clear(dmqlist[numQs-1].q);
    }
    if (!niresult) {
      if (verbose) printf("Giving up, sorry.\n");
      result = 1;
      goto end_ecpp;
    }
    if (verbose > 1) {
      gmp_printf("N[%lu] = %Zd\n", (unsigned long) i, Ni);
      gmp_printf("a = %Zd\n", a);
      gmp_printf("b = %Zd\n", b);
      gmp_printf("m = %Zd\n", m);
      gmp_printf("q = %Zd\n", q);
      gmp_printf("P = (%Zd, %Zd)\n", P.x, P.y);
      fflush(stdout);
    }
    if (prooftextptr != 0) {
      int prooflen = mpz_sizeinbase(Ni, 10) * 7 + 20;
      size_t offset = proofptr - proofstr;
      if ( (offset + 1 + prooflen) > proofstr_size ) {
        proofstr_size += 4*prooflen;
        Renew(proofstr, proofstr_size, char);
        proofptr = proofstr + offset;
      }
      proofptr += gmp_sprintf(proofptr, "%Zd : ECPP : ", Ni);
      proofptr += gmp_sprintf(proofptr,
                              "%Zd %Zd %Zd %Zd (%Zd:%Zd)\n",
                              a, b, m, q, P.x, P.y);
    }
    mpz_set(Ni, q);
  }

end_ecpp:
  mpz_clear(Ni);
  mpz_clear(a); mpz_clear(b); mpz_clear(u); mpz_clear(v);
  mpz_clear(m); mpz_clear(q);
  mpz_clear(mD);
  mpz_clear(t); mpz_clear(t2);
  mpz_clear(minfactor);
  mpz_clear(P.x); mpz_clear(P.y);
  Safefree(dlist);
  for (k = 0; k < 6; k++)
    mpz_clear(mlist[k]);
  for ( ; numQs > 0; numQs--) {
    mpz_clear(dmqlist[numQs-1].m);
    mpz_clear(dmqlist[numQs-1].q);
  }
  Safefree(dmqlist);
  /* Append our proof */
  if (prooftextptr != 0) {
    if (result != 2) {
      Safefree(proofstr);
    } else if (*prooftextptr == 0) {
      *prooftextptr = proofstr;
    } else {
      Renew(*prooftextptr, strlen(*prooftextptr) + strlen(proofstr) + 1, char);
      (void) strcat(*prooftextptr, proofstr);
      Safefree(proofstr);
    }
  }

  return result;
}

#ifdef STANDALONE_ECPP
static void dieusage(char* prog) {
  printf("ECPP-DJ version 1.0.  Dana Jacobsen\n\n");
  printf("Usage: %s [options] <number>\n\n", prog);
  printf("Options:\n");
  printf("   -v     set verbose\n");
  printf("   -V     set extra verbose\n");
  printf("   -c     print certificate\n");
  printf("   -nm1   use n-1 proof only\n");
  printf("   -aks   use AKS for proof\n");
  printf("   -help  this message\n");
  exit(1);
}

int main(int argc, char **argv)
{
  mpz_t n;
  int isprime, i, do_printcert;
  int do_nminus1 = 0;
  int do_aks = 0;
  int do_aprcl = 0;
  int do_bpsw = 0;
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
    } else if (isprime == 1) {
      /* We really shouldn't ever see this. */
      printf("PROBABLY PRIME\n");
    } else {
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
    }
    if (cert != 0) {
      Safefree(cert);
      cert = 0;
    }
  }
  mpz_clear(n);
  _GMP_destroy();
  return 0;
}
#endif
