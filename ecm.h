#ifndef MPU_ECM_H
#define MPU_ECM_H

#include <gmp.h>
#include "ptypes.h"

/* Finite output coordinates are reduced modulo n; infinity is independent
 * of (x,y).  In particular, the finite point (0,1) is not the identity. */
struct ec_affine_point  { mpz_t x, y; int infinity; };
/* Multiply on a short Weierstrass curve modulo odd n with k >= 0.  Consume
 * k to zero on success, or partially on failure, allowing R to alias P.
 * Apart from R=P, writable arguments must not alias any other argument;
 * each point's x and y must be distinct.  Return 1 on a failed inversion, with
 * its gcd in d and R unchanged; for valid curve points this is a factor.
 * Otherwise return 0, set d=1, and inspect R->infinity for the identity. */
extern int ec_affine_multiply(
  const mpz_t a, mpz_t k, const mpz_t n,
  const struct ec_affine_point P, struct ec_affine_point *R,
  mpz_t d);

/* You should use projective instead.  n and f must be distinct. */
extern int  _GMP_ecm_factor_affine(const mpz_t n, mpz_t f, UV BMax, UV ncurves);

/* Main ECM factoring routine.  n and f must be distinct.  Projective scratch
 * is file-global: calls must be serialized (not reentrant or thread-safe).
 * B2<B1 requests the default stage-2 bound; an empty stage 2 is skipped. */
extern int  _GMP_ecm_factor_projective(const mpz_t n, mpz_t f, UV B1, UV B2, UV ncurves);

#endif
