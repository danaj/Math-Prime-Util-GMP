#ifndef MPU_ECM_H
#define MPU_ECM_H

#include <gmp.h>
#include "ptypes.h"

struct ec_affine_point  { mpz_t x, y; };
extern int ec_affine_multiply(
  const mpz_t a, mpz_t k, const mpz_t n,
  const struct ec_affine_point P, struct ec_affine_point *R,
  mpz_t d);

extern int  _GMP_ecm_factor_affine(const mpz_t n, mpz_t f, UV BMax, UV ncurves);
extern int  _GMP_ecm_factor_projective(const mpz_t n, mpz_t f, UV B1, UV B2, UV ncurves);

#endif
