#ifndef MPU_ECM_H
#define MPU_ECM_H

#include <gmp.h>
#include "ptypes.h"

struct ec_affine_point  { mpz_t x, y; };
extern int ec_affine_multiply(
  mpz_t a, mpz_t k, mpz_t n,
  struct ec_affine_point P, struct ec_affine_point *R,
  mpz_t d);

extern int  _GMP_ecm_factor_affine(mpz_t n, mpz_t f, UV BMax, UV ncurves);
extern int  _GMP_ecm_factor_projective(mpz_t n, mpz_t f, UV BMax, UV ncurves);

#endif
