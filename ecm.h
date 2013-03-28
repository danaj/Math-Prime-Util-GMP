#ifndef MPU_ECM_H
#define MPU_ECM_H

#include <gmp.h>
#include "ptypes.h"

extern int  _GMP_ecm_factor_affine(mpz_t n, mpz_t f, UV BMax, UV ncurves);
extern int  _GMP_ecm_factor_projective(mpz_t n, mpz_t f, UV BMax, UV ncurves);

#endif
