#ifndef MPU_ECPP_H
#define MPU_ECPP_H

#include <gmp.h>
#include "ptypes.h"

extern int _GMP_ecpp(mpz_t N, char** prooftextptr);

extern int ecpp_check_point(mpz_t x, mpz_t y, mpz_t m, mpz_t q, mpz_t a,
                            mpz_t N, mpz_t t, mpz_t t2);

#endif
