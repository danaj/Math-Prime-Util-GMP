#ifndef MPU_ECPP_H
#define MPU_ECPP_H

#include <gmp.h>
#include "ptypes.h"

extern void init_ecpp_gcds(UV nsize);
extern void destroy_ecpp_gcds(void);

extern int _GMP_ecpp(const mpz_t N, char** prooftextptr);
extern int _GMP_ecpp_fps(const mpz_t N, char** prooftextptr);

extern int ecpp_check_point(const mpz_t x, const mpz_t y,
                            const mpz_t m, const mpz_t q, mpz_t a,
                            const mpz_t N, mpz_t t, mpz_t t2);

#endif
