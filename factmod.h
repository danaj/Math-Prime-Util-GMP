#ifndef MPU_GMP_FACTMOD_H
#define MPU_GMP_FACTMOD_H

#include <gmp.h>

extern void factorialmod(mpz_t r, const mpz_t n, const mpz_t m);
extern void binomialmod(mpz_t r, const mpz_t n, const mpz_t k, const mpz_t m);

#endif
