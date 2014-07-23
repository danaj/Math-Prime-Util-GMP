#ifndef MPU_FACTOR_H
#define MPU_FACTOR_H

#include <gmp.h>
#include "ptypes.h"

extern int factor(mpz_t n, mpz_t* factors[], int* exponents[]);
extern void clear_factors(int nfactors, mpz_t* pfactors[], int* pexponents[]);

extern int moebius(mpz_t n);
extern int liouville(mpz_t n);
extern void totient(mpz_t totient, mpz_t n);
extern void jordan_totient(mpz_t tot, mpz_t n, unsigned long k);
extern void carmichael_lambda(mpz_t lambda, mpz_t n);
extern void znorder(mpz_t res, mpz_t a, mpz_t n);
extern void znprimroot(mpz_t root, mpz_t n);

#endif
