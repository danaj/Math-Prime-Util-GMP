#ifndef MPU_FACTOR_H
#define MPU_FACTOR_H

#include <gmp.h>
#include "ptypes.h"

extern void _init_factor(void);

extern int factor(mpz_t n, mpz_t* factors[], int* exponents[]);
extern void clear_factors(int nfactors, mpz_t* pfactors[], int* pexponents[]);

extern uint32_t omega(mpz_t n);
extern uint32_t bigomega(mpz_t n);
extern void sigma(mpz_t res, mpz_t n, unsigned long k);
extern int moebius(mpz_t n);
extern int liouville(mpz_t n);
extern int is_semiprime(mpz_t n);
extern void totient(mpz_t totient, mpz_t n);
extern void jordan_totient(mpz_t tot, mpz_t n, unsigned long k);
extern void carmichael_lambda(mpz_t lambda, mpz_t n);
extern void znorder(mpz_t res, mpz_t a, mpz_t n);
extern void znprimroot(mpz_t root, mpz_t n);
extern void ramanujan_tau(mpz_t res, mpz_t n);

extern UV   _GMP_trial_factor(mpz_t n, UV from_n, UV to_n);
extern int  _GMP_prho_factor(mpz_t n, mpz_t f, UV a, UV rounds);
extern int  _GMP_pbrent_factor(mpz_t n, mpz_t f, UV a, UV rounds);
extern int  _GMP_pminus1_factor(mpz_t n, mpz_t f, UV B1, UV B2);
extern int  _GMP_pplus1_factor(mpz_t n, mpz_t f, UV P0, UV B1, UV B2);
extern int  _GMP_holf_factor(mpz_t n, mpz_t f, UV rounds);
extern int  _GMP_squfof_factor(mpz_t n, mpz_t f, UV rounds);
extern int  _GMP_cheb_factor(mpz_t n, mpz_t f, UV B, UV initx);

extern unsigned long   power_factor(mpz_t n, mpz_t f);

extern mpz_t* divisor_list(int* ndivisors, mpz_t n);

extern int is_smooth(mpz_t n, mpz_t k);
extern int is_rough(mpz_t n, mpz_t k);
extern int is_powerful(mpz_t n, uint32_t k);
extern int is_almost_prime(uint32_t k, mpz_t n);

#endif
