#ifndef MPU_FACTOR_H
#define MPU_FACTOR_H

#include <gmp.h>
#include "ptypes.h"

extern void _init_factor(void);

extern int factor(const mpz_t n, mpz_t* factors[], int* exponents[]);
extern void clear_factors(int nfactors, mpz_t* pfactors[], int* pexponents[]);

extern uint32_t omega(const mpz_t n);
extern uint32_t bigomega(const mpz_t n);
extern void sigma(mpz_t res, const mpz_t n, unsigned long k);
extern int moebius(const mpz_t n);
extern int is_square_free_ui(unsigned long n);
extern int is_square_free(const mpz_t n);
extern int liouville(const mpz_t n);
extern int is_semiprime(const mpz_t n);
extern void totient(mpz_t totient, const mpz_t n);
extern void jordan_totient(mpz_t tot, const mpz_t n, unsigned long k);
extern void carmichael_lambda(mpz_t lambda, const mpz_t n);
extern void znorder(mpz_t res, const mpz_t a, const mpz_t n);
extern void znprimroot(mpz_t root, const mpz_t n);
/* extern void ramanujan_tau(mpz_t res, const mpz_t n); */

extern UV   _GMP_trial_factor(const mpz_t n, UV from_n, UV to_n);
extern int  _GMP_prho_factor(const mpz_t n, mpz_t f, UV a, UV rounds);
extern int  _GMP_pbrent_factor(const mpz_t n, mpz_t f, UV a, UV rounds);
extern int  _GMP_pminus1_factor(const mpz_t n, mpz_t f, UV B1, UV B2);
extern int  _GMP_pplus1_factor(const mpz_t n, mpz_t f, UV P0, UV B1, UV B2);
extern int  _GMP_holf_factor(const mpz_t n, mpz_t f, UV rounds);
extern int  _GMP_squfof_factor(const mpz_t n, mpz_t f, UV rounds);
extern int  _GMP_cheb_factor(const mpz_t n, mpz_t f, UV B, UV initx);

extern unsigned long   power_factor(const mpz_t n, mpz_t f);

extern mpz_t* divisor_list(int* ndivisors, const mpz_t n, const mpz_t maxd);

extern int is_smooth(const mpz_t n, const mpz_t k);
extern int is_rough(const mpz_t n, const mpz_t k);
extern int is_powerful(const mpz_t n, uint32_t k);
extern int is_almost_prime(uint32_t k, const mpz_t n);

/* Trial factor iterator.  next returns factors <= B until it returns 0. */
extern void* trial_factor_iterator_create(const mpz_t n, UV B);
extern int   trial_factor_iterator_next(unsigned long *f, uint32_t *e, void* ctx);
extern void  trial_factor_iterator_n(mpz_t n, void* ctx);
extern void  trial_factor_iterator_destroy(void* ctx);

/* mostly for testing, these do trial factoring of n up to B. */
extern int tfall1(mpz_t n, UV B, mpz_t* pfactors[]);
extern int tfall2(mpz_t n, UV B, mpz_t* pfactors[]);

#endif
