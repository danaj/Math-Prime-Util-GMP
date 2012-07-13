#ifndef MPU_GMPMAIN_H
#define MPU_GMPMAIN_H

#include <gmp.h>
#include "ptypes.h"

extern int  _GMP_miller_rabin(mpz_t n, mpz_t a);
extern int  _GMP_is_strong_lucas_pseudoprime(mpz_t n);
extern int  _GMP_trial_div(mpz_t n, UV to_n);

extern int  _GMP_is_prob_prime(mpz_t n);
extern int  _GMP_is_prime(mpz_t n);
extern void _GMP_next_prime(mpz_t n);
extern void _GMP_prev_prime(mpz_t n);

extern int  _GMP_prho_factor(mpz_t n, mpz_t f, UV a, UV rounds);
extern int  _GMP_pbrent_factor(mpz_t n, mpz_t f, UV a, UV rounds);
extern int  _GMP_pminus1_factor(mpz_t n, mpz_t f, UV smoothness_bound, UV B2_multiplier);
extern int  _GMP_holf_factor(mpz_t n, mpz_t f, UV rounds);
extern UV   _GMP_trial_factor(mpz_t n, UV from_n, UV to_n);
extern int  _GMP_squfof_factor(mpz_t n, mpz_t f, UV rounds);

#endif
