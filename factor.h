/* vim: set et ts=2 sw=2 sts=2: */
#ifndef MPU_FACTOR_H
#define MPU_FACTOR_H

#include <gmp.h>
#include "ptypes.h"

typedef enum {
  FS_INIT = 0,
  FS_TRIAL,
  FS_POWER,
  FS_LARGE,
  FS_TERM,
} fs_state_t;

/* Max number of factors on the unfactored stack, not the max total factors.
 * This is used when we split n into two or more composites.  Since we work
 * on the smaller of the composites first, this rarely goes above 10 even
 * with thousands of non-trivial factors. */
#define MAX_FACTORS 128

typedef struct factor_state_s {
  fs_state_t state;
  mpz_t n;  /* remaining number to be factored */
  mpz_t f;  /* new factor found */
  int e;    /* new exponent found */
  int ef;   /* exponent multiplier */
  UV tlim;  /* p^2 limit checked by trial division */

  /* used only for trial division phase */
  UV sp;  /* smallprime index */
  UV un;

  /* used only after trial division phase */
  int log;    /* verbose_level */
  int ntofac; /* number of additional factors in tofac_stack[] */
  mpz_t tofac_stack[MAX_FACTORS];
} factor_state;

extern void _init_factor(void);

extern int factor(mpz_t n, mpz_t* factors[], int* exponents[]);
extern int factor_one(factor_state* fs);
extern void clear_factors(int nfactors, mpz_t* pfactors[], int* pexponents[]);

extern int omega(mpz_t n);
extern int bigomega(mpz_t n);
extern void sigma(mpz_t res, mpz_t n, UV k);
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

extern UV   power_factor(mpz_t n, mpz_t f);

extern mpz_t* divisor_list(int* ndivisors, mpz_t n);

extern int is_smooth(mpz_t n, mpz_t k);
extern int is_rough(mpz_t n, mpz_t k);
extern int is_powerful(mpz_t n, uint32_t k);
extern int is_almost_prime(uint32_t k, mpz_t n);

#endif
