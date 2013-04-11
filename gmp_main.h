#ifndef MPU_GMPMAIN_H
#define MPU_GMPMAIN_H

#include <gmp.h>
#include "ptypes.h"

extern void _GMP_set_verbose(int v);
extern int  _GMP_get_verbose(void);
extern gmp_randstate_t*  _GMP_get_randstate(void);

extern void _GMP_init(void);
extern void _GMP_destroy(void);

extern int  _GMP_miller_rabin(mpz_t n, mpz_t a);
extern int  _GMP_is_strong_lucas_pseudoprime(mpz_t n);

extern UV   _GMP_trial_factor(mpz_t n, UV from_n, UV to_n);

extern int  _GMP_is_prime(mpz_t n);
extern int  _GMP_is_prob_prime(mpz_t n);
extern int  _GMP_is_provable_prime(mpz_t n, char * prooftext);
extern int  _GMP_is_aks_prime(mpz_t n);
extern void _GMP_next_prime(mpz_t n);
extern void _GMP_prev_prime(mpz_t n);

/* extern int _GMP_primality_pocklington(mpz_t n, int do_quick); */
extern int _GMP_primality_bls(mpz_t n, int do_quick, char ** prooftextptr);

extern int  _GMP_prho_factor(mpz_t n, mpz_t f, UV a, UV rounds);
extern int  _GMP_pbrent_factor(mpz_t n, mpz_t f, UV a, UV rounds);
extern int  _GMP_pminus1_factor(mpz_t n, mpz_t f, UV B1, UV B2);
extern int  _GMP_pminus1_factor2(mpz_t n, mpz_t f, UV rounds);
extern int  _GMP_holf_factor(mpz_t n, mpz_t f, UV rounds);
extern int  _GMP_squfof_factor(mpz_t n, mpz_t f, UV rounds);
extern int  _GMP_power_factor(mpz_t n, mpz_t f);

extern void _GMP_pn_primorial(mpz_t prim, UV n);
extern void _GMP_primorial(mpz_t prim, mpz_t n);
extern void _GMP_lcm_of_consecutive_integers(UV B, mpz_t m);

#endif
