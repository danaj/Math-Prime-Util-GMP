#ifndef MPU_GMPPRIMALITY_H
#define MPU_GMPPRIMALITY_H

#include <gmp.h>
#include "ptypes.h"

extern int  _GMP_miller_rabin(mpz_t n, mpz_t a);
extern int  _GMP_is_lucas_pseudoprime(mpz_t n, int strength);
extern int  _GMP_is_almost_extra_strong_lucas_pseudoprime(mpz_t n, UV incr);
extern int  _GMP_is_frobenius_underwood_pseudoprime(mpz_t n);
extern int  _GMP_is_frobenius_khashin_pseudoprime(mpz_t n);
extern int  is_perrin_pseudoprime(mpz_t n);
extern int  is_frobenius_pseudoprime(mpz_t n, IV P, IV Q);
extern int  is_frobenius_cp_pseudoprime(mpz_t n, UV ntests);
extern int  _GMP_miller_rabin_random(mpz_t n, UV numbases, char* seedstr);

extern void lucas_seq(mpz_t U, mpz_t V, mpz_t n, IV P, IV Q, mpz_t k,
                      mpz_t Qk, mpz_t t);
extern void alt_lucas_seq(mpz_t U, mpz_t V, mpz_t n, IV P, IV Q, mpz_t k,
                          mpz_t Qk, mpz_t t);
extern void lucasuv(mpz_t Uh, mpz_t Vl, IV P, IV Q, mpz_t k);
extern int lucas_lehmer(UV p);
extern int llr(mpz_t N);
extern int proth(mpz_t N);
extern int is_proth_form(mpz_t N);

int _GMP_BPSW(mpz_t n);
int is_deterministic_miller_rabin_prime(mpz_t n);  /* assumes n is BPSW */
extern int  is_miller_prime(mpz_t n, int assume_grh);

extern int  _GMP_is_prime(mpz_t n);
extern int  _GMP_is_prob_prime(mpz_t n);
extern int  _GMP_is_provable_prime(mpz_t n, char ** prooftext);

#endif
