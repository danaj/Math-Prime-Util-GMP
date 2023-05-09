#ifndef MPU_GMPPRIMALITY_H
#define MPU_GMPPRIMALITY_H

#include <gmp.h>
#include "ptypes.h"

extern int is_pseudoprime(mpz_t n, mpz_t a);
extern int is_euler_pseudoprime(mpz_t n, mpz_t a);

extern int  miller_rabin(mpz_t n, mpz_t a);
extern int  miller_rabin_ui(mpz_t n, unsigned long a);
extern int  miller_rabin_random(mpz_t n, UV numbases, char* seedstr);

extern int  _GMP_is_lucas_pseudoprime(mpz_t n, int strength);
extern int  _GMP_is_almost_extra_strong_lucas_pseudoprime(mpz_t n, UV incr);
extern int  _GMP_is_frobenius_underwood_pseudoprime(mpz_t n);
extern int  _GMP_is_frobenius_khashin_pseudoprime(mpz_t n);
extern int  is_perrin_pseudoprime(mpz_t n, int restricted);
extern int  is_euler_plumb_pseudoprime(mpz_t n);
extern int  is_frobenius_pseudoprime(mpz_t n, IV P, IV Q);
extern int  is_frobenius_cp_pseudoprime(mpz_t n, UV ntests);

extern int lucas_lehmer(UV p);
extern int llr(mpz_t N);
extern int proth(mpz_t N);
extern int is_proth_form(mpz_t N);

extern int _GMP_BPSW(mpz_t n);
extern int is_deterministic_miller_rabin_prime(mpz_t n);  /* assumes n is BPSW */
extern int  is_miller_prime(mpz_t n, int assume_grh);
extern int is_bpsw_dmr_prime(mpz_t n);

extern int  _GMP_is_prime(mpz_t n);
extern int  _GMP_is_prob_prime(mpz_t n);
extern int  _GMP_is_provable_prime(mpz_t n, char ** prooftext);

#endif
