#ifndef MPU_SIQS_DEP_H
#define MPU_SIQS_DEP_H

#include <gmp.h>
#include "ptypes.h"

/* Host adapters for SIQS and its solver. A standalone C embedding supplies
 * these three definitions (mpu-siqs.c contains the reference implementations).
 * They must not manage the prime cache themselves. For PSIQS/concurrent C
 * calls, use reentrant adapters, caller-local temporaries, and immutable or
 * synchronized configuration; don't enter unaudited Perl error/API machinery.
 * See tools/README-siqs-embedding.txt for the complete host contract.
 * The ECPP standalone build retains the full MPU-GMP host functions. */
#if defined(STANDALONE) && !defined(STANDALONE_ECPP)

/* SIQS starts progress output above 2; return 0 for quiet embedding. */
extern int siqs_verbose_level(void);
/* Nonzero means prime/probably prime; do not modify n. */
extern int siqs_is_prob_prime(const mpz_t n);
/* f is already initialized. Return nonzero only with a proper divisor in f;
 * return 0 when the bounded attempt misses. Do not clear n or f. */
extern int siqs_pbrent_factor(const mpz_t n, mpz_t f, UV a, UV rounds);

#else

#include "utility.h"
#include "primality.h"
#include "factor.h"

#define siqs_verbose_level()         get_verbose_level()
#define siqs_is_prob_prime(n)        _GMP_is_prob_prime(n)
#define siqs_pbrent_factor(n,f,a,r)  _GMP_pbrent_factor(n,f,a,r)

#endif

#endif
