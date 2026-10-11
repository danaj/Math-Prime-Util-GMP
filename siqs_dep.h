#ifndef MPU_SIQS_DEP_H
#define MPU_SIQS_DEP_H

#include <gmp.h>
#include "ptypes.h"

/* Host adapters for SIQS.  A standalone C embedding supplies this definition
 * (mpu-siqs.c contains the reference implementations).
 * They must not manage the prime cache themselves.
 * For PSIQS/concurrent C calls, use reentrant adapters, caller-local
 * temporaries, and immutable or synchronized configuration; don't enter
 * unaudited Perl error/API machinery.
 * See tools/README-siqs-embedding.txt for the complete host contract.
 * The ECPP standalone build retains the full MPU-GMP host functions. */

#if defined(STANDALONE) && !defined(STANDALONE_ECPP)

/* Nonzero means prime/probably prime; do not modify n. */
extern int siqs_is_prob_prime(const mpz_t n);

#else  /* Below means we're in the Perl module or ECPP */

#include "primality.h"
#define siqs_is_prob_prime(n)        _GMP_is_prob_prime(n)

#endif

#endif
