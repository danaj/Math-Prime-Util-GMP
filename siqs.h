#ifndef MPU_SIQS_H
#define MPU_SIQS_H

#include <gmp.h>
#include "ptypes.h"

/* While we work for small sizes, performance under 36-bits is suboptimal. */
#define MPU_SIQS_MIN_BITS   1U
/* 431 is a hard cap for this implementation.  No tuning done > 330 bits. */
#define MPU_SIQS_MAX_BITS 370U

/* Parallel collection is optional; small/inline-solving policies stay serial. */
#define PSIQS_MAX_THREADS 256U
#ifndef PSIQS_SERIAL_MAX_BITS
# define PSIQS_SERIAL_MAX_BITS 64U
#endif

/* n must be positive.  Return an allocated multiplicative partition of n.
 * Partition elements are not necessarily prime, and any partial splitting is
 * retained.  If no split is found, the sole element is n.  The SIQS stage is
 * attempted only when the post-trial cofactor is between MPU_SIQS_MIN_BITS
 * and MPU_SIQS_MAX_BITS, inclusive.
 *
 * trial_start is the first candidate not already checked for small factors;
 * the returned array must be released with gmp_siqs_free. */
extern mpz_t *gmp_siqs(const mpz_t n, uint32_t *nfactors,
                      uint32_t trial_start);
extern void gmp_siqs_free(mpz_t *factors, uint32_t nfactors);

#ifdef PSIQS
/* Same partition and trial_start contract as gmp_siqs.  Release the result
 * with gmp_siqs_free.  nthreads must be in [1, PSIQS_MAX_THREADS]; one worker
 * uses the ordinary serial path. */
extern mpz_t *gmp_psiqs(const mpz_t n, uint32_t *nfactors,
                       uint32_t trial_start, uint32_t nthreads);
#endif

#endif
