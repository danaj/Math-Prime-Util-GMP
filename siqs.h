#ifndef MPU_SIQS_H
#define MPU_SIQS_H

#include <gmp.h>
#include "ptypes.h"

#define MPU_SIQS_MIN_BITS  40U
#define MPU_SIQS_MAX_BITS 350U

/* n must be positive.  Return an allocated multiplicative partition of n.
 * Partition elements are not necessarily prime, and any partial splitting is
 * retained.  If no split is found, the sole element is n.  The SIQS stage is
 * attempted only when the post-trial cofactor is between MPU_SIQS_MIN_BITS
 * and MPU_SIQS_MAX_BITS, inclusive.
 *
 * trial_start is the first candidate not already checked for small factors;
 * the returned array must be released with _GMP_siqs_free. */
extern mpz_t *_GMP_siqs(const mpz_t n, uint32_t *nfactors,
                        uint32_t trial_start);
extern void _GMP_siqs_free(mpz_t *factors, uint32_t nfactors);

#endif
