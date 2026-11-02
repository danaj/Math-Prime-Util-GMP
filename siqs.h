#ifndef MPU_SIQS_H
#define MPU_SIQS_H

#include <gmp.h>
#include "ptypes.h"


/* Lowest post-trial cofactor size accepted by the SIQS stage. */
#define MPU_SIQS_MIN_BITS   1U
/* Arithmetic and policy-coverage limit, not a performance guarantee.
 * We have currently verified good performance through 400 bits. */
#define MPU_SIQS_MAX_BITS 512U


/* Host lifetime for both entry points:
 *   - No prime-cache startup/shutdown is required.  The independent pr32_iter
 *     uses immutable base tables and iterator-owned segment buffers.
 *   - Calls have independent writable state and may run concurrently with
 *     separate output storage. Inputs may be shared read-only, but must not
 *     be mutated or cleared during calls.
 * See README-SIQS.md for module use and background. */

/* n must be positive.  Return an allocated multiplicative partition of n.
 * Partition elements are not necessarily prime, and any partial splitting is
 * retained.  If no split is found, the sole element is n.  The SIQS stage is
 * attempted only when the post-trial cofactor is between MPU_SIQS_MIN_BITS
 * and MPU_SIQS_MAX_BITS, inclusive.
 *
 * trial_start is the first candidate not already checked for small factors;
 * use 0 for ordinary pretests.  nfactors must be non-NULL and receives the
 * array length, including repeated factors.  n is not modified.  The caller
 * owns all returned mpz_t values and the array: release both with
 * gmp_siqs_free, without separately clearing/freeing its elements first.
 * There is no complete-factorization flag; callers must check the partition.
 * verbose is independent of host/module settings; negative values become 0.
 * 0 is quiet, 1 prints setup/final summaries, 2 adds periodic progress,
 * and 3 adds diagnostics.
 * Standalone fatal allocation/invariant errors exit rather than return NULL. */
extern mpz_t *gmp_siqs(const mpz_t n, uint32_t *nfactors,
                      uint32_t trial_start, int verbose);
extern void gmp_siqs_free(mpz_t *factors, uint32_t nfactors);



#endif
