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

/* Host lifetime for both entry points:
 *   - The shared prime-iterator cache must already be initialized. Standalone
 *     C hosts include prime_iterator.h and call prime_iterator_global_startup
 *     once before any calls, then prime_iterator_global_shutdown only after
 *     all calls/iterators finish and all caller threads have joined.
 *   - Startup/shutdown are not reference-counted or safe to race with calls.
 *     Do not repeat startup while the cache is live. These SIQS functions do
 *     not manage it. Existing MPU-GMP hosts/drivers already own this lifetime.
 *   - Standalone calls have independent writable state and may run concurrently
 *     with reentrant siqs_dep.h host adapters and separate output storage.
 *     Inputs may be shared read-only, but must not be mutated/cleared during
 *     calls. Keep adapter configuration immutable or properly synchronized.
 * See tools/README-siqs-embedding.txt for adapters, builds, and an example. */

/* n must be positive.  Return an allocated multiplicative partition of n.
 * Partition elements are not necessarily prime, and any partial splitting is
 * retained.  If no split is found, the sole element is n.  The SIQS stage is
 * attempted only when the post-trial cofactor is between MPU_SIQS_MIN_BITS
 * and MPU_SIQS_MAX_BITS, inclusive.
 *
 * trial_start is the first candidate not already checked for small factors;
 * use 0 for ordinary pretests. nfactors must be non-NULL and receives the
 * array length, including repeated factors. n is not modified. The caller
 * owns all returned mpz_t values and the array: release both with
 * gmp_siqs_free, without separately clearing/freeing its elements first.
 * There is no complete-factorization flag; callers must check the partition.
 * verbose is independent of host/module settings: 0 is quiet, 1 prints
 * setup/final summaries, 2 adds periodic progress, and 3 adds diagnostics.
 * Standalone fatal allocation/invariant errors exit rather than return NULL. */
extern mpz_t *gmp_siqs(const mpz_t n, uint32_t *nfactors,
                      uint32_t trial_start, int verbose);
extern void gmp_siqs_free(mpz_t *factors, uint32_t nfactors);

#ifdef PSIQS
/* Same partition and trial_start contract as gmp_siqs.  Release the result
 * with gmp_siqs_free.  nthreads must be in [1, PSIQS_MAX_THREADS]; one worker
 * uses the ordinary serial path. Requires a STANDALONE/PSIQS pthread build;
 * native workers must not call unaudited Perl host adapters. Concurrent outer
 * callers each create their own pool; nthreads is not a process-wide cap. */
extern mpz_t *gmp_psiqs(const mpz_t n, uint32_t *nfactors,
                       uint32_t trial_start, int verbose, uint32_t nthreads);
#endif

#endif
