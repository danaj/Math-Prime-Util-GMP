Embedding SIQS/PSIQS in a C host
===============================

This describes the current in-tree C interfaces, not a separately packaged
library or a new initialization wrapper. Public declarations are in siqs.h;
standalone host adapters are in siqs_dep.h. mpu-siqs.c is the reference driver
and adapter implementation. Its normal command-line startup already handles
the lifecycle below. The Perl/MPU-GMP host also has its own initialization.

Shared prime-cache lifetime
--------------------------

The prime iterator uses process-global prime/sieve caches. SIQS does not own,
initialize, reference-count, or shut down those caches for each factor call.

A standalone embedding must include prime_iterator.h and have one host owner:

  1. Call prime_iterator_global_startup() once before any SIQS/prime-iterator
     use and before launching caller threads.
  2. Keep that cache alive and unchanged across repeated or concurrent calls.
  3. Finish all calls, join all caller threads, and destroy any other remaining
     prime_iterator instances before prime_iterator_global_shutdown().
  4. Shut down once when that owner's cache lifetime ends.

Do not initialize/shut down around each concurrent call or from an individual
worker, primality callback, or splitter callback. Repeating startup while live
is not safe/idempotent; it replaces shared allocations rather than acquiring
another reference. Shutdown invalidates the cache for every user in that
process. Concurrent startup/shutdown and cache reads are not supported.

If another MPU-GMP host already owns the cache, share its established lifetime
instead of initializing it again. For example, _GMP_init/_GMP_destroy manage
the cache along with other MPU-GMP state; do not independently start/stop it
inside that lifetime. This does not make the entire MPU-GMP API thread-safe.

The standalone API has per-call writable sieve/graph/solver state. Separate
gmp_siqs or gmp_psiqs calls may run concurrently under the stable cache lifetime
with reentrant adapters. Each caller must supply its own count output and own
its returned array. Inputs may be shared read-only, but must not be mutated
or cleared while calls use them. Returned mpz_t values are independent copies,
not borrowed inputs or prime-cache storage.

A gmp_psiqs call joins its own workers before returning. The outer host still
has to wait for its own caller threads before shutdown. Each call's nthreads
is a pool size, not a global scheduler/core limit: simultaneous calls can
oversubscribe the machine and consume the sum of their private scratch.

Standalone host adapters
------------------------

Compile a standalone embedding with STANDALONE, without STANDALONE_ECPP.
Supply these three definitions with the exact siqs_dep.h prototypes:

  int siqs_verbose_level(void);
    Return a stable verbosity value. Zero is quiet; SIQS starts progress
    output above 2. Keep configuration immutable during calls or synchronize
    access. Concurrent progress output can interleave on shared stdout/stderr.

  int siqs_is_prob_prime(const mpz_t n);
    Return nonzero for prime/probably-prime, zero for composite; do not modify n.
    The current standalone driver uses mpz_probab_prime_p(n, 25).

  int siqs_pbrent_factor(const mpz_t n, mpz_t f, UV a, UV rounds);
    n is read-only and f is already initialized by SIQS. Perform a bounded
    Pollard-Brent attempt using a and rounds. Return nonzero only after setting
    a proper divisor 1 < f < n that divides n; return zero on a miss. Do not
    clear n or f. Use local temporaries/RNG state, not shared writable scratch.

The driver contains a portable GMP implementation of the third adapter.
Reuse/adapt its three host functions in your host source; do not link its main
alongside your own main. Include siqs_dep.h in the adapter translation unit.
The UV type comes from ptypes.h; use the same build flags/types across units.

Callbacks may execute in native worker threads. They must be reentrant and
must not call unaudited Perl APIs/error machinery. Normal non-standalone builds
map these names to MPU-GMP functions; STANDALONE_ECPP deliberately retains that
full-host mapping. Those modes are not the lightweight standalone adapter mode,
and PSIQS currently requires the standalone host.

Do not throw/unwind/longjmp out of a worker callback to attempt recovery.
Standalone allocation/invariant errors and invalid API worker counts are fatal:
the current croak prints to stderr and exits the process with status 3. There
is no recoverable allocation-error/NULL-result contract. Partial thread creation
instead uses the successfully created relation workers, or serial collection
if none started. This is not general out-of-memory recovery.

Calling and owning results
--------------------------

  factors = gmp_siqs(n, &count, trial_start);
  factors = gmp_psiqs(n, &count, trial_start, nthreads);  /* PSIQS build */

n must be an initialized positive mpz_t; count must point to caller-owned
uint32_t output storage. n is not modified. Use trial_start=0 for ordinary
pretests. A larger value means smaller trial candidates have already been
checked; do not skip unperformed trial division accidentally.

A normal return gives an allocated multiplicative partition, not a promise
of complete prime factorization. count is its length, including repeated
factors. Do not assume ordering or use count>1 as a completeness test. Check
that the product is n and test each element's primality as appropriate for
your application. n=1 returns the identity partition [1]. An unsplit composite,
unsupported SIQS cofactor size, or missed split can leave composite elements;
partial splits are still retained.

The size limits in siqs.h control the SIQS stage on the post-trial cofactor,
not a blanket guarantee that all original inputs below the limit split, nor
that larger inputs cannot return factors through pretests.

Release the array and every owned mpz_t together:

  gmp_siqs_free(factors, count);

Do not clear its elements first, call plain free on the array, or use it after
release. The same free function applies to both entry points. If a host wants
complete factoring, it must iteratively refine composite elements and enforce
progress/bounded retries; never keep retrying an unchanged composite forever.
The mpu-siqs command-line driver handles this separately.

Minimal caller example
----------------------

Save the following as my-factor.c. It requires the host-adapter source
described above; the callback definitions are intentionally not duplicated
here. It reports an incomplete partition rather than silently labelling it
a complete factorization.

  #include <stdio.h>
  #include <gmp.h>
  #include "siqs.h"
  #include "prime_iterator.h"

  int main(int argc, char **argv) {
    mpz_t n, product;
    mpz_t *factors;
    uint32_t count, i;
    int complete = 1, status = 0;

    if (argc != 2) {
      fprintf(stderr, "usage: my-factor POSITIVE_INTEGER\n");
      return 2;
    }
    mpz_init(n);
    if (mpz_set_str(n, argv[1], 10) != 0 || mpz_sgn(n) <= 0) {
      fprintf(stderr, "invalid positive decimal integer\n");
      mpz_clear(n);
      return 2;
    }

    prime_iterator_global_startup();
  #ifdef PSIQS
    factors = gmp_psiqs(n, &count, 0, 4);
  #else
    factors = gmp_siqs(n, &count, 0);
  #endif
    mpz_init_set_ui(product, 1);
    for (i = 0; i < count; i++) {
      mpz_mul(product, product, factors[i]);
      if (mpz_cmp_ui(factors[i], 1) > 0 &&
          mpz_probab_prime_p(factors[i], 25) == 0)
        complete = 0;
    }
    if (mpz_cmp(product, n) != 0) {
      fprintf(stderr, "invalid returned factor partition\n");
      status = 3;
    } else {
      for (i = 0; i < count; i++)
        gmp_printf("%s%Zd", i ? " " : "", factors[i]);
      printf("%s\n", complete ? "" : " [incomplete]");
      status = complete ? 0 : 1;
    }

    gmp_siqs_free(factors, count);
    mpz_clear(product);
    mpz_clear(n);
    prime_iterator_global_shutdown();
    return status;
  }

Build from the repository root, supplying your three adapters in my-siqs-host.c:

  cc -O3 -DSTANDALONE -I. -o my-factor my-factor.c my-siqs-host.c \
    siqs.c lanczos.c prime_iterator.c squfof126.c pbrent63.c -lgmp -lm
  ./my-factor 22095311209999409685885162322219

For parallel support, add -DPSIQS -pthread to that command. Both public entry
points are then available; the example requests four workers. Valid counts
are 1..PSIQS_MAX_THREADS. One worker and small/inline-solving policies use the
serial path. Do not compile psiqs_inc.c or planczos_inc.c separately; they
are included by siqs.c/lanczos.c. Use -march=native only for a local-machine
binary, not one intended for arbitrary other CPUs.

The example is one synchronous call. For repeated inputs, keep startup before
the whole loop and shutdown after it. For concurrent inputs, put startup before
launching caller threads; each thread uses its own input/count/result storage
(or shares only read-only input), and the host joins them all before shutdown.

Validation and scope
--------------------

  make check-psiqs SIQS_CHECK_ARGS='--suite workers --extended'

The persistent worker suite exercises one cache startup, repeated public
1/2/4-thread calls, two simultaneous callers with independent pools, result
ownership, and one shutdown after joins. It also checks reduced/failed pool
creation and cleanup. It is a regression check, not exhaustive proof of every
schedule or a claim that all other MPU-GMP routines are reentrant.

This guide documents the existing contract. It adds no global reference
counter, hidden per-call initialization, new embedding wrapper, or cleanup/
error-policy changes.
