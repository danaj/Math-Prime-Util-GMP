SIQS standalone health checks
=============================

These optional developer checks exercise the current SIQS implementation.
They do not build/load the Perl XS module, run performance benchmarks, or
attempt complete large factorizations. They are not part of ordinary
"make test". Production sources and benchmark builds are not modified.

After "perl Makefile.PL":

  make check-siqs
  make test-siqs                          # alias
  make check-psiqs                        # adds --threaded
  make check-siqs SIQS_CHECK_ARGS='--extended --verbose'
  make check-siqs SIQS_CHECK_ARGS='--sanitize'
  make check-psiqs SIQS_CHECK_ARGS='--extended'

The Make target uses only a core-Perl (5.10+) build wrapper, the C compiler,
and GMP.
check-psiqs runs the same suites with the optional pthread Lanczos/worker checks
enabled by default. Both targets share SIQS_CHECK_ARGS for additional options.
The wrapper compiles one checker into a temporary directory and removes it
afterward. It does not require Math::Prime::Util or Math::Prime::Util::GMP to
be installed or built. Use --keep to retain the binary.

Direct build without Perl or a generated Makefile, from the repository root:

  cc -O3 -DSTANDALONE -o /tmp/siqs-check tools/siqs-check.c \
    prime_iterator.c squfof126.c pbrent63.c -lgmp -lm
  /tmp/siqs-check --suite sieve
  /tmp/siqs-check --extended --verbose

Add -march=native only for a binary intended for the local machine.
Add -DPSIQS -pthread to include the optional Lanczos/worker pool checks.
Do not compile lanczos.c separately: the matrix suite includes it to inspect
private packing/kernels, just as the main checker includes siqs.c.

Suites and options
------------------

One executable runs the named suites ("sieve", "relations", "matrix", "workers"). The default
"all" selection runs every registered suite. --list describes available
suites; --suite selects one. --extended adds more fixtures, policy endpoints,
and polynomials. --verbose reports individual polynomial/matrix fixtures. New suites can
be registered without creating separate executables or a large framework.

  perl tools/siqs-check.pl --list
  perl tools/siqs-check.pl --suite sieve --extended
  perl tools/siqs-check.pl --suite relations
  perl tools/siqs-check.pl --suite matrix --threaded --extended
  make check-psiqs SIQS_CHECK_ARGS='--suite workers'
  perl tools/siqs-check.pl --block-size 32768
  perl tools/siqs-check.pl --block-size 0
  perl tools/siqs-check.pl --sanitize
  perl tools/siqs-check.pl --sanitize --debug

Block size is a production compile-time setting, so testing another maximum
requires another build. A single build exercises both blocked and unblocked
contexts where its activation rules permit. Zero disables blocking entirely.
The runner supports --cc, --cflags, and --ldflags for other compilers/GMP
installations. No special CPU architecture is required by the checker.

Sanitizers need compiler/runtime support. --sanitize preserves release-style
padded stores; --debug separately enables SIQS_DEBUG's additional assertions
and its bounds-checked padding behavior. Both are useful and are different
checks. A run returns 0 on success, 1 on a failed check, or 2 on invalid checker
arguments. Build errors and fatal engine diagnostics also return nonzero.

Sieve suite
-----------

The independent reference uses uint32_t logical scores and simple bounded
root progressions. It does not reuse production byte kernels, fixed-hit
tiers, block-root advancement, or candidate scanning. Equal roots contribute
once. A checked per-cell upper bound ensures the reference itself cannot
overflow. Every logical cell must agree with production and remain <=255.
An intentional 256-score fixture verifies that the reference detects overflow
that both byte implementations could hide; it is not a real-input failure.

Candidate positions, bias-adjusted scores, initial hit links, and narrow/wide
candidate-map entries are independently checked. Clearing and reusing both
maps is checked, as is preservation of the original polynomial roots during
blocked sieving. Padding is allowed to wrap; it is never treated as a logical
score or scanned for candidates. Sanitizers check the actual padded stores.

Synthetic fixtures cover word-scan tails, every fixed-hit tier boundary,
equal roots, roots at zero/p-1, the localized one-hit crossover, activation
thresholds, balanced blocks, partial final blocks, and workspace reuse.
Periods need not be prime to test the storage kernels. The extended mode
also covers factor-base array counts crossing 65535/65536.

Real fixtures use deterministic semiprimes and production parameter/base/A/B
setup. Their roots are independently checked by evaluating the GMP polynomial.
Coverage includes actual q=1/q=2 sieving, d=1/d=2, primes in A, primes dividing
k, later Gray-code polynomials, later families, and larger blocked intervals.
Small fixtures are required to survive the normal public trial-division bound;
setup must not shortcut on a factor found in the base. The extended grid goes
through 370 bits, plus an intentionally long interval on a smaller input.
No solver, pthread pool, or full-factorization run is used by this suite.

The maximum reported score describes the sampled polynomials, not a proof
for every input in the supported range. A failed check reports its suite,
fixture, first discrepancy, and context geometry for reproduction.

Relations suite
---------------

Small valid congruences modulo 35 are constructed with independently computed
products and exhaustive square-root searches. Every full relation is checked
as y^2 == product(p^exponent) mod N, including the sign row. Tests cover smooth
vector ownership transfer, repeated 1LP pairs sharing their original anchor,
duplicate rejection, two-vertex and longer cycles, self-loops, and overlapping
cycles in forward/reverse/interleaved forest insertion orders. These insertion
orders model different worker merge orders without running polynomial workers.
Long paths exceed the initial 64-entry cycle buffers and check their reuse.

Inverse failures must discover the factor, not fabricate a full relation.
Cycle-closing and rejected raw allocations must rewind, while non-last frees
must not rewind another live allocation. Packed row/exponent storage limits
are checked separately. Extended checks also grow/reallocate anchor hash tables
and raw arena blocks. Sanitizers check cleanup of retained anchors/forest edges.
Composite/repeated divisor insertions must preserve the multiplicative partition
and invalidate cached primality where necessary.

The earlier exponent fixtures are retained: exact uint32_t limits, repeated
overflow rejection without mutation, scratch reset, and a valid later cycle or
dependency after an overflowing one. Artificial wide raw/full exponents test
arithmetic limits; they do not claim such an overflow occurred on a real input.

Matrix suite
------------

Relation exponent parity is checked against the generated columns (including
the sign row and zero-weight columns). An independent fixed-point pruning oracle
checks singleton chains and heavy-column trimming. Retained columns must still
match their original vectors and IDs; dependencies are remapped and verified
against the immutable original matrix, including removed columns as zero lanes.

Scalar original-column multiplication/transposition checks the packed/unpacked
kernels, dense input words, retained/post-row bitmaps, and aliased symmetric
products. Packing tests straddle 1023/1024 active rows and 32768/32769 columns.
Dense/Lanczos tests straddle the 1536/1537-column dispatch boundary; empty,
zero-rank, rank-deficient, and no-kernel dense cases are also checked. Every
returned dependency must be nonempty, linearly independent, and annihilate the
original matrix. Both normal and retain-all-row Lanczos paths are exercised.

--threaded compiles with PSIQS and pthreads. It checks 2/3/4-thread symmetric
products, inner products, and masked vector accumulation against scalar
references on deliberately uneven column weights. Fixed seeds must produce
the same serial/threaded solver output; packed/small matrices also exercise
the serial fallback. The extended suite includes an actually threaded solve
above 32768 columns. These are correctness checks, not scaling benchmarks.

No production code changes are required for these relation/matrix tests.

Workers suite
-------------

Enabled by check-psiqs or --threaded; otherwise it explicitly reports SKIP.
Small genuine polynomial jobs check 1/2/4-worker reuse, independent scratch,
distinct A families, exactly-once merges, and reuse after a forced wide
candidate-map allocation. Buffered GMP records must be cleared before their
arena is reset. Normal blocks are retained; older/oversized blocks are released.
Synthetic buffers check deep-copy ownership, cross-worker cycles, and factor
discovery during a partial merge or without any emitted relation.

Test-only pthread gates create a completed/READY worker, a running worker,
an assigned job not yet started, and an idle worker at shutdown. Gates do not
use sleeps or rely on which core the OS chooses. Workers must join before
their results are drained or scratch freed. A short final polynomial budget
checks reserved versus executed work and the collector's serial fallback.
Another gated case requests a stop after 32 polynomials and checks
that both exit at the next stop poll, with valid partial relation buffers.
Those buffers merge exactly once; retry collection keeps the existing context,
chooses fresh A values instead of replaying the abandoned tails, and uses only
the unspent polynomial budget. This tests retry collection, not an injected
failed matrix solve. No production scheduling gates or timers are added.

Test-only wrappers inject top-level allocation failures, mutex/condition
initialization failures, and selective or total pthread_create failures.
The regular pool must retain successfully started workers, skip failed slots,
and join/clear initialized resources exactly once. The extended suite requests
256 slots but permits only two real workers to start; it does not launch 256
threads. Failed-slot scratch must be released immediately, while its condition
is kept for teardown. Scratch ownership is independent of the started flag;
joining clears thread lifetime but must not lose ownership or cause a double
clear. The initial allocation peak is unchanged, and worker-scratch allocation
failure still follows the existing fatal policy.
Retained Lanczos failure cases instead require serial fallback and fixed-seed
output identical to serial. Resource counters must return to zero.

Repeated public 1/2/4-thread calls and two simultaneous callers verify exact
prime partitions after one host prime-cache startup. Startup/shutdown ownership
is unchanged. On POSIX hosts, child processes verify invalid counts 0/257 and
the existing fatal worker-scratch allocation policy: stderr diagnostics and
exit status 3, with no stdout. Those cases explicitly skip on Windows.

A POSIX test-only 120-second watchdog turns hangs into failures. It is not a
production timer; override SIQS_CHECK_WORKER_TIMEOUT in --cflags for a slow
host or set it to zero when debugging. Fatal-case children have a separate
20-second bound. The wrappers/scheduler are compiled only into the checker;
ordinary production source is not rewritten or exported for testing.

ASan/UBSan can be used as above. Where ThreadSanitizer is supported:

  make check-psiqs SIQS_CFLAGS='-O1 -g -fsanitize=thread' \
    SIQS_CHECK_ARGS='--suite workers'

Do not combine ThreadSanitizer with --sanitize (ASan). Normal worker jobs do
not acquire the allocation-injection lock, to avoid hiding races. These are
regression fixtures, not an exhaustive proof of every possible schedule.
The embedding/prime-cache lifetime contract is documented in siqs.h and
tools/README-siqs-embedding.txt. The simultaneous-call checks use that contract:
one startup, stable cache during all callers, then one shutdown after joins.
