SIQS BYTE-SIEVE BLOCK-SIZE SCREEN
===============================

Build only (regenerate Makefile after pulling the new target if necessary;
Perl module tests/build do not run this tool):
  make siqs-sieve-bench

Or compile without Perl/MakeMaker:
  cc -O3 -march=native -DSTANDALONE -DPSIQS -pthread \
    -o siqs-sieve-bench tools/siqs-sieve-bench.c lanczos.c \
    prime_iterator.c squfof126.c pbrent63.c -lgmp -lm

Run on an idle machine:
  ./siqs-sieve-bench
  ./siqs-sieve-bench --threads 4 --output blocks-4.tsv

For a short staged screen and advisory compiler flags, see siqs-autotune.pl
and README-siqs-autotune.txt. Ordinary sweeps below remain available unchanged.

Default input is RSA-100 (330 bits). An alternative decimal input
can be supplied as the final argument. Default sizes are off, 16/32/64 KiB,
with three alternating-order sweeps and roughly one second per timed batch.
The tool does not sweep thread counts: choose the intended sieve concurrency.
One thread is genuinely serial; at larger counts the caller is one sieve worker
and a persistent pool supplies the others. Thread creation is outside timing.

Custom block maxima need not be powers of two:
  ./siqs-sieve-bench --threads 4 --sizes 0,12,14,16,18,24,32,64,96,128,256

IMPORTANT: The maximum is not the actual balanced block size. Production's
activation rule is retained by default: logical interval >= 2.5 x maximum.
The default input's interval is 582 KiB with current parameters, so 16/32/64 KiB
are all active without any threshold selection by the user. Actual sizes and
dense-prime counts are printed. Geometrically identical configurations reuse
the same measurements; all sizes share one executable, avoiding separate-build
code-layout differences. Off is always included as the control. Block size also
changes the dense/sparse prime split.

For comparison with the smaller 250-bit example:
  ./siqs-sieve-bench --threads 4 \
    1401811817899460116600945074728583412740519573015376930481203561750251051823

Its current interval is 71.5 KiB: the 32/64 KiB maxima fall back to unblocked
under the default rule. Both 12/14 KiB maxima produce six balanced blocks of
at most 12203 bytes. For our later studies only, --min-length K supplies an
advanced override in KiB; zero bypasses the threshold for size-only screening.
Intervals fitting within a single block always stay unblocked. The tool does
not tune the activation threshold automatically. Kernel winners alone are not
prescriptions for that crossover, which also depends on interactions outside
this isolated replay.

Timed work is the actual siqs_run_sieve: physical byte initialization, blocked
root copies/advancement, dense sieving, and the full-interval sparse suffix.
Multiplier selection, factor-base/polynomial preparation, next_B root updates,
candidate finding, resieving, evaluation, relation handling and solving are NOT
timed or performed in the replay. There is no full factorization.

Every worker owns normal production scratch and one distinct, genuinely prepared
A/polynomial. Its roots remain fixed for this initial cache-locality screen.
Repeated full clears prevent accumulation between sieve calls; a rotating byte
read makes every call observable. This is a warmed kernel replay, not a model
of the full factoring cache mix or the full distribution of polynomials. Confirm
promising settings end-to-end on several inputs/bit sizes before changing defaults.

The harness includes siqs.c and uses its allocator/balancing and sieve kernels
unchanged. Runtime block settings exist only in the harness; no production hooks,
source edits or policy changes are needed. Before timing and after every batch,
every logical sieve byte is compared to an unblocked reference for that worker.
Padding bytes are intentionally not compared. New block workspaces, correctness
checks and warmups are outside the timers.

Wall uses CLOCK_MONOTONIC; aggregate process CPU uses getrusage. Dispatch/wait
overhead is amortized over an entire batch, never one barrier per polynomial.
Calibration adjusts calls per worker separately for each distinct geometry, so
compare normalized wall microseconds per completed sieve (throughput), NOT raw
batch times. CPU microseconds per sieve indicates aggregate processor cost.
--polynomials N skips calibration for a reproducible fixed-call comparison;
--seconds changes the default calibration target; --repeat changes sweep count.

--output writes and flushes a TSV after each measured batch and refuses to replace
an existing file. Alias configurations appear in the console summary, not as
duplicate fabricated measurements in the raw TSV. The tool reports the fastest
measured geometry but never changes a build configuration or production setting.
Near ties inside the observed spread should not drive tuning. Run several inputs
and enough repeats; SMT siblings share core/cache resources. Threads are not
pinned, and the tool does not manage NUMA placement. Do not run alongside
the factoring queue or other heavy work when collecting performance results.

Inputs must be odd, non-power composites of 100..MPU_SIQS_MAX_BITS bits without
factor-base divisors. Setup retains polynomial scratch per worker, so large
factor bases and very high thread counts require substantial memory. Requires
GMP, POSIX pthreads, clock_gettime(CLOCK_MONOTONIC), and getrusage; no MPU Perl
installation or external benchmarking dependencies are required.
