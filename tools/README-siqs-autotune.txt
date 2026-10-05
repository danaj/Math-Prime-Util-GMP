OPTIONAL SIQS HARDWARE SCREEN
============================

Build the shared helper, then run the advisory tool:

  make siqs-sieve-bench
  tools/siqs-autotune.pl
  tools/siqs-autotune.pl --flags --threads 8 --output tune-8.tsv

Or measure and rebuild msiqs in one step:

  make siqs-tuned
  make siqs-tuned-8
  make siqs-tuned SIQS_TUNE_ARGS='--threads 8'

The eight-thread forms screen for eight threads; they do not change msiqs's
default thread count. Still pass -threads 8 when running that executable.
The siqs-tuned-8 alias prints its recursive make command as a usage example.

Compilation stays in the Makefile / existing standalone builder. siqs-tuned is
a phony override: it always measures and rebuilds msiqs, even if the executable
already exists. Failed measurements prevent compilation. Ordinary make siqs
does not run tuning. The Perl module is not rebuilt by these targets.
CC/SIQS_CFLAGS/SIQS_ARCH_FLAGS retain their normal meanings; EXTRA_FLAGS adds
manual options to standalone builds. No non-core Perl modules are used.
Flags are not saved: a later ordinary rebuild uses the ordinary build options.

OUTPUT MODES
------------

  tools/siqs-autotune.pl --flags
  tools/siqs-autotune.pl --onlyflags
  make siqs EXTRA_FLAGS="$(tools/siqs-autotune.pl --flags)"

--flags prints only the final compiler flags on stdout, with concise preparation,
stage lists, decision and serial-speed information on stderr. It is the mode
used by make siqs-tuned: progress stays visible during command substitution.
--onlyflags (alias --flags-only) also prints just flags on stdout, but suppresses
normal progress completely. Without either switch, stdout contains full
measurements, aggregate results and compiler flags. Errors still go to stderr.
The two flags-only modes are mutually exclusive.

For the manual make form, remove an existing msiqs first: make does not detect
variable changes. siqs-tuned always forces the rebuild. In shell scripts, check
the autotune command's exit status before invoking make; inline command
substitution alone does not propagate that failure.

The result line contains -DSIQS_SIEVE_BLOCK_SIZE and/or
-DSIQS_PROGRESS_WORK_PER_SEC definitions. No partial flags are emitted after
failed measurements or malformed helper output. --bench supplies another
measurement executable. --output creates an exclusive TSV, never overwriting
an existing file. Its stage column identifies serial/threaded/final/extended;
n identifies the input, and threads is the actual measured worker count.
Extended rows have repeat indices 3 and 4 because they add to the first two
final samples. Geometrically identical measurements are reused, not fabricated
as duplicate measured rows.

SCREEN POLICY
-------------

Default inputs: RSA-100/RSA-110 (330/364 bits).
Default block maxima: 0/16/32/48/64/96/128 KiB; --sizes supplies another list.
Optional decimal INTEGER arguments REPLACE the defaults: supply 1..16 distinct
odd non-power composites of 100..MPU_SIQS_MAX_BITS bits, without FB divisors.
Duplicate integers are rejected so an accidental repeat cannot bias weighting.
Default --threads is 1. No automatic core-count, affinity, NUMA or heterogeneous
core selection. Use an idle machine in a representative power/thermal state.

Each stage measures every surviving geometry on EVERY input. For each input,
use the median wall microseconds per completed sieve, normalize by its fastest
surviving choice, then take an equal-weight geometric mean across inputs.
Thus larger inputs do not dominate simply because their sieves take longer.
Cull/select from that aggregate, not a vote among per-input winners: a
compromise can win without being individually best anywhere.

1. Approximately 200 ms per distinct geometry/input, serial. Drop aggregate
   choices more than 15% slower than the leader; keep at least two.
2. Approximately 400 ms per survivor/input at --threads. Drop aggregate choices
   more than 5% slower; keep at least two.
3. Two approximately 400 ms final sweeps per finalist/input at --threads.
4. If the leading aggregates are within 3%, or leading contenders show >5%
   timing spread, do one longer confirmation: two 1-second sweeps per retained
   contender/input. Keep choices within 5% of the leader (at least two).
   Combine all four final samples using their median. Never extend again.

Choose the smaller maximum among choices within 2% of the fastest aggregate.
A decision line explains either the aggregate near-tie/smaller-block preference,
or the aggregate speed ratio against the runner-up plus per-input win count.
The ratio describes the aggregate; the win count describes consistency, not
statistical significance or end-to-end factoring speed. If the latest final
sweeps still vary by >5%, repeat on an idle machine.

Input and candidate sweep orders reverse to reduce systematic order bias.
A short pilot calibrates stage 1; later stages reuse observed per-worker time.
Durations are approximate, not hard deadlines. There is one persistent pool
of the requested size, shared across input scratch; inactive inputs do not
create extra threads. All prepared input scratch remains resident during the
screen, so large input sets or thread counts can require substantial memory.
Inputs are fully prepared outside the measurement timers.

The loose serial cut is a heuristic: serial losers can win under different
shared-cache/SMT pressure. No eliminated candidate is restored. For a more
exhaustive check, use the ordinary siqs-sieve-bench sweeps with longer/repeated
measurements. Every measured worker must finish its assigned work and every
logical sieve byte must match its input's unblocked reference.

The block maximum is NOT the balanced actual size; both are shown in full
output. Equal geometries on one input reuse that measurement, but candidates
are globally aliased only if their geometries match on ALL inputs.
Off stays a candidate and can win. With no active blocked geometry, no block
recommendation is emitted. Production balancing and activation at 2.5 x maximum
are retained: a different maximum also changes its derived activation length.
This tool does not tune that ratio. Confirm medium-size behavior separately
before making a hardware recommendation a general-purpose default.

This is a warmed byte-sieve replay on one genuine fixed polynomial per worker
and input, not the full factoring cache mix or polynomial distribution. Default
tuning takes roughly 10--30 seconds on the reference machine, plus compilation
when needed; larger lists, high concurrency or slower hardware cost more.
See README-siqs-sieve-bench.txt for details of the measured kernel.

SERIAL SPEED / PROGRESS WORK RATE
--------------------------------

After screening, release the sieve workers and use a fresh serial context on
the FIRST input with the selected geometry. The default two-second probe includes
polynomial/root generation, byte sieving, candidate scan, resieving, evaluation,
graph insertion and normal readiness checks. Setup, cleanup and solving are
excluded. Check time every 32 polynomials and at family boundaries, never in a
kernel. Natural readiness/factor discovery can stop small inputs early; very
short probes produce no work-rate recommendation.

"667M work units per CPU second" means roughly 667 million units, where a unit
is M x polynomials and M is the half-interval. SIQS_PROGRESS_WORK_PER_SEC is
rounded to a million. It is an early-collection estimate, not a throughput
promise: other sizes, large late graphs, background activity or core classes can
change cadence. Threaded reporting retains its existing worker-count scaling.
This flag ONLY changes verbose-output cadence, never factoring parameters.
SIQS_PROGRESS_OUTPUT_EVERY_NSECS remains the independent human-facing knob.

--work-seconds 0 skips speed calibration; a larger value (up to 60) extends it.
No production timers, signal handlers, pretests, band or activation-ratio changes
are installed. Source/build settings and TODO files are never edited by the tool.

VALIDATION
----------

Optional aggregate/alias/pool/protocol tests: make check-siqs-autotune.
Wrapper-only tests without compilation: perl xt/siqs-autotune.t.
Requires the existing GMP/POSIX-pthreads helper, independent of the Perl module's
normal build/test workflow.
