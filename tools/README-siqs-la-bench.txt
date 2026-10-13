Dense / block-Lanczos crossover on real SIQS matrices
===================================================

Build from the repository root:
  make siqs-la-bench
or, without the Perl build setup:
  cc -O3 -march=native -DSTANDALONE -UPSIQS -o siqs-la-bench \
    tools/siqs-la-bench.c siqs.c prime_iterator.c squfof126.c -lgmp -lm

Run on a QUIET machine, with unrelated factoring/benchmark queues paused:
  ./siqs-la-bench --repeats 6 < inputs.txt > crossover.tsv

Inputs are positive decimal composites, one per line; blank lines and comment
lines are allowed. Balanced semiprimes roughly 130..190 bits are useful around
the current 1536-column crossover. The tool rejects inputs above 245 bits to
avoid capturing recursive cofactor SIQS calls from the 2LP policies. It uses
the ordinary serial SIQS entry point with current production parameters and
primality adapter. No factor-base/interval settings or solver policy are changed.

Including lanczos.c with the reduction entry point renamed lets a tool-only
wrapper see each real matrix immediately after reduction. Both solvers see
identical independent copies. The normal SIQS solver then runs and factoring
continues, with its partition checked for the exact product and probable-prime
factors. Input preparation and full factoring are NOT the timings reported.
Failed input validation or incomplete factoring returns nonzero.
An input set with no matrices inside the selected column range also returns
nonzero. This is a POSIX tool (tested on macOS; clock APIs also exist on Linux).

Options:
  --repeats N    Solver repetitions per matrix (default 6, maximum 100).
                Even counts give each seed both dense-first/Lanczos-first orders.
  --min-cols N   Minimum eligible reduced column count (default 384).
  --max-cols N   Maximum eligible reduced column count (default 8192).
                Column-option limits are 32768; large dense solves may be slow.

No synthetic matrix model is used. Timing includes each public solver's normal
setup, allocations, internal retries and production verification. Matrix copy,
extra independent dependency verification and release are outside the timed
call. CPU uses CLOCK_PROCESS_CPUTIME_ID; wall uses CLOCK_MONOTONIC. Both report
per-matrix medians in milliseconds. Order alternates between repetitions and
inputs. A fixed benchmark seed, incremented once per pair of repetitions, is
independent of SIQS's own RNG. Threading is disabled for both solvers here.

Every returned dependency must be nonempty and satisfy the original sparse
matrix; masks must match the reported result bits. Different dependency counts
are legitimate. This is a solver-time comparison, not evidence that one set
of dependencies is more useful for factoring. Dense's production fallback to
Lanczos, if its basis gives trivial factors, is outside these solver timings.

TSV columns record input ordinal/bit size, solve ordinal, rows/columns/nonzeros,
repeat count, dense/Lanczos CPU and wall medians, Lanczos/dense CPU ratio, and
dependency counts. Ratio > 1 means dense was faster. Keep the input file, TSV,
compiler flags, binary/source hashes and machine description with the report.
Use multiple independent numbers near each size, group by ACTUAL reduced column
count, and look for a sustained change rather than fitting one matrix's result.
Do not subtract unrelated sieving work or use total factor time as solver time.

The benchmark adds no hooks to production SIQS/Lanczos and never changes
LA_DENSE_CROSSOVER_COLS. That policy is in lanczos.h, currently 1536 columns.
The older September 3 study is misc/siqs-dense-solver-results.md; its original
/tmp harness/data are no longer available, and it predates today's panel solver.
