BLOCK-LANCZOS SCALING BENCHMARK
==============================

Build only (does not run a large benchmark):
  make lanczos-bench

Or compile without Perl/MakeMaker:
  cc -O3 -march=native -DSTANDALONE -DPSIQS -pthread \
    -o lanczos-bench tools/lanczos-bench.c -lm

GMP headers must be available because the shared host header includes gmp.h;
the benchmark does not call GMP, SIQS, or any factorization routine. Use the
usual -I path if those headers are not in the compiler's search path. Requires
POSIX pthreads, clock_gettime(CLOCK_MONOTONIC), and getrusage.

Example full-size run:
  ./lanczos-bench --max-threads 64 --repeat 3 --output lanczos-110.tsv

Quick functional/small scaling check:
  ./lanczos-bench --rows 33000 --max-threads 8 --output lanczos-small.tsv

Sweep: 1, 2, 4, 8, ... up to the maximum, including the exact maximum if not
a power of two (e.g. --max-threads 10 includes 10). Counts include the caller,
as in the production solver. The first solve is serial. Repeated sweeps
alternate ascending/descending order; the summary gives median, min/max,
median CPU, and wall speedup against the serial median. TSV records every
solve and is flushed after each result. Existing output files are not replaced.
--verbose adds the solver's iteration diagnostics. --help lists all options.

The tool defaults to 183590 rows and 183654 columns: the reduced RSA-110
dimensions recorded in misc/siqs-timing-info.txt. It GENERATES an approximation,
not RSA-110's actual relations. The default 50 +/- 16 sparse-tail entries per
column, plus approximately nine frequent-head entries, are modelling choices,
not measurements of that real matrix. --weight changes the tail density,
--head changes the frequent-head size, and --seed generates another matrix.
For example, compare --weight 35, 50, and 70 to test sensitivity to density.

Tail hits mix uniform and lower-row-biased sampling; head incidence decays
with row number. Two tail entries per column ensure each tail row occurs at
least twice. Entries are distinct within each column; weights vary and columns
are ordered as after production's reduction. There are no planted duplicate
columns or known null vectors. This models a reduced sparse core rather than
the raw relation matrix; generation and pre-solver pruning are not timed.

Timed work includes the normal matrix conversion, pool creation, complete
block-Lanczos iteration, dependency extraction/verification inside the solver,
retries, and pool/storage teardown. It does not include matrix generation or
the harness's independent validation. Each solve uses identical matrix data
and solver seeds. Outside the timer, the harness checks A*X=0, nonempty and
independent dependency lanes, and bit-for-bit equality with the serial result.
This finds dependencies, not integer factors.

The harness includes lanczos.c to observe actual thread count and execute its
unchanged private recurrence. Its small setup/retry wrapper mirrors the public
entry point. No production sources, algorithms, or pool policies are changed.
The BENCHMARK build defaults NLA_MAX_THREADS to 256 so a 128/192-core machine
can explore higher counts; production's provisional 32-thread cap is unchanged.
Override the harness cap at compilation with -DNLA_MAX_THREADS=N if desired.
Out-of-range requests are rejected. Actual pool size is recorded; a packed
small matrix (at most 32768 columns) or failed pool creation is clearly marked
as serial fallback, not misleadingly labelled parallel. A summary actual count
of zero means pool size varied between repeats; consult the per-solve TSV.

These are solver scaling results, not predictions of overall factoring time.
Synthetic matrices cannot capture every QS incidence/rank/cache peculiarity.
For a reliable choice, use several seeds/densities on the target hardware,
with no competing work. Large thread counts allocate one private row buffer
per solver worker: at the default size, roughly 1.4 MiB per worker in addition
to shared matrix/solver storage. Oversubscribed CPUs and NUMA placement can
strongly change results. The tool does not pin threads or select a production
thread limit automatically.
