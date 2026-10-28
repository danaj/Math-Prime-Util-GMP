# MPU-SIQS Implementation details

## Implementation notes

The implementation combines several pieces of engineering rather than one
isolated optimization:

- Self-initializing polynomial families and size-dependent multiplier and
  polynomial selection.
- A cache-blocked byte sieve, inexpensive root updates, and specialized
  one-hit and two-hit kernels for large primes.
- Partial relations with one or two large primes, combined into full relations
  before constructing the matrix.  Collection readiness is based on useful
  relations from the graph, not just the raw number of partials received.
- A bounded incremental dependency solver ([BIDS](#bids)) for small factor bases.
  It can use a dependency as soon as it arrives, without collecting a complete
  batch first.  Larger problems use dense elimination or block Lanczos.
- Asynchronous parallel collection with private worker state, plus a
  pthread block Lanczos implementation for the sparse matrix stage.
- Measured parameter bands and compact storage for relations and matrices.

The tuning balances factor-base size, sieve interval, and large-prime
bounds.  A setting that produces more partials need not finish sooner,
as higher bounds can increase cofactor-splitting costs, storage, graph
processing, and eventual matrix work.
Similarly, saving memory can be worthwhile even when relation
collection takes a little longer.

Profiles and controlled comparisons are needed to separate the effects
of parameter choices, sieving, relation handling, and linear algebra.


## Goals

- Reliable factoring throughout the supported range.  MPU-SIQS should
  strive to split every supported composite, including unusual small
  inputs.  Targeted recovery profiles handle difficult cases without
  making ordinary inputs pay for broader default settings.

- Let the caller choose the factoring strategy.  After modest trial
  division and basic input checks, we use SIQS.  The caller decides
  whether to run SQUFOF, p−1, ECM, or other methods first.

- Find multiple factors when practical.  We use the collected relations
  to split additional composite cofactors where possible.  A call can
  therefore produce more than two factors, although complete
  factorization is not guaranteed.

- Make SIQS fast across its supported range.  Small inputs deserve
  attention as do large ones.  Our goal is to make SIQS as fast as
  practical, while the caller decides when to use another algorithm.

- Make parallelism useful at moderate sizes as well as large.  Aim
  for low parallel overhead and good scaling, while retaining serial
  execution where additional workers would cost more than they help.
  This includes our use of parallel block Lanczos as soon as it is faster.

- Keep memory use reasonable.  Compact storage and modest per-worker
  memory matter, especially when running many workers.

- Use portable C that compilers can optimize.  Structure data and loops
  to support compiler-generated vectorization without maintaining
  hand-written assembly or SIMD implementations.


## Sieving

Our byte sieve applies dense primes in cache-sized blocks, then applies
sparse primes across the full interval.  The kernels group primes by
the number of hits each root can have in the range being sieved, with
specialized one-hit and two-hit kernels for sparse updates.  Careful
organization of these loops improves cache locality and gives the
compiler opportunities for vectorization.

I have worked on optimizing segmented prime sieves for over 10 years, so my
approach to the problem was from that perspective.  Duplicating work, as
long as it keeps the working set in cache, can be *much* faster than
monolithic work even with a higher raw instruction count.  I also started
with modern systems -- the initial development was done on an Apple M1 Pro
(ARM64) with modern compilers.

We implemented and tested a bucket sieve, but the direct kernels performed
better in our measurements.  Buckets avoid repeatedly scanning primes that
have no hits in a block, but introduce work to build, store, and traverse
hit lists.  At our tested sizes, the direct approach was the better tradeoff.
A different bucket structure might reverse that, but we didn't find one.

This is an important design difference from YAFU.  Its x86 implementation
has highly optimized AVX-512 sieve kernels.
When trying to explain why our end-to-end factoring was consistently faster
at all measured sizes, we captured profile data for a 100-digit input.
That showed that YAFU spent most of its time in root updates and bucket
work rather than its fast dense sieving.
With each program’s default parameters, the difference in combined sieving
and root-preparation cost accounted for most of the observed runtime gap.

## Policy bands

Implementations often select parameters using tables indexed by digit
or bit length, or calculate them from formulas.  We chose to select
parameters using a table of policy bands, each covering a range of input
bit lengths.  The bit size is measured after trial division and before
applying the multiplier.  Some policy values are fixed within a band,
others vary linearly with bit length.  Factor-base sizes, sieve intervals,
and other working parameters are derived from these values.
This keeps tuning centralized and expresses broad trends with
relatively few coefficients.

The tool `tools/parse-siqs-bands.pl` will read the source and parse the
bands into an easily readable table of the important parameter settings.

As SIQS with tiny inputs has a number of unique challenges, there are
backup policy bands for inputs up to 64 bits.  These are very rarely
used but essential to allow simultaneous goals of always factoring and
maintaining the best average performance.  They cost nothing for inputs
that don't need them.

## Choosers

The SIQS algorithm uses 'A families' each containing $2^{q-1}$ polynomials.
The A values are constructed using 'q' primes from the factor base, with
a product A that is close to a target value.  This is well documented in
the literature, e.g. Carrier and Wagstaff (2002), Kechlibar (2005),
Crandall and Pomerance (2010), and Milan (2010) among many others.

This still leaves the exact method used up to the implementation, and while
we started with a simple yet effective method similar to
[JML](https://github.com/TilmanNeumann/java-math-library), we have
replaced it with new systems of our own design.

Much of the research literature stresses constructing A such that it is
extremely close to the target, e.g. a distance of 0.1%.  In one experiment
the legacy chooser was measured at medians of 3-5% with 95th percentile
values sometimes exceeding 50%, yet it performed quite well.  The main
problems we encountered were not distance from the target, but trivial
dependencies resulting from families sharing many primes.  This was much
less of a problem at higher q, but it was critical for very small inputs
using q <= 4.

An experiment using 185-218 bit inputs (q=9 for our policies) used an
offline oracle to supply A values within 0.01% of the target.  Even with
selection cost excluded, the tighter values had almost no effect on
end-to-end factoring times.
My reading led me to expect a larger benefit from tighter A selection.
Several parts of the current implementation grew out of this work:

- BIDS, the Bounded Incremental Dependency Solver, described in
  [the next section](#bids).  This takes care of a number of problems
  unique to very small inputs, while improving overall performance.

- Deterministic choosers for q=1, q=2, and q=3.  These have excellent
  median A distances and deal with some of the unique issues of tiny
  inputs.

- SAR, Staged Adaptive Random chooser, reduces median A distance
  and improves the outliers compared to the legacy chooser for ranges
  we measured.  Still quite fast for smaller inputs.

- SAC, Staged Adaptive Cubic chooser.
  Cubic generates a sequence of unique prime prefixes.  Any two prefixes
  share at most three primes within a batch. The $F_{19}$ extension
  provides 19 batches, with at most four shared prefix primes across
  batches.  In all cases, the final prime is fitted from the eligible
  factor base, and the prefix-overlap guarantee does not apply to it.
  It falls back to SAR if the cycle is exhausted or a bounded search
  finds no usable A.

While Cubic was designed and implemented before reading Carrier and
Wagstaff, the concept can be seen as vaguely similar to their idea of
using a deterministic lexicographic map (section 3.1) for the prefix.
The latter has a longer cycle but gives no prime-duplication guarantee,
which is the problem we were trying to solve with Cubic.

The SAR chooser works at all q >= 2, and is an excellent all-around
choice.  It could be tuned to give tighter distances if desired (at
a performance cost for smaller inputs).
SAC (Cubic) works for all 5 <= q <= 16.  After testing, there was no
strong evidence to use Cubic over SAR as they are currently implemented.

The current implementation uses deterministic for `q <= 2`, and SAR
for `q >= 3`.  Cubic is not used by default (it can be enabled with a define).
The staged approach combines widening the accepted A error as well
as widening the factor base window after repeated misses.

## BIDS

Typical operation of a quadratic sieve involves collecting roughly FB
relations plus a surplus, possibly using matrix readiness
tests near the end to decide when to run the GF(2) linear algebra solver.
Solvers include dense solvers (e.g. Gaussian Elimination) and Block Lanczos.
If the resulting dependencies yield no useful factors, collection resumes
to gather more relations, followed by another solve step.
Careful attention can make this work well for large inputs.
With small inputs, an apparently sufficient collection of relations
can still yield only trivial factors.  Investigating this problem,
including the chooser work, led to an incremental solver.

BIDS keeps an elimination basis as full relations arrive.  Each relation
either extends the basis or produces a dependency, which is immediately
tested for factors.  If it produces only trivial factors, collection
continues using the existing basis.  This avoids repeated batch solves.
It can also complete a factorization before the conventional relation
target is reached, saving sieving work as well as solver work.
A relation budget bounds each incremental attempt.

BIDS also works with parallel relation collection.  The workers collect
relations while the coordinator maintains a single elimination basis and
tests new dependencies. Once the factorization is complete, the coordinator
requests a stop -- workers may still have some work in flight.

Since this is a dense solver, both memory and time growth make it
appropriate for smaller factor bases.  We continue to use block
Lanczos for large sparse matrices.  What came as a surprise was
how well this system, meant to fix problems on 64-bit inputs, works
even at larger input sizes.  MPU-SIQS currently uses BIDS when the
factor base is less than 4096 primes, or about 196-bit (59 digit)
inputs with current parameters.

To the best of my knowledge, no other current SIQS implementation uses
something like this.  I later found the concept of an incremental dense
solver for quadratic sieve described in
[*Compact incremental Gaussian elimination over Z/2Z*](https://infoscience.epfl.ch/bitstreams/6551ea6d-64fa-463a-8dd9-6542168ba782/download), Lenstra and Manasse (1988).
We developed BIDS independently, before finding this paper.

## How LP and R are selected

Policy bands use bits(N), where N is the input remaining after
trial division, before the multiplier k is applied.  Each band
contains `K_L` and `K_R` curves, where `P` is the actual largest
factor-base prime.  `LP` limits each individual large prime,
while `R` limits the remaining cofactor after factor-base division.

`LP = floor(K_L × P)`, capped at $2^{38}−1$.  With
`SIQS_R_DERIVED_FROM_LP`, `R` is the final capped `LP` in smooth/1LP mode, or
`LP × P` in 2LP mode.  An explicit `K_R` instead sets `R = floor(K_R × P)` in
smooth/1LP mode, or `floor(K_R × P²)` in 2LP mode.  `R` saturates at 2^64−1.
Smooth-only profiles use `K_L = 1`.

`SIQS_POLICY_LINEAR(base, step, origin)` permits extrapolation from
its origin.  `SIQS_POLICY_RAMP(beg, end, first_bits, last_bits)` stores
endpoint values and interpolates linearly between them; extrapolation
is an error.  The policy table and its comments in [siqs.c](siqs.c)
describe the other parameters.

The sieve-hit exponent sets a separate coarse residual allowance,
raised to at least `2 × R`.  The fine filter uses the final `R`.

The lower bands have been tuned with these separate K curves.  Settings
above the end-to-end tested range remain provisional.  Changing factor-base
size also changes P, and therefore both LP and R.
The four tuning knobs -- factor-base size, sieve interval, `K_L`,
and `K_R` -- must therefore be considered together, even though their
policy curves can be edited separately.

## Other

Before starting SIQS, trial division tests primes below 200 for inputs
up to 40 bits, and below five times the original input bit length for
larger inputs, capped at 5,000.
The `trial_start` argument lets callers skip redundant trial division.
Tests on composites below 2^32 also succeeded with a trial-division
bound of 10.  The balanced-semiprime benchmarks in the README and graphs
have factors above these bounds.

## Build and API notes

Run the commands in this section from the siqs/ directory.

`make siqs` builds `msiqs` with pthread support when a compile/link probe
succeeds, otherwise serial-only.  `-t N`, `--threads N`, and `-threads N`
are equivalent ways to request workers in a threaded build.
`make siqs-serial` always rebuilds `msiqs` without
threading.  Remove `msiqs` before rebuilding to return to
autodetection or after changing compiler settings.

The standalone build normally uses `-O3 -march=native -falign-loops=32`;
compiler flags can be overridden with `CFLAGS`.  The 32-byte loop alignment
avoids placement-sensitive sieve slowdowns seen on x86; retain it with custom
flags when the compiler supports it.  The Makefile also accepts `CC`, `CPPFLAGS`,
`LDFLAGS`, and `LDLIBS`.
There is a build helper using Perl, but the executable itself has nothing
to do with Perl.

```sh
cc -O3 -march=native -falign-loops=32 -DSTANDALONE -DPSIQS -pthread -o msiqs mpu-siqs.c siqs.c lanczos.c pr32_iter.c ../squfof126.c -lgmp -lm
```

Remove `-DPSIQS -pthread` to disable threading.
The default is one worker; use -t N to request more.
The standalone program does not autodetect the number of threads, as this
is problematic on machines with SMT (x86) and on newer processors with
heterogeneous cores ("performance" and "efficiency" cores).

Thread count is a resource choice, not a
promise of linear scaling, particularly during linear algebra.

For use in another C project, see the
[SIQS embedding guide](tools/README-siqs-embedding.txt) and the
[matrix-solver embedding guide](tools/README-lanczos-embedding.md).
The low-level SIQS interface returns an array of factors whose product
equals the input.  Some factors may remain composite.
