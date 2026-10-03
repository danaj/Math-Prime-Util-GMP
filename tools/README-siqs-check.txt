SIQS standalone health checks
=============================

These optional developer checks exercise the current SIQS implementation.
They do not build/load the Perl XS module, run performance benchmarks, or
attempt complete large factorizations. They are not part of ordinary
"make test". Production sources and benchmark builds are not modified.

After "perl Makefile.PL":

  make check-siqs
  make test-siqs                          # alias
  make check-siqs SIQS_CHECK_ARGS='--extended --verbose'
  make check-siqs SIQS_CHECK_ARGS='--sanitize'

The Make target uses only a core-Perl (5.10+) build wrapper, the C compiler,
and GMP.
The wrapper compiles one checker into a temporary directory and removes it
afterward. It does not require Math::Prime::Util or Math::Prime::Util::GMP to
be installed or built. Use --keep to retain the binary.

Direct build without Perl or a generated Makefile, from the repository root:

  cc -O3 -DSTANDALONE -o /tmp/siqs-check tools/siqs-check.c \
    lanczos.c prime_iterator.c squfof126.c pbrent63.c -lgmp -lm
  /tmp/siqs-check --suite sieve
  /tmp/siqs-check --extended --verbose

Add -march=native only for a binary intended for the local machine.

Suites and options
------------------

One executable runs the named suites (currently just "sieve"). The default
"all" selection runs every registered suite. --list describes available
suites; --suite selects one. --extended adds more fixtures, policy endpoints,
and polynomials. --verbose reports each real-polynomial setup. New suites can
be registered without creating separate executables or a large framework.

  perl tools/siqs-check.pl --list
  perl tools/siqs-check.pl --suite sieve --extended
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
