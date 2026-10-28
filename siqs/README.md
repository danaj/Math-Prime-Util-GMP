# SIQS in Math::Prime::Util::GMP

MPU-SIQS is a portable C implementation of the Self-Initializing Quadratic
Sieve for integer factoring, using GMP for large-integer arithmetic.  It is
the default quadratic sieve in Math::Prime::Util::GMP 0.55.
MPU-SIQS was written in late 2026 by Dana Jacobsen, using OpenAI Codex
to write parts of the code, tests, and tools.

This README focuses on building and using the standalone SIQS tool.

Both SIQS and the block Lanczos solver can be embedded in another C
project in addition to their use in the included standalone program.
Optional pthreads accelerate both relation collection and block Lanczos
linear algebra.  The sieve and matrix solvers use portable C without
hand-written assembly or explicit SIMD intrinsics.

## Installation and use

Building requires a C compiler and GMP.  You do *not* need Perl to build
the standalone tool, which can be compiled with a single CC line.

### Standalone program

From the repository root (no module configuration or installation needed):

```sh
cd siqs
make
./msiqs 22095311209999409685885162322219
```

To request eight workers and periodic progress reports:

```sh
./msiqs -t 8 -vv 22095311209999409685885162322219
```

The Makefile detects pthread support automatically.  The default is one
factoring thread, and very small inputs may stay serial even when more
are requested.  `-v` prints setup and final summaries, `-vv` adds progress.
With no integers on the command line, `msiqs` reads whitespace or newline
separated inputs from stdin.
See [build and API notes](IMPLEMENTATION.md#build-and-api-notes) for
compiler options and serial-only builds.

To save relation collection and resume it later, use an optional journal:

```sh
./msiqs -t 8 --journal job.journal 22095311209999409685885162322219
./msiqs -t 8 --journal job.journal
python3 tools/journal-info.py --follow --eta job.journal
```

The [journal guide](JOURNAL.md) describes restart behavior, inspection and
the file format. Journalling is off unless requested.

### Embedding in C

The [SIQS embedding guide](tools/README-siqs-embedding.txt) gives build
instructions and host examples.
The [matrix solvers](tools/README-lanczos-embedding.md) can also be
embedded independently, without GMP or Perl, and are not tied to this
SIQS implementation.

## SIQS Performance at a glance

All settings left at default, no tuning for machine or input.

Representative timing results on an Apple M1 Pro, using multiple
balanced semiprimes and an October 2026 SIQS build.

| Decimal digits | Median time <br> 1 thread | Median time <br> 8 threads | Peak RSS <br> 1 thread |
| ---: | ---: | ---: | ---: |
| 30 | 0.95 ms | 1.08ms | 2.8 MiB |
| 40 | 7.2 ms | 3.0ms | 3.3 MiB |
| 60 | 0.51 s | 0.089s | 11 MiB |
| 80 | 39 s | 5.6s | 66 MiB |
| 100 | 51 min | 7 min | 410 MiB |

With 8 threads on a 2021 MacBook M1 Pro, RSA-100 completed in 4m 58.7s,
RSA-110 in 51m 36.5s.

A 2026 MacBook M5 Pro with 18 threads completed RSA-100 in 2m 16.2s.

Timing results on a 48-core AMD Turin (EPYC 9B45, Zen5), 48 threads:

| Input    | digits | bits | time | CPU time |
| :---:    | ---: | ---: | ---: | ---: |
| RSA-100  | 100  | 330 |      42.7 s | |
| RSA-110  | 110  | 364 |   7m 29s | 6 core hrs |
| RSA-120  | 120  | 397 |  63m 22s | 51 core hrs |
| RSA-130  | 130  | 430 | 9h 28m | 454 core hrs |

A 96-core Xeon 8581C completed RSA-130 in 5h 44m on October 24 2026.

More comparisons can be seen in the [benchmarks document](BENCHMARKS.md)

### Supported range

Parameters cover inputs through 512 bits. Complete factoring runs have
been verified through 430 bits, while the larger range has had component
testing.  The design expectation of MPU-SIQS is to factor all composite
inputs within its size range.  Composites with multiple factors are
usually fully factored, but this is not guaranteed.

MPU-SIQS implements the Self-Initializing Quadratic Sieve.  It does not
implement ECM or NFS, which can be more practical choices for some inputs.
While it works well from 1 to at least 430 bits, typically SIQS is one part
of a generic factoring package, as in MPU's factor routine, or in
[YAFU](http://sourceforge.net/projects/yafu/).

## Further reading

- [implementation details](IMPLEMENTATION.md)
- [evolution and history](EVOLUTION.md)
- [benchmark comparisons and methodology](BENCHMARKS.md)
- [benchmark graphs](graphs/).
