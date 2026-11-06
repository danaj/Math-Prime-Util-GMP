# SIQS in Math::Prime::Util::GMP

This snapshot comes from msiqs revision `755054b302e9`.

The module contains a frozen serial edition of [msiqs](https://github.com/danaj/msiqs),
a portable Self-Initializing Quadratic Sieve implementation using GMP.
It began as part of this module in 2026, written by Dana Jacobsen with
assistance from OpenAI Codex. The standalone project now develops independently.

Install the module as described in [README](README). Use `factor($n)` for the
general factoring strategy, or `qs_factor($n)` to call SIQS directly:

```perl
use Math::Prime::Util::GMP qw(factor qs_factor);
my @factors = factor($n);
my @partition = qs_factor($n);
```

`qs_factor` returns a multiplicative partition whose elements are not
necessarily prime.  `factor` continues factoring composites and reports an
explicit error if it exhausts its methods.

## Performance

SIQS is effective for balanced composites and can find several factors in
one call. Work depends on input size, factor shape, and the machine. The
serial implementation supports inputs through 512 bits; that arithmetic
limit is not a claim that factoring at the upper end is inexpensive.

Representative October 2026 measurements on an Apple M1 Pro, across
multiple balanced semiprimes, were approximately:

| Decimal digits | Median single-thread time |
| ---: | ---: |
| 30 | 0.95 ms |
| 40 | 7.2 ms |
| 60 | 0.51 s |
| 80 | 39 s |
| 100 | 51 min |

These historical measurements describe the tested builds, rather than a
new measurement of every later snapshot. The full [benchmark methodology and results](https://github.com/danaj/msiqs/blob/master/docs/BENCHMARKS.md)
and [evolution study](https://github.com/danaj/msiqs/blob/master/docs/EVOLUTION.md)
are maintained in msiqs.

## Maintenance

This edition has its own frozen sources and receives deliberate, manually
ported updates. It has no parallel workers, journal code, 3LP, or timing
instrumentation. No standalone SIQS program is built from this package.
The public Perl factoring interfaces and their behavior are unchanged.

See msiqs's [embedded package](https://github.com/danaj/msiqs/tree/755054b302e9edc55465f60dfacff7a8b14ffbe9/embedded)
for the source snapshot, native TinySIQS128, and embedding instructions.
TinySIQS128 is independent of GMP and is not used by this GMP module.
The embedded source license is included as [LICENSE-SIQS](LICENSE-SIQS);
the module's existing licensing remains in [LICENSE](LICENSE).
