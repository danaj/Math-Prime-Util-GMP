# Embedding the binary matrix solvers in C

The solvers find right-nullspace dependencies over GF(2): subsets of matrix
columns whose XOR is zero. They work independently of SIQS and return
dependencies, not integer factors. A factoring host supplies its own relation
matrix and uses the dependencies in its square-root step.

Copy `lanczos.c`, `lanczos.h`, and `ptypes.h` into your project. For parallel
Lanczos, also copy `planczos_inc.c`. Compile with `STANDALONE`; no GMP, Perl,
SIQS adapters, prime caches, startup, or shutdown calls are needed. Keep the
copyright notices and consult the repository's `LICENSE` for redistribution
terms.

## Simple solver calls

All three solvers use the same column array and result format. Given
`unsigned long nrows, ncols`, a `la_col_t *cols`, and caller-owned
`uint64_t mask`, the basic calls are:

```c
/* Exact dense elimination, usually the best choice for small matrices. */
uint64_t *deps = la_dense_nullspace(nrows, ncols, cols, &mask);
```

```c
/* Serial block Lanczos. Zero means no externally packed dense rows.
 * Seed is a caller-owned value; the final zero requests quiet output. */
uint64_t *deps = la_block_lanczos(nrows, 0, ncols, cols,
                                 UINT64_C(0x83d2e5b79a4c610f), &mask, 0);
```

```c
/* Requires a PSIQS build. Request four total threads, including the caller.
 * retain_all_rows = 0, verbose = 0, nthreads = 4. */
uint64_t *deps = la_block_lanczos_threaded(nrows, 0, ncols, cols,
                                          UINT64_C(0x83d2e5b79a4c610f),
                                          &mask, 0, 0, 4);
```

These are alternatives, not three allocations to put into one variable
without freeing the earlier result. Check `deps != NULL && mask != 0`, use
the dependencies, then `free(deps)`. The solvers do not free or modify the
supplied columns. The optional matrix reducer described below is different:
it reorders columns and frees discarded column data.

Lanczos takes one `uint64_t` seed. Each call uses standard SplitMix64 to
initialize its random blocks, advancing one local state across internal
retries. Any seed is valid, including zero. The examples use an arbitrary
fixed seed for reproducibility; a host can supply different seeds for fresh
dependency samples. There is no shared RNG state, and serial and threaded
calls with the same matrix, seed, and row-retention setting agree.

The host chooses the solver; there is no automatic public dispatcher. SIQS
uses dense elimination through `LA_DENSE_CROSSOVER_COLS` (currently 1536
columns after reduction), then Lanczos. This is a measured SIQS crossover,
not a correctness limit or a universal optimum. Calling the threaded entry
point does not select dense elimination automatically.

## Matrix representation

Include `lanczos.h`. Its column type is:

```c
typedef struct {
  uint32_t *data;
  uint32_t weight;
  uint32_t orig;
} la_col_t;
```

For the ordinary format, each column's `data[0..weight-1]` lists the
zero-based row indices where that column has a one. Reduce integer entries
modulo two first: even entries disappear, and each remaining row occurs
once per column. Sorting row indices is useful for a canonical representation
but is not required by the solvers. Every index must be less than `nrows`.
An all-zero column has weight zero and may have `data = NULL`.
`mask` must point to valid output storage; a nonempty matrix needs `ncols`
valid column descriptors and the data described by each descriptor.

`orig` is a host column identifier, normally its position in the original
matrix. The solvers return results in the supplied array's order, not `orig`
order. The reducer preserves `orig` so a host can map retained columns back
to its original relations.

Column data may be static, on the stack, or dynamically allocated when
calling a solver directly. Keep it alive and unchanged until the call
returns. If using `la_reduce_matrix`, each column's data must instead be a
separate allocation that can be passed to `free`.

The reducer and Lanczos reject row or column dimensions above `UINT32_MAX`;
row indices and column weights are 32-bit. Available address space and
memory can impose much smaller practical limits. The internal 16-bit packed
column representation is selected only when the matrix fits it; it does
not impose a 65536-column limit on the public Lanczos interface.

## Reading and freeing dependencies

A successful result is an allocated array of `ncols` `uint64_t` values.
Bit `d` of `deps[c]` says whether column `c` belongs to dependency `d`.
The set bits in `mask` identify the returned dependencies; inspect bits
0 through 63 rather than assuming a fixed count or a particular basis.

```c
unsigned int d;
unsigned long c;

for (d = 0; d < 64; d++) {
  uint64_t bit = (uint64_t)1 << d;
  if (!(mask & bit))
    continue;
  for (c = 0; c < ncols; c++)
    if (la_get_null_entry(deps, c, d)) {
      /* Column c participates in dependency d. */
    }
}
free(deps);
```

Dense elimination returns up to 64 independent dependencies. Lanczos returns
a sample of up to 64, not a guarantee of the complete nullspace. Both verify
returned dependencies against the supplied matrix. Different seeds or solver
choices can produce different valid bases. A valid dependency does not
necessarily produce a nontrivial factor in a factoring application.

`NULL` with `mask = 0` means no dependencies were returned, not proof that
the matrix has no nullspace. Dense elimination can also return this on an
allocation/size failure; Lanczos can exhaust its internal attempts. A host
can try a different solver or seed, or collect more relations. Bound
host retries rather than retrying the same unsuccessful work indefinitely.

Use ordinary `free` for the returned array, with a compatible C runtime.
The column array and its data remain the host's responsibility.

## Complete small example

Save this as `my-nullspace.c`. It verifies and prints dependencies for a
three-row, four-column matrix. Its columns 0, 1, and 2 XOR to zero; columns
0 and 3 are identical, giving another dependency. No reducer is used, so
the stack-allocated column data needs no heap allocation or cleanup.

```c
/* my-nullspace.c */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lanczos.h"

int main(int argc, char **argv) {
  uint32_t c0[] = {0, 1}, c1[] = {1, 2};
  uint32_t c2[] = {0, 2}, c3[] = {0, 1};
  la_col_t cols[] = {{c0, 2, 0}, {c1, 2, 1},
                     {c2, 2, 2}, {c3, 2, 3}};
  unsigned long nrows = 3, ncols = 4, c, i;
  unsigned int d, selected;
  unsigned char parity[3];
  uint64_t mask = 0, *deps;

  if (argc != 2) {
    fprintf(stderr, "usage: my-nullspace dense|lanczos|wide|threaded\n");
    return 2;
  }
  if (strcmp(argv[1], "dense") == 0) {
    deps = la_dense_nullspace(nrows, ncols, cols, &mask);
  } else if (strcmp(argv[1], "lanczos") == 0) {
    deps = la_block_lanczos(nrows, 0, ncols, cols,
                            UINT64_C(0x83d2e5b79a4c610f), &mask, 0);
  } else if (strcmp(argv[1], "wide") == 0) {
    deps = la_block_lanczos_wide(nrows, 0, ncols, cols,
                                 UINT64_C(0x83d2e5b79a4c610f), &mask, 0);
  } else if (strcmp(argv[1], "threaded") == 0) {
#ifdef PSIQS
    deps = la_block_lanczos_threaded(nrows, 0, ncols, cols,
                                    UINT64_C(0x83d2e5b79a4c610f), &mask, 0, 0, 4);
#else
    fprintf(stderr, "rebuild with -DPSIQS -pthread for this entry point\n");
    return 2;
#endif
  } else {
    fprintf(stderr, "unknown solver: %s\n", argv[1]);
    return 2;
  }
  if (deps == NULL || mask == 0) {
    fprintf(stderr, "no dependencies returned\n");
    free(deps);
    return 1;
  }

  for (d = 0; d < 64; d++) {
    uint64_t bit = (uint64_t)1 << d;
    if (!(mask & bit))
      continue;
    memset(parity, 0, sizeof(parity));
    selected = 0;
    for (c = 0; c < ncols; c++) {
      if (deps[c] & bit) {
        selected++;
        for (i = 0; i < cols[c].weight; i++)
          parity[cols[c].data[i]] ^= 1;
      }
    }
    if (selected == 0 || parity[0] || parity[1] || parity[2]) {
      fprintf(stderr, "invalid dependency\n");
      free(deps);
      return 3;
    }
    printf("dependency %u:", d);
    for (c = 0; c < ncols; c++)
      if (deps[c] & bit)
        printf(" %lu", c);
    putchar('\n');
  }
  free(deps);
  return 0;
}
```

Build from the repository root, or adapt the paths to your copied files:

```sh
cc -O3 -DSTANDALONE -I. -o my-nullspace my-nullspace.c lanczos.c
./my-nullspace dense
./my-nullspace lanczos
./my-nullspace wide
```

For parallel support:

```sh
cc -O3 -DSTANDALONE -DPSIQS -pthread -I. \
  -o my-nullspace my-nullspace.c lanczos.c
./my-nullspace threaded
```

The tiny example exercises the threaded entry point's serial fallback,
not actual parallel execution. Real worker creation currently requires more
than 32768 columns. Neither build needs `-lgmp` or `-lm`. Do not compile
`planczos_inc.c` separately: `lanczos.c` includes it. Use matching build
defines for the caller and solver; compiling with `PSIQS` preserves the
serial entry points as well.

## Optional matrix reduction

For a factoring matrix, the host can peel singleton rows and trim heavy
excess columns before solving:

```c
unsigned long original_ncols = ncols;

/* Here every cols[c].data must be separately malloc-allocated. */
la_reduce_matrix(&nrows, &ncols, cols, 0);
deps = la_dense_nullspace(nrows, ncols, cols, &mask);
/* Or call Lanczos on the reduced matrix. */
```

Reduction is optional and tailored to factoring: it can discard valid
dependencies when trimming excess columns. Skip it if you need to preserve
the nullspace of a general matrix. It sorts and compacts `cols` in place,
updates `ncols`, and frees discarded columns' `data`. Surviving entries
occupy `cols[0..ncols-1]`; trailing `data` pointers are set to NULL.
`nrows` remains the original row-index bound, not the number of active rows;
the surviving row indices are not renumbered. Do not lower it to the printed
active-row count. The solvers perform their own internal row compaction.

Initialize distinct `orig` values before reduction. To expand successful
dependencies into the original column order:

```c
uint64_t *original_deps;
unsigned long c;

original_deps = (uint64_t *)calloc(original_ncols, sizeof(*original_deps));
if (original_deps == NULL) {
  /* Handle allocation failure before using original_deps. */
} else {
  for (c = 0; c < ncols; c++)
    original_deps[cols[c].orig] = deps[c];
  /* Use original_deps, then free it. Discarded columns stay zero. */
  free(original_deps);
}
free(deps);
for (c = 0; c < original_ncols; c++)
  free(cols[c].data);
/* Also free cols if the descriptor array itself was allocated. */
```

This mapping assumes `orig` was initialized to the original column index
and the solve succeeded. Do not retain owning aliases to the old `data`
pointers: reduction can free or move them. Do not pass stack/static column
data, slices of one shared allocation, or pointers requiring another
deallocator to the reducer. For such input, make independent heap copies
first, or skip reduction.

## Packed input rows and the wider Lanczos sample

The easiest interface uses `dense_rows = 0` and lists every one explicitly.
For an existing host format with a packed dense head, Lanczos also accepts
`dense_rows > 0`: after each column's `weight` sparse indices, append
`ceil(dense_rows / 32)` `uint32_t` words. Bit `r % 32` of word `r / 32`
encodes row `r`, for `0 <= r < dense_rows <= nrows`. Sparse indices should
then refer only to rows at or above `dense_rows`, without duplicating the
packed head. Clear unused high bits of the last packed word.

This suffix format is supported by serial, wide, and threaded Lanczos.
It is not supported by dense elimination or `la_reduce_matrix`; expand
the head into explicit row indices before using either of those functions.
It is separate from the solver's automatic internal packing optimization.

For another dependency sample, use:

```c
deps = la_block_lanczos_wide(nrows, 0, ncols, cols,
                            UINT64_C(0x83d2e5b79a4c610f), &mask, 0);
/* Or use retain_all_rows = 1 in la_block_lanczos_threaded. */
```

The ordinary solver can postpone 48 frequent rows until dependency
extraction when it packs a matrix. The wide variant retains those rows
during iteration instead. It uses the same 64-bit result format and is
not a larger-index or larger-block-width interface. This is an optional
retry choice, not a guarantee of more dependencies. Free an earlier result
before replacing it.

## Parallel behavior and errors

Thread counts are caller-selected; use 1 for serial execution. The count
includes the calling thread, which participates in the kernels. A pool
with N total threads creates N-1 pthreads, lives through the solve and its
internal retries, and joins every worker before returning.

The current implementation falls back to serial for a request below two
threads, a packed matrix, at most 32768 columns, or zero input rows. The
default `NLA_MAX_THREADS` is 32; larger requests are clamped to that cap.
Override it when compiling `lanczos.c` with `-DNLA_MAX_THREADS=N` if needed.
This is independent of SIQS worker limits. With verbosity above zero,
`Lanczos using N threads` reports an actually created pool. If optional
pool allocation, initialization, or thread creation fails, all started
workers are joined and the solve falls back to serial.

Each pool thread needs roughly `8 * nrows` bytes for a private row buffer,
plus smaller kernel scratch. More threads are not necessarily faster;
benchmark representative matrices on the target machine. Threads are not
pinned, and the solver does not detect performance versus efficiency cores,
NUMA topology, or a machine-wide concurrency budget.

Separate standalone solves may run concurrently, with private scratch,
seed state, result storage, and pools. Immutable column data can be shared;
each call needs its own mask and result. Never run the destructive reducer
concurrently with any use of that matrix. Concurrent verbose output can
interleave, and concurrent pools can oversubscribe CPU and memory resources.

Verbosity is per call: 0 quiet, 1 summaries, 2 the progress level, and 3
detailed diagnostics. The current solver's level 2 has the same messages
as level 1; it does not provide periodic progress or heartbeat output.
Negative values are clamped to zero. Dense elimination has no verbosity
argument and does not print normal diagnostics.

The standalone `croak` in `ptypes.h` prints a fatal diagnostic to stderr
and exits with status 3. Invalid matrix dimensions/indices and most
reducer/Lanczos allocation failures use that path; not every failure is
a recoverable NULL return. There is no host error callback or cancellation
API. To change fatal error policy, adapt the standalone `croak` definition;
do not unwind or longjmp past active solver workers.

## Validation and scaling

The matrix health checks cover reduction, original-column mapping,
dense and Lanczos nullspaces, packed input rows, and serial/threaded agreement:

```sh
perl tools/siqs-check.pl --suite matrix --threaded
perl tools/siqs-check.pl --suite matrix --threaded --extended
```

Those repository checks include SIQS fixtures and need GMP, unlike the
embedded solver itself. The extended suite includes an actual parallel
solve above the worker threshold. For scaling measurements, see
[README-lanczos-bench.txt](README-lanczos-bench.txt); its synthetic matrices
help compare thread counts without collecting a large factorization first.
