/*
 * Standalone block-Lanczos scaling benchmark.  See README-lanczos-bench.txt.
 * Uses synthetic QS-like columns, not a factorization or a stored QS matrix.
 * Including the solver exposes its actual pool size and unchanged recurrence;
 * no instrumentation or benchmark hooks are needed in production sources.
 *
 * cc -O3 -march=native -DSTANDALONE -DPSIQS -pthread \
 *    -o lanczos-bench tools/lanczos-bench.c -lm
 * Copyright (c) 2026 Dana Jacobsen. See LICENSE for redistribution terms.
 */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#include <errno.h>
#include <fcntl.h>
#include <sys/resource.h>
#include <time.h>
#include <unistd.h>

#if !defined(STANDALONE) || !defined(PSIQS)
#error "Build lanczos-bench with -DSTANDALONE -DPSIQS -pthread"
#endif
/* The benchmark may explore beyond production's provisional 32-thread cap.
 * This does not change the solver cap in siqs/psiqs builds. */
#ifndef NLA_MAX_THREADS
#define NLA_MAX_THREADS 256U
#endif
#include "../lanczos.c"

static int bench_verbose;

typedef struct {
  uint32_t requested, effective, attempts, dependencies;
  int packed;
  double wall, cpu;
} bench_result_t;

static uint32_t bench_random(uint32_t *state) {
  uint32_t x = *state;
  x ^= x << 13; x ^= x >> 17; x ^= x << 5;
  *state = x;
  return x;
}

static int bench_compare_u32(const void *a, const void *b) {
  uint32_t x = *(const uint32_t *)a, y = *(const uint32_t *)b;
  return x < y ? -1 : x != y;
}
static int bench_compare_double(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return x < y ? -1 : x != y;
}
static double bench_median(double *values, uint32_t count) {
  qsort(values, count, sizeof(*values), bench_compare_double);
  return count % 2U ? values[count / 2U]
       : 0.5 * (values[count / 2U - 1U] + values[count / 2U]);
}
static double bench_wall(void) {
  struct timespec now;
  if (clock_gettime(CLOCK_MONOTONIC, &now) != 0)
    croak("lanczos-bench: clock_gettime failed: %s", strerror(errno));
  return (double)now.tv_sec + (double)now.tv_nsec / 1e9;
}
static double bench_cpu(void) {
  struct rusage usage;
  if (getrusage(RUSAGE_SELF, &usage) != 0)
    croak("lanczos-bench: getrusage failed: %s", strerror(errno));
  return (double)usage.ru_utime.tv_sec + (double)usage.ru_stime.tv_sec
       + ((double)usage.ru_utime.tv_usec + usage.ru_stime.tv_usec) / 1e6;
}

/* A reduced-core approximation: varying column weights; frequent small-prime
 * rows; a mostly uniform sparse tail with a minority of lower-row-biased hits.
 * Two distinct tail rows per column ensure that every tail row occurs at least
 * twice, without planting dependencies, duplicate columns, or an easy block
 * decomposition.  All entries (including the head) use ordinary sparse indices,
 * just as siqs_build_matrix does; input_dense_rows is zero. */
static la_col_t *bench_matrix(uint32_t rows, uint32_t count, uint32_t weight,
                              uint32_t head, uint32_t seed, size_t *entries) {
  uint32_t c, row, width = weight / 3U, tail = rows - head;
  uint32_t random_state = seed;
  uint32_t *marks = (uint32_t *)nla_calloc(rows, sizeof(*marks));
  uint32_t *counts = (uint32_t *)nla_calloc(rows, sizeof(*counts));
  la_col_t *cols = (la_col_t *)nla_calloc(count, sizeof(*cols));
  uint32_t min_count = UINT32_MAX, max_count = 0, inactive = 0;
  *entries = 0;
  for (c = 0; c < count; c++) {
    uint32_t i, n = 0, token = c + 1U;
    uint32_t target = weight - width + bench_random(&random_state) % (2U * width + 1U);
    uint32_t *data = (uint32_t *)nla_malloc(target + head, sizeof(*data));
    data[n++] = head + c % tail;
    data[n++] = head + (uint32_t)(((uint64_t)c + tail / 2U) % tail);
    marks[data[0]] = marks[data[1]] = token;
    while (n < target) {
      uint32_t offset = bench_random(&random_state) % tail;
      if ((bench_random(&random_state) & 3U) == 0)
        offset = (uint32_t)((uint64_t)offset * (bench_random(&random_state) % tail) / tail);
      row = head + offset;
      if (marks[row] != token) {
        marks[row] = token;
        data[n++] = row;
      }
    }
    for (row = 0; row < head; row++)
      if (bench_random(&random_state) % (2U + row / 4U) == 0)
        data[n++] = row;
    qsort(data, n, sizeof(*data), bench_compare_u32);
    cols[c].orig = c;
    cols[c].weight = n;
    cols[c].data = data;
    if ((size_t)n > (size_t)-1 - *entries)
      croak("lanczos-bench: matrix weight overflow");
    *entries += n;
    for (i = 0; i < n; i++) counts[data[i]]++;
  }
  for (row = 0; row < rows; row++) {
    if (counts[row] == 0) inactive++;
    if (counts[row] < min_count) min_count = counts[row];
    if (counts[row] > max_count) max_count = counts[row];
  }
  printf("Synthetic matrix: %u rows x %u columns, %lu entries (%.2f/column)\n",
         rows, count, (unsigned long)*entries, (double)*entries / count);
  printf("Row incidence: min %u, max %u, unused %u; head %u rows; seed %u\n",
         min_count, max_count, inactive, head, seed);
  /* This is the column ordering passed to Lanczos after la_reduce_matrix. */
  qsort(cols, count, sizeof(*cols), nla_compare_columns);
  free(counts);
  free(marks);
  return cols;
}

/* Mirrors only the public entry point's setup/retry/cleanup.  The matrix
 * conversion, pool, kernels, iteration, and dependency extraction are the
 * actual solver.  Private access lets us report serial fallback, rather than
 * silently labelling a capped or failed pool as the requested thread count. */
static uint64_t *bench_solve(uint32_t rows, uint32_t count, la_col_t *cols,
                             uint32_t threads, bench_result_t *run,
                             uint64_t *mask) {
  nla_matrix_t matrix;
  uint64_t *result = NULL;
  uint64_t rng_state = UINT64_C(0x83d2e5b79a4c610f);
  uint32_t attempt;
  double wall = bench_wall(), cpu = bench_cpu();
  memset(run, 0, sizeof(*run));
  run->requested = threads;
  run->effective = 1U;
  *mask = 0;
  nla_matrix_init(&matrix, rows, 0, count, cols, NLA_POST_ROWS, bench_verbose ? 3 : 0);
  run->packed = matrix.packed;
  matrix.pool = nla_pool_create(&matrix, threads);
  if (matrix.pool != NULL) run->effective = matrix.pool->nthreads;
  for (attempt = 0; attempt < NLA_MAX_ATTEMPTS; attempt++) {
    run->attempts++;
    result = nla_block_lanczos_once(&matrix, &rng_state, mask, bench_verbose ? 3 : 0);
    if (result != NULL && *mask != 0) break;
    free(result);
    result = NULL;
  }
  nla_pool_destroy(matrix.pool);
  nla_matrix_clear(&matrix);
  run->cpu = bench_cpu() - cpu;
  run->wall = bench_wall() - wall;
  return result;
}

/* Independent scalar check, outside the timer.  Check A*X=0, nonempty lanes,
 * and lane independence; also compare bit-for-bit with the serial solve. */
static uint32_t bench_verify(uint32_t rows, uint32_t count, const la_col_t *cols,
                             const uint64_t *result, uint64_t mask) {
  uint32_t c, i, dependencies = 0;
  uint64_t used = 0, basis[64] = {0};
  uint64_t *parity = (uint64_t *)nla_calloc(rows, sizeof(*parity));
  if (result == NULL || mask == 0)
    croak("lanczos-bench: solver returned no dependencies");
  for (c = 0; c < count; c++) {
    uint64_t bits = result[c];
    if ((bits & ~mask) != 0) croak("lanczos-bench: dependency mask mismatch");
    used |= bits;
    for (i = 0; i < cols[c].weight; i++) parity[cols[c].data[i]] ^= bits;
    for (i = 0; i < 64U && bits != 0; i++) {
      if (bits & (UINT64_C(1) << i)) {
        if (basis[i]) bits ^= basis[i];
        else { basis[i] = bits; break; }
      }
    }
  }
  for (i = 0; i < rows; i++)
    if (parity[i] != 0) croak("lanczos-bench: invalid dependency");
  if (used != mask) croak("lanczos-bench: empty dependency");
  for (i = 0; i < 64U; i++) {
    if ((basis[i] != 0) != ((mask & (UINT64_C(1) << i)) != 0))
      croak("lanczos-bench: dependent dependency lanes");
    if (basis[i]) dependencies++;
  }
  free(parity);
  return dependencies;
}

static uint32_t bench_number(const char *text) {
  char *end;
  unsigned long value;
  if (text[0] < '0' || text[0] > '9')
    croak("lanczos-bench: invalid number: %s", text);
  errno = 0;
  value = strtoul(text, &end, 10);
  if (errno != 0 || *end != '\0' || value > UINT32_MAX)
    croak("lanczos-bench: invalid number: %s", text);
  return (uint32_t)value;
}

static void bench_usage(void) {
  printf("Usage: lanczos-bench [options]\n"
    "  --max-threads N  Sweep 1,2,4,8,...,N (default 8; includes exact N).\n"
    "  --repeat N       Runs per count (default 1; alternate sweep direction).\n"
    "  --rows N         Matrix rows (default 183590, observed RSA-110 size).\n"
    "  --weight N       Mean sparse-tail entries per column (default 50).\n"
    "  --head N         Frequent small-prime rows (default 64).\n"
    "  --seed N         Nonzero matrix seed (default 1).\n"
    "  --output FILE    Write per-solve TSV; refuses to overwrite a file.\n"
    "  --verbose        Also print the solver's iteration diagnostics.\n"
    "Build thread cap: %u (production cap is unchanged).\n", (unsigned int)NLA_MAX_THREADS);
}

int main(int argc, char **argv) {
  uint32_t rows = 183590U, weight = 50U, head = 64U, seed = 1U;
  uint32_t maximum = 8U, repeats = 1U, counts[33], ncounts = 0, count;
  uint32_t repeat, point, c;
  size_t entries;
  const char *output = NULL;
  FILE *tsv = NULL;
  la_col_t *cols;
  bench_result_t *runs;
  uint64_t *reference = NULL, reference_mask = 0;
  double *times, serial_median;
  int arg;
  for (arg = 1; arg < argc; arg++) {
    const char *option = argv[arg];
    uint32_t value;
    if (strcmp(option, "--help") == 0) { bench_usage(); return 0; }
    if (strcmp(option, "--verbose") == 0) { bench_verbose = 1; continue; }
    if (++arg == argc) croak("lanczos-bench: missing value for %s", option);
    if (strcmp(option, "--output") == 0) { output = argv[arg]; continue; }
    value = bench_number(argv[arg]);
    if      (strcmp(option, "--max-threads") == 0) maximum = value;
    else if (strcmp(option, "--repeat") == 0) repeats = value;
    else if (strcmp(option, "--rows") == 0) rows = value;
    else if (strcmp(option, "--weight") == 0) weight = value;
    else if (strcmp(option, "--head") == 0) head = value;
    else if (strcmp(option, "--seed") == 0) seed = value;
    else croak("lanczos-bench: unknown option: %s", option);
  }
  if (maximum == 0 || maximum > NLA_MAX_THREADS)
    croak("lanczos-bench: max threads must be 1..%u", (unsigned int)NLA_MAX_THREADS);
  if (repeats == 0 || seed == 0 || weight < 2U || rows < 128U ||
      rows > UINT32_MAX - 64U || head >= rows ||
      (uint64_t)weight + weight / 3U > rows - head)
    croak("lanczos-bench: invalid size, weight, head, seed, or repeat count");
  count = rows + 64U;
  counts[ncounts++] = 1U;
  for (c = 2U; c <= maximum; c *= 2U) {
    counts[ncounts++] = c;
    if (c > maximum / 2U) break;
  }
  if (counts[ncounts - 1U] != maximum) counts[ncounts++] = maximum;
  if ((size_t)repeats > (size_t)-1 / ncounts)
    croak("lanczos-bench: repeat count exceeds allocation range");
  if (output != NULL) {
    int descriptor = open(output, O_WRONLY | O_CREAT | O_EXCL, 0666);
    if (descriptor < 0) croak("lanczos-bench: cannot create %s: %s", output, strerror(errno));
    tsv = fdopen(descriptor, "w");
    if (tsv == NULL) {
      close(descriptor);
      croak("lanczos-bench: fdopen failed: %s", strerror(errno));
    }
    fprintf(tsv, "rows\tcolumns\tentries\ttail_weight\thead\tseed\trepeat\trequested_threads\teffective_threads\tpacked\twall_s\tcpu_s\tdependencies\tattempts\tverified\n");
    fflush(tsv);
  }
  setvbuf(stdout, NULL, _IOLBF, 0);
  cols = bench_matrix(rows, count, weight, head, seed, &entries);
  printf("Full solves, excluding generation/independent verification; cap %u.\n",
         (unsigned int)NLA_MAX_THREADS);
  runs = (bench_result_t *)nla_calloc((size_t)ncounts * repeats, sizeof(*runs));
  times = (double *)nla_malloc(repeats, sizeof(*times));
  for (repeat = 0; repeat < repeats; repeat++) {
    for (point = 0; point < ncounts; point++) {
      uint32_t index = repeat % 2U ? ncounts - 1U - point : point;
      bench_result_t *run = runs + (size_t)index * repeats + repeat;
      uint64_t mask, *result;
      printf("BEGIN repeat %u/%u, threads %u\n", repeat + 1U, repeats, counts[index]);
      result = bench_solve(rows, count, cols, counts[index], run, &mask);
      run->dependencies = bench_verify(rows, count, cols, result, mask);
      if (reference == NULL) { reference = result; reference_mask = mask; }
      else {
        if (mask != reference_mask ||
            memcmp(result, reference, (size_t)count * sizeof(*result)) != 0)
          croak("lanczos-bench: solve differs from the serial reference");
        free(result);
      }
      printf("DONE  threads %u (actual %u): %.3f wall, %.3f CPU, %u dependencies, %u attempt(s)%s\n",
        run->requested, run->effective, run->wall, run->cpu,
        run->dependencies, run->attempts,
        run->effective != run->requested ? " [SERIAL FALLBACK]" : "");
      if (tsv != NULL) {
        fprintf(tsv, "%u\t%u\t%lu\t%u\t%u\t%u\t%u\t%u\t%u\t%d\t%.9f\t%.9f\t%u\t%u\t1\n",
          rows, count, (unsigned long)entries, weight, head, seed, repeat + 1U,
          run->requested, run->effective, run->packed, run->wall, run->cpu,
          run->dependencies, run->attempts);
        if (fflush(tsv) != 0) croak("lanczos-bench: TSV write failed");
      }
    }
  }
  for (repeat = 0; repeat < repeats; repeat++) times[repeat] = runs[repeat].wall;
  serial_median = bench_median(times, repeats);
  printf("\n threads  actual   median wall   min..max wall      speedup   median CPU\n");
  for (point = 0; point < ncounts; point++) {
    double median, low, high, cpu;
    uint32_t effective = runs[(size_t)point * repeats].effective;
    for (repeat = 0; repeat < repeats; repeat++)
      times[repeat] = runs[(size_t)point * repeats + repeat].wall;
    median = bench_median(times, repeats);
    low = times[0]; high = times[repeats - 1U];
    for (repeat = 0; repeat < repeats; repeat++) {
      bench_result_t *run = runs + (size_t)point * repeats + repeat;
      times[repeat] = run->cpu;
      if (run->effective != effective) effective = 0U;
    }
    cpu = bench_median(times, repeats);
    printf(" %7u  %6u  %10.3fs  %7.3f..%-7.3f  %7.2fx  %10.3fs\n",
      counts[point], effective, median, low, high, serial_median / median, cpu);
  }
  if (tsv != NULL && fclose(tsv) != 0) croak("lanczos-bench: TSV close failed");
  free(reference); free(runs); free(times);
  for (c = 0; c < count; c++) free(cols[c].data);
  free(cols);
  return 0;
}
