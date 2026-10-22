/* Real reduced SIQS matrices: dense versus serial block Lanczos.
 * See README-siqs-la-bench.txt. No production timing/capture hooks.
 * cc -O3 -march=native -DSTANDALONE -UPSIQS -o siqs-la-bench \
 *   tools/siqs-la-bench.c siqs.c ../prime_iterator.c ../squfof126.c -lgmp -lm
 * Copyright (c) 2026 Dana Jacobsen. See LICENSE for redistribution terms.
 */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#include <errno.h>
#include <inttypes.h>
#include <time.h>

#if !defined(STANDALONE) || defined(PSIQS)
#error "Build siqs-la-bench with -DSTANDALONE -UPSIQS (serial)"
#endif
#define la_reduce_matrix bench_original_reduce
#include "../lanczos.c"
#undef la_reduce_matrix
#include "../siqs.h"
#include "../siqs_dep.h"
#include "../../prime_iterator.h"

int siqs_is_prob_prime(const mpz_t n) { return mpz_probab_prime_p(n, 25); }

typedef struct {
  double cpu, wall;
  unsigned int dependencies;
} bench_time_t;

static unsigned int bench_repeats = 6, bench_input_bits;
static unsigned long bench_min_cols = 384, bench_max_cols = 8192;
static unsigned long bench_input_id, bench_matrix_id, bench_total_matrices;
static uint64_t bench_seed = UINT64_C(0x6ad37249f08e15cb);

static double bench_clock(clockid_t clock_id) {
  struct timespec value;
  if (clock_gettime(clock_id, &value) != 0)
    croak("siqs-la-bench: clock_gettime failed\n");
  return (double)value.tv_sec + 1e-9 * (double)value.tv_nsec;
}

static unsigned int bench_bit_count(uint64_t mask) {
  unsigned int count = 0;
  while (mask != 0) { count += (unsigned int)(mask & 1U); mask >>= 1; }
  return count;
}

static int bench_double_compare(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return x < y ? -1 : x > y;
}

static double bench_median(double *values) {
  qsort(values, bench_repeats, sizeof(*values), bench_double_compare);
  return bench_repeats & 1U ? values[bench_repeats / 2U]
       : (values[bench_repeats / 2U - 1U] + values[bench_repeats / 2U]) / 2.0;
}

static la_col_t *bench_copy(unsigned long ncols, const la_col_t *columns) {
  la_col_t *copy;
  unsigned long c;
  New(0, copy, (size_t)ncols, la_col_t);
  for (c = 0; c < ncols; c++) {
    copy[c] = columns[c];
    New(0, copy[c].data, (size_t)columns[c].weight, uint32_t);
    if (columns[c].weight != 0)
      memcpy(copy[c].data, columns[c].data, (size_t)columns[c].weight * sizeof(uint32_t));
  }
  return copy;
}

static void bench_free_columns(unsigned long ncols, la_col_t *columns) {
  unsigned long c;
  for (c = 0; c < ncols; c++) free(columns[c].data);
  free(columns);
}

/* Separate from production's verification; outside the timed solver call. */
static void bench_verify(unsigned long nrows, unsigned long ncols,
                         const la_col_t *columns, const uint64_t *deps, uint64_t mask) {
  uint64_t *parity, present = 0;
  unsigned long c, row;
  if (deps == NULL || mask == 0)
    croak("siqs-la-bench: solver returned no dependencies\n");
  Newz(0, parity, nrows != 0 ? (size_t)nrows : 1U, uint64_t);
  for (c = 0; c < ncols; c++) {
    uint32_t i;
    if ((deps[c] & ~mask) != 0)
      croak("siqs-la-bench: dependency outside result mask\n");
    present |= deps[c];
    for (i = 0; i < columns[c].weight; i++) {
      row = columns[c].data[i];
      if (row >= nrows) croak("siqs-la-bench: row outside matrix\n");
      parity[row] ^= deps[c];
    }
  }
  if (present != mask) croak("siqs-la-bench: empty dependency reported\n");
  for (row = 0; row < nrows; row++)
    if (parity[row] != 0) croak("siqs-la-bench: invalid dependency\n");
  free(parity);
}

static bench_time_t bench_solve(int dense, unsigned long nrows, unsigned long ncols,
                                const la_col_t *columns, uint64_t seed) {
  bench_time_t result;
  la_col_t *copy = bench_copy(ncols, columns);
  uint64_t *dependencies, mask;
  double wall_start = bench_clock(CLOCK_MONOTONIC);
  double cpu_start = bench_clock(CLOCK_PROCESS_CPUTIME_ID);
  if (dense)
    dependencies = la_dense_nullspace(nrows, ncols, copy, &mask);
  else
    dependencies = la_block_lanczos(nrows, 0, ncols, copy, seed, &mask, 0);
  result.cpu = bench_clock(CLOCK_PROCESS_CPUTIME_ID) - cpu_start;
  result.wall = bench_clock(CLOCK_MONOTONIC) - wall_start;
  bench_verify(nrows, ncols, columns, dependencies, mask);
  result.dependencies = bench_bit_count(mask);
  free(dependencies);
  bench_free_columns(ncols, copy);
  return result;
}

void la_reduce_matrix(unsigned long *nrows, unsigned long *ncols,
                      la_col_t *columns, int verbose) {
  double cpu[2][100], wall[2][100], median_cpu[2], median_wall[2];
  uint64_t nnz = 0;
  unsigned int repeat, dense_deps = 0, lanczos_min_deps = 64, lanczos_max_deps = 0;
  unsigned long c;
  bench_original_reduce(nrows, ncols, columns, verbose);
  bench_matrix_id++;
  if (*ncols < bench_min_cols || *ncols > bench_max_cols) return;
  for (c = 0; c < *ncols; c++) nnz += columns[c].weight;
  for (repeat = 0; repeat < bench_repeats; repeat++) {
    unsigned int first = (unsigned int)((bench_input_id + bench_matrix_id + repeat) & 1U);
    unsigned int order;
    /* Each seed gets both launch orientations with even repeat counts. */
    uint64_t seed = bench_seed + (uint64_t)(repeat / 2U);
    for (order = 0; order < 2; order++) {
      unsigned int solver = first ^ order;
      bench_time_t value = bench_solve(solver == 0, *nrows, *ncols, columns, seed);
      cpu[solver][repeat] = value.cpu;
      wall[solver][repeat] = value.wall;
      if (solver == 0) dense_deps = value.dependencies;
      else {
        if (value.dependencies < lanczos_min_deps) lanczos_min_deps = value.dependencies;
        if (value.dependencies > lanczos_max_deps) lanczos_max_deps = value.dependencies;
      }
    }
  }
  for (c = 0; c < 2; c++) {
    median_cpu[c] = bench_median(cpu[c]);
    median_wall[c] = bench_median(wall[c]);
  }
  printf("%lu\t%u\t%lu\t%lu\t%lu\t%" PRIu64 "\t%u\t%.6f\t%.6f\t%.6f\t%.6f\t%.5f\t%u\t%u\t%u\n",
         bench_input_id, bench_input_bits, bench_matrix_id, *nrows, *ncols, nnz,
         bench_repeats, median_cpu[0] * 1000, median_cpu[1] * 1000,
         median_wall[0] * 1000, median_wall[1] * 1000,
         median_cpu[1] / median_cpu[0], dense_deps, lanczos_min_deps, lanczos_max_deps);
  fflush(stdout);
  bench_total_matrices++;
}

static unsigned long bench_number(const char *text, unsigned long maximum) {
  char *end;
  unsigned long value;
  if (text[0] < '0' || text[0] > '9') croak("siqs-la-bench: expected positive integer\n");
  errno = 0;
  value = strtoul(text, &end, 10);
  if (errno || *end != '\0' || value == 0 || value > maximum)
    croak("siqs-la-bench: integer outside option range\n");
  return value;
}

int main(int argc, char **argv) {
  char *line = NULL;
  size_t allocated = 0;
  mpz_t n, product;
  int argument, status = 0;
  for (argument = 1; argument < argc; argument++) {
    if (strcmp(argv[argument], "--help") == 0) {
      puts("siqs-la-bench [--repeats 6] [--min-cols 384] [--max-cols 8192]");
      puts("Read one positive decimal input per line on stdin (at most 245 bits).");
      puts("TSV timings compare current solvers on identical real reduced matrices.");
      return 0;
    }
    if (argument + 1 >= argc) croak("siqs-la-bench: missing option value\n");
    if (strcmp(argv[argument], "--repeats") == 0)
      bench_repeats = (unsigned int)bench_number(argv[++argument], 100);
    else if (strcmp(argv[argument], "--min-cols") == 0)
      bench_min_cols = bench_number(argv[++argument], 32768);
    else if (strcmp(argv[argument], "--max-cols") == 0)
      bench_max_cols = bench_number(argv[++argument], 32768);
    else croak("siqs-la-bench: unknown option\n");
  }
  if (bench_min_cols > bench_max_cols) croak("siqs-la-bench: invalid column range\n");
  fprintf(stderr, "Current dense cutoff: %lu reduced columns; %u repeats, paired seed/order rotation.\n",
          LA_DENSE_CROSSOVER_COLS, bench_repeats);
  puts("input\tbits\tmatrix\trows\tcols\tnnz\trepeats\tdense_cpu_ms\tlanczos_cpu_ms\tdense_wall_ms\tlanczos_wall_ms\tlanczos_over_dense\tdense_deps\tlanczos_min_deps\tlanczos_max_deps");
  prime_iterator_global_startup();
  mpz_inits(n, product, NULL);
  while (getline(&line, &allocated, stdin) != -1) {
    char *start = line, *end;
    mpz_t *factors;
    uint32_t nfactors, i;
    unsigned long matrices_before = bench_total_matrices;
    while (*start == ' ' || *start == '\t' || *start == '\r' || *start == '\n') start++;
    if (*start == '\0' || *start == '#') continue;
    end = start;
    while (*end >= '0' && *end <= '9') end++;
    while (*end == ' ' || *end == '\t' || *end == '\r' || *end == '\n') *end++ = '\0';
    if (*end != '\0' || end == start || mpz_set_str(n, start, 10) != 0 || mpz_cmp_ui(n, 4) < 0) {
      fprintf(stderr, "siqs-la-bench: invalid decimal input\n"); status = 2; break;
    }
    bench_input_bits = (unsigned int)mpz_sizeinbase(n, 2);
    /* Keep this small-matrix tool out of 2LP's nested SIQS cofactor calls. */
    if (bench_input_bits >= 246U) {
      fprintf(stderr, "siqs-la-bench: input exceeds 245 bits\n"); status = 2; break;
    }
    bench_input_id++;
    bench_matrix_id = 0;
    factors = gmp_siqs(n, &nfactors, 0, 0);
    mpz_set_ui(product, 1);
    for (i = 0; i < nfactors; i++) {
      if (mpz_cmp_ui(factors[i], 1) <= 0 || !siqs_is_prob_prime(factors[i])) status = 1;
      mpz_mul(product, product, factors[i]);
    }
    if (nfactors < 2 || mpz_cmp(product, n) != 0) status = 1;
    gmp_siqs_free(factors, nfactors);
    fprintf(stderr, "Input %lu (%u bits): %lu measured matrices.\n", bench_input_id,
            bench_input_bits, bench_total_matrices - matrices_before);
  }
  if (ferror(stdin)) status = 2;
  if (bench_total_matrices == 0 && status == 0) status = 1;
  free(line);
  mpz_clears(n, product, NULL);
  prime_iterator_global_shutdown();
  return status;
}
