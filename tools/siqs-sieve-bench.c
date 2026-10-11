/*
 * Fixed-work SIQS byte-sieve block-size screen. No full factorization.
 * See README-siqs-sieve-bench.txt. Includes the production private interfaces,
 * as siqs-check does; no production hooks, policies or binaries are changed.
 *
 * cc -O3 -march=native -DSTANDALONE -DPSIQS -pthread \
 *   -o siqs-sieve-bench tools/siqs-sieve-bench.c lanczos.c \
 *   prime_iterator.c squfof126.c -lgmp -lm
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
#error "Build siqs-sieve-bench with -DSTANDALONE -DPSIQS -pthread"
#endif
#include "../ptypes.h"
/* The production allocator is the sole source of balancing/activation math.
 * These are runtime inputs only in this tool. Its parser enforces the size
 * bounds normally enforced by siqs.c's compile-time check (identifiers in
 * that #if evaluate to zero). Change them only while all workers are idle. */
#ifdef SIQS_SIEVE_BLOCK_SIZE
#error "Use --sizes rather than a compile-time block-size override"
#endif
#ifdef SIQS_SIEVE_BLOCK_MIN_LENGTH
#error "Use --min-length rather than a compile-time threshold override"
#endif
static uint32_t bench_block_maximum, bench_block_minimum;
#define SIQS_SIEVE_BLOCK_SIZE bench_block_maximum
#define SIQS_SIEVE_BLOCK_MIN_LENGTH bench_block_minimum
#define main siqs_sieve_bench_unused_driver_main
#include "../mpu-siqs.c"
#undef main
#include "../siqs.c"

#define BENCH_MAX_SIZES 64U
#define BENCH_MAX_ITERATIONS UINT64_C(100000000)
static const char *bench_default_n =
  "1522605027922533360535618378132637429718068114961380688657908494580122963258952897654000350692006139";

typedef struct bench_pool_t bench_pool_t;
typedef struct {
  psiqs_worker_t prepared;
  bench_pool_t *pool;
  pthread_t thread;
  uint32_t index;
  uint8_t *reference;
  uint64_t checksum;
  uint64_t completed_iterations;
} bench_worker_t;

struct bench_pool_t {
  bench_worker_t *workers, *anchors;
  uint32_t count, pending;
  uint64_t generation, iterations;
  pthread_mutex_t mutex;
  pthread_cond_t work, done;
  int stop;
};

typedef struct {
  uint32_t maximum, minimum, block, blocks, dense, alias;
  uint64_t iterations;
  double *wall_per_poly, *cpu_per_poly;
  double median, minimum_time, maximum_time, median_cpu;
  double last_wall;
  int retained;
} bench_case_t;

static double bench_wall(void) {
  struct timespec now;
  if (clock_gettime(CLOCK_MONOTONIC, &now) != 0)
    croak("siqs-sieve-bench: clock_gettime failed: %s", strerror(errno));
  return (double)now.tv_sec + (double)now.tv_nsec / 1e9;
}

static double bench_cpu(void) {
  struct rusage usage;
  if (getrusage(RUSAGE_SELF, &usage) != 0)
    croak("siqs-sieve-bench: getrusage failed: %s", strerror(errno));
  return (double)usage.ru_utime.tv_sec + (double)usage.ru_stime.tv_sec
       + ((double)usage.ru_utime.tv_usec + usage.ru_stime.tv_usec) / 1e6;
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

/* The prepared real polynomial remains fixed. Reading a rotating byte after
 * every call makes every completed sieve observable, not just the last one.
 * No polynomial/root generation, candidate scan, resieve or evaluation here. */
static void bench_sieve(bench_worker_t *worker, uint64_t iterations) {
  siqs_ctx_t *ctx = &worker->prepared.ctx;
  uint32_t pos = worker->index % ctx->sieve_length;
  uint64_t i, checksum = 0;
  for (i = 0; i < iterations; i++) {
    siqs_run_sieve(ctx);
    checksum += ctx->sieve[pos];
    pos += 7919U;
    if (pos >= ctx->sieve_length) pos -= ctx->sieve_length;
  }
  worker->checksum = checksum;
  worker->completed_iterations = iterations;
}

static void *bench_thread(void *argument) {
  bench_worker_t *worker = (bench_worker_t *)argument;
  bench_pool_t *pool = worker->pool;
  uint64_t generation = 0;
  pthread_mutex_lock(&pool->mutex);
  for (;;) {
    uint64_t iterations;
    bench_worker_t *job;
    while (!pool->stop && generation == pool->generation)
      pthread_cond_wait(&pool->work, &pool->mutex);
    if (pool->stop) break;
    generation = pool->generation;
    iterations = pool->iterations;
    /* Scratch may change only while the pool is parked. Anchors live until
     * every thread has joined; copy the active job while holding the mutex. */
    job = &pool->workers[worker->index];
    pthread_mutex_unlock(&pool->mutex);
    bench_sieve(job, iterations);
    pthread_mutex_lock(&pool->mutex);
    if (--pool->pending == 0) pthread_cond_signal(&pool->done);
  }
  pthread_mutex_unlock(&pool->mutex);
  return NULL;
}

/* Caller is worker zero: --threads 1 is truly serial, with no worker thread.
 * The persistent pool dispatches once per timed batch, never per polynomial. */
static void bench_run_count(bench_pool_t *pool, uint32_t active,
                            uint64_t iterations, double *wall, double *cpu,
                            uint64_t *checksum) {
  double wall_start, cpu_start;
  uint32_t i;
  /* Serial screening does not wake the parked pool.  Threaded batches always
   * use its entire requested count; no partial-pool scheduling is needed. */
  if (active != 1U && active != pool->count)
    croak("siqs-sieve-bench: invalid active worker count");
  if (active > 1U) {
    pthread_mutex_lock(&pool->mutex);
    pool->iterations = iterations;
    pool->pending = pool->count - 1U;
    pool->generation++;
  }
  for (i = 0; i < active; i++) pool->workers[i].completed_iterations = 0;
  cpu_start = bench_cpu();
  wall_start = bench_wall();
  if (active > 1U) {
    pthread_cond_broadcast(&pool->work);
    pthread_mutex_unlock(&pool->mutex);
  }
  bench_sieve(&pool->workers[0], iterations);
  if (active > 1U) {
    pthread_mutex_lock(&pool->mutex);
    while (pool->pending != 0)
      pthread_cond_wait(&pool->done, &pool->mutex);
  }
  *wall = bench_wall() - wall_start;
  *cpu = bench_cpu() - cpu_start;
  if (active > 1U) pthread_mutex_unlock(&pool->mutex);
  *checksum = 0;
  for (i = 0; i < active; i++) {
    siqs_ctx_t *ctx = &pool->workers[i].prepared.ctx;
    if (pool->workers[i].completed_iterations != iterations)
      croak("siqs-sieve-bench: wrong input scratch used by worker %u", i);
    if (memcmp(ctx->sieve, pool->workers[i].reference, ctx->sieve_length) != 0)
      croak("siqs-sieve-bench: sieve mismatch after timed run (worker %u)", i);
    *checksum += pool->workers[i].checksum;
  }
}

static void bench_run(bench_pool_t *pool, uint64_t iterations,
                       double *wall, double *cpu, uint64_t *checksum) {
  bench_run_count(pool, pool->count, iterations, wall, cpu, checksum);
}

static void bench_configure(bench_pool_t *pool, bench_case_t *test) {
  uint32_t i;
  bench_block_maximum = test->maximum;
  bench_block_minimum = test->minimum;
  for (i = 0; i < pool->count; i++) {
    bench_worker_t *worker = &pool->workers[i];
    siqs_ctx_t *ctx = &worker->prepared.ctx;
    if (ctx->block_prime_count != 0) {
      free(ctx->block_root1); free(ctx->block_root2); free(ctx->block_step);
    }
    ctx->block_root1 = ctx->block_root2 = ctx->block_step = NULL;
    ctx->block_prime_count = ctx->block_length = 0;
    siqs_block_workspace_allocate(ctx);
    siqs_run_sieve(ctx);
    if (memcmp(ctx->sieve, worker->reference, ctx->sieve_length) != 0)
      croak("siqs-sieve-bench: blocked/unblocked mismatch (worker %u)", i);
    if (i == 0) {
      test->dense = ctx->block_prime_count;
      test->block = test->dense ? ctx->block_length : 0;
      test->blocks = test->block
        ? (ctx->sieve_length - 1U) / test->block + 1U : 0;
    } else if (test->dense != ctx->block_prime_count ||
               test->block != (ctx->block_prime_count ? ctx->block_length : 0)) {
      croak("siqs-sieve-bench: workers have different block geometry");
    }
  }
}

static uint32_t bench_unsigned(const char *text, uint32_t maximum) {
  unsigned long value;
  char *end;
  errno = 0;
  if (*text < '0' || *text > '9')
    croak("siqs-sieve-bench: invalid unsigned value %s", text);
  value = strtoul(text, &end, 10);
  if (errno || *end || value > maximum)
    croak("siqs-sieve-bench: invalid unsigned value %s", text);
  return (uint32_t)value;
}

static uint32_t bench_sizes(const char *text, bench_case_t *cases) {
  uint32_t count = 1;
  /* Always include the unblocked control first. */
  cases[0].maximum = 0;
  while (*text) {
    unsigned long kib;
    char *end;
    uint32_t i, bytes;
    errno = 0;
    if (*text < '0' || *text > '9') croak("siqs-sieve-bench: invalid size list");
    kib = strtoul(text, &end, 10);
    if (errno || (*end && *end != ',') ||
        (kib != 0 && (kib < 4 || kib > 1024)))
      croak("siqs-sieve-bench: sizes must be 0 or 4..1024 KiB");
    bytes = (uint32_t)kib * 1024U;
    for (i = 0; i < count && cases[i].maximum != bytes; i++) { }
    if (i == count) {
      if (count == BENCH_MAX_SIZES) croak("siqs-sieve-bench: too many sizes");
      cases[count++].maximum = bytes;
    }
    if (*end == ',' && end[1] == '\0') croak("siqs-sieve-bench: empty final size");
    text = *end ? end + 1 : end;
  }
  return count;
}

static void bench_write_result(FILE *tsv, const bench_pool_t *pool,
                                const bench_case_t *test, uint32_t round,
                                uint32_t active, double wall, double cpu,
                                uint64_t checksum, const char *stage) {
  const siqs_ctx_t *ctx = &pool->workers[0].prepared.ctx;
  double total = (double)test->iterations * active;
  if (!tsv) return;
  if (stage) fprintf(tsv, "%s\t", stage);
  gmp_fprintf(tsv, "%Zd\t%u\t%lu\t%u\t%u\t%u\t%u\t", ctx->n,
    ctx->params.bits, ctx->multiplier, ctx->params.poly_d,
    ctx->params.q_count, ctx->params.fb_size, ctx->sieve_length);
  fprintf(tsv, "%u\t%u\t%u\t%u\t%u\t%u\t%u\t%llu\t%.9f\t%.9f\t%.6f\t%.6f\t%llu\n",
    round + 1U, active, test->maximum, test->minimum, test->block,
    test->blocks, test->dense, (unsigned long long)test->iterations,
    wall, cpu, wall / total * 1e6, cpu / total * 1e6,
    (unsigned long long)checksum);
  if (fflush(tsv)) croak("siqs-sieve-bench: could not write results");
}

static void bench_validate_input(mpz_t n, const char *number) {
  const char *p;
  for (p = number; *p; p++)
    if (*p < '0' || *p > '9')
      croak("siqs-sieve-bench: input must be a positive decimal integer");
  if (mpz_set_str(n, number, 10) != 0 || mpz_sizeinbase(n, 2) < 100U ||
      mpz_sizeinbase(n, 2) > MPU_SIQS_MAX_BITS || mpz_even_p(n) ||
      siqs_is_prob_prime(n) || mpz_perfect_power_p(n))
    croak("siqs-sieve-bench: use an odd, non-power composite of 100..%u bits",
          (unsigned)MPU_SIQS_MAX_BITS);
}

static bench_worker_t *bench_prepare_workers(bench_pool_t *pool, const mpz_t n) {
  siqs_ctx_t master;
  siqs_poly_t dispatch;
  siqs_factor_array_t factors;
  bench_worker_t *workers;
  uint32_t i;
  /* References are unblocked regardless of the previously measured input. */
  bench_block_maximum = bench_block_minimum = 0;
  siqs_factor_array_init(&factors, n);
  siqs_ctx_init(&master, n, n, &factors, NULL, 0);
  if (!siqs_ctx_allocate(&master))
    croak("siqs-sieve-bench: input has a factor in its factor base; choose another input");
  siqs_poly_init(&master, &dispatch);
  workers = (bench_worker_t *)siqs_calloc(pool->count, sizeof(*workers));
  gmp_printf("Input %Zd (%u bits)\n", n, master.params.bits);
  printf("k=%lu d=%u q=%u FB=%u sieve=%u B (%.2f KiB), first prime=%u; workers=%u\n",
    master.multiplier, master.params.poly_d, master.params.q_count,
    master.params.fb_size, master.sieve_length, master.sieve_length / 1024.0,
    master.prime[master.params.sieve_start], pool->count);
  puts("Preparing one distinct real A/polynomial per worker; setup is not timed.");
  fflush(stdout);
  for (i = 0; i < pool->count; i++) {
    bench_worker_t *worker = &workers[i];
    siqs_ctx_t *ctx = &worker->prepared.ctx;
    worker->pool = pool;
    worker->index = i;
    psiqs_worker_init(&worker->prepared, &master, i);
    if (!psiqs_assign_family(&master, &dispatch, &worker->prepared, UINT32_MAX))
      croak("siqs-sieve-bench: could not prepare distinct A families");
    siqs_first_B_and_roots(ctx, &worker->prepared.poly);
    siqs_set_family_sieve_initial(ctx);
    siqs_run_sieve(ctx);
    worker->reference = (uint8_t *)siqs_malloc(ctx->sieve_length);
    memcpy(worker->reference, ctx->sieve, ctx->sieve_length);
  }
  siqs_poly_clear(&master, &dispatch);
  siqs_ctx_clear(&master);
  free(factors.primality);
  gmp_siqs_free(factors.values, factors.count);
  return workers;
}

static void bench_pool_start(bench_pool_t *pool) {
  uint32_t i;
  pool->anchors = pool->workers;
  if (pthread_mutex_init(&pool->mutex, NULL) ||
      pthread_cond_init(&pool->work, NULL) || pthread_cond_init(&pool->done, NULL))
    croak("siqs-sieve-bench: cannot initialize worker synchronization");
  for (i = 1; i < pool->count; i++)
    if (pthread_create(&pool->anchors[i].thread, NULL, bench_thread, &pool->anchors[i]))
      croak("siqs-sieve-bench: cannot create all requested threads");
}

static void bench_pool_select(bench_pool_t *pool, bench_worker_t *workers) {
  pthread_mutex_lock(&pool->mutex);
  if (pool->pending != 0) croak("siqs-sieve-bench: cannot switch a running pool");
  pool->workers = workers;
  pthread_mutex_unlock(&pool->mutex);
}

static void bench_pool_stop(bench_pool_t *pool) {
  uint32_t i;
  pthread_mutex_lock(&pool->mutex);
  pool->stop = 1;
  pthread_cond_broadcast(&pool->work);
  pthread_mutex_unlock(&pool->mutex);
  for (i = 1; i < pool->count; i++)
    if (pthread_join(pool->anchors[i].thread, NULL))
      croak("siqs-sieve-bench: join failed");
  pthread_cond_destroy(&pool->done);
  pthread_cond_destroy(&pool->work);
  pthread_mutex_destroy(&pool->mutex);
}

static void bench_clear_workers(bench_worker_t *workers, uint32_t count) {
  uint32_t i;
  for (i = 0; i < count; i++) {
    free(workers[i].reference);
    psiqs_worker_clear(&workers[i].prepared);
  }
  free(workers);
}

#include "siqs-autotune.inc.c"

static void bench_usage(void) {
  puts("usage: siqs-sieve-bench [options] [INTEGER ...]\n"
    "  --threads N     Sieve workers (default 1; no automatic thread sweep).\n"
    "  --sizes LIST    Maxima in KiB (default 0,16,32,64; off always included).\n"
    "  --seconds S     Approximate seconds per measured batch (default 1).\n"
    "  --repeat N      Sweeps (default 3; alternating direction).\n"
    "  --polynomials N Fixed calls per worker instead of time calibration.\n"
    "  --output FILE   Per-run TSV; existing files are never replaced.\n"
    "  --min-length K  Advanced: override activation threshold in KiB.\n"
    "                  Default is production's 2.5 x each maximum; 0 bypasses it.\n"
    "  --autotune      Staged serial/threaded screen; see tools/siqs-autotune.pl.\n"
    "  --work-seconds S Collection probe for autotune (default 2; 0 skips it).\n"
    "Normal mode: one input, default RSA-100 (330 bits); only byte-sieving is timed.\n"
    "Autotune: defaults RSA-100/RSA-110 (330/364 bits), or supply distinct inputs.\n"
    "Reports actual balanced blocks and reuses geometrically identical cases.");
}

#ifndef SIQS_SIEVE_BENCH_MAIN
# define SIQS_SIEVE_BENCH_MAIN main
#endif
int SIQS_SIEVE_BENCH_MAIN(int argc, char **argv) {
  bench_case_t cases[BENCH_MAX_SIZES];
  bench_pool_t pool;
  mpz_t n;
  const char *number = bench_default_n, *size_list = "0,16,32,64", *output = NULL;
  const char *numbers[AUTOTUNE_MAX_INPUTS];
  uint32_t threads = 1, repeats = 3, minimum = UINT32_MAX;
  uint32_t count, i, j, round, best = 0;
  uint64_t fixed_iterations = 0;
  double seconds = 1.0, work_seconds = 2.0;
  FILE *tsv = NULL;
  int arg, number_seen = 0, short_runs = 0, autotune = 0, fixed_options = 0;
  memset(cases, 0, sizeof(cases));
  memset(&pool, 0, sizeof(pool));
  for (arg = 1; arg < argc; arg++) {
    const char *option = argv[arg], *value;
    if (strcmp(option, "--help") == 0) { bench_usage(); return 0; }
    if (strcmp(option, "--autotune") == 0) { autotune = 1; continue; }
    if (*option != '-') {
      if ((uint32_t)number_seen == AUTOTUNE_MAX_INPUTS)
        croak("siqs-sieve-bench: too many inputs");
      numbers[number_seen++] = option;
      continue;
    }
    if (++arg == argc) croak("siqs-sieve-bench: missing value for %s", option);
    value = argv[arg];
    if (!strcmp(option, "--threads") || !strcmp(option, "-threads"))
      threads = bench_unsigned(value, PSIQS_MAX_THREADS);
    else if (!strcmp(option, "--repeat")) {
      repeats = bench_unsigned(value, 1000U); fixed_options = 1;
    }
    else if (!strcmp(option, "--sizes")) size_list = value;
    else if (!strcmp(option, "--output")) output = value;
    else if (!strcmp(option, "--min-length"))
      minimum = bench_unsigned(value, 1048576U) * 1024U;
    else if (!strcmp(option, "--polynomials")) {
      fixed_options = 1;
      fixed_iterations = bench_unsigned(value, (uint32_t)BENCH_MAX_ITERATIONS);
      if (!fixed_iterations) croak("siqs-sieve-bench: polynomials must be positive");
    } else if (!strcmp(option, "--seconds")) {
      char *end;
      fixed_options = 1;
      errno = 0;
      seconds = strtod(value, &end);
      if (errno || *end || !(seconds > 0.0 && seconds <= 3600.0))
        croak("siqs-sieve-bench: seconds must be in (0,3600]");
    } else if (!strcmp(option, "--work-seconds")) {
      char *end;
      errno = 0;
      work_seconds = strtod(value, &end);
      if (errno || *end || !(work_seconds >= 0.0 && work_seconds <= 60.0))
        croak("siqs-sieve-bench: work-seconds must be in [0,60]");
    } else croak("siqs-sieve-bench: unknown option %s", option);
  }
  if (autotune && fixed_options)
    croak("siqs-sieve-bench: autotune supplies its own seconds/repeat/polynomials");
  if (autotune && minimum != UINT32_MAX)
    croak("siqs-sieve-bench: autotune retains the production activation threshold");
  if (!autotune && number_seen > 1)
    croak("siqs-sieve-bench: supply only one input outside autotune mode");
  if (!threads || !repeats) croak("siqs-sieve-bench: threads/repeat must be positive");
  if (!*size_list) croak("siqs-sieve-bench: empty size list");
  count = bench_sizes(size_list, cases);
  for (i = 0; i < count; i++) {
    cases[i].minimum = minimum == UINT32_MAX ? 5U * cases[i].maximum / 2U : minimum;
    cases[i].alias = i;
    if (!autotune) {
      cases[i].wall_per_poly = (double *)siqs_calloc(repeats, sizeof(double));
      cases[i].cpu_per_poly = (double *)siqs_calloc(repeats, sizeof(double));
    }
  }
  if (output) {
    int fd = open(output, O_WRONLY | O_CREAT | O_EXCL, 0644);
    if (fd < 0) croak("siqs-sieve-bench: cannot create %s: %s", output, strerror(errno));
    tsv = fdopen(fd, "w");
    if (!tsv) croak("siqs-sieve-bench: fdopen failed: %s", strerror(errno));
    if (autotune) fprintf(tsv, "stage\t");
    fprintf(tsv, "n\tbits\tmultiplier\td\tq\tfb_size\tsieve_bytes\trepeat\tthreads\tmax_bytes\tmin_length_bytes\tblock_bytes\tblocks\tdense_primes\titerations_per_worker\twall_seconds\tcpu_seconds\twall_us_per_poly\tcpu_us_per_poly\tchecksum\n");
    fflush(tsv);
  }

  verbose_level = 0;
  prime_iterator_global_startup();
  if (autotune) {
    autotune_run(numbers, (uint32_t)number_seen, threads, cases, count, tsv, work_seconds);
    if (tsv && fclose(tsv)) croak("siqs-sieve-bench: could not close results");
    prime_iterator_global_shutdown();
    return 0;
  }
  if (number_seen) number = numbers[0];
  mpz_init(n);
  bench_validate_input(n, number);
  pool.count = threads;
  pool.workers = bench_prepare_workers(&pool, n);
  bench_pool_start(&pool);
  /* Calibrate once per distinct geometry, not once per repeated sweep. */
  for (i = 0; i < count; i++) {
    double wall, cpu;
    uint64_t checksum, iterations = 8;
    bench_configure(&pool, &cases[i]);
    for (j = 0; j < i; j++)
      if (cases[j].block == cases[i].block) { cases[i].alias = cases[j].alias; break; }
    printf("max=%u KiB: ", cases[i].maximum / 1024U);
    if (cases[i].block)
      printf("%u balanced blocks <= %u B (%.2f KiB), %u dense primes",
        cases[i].blocks, cases[i].block, cases[i].block / 1024.0, cases[i].dense);
    else printf("unblocked");
    if (cases[i].alias != i) {
      printf("; identical to max=%u KiB, results reused\n",
             cases[cases[i].alias].maximum / 1024U);
      continue;
    }
    putchar('\n');
    bench_run(&pool, iterations, &wall, &cpu, &checksum); /* Warmup. */
    if (fixed_iterations) cases[i].iterations = fixed_iterations;
    else {
      do {
        bench_run(&pool, iterations, &wall, &cpu, &checksum);
        if (wall >= 0.05 || iterations >= BENCH_MAX_ITERATIONS / 2U) break;
        iterations *= 2U;
      } while (1);
      {
        double estimate = ceil((double)iterations * seconds / wall);
        if (!(estimate >= 1.0 && estimate <= (double)BENCH_MAX_ITERATIONS))
          croak("siqs-sieve-bench: calibrated work outside supported range");
        cases[i].iterations = (uint64_t)estimate;
      }
    }
    fflush(stdout);
  }

  for (round = 0; round < repeats; round++) {
    for (j = 0; j < count; j++) {
      bench_case_t *test;
      double wall, cpu, total;
      uint64_t checksum;
      i = round % 2U ? count - 1U - j : j;
      test = &cases[i];
      if (test->alias != i) continue;
      bench_configure(&pool, test);
      bench_run(&pool, 8, &wall, &cpu, &checksum); /* Untimed warmup batch. */
      bench_run(&pool, test->iterations, &wall, &cpu, &checksum);
      if (wall < 0.05) short_runs = 1;
      total = (double)test->iterations * threads;
      test->wall_per_poly[round] = wall / total;
      test->cpu_per_poly[round] = cpu / total;
      printf("repeat=%u/%u max=%u KiB calls/worker=%llu wall=%.6fs CPU=%.6fs (%.3f us/sieve)\n",
        round + 1U, repeats, test->maximum / 1024U,
        (unsigned long long)test->iterations, wall, cpu, wall / total * 1e6);
      fflush(stdout);
      bench_write_result(tsv, &pool, test, round, threads, wall, cpu,
                          checksum, NULL);
    }
  }
  for (i = 0; i < count; i++) {
    bench_case_t *test = &cases[i];
    if (test->alias != i) continue;
    test->median = bench_median(test->wall_per_poly, repeats);
    test->minimum_time = test->wall_per_poly[0];
    test->maximum_time = test->wall_per_poly[repeats - 1U];
    test->median_cpu = bench_median(test->cpu_per_poly, repeats);
    if (test->median < cases[best].median) best = i;
  }
  puts("\n max KiB  actual B  blocks   median us/sieve      min..max us    speedup   CPU us/sieve");
  for (i = 0; i < count; i++) {
    bench_case_t *test = &cases[i], *timing = &cases[test->alias];
    printf("%8u %9u %7u %17.3f %9.3f..%-9.3f %7.3fx %14.3f%s\n",
      test->maximum / 1024U, test->block, test->blocks, timing->median * 1e6,
      timing->minimum_time * 1e6, timing->maximum_time * 1e6,
      cases[0].median / timing->median, timing->median_cpu * 1e6,
      test->alias != i ? " (identical geometry; reused)" : "");
  }
  if (cases[best].block)
    printf("Fastest measured geometry: max=%u KiB, actual=%u B.\n",
      cases[best].maximum / 1024U, cases[best].block);
  else puts("Fastest measured geometry: unblocked.");
if (short_runs)
    puts("WARNING: batches below 50 ms are smoke tests, not reliable timing comparisons.");
  puts("All logical sieve bytes matched the unblocked reference. No policy changes applied.");
  bench_pool_stop(&pool);
  bench_clear_workers(pool.workers, pool.count);
  for (i = 0; i < count; i++) { free(cases[i].wall_per_poly); free(cases[i].cpu_per_poly); }
  if (tsv && fclose(tsv)) croak("siqs-sieve-bench: could not close results");
  prime_iterator_global_shutdown();
  mpz_clear(n);
  return 0;
}
