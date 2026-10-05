/* Included by siqs-sieve-bench.c, not a separate compilation unit.
 * Advisory hardware screening only: no factoring-policy or source edits.
 * Copyright (c) 2026 Dana Jacobsen. See LICENSE for redistribution terms. */
#include <stdarg.h>

#define AUTOTUNE_MAX_INPUTS 16U
#define AUTOTUNE_MAX_SAMPLES 4U
static const char *autotune_large_n =
  "35794234179725868774991807832568455403003778024228226193532908190484670252364677411513516111204504060317568667";

typedef struct {
  mpz_t n;
  bench_worker_t *workers;
  bench_case_t cases[BENCH_MAX_SIZES];
  double wall[BENCH_MAX_SIZES][AUTOTUNE_MAX_SAMPLES];
  double cpu[BENCH_MAX_SIZES][AUTOTUNE_MAX_SAMPLES];
} autotune_input_t;

/* The wrapper routes these concise lines to stderr in --flags mode, while
 * keeping numeric compiler flags alone on stdout. Detailed rows remain
 * available in normal mode and the optional TSV. Flush before doing work. */
#if defined(__GNUC__) || defined(__clang__)
static void autotune_progress(const char *format, ...)
  __attribute__((format(printf, 1, 2)));
#endif
static void autotune_progress(const char *format, ...) {
  va_list args;
  printf("SIQS_AUTOTUNE_PROGRESS\t");
  va_start(args, format);
  vprintf(format, args);
  va_end(args);
  putchar('\n');
  fflush(stdout);
}

/* Always retain the two fastest surviving distinct geometries. Serial culling
 * remains deliberately loose because concurrency can change the ranking. */
static void autotune_cull(bench_case_t *cases, uint32_t count, double margin) {
  uint32_t i, first = UINT32_MAX, second = UINT32_MAX;
  double best;
  for (i = 0; i < count; i++) {
    if (!cases[i].retained) continue;
    if (first == UINT32_MAX || cases[i].median < cases[first].median) {
      second = first;
      first = i;
    } else if (second == UINT32_MAX || cases[i].median < cases[second].median) {
      second = i;
    }
  }
  best = cases[first].median;
  for (i = 0; i < count; i++)
    if (cases[i].retained && i != first && i != second &&
        cases[i].median > best * (1.0 + margin))
      cases[i].retained = 0;
}

/* A maximum can alias another on one input but not another. Only collapse a
 * candidate globally when its actual geometry matches on EVERY input. */
static void autotune_aliases(autotune_input_t *inputs, uint32_t ninputs,
                              bench_case_t *choices, uint32_t count) {
  uint32_t i, j, p;
  for (i = 0; i < count; i++) {
    choices[i].alias = i;
    for (j = 0; j < i; j++) {
      for (p = 0; p < ninputs; p++)
        if (inputs[p].cases[i].block != inputs[p].cases[j].block) break;
      if (p == ninputs) { choices[i].alias = choices[j].alias; break; }
    }
    choices[i].retained = choices[i].alias == i;
  }
}

static uint64_t autotune_iterations(double estimate) {
  if (!(estimate >= 1.0 && estimate <= (double)BENCH_MAX_ITERATIONS))
    croak("siqs-autotune: calibrated work outside supported range");
  return (uint64_t)ceil(estimate);
}

/* Each input gets equal weight, not weight proportional to its sieve time.
 * Aggregate relative times, not independently selected winners: a compromise
 * can win overall even when it is not the best on either input. */
static void autotune_scores(autotune_input_t *inputs, uint32_t ninputs,
                             bench_case_t *choices, uint32_t count,
                             uint32_t samples, uint32_t fresh_start) {
  uint32_t p, i, s;
  for (i = 0; i < count; i++) choices[i].median = 0.0;
  for (p = 0; p < ninputs; p++) {
    double best = HUGE_VAL;
    for (i = 0; i < count; i++) {
      bench_case_t *test = &inputs[p].cases[i];
      if (!choices[i].retained) continue;
      test->minimum_time = test->maximum_time = test->wall_per_poly[fresh_start];
      for (s = fresh_start + 1U; s < samples; s++) {
        if (test->wall_per_poly[s] < test->minimum_time)
          test->minimum_time = test->wall_per_poly[s];
        if (test->wall_per_poly[s] > test->maximum_time)
          test->maximum_time = test->wall_per_poly[s];
      }
      test->median = bench_median(test->wall_per_poly, samples);
      test->median_cpu = bench_median(test->cpu_per_poly, samples);
      if (!(test->median > 0.0 && test->median < HUGE_VAL))
        croak("siqs-autotune: invalid timing result");
      if (test->median < best) best = test->median;
    }
    for (i = 0; i < count; i++)
      if (choices[i].retained)
        choices[i].median += log(inputs[p].cases[i].median / best);
  }
  for (i = 0; i < count; i++)
    if (choices[i].retained)
      choices[i].median = exp(choices[i].median / ninputs);
}

static void autotune_stage(bench_pool_t *pool, autotune_input_t *inputs,
                            uint32_t ninputs, bench_case_t *choices,
                            uint32_t count, uint32_t stage, double seconds,
                            uint32_t passes, uint32_t first_sample, FILE *tsv) {
  static const char *names[4] = { "serial", "threaded", "final", "extended" };
  uint32_t active = stage == 0 ? 1U : pool->count;
  uint32_t pass, a, p, j, i, k;
  printf("SIQS_AUTOTUNE_PROGRESS\tTesting block stage %u%s (%u inputs, %u thread%s):",
    stage + 1U, stage == 3U ? " (longer confirmation)" : "", ninputs,
    active, active == 1U ? "" : "s");
  for (i = 0; i < count; i++)
    if (choices[i].retained) printf(" %u", choices[i].maximum / 1024U);
  puts(" KiB");
  fflush(stdout);
  for (pass = 0; pass < passes; pass++) {
    uint32_t sample = first_sample + pass;
    for (a = 0; a < ninputs; a++) {
      int measured[BENCH_MAX_SIZES];
      memset(measured, 0, sizeof(measured));
      p = (stage + pass) % 2U ? ninputs - 1U - a : a;
      bench_pool_select(pool, inputs[p].workers);
      for (j = 0; j < count; j++) {
        bench_case_t *test;
        double wall, cpu;
        uint64_t checksum;
        i = (stage + pass + p) % 2U ? count - 1U - j : j;
        if (!choices[i].retained) continue;
        test = &inputs[p].cases[i];
        /* Reuse equal geometry locally, even if the global candidates differ
         * on another input. Never fabricate a duplicate measured TSV row. */
        for (k = 0; k < count; k++)
          if (measured[k] && inputs[p].cases[k].block == test->block) break;
        if (k < count) {
          bench_case_t *source = &inputs[p].cases[k];
          test->iterations = source->iterations;
          test->last_wall = source->last_wall;
          test->wall_per_poly[sample] = source->wall_per_poly[sample];
          test->cpu_per_poly[sample] = source->cpu_per_poly[sample];
          measured[i] = 1;
          printf("  input=%u max=%u KiB: identical to max=%u KiB; reused\n",
            p + 1U, test->maximum / 1024U, source->maximum / 1024U);
          continue;
        }
        bench_configure(pool, test);
        bench_run_count(pool, active, 8, &wall, &cpu, &checksum);
        if (stage == 0) {
          test->iterations = 8;
          do {
            bench_run_count(pool, active, test->iterations, &wall, &cpu, &checksum);
            if (wall >= 0.010 || test->iterations >= BENCH_MAX_ITERATIONS/2U) break;
            test->iterations *= 2U;
          } while (1);
          test->last_wall = wall;
        }
        test->iterations = autotune_iterations(
            (double)test->iterations * seconds / test->last_wall);
        bench_run_count(pool, active, test->iterations, &wall, &cpu, &checksum);
        test->last_wall = wall;
        test->wall_per_poly[sample] = wall / ((double)test->iterations * active);
        test->cpu_per_poly[sample] = cpu / ((double)test->iterations * active);
        measured[i] = 1;
        bench_write_result(tsv, pool, test, sample, active, wall, cpu, checksum,
                            names[stage]);
        printf("  input=%u max=%4u KiB actual=%7u B: %8.3f us/sieve (%.3f s)\n",
          p + 1U, test->maximum / 1024U, test->block,
          test->wall_per_poly[sample] * 1e6, wall);
        fflush(stdout);
      }
    }
  }
  autotune_scores(inputs, ninputs, choices, count, first_sample + passes, first_sample);
}

/* Close aggregate scores or visibly noisy contenders get exactly one longer
 * stage. No timers/confidence machinery in production; this is a bounded tool. */
static int autotune_uncertain(const autotune_input_t *inputs, uint32_t ninputs,
                                const bench_case_t *choices, uint32_t count) {
  uint32_t i, p;
  double best = HUGE_VAL, second = HUGE_VAL;
  for (i = 0; i < count; i++)
    if (choices[i].retained) {
      double t = choices[i].median;
      if (t < best) { second = best; best = t; }
      else if (t < second) second = t;
    }
  if (second == HUGE_VAL) return 0;
  if (second <= best * 1.03) return 1;
  for (i = 0; i < count; i++)
    if (choices[i].retained && choices[i].median <= second)
      for (p = 0; p < ninputs; p++)
        if (inputs[p].cases[i].maximum_time >
            inputs[p].cases[i].minimum_time * 1.05) return 1;
  return 0;
}

static uint32_t autotune_choose(const bench_case_t *choices, uint32_t count) {
  uint32_t i, best = 0, chosen;
  double best_time = HUGE_VAL;
  for (i = 0; i < count; i++)
    if (choices[i].retained && choices[i].median < best_time) {
      best = i;
      best_time = choices[i].median;
    }
  chosen = best;
  for (i = 0; i < count; i++)
    if (choices[choices[i].alias].retained &&
        choices[choices[i].alias].median <= best_time * 1.02 &&
        choices[i].maximum < choices[chosen].maximum) chosen = i;
  return chosen;
}

static void autotune_block_label(char *text, uint32_t maximum) {
  if (maximum) sprintf(text, "%uk", maximum / 1024U);
  else strcpy(text, "unblocked");
}

static void autotune_decision(const autotune_input_t *inputs,
                               const bench_case_t *choices, uint32_t count,
                               uint32_t chosen, uint32_t ninputs) {
  uint32_t i, p, wins = 0, best = UINT32_MAX, second = UINT32_MAX, other;
  char selected[24], compared[24];
  for (i = 0; i < count; i++)
    if (choices[i].retained) {
      if (best == UINT32_MAX || choices[i].median < choices[best].median) {
        second = best; best = i;
      } else if (second == UINT32_MAX || choices[i].median < choices[second].median) {
        second = i;
      }
    }
  autotune_block_label(selected, choices[chosen].maximum);
  if (chosen != best ||
      (second != UINT32_MAX && choices[second].median <= choices[best].median * 1.02)) {
    double selected_time, other_time, gap;
    other = chosen != best ? best : second;
    autotune_block_label(compared, choices[other].maximum);
    selected_time = choices[choices[chosen].alias].median;
    other_time = choices[choices[other].alias].median;
    gap = selected_time > other_time ? selected_time / other_time - 1.0
                                    : other_time / selected_time - 1.0;
    autotune_progress("Tied: choosing %s over %s (%.2f%% aggregate gap; prefer smaller)",
      selected, compared, gap * 100.0);
  } else if (second != UINT32_MAX) {
    autotune_block_label(compared, choices[second].maximum);
    for (p = 0; p < ninputs; p++)
      if (inputs[p].cases[choices[chosen].alias].median < inputs[p].cases[second].median)
        wins++;
    autotune_progress("%s: aggregate %.3fx faster than %s across %u input%s (%u/%u wins)",
      selected, choices[second].median / choices[best].median, compared,
      ninputs, ninputs == 1U ? "" : "s", wins, ninputs);
  } else {
    autotune_progress("Selected block size %s across %u input%s (one distinct geometry)",
      selected, ninputs, ninputs == 1U ? "" : "s");
  }
}

/* Probe the normal serial collection pipeline with a cooperative time limit.
 * Unlike byte-sieve screening this includes A/root generation, candidate
 * finding, resieving, evaluation, graph insertion and readiness checks. Setup,
 * solving and cleanup are excluded. This duplicates only the small driver loop,
 * not the production kernels; there is no timer hook in production siqs.c. */
static uint64_t autotune_work_rate(const mpz_t n, double seconds) {
  siqs_ctx_t ctx;
  siqs_poly_t poly;
  siqs_factor_array_t factors;
  uint32_t families = 0, polys = 0, next_check, check_interval;
  double wall_start, cpu_start, wall, cpu, rate;
  int done = 0;
  uint64_t rounded = 0;
  if (seconds == 0.0) return 0;
  autotune_progress("Starting serial speed calculation (%.1f seconds)", seconds);
  printf("\nCalibrating serial collection work rate for about %.2f seconds...\n", seconds);
  fflush(stdout);
  siqs_factor_array_init(&factors, n);
  siqs_ctx_init(&ctx, n, n, &factors, NULL);
  if (!siqs_ctx_allocate(&ctx) || ctx.factor_found)
    croak("siqs-autotune: work-rate input has a factor in its factor base");
  siqs_poly_init(&ctx, &poly);
  check_interval = ctx.params.fb_size / 128U;
  if (check_interval < SIQS_MATRIX_CHECK_MIN) check_interval = SIQS_MATRIX_CHECK_MIN;
  if (check_interval > SIQS_MATRIX_CHECK_MAX) check_interval = SIQS_MATRIX_CHECK_MAX;
  next_check = ctx.params.fb_size - ctx.params.fb_size / 4U;
  if (next_check < 256U) next_check = 256U;
  cpu_start = bench_cpu();
  wall_start = bench_wall();
  while (!done) {
    if (!siqs_new_family(&ctx, &poly)) break;
    families++;
    for (;;) {
      siqs_sieve_polynomial(&ctx, &poly);
      polys++;
      if (ctx.factor_found || ctx.full_count >= ctx.params.target_relations)
        done = 1;
      if (!done && ctx.full_count >= next_check) {
        uint32_t rows, cols;
        done = siqs_matrix_ready(&ctx, &rows, &cols);
        next_check = ctx.full_count + check_interval;
      }
      if (!done && (polys & 31U) == 0 && bench_wall() - wall_start >= seconds)
        done = 1;
      if (done || !siqs_next_B(&ctx, &poly)) break;
    }
    if (bench_wall() - wall_start >= seconds) done = 1;
  }
  wall = bench_wall() - wall_start;
  cpu = bench_cpu() - cpu_start;
  printf("Collection: %u polynomials, %u families, %.3f wall / %.3f CPU seconds.\n",
    polys, families, wall, cpu);
  if (cpu >= 0.05 && polys >= 32U) {
    rate = (double)ctx.params.half_interval * polys / cpu;
    /* A coarse calibration, not an assertion of seven significant digits. */
    if (rate >= 1000000.0 && rate <= 1e15)
      rounded = (uint64_t)(rate / 1000000.0 + 0.5) * UINT64_C(1000000);
  }
  if (rounded) {
    autotune_progress("Serial speed: %.0fM work units per CPU second", (double)rounded/1000000.0);
    printf("Work rate: approximately %llu M*polynomials per CPU second.\n",
      (unsigned long long)rounded);
  } else {
    autotune_progress("Serial speed: probe too short; keeping the default");
    puts("Work-rate probe too short for a useful recommendation; keep the default.");
  }
  siqs_poly_clear(&ctx, &poly);
  siqs_ctx_clear(&ctx);
  free(factors.primality);
  gmp_siqs_free(factors.values, factors.count);
  return rounded;
}


static void autotune_run(const char **numbers, uint32_t ninputs, uint32_t threads,
                           bench_case_t *choices, uint32_t count, FILE *tsv,
                           double work_seconds) {
  autotune_input_t *inputs;
  bench_pool_t pool;
  uint32_t p, i, j, chosen, total_samples = 2U;
  uint64_t rate;
  int have_blocked = 0, noisy = 0;
  if (ninputs == 0) {
    numbers[0] = bench_default_n;
    numbers[1] = autotune_large_n;
    ninputs = 2;
  }
  inputs = (autotune_input_t *)siqs_calloc(ninputs, sizeof(*inputs));
  memset(&pool, 0, sizeof(pool));
  pool.count = threads;
  for (p = 0; p < ninputs; p++) {
    mpz_init(inputs[p].n);
    bench_validate_input(inputs[p].n, numbers[p]);
    for (j = 0; j < p; j++)
      if (mpz_cmp(inputs[p].n, inputs[j].n) == 0)
        croak("siqs-autotune: duplicate inputs would bias the aggregate");
  }
  autotune_progress("Preparing %u input%s for %u thread%s",
    ninputs, ninputs == 1U ? "" : "s", threads, threads == 1U ? "" : "s");
  for (p = 0; p < ninputs; p++) {
    inputs[p].workers = bench_prepare_workers(&pool, inputs[p].n);
    pool.workers = inputs[p].workers; /* Threads have not yet started. */
    for (i = 0; i < count; i++) {
      bench_case_t *test = &inputs[p].cases[i];
      *test = choices[i];
      test->wall_per_poly = inputs[p].wall[i];
      test->cpu_per_poly = inputs[p].cpu[i];
      bench_configure(&pool, test);
      if (test->block) have_blocked = 1;
      printf("  input=%u max=%u KiB: actual=%u B, blocks=%u\n",
        p + 1U, test->maximum / 1024U, test->block, test->blocks);
    }
  }
  pool.workers = inputs[0].workers;
  bench_pool_start(&pool);
  autotune_aliases(inputs, ninputs, choices, count);
  puts("Equal-weight geometric mean of per-input relative times; cuts +15%, then +5%.");
  autotune_stage(&pool, inputs, ninputs, choices, count, 0, 0.200, 1, 0, tsv);
  autotune_cull(choices, count, 0.15);
  autotune_stage(&pool, inputs, ninputs, choices, count, 1, 0.400, 1, 0, tsv);
  autotune_cull(choices, count, 0.05);
  autotune_stage(&pool, inputs, ninputs, choices, count, 2, 0.400, 2, 0, tsv);
  if (autotune_uncertain(inputs, ninputs, choices, count)) {
    autotune_cull(choices, count, 0.05);
    autotune_stage(&pool, inputs, ninputs, choices, count, 3, 1.000, 2, 2, tsv);
    total_samples = 4U;
  }
  chosen = autotune_choose(choices, count);
  puts("\nAggregate results (relative to the fastest aggregate):");
  {
    double best = HUGE_VAL;
    for (i = 0; i < count; i++)
      if (choices[i].retained && choices[i].median < best) best = choices[i].median;
    for (i = 0; i < count; i++)
      if (choices[i].retained) {
        printf("  max=%4u KiB: +%.2f%%", choices[i].maximum / 1024U,
          (choices[i].median / best - 1.0) * 100.0);
        for (p = 0; p < ninputs; p++) {
          const bench_case_t *test = &inputs[p].cases[i];
          printf(", input %u %.3f us/sieve", p + 1U, test->median * 1e6);
          if (test->maximum_time > test->minimum_time * 1.05) noisy = 1;
        }
        printf(" (%u samples/input)\n", total_samples);
      }
  }
  autotune_decision(inputs, choices, count, chosen, ninputs);
  if (noisy) autotune_progress("Warning: final timing spread exceeds 5%%; repeat on an idle machine");
  puts("All logical sieve bytes matched the unblocked references. No policy changes applied.");
  bench_pool_stop(&pool);
  for (p = 0; p < ninputs; p++) bench_clear_workers(inputs[p].workers, threads);
  bench_block_maximum = choices[chosen].maximum;
  bench_block_minimum = choices[chosen].minimum;
  rate = autotune_work_rate(inputs[0].n, work_seconds);
  if (!have_blocked)
    autotune_progress("No active blocked geometry: keeping the default block policy");
  puts("Advisory aggregate for these inputs/concurrency; speed calibration affects output only.");
  printf("SIQS_AUTOTUNE_FLAGS\t");
  if (have_blocked) printf("-DSIQS_SIEVE_BLOCK_SIZE=%uU", choices[chosen].maximum);
  if (rate) printf("%s-DSIQS_PROGRESS_WORK_PER_SEC=%lluULL",
    have_blocked ? " " : "", (unsigned long long)rate);
  putchar('\n');
  for (p = 0; p < ninputs; p++) mpz_clear(inputs[p].n);
  free(inputs);
}
