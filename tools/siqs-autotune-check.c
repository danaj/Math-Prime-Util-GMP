/* Optional aggregate/pool tests, not part of the Perl module's test suite.
 * Build with the same dependencies as siqs-sieve-bench.
 * Copyright (c) 2026 Dana Jacobsen. See LICENSE for redistribution terms. */
#define SIQS_SIEVE_BENCH_MAIN siqs_autotune_unused_benchmark_main
#include "siqs-sieve-bench.c"

static void check(int condition, const char *message) {
  if (!condition) croak("autotune check: %s", message);
}

static void check_decision(const autotune_input_t *inputs,
                            const bench_case_t *choices, uint32_t chosen,
                            const char *expected) {
  FILE *output = tmpfile();
  int saved_stdout;
  char text[512];
  size_t length;
  check(output != NULL, "decision capture file");
  fflush(stdout);
  saved_stdout = dup(STDOUT_FILENO);
  check(saved_stdout >= 0, "decision capture descriptor");
  check(dup2(fileno(output), STDOUT_FILENO) >= 0, "redirect decision output");
  autotune_decision(inputs, choices, 3, chosen, 2);
  fflush(stdout);
  check(dup2(saved_stdout, STDOUT_FILENO) >= 0, "restore decision output");
  close(saved_stdout);
  rewind(output);
  length = fread(text, 1, sizeof(text)-1U, output);
  text[length] = '\0';
  fclose(output);
  check(strstr(text, expected) != NULL, "decision wording/strength");
}

static void check_scores(void) {
  autotune_input_t *inputs = (autotune_input_t *)siqs_calloc(2, sizeof(*inputs));
  bench_case_t choices[3];
  double saved[3];
  static const double times[2][3] = { { 100, 160, 120 }, { 400, 240, 300 } };
  uint32_t p, i;
  memset(choices, 0, sizeof(choices));
  for (i = 0; i < 3; i++) {
    choices[i].maximum = (64U + i*32U)*1024U;
    choices[i].alias = i;
    choices[i].retained = 1;
    for (p = 0; p < 2; p++) {
      inputs[p].cases[i].wall_per_poly = inputs[p].wall[i];
      inputs[p].cases[i].cpu_per_poly = inputs[p].cpu[i];
      inputs[p].wall[i][0] = inputs[p].cpu[i][0] = times[p][i];
    }
  }
  autotune_scores(inputs, 2, choices, 3, 1, 0);
  check(autotune_choose(choices, 3) == 2, "aggregate must find compromise winner");
  for (i = 0; i < 3; i++) {
    saved[i] = choices[i].median;
    inputs[1].wall[i][0] *= 1000.0;
  }
  autotune_scores(inputs, 2, choices, 3, 1, 0);
  for (i = 0; i < 3; i++)
    check(fabs(saved[i] - choices[i].median) < 1e-12,
          "slower input must not acquire extra aggregate weight");

  choices[0].median = 1.6;
  choices[1].median = 1.015;
  choices[2].median = 1.0;
  check(autotune_choose(choices, 3) == 1, "prefer smaller within 2 percent");
  check_decision(inputs, choices, 1,
    "Tied: choosing 96k over 128k (1.50% aggregate gap; prefer smaller)");
  check(autotune_uncertain(inputs, 2, choices, 3), "close scores extend");
  choices[1].median = 1.0;
  choices[2].median = 1.015;
  check_decision(inputs, choices, 1,
    "Tied: choosing 96k over 128k (1.50% aggregate gap; prefer smaller)");
  choices[1].median = 1.1;
  choices[2].median = 1.0;
  check(autotune_choose(choices, 3) == 2, "clear winner retained");
  check_decision(inputs, choices, 2,
    "128k: aggregate 1.100x faster than 96k across 2 inputs (1/2 wins)");
  check(!autotune_uncertain(inputs, 2, choices, 3), "clear/stable scores do not extend");
  inputs[1].cases[2].maximum_time *= 1.10;
  check(autotune_uncertain(inputs, 2, choices, 3), "noisy contender extends");
  autotune_cull(choices, 3, 0.05);
  check(!choices[0].retained && choices[1].retained && choices[2].retained,
        "culling retains at least two");
  check(!autotune_uncertain(inputs, 2, choices + 2, 1), "singleton does not extend");

  for (p = 0; p < 2; p++) {
    inputs[p].cases[0].block = 0;
    inputs[p].cases[1].block = 100;
    inputs[p].cases[2].block = 100;
  }
  inputs[1].cases[2].block = 200;
  autotune_aliases(inputs, 2, choices, 3);
  check(choices[2].retained, "one-input alias is not a global alias");
  inputs[1].cases[2].block = 100;
  autotune_aliases(inputs, 2, choices, 3);
  check(!choices[2].retained && choices[2].alias == 1, "all-input alias is reused");
  choices[1].median = 1.0;
  choices[0].median = 1.6;
  choices[2].maximum = choices[1].maximum - 1024U;
  check(autotune_choose(choices, 3) == 2, "smaller alias may represent winner");
  free(inputs);
  puts("PASS aggregate: compromise, scale invariance, ties, culling, extension and aliases");
}

static void check_pool_inputs(void) {
  bench_pool_t pool;
  bench_worker_t *workers[2];
  bench_case_t geometry;
  mpz_t n[2];
  uint32_t p, repeat;
  double wall, cpu;
  uint64_t checksum;
  memset(&pool, 0, sizeof(pool));
  memset(&geometry, 0, sizeof(geometry));
  pool.count = 3;
  geometry.maximum = 64U*1024U;
  geometry.minimum = 5U*geometry.maximum/2U;
  prime_iterator_global_startup();
  verbose_level = 0;
  for (p = 0; p < 2; p++) {
    mpz_init(n[p]);
    bench_validate_input(n[p], p ? autotune_large_n : bench_default_n);
    workers[p] = bench_prepare_workers(&pool, n[p]);
  }
  pool.workers = workers[0];
  bench_pool_start(&pool);
  for (repeat = 0; repeat < 3; repeat++)
    for (p = 0; p < 2; p++) {
      bench_pool_select(&pool, workers[p]);
      bench_configure(&pool, &geometry);
      bench_run_count(&pool, 1, 4, &wall, &cpu, &checksum);
      bench_run(&pool, 4, &wall, &cpu, &checksum);
    }
  bench_pool_stop(&pool);
  for (p = 0; p < 2; p++) {
    bench_clear_workers(workers[p], pool.count);
    mpz_clear(n[p]);
  }
  prime_iterator_global_shutdown();
  puts("PASS pool: repeated serial/threaded switches, exact work and byte references");
}

int main(void) {
  check_scores();
  check_pool_inputs();
  return 0;
}
