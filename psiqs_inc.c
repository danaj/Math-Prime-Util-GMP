/*============================================================================
  Parallel family collection for SIQS.

  Included only by siqs.c when PSIQS is defined.  Do not compile this file
  separately: it uses the shared private SIQS types and static helpers.

  Workers own polynomial/sieve scratch and buffer raw relations.  The caller
  assigns distinct A families and merges each joined batch into one relation
  graph; factor partition updates and relation merging remain serial.

  Copyright (c) 2026 Dana Jacobsen
============================================================================*/

/* Native worker calls into Perl's error/primality adapters need a separate
 * audit.  For now, enable this implementation only in the standalone host. */
#ifndef STANDALONE
# error "PSIQS currently requires the standalone host"
#endif

#include <pthread.h>

typedef struct {
  siqs_ctx_t ctx;
  siqs_poly_t poly;
  siqs_factor_array_t result;
  uint32_t limit;
  uint32_t polynomials;
} psiqs_worker_t;

/* Workers own every writable field.  Reuse scratch allocation and cleanup,
 * but never allocate a relation graph or rebuild the common factor base. */
static void psiqs_worker_init(psiqs_worker_t *worker,
                              const siqs_ctx_t *master, uint32_t index) {
  siqs_ctx_t *ctx = &worker->ctx;
  uint32_t i;
  memset(worker, 0, sizeof(*worker));
  ctx->original_n = master->original_n;
  ctx->params = master->params;
  ctx->multiplier = master->multiplier;
  ctx->largest_fb_prime = master->largest_fb_prime;
  ctx->buffer_relations = 1;
  mpz_init_set(ctx->n, master->n);
  mpz_init_set(ctx->kn, master->kn);
  mpz_init(ctx->eval.y);
  mpz_init(ctx->eval.q);
  mpz_init(ctx->eval.rest);
  ctx->fb = (siqs_fb_t *)siqs_malloc(
      (size_t)ctx->params.fb_size * sizeof(*ctx->fb));
  memcpy(ctx->fb, master->fb,
         (size_t)ctx->params.fb_size * sizeof(*ctx->fb));
  for (i = 0; i < ctx->params.fb_size; i++)
    ctx->fb[i].in_a = 0;
  ctx->cofactor_rng.state = siqs_mix64(
      master->cofactor_rng.state ^
      ((uint64_t)index + 1U) * UINT64_C(0x9e3779b97f4a7c15));
  siqs_factor_array_init(&worker->result, ctx->n);
  ctx->result = &worker->result;
  siqs_workspace_allocate(ctx);
  siqs_poly_init(ctx, &worker->poly);
}

/* Select distinct A values serially, then let each worker initialize roots. */
static int psiqs_assign_family(siqs_ctx_t *master, siqs_poly_t *dispatch,
                               psiqs_worker_t *worker, uint32_t limit) {
  siqs_ctx_t *ctx = &worker->ctx;
  siqs_poly_t *poly = &worker->poly;
  uint32_t i;
  if (!siqs_choose_A(master, dispatch))
    return 0;
  for (i = 0; i < poly->q_count; i++)
    if (poly->a_index[i] < ctx->params.fb_size)
      ctx->fb[poly->a_index[i]].in_a = 0;
  memcpy(poly->a_index, dispatch->a_index,
         (size_t)poly->q_count * sizeof(*poly->a_index));
  for (i = 0; i < poly->q_count; i++)
    ctx->fb[poly->a_index[i]].in_a = 1;
  mpz_set(poly->A, dispatch->A);
  mpz_set(poly->DA, dispatch->DA);
  worker->limit = limit;
  worker->polynomials = 0;
  return 1;
}

/* No shared mutation, logging, graph insertion or solving occurs here. */
static void *psiqs_sieve_family(void *argument) {
  psiqs_worker_t *worker = (psiqs_worker_t *)argument;
  siqs_ctx_t *ctx = &worker->ctx;
  siqs_poly_t *poly = &worker->poly;
  siqs_first_B_and_roots(ctx, poly);
  siqs_set_family_sieve_initial(ctx);
  do {
    siqs_sieve_polynomial(ctx, poly);
    worker->polynomials++;
    if (ctx->factor_found || worker->polynomials >= worker->limit)
      break;
  } while (siqs_next_B(ctx, poly));
  return NULL;
}

/* Copy buffered records into the collector's arena; never transfer pointers
 * between arenas.  All partials enter the same graph, including cross-worker
 * cycles.  The worker buffer is discarded only after its GMP values clear. */
static void psiqs_merge_worker(siqs_ctx_t *master, psiqs_worker_t *worker) {
  siqs_ctx_t *ctx = &worker->ctx;
  uint32_t i;
  master->total_candidates += ctx->total_candidates;
  master->split_attempts += ctx->split_attempts;
  master->split_squfof_or_square += ctx->split_squfof_or_square;
  master->split_rho += ctx->split_rho;
  master->split_failures += ctx->split_failures;
  ctx->total_candidates = 0;
  ctx->split_attempts = 0;
  ctx->split_squfof_or_square = 0;
  ctx->split_rho = 0;
  ctx->split_failures = 0;
  /* A polynomial zero can discover a divisor without emitting a relation.
   * Refine the parent's partition only here, after every worker has joined. */
  if (ctx->factor_found) {
    for (i = 0; i < worker->result.count; i++)
      if (siqs_insert_divisor(master->result, worker->result.values[i]))
        master->factor_found = 1;
  }
  for (i = 0; i < ctx->raw_count; i++) {
    const siqs_raw_relation_t *raw = ctx->raw[i];
    if (!master->factor_found) {
      siqs_factor_t *factors = siqs_eval_reserve_factors(
          master, raw->nfactors);
      siqs_raw_relation_t *copy;
      uint32_t j;
      for (j = 0; j < raw->nfactors; j++) {
        factors[j].row = siqs_raw_factor_row(raw, j);
        factors[j].exponent = siqs_raw_factor_exponent(raw, j);
      }
      copy = siqs_raw_relation_new(master, raw->y, factors,
                                   raw->nfactors, raw->lp1, raw->lp2);
      siqs_accept_raw_relation(master, copy);
    }
  }
  while (ctx->raw_count != 0)
    siqs_raw_relation_free(ctx, ctx->raw[--ctx->raw_count]);
  siqs_raw_arena_clear(&ctx->raw_arena);
}

static void psiqs_worker_clear(psiqs_worker_t *worker) {
  siqs_poly_clear(&worker->ctx, &worker->poly);
  siqs_ctx_clear(&worker->ctx);
  free(worker->result.primality);
  gmp_siqs_free(worker->result.values, worker->result.count);
}

/* Batch-and-join collection intentionally sacrifices load balancing and some
 * stop precision.  It needs no locks or concurrent graph/solver operations.
 * Readiness, progress and the aggregate polynomial guard are checked between
 * batches.  Unfinished final families are bounded by that same global guard. */
static int psiqs_collect_relations(siqs_ctx_t *ctx, siqs_poly_t *dispatch,
                                   uint32_t target,
                                   uint32_t *next_matrix_check,
                                   uint32_t *family_count,
                                   uint32_t *poly_count) {
  psiqs_worker_t *workers = (psiqs_worker_t *)siqs_calloc(
      ctx->nthreads, sizeof(*workers));
  pthread_t *threads = (pthread_t *)siqs_malloc(
      (size_t)ctx->nthreads * sizeof(*threads));
  int *started = (int *)siqs_malloc(
      (size_t)ctx->nthreads * sizeof(*started));
  uint32_t i, check_interval = ctx->params.fb_size / 128;
  uint32_t max_polynomials =
      ctx->params.fb_size > UINT32_MAX / SIQS_MAX_POLYNOMIALS_PER_FB
        ? UINT32_MAX : ctx->params.fb_size * SIQS_MAX_POLYNOMIALS_PER_FB;
  uint64_t report_step =
      (SIQS_PROGRESS_WORK_INTERVAL + ctx->params.half_interval - 1U) /
      ctx->params.half_interval;
  uint64_t next_report = (uint64_t)*poly_count + report_step;
  uint32_t last_report_count = UINT32_MAX, last_report_polys = UINT32_MAX;
  int verbose = siqs_verbose_level(), complete = 0;
  if (max_polynomials < 1000000U)
    max_polynomials = 1000000U;
  if (check_interval < SIQS_MATRIX_CHECK_MIN)
    check_interval = SIQS_MATRIX_CHECK_MIN;
  if (check_interval > SIQS_MATRIX_CHECK_MAX)
    check_interval = SIQS_MATRIX_CHECK_MAX;
  for (i = 0; i < ctx->nthreads; i++)
    psiqs_worker_init(&workers[i], ctx, i);

  for (;;) {
    uint32_t count = 0, remaining = max_polynomials - *poly_count;
    int exhausted = 0;
    if (ctx->factor_found || ctx->full_count >= target) {
      complete = 1;
      break;
    }
    while (count < ctx->nthreads && remaining != 0) {
      uint32_t limit = dispatch->b_limit < remaining
                     ? dispatch->b_limit : remaining;
      if (!psiqs_assign_family(ctx, dispatch, &workers[count], limit)) {
        exhausted = 1;
        break;
      }
      remaining -= limit;
      count++;
    }
    if (count == 0)
      break;
    *family_count += count;
    for (i = 0; i < count; i++) {
      int error = pthread_create(&threads[i], NULL,
                                 psiqs_sieve_family, &workers[i]);
      started[i] = error == 0;
      if (error != 0) {
        fprintf(stderr, "PSIQS: pthread_create: %s; running family serially\n",
                strerror(error));
        psiqs_sieve_family(&workers[i]);
      }
    }
    for (i = 0; i < count; i++) {
      if (started[i]) {
        int error = pthread_join(threads[i], NULL);
        if (error != 0)
          croak("PSIQS: pthread_join failed: %s", strerror(error));
      }
    }
    for (i = 0; i < count; i++) {
      *poly_count += workers[i].polynomials;
      psiqs_merge_worker(ctx, &workers[i]);
    }
    if (ctx->full_count >= *next_matrix_check && !ctx->factor_found) {
      uint32_t rows, columns;
      int ready = siqs_matrix_ready(ctx, &rows, &columns);
      *next_matrix_check = ctx->full_count + check_interval;
      if (ready && verbose > 3) {
        siqs_print_relation_report(ctx, target, *poly_count);
        last_report_count = ctx->full_count;
        last_report_polys = *poly_count;
      }
      if ((ready && verbose > 3) || verbose > 4) {
        printf("# siqs matrix core %u columns, %u rows%s\n",
               columns, rows, ready ? ", ready" : "");
        fflush(stdout);
      }
      if (ready) {
        complete = 1;
        break;
      }
    }
    if (verbose > 3 && (uint64_t)*poly_count >= next_report) {
      siqs_print_relation_report(ctx, target, *poly_count);
      last_report_count = ctx->full_count;
      last_report_polys = *poly_count;
      next_report = (uint64_t)*poly_count + report_step;
    }
    if (ctx->factor_found || ctx->full_count >= target) {
      complete = 1;
      break;
    }
    if (exhausted || *poly_count >= max_polynomials)
      break;
  }
  if (verbose > 3 &&
      (last_report_count != ctx->full_count || last_report_polys != *poly_count))
    siqs_print_relation_report(ctx, target, *poly_count);
  for (i = 0; i < ctx->nthreads; i++)
    psiqs_worker_clear(&workers[i]);
  free(started);
  free(threads);
  free(workers);
  return complete;
}

/* Same partition contract as gmp_siqs; worker count belongs to this call. */
mpz_t *gmp_psiqs(const mpz_t n, uint32_t *nfactors,
                 uint32_t trial_start, uint32_t nthreads) {
  if (nthreads == 0 || nthreads > PSIQS_MAX_THREADS)
    croak("PSIQS: worker count must be between 1 and %u",
          (unsigned)PSIQS_MAX_THREADS);
  return siqs_factor(n, nfactors, trial_start, nthreads);
}
