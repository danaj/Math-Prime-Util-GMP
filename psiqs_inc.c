/*============================================================================
  Parallel family collection for SIQS.

  Included only by siqs.c when PSIQS is defined.  Do not compile this file
  separately: it uses the shared private SIQS types and static helpers.

  Workers own polynomial/sieve scratch and buffer raw relations.  The caller
  assigns distinct A families and merges each completed buffer into one relation
  graph; factor partition updates and relation merging remain serial.

  Copyright (c) 2026 Dana Jacobsen
============================================================================*/

/* Native worker calls into Perl's error/primality adapters need a separate
 * audit.  For now, enable this implementation only in the standalone host. */
#ifndef STANDALONE
# error "PSIQS currently requires the standalone host"
#endif

#include <pthread.h>

typedef struct psiqs_pool_t psiqs_pool_t;

typedef enum { PSIQS_IDLE, PSIQS_WORK, PSIQS_READY } psiqs_worker_state_t;

typedef struct {
  siqs_ctx_t ctx;
  siqs_poly_t poly;
  siqs_factor_array_t result;
  uint32_t limit;
  uint32_t polynomials;
  psiqs_pool_t *pool;
  pthread_t thread;
  pthread_cond_t work;
  uint32_t index;
  psiqs_worker_state_t state;
  int started;
  int scratch_initialized; /* Independent of thread lifetime and condition. */
} psiqs_worker_t;

struct psiqs_pool_t {
  psiqs_worker_t *workers;
  uint32_t *ready;
  uint32_t count, initialized, conditions, live;
  uint32_t head, tail, queued;
  pthread_mutex_t mutex;
  pthread_cond_t done;
  int stop;
};

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
  worker->scratch_initialized = 1;
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
    /* A q=12 family has 2048 polynomials.  At readiness, finish only
     * a short prefix instead of delaying join for the rest of that family.
     * Keep every buffered relation; retries select fresh A values.  Read the
     * existing stop flag under its mutex, outside the sieve/root kernels. */
    if ((worker->polynomials & 31U) == 0) {
      int stop;
      pthread_mutex_lock(&worker->pool->mutex);
      stop = worker->pool->stop;
      pthread_mutex_unlock(&worker->pool->mutex);
      if (stop)
        break;
    }
  } while (siqs_next_B(ctx, poly));
  return NULL;
}

/* After every buffered GMP value is cleared, retain one normal-sized block
 * for the next family.  Older blocks and oversized one-offs are released;
 * final worker cleanup still frees the retained block. */
static void psiqs_raw_arena_reset(siqs_raw_arena_t *arena) {
  siqs_raw_block_t *block = arena->current;
  if (block != NULL &&
      block->capacity <= SIQS_RAW_BLOCK_MAX - sizeof(*block)) {
    arena->current = block->previous;
    siqs_raw_arena_clear(arena);
    block->previous = NULL;
    block->used = 0;
    arena->current = block;
  } else {
    siqs_raw_arena_clear(arena);
  }
}

/* Copy buffered records into the collector's arena; never transfer pointers
 * between arenas.  All partials enter the same graph, including cross-worker
 * cycles.  Clear worker GMP values before resetting its buffer storage. */
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
   * Refine the parent's partition only here, while this worker is parked. */
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
  psiqs_raw_arena_reset(&ctx->raw_arena);
}

static void psiqs_worker_clear(psiqs_worker_t *worker) {
  if (!worker->scratch_initialized)
    return;
  siqs_poly_clear(&worker->ctx, &worker->poly);
  siqs_ctx_clear(&worker->ctx);
  free(worker->result.primality);
  gmp_siqs_free(worker->result.values, worker->result.count);
  memset(&worker->result, 0, sizeof(worker->result));
  worker->scratch_initialized = 0;
}

/* Each worker has one result buffer.  Publishing it parks that worker until
 * the caller finishes merging and assigns another family.  Other workers
 * continue independently; the queue never holds more than one result each. */
static void *psiqs_pool_worker(void *argument) {
  psiqs_worker_t *worker = (psiqs_worker_t *)argument;
  psiqs_pool_t *pool = worker->pool;
  pthread_mutex_lock(&pool->mutex);
  for (;;) {
    while (!pool->stop && worker->state != PSIQS_WORK)
      pthread_cond_wait(&worker->work, &pool->mutex);
    if (pool->stop)
      break;
    pthread_mutex_unlock(&pool->mutex);
    psiqs_sieve_family(worker);
    pthread_mutex_lock(&pool->mutex);
    worker->state = PSIQS_READY;
    if (pool->queued >= pool->count)
      croak("PSIQS: completed-family queue overflow");
    pool->ready[pool->tail] = worker->index;
    if (++pool->tail == pool->count)
      pool->tail = 0;
    pool->queued++;
    pthread_cond_signal(&pool->done);
  }
  pthread_mutex_unlock(&pool->mutex);
  return NULL;
}

/* Stop assigning work, finish short polynomial prefixes of running families,
 * and cancel jobs not yet begun.  Unrun family tails are not resumed on retry.
 * The caller drains completed buffers after joining, before freeing scratch. */
static void psiqs_pool_join(psiqs_pool_t *pool) {
  uint32_t i;
  if (pool->live == 0)
    return;
  pthread_mutex_lock(&pool->mutex);
  pool->stop = 1;
  for (i = 0; i < pool->conditions; i++)
    if (pool->workers[i].started)
      pthread_cond_signal(&pool->workers[i].work);
  pthread_mutex_unlock(&pool->mutex);
  for (i = 0; i < pool->count; i++) {
    if (pool->workers[i].started) {
      int error = pthread_join(pool->workers[i].thread, NULL);
      if (error != 0)
        croak("PSIQS: pthread_join failed: %s", strerror(error));
      pool->workers[i].started = 0;
    }
  }
  pool->live = 0;
}

static void psiqs_pool_destroy(psiqs_pool_t *pool) {
  uint32_t i;
  if (pool == NULL)
    return;
  psiqs_pool_join(pool);
  for (i = 0; i < pool->conditions; i++)
    pthread_cond_destroy(&pool->workers[i].work);
  pthread_cond_destroy(&pool->done);
  pthread_mutex_destroy(&pool->mutex);
  for (i = 0; i < pool->initialized; i++)
    psiqs_worker_clear(pool->workers + i);
  free(pool->ready);
  free(pool->workers);
  free(pool);
}

/* Creation failures reduce the pool; no jobs run until initialization ends.
 * If no worker can start, the caller uses the ordinary serial collector. */
static psiqs_pool_t *psiqs_pool_create(siqs_ctx_t *ctx) {
  psiqs_pool_t *pool = (psiqs_pool_t *)calloc(1, sizeof(*pool));
  uint32_t i;
  if (pool == NULL)
    return NULL;
  pool->count = ctx->nthreads;
  pool->workers = (psiqs_worker_t *)calloc(pool->count, sizeof(*pool->workers));
  pool->ready = (uint32_t *)malloc((size_t)pool->count * sizeof(*pool->ready));
  if (pool->workers == NULL || pool->ready == NULL)
    goto fail_mutex;
  if (pthread_mutex_init(&pool->mutex, NULL) != 0)
    goto fail_mutex;
  if (pthread_cond_init(&pool->done, NULL) != 0)
    goto fail_done;
  for (i = 0; i < pool->count; i++) {
    psiqs_worker_t *worker = pool->workers + i;
    psiqs_worker_init(worker, ctx, i);
    pool->initialized++;
    worker->pool = pool;
    worker->index = i;
    if (pthread_cond_init(&worker->work, NULL) != 0) {
      psiqs_pool_destroy(pool);
      return NULL;
    }
    pool->conditions++;
  }
  for (i = 0; i < pool->count; i++) {
    psiqs_worker_t *worker = pool->workers + i;
    int error = pthread_create(&worker->thread, NULL, psiqs_pool_worker, worker);
    if (error == 0) {
      worker->started = 1;
      pool->live++;
    } else {
      fprintf(stderr, "PSIQS: pthread_create: %s; reducing worker pool\n",
              strerror(error));
      /* No thread owns this slot. Keep its condition for pool teardown, but
       * release scratch now; the ownership flag prevents a second clear. */
      psiqs_worker_clear(worker);
    }
  }
  if (pool->live == 0) {
    psiqs_pool_destroy(pool);
    return NULL;
  }
  return pool;

fail_done:
  pthread_mutex_destroy(&pool->mutex);
fail_mutex:
  free(pool->ready);
  free(pool->workers);
  free(pool);
  return NULL;
}

/* Only the caller assigns families, preserving the single-owner A hash/RNG. */
static int psiqs_pool_assign(psiqs_pool_t *pool, siqs_ctx_t *ctx,
                              siqs_poly_t *dispatch, psiqs_worker_t *worker,
                              uint32_t limit) {
  if (!psiqs_assign_family(ctx, dispatch, worker, limit))
    return 0;
  pthread_mutex_lock(&pool->mutex);
  worker->state = PSIQS_WORK;
  pthread_cond_signal(&worker->work);
  pthread_mutex_unlock(&pool->mutex);
  return 1;
}

/* Taking a result acquires every write made by that worker before publication. */
static psiqs_worker_t *psiqs_pool_take(psiqs_pool_t *pool) {
  psiqs_worker_t *worker;
  pthread_mutex_lock(&pool->mutex);
  while (pool->queued == 0)
    pthread_cond_wait(&pool->done, &pool->mutex);
  worker = pool->workers + pool->ready[pool->head];
  if (++pool->head == pool->count)
    pool->head = 0;
  pool->queued--;
  worker->state = PSIQS_IDLE;
  pthread_mutex_unlock(&pool->mutex);
  return worker;
}

/* A bounded, asynchronous A-family pool: no per-polynomial locks or graph
 * sharing.  Poll for stop every 32 polynomials under the existing mutex.
 * Reserve polynomial budgets on assignment, refund unused work on
 * completion, and join every thread before the matrix solver can run. */
static int psiqs_collect_relations(siqs_ctx_t *ctx, siqs_poly_t *dispatch,
                                   uint32_t target,
                                   uint32_t *next_matrix_check,
                                   uint32_t *family_count,
                                   uint32_t *poly_count) {
  psiqs_pool_t *pool;
  uint32_t i, check_interval = ctx->params.fb_size / 128;
  uint32_t max_polynomials =
      ctx->params.fb_size > UINT32_MAX / SIQS_MAX_POLYNOMIALS_PER_FB
        ? UINT32_MAX : ctx->params.fb_size * SIQS_MAX_POLYNOMIALS_PER_FB;
  uint64_t report_step =
      (SIQS_PROGRESS_WORK_INTERVAL + ctx->params.half_interval - 1U) /
      ctx->params.half_interval;
  uint64_t next_report = (uint64_t)*poly_count + report_step;
  uint32_t last_report_count = UINT32_MAX, last_report_polys = UINT32_MAX;
  uint32_t remaining, active = 0;
  int verbose = siqs_verbose_level(), complete = 0, exhausted = 0, ready = 0;
  if (max_polynomials < 1000000U)
    max_polynomials = 1000000U;
  if (check_interval < SIQS_MATRIX_CHECK_MIN)
    check_interval = SIQS_MATRIX_CHECK_MIN;
  if (check_interval > SIQS_MATRIX_CHECK_MAX)
    check_interval = SIQS_MATRIX_CHECK_MAX;
  if (ctx->factor_found || ctx->full_count >= target)
    return 1;
  if (*poly_count >= max_polynomials)
    return 0;
  pool = psiqs_pool_create(ctx);
  if (pool == NULL) {
    fprintf(stderr, "PSIQS: worker pool unavailable; using serial collection\n");
    return siqs_collect_relations(ctx, dispatch, target, next_matrix_check,
                                  family_count, poly_count);
  }
  remaining = max_polynomials - *poly_count;
  for (i = 0; i < pool->count && remaining != 0; i++) {
    if (pool->workers[i].started) {
      uint32_t limit = dispatch->b_limit < remaining
                     ? dispatch->b_limit : remaining;
      if (!psiqs_pool_assign(pool, ctx, dispatch, pool->workers + i, limit)) {
        exhausted = 1;
        break;
      }
      remaining -= limit;
      active++;
    }
  }
  while (active != 0) {
    psiqs_worker_t *worker = psiqs_pool_take(pool);
    active--;
    remaining += worker->limit - worker->polynomials;
    (*family_count)++;
    *poly_count += worker->polynomials;
    psiqs_merge_worker(ctx, worker);
    if (ctx->full_count >= *next_matrix_check && !ctx->factor_found) {
      uint32_t rows, columns;
      ready = siqs_matrix_ready(ctx, &rows, &columns);
      *next_matrix_check = ctx->full_count + check_interval;
      if (!ready && verbose > 4) {
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
    if (!exhausted && remaining != 0) {
      uint32_t limit = dispatch->b_limit < remaining
                     ? dispatch->b_limit : remaining;
      if (psiqs_pool_assign(pool, ctx, dispatch, worker, limit)) {
        remaining -= limit;
        active++;
      } else {
        exhausted = 1;
      }
    }
  }
  psiqs_pool_join(pool);
  /* No worker can still mutate these buffers.  Include every family that ran,
   * even if readiness was detected before its completion; canceled jobs have
   * no results.  Total executed polynomials cannot exceed reserved budgets. */
  for (i = 0; i < pool->count; i++) {
    psiqs_worker_t *worker = pool->workers + i;
    if (worker->state == PSIQS_READY) {
      (*family_count)++;
      *poly_count += worker->polynomials;
      psiqs_merge_worker(ctx, worker);
    }
  }
  if (ctx->factor_found || ctx->full_count >= target)
    complete = 1;
  if (ready && verbose > 3) {
    uint32_t rows, columns;
    (void)siqs_matrix_ready(ctx, &rows, &columns);
    siqs_print_relation_report(ctx, target, *poly_count);
    last_report_count = ctx->full_count;
    last_report_polys = *poly_count;
    printf("# siqs matrix core %u columns, %u rows, ready\n", columns, rows);
    fflush(stdout);
  }
  if (verbose > 3 &&
      (last_report_count != ctx->full_count || last_report_polys != *poly_count))
    siqs_print_relation_report(ctx, target, *poly_count);
  psiqs_pool_destroy(pool);
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
