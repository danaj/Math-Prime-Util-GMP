/* Worker lifecycle suite, included by siqs-check.c, not a compilation unit.
 * Only test source is interposed. Production code and fallback policy remain
 * unchanged. Scheduling gates use native pthread operations, never sleeps. */
#ifdef PSIQS
#ifndef SIQS_CHECK_WORKER_TIMEOUT
# define SIQS_CHECK_WORKER_TIMEOUT 120U
#endif
typedef struct {
  siqs_ctx_t ctx;
  siqs_factor_array_t result;
  siqs_poly_t dispatch;
  mpz_t n, prime;
} worker_fixture_t;

static void worker_resources_clear(void) {
  CHECK(seam_threads_live == 0 && seam_mutexes_live == 0 && seam_conditions_live == 0);
  CHECK(seam_created == seam_joined);
}

static void worker_fault_reset(void) {
  worker_resources_clear();
  seam_alloc_after = seam_mutex_after = seam_cond_after = seam_thread_after = -1;
  seam_thread_alternate = seam_quiet = seam_exit_status = 0;
  seam_thread_calls = seam_created = seam_joined = seam_warnings = 0;
}

static void worker_fixture_open(worker_fixture_t *f, uint32_t threads, int wide) {
  gmp_randstate_t random;
  mpz_t q;
  gmp_randinit_default(random); gmp_randseed_ui(random, 202610036U);
  mpz_init(f->n); mpz_init(f->prime); mpz_init(q);
  do {
    mpz_urandomb(f->prime, random, 77); mpz_setbit(f->prime, 76);
    mpz_nextprime(f->prime, f->prime);
    mpz_urandomb(q, random, 78); mpz_setbit(q, 77); mpz_nextprime(q, q);
    mpz_mul(f->n, f->prime, q);
  } while (mpz_sizeinbase(f->n, 2) != 155 || mpz_cmp(f->prime, q) == 0);
  mpz_clear(q); gmp_randclear(random);
  siqs_factor_array_init(&f->result, f->n);
  siqs_ctx_init(&f->ctx, f->n, f->n, &f->result, NULL);
  if (wide) f->ctx.params.half_interval = 40000;
  f->ctx.nthreads = threads;
  CHECK(siqs_ctx_allocate(&f->ctx));
  CHECK(!f->ctx.factor_found && f->ctx.params.q_count > 2);
  siqs_poly_init(&f->ctx, &f->dispatch);
  case_ctx = &f->ctx;
}

static void worker_fixture_close(worker_fixture_t *f) {
  uint32_t i, count;
  mpz_t product, *values;
  for (i = 0; i < f->ctx.full_count; i++) relation_verify_full(&f->ctx, f->ctx.full[i]);
  CHECK(f->ctx.factor_touched_count == 0);
  mpz_init_set_ui(product, 1);
  for (i = 0; i < f->result.count; i++) mpz_mul(product, product, f->result.values[i]);
  CHECK(mpz_cmp(product, f->n) == 0); mpz_clear(product);
  siqs_poly_clear(&f->ctx, &f->dispatch); siqs_ctx_clear(&f->ctx);
  values = siqs_factor_array_release(&f->result, &count); gmp_siqs_free(values, count);
  mpz_clear(f->n); mpz_clear(f->prime); case_ctx = NULL;
}

static void worker_assert_buffer_cleared(const psiqs_worker_t *worker) {
  const siqs_ctx_t *ctx = &worker->ctx;
  CHECK(ctx->raw_count == 0 && ctx->total_candidates == 0 && ctx->split_attempts == 0);
  CHECK(ctx->split_squfof_or_square == 0 && ctx->split_rho == 0 && ctx->split_failures == 0);
  if (ctx->raw_arena.current != NULL) {
    CHECK(ctx->raw_arena.current->used == 0 && ctx->raw_arena.current->previous == NULL);
  }
}

static void worker_reuse(void) {
  static const uint32_t counts[] = {1, 2, 4};
  uint32_t t, round, i, used;
  for (t = 0; t < sizeof(counts) / sizeof(*counts); t++) {
    worker_fixture_t f;
    psiqs_pool_t *pool;
    mpz_t assigned[32];
    uint32_t *wide_map;
    case_name = "distinct-families-and-wide-map-reuse";
    worker_fault_reset(); worker_fixture_open(&f, counts[t], 1);
    pool = psiqs_pool_create(&f.ctx);
    CHECK(pool != NULL && pool->live == counts[t]);
    CHECK(pool->initialized == counts[t] && pool->conditions == counts[t]);
    for (i = 0; i < counts[t]; i++) CHECK(pool->workers[i].scratch_initialized);
    for (i = 0; i < 32; i++) mpz_init(assigned[i]);
    /* Deliberately migrate one parked worker's map beyond 16-bit entries. */
    memset(pool->workers[0].ctx.sieve, 128, pool->workers[0].ctx.sieve_length);
    for (i = 0; i < 65536; i++) siqs_add_candidate(&pool->workers[0].ctx, i);
    CHECK(pool->workers[0].ctx.candidate_wide);
    wide_map = pool->workers[0].ctx.candidate_at_wide;
    siqs_clear_candidate_map(&pool->workers[0].ctx);
    CHECK(!pool->workers[0].ctx.candidate_wide && wide_map != NULL);
    used = 0;
    for (round = 0; round < (extended ? 6U : 3U); round++) {
      unsigned char seen[4] = {0};
      for (i = 0; i < counts[t]; i++) {
        uint32_t j;
        psiqs_worker_t *w = pool->workers + i;
        CHECK(w->ctx.fb != f.ctx.fb && w->ctx.sieve != f.ctx.sieve);
        CHECK(psiqs_pool_assign(pool, &f.ctx, &f.dispatch, w, 1U + round % 3U));
        CHECK(used < 32);
        mpz_set(assigned[used], w->poly.A);
        for (j = 0; j < used; j++) CHECK(mpz_cmp(assigned[j], assigned[used]) != 0);
        used++;
      }
      for (i = 0; i < counts[t]; i++) {
        psiqs_worker_t *w = psiqs_pool_take(pool);
        uint64_t candidates = f.ctx.total_candidates + w->ctx.total_candidates;
        uint32_t full;
        CHECK(w->index < counts[t] && !seen[w->index]); seen[w->index] = 1;
        CHECK(w->state == PSIQS_IDLE && w->polynomials == 1U + round % 3U);
        psiqs_merge_worker(&f.ctx, w);
        CHECK(f.ctx.total_candidates == candidates); worker_assert_buffer_cleared(w);
        /* A second merge of the now-empty parked buffer must add nothing. */
        full = f.ctx.full_count;
        psiqs_merge_worker(&f.ctx, w);
        CHECK(f.ctx.total_candidates == candidates && f.ctx.full_count == full);
      }
      CHECK(pool->workers[0].ctx.candidate_at_wide == wide_map);
      for (i = 0; i < 65536; i++) CHECK(wide_map[i] == 0);
      CHECK(!f.ctx.factor_found);
    }
    psiqs_pool_destroy(pool); worker_resources_clear();
    for (i = 0; i < 32; i++) mpz_clear(assigned[i]);
    worker_fixture_close(&f);
  }
  puts("PASS workers: 1/2/4-worker reuse, distinct A families, exactly-once merges, wide candidate scratch");
}

typedef struct {
  pthread_mutex_t mutex;
  pthread_cond_t changed;
  psiqs_pool_t *pool;
  unsigned char entered[4], running[4], ready[4];
  int release_start, release_work, joining;
} worker_schedule_t;
static worker_schedule_t worker_schedule;

static void worker_schedule_enter(void *argument) {
  psiqs_worker_t *w = (psiqs_worker_t *)argument;
  worker_schedule_t *s = &worker_schedule;
  pthread_mutex_lock(&s->mutex);
  CHECK(w->index < 4);
  if (s->pool == NULL) s->pool = w->pool;
  CHECK(s->pool == w->pool);
  s->entered[w->index] = 1; pthread_cond_broadcast(&s->changed);
  while (w->index == 2 && !s->release_start) pthread_cond_wait(&s->changed, &s->mutex);
  pthread_mutex_unlock(&s->mutex);
}
static int worker_schedule_before_unlock(pthread_mutex_t *mutex, void *argument) {
  psiqs_worker_t *w = (psiqs_worker_t *)argument;
  worker_schedule_t *s = &worker_schedule;
  int pause = 0;
  /* Caller still owns pool->mutex; read state only under that lock. */
  if (mutex == &w->pool->mutex && w->state == PSIQS_WORK && !w->pool->stop) {
    pthread_mutex_lock(&s->mutex);
    s->running[w->index] = 1; pthread_cond_broadcast(&s->changed);
    pause = w->index == 1;
    pthread_mutex_unlock(&s->mutex);
  }
  return pause;
}
static void worker_schedule_after_unlock(void *argument) {
  worker_schedule_t *s = &worker_schedule;
  (void)argument;
  pthread_mutex_lock(&s->mutex);
  while (!s->release_work) pthread_cond_wait(&s->changed, &s->mutex);
  pthread_mutex_unlock(&s->mutex);
}
static void worker_schedule_signal(pthread_cond_t *cond, void *argument) {
  worker_schedule_t *s = &worker_schedule;
  psiqs_worker_t *w = (psiqs_worker_t *)argument;
  pthread_mutex_lock(&s->mutex);
  if (w != NULL && cond == &w->pool->done) s->ready[w->index] = 1;
  /* join signals worker 0 only after publishing stop under pool->mutex. */
  if (s->pool != NULL && cond == &s->pool->workers[0].work && s->pool->stop) {
    s->joining = s->release_start = s->release_work = 1;
  }
  pthread_cond_broadcast(&s->changed);
  pthread_mutex_unlock(&s->mutex);
}

static void worker_shutdown(void) {
  worker_fixture_t f;
  worker_schedule_t *s = &worker_schedule;
  psiqs_pool_t *pool;
  uint32_t i, families = 0, polynomials = 0;
  case_name = "READY-running-canceled-and-idle-shutdown";
  worker_fault_reset(); worker_fixture_open(&f, 4, 0);
  memset(s, 0, sizeof(*s));
  CHECK(pthread_mutex_init(&s->mutex, NULL) == 0);
  CHECK(pthread_cond_init(&s->changed, NULL) == 0);
  seam_controlled_start = psiqs_pool_worker; seam_enter = worker_schedule_enter;
  seam_unlock_before = worker_schedule_before_unlock;
  seam_unlock_after = worker_schedule_after_unlock; seam_signal = worker_schedule_signal;
  pool = psiqs_pool_create(&f.ctx); CHECK(pool != NULL);
  pthread_mutex_lock(&s->mutex);
  while (!s->entered[0] || !s->entered[1] || !s->entered[2] || !s->entered[3])
    pthread_cond_wait(&s->changed, &s->mutex);
  pthread_mutex_unlock(&s->mutex);
  for (i = 0; i < 3; i++) CHECK(psiqs_pool_assign(pool, &f.ctx, &f.dispatch, pool->workers + i, i == 0 ? 2U : 1U));
  pthread_mutex_lock(&s->mutex);
  while (!s->ready[0] || !s->running[1]) pthread_cond_wait(&s->changed, &s->mutex);
  pthread_mutex_unlock(&s->mutex);
  pthread_mutex_lock(&pool->mutex);
  CHECK(pool->workers[0].state == PSIQS_READY && pool->workers[1].state == PSIQS_WORK);
  CHECK(pool->workers[2].state == PSIQS_WORK && pool->workers[3].state == PSIQS_IDLE);
  CHECK(pool->queued == 1);
  pthread_mutex_unlock(&pool->mutex);
  psiqs_pool_join(pool);
  CHECK(s->joining && pool->live == 0 && seam_threads_live == 0);
  CHECK(pool->workers[0].state == PSIQS_READY && pool->workers[1].state == PSIQS_READY);
  CHECK(pool->workers[2].polynomials == 0 && !s->running[2] && !s->ready[2]);
  CHECK(pool->workers[3].polynomials == 0 && pool->queued == 2);
  for (i = 0; i < pool->count; i++) {
    CHECK(!pool->workers[i].started);
    if (pool->workers[i].state == PSIQS_READY) {
      families++; polynomials += pool->workers[i].polynomials;
      psiqs_merge_worker(&f.ctx, pool->workers + i);
      worker_assert_buffer_cleared(pool->workers + i);
    }
  }
  CHECK(families == 2 && polynomials == 3);
  psiqs_pool_destroy(pool); worker_resources_clear();
  seam_controlled_start = NULL; seam_enter = NULL; seam_unlock_before = NULL;
  seam_unlock_after = NULL; seam_signal = NULL;
  CHECK(pthread_cond_destroy(&s->changed) == 0);
  CHECK(pthread_mutex_destroy(&s->mutex) == 0);
  worker_fixture_close(&f);
  puts("PASS workers: controlled READY/running/not-started/idle shutdown, join-before-drain accounting");
}

static int worker_early_stop_before_unlock(pthread_mutex_t *mutex, void *argument) {
  psiqs_worker_t *w = (psiqs_worker_t *)argument;
  worker_schedule_t *s = &worker_schedule;
  /* Pause after each worker's first stop poll, while its cached flag is false.
   * The coordinator's join publishes stop and releases both gates; the next
   * poll must finish these incomplete families at exactly 64 polynomials. */
  if (mutex == &w->pool->mutex && w->state == PSIQS_WORK &&
      !w->pool->stop && w->polynomials == 32U) {
    pthread_mutex_lock(&s->mutex);
    s->running[w->index] = 1; pthread_cond_broadcast(&s->changed);
    pthread_mutex_unlock(&s->mutex);
    return 1;
  }
  return 0;
}

static void worker_verify_raw(const siqs_ctx_t *ctx) {
  uint32_t i, j;
  mpz_t lhs, rhs, lp;
  mpz_init(lhs); mpz_init(rhs); mpz_init(lp);
  for (i = 0; i < ctx->raw_count; i++) {
    const siqs_raw_relation_t *raw = ctx->raw[i];
    siqs_factor_t *factors = (siqs_factor_t *)check_allocate(raw->nfactors, sizeof(*factors));
    for (j = 0; j < raw->nfactors; j++) {
      factors[j].row = siqs_raw_factor_row(raw, j);
      factors[j].exponent = siqs_raw_factor_exponent(raw, j);
    }
    relation_product(ctx, rhs, factors, raw->nfactors); free(factors);
    siqs_mpz_set_u64(lp, raw->lp1); mpz_mul(rhs, rhs, lp);
    siqs_mpz_set_u64(lp, raw->lp2); mpz_mul(rhs, rhs, lp);
    mpz_mod(rhs, rhs, ctx->n);
    mpz_mul(lhs, raw->y, raw->y); mpz_mod(lhs, lhs, ctx->n);
    CHECK(mpz_cmp(lhs, rhs) == 0);
  }
  mpz_clear(lhs); mpz_clear(rhs); mpz_clear(lp);
}

static void worker_early_shutdown(void) {
  worker_fixture_t f;
  worker_schedule_t *s = &worker_schedule;
  psiqs_pool_t *pool;
  mpz_t abandoned[2];
  uint32_t i, max_polys, families = 0, polys, next = UINT32_MAX;
  uint32_t buffered = 0;
  case_name = "partial-family-stop-and-fresh-family-retry";
  worker_fault_reset(); worker_fixture_open(&f, 2, 0);
  CHECK(f.dispatch.b_limit >= 128U);
  max_polys = f.ctx.params.fb_size * SIQS_MAX_POLYNOMIALS_PER_FB;
  if (max_polys < 1000000U) max_polys = 1000000U;
  polys = max_polys - 196U; /* Reserve 128 + 68 before either worker stops. */
  memset(s, 0, sizeof(*s));
  CHECK(pthread_mutex_init(&s->mutex, NULL) == 0);
  CHECK(pthread_cond_init(&s->changed, NULL) == 0);
  seam_controlled_start = psiqs_pool_worker; seam_enter = worker_schedule_enter;
  seam_unlock_before = worker_early_stop_before_unlock;
  seam_unlock_after = worker_schedule_after_unlock; seam_signal = worker_schedule_signal;
  pool = psiqs_pool_create(&f.ctx); CHECK(pool != NULL && pool->live == 2);
  for (i = 0; i < 2; i++) {
    mpz_init(abandoned[i]);
    CHECK(psiqs_pool_assign(pool, &f.ctx, &f.dispatch, pool->workers + i, i == 0 ? 128U : 68U));
    mpz_set(abandoned[i], pool->workers[i].poly.A);
  }
  pthread_mutex_lock(&s->mutex);
  while (!s->running[0] || !s->running[1]) pthread_cond_wait(&s->changed, &s->mutex);
  pthread_mutex_unlock(&s->mutex);
  psiqs_pool_join(pool);
  CHECK(s->joining && pool->live == 0 && pool->queued == 2);
  for (i = 0; i < 2; i++) {
    psiqs_worker_t *w = pool->workers + i;
    uint64_t candidates = f.ctx.total_candidates + w->ctx.total_candidates;
    uint32_t full;
    CHECK(w->state == PSIQS_READY && w->polynomials == 64U);
    CHECK(w->polynomials < w->limit);
    CHECK(w->ctx.factor_touched_count == 0);
    worker_verify_raw(&w->ctx); buffered += w->ctx.raw_count;
    families++; polys += w->polynomials;
    psiqs_merge_worker(&f.ctx, w); worker_assert_buffer_cleared(w);
    CHECK(f.ctx.total_candidates == candidates);
    full = f.ctx.full_count;
    psiqs_merge_worker(&f.ctx, w);
    CHECK(f.ctx.total_candidates == candidates && f.ctx.full_count == full);
  }
  CHECK(buffered != 0 && families == 2 && max_polys - polys == 68U);
  psiqs_pool_destroy(pool); worker_resources_clear();
  seam_controlled_start = NULL; seam_enter = NULL; seam_unlock_before = NULL;
  seam_unlock_after = NULL; seam_signal = NULL;
  CHECK(pthread_cond_destroy(&s->changed) == 0);
  CHECK(pthread_mutex_destroy(&s->mutex) == 0);
  CHECK(!f.ctx.factor_found);
  /* Model collection after an unsuccessful attempt, retaining the context and
   * partial-family relations.  The unused 68 polynomials are still available;
   * do not replay either abandoned A.  No solver outcome is injected here. */
  CHECK(!psiqs_collect_relations(&f.ctx, &f.dispatch, UINT32_MAX,
                                 &next, &families, &polys));
  CHECK(polys == max_polys && families == 3 && !f.ctx.factor_found);
  for (i = 0; i < 2; i++) {
    CHECK(mpz_cmp(f.dispatch.A, abandoned[i]) != 0); mpz_clear(abandoned[i]);
  }
  worker_resources_clear(); worker_fixture_close(&f);
  puts("PASS workers: stop within 32 polynomials, valid partial buffers, exactly-once merge and fresh-A retry budget");
}

static void worker_merge_ownership(void) {
  siqs_ctx_t master;
  siqs_factor_array_t result;
  psiqs_worker_t a, b;
  siqs_factor_t f = {1, 2}, inverse = {2, 1};
  siqs_raw_relation_t *source;
  siqs_raw_block_t *block;
  uint32_t full;
  case_name = "cross-worker-cycles-and-factor-merge";
  memset(&a, 0, sizeof(a)); memset(&b, 0, sizeof(b));
  relation_context(&master, &result, 2);
  relation_context(&a.ctx, &a.result, 2); relation_context(&b.ctx, &b.result, 2);
  source = relation_raw(&a.ctx, 1, 11, &f, 1); block = a.ctx.raw_arena.current;
  siqs_store_raw(&a.ctx, source);
  a.ctx.total_candidates = 7; a.ctx.split_attempts = 3;
  a.ctx.split_squfof_or_square = 2; a.ctx.split_rho = 1;
  psiqs_merge_worker(&master, &a);
  CHECK(master.raw_count == 1 && master.raw[0] != source);
  CHECK(a.ctx.raw_arena.current == block); worker_assert_buffer_cleared(&a);
  CHECK(master.total_candidates == 7 && master.split_attempts == 3);
  f.row = 2;
  siqs_store_raw(&b.ctx, relation_raw(&b.ctx, 1, 11, &f, 1));
  psiqs_merge_worker(&master, &b);
  CHECK(master.full_count == 1); relation_verify_full(&master, master.full[0]);
  /* Reuse a's retained block, then close a cycle whose LP is a factor. */
  siqs_store_raw(&a.ctx, relation_raw(&a.ctx, 1, 5, &inverse, 1));
  psiqs_merge_worker(&master, &a);
  source = relation_raw(&b.ctx, 1, 5, &inverse, 1);
  mpz_sub(source->y, b.ctx.n, source->y); siqs_store_raw(&b.ctx, source);
  siqs_store_raw(&b.ctx, relation_raw(&b.ctx, 1, 1, &f, 1));
  full = master.full_count;
  psiqs_merge_worker(&master, &b);
  CHECK(master.factor_found && result.count == 2 && master.full_count == full);
  CHECK(siqs_all_factors_prime(&result));
  worker_assert_buffer_cleared(&b);
  relation_finish(&a.ctx, &a.result); relation_finish(&b.ctx, &b.result);
  relation_finish(&master, &result);
  /* Also cover a worker-reported divisor with no emitted relation. */
  relation_context(&master, &result, 2);
  relation_context(&a.ctx, &a.result, 2);
  mpz_set_ui(a.ctx.eval.y, 5);
  CHECK(siqs_insert_divisor(&a.result, a.ctx.eval.y)); a.ctx.factor_found = 1;
  psiqs_merge_worker(&master, &a);
  CHECK(master.factor_found && result.count == 2);
  relation_finish(&a.ctx, &a.result); relation_finish(&master, &result);
  /* Arena reset retains a normal latest block but releases an oversized
   * one-off. There are no GMP records in these storage-only allocations. */
  memset(&a.ctx.raw_arena, 0, sizeof(a.ctx.raw_arena));
  (void)siqs_raw_arena_alloc(&a.ctx.raw_arena, 64);
  (void)siqs_raw_arena_alloc(&a.ctx.raw_arena, SIQS_RAW_BLOCK_INITIAL * 2U);
  block = a.ctx.raw_arena.current; CHECK(block->previous != NULL);
  psiqs_raw_arena_reset(&a.ctx.raw_arena);
  CHECK(a.ctx.raw_arena.current == block && block->used == 0 && block->previous == NULL);
  (void)siqs_raw_arena_alloc(&a.ctx.raw_arena, SIQS_RAW_BLOCK_MAX + 64U);
  psiqs_raw_arena_reset(&a.ctx.raw_arena); CHECK(a.ctx.raw_arena.current == NULL);
  puts("PASS workers: deep-copy ownership, cross-buffer cycles, factor discovery mid-merge and no-relation divisors");
}

static void worker_creation_failures(void) {
  worker_fixture_t f;
  psiqs_pool_t *pool;
  uint32_t i, seen = 0;
  case_name = "pool-creation-failures";
  worker_fault_reset(); worker_fixture_open(&f, 4, 0);
  for (i = 0; i < 3; i++) {
    worker_fault_reset(); seam_alloc_after = (int)i;
    CHECK(psiqs_pool_create(&f.ctx) == NULL); worker_resources_clear();
  }
  worker_fault_reset(); seam_mutex_after = 0;
  CHECK(psiqs_pool_create(&f.ctx) == NULL); worker_resources_clear();
  for (i = 0; i <= 4; i++) {
    worker_fault_reset(); seam_cond_after = (int)i;
    CHECK(psiqs_pool_create(&f.ctx) == NULL); worker_resources_clear();
  }
  worker_fault_reset(); seam_thread_after = 0; seam_quiet = 1;
  CHECK(psiqs_pool_create(&f.ctx) == NULL);
  CHECK(seam_thread_calls == 4 && seam_created == 0 && seam_warnings == 4);
  worker_resources_clear();
  worker_fault_reset(); seam_thread_alternate = seam_quiet = 1;
  pool = psiqs_pool_create(&f.ctx);
  CHECK(pool != NULL && pool->live == 2 && pool->initialized == 4);
  CHECK(pool->workers[0].started && !pool->workers[1].started);
  CHECK(pool->workers[2].started && !pool->workers[3].started);
  CHECK(seam_warnings == 2);
  for (i = 0; i < pool->count; i++) {
    psiqs_worker_t *w = pool->workers + i;
    CHECK(w->scratch_initialized == w->started);
    if (w->started) CHECK(w->ctx.fb != NULL && w->ctx.sieve != NULL);
    else {
      CHECK(w->ctx.fb == NULL && w->ctx.sieve == NULL && w->poly.a_index == NULL);
      CHECK(w->result.values == NULL && w->result.primality == NULL && w->result.count == 0);
      psiqs_worker_clear(w); /* Re-clearing an unowned slot must be harmless. */
      CHECK(!w->scratch_initialized);
    }
  }
  for (i = 0; i < pool->count; i++) if (pool->workers[i].started)
    CHECK(psiqs_pool_assign(pool, &f.ctx, &f.dispatch, pool->workers + i, 1));
  for (i = 0; i < 2; i++) {
    psiqs_worker_t *w = psiqs_pool_take(pool);
    CHECK(w->index == 0 || w->index == 2);
    CHECK(!(seen & (1U << w->index))); seen |= 1U << w->index;
    psiqs_merge_worker(&f.ctx, w);
  }
  CHECK(seen == 5U);
  CHECK(pool->workers[1].polynomials == 0 && pool->workers[3].polynomials == 0);
  psiqs_pool_join(pool);
  for (i = 0; i < pool->count; i++) {
    CHECK(!pool->workers[i].started);
    /* Joining changes thread lifetime, not scratch ownership. */
    CHECK(pool->workers[i].scratch_initialized == (i == 0 || i == 2));
  }
  psiqs_pool_destroy(pool); worker_resources_clear(); worker_fault_reset();
  if (extended) {
    f.ctx.nthreads = PSIQS_MAX_THREADS;
    seam_thread_after = 2; seam_quiet = 1;
    pool = psiqs_pool_create(&f.ctx);
    CHECK(pool != NULL && pool->count == PSIQS_MAX_THREADS && pool->live == 2);
    CHECK(seam_thread_calls == PSIQS_MAX_THREADS && seam_created == 2);
    CHECK(seam_warnings == PSIQS_MAX_THREADS - 2U);
    for (i = 0; i < pool->count; i++) {
      psiqs_worker_t *w = pool->workers + i;
      CHECK(w->scratch_initialized == (i < 2));
      CHECK((w->ctx.fb != NULL) == (i < 2));
      CHECK((w->result.values != NULL) == (i < 2));
    }
    psiqs_pool_destroy(pool); worker_resources_clear(); worker_fault_reset();
  }
  worker_fixture_close(&f);
  puts("PASS workers: allocation/init failures, zero-live fallback, noncontiguous reduced pools and cleanup");
}

static void worker_budget_and_fallback(void) {
  uint32_t test;
  for (test = 0; test < 2; test++) {
    worker_fixture_t f;
    uint32_t max_polys, families = 0, polys, next = UINT32_MAX;
    case_name = "shortened-budget-and-serial-fallback";
    worker_fault_reset(); worker_fixture_open(&f, 4, 0);
    max_polys = f.ctx.params.fb_size * SIQS_MAX_POLYNOMIALS_PER_FB;
    if (max_polys < 1000000U) max_polys = 1000000U;
    polys = max_polys - 5U;
    if (test == 0) f.dispatch.b_limit = 3; /* Reserve 3, then the last 2. */
    else { seam_thread_after = 0; seam_quiet = 1; }
    CHECK(!psiqs_collect_relations(&f.ctx, &f.dispatch, UINT32_MAX,
                                   &next, &families, &polys));
    CHECK(polys == max_polys && families == (test == 0 ? 2U : 1U));
    CHECK(!f.ctx.factor_found);
    if (test == 1) CHECK(seam_created == 0 && seam_warnings == 5);
    worker_resources_clear(); worker_fault_reset(); worker_fixture_close(&f);
  }
  puts("PASS workers: shortened last-family budget, executed-work accounting, total-create serial fallback");
}

/* Reuse the earlier Lanczos resource-failure cases, using the existing matrix
 * fixture generator and deterministic output checks rather than another tool. */
static void worker_lanczos_failures(void) {
  la_col_t *cols;
  nla_matrix_t matrix;
  uint32_t i;
  uint64_t mask, other_mask, *serial, *fallback;
  case_name = "Lanczos-pool-failure-serial-fallback";
  worker_fault_reset(); cols = matrix_fixture(129, 32769, 37);
  nla_matrix_init(&matrix, 129, 37, 32769, cols, 48);
  CHECK(!matrix.packed);
  for (i = 0; i < 6; i++) {
    worker_fault_reset(); seam_alloc_after = (int)i;
    CHECK(nla_pool_create(&matrix, 4) == NULL); worker_resources_clear();
  }
  worker_fault_reset(); seam_mutex_after = 0;
  CHECK(nla_pool_create(&matrix, 4) == NULL); worker_resources_clear();
  for (i = 0; i < 2; i++) {
    worker_fault_reset(); seam_cond_after = (int)i;
    CHECK(nla_pool_create(&matrix, 4) == NULL); worker_resources_clear();
  }
  for (i = 0; i < 3; i++) {
    worker_fault_reset(); seam_thread_after = (int)i;
    CHECK(nla_pool_create(&matrix, 4) == NULL);
    CHECK(seam_created == i); worker_resources_clear();
  }
  worker_fault_reset();
  CHECK(nla_pool_create(&matrix, 0) == NULL && nla_pool_create(&matrix, 1) == NULL);
  serial = la_block_lanczos(129, 37, 32769, cols, 31, 47, &mask);
  matrix_verify_dependencies(129, 32769, 37, cols, serial, mask);
  seam_thread_after = 1;
  fallback = la_block_lanczos_threaded(129, 37, 32769, cols, 31, 47, &other_mask, 4, 0);
  CHECK(mask == other_mask && memcmp(serial, fallback, 32769U * sizeof(uint64_t)) == 0);
  matrix_verify_dependencies(129, 32769, 37, cols, fallback, other_mask);
  free(serial); free(fallback); nla_matrix_clear(&matrix); matrix_free(cols, 32769);
  worker_resources_clear(); worker_fault_reset();
  puts("PASS workers: retained Lanczos allocation/init/partial-create failures and fixed-seed serial fallback");
}

#ifndef _WIN32
static void worker_timeout(int signal_number) {
  static const char message[] = "FAIL workers: watchdog expired (possible hang)\n";
  ssize_t written;
  (void)signal_number;
  written = write(STDERR_FILENO, message, sizeof(message) - 1U);
  (void)written; /* Best-effort diagnostic; exit even if stderr is unwritable. */
  _exit(1);
}
static void worker_fatal_cases(void) {
  worker_fixture_t f;
  uint32_t test;
  case_name = "invalid-counts-and-fatal-scratch-OOM";
  worker_fault_reset(); worker_fixture_open(&f, 2, 0);
  for (test = 0; test < 3; test++) {
    int errpipe[2], outpipe[2], status;
    char diagnostic[512], output[8];
    ssize_t bytes;
    pid_t child;
    CHECK(pipe(errpipe) == 0 && pipe(outpipe) == 0); fflush(NULL);
    child = fork(); CHECK(child >= 0);
    if (child == 0) {
      uint32_t count;
      close(errpipe[0]); close(outpipe[0]);
      if (dup2(errpipe[1], STDERR_FILENO) < 0 || dup2(outpipe[1], STDOUT_FILENO) < 0) _exit(1);
      close(errpipe[1]); close(outpipe[1]);
      alarm(20); seam_exit_status = 3;
      if (test < 2) (void)gmp_psiqs(f.n, &count, 0, test == 0 ? 0 : PSIQS_MAX_THREADS + 1U);
      else { seam_alloc_after = 3; (void)psiqs_pool_create(&f.ctx); }
      _exit(1); /* Returning instead of croaking is a failed check. */
    }
    close(errpipe[1]); close(outpipe[1]);
    CHECK(waitpid(child, &status, 0) == child);
    bytes = read(errpipe[0], diagnostic, sizeof(diagnostic) - 1U);
    CHECK(bytes > 0); diagnostic[bytes] = '\0';
    CHECK(read(outpipe[0], output, sizeof(output)) == 0);
    CHECK(WIFEXITED(status) && WEXITSTATUS(status) == 3);
    CHECK(strstr(diagnostic, test < 2 ? "worker count must be between" : "unable to allocate memory") != NULL);
    CHECK(diagnostic[bytes - 1] == '\n'); close(errpipe[0]); close(outpipe[0]);
  }
  worker_fixture_close(&f);
  puts("PASS workers: invalid 0/257 requests and fatal worker-scratch OOM, stderr/status verified in children");
}
#endif

static void worker_public_calls(void) {
  worker_fixture_t f;
  uint32_t threads;
  case_name = "repeated-public-factor-calls";
  worker_fault_reset(); worker_fixture_open(&f, 2, 0);
  for (threads = 1; threads <= 4; threads *= 2) {
    uint32_t count, i;
    mpz_t *values, product;
    values = gmp_psiqs(f.n, &count, 0, threads); CHECK(values != NULL && count == 2);
    mpz_init_set_ui(product, 1);
    for (i = 0; i < count; i++) {
      CHECK(mpz_probab_prime_p(values[i], 25)); mpz_mul(product, product, values[i]);
    }
    CHECK(mpz_cmp(product, f.n) == 0);
    mpz_clear(product); gmp_siqs_free(values, count); worker_resources_clear();
  }
  worker_fixture_close(&f);
  puts("PASS workers: repeated public 1/2/4-thread calls return exact prime partitions");
}

typedef struct {
  mpz_srcptr n;
  mpz_t *values;
  uint32_t count;
} worker_public_call_t;
static void *worker_public_call(void *argument) {
  worker_public_call_t *call = (worker_public_call_t *)argument;
  call->values = gmp_psiqs(call->n, &call->count, 0, 2);
  return NULL;
}
static void worker_simultaneous_calls(void) {
  worker_fixture_t f;
  worker_public_call_t calls[2];
  pthread_t threads[2];
  uint32_t i, j;
  mpz_t product;
  case_name = "simultaneous-independent-caller-pools";
  worker_fault_reset(); worker_fixture_open(&f, 2, 0);
  memset(calls, 0, sizeof(calls));
  for (i = 0; i < 2; i++) {
    calls[i].n = f.n;
    CHECK(pthread_create(threads + i, NULL, worker_public_call, calls + i) == 0);
  }
  for (i = 0; i < 2; i++) CHECK(pthread_join(threads[i], NULL) == 0);
  mpz_init(product);
  for (i = 0; i < 2; i++) {
    CHECK(calls[i].values != NULL && calls[i].count == 2);
    mpz_set_ui(product, 1);
    for (j = 0; j < calls[i].count; j++) {
      CHECK(mpz_probab_prime_p(calls[i].values[j], 25));
      mpz_mul(product, product, calls[i].values[j]);
    }
    CHECK(mpz_cmp(product, f.n) == 0);
    gmp_siqs_free(calls[i].values, calls[i].count);
  }
  mpz_clear(product); worker_resources_clear(); worker_fixture_close(&f);
  puts("PASS workers: two simultaneous callers own independent pools after one host startup");
}

static void suite_workers(void) {
  worker_fault_reset();
  CHECK(pthread_key_create(&seam_key, NULL) == 0); seam_key_ready = 1;
#ifndef _WIN32
  CHECK(signal(SIGALRM, worker_timeout) != SIG_ERR); alarm(SIQS_CHECK_WORKER_TIMEOUT);
#endif
  prime_iterator_global_startup();
  worker_reuse(); worker_shutdown(); worker_early_shutdown();
  worker_merge_ownership(); worker_creation_failures();
  worker_budget_and_fallback(); worker_lanczos_failures(); worker_public_calls();
  worker_simultaneous_calls();
#ifndef _WIN32
  worker_fatal_cases(); alarm(0); CHECK(signal(SIGALRM, SIG_DFL) != SIG_ERR);
#else
  puts("SKIP workers: child-process fatal/watchdog checks need POSIX fork/alarm");
#endif
  prime_iterator_global_shutdown();
  worker_fault_reset(); CHECK(pthread_key_delete(seam_key) == 0); seam_key_ready = 0;
  puts("PASS workers: polynomial pool lifecycle/fault coverage and Lanczos fallback"); fflush(stdout);
}
#else
static void suite_workers(void) {
  puts("SKIP workers: use make check-psiqs or --threaded");
}
#endif
