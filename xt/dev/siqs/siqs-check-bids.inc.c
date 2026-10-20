/* BIDS suite, included by siqs-check.c. Valid congruences exercise the
 * scalar/packed bases, continued partition refinement and collector stops. */
static void bids_fixture_open(siqs_ctx_t *ctx, siqs_factor_array_t *result,
                               uint32_t n, uint32_t fb) {
  memset(ctx, 0, sizeof(*ctx));
  mpz_init_set_ui(ctx->n, n); mpz_init_set(ctx->kn, ctx->n);
  mpz_init(ctx->eval.y); mpz_init(ctx->eval.q); mpz_init(ctx->eval.rest);
  ctx->original_n = ctx->n;
  ctx->params.fb_size = fb;
  ctx->params.max_large_primes = 1;
  ctx->fb = (siqs_fb_t *)check_allocate(fb, sizeof(*ctx->fb));
  ctx->factor_counts = (uint32_t *)check_allocate(fb + 1U, sizeof(uint32_t));
  ctx->factor_touched = (uint32_t *)check_allocate(fb + 1U, sizeof(uint32_t));
  siqs_factor_array_init(result, ctx->n); ctx->result = result;
  ctx->use_bids = 1;
  ctx->bids_cap = fb + 1U + SIQS_BIDS_EXTRA_RELS;
  case_ctx = ctx;
}

static void bids_fixture_add(siqs_ctx_t *ctx, uint32_t y,
                              const siqs_factor_t *factors, uint32_t count) {
  siqs_full_relation_t *r = (siqs_full_relation_t *)check_allocate(1, sizeof(*r));
  mpz_init_set_ui(r->y, y);
  r->nfactors = count;
  r->factors = (siqs_factor_t *)check_allocate(count, sizeof(*r->factors));
  if (count) memcpy(r->factors, factors, (size_t)count * sizeof(*factors));
  relation_verify_full(ctx, r);
  siqs_store_full(ctx, r);
}

static void bids_fixture_close(siqs_ctx_t *ctx, siqs_factor_array_t *result) {
  uint32_t i, count;
  mpz_t product, *values;
  CHECK(ctx->factor_touched_count == 0);
  for (i = 0; i <= ctx->params.fb_size; i++) CHECK(ctx->factor_counts[i] == 0);
  mpz_init_set_ui(product, 1);
  for (i = 0; i < result->count; i++) mpz_mul(product, product, result->values[i]);
  CHECK(mpz_cmp(product, ctx->n) == 0);
  mpz_clear(product);
  siqs_ctx_clear(ctx);
  values = siqs_factor_array_release(result, &count); gmp_siqs_free(values, count);
  case_ctx = NULL;
}

/* Labelled rows all use the prime 131 == 1 (mod 65); the sign root is 8.
 * This gives valid congruences with any desired parity. It is an artificial
 * matrix fixture, not a real factor base. The column transform is triangular
 * in a shuffled row order, so rank is known independently of elimination. */
static void bids_layout(uint32_t rows) {
  siqs_ctx_t ctx;
  siqs_factor_array_t result;
  siqs_factor_t *factors = (siqs_factor_t *)check_allocate(rows, sizeof(*factors));
  uint32_t *order = (uint32_t *)check_allocate(rows, sizeof(*order));
  uint32_t i;
  case_name = "shuffled-pivots-and-provenance-prefixes";
  bids_fixture_open(&ctx, &result, 65, rows - 1U);
  for (i = 0; i < rows; i++) {
    order[i] = i;
    factors[i].row = i; factors[i].exponent = 1;
    if (i) ctx.fb[i - 1U].p = 131;
  }
  for (i = rows; i > 1U; i--) {
    uint32_t j = matrix_random() % i, swap = order[i - 1U];
    order[i - 1U] = order[j]; order[j] = swap;
  }
  for (i = 0; i < rows; i++) {
    siqs_factor_t f[2] = {{order[i], 1}, {i ? order[i - 1U] : 0, 1}};
    uint32_t count = i ? 2 : 1;
    int sign = order[i] == 0 || (i && order[i - 1U] == 0);
    if (i && f[0].row > f[1].row) {
      siqs_factor_t swap = f[0]; f[0] = f[1]; f[1] = swap;
    }
    bids_fixture_add(&ctx, sign ? 8 : 1, f, count);
    CHECK(ctx.bids->rank == i + 1U && !ctx.bids->failed);
  }
  CHECK(ctx.bids->dependencies == 0 && !ctx.bids_stop && !ctx.factor_found);
  CHECK((ctx.bids->words == 1U) == (rows <= 64U));
  CHECK((ctx.bids->storage == NULL) == (rows <= 64U));
  bids_fixture_add(&ctx, 8, factors, rows);
  bids_fixture_add(&ctx, order[0] == 0 ? 8 : 1, factors + order[0], 1);
  CHECK(ctx.bids->rank == rows && ctx.bids->processed == rows + 2U);
  CHECK(ctx.bids->dependencies == 2 && ctx.bids->trivial == 2);
  CHECK(!ctx.factor_found && !ctx.bids->failed && !ctx.bids_stop);
  bids_fixture_close(&ctx, &result); free(factors); free(order);
}

static void bids_soft_budget(uint32_t fb) {
  siqs_ctx_t ctx;
  siqs_factor_array_t result;
  siqs_factor_t f = {1, 1};
  uint32_t i, cap;
  case_name = "partial-split-soft-stop-and-late-completion";
  bids_fixture_open(&ctx, &result, 105, fb);
  ctx.fb[0].p = 211; /* 1 modulo 105. */
  cap = ctx.bids_cap;
  bids_fixture_add(&ctx, 1, &f, 1);
  bids_fixture_add(&ctx, 41, &f, 1); /* 5 * 21: continue with the existing basis. */
  CHECK(ctx.bids->rank == 1 && ctx.bids->progress && ctx.bids->partial == 1);
  CHECK(!ctx.factor_found && !siqs_all_factors_prime(&result));
  for (i = ctx.full_count; i < cap; i++) {
    CHECK(!ctx.bids_stop); bids_fixture_add(&ctx, 1, NULL, 0);
  }
  CHECK(ctx.bids_stop == SIQS_BIDS_BUDGET && !ctx.factor_found);
  /* A new candidate list must not start after the stop request. */
  ctx.candidate_count = 1; siqs_evaluate_candidates(&ctx, NULL);
  CHECK(ctx.full_count == cap);
  for (i = 0; i < 10; i++) bids_fixture_add(&ctx, 1, NULL, 0);
  CHECK(ctx.full_count == cap + 10U && !ctx.factor_found);
  bids_fixture_add(&ctx, 34, &f, 1); /* A queued relation splits 21 into 3 * 7. */
  CHECK(ctx.factor_found && siqs_all_factors_prime(&result));
  CHECK(ctx.bids->rank == 1 && ctx.bids->original[0] == 0);
  CHECK(ctx.bids->processed == ctx.full_count && ctx.full_count == cap + 11U);
  bids_fixture_close(&ctx, &result);

  case_name = "nondependency-partial-split";
  bids_fixture_open(&ctx, &result, 105, fb);
  siqs_one_lp_init(&ctx.one_lp);
  CHECK(!siqs_materialize_one_lp_pair(&ctx, NULL, NULL, 5));
  CHECK(!ctx.factor_found && ctx.bids->progress && ctx.bids->partial == 1);
  CHECK(ctx.full_count == 0 && !siqs_all_factors_prime(&result));
  bids_fixture_add(&ctx, 34, NULL, 0);
  CHECK(ctx.factor_found && siqs_all_factors_prime(&result));
  bids_fixture_close(&ctx, &result);
}

static void bids_uncertified(uint32_t fb) {
  siqs_ctx_t ctx;
  siqs_factor_array_t result;
  siqs_factor_t f = {1, UINT32_MAX};
  case_name = "overflow-does-not-certify-a-dependency";
  bids_fixture_open(&ctx, &result, 65, fb);
  ctx.fb[0].p = 131;
  bids_fixture_add(&ctx, 1, &f, 1);
  f.exponent = 1; bids_fixture_add(&ctx, 1, &f, 1);
  CHECK(ctx.bids->failed && ctx.bids_stop == SIQS_BIDS_UNCERTIFIED);
  CHECK(ctx.bids->unknown == 1 && !ctx.factor_found);
  /* Pending records remain owned even after the solver becomes unusable. */
  bids_fixture_add(&ctx, 1, NULL, 0);
  CHECK(ctx.full_count == 3 && ctx.bids->processed == 2);
  bids_fixture_close(&ctx, &result);
}

static void bids_perfect_powers(void) {
  siqs_factor_array_t result;
  mpz_t n, d, product, *values;
  uint32_t i, count, sevens = 0, elevens = 0;
  case_name = "prime-and-composite-power-refinement";
  mpz_init_set_ui(n, 41503); mpz_init_set_ui(d, 343); mpz_init_set_ui(product, 1);
  siqs_factor_array_init(&result, n);
  CHECK(siqs_insert_divisor(&result, d)); /* 7^3 * 11^2. */
  CHECK(siqs_all_factors_prime(&result) && result.count == 5);
  for (i = 0; i < result.count; i++) {
    sevens += mpz_cmp_ui(result.values[i], 7) == 0;
    elevens += mpz_cmp_ui(result.values[i], 11) == 0;
  }
  CHECK(sevens == 3 && elevens == 2);
  values = siqs_factor_array_release(&result, &count); gmp_siqs_free(values, count);
  mpz_ui_pow_ui(n, 35, 6);
  siqs_factor_array_init(&result, n);
  CHECK(!siqs_all_factors_prime(&result) && result.count == 6);
  for (i = 0; i < result.count; i++) {
    CHECK(mpz_cmp_ui(result.values[i], 35) == 0);
    mpz_mul(product, product, result.values[i]);
  }
  CHECK(mpz_cmp(product, n) == 0);
  CHECK(!siqs_all_factors_prime(&result) && result.count == 6);
  values = siqs_factor_array_release(&result, &count); gmp_siqs_free(values, count);
  mpz_clear(n); mpz_clear(d); mpz_clear(product);
}

#ifdef PSIQS
static void bids_worker_late_relation(void) {
  siqs_ctx_t master;
  siqs_factor_array_t result;
  psiqs_worker_t worker;
  siqs_factor_t f = {1, 1};
  mpz_t y;
  case_name = "queued-worker-completes-after-soft-stop";
  memset(&worker, 0, sizeof(worker));
  bids_fixture_open(&master, &result, 105, 2);
  siqs_hashset_init(&master.relation_hashes, 16);
  master.fb[0].p = 211; master.bids_cap = 2;
  bids_fixture_add(&master, 1, &f, 1);
  bids_fixture_add(&master, 41, &f, 1);
  CHECK(master.bids_stop == SIQS_BIDS_BUDGET && !master.factor_found);
  bids_fixture_open(&worker.ctx, &worker.result, 105, 2);
  worker.ctx.use_bids = 0; worker.ctx.buffer_relations = 1;
  worker.ctx.fb[0].p = 211;
  mpz_init_set_ui(y, 34);
  siqs_store_raw(&worker.ctx, siqs_raw_relation_new(&worker.ctx, y, &f, 1, 1, 1));
  psiqs_merge_worker(&master, &worker);
  CHECK(master.full_count == 3 && master.bids->processed == 3);
  CHECK(master.factor_found && siqs_all_factors_prime(&result));
  CHECK(worker.ctx.raw_count == 0);
  mpz_clear(y);
  bids_fixture_close(&worker.ctx, &worker.result);
  bids_fixture_close(&master, &result);
}

static void bids_workers(void) {
  uint32_t test;
  for (test = 0; test < 3; test++) {
    worker_fixture_t f;
    uint32_t next = 1, families = 0, polynomials = 0;
    case_name = "worker-budget-drain-and-serial-fallback";
    worker_fault_reset(); worker_fixture_open(&f, test == 0 ? 2 : 4, 0);
    f.ctx.use_bids = 1; f.ctx.bids_cap = 1;
    if (test == 2) { seam_thread_after = 0; seam_quiet = 1; }
    CHECK(psiqs_collect_relations(&f.ctx, &f.dispatch, 1,
                                   &next, &families, &polynomials));
    if (detailed) printf("  bids workers %u: %u families, %u polynomials, %u relations, stop=%u\n",
                         test, families, polynomials, f.ctx.full_count, (unsigned)f.ctx.bids_stop);
    CHECK(f.ctx.bids != NULL && f.ctx.full_count >= f.ctx.bids_cap);
    CHECK(f.ctx.bids->processed == f.ctx.full_count && !f.ctx.bids->failed);
    CHECK(f.ctx.bids_stop == SIQS_BIDS_BUDGET && families > 0 && polynomials > 0);
    CHECK(f.ctx.matrix_ready_incidence == NULL && next == 1);
    if (test == 2) CHECK(seam_created == 0 && seam_warnings == 5);
    worker_resources_clear(); worker_fault_reset(); worker_fixture_close(&f);
  }
  puts("PASS bids: two/four-worker backlog drainage and serial fallback without readiness scans");
}
#endif

static void suite_bids(void) {
  static const uint32_t sizes[] = {2, 63, 64, 65, 66, 127, 128, 129,
                                   511, 512, 513, 4095, 4096};
  siqs_ctx_t ctx;
  uint32_t i;
  case_name = "eligibility";
  memset(&ctx, 0, sizeof(ctx)); ctx.params.fb_size = SIQS_BIDS_MAX_FB;
  CHECK(siqs_bids_enabled(&ctx));
  ctx.params.fb_size++; CHECK(!siqs_bids_enabled(&ctx));
  ctx.params.fb_size--; ctx.buffer_relations = 1; CHECK(!siqs_bids_enabled(&ctx));
  for (i = 0; i < sizeof(sizes) / sizeof(*sizes); i++) bids_layout(sizes[i]);
  bids_soft_budget(2); bids_soft_budget(65);
  bids_uncertified(2); bids_uncertified(65);
  bids_perfect_powers();
  puts("PASS bids: scalar/packed boundaries, shuffled pivots, partial refinement, soft budget and uncertified stops");
#ifdef PSIQS
  prime_iterator_global_startup();
  bids_worker_late_relation();
  bids_workers();
  prime_iterator_global_shutdown();
#endif
}
