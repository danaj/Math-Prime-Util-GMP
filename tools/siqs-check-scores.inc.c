/* Exact byte-score safety checks; included by siqs-check.c, not a unit. */
static void score_reference_upper(mpz_t upper, const siqs_ctx_t *ctx) {
  mpz_t target, a, da, b, y, value, t;
  uint32_t endpoint;
  mpz_inits(target, a, da, b, y, value, t, NULL);
  mpz_mul_ui(target, ctx->kn, 2); mpz_sqrt(target, target);
  mpz_fdiv_q_ui(target, target, ctx->params.half_interval);
  mpz_fdiv_q_ui(target, target, ctx->params.poly_d);
  mpz_set_ui(upper, 0);
  for (endpoint = 0; endpoint < 2; endpoint++) {
    if (endpoint) mpz_mul_ui(a, target, 4);
    else mpz_cdiv_q_ui(a, target, 4);
    mpz_mul_ui(da, a, ctx->params.poly_d);
    mpz_mul_ui(b, a, ctx->params.q_count + 2U); mpz_fdiv_q_2exp(b, b, 1);
    mpz_mul_ui(y, da, ctx->params.half_interval); mpz_add(y, y, b);
    mpz_mul(value, y, y); mpz_sub(value, value, ctx->kn);
    mpz_abs(value, value); mpz_cdiv_q(value, value, da);
    if (mpz_cmp(value, upper) > 0) mpz_set(upper, value);
    mpz_cdiv_q(t, ctx->kn, da);
    if (mpz_cmp(t, upper) > 0) mpz_set(upper, t);
  }
  if (ctx->params.poly_d == 2) mpz_fdiv_q_2exp(upper, upper, 2);
  mpz_clears(target, a, da, b, y, value, t, NULL);
}

/* Independent Gray-code enumeration, not another dynamic program. */
static void check_score_certificates(void) {
  enum { count = 16 };
  siqs_ctx_t ctx;
  siqs_fb_t fb[count + 1];
  mpz_t product, upper, prime;
  uint32_t round, i, m, score, maximum, bound;
  case_name = "exact-score-certificate-subsets";
  memset(&ctx, 0, sizeof(ctx)); memset(fb, 0, sizeof(fb));
  ctx.fb = fb; ctx.params.fb_size = count + 1; ctx.params.sieve_start = 1;
  ctx.params.q_count = 4; ctx.params.half_interval = 4096;
  fb[0].p = 2;
  mpz_init(ctx.kn); mpz_inits(product, upper, prime, NULL);
  for (round = 0; round < 18; round++) {
    uint32_t previous = 0;
    int certified;
    ctx.params.bits = 129; ctx.params.sieve_start = 1;
    ctx.params.sieve_score_bound = 0;
    ctx.params.poly_d = round & 1U ? 2 : 1;
    mpz_set_ui(ctx.kn, 1); mpz_mul_2exp(ctx.kn, ctx.kn, 40 + 2 * round);
    mpz_add_ui(ctx.kn, ctx.kn, 1);
    mpz_set_ui(prime, round < 2 ? 4090 : 3);
    for (i = 1; i <= count; i++) {
      mpz_nextprime(prime, prime); fb[i].p = (uint32_t)mpz_get_ui(prime);
      fb[i].logp = (uint8_t)(round < 2 ? SIQS_SIEVE_SCORE_TOP - 1
                          : round < 4 ? 1 + (i % 3) * (SIQS_SIEVE_SCORE_TOP / 8)
                          : SIQS_SIEVE_SCORE_TOP / 4);
    }
    if (round < 2) {
      /* These weights are constant.  Mark their trivial capped bound to
       * exercise skipping a same-weight tail after the retention limit. */
      ctx.params.sieve_score_bound = SIQS_SIEVE_SCORE_MAX + 1U;
    } else if (round >= 6) {
      /* Exercise the cheap bound with real rounded weights, both d values,
       * and different scales; the same independent enumeration checks it. */
      ctx.params.sieve_free_units = SIQS_SIEVE_SCORE_TOP / 8U;
      ctx.params.sieve_hit_bound = round < 12 ? 16.0 : 1048576.0;
      siqs_set_log_weights(&ctx);
      CHECK(ctx.params.sieve_score_bound != 0);
    }
    score_reference_upper(upper, &ctx);
    maximum = score = 0; mpz_set_ui(product, 1);
    for (m = 1; m < (1U << count); m++) {
      uint32_t gray = m ^ (m >> 1), changed = gray ^ previous;
      uint32_t index = siqs_ctz32(changed) + 1;
      if (index >= ctx.params.sieve_start) {
        if (gray & changed) { mpz_mul_ui(product, product, fb[index].p); score += fb[index].logp; }
        else { mpz_divexact_ui(product, product, fb[index].p); score -= fb[index].logp; }
      }
      if (mpz_cmp(product, upper) <= 0 && score > maximum) maximum = score;
      previous = gray;
    }
    certified = siqs_certify_sieve_scores(&ctx, &bound);
    CHECK(certified == (maximum <= SIQS_SIEVE_SCORE_MAX));
    CHECK(bound == (maximum > SIQS_SIEVE_SCORE_MAX ? SIQS_SIEVE_SCORE_MAX + 1U : maximum));
    if (ctx.params.sieve_score_bound != 0) {
      CHECK(bound <= ctx.params.sieve_score_bound);
      CHECK(mpz_sizeinbase(upper, 2) <= siqs_sieve_score_envelope_bits(&ctx));
      if (ctx.params.sieve_score_bound <= SIQS_SIEVE_SCORE_MAX) {
        ctx.params.sieve_initial = (uint8_t)(SIQS_SIEVE_SCORE_MAX
                                             - ctx.params.sieve_score_bound);
        ctx.params.stage1_bias = 0;
        siqs_select_sieve_score_mode(&ctx);
        CHECK(!ctx.params.sieve_raw_scores && !ctx.params.sieve_saturating);
      }
    }
    ctx.params.bits = 129; ctx.params.sieve_initial = 0; ctx.params.stage1_bias = 0;
    siqs_select_sieve_score_mode(&ctx);
    CHECK(ctx.params.sieve_saturating == !certified);
    CHECK(ctx.params.sieve_raw_scores == !certified);
    if (certified) {
      /* Exact physical-byte boundary, independent of the usual initializer. */
      ctx.params.sieve_initial = (uint8_t)(SIQS_SIEVE_SCORE_MAX - bound);
      siqs_select_sieve_score_mode(&ctx); CHECK(!ctx.params.sieve_raw_scores);
      ctx.params.stage1_bias = 1;
      siqs_select_sieve_score_mode(&ctx);
      CHECK(ctx.params.sieve_raw_scores && !ctx.params.sieve_saturating);
    }
  }
  ctx.params.q_count = 3; CHECK(!siqs_certify_sieve_scores(&ctx, &bound));
  CHECK(siqs_sieve_score_envelope_bits(&ctx) == 0);
  ctx.params.q_count = 4; mpz_set_ui(ctx.kn, 49);
  CHECK(!siqs_certify_sieve_scores(&ctx, &bound));
  CHECK(siqs_sieve_score_envelope_bits(&ctx) == 0);
  mpz_set_ui(ctx.kn, 57); ctx.params.poly_d = 1;
  siqs_set_log_weights(&ctx);
  CHECK(ctx.params.sieve_score_bound == 0);
  CHECK(!siqs_certify_sieve_scores(&ctx, &bound)); /* rounded target is zero */
  siqs_select_sieve_score_mode(&ctx);
  CHECK(ctx.params.sieve_raw_scores && ctx.params.sieve_saturating);
  mpz_set_ui(ctx.kn, 0); ctx.params.bits = 128;
  ctx.params.sieve_raw_scores = ctx.params.sieve_saturating = 1;
  siqs_select_sieve_score_mode(&ctx);
  CHECK(!ctx.params.sieve_raw_scores && !ctx.params.sieve_saturating);
  mpz_clears(ctx.kn, product, upper, prime, NULL);
  puts("PASS scores: subset certificates, cheap headroom proof, exact byte boundaries, uncertified cases, small-input bypass");
}

static void check_raw_score_scan(void) {
  siqs_ctx_t ctx;
  uint32_t t, lane, a, b, i, cells[] = {2, 7};
  case_name = "raw-score-scanner";
  for (t = 1; t <= 128; t++) {
    uint64_t delta = UINT64_C(0x0101010101010101) * (128U - t);
    for (lane = 0; lane < 8; lane++) for (a = 0; a < 256; a++) {
      for (b = 0; b < (extended ? 256U : 4U); b++) {
        uint32_t other = extended ? b : (b == 0 ? 0 : b == 1 ? 127 : b == 2 ? 128 : 255);
        uint64_t word = ((uint64_t)a << (8 * lane)) |
                        ((uint64_t)other << (8 * ((lane + 1) & 7)));
        CHECK((siqs_raw_candidate_word(word, delta) != 0) == (a >= t || other >= t));
      }
    }
  }
  synthetic_init(&ctx, 65560, cells, 2);
  ctx.params.sieve_raw_scores = 1;
  for (t = 0; t <= SIQS_SIEVE_SCORE_TOP; t++) {
    uint32_t count = 0;
    ctx.sieve_length = t == 0 ? 65560 : 65;
    ctx.params.half_interval = ctx.sieve_length / 2;
    ctx.active_sieve_initial = (uint8_t)(SIQS_SIEVE_SCORE_TOP - t);
    for (i = 0; i < ctx.sieve_length; i++) ctx.sieve[i] = (uint8_t)(siqs_mix64(i + t) & SIQS_SIEVE_SCORE_MAX);
    siqs_find_candidates(&ctx);
    for (i = 0; i < ctx.sieve_length; i++) if (ctx.sieve[i] >= t) {
      CHECK(count < ctx.candidate_count);
      CHECK(ctx.candidates[count].x == (int32_t)i - (int32_t)ctx.params.half_interval);
      CHECK(ctx.candidates[count].sieve_score == ctx.sieve[i]); count++;
    }
    CHECK(count == ctx.candidate_count);
    if (!t) CHECK(ctx.candidate_wide);
    siqs_clear_candidate_map(&ctx);
  }
  ctx.active_sieve_initial = (uint8_t)(SIQS_SIEVE_SCORE_TOP + 1U);
  siqs_find_candidates(&ctx); CHECK(ctx.candidate_count == ctx.sieve_length);
  siqs_clear_candidate_map(&ctx); synthetic_clear(&ctx);
  puts("PASS scores: adjustable word thresholds, byte tails, all-candidate scans and wide maps");
}

static void check_saturating_scores(void) {
  const uint32_t lengths[] = {33, 4096, 65537, 164097};
  uint32_t a, b, i, k, round, count = extended ? 24000 : 6000;
  uint32_t *periods = (uint32_t *)check_allocate(count, sizeof(uint32_t));
  mpz_t prime;
  case_name = "saturating-kernel-oracle";
  for (a = 0; a <= SIQS_SIEVE_SCORE_MAX; a++) for (b = 0; b < SIQS_SIEVE_SCORE_TOP; b++) {
    uint8_t cell = (uint8_t)a;
    siqs_sieve_add_sat(&cell, (uint8_t)b);
    CHECK(cell == (a + b > SIQS_SIEVE_SCORE_MAX ? SIQS_SIEVE_SCORE_MAX : a + b));
  }
  mpz_init_set_ui(prime, 1);
  for (i = 0; i < count; i++) { mpz_nextprime(prime, prime); periods[i] = (uint32_t)mpz_get_ui(prime); }
  for (k = 0; k < sizeof(lengths) / sizeof(*lengths); k++) {
    siqs_ctx_t ctx;
    synthetic_init(&ctx, lengths[k], periods, count);
    ctx.params.sieve_raw_scores = ctx.params.sieve_saturating = 1;
    ctx.active_sieve_initial = SIQS_SIEVE_SCORE_TOP / 2;
    for (round = 0; round < 3; round++) {
      synthetic_roots(&ctx, round);
      for (i = 0; i < count; i++) ctx.sieve_logp[i] = (uint8_t)(1 + i % (SIQS_SIEVE_SCORE_TOP - 1));
      compare_sieve(&ctx);
    }
    synthetic_clear(&ctx);
  }
  free(periods); mpz_clear(prime);
  puts("PASS scores: exhaustive saturation, all hit tiers, equal roots, blocking and padding");
}

/* A genuine polynomial and FB-smooth Q(0), not an accepted-A fixture.
 * Force the saturating path to check the filter/evaluator independently. */
static void check_saturated_fine_filter(void) {
  enum { fb_count = 48 };
  siqs_ctx_t ctx;
  siqs_poly_t poly;
  siqs_factor_array_t result;
  mpz_t n, A, q0, prime;
  uint32_t i, x, roots[2], found, row;
  case_name = "saturated-fine-filter-exact-relation";
  memset(&ctx, 0, sizeof(ctx)); mpz_inits(n, A, q0, prime, NULL);
  ctx.params.fb_size = fb_count; ctx.params.q_count = 4;
  ctx.params.poly_d = 1; ctx.params.half_interval = 4096; ctx.sieve_length = 8192;
  ctx.params.max_large_primes = 2; ctx.params.sieve_start = 1;
  ctx.params.sieve_raw_scores = ctx.params.sieve_saturating = 1;
  ctx.params.smooth_bound = ctx.params.large_prime_bound = 1;
  ctx.params.log_scale = 100.0; ctx.active_sieve_initial = 1;
  ctx.fb = (siqs_fb_t *)check_allocate(fb_count, sizeof(*ctx.fb));
  ctx.prime = (uint32_t *)check_allocate(fb_count, sizeof(uint32_t));
  ctx.root1 = (uint32_t *)check_allocate(fb_count, sizeof(uint32_t));
  ctx.root2 = (uint32_t *)check_allocate(fb_count, sizeof(uint32_t));
  ctx.fb_reciprocal = (uint32_t *)check_allocate(fb_count, sizeof(uint32_t));
  ctx.sieve_logp = (uint8_t *)check_allocate(fb_count, 1);
  ctx.fb_in_a = (uint8_t *)check_allocate(fb_count, 1);
  mpz_set_ui(prime, 1); mpz_set_ui(A, 1); mpz_set_ui(q0, 1);
  for (i = 0; i < fb_count; i++) {
    uint32_t p;
    mpz_nextprime(prime, prime); p = (uint32_t)mpz_get_ui(prime);
    ctx.prime[i] = ctx.fb[i].p = p; ctx.fb[i].sqrt_kn = 1;
    ctx.sieve_logp[i] = ctx.fb[i].logp = SIQS_SIEVE_SCORE_TOP - 1U;
    ctx.fb_reciprocal[i] = (uint32_t)(UINT64_C(0x100000000) / p);
    if (p == 17 || p == 19 || p == 23 || p == 29) { ctx.fb_in_a[i] = 1; mpz_mul_ui(A, A, p); }
    else mpz_mul_ui(q0, q0, p);
  }
  mpz_mul(n, A, q0); mpz_add_ui(n, n, 1);
  siqs_factor_array_init(&result, n); ctx.result = &result; ctx.original_n = n;
  mpz_init_set(ctx.n, n); mpz_init_set(ctx.kn, n);
  mpz_init(ctx.eval.q); mpz_init(ctx.eval.y); mpz_init(ctx.eval.rest);
  ctx.largest_fb_prime = ctx.prime[fb_count - 1];
  siqs_poly_init(&ctx, &poly); mpz_set(poly.A, A); mpz_set(poly.DA, A);
  mpz_set_ui(poly.B, 1); mpz_neg(poly.C, q0); row = 0;
  for (i = 0; i < fb_count; i++) {
    uint32_t p = ctx.prime[i], am = (uint32_t)mpz_fdiv_ui(A, p), cm = (uint32_t)mpz_fdiv_ui(poly.C, p);
    found = 0;
    if (ctx.fb_in_a[i]) poly.a_index[row++] = i;
    for (x = 0; x < p; x++) if (((uint64_t)am * x * x + 2U * x + cm) % p == 0) {
      CHECK(found < 2); roots[found++] = (x + ctx.params.half_interval) % p;
    }
    CHECK(found); if (found == 1) roots[1] = roots[0];
    ctx.root1[i] = roots[0] < roots[1] ? roots[0] : roots[1];
    ctx.root2[i] = roots[0] < roots[1] ? roots[1] : roots[0];
  }
  CHECK(row == 4); ctx.resieve_identity_index = ctx.resieve_one_subtract_index = fb_count;
  ctx.sieve = (uint8_t *)check_allocate(check_sieve_allocation(&ctx), 1);
  ctx.candidate_at = (uint16_t *)check_allocate(ctx.sieve_length, sizeof(uint16_t));
  ctx.factor_counts = (uint32_t *)check_allocate(fb_count + 1, sizeof(uint32_t));
  ctx.factor_touched = (uint32_t *)check_allocate(fb_count + 1, sizeof(uint32_t));
  siqs_graph_init(&ctx.graph); siqs_hashset_init(&ctx.relation_hashes, 128);
  siqs_run_sieve(&ctx); CHECK(ctx.sieve[4096] == SIQS_SIEVE_SCORE_MAX);
  siqs_find_candidates(&ctx); CHECK(ctx.candidate_at[4096] || ctx.candidate_wide);
  siqs_clear_candidate_map(&ctx); siqs_add_candidate(&ctx, 4096);
  ctx.params.sieve_saturating = 0; siqs_filter_candidates(&ctx, &poly);
  CHECK(ctx.candidate_count == 0);
  ctx.params.sieve_saturating = 1; siqs_clear_candidate_map(&ctx); siqs_add_candidate(&ctx, 4096);
  siqs_filter_candidates(&ctx, &poly); CHECK(ctx.candidate_count == 1 && ctx.hit_count == 1);
  siqs_resieve_candidates(&ctx, 1, fb_count); siqs_evaluate_candidate(&ctx, &poly, 0);
  CHECK(ctx.accepted_smooth == 1 && ctx.full_count == 1);
  relation_verify_full(&ctx, ctx.full[0]);
  siqs_poly_clear(&ctx, &poly); siqs_ctx_clear(&ctx);
  for (i = 0; i < result.count; i++) mpz_clear(result.values[i]);
  free(result.values); free(result.primality); mpz_clears(n, A, q0, prime, NULL);
  puts("PASS scores: saturation bypass preserves omitted hits and a verified smooth relation");
}

/* Extended structural checks, not factoring benchmarks or new supported
 * policies. Preserve the worst-case kN exactly; suppress only the FB builder's
 * early-divisor shortcut for the artificial input 2^bits-1. No sieve-length
 * arrays, relation graph, or matrix are allocated. */
static void check_score_upper_envelope(void) {
  static const struct {
    uint32_t bits, k, q, count, half;
  } cases[] = {
    {400, 2047, 12, 60000, 2097152},
    {432, 2047, 12, 60000, 4194304},
    {512,    1, 12, 60000, 16777216},
    {512, 2047, 12, 60000, 16777216},
    {512, 2047, 13, 60000, 16777216},
    {512, 2047, 13, 600000, 16777216}
  };
  mpz_t n, y, value, magnitude, upper, t;
  uint32_t c;
  case_name = "upper-score-envelope";
  prime_iterator_global_startup();
  mpz_inits(n, y, value, magnitude, upper, t, NULL);
  for (c = 0; c < sizeof(cases) / sizeof(*cases); c++) {
    siqs_ctx_t ctx;
    siqs_poly_t poly;
    siqs_factor_array_t result;
    siqs_policy_band_t band = *siqs_policy_band(400);
    double ln_n = cases[c].bits * M_LN2;
    uint32_t i, polynomial, bound, observed = 0;
    int certified;
    band.first_bits = band.last_bits = (uint16_t)cases[c].bits;
    band.q_count = (uint8_t)cases[c].q;
    band.k_l = (siqs_policy_curve_t){384.0, 0.0, band.first_bits, UINT16_MAX};
    band.k_r = (siqs_policy_curve_t){768.0, 0.0, band.first_bits, UINT16_MAX};
    band.fixed_half_interval = cases[c].half;
    band.fb_coefficient = (siqs_policy_curve_t){
      log(cases[c].count + 0.5) / sqrt(ln_n * log(ln_n)),
      0.0, band.first_bits, UINT16_MAX};
    mpz_set_ui(n, 1); mpz_mul_2exp(n, n, cases[c].bits); mpz_sub_ui(n, n, 1);
    siqs_factor_array_init(&result, n);
    siqs_ctx_init(&ctx, n, n, &result, &band, 0); case_ctx = &ctx;
    ctx.multiplier = cases[c].k; mpz_mul_ui(ctx.kn, n, ctx.multiplier);
    ctx.params.poly_d = mpz_fdiv_ui(ctx.kn, 8) == 1 ? 2U : 1U;
    mpz_set_ui(ctx.n, 1); CHECK(siqs_build_factor_base(&ctx));
    CHECK(ctx.params.fb_size == cases[c].count);
    siqs_set_large_prime_bounds(&ctx); siqs_set_log_weights(&ctx);
    certified = siqs_certify_sieve_scores(&ctx, &bound);
    CHECK(ctx.params.sieve_score_bound != 0);
    CHECK(bound <= ctx.params.sieve_score_bound);
    score_reference_upper(upper, &ctx);
    CHECK(mpz_sizeinbase(upper, 2) <= siqs_sieve_score_envelope_bits(&ctx));
    siqs_select_sieve_score_mode(&ctx);
    CHECK(ctx.params.sieve_saturating == !certified);
    CHECK(!ctx.params.sieve_saturating || ctx.params.sieve_raw_scores);
    ctx.prime = siqs_malloc(siqs_array_size(ctx.params.fb_size, sizeof(uint32_t)));
    ctx.fb_in_a = siqs_calloc(ctx.params.fb_size, 1);
    ctx.root1 = siqs_malloc(siqs_array_size(ctx.params.fb_size, sizeof(uint32_t)));
    ctx.root2 = siqs_malloc(siqs_array_size(ctx.params.fb_size, sizeof(uint32_t)));
    for (i = 0; i < ctx.params.fb_size; i++) ctx.prime[i] = ctx.fb[i].p;
    siqs_hashset_init(&ctx.a_hashes, 1024); siqs_poly_init(&ctx, &poly);
    CHECK(siqs_new_family(&ctx, &poly)); score_reference_upper(upper, &ctx);
    for (polynomial = 0; polynomial < 3; polynomial++) {
      uint32_t sample;
      mpz_mul_ui(t, poly.A, 4); CHECK(mpz_cmp(t, poly.target_A) >= 0);
      mpz_mul_ui(t, poly.target_A, 4); CHECK(mpz_cmp(poly.A, t) <= 0);
      mpz_mul_ui(t, poly.A, poly.q_count + 2U); mpz_fdiv_q_2exp(t, t, 1);
      mpz_abs(magnitude, poly.B); CHECK(mpz_cmp(magnitude, t) <= 0);
      for (sample = 0; sample < 5; sample++) {
        uint32_t pos = sample * (cases[c].half / 2U), score = 0;
        mpz_set_si(y, (int32_t)pos - (int32_t)cases[c].half);
        mpz_mul(y, y, poly.DA); mpz_add(y, y, poly.B);
        mpz_mul(value, y, y); mpz_sub(value, value, ctx.kn);
        CHECK(mpz_divisible_p(value, poly.DA)); mpz_divexact(value, value, poly.DA);
        CHECK(mpz_sgn(value) != 0); mpz_abs(magnitude, value);
        if (ctx.params.poly_d == 2) {
          CHECK(mpz_divisible_2exp_p(magnitude, 2)); mpz_fdiv_q_2exp(magnitude, magnitude, 2);
        }
        CHECK(mpz_cmp(magnitude, upper) <= 0);
        for (i = ctx.params.sieve_start; i < ctx.params.fb_size; i++) {
          uint32_t residue = pos % ctx.prime[i];
          int hit = residue == ctx.root1[i] || residue == ctx.root2[i];
          CHECK(hit == (mpz_fdiv_ui(value, ctx.prime[i]) == 0));
          if (hit) score += ctx.fb[i].logp;
        }
        if (certified) CHECK(score <= bound);
        if (!ctx.params.sieve_raw_scores)
          CHECK(score + siqs_physical_sieve_initial(&ctx) <= SIQS_SIEVE_SCORE_MAX);
        if (score > observed) observed = score;
      }
      if (polynomial < 2) CHECK(siqs_next_B(&ctx, &poly));
    }
    if (detailed) {
      printf("  ENVELOPE bits=%u k=%u q=%u FB=%u bound=%u observed=%u raw=%u sat=%u\n",
        cases[c].bits, cases[c].k, cases[c].q, ctx.params.fb_size, bound, observed,
        ctx.params.sieve_raw_scores, ctx.params.sieve_saturating);
      fflush(stdout);
    }
    siqs_poly_clear(&ctx, &poly); siqs_ctx_clear(&ctx); case_ctx = NULL;
    for (i = 0; i < result.count; i++) mpz_clear(result.values[i]);
    free(result.values); free(result.primality);
  }
  mpz_clears(n, y, value, magnitude, upper, t, NULL);
  prime_iterator_global_shutdown();
  puts("PASS scores: 400/432/512-bit setup and polynomial envelope, k=1/2047, q=12/13");
}

#ifdef PSIQS
static void check_score_worker_inheritance(void) {
  worker_fixture_t fixture;
  uint32_t mode;
  case_name = "worker-score-mode-inheritance";
  prime_iterator_global_startup();
  worker_fault_reset(); worker_fixture_open(&fixture, 1, 0);
  for (mode = SIQS_SIEVE_SCORE_BITS == 8U ? 0U : 2U; mode < 3; mode++) {
    psiqs_worker_t worker;
    fixture.ctx.params.sieve_raw_scores = mode != 0;
    fixture.ctx.params.sieve_saturating = mode == 2;
    psiqs_worker_init(&worker, &fixture.ctx, 0);
    CHECK(worker.ctx.params.sieve_raw_scores == fixture.ctx.params.sieve_raw_scores);
    CHECK(worker.ctx.params.sieve_saturating == fixture.ctx.params.sieve_saturating);
    CHECK(psiqs_assign_family(&fixture.ctx, &fixture.dispatch, &worker, 1));
    siqs_first_B_and_roots(&worker.ctx, &worker.poly); siqs_set_family_sieve_initial(&worker.ctx);
    compare_sieve(&worker.ctx); psiqs_worker_clear(&worker);
  }
  worker_fixture_close(&fixture); case_ctx = NULL;
  prime_iterator_global_shutdown();
  puts("PASS scores: worker inheritance of legacy/raw/saturating modes and immutable weights");
}
#endif

static void suite_scores(void) {
  check_score_certificates(); check_raw_score_scan();
  check_saturating_scores(); check_saturated_fine_filter();
  if (extended) check_score_upper_envelope();
#ifdef PSIQS
  check_score_worker_inheritance();
#endif
}
