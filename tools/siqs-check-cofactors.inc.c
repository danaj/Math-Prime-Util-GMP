/* Cofactor acceptance/counter and nested-call checks, not a compilation unit. */
static void cofactor_context(siqs_ctx_t *ctx, siqs_factor_array_t *result) {
  memset(ctx, 0, sizeof(*ctx));
  memset(result, 0, sizeof(*result));
  ctx->result = result;
  ctx->params.bits = 380;
  ctx->params.max_large_primes = 2;
  ctx->params.smooth_bound = UINT64_MAX;
  ctx->params.large_prime_bound = SIQS_LP_MAX;
  ctx->largest_fb_prime = 100;
  ctx->cofactor_rng.state = UINT64_C(0x392d0af6d352bbad);
}
static void cofactor_counts(const siqs_ctx_t *ctx) {
  CHECK(ctx->split_attempts == ctx->split_square + ctx->split_siqs +
        ctx->split_squfof + ctx->split_rho + ctx->split_fail + ctx->split_rejected);
#ifdef SIQS_TIMING
  CHECK(ctx->result->cofactor_calls == ctx->split_attempts);
  CHECK(ctx->result->primality_composite >= ctx->split_attempts);
#endif
}
static int cofactor_resolve(siqs_ctx_t *ctx, const char *value) {
  mpz_t n;
  uint64_t a, b, input;
  int accepted;
  mpz_init_set_str(n, value, 10);
  accepted = siqs_resolve_cofactor(ctx, n, &a, &b);
  if (accepted) {
    CHECK(siqs_mpz_to_u64(n, &input));
    CHECK(a != 0 && input % a == 0 && input / a == b);
    CHECK((a == 1 || (a > ctx->largest_fb_prime && a <= ctx->params.large_prime_bound)) &&
          (b == 1 || (b > ctx->largest_fb_prime && b <= ctx->params.large_prime_bound)));
  } else CHECK(a == 1 && b == 1);
  mpz_clear(n);
  cofactor_counts(ctx);
  return accepted;
}
static void cofactor_classification(void) {
  siqs_ctx_t ctx;
  siqs_factor_array_t result;
  mpz_t divisor;
  uint64_t a, b;
  case_name = "split-methods-versus-failures-and-rejections";
  mpz_init_set_ui(divisor, 0);
  CHECK(!siqs_u64_split_pair(divisor, 10403, &a, &b));
  mpz_set_ui(divisor, 10403);
  CHECK(!siqs_u64_split_pair(divisor, 10403, &a, &b));
  mpz_set_ui(divisor, 2);
  CHECK(!siqs_u64_split_pair(divisor, 10403, &a, &b));
  mpz_set_ui(divisor, 101);
  CHECK(siqs_u64_split_pair(divisor, 10403, &a, &b) && a == 101 && b == 103);
  mpz_set_ui(divisor, 1); mpz_mul_2exp(divisor, divisor, 64);
  CHECK(!siqs_u64_split_pair(divisor, 10403, &a, &b));
  mpz_clear(divisor);
  cofactor_context(&ctx, &result);
  CHECK(cofactor_resolve(&ctx, "10201")); /* 101^2 */
  CHECK(ctx.split_square == 1);
  cofactor_context(&ctx, &result);
  CHECK(cofactor_resolve(&ctx, "10403")); /* 101*103 */
  CHECK(ctx.split_squfof == 1);
  cofactor_context(&ctx, &result); ctx.params.large_prime_bound = 103;
  CHECK(!cofactor_resolve(&ctx, "10807")); /* 101*107: split, LP-invalid */
  CHECK(ctx.split_rejected == 1 && ctx.split_fail == 0 && ctx.split_squfof == 0);
  cofactor_context(&ctx, &result); ctx.params.large_prime_bound = 100;
  CHECK(!cofactor_resolve(&ctx, "10201"));
  CHECK(ctx.split_rejected == 1 && ctx.split_square == 0);
  cofactor_context(&ctx, &result);
  CHECK(!cofactor_resolve(&ctx, "1113121")); /* 101*103*107, not a prime pair */
  CHECK(ctx.split_rejected == 1);
  cofactor_miss_squfof = 1;
  cofactor_miss_siqs = 1;
  cofactor_context(&ctx, &result);
  CHECK(!cofactor_resolve(&ctx, "47053")); /* 211*223, beyond inner trial limit */
  CHECK(ctx.split_fail == 1 && ctx.split_rejected == 0 && ctx.split_rho == 0);
  cofactor_miss_squfof = cofactor_miss_siqs = 0;
  cofactor_context(&ctx, &result);
  CHECK(cofactor_resolve(&ctx, "1") && ctx.split_attempts == 0);
  CHECK(!cofactor_resolve(&ctx, "18446744073709551616") && ctx.split_attempts == 0);
  ctx.params.large_prime_bound = 100;
  CHECK(!cofactor_resolve(&ctx, "107") && ctx.split_attempts == 0);
  ctx.params.smooth_bound = 5000;
  CHECK(!cofactor_resolve(&ctx, "10403") && ctx.split_attempts == 0);
  puts("PASS cofactors: accepted methods, bounded misses, invalid pairs, and counter conservation");
}

static void cofactor_trace_reset(void) {
  cofactor_squfof_calls = cofactor_prime_calls = 0;
}

static void cofactor_lp_ceiling(void) {
  siqs_ctx_t ctx;
  siqs_factor_array_t result;
  case_name = "38-bit-LP-policy-ceiling";
  cofactor_context(&ctx, &result);
  /* Explicit independent R preserves this fixture's maximum residual range
   * even when LP hits its storage ceiling.  Coupled R is tested separately. */
  ctx.params.lp_multiplier = DBL_MAX;
  ctx.params.residual_multiplier = DBL_MAX;
  siqs_set_large_prime_bounds(&ctx);
  CHECK(ctx.params.large_prime_bound == SIQS_LP_MAX);
  CHECK(ctx.params.smooth_bound == UINT64_MAX);
  /* Known primes above the former 36-bit ceiling and just below/above
   * the new ceiling.  Prime residuals need no cofactor split attempt. */
  CHECK(cofactor_resolve(&ctx, "68719476767"));
  CHECK(cofactor_resolve(&ctx, "274877906899"));
  CHECK(!cofactor_resolve(&ctx, "274877906951"));
  CHECK(ctx.split_attempts == 0);
  puts("PASS cofactors: 38-bit LP ceiling, wider accepted labels, unchanged residual bound");
}

static void cofactor_a_search_stages(void) {
  static const struct {
    const char *n;
    uint32_t q, stage;
  } cases[] = {
    { "137438959313", 2, 0 }, /* No tolerance window can contain an A. */
    { "68719477097", 2, 0 },
    { "104339829049", 2, 1 },
    { "350190603377", 2, 2 },
    { "39586268787172817", 3, 2 },
    { "597107386758137", 3, 3 }, /* Exactly one triple in the final window. */
    { "610235247872417", 3, 0 }, /* Smallest triple exceeds every window. */
    { "68719477433", 2, 3 },
    { "137438955233", 2, 3 },
    { "37853468033", 1, 1 },
    { "1479434068249753096553621", 4, 1 },
    { "1479434068249753096553621", 4, 2 },
    { "1479434068249753096553621", 4, 0 }
  };
  siqs_ctx_t ctx;
  siqs_poly_t poly;
  siqs_factor_array_t partition;
  siqs_rng_t rng;
  mpz_t n, scaled;
  mpz_t *factors;
  uint32_t i, j, count, draws, limit, tolerance;
  int success;
  case_name = "A-search-stage-skipping-and-cleanup";
  mpz_init(n); mpz_init(scaled);
  for (i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
    mpz_set_str(n, cases[i].n, 10);
    siqs_factor_array_init(&partition, n);
    siqs_ctx_init(&ctx, n, n, &partition, NULL, 0);
    CHECK(siqs_ctx_allocate(&ctx));
    siqs_poly_init(&ctx, &poly);
    CHECK(poly.q_count == cases[i].q);
    /* q=1 must ignore its tolerance window.  Give the other q=4 fixtures
     * exactly four eligible primes: their sole A fits tolerance 4 but not
     * 2 for stage 2, or only the forbidden final tolerance for exhaustion. */
    if (poly.q_count == 1)
      mpz_set_ui(poly.target_A, 1);
    if (poly.q_count == 4 && cases[i].stage != 1) {
      mpz_set_ui(scaled, 1);
      for (j = 1; j < ctx.params.fb_size; j++) {
        ctx.fb[j].sqrt_kn = j <= 4 ? 1U : 0U;
        if (j <= 4) mpz_mul_ui(scaled, scaled, ctx.fb[j].p);
      }
      mpz_fdiv_q_ui(poly.target_A, scaled, cases[i].stage == 2 ? 3U : 6U);
      CHECK(ctx.params.a_final_tolerance >= 8);
    }
    rng = ctx.poly_rng;
    success = siqs_choose_A(&ctx, &poly);
    CHECK(success == (cases[i].stage != 0));
    /* Each sampled prefix prime consumes exactly one RNG draw.  Small-q
     * fixtures skip impossible stages.  q=4 spends 5K attempts locally
     * before stage 2, and stops after its total 10K budget if that fails. */
    limit = success || poly.q_count == 4 ? 10000U * (poly.q_count - 1U) : 0;
    for (draws = 0; rng.state != ctx.poly_rng.state && draws < limit; draws++)
      (void)siqs_rand64(&rng);
    CHECK(rng.state == ctx.poly_rng.state);
    if (poly.q_count == 4) {
      uint32_t prefix_draws = poly.q_count - 1U;
      if (cases[i].stage == 1)
        CHECK(draws > 0 && draws <= 5000U * prefix_draws);
      else if (cases[i].stage == 2)
        CHECK(draws > 5000U * prefix_draws);
      else
        CHECK(draws == 10000U * prefix_draws);
    }
    if (success) {
      CHECK(ctx.a_hashes.count == 1);
      for (j = 0; j < poly.q_count; j++) {
        CHECK(poly.a_index[j] > 0 && poly.a_index[j] < ctx.params.fb_size);
        CHECK(ctx.fb[poly.a_index[j]].in_a && ctx.fb[poly.a_index[j]].sqrt_kn != 0);
        if (j) CHECK(poly.a_index[j - 1] < poly.a_index[j]);
      }
      if (poly.q_count != 1) {
        tolerance = cases[i].stage == 1 ? 2U : cases[i].stage == 2 ? 4U
                                      : ctx.params.a_final_tolerance;
        mpz_mul_ui(scaled, poly.target_A, tolerance);
        CHECK(mpz_cmp(poly.A, scaled) <= 0);
        mpz_mul_ui(scaled, poly.A, tolerance);
        CHECK(mpz_cmp(scaled, poly.target_A) >= 0);
      } else {
        rng = ctx.poly_rng;
        CHECK(!siqs_choose_A(&ctx, &poly) && ctx.poly_rng.state == rng.state);
      }
    } else {
      CHECK(ctx.a_hashes.count == 0);
      for (j = 0; j < ctx.params.fb_size; j++) CHECK(!ctx.fb[j].in_a);
    }
    if (i == 0) {
      /* Fewer than q eligible odd primes must fail without entering the
       * nearest-available search, regardless of the size of the target. */
      for (j = 1; j < ctx.params.fb_size; j++) ctx.fb[j].sqrt_kn = 0;
      ctx.fb[1].sqrt_kn = 1;
      mpz_set_ui(poly.target_A, 1000000);
      rng = ctx.poly_rng;
      CHECK(!siqs_choose_A(&ctx, &poly) && ctx.poly_rng.state == rng.state);
    }
    if (i == 1) {
      /* The largest eligible product is below even the final window.
       * Reject without consuming random draws, just like a too-small target. */
      mpz_set_ui(poly.target_A, 1);
      mpz_mul_2exp(poly.target_A, poly.target_A, 64);
      rng = ctx.poly_rng;
      CHECK(!siqs_choose_A(&ctx, &poly) && ctx.poly_rng.state == rng.state);
      CHECK(ctx.a_hashes.count == 0);
      for (j = 0; j < ctx.params.fb_size; j++) CHECK(!ctx.fb[j].in_a);
    }
    siqs_poly_clear(&ctx, &poly);
    for (j = 0; j < ctx.params.fb_size; j++) CHECK(!ctx.fb[j].in_a);
    siqs_ctx_clear(&ctx);
    factors = siqs_factor_array_release(&partition, &count);
    gmp_siqs_free(factors, count);
  }
  /* Exercise both late-stage examples through the public entry. */
  mpz_set_str(n, "68719477433", 10);
  factors = gmp_siqs(n, &count, 0, 0);
  CHECK(count == 2 &&
        ((mpz_cmp_ui(factors[0], 431) == 0 && mpz_cmp_ui(factors[1], 159441943) == 0) ||
         (mpz_cmp_ui(factors[1], 431) == 0 && mpz_cmp_ui(factors[0], 159441943) == 0)));
  gmp_siqs_free(factors, count);
  /* This fixture has three prime factors, not a semiprime. */
  mpz_set_str(n, "137438955233", 10);
  factors = gmp_siqs(n, &count, 0, 0);
  CHECK(count == 3); mpz_set_ui(scaled, 1);
  for (i = 0; i < count; i++) {
    CHECK(mpz_cmp_ui(factors[i], 2797) == 0 || mpz_cmp_ui(factors[i], 2819) == 0 ||
          mpz_cmp_ui(factors[i], 17431) == 0);
    mpz_mul(scaled, scaled, factors[i]);
  }
  CHECK(mpz_cmp(scaled, n) == 0); gmp_siqs_free(factors, count);
  mpz_clear(n); mpz_clear(scaled);
  puts("PASS cofactors: A-search stages skip impossible low-q windows, widen q4 within 10K attempts, preserve q1, and clear flags");
}

static void cofactor_low_smooth_recovery(void) {
  static const char *names[] = {
    "smooth_k1_q1_low_recovery_8k",
    "smooth_k1_q1_low_recovery_16k",
    "smooth_k1_q1_low_recovery_legacy"
  };
  static const char *values[][3] = {
    { "84098302697", "203653", "412949" },
    { "81126511553", "8087", "10031719" },
    { "89409580193", "6761", "13224313" },
    { "99430906937", "239689", "414833" },
    { "178537103153", "509", "350760517" },
    { "137438959313", "46099", "2981387" },
    { "68719477097", "24793", "2771729" }
  };
  const siqs_policy_band_t *profiles[3] = { NULL, NULL, NULL };
  siqs_factor_array_t partition;
  siqs_parameters_t primary, recovery;
  mpz_t n, divisor, root, p, q;
  mpz_t *factors;
  uint32_t i, j, bits, count;
  case_name = "37-41-q1-smooth-recovery";
  for (i = 0; i < SIQS_RECOVERY_POLICY_COUNT; i++)
    for (j = 0; j < 3; j++)
      if (strcmp(siqs_recovery_policies[i].name, names[j]) == 0)
        profiles[j] = &siqs_recovery_policies[i];
  mpz_init(n); mpz_init(divisor); mpz_init(root); mpz_init(p); mpz_init(q);
  for (j = 0; j < 3; j++) {
    CHECK(profiles[j] != NULL &&
          profiles[j]->first_bits == MPU_SIQS_MIN_BITS &&
          profiles[j]->last_bits == 41);
    for (bits = 37; bits <= 41; bits++) {
      mpz_set_ui(n, 1); mpz_mul_2exp(n, n, bits - 1U); mpz_add_ui(n, n, 1);
      siqs_select_parameters(&primary, n, NULL);
      siqs_select_parameters(&recovery, n, profiles[j]);
      CHECK(primary.q_count == 2 && recovery.q_count == 1 &&
            recovery.max_large_primes == 1 && recovery.lp_multiplier == 1.0 &&
            recovery.residual_multiplier == 0.0);
      CHECK(primary.fb_size == recovery.fb_size &&
            primary.stage1_bias == recovery.stage1_bias &&
            primary.relation_extra == recovery.relation_extra &&
            primary.sieve_free_units == recovery.sieve_free_units &&
            primary.a_final_tolerance == recovery.a_final_tolerance);
      if (j < 2) CHECK(recovery.half_interval == (8192U << j));
    }
  }
  for (i = 0; i < sizeof(values) / sizeof(values[0]); i++) {
    mpz_set_str(n, values[i][0], 10);
    mpz_set_str(p, values[i][1], 10); mpz_set_str(q, values[i][2], 10);
    CHECK(mpz_probab_prime_p(p, 25) && mpz_probab_prime_p(q, 25));
    mpz_mul(divisor, p, q); CHECK(mpz_cmp(n, divisor) == 0);
    /* Force each profile on the original fixture, even if future primary
     * tuning fixes it; also force the first recovery on the fresh misses. */
    for (j = 0; j < (i == 0 ? 3U : 1U); j++) {
      siqs_factor_array_init(&partition, n);
      CHECK(siqs_try_policy(n, n, &partition, profiles[j], divisor, root, 0, 1));
      factors = siqs_factor_array_release(&partition, &count);
      CHECK(count == 2 &&
            ((mpz_cmp(factors[0], p) == 0 && mpz_cmp(factors[1], q) == 0) ||
             (mpz_cmp(factors[0], q) == 0 && mpz_cmp(factors[1], p) == 0)));
      gmp_siqs_free(factors, count);
    }
    factors = gmp_siqs(n, &count, 0, 0);
    CHECK(count == 2 &&
          ((mpz_cmp(factors[0], p) == 0 && mpz_cmp(factors[1], q) == 0) ||
           (mpz_cmp(factors[0], q) == 0 && mpz_cmp(factors[1], p) == 0)));
    gmp_siqs_free(factors, count);
  }
  mpz_clear(n); mpz_clear(divisor); mpz_clear(root); mpz_clear(p); mpz_clear(q);
  puts("PASS cofactors: smooth q1 recovery covers 37-41 bits and splits the A-exhaustion fixtures");
}

/* Exercise the recovery itself even if future primary tuning also fixes the
 * fixture.  Its ordinary public and nested entry paths are checked below. */
static void cofactor_q3_recovery(void) {
  const siqs_policy_band_t *profile = NULL;
  siqs_factor_array_t partition;
  siqs_parameters_t primary, recovery;
  mpz_t n, divisor, root, p, q;
  mpz_t *factors;
  uint32_t i, count;
  case_name = "50-64-q3-one-lp-recovery";
  for (i = 0; i < SIQS_RECOVERY_POLICY_COUNT; i++)
    if (strcmp(siqs_recovery_policies[i].name, "one_lp_k60_q3_recovery") == 0)
      profile = &siqs_recovery_policies[i];
  CHECK(profile != NULL && profile->first_bits == 50 && profile->last_bits == 64);
  mpz_init(n); mpz_init(divisor); mpz_init(root); mpz_init(p); mpz_init(q);
  for (i = 50; i <= 64; i++) {
    mpz_set_ui(n, 1); mpz_mul_2exp(n, n, i - 1U); mpz_add_ui(n, n, 1);
    siqs_select_parameters(&primary, n, NULL);
    siqs_select_parameters(&recovery, n, profile);
    CHECK(recovery.q_count == 3 && recovery.max_large_primes == 1 &&
          recovery.lp_multiplier == 60.0 && recovery.residual_multiplier == 0.0);
    CHECK(primary.fb_size == recovery.fb_size &&
          primary.half_interval == recovery.half_interval &&
          primary.stage1_bias == recovery.stage1_bias &&
          primary.relation_extra == recovery.relation_extra &&
          primary.sieve_free_units == recovery.sieve_free_units &&
          primary.a_final_tolerance == recovery.a_final_tolerance);
  }
  mpz_set_str(n, "39586268787172817", 10);
  mpz_set_str(p, "4728917", 10); mpz_set_str(q, "8371106701", 10);
  mpz_mul(divisor, p, q); CHECK(mpz_cmp(n, divisor) == 0);
  siqs_factor_array_init(&partition, n);
  CHECK(siqs_try_policy(n, n, &partition, profile, divisor, root, 0, 1));
  factors = siqs_factor_array_release(&partition, &count);
  CHECK(count == 2 &&
        ((mpz_cmp(factors[0], p) == 0 && mpz_cmp(factors[1], q) == 0) ||
         (mpz_cmp(factors[0], q) == 0 && mpz_cmp(factors[1], p) == 0)));
  gmp_siqs_free(factors, count);
  factors = gmp_siqs(n, &count, 0, 0);
  CHECK(count == 2 &&
        ((mpz_cmp(factors[0], p) == 0 && mpz_cmp(factors[1], q) == 0) ||
         (mpz_cmp(factors[0], q) == 0 && mpz_cmp(factors[1], p) == 0)));
  gmp_siqs_free(factors, count);
  mpz_clear(n); mpz_clear(divisor); mpz_clear(root); mpz_clear(p); mpz_clear(q);
  puts("PASS cofactors: 50-64-bit q3 recovery preserves primary geometry and splits the exhaustion fixture");
}

static void cofactor_q2_geometry_recovery(void) {
  static const struct { const char *n, *p, *q; } cases[] = {
    { "597107386758137", "31287313", "19084649" },
    { "610235247872417", "19489627", "31310771" },
    { "1330578890360873", "29816789", "44625157" },
    { "2941303423812353", "48138007", "61101479" },
    { "39586268787172817", "4728917", "8371106701" }
  };
  const siqs_policy_band_t *profile = NULL;
  siqs_factor_array_t partition;
  siqs_parameters_t primary, recovery;
  mpz_t n, divisor, root, p, q;
  mpz_t *factors;
  uint32_t i, count;
  case_name = "50-64-q3-to-q2-geometry-recovery";
  for (i = 0; i < SIQS_RECOVERY_POLICY_COUNT; i++)
    if (strcmp(siqs_recovery_policies[i].name,
               "one_lp_k60_q2_geometry_recovery") == 0)
      profile = &siqs_recovery_policies[i];
  CHECK(profile != NULL && profile->first_bits == 50 && profile->last_bits == 64);
  mpz_init(n); mpz_init(divisor); mpz_init(root); mpz_init(p); mpz_init(q);
  for (i = 50; i <= 64; i++) {
    mpz_set_ui(n, 1); mpz_mul_2exp(n, n, i - 1U); mpz_add_ui(n, n, 1);
    siqs_select_parameters(&primary, n, NULL);
    siqs_select_parameters(&recovery, n, profile);
    CHECK(primary.q_count == 3 && recovery.q_count == 2 &&
          recovery.max_large_primes == 1 && recovery.lp_multiplier == 60.0 &&
          recovery.residual_multiplier == 0.0);
    CHECK(primary.fb_size == recovery.fb_size &&
          primary.half_interval == recovery.half_interval &&
          primary.stage1_bias == recovery.stage1_bias &&
          primary.relation_extra == recovery.relation_extra &&
          primary.sieve_free_units == recovery.sieve_free_units &&
          primary.a_final_tolerance == recovery.a_final_tolerance &&
          primary.multiplier_refine_divisor == recovery.multiplier_refine_divisor);
  }
  for (i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
    mpz_set_str(n, cases[i].n, 10);
    mpz_set_str(p, cases[i].p, 10); mpz_set_str(q, cases[i].q, 10);
    mpz_mul(divisor, p, q); CHECK(mpz_cmp(n, divisor) == 0);
    siqs_factor_array_init(&partition, n);
    CHECK(siqs_try_policy(n, n, &partition, profile, divisor, root, 0, 1));
    factors = siqs_factor_array_release(&partition, &count);
    CHECK(count == 2 &&
          ((mpz_cmp(factors[0], p) == 0 && mpz_cmp(factors[1], q) == 0) ||
           (mpz_cmp(factors[0], q) == 0 && mpz_cmp(factors[1], p) == 0)));
    gmp_siqs_free(factors, count);
    factors = gmp_siqs(n, &count, 0, 0);
    CHECK(count == 2 &&
          ((mpz_cmp(factors[0], p) == 0 && mpz_cmp(factors[1], q) == 0) ||
           (mpz_cmp(factors[0], q) == 0 && mpz_cmp(factors[1], p) == 0)));
    gmp_siqs_free(factors, count);
  }
  mpz_clear(n); mpz_clear(divisor); mpz_clear(root); mpz_clear(p); mpz_clear(q);
  puts("PASS cofactors: terminal q2/K60 recovery changes A geometry and splits all five scarcity fixtures");
}

static void cofactor_cascade(void) {
  siqs_ctx_t ctx;
  siqs_factor_array_t result;
  case_name = "cofactor-preferred-method-and-fallback-order";
  cofactor_trace = 1;
  cofactor_trace_reset(); cofactor_context(&ctx, &result);
  CHECK(cofactor_resolve(&ctx, "47053"));
  CHECK(ctx.split_squfof == 1 && cofactor_squfof_calls == 1 &&
        cofactor_prime_calls == 0);
  cofactor_miss_squfof = 1;
  cofactor_trace_reset(); cofactor_context(&ctx, &result);
  CHECK(cofactor_resolve(&ctx, "47053"));
  CHECK(cofactor_squfof_calls == 1);
  CHECK(ctx.split_siqs == 1 && cofactor_prime_calls != 0);
  cofactor_miss_siqs = 1;
  cofactor_trace_reset(); cofactor_context(&ctx, &result);
  CHECK(!cofactor_resolve(&ctx, "47053"));
  CHECK(ctx.split_fail == 1 && ctx.split_rejected == 0 && ctx.split_rho == 0 &&
        cofactor_squfof_calls == 1);
  CHECK(cofactor_prime_calls == 1);

  /* Above the crossover, a SIQS miss is terminal: no SQUFOF retry. */
  cofactor_miss_squfof = 0;
  cofactor_trace_reset(); cofactor_context(&ctx, &result);
  ctx.largest_fb_prime = 397; ctx.params.large_prime_bound = UINT64_MAX;
  /* 1009 * 9007199254740881: both prime; skip the inner trial limit. */
  CHECK(!cofactor_resolve(&ctx, "9088264048033548929"));
  CHECK(ctx.split_fail == 1 && ctx.split_rejected == 0 && ctx.split_rho == 0 &&
        cofactor_squfof_calls == 0 && cofactor_prime_calls == 1);
  cofactor_trace_reset(); cofactor_context(&ctx, &result); ctx.params.bits = 64;
  ctx.largest_fb_prime = 397; ctx.params.large_prime_bound = UINT64_MAX;
  CHECK(cofactor_resolve(&ctx, "9088264048033548929"));
  CHECK(ctx.split_siqs == 0 && cofactor_squfof_calls == 1 && cofactor_prime_calls == 0);
  /* A nested-size input cannot re-enter SIQS even when SQUFOF misses. */
  cofactor_miss_squfof = 1;
  cofactor_trace_reset(); cofactor_context(&ctx, &result); ctx.params.bits = 64;
  CHECK(!cofactor_resolve(&ctx, "47053"));
  CHECK(ctx.split_fail == 1 && ctx.split_rejected == 0 &&
        cofactor_squfof_calls == 1 && cofactor_prime_calls == 0);
  cofactor_miss_siqs = cofactor_miss_squfof = cofactor_trace = 0;
  puts("PASS cofactors: SQUFOF/SIQS crossover, low-bit SIQS recovery, terminal misses, and recursion guard");
}
#ifndef _WIN32
static void cofactor_miss_notice(void) {
  unsigned int verbose;
  siqs_ctx_t ctx;
  siqs_factor_array_t result;
  mpz_t n;
  case_name = "cofactor-siqs-miss-diagnostic";
  cofactor_miss_squfof = cofactor_miss_siqs = 1;
  mpz_init_set_ui(n, 47053);
  for (verbose = 0; verbose <= 1; verbose++) {
    FILE *capture = tmpfile();
    char message[256];
    int saved_stderr, accepted;
    uint64_t a, b;
    CHECK(capture != NULL); fflush(stderr);
    saved_stderr = dup(STDERR_FILENO); CHECK(saved_stderr >= 0);
    CHECK(dup2(fileno(capture), STDERR_FILENO) >= 0);
    cofactor_context(&ctx, &result); ctx.verbose = (int)verbose;
    accepted = siqs_resolve_cofactor(&ctx, n, &a, &b);
    fflush(stderr);
    CHECK(dup2(saved_stderr, STDERR_FILENO) >= 0); close(saved_stderr);
    CHECK(!accepted && ctx.split_fail == 1 && ctx.split_rejected == 0 &&
          ctx.split_rho == 0 && a == 1 && b == 1);
    cofactor_counts(&ctx);
    rewind(capture);
    if (verbose) {
      CHECK(fgets(message, sizeof(message), capture) != NULL);
      CHECK(strstr(message, "failed to split 47053 (16 bits)") != NULL);
    }
    CHECK(fgets(message, sizeof(message), capture) == NULL);
    fclose(capture);
  }
  mpz_clear(n);
  cofactor_miss_squfof = cofactor_miss_siqs = 0;
  puts("PASS cofactors: SIQS miss drops the relation, reports residual at verbose level 1+, and stays quiet at level 0");
}
#endif

/* Host verbosity is deliberately high: the complete inner call, including
 * reduction and any Lanczos fallback, must remain quiet without mutating it. */
static void cofactor_nested_case(uint32_t bits) {
  siqs_ctx_t ctx;
  siqs_factor_array_t result;
  mpz_t n;
  uint64_t a, b, input;
  uint32_t i;
  static const char *values[] = {
    "36028778899572911", /* 134217689 * 268435399, 55 bits */
    "72057554846356433", /* 268435399 * 268435367, 56 bits */
    "39586268787172817", /* 4728917 * 8371106701: primary A exhaustion, 56 bits */
    "144115156668907691", /* 268435399 * 536870909, 57 bits */
    "288230356824359011", /* 536870909 * 536870879, 58 bits */
    "576460727070490747", /* 536870909 * 1073741783, 59 bits */
    "18446743979220271189" /* 4294967291 * 4294967279, 64 bits */
  };
  for (i = 0; i < sizeof(values) / sizeof(values[0]); i++) {
    cofactor_context(&ctx, &result); ctx.params.bits = bits;
    mpz_init_set_str(n, values[i], 10);
    CHECK(siqs_mpz_to_u64(n, &input));
    CHECK(siqs_resolve_cofactor(&ctx, n, &a, &b));
    CHECK(a > ctx.largest_fb_prime && b > ctx.largest_fb_prime);
    CHECK(input % a == 0 && input / a == b);
    cofactor_counts(&ctx);
    if (bits > 64 && mpz_sizeinbase(n, 2) >= 56U)
      CHECK(ctx.split_siqs == 1);
    else if (bits <= 64) CHECK(ctx.split_siqs == 0);
    mpz_clear(n);
  }
}
static void cofactor_siqs_helper(void) {
  mpz_t n, p, q;
  uint64_t a = 17, b = 19;
  case_name = "siqs-largest-two-u64-partition-values";
  mpz_init(n); mpz_init(p); mpz_init(q);
  CHECK(!siqs_split_to_u64(n, 0, &a, &b) && a == 0 && b == 0);
  mpz_set_si(n, -7);
  CHECK(!siqs_split_to_u64(n, 0, &a, &b) && a == 0 && b == 0);
  mpz_set_ui(n, 1);
  CHECK(!siqs_split_to_u64(n, 0, &a, &b) && a == 0 && b == 0);
  mpz_set_ui(n, 101);
  CHECK(!siqs_split_to_u64(n, 0, &a, &b) && a == 0 && b == 0);
  mpz_set_ui(n, 9);
  CHECK(siqs_split_to_u64(n, 0, &a, &b) && a == 3 && b == 3);
  mpz_set_ui(n, 1155); /* 3*5*7*11: largest two, not a complete pair. */
  CHECK(siqs_split_to_u64(n, 0, &a, &b) && a == 7 && b == 11);
  mpz_set_ui(n, 27); /* Three repeated entries: retain the largest two. */
  CHECK(siqs_split_to_u64(n, 0, &a, &b) && a == 3 && b == 3);
  mpz_set_str(p, "4294967311", 10); mpz_set_str(q, "4294967357", 10);
  CHECK(mpz_probab_prime_p(p, 25) && mpz_probab_prime_p(q, 25));
  mpz_mul(n, p, q); CHECK(mpz_sizeinbase(n, 2) == 65);
  CHECK(siqs_split_to_u64(n, 0, &a, &b) &&
        a == UINT64_C(4294967311) && b == UINT64_C(4294967357));
  mpz_set_str(p, "18446744073709551629", 10);
  CHECK(mpz_probab_prime_p(p, 25)); mpz_mul_ui(n, p, 3);
  CHECK(!siqs_split_to_u64(n, 0, &a, &b) && a == 0 && b == 0);
  mpz_clear(n); mpz_clear(p); mpz_clear(q);
}
#ifdef PSIQS
static void *cofactor_thread(void *unused) {
  unsigned int repeat;
  (void)unused;
  for (repeat = 0; repeat < 16; repeat++) cofactor_nested_case(380);
  return NULL;
}
#endif
#ifdef SIQS_TIMING
static void cofactor_timing(void) {
  siqs_ctx_t ctx, other;
  siqs_factor_array_t first, second;
  mpz_t n;
  uint64_t a, b, elapsed;
  uint32_t count;
  mpz_t *factors;
  case_name = "splitter-timing-and-per-call-ownership";
  mpz_init_set_ui(n, 103);
  cofactor_context(&ctx, &first);
  cofactor_context(&other, &second);
  siqs_factor_array_init(&first, n);
  siqs_factor_array_init(&second, n);
  CHECK(siqs_resolve_cofactor(&ctx, n, &a, &b)); /* Prime: no splitter. */
  CHECK(first.cofactor_time == 0 && first.cofactor_calls == 0);
  CHECK(first.primality_time == 0 && first.primality_prime == 0 &&
        first.primality_composite == 0); /* Below pmax^2: no pretest. */
  ctx.params.bits = 64;
  mpz_set_ui(n, 65537);
  CHECK(siqs_resolve_cofactor(&ctx, n, &a, &b)); /* Prime above pmax^2. */
  CHECK(first.primality_prime == 1 && first.primality_composite == 0 &&
        first.cofactor_calls == 0);
  mpz_set_ui(n, 10201); /* 101^2: always a split attempt, regardless of size. */
  CHECK(siqs_resolve_cofactor(&ctx, n, &a, &b));
  CHECK(first.cofactor_calls == 1);
  elapsed = first.cofactor_time;
  ctx.params.large_prime_bound = 100;
  CHECK(!siqs_resolve_cofactor(&ctx, n, &a, &b)); /* Split, but policy-invalid. */
  CHECK(first.cofactor_calls == 2 && first.cofactor_time >= elapsed);
  mpz_set_ui(n, 65537);
  CHECK(!siqs_resolve_cofactor(&ctx, n, &a, &b)); /* Prime, but above LP. */
  CHECK(first.primality_prime == 2 && first.primality_composite == 2 &&
        first.cofactor_calls == 2);
  elapsed = first.primality_time;
  mpz_set_ui(n, 10201);
  ctx.params.smooth_bound = 100;
  CHECK(!siqs_resolve_cofactor(&ctx, n, &a, &b)); /* Early R rejection. */
  CHECK(first.cofactor_calls == 2);
  CHECK(first.primality_time == elapsed && first.primality_prime == 2 &&
        first.primality_composite == 2);
  CHECK(second.cofactor_time == 0 && second.cofactor_calls == 0);
  CHECK(second.primality_time == 0 && second.primality_prime == 0 &&
        second.primality_composite == 0);
  CHECK(siqs_resolve_cofactor(&other, n, &a, &b));
  CHECK(second.cofactor_calls == 1 && first.cofactor_calls == 2);
  CHECK(second.primality_prime == 0 && second.primality_composite == 1 &&
        first.primality_prime == 2 && first.primality_composite == 2);
  cofactor_counts(&ctx); cofactor_counts(&other);
  factors = siqs_factor_array_release(&first, &count); gmp_siqs_free(factors, count);
  factors = siqs_factor_array_release(&second, &count); gmp_siqs_free(factors, count);
  mpz_clear(n);
  puts("PASS cofactors: split-only timing, no size gate, early exits excluded, private totals");
  puts("PASS cofactors: primality timing, prime/composite answers, shortcut exclusion, private totals");
}
#endif

static void cofactor_mont64_case(uint64_t n, uint64_t a, uint64_t b) {
  mont64_t ctx;
  mpz_t zn, za, zb, expected;
  uint64_t ma, mb, value;
  mpz_init(zn); mpz_init(za); mpz_init(zb); mpz_init(expected);
  mpz_import(zn, 1, -1, sizeof(n), 0, 0, &n);
  mpz_import(za, 1, -1, sizeof(a), 0, 0, &a);
  mpz_import(zb, 1, -1, sizeof(b), 0, 0, &b);
  mont64_init(&ctx, n);
  CHECK(ctx.one > 0 && ctx.one < n && mont64_enter(1, &ctx) == ctx.one);
#if MONT64_HAVE_UINT128
  CHECK(n * ctx.ninv == UINT64_MAX);
  mpz_set_ui(expected, 1); mpz_mul_2exp(expected, expected, 64);
  mpz_mod(expected, expected, zn);
  CHECK(siqs_mpz_to_u64(expected, &value) && value == ctx.one);
  mpz_set_ui(expected, 1); mpz_mul_2exp(expected, expected, 128);
  mpz_mod(expected, expected, zn);
  CHECK(siqs_mpz_to_u64(expected, &value) && value == ctx.r2);
#else
  CHECK(ctx.one == 1);
#endif
  ma = mont64_enter(a, &ctx); mb = mont64_enter(b, &ctx);
  CHECK(ma < n && mb < n);
  CHECK(mont64_exit(ma, &ctx) == a % n && mont64_exit(mb, &ctx) == b % n);
  mpz_mul(expected, za, zb); mpz_mod(expected, expected, zn);
  CHECK(siqs_mpz_to_u64(expected, &value));
  CHECK(mont64_mulmod(a, b, n) == value);
  CHECK(mont64_exit(mont64_mul(ma, mb, &ctx), &ctx) == value);
  mpz_add(expected, za, zb); mpz_mod(expected, expected, zn);
  CHECK(siqs_mpz_to_u64(expected, &value));
  CHECK(mont64_exit(mont64_add(ma, mb, n), &ctx) == value);
  mpz_sub(expected, za, zb); mpz_mod(expected, expected, zn);
  CHECK(siqs_mpz_to_u64(expected, &value));
  CHECK(mont64_exit(mont64_sub(ma, mb, n), &ctx) == value);
  mpz_clear(zn); mpz_clear(za); mpz_clear(zb); mpz_clear(expected);
}
static void cofactor_mont64(void) {
  static const uint64_t moduli[] = {
    3, 5, 7, 65537, UINT64_C(4294967291),
    UINT64_C(9223372036854775807), UINT64_C(9223372036854775809),
    UINT64_C(18446744073709551557), UINT64_MAX - 2, UINT64_MAX
  };
  uint64_t n, values[5];
  uint32_t i, j, k;
  siqs_rng_t rng;
  case_name = "mont64-arithmetic-and-full-width-carries";
  for (i = 0; i < sizeof(moduli) / sizeof(moduli[0]); i++) {
    n = moduli[i];
    values[0] = 0; values[1] = 1; values[2] = n - 2;
    values[3] = n - 1; values[4] = UINT64_MAX;
    for (j = 0; j < 5; j++) for (k = 0; k < 5; k++)
      cofactor_mont64_case(n, values[j], values[k]);
  }
  rng.state = UINT64_C(0x1705ca9eb82d643f);
  for (i = 0; i < 2048; i++) {
    uint64_t a, b;
    n = siqs_rand64(&rng) | 1U;
    if (n < 3) n = 3;
    a = siqs_rand64(&rng); b = siqs_rand64(&rng);
    cofactor_mont64_case(n, a, b);
  }
  printf("PASS cofactors: mont64 %s arithmetic matches GMP, conversions, carries and full-width moduli\n",
         MONT64_HAVE_UINT128 ? "128-bit" : "portable");
}

/* Exercise the MR2 stage independently of the combined BPSW test. */
static int cofactor_mr2_native(uint64_t n) {
  mont64_t ctx;
  if (n < 2 || (n & 1U) == 0)
    return n == 2;
  mont64_init(&ctx, n);
  return siqs_miller_rabin_base_2_mont(&ctx);
}
/* Independent GMP-arithmetic reference for the base-2 stage. */
static int cofactor_mr2_reference(uint64_t n) {
  mpz_t z, d, x, minus_one;
  mp_bitcnt_t s, i;
  int pass;
  if (n < 2 || (n & 1U) == 0)
    return n == 2;
  mpz_init(z); mpz_init(d); mpz_init(x); mpz_init(minus_one);
  mpz_import(z, 1, -1, sizeof(n), 0, 0, &n);
  mpz_sub_ui(minus_one, z, 1);
  s = mpz_scan1(minus_one, 0);
  mpz_fdiv_q_2exp(d, minus_one, s);
  mpz_set_ui(x, 2);
  mpz_powm(x, x, d, z);
  pass = mpz_cmp_ui(x, 1) == 0 || mpz_cmp(x, minus_one) == 0;
  for (i = 1; !pass && i < s; i++) {
    mpz_mul(x, x, x);
    mpz_mod(x, x, z);
    pass = mpz_cmp(x, minus_one) == 0;
  }
  mpz_clear(z); mpz_clear(d); mpz_clear(x); mpz_clear(minus_one);
  return pass;
}
static void cofactor_mr2(void) {
  static const uint64_t edges[] = {
    UINT64_C(1373653), UINT64_C(3215031751), UINT64_C(341550071728321),
    UINT64_C(3825123056546413051), UINT64_C(18446744073709551557), UINT64_MAX
  };
  siqs_rng_t rng;
  siqs_ctx_t ctx;
  siqs_factor_array_t result;
  uint64_t n;
  uint32_t i;
  case_name = "native-base-2-miller-rabin";
  for (n = 0; n < 4096; n++)
    CHECK(cofactor_mr2_native(n) == cofactor_mr2_reference(n));
  for (i = 0; i < sizeof(edges) / sizeof(edges[0]); i++)
    CHECK(cofactor_mr2_native(edges[i]) == cofactor_mr2_reference(edges[i]));
  for (i = 1; i < 64; i++) {
    n = UINT64_C(1) << i;
    CHECK(cofactor_mr2_native(n - 1) == cofactor_mr2_reference(n - 1));
    CHECK(cofactor_mr2_native(n + 1) == cofactor_mr2_reference(n + 1));
  }
  rng.state = UINT64_C(0x98374a52d1cf608b);
  for (i = 0; i < 2048; i++) {
    n = siqs_rand64(&rng);
    if (i & 1U) n |= 1U;
    if (i & 2U) n >>= 32;
    CHECK(cofactor_mr2_native(n) == cofactor_mr2_reference(n));
  }
  /* 2047 = 23*89 passes base 2, but must fail the combined test. */
  CHECK(cofactor_mr2_native(2047));
  cofactor_context(&ctx, &result);
  ctx.largest_fb_prime = 19; ctx.params.large_prime_bound = 100;
  CHECK(cofactor_resolve(&ctx, "2047") && ctx.split_attempts == 1);
#ifdef SIQS_TIMING
  CHECK(result.primality_prime == 0 && result.primality_composite == 1);
#endif
  puts("PASS cofactors: native MR2 matches GMP base-2 reference, edges and pseudoprimes");
}

static void cofactor_bpsw_check(uint64_t n, mpz_t z) {
  mpz_import(z, 1, -1, sizeof(n), 0, 0, &n);
  CHECK(siqs_bpsw_u64(z, n) == (mpz_probab_prime_p(z, 25) != 0));
}
static void cofactor_bpsw(void) {
  static const uint64_t edges[] = {
    UINT64_C(2047), UINT64_C(1373653), UINT64_C(3215031751),
    UINT64_C(341550071728321), UINT64_C(3825123056546413051),
    UINT64_C(18446744073709551557), UINT64_MAX
  };
  static const uint64_t lucas_pseudoprimes[] = {989,3239,5777,10877,27971,29681};
  mpz_t z, a;
  siqs_rng_t rng;
  uint64_t n, value;
  uint32_t i;
  case_name = "native-mr2-and-almost-extra-strong-lucas";
  mpz_init(z); mpz_init(a);
  for (n = 0; n < 65536; n++) cofactor_bpsw_check(n, z);
  for (i = 0; i < sizeof(edges) / sizeof(edges[0]); i++)
    cofactor_bpsw_check(edges[i], z);
  for (i = 0; i < sizeof(lucas_pseudoprimes) / sizeof(lucas_pseudoprimes[0]); i++) {
    mont64_t ctx;
    mont64_init(&ctx, lucas_pseudoprimes[i]);
    CHECK(siqs_lucas_aes_mont(&ctx));
    cofactor_bpsw_check(lucas_pseudoprimes[i], z);
  }
  for (i = 1; i < 64; i++) {
    n = UINT64_C(1) << i;
    cofactor_bpsw_check(n - 1, z); cofactor_bpsw_check(n + 1, z);
  }
  /* Base-2 Wieferich-prime squares exercise the explicit square guard. */
  cofactor_bpsw_check(UINT64_C(1093) * 1093, z);
  cofactor_bpsw_check(UINT64_C(3511) * 3511, z);
  rng.state = UINT64_C(0x29c514e670bd38a1);
  for (i = 0; i < 2048; i++) {
    n = siqs_rand64(&rng) | 1U;
    if (i & 1U) n >>= 32;
    cofactor_bpsw_check(n, z);
    if (n >= 3 && (n & 1U)) {
      value = siqs_rand64(&rng);
      mpz_import(a, 1, -1, sizeof(value), 0, 0, &value);
      CHECK(siqs_jacobi_u64(value, n) == mpz_jacobi(a, z));
    }
    if (i < 256) {
      mpz_nextprime(z, z);
      if (siqs_mpz_to_u64(z, &value)) CHECK(siqs_bpsw_u64(z, value));
      value = siqs_rand64(&rng) & UINT32_MAX;
      cofactor_bpsw_check(value * value, z);
    }
  }
  mpz_clear(z); mpz_clear(a);
  puts("PASS cofactors: native MR2/AES matches GMP, small exhaustive range, 64-bit primes/squares, Jacobi and pseudoprimes");
}

static void suite_cofactors(void) {
  int saved_verbose = verbose_level;
#ifndef _WIN32
  FILE *capture;
  int saved_stdout;
#endif
  case_name = "nested-siqs-ownership-reentrancy-and-quiet-output";
  prime_iterator_global_startup();
  cofactor_mont64();
  cofactor_mr2();
  cofactor_bpsw();
#ifdef SIQS_TIMING
  cofactor_timing();
#endif
  cofactor_a_search_stages();
  cofactor_low_smooth_recovery();
  cofactor_q3_recovery();
  cofactor_q2_geometry_recovery();
  cofactor_classification();
  cofactor_lp_ceiling();
  cofactor_cascade();
#ifndef _WIN32
  cofactor_miss_notice();
#endif
  case_name = "nested-siqs-ownership-reentrancy-and-quiet-output";
  verbose_level = 5;
#ifndef _WIN32
  capture = tmpfile(); CHECK(capture != NULL);
  fflush(stdout); saved_stdout = dup(STDOUT_FILENO); CHECK(saved_stdout >= 0);
  CHECK(dup2(fileno(capture), STDOUT_FILENO) >= 0);
#endif
  cofactor_siqs_helper();
  case_name = "nested-siqs-ownership-reentrancy-and-quiet-output";
  cofactor_nested_case(380);
  cofactor_nested_case(64); /* Outer-input guard prohibits a second SIQS level. */
#ifdef PSIQS
  {
    pthread_t threads[4];
    unsigned int i;
    for (i = 0; i < 4; i++) CHECK(pthread_create(&threads[i], NULL, cofactor_thread, NULL) == 0);
    for (i = 0; i < 4; i++) CHECK(pthread_join(threads[i], NULL) == 0);
  }
#endif
  CHECK(verbose_level == 5);
#ifndef _WIN32
  fflush(stdout); CHECK(ftell(capture) == 0);
  CHECK(dup2(saved_stdout, STDOUT_FILENO) >= 0);
  close(saved_stdout); fclose(capture);
#endif
  verbose_level = saved_verbose;
  prime_iterator_global_shutdown();
  puts("PASS cofactors: largest-two helper, 55-59/64-bit cases, wide input, overflow, and quiet/reentrant calls");
}
