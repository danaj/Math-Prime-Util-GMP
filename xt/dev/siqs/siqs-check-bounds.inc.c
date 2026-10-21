/* Cold-path arithmetic and allocation-boundary checks.  Huge sizes are
 * descriptors only: fatal fixtures must reject them before allocation. */
static void check_arithmetic_bounds(void) {
  const size_t maximum = (size_t)-1;
  const uint64_t maximum_polynomials = UINT32_MAX / 2U;
  uint32_t capacity, i;
  siqs_parameters_t parameters;
  siqs_hashset_t set;
  siqs_graph_t graph;
  case_name = "checked-arithmetic-and-hash-thresholds";
  CHECK(siqs_array_size(0, maximum) == 0);
  CHECK(siqs_array_size(maximum, 1) == maximum);
  CHECK(siqs_array_size(maximum / 8, 8) == maximum - maximum % 8);
  CHECK(siqs_size_add(maximum - 8, 8) == maximum);
  CHECK(siqs_size_add(maximum, 0) == maximum);
  CHECK(siqs_capacity_double(UINT32_C(1) << 30) == (UINT32_C(1) << 31));
  for (capacity = 16;; capacity *= 2U) {
    uint32_t threshold = siqs_hash_grow_at(capacity);
    CHECK((uint64_t)threshold * 10U >= (uint64_t)capacity * 7U);
    CHECK((uint64_t)(threshold - 1U) * 10U < (uint64_t)capacity * 7U);
    CHECK(threshold < capacity);
    if (capacity == (UINT32_C(1) << 31)) break;
  }
  /* Setup bounds suffice for target arithmetic and run-wide u64 statistics. */
  CHECK((uint64_t)SIQS_PACKED_FACTOR_ROW_MASK + 1U + UINT16_MAX < UINT32_MAX);
  CHECK((uint64_t)SIQS_PACKED_FACTOR_ROW_MASK * SIQS_MAX_POLYNOMIALS_PER_FB
        < (UINT64_C(1) << 30));
  CHECK(maximum_polynomials < (UINT64_C(1) << 31));
  /* Even the loose upper bound of one solve retry per newly added relation
   * leaves room for a stopped, zero-polynomial family per worker per retry. */
  CHECK(maximum_polynomials
        + (SIQS_PACKED_FACTOR_ROW_MASK + 1U + SIQS_MAX_EXTRA_RELS)
        < (UINT64_C(1) << 32));
  CHECK(maximum_polynomials
        + (uint64_t)(SIQS_PACKED_FACTOR_ROW_MASK + 1U + SIQS_MAX_EXTRA_RELS)
          * PSIQS_MAX_THREADS < (UINT64_C(1) << 32));
  CHECK((uint64_t)SIQS_HALF_INTERVAL_MAX * 2U == (UINT64_C(1) << 31));
  CHECK(maximum_polynomials * ((uint64_t)SIQS_HALF_INTERVAL_MAX * 2U)
        < (UINT64_C(1) << 62));
  memset(&parameters, 0, sizeof(parameters));
  parameters.fb_size = SIQS_PACKED_FACTOR_ROW_MASK;
  parameters.bits = 400;
  CHECK(siqs_polynomial_limit(&parameters) ==
        SIQS_PACKED_FACTOR_ROW_MASK * SIQS_MAX_POLYNOMIALS_PER_FB);
  parameters.bits = 401;
  CHECK(siqs_polynomial_limit(&parameters) == maximum_polynomials);
  parameters.bits = 512;
  CHECK(siqs_polynomial_limit(&parameters) == maximum_polynomials);
  parameters.fb_size = 1000;
  parameters.bits = 400;
  CHECK(siqs_polynomial_limit(&parameters) == 1000000U);
  parameters.bits = 401;
  CHECK(siqs_polynomial_limit(&parameters) == maximum_polynomials);
  parameters.fb_size = 1;
  CHECK(siqs_polynomial_limit(&parameters) == maximum_polynomials);
  parameters.bits = 400;
  CHECK(siqs_polynomial_limit(&parameters) == 1000000U);
  CHECK(siqs_array_size(SIQS_MAX_A_PRIMES - 1U,
        siqs_array_size(SIQS_PACKED_FACTOR_ROW_MASK, sizeof(uint32_t)))
        == (size_t)(SIQS_MAX_A_PRIMES - 1U) * SIQS_PACKED_FACTOR_ROW_MASK * 4U);
  siqs_hashset_init(&set, 16);
  siqs_graph_init(&graph);
  for (i = 0; i < 3000; i++) {
    CHECK(siqs_hashset_insert(&set, i + 1U));
    CHECK(!siqs_hashset_insert(&set, i + 1U));
    CHECK(set.count == i + 1U && set.count <= set.grow_at);
    CHECK(set.grow_at == siqs_hash_grow_at(set.alloc));
    CHECK(siqs_graph_vertex(&graph, i + 2U) == i + 1U);
    CHECK(siqs_graph_vertex(&graph, i + 2U) == i + 1U);
    CHECK(graph.vertex_count == i + 2U &&
          graph.vertex_count <= graph.table_grow_at);
    CHECK(graph.table_grow_at == siqs_hash_grow_at(graph.table_alloc));
  }
  siqs_hashset_clear(&set);
  siqs_graph_clear(&graph);
  puts("PASS policies: checked byte/count arithmetic, 70% thresholds through 2^31 and real hash/graph growth");
}

static void check_parameter_envelope(void) {
  siqs_policy_band_t band = *siqs_policy_band(300);
  siqs_parameters_t p;
  mpz_t n;
  double ln_n = 300.0 * M_LN2, ln_term = sqrt(ln_n * log(ln_n));
  siqs_ctx_t ctx;
  case_name = "parameter-conversion-envelope";
  mpz_init_set_ui(n, 1); mpz_mul_2exp(n, n, 299);
  band.interval_base_scale = (siqs_policy_curve_t){0.0, 0.0, 300, UINT16_MAX};
  band.interval_lp_scale.amount = 0.0;
  siqs_select_parameters(&p, n, &band);
  CHECK(p.half_interval == 4096);
  band.half_interval_floor = SIQS_SIEVE_ALIGN;
  band.interval_base_scale.begin = -0.01;
  siqs_select_parameters(&p, n, &band);
  CHECK(p.half_interval == SIQS_SIEVE_ALIGN);
  band.interval_base_scale.begin = 0.0;
  band.half_interval_floor = 1792;
  siqs_select_parameters(&p, n, &band);
  CHECK(p.half_interval == 1792);
  band.fixed_half_interval = SIQS_SIEVE_ALIGN;
  siqs_select_parameters(&p, n, &band);
  CHECK(p.half_interval == SIQS_SIEVE_ALIGN);
  band.fixed_half_interval = band.half_interval_floor = 0;
  band.fb_floor = 0;
  band.fb_coefficient = (siqs_policy_curve_t){
      log((double)SIQS_PACKED_FACTOR_ROW_MASK + 0.5) / ln_term,
      0.0, 300, UINT16_MAX};
  band.interval_base_scale.begin = ((double)SIQS_HALF_INTERVAL_MAX + 0.5)
                                  / (6144.0 + exp(0.37 * ln_term));
  siqs_select_parameters(&p, n, &band);
  CHECK(p.fb_size == SIQS_PACKED_FACTOR_ROW_MASK);
  CHECK(p.half_interval == SIQS_HALF_INTERVAL_MAX);
  band.fixed_half_interval = SIQS_HALF_INTERVAL_MAX;
  siqs_select_parameters(&p, n, &band);
  CHECK(p.half_interval == SIQS_HALF_INTERVAL_MAX);
  mpz_clear(n);
  /* A synthetic iterator checks the largest permitted prime without walking
   * there or allocating a huge factor base.  It divides kN, not N. */
  memset(&ctx, 0, sizeof(ctx));
  ctx.params.fb_size = 2;
  mpz_init_set_ui(ctx.n, 323);
  mpz_init_set_ui(ctx.kn, UINT32_C(0x7fffffff));
  mpz_mul_ui(ctx.kn, ctx.kn, 323);
  prime_iterator_global_startup();
  forced_fb_prime = (UV)UINT32_C(0x7fffffff);
  CHECK(siqs_build_factor_base(&ctx));
  CHECK(ctx.largest_fb_prime == UINT32_C(0x7fffffff));
  forced_fb_prime = 0;
  prime_iterator_global_shutdown();
  free(ctx.fb);
  mpz_clear(ctx.n); mpz_clear(ctx.kn);
  puts("PASS policies: valid fractional FB/M limits, zero-scale minimum, fixed M limit and P=2^31-1");
}

/* Independent GMP products check the effective, capped scalar bounds. */
static void check_upper_bound_value(uint64_t actual, uint64_t a, uint64_t b,
                                    double multiplier, uint64_t maximum) {
  mpz_t expected, factor, limit, observed;
  mpz_inits(expected, factor, limit, observed, NULL);
  siqs_mpz_set_u64(expected, a);
  siqs_mpz_set_u64(factor, b);
  mpz_mul(expected, expected, factor);
  siqs_mpz_set_u64(limit, maximum);
  if (multiplier == floor(multiplier)) {
    mpz_set_d(factor, multiplier);
    mpz_mul(expected, expected, factor);
  } else {
    double value;
    CHECK(b == 0 || a <= UINT64_MAX / b);
    value = (double)(a * b) * multiplier;
    if (value >= (double)maximum) mpz_set(expected, limit);
    else mpz_set_d(expected, value);
  }
  if (mpz_cmp(expected, limit) > 0) mpz_set(expected, limit);
  siqs_mpz_set_u64(observed, actual);
  CHECK(mpz_cmp(expected, observed) == 0);
  mpz_clears(expected, factor, limit, observed, NULL);
}

/* Resolve actual upper policies, but never build an FB or worker workspace.
 * P values below are synthetic, including nonprime cap-boundary fixtures. */
static void check_upper_policy_sweep(void) {
  static const unsigned long multipliers[] = {
    1, 2, 3, 5, 7, 19, 41, 43, 2047, 4095
  };
  static const uint32_t p_values[] = {
    1009, 993481, 25695647, 154981283, 154981284,
    715827882, 715827883, UINT32_C(0x7fffffff)
  };
  siqs_ctx_t ctx;
  uint32_t bits, edge, residue, i, j, setups = 0, bounds = 0, d_modes = 0;
  uint32_t maximum_fb = 0, maximum_m = 0, lp_clips = 0, r_clips = 0;
  char name[80];
  memset(&ctx, 0, sizeof(ctx));
  mpz_inits(ctx.n, ctx.kn, NULL);
  case_ctx = &ctx;
  case_name = "upper-bound-GMP-reference";
  check_upper_bound_value(1261, 1009, 1, 1.25, SIQS_LP_MAX);
  check_upper_bound_value(2545202, 1009, 1009, 2.5, SIQS_RESIDUAL_PRODUCT_MAX);
  for (bits = 400; bits <= 512; bits++) {
    const siqs_policy_band_t *band = siqs_policy_band(bits);
    uint32_t owners = 0;
    case_name = "upper-policy-ownership";
    CHECK(band != NULL);
    for (i = 0; i < SIQS_POLICY_BAND_COUNT; i++)
      if (bits >= siqs_policy_bands[i].first_bits &&
          bits <= siqs_policy_bands[i].last_bits) owners++;
    CHECK(owners == 1);
    if (bits == 401U || bits == 431U) {
      const siqs_policy_band_t *previous = siqs_policy_band(bits - 1U);
      case_name = "upper-K-ramp-seam";
      CHECK(siqs_policy_curve_value(&previous->k_l, bits - 1U)
            == siqs_policy_curve_value(&band->k_l, bits));
      CHECK(siqs_policy_curve_value(&previous->k_r, bits - 1U)
            == siqs_policy_curve_value(&band->k_r, bits));
    }
    if (bits == 431U) {
      const siqs_policy_curve_t *previous =
          &siqs_policy_band(430)->fb_coefficient;
      double continued = 2.0 * siqs_policy_curve_value(previous, 430)
                       - siqs_policy_curve_value(previous, 429);
      case_name = "430-431-FB-curve-continuation";
      CHECK(fabs(siqs_policy_curve_value(&band->fb_coefficient, 431)
                 - continued) < 1e-12);
    }
    for (edge = 0; edge < 2; edge++) {
      for (residue = 1; residue < 8; residue += 2) {
        siqs_parameters_t selected;
        mpz_set_ui(ctx.n, 1);
        mpz_mul_2exp(ctx.n, ctx.n, edge ? bits : bits - 1U);
        if (edge) mpz_sub_ui(ctx.n, ctx.n, 8U - residue);
        else mpz_add_ui(ctx.n, ctx.n, residue);
        sprintf(name, "upper-policy-b%u-edge%u-mod8-%u", bits, edge, residue);
        case_name = name;
        CHECK(mpz_sizeinbase(ctx.n, 2) == bits);
        CHECK(mpz_fdiv_ui(ctx.n, 8) == residue);
        memset(&selected, 0, sizeof(selected));
        siqs_select_parameters(&selected, ctx.n, NULL);
        CHECK(selected.bits == bits);
        CHECK(selected.policy_first_bits == band->first_bits);
        CHECK(selected.policy_last_bits == band->last_bits);
        CHECK(selected.max_large_primes == 2);
        CHECK(selected.q_count >= 1 && selected.q_count <= SIQS_MAX_A_PRIMES);
        CHECK(selected.fb_size >= 1 && selected.fb_size <= SIQS_PACKED_FACTOR_ROW_MASK);
        CHECK(selected.half_interval >= 4096 && selected.half_interval <= SIQS_HALF_INTERVAL_MAX);
        CHECK((selected.half_interval & (SIQS_SIEVE_ALIGN - 1U)) == 0);
        CHECK(selected.target_relations == selected.fb_size + 1U + selected.relation_extra);
        CHECK((uint64_t)siqs_polynomial_limit(&selected) < (UINT64_C(1) << 31));
        CHECK(siqs_polynomial_limit(&selected) == (bits > 400U
              ? UINT32_MAX / 2U
              : selected.fb_size * SIQS_MAX_POLYNOMIALS_PER_FB));
        CHECK(selected.sieve_hit_bound_nominal > 0.0 && selected.sieve_hit_bound_nominal <= DBL_MAX);
        CHECK(selected.sieve_free_units < SIQS_SIEVE_SCORE_TOP);
        if (selected.fb_size > maximum_fb) maximum_fb = selected.fb_size;
        if (selected.half_interval > maximum_m) maximum_m = selected.half_interval;
        for (i = 0; i < sizeof(multipliers) / sizeof(multipliers[0]); i++) {
          ctx.params = selected;
          ctx.multiplier = multipliers[i];
          mpz_mul_ui(ctx.kn, ctx.n, ctx.multiplier);
          ctx.params.poly_d = mpz_fdiv_ui(ctx.kn, 8) == 1 ? 2U : 1U;
          d_modes |= 1U << ctx.params.poly_d;
          CHECK(mpz_sizeinbase(ctx.kn, 2) <= bits + 12U);
          CHECK(ctx.params.bits == mpz_sizeinbase(ctx.n, 2));
          setups++;
          for (j = 0; j < sizeof(p_values) / sizeof(p_values[0]); j++) {
            uint64_t pmax = p_values[j];
            sprintf(name, "upper-policy-b%u-edge%u-mod8-%u-k%lu-P%u",
                    bits, edge, residue, ctx.multiplier, p_values[j]);
            ctx.largest_fb_prime = p_values[j];
            siqs_set_large_prime_bounds(&ctx);
            check_upper_bound_value(ctx.params.large_prime_bound, pmax, 1,
                                    ctx.params.lp_multiplier, SIQS_LP_MAX);
            if (ctx.params.residual_multiplier != 0.0)
              check_upper_bound_value(ctx.params.smooth_bound, pmax, pmax,
                  ctx.params.residual_multiplier, SIQS_RESIDUAL_PRODUCT_MAX);
            else
              check_upper_bound_value(ctx.params.smooth_bound,
                  ctx.params.large_prime_bound, pmax, 1.0, SIQS_RESIDUAL_PRODUCT_MAX);
            CHECK(ctx.params.sieve_hit_bound >= ctx.params.sieve_hit_bound_nominal);
            CHECK(ctx.params.sieve_hit_bound >=
                  ctx.params.sieve_hit_residual_multiplier * (double)ctx.params.smooth_bound);
            CHECK(ctx.params.sieve_hit_bound > 0.0 && ctx.params.sieve_hit_bound <= DBL_MAX);
            if (ctx.params.large_prime_bound == SIQS_LP_MAX) lp_clips++;
            if (ctx.params.smooth_bound == SIQS_RESIDUAL_PRODUCT_MAX) r_clips++;
            bounds++;
          }
        }
      }
    }
  }
  case_ctx = NULL;
  mpz_clears(ctx.n, ctx.kn, NULL);
  case_name = "upper-production-policy-sweep";
  CHECK(d_modes == ((1U << 1) | (1U << 2)));
  printf("PASS policies: 400-512, %u setup classes, %u synthetic-P bounds, "
         "both d modes, peak FB=%u M=%u, LP/R clips=%u/%u; no sieving\n",
         setups, bounds, maximum_fb, maximum_m, lp_clips, r_clips);
}

/* Only the cutoff arithmetic needs testing here: synthetic increasing
 * entries avoid building a factor base or allocating a large sieve. */
static uint32_t check_sieve_start_fixture(double exponent, uint32_t floor) {
  enum { count = 512 };
  siqs_ctx_t ctx;
  siqs_fb_t fb[count];
  uint32_t i, result;
  memset(&ctx, 0, sizeof(ctx)); memset(fb, 0, sizeof(fb));
  for (i = 0; i < count; i++) { fb[i].p = i + 2U; fb[i].sqrt_kn = 1; }
  ctx.fb = fb; ctx.params.fb_size = count;
  ctx.params.half_interval = 4096; ctx.params.poly_d = 1;
  ctx.params.sieve_free_units = 16; ctx.params.sieve_hit_bound = 1.0e30;
  ctx.params.sieve_start_index_exponent = exponent;
  ctx.params.sieve_start_prime_floor = floor;
  mpz_init_set_ui(ctx.kn, 1); mpz_mul_2exp(ctx.kn, ctx.kn, 512);
  mpz_sub_ui(ctx.kn, ctx.kn, 1); mpz_mul_ui(ctx.kn, ctx.kn, 2047);
  siqs_set_log_weights(&ctx);
  result = ctx.params.sieve_start; mpz_clear(ctx.kn);
  return result;
}

static void check_first_sieved_bounds(void) {
  case_name = "omitted-prime-setup-bound";
  CHECK(siqs_icbrt_u64(SIQS_PACKED_FACTOR_ROW_MASK) == 127U);
  CHECK(check_sieve_start_fixture(0.0, 0) == 8U);
  CHECK(check_sieve_start_fixture(0.0, 258) == SIQS_POSTFILTER_MAX_SMALL);
  CHECK(check_sieve_start_fixture(0.5, 0) == 22U);
  puts("PASS policies: maximum-FB cbrt and optional first-sieved overrides fit the checked postfilter boundary");
}

static void check_allocation_errors(void) {
#ifndef _WIN32
  unsigned int fixture;
  case_name = "cold-overflow-diagnostics";
  fflush(NULL);
  for (fixture = 0; fixture < 38; fixture++) {
    int status;
    pid_t child = fork();
    CHECK(child >= 0);
    if (child == 0) {
      const uint32_t limit = UINT32_C(1) << 31;
      siqs_policy_band_t band = *siqs_policy_band(300);
      siqs_parameters_t p;
      siqs_ctx_t ctx;
      siqs_hashset_t set;
      siqs_one_lp_state_t one;
      siqs_graph_t graph;
      siqs_factor_array_t factors;
      siqs_candidate_t candidate;
      siqs_fb_t fb[2];
      mpz_t n;
      uint8_t sieve = 128;
      uint32_t *path = NULL, path_alloc = limit;
      memset(&ctx, 0, sizeof(ctx));
      memset(&set, 0, sizeof(set)); memset(&one, 0, sizeof(one));
      memset(&graph, 0, sizeof(graph)); memset(&factors, 0, sizeof(factors));
      memset(&candidate, 0, sizeof(candidate)); memset(fb, 0, sizeof(fb));
      mpz_init_set_ui(n, 1); mpz_mul_2exp(n, n, 299);
      if (freopen("/dev/null", "w", stderr) == NULL) _exit(1);
      switch (fixture) {
        case 0: (void)siqs_array_size((size_t)-1, 2); break;
        case 1: (void)siqs_size_add((size_t)-1, 1); break;
        case 2: (void)siqs_calloc((size_t)-1, 2); break;
        case 3: (void)siqs_capacity_double(limit); break;
        case 4: siqs_hashset_init(&set, limit + 1U); break;
        case 5: set.alloc = limit; siqs_hashset_grow(&set); break;
        case 6: one.alloc = limit; siqs_one_lp_rehash(&one); break;
        case 7:
          graph.table_alloc = limit;
          graph.vertex_count = graph.table_grow_at = siqs_hash_grow_at(limit);
          (void)siqs_graph_vertex(&graph, 17); break;
        case 8:
          ctx.raw_alloc = ctx.raw_count = limit;
          (void)siqs_store_raw(&ctx, NULL); break;
        case 9:
          ctx.full_alloc = ctx.full_count = limit;
          siqs_store_full(&ctx, NULL); break;
        case 10:
          ctx.candidate_alloc = ctx.candidate_count = limit; ctx.sieve = &sieve;
          siqs_add_candidate(&ctx, 0); break;
        case 11:
          ctx.hit_alloc = ctx.hit_count = limit; ctx.candidates = &candidate;
          siqs_add_hit(&ctx, 0, 0); break;
        case 12: siqs_graph_reserve_cycle_path(&path, &path_alloc, limit + 1U); break;
        case 13: band.fb_coefficient.begin = 1.0; break; /* finite FB > u32 */
        case 14: band.fb_coefficient.begin = 1000.0; break; /* exp overflow */
        case 15: band.fb_floor = 0; band.fb_coefficient.begin = -1000.0; break;
        case 16: band.q_count = 0; break;
        case 17: band.q_count = SIQS_MAX_A_PRIMES + 1U; break;
        case 18: band.fixed_half_interval = SIQS_HALF_INTERVAL_MAX + SIQS_SIEVE_ALIGN; break;
        case 19: band.interval_base_scale.begin = -1.0; break;
        case 20: band.interval_base_scale.begin = 1.0e100; break;
        case 21: band.interval_base_scale.begin = 1.0e6; break;
        case 22: band.sieve_hit_exponent.begin = 1000.0; break;
        case 23: band.sieve_free_units = 128; break;
        case 24:
          ctx.params.fb_size = 2; ctx.fb = fb;
          ctx.params.half_interval = 4096; ctx.params.sieve_hit_bound = 1;
          ctx.params.poly_d = 1; ctx.params.sieve_start_index_exponent = 1000;
          fb[0].p = 2; fb[1].p = 3; mpz_init_set(ctx.kn, n);
          siqs_set_log_weights(&ctx); break;
        case 25: case 26: case 27:
          prime_iterator_global_startup();
          ctx.params.fb_size = 2;
          mpz_init_set_ui(ctx.n, 323); mpz_init_set_ui(ctx.kn, 323);
          forced_fb_prime = fixture == 25 ? (UV)limit : UV_MAX;
#if BITS_PER_WORD == 64
          if (fixture == 27) forced_fb_prime = (UV)(UINT64_C(1) << 32) + 9U;
#endif
          (void)siqs_build_factor_base(&ctx); break;
        case 28: (void)siqs_size_add((size_t)-1 - 7U, 8); break;
        case 29:
          factors.alloc = factors.count = limit;
          siqs_factor_array_append(&factors, n); break;
        case 30: (void)check_sieve_start_fixture(0.0, 259); break;
        case 31: (void)check_sieve_start_fixture(0.9, 0); break;
        case 32: band.half_interval_floor = SIQS_SIEVE_ALIGN - 1; break;
        case 33: band.half_interval_floor = SIQS_SIEVE_ALIGN + 1; break;
        case 34: band.half_interval_floor = SIQS_HALF_INTERVAL_MAX + SIQS_SIEVE_ALIGN; break;
        case 35: band.half_interval_floor = SIQS_SIEVE_ALIGN;
                 band.interval_base_scale.begin = -HUGE_VAL; break;
        case 36: band.half_interval_floor = SIQS_SIEVE_ALIGN;
                 band.interval_base_scale.begin = NAN; break;
        case 37: band.fixed_half_interval = SIQS_SIEVE_ALIGN - 1; break;
      }
      if ((fixture >= 13 && fixture <= 23) || fixture >= 32)
        siqs_select_parameters(&p, n, &band);
      _exit(0);
    }
    CHECK(waitpid(child, &status, 0) == child);
    CHECK(WIFEXITED(status) && WEXITSTATUS(status) == 3);
  }
  puts("PASS policies: 38 preallocation/conversion/growth/cutoff failures rejected without huge allocations");
#else
  puts("SKIP policies: fatal allocation-boundary checks require POSIX fork");
#endif
}
