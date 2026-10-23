/* Production-policy setup only: no sieve, relation collection or solver. */
#include <time.h>

#ifndef SIQS_CHECK_SETUP_TIMEOUT
# define SIQS_CHECK_SETUP_TIMEOUT 120U
#endif
/* Major fixed-array payload only; not an RSS/allocator/sanitizer limit. */
#define SIQS_CHECK_SETUP_BYTES (128U * 1024U * 1024U)

#ifndef _WIN32
static void setup_timeout(int signal_number) {
  static const char message[] = "FAIL setup: per-case wall-time limit exceeded\n";
  ssize_t written;
  (void)signal_number;
  written = write(STDERR_FILENO, message, sizeof(message) - 1U);
  (void)written;
  _exit(1);
}
#endif

static void setup_input(mpz_t n, uint32_t bits, uint32_t residue) {
  mpz_t p, q, limit;
  mpz_inits(p, q, limit, NULL);
  mpz_set_ui(p, 1); mpz_mul_2exp(p, p, bits / 2U);
  mpz_sub_ui(p, p, 1048576); mpz_nextprime(p, p);
  CHECK(mpz_sizeinbase(p, 2) == bits / 2U);
  mpz_set_ui(limit, 1); mpz_mul_2exp(limit, limit, (bits + 1U) / 2U);
  mpz_sub_ui(q, limit, 2097152);
  do { mpz_nextprime(q, q); }
  while ((mpz_fdiv_ui(p, 8) * mpz_fdiv_ui(q, 8)) % 8U != residue);
  CHECK(mpz_cmp(q, limit) < 0 && mpz_cmp(p, q) != 0);
  mpz_mul(n, p, q);
  CHECK(mpz_sizeinbase(n, 2) == bits && mpz_fdiv_ui(n, 8) == residue);
  mpz_clears(p, q, limit, NULL);
}

static void setup_polynomial_value(mpz_t value, const siqs_ctx_t *ctx,
                                   const siqs_poly_t *poly, uint32_t pos) {
  mpz_set_si(value, (int32_t)pos - (int32_t)ctx->params.half_interval);
  mpz_mul(value, value, poly->DA); mpz_add(value, value, poly->B);
  mpz_mul(value, value, value); mpz_sub(value, value, ctx->kn);
  CHECK(mpz_divisible_p(value, poly->DA));
  mpz_divexact(value, value, poly->DA);
}

static void setup_check_polynomial(const siqs_ctx_t *ctx,
    const siqs_poly_t *poly, const mpz_t upper, int certified, uint32_t bound) {
  mpz_t value, product;
  uint32_t i, sample;
  mpz_inits(value, product, NULL);
  mpz_set_ui(product, 1);
  for (i = 0; i < poly->q_count; i++) {
    CHECK(poly->a_index[i] < ctx->params.fb_size);
    CHECK(i == 0 || poly->a_index[i - 1U] < poly->a_index[i]);
    CHECK(ctx->fb_in_a[poly->a_index[i]]);
    mpz_mul_ui(product, product, ctx->prime[poly->a_index[i]]);
  }
  CHECK(mpz_cmp(product, poly->A) == 0);
  mpz_mul_ui(product, poly->A, ctx->params.poly_d);
  CHECK(mpz_cmp(product, poly->DA) == 0);
  mpz_mul_ui(product, poly->A, 4);
  CHECK(mpz_cmp(product, poly->target_A) >= 0);
  mpz_mul_ui(product, poly->target_A, 4);
  CHECK(mpz_cmp(poly->A, product) <= 0);
  mpz_mul_ui(product, poly->A, poly->q_count + 2U);
  mpz_fdiv_q_2exp(product, product, 1);
  mpz_abs(value, poly->B); CHECK(mpz_cmp(value, product) <= 0);
  setup_polynomial_value(value, ctx, poly, ctx->params.half_interval);
  CHECK(mpz_cmp(value, poly->C) == 0);
  for (i = 0; i < ctx->params.fb_size; i++)
    CHECK(ctx->root1[i] <= ctx->root2[i] && ctx->root2[i] < ctx->prime[i]);
  /* Check small/special primes, evenly spaced primes, and every A prime. */
  for (i = 0; i < 32U + poly->q_count; i++) {
    uint32_t j = i < 16U ? i : i < 32U
      ? (uint32_t)((uint64_t)(i - 16U) * (ctx->params.fb_size - 1U) / 15U)
      : poly->a_index[i - 32U];
    uint32_t r;
    for (r = 0; r < (ctx->root1[j] == ctx->root2[j] ? 1U : 2U); r++) {
      setup_polynomial_value(value, ctx, poly, r ? ctx->root2[j] : ctx->root1[j]);
      CHECK(mpz_divisible_ui_p(value, ctx->prime[j]));
    }
  }
  for (sample = 0; sample < 4; sample++) {
    uint32_t m = ctx->params.half_interval;
    uint32_t pos = sample == 0 ? 0 : sample == 1 ? m - 1U
                 : sample == 2 ? m : 2U * m - 1U;
    uint32_t score = 0;
    setup_polynomial_value(value, ctx, poly, pos);
    CHECK(mpz_sgn(value) != 0);
    mpz_abs(value, value);
    if (ctx->params.poly_d == 2) {
      CHECK(mpz_divisible_2exp_p(value, 2));
      mpz_fdiv_q_2exp(value, value, 2);
    }
    CHECK(mpz_cmp(value, upper) <= 0);
    for (i = ctx->params.sieve_start; i < ctx->params.fb_size; i++) {
      uint32_t residue = pos % ctx->prime[i];
      if (residue == ctx->root1[i] || residue == ctx->root2[i])
        score += ctx->fb[i].logp;
    }
    if (certified) CHECK(score <= bound);
    if (!ctx->params.sieve_raw_scores)
      CHECK(score + siqs_physical_sieve_initial(ctx) <= SIQS_SIEVE_SCORE_MAX);
  }
  mpz_clears(value, product, NULL);
}

static void setup_real_case(uint32_t bits, uint32_t k, uint32_t residue,
                            uint32_t *d_modes, uint32_t *score_modes) {
  siqs_ctx_t ctx;
  siqs_poly_t poly;
  siqs_factor_array_t result;
  mpz_t n, upper;
  uint32_t i, family, polynomial, bound = 0;
  uint64_t payload;
  int certified;
  char name[80];
  clock_t begin = clock();
  sprintf(name, "production-b%u-k%u-mod8-%u", bits, k, residue);
  case_name = name;
#ifndef _WIN32
  alarm(SIQS_CHECK_SETUP_TIMEOUT);
#endif
  mpz_inits(n, upper, NULL);
  setup_input(n, bits, residue);
  siqs_factor_array_init(&result, n);
  siqs_ctx_init(&ctx, n, n, &result, NULL, 0);
  case_ctx = &ctx;
  ctx.multiplier = k; mpz_mul_ui(ctx.kn, n, k);
  ctx.params.poly_d = mpz_fdiv_ui(ctx.kn, 8) == 1 ? 2U : 1U;
  *d_modes |= 1U << ctx.params.poly_d;
  CHECK(!mpz_perfect_square_p(ctx.kn));
  payload = (uint64_t)ctx.params.fb_size *
      (sizeof(siqs_fb_t) + 3U * sizeof(uint32_t) + sizeof(uint8_t)
       + (ctx.params.q_count - 1U) * sizeof(uint32_t));
  CHECK(payload <= SIQS_CHECK_SETUP_BYTES);
  CHECK(siqs_build_factor_base(&ctx) && !ctx.factor_found);
  CHECK(ctx.largest_fb_prime < (UINT32_C(1) << 31));
  siqs_set_large_prime_bounds(&ctx);
  check_upper_bound_value(ctx.params.large_prime_bound, ctx.largest_fb_prime, 1,
                         ctx.params.lp_multiplier, SIQS_LP_MAX);
  if (ctx.params.residual_multiplier != 0.0)
    check_upper_bound_value(ctx.params.smooth_bound, ctx.largest_fb_prime,
        ctx.largest_fb_prime, ctx.params.residual_multiplier, SIQS_RESIDUAL_PRODUCT_MAX);
  else
    check_upper_bound_value(ctx.params.smooth_bound, ctx.params.large_prime_bound,
                            ctx.largest_fb_prime, 1.0, SIQS_RESIDUAL_PRODUCT_MAX);
  siqs_set_log_weights(&ctx);
  certified = siqs_certify_sieve_scores(&ctx, &bound);
  siqs_select_sieve_score_mode(&ctx);
  CHECK(ctx.params.sieve_saturating == !certified);
  CHECK(!ctx.params.sieve_saturating || ctx.params.sieve_raw_scores);
  *score_modes |= 1U << (ctx.params.sieve_saturating ? 2U
                        : ctx.params.sieve_raw_scores ? 1U : 0U);
  ctx.prime = siqs_malloc(siqs_array_size(ctx.params.fb_size, sizeof(uint32_t)));
  ctx.fb_in_a = siqs_calloc(ctx.params.fb_size, 1);
  ctx.root1 = siqs_malloc(siqs_array_size(ctx.params.fb_size, sizeof(uint32_t)));
  ctx.root2 = siqs_malloc(siqs_array_size(ctx.params.fb_size, sizeof(uint32_t)));
  ctx.sieve_length = 2U * ctx.params.half_interval;
  for (i = 0; i < ctx.params.fb_size; i++) ctx.prime[i] = ctx.fb[i].p;
  siqs_hashset_init(&ctx.a_hashes, 1024);
  siqs_poly_init(&ctx, &poly);
  score_reference_upper(upper, &ctx);
  for (family = 0; family < 2; family++) {
    CHECK(siqs_new_family(&ctx, &poly));
    for (polynomial = 0; polynomial < 3; polynomial++) {
      setup_check_polynomial(&ctx, &poly, upper, certified, bound);
      if (polynomial < 2) CHECK(siqs_next_B(&ctx, &poly));
    }
  }
  CHECK(ctx.sieve == NULL && ctx.candidate_at == NULL && ctx.candidate_at_wide == NULL);
  CHECK(ctx.raw_count == 0 && ctx.full_count == 0);
  if (detailed) {
    printf("PASS setup: bits=%u k=%u d=%u q=%u FB=%u P=%u M=%u\n",
        bits, k, ctx.params.poly_d, ctx.params.q_count, ctx.params.fb_size,
        ctx.largest_fb_prime, ctx.params.half_interval);
    printf("  LP=%llu R=%llu score=%s bound=%u payload=%.1f MiB CPU=%.2fs\n",
        (unsigned long long)ctx.params.large_prime_bound,
        (unsigned long long)ctx.params.smooth_bound,
        ctx.params.sieve_saturating ? "saturating" : ctx.params.sieve_raw_scores ? "raw" : "legacy",
        bound, (double)payload / (1024.0 * 1024.0),
        (double)(clock() - begin) / CLOCKS_PER_SEC);
    fflush(stdout);
  }
  siqs_poly_clear(&ctx, &poly);
  siqs_ctx_clear(&ctx); case_ctx = NULL;
  for (i = 0; i < result.count; i++) mpz_clear(result.values[i]);
  free(result.values); free(result.primality);
  mpz_clears(n, upper, NULL);
#ifndef _WIN32
  alarm(0);
#endif
}

/* Exhaust the raw streams, not millions of actual polynomial families. Five
 * finite evaluations uniquely identify an F_19 quartic; four identify a cubic.
 * Compare against explicit powers, independently of production Horner fitting. */
static void setup_check_cubic_stream(void) {
#if USE_CUBIC_CHOOSER
  unsigned test;
  for (test = 0; test < 2; test++) {
    siqs_ctx_t ctx;
    siqs_poly_t poly;
    siqs_factor_array_t result;
    siqs_A_search_t *g;
    siqs_A_cubic_t *c;
    uint32_t samples[19][4][SIQS_MAX_A_PRIMES];
    uint32_t frequency[SIQS_MAX_A_PRIMES][19] = {{0}};
    uint32_t prefix[SIQS_MAX_A_PRIMES], state, initial_state, draw, q, batches, key_words;
    unsigned char *seen;
    int have_best = 0;
    mpz_t n;
    case_name = "cubic-first-batch-compatibility-and-quartic-extension";
    mpz_init(n); setup_input(n, test == 0 ? 145 : 280, 1);
    siqs_factor_array_init(&result, n);
    siqs_ctx_init(&ctx, n, n, &result, NULL, 0); case_ctx = &ctx;
    ctx.multiplier = 1; ctx.params.poly_d = 2;
    mpz_set(ctx.kn, n);
    CHECK(siqs_ctx_allocate(&ctx) && !ctx.factor_found);
    siqs_poly_init(&ctx, &poly);
    poly.A_search = g = siqs_choose_A_prepare(&ctx, &poly, 1);
    /* These are raw codewords, not prime lookups: exercise F_7's infinity
     * coordinate and the maximum supported number of F_19 coordinates. */
    c = g->cubic; q = test == 0 ? 9U : SIQS_MAX_A_PRIMES;
    CHECK(c != NULL && !g->random_mode && q >= 5);
    CHECK(c->field == (test == 0 ? 7U : 19U));
    CHECK(c->period == (test == 0 ? 2401U : 130321U));
    CHECK(c->limit == (test == 0 ? 2401U : 2476099U));
    CHECK(c->batch == 0 && c->batch_end == c->period);
    batches = test == 0 ? 1U : 19U;
    key_words = test == 0 ? 4U : 5U;
    seen = (unsigned char *)check_allocate(c->limit, 1);
    initial_state = state = c->state;
    for (draw = 0; draw < c->limit; draw++) {
      uint32_t coefficients[4], value, key = 0, batch = draw/c->period, x;
      state = ((c->field+1U)*state+1U)%c->period;
      value = state;
      for (x = 0; x < 4; x++) {
        coefficients[x] = value%c->field; value /= c->field;
      }
      CHECK(siqs_A_cubic_prefix(g, q, prefix));
      CHECK(c->state == state && c->draws == draw+1U);
      CHECK(c->batch == batch && c->batch_end == (batch+1U)*c->period);
      for (x = 0; x < q-1U; x++) {
        uint32_t low = c->first+c->field*x, slot;
        uint32_t expected = x == c->field ? coefficients[3]
          : (coefficients[0]+coefficients[1]*x+coefficients[2]*x*x
             + coefficients[3]*x*x*x+batch*x*x*x*x)%c->field;
        CHECK(prefix[x] >= low && prefix[x] < low+c->field);
        slot = prefix[x]-low;
        CHECK(slot == expected);
        frequency[x][slot]++;
        if (x < key_words) key = key*c->field+slot;
      }
      CHECK(key < c->limit && !seen[key]); seen[key] = 1;
      if (draw%c->period < 4U)
        memcpy(samples[batch][draw%c->period], prefix, (q-1U)*sizeof(*prefix));
    }
    CHECK(!siqs_A_cubic_prefix(g, q, prefix));
    CHECK(!siqs_A_cubic_prefix(g, q, prefix) && c->draws == c->limit);
    for (uint32_t x = 0; x < q-1U; x++)
      for (uint32_t slot = 0; slot < c->field; slot++)
        CHECK(frequency[x][slot] == c->limit/c->field);
    for (uint32_t a = 0; a < batches*4U; a++)
      for (uint32_t b = a+1U; b < batches*4U; b++) {
        uint32_t shared = 0;
        for (uint32_t x = 0; x < q-1U; x++)
          shared += samples[a/4U][a%4U][x] == samples[b/4U][b%4U][x];
        CHECK(shared <= (a/4U == b/4U ? 3U : 4U));
      }
    /* Exhaustion still hands control to the normal fitting/SAR machinery. */
    (void)siqs_adaptive_A_candidate(&ctx, &poly, &have_best);
    CHECK(g->random_mode && c->draws == c->limit);
    /* Fit real A families at every added-batch boundary, retaining the
     * actual policy's q and geometry. No huge factoring run is required. */
    if (test != 0) {
      siqs_A_search_clear(&poly);
      poly.A_search = g = siqs_choose_A_prepare(&ctx, &poly, 1);
      c = g->cubic;
      for (uint32_t batch = 1; batch < 19U; batch++) {
        mpz_t product;
        c->draws = c->batch_end = batch*c->period;
        c->batch = batch-1U; c->state = initial_state;
        CHECK(siqs_new_family(&ctx, &poly));
        CHECK(!g->random_mode && c->batch == batch);
        mpz_init_set_ui(product, 1);
        for (uint32_t x = 0; x < poly.q_count; x++) {
          CHECK(poly.a_index[x] < ctx.params.fb_size);
          for (uint32_t y = 0; y < x; y++) CHECK(poly.a_index[x] != poly.a_index[y]);
          mpz_mul_ui(product, product, ctx.fb[poly.a_index[x]].p);
        }
        CHECK(mpz_cmp(product, poly.A) == 0);
        CHECK(mpz_cmp(poly.A, g->hard_lower) >= 0 &&
              mpz_cmp(poly.A, g->hard_upper) <= 0);
        mpz_mul(product, poly.B, poly.B); mpz_sub(product, product, ctx.kn);
        CHECK(mpz_divisible_p(product, poly.DA));
        mpz_divexact(product, product, poly.DA);
        CHECK(mpz_cmp(product, poly.C) == 0);
        mpz_clear(product);
      }
    }
    free(seen);
    siqs_poly_clear(&ctx, &poly); siqs_ctx_clear(&ctx);
    gmp_siqs_free(result.values, result.count); free(result.primality);
    mpz_clear(n); case_ctx = NULL;
  }
  puts("PASS setup: complete F_7/F_19 streams, first-batch compatibility, uniqueness/balance/overlap, SAR handoff and 18 real extension families");
#endif
}

static void setup_check_cubic_diagnostics(void) {
#if USE_CUBIC_CHOOSER && !defined(_WIN32)
  siqs_ctx_t ctx;
  siqs_poly_t poly;
  siqs_factor_array_t result;
  mpz_t n;
  unsigned reason, verbosity;
  case_name = "cubic-to-sar-diagnostics";
  mpz_init(n); setup_input(n, 145, 1);
  siqs_factor_array_init(&result, n);
  siqs_ctx_init(&ctx, n, n, &result, NULL, 0); case_ctx = &ctx;
  ctx.multiplier = 1; ctx.params.poly_d = 2;
  CHECK(siqs_build_factor_base(&ctx) && !ctx.factor_found);
  siqs_poly_init(&ctx, &poly);
  for (reason = 0; reason < 2; reason++) {
    for (verbosity = 0; verbosity < 2; verbosity++) {
      siqs_A_search_t *g;
      FILE *capture = tmpfile();
      char message[256];
      int saved_stdout, have_best = 0;
      CHECK(capture != NULL);
      poly.A_search = g = siqs_choose_A_prepare(&ctx, &poly, 1);
      CHECK(g->cubic != NULL && !g->random_mode);
      /* Reject every candidate without depending on a particular prime set. */
      mpz_set_ui(g->hard_lower, 1); mpz_set_ui(g->hard_upper, 0);
#if HAVE_UINT128 && ULONG_MAX > UINT32_MAX
      siqs_A_native_guard(g);
#endif
      if (reason == 0) g->cubic->draws = g->cubic->limit;
      ctx.verbose = (int)verbosity;
      fflush(stdout); saved_stdout = dup(STDOUT_FILENO);
      CHECK(saved_stdout >= 0 && dup2(fileno(capture), STDOUT_FILENO) >= 0);
      CHECK(!siqs_adaptive_A_candidate(&ctx, &poly, &have_best));
      CHECK(g->random_mode && !have_best);
      CHECK(g->cubic->draws == (reason == 0 ? g->cubic->limit
                                          : SIQS_A_CUBIC_RETRIES));
      /* Repeated SAR searches must not repeat the transition message. */
      CHECK(!siqs_adaptive_A_candidate(&ctx, &poly, &have_best));
      fflush(stdout); CHECK(dup2(saved_stdout, STDOUT_FILENO) >= 0);
      close(saved_stdout); rewind(capture);
      if (verbosity) {
        CHECK(fgets(message, sizeof(message), capture) != NULL);
        CHECK(strstr(message, "Cubic -> SAR:") != NULL);
        CHECK(strstr(message, reason == 0 ? "prefix period exhausted"
                        : "bounded search found no usable A") != NULL);
        CHECK(strstr(message, "prefix proposals used") != NULL);
      }
      CHECK(fgets(message, sizeof(message), capture) == NULL);
      fclose(capture); siqs_A_search_clear(&poly);
    }
  }
  ctx.verbose = 0;
  siqs_poly_clear(&ctx, &poly); siqs_ctx_clear(&ctx);
  gmp_siqs_free(result.values, result.count); free(result.primality);
  mpz_clear(n); case_ctx = NULL;
  puts("PASS setup: both Cubic-to-SAR reasons, proposal counts, one-time reporting and quiet mode");
#else
  puts("SKIP setup: Cubic diagnostic capture requires Cubic and POSIX descriptors");
#endif
}

static void suite_setup(void) {
  static const uint32_t sizes[] = {400, 401, 430, 431, 463, 496, 512};
  static const struct { uint32_t k, residue; } classes[] = {
    {1, 1}, {2, 3}, {2047, 1}, {2047, 7}
  };
  uint32_t i, j, d_modes = 0, score_modes = 0;
#ifndef _WIN32
  void (*previous_handler)(int);
#endif
  setup_check_cubic_stream();
  setup_check_cubic_diagnostics();
  if (!extended) {
    puts("SKIP setup: real upper-policy FB/A fixtures require --extended");
    return;
  }
  case_name = "production-upper-policy-setup";
#ifndef _WIN32
  previous_handler = signal(SIGALRM, setup_timeout);
  CHECK(previous_handler != SIG_ERR);
#else
  puts("SKIP setup timeout: per-case wall limit requires POSIX alarm");
#endif

  for (i = 0; i < sizeof(sizes) / sizeof(*sizes); i++)
    for (j = 0; j < sizeof(classes) / sizeof(*classes); j++)
      setup_real_case(sizes[i], classes[j].k, classes[j].residue,
                      &d_modes, &score_modes);
  setup_real_case(512, 4095, 7, &d_modes, &score_modes);

#ifndef _WIN32
  CHECK(signal(SIGALRM, previous_handler) != SIG_ERR);
#endif
  case_name = "production-upper-policy-setup";
  CHECK(d_modes == ((1U << 1) | (1U << 2)));
  printf("PASS setup: 29 production-policy fixtures, 58 A families, "
         "174 polynomials; both d modes, score-mode mask=%u; no sieving\n",
         score_modes);
}
