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

static void suite_setup(void) {
  static const uint32_t sizes[] = {400, 401, 430, 431, 463, 496, 512};
  static const struct { uint32_t k, residue; } classes[] = {
    {1, 1}, {2, 3}, {2047, 1}, {2047, 7}
  };
  uint32_t i, j, d_modes = 0, score_modes = 0;
#ifndef _WIN32
  void (*previous_handler)(int);
#endif
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
