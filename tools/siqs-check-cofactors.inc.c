/* Cofactor acceptance/counter and nested-call checks, not a compilation unit. */
static void cofactor_context(siqs_ctx_t *ctx) {
  memset(ctx, 0, sizeof(*ctx));
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
  mpz_t divisor;
  uint64_t a, b;
  case_name = "split-methods-versus-failures-and-rejections";
  cofactor_miss_native = 1;
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
  cofactor_context(&ctx);
  CHECK(cofactor_resolve(&ctx, "10201")); /* 101^2 */
  CHECK(ctx.split_square == 1);
  cofactor_context(&ctx);
  CHECK(cofactor_resolve(&ctx, "10403")); /* 101*103 */
  CHECK(ctx.split_squfof == 1);
  cofactor_context(&ctx); ctx.params.large_prime_bound = 103;
  CHECK(!cofactor_resolve(&ctx, "10807")); /* 101*107: split, LP-invalid */
  CHECK(ctx.split_rejected == 1 && ctx.split_fail == 0 && ctx.split_squfof == 0);
  cofactor_context(&ctx); ctx.params.large_prime_bound = 100;
  CHECK(!cofactor_resolve(&ctx, "10201"));
  CHECK(ctx.split_rejected == 1 && ctx.split_square == 0);
  cofactor_context(&ctx);
  CHECK(!cofactor_resolve(&ctx, "1113121")); /* 101*103*107, not a prime pair */
  CHECK(ctx.split_rejected == 1);
  cofactor_miss_squfof = 1;
  cofactor_miss_siqs = 1;
  cofactor_context(&ctx);
  CHECK(cofactor_resolve(&ctx, "47053")); /* 211*223, beyond inner trial limit */
  CHECK(ctx.split_rho == 1);
  cofactor_miss_rho = 1;
  cofactor_context(&ctx);
  CHECK(!cofactor_resolve(&ctx, "47053"));
  CHECK(ctx.split_fail == 1 && ctx.split_rejected == 0);
  cofactor_miss_rho = cofactor_miss_squfof = cofactor_miss_siqs = 0;
  cofactor_context(&ctx);
  CHECK(cofactor_resolve(&ctx, "1") && ctx.split_attempts == 0);
  CHECK(!cofactor_resolve(&ctx, "18446744073709551616") && ctx.split_attempts == 0);
  ctx.params.large_prime_bound = 100;
  CHECK(!cofactor_resolve(&ctx, "107") && ctx.split_attempts == 0);
  ctx.params.smooth_bound = 5000;
  CHECK(!cofactor_resolve(&ctx, "10403") && ctx.split_attempts == 0);
  cofactor_miss_native = 0;
  puts("PASS cofactors: accepted methods, bounded misses, invalid pairs, and counter conservation");
}

static void cofactor_trace_reset(void) {
  cofactor_squfof_calls = cofactor_prime_calls = cofactor_rho_calls = 0;
}
static void cofactor_cascade(void) {
  siqs_ctx_t ctx;
  case_name = "cofactor-preferred-method-and-fallback-order";
  cofactor_miss_native = cofactor_trace = 1;
  cofactor_trace_reset(); cofactor_context(&ctx);
  CHECK(cofactor_resolve(&ctx, "47053"));
  CHECK(ctx.split_squfof == 1 && cofactor_squfof_calls == 1 &&
        cofactor_prime_calls == 0 && cofactor_rho_calls == 0);
  cofactor_miss_squfof = 1;
  cofactor_trace_reset(); cofactor_context(&ctx);
  CHECK(cofactor_resolve(&ctx, "47053"));
  CHECK(cofactor_squfof_calls == 1);
  CHECK(ctx.split_siqs == 1 && cofactor_prime_calls != 0 && cofactor_rho_calls == 0);
  cofactor_miss_siqs = 1;
  cofactor_trace_reset(); cofactor_context(&ctx);
  CHECK(cofactor_resolve(&ctx, "47053"));
  CHECK(ctx.split_rho == 1 && cofactor_squfof_calls == 1 && cofactor_rho_calls == 1);
  CHECK(cofactor_prime_calls == 1);

  /* Easy rho recovery, but large enough that enabled SIQS skips SQUFOF. */
  cofactor_miss_squfof = 0;
  cofactor_trace_reset(); cofactor_context(&ctx);
  ctx.largest_fb_prime = 397; ctx.params.large_prime_bound = UINT64_MAX;
  /* 1009 * 9007199254740881: both prime; skip the inner trial limit. */
  CHECK(cofactor_resolve(&ctx, "9088264048033548929"));
  CHECK(ctx.split_rho == 1 && cofactor_squfof_calls == 0 &&
        cofactor_prime_calls == 1 && cofactor_rho_calls == 1);
  cofactor_trace_reset(); cofactor_context(&ctx); ctx.params.bits = 64;
  ctx.largest_fb_prime = 397; ctx.params.large_prime_bound = UINT64_MAX;
  CHECK(cofactor_resolve(&ctx, "9088264048033548929"));
  CHECK(ctx.split_siqs == 0 && cofactor_squfof_calls == 1 && cofactor_prime_calls == 0);
  cofactor_miss_siqs = cofactor_miss_native = cofactor_trace = 0;
  puts("PASS cofactors: SQUFOF/SIQS crossover, low-bit SIQS recovery, rho fallback, and recursion guard");
}
#ifndef _WIN32
static void cofactor_miss_notice(void) {
  unsigned int verbose;
  siqs_ctx_t ctx;
  mpz_t n;
  case_name = "cofactor-siqs-miss-diagnostic";
  cofactor_miss_native = cofactor_miss_squfof = cofactor_miss_siqs = 1;
  mpz_init_set_ui(n, 47053);
  for (verbose = 0; verbose <= 1; verbose++) {
    FILE *capture = tmpfile();
    char message[256];
    int saved_stderr, accepted;
    uint64_t a, b;
    CHECK(capture != NULL); fflush(stderr);
    saved_stderr = dup(STDERR_FILENO); CHECK(saved_stderr >= 0);
    CHECK(dup2(fileno(capture), STDERR_FILENO) >= 0);
    cofactor_context(&ctx); ctx.verbose = (int)verbose;
    accepted = siqs_resolve_cofactor(&ctx, n, &a, &b);
    fflush(stderr);
    CHECK(dup2(saved_stderr, STDERR_FILENO) >= 0); close(saved_stderr);
    CHECK(accepted && ctx.split_rho == 1 && a == 211 && b == 223);
    cofactor_counts(&ctx);
    rewind(capture);
    if (verbose) {
      CHECK(fgets(message, sizeof(message), capture) != NULL);
      CHECK(strstr(message, "failed to split 47053 (16 bits); trying rho") != NULL);
    }
    CHECK(fgets(message, sizeof(message), capture) == NULL);
    fclose(capture);
  }
  mpz_clear(n);
  cofactor_miss_native = cofactor_miss_squfof = cofactor_miss_siqs = 0;
  puts("PASS cofactors: SIQS miss reports residual at verbose level 1+, quiet level stays silent");
}
#endif

/* Host verbosity is deliberately high: the complete inner call, including
 * reduction and any Lanczos fallback, must remain quiet without mutating it. */
static void cofactor_nested_case(uint32_t bits) {
  siqs_ctx_t ctx;
  mpz_t n;
  uint64_t a, b, input;
  uint32_t i;
  static const char *values[] = {
    "36028778899572911", /* 134217689 * 268435399, 55 bits */
    "72057554846356433", /* 268435399 * 268435367, 56 bits */
    "144115156668907691", /* 268435399 * 536870909, 57 bits */
    "288230356824359011", /* 536870909 * 536870879, 58 bits */
    "576460727070490747", /* 536870909 * 1073741783, 59 bits */
    "18446743979220271189" /* 4294967291 * 4294967279, 64 bits */
  };
  for (i = 0; i < sizeof(values) / sizeof(values[0]); i++) {
    cofactor_context(&ctx); ctx.params.bits = bits;
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
static void suite_cofactors(void) {
  int saved_verbose = verbose_level;
#ifndef _WIN32
  FILE *capture;
  int saved_stdout;
#endif
  case_name = "nested-siqs-ownership-reentrancy-and-quiet-output";
  prime_iterator_global_startup();
  cofactor_classification();
  cofactor_cascade();
#ifndef _WIN32
  cofactor_miss_notice();
#endif
  case_name = "nested-siqs-ownership-reentrancy-and-quiet-output";
  cofactor_miss_native = 1;
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
  verbose_level = saved_verbose; cofactor_miss_native = 0;
  prime_iterator_global_shutdown();
  puts("PASS cofactors: largest-two helper, 55-59/64-bit cases, wide input, overflow, and quiet/reentrant calls");
}
