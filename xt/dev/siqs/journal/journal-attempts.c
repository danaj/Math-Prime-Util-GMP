/* Exercise completed unsuccessful/recovery attempts without a long failure.
 * This is a private integration test, never linked into the standalone tool. */
#define main journal_frontend_main
#include "../../../../siqs/mpu-siqs.c"
#undef main
#include "../../../../siqs/siqs.c"

static int check_codecs(void) {
  static const struct {
    uint32_t value;
    unsigned count;
    unsigned char bytes[5];
  } cases[] = {
    {0,1,{0}}, {1,1,{1}}, {127,1,{127}}, {128,2,{128,1}},
    {255,2,{255,1}}, {16383,2,{255,127}}, {16384,3,{128,128,1}},
    {2047,2,{255,15}}, {UINT32_C(4194303),4,{255,255,255,1}},
    {UINT32_MAX,5,{255,255,255,255,15}}
  };
  unsigned char encoded[8], *end;
  unsigned i, width;
  sj_reader r;
  for (i = 0; i < sizeof(cases)/sizeof(cases[0]); i++) {
    end = sj_put_varint(encoded, cases[i].value);
    if ((size_t)(end-encoded) != cases[i].count ||
        memcmp(encoded, cases[i].bytes, cases[i].count) != 0) return 10;
    r.data = encoded; r.left = cases[i].count;
    if (sj_get_varint(&r) != cases[i].value || r.left != 0) return 11;
  }
  for (width = 1; width <= 8; width++) {
    uint64_t value = UINT64_C(0xfedcba9876543210);
    uint64_t mask = UINT64_MAX >> (8U * (8U-width));
    end = sj_put_uint(encoded, value, width);
    if ((size_t)(end-encoded) != width) return 12;
    for (i = 0; i < width; i++)
      if (encoded[i] != (unsigned char)(value >> (8U*i))) return 13;
    r.data = encoded; r.left = width;
    if (sj_get_uint(&r,width) != (value&mask) || r.left != 0) return 14;
  }
  puts("PASS canonical ULEB128 boundary encodings and 1-8-byte LP fields");
  return 0;
}

/* A tiny valid relation stream exercises the production writer, reader and
 * graph acceptance without running a full factorization or linear algebra. */
static int check_events(const char *path, int solve) {
  siqs_journal_t *journal;
  siqs_factor_array_t result;
  siqs_ctx_t ctx;
  siqs_poly_t poly;
  siqs_factor_t factor = {1, 1};
  mpz_t n, y;
  uint64_t families = 0, polynomials = 0, p;
  int replay;
  mpz_init_set_ui(n, 1); mpz_mul_2exp(n, n, 89); mpz_sub_ui(n, n, 1);
  mpz_init_set_ui(y, 1); mpz_mul_2exp(y, y, 107); mpz_sub_ui(y, y, 1);
  mpz_mul(n, n, y);
  journal = gmp_siqs_journal_open(path, n); replay = journal->reading;
  siqs_factor_array_init(&result, n); result.journal = journal;
  siqs_ctx_init(&ctx, n, n, &result, NULL, 0);
  ctx.params.policy_name = "journal_event_codec";
  ctx.params.max_large_primes = 2;
  if (ctx.params.large_prime_bound > UINT32_MAX) return 20;
  ctx.params.smooth_bound = ctx.params.large_prime_bound * ctx.params.large_prime_bound;
  ctx.params.residual_multiplier = ctx.params.lp_multiplier * ctx.params.lp_multiplier;
  if (solve) { ctx.use_bids = 0; ctx.inline_matrix_solves = 0; }
  if (!siqs_ctx_allocate(&ctx)) return 21;
  siqs_poly_init(&ctx, &poly);
  ctx.matrix_next_target = ctx.params.target_relations;
  if (sj_attach(&ctx, &poly, &families, &polynomials)) return 22;
  if (!replay) {
    sj_begin_attempt(&ctx);
    if (!siqs_choose_A(&ctx, &poly)) return 23;
    siqs_journal_family(&ctx, &poly);
    mpz_set_ui(y, 1);
    siqs_accept_raw_relation(&ctx, siqs_raw_relation_new(&ctx, y, NULL, 0, 1, 1));
    mpz_set_ui(y, 2);
    siqs_accept_raw_relation(&ctx, siqs_raw_relation_new(&ctx, y, &factor, 1, 1, 2));
    p = (uint64_t)ctx.fb[ctx.params.fb_size - 1].p + 1U;
    while (!siqs_bpsw_u64(p)) p++;
    if (p > ctx.params.large_prime_bound || p > ULONG_MAX) return 24;
    mpz_set_ui(y, (unsigned long)p);
    siqs_accept_raw_relation(&ctx, siqs_raw_relation_new(&ctx, y, NULL, 0, p, p));
    ctx.total_candidates = 3;
    if (solve) {
      uint64_t sequence = journal->sequence;
      uint32_t k;
      for (k = 1; ctx.full_count < SIQS_MATRIX_EXTRA_RELS(&ctx); k++) {
        factor.exponent = 2U*k;
        mpz_set_ui(y, 1); mpz_mul_2exp(y, y, k);
        siqs_accept_raw_relation(&ctx, siqs_raw_relation_new(&ctx, y, &factor, 1, 1, 1));
        ctx.total_candidates++;
      }
      if (journal->batch.size == 0 || journal->sequence != sequence) return 25;
      (void)siqs_solve(&ctx);
      if (journal->batch.size != 0 || journal->sequence != sequence + 1U) return 26;
    } else sj_flush(&ctx);
  }
  printf("{\"full\":%u,\"smooth\":%llu,\"one_lp\":%llu,\"two_lp\":%llu,\"families\":%llu}\n",
    ctx.full_count, (unsigned long long)ctx.accepted_smooth,
    (unsigned long long)ctx.accepted_one_lp, (unsigned long long)ctx.accepted_two_lp,
    (unsigned long long)journal->families);
  siqs_poly_clear(&ctx, &poly); siqs_ctx_clear(&ctx);
  gmp_siqs_journal_close(journal);
  gmp_siqs_free(result.values, result.count); free(result.primality);
  mpz_clear(n); mpz_clear(y);
  return 0;
}

static int check_clock(const char *path) {
  siqs_journal_t *journal;
  siqs_ctx_t ctx;
  mpz_t n, values[2], *cached;
  uint32_t count;
  uint64_t origin;
  mpz_init_set_ui(n, 77);
  mpz_init_set_ui(values[0], 7); mpz_init_set_ui(values[1], 11);
  journal = gmp_siqs_journal_open(path, n);
  gmp_siqs_journal_close(journal);
  journal = gmp_siqs_journal_open(path, n);
  if (sj_peek(journal) || journal->writing_run || journal->clock_started) return 40;
  memset(&ctx, 0, sizeof(ctx)); ctx.journal = journal;
  sj_begin_attempt(&ctx);
  if (!journal->clock_started || journal->writing_run || journal->run != 0) return 41;
  /* Simulate elapsed work without sleeping or running a factor job. */
  journal->mono_start /= 2;
  origin = journal->mono_start;
  sj_begin_attempt(&ctx);
  if (journal->mono_start != origin) return 42;
  sj_finish_call(journal, n, values, 2);
  if (!journal->writing_run || journal->run != 1 || journal->mono_start != origin) return 43;
  sj_begin_attempt(&ctx);
  if (journal->mono_start != origin) return 44;
  gmp_siqs_journal_finish(journal, values, 2, 1);
  gmp_siqs_journal_close(journal);
  journal = gmp_siqs_journal_open(path, n);
  cached = sj_cached_call(journal, n, &count);
  if (cached == NULL || count != 2) return 45;
  gmp_siqs_free(cached, count);
  gmp_siqs_journal_finish(journal, values, 2, 1);
  if (journal->clock_started || journal->writing_run || journal->run != 1) return 46;
  gmp_siqs_journal_close(journal);
  mpz_clear(n); mpz_clear(values[0]); mpz_clear(values[1]);
  puts("PASS resumed clock starts before writing, retains its origin, and replay starts no run");
  return 0;
}

int main(int argc, char **argv) {
  siqs_journal_t *journal;
  siqs_factor_array_t result;
  mpz_t n, p;
  int attempt, was_replay, wide;
  if (argc == 3 && strcmp(argv[1], "--events") == 0) return check_events(argv[2], 0);
  if (argc == 3 && strcmp(argv[1], "--solve-events") == 0) return check_events(argv[2], 1);
  if (argc == 3 && strcmp(argv[1], "--clock") == 0) return check_clock(argv[2]);
  if (argc == 2 && strcmp(argv[1], "--codecs") == 0) return check_codecs();
  if (argc != 2 && (argc != 3 || strcmp(argv[2], "--wide") != 0)) return 2;
  wide = argc == 3;
  mpz_init_set_str(n, "1000000016000000063", 10);
  mpz_init_set_ui(p, 1000000007UL);
  if (wide) {
    /* Two known 256-bit primes: largest 256-bit prime and secp256k1's p.
     * Exercise real 512-bit setup/header/replay, without any sieving. */
    mpz_set_ui(n, 1); mpz_mul_2exp(n, n, 256); mpz_sub_ui(n, n, 189);
    mpz_set_ui(p, 1); mpz_mul_2exp(p, p, 256);
    mpz_sub_ui(p, p, 4294967295UL); mpz_sub_ui(p, p, 978);
    mpz_mul(n, n, p);
    if (mpz_sizeinbase(n, 2) != 512) return 7;
  }
  journal = gmp_siqs_journal_open(argv[1], n);
  was_replay = journal->reading;
  siqs_factor_array_init(&result, n); result.journal = journal;
  for (attempt = 0; attempt < (wide ? 1 : 2); attempt++) {
    siqs_ctx_t ctx;
    siqs_poly_t poly;
    uint64_t families = 0, polys = 0;
    int completed;
    siqs_ctx_init(&ctx, n, n, &result, NULL, 0);
    /* Give the second attempt different resolved interval geometry. */
    if (attempt) ctx.params.half_interval += 64;
    if (!siqs_ctx_allocate(&ctx)) return 3;
    siqs_poly_init(&ctx, &poly);
    completed = sj_attach(&ctx, &poly, &families, &polys);
    if (completed != was_replay) return 4;
    if (!completed) {
      sj_begin_attempt(&ctx);
      if (!wide) {
        if (!siqs_new_family(&ctx, &poly)) return 5;
        siqs_sieve_polynomial(&ctx, &poly); polys++;
      }
      if (attempt) {
        siqs_insert_divisor(&result, p);
        ctx.factor_found = 1;
      }
      sj_end_attempt(&ctx);
    }
    if (!!ctx.factor_found != !!attempt) return 6;
    siqs_poly_clear(&ctx, &poly); siqs_ctx_clear(&ctx);
  }
  sj_finish_call(journal, n, result.values, result.count);
  gmp_siqs_journal_close(journal);
  gmp_siqs_free(result.values, result.count); free(result.primality);
  mpz_clear(n); mpz_clear(p);
  puts(wide ? "PASS real 512-bit resolved parameters and ordered factor-base header" :
       "PASS unsuccessful attempt followed by different resolved recovery context");
  return 0;
}
