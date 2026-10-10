/* Standalone SIQS sanity/regression checker, not part of Perl's test suite.
 *
 * Direct build from the repository root (GMP required, Perl not required):
 *   cc -O3 -DSTANDALONE -o /tmp/siqs-check tools/siqs-check.c \
 *     prime_iterator.c squfof126.c pbrent63.c -lgmp -lm
 *   /tmp/siqs-check --suite sieve
 *
 * Include the current implementation to test its private interfaces without
 * exporting them or adding production testing hooks.  Keep reference code
 * independent of production sieve/scan kernels.  New suites can be added to
 * the small dispatcher below; no full factorization is needed here.
 *
 * Copyright (c) 2026 Dana Jacobsen.
 */
#include <stdarg.h>
#ifndef _WIN32
# include <unistd.h>
# include <sys/wait.h>
#endif
#define main siqs_check_unused_driver_main
#include "../mpu-siqs.c"
#undef main
#include "siqs-check-thread-seams.h"
#include "siqs-check-cofactor-seams.h"
#ifndef SIQS_CHECK_SOURCE
# define SIQS_CHECK_SOURCE "../siqs.c"
#endif
#include SIQS_CHECK_SOURCE
#include "../mont64.h"
#undef squfof126
#undef uvpbrent63
#undef siqs_pbrent_factor
#undef siqs_is_prob_prime
#undef malloc
#undef calloc
#undef pthread_create
#undef pthread_join
#undef pthread_mutex_init
#undef pthread_mutex_destroy
#undef pthread_cond_init
#undef pthread_cond_destroy
#undef pthread_mutex_unlock
#undef pthread_cond_signal
#undef fprintf
#undef exit

static int extended, detailed;
static const char *suite_name = "startup";
static const char *case_name = "startup";
static const siqs_ctx_t *case_ctx;
static uint32_t comparisons, maximum_score, d_seen, special_a, special_k;

static void check_fail(const char *format, ...) {
  va_list args;
  fprintf(stderr, "FAIL %s/%s: ", suite_name, case_name);
  va_start(args, format);
  vfprintf(stderr, format, args);
  va_end(args);
  fputc('\n', stderr);
  if (case_ctx != NULL)
    fprintf(stderr, "  bits=%u k=%lu d=%u FB=%u length=%u first=%u "
            "block=%u dense-primes=%u\n", case_ctx->params.bits,
            case_ctx->multiplier, case_ctx->params.poly_d,
            case_ctx->params.fb_size, case_ctx->sieve_length,
            case_ctx->params.sieve_start, case_ctx->block_length,
            case_ctx->block_prime_count);
  exit(1);
}

#define CHECK(condition) do { \
  if (!(condition)) check_fail("%s:%u: %s", __FILE__, \
                               (unsigned)__LINE__, #condition); \
} while (0)

static void *check_allocate(size_t count, size_t size) {
  void *p;
  CHECK(size == 0 || count <= (size_t)-1 / size);
  p = calloc(count != 0 ? count : 1, size);
  CHECK(p != NULL);
  return p;
}

/* Mirrors the allocation contract, not the store kernels.  Only logical
 * scores are compared: padding may intentionally wrap and is never scanned. */
static size_t check_sieve_allocation(const siqs_ctx_t *ctx) {
  size_t length = 2U * (size_t)ctx->sieve_length;
  size_t pmax = (size_t)ctx->prime[ctx->params.fb_size - 1U] + 1U;
  return (length > pmax ? length : pmax) + 8U;
}

/* A deliberately boring reference: walk each canonical progression and add
 * to a wide logical cell.  Do not reuse production tier dispatch, padded
 * stores, block-root advancement, byte addition, or candidate scanning. */
static uint32_t *reference_sieve(const siqs_ctx_t *ctx, uint32_t *maximum) {
  uint32_t i, pos, initial = (uint32_t)ctx->active_sieve_initial
                          + ctx->params.stage1_bias;
  uint64_t ceiling = initial;
  uint32_t *score = (uint32_t *)check_allocate(ctx->sieve_length,
                                               sizeof(uint32_t));
  for (pos = 0; pos < ctx->sieve_length; pos++)
    score[pos] = initial;
  for (i = ctx->params.sieve_start; i < ctx->params.fb_size; i++) {
    size_t at;
    uint32_t p = ctx->prime[i], a = ctx->root1[i], b = ctx->root2[i];
    CHECK(p != 0 && a <= b && b < p);
    /* Each prime contributes at most once to a particular logical cell. */
    ceiling += ctx->sieve_logp[i];
    CHECK(ceiling <= UINT32_MAX);
    for (at = a; at < ctx->sieve_length; at += p)
      score[at] += ctx->sieve_logp[i];
    if (b != a)
      for (at = b; at < ctx->sieve_length; at += p)
        score[at] += ctx->sieve_logp[i];
  }
  *maximum = initial;
  for (pos = 0; pos < ctx->sieve_length; pos++)
    if (score[pos] > *maximum)
      *maximum = score[pos];
  return score;
}

static void compare_sieve(siqs_ctx_t *ctx) {
  uint32_t pos, count = 0, peak, i;
  size_t roots_size = (size_t)ctx->params.fb_size * sizeof(uint32_t);
  uint32_t *roots = (uint32_t *)check_allocate(2, roots_size);
  uint32_t *score;
  case_ctx = ctx;
  CHECK(ctx->params.sieve_start <= ctx->params.fb_size);
  CHECK(ctx->block_prime_count <= ctx->params.fb_size - ctx->params.sieve_start);
  for (i = 0; i < ctx->block_prime_count; i++)
    CHECK(ctx->block_step[i] ==
          ctx->block_length % ctx->prime[ctx->params.sieve_start + i]);
  score = reference_sieve(ctx, &peak);
  if (peak > UINT8_MAX)
    for (pos = 0; pos < ctx->sieve_length; pos++)
      if (score[pos] > UINT8_MAX)
        check_fail("reference score overflow at pos=%u: score=%u (>255)",
                   pos, score[pos]);
  if (peak > maximum_score) maximum_score = peak;
  memcpy(roots, ctx->root1, roots_size);
  memcpy((unsigned char *)roots + roots_size, ctx->root2, roots_size);
  /* The initializer must erase stale bytes; padding values must not matter. */
  memset(ctx->sieve, 0xd3, check_sieve_allocation(ctx));
  siqs_run_sieve(ctx);
  for (pos = 0; pos < ctx->sieve_length; pos++)
    if (score[pos] != ctx->sieve[pos])
      check_fail("byte mismatch at pos=%u: reference=%u production=%u",
                 pos, score[pos], ctx->sieve[pos]);
  for (i = 0; i < ctx->params.fb_size; i++)
    if (roots[i] != ctx->root1[i] ||
        roots[ctx->params.fb_size + i] != ctx->root2[i])
      check_fail("original roots changed at FB index=%u: (%u,%u) -> (%u,%u)",
                 i, roots[i], roots[ctx->params.fb_size + i],
                 ctx->root1[i], ctx->root2[i]);

  siqs_find_candidates(ctx);
  for (pos = 0; pos < ctx->sieve_length; pos++) {
    uint32_t map = ctx->candidate_wide ? ctx->candidate_at_wide[pos]
                                     : ctx->candidate_at[pos];
    if (score[pos] >= 128U) {
      if (count >= ctx->candidate_count)
        check_fail("missing candidate at pos=%u index=%u score=%u", pos,
                   count, score[pos] - ctx->params.stage1_bias);
      if (ctx->candidates[count].x !=
          (int32_t)pos - (int32_t)ctx->params.half_interval ||
          ctx->candidates[count].sieve_score !=
          score[pos] - ctx->params.stage1_bias)
        check_fail("candidate mismatch at pos=%u index=%u: x=%d score=%u",
                   pos, count, ctx->candidates[count].x,
                   (unsigned)ctx->candidates[count].sieve_score);
      CHECK(ctx->candidates[count].first_hit == SIQS_NO_INDEX);
      if (map != count + 1U)
        check_fail("candidate map mismatch at pos=%u: reference=%u actual=%u",
                   pos, count + 1U, map);
      count++;
    } else {
      if (map != 0)
        check_fail("unexpected candidate map entry at pos=%u: %u", pos, map);
    }
  }
  if (count != ctx->candidate_count)
    check_fail("candidate count mismatch: reference=%u actual=%u", count,
               ctx->candidate_count);
  CHECK(ctx->candidate_wide == (count > UINT16_MAX));
  siqs_clear_candidate_map(ctx);
  CHECK(ctx->candidate_count == 0 && !ctx->candidate_wide);
  for (pos = 0; pos < ctx->sieve_length; pos++) {
    CHECK(ctx->candidate_at[pos] == 0);
    if (ctx->candidate_at_wide != NULL)
      CHECK(ctx->candidate_at_wide[pos] == 0);
  }
  comparisons++;
  free(roots);
  free(score);
  case_ctx = NULL;
}

static int period_compare(const void *va, const void *vb) {
  uint32_t a = *(const uint32_t *)va, b = *(const uint32_t *)vb;
  return a < b ? -1 : a > b;
}

static void synthetic_init(siqs_ctx_t *ctx, uint32_t length,
                            const uint32_t *periods, uint32_t count) {
  uint32_t i;
  memset(ctx, 0, sizeof(*ctx));
  ctx->sieve_length = length;
  ctx->params.half_interval = length / 2;
  ctx->params.fb_size = count;
  ctx->prime = (uint32_t *)check_allocate(count, sizeof(uint32_t));
  memcpy(ctx->prime, periods, (size_t)count * sizeof(uint32_t));
  qsort(ctx->prime, count, sizeof(uint32_t), period_compare);
  for (i = 0; i < count; i++) CHECK(ctx->prime[i] != 0);
  ctx->root1 = (uint32_t *)check_allocate(count, sizeof(uint32_t));
  ctx->root2 = (uint32_t *)check_allocate(count, sizeof(uint32_t));
  ctx->sieve_logp = (uint8_t *)check_allocate(count, 1);
  ctx->sieve = (uint8_t *)check_allocate(check_sieve_allocation(ctx), 1);
  ctx->candidate_at = (uint16_t *)check_allocate(length, sizeof(uint16_t));
  siqs_block_workspace_allocate(ctx);
}

static void synthetic_clear(siqs_ctx_t *ctx) {
  free(ctx->prime); free(ctx->root1); free(ctx->root2); free(ctx->sieve_logp);
  free(ctx->sieve); free(ctx->candidate_at); free(ctx->candidate_at_wide);
  free(ctx->candidates); free(ctx->block_root1); free(ctx->block_root2);
  free(ctx->block_step);
}

static void synthetic_roots(siqs_ctx_t *ctx, uint32_t round) {
  uint32_t i;
  for (i = 0; i < ctx->params.fb_size; i++) {
    uint32_t p = ctx->prime[i];
    uint32_t a = (uint32_t)((UINT64_C(1729) * (i + 1U) + 17U * round) % p);
    uint32_t b = (uint32_t)((UINT64_C(7919) * (i + 3U) + 31U * round) % p);
    if ((i + round) % 4U == 0) b = a;
    ctx->root1[i] = a < b ? a : b;
    ctx->root2[i] = a < b ? b : a;
    ctx->sieve_logp[i] = 1;
  }
}

static void check_tier_boundaries(void) {
  uint32_t block = SIQS_SIEVE_BLOCK_SIZE ? SIQS_SIEVE_BLOCK_SIZE : 65536U;
  uint32_t minimum = SIQS_SIEVE_BLOCK_SIZE ? SIQS_SIEVE_BLOCK_MIN_LENGTH
                                        : 5U * block / 2U;
  uint32_t lengths[] = {31, 32, 33, 63, 64, 65, 8191, 8192, 8193,
    65535, 65536, 65537, minimum - 1U, minimum, minimum + 1U,
    3U * block + 17U};
  uint32_t n;
  case_name = "tier-boundaries";
  for (n = 0; n < sizeof(lengths) / sizeof(*lengths); n++) {
    uint32_t length = lengths[n], periods[100], count = 0, divisor, round;
    uint32_t span = length;
    siqs_ctx_t ctx;
    if (SIQS_SIEVE_BLOCK_SIZE && length > block) {
      uint32_t blocks = (length - 1U) / block + 1U;
      span = length / blocks + (length % blocks != 0);
    }
    periods[count++] = 2; periods[count++] = 3; periods[count++] = 7;
    /* Both whole-interval and balanced-block fixed-hit cutoffs. */
    for (divisor = 1; divisor <= 6; divisor++) {
      uint32_t cutoff[2], j;
      cutoff[0] = length / divisor; cutoff[1] = span / divisor;
      for (j = 0; j < 2; j++) {
        if (cutoff[j] > 1) periods[count++] = cutoff[j] - 1U;
        if (cutoff[j] != 0) periods[count++] = cutoff[j];
        periods[count++] = cutoff[j] + 1U;
      }
    }
    periods[count++] = 2U * length + 1U;
    synthetic_init(&ctx, length, periods, count);
    for (round = 0; round < 5; round++) {
      uint32_t i;
      synthetic_roots(&ctx, round);
      if (round >= 3)
        for (i = 0; i < count; i++) {
          ctx.root1[i] = round == 3 ? 0 : ctx.prime[i] - 1U;
          ctx.root2[i] = i % 4U == 0 ? ctx.root1[i] : ctx.prime[i] - 1U;
        }
      ctx.active_sieve_initial = round == 1 ? 120 : 100;
      ctx.params.stage1_bias = 8;
      ctx.params.sieve_start = round >= 2 ? count / 3 : 0;
      /* Rebuild the block split after changing the first-sieved index. */
      if (round == 2) {
        free(ctx.block_root1); free(ctx.block_root2); free(ctx.block_step);
        ctx.block_root1 = ctx.block_root2 = ctx.block_step = NULL;
        ctx.block_prime_count = 0;
        siqs_block_workspace_allocate(&ctx);
      }
      compare_sieve(&ctx);
    }
    synthetic_clear(&ctx);
  }
  puts("PASS sieve: tier boundaries, scan tails, balanced blocks, scratch reuse");
  fflush(stdout);
}

static void check_two_hit_local(void) {
  static const uint32_t lengths[] = {65535, 65536, 65537, 131073, 196625};
  uint32_t n;
  case_name = "two-hit/local";
  for (n = 0; n < sizeof(lengths) / sizeof(*lengths); n++) {
    uint32_t length = lengths[n], count = length - length / 2U;
    uint32_t *periods, i, round;
    siqs_ctx_t ctx;
    if (count > 4096U) count = 4096U;
    periods = (uint32_t *)check_allocate(count + 1U, sizeof(uint32_t));
    periods[0] = 7; /* Exercise a dense prefix when blocking is active. */
    for (i = 0; i < count; i++) periods[i + 1U] = length / 2U + 1U + i;
    synthetic_init(&ctx, length, periods, count + 1U);
    free(periods);
    for (round = 0; round < 3U; round++) {
      synthetic_roots(&ctx, round);
      if (round == 1U) {
        /* Mostly missing final hits, with both distinct and equal roots.
         * More than 2048 primes cycles the 4 KiB sink stripe. */
        for (i = 1; i <= count; i++) {
          uint32_t p = ctx.prime[i];
          ctx.root1[i] = i % 4U == 0 ? p - 1U : p - 2U;
          ctx.root2[i] = p - 1U;
        }
      }
      if (round != 0) {
        /* Adjacent final hits at length-1 and length: one logical, one
         * redirected.  Only one prime uses these roots to avoid overflow. */
        uint32_t root = length - ctx.prime[1];
        ctx.root1[1] = root - 1U;
        ctx.root2[1] = root;
      }
      ctx.active_sieve_initial = 96;
      ctx.params.stage1_bias = 8;
      compare_sieve(&ctx);
    }
    synthetic_clear(&ctx);
  }
  puts("PASS sieve: two-hit localization, exact-end hits, duplicates, sink wrap");
  fflush(stdout);
}

static void check_one_hit_and_maps(void) {
  uint32_t counts[] = {SIQS_ONE_HIT_LOCAL_MIN_PRIMES - 1U,
    SIQS_ONE_HIT_LOCAL_MIN_PRIMES, SIQS_ONE_HIT_LOCAL_MIN_PRIMES + 1U,
    8193, 65535, 65536, 65537};
  uint32_t n, limit = extended ? 7U : 4U;
  case_name = "one-hit/maps";
  for (n = 0; n < limit; n++) {
    uint32_t count = counts[n], length = n < 3 ? 8192U : 131073U;
    uint32_t *periods = (uint32_t *)check_allocate(count, sizeof(uint32_t));
    uint32_t i, round;
    siqs_ctx_t ctx;
    for (i = 0; i < count; i++) periods[i] = length + 1U + i;
    synthetic_init(&ctx, length, periods, count);
    free(periods);
    for (round = 0; round < 3; round++) {
      synthetic_roots(&ctx, round);
      ctx.params.stage1_bias = 8;
      ctx.active_sieve_initial = round == 1 ? 120 : 96;
      compare_sieve(&ctx);
    }
    synthetic_clear(&ctx);
  }
  puts("PASS sieve: localized one-hit threshold and narrow/wide candidate maps");
  fflush(stdout);
}

static void check_oracle_overflow(void) {
  siqs_ctx_t ctx;
  uint32_t p = 7, peak, *score;
  case_name = "oracle-overflow-self-check";
  synthetic_init(&ctx, 33, &p, 1);
  ctx.active_sieve_initial = 127;
  ctx.params.stage1_bias = 127;
  ctx.sieve_logp[0] = 2;
  /* Equal roots must contribute once: 254+2 == 256, not 258 or zero. */
  score = reference_sieve(&ctx, &peak);
  CHECK(peak == 256 && score[0] == 256);
  CHECK(score[1] == 254);
  CHECK(score[0] >= 128 && !((uint8_t)score[0] & 0x80U));
  free(score);
  synthetic_clear(&ctx);
  puts("PASS sieve: wide reference detects overflow hidden by byte wrapping");
}

/* These helpers build genuine polynomials but never collect or solve a full
 * factorization.  Small fixtures must survive the public trial-division
 * pretest and must not find their factor while constructing the base. */
static void make_semiprime(mpz_t n, gmp_randstate_t random, uint32_t bits) {
  mpz_t p, q;
  mpz_init(p); mpz_init(q);
  do {
    mpz_urandomb(p, random, bits / 2);
    mpz_setbit(p, bits / 2 - 1U); mpz_nextprime(p, p);
    mpz_urandomb(q, random, bits - bits / 2);
    mpz_setbit(q, bits - bits / 2 - 1U); mpz_nextprime(q, q);
    mpz_mul(n, p, q);
  } while (mpz_sizeinbase(n, 2) != bits || mpz_cmp(p, q) == 0);
  mpz_clear(p); mpz_clear(q);
}

static void check_trial_survivor(const mpz_t n) {
  PRIME_ITERATOR(iter);
  UV p;
  uint32_t bits = (uint32_t)mpz_sizeinbase(n, 2);
  uint32_t limit = bits <= 40 ? 200U : bits >= 1000 ? 5000U : 5U * bits;
  prime_iterator_setprime(&iter, 1);
  for (p = prime_iterator_next(&iter); p < limit;
       p = prime_iterator_next(&iter))
    CHECK(!mpz_divisible_ui_p(n, p));
  prime_iterator_destroy(&iter);
}

static void check_polynomial_roots(siqs_ctx_t *ctx, const siqs_poly_t *poly) {
  uint32_t i, r;
  mpz_t x, value;
  mpz_init(x); mpz_init(value);
  for (i = 0; i < ctx->params.fb_size; i++) {
    CHECK(ctx->root1[i] <= ctx->root2[i] && ctx->root2[i] < ctx->prime[i]);
    if (ctx->fb[i].in_a) special_a++;
    if (i != 0 && ctx->fb[i].sqrt_kn == 0) special_k++;
    for (r = 0; r < (ctx->root1[i] == ctx->root2[i] ? 1U : 2U); r++) {
      uint32_t root = r == 0 ? ctx->root1[i] : ctx->root2[i];
      CHECK(root <= INT32_MAX);
      mpz_set_si(x, (int32_t)root - (int32_t)ctx->params.half_interval);
      mpz_mul(value, poly->DA, x);
      mpz_add(value, value, poly->B);
      mpz_mul(value, value, value);
      mpz_sub(value, value, ctx->kn);
      CHECK(mpz_divisible_p(value, poly->DA));
      mpz_divexact(value, value, poly->DA);
      if (!mpz_divisible_ui_p(value, ctx->prime[i]))
        check_fail("polynomial root mismatch: FB index=%u prime=%u root=%u",
                   i, ctx->prime[i], root);
    }
  }
  mpz_clear(x); mpz_clear(value);
}

static void check_A_search_cache(void) {
  static const uint32_t primes[] = {2, 3, 5, 7, 11, 13};
  static const struct {
    uint32_t q, target, stage, tolerance;
  } windows[] = {
    /* With eligible 3/5/7: q2 products span 15..35; q3 is exactly 105. */
    {2, 8, 1, 8}, {2, 4, 2, 8}, {2, 2, 3, 8}, {2, 1, 0, 8},
    {2, 70, 1, 8}, {2, 71, 2, 8}, {2, 140, 2, 8},
    {2, 141, 3, 8}, {2, 280, 3, 8}, {2, 281, 0, 8},
    {3, 53, 1, 8}, {3, 52, 2, 8}, {3, 27, 2, 8},
    {3, 26, 3, 8}, {3, 14, 3, 8}, {3, 13, 0, 8},
    {3, 210, 1, 8}, {3, 211, 2, 8}, {3, 420, 2, 8},
    {3, 421, 3, 8}, {3, 840, 3, 8}, {3, 841, 0, 8},
    {1, 0, 1, 8}, {4, 0, 1, 8}, {2, 2, 0, 4}, {3, 14, 0, 4}
  };
  siqs_ctx_t ctx;
  siqs_poly_t poly;
  siqs_fb_t fb[6];
  uint32_t eligible, occupied, center, i;
  size_t fixture;
  case_name = "A-search-cache";
  memset(&ctx, 0, sizeof(ctx)); memset(fb, 0, sizeof(fb));
  ctx.fb = fb; ctx.params.fb_size = 6U;
  case_ctx = &ctx;
  for (i = 0; i < 6U; i++) fb[i].p = primes[i];
  for (eligible = 1U; eligible < 32U; eligible++) {
    for (occupied = 0U; occupied < 32U; occupied++) {
      if ((eligible & ~occupied) == 0U) continue;
      for (i = 1U; i < 6U; i++) {
        fb[i].sqrt_kn = (eligible >> (i - 1U)) & 1U;
        fb[i].in_a = (occupied >> (i - 1U)) & 1U;
      }
      for (center = 1U; center < 6U; center++) {
        uint32_t step, expected = SIQS_NO_INDEX;
        /* Independent old scan: upward wins every equal-distance tie. */
        for (step = 0U; step < 6U; step++) {
          uint32_t up = center + step, down = center >= step ? center - step : 0U;
          if (up < 6U && fb[up].sqrt_kn && !fb[up].in_a) {
            expected = up; break;
          }
          if (step && down && fb[down].sqrt_kn && !fb[down].in_a) {
            expected = down; break;
          }
        }
        CHECK(expected != SIQS_NO_INDEX);
        CHECK(siqs_nearest_available_fb_from_index(&ctx, center) == expected);
        CHECK(siqs_nearest_available_fb(&ctx, fb[center].p) == expected);
      }
    }
  }
  ctx.params.fb_size = 4U;
  ctx.params.half_interval = ctx.params.poly_d = 1U;
  mpz_init_set_ui(ctx.kn, 100U);
  for (fixture = 0; fixture < sizeof(windows) / sizeof(*windows); fixture++) {
    uint32_t repeat;
    for (i = 1U; i < 4U; i++) { fb[i].sqrt_kn = 1U; fb[i].in_a = 0U; }
    ctx.params.q_count = windows[fixture].q;
    ctx.params.a_final_tolerance = windows[fixture].tolerance;
    siqs_poly_init(&ctx, &poly);
    CHECK(!poly.a_search_ready && !poly.a_search_stage &&
          !poly.a_search_center && !poly.a_search_variance);
    /* Target and policy are finalized before the first lazy preparation. */
    mpz_set_ui(poly.target_A, windows[fixture].target);
    if (windows[fixture].stage != 0U)
      siqs_prepare_A_search(&ctx, &poly);
    else {
      poly.a_index[0] = 1U;
      for (repeat = 0U; repeat < 2U; repeat++) {
        fb[1].in_a = 1U;
        CHECK(!siqs_choose_A(&ctx, &poly));
        CHECK(!fb[1].in_a && poly.a_search_ready && !poly.a_search_stage &&
              !poly.a_search_center && !poly.a_search_variance);
      }
    }
    CHECK(poly.a_search_ready && poly.a_search_stage == windows[fixture].stage);
    if (poly.a_search_stage != 0U)
      CHECK(poly.a_search_center > 0U && poly.a_search_center < 4U &&
            poly.a_search_variance >= 8U);
    siqs_poly_clear(&ctx, &poly);
  }
  mpz_clear(ctx.kn);
  case_ctx = NULL;
  puts("PASS sieve: upward-first FB search, exact A windows, lazy cache and reset");
}

static void check_real_case(const mpz_t n, unsigned long forced_k,
                            uint32_t forced_d, uint32_t half) {
  siqs_ctx_t ctx;
  siqs_poly_t poly;
  siqs_factor_array_t result;
  uint32_t family, polynomial, total = 0;
  uint32_t cached_stage = 0U, cached_center = 0U, cached_variance = 0U;
  case_name = "real-polynomials";
  check_trial_survivor(n);
  siqs_factor_array_init(&result, n);
  siqs_ctx_init(&ctx, n, n, &result, NULL, 0);
  case_ctx = &ctx;
  if (forced_k != 0) {
    ctx.multiplier = forced_k;
    mpz_mul_ui(ctx.kn, n, forced_k);
    ctx.params.poly_d = mpz_fdiv_ui(ctx.kn, 8) == 1 ? 2U : 1U;
  }
  if (forced_d != 0) {
    CHECK(forced_d == 1 || (forced_d == 2 && mpz_fdiv_ui(ctx.kn, 8) == 1));
    ctx.params.poly_d = forced_d;
  }
  if (half != 0) ctx.params.half_interval = half;
  CHECK(siqs_ctx_allocate(&ctx));
  siqs_poly_init(&ctx, &poly);
  CHECK(!poly.a_search_ready && !poly.a_search_stage &&
        !poly.a_search_center && !poly.a_search_variance);
  d_seen |= 1U << ctx.params.poly_d;
  for (family = 0; family < 2; family++) {
    uint32_t i, marked = 0U;
    int found = siqs_new_family(&ctx, &poly);
    CHECK(poly.a_search_ready);
    if (family == 0U) {
      cached_stage = poly.a_search_stage;
      cached_center = poly.a_search_center;
      cached_variance = poly.a_search_variance;
    }
    CHECK(poly.a_search_stage == cached_stage &&
          poly.a_search_center == cached_center &&
          poly.a_search_variance == cached_variance);
    for (i = 1U; i < ctx.params.fb_size; i++) marked += ctx.fb[i].in_a != 0;
    CHECK(marked == (found ? poly.q_count : 0U));
    if (!found) break;
    for (polynomial = 0; polynomial < (extended ? 8U : 3U); polynomial++) {
      case_ctx = &ctx;
      check_polynomial_roots(&ctx, &poly);
      compare_sieve(&ctx);
      total++;
      if (!siqs_next_B(&ctx, &poly)) break;
    }
  }
  CHECK(total != 0);
  if (detailed) {
    printf("  CHECK bits=%u k=%lu d=%u q=%u FB=%u M=%u block=%u polys=%u\n",
           ctx.params.bits, ctx.multiplier, ctx.params.poly_d,
           ctx.params.q_count, ctx.params.fb_size, ctx.params.half_interval,
           ctx.block_length, total);
    fflush(stdout);
  }
  siqs_poly_clear(&ctx, &poly);
  siqs_ctx_clear(&ctx);
  {
    uint32_t count;
    mpz_t *values = siqs_factor_array_release(&result, &count);
    gmp_siqs_free(values, count);
  }
  case_ctx = NULL;
}

static void check_real_polynomials(void) {
  static const uint32_t quick_bits[] = {49, 81, 130, 193, 246, 311, 330};
  static const uint32_t more_bits[] = {33, 36, 37, 43, 65, 96, 104, 114,
    144, 145, 167, 184, 185, 218, 219, 237, 245, 269, 270, 299, 300, 310, 364, 370};
  gmp_randstate_t random;
  mpz_t n;
  uint32_t i;
  case_name = "real-polynomials";
  prime_iterator_global_startup();
  gmp_randinit_default(random); gmp_randseed_ui(random, 20261003U);
  mpz_init_set_str(n, "100160063", 10);
  check_real_case(n, 0, 0, 0); /* Known q=1 actual-sieve fixture. */
  for (i = 0; i < sizeof(quick_bits) / sizeof(*quick_bits); i++) {
    make_semiprime(n, random, quick_bits[i]);
    check_real_case(n, 0, 0, 0);
  }
  /* Cover d=1 and d=2 explicitly, plus one-root primes dividing k. */
  do { make_semiprime(n, random, 130); } while (mpz_fdiv_ui(n, 8) != 3);
  check_real_case(n, 3, 1, 0);
  check_real_case(n, 3, 2, 0);
  CHECK(d_seen == ((1U << 1) | (1U << 2)) && special_a && special_k);
  if (extended) {
    for (i = 0; i < sizeof(more_bits) / sizeof(*more_bits); i++) {
      make_semiprime(n, random, more_bits[i]);
      check_real_case(n, 0, 0, 0);
    }
    make_semiprime(n, random, 240);
    check_real_case(n, 0, 0, 524544U); /* 1 MiB plus a partial last block. */
  }
  mpz_clear(n); gmp_randclear(random); prime_iterator_global_shutdown();
  puts("PASS sieve: genuine q=1/q=2 and larger polynomials, d=1/d=2, special roots");
  fflush(stdout);
}

static void suite_sieve(void) {
  check_tier_boundaries();
  check_two_hit_local();
  check_one_hit_and_maps();
  check_oracle_overflow();
  check_A_search_cache();
  check_real_polynomials();
  printf("PASS sieve: %u comparisons, maximum observed physical score %u/255\n",
         comparisons, maximum_score);
}

static void check_policy_curves(void) {
  static const siqs_policy_curve_t line = {60.0, 1.25, 167, UINT16_MAX};
  static const siqs_policy_curve_t ramp = {60.0, 80.0, 167, 177};
  static const siqs_policy_curve_t descending = {120.0, 60.0, 167, 177};
  static const siqs_policy_curve_t exact = {1.0e16, 0.125, 167, 177};
  static const siqs_policy_curve_t single = {7.25, 7.25, 167, 167};
  static const siqs_policy_curve_t derived = {0.0, 0.0, UINT16_MAX, UINT16_MAX};
  case_name = "linear-and-endpoint-curves";
  CHECK(siqs_policy_curve_value(&line, 160) == 51.25);
  CHECK(siqs_policy_curve_value(&line, 180) == 76.25);
  CHECK(siqs_policy_curve_value(&ramp, 167) == ramp.begin);
  CHECK(siqs_policy_curve_value(&ramp, 172) == 70.0);
  CHECK(siqs_policy_curve_value(&ramp, 177) == ramp.step_or_end);
  CHECK(siqs_policy_curve_value(&descending, 172) == 90.0);
  CHECK(siqs_policy_curve_value(&exact, 167) == exact.begin);
  CHECK(siqs_policy_curve_value(&exact, 177) == exact.step_or_end);
  CHECK(siqs_policy_curve_value(&single, 167) == single.begin);
  CHECK(siqs_policy_derived_r(&derived) && !siqs_policy_derived_r(&line));
  puts("PASS policies: linear extrapolation, interpolation, exact endpoints, one-bit ramps and R sentinel");
}

static void check_policy_bounds(void) {
  siqs_ctx_t ctx;
  uint64_t expected;
  case_name = "factor-base-relative-bounds";
  memset(&ctx, 0, sizeof(ctx));
  ctx.largest_fb_prime = 1009;
  ctx.params.lp_multiplier = 1.25;
  ctx.params.sieve_hit_bound_nominal = 10.0;
  ctx.params.sieve_hit_residual_multiplier = 2.0;
  ctx.params.max_large_primes = 1;
  siqs_set_large_prime_bounds(&ctx);
  CHECK(ctx.params.large_prime_bound == 1261 && ctx.params.smooth_bound == 1261);
  CHECK(ctx.params.sieve_hit_bound == 2522.0);
  ctx.params.max_large_primes = 2;
  siqs_set_large_prime_bounds(&ctx);
  expected = UINT64_C(1261) * UINT64_C(1009);
  CHECK(ctx.params.large_prime_bound == 1261 && ctx.params.smooth_bound == expected);
  CHECK(ctx.params.sieve_hit_bound == 2.0 * (double)expected);
  ctx.params.residual_multiplier = 2.5;
  siqs_set_large_prime_bounds(&ctx);
  CHECK(ctx.params.large_prime_bound == 1261 && ctx.params.smooth_bound == 2545202);
  ctx.params.max_large_primes = 1;
  siqs_set_large_prime_bounds(&ctx);
  CHECK(ctx.params.smooth_bound == 2522);
  ctx.params.sieve_hit_bound_nominal = 1.0e8;
  siqs_set_large_prime_bounds(&ctx);
  CHECK(ctx.params.sieve_hit_bound == 1.0e8);
  ctx.params.lp_multiplier = DBL_MAX;
  ctx.params.residual_multiplier = 0.0;
  siqs_set_large_prime_bounds(&ctx);
  CHECK(ctx.params.large_prime_bound == SIQS_LP_MAX && ctx.params.smooth_bound == SIQS_LP_MAX);
  ctx.params.max_large_primes = 2;
  siqs_set_large_prime_bounds(&ctx);
  CHECK(ctx.params.smooth_bound == SIQS_LP_MAX * UINT64_C(1009));
  ctx.largest_fb_prime = UINT32_MAX;
  ctx.params.residual_multiplier = DBL_MAX;
  siqs_set_large_prime_bounds(&ctx);
  CHECK(ctx.params.smooth_bound == UINT64_MAX);
  ctx.params.residual_multiplier = 0.0;
  siqs_set_large_prime_bounds(&ctx);
  CHECK(ctx.params.smooth_bound == UINT64_MAX);
  CHECK(siqs_bound_product(UINT64_MAX, 2, UINT64_MAX) == UINT64_MAX);
  CHECK(siqs_bound_product(UINT64_MAX, 0, UINT64_MAX) == 0);
  CHECK(siqs_scaled_bound(1009, 1.25, UINT64_MAX) == 1261);
  CHECK(siqs_scaled_bound(UINT64_C(9007199254740993), 1.0, UINT64_MAX) ==
        UINT64_C(9007199254740993));
  CHECK(siqs_scaled_bound(UINT64_MAX, 1.0, UINT64_MAX) == UINT64_MAX);
  CHECK(siqs_scaled_bound(UINT64_MAX, 1.5, UINT64_MAX) == UINT64_MAX);
  CHECK(siqs_scaled_bound(1, (double)UINT64_MAX, UINT64_MAX) == UINT64_MAX);
  CHECK(siqs_scaled_bound(1, DBL_MAX, UINT64_MAX) == UINT64_MAX);
  CHECK(siqs_scaled_bound(UINT64_C(4294967297), 4294967297.0, UINT64_MAX) == UINT64_MAX);
  CHECK(siqs_scaled_bound(0, DBL_MAX, UINT64_MAX) == 0);
  puts("PASS policies: fractional/integral K, capped LP, coupled/independent R, final sieve bounds and saturation");
}

static void check_policy_profiles(void) {
  static const struct { uint32_t last; double k; } original[] = {
    {95, 1.0}, {103, 2.0}, {144, 4.0}, {166, 8.0}, {177, 16.0},
    {192, 20.0}, {218, 32.0}, {236, 48.0}, {245, 72.0}
  };
  siqs_policy_t policy;
  siqs_ctx_t ctx;
  uint32_t bits, i, index = 0;
  case_name = "primary-and-recovery-K-values";
  memset(&ctx, 0, sizeof(ctx));
  ctx.largest_fb_prime = 1009;
  for (bits = MPU_SIQS_MIN_BITS; bits <= MPU_SIQS_MAX_BITS; bits++) {
    siqs_resolve_policy(&policy, bits, NULL);
    CHECK(policy.lp_multiplier >= 1.0);
    if (bits >= 246 && bits <= 269) {
      CHECK(policy.max_large_primes == 2 && policy.lp_multiplier == 112.0);
      CHECK(policy.residual_multiplier == 168.0);
      ctx.params.lp_multiplier = policy.lp_multiplier;
      ctx.params.residual_multiplier = policy.residual_multiplier;
      ctx.params.max_large_primes = policy.max_large_primes;
      siqs_set_large_prime_bounds(&ctx);
      CHECK(ctx.params.large_prime_bound == UINT64_C(112) * 1009U);
      CHECK(ctx.params.smooth_bound == UINT64_C(168) * 1009U * 1009U);
      continue;
    }
    if (bits >= 270 && bits <= 299) {
      double t = (double)(bits - 270U) / 29.0;
      double k_l = 112.0 + 54.0 * t;
      double k_r = 168.0 + 82.0 * t;
      CHECK(policy.max_large_primes == 2);
      CHECK(fabs(policy.lp_multiplier - k_l) < 1e-12);
      CHECK(fabs(policy.residual_multiplier - k_r) < 1e-12);
      ctx.params.lp_multiplier = policy.lp_multiplier;
      ctx.params.residual_multiplier = policy.residual_multiplier;
      ctx.params.max_large_primes = policy.max_large_primes;
      siqs_set_large_prime_bounds(&ctx);
      CHECK(ctx.params.large_prime_bound == (uint64_t)(k_l * 1009.0));
      CHECK(ctx.params.smooth_bound == (uint64_t)(k_r * (1009.0 * 1009.0)));
      continue;
    }
    if (bits > 245) {
      CHECK(policy.max_large_primes == 2);
      ctx.params.lp_multiplier = policy.lp_multiplier;
      ctx.params.residual_multiplier = policy.residual_multiplier;
      ctx.params.max_large_primes = policy.max_large_primes;
      siqs_set_large_prime_bounds(&ctx);
      CHECK(ctx.params.large_prime_bound == (uint64_t)(policy.lp_multiplier * 1009.0));
      if (policy.residual_multiplier == 0.0)
        CHECK(ctx.params.smooth_bound == ctx.params.large_prime_bound * 1009U);
      else
        CHECK(ctx.params.smooth_bound ==
              (uint64_t)(policy.residual_multiplier * (1009.0 * 1009.0)));
      continue;
    }
    CHECK(policy.residual_multiplier == 0.0);
    while (bits > original[index].last) index++;
    CHECK(policy.max_large_primes == 1 && policy.lp_multiplier == original[index].k);
    ctx.params.lp_multiplier = policy.lp_multiplier;
    ctx.params.residual_multiplier = policy.residual_multiplier;
    ctx.params.max_large_primes = policy.max_large_primes;
    siqs_set_large_prime_bounds(&ctx);
    CHECK(ctx.params.large_prime_bound == (uint64_t)original[index].k * 1009U);
    CHECK(ctx.params.smooth_bound == ctx.params.large_prime_bound);
  }
  for (i = 0; i < SIQS_RECOVERY_POLICY_COUNT; i++) {
    const siqs_policy_band_t *band = &siqs_recovery_policies[i];
    for (bits = band->first_bits; bits <= band->last_bits; bits++) {
      siqs_resolve_policy(&policy, bits, band);
      CHECK(policy.max_large_primes == 1 && policy.residual_multiplier == 0.0);
      CHECK(policy.lp_multiplier == (i < 3 ? 1.0 : 60.0));
    }
  }
  puts("PASS policies: every primary bit and recovery profile, unchanged smooth/1LP K values");
}

static void check_policy_errors(void) {
#ifndef _WIN32
  unsigned int fixture;
  case_name = "invalid-policy-errors";
  fflush(NULL);
  for (fixture = 0; fixture < 11; fixture++) {
    int status;
    pid_t child = fork();
    CHECK(child >= 0);
    if (child == 0) {
      siqs_policy_curve_t curve = {1.0, 1.0, 300, 310};
      siqs_policy_band_t band = *siqs_policy_band(300);
      siqs_policy_t policy;
      uint32_t bits = 300;
      if (freopen("/dev/null", "w", stderr) == NULL) _exit(1);
      switch (fixture) {
        case 0: bits = 299; break; /* even a constant RAMP cannot extrapolate */
        case 1: bits = 311; break;
        case 2: curve.first_bits = 311; break;
        case 3: curve.last_bits = 300; curve.step_or_end = 2.0; break;
        case 4: curve.first_bits = curve.last_bits = UINT16_MAX; break;
        case 5: curve.first_bits = UINT16_MAX; break;
        case 6: curve.begin = HUGE_VAL; break;
        case 7: band.k_l.begin = band.k_l.step_or_end = 0.0; break;
        case 8: band.k_r = curve; band.k_r.begin = band.k_r.step_or_end = 0.0; break;
        case 9: band.k_r = curve; band.k_r.begin = band.k_r.step_or_end = -1.0; break;
        case 10: curve.begin = DBL_MAX; curve.step_or_end = DBL_MAX;
                 curve.first_bits = 1; curve.last_bits = UINT16_MAX; break;
      }
      if (fixture >= 7 && fixture <= 9) siqs_resolve_policy(&policy, 300, &band);
      else (void)siqs_policy_curve_value(&curve, bits);
      _exit(0);
    }
    CHECK(waitpid(child, &status, 0) == child);
    CHECK(WIFEXITED(status) && WEXITSTATUS(status) == 3);
  }
  puts("PASS policies: invalid ranges, extrapolation, numeric sentinel misuse, nonfinite values and invalid K rejected");
#else
  puts("SKIP policies: fatal-diagnostic checks require POSIX fork");
#endif
}

static void suite_policies(void) {
  check_policy_curves();
  check_policy_bounds();
  check_policy_profiles();
  check_policy_errors();
}

#include "siqs-check-relations.inc.c"
#include "siqs-check-matrix.inc.c"
#include "siqs-check-workers.inc.c"
#include "siqs-check-cofactors.inc.c"

typedef struct {
  const char *name;
  const char *description;
  void (*run)(void);
} check_suite_t;

static const check_suite_t suites[] = {
  {"policies", "K curves, exact endpoints, LP/R bounds, saturation and recovery policies",
   suite_policies},
  {"sieve", "wide-score oracle, byte kernels, blocking, roots, candidate maps",
   suite_sieve},
  {"relations", "congruences, cycle paths, ownership, partition and exponent limits",
   suite_relations},
  {"matrix", "original-column oracles, reduction, packing and nullspace solvers",
   suite_matrix},
  {"workers", "pthread pool lifecycle, reuse, shutdown and injected failures",
   suite_workers},
  {"cofactors", "split counters, acceptance, bounded nested SIQS and quiet/reentrant calls",
   suite_cofactors}
};

static void usage(void) {
  puts("usage: siqs-check [--suite all|policies|sieve|relations|matrix|workers|cofactors] [--extended] [--verbose] [--list]");
}

int main(int argc, char **argv) {
  const char *selected = "all";
  uint32_t i, ran = 0;
  int argument;
  for (argument = 1; argument < argc; argument++) {
    if (strcmp(argv[argument], "--suite") == 0 && argument + 1 < argc)
      selected = argv[++argument];
    else if (strcmp(argv[argument], "--extended") == 0) extended = 1;
    else if (strcmp(argv[argument], "--verbose") == 0) detailed = 1;
    else if (strcmp(argv[argument], "--help") == 0) { usage(); return 0; }
    else if (strcmp(argv[argument], "--list") == 0) {
      for (i = 0; i < sizeof(suites) / sizeof(*suites); i++)
        printf("%s: %s\n", suites[i].name, suites[i].description);
      return 0;
    } else { usage(); return 2; }
  }
  verbose_level = 0;
  printf("SIQS checks: block maximum %u bytes, %s suite\n",
         (unsigned)SIQS_SIEVE_BLOCK_SIZE, extended ? "extended" : "quick");
  fflush(stdout);
  for (i = 0; i < sizeof(suites) / sizeof(*suites); i++)
    if (strcmp(selected, "all") == 0 || strcmp(selected, suites[i].name) == 0) {
      suite_name = suites[i].name;
      suites[i].run();
      ran++;
    }
  if (ran == 0) { fprintf(stderr, "unknown suite: %s\n", selected); return 2; }
  return 0;
}
