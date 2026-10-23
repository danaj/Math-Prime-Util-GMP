/*============================================================================
 * siqs_choose_inc.c -- included by siqs.c after its private types and
 * allocation, RNG and hash-set helpers have been defined.
 *
 * Entry points used by siqs.c:
 *   siqs_choose_A(ctx, poly)  -- select a family; zero on bounded failure.
 *   siqs_choose_A_clear(poly) -- release the chooser's private state.
 *
 * The caller initializes q_count and target_A before the first selection.
 * target_A already includes the multiplier, interval and poly_d divisor.
 * Chooser preparation is deferred until first use, so worker polynomial
 * generators do not allocate selection state. Geometry and arithmetic mode
 * remain fixed for that generator. The main file owns polynomial B/root
 * construction, collection, recovery, and the scorer/certifier policy.
 *
 * Successful selection installs sorted distinct eligible FB indices, A,
 * DA=d*A and fb_in_a marks, and records the family fingerprint. The caller's
 * existing recovery policies handle terminal failure. No family is put
 * into history merely because it was considered during a search.
 *
 * This is more complicated than it needs to be, as it became as much a
 * research project as it is a minimal practical implementation.
 *
 * Copyright (c) 2026 Dana Jacobsen. See siqs.c for authorship and credits.
 *==========================================================================*/

/* See siqs_choose_A() at the end for what is used.  See comments closer to
 * methods for more technical descriptions.
 *
 * Methods available here:
 *
 *  deterministic, works for q=1, q=2, q=3
 *
 *      Custom written to walk the FB in a good order.  Great A values,
 *      but q=3 can be slower than we'd like at that size.
 *
 *  SAR, works for q >= 2
 *
 *      Staged Adaptive Random selects q-1 distinct eligible primes from a
 *      root-centered window, then fits the final prime from the full eligible
 *      FB. It keeps the best unused family while progressively relaxing the
 *      accepted A error. Repeated misses modestly widen the sampling window.
 *
 *  SAC, works for q >= 5
 *
 *      Staged Adaptive Cubic shares SAR's final-prime fitting and acceptance
 *      logic. A polynomial of degree at most three modulo a prime selects
 *      one prime from each of q-1 disjoint pools. Generated prefixes are
 *      unique throughout the cycle; any two share at most three primes.
 *      Each pool prime occurs equally often over the full generated cycle.
 *      Final-prime fitting can add overlap; SAR fallback has no such bound.
 *      The motivation was to reduce redundant relations.
 *
 *  References:
 *     Carrier and Wagstaff, 2002.
 *     Kechlibar, 2005.
 *     Milan, 2010 (TIFA uses random selection prefix).
 */

/*----------------------------------------------------------------------------
 * Common ordering and family history
 *----------------------------------------------------------------------------*/

/* Largest A/target or target/A allowed by the active q>=2 policies.
 * Keep the chooser and score proofs on this same contract, including SAR
 * guard widening and deterministic-pair recovery. q1 has no relative
 * bound; its score proof uses the actual family's polynomial instead. */
static uint32_t siqs_A_maximum_factor(const siqs_ctx_t *ctx, uint32_t q) {
  if (q == 2U) return ctx->params.bits <= 49U ? 32U : 8U;
  if (q == 3U) return 16U;
  return q >= 4U && q <= SIQS_MAX_A_PRIMES ? 2U : 0U;
}

static void siqs_sort_A_indices(uint32_t *indices, uint32_t count) {
  uint32_t i;
  for (i = 1U; i < count; i++) {
    uint32_t key = indices[i], j = i;
    while (j != 0U && indices[j - 1U] > key) {
      indices[j] = indices[j - 1U];
      j--;
    }
    indices[j] = key;
  }
}

static uint64_t siqs_a_fingerprint(const uint32_t *indices, uint32_t count) {
  uint32_t i;
  uint64_t h = UINT64_C(0x6a09e667f3bcc909);
  for (i = 0; i < count; i++)
    h = siqs_mix64(h ^ ((uint64_t)indices[i] +
                        UINT64_C(0x9e3779b97f4a7c15)));
  return h;
}


/*----------------------------------------------------------------------------
 * q=1: nearest eligible prime, with one backup
 *----------------------------------------------------------------------------*/
/* At q=1 the entire product is one prime. Select by absolute distance from
 * target_A, with the smaller prime winning ties. Keep one backup, then let
 * recovery change the geometry. No A tolerance is imposed for this method.
 * With BIDS and the corrected retry policy this small deterministic chooser
 * was competitive with the old random chooser, without its search loops. */
static int siqs_choose_A_q1(siqs_ctx_t *ctx, siqs_poly_t *poly) {
  uint32_t size = ctx->params.fb_size, index;
  if (poly->a_index[0] < size)
    ctx->fb_in_a[poly->a_index[0]] = 0;
  poly->a_index[0] = SIQS_NO_INDEX;
  if (!poly->q1_ready) {
    uint32_t lo = 1, hi = size;
    while (lo < hi) {
      uint32_t mid = lo + (hi - lo) / 2;
      if (mpz_cmp_ui(poly->target_A, ctx->fb[mid].p) > 0)
        lo = mid + 1;
      else
        hi = mid;
    }
    poly->q1_lower = lo - 1;
    poly->q1_upper = lo;
    poly->q1_ready = 1;
  }
  for (;;) {
    int take_lower;
    while (poly->q1_lower && ctx->fb[poly->q1_lower].sqrt_kn == 0)
      poly->q1_lower--;
    while (poly->q1_upper < size && ctx->fb[poly->q1_upper].sqrt_kn == 0)
      poly->q1_upper++;
    if (!poly->q1_lower && poly->q1_upper == size)
      return 0;
    if (!poly->q1_lower)
      take_lower = 0;
    else if (poly->q1_upper == size)
      take_lower = 1;
    else {
      /* Both sides exist, so target_A fits uint32.  All sums fit uint64. */
      uint64_t target = mpz_get_ui(poly->target_A);
      uint64_t lower = ctx->fb[poly->q1_lower].p;
      uint64_t upper = ctx->fb[poly->q1_upper].p;
      take_lower = 2 * target <= lower + upper;
    }
    index = take_lower ? poly->q1_lower-- : poly->q1_upper++;
    if (++poly->q1_count >= 2U) {
      poly->q1_lower = 0;
      poly->q1_upper = size;
    }
    if (!siqs_hashset_insert(&ctx->a_hashes, siqs_a_fingerprint(&index, 1)))
      continue;
    poly->a_index[0] = index;
    ctx->fb_in_a[index] = 1;
    mpz_set_ui(poly->A, ctx->fb[index].p);
    mpz_mul_ui(poly->DA, poly->A, ctx->params.poly_d);
    return 1;
  }
}


/*----------------------------------------------------------------------------
 * q=2: deterministic pairs and disjoint-prime passes
 *----------------------------------------------------------------------------*/
/* Most q=2 inputs finish in the first few families, so preparation and the
 * first choices matter. We walk canonical prefixes and fit the final prime
 * by absolute distance. Within a pass, primes already emitted are skipped;
 * another pass may reuse primes, while family history prevents exact reuse.
 *
 * Start with narrow percentage bands and widen only after an empty pass.
 * The primary q=2 range uses final factors 2,4,8,16,32; the factor-32 tail
 * avoids rare cases which otherwise immediately need q=1 recovery. Existing
 * higher-bit q=2 recovery retains factor 8. These are collection policies,
 * not proofs of scorer safety.
 *
 * q=2 supplies a constant row size of two.
 * The row engine also supports q=3. */
typedef struct {
  uint8_t *used;
  uint32_t i, j, emitted, stage;
  uint64_t pass;             /* For review/observability; not an epoch stamp. */
  mpz_t lower, upper, prefix, scratch;
  uint64_t native[3];        /* target, lower, upper when all fit uint64 */
  int use_native;
  uint32_t stage_count;
} siqs_A_pairs_t;

static void siqs_A_pairs_band(siqs_A_pairs_t *s, const siqs_poly_t *poly) {
  static const unsigned permille[10] = {2, 4, 8, 16, 32, 50, 100, 200, 400, 500};
  if (s->stage < 10) {
    mpz_mul_ui(s->scratch, poly->target_A, permille[s->stage]);
    mpz_fdiv_q_ui(s->scratch, s->scratch, 1000);
    mpz_sub(s->lower, poly->target_A, s->scratch);
    mpz_add(s->upper, poly->target_A, s->scratch);
  } else {
    unsigned factor = 1U << (s->stage - 9U);
    mpz_cdiv_q_ui(s->lower, poly->target_A, factor);
    mpz_mul_ui(s->upper, poly->target_A, factor);
  }
  s->use_native = mpz_sizeinbase(s->upper, 2) <= 64;
  if (s->use_native) {
    memset(s->native, 0, sizeof(s->native));
    mpz_export(s->native, NULL, -1, sizeof(uint64_t), 0, 0, poly->target_A);
    mpz_export(s->native+1, NULL, -1, sizeof(uint64_t), 0, 0, s->lower);
    mpz_export(s->native+2, NULL, -1, sizeof(uint64_t), 0, 0, s->upper);
  }
}

/* q must be 2 or 3; FB index zero is prime 2. maximum_factor is a
 * power of two in [4,2^30].
 * The FB and target remain fixed for this state's lifetime. */
static void siqs_A_pairs_init(siqs_A_pairs_t *s, const siqs_ctx_t *ctx,
                        const siqs_poly_t *poly, uint32_t maximum_factor) {
  uint32_t factor;
  if (maximum_factor < 4U || maximum_factor > (1U << 30) ||
      (maximum_factor & (maximum_factor - 1U)) != 0U)
    croak("SIQS: invalid small-q maximum factor");
  memset(s, 0, sizeof(*s));
  s->stage_count = 12U;
  for (factor = 4U; factor < maximum_factor; factor <<= 1)
    s->stage_count++;
  s->used = siqs_calloc(ctx->params.fb_size, sizeof(*s->used));
  mpz_inits(s->lower, s->upper, s->prefix, s->scratch, NULL);
  s->pass = 1;
  s->i = 1;
  s->j = 2;
  siqs_A_pairs_band(s, poly);
}

static void siqs_A_pairs_clear(siqs_A_pairs_t *s) {
  free(s->used);
  mpz_clears(s->lower, s->upper, s->prefix, s->scratch, NULL);
}

/* Values above UINT32_MAX need only the sentinel UINT32_MAX+1. */
static uint64_t siqs_A_pairs_project(siqs_A_pairs_t *s, const mpz_t value,
                              unsigned which, uint64_t prefix, int ceiling) {
  if (s->use_native)
    return s->native[which] / prefix +
           (ceiling && s->native[which] % prefix != 0);
  if (ceiling) mpz_cdiv_q(s->scratch, value, s->prefix);
  else         mpz_fdiv_q(s->scratch, value, s->prefix);
  return mpz_cmp_ui(s->scratch, UINT32_MAX) > 0 ? UINT64_C(1) << 32
                                              : mpz_get_ui(s->scratch);
}

/* Returns 0 after every fingerprint-available family in the final band has
 * been emitted.  This is not exhaustion outside the configured guard.
 */
static int siqs_A_pairs_next(siqs_A_pairs_t *s, siqs_ctx_t *ctx,
                              siqs_poly_t *poly, uint32_t q) {
  uint32_t size = ctx->params.fb_size, at;
  for (at = 0; at < q; at++) {
    if (poly->a_index[at] < size) ctx->fb_in_a[poly->a_index[at]] = 0;
    poly->a_index[at] = SIQS_NO_INDEX;
  }
  while (s->stage < s->stage_count) {
    uint32_t index[3], last, low, high, end;
    uint64_t prefix, wanted, minimum, maximum;
    if ((q == 2 ? s->i : s->j) + 1 >= size) {
      if (s->emitted == 0) {
        if (++s->stage == s->stage_count) return 0;
        siqs_A_pairs_band(s, poly);
      }
      memset(s->used, 0, size * sizeof(*s->used));
      s->emitted = 0;
      s->pass++;
      s->i = 1;
      s->j = 2;
      if (size <= q) continue;
    }
    /* q2 rows ascend; q3 prefixes visit (1,2), (1,3), (2,3), ... */
    index[0] = s->i;
    index[1] = s->j;
    last = index[q-2];
    if (q == 2) s->i++;
    else if (++s->i == s->j) { s->i = 1; s->j++; }
    for (at = 0; at + 1 < q; at++)
      if (s->used[index[at]] || ctx->fb[index[at]].sqrt_kn == 0) break;
    if (at + 1 < q) continue;
    prefix = ctx->fb[index[0]].p;
    if (q == 3) prefix *= ctx->fb[index[1]].p;
    if (!s->use_native) siqs_mpz_set_u64(s->prefix, prefix);
    minimum = siqs_A_pairs_project(s, s->lower, 1, prefix, 1);
    maximum = siqs_A_pairs_project(s, s->upper, 2, prefix, 0);
    if (minimum > maximum || maximum <= ctx->fb[last].p ||
        minimum > ctx->fb[size-1].p) continue;
    wanted = siqs_A_pairs_project(s, poly->target_A, 0, prefix, 1);
    high = last + 1;
    end = size;
    while (high < end) {
      uint32_t mid = high + (end - high) / 2;
      if (ctx->fb[mid].p < wanted) high = mid + 1;
      else end = mid;
    }
    low = high - 1;
    for (;;) {
      uint32_t pick;
      while (low > last &&
             (s->used[low] || ctx->fb[low].sqrt_kn == 0)) low--;
      while (high < size &&
             (s->used[high] || ctx->fb[high].sqrt_kn == 0)) high++;
      if (low > last && ctx->fb[low].p < minimum) low = last;
      if (high < size && ctx->fb[high].p > maximum) high = size;
      if (low == last && high == size) break;
      /* Odd closing primes have integer midpoints.  Thus ceil(T/prefix)
       * preserves exact absolute-distance comparisons, including ties. */
      pick = high == size || (low > last &&
             2 * wanted <= (uint64_t)ctx->fb[low].p + ctx->fb[high].p)
           ? low-- : high++;
      index[q-1] = pick;
      if (!siqs_hashset_insert(&ctx->a_hashes, siqs_a_fingerprint(index, q)))
        continue;
      for (at = 0; at < q; at++) {
        poly->a_index[at] = index[at];
        s->used[index[at]] = ctx->fb_in_a[index[at]] = 1;
      }
      siqs_mpz_set_u64(poly->A, prefix);
      mpz_mul_ui(poly->A, poly->A, ctx->fb[pick].p);
      mpz_mul_ui(poly->DA, poly->A, ctx->params.poly_d);
      s->emitted++;
      return 1;
    }
  }
  return 0;
}
static int siqs_choose_A_pairs(siqs_ctx_t *ctx, siqs_poly_t *poly,
                                uint32_t q, uint32_t maximum_factor) {
  if (poly->q2_search == NULL) {
    siqs_A_pairs_t *s = siqs_malloc(sizeof(*s));
    siqs_A_pairs_init(s, ctx, poly, maximum_factor);
    poly->q2_search = s;
  }
  return siqs_A_pairs_next(poly->q2_search, ctx, poly, q);
}


static int siqs_choose_A_deterministic(siqs_ctx_t *ctx, siqs_poly_t *poly) {
  if (poly->q_count == 1)
    return siqs_choose_A_q1(ctx, poly);

  if (poly->q_count == 2)
    return siqs_choose_A_pairs(ctx, poly, 2U,
                               siqs_A_maximum_factor(ctx, 2U));

  if (poly->q_count == 3)
    return siqs_choose_A_pairs(ctx, poly, 3U, 8U);

  croak("SIQS: Only q = {1, 2, 3} have deterministic A choosers");
  return 0;
}


/*----------------------------------------------------------------------------
 * SAR and Cubic: shared fitting and adaptive acceptance
 *----------------------------------------------------------------------------*/
/* Staged Adaptive RANDOM (SAR) samples q-1 distinct primes from a prepared
 * root-centered interval and fits the final prime. Cubic supplies structured
 * prefixes to the same fitting and acceptance machinery.
 *
 * Keep the best unused family while relaxing the accepted multiplicative
 * error through the table below. Start with one valid proposal, then
 * try one more per failed error level. A later failed search leaves the
 * retained best intact. Very tight distances gave attractive statistics
 * without paying for themselves in full factoring, so we retain this
 * inexpensive bounded schedule.
 *
 * SAR gets 128 raw proposals per candidate. Every 32 misses widen its local
 * radius by ceil(starting_radius/16), clipped at the eligible FB edges. Each
 * new candidate search starts from the prepared radius. q=3 alone has the
 * tested outer guard stages 2,4,8,16; higher q keeps its factor-two guard.
 *
 * With 64-bit unsigned long and native uint128_t, setup checks that target_A
 * and even the largest possible FB-prime product fit in 64 bits. Native
 * fitting uses uint64_t products and uint128_t cross-products for exact
 * comparisons.
 * Otherwise the same selection and rounding rules use GMP. This proof is
 * independent of the input's bit size and is done once during preparation. */
#ifndef USE_CUBIC_CHOOSER
# define USE_CUBIC_CHOOSER 1
#endif
#ifndef SIQS_A_CUBIC_MIN_Q
# define SIQS_A_CUBIC_MIN_Q 8U
#endif

typedef struct {
  uint32_t first, field, period, state, draws;
} siqs_A_cubic_t;

enum { SIQS_A_RANDOM_RETRIES = 128, SIQS_A_RANDOM_BLOCK = 32,
       SIQS_A_CUBIC_RETRIES = 32 };

/* Accepted multiplicative error, in hundredths of one percent (100 = 1%).
 * Edit this table to change the values or the number of error levels.
 * For error e and target T, accept T/(1+e) <= A <= T*(1+e). */
static const uint32_t siqs_A_error_basis_points[] = {
   50U,  /* 0.5% */
  100U,  /* 1%   */
  200U,  /* 2%   */
  400U,  /* 4%   */
  800U   /* 8%   */
};
enum { SIQS_A_ERROR_SCALE = 10000U,
       SIQS_A_ERROR_LEVELS = sizeof(siqs_A_error_basis_points) /
                            sizeof(siqs_A_error_basis_points[0]) };

struct siqs_A_search {
  const siqs_ctx_t *ctx;
  uint32_t *eligible;
  uint32_t *native_primes;
  uint32_t size, center, radius;
  siqs_A_cubic_t *cubic;
  int random_mode;
  mpz_t target, squared, wanted, candidate, scratch, base, best, selected;
  mpz_t hard_lower, hard_upper;
  mpz_t lower[SIQS_A_ERROR_LEVELS], upper[SIQS_A_ERROR_LEVELS];
  uint32_t selected_indices[SIQS_MAX_A_PRIMES];
#if HAVE_UINT128 && ULONG_MAX > UINT32_MAX
  int native;
  unsigned long native_target;
  uint64_t native_best, native_selected, native_hard_lower, native_hard_upper;
  uint64_t native_lower[SIQS_A_ERROR_LEVELS], native_upper[SIQS_A_ERROR_LEVELS];
  uint128_t native_squared;
#endif
};

/* Both transitions are permanent for this prepared search; report only once. */
static void siqs_A_switch_to_sar(siqs_A_search_t *g, const char *reason) {
  if (g->random_mode) return;
  g->random_mode = 1;
  if (g->ctx->verbose > 0) {
    printf("# siqs A chooser Cubic -> SAR: %s (%u/%u prefix proposals used)\n",
           reason, g->cubic->draws, g->cubic->period);
    fflush(stdout);
  }
}

static void siqs_A_search_clear(siqs_poly_t *poly) {
  siqs_A_search_t *g = poly->A_search;
  uint32_t i;
  if (g == NULL)
    return;
  free(g->eligible);
  free(g->native_primes);
  free(g->cubic);
  for (i = 0; i < SIQS_A_ERROR_LEVELS; i++)
    mpz_clears(g->lower[i], g->upper[i], NULL);
  mpz_clears(g->target, g->squared, g->wanted, g->candidate, g->scratch,
             g->base, g->best, g->selected, g->hard_lower, g->hard_upper, NULL);
  free(g);
  poly->A_search = NULL;
}

static uint32_t siqs_A_prime(const siqs_A_search_t *g, uint32_t position) {
  return g->ctx->fb[g->eligible[position]].p;
}

static int siqs_A_available(uint32_t position, const uint32_t *prefix,
                            uint32_t count) {
  uint32_t i;
  for (i = 0; i < count; i++)
    if (prefix[i] == position)
      return 0;
  return 1;
}

#if HAVE_UINT128 && ULONG_MAX > UINT32_MAX
/* All candidate products fit uint64. Upper bounds may exceed it: clamp
 * those bounds, preserving every comparison with a representable product. */
static uint64_t siqs_A_native_bound(mpz_srcptr bound) {
  return mpz_fits_ulong_p(bound) ? mpz_get_ui(bound) : ULONG_MAX;
}

static void siqs_A_native_guard(siqs_A_search_t *g) {
  if (g->native) {
    g->native_hard_lower = mpz_get_ui(g->hard_lower);
    g->native_hard_upper = siqs_A_native_bound(g->hard_upper);
  }
}

/* Minimize max(A/T,T/A), not absolute distance. When candidates straddle T,
 * their product compared with T^2 selects the closer multiplicative ratio.
 * Strict comparisons preserve the first candidate on an exact tie. */
static int siqs_A_closer64(const siqs_A_search_t *g, uint64_t a, uint64_t b) {
  uint64_t t = g->native_target;
  if (a <= t && b <= t) return a > b;
  if (a >= t && b >= t) return a < b;
  uint128_t product = (uint128_t)a * b;
  return a < t ? product > g->native_squared : product < g->native_squared;
}
#endif

static int siqs_A_closer(siqs_A_search_t *g, mpz_srcptr a, mpz_srcptr b) {
#if HAVE_UINT128 && ULONG_MAX > UINT32_MAX
  if (g->native) {
    return siqs_A_closer64(g, mpz_get_ui(a), mpz_get_ui(b));
  }
#endif
  int as = mpz_cmp(a, g->target), bs = mpz_cmp(b, g->target);
  if (as <= 0 && bs <= 0)
    return mpz_cmp(a, b) > 0;
  if (as >= 0 && bs >= 0)
    return mpz_cmp(a, b) < 0;
  mpz_mul(g->scratch, a, b);
  return as < 0 ? mpz_cmp(g->scratch, g->squared) > 0
                : mpz_cmp(g->scratch, g->squared) < 0;
}

static uint32_t siqs_A_lower_bound(const siqs_A_search_t *g, mpz_srcptr wanted) {
  uint32_t lo = 0, hi = g->size;
  while (lo < hi) {
    uint32_t mid = lo + (hi - lo) / 2;
    if (mpz_cmp_ui(wanted, siqs_A_prime(g, mid)) > 0)
      lo = mid + 1;
    else
      hi = mid;
  }
  return lo;
}

static uint32_t siqs_A_last_prime(siqs_A_search_t *g, const uint32_t *prefix,
                                  mpz_srcptr base) {
  uint32_t lo, side, best = SIQS_NO_INDEX;
  uint32_t count = g->ctx->params.q_count - 1;
  mpz_cdiv_q(g->wanted, g->target, base);
  lo = siqs_A_lower_bound(g, g->wanted);
  for (side = 0; side < 2; side++) {
    int64_t at = side ? (int64_t)lo : (int64_t)lo - 1;
    while (at >= 0 && at < g->size &&
           !siqs_A_available((uint32_t)at, prefix, count))
      at += side ? 1 : -1;
    if (at < 0 || at >= g->size)
      continue;
    mpz_mul_ui(g->candidate, base, siqs_A_prime(g, (uint32_t)at));
    if (best == SIQS_NO_INDEX || siqs_A_closer(g, g->candidate, g->best)) {
      mpz_set(g->best, g->candidate);
      best = (uint32_t)at;
    }
  }
  if (best == SIQS_NO_INDEX)
    croak("SIQS: no eligible final A prime");
  return best;
}

#if HAVE_UINT128 && ULONG_MAX > UINT32_MAX
/* The caller proves that even the largest possible product fits uint64_t.
 * Cross-products use 128 bits, preserving the exact geometric comparison. */
static uint32_t siqs_A_last_prime64(siqs_A_search_t *g, const uint32_t *prefix) {
  uint32_t count = g->ctx->params.q_count-1U, lo = 0, hi = g->size;
  uint32_t best = SIQS_NO_INDEX;
  uint64_t base = 1, best_a = 0;
  for (uint32_t i = 0; i < count; i++) base *= g->native_primes[prefix[i]];
  uint64_t wanted = g->native_target/base + (g->native_target%base != 0);
  while (lo < hi) {
    uint32_t mid = lo + (hi-lo)/2;
    if (g->native_primes[mid] < wanted) lo = mid+1;
    else hi = mid;
  }
  for (unsigned side = 0; side < 2; side++) {
    int64_t at = side ? (int64_t)lo : (int64_t)lo-1;
    while (at >= 0 && at < g->size &&
           !siqs_A_available((uint32_t)at,prefix,count)) at += side ? 1 : -1;
    if (at < 0 || at >= g->size) continue;
    uint64_t a = base * g->native_primes[(uint32_t)at];
    if (best == SIQS_NO_INDEX || siqs_A_closer64(g,a,best_a)) {
      best_a = a; best = (uint32_t)at;
    }
  }
  if (best == SIQS_NO_INDEX) croak("SIQS: no eligible final A prime");
  g->native_best = best_a;
  return best;
}
#endif

#if USE_CUBIC_CHOOSER
/*----------------------------------------------------------------------------
 * Cubic prefixes over disjoint prime pools
 *----------------------------------------------------------------------------
 * Each coordinate owns a disjoint pool of field_size primes. A degree-three
 * polynomial over F_7 or F_19 selects one slot per coordinate; the leading
 * coefficient supplies the extra coordinate when needed. Distinct codewords
 * share at most three prefix primes. Fitting the final prime can add overlap,
 * and the RANDOM fallback has no such bound.
 *
 * Cubic is supported from q=5 through q=16; policy currently starts at q=8.
 * F_7 serves q<=9, F_19 the rest. The complete period is field_size^4, and
 * rejected proposals also consume states. On exhaustion, or an initial
 * bounded miss, use the already prepared RANDOM fallback. No batch is used.
 * TODO: a batch coordinate could extend the stream; overlap bounds between
 * different batches would need a separate analysis.
 *
 * F_11 may deserve another look in selected ranges. Our field experiments
 * did not justify adding it to this policy. Pool-mean balancing could center
 * the remaining prime near a consistent 1x--2x root; q=8 factoring tests did
 * not find an overall gain, so retain the measured window.
 *--------------------------------------------------------------------------*/
static void siqs_A_cubic_init(siqs_A_search_t *g, const siqs_poly_t *poly,
                               double root) {
  uint32_t q = poly->q_count, field = q <= 9U ? 7U : 19U;
  uint32_t window = (q-1U)*field;
  if (q < SIQS_A_CUBIC_MIN_Q || q < 5U || q-1U > field+1U || g->size < window)
    return;
  /* Replacing high slots by each pool's first prime lowers the prefix mean
   * while keeping pools disjoint, but can repeat prefixes and weaken the
   * overlap bound. q6/q7 tests found no sieve-work gain. Use RANDOM instead. */
  if (q == 6U || q == 7U) {
    double minimum = 1.0, desired = mpz_get_d(g->target)/root;
    for (uint32_t j = 0; j < q-1U; j++) {
      double sum = 0.0;
      for (uint32_t k = 0; k < field; k++) sum += siqs_A_prime(g,j*field+k);
      minimum *= sum/field;
    }
    if (minimum > desired) return;
  }
  siqs_A_cubic_t *c = siqs_calloc(1,sizeof(*c));
  c->field = field;
  c->period = field*field*field*field;
  c->state = ((uint64_t)(uint32_t)g->ctx->poly_rng.state+77U)%c->period;
  c->first = g->center > window/2U ? g->center-window/2U : 0;
  if (c->first > g->size-window) c->first = g->size-window;
  g->cubic = c;
  g->random_mode = 0;
  g->radius = window; /* Twice the nominal Cubic radius, ready for RANDOM. */
}

/* For either odd prime p=field, state -> (p+1)*state+1 has a full period
 * modulo p^4. Its four base-p digits therefore visit every coefficient tuple
 * once; the seed changes only the starting point. State persists across
 * families, and even rejected proposals consume a draw. Report exhaustion
 * instead of wrapping and repeating the prefix stream.
 *
 * Coordinates x<p evaluate that cubic; x=p, when present, is its leading
 * coefficient (the projective point at infinity). Two distinct cubics agree
 * at at most three coordinates, including infinity. Disjoint coordinate
 * pools turn that into at most three shared prefix primes. This bound does
 * not cover the fitted final prime or the caller's RANDOM fallback. */
static int siqs_A_cubic_prefix(siqs_A_search_t *g, uint32_t q,
                                uint32_t *prefix) {
  siqs_A_cubic_t *c = g->cubic;
  uint32_t value, coefficients[4];
  if (c->draws == c->period) return 0;
  c->state = ((c->field+1U)*c->state+1U)%c->period;
  c->draws++;
  value = c->state;
  for (unsigned i = 0; i < 4; i++) {
    coefficients[i] = value%c->field; value /= c->field;
  }
  for (uint32_t x = 0; x < q-1U; x++) {
    value = 0;
    if (x == c->field) value = coefficients[3];
    else for (unsigned i = 4; i > 0; i--)
      value = (value*x+coefficients[i-1U])%c->field;
    prefix[x] = c->first+c->field*x+value;
  }
  return 1;
}
#endif

static siqs_A_search_t *siqs_choose_A_prepare(const siqs_ctx_t *ctx,
                                          const siqs_poly_t *poly, int cubic) {
  siqs_A_search_t *g = (siqs_A_search_t *)siqs_calloc(1, sizeof(*g));
  uint32_t i, lo, hi, wanted;
  double root;
  g->ctx = ctx;
  g->random_mode = 1;
  mpz_init_set(g->target, poly->target_A);
  mpz_inits(g->squared, g->wanted, g->candidate, g->scratch,
             g->base, g->best, g->selected, g->hard_lower, g->hard_upper, NULL);
  mpz_mul(g->squared, g->target, g->target);
  mpz_cdiv_q_ui(g->hard_lower, g->target, 2U);
  mpz_mul_ui(g->hard_upper, g->target, 2U);
  /* Cache exact integer bounds once: ceil(T/(1+e)), floor(T*(1+e)). */
  for (i = 0; i < SIQS_A_ERROR_LEVELS; i++) {
    uint32_t numerator = SIQS_A_ERROR_SCALE + siqs_A_error_basis_points[i];
    mpz_inits(g->lower[i], g->upper[i], NULL);
    mpz_mul_ui(g->lower[i], g->target, SIQS_A_ERROR_SCALE);
    mpz_cdiv_q_ui(g->lower[i], g->lower[i], numerator);
    mpz_mul_ui(g->upper[i], g->target, numerator);
    mpz_fdiv_q_ui(g->upper[i], g->upper[i], SIQS_A_ERROR_SCALE);
  }
  g->eligible = (uint32_t *)siqs_malloc(
      siqs_array_size(ctx->params.fb_size, sizeof(*g->eligible)));
  for (i = 1; i < ctx->params.fb_size; i++)
    if (ctx->fb[i].sqrt_kn != 0)
      g->eligible[g->size++] = i;
  if (g->size == 0 || mpz_sgn(g->target) <= 0)
    return g;
  root = pow(mpz_get_d(g->target), 1.0 / poly->q_count);
  wanted = root < 3.0 ? 3U : root >= UINT32_MAX ? UINT32_MAX : (uint32_t)root;
  lo = 0;
  hi = g->size;
  while (lo < hi) {
    uint32_t mid = lo + (hi - lo) / 2U;
    if (siqs_A_prime(g, mid) < wanted)
      lo = mid + 1U;
    else
      hi = mid;
  }
  g->center = lo == g->size ? g->size - 1U : lo;
  if (lo > 0 && lo < g->size &&
      wanted - siqs_A_prime(g, lo-1U) <= siqs_A_prime(g, lo) - wanted)
    g->center = lo - 1U;

  /* SAR draws random prefix primes from a window around the root center.
   * Kechlibar calls for an appropriately long interval in the factor base;
   * "jasonp" gave roughly 50 candidate primes as an example on MathOverflow
   * in 2011.
   *
   * We follow Tilman Neumann's idea in JML of scaling sqrt(FB) to set the
   * radius. JML uses 0.75. Our 51-144-bit SAR tests worked well with
   * coefficients from 0.2 to 1.0 (with the minimum radius below). */
  /* TODO: More benchmarks, 0.5 might be better for us */
  g->radius = (uint32_t)(0.75 * sqrt((double)ctx->params.fb_size));
  if (g->radius < 8U)
    g->radius = 8U;

#if HAVE_UINT128 && ULONG_MAX > UINT32_MAX
  /* Prove that every possible family product fits before using native
   * fitting, shared by SAR and Cubic. */
  unsigned long largest = siqs_A_prime(g, g->size-1U), product = 1;
  g->native = mpz_fits_ulong_p(g->target);
  for (i = 0; g->native && i < poly->q_count; i++) {
    if (product > ULONG_MAX / largest) g->native = 0;
    else product *= largest;
  }
  g->native_target = mpz_get_ui(g->target);
  g->native_squared = (uint128_t)g->native_target * g->native_target;
  if (g->native) {
    /* Only native fitting needs this cache; GMP searches skip it. */
    g->native_primes = siqs_malloc(siqs_array_size(g->size, sizeof(*g->native_primes)));
    for (i = 0; i < g->size; i++) g->native_primes[i] = siqs_A_prime(g, i);
    siqs_A_native_guard(g);
    for (i = 0; i < SIQS_A_ERROR_LEVELS; i++) {
      g->native_lower[i] = mpz_get_ui(g->lower[i]);
      g->native_upper[i] = siqs_A_native_bound(g->upper[i]);
    }
  }
#endif

#if USE_CUBIC_CHOOSER
  if (cubic) siqs_A_cubic_init(g, poly, root);
#endif

  return g;
}

/* Sample distinct positions from [center-radius, center+radius), clipped
 * to [0, limit). The caller supplies all policy and geometry. */
static int siqs_random_positions(siqs_rng_t *rng, uint32_t limit,
                                  uint32_t center, uint32_t radius,
                                  uint32_t count, uint32_t *positions) {
  uint32_t lo, hi, span;
  uint64_t upper;
  if (!count) return 1;
  if (center >= limit) return 0;
  lo = center > radius ? center - radius : 0;
  upper = (uint64_t)center + radius;
  hi = upper < limit ? (uint32_t)upper : limit;
  span = hi - lo;
  if (span < count) return 0;

  /* Floyd sampling: exactly count draws, with no duplicate-retry loop. */
  for (uint32_t i = 0, j = span - count; i < count; i++, j++) {
    uint32_t pick = lo + siqs_rand_range(rng, j + 1U);
    for (uint32_t k = 0; k < i; k++)
      if (positions[k] == pick) { pick = lo + j; break; }
    positions[i] = pick;
  }
  return 1;
}

static int siqs_A_history_contains(const siqs_hashset_t *set, uint64_t h) {
  uint32_t pos;
  h = siqs_mix64(h);
  pos = (uint32_t)h & (set->alloc - 1U);
  while (set->slots[pos] != SIQS_HASH_EMPTY) {
    if (set->slots[pos] == h)
      return 1;
    pos = (pos + 1U) & (set->alloc - 1U);
  }
  return 0;
}

/* Find one unused guarded candidate. Later misses leave the retained best
 * intact. Only the final accepted family is inserted into the history. */
static int siqs_adaptive_A_candidate(siqs_ctx_t *ctx, siqs_poly_t *poly,
                                    int *have_best) {
  siqs_A_search_t *g = poly->A_search;
  uint32_t q = poly->q_count, radius = g->radius;
  uint32_t step = radius / 16U + (radius % 16U != 0U);
  uint32_t maximum = g->center, attempt, i;
  uint32_t retry_limit = g->random_mode ? SIQS_A_RANDOM_RETRIES : SIQS_A_CUBIC_RETRIES;
  uint32_t prefix[SIQS_MAX_A_PRIMES], indices[SIQS_MAX_A_PRIMES];
  if (g->center < g->size && g->size - g->center > maximum)
    maximum = g->size - g->center;
  /* Every 32 misses add ceil(starting radius/16); stop expanding at the FB
   * edges. The next candidate search starts from the prepared radius. */
  for (attempt = 0; attempt < retry_limit; attempt++) {
#if USE_CUBIC_CHOOSER
    if (!g->random_mode && !siqs_A_cubic_prefix(g, q, prefix)) {
      siqs_A_switch_to_sar(g, "prefix period exhausted");
      retry_limit = SIQS_A_RANDOM_RETRIES;
      attempt = 0;
    }
#endif
    if (g->random_mode) {
      if (attempt != 0 && attempt % SIQS_A_RANDOM_BLOCK == 0 && radius < maximum) {
        uint32_t room = maximum - radius;
        radius += step < room ? step : room;
      }
      if (!siqs_random_positions(&ctx->poly_rng, g->size, g->center,
                                   radius, q-1U, prefix))
        return 0;
    }
#if HAVE_UINT128 && ULONG_MAX > UINT32_MAX
    if (g->native) prefix[q-1U] = siqs_A_last_prime64(g, prefix);
    else
#endif
    {
      mpz_set_ui(g->base, 1U);
      for (i = 0; i < q-1U; i++)
        mpz_mul_ui(g->base, g->base, siqs_A_prime(g, prefix[i]));
      prefix[q-1U] = siqs_A_last_prime(g, prefix, g->base);
    }
    int guarded;
#if HAVE_UINT128 && ULONG_MAX > UINT32_MAX
    if (g->native)
      guarded = g->native_best >= g->native_hard_lower &&
                g->native_best <= g->native_hard_upper;
    else
#endif
      guarded = mpz_cmp(g->best, g->hard_lower) >= 0 &&
                mpz_cmp(g->best, g->hard_upper) <= 0;
    if (!guarded) {
      continue;
    }
    for (i = 0; i < q; i++)
      indices[i] = g->eligible[prefix[i]];
    siqs_sort_A_indices(indices, q);
    uint64_t fingerprint = siqs_a_fingerprint(indices, q);
    if (siqs_A_history_contains(&ctx->a_hashes, fingerprint)) {
      continue;
    }
    int improved;
#if HAVE_UINT128 && ULONG_MAX > UINT32_MAX
    if (g->native)
      improved = !*have_best || siqs_A_closer64(g, g->native_best, g->native_selected);
    else
#endif
      improved = !*have_best || siqs_A_closer(g, g->best, g->selected);
    if (improved) {
#if HAVE_UINT128 && ULONG_MAX > UINT32_MAX
      if (g->native) g->native_selected = g->native_best;
      else
#endif
      mpz_set(g->selected, g->best);
      memcpy(g->selected_indices, indices, q * sizeof(*indices));
      *have_best = 1;
    }
    return 1;
  }
  /* An initial bounded Cubic miss gets one bounded RANDOM rescue. Later
   * misses preserve the best candidate already retained by the caller. */
  if (!g->random_mode && !*have_best) {
    siqs_A_switch_to_sar(g, "bounded search found no usable A");
    return siqs_adaptive_A_candidate(ctx, poly, have_best);
  }
  return 0;
}

static int siqs_adaptive_choose_A(siqs_ctx_t *ctx, siqs_poly_t *poly) {
  siqs_A_search_t *g;
  uint32_t i, q = poly->q_count;
  int have_best = 0;
  if (q < 2U || q > SIQS_MAX_A_PRIMES)
    croak("SIQS: unsupported RANDOM geometry");
  for (i = 0; i < q; i++)
    if (poly->a_index[i] < ctx->params.fb_size)
      ctx->fb_in_a[poly->a_index[i]] = 0;
  if (poly->A_search == NULL)
    poly->A_search = siqs_choose_A_prepare(ctx, poly, 0);
  g = poly->A_search;
  if (g->size < q || !siqs_adaptive_A_candidate(ctx, poly, &have_best))
    return 0;
  /* Check the retained best against each precomputed error band. If
   * outside, try one new valid proposal before advancing to the next band.
   * After the last band (or a failed search), return the best we have. */
  for (i = 0; i < SIQS_A_ERROR_LEVELS; i++) {
    int within;
#if HAVE_UINT128 && ULONG_MAX > UINT32_MAX
    if (g->native)
      within = g->native_selected >= g->native_lower[i] &&
               g->native_selected <= g->native_upper[i];
    else
#endif
      within = mpz_cmp(g->selected, g->lower[i]) >= 0 &&
               mpz_cmp(g->selected, g->upper[i]) <= 0;
    if (within) break;
    if (!siqs_adaptive_A_candidate(ctx, poly, &have_best))
      break;
  }
  if (!siqs_hashset_insert(&ctx->a_hashes,
                            siqs_a_fingerprint(g->selected_indices, q)))
    croak("SIQS: accepted an already used A family");
  memcpy(poly->a_index, g->selected_indices, q * sizeof(*poly->a_index));
  for (i = 0; i < q; i++)
    ctx->fb_in_a[poly->a_index[i]] = 1;
#if HAVE_UINT128 && ULONG_MAX > UINT32_MAX
  if (g->native) mpz_set_ui(poly->A, g->native_selected);
  else
#endif
  mpz_set(poly->A, g->selected);
  mpz_mul_ui(poly->DA, poly->A, ctx->params.poly_d);
  return 1;
}



/* Keep the usual factor-two path intact. Relax only after its initial
 * bounded search fails. Restore factor two for the next family. */
static int siqs_choose_A_sar(siqs_ctx_t *ctx, siqs_poly_t *poly) {
  uint32_t maximum_factor = siqs_A_maximum_factor(ctx, poly->q_count);
  if (siqs_adaptive_choose_A(ctx, poly)) return 1;
  siqs_A_search_t *g = poly->A_search;
  int found = 0;
  for (uint32_t factor = 4; factor <= maximum_factor; factor *= 2) {
    mpz_cdiv_q_ui(g->hard_lower, g->target, factor);
    mpz_mul_ui(g->hard_upper, g->target, factor);
#if HAVE_UINT128 && ULONG_MAX > UINT32_MAX
    siqs_A_native_guard(g);
#endif
    found = siqs_adaptive_choose_A(ctx, poly);
    if (found) break;
  }
  mpz_cdiv_q_ui(g->hard_lower, g->target, 2);
  mpz_mul_ui(g->hard_upper, g->target, 2);
#if HAVE_UINT128 && ULONG_MAX > UINT32_MAX
  siqs_A_native_guard(g);
#endif
  return found;
}

#if USE_CUBIC_CHOOSER
static int siqs_choose_A_cubic(siqs_ctx_t *ctx, siqs_poly_t *poly) {
  if (poly->A_search == NULL)
    poly->A_search = siqs_choose_A_prepare(ctx, poly, 1);
  return siqs_adaptive_choose_A(ctx, poly);
}
#endif

/*----------------------------------------------------------------------------
 * Policy and lifecycle: entry points used by siqs.c
 *----------------------------------------------------------------------------*/

static int siqs_choose_A(siqs_ctx_t *ctx, siqs_poly_t *poly) {

  /* At q=1 we're choosing only one prime. */
  if (poly->q_count == 1U)
    return siqs_choose_A_deterministic(ctx, poly);

  /* SAR also works at q=2. This policy selects deterministic q2 for its
   * tighter A distances and small measured factoring-time advantage. */
  if (poly->q_count == 2U)
    return siqs_choose_A_deterministic(ctx, poly);

  /* The shared deterministic search gave excellent A distances in
   * 50- to 80-bit tests.  That didn't translate into a performance or
   * recovery advantage over SAR, and it is slower to run. */
  if (poly->q_count == 3U)
    return siqs_choose_A_sar(ctx, poly);

  if (poly->q_count == 4U)
    return siqs_choose_A_sar(ctx, poly);

  /* q >= 5 use either Cubic or SAR */
#if USE_CUBIC_CHOOSER
  if (poly->q_count >= SIQS_A_CUBIC_MIN_Q)
    return siqs_choose_A_cubic(ctx, poly);
#endif
  return siqs_choose_A_sar(ctx, poly);
}

static void siqs_choose_A_clear(siqs_poly_t *poly) {
  if (poly->q2_search != NULL) {
    siqs_A_pairs_clear(poly->q2_search);
    free(poly->q2_search);
    poly->q2_search = NULL;
  }
  siqs_A_search_clear(poly);
}
