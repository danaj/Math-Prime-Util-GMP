/* Relation suite, included by siqs-check.c; not a separate compilation unit.
 * Small valid congruences exercise ownership, cycle paths, and partitioning.
 * Artificial large exponents separately exercise arithmetic limits. */
#define RELATION_CANARY UINT32_C(0x91abcdef)

static void relation_context(siqs_ctx_t *ctx, siqs_factor_array_t *result,
                             uint32_t mode) {
  memset(ctx, 0, sizeof(*ctx));
  mpz_init_set_ui(ctx->n, 35);
  mpz_init_set(ctx->kn, ctx->n);
  mpz_init(ctx->eval.y); mpz_init(ctx->eval.q); mpz_init(ctx->eval.rest);
  ctx->original_n = ctx->n;
  ctx->params.fb_size = 2;
  ctx->params.max_large_primes = mode;
  ctx->fb = (siqs_fb_t *)check_allocate(2, sizeof(*ctx->fb));
  ctx->fb[0].p = 2; ctx->fb[1].p = 3;
  ctx->factor_counts = (uint32_t *)check_allocate(4, sizeof(uint32_t));
  ctx->factor_touched = (uint32_t *)check_allocate(4, sizeof(uint32_t));
  ctx->factor_counts[3] = ctx->factor_touched[3] = RELATION_CANARY;
  siqs_factor_array_init(result, ctx->n);
  ctx->result = result;
  siqs_graph_init(&ctx->graph);
  if (mode == 1) siqs_one_lp_init(&ctx->one_lp);
  siqs_hashset_init(&ctx->relation_hashes, 16);
  case_ctx = ctx;
}

static void relation_scratch(const siqs_ctx_t *ctx) {
  uint32_t i;
  CHECK(ctx->factor_touched_count == 0);
  for (i = 0; i < 3; i++) CHECK(ctx->factor_counts[i] == 0);
  CHECK(ctx->factor_counts[3] == RELATION_CANARY);
  CHECK(ctx->factor_touched[3] == RELATION_CANARY);
}

/* Independent congruence arithmetic, including the sign row. */
static void relation_product(const siqs_ctx_t *ctx, mpz_t product,
                             const siqs_factor_t *factors, uint32_t count) {
  uint32_t i;
  mpz_t power;
  mpz_init(power);
  mpz_set_ui(product, 1);
  for (i = 0; i < count; i++) {
    CHECK(factors[i].row <= ctx->params.fb_size);
    if (i != 0) CHECK(factors[i - 1U].row < factors[i].row);
    if (factors[i].row == 0) mpz_set_si(power, -1);
    else mpz_set_ui(power, ctx->fb[factors[i].row - 1U].p);
    mpz_powm_ui(power, power, factors[i].exponent, ctx->n);
    mpz_mul(product, product, power);
    mpz_mod(product, product, ctx->n);
  }
  mpz_clear(power);
}

static void relation_verify_full(const siqs_ctx_t *ctx,
                                 const siqs_full_relation_t *r) {
  mpz_t lhs, rhs;
  mpz_init(lhs); mpz_init(rhs);
  relation_product(ctx, rhs, r->factors, r->nfactors);
  mpz_mul(lhs, r->y, r->y); mpz_mod(lhs, lhs, ctx->n);
  CHECK(mpz_cmp(lhs, rhs) == 0);
  mpz_clear(lhs); mpz_clear(rhs);
}

/* Find a square root by exhaustive search modulo 35, not by using SIQS.
 * Also check the packed/wide accessor round trip before accepting the raw. */
static siqs_raw_relation_t *relation_raw(siqs_ctx_t *ctx, uint64_t lp1,
    uint64_t lp2, const siqs_factor_t *factors, uint32_t count) {
  siqs_raw_relation_t *raw;
  uint32_t i, root;
  mpz_t rhs, lp, y;
  mpz_init(rhs); mpz_init(lp); mpz_init(y);
  relation_product(ctx, rhs, factors, count);
  siqs_mpz_set_u64(lp, lp1); mpz_mul(rhs, rhs, lp);
  siqs_mpz_set_u64(lp, lp2); mpz_mul(rhs, rhs, lp);
  mpz_mod(rhs, rhs, ctx->n);
  for (root = 0; root < 35; root++)
    if (root * root % 35U == mpz_get_ui(rhs)) break;
  CHECK(root < 35);
  mpz_set_ui(y, root);
  raw = siqs_raw_relation_new(ctx, y, factors, count, lp1, lp2);
  CHECK(raw->lp1 == (lp1 < lp2 ? lp1 : lp2));
  CHECK(raw->lp2 == (lp1 < lp2 ? lp2 : lp1));
  for (i = 0; i < count; i++) {
    CHECK(siqs_raw_factor_row(raw, i) == factors[i].row);
    CHECK(siqs_raw_factor_exponent(raw, i) == factors[i].exponent);
  }
  mpz_clear(rhs); mpz_clear(lp); mpz_clear(y);
  return raw;
}

static void relation_finish(siqs_ctx_t *ctx, siqs_factor_array_t *result) {
  uint32_t i, count;
  mpz_t product, *values;
  relation_scratch(ctx);
  for (i = 0; i < ctx->full_count; i++) relation_verify_full(ctx, ctx->full[i]);
  /* Check the partition independently, before releasing its owning context. */
  mpz_init_set_ui(product, 1);
  for (i = 0; i < result->count; i++) {
    CHECK(mpz_cmp_ui(result->values[i], 1) > 0);
    mpz_mul(product, product, result->values[i]);
  }
  CHECK(mpz_cmp(product, ctx->n) == 0);
  mpz_clear(product);
  siqs_ctx_clear(ctx);
  values = siqs_factor_array_release(result, &count);
  gmp_siqs_free(values, count);
  case_ctx = NULL;
}

static void relation_smooth_and_pairs(void) {
  siqs_ctx_t ctx;
  siqs_factor_array_t result;
  siqs_raw_relation_t *r, *anchor;
  siqs_factor_t smooth[3] = {{0, 2}, {1, 2}, {2, 2}};
  siqs_factor_t left[1] = {{1, 2}}, right[2] = {{1, 2}, {2, 2}};
  siqs_factor_t *moved;
  size_t used;
  case_name = "smooth-and-one-LP";
  relation_context(&ctx, &result, 1);
  r = relation_raw(&ctx, 1, 1, smooth, 3);
  moved = r->factors.wide;
  siqs_accept_raw_relation(&ctx, r);
  CHECK(ctx.full_count == 1 && ctx.full[0]->factors == moved);
  CHECK(ctx.raw_count == 0 && ctx.accepted_smooth == 1);
  anchor = relation_raw(&ctx, 11, 1, left, 1);
  siqs_accept_raw_relation(&ctx, anchor);
  used = ctx.raw_arena.current->used;
  /* Exact duplicate rejection must neither add an anchor nor leak a slot. */
  siqs_accept_raw_relation(&ctx, relation_raw(&ctx, 1, 11, left, 1));
  CHECK(ctx.raw_arena.current->used == used);
  CHECK(ctx.one_lp.count == 1 && ctx.accepted_one_lp == 1);
  siqs_accept_raw_relation(&ctx, relation_raw(&ctx, 1, 11, right, 2));
  CHECK(ctx.raw_arena.current->used == used && ctx.full_count == 2);
  CHECK(ctx.full[1]->nfactors == 2);
  CHECK(ctx.full[1]->factors[0].exponent == 4);
  CHECK(ctx.full[1]->factors[1].exponent == 2);
  /* A third distinct partial still pairs with the original anchor. */
  siqs_accept_raw_relation(&ctx, relation_raw(&ctx, 1, 11, right + 1, 1));
  CHECK(ctx.full_count == 3 && ctx.raw_arena.current->used == used);
  CHECK(ctx.full[2]->factors[0].exponent == 2);
  CHECK(ctx.full[2]->factors[1].exponent == 2);
  CHECK(siqs_one_lp_anchor_relation(&ctx.one_lp, 11, NULL) == anchor);
  relation_finish(&ctx, &result);
}

static void relation_inverse_failures(void) {
  uint32_t mode;
  for (mode = 1; mode <= 2; mode++) {
    siqs_ctx_t ctx;
    siqs_factor_array_t result;
    siqs_factor_t factor = {2, 1};
    case_name = mode == 1 ? "one-LP-inverse-factor" : "cycle-inverse-factor";
    relation_context(&ctx, &result, mode);
    if (mode == 1) {
      siqs_raw_relation_t *r = relation_raw(&ctx, 1, 5, &factor, 1);
      siqs_accept_raw_relation(&ctx, r);
      r = relation_raw(&ctx, 1, 5, &factor, 1);
      mpz_sub(r->y, ctx.n, r->y); /* Other valid root, not a duplicate. */
      siqs_accept_raw_relation(&ctx, r);
    } else {
      siqs_accept_raw_relation(&ctx, relation_raw(&ctx, 5, 5, NULL, 0));
      CHECK(ctx.raw_count == 0);
    }
    CHECK(ctx.factor_found && ctx.full_count == 0 && result.count == 2);
    CHECK(siqs_all_factors_prime(&result));
    CHECK(mpz_cmp_ui(result.values[0], 5) == 0);
    CHECK(mpz_cmp_ui(result.values[1], 7) == 0);
    relation_finish(&ctx, &result);
  }
}

static void relation_general_pairs(void) {
  siqs_ctx_t ctx;
  siqs_factor_array_t result;
  siqs_factor_t left = {1, 2}, right[2] = {{1, 2}, {2, 2}};
  case_name = "two-vertex-cycle-and-smooth-self-loop";
  relation_context(&ctx, &result, 2);
  siqs_accept_raw_relation(&ctx, relation_raw(&ctx, 1, 1, &left, 1));
  CHECK(ctx.raw_count == 0 && ctx.full_count == 1);
  siqs_accept_raw_relation(&ctx, relation_raw(&ctx, 1, 11, &left, 1));
  siqs_accept_raw_relation(&ctx, relation_raw(&ctx, 1, 11, right, 2));
  CHECK(ctx.raw_count == 1 && ctx.full_count == 2);
  CHECK(ctx.full[1]->nfactors == 2 && ctx.full[1]->factors[0].exponent == 4);
  CHECK(ctx.full[1]->factors[1].exponent == 2);
  siqs_accept_raw_relation(&ctx, relation_raw(&ctx, 1, 11, right + 1, 1));
  CHECK(ctx.raw_count == 1 && ctx.full_count == 3);
  CHECK(ctx.full[2]->nfactors == 2 && ctx.full[2]->factors[0].exponent == 2);
  CHECK(ctx.full[2]->factors[1].exponent == 2);
  relation_finish(&ctx, &result);
}

static void relation_graph_paths(void) {
  uint32_t pass, i, length = extended ? 160U : 80U;
  uint64_t *labels = (uint64_t *)check_allocate(length + 1U, sizeof(uint64_t));
  mpz_t prime;
  mpz_init_set_ui(prime, 35);
  labels[0] = 1;
  for (i = 1; i <= length; i++) {
    unsigned long residue;
    /* Unit squares modulo 35, not just residue 1: division by the LP
     * product must actually work for the congruence checks to pass. */
    do {
      mpz_nextprime(prime, prime);
      residue = mpz_fdiv_ui(prime, 35);
    } while (residue != 1 && residue != 4 && residue != 9 &&
             residue != 11 && residue != 16 && residue != 29);
    labels[i] = (uint64_t)mpz_get_ui(prime);
  }
  mpz_clear(prime);
  /* Reverse/interleaved forest insertion changes rerooting, not the cycles.
   * This models relations arriving in different worker merge orders. */
  for (pass = 0; pass < 3; pass++) {
    siqs_ctx_t ctx;
    siqs_factor_array_t result;
    siqs_factor_t factor = {1, 2};
    uint32_t left[4] = {0, 10, 0, 10};
    uint32_t right[4];
    right[0] = length; right[1] = length - 10U;
    right[2] = length - 1U; right[3] = 10;
    case_name = "overlapping-rerooted-cycles";
    relation_context(&ctx, &result, 2);
    for (i = 0; i < length; i++) {
      uint32_t edge = pass == 0 ? i : pass == 1 ? length - 1U - i
                      : i < length / 2U ? 2U * i : 2U * (i - length / 2U) + 1U;
      siqs_accept_raw_relation(&ctx,
          relation_raw(&ctx, labels[edge], labels[edge + 1U], &factor, 1));
    }
    CHECK(ctx.raw_count == length && ctx.full_count == 0);
    for (i = 0; i < 4; i++) {
      siqs_raw_block_t *block = ctx.raw_arena.current;
      size_t used = block->used;
      if (i == 3) ctx.graph.mark = UINT32_MAX; /* Exercise cycle-mark wrap. */
      siqs_accept_raw_relation(&ctx,
          relation_raw(&ctx, labels[left[i]], labels[right[i]], &factor, 1));
      CHECK(ctx.raw_count == length && ctx.full_count == i + 1U);
      CHECK(ctx.raw[length] == NULL);
      CHECK(ctx.raw_arena.current == block && block->used == used);
      CHECK(ctx.full[i]->nfactors == 1);
      CHECK(ctx.full[i]->factors[0].exponent == 2U * (right[i] - left[i] + 1U));
      relation_verify_full(&ctx, ctx.full[i]);
      relation_scratch(&ctx);
    }
    CHECK(ctx.graph.cycle_edge_alloc > 64 && ctx.graph.cycle_vertex_alloc > 64);
    relation_finish(&ctx, &result);
  }
  free(labels);
}

static void relation_arena_and_partition(void) {
  siqs_ctx_t ctx;
  siqs_factor_array_t result;
  siqs_factor_t f = {1, 2};
  siqs_raw_relation_t *first, *second;
  size_t first_end, both_end;
  uint32_t count;
  siqs_factor_t limit = {SIQS_PACKED_FACTOR_ROW_MASK, SIQS_PACKED_FACTOR_EXP_MAX};
  mpz_t n, divisor, product, *values;
  case_name = "arena-non-LIFO-and-partition";
  relation_context(&ctx, &result, 2);
  first = relation_raw(&ctx, 1, 71, &f, 1);
  first_end = ctx.raw_arena.current->used;
  second = relation_raw(&ctx, 1, 211, &f, 1);
  both_end = ctx.raw_arena.current->used;
  siqs_raw_relation_free(&ctx, first);
  CHECK(ctx.raw_arena.current->used == both_end);
  siqs_raw_relation_free(&ctx, second);
  CHECK(ctx.raw_arena.current->used == first_end);
  /* Storage-only endpoint fixture, not an accepted relation congruence. */
  first = siqs_raw_relation_new(&ctx, ctx.eval.y, &limit, 1, 1, 71);
  CHECK(siqs_raw_factor_row(first, 0) == SIQS_PACKED_FACTOR_ROW_MASK);
  CHECK(siqs_raw_factor_exponent(first, 0) == SIQS_PACKED_FACTOR_EXP_MAX);
  siqs_raw_relation_free(&ctx, first);
  CHECK(ctx.raw_arena.current->used == first_end);
  relation_finish(&ctx, &result);
  mpz_init_set_ui(n, 1575); /* 3^2 * 5^2 * 7: repeated/composite divisors. */
  mpz_init(divisor); mpz_init(product);
  siqs_factor_array_init(&result, n);
  mpz_set_ui(divisor, 15); CHECK(siqs_insert_divisor(&result, divisor));
  CHECK(!siqs_insert_divisor(&result, divisor));
  CHECK(!siqs_all_factors_prime(&result)); /* Cache composite status before refinement. */
  mpz_set_ui(divisor, 3); CHECK(siqs_insert_divisor(&result, divisor));
  CHECK(!siqs_insert_divisor(&result, divisor));
  CHECK(siqs_all_factors_prime(&result));
  CHECK(siqs_all_factors_prime(&result)); /* Repeated cache use. */
  CHECK(result.count == 5);
  mpz_set_ui(product, 1);
  for (count = 0; count < result.count; count++) {
    CHECK(mpz_probab_prime_p(result.values[count], 25) != 0);
    mpz_mul(product, product, result.values[count]);
  }
  CHECK(mpz_cmp(product, n) == 0);
  values = siqs_factor_array_release(&result, &count); gmp_siqs_free(values, count);
  mpz_clear(n); mpz_clear(divisor); mpz_clear(product);
}

static void relation_anchor_growth(void) {
  siqs_ctx_t ctx;
  siqs_factor_array_t result;
  siqs_factor_t left = {1, 2}, right = {2, 2};
  mpz_t prime;
  uint32_t i;
  uint64_t last = 0;
  case_name = "anchor-rehash-and-arena-growth";
  relation_context(&ctx, &result, 1);
  mpz_init_set_ui(prime, 35);
  for (i = 0; i < 900; i++) {
    do { mpz_nextprime(prime, prime); } while (mpz_fdiv_ui(prime, 35) != 1);
    last = (uint64_t)mpz_get_ui(prime);
    siqs_accept_raw_relation(&ctx, relation_raw(&ctx, 1, last, &left, 1));
  }
  CHECK(ctx.one_lp.count == 900 && ctx.one_lp.alloc > 1024);
  CHECK(ctx.raw_arena.current->previous != NULL);
  siqs_accept_raw_relation(&ctx, relation_raw(&ctx, 1, last, &right, 1));
  CHECK(ctx.full_count == 1 && ctx.full[0]->nfactors == 2);
  mpz_clear(prime);
  relation_finish(&ctx, &result);
}

static siqs_full_relation_t *relation_full_power(siqs_ctx_t *ctx,
                                                uint32_t row, uint32_t exp) {
  siqs_full_relation_t *r = (siqs_full_relation_t *)check_allocate(1, sizeof(*r));
  CHECK((exp & 1U) == 0 && row != 0);
  mpz_init_set_ui(r->y, ctx->fb[row - 1U].p);
  mpz_powm_ui(r->y, r->y, exp / 2U, ctx->n);
  r->factors = (siqs_factor_t *)check_allocate(1, sizeof(*r->factors));
  r->nfactors = 1; r->factors[0].row = row; r->factors[0].exponent = exp;
  siqs_store_full(ctx, r);
  return r;
}

static void relation_exponent_limits(void) {
  siqs_ctx_t ctx;
  siqs_factor_array_t result;
  siqs_factor_t f = {1, 2};
  siqs_raw_relation_t *raw[2];
  uint32_t edges[2] = {0, 1}, i;
  la_col_t columns[3];
  uint64_t nullrows[3] = {UINT64_C(1), UINT64_C(1), UINT64_C(2)};
  case_name = "exponent-overflow-and-recovery";
  relation_context(&ctx, &result, 2);
  CHECK(siqs_touch_factor(&ctx, 1, 0) && ctx.factor_touched_count == 0);
  CHECK(siqs_touch_factor(&ctx, 1, UINT32_MAX - 1U));
  CHECK(siqs_touch_factor(&ctx, 1, 1));
  CHECK(!siqs_touch_factor(&ctx, 1, 1));
  CHECK(!siqs_touch_factor(&ctx, 1, 1));
  CHECK(ctx.factor_counts[1] == UINT32_MAX && ctx.factor_touched_count == 1);
  siqs_reset_touched_factors(&ctx);
  /* Synthetic wide vectors bypass packed raw's 12-bit exponent limit. */
  for (i = 0; i < 2; i++) {
    raw[i] = relation_raw(&ctx, 1, 1, &f, 1);
    siqs_store_raw(&ctx, raw[i]);
    raw[i]->factors.wide[0].exponent = UINT32_C(2147483648) - (i == 0 ? 2U : 0U);
    mpz_set_ui(raw[i]->y, 2);
    mpz_powm_ui(raw[i]->y, raw[i]->y, raw[i]->factors.wide[0].exponent / 2U, ctx.n);
  }
  siqs_materialize_cycle(&ctx, edges, 2, NULL, 0);
  CHECK(ctx.full_count == 1 && ctx.full[0]->factors[0].exponent == UINT32_MAX - 1U);
  raw[0]->factors.wide[0].exponent += 2U;
  mpz_set(raw[0]->y, raw[1]->y);
  siqs_materialize_cycle(&ctx, edges, 2, NULL, 0);
  CHECK(ctx.full_count == 1 && !ctx.factor_found); relation_scratch(&ctx);
  for (i = 0; i < 2; i++) {
    raw[i]->factors.wide[0].exponent = 2;
    mpz_set_ui(raw[i]->y, 2);
  }
  siqs_materialize_cycle(&ctx, edges, 2, NULL, 0);
  CHECK(ctx.full_count == 2 && ctx.full[1]->factors[0].exponent == 4);
  relation_finish(&ctx, &result);
  relation_context(&ctx, &result, 2);
  memset(columns, 0, sizeof(columns));
  for (i = 0; i < 3; i++) columns[i].orig = i;
  (void)relation_full_power(&ctx, 1, UINT32_C(2147483648));
  (void)relation_full_power(&ctx, 1, UINT32_C(2147483648));
  mpz_set_ui(relation_full_power(&ctx, 2, 2)->y, 18); /* 18^2 == 3^2 mod 35. */
  CHECK(!siqs_test_dependencies(&ctx, columns, 3, nullrows, UINT64_C(1)));
  CHECK(!ctx.factor_found && result.count == 1); relation_scratch(&ctx);
  CHECK(siqs_test_dependencies(&ctx, columns, 3, nullrows, UINT64_C(3)));
  CHECK(ctx.factor_found && result.count == 2 && siqs_all_factors_prime(&result));
  relation_finish(&ctx, &result);
}

static void suite_relations(void) {
  relation_smooth_and_pairs();
  relation_general_pairs();
  relation_inverse_failures();
  relation_graph_paths();
  relation_arena_and_partition();
  relation_exponent_limits();
  if (extended) relation_anchor_growth();
  puts("PASS relations: congruences, shared/rerooted cycles, ownership, inverse factors, exponent limits");
  fflush(stdout);
}
