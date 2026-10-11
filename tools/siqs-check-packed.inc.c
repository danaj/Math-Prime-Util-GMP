/* Packed-factor checks, included by the relation suite.
 *
 * The envelope is deliberately independent of tuning:
 *   N < 2^512, k < 2^12, q <= 16, p_i < 2^32, |x| <= M <= 2^30.
 *   A < 2^512; H_i <= A/2; |B| <= (q/2+1)A <= 9A.
 *   |C| <= 81A+kN and |Q| <= 2AM^2+18AM+81A+kN < 2^574.
 * For nonzero Q, v_p(Q) <= 573. Even assigning all q+1 A/d extras to
 * one prime gives <=590. This is a bound on one raw entry, not cycle sums.
 *
 * Large-N fixtures below construct genuine integral polynomials and check
 * their identities/roots, but use tiny factor bases and do no collection.
 * They do not validate large-input factoring performance or reliability.
 */

static void relation_packed_envelope(void) {
  static const uint32_t multipliers[] = {1U, 2U, 2047U, 4095U};
  mpz_t n, kn, a, bound, term, limit;
  uint32_t bits, q, k, d, cases = 0, maximum = 0;
  case_name = "512-bit-raw-exponent-envelope";
  CHECK(SIQS_PACKED_FACTOR_ROW_BITS == 21U);
  CHECK(SIQS_PACKED_FACTOR_ROW_MASK == 2097151U);
  CHECK(SIQS_PACKED_FACTOR_EXP_MAX == 2047U);
  CHECK(SIQS_MAX_A_PRIMES == 16U);
  CHECK(SIQS_EVAL_MAX_EXTRA_FACTORS == SIQS_MAX_A_PRIMES + 2U);
  CHECK(SIQS_HALF_INTERVAL_MAX <= (UINT32_C(1) << 30));
  CHECK(MPU_SIQS_MAX_BITS <= 512U);
  CHECK(SIQS_MULTIPLIER_MAX <= 4095U && SIQS_MULTIPLIER_EVEN_MAX <= 4094U);
  mpz_inits(n, kn, a, bound, term, limit, NULL);
  mpz_set_ui(limit, 1); mpz_mul_2exp(limit, limit, 574);
  for (bits = 1; bits <= 512; bits++) {
    mpz_set_ui(n, 1); mpz_mul_2exp(n, n, bits); mpz_sub_ui(n, n, 1);
    for (q = 1; q <= SIQS_MAX_A_PRIMES; q++) {
      mpz_set_ui(a, UINT32_MAX); mpz_pow_ui(a, a, q);
      for (k = 0; k < sizeof(multipliers) / sizeof(*multipliers); k++) {
        mpz_mul_ui(kn, n, multipliers[k]);
        if (bits == 512 && multipliers[k] == 2047U)
          CHECK(mpz_sizeinbase(kn, 2) == 523U);
        if (bits == 512 && multipliers[k] == 4095U)
          CHECK(mpz_sizeinbase(kn, 2) == 524U);
        for (d = 1; d <= 2; d++) {
          uint32_t exponent;
          mpz_mul_ui(bound, a, d);
          mpz_mul_ui(bound, bound, SIQS_HALF_INTERVAL_MAX);
          mpz_mul_ui(bound, bound, SIQS_HALF_INTERVAL_MAX);
          mpz_mul_ui(term, a, 18U);
          mpz_mul_ui(term, term, SIQS_HALF_INTERVAL_MAX);
          mpz_add(bound, bound, term);
          mpz_mul_ui(term, a, 81U);
          mpz_add(bound, bound, term); mpz_add(bound, bound, kn);
          CHECK(mpz_cmp(bound, limit) < 0);
          exponent = (uint32_t)mpz_sizeinbase(bound, 2) - 1U + q + 1U;
          CHECK(exponent <= 590U && exponent <= SIQS_PACKED_FACTOR_EXP_MAX);
          if (exponent > maximum) maximum = exponent;
          cases++;
        }
      }
    }
  }
  CHECK(maximum == 590U);
  mpz_clears(n, kn, a, bound, term, limit, NULL);
  printf("PASS relations: %u exact 1-512-bit envelope cases; raw exponent <=%u <2048\n",
         cases, maximum);
}

static void relation_packed_boundaries(void) {
  static const uint32_t rows[] = {0U, 1U, 2097150U, 2097151U};
  static const uint32_t exponents[] = {0U, 1U, 590U, 2046U, 2047U};
  siqs_raw_relation_t raw;
  mpz_t power;
  uint32_t word, i, j;
  case_name = "21-11-storage-boundaries";
  memset(&raw, 0, sizeof(raw));
  raw.lp2 = 71U; raw.nfactors = 1U; raw.factors.packed = &word;
  for (i = 0; i < sizeof(rows) / sizeof(*rows); i++) {
    for (j = 0; j < sizeof(exponents) / sizeof(*exponents); j++) {
      word = siqs_pack_raw_factor(rows[i], exponents[j]);
      CHECK(siqs_raw_factor_row(&raw, 0) == rows[i]);
      CHECK(siqs_raw_factor_exponent(&raw, 0) == exponents[j]);
    }
  }
  CHECK(siqs_pack_raw_factor(2097151U, 2047U) == UINT32_MAX);
  /* Synthetic valuation endpoint, not a claimed polynomial relation. */
  mpz_init_set_ui(power, 1); mpz_mul_2exp(power, power, 573U);
  CHECK(siqs_remove_ui(power, 2U) == 573U && mpz_cmp_ui(power, 1U) == 0);
  mpz_clear(power);
  puts("PASS relations: 21/11 row/exponent endpoints, combined UINT32_MAX word and v2=573");
}

static void relation_large_polynomial(uint32_t q, uint32_t d,
                                      uint64_t *count, uint32_t *max_bits) {
  siqs_ctx_t ctx;
  siqs_poly_t poly;
  siqs_fb_t fb[SIQS_MAX_A_PRIMES + 2U];
  uint32_t primes[SIQS_MAX_A_PRIMES + 2U];
  uint32_t root1[SIQS_MAX_A_PRIMES + 2U], root2[SIQS_MAX_A_PRIMES + 2U];
  uint8_t in_a[SIQS_MAX_A_PRIMES + 2U];
  mpz_t prime, expected_b, offset, bound, y, value, limit;
  uint32_t i;
  int point;
  memset(&ctx, 0, sizeof(ctx)); memset(fb, 0, sizeof(fb));
  memset(in_a, 0, sizeof(in_a));
  ctx.params.bits = 512U; ctx.params.q_count = q; ctx.params.poly_d = d;
  ctx.params.fb_size = q + 2U;
  ctx.params.half_interval = SIQS_HALF_INTERVAL_MAX;
  ctx.sieve_length = 2U * ctx.params.half_interval;
  ctx.multiplier = 2047U;
  ctx.fb = fb; ctx.prime = primes; ctx.fb_in_a = in_a;
  ctx.root1 = root1; ctx.root2 = root2;
  mpz_init_set_ui(ctx.n, 1); mpz_mul_2exp(ctx.n, ctx.n, 512U);
  mpz_sub_ui(ctx.n, ctx.n, d == 2U ? 1U : 3U);
  mpz_init(ctx.kn); mpz_mul_ui(ctx.kn, ctx.n, ctx.multiplier);
  CHECK((mpz_fdiv_ui(ctx.kn, 8U) == 1U) == (d == 2U));
  mpz_inits(prime, expected_b, offset, bound, y, value, limit, NULL);
  mpz_set_ui(prime, (UINT32_C(1) << 30) - 1024U);
  fb[0].p = primes[0] = 2U;
  /* Include one ordinary two-root prime as well as 2 and all A primes. */
  for (i = 1; i <= q + 1U; i++) {
    uint32_t p, residue;
    do {
      mpz_nextprime(prime, prime);
      CHECK(mpz_cmp_ui(prime, UINT32_C(1) << 31) < 0);
      p = (uint32_t)mpz_get_ui(prime);
      residue = (uint32_t)mpz_fdiv_ui(ctx.kn, p);
    } while (residue == 0 || siqs_legendre_u32(residue, p) != 1);
    fb[i].p = primes[i] = p;
    fb[i].sqrt_kn = siqs_sqrtmod_u32(residue, p);
    CHECK(fb[i].sqrt_kn != 0);
    in_a[i] = i <= q;
  }
  case_ctx = &ctx;
  siqs_poly_init(&ctx, &poly);
  mpz_set_ui(poly.A, 1U);
  for (i = 0; i < q; i++) {
    poly.a_index[i] = i + 1U;
    mpz_mul_ui(poly.A, poly.A, primes[i + 1U]);
  }
  mpz_mul_ui(poly.DA, poly.A, d);
  siqs_first_B_and_roots(&ctx, &poly);
  mpz_set(offset, poly.B);
  for (i = 0; i < q; i++) {
    CHECK(mpz_sgn(poly.H[i]) >= 0);
    mpz_mul_ui(value, poly.H[i], 2U);
    CHECK(mpz_cmp(value, poly.A) <= 0);
    mpz_sub(offset, offset, poly.H[i]);
  }
  CHECK(mpz_sgn(offset) == 0 || mpz_cmp(offset, poly.A) == 0);
  mpz_mul_ui(bound, poly.A, 9U);
  mpz_set_ui(limit, 1U); mpz_mul_2exp(limit, limit, 574U);
  do {
    uint32_t gray = poly.b_index ^ (poly.b_index >> 1);
    mpz_set(expected_b, offset);
    for (i = 0; i < q; i++) {
      if (gray & (1U << i)) mpz_sub(expected_b, expected_b, poly.H[i]);
      else mpz_add(expected_b, expected_b, poly.H[i]);
    }
    CHECK(mpz_cmp(expected_b, poly.B) == 0);
    mpz_abs(value, poly.B); CHECK(mpz_cmp(value, bound) <= 0);
    for (point = -1; point <= 1; point++) {
      int32_t x = point * (int32_t)ctx.params.half_interval;
      uint32_t width;
      /* Independent squared-numerator identity, not Horner evaluation. */
      mpz_mul_si(y, poly.DA, x); mpz_add(y, y, poly.B);
      mpz_mul(value, y, y); mpz_sub(value, value, ctx.kn);
      CHECK(mpz_divisible_p(value, poly.DA));
      mpz_divexact(value, value, poly.DA);
      if (point == 0) CHECK(mpz_cmp(value, poly.C) == 0);
      mpz_abs(value, value);
      CHECK(mpz_cmp(value, limit) < 0);
      width = (uint32_t)mpz_sizeinbase(value, 2);
      if (width > *max_bits) *max_bits = width;
    }
    /* Independent root spot checks also run without SIQS_DEBUG. */
    if ((poly.b_index & 1023U) == 0U || poly.b_index + 1U == poly.b_limit) {
      uint32_t which;
      for (i = 0; i < ctx.params.fb_size; i++) {
        uint32_t p = primes[i];
        for (which = 0; which < 2U; which++) {
          uint32_t root = which == 0U ? root1[i] : root2[i];
          uint32_t x;
          CHECK(root < p);
          x = (uint32_t)(((uint64_t)root + p - ctx.params.half_interval % p) % p);
          mpz_mul_ui(value, poly.DA, x);
          mpz_addmul_ui(value, poly.B, 2U);
          mpz_mul_ui(value, value, x); mpz_add(value, value, poly.C);
          CHECK(mpz_divisible_ui_p(value, p));
        }
      }
    }
    (*count)++;
  } while (siqs_next_B(&ctx, &poly));
  CHECK(poly.b_index + 1U == poly.b_limit);
  CHECK(!siqs_next_B(&ctx, &poly));
  siqs_poly_clear(&ctx, &poly);
  mpz_clears(ctx.n, ctx.kn, prime, expected_b, offset, bound, y, value, limit, NULL);
  case_ctx = NULL;
}

static void relation_large_polynomials(void) {
  static const uint32_t qs[] = {1U, 2U, 3U, 13U, 16U};
  uint64_t count = 0;
  uint32_t i, q, d, maximum = 0;
  case_name = "512-bit-real-polynomial-helpers";
  if (extended) {
    for (q = 1; q <= SIQS_MAX_A_PRIMES; q++)
      for (d = 1; d <= 2; d++) relation_large_polynomial(q, d, &count, &maximum);
  } else {
    for (i = 0; i < sizeof(qs) / sizeof(*qs); i++)
      for (d = 1; d <= 2; d++) relation_large_polynomial(qs[i], d, &count, &maximum);
  }
  printf("PASS relations: 512-bit tiny-FB fixtures, %llu Gray polynomials, both d modes, max |Q| %u bits\n",
         (unsigned long long)count, maximum);
}

static void relation_zero_polynomial(void) {
  siqs_ctx_t ctx;
  siqs_factor_array_t result;
  siqs_poly_t poly;
  case_name = "zero-Q-bypasses-removal-and-packing";
  relation_context(&ctx, &result, 2U);
  ctx.params.q_count = ctx.params.poly_d = 1U;
  ctx.params.half_interval = 16U;
  ctx.fb_in_a = (uint8_t *)check_allocate(2U, sizeof(*ctx.fb_in_a));
  ctx.fb_in_a[1] = 1U;
  ctx.multiplier = 35U; mpz_mul_ui(ctx.kn, ctx.n, ctx.multiplier);
  ctx.candidates = (siqs_candidate_t *)check_allocate(1U, sizeof(*ctx.candidates));
  ctx.candidates[0].x = -12;
  ctx.candidates[0].first_hit = SIQS_NO_INDEX;
  siqs_poly_init(&ctx, &poly);
  poly.a_index[0] = 1U;
  mpz_set_ui(poly.A, 3U); mpz_set_ui(poly.DA, 3U); mpz_set_ui(poly.B, 1U);
  siqs_compute_C(&ctx, &poly); CHECK(mpz_cmp_si(poly.C, -408) == 0);
  siqs_evaluate_candidate(&ctx, &poly, 0U);
  CHECK(mpz_sgn(ctx.eval.q) == 0 && ctx.split_attempts == 0);
  CHECK(ctx.full_count == 0 && !ctx.factor_found && result.count == 1U);
  siqs_poly_clear(&ctx, &poly);
  relation_finish(&ctx, &result);
  puts("PASS relations: zero Q exits before factor removal, splitting or packed storage");
}

static void relation_packed_failures(void) {
#ifndef _WIN32
  uint32_t fixture;
  case_name = "packed-and-polynomial-envelope-rejections";
  fflush(NULL);
  for (fixture = 0; fixture < 11U; fixture++) {
    int status;
    pid_t child = fork();
    CHECK(child >= 0);
    if (child == 0) {
      siqs_ctx_t ctx;
      siqs_poly_t poly;
      uint32_t row = 1U, exponent = 1U;
      if (freopen("/dev/null", "w", stderr) == NULL) _exit(1);
      memset(&ctx, 0, sizeof(ctx));
      ctx.params.fb_size = ctx.params.q_count = ctx.params.poly_d = 1U;
      ctx.params.half_interval = 1U;
      mpz_init_set_ui(ctx.kn, 1U);
      switch (fixture) {
        case 0: row = SIQS_PACKED_FACTOR_ROW_MASK + 1U; break;
        case 1: exponent = SIQS_PACKED_FACTOR_EXP_MAX + 1U; break;
        case 2: row = SIQS_PACKED_FACTOR_ROW_MASK + 1U;
                exponent = SIQS_PACKED_FACTOR_EXP_MAX + 1U; break;
        case 3: row = UINT32_MAX; break;
        case 4: exponent = UINT32_MAX; break;
        case 5: ctx.params.q_count = 0U; break;
        case 6: ctx.params.q_count = SIQS_MAX_A_PRIMES + 1U; break;
        case 7: ctx.params.half_interval = 0U; break;
        case 8: ctx.params.half_interval = SIQS_HALF_INTERVAL_MAX + 1U; break;
        case 9: ctx.params.poly_d = 0U; break;
        case 10: ctx.params.poly_d = 3U; break;
      }
      if (fixture < 5U) (void)siqs_pack_raw_factor(row, exponent);
      else siqs_poly_init(&ctx, &poly);
      _exit(0);
    }
    CHECK(waitpid(child, &status, 0) == child);
    CHECK(WIFEXITED(status) && WEXITSTATUS(status) == 3);
  }
  puts("PASS relations: row 2097152, exponent 2048, oversized words and invalid q/M/d rejected");
#else
  puts("SKIP relations: fatal packed/envelope diagnostics require POSIX fork");
#endif
}

