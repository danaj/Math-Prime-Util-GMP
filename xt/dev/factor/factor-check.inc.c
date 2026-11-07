/* Included after the unchanged factor.c prefix by factor-check.pl.
 * Factor-loop methods are scripted; separate p-1 checks use the real code.
 * GMP primality, factor bookkeeping and the prime iterator remain real. */
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

enum { CHECK_ECM, CHECK_PM1, CHECK_PBRENT, CHECK_HOLF, CHECK_CHEB,
       CHECK_SIQS, CHECK_TINYECM64 };
typedef struct { int kind; UV a, b, b2; } check_event_t;
static check_event_t check_events[256];
static size_t check_nevents;
static unsigned check_siqs_calls, check_releases, check_assertions, check_fail_after;
static mpz_t check_input, check_parts[8];
static unsigned check_nparts, check_primality_calls;
static int check_split, check_verbose, check_pm1_split;
static jmp_buf check_failure_jump;
static int check_expecting_failure;
static char check_failure_message[128];
static size_t check_allocations;
static mpz_t *check_failure_factors;
static int *check_failure_exponents;

#define CHECK(c) do { check_assertions++; if (!(c)) { \
  fprintf(stderr, "FAIL line %d: %s\n", __LINE__, #c); exit(1); \
} } while (0)

static void check_croak(const char *format, ...) {
  va_list args;
  va_start(args, format);
  vsnprintf(check_failure_message, sizeof(check_failure_message), format, args);
  va_end(args);
  if (check_expecting_failure) longjmp(check_failure_jump, 1);
  fprintf(stderr, "Unexpected croak: %s\n", check_failure_message);
  exit(1);
}

static int check_printf(const char *format, ...) {
  va_list args;
  int result;
  if (check_expecting_failure && check_verbose < 2) return 0;
  va_start(args, format);
  result = gmp_vprintf(format, args);
  va_end(args);
  return result;
}

/* Count GMP limbs, including on systems whose address sanitizer does not
 * provide leak detection. */
static void *check_malloc(size_t count, size_t size) {
  void *mem = mpu_malloc(count, size);
  if (mem != NULL) check_allocations++;
  return mem;
}
static void *check_realloc(void *mem, size_t count, size_t size) {
  int had_mem = mem != NULL;
  void *result = mpu_realloc(mem, count, size);
  if (!had_mem && result != NULL) check_allocations++;
  if (had_mem && result == NULL) check_allocations--;
  return result;
}
static void check_free(void *mem) {
  if (mem != NULL) {
    CHECK(check_allocations > 0);
    check_allocations--;
  }
  free(mem);
}
static void *check_gmp_malloc(size_t size) { return check_malloc(size, 1); }
static void *check_gmp_realloc(void *mem, size_t old_size, size_t new_size) {
  (void)old_size;
  return check_realloc(mem, new_size, 1);
}
static void check_gmp_free(void *mem, size_t size) {
  (void)size;
  check_free(mem);
}

static void check_terminal_failure(void) {
  size_t allocations = check_allocations;
  check_failure_factors = NULL;
  check_failure_exponents = NULL;
  check_failure_message[0] = '\0';
  check_expecting_failure = 1;
  if (setjmp(check_failure_jump) == 0) {
    (void)factor(check_input, &check_failure_factors, &check_failure_exponents);
    CHECK(0); /* Exhaustion must throw, never return an unsplit composite. */
  }
  check_expecting_failure = 0;
  CHECK(strcmp(check_failure_message, "factor() unable to factor") == 0);
  CHECK(check_failure_factors == NULL && check_failure_exponents == NULL);
  CHECK(check_allocations == allocations);
}

static void check_record(int kind, UV a, UV b) {
  CHECK(check_nevents < sizeof(check_events) / sizeof(check_events[0]));
  check_events[check_nevents].kind = kind;
  check_events[check_nevents].a = a;
  check_events[check_nevents].b2 = 0;
  check_events[check_nevents++].b = b;
}

int get_verbose_level(void) { return check_verbose; }
uint64_t irand_u64(void) { return UINT64_C(0x123456789abcdef0); }
int _GMP_is_prob_prime(const mpz_t n) {
  check_primality_calls++;
  return mpz_probab_prime_p(n, 20) != 0;
}
void mpz_set_uv(mpz_t n, UV v) { mpz_import(n, 1, 1, sizeof(v), 0, 0, &v); }
void mpz_product_ui(mpz_t n, unsigned long *v, unsigned long count) {
  unsigned long i;
  mpz_set_ui(n, 1);
  for (i = 0; i < count; i++) mpz_mul_ui(n, n, v[i]);
}
int _GMP_tinyecm64_factor(const mpz_t n, mpz_t f, int try_small) {
  (void)n; (void)f; check_record(CHECK_TINYECM64, (UV)try_small, 0); return 0;
}
int _GMP_pminus1_factor(const mpz_t n, mpz_t f, UV b1, UV b2) {
  (void)n;
  check_record(CHECK_PM1, b1, b2);
  if (check_pm1_split) {
    unsigned i;
    check_pm1_split--;
    for (i = 3; i < 8; i++)
      if (mpz_cmp(check_parts[i], n) < 0 && mpz_divisible_p(n, check_parts[i])) break;
    CHECK(i < 8);
    mpz_set(f, check_parts[i]);
    return 1;
  }
  return 0;
}
int _GMP_ecm_factor_projective(const mpz_t n, mpz_t f, UV b1, UV b2, UV curves,
                              uint64_t seed, int verbose) {
  (void)n; (void)f;
  CHECK(seed == UINT64_C(0x123456789abcdef0) && verbose == check_verbose);
  check_record(CHECK_ECM, b1, curves);
  check_events[check_nevents-1].b2 = b2; return 0;
}
int _GMP_pbrent_factor(const mpz_t n, mpz_t f, UV a, UV rounds) {
  (void)n; (void)f; check_record(CHECK_PBRENT, a, rounds); return 0;
}
int _GMP_holf_factor(const mpz_t n, mpz_t f, UV rounds) {
  (void)n; (void)f; check_record(CHECK_HOLF, rounds, 0); return 0;
}
int _GMP_cheb_factor(const mpz_t n, mpz_t f, UV bound, UV x) {
  (void)n; (void)f; check_record(CHECK_CHEB, bound, x); return 0;
}

mpz_t *gmp_siqs(const mpz_t n, uint32_t *count, uint32_t trial_start, int verbose) {
  mpz_t *parts;
  unsigned i;
  int split;
  CHECK(trial_start == 64007);
  CHECK(verbose == check_verbose - 2);
  check_siqs_calls++;
  check_record(CHECK_SIQS, (UV)mpz_sizeinbase(n, 2), 0);
  split = check_split && (!check_fail_after || check_siqs_calls <= check_fail_after);
  *count = split ? check_nparts : 1;
  if (*count == 0) return NULL;
  New(0, parts, *count, mpz_t);
  if (!split) {
    mpz_init_set(parts[0], n);
  } else if (mpz_cmp(n, check_input) == 0) {
    for (i = 0; i < *count; i++) mpz_init_set(parts[i], check_parts[i]);
  } else {
    /* Follow-up composite work: split by the first scripted proper divisor. */
    *count = 2;
    for (i = 0; i < 8; i++)
      if (mpz_cmp_ui(check_parts[i], 1) > 0 && mpz_cmp(check_parts[i], n) < 0 &&
          mpz_divisible_p(n, check_parts[i])) break;
    CHECK(i < 8);
    mpz_init_set(parts[0], check_parts[i]);
    mpz_init(parts[1]);
    mpz_divexact(parts[1], n, parts[0]);
  }
  return parts;
}
void gmp_siqs_free(mpz_t *parts, uint32_t count) {
  uint32_t i;
  for (i = 0; i < count; i++) mpz_clear(parts[i]);
  Safefree(parts);
  check_releases++;
}

static void check_reset(void) {
  check_nevents = 0;
  check_siqs_calls = check_releases = 0;
  check_fail_after = 0;
  check_split = 0;
  check_pm1_split = 0;
  check_nparts = 2;
}
static void check_semiprime(unsigned bits) {
  mpz_t q;
  mpz_init(q);
  mpz_set_ui(check_parts[0], 1);
  mpz_mul_2exp(check_parts[0], check_parts[0], bits/2 - 1);
  mpz_nextprime(check_parts[0], check_parts[0]);
  mpz_set_ui(q, 1);
  mpz_mul_2exp(q, q, bits - 1);
  mpz_cdiv_q(q, q, check_parts[0]);
  mpz_nextprime(check_parts[1], q);
  mpz_mul(check_input, check_parts[0], check_parts[1]);
  CHECK(mpz_sizeinbase(check_input, 2) == bits);
  CHECK(mpz_cmp_ui(check_parts[0], 64007) > 0);
  CHECK(mpz_cmp_ui(check_parts[1], 64007) > 0);
  mpz_clear(q);
}

/* Check the intended non-SIQS sequence independently of the policy tables.
 * A SIQS miss must leave the complete five-stage fallback intact. */
static void check_expected(size_t *index, int kind, UV a, UV b) {
  while (*index < check_nevents && check_events[*index].kind == CHECK_SIQS) (*index)++;
  CHECK(*index < check_nevents);
  CHECK(check_events[*index].kind == kind);
  CHECK(check_events[*index].a == a);
  CHECK(check_events[*index].b == b);
  (*index)++;
}
static void check_expected_ecm(size_t *index, UV b1, UV curves, UV b2) {
  check_expected(index, CHECK_ECM, b1, curves);
  CHECK(check_events[*index-1].b2 == b2);
}
static void check_fallback_recipe(unsigned bits) {
  size_t index = 0;
  UV base, early, m, pretest_b2_ratio = bits <= 224 ? 40 : 0;
  if (bits <= 64) check_expected(&index, CHECK_TINYECM64, 0, 0);
  early = bits < 88 ? 0 : bits < 112 ? 1000 : bits < 176 ? 2000 : 20000;
  if (early) check_expected(&index, CHECK_PM1, early, (bits <= 224 ? 5 : 10)*early);
  if (bits >= 128) check_expected_ecm(&index, 200, 4, 200*pretest_b2_ratio);
  if (bits >= 160) check_expected_ecm(&index, 600, 4, 600*pretest_b2_ratio);
  if (bits >= 176) check_expected_ecm(&index, 600, 12, 600*pretest_b2_ratio);
  if (bits >= 192) check_expected_ecm(&index, 2000, 10, 2000*pretest_b2_ratio);
  if (bits >= 208) {
    check_expected(&index, CHECK_PM1, 200000, 3000000);
    check_expected_ecm(&index, 20000, 2, 20000*pretest_b2_ratio);
  }
  if (bits >= 240) check_expected_ecm(&index, 30000, 20, 0);
  if (bits >= 256) check_expected_ecm(&index, 40000, 40, 0);
  if (bits >= 270) check_expected_ecm(&index, 80000, 40, 0);
  if (bits >= 288) check_expected_ecm(&index, 160000, 80, 0);
  base = bits < 65 ? 5000 : bits < 160 ? 20000 : bits < 192 ? 30000
       : bits < 224 ? 40000 : bits < 256 ? 80000 : bits < 512 ? 160000 : 320000;
  CHECK(factor_ecm_base_bound(bits) == base);
  check_expected_ecm(&index, 2*base, 20, 0);
  check_expected(&index, CHECK_PBRENT, 1, 1048576);
  check_expected_ecm(&index, 4*base, 20, 0);
  check_expected_ecm(&index, 8*base, 20, 0);
  check_expected(&index, CHECK_HOLF, 1048576, 0);
  check_expected(&index, CHECK_CHEB, 50000, 0);
  check_expected(&index, CHECK_PM1, 5000000, 100000000);
  check_expected_ecm(&index, 32*base, 40, 0);
  for (m = 8; m <= 128; m *= 2) check_expected_ecm(&index, m*base, 100, 0);
  while (index < check_nevents && check_events[index].kind == CHECK_SIQS) index++;
  CHECK(index == check_nevents);
}

static uint64_t check_ecm_b2_32(uint64_t b1, uint64_t b2) {
  return ecm_resolve_b2_32((uint32_t)b1, (uint32_t)b2);
}

static void check_ecm_bounds(uint64_t (*resolve)(uint64_t, uint64_t),
                             uint64_t maximum, uint64_t limit) {
  uint64_t b1, want, d, last;
  unsigned i;
  CHECK(resolve(0, 0) == 0);
  CHECK(resolve(7, 6) == 700);
  CHECK(resolve(7, 7) == 7); /* Explicit stage-1-only request. */
  CHECK(resolve(1000, 100000) == 100000);
  CHECK(resolve(40960000, 0) == UINT64_C(4096000000));
  CHECK(resolve(maximum, 0) == limit);
  CHECK(resolve(maximum, maximum) == limit);
  CHECK(resolve(1, maximum) == limit);
  CHECK(resolve(1, limit - 1) == limit - 1);
  CHECK(resolve(1, limit) == limit);
  CHECK(resolve(1, limit + 1) == limit);
  for (i = 0; i <= 1000; i++) {
    b1 = limit/100 - 500 + i;
    want = b1 > limit/100 ? limit : 100*b1;
    CHECK(resolve(b1, 0) == want);
    CHECK(resolve(b1, b1 - 1) == want);
  }
  d = (uint64_t)sqrt((double)limit / 2.0);
  d += d & 1;
  CHECK(limit <= maximum - 3*d);
  last = 1 + ((limit + d - 2) / (2*d)) * (2*d);
  CHECK(last <= maximum - 2*d);
  CHECK(last + 2*d >= limit + d); /* Final increment exits, without wrap. */
}

static void check_names(void) {
  factor_step_t step = {0, FACTOR_METHOD_LARGE_ECM, 0, 0, NULL};
  char first[64], second[64];
  const char *name;
  name = factor_step_name(&factor_pretests[0], 20000, 200000,
                          first, sizeof(first));
  CHECK(name == first && strcmp(name, "p-1 (20k)") == 0);
  name = factor_step_name(&factor_ladder[1], 320000, 20,
                          second, sizeof(second));
  CHECK(name == second && strcmp(name, "ecm (320k,20)") == 0);
  CHECK(strcmp(first, "p-1 (20k)") == 0); /* Caller buffers are independent. */
  name = factor_step_name(&step, 1310720000, 100, first, sizeof(first));
  CHECK(strcmp(name, "ecm (1310720k,100)") == 0);
  name = factor_step_name(&factor_ladder[0], 0, 0, first, sizeof(first));
  CHECK(name == factor_ladder[0].description && strcmp(name, "SIQS") == 0);
  name = factor_step_name(&factor_pretests[1], 200, 4, first, sizeof(first));
  CHECK(name == factor_pretests[1].description && strcmp(name, "ecm (200)") == 0);
}

static void check_schedule(void) {
  static const unsigned sizes[] = {36,52,53,55,56,59,63,64,65,87,88,89,95,96,111,112,127,128,129,
    143,144,159,160,161,175,176,177,179,180,181,191,192,193,199,200,201,207,208,209,223,224,225,239,240,255,256,269,270,
    287,288,319,320,321,359,360,361,399,400,401,511,512,513,600};
  size_t i, j;
  for (i = 0; i < sizeof(sizes)/sizeof(sizes[0]); i++) {
    unsigned bits = sizes[i];
    UV base = bits < 512 ? 160000 : 320000;
    check_reset();
    check_semiprime(bits);
    check_terminal_failure();
    CHECK(check_siqs_calls == (bits >= MPU_SIQS_MIN_BITS && bits <= MPU_SIQS_MAX_BITS));
    CHECK(check_releases == check_siqs_calls);
    check_fallback_recipe(bits);
    if (bits > 320 && bits <= MPU_SIQS_MAX_BITS) {
      for (j = 0; j < check_nevents && check_events[j].kind != CHECK_SIQS; j++) ;
      CHECK(j > 0 && j < check_nevents);
      CHECK(check_events[j-1].kind == CHECK_ECM);
      CHECK(check_events[j-1].a == (bits <= 360 ? 2*base : bits <= 400 ? 8*base : 32*base));
      CHECK(check_events[j-1].b == (bits <= 400 ? 20 : 100));
    }
  }
}

static void check_partition_value(factor_partition_t *part, const mpz_t expected)
{
  mpz_t product, power, g;
  size_t i, j;
  mpz_inits(product, power, g, NULL);
  mpz_set_ui(product, 1);
  for (i = 0; i < part->count; i++) {
    CHECK(part->pieces[i].exponent > 0 && mpz_cmp_ui(part->pieces[i].value, 1) > 0);
    mpz_pow_ui(power, part->pieces[i].value, part->pieces[i].exponent);
    mpz_mul(product, product, power);
    for (j = 0; j < i; j++) {
      mpz_gcd(g, part->pieces[i].value, part->pieces[j].value);
      CHECK(mpz_cmp_ui(g, 1) == 0);
    }
  }
  CHECK(mpz_cmp(product, expected) == 0);
  mpz_clears(product, power, g, NULL);
}

static void check_partition(void)
{
  factor_partition_t part;
  mpz_t n, d, root, product, power;
  size_t index, i;
  unsigned seed = 7, repeat, j;
  int result;
  mpz_inits(n, d, root, product, power, NULL);

  factor_partition_init(&part);
  mpz_set_ui(n, 12);
  CHECK(factor_partition_append(&part, n, 1, FACTOR_COMPOSITE));
  mpz_set_ui(d, 6);
  CHECK(factor_partition_split(&part, 0, d));
  check_partition_value(&part, n);
  CHECK(part.count == 2);
  CHECK(factor_partition_prepare(&part, root, &index, 0) && index == part.count);
  factor_partition_clear(&part);

  /* Overlapping composite pieces must become a coprime basis. */
  factor_partition_init(&part);
  mpz_set_ui(n, 12); mpz_set_ui(d, 18);
  CHECK(factor_partition_append(&part, n, 2, FACTOR_UNKNOWN));
  CHECK(factor_partition_append(&part, d, 3, FACTOR_UNKNOWN));
  mpz_pow_ui(product, n, 2); mpz_pow_ui(power, d, 3); mpz_mul(product, product, power);
  CHECK(factor_partition_refine(&part));
  check_partition_value(&part, product);
  CHECK(part.count == 2);
  CHECK(factor_partition_prepare(&part, root, &index, 0) && index == part.count);
  factor_partition_clear(&part);

  /* A proven lower factor bound can settle prime pieces without PRP.
   * Equality is not enough: 7^2 must still be recognized as a power. */
  factor_partition_init(&part);
  mpz_set_ui(n, 11);
  CHECK(factor_partition_append(&part, n, 1, FACTOR_UNKNOWN));
  mpz_set_ui(n, 47);
  CHECK(factor_partition_append(&part, n, 1, FACTOR_UNKNOWN));
  mpz_set_ui(n, 169);
  CHECK(factor_partition_append(&part, n, 1, FACTOR_UNKNOWN));
  mpz_set_ui(product, 11UL*47UL*169UL);
  check_primality_calls = 0;
  CHECK(factor_partition_prepare(&part, root, &index, 7UL*7UL));
  CHECK(index == part.count && check_primality_calls == 0);
  CHECK(mpz_cmp_ui(part.pieces[2].value, 13) == 0 && part.pieces[2].exponent == 2);
  check_partition_value(&part, product);
  factor_partition_clear(&part);

  factor_partition_init(&part);
  mpz_set_ui(n, 49);
  CHECK(factor_partition_append(&part, n, 1, FACTOR_UNKNOWN));
  CHECK(factor_partition_prepare(&part, root, &index, 49));
  CHECK(part.count == 1 && mpz_cmp_ui(part.pieces[0].value, 7) == 0);
  CHECK(part.pieces[0].exponent == 2 && part.pieces[0].state == FACTOR_PRIME);
  check_partition_value(&part, n);
  factor_partition_clear(&part);

  /* Multiplicity avoids a power-sized work array or repeated root work. */
  factor_partition_init(&part);
  mpz_set_ui(n, 2); mpz_pow_ui(n, n, 100000);
  CHECK(factor_partition_append(&part, n, 1, FACTOR_UNKNOWN));
  CHECK(factor_partition_prepare(&part, root, &index, 0));
  CHECK(index == part.count && part.count == 1);
  CHECK(mpz_cmp_ui(part.pieces[0].value, 2) == 0 && part.pieces[0].exponent == 100000);
  check_partition_value(&part, n);
  factor_partition_clear(&part);

  /* Removing all copies through gcd refinement must also be bounded. */
  factor_partition_init(&part);
  CHECK(factor_partition_append(&part, n, 1, FACTOR_COMPOSITE));
  mpz_set_ui(d, 2);
  CHECK(factor_partition_split(&part, 0, d));
  CHECK(part.count == 1 && part.pieces[0].exponent == 100000);
  check_partition_value(&part, n);
  factor_partition_clear(&part);

  factor_partition_init(&part);
  mpz_set_ui(n, 6);
  CHECK(factor_partition_append(&part, n, INT_MAX, FACTOR_COMPOSITE));
  CHECK(factor_partition_append(&part, n, 1, FACTOR_UNKNOWN));
  CHECK(!factor_partition_refine(&part));
  CHECK(strcmp(part.error, "factor() exponent overflow") == 0);
  factor_partition_clear(&part);

  /* Each overlap path must check weighted exponents before multiplying. */
  for (j = 0; j < 3; j++) {
    factor_partition_init(&part);
    mpz_set_ui(n, j == 0 ? 6 : j == 1 ? 2 : 8);
    CHECK(factor_partition_append(&part, n, j == 1 ? 1 : INT_MAX, FACTOR_COMPOSITE));
    mpz_set_ui(n, j == 0 ? 10 : j == 1 ? 8 : 2);
    CHECK(factor_partition_append(&part, n, j == 1 ? INT_MAX : 1, FACTOR_COMPOSITE));
    CHECK(!factor_partition_refine(&part));
    CHECK(strcmp(part.error, "factor() exponent overflow") == 0);
    factor_partition_clear(&part);
  }
  factor_partition_init(&part);
  mpz_set_ui(n, 4);
  CHECK(factor_partition_append(&part, n, INT_MAX, FACTOR_UNKNOWN));
  CHECK(!factor_partition_prepare(&part, root, &index, 0));
  CHECK(strcmp(part.error, "factor() exponent overflow") == 0);
  factor_partition_clear(&part);

  factor_partition_init(&part);
  mpz_set_ui(n, 6);
  CHECK(factor_partition_append(&part, n, 1, FACTOR_COMPOSITE));
  for (j = 0; j < 4; j++) {
    mpz_set_ui(d, j == 3 ? 4 : j == 2 ? 6 : j);
    CHECK(!factor_partition_split(&part, 0, d));
    CHECK(mpz_cmp_ui(part.pieces[0].value, 6) == 0 && part.count == 1);
  }
  factor_partition_clear(&part);

  /* Every intermediate refinement preserves the independent product. */
  for (repeat = 0; repeat < 200; repeat++) {
    factor_partition_init(&part);
    mpz_set_ui(product, 1);
    for (j = 0; j < 6; j++) {
      unsigned value, exponent;
      seed = seed*1664525U + 1013904223U; value = 2 + seed%499;
      seed = seed*1664525U + 1013904223U; exponent = 1 + seed%7;
      mpz_set_ui(n, value);
      CHECK(factor_partition_append(&part, n, (int)exponent, FACTOR_UNKNOWN));
      mpz_pow_ui(power, n, exponent); mpz_mul(product, product, power);
    }
    CHECK(factor_partition_refine(&part));
    check_partition_value(&part, product);
    while (1) {
      unsigned long divisor;
      CHECK(factor_partition_prepare(&part, root, &index, 0));
      check_partition_value(&part, product);
      if (index == part.count) break;
      for (divisor = 2; !mpz_divisible_ui_p(part.pieces[index].value, divisor); divisor++) ;
      mpz_set_ui(d, divisor);
      CHECK(factor_partition_split(&part, index, d));
      check_partition_value(&part, product);
    }
    for (i = 0; i < part.count; i++) CHECK(part.pieces[i].state == FACTOR_PRIME);
    factor_partition_clear(&part);
  }

  /* Low-level SIQS misses and invalid partitions must release their arrays. */
  check_reset(); check_semiprime(100);
  factor_partition_init(&part);
  CHECK(factor_partition_append(&part, check_input, 1, FACTOR_COMPOSITE));
  CHECK(factor_partition_siqs(&part, 0, check_input, 0) == 0);
  check_partition_value(&part, check_input);
  check_split = 1; check_nparts = 0;
  CHECK(factor_partition_siqs(&part, 0, check_input, 0) == 0);
  check_nparts = 2; mpz_add_ui(check_parts[0], check_parts[0], 1);
  result = factor_partition_siqs(&part, 0, check_input, 0);
  CHECK(result == -1 && strcmp(part.error, "SIQS: invalid result partition") == 0);
  check_partition_value(&part, check_input);
  CHECK(check_releases == check_siqs_calls);
  factor_partition_clear(&part);
  mpz_clears(n, d, root, product, power, NULL);
}

static void check_output(const mpz_t n, int count, mpz_t *factors, int *exponents)
{
  mpz_t product, power;
  int i;
  mpz_inits(product, power, NULL);
  mpz_set_ui(product, 1);
  for (i = 0; i < count; i++) {
    CHECK(_GMP_is_prob_prime(factors[i]) && exponents[i] > 0);
    if (i) CHECK(mpz_cmp(factors[i-1], factors[i]) < 0);
    mpz_pow_ui(power, factors[i], exponents[i]);
    mpz_mul(product, product, power);
  }
  CHECK(mpz_cmp(product, n) == 0);
  mpz_clears(product, power, NULL);
}

static void check_real_loop(void)
{
  mpz_t *factors;
  mpz_t product;
  int *exponents, nfactors, i;
  check_reset(); check_semiprime(100); check_split = 1;
  nfactors = factor(check_input, &factors, &exponents);
  CHECK(nfactors == 2 && check_siqs_calls == 1);
  check_output(check_input, nfactors, factors, exponents);
  clear_factors(nfactors, &factors, &exponents);

  /* SIQS returns p^2 and q. Its composite power needs no second SIQS. */
  check_reset(); check_semiprime(100);
  mpz_set(check_parts[3], check_parts[0]);
  mpz_mul(check_parts[0], check_parts[0], check_parts[0]);
  mpz_mul(check_input, check_parts[0], check_parts[1]);
  check_split = 1;
  nfactors = factor(check_input, &factors, &exponents);
  CHECK(nfactors == 2 && check_siqs_calls == 1);
  check_output(check_input, nfactors, factors, exponents);
  clear_factors(nfactors, &factors, &exponents);

  /* pq, qr, pr, s: one SIQS partition resolves all four distinct primes. */
  check_reset(); check_split = 1; check_nparts = 4;
  for (i = 3; i < 7; i++) {
    mpz_set_ui(check_parts[i], 1000000 + 100*i);
    mpz_nextprime(check_parts[i], check_parts[i]);
  }
  mpz_mul(check_parts[0], check_parts[4], check_parts[5]);
  mpz_mul(check_parts[1], check_parts[5], check_parts[6]);
  mpz_mul(check_parts[2], check_parts[4], check_parts[6]);
  mpz_set_ui(check_input, 1);
  for (i = 0; i < 4; i++) mpz_mul(check_input, check_input, check_parts[i]);
  nfactors = factor(check_input, &factors, &exponents);
  CHECK(nfactors == 4 && check_siqs_calls == 1);
  check_output(check_input, nfactors, factors, exponents);
  clear_factors(nfactors, &factors, &exponents);

  /* Exceed the previous fixed stack limit with real small prime factors. */
  check_reset(); mpz_init_set_ui(product, 1); mpz_set_ui(check_parts[0], 1);
  for (i = 0; i < 200; i++) {
    mpz_nextprime(check_parts[0], check_parts[0]);
    mpz_mul(product, product, check_parts[0]);
  }
  nfactors = factor(product, &factors, &exponents);
  CHECK(nfactors == 200 && check_siqs_calls == 0);
  check_output(product, nfactors, factors, exponents);
  clear_factors(nfactors, &factors, &exponents);
  mpz_clear(product);
}

static void check_pm1_boundaries(void)
{
  int (*methods[2])(const mpz_t, mpz_t, UV, UV) =
    { check_pm1_list, check_pm1_iterator };
  mpz_t exponent, a, t, n, pz, qz, f;
  unsigned long bound, q, p, k;
  unsigned method;
  mpz_inits(exponent, a, t, n, pz, qz, f, NULL);

  /* Independent lcm/powm arithmetic supplies a factor requiring precisely
   * the first prime after B1. Exercise prime and nonprime bounds. */
  for (bound = 7; bound <= 100; bound++) {
    mpz_set_ui(exponent, 1);
    for (k = 2; k <= bound; k++) mpz_lcm_ui(exponent, exponent, k);
    mpz_set_ui(qz, bound); mpz_nextprime(qz, qz); q = mpz_get_ui(qz);
    for (k = 2; k <= 8*bound; k += 2) {
      p = k*q+1;
      mpz_set_ui(pz, p);
      if (!mpz_probab_prime_p(pz, 20)) continue;
      mpz_mul_ui(n, pz, 1000000007);
      mpz_set_ui(a, 2); mpz_powm(a, a, exponent, n);
      mpz_sub_ui(t, a, 1); mpz_gcd(f, t, n);
      if (mpz_cmp_ui(f, 1) != 0) continue;
      mpz_powm_ui(t, a, q, n); mpz_sub_ui(t, t, 1); mpz_gcd(f, t, n);
      if (mpz_cmp_ui(f, p) == 0) break;
    }
    CHECK(k <= 8*bound);
    for (method = 0; method < 2; method++) {
      CHECK(!methods[method](n, f, bound, q-1));
      CHECK(mpz_cmp(f, n) == 0);
      CHECK(methods[method](n, f, bound, q));
      CHECK(mpz_cmp_ui(f, p) == 0);
    }
  }
  /* A later prime is included at B2, but never when it lies beyond B2. */
  mpz_set_ui(n, 84463); mpz_mul_ui(n, n, 1000000007);
  for (method = 0; method < 2; method++) {
    CHECK(!methods[method](n, f, 2000, 2010));
    CHECK(mpz_cmp(f, n) == 0);
    CHECK(methods[method](n, f, 2000, 2011));
    CHECK(mpz_cmp_ui(f, 84463) == 0);
  }
  /* A whole-input GCD must terminate the non-backtracking second stage. */
  mpz_set_ui(n, 72109); mpz_mul_ui(n, n, 84127);
  for (method = 0; method < 2; method++) {
    check_pm1_stage2_terms = 0;
    CHECK(!methods[method](n, f, 2000, 1000000));
    CHECK(mpz_cmp(f, n) == 0);
    CHECK(check_pm1_stage2_terms == 34);
  }
  /* 4547-1 = 2*2273; prime 2273 is stage-2 term 35. The second GCD
   * must be at term 162, retaining the early first check at term 34. */
  mpz_set_ui(n, 4547); mpz_mul_ui(n, n, 1000000007);
  for (method = 0; method < 2; method++) {
    check_pm1_stage2_terms = 0;
    CHECK(methods[method](n, f, 2000, 20000));
    CHECK(mpz_cmp_ui(f, 4547) == 0);
    CHECK(check_pm1_stage2_terms == 162);
  }
#if BITS_PER_WORD == 64
  /* 2^32 is an actual stage-1 prime power, not just an oversized bound. */
  mpz_set_ui(n, 10403);
  CHECK(check_pm1_ui32(n, f, UINT64_C(4294967296), UINT64_C(4294967296)));
  CHECK(mpz_cmp_ui(f, 101) == 0);
#endif
  CHECK(check_pm1_b2_32(2000) == 20000);
  CHECK(check_pm1_b2_32(UINT32_MAX/10) == UINT32_MAX/10*10);
  CHECK(check_pm1_b2_32(UINT32_MAX/10+1) == UINT32_MAX);
  CHECK(check_pm1_b2_32(UINT32_MAX) == UINT32_MAX);
  CHECK(check_pm1_b2_64(2000) == 20000);
  CHECK(check_pm1_b2_64(UINT64_MAX/10) == UINT64_MAX/10*10);
  CHECK(check_pm1_b2_64(UINT64_MAX/10+1) == UINT64_MAX);
  CHECK(check_pm1_b2_64(UINT64_MAX) == UINT64_MAX);
  mpz_clears(exponent, a, t, n, pz, qz, f, NULL);
}

static void check_early_pm1_until_miss(void)
{
  mpz_t *factors;
  int *exponents, count;
  unsigned hit, i, early, large;

  for (hit = 0; hit <= 4; hit++) {
    check_reset(); check_split = 1; check_pm1_split = hit;
    for (i = 0; i < 8; i++) mpz_set_ui(check_parts[i], 1);
    mpz_set_ui(check_parts[3], 1);
    mpz_mul_2exp(check_parts[3], check_parts[3], 63);
    mpz_nextprime(check_parts[3], check_parts[3]);
    for (i = 4; i < 8; i++) mpz_nextprime(check_parts[i], check_parts[i-1]);
    mpz_set(check_parts[0], check_parts[3]);
    for (i = 4; i < 8; i++) mpz_mul(check_parts[1], check_parts[1], check_parts[i]);
    mpz_mul(check_input, check_parts[0], check_parts[1]);
    count = factor(check_input, &factors, &exponents);
    CHECK(count == 5 && check_siqs_calls == 4-hit);
    check_output(check_input, count, factors, exponents);
    early = large = 0;
    for (i = 0; i < check_nevents; i++) {
      if (check_events[i].kind != CHECK_PM1) continue;
      if (check_events[i].a == 20000 || check_events[i].a == 2000) early++;
      else {
        CHECK(check_events[i].a == 200000 && check_events[i].b == 3000000);
        large++;
      }
    }
    CHECK(early == (hit < 4 ? hit+1 : 4));
    CHECK(large == (hit < 2 ? 2-hit : 0));
    clear_factors(count, &factors, &exponents);
  }

  /* A prime and a prime power finish classification before early p-1. */
  for (i = 1; i <= 3; i += 2) {
    check_reset();
    mpz_pow_ui(check_input, check_parts[3], i);
    count = factor(check_input, &factors, &exponents);
    CHECK(count == 1 && exponents[0] == (int)i && check_nevents == 0);
    check_output(check_input, count, factors, exponents);
    clear_factors(count, &factors, &exponents);
  }
}

static void check_failure_cleanup(void)
{
  mpz_t original;
  int i;
  check_reset();
  for (i = 3; i < 8; i++) {
    mpz_set_ui(check_parts[i], 1000000 + 100*i);
    mpz_nextprime(check_parts[i], check_parts[i]);
  }
  mpz_set(check_parts[0], check_parts[7]);
  mpz_mul(check_parts[1], check_parts[3], check_parts[4]);
  mpz_mul(check_parts[2], check_parts[5], check_parts[6]);
  mpz_set_ui(check_input, 1);
  for (i = 0; i < 3; i++) mpz_mul(check_input, check_input, check_parts[i]);
  mpz_init_set(original, check_input);
  check_split = 1; check_nparts = 3; check_fail_after = 1;
  check_terminal_failure();
  CHECK(check_siqs_calls == 2 && check_releases == 2);
  CHECK(mpz_cmp(check_input, original) == 0);
  mpz_clear(original);
}


int main(void) {
  unsigned i;
  mp_set_memory_functions(check_gmp_malloc, check_gmp_realloc, check_gmp_free);
  mpz_init(check_input);
  for (i = 0; i < 8; i++) mpz_init(check_parts[i]);
  prime_iterator_global_startup();
  _init_factor();
  check_ecm_bounds(check_ecm_b2_32, UINT32_MAX, UINT64_C(4294770687));
  check_ecm_bounds(ecm_resolve_b2_64, UINT64_MAX,
                   UINT64_MAX - UINT64_C(12884901888));
  check_names();
#ifdef FACTOR_CHECK_TRACE
  /* Optional verbose scheduling trace; no real high-bit factoring is done. */
  check_verbose = 2;
#endif
  check_schedule();
  check_verbose = 0;
  check_pm1_boundaries();
  check_partition();
  check_early_pm1_until_miss();
  check_failure_cleanup();
  check_real_loop();
  _destroy_factor();
  prime_iterator_global_shutdown();
  for (i = 0; i < 8; i++) mpz_clear(check_parts[i]);
  mpz_clear(check_input);
  CHECK(check_allocations == 0);
  printf("PASS factor partition/ladder: %u checks; p-1 boundaries, "
         "scheduling through 600 bits, 32/64-bit ECM bounds, failure cleanup\n",
         check_assertions);
  return 0;
}
