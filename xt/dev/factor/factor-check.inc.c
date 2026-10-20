/* Included after the unchanged factor.c prefix by factor-check.pl.
 * All factoring methods are scripted: no high-bit input is actually split.
 * GMP primality, factor bookkeeping and the prime iterator remain real. */
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

enum { CHECK_ECM, CHECK_PM1, CHECK_PBRENT, CHECK_HOLF, CHECK_CHEB,
       CHECK_SIQS, CHECK_NATIVE, CHECK_SQUFOF };
typedef struct { int kind; UV a, b; } check_event_t;
static check_event_t check_events[256];
static size_t check_nevents;
static unsigned check_siqs_calls, check_releases, check_assertions, check_fail_after;
static mpz_t check_input, check_parts[8];
static unsigned check_nparts;
static int check_split, check_verbose;
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
  check_events[check_nevents++].b = b;
}

int get_verbose_level(void) { return check_verbose; }
int _GMP_is_prob_prime(const mpz_t n) { return mpz_probab_prime_p(n, 20) != 0; }
unsigned long power_factor(const mpz_t n, mpz_t f) { (void)n; (void)f; return 0; }
void mpz_set_uv(mpz_t n, UV v) { mpz_import(n, 1, 1, sizeof(v), 0, 0, &v); }
void mpz_product_ui(mpz_t n, unsigned long *v, unsigned long count) {
  unsigned long i;
  mpz_set_ui(n, 1);
  for (i = 0; i < count; i++) mpz_mul_ui(n, n, v[i]);
}
int pbrent63(const mpz_t n, mpz_t f, UV rounds) {
  (void)n; (void)f; check_record(CHECK_NATIVE, rounds, 0); return 0;
}
int squfof126(const mpz_t n, mpz_t f, UV rounds) {
  (void)n; (void)f; check_record(CHECK_SQUFOF, rounds, 0); return 0;
}
int _GMP_pminus1_factor(const mpz_t n, mpz_t f, UV b1, UV b2) {
  (void)n; (void)f; check_record(CHECK_PM1, b1, b2); return 0;
}
int _GMP_ecm_factor_projective(const mpz_t n, mpz_t f, UV b1, UV b2, UV curves) {
  (void)n; (void)f; CHECK(b2 == 0); check_record(CHECK_ECM, b1, curves); return 0;
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

/* Compare all non-SIQS calls with the old recipe, independently of the new
 * tables.  A SIQS miss must leave the complete fallback sequence intact. */
static void check_expected(size_t *index, int kind, UV a, UV b) {
  while (*index < check_nevents && check_events[*index].kind == CHECK_SIQS) (*index)++;
  CHECK(*index < check_nevents);
  CHECK(check_events[*index].kind == kind);
  CHECK(check_events[*index].a == a);
  CHECK(check_events[*index].b == b);
  (*index)++;
}
static void check_old_recipe(unsigned bits) {
  size_t index = 0;
  UV base, early, m;
  if (bits <= 63) check_expected(&index, CHECK_NATIVE, 400000, 0);
  if (bits >= 53 && bits <= 59) check_expected(&index, CHECK_PM1, 1000, 10000);
  if (bits < 56) check_expected(&index, CHECK_SQUFOF, bits < 53 ? 400000 : 1000000, 0);
  early = bits < 96 ? 0 : bits < 112 ? 1000 : bits < 144 ? 2000 : 20000;
  if (early) check_expected(&index, CHECK_PM1, early, 10*early);
  if (bits >= 160) check_expected(&index, CHECK_ECM, 200, 4);
  if (bits >= 200) {
    check_expected(&index, CHECK_ECM, 600, 20);
    check_expected(&index, CHECK_ECM, 2000, 10);
  }
  if (bits >= 208) {
    check_expected(&index, CHECK_PM1, 200000, 3000000);
    check_expected(&index, CHECK_ECM, 20000, 2);
  }
  if (bits >= 240) check_expected(&index, CHECK_ECM, 30000, 20);
  if (bits >= 256) check_expected(&index, CHECK_ECM, 40000, 40);
  if (bits >= 270) check_expected(&index, CHECK_ECM, 80000, 40);
  if (bits >= 288) check_expected(&index, CHECK_ECM, 160000, 80);
  base = bits < 65 ? 5000 : bits < 160 ? 20000 : bits < 192 ? 30000
       : bits < 224 ? 40000 : bits < 256 ? 80000 : bits < 512 ? 160000 : 320000;
  check_expected(&index, CHECK_ECM, 2*base, 20);
  check_expected(&index, CHECK_PBRENT, 1, 1048576);
  check_expected(&index, CHECK_ECM, 4*base, 20);
  check_expected(&index, CHECK_ECM, 8*base, 20);
  check_expected(&index, CHECK_HOLF, 1048576, 0);
  check_expected(&index, CHECK_CHEB, 50000, 0);
  check_expected(&index, CHECK_PM1, 5000000, 100000000);
  check_expected(&index, CHECK_ECM, 32*base, 40);
  for (m = 8; m <= 4096; m *= 2) check_expected(&index, CHECK_ECM, m*base, 100);
  while (index < check_nevents && check_events[index].kind == CHECK_SIQS) index++;
  CHECK(index == check_nevents);
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
  static const unsigned sizes[] = {36,52,53,55,56,59,63,64,65,95,96,111,112,
    143,144,159,160,191,192,199,200,207,208,223,224,239,240,255,256,269,270,
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
    check_old_recipe(bits);
    if (bits > 320 && bits <= MPU_SIQS_MAX_BITS) {
      for (j = 0; j < check_nevents && check_events[j].kind != CHECK_SIQS; j++) ;
      CHECK(j > 0 && j < check_nevents);
      CHECK(check_events[j-1].kind == CHECK_ECM);
      CHECK(check_events[j-1].a == (bits <= 360 ? 2*base : bits <= 400 ? 8*base : 32*base));
      CHECK(check_events[j-1].b == (bits <= 400 ? 20 : 100));
    }
  }
}

static void check_helper(void) {
  mpz_t n, f, product, stack[MAX_FACTORS];
  mpz_t *factors = NULL;
  int *exponents = NULL, nfactors = 0, nstack = 0;
  unsigned i;
  factor_siqs_result_t result;
  mpz_inits(n, f, product, NULL);
  check_reset();
  check_semiprime(100);
  mpz_set(n, check_input);
  mpz_set_ui(f, 17);
  result = factor_try_siqs(n, f, stack, &nstack, &nfactors, &factors, &exponents, 0);
  CHECK(result == FACTOR_SIQS_NO_SPLIT && nstack == 0 && nfactors == 0);
  CHECK(mpz_cmp(n, check_input) == 0 && mpz_cmp_ui(f, 17) == 0);
  check_split = 1;
  check_nparts = 0;
  result = factor_try_siqs(n, f, stack, &nstack, &nfactors, &factors, &exponents, 0);
  CHECK(result == FACTOR_SIQS_NO_SPLIT && mpz_cmp(n, check_input) == 0);
  check_nparts = 2;
  result = factor_try_siqs(n, f, stack, &nstack, &nfactors, &factors, &exponents, 0);
  CHECK(result == FACTOR_SIQS_FACTOR && nstack == 0 && nfactors == 0);
  CHECK(mpz_cmp(n, check_input) == 0 && mpz_cmp(f, check_parts[0]) == 0);
  mpz_set(check_parts[1], check_parts[0]);
  mpz_mul(check_input, check_parts[0], check_parts[1]);
  mpz_set(n, check_input);
  check_verbose = 3;
  result = factor_try_siqs(n, f, stack, &nstack, &nfactors, &factors, &exponents, 3);
  check_verbose = 0;
  CHECK(result == FACTOR_SIQS_FACTOR && nstack == 0 && nfactors == 0);
  CHECK(mpz_cmp(n, check_input) == 0 && mpz_cmp(f, check_parts[0]) == 0);

  /* Repeated primes and one remaining composite: account for multiplicity. */
  check_nparts = 4;
  mpz_set_ui(check_parts[0], 1000003);
  mpz_set(check_parts[1], check_parts[0]);
  mpz_set_ui(check_parts[2], 1000033);
  mpz_set_ui(check_parts[3], 1000037);
  mpz_mul(check_parts[3], check_parts[3], check_parts[2]);
  mpz_set_ui(check_input, 1);
  for (i = 0; i < check_nparts; i++) mpz_mul(check_input, check_input, check_parts[i]);
  mpz_set(n, check_input);
  result = factor_try_siqs(n, f, stack, &nstack, &nfactors, &factors, &exponents, 0);
  CHECK(result == FACTOR_SIQS_WORK_UPDATED && nfactors == 2 && nstack == 0);
  CHECK(exponents[0] == 2 && exponents[1] == 1);
  CHECK(mpz_cmp(n, check_parts[3]) == 0);
  clear_factors(nfactors, &factors, &exponents);
  factors = NULL; exponents = NULL; nfactors = 0;

  /* Composite-only partitions preserve unpushed products under stack pressure. */
  check_nparts = 3;
  for (i = 0; i < check_nparts; i++) {
    mpz_set_ui(check_parts[i], 1000003 + 30*i);
    mpz_mul_ui(check_parts[i], check_parts[i], 1000037);
  }
  mpz_set_ui(check_input, 1);
  for (i = 0; i < check_nparts; i++) mpz_mul(check_input, check_input, check_parts[i]);
  mpz_set(n, check_input);
  nstack = MAX_FACTORS-2;
  result = factor_try_siqs(n, f, stack, &nstack, &nfactors, &factors, &exponents, 0);
  CHECK(result == FACTOR_SIQS_WORK_UPDATED && nstack == MAX_FACTORS-1 && nfactors == 0);
  mpz_mul(product, n, stack[MAX_FACTORS-2]);
  CHECK(mpz_cmp(product, check_input) == 0);
  mpz_clear(stack[MAX_FACTORS-2]);
  mpz_set(n, check_input);
  result = factor_try_siqs(n, f, stack, &nstack, &nfactors, &factors, &exponents, 0);
  CHECK(result == FACTOR_SIQS_NO_SPLIT && mpz_cmp(n, check_input) == 0);
  CHECK(check_releases == check_siqs_calls);
  mpz_clears(n, f, product, NULL);
}

static void check_real_loop(void) {
  mpz_t *factors;
  mpz_t product;
  int *exponents, nfactors, i;
  check_reset();
  check_semiprime(100);
  check_split = 1;
  nfactors = factor(check_input, &factors, &exponents);
  CHECK(nfactors == 2 && check_siqs_calls == 1);
  clear_factors(nfactors, &factors, &exponents);

  /* A composite returned in a pair is queued normally, then gets a fresh
   * SIQS attempt.  The second reply consists of the two equal prime factors. */
  check_reset();
  check_semiprime(100);
  mpz_set(check_parts[3], check_parts[0]);
  mpz_mul(check_parts[0], check_parts[0], check_parts[0]);
  mpz_mul(check_input, check_parts[0], check_parts[1]);
  check_nparts = 2;
  check_split = 1;
  /* Return composite + prime at the root, then prime + prime. */
  nfactors = factor(check_input, &factors, &exponents);
  CHECK(nfactors == 2 && check_siqs_calls == 2);
  mpz_init_set_ui(product, 1);
  for (i = 0; i < nfactors; i++) {
    CHECK(_GMP_is_prob_prime(factors[i]));
    while (exponents[i]-- > 0) mpz_mul(product, product, factors[i]);
  }
  CHECK(mpz_cmp(product, check_input) == 0);
  mpz_clear(product);
  clear_factors(nfactors, &factors, &exponents);
  CHECK(check_releases == check_siqs_calls);

  /* A three-composite partition takes the WORK_UPDATED path; all queued
   * and retained cofactors get a fresh attempt and conserve the full input. */
  check_reset();
  check_split = 1;
  check_nparts = 3;
  mpz_init(product);
  mpz_set_ui(check_input, 1);
  for (i = 0; i < 3; i++) {
    mpz_set_ui(check_parts[3+i], 1000000 + 100*i);
    mpz_nextprime(check_parts[3+i], check_parts[3+i]);
    mpz_nextprime(product, check_parts[3+i]);
    mpz_mul(check_parts[i], check_parts[3+i], product);
    mpz_mul(check_input, check_input, check_parts[i]);
  }
  nfactors = factor(check_input, &factors, &exponents);
  CHECK(nfactors == 6 && check_siqs_calls == 4);
  mpz_set_ui(product, 1);
  for (i = 0; i < nfactors; i++) {
    CHECK(_GMP_is_prob_prime(factors[i]) && exponents[i] == 1);
    mpz_mul(product, product, factors[i]);
  }
  CHECK(mpz_cmp(product, check_input) == 0);
  CHECK(check_releases == check_siqs_calls);
  clear_factors(nfactors, &factors, &exponents);
  mpz_clear(product);
}

static void check_failure_cleanup(void) {
  mpz_t original;
  check_reset();
  check_semiprime(100);
  mpz_mul(check_parts[2], check_parts[0], check_parts[1]);
  mpz_mul(check_parts[1], check_parts[1], check_parts[1]);
  mpz_set_ui(check_parts[0], 1000003);
  mpz_mul(check_input, check_parts[0], check_parts[1]);
  mpz_mul(check_input, check_input, check_parts[2]);
  mpz_init_set(original, check_input);
  check_split = 1;
  check_nparts = 3;
  check_fail_after = 1;
  /* Accept one prime and queue a composite, then exhaust the retained
   * composite.  The croak must release partial results and pending work. */
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
  _init_factor();
  check_names();
#ifdef FACTOR_CHECK_TRACE
  /* Optional verbose scheduling trace; no real high-bit factoring is done. */
  check_verbose = 2;
#endif
  check_schedule();
  check_verbose = 0;
  check_helper();
  check_failure_cleanup();
  check_real_loop();
  _destroy_factor();
  for (i = 0; i < 8; i++) mpz_clear(check_parts[i]);
  mpz_clear(check_input);
  CHECK(check_allocations == 0);
  printf("PASS factor ladder/helper: %u checks; boundaries through 600 bits, failure cleanup\n", check_assertions);
  return 0;
}
