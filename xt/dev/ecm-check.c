/* Private ECM arithmetic/allocation regressions, not shipped with the package.
 * Run with: perl xt/dev/ecm-check.pl [--sanitize] [--threads]
 * Includes the production implementation to exercise its internal chains.
 * The independent oracle uses ordinary affine arithmetic over a prime field.
 * Copyright (c) 2026 Dana Jacobsen. See LICENSE for redistribution terms. */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <setjmp.h>
#include <string.h>
#include <gmp.h>
#include "../../ptypes.h"
#ifdef ECM_CHECK_THREADS
#include <pthread.h>
#endif

#ifdef ECM_CHECK_AFFINE_ONLY
#include "verifier-affine.h"
#else
/* Trap before any huge test allocation.  Croaks and reaching Newz are
 * distinguishable, so the overflow regression cannot accidentally allocate. */
static jmp_buf allocation_trap;
static int allocation_probe;
static UV probed_count;
static unsigned allocation_calls;
static char error_text[128];
#ifdef ECM_CHECK_THREADS
static pthread_mutex_t allocation_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif

static void check_croak(const char *format, ...)
{
  va_list args;
  va_start(args, format);
  vsnprintf(error_text, sizeof(error_text), format, args);
  va_end(args);
  if (allocation_probe) longjmp(allocation_trap, 1);
  fprintf(stderr, "Unexpected croak: %s\n", error_text);
  abort();
}

static void *check_calloc(UV count, size_t size)
{
#ifdef ECM_CHECK_THREADS
  pthread_mutex_lock(&allocation_mutex);
#endif
  allocation_calls++;
#ifdef ECM_CHECK_THREADS
  pthread_mutex_unlock(&allocation_mutex);
#endif
  if (allocation_probe) {
    probed_count = count;
    longjmp(allocation_trap, 2);
  }
  return mpu_calloc((size_t)count, size);
}

#undef croak
#define croak check_croak
#undef Newz
#define Newz(id, mem, count, type) \
  ((mem) = (type*)check_calloc((count), sizeof(type)))
#include "../../ecm.c"

#endif

void mpz_set_uv(mpz_t n, UV value)
{ mpz_import(n, 1, 1, sizeof(value), 0, 0, &value); }

static unsigned long checks;
#define CHECK(c) do { \
  checks++; \
  if (!(c)) { \
    fprintf(stderr, "FAIL line %d: %s\n", __LINE__, #c); \
    exit(1); \
  } \
} while (0)

typedef struct { int64_t x, y; int infinity; } reference_point;

static int64_t mod(int64_t a, int64_t p)
{ a %= p; return a < 0 ? a+p : a; }

static int64_t power(int64_t a, int64_t e, int64_t p)
{
  int64_t r = 1;
  for (; e != 0; e >>= 1, a = a*a % p)
    if (e & 1) r = r*a % p;
  return r;
}

/* y^2 = x^3 + a2*x^2 + a4*x + a6.  All test moduli are <= 10^9+7,
 * so products and the doubling numerator fit in signed 64 bits. */
static reference_point add(reference_point P, reference_point Q,
                           int64_t a2, int64_t a4, int64_t p)
{
  reference_point R = {0, 0, 1};
  int64_t m, numerator, denominator;
  if (P.infinity) return Q;
  if (Q.infinity) return P;
  if (P.x == Q.x) {
    if (mod(P.y+Q.y, p) == 0) return R;
    CHECK(P.y == Q.y);
    numerator = mod(3*P.x*P.x + 2*a2*P.x + a4, p);
    denominator = mod(2*P.y, p);
  } else {
    numerator = mod(Q.y-P.y, p);
    denominator = mod(Q.x-P.x, p);
  }
  CHECK(denominator != 0);
  m = numerator * power(denominator, p-2, p) % p;
  R.x = mod(m*m - a2 - P.x - Q.x, p);
  R.y = mod(m*(P.x-R.x) - P.y, p);
  R.infinity = 0;
  return R;
}

static reference_point multiply(reference_point P, UV k,
                                int64_t a2, int64_t a4, int64_t p)
{
  reference_point R = {0, 0, 1};
  while (k != 0) {
    if (k & 1) R = add(R, P, a2, a4, p);
    k >>= 1;
    if (k != 0) P = add(P, P, a2, a4, p);
  }
  return R;
}

static void check_affine_point(reference_point P, unsigned long a,
                               unsigned long p, UV limit)
{
  struct ec_affine_point input, output;
  mpz_t A, N, k, d;
  UV scalar;
  mpz_inits(A, N, k, d, input.x, input.y, output.x, output.y, NULL);
  mpz_set_ui(A, a); mpz_set_ui(N, p);
  for (scalar = 0; scalar <= limit; scalar++) {
    reference_point expected = multiply(P, scalar, 0, a, p);
    int alias;
    for (alias = 0; alias <= 1; alias++) {
      struct ec_affine_point *R = alias ? &input : &output;
      /* Noncanonical input coordinates must also be handled. */
      mpz_set_si(input.x, P.x - (int64_t)p);
      mpz_set_si(input.y, P.y + (int64_t)p);
      input.infinity = P.infinity;
      mpz_set_uv(k, scalar);
      CHECK(!ec_affine_multiply(A, k, N, input, R, d));
      CHECK(mpz_cmp_ui(d, 1) == 0 && mpz_sgn(k) == 0);
      CHECK(R->infinity == expected.infinity);
      if (!expected.infinity) {
        CHECK(mpz_cmp_ui(R->x, (unsigned long)expected.x) == 0);
        CHECK(mpz_cmp_ui(R->y, (unsigned long)expected.y) == 0);
      }
    }
  }
  mpz_clears(A, N, k, d, input.x, input.y, output.x, output.y, NULL);
}

static void check_affine(void)
{
  const reference_point P = {3, 6, 0};    /* order 5, a=2, b=3, p=97 */
  const reference_point finite01 = {0, 1, 0}; /* a=2, b=1, p=97 */
  const reference_point torsion = {96, 0, 0}; /* order 2, a=2, b=3 */
  const reference_point infinity = {0, 1, 1};
  struct ec_affine_point input, other, output;
  mpz_t a, n, k, d, m, t1, t2;
  check_affine_point(P, 2, 97, 512);
  check_affine_point(finite01, 2, 97, 128);
  check_affine_point(torsion, 2, 97, 16);
  check_affine_point(infinity, 2, 97, 16);
  mpz_inits(a,n,k,d,m,t1,t2,input.x,input.y,other.x,other.y,
            output.x,output.y,NULL);
  mpz_set_ui(a, 2); mpz_set_ui(n, 35);
  input.infinity = other.infinity = 0;
  /* Failed doubling inversion: gcd(2*y,n)=5.  Do not report infinity. */
  mpz_set_ui(input.x, 3); mpz_set_ui(input.y, 5); mpz_set_ui(k, 2);
  mpz_set_ui(output.x, 42); mpz_set_ui(output.y, 43); output.infinity = 0;
  CHECK(ec_affine_multiply(a, k, n, input, &output, d));
  CHECK(mpz_cmp_ui(d, 5) == 0);
  CHECK(mpz_cmp_ui(output.x, 42) == 0 && !output.infinity);
  /* Equal x with CRT-mixed signs on y^2=x^3+2*x+1 modulo 35. */
  mpz_set_ui(input.x, 0); mpz_set_ui(input.y, 1);
  mpz_set_ui(other.x, 0); mpz_set_ui(other.y, 6);
  CHECK(_ec_add_AB(a, n, input, other, &output, d, m, t1, t2));
  CHECK(mpz_cmp_ui(d, 5) == 0);
  /* Ordinary addition with a nonunit x difference, on b=3. */
  mpz_set_ui(input.x, 3); mpz_set_ui(input.y, 6);
  mpz_set_ui(other.x, 13); mpz_set_ui(other.y, 14);
  CHECK(_ec_add_AB(a, n, input, other, &output, d, m, t1, t2));
  CHECK(mpz_cmp_ui(d, 5) == 0);
  mpz_clears(a,n,k,d,m,t1,t2,input.x,input.y,other.x,other.y,
             output.x,output.y,NULL);
}

#ifndef ECM_CHECK_AFFINE_ONLY
static void check_random(void)
{
  uint64_t state = 0, copy;
  mpz_t n, r, temporary, expected, repeated;
  unsigned i;
  CHECK(ec_rand64(&state) == UINT64_C(0xe220a8397b1dcdaf));
  CHECK(ec_rand64(&state) == UINT64_C(0x6e789e6aa1b965f4));
  CHECK(ec_rand64(&state) == UINT64_C(0x06c45d188009454f));
  mpz_inits(n,r,temporary,expected,repeated,NULL);
  mpz_set_str(n, "ffffffffffffffffffffffffffffffff", 16);
  mpz_set_str(expected, "e220a8397b1dcdaf6e789e6aa1b965f4", 16);
  state = 0;
  ec_random_mod(r, temporary, n, &state);
  CHECK(mpz_cmp(r, expected) == 0); /* Multiword, endian/UV-independent. */
  mpz_set_ui(n, 97);
  for (i = 0; i < 256; i++) {
    copy = state;
    ec_random_mod(r, temporary, n, &state);
    ec_random_mod(repeated, temporary, n, &copy);
    CHECK(mpz_cmp(r, repeated) == 0 && state == copy);
    CHECK(mpz_sgn(r) >= 0 && mpz_cmp(r, n) < 0);
  }
  mpz_clears(n,r,temporary,expected,repeated,NULL);
}

static UV gcd_uv(UV a, UV c)
{ while (c != 0) { UV t = a%c; a = c; c = t; } return a; }

static void check_scalar(ecm_ctx *ctx, reference_point P, unsigned long p, UV k,
                         unsigned long scale)
{
  reference_point expected = multiply(P, k, 6, 1, p);
  mpz_t x, z;
  int binary;
  mpz_inits(x, z, NULL);
  mpz_set_ui(ctx->n, p); mpz_set_ui(ctx->b, 2); /* (Montgomery A+2)/4 */
  for (binary = 0; binary <= 1; binary++) {
    mpz_set_ui(x, (unsigned long)((uint64_t)P.x*scale % p));
    mpz_set_ui(z, scale);
    if (binary) ec_mult_binary(ctx, k, x, z);
    else        ec_mult(ctx, k, x, z);
    if (expected.infinity) {
      CHECK(mpz_sgn(z) == 0 && mpz_sgn(x) != 0);
    } else {
      int ok = mpz_sgn(z) != 0 && mpz_cmp_ui(x,
        (unsigned long)((uint64_t)expected.x*mpz_get_ui(z) % p)) == 0;
      if (!ok)
        fprintf(stderr, "scalar=%"UVuf" binary=%d p=%lu\n", k, binary, p);
      CHECK(ok);
    }
  }
  mpz_clears(x, z, NULL);
}

static void check_projective(ecm_ctx *ctx)
{
  const double ratios[] = {1.61803398875, 1.72360679775, 1.618347119656,
                           1.617914406529, 1.58017872826};
  const reference_point old_case = {3, 293, 0};
  const UV powers[] = {3,9,25,27,49,81,121,125,169,729,2401,19683};
  const UV boundaries[] = {UV_MAX/6-2, UV_MAX/6-1, UV_MAX/6,
                            UV_MAX/6+1, UV_MAX/6+2, UV_MAX-2, UV_MAX-1, UV_MAX};
  reference_point P = {0,0,0}, torsion = {0,0,0};
  const unsigned long p = 1000000007UL;
  UV k;
  unsigned i;
  for (P.x = 2; ; P.x++) {
    int64_t rhs = mod(P.x*P.x*P.x + 6*P.x*P.x + P.x, p);
    P.y = power(rhs, (p+1)/4, p);
    if (P.y != 0 && P.y*P.y % p == rhs) break;
  }
  for (i = 0; i < sizeof(powers)/sizeof(powers[0]); i++)
    check_scalar(ctx, old_case, 1009, powers[i], 1);
  for (k = 0; k <= 4096; k++) {
    check_scalar(ctx, P, p, k, k%17+1);
    if (k > 2 && (k & 1)) {
      for (i = 0; i < sizeof(ratios)/sizeof(ratios[0]); i++) {
        UV r = (UV)((double)k/ratios[i]+0.5);
        /* Coprime starts must work; some noncoprime starts also finish
         * at one after absorbing the shared factor in a PRAC step. */
        CHECK(gcd_uv(k,r) != 1 || lucas_cost(k,r) != UV_MAX);
      }
    }
  }
  /* All five original choices for 7^4 are invalid. */
  for (i = 0; i < sizeof(ratios)/sizeof(ratios[0]); i++)
    CHECK(lucas_cost(2401, (UV)(2401/ratios[i]+0.5)) == UV_MAX);
  CHECK(lucas_cost(2401, 2401/2+1) != UV_MAX);
  CHECK(lucas_cost(9, 6) == UV_MAX);
  CHECK(lucas_cost(189, 120) != UV_MAX);
  for (i = 0; i < sizeof(boundaries)/sizeof(boundaries[0]); i++)
    check_scalar(ctx, P, p, boundaries[i], 3);
#if UV_MAX > UINT32_MAX
  check_scalar(ctx, P, p, (UV)UINT32_MAX, 5);
  check_scalar(ctx, P, p, (UV)UINT32_MAX+2, 5);
  check_scalar(ctx, P, p, UVCONST(10000000019), 5);
#endif
  for (k = 0; k < 6; k++) check_scalar(ctx, torsion, p, k, 1);
}

static void check_stage2(ecm_ctx *ctx)
{
  const UV limits[] = {13,19,31,59,127};
  unsigned i;
  mpz_t x,z,f;
  mpz_inits(x,z,f,NULL);
  /* Valid Q with order 13 modulo 1009 and order >100 modulo 1051.
   * The left half of the giant-step interval must find factor 1009. */
  mpz_set_ui(ctx->n, 1060459); mpz_set_ui(ctx->b, 2);
  for (i = 0; i < sizeof(limits)/sizeof(limits[0]); i++) {
    mpz_set_ui(x, 972917); mpz_set_ui(z, 1);
    CHECK(ec_stage2(ctx, 7, limits[i], x, z, f) != 0);
    CHECK(mpz_cmp_ui(f, 1009) == 0);
  }
  /* Empty ranges do no allocation; B2=2 exercises the minimum even D. */
  i = allocation_calls;
  CHECK(ec_stage2(ctx, 0, 1, x, z, f) == 0);
  CHECK(ec_stage2(ctx, 7, 7, x, z, f) == 0);
  CHECK(allocation_calls == i);
  mpz_set_ui(x, 972917); mpz_set_ui(z, 1);
  CHECK(ec_stage2(ctx, 0, 2, x, z, f) == 0);

  /* Normalization reports a factor or a lost curve explicitly, before
   * allocating a table.  Neither is a successful normalization. */
  mpz_set_ui(ctx->n, 35);
  mpz_set_ui(x, 6); mpz_set_ui(z, 2);
  CHECK(ec_normalize(ctx, x, z, f));
  CHECK(mpz_cmp_ui(x, 3) == 0 && mpz_cmp_ui(z, 1) == 0);
  CHECK(mpz_cmp_ui(f, 1) == 0);
  i = allocation_calls;
  mpz_set_ui(z, 5);
  CHECK(ec_stage2(ctx, 7, 19, x, z, f) == 2);
  CHECK(mpz_cmp_ui(f, 5) == 0);
  mpz_set_ui(z, 35);
  CHECK(ec_stage2(ctx, 7, 19, x, z, f) == 0);
  CHECK(mpz_cmp_ui(f, 35) == 0);
  CHECK(allocation_calls == i);
  mpz_clears(x,z,f,NULL);
}

static void check_public_tiny_bounds(void)
{
  mpz_t n, f;
  UV b1, b2;
  mpz_init_set_str(n, "1000036000099", 10); mpz_init(f);
  for (b1 = 0; b1 <= 3; b1++) {
    for (b2 = 0; b2 <= 3; b2++) {
      int found = _GMP_ecm_factor_projective(n, f, b1, b2, 1, UINT64_C(11), 0);
      if (found)
        CHECK(mpz_cmp_ui(f, 1) > 0 && mpz_cmp(f, n) < 0 &&
              mpz_divisible_p(n, f));
      else CHECK(found == 0);
    }
  }
  mpz_clears(n,f,NULL);
}

#ifdef ECM_CHECK_THREADS
#define THREAD_JOBS 12
typedef struct { int status, valid; unsigned long factor; } ecm_result;
typedef struct { unsigned id; ecm_result results[THREAD_JOBS]; } ecm_task;

static void *run_ecm_task(void *argument)
{
  const char *numbers[] = {"1060459", "100160063", "1000036000099"};
  ecm_task *task = (ecm_task*)argument;
  unsigned i;
  mpz_t n,f;
  mpz_inits(n,f,NULL);
  for (i = 0; i < THREAD_JOBS; i++) {
    ecm_result *result = &task->results[i];
    UV b1 = 40+5*i;
    uint64_t seed = ((uint64_t)task->id << 32) | i;
    mpz_set_str(n, numbers[(task->id+i)%3], 10);
    result->status = _GMP_ecm_factor_projective(n,f,b1,10*b1,2,seed,0);
    result->valid = result->status == 0 ||
      ((result->status == 1 || result->status == 2) &&
       mpz_cmp_ui(f,1) > 0 && mpz_cmp(f,n) < 0 && mpz_divisible_p(n,f) &&
       mpz_fits_ulong_p(f));
    result->factor = result->status ? mpz_get_ui(f) : 0;
  }
  mpz_clears(n,f,NULL);
  return NULL;
}

static void check_concurrent(void)
{
  pthread_t threads[4];
  ecm_task serial[4], parallel[4];
  unsigned i,j;
  for (i = 0; i < 4; i++) {
    serial[i].id = parallel[i].id = i;
    run_ecm_task(&serial[i]);
  }
  /* The prime cache is already initialized, and each call owns its numbers,
   * scratch and seed.  Check full-driver results against matching serial calls. */
  for (i = 0; i < 4; i++)
    CHECK(pthread_create(&threads[i],NULL,run_ecm_task,&parallel[i]) == 0);
  for (i = 0; i < 4; i++) {
    CHECK(pthread_join(threads[i],NULL) == 0);
    for (j = 0; j < THREAD_JOBS; j++) {
      CHECK(serial[i].results[j].valid && parallel[i].results[j].valid);
      CHECK(serial[i].results[j].status == parallel[i].results[j].status);
      CHECK(serial[i].results[j].factor == parallel[i].results[j].factor);
    }
  }
  CHECK(pthread_mutex_destroy(&allocation_mutex) == 0);
}
#endif

#if defined(ECM_CHECK_ALLOC32) && UV_MAX > UINT32_MAX
static void probe_stage2_allocation(ecm_ctx *ctx, UV bound, int expected, UV count)
{
  /* Static GMP objects remain defined across longjmp after normalization
   * has modified them.  A nonvolatile automatic object would not. */
  static mpz_t x,z,f;
  int status;
  mpz_inits(x,z,f,NULL);
  mpz_set_ui(ctx->n, 1060459); mpz_set_ui(ctx->b, 2);
  mpz_set_ui(x, 972917); mpz_set_ui(z, 1);
  probed_count = 0;
  allocation_probe = 1;
  switch (setjmp(allocation_trap)) {
    case 0:
      ec_stage2(ctx, 7, bound, x, z, f);
      CHECK(0); /* Must reject or reach the nonallocating hook. */
      status = 0;
      break;
    case 1: status = 1; break;
    case 2: status = 2; break;
    default: status = -1; break;
  }
  allocation_probe = 0;
  CHECK(status == expected);
  if (status == 1)
    CHECK(strcmp(error_text, "ECM stage-2 allocation size overflow") == 0);
  else
    CHECK(probed_count == count);
  mpz_clears(x,z,f,NULL);
}

static void check_allocation_limit(ecm_ctx *ctx)
{
  UV d = ((UV)(SIZE_MAX / sizeof(mpz_t)) - 1)/2;
  d -= d%2;
  /* Last feasible even width, then the next even width. */
  probe_stage2_allocation(ctx, 2*d*d, 2, 2*d+1);
  d += 2;
  probe_stage2_allocation(ctx, 2*d*d, 1, 0);
  /* Old 64-to-32 count truncation case: D=2^31, count=2^32+1. */
  probe_stage2_allocation(ctx, UVCONST(9223372036854775808), 1, 0);
}
#endif
#endif

int main(void)
{
  check_affine();
#ifndef ECM_CHECK_AFFINE_ONLY
  ecm_ctx context, *ctx = &context;
  mpz_t modulus;
  mpz_init_set_ui(modulus, 1060459);
  ec_ctx_init(ctx, modulus, 0);
  prime_iterator_global_startup();
  check_random();
  check_projective(ctx);
  check_stage2(ctx);
#if defined(ECM_CHECK_ALLOC32) && UV_MAX > UINT32_MAX
  check_allocation_limit(ctx);
#endif
  ec_ctx_clear(ctx);
  mpz_clear(modulus);
  check_public_tiny_bounds();
#ifdef ECM_CHECK_THREADS
  check_concurrent();
#endif
  prime_iterator_global_shutdown();
#ifdef ECM_CHECK_ALLOC32
  printf("PASS ECM arithmetic (%u-bit UV, modeled 32-bit allocation ceiling): %lu checks\n",
#else
  printf("PASS ECM arithmetic (%u-bit UV): %lu checks\n",
#endif
         (unsigned)(sizeof(UV)*8), checks);
#else
  printf("PASS certificate-verifier affine arithmetic: %lu checks\n", checks);
#endif
  return 0;
}
