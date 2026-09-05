/*============================================================================

  Self-initializing quadratic sieve

  This is an independent implementation for Math::Prime::Util::GMP.  Its
  polynomial, sieve, relation graph, and matrix interfaces are deliberately
  local to this file.  The linear algebra is supplied by lanczos.c.

  For a polynomial family we choose square-free a, d in {1,2}, D = d*a,
  and B so that the following normalized polynomial is integral:

      q(x) = ((D*x + B)^2 - kN) / D.

  A relation therefore gives

      (D*x + B)^2 = D*q(x) (mod N).

  The prime factors of D are included in every relation.  One-large-prime
  partials use a dedicated anchor table; two-large-prime partials use a
  general graph (the value 1 is a distinguished graph vertex).  Each repeated
  1LP label or fundamental graph cycle is turned into a full relation before
  linear algebra, so the matrix contains factor-base rows only.

  Copyright (c) 2026 Dana Jacobsen
  Written by Dana Jacobsen, September 2026, with assistance from OpenAI Codex.

============================================================================*/

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>

#include "ptypes.h"
#include "siqs.h"
#include "factor.h"
#include "lanczos.h"
#include "pbrent63.h"
#include "primality.h"
#include "prime_iterator.h"
#include "squfof126.h"
#include "utility.h"

#ifndef UINT8_MAX
# define UINT8_MAX 255U
#endif
#ifndef UINT16_MAX
# define UINT16_MAX 65535U
#endif
#ifndef UINT32_MAX
# define UINT32_MAX 4294967295U
#endif
#ifndef UINT64_C
# define UINT64_C(value) ((uint64_t)(value))
#endif

#ifndef M_LN2
# define M_LN2 0.69314718055994530942
#endif

#define SIQS_TRIAL_LIMIT        1000U
#define SIQS_MAX_EXTRA_RELS      512U
/* Exact and randomized solving both had ample dependency yield with a
 * 32-column reduced-core surplus throughout the tuned 1LP range.  Keep the
 * older conservative surplus for the independently tuned 2LP path. */
#define SIQS_MATRIX_EXTRA_RELS_1LP_DEFAULT 32U
#define SIQS_MATRIX_EXTRA_RELS_2LP_DEFAULT 64U
#define SIQS_MATRIX_CHECK_MIN     32U
#define SIQS_MATRIX_CHECK_MAX    512U
#define SIQS_MATRIX_RETRY_BATCH_MAX 128U
#define SIQS_EVAL_MAX_EXTRA_FACTORS 18U
#define SIQS_EVAL_INITIAL_FACTORS   64U
/* The packed dense solver wins consistently through about 1200 reduced
 * columns and reaches parity before 1300.  Select by the matrix we actually
 * have, rather than using N bits as a proxy; larger cores use block Lanczos. */
#define SIQS_DENSE_SOLVER_MAX_COLS 1280U
#define SIQS_LP_MAX UINT64_C(0x0000000fffffffff)
#define SIQS_RESIDUAL_PRODUCT_MAX UINT64_C(0x7fffffffffffffff)
#define SIQS_NO_ROOT       UINT32_MAX
#define SIQS_NO_INDEX      UINT32_MAX
#define SIQS_SIEVE_ALIGN          256U
#define SIQS_SIEVE_BLOCK_BITS      15U
#define SIQS_SIEVE_BLOCK_SIZE  (1U << SIQS_SIEVE_BLOCK_BITS)
#define SIQS_SIEVE_BLOCK_MASK  (SIQS_SIEVE_BLOCK_SIZE - 1U)
/* Direct candidate tests are substantially cheaper than sparse candidate-map
 * walks on the measured M1 Pro.  Fixed-work sweeps found a broad plateau from
 * about 6 through the all-direct limit; timer-free full factors selected 8 as
 * the smallest simple value on that plateau. */
#define SIQS_RESIEVE_CUTOFF_COEFFICIENT       8U
#define SIQS_RESIEVE_CUTOFF_FLOOR          1024U
#define SIQS_MATRIX_EXTRA_RELS(ctx) \
    ((ctx)->params.max_large_primes == 1 \
      ? SIQS_MATRIX_EXTRA_RELS_1LP_DEFAULT \
      : SIQS_MATRIX_EXTRA_RELS_2LP_DEFAULT)
/* The bucket sieve remains available for substantially larger intervals.
 * Contiguous sieving wins throughout the production-tuned bit range. */
#define SIQS_SIEVE_BLOCK_MIN UINT32_MAX
/* Automatic d retains d=1 unless the selected kN is 1 modulo 8.  A
 * boundary-enriched full-run sweep favored a 3/16*ln(2) preference for the
 * eligible multipliers; ordinary 160-220-bit samples confirmed that allowing
 * d=2 is healthy below the range where that preference most often matters. */
#define SIQS_MULTIPLIER_D2_BONUS_SIXTEENTHS 3U
/* Multiplier selection is worth shortening while it is a material part of a
 * full factorization.  Fresh full-factor sweeps selected all odd square-free
 * k <= 255 scored through FB/20 at 100--177 bits.  From 178 bits onward the
 * smaller k <= 127 set scored through min(FB,1000) was more stable; extending
 * the shallow search through 184 gained less than 1% and lost most pairwise
 * comparisons.  The 177/178 boundary also coincides with an existing policy
 * transition rather than introducing a new one solely for this selector. */
#ifndef SIQS_MULTIPLIER_MAX
# define SIQS_MULTIPLIER_MAX 255U
#endif
#ifndef SIQS_MULTIPLIER_FULL_MAX
# define SIQS_MULTIPLIER_FULL_MAX 127U
#endif
#ifndef SIQS_MULTIPLIER_SHALLOW_LAST_BITS
# define SIQS_MULTIPLIER_SHALLOW_LAST_BITS 177U
#endif
#ifndef SIQS_MULTIPLIER_SEARCH_DEPTH
# define SIQS_MULTIPLIER_SEARCH_DEPTH 1000U
#endif
#ifndef SIQS_MULTIPLIER_SEARCH_DIVISOR
# define SIQS_MULTIPLIER_SEARCH_DIVISOR 20U
#endif
#ifndef SIQS_MULTIPLIER_SEARCH_FLOOR
# define SIQS_MULTIPLIER_SEARCH_FLOOR 1U
#endif
/* The standard Knuth--Schroeppel size term is 0.5*ln(k).  Our smaller
 * historical 0.5*ln(2)*ln(k) penalty performed better in full-factor sweeps,
 * despite favoring larger k while the parameter policy is indexed by N rather
 * than kN, so retain it deliberately.
 * TODO: Retest the standard coefficient together with retuned FB, interval,
 * and large-prime parameters; changing this term alone is not a fair test of
 * a policy whose other parameters were tuned around the current k choices. */
#ifndef SIQS_MULTIPLIER_SIZE_PENALTY
# define SIQS_MULTIPLIER_SIZE_PENALTY (0.5 * M_LN2)
#endif
#if SIQS_MULTIPLIER_MAX < 1 || SIQS_MULTIPLIER_MAX > 4095 || \
    !(SIQS_MULTIPLIER_MAX & 1)
# error "SIQS_MULTIPLIER_MAX must be an odd value from 1 through 4095"
#endif
#if SIQS_MULTIPLIER_FULL_MAX < 1 || \
    SIQS_MULTIPLIER_FULL_MAX > SIQS_MULTIPLIER_MAX || \
    !(SIQS_MULTIPLIER_FULL_MAX & 1)
# error "SIQS_MULTIPLIER_FULL_MAX must be odd and no larger than MAX"
#endif
#if SIQS_MULTIPLIER_SEARCH_DEPTH < 1
# error "SIQS_MULTIPLIER_SEARCH_DEPTH must be positive"
#endif
#if SIQS_MULTIPLIER_SEARCH_FLOOR < 1
# error "SIQS_MULTIPLIER_SEARCH_FLOOR must be positive"
#endif
#define SIQS_MULTIPLIER_CAPACITY ((SIQS_MULTIPLIER_MAX + 1U) / 2U)
#define SIQS_BUCKET_FB_LIMIT   (1U << (32U - SIQS_SIEVE_BLOCK_BITS))
#define SIQS_POSTFILTER_MAX_SMALL  256U
#define SIQS_HASH_EMPTY        UINT64_C(0)

typedef struct {
  uint32_t p;
  uint32_t sqrt_kn;
  uint8_t logp;
  uint8_t in_a;
} siqs_fb_t;

typedef struct {
  uint32_t row;
  uint32_t exponent;
} siqs_factor_t;

typedef struct {
  mpz_t y;
  mpz_t q;
  mpz_t rest;
  siqs_factor_t *factors;
  uint32_t factor_alloc;
} siqs_eval_workspace_t;

typedef struct {
  mpz_t y;
  siqs_factor_t *factors;
  uint32_t nfactors;
  uint64_t lp1;
  uint64_t lp2;
  uint64_t fingerprint;
} siqs_raw_relation_t;

typedef struct {
  mpz_t y;
  siqs_factor_t *factors;
  uint32_t nfactors;
} siqs_full_relation_t;

typedef struct {
  uint64_t label;
  siqs_raw_relation_t *relation;
} siqs_one_lp_anchor_t;

typedef struct {
  siqs_one_lp_anchor_t *anchors;
  uint32_t count;
  uint32_t alloc;
  mpz_t y;
  mpz_t large_prime;
  mpz_t inverse;
  int initialized;
} siqs_one_lp_state_t;

typedef struct {
  mpz_t *values;
  uint32_t count;
  uint32_t alloc;
} siqs_factor_array_t;

typedef struct {
  const char *policy_name;
  uint32_t bits;
  uint32_t policy_first_bits;
  uint32_t policy_last_bits;
  uint32_t max_large_primes;
  uint32_t one_lp_policy_multiplier;
  uint32_t lp_policy_multiplier_floor;
  uint32_t lp_policy_product_multiplier_floor;
  uint32_t fb_size;
  uint32_t half_interval;
  uint32_t q_count;
  uint32_t poly_d;
  uint32_t target_relations;
  uint32_t sieve_start;
  uint32_t sieve_free_units;
  uint32_t sieve_start_prime_floor;
  uint32_t fb_floor;
  uint32_t relation_extra;
  uint8_t sieve_initial;
  uint8_t stage1_bias;
  uint64_t large_prime_bound;
  uint64_t smooth_bound;
  double fb_coefficient;
  double interval_scale;
  double smooth_bound_exponent;
  double sieve_hit_exponent;
  double sieve_hit_bound_nominal;
  double sieve_hit_residual_multiplier;
  double sieve_hit_bound;
  double log_scale;
} siqs_parameters_t;

typedef struct {
  uint64_t state;
} siqs_rng_t;

typedef struct {
  uint64_t *slots;
  uint32_t alloc;
  uint32_t count;
} siqs_hashset_t;

typedef struct {
  uint64_t value;
  uint32_t dsu_parent;
  uint32_t dsu_size;
  uint32_t tree_parent;
  uint32_t parent_edge;
  uint32_t mark;
} siqs_vertex_t;

typedef struct {
  siqs_vertex_t *vertices;
  uint32_t vertex_count;
  uint32_t vertex_alloc;
  uint32_t *table;
  uint32_t table_alloc;
  uint32_t mark;
} siqs_graph_t;

typedef struct {
  /* A is the square-free factor-base product; DA = d*A is both the linear
   * coefficient and the divisor of the normalized polynomial. */
  mpz_t A;
  mpz_t DA;
  mpz_t B;
  mpz_t C;
  mpz_t target_A;
  mpz_t *H;
  uint32_t *a_index;
  uint32_t *corrections;
  uint32_t q_count;
  uint32_t b_index;
  uint32_t b_limit;
} siqs_poly_t;

typedef struct {
  int32_t x;
  uint32_t first_hit;
  uint8_t sieve_score;
} siqs_candidate_t;

typedef struct {
  uint32_t fb_index;
  uint32_t next;
} siqs_hit_t;

typedef struct {
  mpz_srcptr original_n;
  mpz_t n;
  mpz_t kn;
  unsigned long multiplier;
  siqs_parameters_t params;
  siqs_fb_t *fb;
  uint32_t *prime;
  uint32_t *root1;
  uint32_t *root2;
  uint32_t *fb_reciprocal;
  uint32_t largest_fb_prime;
  uint8_t *sieve;
  uint32_t sieve_length;
  uint32_t sieve_offset;
  uint32_t active_sieve_length;
  uint8_t active_sieve_initial;
  uint32_t block_count;
  uint32_t block_large_index;
  uint32_t *bucket_bounds;
  uint32_t *bucket_fill;
  uint32_t *bucket_events;
  uint32_t bucket_event_count;
  uint32_t bucket_event_alloc;
  int use_block_sieve;
  uint16_t *candidate_at;
  uint32_t *candidate_at_wide;
  int candidate_wide;
  siqs_candidate_t *candidates;
  uint32_t candidate_count;
  uint32_t candidate_alloc;
  siqs_hit_t *hits;
  uint32_t hit_count;
  uint32_t hit_alloc;
  siqs_raw_relation_t **raw;
  uint32_t raw_count;
  uint32_t raw_alloc;
  siqs_full_relation_t **full;
  uint32_t full_count;
  uint32_t full_alloc;
  siqs_graph_t graph;
  siqs_one_lp_state_t one_lp;
  siqs_hashset_t relation_hashes;
  siqs_hashset_t a_hashes;
  siqs_factor_array_t *result;
  siqs_rng_t poly_rng;
  siqs_rng_t cofactor_rng;
  siqs_rng_t la_rng;
  siqs_eval_workspace_t eval;
  uint32_t *factor_counts;
  uint32_t *factor_touched;
  uint32_t factor_touched_count;
  uint64_t total_candidates;
  uint64_t accepted_smooth;
  uint64_t accepted_one_lp;
  uint64_t accepted_two_lp;
  uint64_t split_attempts;
  uint64_t split_squares;
  uint64_t split_squfof;
  uint64_t split_rho;
  uint64_t split_failures;
  int factor_found;
} siqs_ctx_t;

static int siqs_use_one_lp_relation_path(const siqs_ctx_t *ctx) {
  return ctx->params.max_large_primes == 1;
}

/*----------------------------------------------------------------------------
 * Basic allocation, random, and integer helpers
 *----------------------------------------------------------------------------*/

static void *siqs_malloc(size_t size) {
  void *p = malloc(size ? size : 1);
  if (p == NULL)
    croak("SIQS: unable to allocate memory");
  return p;
}

static void *siqs_calloc(size_t count, size_t size) {
  void *p = calloc(count ? count : 1, size ? size : 1);
  if (p == NULL)
    croak("SIQS: unable to allocate memory");
  return p;
}

static void *siqs_realloc(void *old, size_t size) {
  void *p = realloc(old, size ? size : 1);
  if (p == NULL)
    croak("SIQS: unable to grow allocation");
  return p;
}

static uint64_t siqs_mix64(uint64_t x) {
  x ^= x >> 30;
  x *= UINT64_C(0xbf58476d1ce4e5b9);
  x ^= x >> 27;
  x *= UINT64_C(0x94d049bb133111eb);
  x ^= x >> 31;
  return x ? x : UINT64_C(0x9e3779b97f4a7c15);
}

static uint64_t siqs_rand64(siqs_rng_t *rng) {
  uint64_t x = rng->state;
  x ^= x >> 12;
  x ^= x << 25;
  x ^= x >> 27;
  rng->state = x;
  return x * UINT64_C(2685821657736338717);
}

static uint32_t siqs_rand_range(siqs_rng_t *rng, uint32_t limit) {
  return limit ? (uint32_t)(siqs_rand64(rng) % limit) : 0;
}


static void siqs_mpz_set_u64(mpz_t z, uint64_t value) {
  mpz_import(z, 1, -1, sizeof(value), 0, 0, &value);
}

static uint32_t siqs_inverse_u32(uint32_t a, uint32_t p) {
  int64_t t = 0, nt = 1;
  int64_t r = p, nr = a % p;
  while (nr != 0) {
    int64_t q = r / nr;
    int64_t z = t - q * nt;
    t = nt;
    nt = z;
    z = r - q * nr;
    r = nr;
    nr = z;
  }
  if (r != 1)
    return 0;
  if (t < 0)
    t += p;
  return (uint32_t)t;
}

static uint32_t siqs_powmod_u32(uint32_t a, uint32_t e, uint32_t p) {
  uint64_t base = a % p;
  uint64_t result = 1;
  while (e != 0) {
    if (e & 1U)
      result = result * base % p;
    e >>= 1;
    if (e != 0)
      base = base * base % p;
  }
  return (uint32_t)result;
}

static int siqs_legendre_u32(uint32_t a, uint32_t p) {
  uint32_t r;
  if (a == 0)
    return 0;
  r = siqs_powmod_u32(a, (p - 1) / 2, p);
  return r == 1 ? 1 : r == p - 1 ? -1 : 0;
}

/* Jacobi(a,n) for an odd n.  The multiplier search uses this only with a
 * prime n, where it is exactly the Legendre symbol, but the binary form is
 * substantially cheaper than modular exponentiation for every candidate. */
static int siqs_jacobi_odd_u32(uint32_t a, uint32_t n) {
  int symbol = 1;
  a %= n;
  while (a != 0) {
    while (!(a & 1U)) {
      a >>= 1;
      if ((n & 7U) == 3U || (n & 7U) == 5U)
        symbol = -symbol;
    }
    {
      uint32_t t = a;
      a = n;
      n = t;
    }
    if ((a & 3U) == 3U && (n & 3U) == 3U)
      symbol = -symbol;
    a %= n;
  }
  return n == 1 ? symbol : 0;
}

static uint32_t siqs_sqrtmod_u32(uint32_t n, uint32_t p) {
  uint32_t q, s, z, c, x, t, m;
  if (p == 2)
    return n & 1U;
  n %= p;
  if (n == 0)
    return 0;
  if ((p & 3U) == 3U)
    return siqs_powmod_u32(n, (p + 1) / 4, p);

  q = p - 1;
  s = 0;
  while (!(q & 1U)) {
    q >>= 1;
    s++;
  }
  for (z = 2; siqs_legendre_u32(z, p) != -1; z++)
    ;
  c = siqs_powmod_u32(z, q, p);
  x = siqs_powmod_u32(n, (q + 1) / 2, p);
  t = siqs_powmod_u32(n, q, p);
  m = s;
  while (t != 1) {
    uint32_t i, b, tt = t;
    for (i = 1; i < m; i++) {
      tt = (uint32_t)((uint64_t)tt * tt % p);
      if (tt == 1)
        break;
    }
    if (i == m)
      return 0;
    b = siqs_powmod_u32(c, 1U << (m - i - 1), p);
    x = (uint32_t)((uint64_t)x * b % p);
    c = (uint32_t)((uint64_t)b * b % p);
    t = (uint32_t)((uint64_t)t * c % p);
    m = i;
  }
  return x;
}

static uint32_t siqs_ctz32(uint32_t n) {
  uint32_t c = 0;
  while (!(n & 1U)) {
    n >>= 1;
    c++;
  }
  return c;
}

/*----------------------------------------------------------------------------
 * Result factor partition
 *----------------------------------------------------------------------------*/

static void siqs_factor_array_init(siqs_factor_array_t *fa, const mpz_t n) {
  fa->alloc = 16;
  fa->count = 1;
  fa->values = (mpz_t *)siqs_malloc(fa->alloc * sizeof(mpz_t));
  mpz_init_set(fa->values[0], n);
}

static void siqs_factor_array_append(siqs_factor_array_t *fa, const mpz_t n) {
  if (fa->count == fa->alloc) {
    fa->alloc *= 2;
    fa->values = (mpz_t *)siqs_realloc(fa->values,
                                       fa->alloc * sizeof(mpz_t));
  }
  mpz_init_set(fa->values[fa->count++], n);
}

/* Refine the represented multiplicative partition by every copy of d. */
static int siqs_insert_divisor(siqs_factor_array_t *fa, const mpz_t d) {
  uint32_t i;
  int changed = 0;
  mpz_t g, q;
  if (mpz_cmp_ui(d, 1) <= 0)
    return 0;
  mpz_init(g);
  mpz_init(q);
  for (i = 0; i < fa->count; i++) {
    mpz_gcd(g, fa->values[i], d);
    if (mpz_cmp_ui(g, 1) <= 0 || mpz_cmp(g, fa->values[i]) == 0)
      continue;
    mpz_divexact(q, fa->values[i], g);
    mpz_set(fa->values[i], g);
    siqs_factor_array_append(fa, q);
    changed = 1;
  }
  mpz_clear(g);
  mpz_clear(q);
  return changed;
}

static int siqs_all_factors_prime(const siqs_factor_array_t *fa) {
  uint32_t i;
  for (i = 0; i < fa->count; i++)
    if (!_GMP_is_prob_prime(fa->values[i]))
      return 0;
  return 1;
}

static void siqs_verify_partition(const mpz_t n,
                                  const siqs_factor_array_t *fa) {
  uint32_t i;
  mpz_t product;
  mpz_init_set_ui(product, 1);
  for (i = 0; i < fa->count; i++) {
    if (mpz_cmp_ui(fa->values[i], 1) <= 0)
      croak("SIQS: invalid result factor");
    mpz_mul(product, product, fa->values[i]);
  }
  if (mpz_cmp(product, n) != 0)
    croak("SIQS: result factors do not multiply to the input");
  mpz_clear(product);
}

void _GMP_siqs_free(mpz_t *factors, uint32_t nfactors) {
  uint32_t i;
  for (i = 0; i < nfactors; i++)
    mpz_clear(factors[i]);
  free(factors);
}

/*----------------------------------------------------------------------------
 * Small open-addressed set used for relation and polynomial fingerprints
 *----------------------------------------------------------------------------*/

static void siqs_hashset_init(siqs_hashset_t *set, uint32_t initial) {
  uint32_t n = 16;
  while (n < initial)
    n <<= 1;
  set->slots = (uint64_t *)siqs_calloc(n, sizeof(uint64_t));
  set->alloc = n;
  set->count = 0;
}

static void siqs_hashset_clear(siqs_hashset_t *set) {
  free(set->slots);
  memset(set, 0, sizeof(*set));
}

static void siqs_hashset_grow(siqs_hashset_t *set) {
  uint64_t *old = set->slots;
  uint32_t oldn = set->alloc, i;
  set->alloc *= 2;
  set->slots = (uint64_t *)siqs_calloc(set->alloc, sizeof(uint64_t));
  set->count = 0;
  for (i = 0; i < oldn; i++) {
    uint64_t h = old[i];
    if (h != SIQS_HASH_EMPTY) {
      uint32_t pos = (uint32_t)h & (set->alloc - 1);
      while (set->slots[pos] != SIQS_HASH_EMPTY)
        pos = (pos + 1) & (set->alloc - 1);
      set->slots[pos] = h;
      set->count++;
    }
  }
  free(old);
}

/* Return 1 for a newly inserted value and 0 when it was already present. */
static int siqs_hashset_insert(siqs_hashset_t *set, uint64_t h) {
  uint32_t pos;
  h = siqs_mix64(h);
  if (set->count * 10 >= set->alloc * 7)
    siqs_hashset_grow(set);
  pos = (uint32_t)h & (set->alloc - 1);
  while (set->slots[pos] != SIQS_HASH_EMPTY) {
    if (set->slots[pos] == h)
      return 0;
    pos = (pos + 1) & (set->alloc - 1);
  }
  set->slots[pos] = h;
  set->count++;
  return 1;
}

/*----------------------------------------------------------------------------
 * Parameter selection, multiplier, and factor base
 *----------------------------------------------------------------------------*/

typedef struct {
  double base;
  double step;
  uint16_t origin_bits;
  uint8_t staged;
} siqs_policy_linear_t;

typedef struct {
  double amount;
  int16_t numerator_at_first;
  int8_t numerator_step;
  uint8_t denominator;
} siqs_policy_ratio_t;

typedef struct {
  const char *name;
  uint16_t first_bits;
  uint16_t last_bits;
  uint8_t production_large_primes;
  uint8_t q_count;
  uint8_t stage1_bias_base;
  uint8_t stage1_bias_step_bits;
  uint16_t stage1_bias_origin_bits;
  uint8_t one_lp_multiplier;
  uint8_t two_lp_multiplier_floor;
  uint8_t two_lp_product_floor;
  uint8_t sieve_free_units;
  uint16_t sieve_start_prime_floor;
  uint16_t fb_floor;
  uint16_t relation_extra;
  siqs_policy_linear_t fb_coefficient;
  siqs_policy_linear_t interval_base_scale;
  siqs_policy_ratio_t interval_lp_scale;
  double smooth_bound_exponent;
  siqs_policy_linear_t sieve_hit_exponent;
} siqs_policy_band_t;

typedef struct {
  const siqs_policy_band_t *band;
  uint32_t max_large_primes;
  uint32_t q_count;
  uint32_t one_lp_multiplier;
  uint32_t lp_multiplier_floor;
  uint32_t lp_product_floor;
  uint32_t sieve_free_units;
  uint32_t sieve_start_prime_floor;
  uint32_t fb_floor;
  uint32_t relation_extra;
  uint8_t stage1_bias;
  double fb_coefficient;
  double interval_scale;
  double smooth_bound_exponent;
  double sieve_hit_exponent;
} siqs_policy_t;

#define SIQS_POLICY_LINEAR(base, step, origin) \
  { (base), (step), (origin), 0 }
#define SIQS_POLICY_STAGED_LINEAR(base, step, origin) \
  { (base), (step), (origin), 1 }
#define SIQS_POLICY_RATIO(amount, numerator, step, denominator) \
  { (amount), (numerator), (step), (denominator) }

/*
 * This is the single production bit-size policy.  Rows are meaningful bands
 * formed by actual structural changes or interpolation anchors; they are not
 * independently tuned bins.  The interval has two factors because the old
 * model's N-size taper and its 2LP pre-ramp overlap from 231 through 239 bits.
 * Keeping them separate preserves that established calculation exactly.
 * Integer columns after the bit range are: LP count, q count, bias base,
 * bias step width/origin, 1LP K, conditional 2LP K/R floors, sieve byte
 * headroom, the first-sieved-prime floor, the factor-base floor, and the
 * initial relation surplus.
 *
 * The LP count is part of each complete policy row rather than an independent
 * crossover knob: changing it also requires changing the smooth exponent,
 * bounds, interval, and sieve-depth policy.  A crossover near 250 bits is
 * typical, but tests after adding the specialized 1LP path kept 1LP faster
 * through 228, found a CPU tie at 229--230 with about half the relation/graph
 * storage, and put the first repeatable 2LP CPU win at 231.  The established
 * two-LP parameter ramps still begin independently at 250 bits.
 *
 * The 1LP factor-base coefficients are joint collection/matrix choices, not
 * smooth-yield targets.  The lower schedule rises from 0.315 to 0.320 before
 * the 145-bit bridge, then the real K change at 185 starts a second schedule
 * which rises from 0.325 to 0.330 at the 210/211 boundary.  Joint geometry
 * tests keep interval scale 2.5 through 130 bits, ramp it to 3.0 at 140, hold
 * through 200, taper it back to 2.5 at 206, and then rejoin the established
 * 211--217 taper.  At the low end, q=2 has the healthiest A supply through
 * 49 bits, q=3 is faster from 50 through 80, and q=4 takes over at 81.  The
 * remaining adjacent full-factor tests put the q-count changes at 96, 117,
 * 145, 178, 206, and 218 bits; those choices remain independent of the
 * 231-bit LP-mode boundary.
 */
static const siqs_policy_band_t siqs_policy_bands[] = {
  /* These low rows remove the old 160-prime and 96-relation fixed-work floors.
   * K=1 sets the large-prime bound to pmax, so they deliberately collect
   * smooth relations only.  q=2 keeps A construction healthy at the new
   * lower limit; q=3 then wins until the measured 80/81 crossover.  Its
   * interval ramp begins at the fixed 4096 floor and meets the established
   * q=4 curve at that boundary.  Their matrices are too small for an early
   * readiness check.  Although direct target sweeps gave minimum CPU at
   * 0/4/2 extra relations, using 2/4/4 cost only 0.3% over 90,000 inputs,
   * cut matrix retries by 39%, and substantially reduced timing variance.
   * Full-factor sweeps put the return to q=5 and ordinary 1LP collection at
   * the 95/96 boundary. */
  { "smooth_k1_q2_low", MPU_SIQS_MIN_BITS, 49, 1, 2, 0, 0, 0,
    1, 60, 60, 8, 0, 48, 2,
    SIQS_POLICY_LINEAR(0.315, 0.0, 65),
    SIQS_POLICY_LINEAR(0.5, 0.0, 65),
    SIQS_POLICY_RATIO(0.0, 0, 0, 0), 0.12,
    SIQS_POLICY_LINEAR(0.15, 0.0, 65) },
  { "smooth_k1_q3_low", 50, 80, 1, 3, 0, 0, 0,
    1, 60, 60, 8, 0, 48, 4,
    SIQS_POLICY_LINEAR(0.315, 0.0, 65),
    SIQS_POLICY_LINEAR(0.0, 0.041666666666666667, 50),
    SIQS_POLICY_RATIO(0.0, 0, 0, 0), 0.12,
    SIQS_POLICY_LINEAR(0.15, 0.0, 65) },
  { "smooth_k1_q4_low", 81, 95, 1, 4, 0, 0, 0,
    1, 60, 60, 8, 0, 48, 4,
    SIQS_POLICY_LINEAR(0.315, 0.0, 65),
    SIQS_POLICY_LINEAR(0.5, 0.05, 65),
    SIQS_POLICY_RATIO(0.0, 0, 0, 0), 0.12,
    SIQS_POLICY_LINEAR(0.15, 0.0, 65) },
  /* Two extra relations cost 0.1% over this band while cutting matrix retries
   * by 68% and reducing timing variance.  Four cost 0.8% and added little. */
  { "one_lp_k2_q5", 96, 116, 1, 5, 0, 0, 0, 2, 60, 60, 8, 0, 160, 2,
    SIQS_POLICY_LINEAR(0.315, 0.0, 100),
    SIQS_POLICY_LINEAR(1.7, 0.03, 96),
    SIQS_POLICY_RATIO(0.0, 0, 0, 0), 0.12,
    SIQS_POLICY_LINEAR(0.15, 0.0, 100) },
  /* Exact dense elimination normally succeeds at the minimum relation count
   * through this band.  Let it try: a 15,500-input audit needed four second
   * matrix attempts, while avoiding the fixed readiness surplus reduced
   * total CPU by about 7%.  The ordinary retry loop handles that rare tail. */
  { "one_lp_k2_q6", 117, 129, 1, 6, 0, 0, 0, 2, 60, 60, 8, 0, 160, 0,
    SIQS_POLICY_LINEAR(0.315, 0.0, 117),
    SIQS_POLICY_LINEAR(2.5, 0.0, 117),
    SIQS_POLICY_RATIO(0.0, 0, 0, 0), 0.12,
    SIQS_POLICY_LINEAR(0.15, 0.0, 117) },
  /* Stopping at the first full-rank-sized matrix remained healthy through
   * 158 bits and saved about 1--5% across these bands.  The 151-bit score
   * transition is also the measured K3/K4 boundary, allowing one less row. */
  { "one_lp_k3_q6_interval_ramp", 130, 139, 1, 6, 0, 0, 0,
    3, 60, 60, 8, 0, 160, 0,
    SIQS_POLICY_LINEAR(0.315, 0.00033333333333333333, 130),
    SIQS_POLICY_LINEAR(2.5, 0.05, 130),
    SIQS_POLICY_RATIO(0.0, 0, 0, 0), 0.12,
    SIQS_POLICY_LINEAR(0.15, 0.0, 130) },
  { "one_lp_k3_q6", 140, 144, 1, 6, 0, 0, 0, 3, 60, 60, 8, 0, 160, 0,
    SIQS_POLICY_LINEAR(0.315, 0.00033333333333333333, 130),
    SIQS_POLICY_LINEAR(3.0, 0.0, 140),
    SIQS_POLICY_RATIO(0.0, 0, 0, 0), 0.12,
    SIQS_POLICY_LINEAR(0.15, 0.0, 140) },
  { "one_lp_k3_q7_bias_ramp", 145, 150, 1, 7, 0, 2, 143,
    3, 60, 60, 8, 0, 160, 0,
    SIQS_POLICY_LINEAR(0.32, 0.0, 146),
    SIQS_POLICY_LINEAR(3.0, 0.0, 146),
    SIQS_POLICY_RATIO(0.0, 0, 0, 0), 0.12,
    SIQS_POLICY_LINEAR(0.15, 0.0, 146) },
  { "one_lp_k4_q7_score_bias_ramp", 151, 158, 1, 7, 0, 2, 143,
    4, 60, 60, 8, 0, 160, 0,
    SIQS_POLICY_LINEAR(0.32, 0.0, 151),
    SIQS_POLICY_LINEAR(3.0, 0.0, 151),
    SIQS_POLICY_RATIO(0.0, 0, 0, 0), 0.12,
    SIQS_POLICY_STAGED_LINEAR(0.15, 0.0003, 150) },
  /* Zero surplus was 1.5% faster over this band and solved all 1,200 audited
   * inputs on the first matrix attempt. */
  { "one_lp_k4_q7", 159, 166, 1, 7, 8, 0, 0, 4, 60, 60, 8, 0, 160, 0,
    SIQS_POLICY_LINEAR(0.32, 0.0, 159),
    SIQS_POLICY_LINEAR(3.0, 0.0, 159),
    SIQS_POLICY_RATIO(0.0, 0, 0, 0), 0.12,
    SIQS_POLICY_STAGED_LINEAR(0.15, 0.0003, 150) },
  /* The readiness check makes the nominal 96 surplus nearly free here:
   * +32 and +96 produced identical work throughout a coarse 167--177
   * sample.  At the 167-bit lower edge, zero was 0.15% slower than +96 in
   * a fresh 600-input order-balanced confirmation, so retain +96. */
  { "one_lp_k5_q7", 167, 177, 1, 7, 8, 0, 0, 5, 60, 60, 8, 0, 160, 96,
    SIQS_POLICY_LINEAR(0.32, 0.0, 167),
    SIQS_POLICY_LINEAR(3.0, 0.0, 167),
    SIQS_POLICY_RATIO(0.0, 0, 0, 0), 0.12,
    SIQS_POLICY_STAGED_LINEAR(0.15, 0.0003, 150) },
  { "one_lp_k5_q8", 178, 184, 1, 8, 8, 0, 0, 5, 60, 60, 8, 0, 160, 96,
    SIQS_POLICY_LINEAR(0.32, 0.0, 178),
    SIQS_POLICY_LINEAR(3.0, 0.0, 178),
    SIQS_POLICY_RATIO(0.0, 0, 0, 0), 0.12,
    SIQS_POLICY_STAGED_LINEAR(0.15, 0.0003, 150) },
  { "one_lp_k8_q8", 185, 200, 1, 8, 8, 0, 0, 8, 60, 60, 8, 0, 160, 96,
    SIQS_POLICY_LINEAR(0.325, 0.0002, 185),
    SIQS_POLICY_LINEAR(3.0, 0.0, 185),
    SIQS_POLICY_RATIO(0.0, 0, 0, 0), 0.12,
    SIQS_POLICY_STAGED_LINEAR(0.15, 0.0003, 150) },
  { "one_lp_k8_q8_interval_taper", 201, 205, 1, 8, 8, 0, 0,
    8, 60, 60, 8, 0, 160, 96,
    SIQS_POLICY_LINEAR(0.325, 0.0002, 185),
    SIQS_POLICY_LINEAR(3.0, -0.0625, 200),
    SIQS_POLICY_RATIO(0.0, 0, 0, 0), 0.12,
    SIQS_POLICY_STAGED_LINEAR(0.15, 0.0003, 150) },
  { "one_lp_k8_q9", 206, 210, 1, 9, 8, 0, 0,
    8, 60, 60, 8, 0, 160, 96,
    SIQS_POLICY_LINEAR(0.325, 0.0002, 185),
    SIQS_POLICY_LINEAR(2.5, 0.0, 208),
    SIQS_POLICY_RATIO(0.0, 0, 0, 0), 0.12,
    SIQS_POLICY_STAGED_LINEAR(0.15, 0.0003, 150) },
  { "one_lp_k12_q9_taper", 211, 217, 1, 9, 8, 0, 0,
    12, 60, 60, 8, 0, 160, 96,
    SIQS_POLICY_LINEAR(0.33, 0.0, 211),
    SIQS_POLICY_LINEAR(3.0, -0.05, 200),
    SIQS_POLICY_RATIO(0.0, 0, 0, 0), 0.12,
    SIQS_POLICY_STAGED_LINEAR(0.15, 0.0003, 150) },
  { "one_lp_k12_q10_interval_taper", 218, 230, 1, 10, 12, 0, 0,
    12, 60, 60, 8, 0, 160, 96,
    SIQS_POLICY_LINEAR(0.332, -0.001, 218),
    SIQS_POLICY_LINEAR(2.5, -0.05, 210),
    SIQS_POLICY_RATIO(0.0, 0, 0, 0), 0.12,
    SIQS_POLICY_STAGED_LINEAR(0.15, 0.0003, 150) },
  { "two_lp_early_interval_taper", 231, 239, 2, 10, 12, 0, 0,
    0, 60, 60, 8, 0, 160, 96,
    SIQS_POLICY_LINEAR(0.30069720, 0.0, 231),
    SIQS_POLICY_LINEAR(2.5, -0.05, 210),
    SIQS_POLICY_RATIO(0.00537337256, 20, 0, 20), 0.16,
    SIQS_POLICY_LINEAR(0.205, 0.0, 231) },
  { "two_lp_early", 240, 249, 2, 10, 12, 0, 0, 0, 60, 60, 8, 0, 160, 96,
    SIQS_POLICY_LINEAR(0.30069720, 0.0, 240),
    SIQS_POLICY_LINEAR(1.0, 0.0, 240),
    SIQS_POLICY_RATIO(0.00537337256, 20, 0, 20), 0.16,
    SIQS_POLICY_LINEAR(0.205, 0.0, 240) },
  { "two_lp_fb_ramp_q10", 250, 259, 2, 10, 12, 0, 0,
    0, 60, 60, 8, 0, 160, 96,
    SIQS_POLICY_LINEAR(0.30069720, 0.000275913, 250),
    SIQS_POLICY_LINEAR(1.0, 0.0, 251),
    SIQS_POLICY_RATIO(0.00537337256, 20, -1, 20), 0.16,
    SIQS_POLICY_LINEAR(0.205, 0.0, 251) },
  { "two_lp_fb_ramp_q11", 260, 266, 2, 11, 12, 0, 0,
    0, 60, 60, 8, 0, 160, 96,
    SIQS_POLICY_LINEAR(0.30069720, 0.000275913, 250),
    SIQS_POLICY_LINEAR(1.0, 0.0, 260),
    SIQS_POLICY_RATIO(0.00537337256, 10, -1, 20), 0.16,
    SIQS_POLICY_LINEAR(0.205, 0.0, 260) },
  { "two_lp_score_release", 267, 269, 2, 11, 12, 0, 0,
    0, 60, 60, 8, 0, 160, 96,
    SIQS_POLICY_LINEAR(0.30069720, 0.000275913, 250),
    SIQS_POLICY_LINEAR(1.0, 0.0, 267),
    SIQS_POLICY_RATIO(0.00537337256, 3, -1, 20), 0.16,
    SIQS_POLICY_STAGED_LINEAR(0.17, 0.0003, 150) },
  { "two_lp_mid_ramp", 270, 299, 2, 11, 12, 0, 0,
    0, 60, 60, 8, 0, 160, 96,
    SIQS_POLICY_LINEAR(0.30621546, 0.000126151333, 270),
    SIQS_POLICY_LINEAR(1.0, 0.0, 270),
    SIQS_POLICY_RATIO(0.15231778066, 0, 1, 30), 0.16,
    SIQS_POLICY_STAGED_LINEAR(0.18, 0.0003, 150) },
  { "two_lp_high_q11", 300, 304, 2, 11, 12, 0, 0,
    0, 0, 0, 16, 384, 160, 96,
    SIQS_POLICY_LINEAR(0.31, 0.0, 301),
    SIQS_POLICY_LINEAR(1.0, 0.0, 301),
    SIQS_POLICY_RATIO(0.15231778066, 50, -1, 50), 0.16,
    SIQS_POLICY_STAGED_LINEAR(0.18, 0.0003, 150) },
  { "two_lp_high_q12", 305, MPU_SIQS_MAX_BITS, 2, 12, 12, 0, 0,
    0, 0, 0, 16, 384, 160, 96,
    SIQS_POLICY_LINEAR(0.31, 0.0, 305),
    SIQS_POLICY_LINEAR(1.0, 0.0, 305),
    SIQS_POLICY_RATIO(0.15231778066, 45, -1, 50), 0.16,
    SIQS_POLICY_STAGED_LINEAR(0.18, 0.0003, 150) }
};

#undef SIQS_POLICY_LINEAR
#undef SIQS_POLICY_STAGED_LINEAR
#undef SIQS_POLICY_RATIO

#define SIQS_POLICY_BAND_COUNT \
  ((uint32_t)(sizeof(siqs_policy_bands) / sizeof(siqs_policy_bands[0])))

static double siqs_policy_linear_value(const siqs_policy_linear_t *line,
                                       uint32_t bits) {
  int32_t distance = (int32_t)bits - (int32_t)line->origin_bits;
  volatile double delta;
  if (line->step == 0.0)
    return line->base;
  /* Most legacy ramps were direct expressions; leave contraction to the
   * compiler just as before.  The sieve increment was first assigned
   * separately, so preserve that distinct rounding path explicitly. */
  if (!line->staged)
    return line->base + line->step * (double)distance;
  delta = line->step * (double)distance;
  return line->base + delta;
}

static double siqs_policy_ratio_value(const siqs_policy_band_t *band,
                                      uint32_t bits) {
  const siqs_policy_ratio_t *ratio = &band->interval_lp_scale;
  int32_t numerator;
  volatile double delta;
  double value = 1.0;
  if (ratio->denominator == 0)
    return value;
  numerator = ratio->numerator_at_first +
      ratio->numerator_step *
        ((int32_t)bits - (int32_t)band->first_bits);
  delta = ratio->amount * (double)numerator / ratio->denominator;
  value += delta;
  return value;
}

static const siqs_policy_band_t *siqs_policy_band(uint32_t bits) {
  uint32_t index;
  for (index = 0; index < SIQS_POLICY_BAND_COUNT; index++) {
    const siqs_policy_band_t *band = &siqs_policy_bands[index];
    if (bits >= band->first_bits && bits <= band->last_bits)
      return band;
  }
  croak("SIQS: no bit-size policy for %u bits", bits);
  return NULL;
}


static void siqs_resolve_policy(siqs_policy_t *policy, uint32_t bits) {
  const siqs_policy_band_t *band;
  memset(policy, 0, sizeof(*policy));
  band = siqs_policy_band(bits);
  policy->band = band;
  policy->max_large_primes = band->production_large_primes;
  policy->q_count = band->q_count;
  policy->stage1_bias = band->stage1_bias_base;
  if (band->stage1_bias_step_bits != 0)
    policy->stage1_bias += (uint8_t)(
        (bits - band->stage1_bias_origin_bits) /
        band->stage1_bias_step_bits);
  policy->one_lp_multiplier = policy->max_large_primes == 1
                            ? band->one_lp_multiplier : 0U;
  policy->lp_multiplier_floor = policy->max_large_primes == 2
                              ? band->two_lp_multiplier_floor : 0U;
  policy->lp_product_floor = policy->max_large_primes == 2
                           ? band->two_lp_product_floor : 0U;
  policy->sieve_free_units = band->sieve_free_units;
  policy->sieve_start_prime_floor = band->sieve_start_prime_floor;
  policy->fb_floor = band->fb_floor;
  policy->relation_extra = band->relation_extra;
  policy->fb_coefficient =
      siqs_policy_linear_value(&band->fb_coefficient, bits);
  policy->interval_scale =
      siqs_policy_linear_value(&band->interval_base_scale, bits) *
      siqs_policy_ratio_value(band, bits);
  policy->smooth_bound_exponent = band->smooth_bound_exponent;
  policy->sieve_hit_exponent =
      siqs_policy_linear_value(&band->sieve_hit_exponent, bits);
}

static void siqs_select_parameters(siqs_parameters_t *p, const mpz_t n) {
  siqs_policy_t policy;
  double ln_n, ln_term, fb, interval, smooth;
  double sieve_hit_residual_multiplier;
  /* n is the post-trial cofactor handed to SIQS.  Size policies deliberately
   * use this N rather than the later multiplier product kN. */
  p->bits = (uint32_t)mpz_sizeinbase(n, 2);
  siqs_resolve_policy(&policy, p->bits);
  p->policy_name = policy.band->name;
  p->policy_first_bits = policy.band->first_bits;
  p->policy_last_bits = policy.band->last_bits;
  p->max_large_primes = policy.max_large_primes;
  p->one_lp_policy_multiplier = policy.one_lp_multiplier;
  p->lp_policy_multiplier_floor = policy.lp_multiplier_floor;
  p->lp_policy_product_multiplier_floor = policy.lp_product_floor;
  p->q_count = policy.q_count;
  p->sieve_free_units = policy.sieve_free_units;
  p->sieve_start_prime_floor = policy.sieve_start_prime_floor;
  p->fb_floor = policy.fb_floor;
  p->relation_extra = policy.relation_extra;
  p->fb_coefficient = policy.fb_coefficient;
  p->interval_scale = policy.interval_scale;
  p->smooth_bound_exponent = policy.smooth_bound_exponent;
  p->sieve_hit_exponent = policy.sieve_hit_exponent;
  p->stage1_bias = policy.stage1_bias;
  ln_n = p->bits * M_LN2;
  ln_term = sqrt(ln_n * log(ln_n));

  fb = exp(p->fb_coefficient * ln_term);
  if (fb < (double)p->fb_floor) fb = (double)p->fb_floor;
  p->fb_size = (uint32_t)fb;

  interval = 6144.0 + exp(0.37 * ln_term);
  interval *= p->interval_scale;
  p->half_interval = (uint32_t)interval;
  p->half_interval = (p->half_interval + SIQS_SIEVE_ALIGN - 1)
                   & ~(SIQS_SIEVE_ALIGN - 1);
  if (p->half_interval < 4096)
    p->half_interval = 4096;

  p->target_relations = p->fb_size + 1 + p->relation_extra;

  /* One-LP mode uses N^0.12 as its individual residual limit.  Two-LP mode
   * uses N^0.16 so the product bound admits two factor-base-sized primes.
   * After the factor base is built, the early 2LP policy may raise K and R
   * to their measured floors; their upper tails remain automatic. */
  smooth = exp(p->smooth_bound_exponent * ln_n);
  p->smooth_bound = (uint64_t)smooth;
  if (p->smooth_bound < UINT64_C(1000000))
    p->smooth_bound = UINT64_C(1000000);

  sieve_hit_residual_multiplier = 2.0;
  p->sieve_hit_residual_multiplier = sieve_hit_residual_multiplier;
  p->sieve_hit_bound_nominal = exp(p->sieve_hit_exponent * ln_n);
  p->sieve_hit_bound = p->sieve_hit_bound_nominal;
  if (p->sieve_hit_bound < sieve_hit_residual_multiplier *
                           (double)p->smooth_bound)
    p->sieve_hit_bound = sieve_hit_residual_multiplier *
                         (double)p->smooth_bound;
  p->log_scale = 0.0; /* Set after the factor base and kN are known. */
}

static int siqs_squarefree_small(uint32_t n) {
  uint32_t p;
  for (p = 3; p * p <= n; p += 2)
    if (n % (p * p) == 0)
      return 0;
  return 1;
}

static unsigned long siqs_choose_multiplier(const mpz_t n,
                                            uint32_t fb_size) {
  uint32_t kval[SIQS_MULTIPLIER_CAPACITY];
  uint32_t accepted[SIQS_MULTIPLIER_CAPACITY] = { 0 };
  uint32_t bits = (uint32_t)mpz_sizeinbase(n, 2);
  uint32_t max_multiplier = bits <= SIQS_MULTIPLIER_SHALLOW_LAST_BITS
                          ? SIQS_MULTIPLIER_MAX
                          : SIQS_MULTIPLIER_FULL_MAX;
  uint32_t kcount = 0, k, i, unfinished, wanted;
#if SIQS_MULTIPLIER_SEARCH_DIVISOR > 0
  if (bits <= SIQS_MULTIPLIER_SHALLOW_LAST_BITS)
    wanted = fb_size / SIQS_MULTIPLIER_SEARCH_DIVISOR;
  else
    wanted = fb_size < SIQS_MULTIPLIER_SEARCH_DEPTH
           ? fb_size : SIQS_MULTIPLIER_SEARCH_DEPTH;
#else
  wanted = fb_size < SIQS_MULTIPLIER_SEARCH_DEPTH
         ? fb_size : SIQS_MULTIPLIER_SEARCH_DEPTH;
#endif
  double score[SIQS_MULTIPLIER_CAPACITY], best_score;
  const uint32_t d2_bonus_sixteenths =
      (uint32_t)SIQS_MULTIPLIER_D2_BONUS_SIXTEENTHS;
  unsigned long best;
  PRIME_ITERATOR(iter);

  if (wanted < SIQS_MULTIPLIER_SEARCH_FLOOR)
    wanted = SIQS_MULTIPLIER_SEARCH_FLOOR;
  if (wanted > fb_size)
    wanted = fb_size;

  for (k = 1; k <= max_multiplier; k += 2) {
    unsigned long mod8;
    if (!siqs_squarefree_small(k))
      continue;
    kval[kcount] = k;
    mod8 = (mpz_fdiv_ui(n, 8) * k) & 7UL;
    score[kcount] = -SIQS_MULTIPLIER_SIZE_PENALTY * log((double)k);
    score[kcount] += mod8 == 1 ? 2.0 * M_LN2
                     : mod8 == 5 ? M_LN2 : 0.5 * M_LN2;
    /* Sweep a small extra preference for the residue class that permits the
     * d=2 polynomial.  Integer sixteenths express the planned fractions of
     * ln2 exactly at compile time. */
    if (mod8 == 1)
      score[kcount] += d2_bonus_sixteenths * M_LN2 / 16.0;
    kcount++;
  }

  unfinished = kcount;
  prime_iterator_setprime(&iter, 2);
  while (unfinished != 0) {
    uint32_t p = (uint32_t)prime_iterator_next(&iter);
    uint32_t nmod = (uint32_t)mpz_fdiv_ui(n, p);
    int nsymbol = nmod == 0 ? 0 : siqs_jacobi_odd_u32(nmod, p);
    double logp = log((double)p);
    for (i = 0; i < kcount; i++) {
      uint32_t km = kval[i] % p;
      if (accepted[i] == wanted)
        continue;
      if (km == 0) {
        score[i] += logp / p;
        accepted[i]++;
      } else if (nsymbol != 0 &&
                 siqs_jacobi_odd_u32(km, p) == nsymbol) {
        score[i] += 2.0 * logp / (p - 1);
        accepted[i]++;
      }
      if (accepted[i] == wanted)
        unfinished--;
    }
  }
  prime_iterator_destroy(&iter);

  best = kval[0];
  best_score = score[0];
  for (i = 1; i < kcount; i++)
    if (score[i] > best_score) {
      best_score = score[i];
      best = kval[i];
    }
  return best;
}

/* Returns zero when a factor of N was encountered during base construction. */
static int siqs_build_factor_base(siqs_ctx_t *ctx) {
  uint32_t count = 0;
  PRIME_ITERATOR(iter);
  mpz_t divisor;
  mpz_init(divisor);
  ctx->fb = (siqs_fb_t *)siqs_calloc(ctx->params.fb_size,
                                     sizeof(siqs_fb_t));

  ctx->fb[count].p = 2;
  ctx->fb[count].sqrt_kn = (uint32_t)mpz_fdiv_ui(ctx->kn, 2);
  count++;
  prime_iterator_setprime(&iter, 2);
  while (count < ctx->params.fb_size) {
    uint32_t p = (uint32_t)prime_iterator_next(&iter);
    uint32_t r = (uint32_t)mpz_fdiv_ui(ctx->kn, p);
    int symbol = siqs_legendre_u32(r, p);
    if (symbol < 0)
      continue;
    if (r == 0 && mpz_divisible_ui_p(ctx->n, p)) {
      mpz_set_ui(divisor, p);
      if (siqs_insert_divisor(ctx->result, divisor))
        ctx->factor_found = 1;
      prime_iterator_destroy(&iter);
      mpz_clear(divisor);
      return 0;
    }
    ctx->fb[count].p = p;
    ctx->fb[count].sqrt_kn = r == 0 ? 0 : siqs_sqrtmod_u32(r, p);
    count++;
  }
  prime_iterator_destroy(&iter);
  mpz_clear(divisor);
  ctx->largest_fb_prime = ctx->fb[ctx->params.fb_size - 1].p;
  return 1;
}

/*----------------------------------------------------------------------------
 * Relations and the large-prime graph
 *----------------------------------------------------------------------------*/

static int siqs_u32_cmp(const void *va, const void *vb) {
  uint32_t a = *(const uint32_t *)va;
  uint32_t b = *(const uint32_t *)vb;
  return a < b ? -1 : a > b ? 1 : 0;
}

static siqs_raw_relation_t *siqs_raw_relation_new(const mpz_t y,
    const siqs_factor_t *factors, uint32_t nfactors,
    uint64_t lp1, uint64_t lp2) {
  siqs_raw_relation_t *r =
      (siqs_raw_relation_t *)siqs_malloc(sizeof(*r));
  mpz_init_set(r->y, y);
  r->factors = (siqs_factor_t *)siqs_malloc(
      (size_t)nfactors * sizeof(*r->factors));
  memcpy(r->factors, factors, (size_t)nfactors * sizeof(*r->factors));
  r->nfactors = nfactors;
  if (lp1 > lp2) {
    uint64_t t = lp1;
    lp1 = lp2;
    lp2 = t;
  }
  r->lp1 = lp1;
  r->lp2 = lp2;
  r->fingerprint = 0;
  return r;
}

static void siqs_raw_relation_free(siqs_raw_relation_t *r) {
  if (r == NULL)
    return;
  mpz_clear(r->y);
  free(r->factors);
  free(r);
}

static void siqs_full_relation_free(siqs_full_relation_t *r) {
  if (r == NULL)
    return;
  mpz_clear(r->y);
  free(r->factors);
  free(r);
}

static uint64_t siqs_relation_fingerprint(const siqs_raw_relation_t *r) {
  uint64_t h;
  uint32_t i;
  h = (uint64_t)mpz_fdiv_ui(r->y, 4294967291UL);
  h ^= (uint64_t)mpz_fdiv_ui(r->y, 4294967279UL) << 32;
  h ^= siqs_mix64(r->lp1 + UINT64_C(0x243f6a8885a308d3));
  h ^= siqs_mix64(r->lp2 + UINT64_C(0x13198a2e03707344));
  for (i = 0; i < r->nfactors; i++)
    h = siqs_mix64(h ^ ((uint64_t)r->factors[i].row << 32)
                     ^ r->factors[i].exponent);
  return h;
}

#ifdef SIQS_DEBUG
static void siqs_verify_raw_relation(const siqs_ctx_t *ctx,
                                     const siqs_raw_relation_t *r) {
  uint32_t i;
  mpz_t lhs, rhs, t;
  mpz_init(lhs);
  mpz_init_set_ui(rhs, 1);
  mpz_init(t);
  mpz_mul(lhs, r->y, r->y);
  mpz_mod(lhs, lhs, ctx->n);
  for (i = 0; i < r->nfactors; i++) {
    uint32_t row = r->factors[i].row;
    if (row == 0) {
      if (r->factors[i].exponent & 1U)
        mpz_neg(rhs, rhs);
    } else {
      mpz_set_ui(t, ctx->fb[row - 1].p);
      mpz_powm_ui(t, t, r->factors[i].exponent, ctx->n);
      mpz_mul(rhs, rhs, t);
      mpz_mod(rhs, rhs, ctx->n);
    }
  }
  if (r->lp1 != 1) {
    siqs_mpz_set_u64(t, r->lp1);
    mpz_mul(rhs, rhs, t);
    mpz_mod(rhs, rhs, ctx->n);
  }
  if (r->lp2 != 1) {
    siqs_mpz_set_u64(t, r->lp2);
    mpz_mul(rhs, rhs, t);
    mpz_mod(rhs, rhs, ctx->n);
  }
  mpz_mod(rhs, rhs, ctx->n);
  if (mpz_cmp(lhs, rhs) != 0)
    croak("SIQS: raw relation invariant failed");
  mpz_clear(lhs);
  mpz_clear(rhs);
  mpz_clear(t);
}
#endif

static uint32_t siqs_store_raw(siqs_ctx_t *ctx, siqs_raw_relation_t *r) {
  if (ctx->raw_count == ctx->raw_alloc) {
    ctx->raw_alloc = ctx->raw_alloc ? ctx->raw_alloc * 2 : 1024;
    ctx->raw = (siqs_raw_relation_t **)siqs_realloc(
        ctx->raw, (size_t)ctx->raw_alloc * sizeof(*ctx->raw));
  }
  ctx->raw[ctx->raw_count] = r;
  return ctx->raw_count++;
}

static void siqs_store_full(siqs_ctx_t *ctx, siqs_full_relation_t *r) {
  if (ctx->full_count == ctx->full_alloc) {
    ctx->full_alloc = ctx->full_alloc ? ctx->full_alloc * 2 : 1024;
    ctx->full = (siqs_full_relation_t **)siqs_realloc(
        ctx->full, (size_t)ctx->full_alloc * sizeof(*ctx->full));
  }
  ctx->full[ctx->full_count++] = r;
}

/* A smooth raw relation is already a factor-base relation.  Move its GMP and
 * factor-vector ownership into the full-relation array instead of retaining a
 * second raw copy. */
static void siqs_materialize_smooth_relation(siqs_ctx_t *ctx,
                                              siqs_raw_relation_t *raw) {
  siqs_full_relation_t *full =
      (siqs_full_relation_t *)siqs_malloc(sizeof(*full));
  mpz_init(full->y);
  mpz_swap(full->y, raw->y);
  full->factors = raw->factors;
  full->nfactors = raw->nfactors;
  raw->factors = NULL;
  raw->nfactors = 0;
  siqs_raw_relation_free(raw);
  siqs_store_full(ctx, full);
}

/* Both input vectors are strictly row-ordered.  Merge them directly, adding
 * exponents for rows shared by the two partial relations. */
static siqs_factor_t *siqs_merge_relation_factors(
    const siqs_raw_relation_t *left, const siqs_raw_relation_t *right,
    uint32_t *merged_count) {
  uint32_t li = 0, ri = 0, count = 0, unique = 0;
  siqs_factor_t *merged;

  while (li < left->nfactors || ri < right->nfactors) {
    if (ri == right->nfactors ||
        (li < left->nfactors &&
         left->factors[li].row < right->factors[ri].row)) {
      li++;
    } else if (li == left->nfactors ||
               right->factors[ri].row < left->factors[li].row) {
      ri++;
    } else {
      li++;
      ri++;
    }
    unique++;
  }
  li = ri = 0;
  merged = (siqs_factor_t *)siqs_malloc(
      (size_t)unique * sizeof(*merged));

  while (li < left->nfactors || ri < right->nfactors) {
    if (ri == right->nfactors ||
        (li < left->nfactors &&
         left->factors[li].row < right->factors[ri].row)) {
      merged[count++] = left->factors[li++];
    } else if (li == left->nfactors ||
               right->factors[ri].row < left->factors[li].row) {
      merged[count++] = right->factors[ri++];
    } else {
      uint32_t exponent = left->factors[li].exponent
                        + right->factors[ri].exponent;
      if (exponent < left->factors[li].exponent)
        croak("SIQS: relation factor exponent overflow");
      merged[count].row = left->factors[li].row;
      merged[count].exponent = exponent;
      count++;
      li++;
      ri++;
    }
  }
  if (count != unique)
    croak("SIQS: relation factor merge count changed");
  *merged_count = count;
  return merged;
}

static uint32_t siqs_graph_hash(uint64_t value, uint32_t mask) {
  return (uint32_t)siqs_mix64(value) & mask;
}

static void siqs_one_lp_init(siqs_one_lp_state_t *state) {
  memset(state, 0, sizeof(*state));
  state->alloc = 1024;
  state->anchors = (siqs_one_lp_anchor_t *)siqs_calloc(
      state->alloc, sizeof(*state->anchors));
  mpz_init(state->y);
  mpz_init(state->large_prime);
  mpz_init(state->inverse);
  state->initialized = 1;
}

static void siqs_one_lp_clear(siqs_one_lp_state_t *state) {
  uint32_t i;
  if (!state->initialized)
    return;
  for (i = 0; i < state->alloc; i++)
    if (state->anchors[i].label != 0)
      siqs_raw_relation_free(state->anchors[i].relation);
  free(state->anchors);
  mpz_clear(state->y);
  mpz_clear(state->large_prime);
  mpz_clear(state->inverse);
  memset(state, 0, sizeof(*state));
}

static void siqs_one_lp_rehash(siqs_one_lp_state_t *state) {
  siqs_one_lp_anchor_t *old = state->anchors;
  uint32_t old_alloc = state->alloc, i;
  if (old_alloc > UINT32_MAX / 2)
    croak("SIQS: one-large-prime anchor table is too large");
  state->alloc *= 2;
  state->anchors = (siqs_one_lp_anchor_t *)siqs_calloc(
      state->alloc, sizeof(*state->anchors));
  for (i = 0; i < old_alloc; i++) {
    if (old[i].label != 0) {
      uint32_t pos = siqs_graph_hash(old[i].label, state->alloc - 1);
      while (state->anchors[pos].label != 0)
        pos = (pos + 1) & (state->alloc - 1);
      state->anchors[pos] = old[i];
    }
  }
  free(old);
}

/* Return the existing anchor for label, or retain relation as its first
 * anchor and return NULL. */
static siqs_raw_relation_t *siqs_one_lp_anchor_relation(
    siqs_one_lp_state_t *state, uint64_t label,
    siqs_raw_relation_t *relation) {
  uint32_t pos;
  if (label <= 1)
    croak("SIQS: invalid one-large-prime anchor label");
  if (state->count * 10 >= state->alloc * 7)
    siqs_one_lp_rehash(state);
  pos = siqs_graph_hash(label, state->alloc - 1);
  while (state->anchors[pos].label != 0) {
    if (state->anchors[pos].label == label)
      return state->anchors[pos].relation;
    pos = (pos + 1) & (state->alloc - 1);
  }
  state->anchors[pos].label = label;
  state->anchors[pos].relation = relation;
  state->count++;
  return NULL;
}

/* Materialize the two-edge cycle (1,p),(1,p).  Dividing the product of the
 * two square roots by p removes the repeated large prime.  A failed inverse
 * is retained as the same factor-discovery opportunity as the general graph.
 */
static int siqs_materialize_one_lp_pair(siqs_ctx_t *ctx,
    const siqs_raw_relation_t *anchor,
    const siqs_raw_relation_t *relation, uint64_t label) {
  siqs_one_lp_state_t *state = &ctx->one_lp;
  siqs_full_relation_t *full;

  siqs_mpz_set_u64(state->large_prime, label);
  if (!mpz_invert(state->inverse, state->large_prime, ctx->n)) {
    mpz_gcd(state->y, state->large_prime, ctx->n);
    if (mpz_cmp_ui(state->y, 1) > 0 && mpz_cmp(state->y, ctx->n) < 0 &&
        siqs_insert_divisor(ctx->result, state->y)) {
      ctx->factor_found = 1;
      siqs_verify_partition(ctx->original_n, ctx->result);
    }
    return 0;
  }

  mpz_mul(state->y, anchor->y, relation->y);
  mpz_mod(state->y, state->y, ctx->n);
  mpz_mul(state->y, state->y, state->inverse);
  mpz_mod(state->y, state->y, ctx->n);

  full = (siqs_full_relation_t *)siqs_malloc(sizeof(*full));
  mpz_init_set(full->y, state->y);
  full->factors = siqs_merge_relation_factors(
      anchor, relation, &full->nfactors);
  siqs_store_full(ctx, full);
  return 1;
}

static void siqs_graph_rehash(siqs_graph_t *g, uint32_t new_size) {
  uint32_t *old = g->table;
  uint32_t i;
  g->table = (uint32_t *)siqs_calloc(new_size, sizeof(uint32_t));
  g->table_alloc = new_size;
  for (i = 0; i < g->vertex_count; i++) {
    uint32_t pos = siqs_graph_hash(g->vertices[i].value, new_size - 1);
    while (g->table[pos] != 0)
      pos = (pos + 1) & (new_size - 1);
    g->table[pos] = i + 1;
  }
  free(old);
}

static void siqs_graph_init(siqs_graph_t *g) {
  memset(g, 0, sizeof(*g));
  g->vertex_alloc = 1024;
  g->vertices = (siqs_vertex_t *)siqs_malloc(
      (size_t)g->vertex_alloc * sizeof(*g->vertices));
  g->table_alloc = 2048;
  g->table = (uint32_t *)siqs_calloc(g->table_alloc, sizeof(uint32_t));
  /* Vertex zero is the distinguished endpoint for one-large-prime edges. */
  g->vertices[0].value = 1;
  g->vertices[0].dsu_parent = 0;
  g->vertices[0].dsu_size = 1;
  g->vertices[0].tree_parent = SIQS_NO_INDEX;
  g->vertices[0].parent_edge = SIQS_NO_INDEX;
  g->vertices[0].mark = 0;
  g->vertex_count = 1;
  g->table[siqs_graph_hash(1, g->table_alloc - 1)] = 1;
  g->mark = 0;
}

static void siqs_graph_clear(siqs_graph_t *g) {
  free(g->vertices);
  free(g->table);
  memset(g, 0, sizeof(*g));
}

static uint32_t siqs_graph_vertex(siqs_graph_t *g, uint64_t value) {
  uint32_t pos;
  if (g->vertex_count * 10 >= g->table_alloc * 7)
    siqs_graph_rehash(g, g->table_alloc * 2);
  pos = siqs_graph_hash(value, g->table_alloc - 1);
  while (g->table[pos] != 0) {
    uint32_t id = g->table[pos] - 1;
    if (g->vertices[id].value == value)
      return id;
    pos = (pos + 1) & (g->table_alloc - 1);
  }
  if (g->vertex_count == g->vertex_alloc) {
    g->vertex_alloc *= 2;
    g->vertices = (siqs_vertex_t *)siqs_realloc(
        g->vertices, (size_t)g->vertex_alloc * sizeof(*g->vertices));
  }
  {
    uint32_t id = g->vertex_count++;
    siqs_vertex_t *v = &g->vertices[id];
    v->value = value;
    v->dsu_parent = id;
    v->dsu_size = 1;
    v->tree_parent = SIQS_NO_INDEX;
    v->parent_edge = SIQS_NO_INDEX;
    v->mark = 0;
    g->table[pos] = id + 1;
    return id;
  }
}

static uint32_t siqs_graph_find(siqs_graph_t *g, uint32_t v) {
  uint32_t root = v;
  while (g->vertices[root].dsu_parent != root)
    root = g->vertices[root].dsu_parent;
  while (g->vertices[v].dsu_parent != v) {
    uint32_t next = g->vertices[v].dsu_parent;
    g->vertices[v].dsu_parent = root;
    v = next;
  }
  return root;
}

/* Reverse parent pointers so that endpoint becomes the root of its tree. */
static void siqs_graph_reroot(siqs_graph_t *g, uint32_t endpoint) {
  uint32_t previous = SIQS_NO_INDEX;
  uint32_t previous_edge = SIQS_NO_INDEX;
  uint32_t current = endpoint;
  while (current != SIQS_NO_INDEX) {
    uint32_t next = g->vertices[current].tree_parent;
    uint32_t next_edge = g->vertices[current].parent_edge;
    g->vertices[current].tree_parent = previous;
    g->vertices[current].parent_edge = previous_edge;
    previous = current;
    previous_edge = next_edge;
    current = next;
  }
}

static void siqs_touch_factor(siqs_ctx_t *ctx, uint32_t row,
                              uint32_t exponent) {
  if (ctx->factor_counts[row] == 0)
    ctx->factor_touched[ctx->factor_touched_count++] = row;
  ctx->factor_counts[row] += exponent;
}

static void siqs_reset_touched_factors(siqs_ctx_t *ctx) {
  uint32_t i;
  for (i = 0; i < ctx->factor_touched_count; i++)
    ctx->factor_counts[ctx->factor_touched[i]] = 0;
  ctx->factor_touched_count = 0;
}

/* Convert a fundamental cycle of raw partials into a factor-base relation. */
static void siqs_materialize_cycle(siqs_ctx_t *ctx,
    const uint32_t *edges, uint32_t edge_count,
    const uint32_t *vertices, uint32_t vertex_count) {
  uint32_t i, j, nfactors;
  siqs_full_relation_t *full;
  mpz_t y, lp_product, lp_value, inverse, gcd;

  mpz_init_set_ui(y, 1);
  mpz_init_set_ui(lp_product, 1);
  mpz_init(lp_value);
  mpz_init(inverse);
  mpz_init(gcd);
  siqs_reset_touched_factors(ctx);

  for (i = 0; i < edge_count; i++) {
    const siqs_raw_relation_t *r = ctx->raw[edges[i]];
    mpz_mul(y, y, r->y);
    mpz_mod(y, y, ctx->n);
    for (j = 0; j < r->nfactors; j++)
      siqs_touch_factor(ctx, r->factors[j].row,
                       r->factors[j].exponent);
  }
  for (i = 0; i < vertex_count; i++) {
    uint64_t lp = ctx->graph.vertices[vertices[i]].value;
    if (lp != 1) {
      siqs_mpz_set_u64(lp_value, lp);
      mpz_mul(lp_product, lp_product, lp_value);
      mpz_mod(lp_product, lp_product, ctx->n);
    }
  }
  if (!mpz_invert(inverse, lp_product, ctx->n)) {
    mpz_gcd(gcd, lp_product, ctx->n);
    if (mpz_cmp_ui(gcd, 1) > 0 && mpz_cmp(gcd, ctx->n) < 0 &&
        siqs_insert_divisor(ctx->result, gcd)) {
      ctx->factor_found = 1;
      siqs_verify_partition(ctx->original_n, ctx->result);
    }
    siqs_reset_touched_factors(ctx);
    mpz_clear(y);
    mpz_clear(lp_product);
    mpz_clear(lp_value);
    mpz_clear(inverse);
    mpz_clear(gcd);
    return;
  }
  mpz_mul(y, y, inverse);
  mpz_mod(y, y, ctx->n);

  qsort(ctx->factor_touched, ctx->factor_touched_count,
        sizeof(uint32_t), siqs_u32_cmp);
  nfactors = ctx->factor_touched_count;
  full = (siqs_full_relation_t *)siqs_malloc(sizeof(*full));
  mpz_init_set(full->y, y);
  full->nfactors = nfactors;
  full->factors = (siqs_factor_t *)siqs_malloc(
      (size_t)nfactors * sizeof(*full->factors));
  for (i = 0; i < nfactors; i++) {
    uint32_t row = ctx->factor_touched[i];
    full->factors[i].row = row;
    full->factors[i].exponent = ctx->factor_counts[row];
  }
  siqs_store_full(ctx, full);
  siqs_reset_touched_factors(ctx);
  mpz_clear(y);
  mpz_clear(lp_product);
  mpz_clear(lp_value);
  mpz_clear(inverse);
  mpz_clear(gcd);
}

/* Add one graph edge.  A cycle is materialized immediately. */
static void siqs_graph_add_relation(siqs_ctx_t *ctx, uint32_t edge) {
  siqs_raw_relation_t *rel = ctx->raw[edge];
  siqs_graph_t *g = &ctx->graph;
  uint32_t u = siqs_graph_vertex(g, rel->lp1);
  uint32_t v = siqs_graph_vertex(g, rel->lp2);
  uint32_t ru = siqs_graph_find(g, u);
  uint32_t rv = siqs_graph_find(g, v);

  if (ru != rv) {
    /* Re-root only the smaller component; total re-rooting work is O(n log n). */
    if (g->vertices[ru].dsu_size > g->vertices[rv].dsu_size) {
      uint32_t t = u; u = v; v = t;
      t = ru; ru = rv; rv = t;
    }
    siqs_graph_reroot(g, u);
    g->vertices[u].tree_parent = v;
    g->vertices[u].parent_edge = edge;
    g->vertices[ru].dsu_parent = rv;
    g->vertices[rv].dsu_size += g->vertices[ru].dsu_size;
    return;
  }

  {
    uint32_t current, lca, ec = 0, vc = 0, ealloc = 64, valloc = 64;
    uint32_t *edges = (uint32_t *)siqs_malloc(
        (size_t)ealloc * sizeof(uint32_t));
    uint32_t *vertices = (uint32_t *)siqs_malloc(
        (size_t)valloc * sizeof(uint32_t));

    if (++g->mark == 0) {
      uint32_t i;
      for (i = 0; i < g->vertex_count; i++)
        g->vertices[i].mark = 0;
      g->mark = 1;
    }
    current = u;
    while (current != SIQS_NO_INDEX) {
      g->vertices[current].mark = g->mark;
      current = g->vertices[current].tree_parent;
    }
    current = v;
    while (g->vertices[current].mark != g->mark)
      current = g->vertices[current].tree_parent;
    lca = current;

    edges[ec++] = edge;
    current = u;
    vertices[vc++] = current;
    while (current != lca) {
      if (ec == ealloc) {
        ealloc *= 2;
        edges = (uint32_t *)siqs_realloc(
            edges, (size_t)ealloc * sizeof(uint32_t));
      }
      edges[ec++] = g->vertices[current].parent_edge;
      current = g->vertices[current].tree_parent;
      if (vc == valloc) {
        valloc *= 2;
        vertices = (uint32_t *)siqs_realloc(
            vertices, (size_t)valloc * sizeof(uint32_t));
      }
      vertices[vc++] = current;
    }
    current = v;
    while (current != lca) {
      if (vc == valloc) {
        valloc *= 2;
        vertices = (uint32_t *)siqs_realloc(
            vertices, (size_t)valloc * sizeof(uint32_t));
      }
      vertices[vc++] = current;
      if (ec == ealloc) {
        ealloc *= 2;
        edges = (uint32_t *)siqs_realloc(
            edges, (size_t)ealloc * sizeof(uint32_t));
      }
      edges[ec++] = g->vertices[current].parent_edge;
      current = g->vertices[current].tree_parent;
    }
    siqs_materialize_cycle(ctx, edges, ec, vertices, vc);
    free(edges);
    free(vertices);
  }
}

/* In one-large-prime mode every partial edge is (1,p).  The first edge for p
 * is its permanent anchor and each later edge closes the same two-edge
 * fundamental cycle that the general graph would materialize. */
static void siqs_accept_one_lp_relation(siqs_ctx_t *ctx,
                                        siqs_raw_relation_t *relation) {
  siqs_raw_relation_t *anchor;

  if (relation->lp1 == 1 && relation->lp2 == 1) {
    siqs_materialize_smooth_relation(ctx, relation);
    return;
  }
  if (relation->lp1 != 1 || relation->lp2 == 1)
    croak("SIQS: one-large-prime path received an invalid relation");

  anchor = siqs_one_lp_anchor_relation(
      &ctx->one_lp, relation->lp2, relation);
  if (anchor == NULL) {
    return;
  }

  (void)siqs_materialize_one_lp_pair(
      ctx, anchor, relation, relation->lp2);
  siqs_raw_relation_free(relation);
}

static void siqs_accept_raw_relation(siqs_ctx_t *ctx,
                                     siqs_raw_relation_t *r) {
  r->fingerprint = siqs_relation_fingerprint(r);
  if (!siqs_hashset_insert(&ctx->relation_hashes, r->fingerprint)) {
    siqs_raw_relation_free(r);
    return;
  }
#ifdef SIQS_DEBUG
  siqs_verify_raw_relation(ctx, r);
#endif
  if (r->lp1 == 1 && r->lp2 == 1)
    ctx->accepted_smooth++;
  else if (r->lp1 == 1)
    ctx->accepted_one_lp++;
  else
    ctx->accepted_two_lp++;
  if (siqs_use_one_lp_relation_path(ctx)) {
    siqs_accept_one_lp_relation(ctx, r);
    return;
  }
  {
    uint32_t edge = siqs_store_raw(ctx, r);
    siqs_graph_add_relation(ctx, edge);
  }
}

/*----------------------------------------------------------------------------
 * Polynomial families
 *----------------------------------------------------------------------------*/

static uint32_t siqs_nearest_fb_index(const siqs_ctx_t *ctx,
                                      uint32_t wanted) {
  uint32_t lo = 1, hi = ctx->params.fb_size;
  while (lo < hi) {
    uint32_t mid = lo + (hi - lo) / 2;
    if (ctx->fb[mid].p < wanted)
      lo = mid + 1;
    else
      hi = mid;
  }
  if (lo >= ctx->params.fb_size)
    return ctx->params.fb_size - 1;
  if (lo > 1 && wanted - ctx->fb[lo - 1].p <= ctx->fb[lo].p - wanted)
    return lo - 1;
  return lo;
}

static uint32_t siqs_nearest_available_fb(siqs_ctx_t *ctx,
                                          uint32_t wanted) {
  uint32_t center = siqs_nearest_fb_index(ctx, wanted);
  uint32_t step;
  for (step = 0; step < ctx->params.fb_size; step++) {
    uint32_t up = center + step;
    if (up < ctx->params.fb_size && !ctx->fb[up].in_a &&
        ctx->fb[up].sqrt_kn != 0)
      return up;
    if (step != 0 && center >= step) {
      uint32_t down = center - step;
      if (down > 0 && !ctx->fb[down].in_a &&
          ctx->fb[down].sqrt_kn != 0)
        return down;
    }
  }
  croak("SIQS: factor base has no available polynomial prime");
  return 0;
}

static uint64_t siqs_a_fingerprint(const uint32_t *indices, uint32_t count) {
  uint32_t i;
  uint64_t h = UINT64_C(0x6a09e667f3bcc909);
  for (i = 0; i < count; i++)
    h = siqs_mix64(h ^ ((uint64_t)indices[i] +
                        UINT64_C(0x9e3779b97f4a7c15)));
  return h;
}

static void siqs_poly_init(siqs_ctx_t *ctx, siqs_poly_t *poly) {
  uint32_t i;
  memset(poly, 0, sizeof(*poly));
  poly->q_count = ctx->params.q_count;
  poly->b_limit = 1U << (poly->q_count - 1);
  poly->a_index = (uint32_t *)siqs_malloc(
      (size_t)poly->q_count * sizeof(uint32_t));
  poly->H = (mpz_t *)siqs_malloc((size_t)poly->q_count * sizeof(mpz_t));
  for (i = 0; i < poly->q_count; i++) {
    poly->a_index[i] = SIQS_NO_INDEX;
    mpz_init(poly->H[i]);
  }
  poly->corrections = (uint32_t *)siqs_malloc(
      (size_t)(poly->q_count - 1) * ctx->params.fb_size * sizeof(uint32_t));
  mpz_init(poly->A);
  mpz_init(poly->DA);
  mpz_init(poly->B);
  mpz_init(poly->C);
  mpz_init(poly->target_A);
  mpz_mul_ui(poly->target_A, ctx->kn, 2);
  mpz_sqrt(poly->target_A, poly->target_A);
  mpz_fdiv_q_ui(poly->target_A, poly->target_A,
                ctx->params.half_interval);
  mpz_fdiv_q_ui(poly->target_A, poly->target_A, ctx->params.poly_d);
}

static void siqs_poly_clear(siqs_ctx_t *ctx, siqs_poly_t *poly) {
  uint32_t i;
  for (i = 0; i < poly->q_count; i++) {
    if (poly->a_index[i] < ctx->params.fb_size)
      ctx->fb[poly->a_index[i]].in_a = 0;
    mpz_clear(poly->H[i]);
  }
  free(poly->H);
  free(poly->a_index);
  free(poly->corrections);
  mpz_clear(poly->A);
  mpz_clear(poly->DA);
  mpz_clear(poly->B);
  mpz_clear(poly->C);
  mpz_clear(poly->target_A);
  memset(poly, 0, sizeof(*poly));
}

static int siqs_choose_A(siqs_ctx_t *ctx, siqs_poly_t *poly) {
  uint32_t attempt, i;
  uint32_t attempt_limit = poly->q_count <= 3 ? 30000U : 10000U;
  double ideal_d = pow(mpz_get_d(poly->target_A),
                       1.0 / poly->q_count);
  uint32_t ideal = ideal_d < 3.0 ? 3U : (uint32_t)ideal_d;
  uint32_t center = siqs_nearest_fb_index(ctx, ideal);
  uint32_t variance = (uint32_t)(0.75 * sqrt((double)ctx->params.fb_size));
  mpz_t product, remaining, scaled;
  if (variance < 8)
    variance = 8;
  mpz_init(product);
  mpz_init(remaining);
  mpz_init(scaled);

  for (i = 0; i < poly->q_count; i++)
    if (poly->a_index[i] < ctx->params.fb_size)
      ctx->fb[poly->a_index[i]].in_a = 0;

  /* The short low-end polynomials have fewer distinct A products.  Preserve
   * the normal near-optimal search, then widen only an otherwise exhausted
   * q=2 or q=3 family rather than failing an otherwise healthy small input. */
  for (attempt = 0; attempt < attempt_limit; attempt++) {
    uint32_t tolerance = attempt < 10000U ? 2U
                       : attempt < 20000U ? 4U : 8U;
    uint64_t fingerprint;
    mpz_set_ui(product, 1);
    for (i = 0; i + 1 < poly->q_count; i++) {
      int64_t wanted_index;
      uint32_t wanted, index;
      if (attempt >= 20000U) {
        wanted_index = 1 + (int64_t)siqs_rand_range(
            &ctx->poly_rng, ctx->params.fb_size - 1);
      } else {
        uint32_t search_variance = attempt < 10000U
                                 ? variance : 2U * variance;
        int64_t offset = (int64_t)siqs_rand_range(
            &ctx->poly_rng, 2 * search_variance + 1) - search_variance;
        wanted_index = (int64_t)center + offset;
        if (wanted_index < 1)
          wanted_index = 1;
        if (wanted_index >= (int64_t)ctx->params.fb_size)
          wanted_index = ctx->params.fb_size - 1;
      }
      wanted = ctx->fb[wanted_index].p;
      index = siqs_nearest_available_fb(ctx, wanted);
      poly->a_index[i] = index;
      ctx->fb[index].in_a = 1;
      mpz_mul_ui(product, product, ctx->fb[index].p);
    }
    mpz_fdiv_q(remaining, poly->target_A, product);
    {
      unsigned long wanted = mpz_fits_ulong_p(remaining)
                           ? mpz_get_ui(remaining) : ULONG_MAX;
      uint32_t index = siqs_nearest_available_fb(
          ctx, wanted > UINT32_MAX ? UINT32_MAX : (uint32_t)wanted);
      poly->a_index[poly->q_count - 1] = index;
      ctx->fb[index].in_a = 1;
      mpz_mul_ui(product, product, ctx->fb[index].p);
    }

    qsort(poly->a_index, poly->q_count, sizeof(uint32_t), siqs_u32_cmp);
    mpz_mul_ui(scaled, product, tolerance);
    if (mpz_cmp(scaled, poly->target_A) >= 0) {
      mpz_mul_ui(scaled, poly->target_A, tolerance);
      if (mpz_cmp(product, scaled) <= 0) {
        fingerprint = siqs_a_fingerprint(poly->a_index, poly->q_count);
        if (siqs_hashset_insert(&ctx->a_hashes, fingerprint)) {
          mpz_set(poly->A, product);
          mpz_mul_ui(poly->DA, poly->A, ctx->params.poly_d);
          mpz_clear(product);
          mpz_clear(remaining);
          mpz_clear(scaled);
          return 1;
        }
      }
    }
    for (i = 0; i < poly->q_count; i++)
      ctx->fb[poly->a_index[i]].in_a = 0;
  }
  mpz_clear(product);
  mpz_clear(remaining);
  mpz_clear(scaled);
  return 0;
}

static void siqs_compute_C(siqs_ctx_t *ctx, siqs_poly_t *poly) {
  mpz_mul(poly->C, poly->B, poly->B);
  mpz_sub(poly->C, poly->C, ctx->kn);
  if (!mpz_divisible_p(poly->C, poly->DA))
    croak("SIQS: polynomial is not integral");
  mpz_divexact(poly->C, poly->C, poly->DA);
}

static void siqs_set_special_roots(siqs_ctx_t *ctx, siqs_poly_t *poly) {
  uint32_t i;
  uint32_t m = ctx->params.half_interval;
  /* With d=2, B is odd and B^2-kN is divisible by 8*A.  Thus
   * q(x) = 2*A*x^2 + 2*B*x + C is divisible by four for every x.
   * The ordinary d=1 polynomial has exactly one root modulo 2. */
  if (ctx->params.poly_d == 2) {
    ctx->root1[0] = 0;
    ctx->root2[0] = 1;
  } else {
    ctx->root1[0] = ((uint32_t)mpz_fdiv_ui(poly->C, 2) + (m & 1U)) & 1U;
    ctx->root2[0] = SIQS_NO_ROOT;
  }
  for (i = 0; i < poly->q_count; i++) {
    uint32_t index = poly->a_index[i];
    uint32_t p = ctx->fb[index].p;
    uint32_t b = (uint32_t)mpz_fdiv_ui(poly->B, p);
    uint32_t c = (uint32_t)mpz_fdiv_ui(poly->C, p);
    uint32_t denom = (uint32_t)((2ULL * b) % p);
    uint32_t inv = siqs_inverse_u32(denom, p);
    uint32_t x;
    if (inv == 0)
      croak("SIQS: singular polynomial root");
    x = (uint32_t)((uint64_t)(c ? p - c : 0) * inv % p);
    ctx->root1[index] = (x + m % p) % p;
    ctx->root2[index] = SIQS_NO_ROOT;
  }
}

#ifdef SIQS_DEBUG
static uint32_t siqs_debug_polynomial_at_root(const siqs_ctx_t *ctx,
                                              const siqs_poly_t *poly,
                                              uint32_t fb_index,
                                              uint32_t root) {
  uint32_t p = ctx->fb[fb_index].p;
  uint32_t mmod = ctx->params.half_interval % p;
  uint32_t x = root >= mmod ? root - mmod : root + p - mmod;
  uint32_t da = (uint32_t)mpz_fdiv_ui(poly->DA, p);
  uint32_t b = (uint32_t)mpz_fdiv_ui(poly->B, p);
  uint32_t c = (uint32_t)mpz_fdiv_ui(poly->C, p);
  uint32_t linear = (uint32_t)(((uint64_t)2 * b) % p);
  uint32_t value = (uint32_t)(((uint64_t)da * x + linear) % p);
  return (uint32_t)(((uint64_t)value * x + c) % p);
}

static void siqs_verify_polynomial(const siqs_ctx_t *ctx,
                                   const siqs_poly_t *poly) {
  uint32_t i;
  mpz_t product, expected_da, difference, quotient, modulus;
  mpz_init_set_ui(product, 1);
  mpz_init(expected_da);
  mpz_init(difference);
  mpz_init(quotient);
  mpz_init(modulus);

  for (i = 0; i < poly->q_count; i++)
    mpz_mul_ui(product, product, ctx->fb[poly->a_index[i]].p);
  if (mpz_cmp(product, poly->A) != 0)
    croak("SIQS: polynomial A does not match its selected primes");
  mpz_mul_ui(expected_da, poly->A, ctx->params.poly_d);
  if (mpz_cmp(expected_da, poly->DA) != 0)
    croak("SIQS: polynomial d*A coefficient is inconsistent");

  mpz_mul(difference, poly->B, poly->B);
  mpz_sub(difference, difference, ctx->kn);
  if (!mpz_divisible_p(difference, poly->DA))
    croak("SIQS: polynomial numerator is not divisible by d*A");
  mpz_divexact(quotient, difference, poly->DA);
  if (mpz_cmp(quotient, poly->C) != 0)
    croak("SIQS: polynomial constant is inconsistent");

  if (ctx->params.poly_d == 2) {
    if (mpz_fdiv_ui(ctx->kn, 8) != 1 || mpz_even_p(poly->B))
      croak("SIQS: invalid d=2 parity invariant");
    mpz_mul_ui(modulus, poly->A, 8);
    if (!mpz_divisible_p(difference, modulus))
      croak("SIQS: d=2 polynomial lacks the congruence modulo 8*A");
    if (mpz_fdiv_ui(poly->C, 4) != 0)
      croak("SIQS: d=2 normalized polynomial is not always divisible by four");
  }

  for (i = 0; i < ctx->params.fb_size; i++) {
    uint32_t p = ctx->fb[i].p;
    int expect_two = i == 0
                   ? ctx->params.poly_d == 2
                   : !ctx->fb[i].in_a && ctx->fb[i].sqrt_kn != 0;
    if (ctx->root1[i] >= p ||
        siqs_debug_polynomial_at_root(ctx, poly, i, ctx->root1[i]) != 0)
      croak("SIQS: first polynomial root invariant failed");
    if (expect_two) {
      if (ctx->root2[i] >= p || ctx->root2[i] == ctx->root1[i] ||
          siqs_debug_polynomial_at_root(ctx, poly, i, ctx->root2[i]) != 0)
        croak("SIQS: second polynomial root invariant failed");
    } else if (ctx->root2[i] != SIQS_NO_ROOT) {
      croak("SIQS: polynomial unexpectedly has a second root");
    }
  }

  mpz_clear(product);
  mpz_clear(expected_da);
  mpz_clear(difference);
  mpz_clear(quotient);
  mpz_clear(modulus);
}
#endif

static void siqs_first_B_and_roots(siqs_ctx_t *ctx, siqs_poly_t *poly) {
  uint32_t i, j;
  mpz_t adiv;
  mpz_init(adiv);
  mpz_set_ui(poly->B, 0);

  for (i = 0; i < poly->q_count; i++) {
    uint32_t index = poly->a_index[i];
    uint32_t p = ctx->fb[index].p;
    uint32_t amod, inv, gamma;
    mpz_divexact_ui(adiv, poly->A, p);
    amod = (uint32_t)mpz_fdiv_ui(adiv, p);
    inv = siqs_inverse_u32(amod, p);
    gamma = (uint32_t)((uint64_t)ctx->fb[index].sqrt_kn * inv % p);
    if (gamma > p / 2)
      gamma = p - gamma;
    mpz_mul_ui(poly->H[i], adiv, gamma);
    mpz_add(poly->B, poly->B, poly->H[i]);
  }
  /* For kN == 1 (mod 8), selecting the odd representative makes
   * B^2-kN divisible by 8*A.  Adding A changes no odd CRT residue, and every
   * Gray-code correction is even, so the stronger congruence persists for
   * the whole family. */
  if (ctx->params.poly_d == 2 && mpz_even_p(poly->B))
    mpz_add(poly->B, poly->B, poly->A);
  siqs_compute_C(ctx, poly);

  for (j = 1; j < ctx->params.fb_size; j++) {
    uint32_t p = ctx->fb[j].p;
    if (ctx->fb[j].in_a)
      continue;
    {
      uint32_t amod = (uint32_t)mpz_fdiv_ui(poly->DA, p);
      uint32_t inva = siqs_inverse_u32(amod, p);
      uint32_t bmod = (uint32_t)mpz_fdiv_ui(poly->B, p);
      uint32_t negb = bmod ? p - bmod : 0;
      uint32_t s = ctx->fb[j].sqrt_kn;
      uint32_t x1 = (uint32_t)((uint64_t)((negb + s) % p) * inva % p);
      uint32_t x2 = (uint32_t)((uint64_t)((negb + p - s) % p) * inva % p);
      ctx->root1[j] = (x1 + ctx->params.half_interval % p) % p;
      if (x1 == x2)
        ctx->root2[j] = SIQS_NO_ROOT;
      else
        ctx->root2[j] = (x2 + ctx->params.half_interval % p) % p;
      for (i = 0; i + 1 < poly->q_count; i++) {
        uint32_t hmod = (uint32_t)mpz_fdiv_ui(poly->H[i], p);
        poly->corrections[(size_t)i * ctx->params.fb_size + j] =
            (uint32_t)((uint64_t)((2ULL * hmod) % p) * inva % p);
      }
    }
  }
  siqs_set_special_roots(ctx, poly);
#ifdef SIQS_DEBUG
  siqs_verify_polynomial(ctx, poly);
#endif
  poly->b_index = 0;
  mpz_clear(adiv);
}

/* Set the expected score contribution of factor-base primes omitted from
 * the byte sieve.  Use the actual roots of this family: ordinary primes
 * have two, while 2, primes dividing kN, and primes in A have one.  Recompute
 * the sum so rounding cannot accumulate across A's. */
static void siqs_set_family_sieve_initial(siqs_ctx_t *ctx) {
  uint32_t i;
  uint32_t initial;
  uint32_t free_units = ctx->params.sieve_free_units;
  double expected = (double)free_units;
  /* Two residue classes modulo 2 account for one guaranteed factor 2 in the
   * usual root sum.  A d=2 normalized polynomial is divisible by 4, so add
   * the second guaranteed factor explicitly. */
  if (ctx->params.poly_d == 2)
    expected += ctx->fb[0].logp;
  for (i = 0; i < ctx->params.sieve_start; i++) {
    uint32_t roots = ctx->root2[i] == SIQS_NO_ROOT ? 1U : 2U;
    expected += roots * (double)ctx->fb[i].logp / ctx->fb[i].p;
  }
  initial = (uint32_t)(expected + 0.5);
  if (initial > 127)
    initial = 127;
  ctx->active_sieve_initial = (uint8_t)initial;
}

static int siqs_new_family(siqs_ctx_t *ctx, siqs_poly_t *poly) {
  uint32_t i;
  for (i = 0; i < poly->q_count; i++)
    if (poly->a_index[i] < ctx->params.fb_size)
      ctx->fb[poly->a_index[i]].in_a = 0;
  if (!siqs_choose_A(ctx, poly))
    return 0;
  siqs_first_B_and_roots(ctx, poly);
  siqs_set_family_sieve_initial(ctx);
  return 1;
}

static INLINE void siqs_update_roots(siqs_ctx_t *ctx,
                                     const siqs_poly_t *poly,
                                     uint32_t bit, int subtract,
                                     uint32_t first, uint32_t end) {
  const uint32_t *corrections = poly->corrections
                              + (size_t)bit * ctx->params.fb_size;
  const uint32_t *primes = ctx->prime;
  uint32_t *root1 = ctx->root1;
  uint32_t *root2 = ctx->root2;
  uint32_t j;
  if (subtract) {
#if defined(__clang__)
# pragma clang loop vectorize(enable)
#endif
    for (j = first; j < end; j++) {
      uint32_t p = primes[j];
      uint32_t corr = corrections[j];
      root1[j] += corr;
      if (root1[j] >= p) root1[j] -= p;
      if (root2[j] != SIQS_NO_ROOT) {
        root2[j] += corr;
        if (root2[j] >= p) root2[j] -= p;
      }
    }
  } else {
    for (j = first; j < end; j++) {
      uint32_t p = primes[j];
      uint32_t corr = corrections[j];
      root1[j] = root1[j] >= corr
               ? root1[j] - corr
               : root1[j] + p - corr;
      if (root2[j] != SIQS_NO_ROOT)
        root2[j] = root2[j] >= corr
                 ? root2[j] - corr
                 : root2[j] + p - corr;
    }
  }
}

static int siqs_next_B(siqs_ctx_t *ctx, siqs_poly_t *poly) {
  uint32_t bit, gray, i, first;
  int subtract;
  if (poly->b_index + 1 >= poly->b_limit)
    return 0;
  poly->b_index++;
  bit = siqs_ctz32(poly->b_index);
  gray = poly->b_index ^ (poly->b_index >> 1);
  subtract = (gray & (1U << bit)) != 0;
  if (subtract) {
    mpz_sub(poly->B, poly->B, poly->H[bit]);
    mpz_sub(poly->B, poly->B, poly->H[bit]);
  } else {
    mpz_add(poly->B, poly->B, poly->H[bit]);
    mpz_add(poly->B, poly->B, poly->H[bit]);
  }
  /* a_index is sorted.  Updating the ranges between its entries avoids an
   * in-A test for every factor-base entry on every polynomial. */
  first = 1;
  for (i = 0; i < poly->q_count; i++) {
    uint32_t a_index = poly->a_index[i];
    siqs_update_roots(ctx, poly, bit, subtract, first, a_index);
    first = a_index + 1;
  }
  siqs_update_roots(ctx, poly, bit, subtract, first,
                    ctx->params.fb_size);
  siqs_compute_C(ctx, poly);
  siqs_set_special_roots(ctx, poly);
#ifdef SIQS_DEBUG
  siqs_verify_polynomial(ctx, poly);
#endif
  return 1;
}

/*----------------------------------------------------------------------------
 * Log sieve, resieving, and relation construction
 *----------------------------------------------------------------------------*/

static double siqs_log_mpz(const mpz_t n) {
  signed long exponent;
  double mantissa = mpz_get_d_2exp(&exponent, n);
  return log(fabs(mantissa)) + (double)exponent * M_LN2;
}

static void siqs_set_log_weights(siqs_ctx_t *ctx) {
  uint32_t i;
  double root_weighted_unsieved = 0.0;
  double log_q = 0.5 * siqs_log_mpz(ctx->kn)
               + log((double)ctx->params.half_interval)
               + 0.5 * M_LN2;
  double denom = log_q - log(ctx->params.sieve_hit_bound);
  double scale;
  uint32_t free_units = ctx->params.sieve_free_units;
  if (denom < 8.0)
    denom = 8.0;
  /* The fixed free units are added back in the initializer, so this leaves
   * the real-log threshold unchanged while reserving enough byte headroom
   * for rounding error at the upper end of the supported range. */
  scale = (128.0 - free_units) / denom;
  ctx->params.log_scale = scale;
  for (i = 0; i < ctx->params.fb_size; i++) {
    uint32_t weight = (uint32_t)(log((double)ctx->fb[i].p) * scale + 0.5);
    if (weight == 0) weight = 1;
    if (weight > 127) weight = 127;
    ctx->fb[i].logp = (uint8_t)weight;
  }

  /* Dense progressions for the smallest primes dominate sieve traffic.  Use
   * their expected log contribution as the byte initializer, start at about
   * cbrt(FB), and replace that estimate with exact division in the candidate
   * postfilter. */
  ctx->params.sieve_start = (uint32_t)cbrt((double)ctx->params.fb_size);
  if (ctx->params.sieve_start < 1)
    ctx->params.sieve_start = 1;
  /* At 300 bits cbrt(FB) still lands among very dense progressions.  A full
   * run improved by 6.8% when the exact postfilter handled primes below 384.
   * Treat this as a floor, not a replacement for cbrt(FB), so naturally
   * larger cutoffs at the upper end of the supported range never move back. */
  if (ctx->params.sieve_start_prime_floor != 0) {
    while (ctx->params.sieve_start < ctx->params.fb_size &&
           ctx->fb[ctx->params.sieve_start].p <
               ctx->params.sieve_start_prime_floor)
      ctx->params.sieve_start++;
    if (ctx->params.sieve_start >= ctx->params.fb_size)
      croak("SIQS: first-sieved-prime floor exceeds the factor base");
  }
  for (i = 0; i < ctx->params.sieve_start; i++) {
    double contribution = (double)ctx->fb[i].logp / ctx->fb[i].p;
    /* The d=2 normalized polynomial vanishes in both residue classes mod 2;
     * d=1 and odd primes dividing kN have one root.  Other factor-base primes
     * have two distinct roots before a family possibly puts one in A. */
    root_weighted_unsieved += contribution *
        (i == 0 ? (double)ctx->params.poly_d
                : ctx->fb[i].sqrt_kn == 0 ? 1.0 : 2.0);
  }
  if (ctx->params.poly_d == 2)
    root_weighted_unsieved += ctx->fb[0].logp;
  i = (uint32_t)(free_units + root_weighted_unsieved + 0.5);
  if (i > 127)
    i = 127;
  ctx->params.sieve_initial = (uint8_t)i;
  ctx->active_sieve_initial = ctx->params.sieve_initial;
}

static INLINE void siqs_sieve_add(uint8_t *cell, uint8_t logp) {
#ifdef SIQS_DEBUG
  uint32_t value = (uint32_t)*cell + logp;
  if (value > UINT8_MAX)
    croak("SIQS: byte sieve score overflow");
  *cell = (uint8_t)value;
#else
  *cell += logp;
#endif
}

static INLINE void siqs_sieve_one_root(uint8_t *sieve, uint32_t length,
                                       uint32_t root, uint32_t p,
                                       uint8_t logp) {
  uint32_t pos = root;
  for (; pos + 4 * p < length; pos += 4 * p) {
    siqs_sieve_add(sieve + pos, logp);
    siqs_sieve_add(sieve + pos + p, logp);
    siqs_sieve_add(sieve + pos + 2 * p, logp);
    siqs_sieve_add(sieve + pos + 3 * p, logp);
  }
  for (; pos < length; pos += p)
    siqs_sieve_add(sieve + pos, logp);
}

/* Alternating-root-gap sieve structure inspired by PARI/GP's MPQS sieve. */
static INLINE void siqs_sieve_two_roots(uint8_t *sieve, uint32_t length,
                                        uint32_t root1, uint32_t root2,
                                        uint32_t p, uint8_t logp) {
  uint32_t pos, gap1, gap2;
  uint32_t p4 = 4 * p;

  if (root1 > root2) {
    uint32_t tmp = root1;
    root1 = root2;
    root2 = tmp;
  }

  pos = root1;
  gap1 = root2 - root1;
  gap2 = p - gap1;
  while (pos + p4 < length) {
    siqs_sieve_add(sieve + pos, logp); pos += gap1;
    siqs_sieve_add(sieve + pos, logp); pos += gap2;
    siqs_sieve_add(sieve + pos, logp); pos += gap1;
    siqs_sieve_add(sieve + pos, logp); pos += gap2;
    siqs_sieve_add(sieve + pos, logp); pos += gap1;
    siqs_sieve_add(sieve + pos, logp); pos += gap2;
    siqs_sieve_add(sieve + pos, logp); pos += gap1;
    siqs_sieve_add(sieve + pos, logp); pos += gap2;
  }
  while (pos < length) {
    siqs_sieve_add(sieve + pos, logp);
    pos += gap1;
    if (pos >= length)
      break;
    siqs_sieve_add(sieve + pos, logp);
    pos += gap2;
  }
}

static uint8_t siqs_physical_sieve_initial(const siqs_ctx_t *ctx) {
  uint32_t initial = (uint32_t)ctx->active_sieve_initial
                   + ctx->params.stage1_bias;
  if (initial > UINT8_MAX)
    croak("SIQS: stage-1 sieve initializer is out of range");
  return (uint8_t)initial;
}

static void siqs_run_sieve(siqs_ctx_t *ctx) {
  uint32_t i;
  ctx->sieve_offset = 0;
  ctx->active_sieve_length = ctx->sieve_length;
  memset(ctx->sieve, siqs_physical_sieve_initial(ctx), ctx->sieve_length);
  for (i = ctx->params.sieve_start; i < ctx->params.fb_size; i++) {
    siqs_fb_t *fb = &ctx->fb[i];
    if (ctx->root2[i] == SIQS_NO_ROOT) {
      siqs_sieve_one_root(ctx->sieve, ctx->sieve_length,
                          ctx->root1[i], ctx->prime[i], fb->logp);
    } else {
      siqs_sieve_two_roots(ctx->sieve, ctx->sieve_length,
                           ctx->root1[i], ctx->root2[i], ctx->prime[i],
                           fb->logp);
    }
  }
}

static void siqs_clear_candidate_map(siqs_ctx_t *ctx) {
  uint32_t i;
  for (i = 0; i < ctx->candidate_count; i++) {
    uint32_t absolute = (uint32_t)(ctx->candidates[i].x
                                 + (int32_t)ctx->params.half_interval);
    uint32_t pos = absolute - ctx->sieve_offset;
#ifdef SIQS_DEBUG
    if (absolute < ctx->sieve_offset ||
        pos >= ctx->active_sieve_length)
      croak("SIQS: candidate lies outside its sieve block");
#endif
    ctx->candidate_at[pos] = 0;
    if (ctx->candidate_wide)
      ctx->candidate_at_wide[pos] = 0;
  }
  ctx->candidate_wide = 0;
  ctx->candidate_count = 0;
  ctx->hit_count = 0;
}

static void siqs_add_candidate(siqs_ctx_t *ctx, uint32_t pos) {
  uint32_t physical_score = ctx->sieve[pos];
  uint32_t logical_score;
  if (physical_score < 128U ||
      physical_score < ctx->params.stage1_bias)
    croak("SIQS: invalid stage-1 candidate score");
  logical_score = physical_score - ctx->params.stage1_bias;
  if (ctx->candidate_count == UINT16_MAX) {
    uint32_t i;
    if (ctx->candidate_at_wide == NULL)
      ctx->candidate_at_wide = (uint32_t *)siqs_calloc(
          ctx->sieve_length, sizeof(uint32_t));
    for (i = 0; i < ctx->candidate_count; i++) {
      uint32_t old_absolute = (uint32_t)(ctx->candidates[i].x
                              + (int32_t)ctx->params.half_interval);
      uint32_t old_pos = old_absolute - ctx->sieve_offset;
#ifdef SIQS_DEBUG
      if (old_absolute < ctx->sieve_offset ||
          old_pos >= ctx->active_sieve_length)
        croak("SIQS: candidate migration crossed a sieve block");
#endif
      ctx->candidate_at_wide[old_pos] = i + 1;
    }
    ctx->candidate_wide = 1;
  }
  if (ctx->candidate_count == ctx->candidate_alloc) {
    ctx->candidate_alloc = ctx->candidate_alloc
                         ? ctx->candidate_alloc * 2 : 1024;
    ctx->candidates = (siqs_candidate_t *)siqs_realloc(
        ctx->candidates,
        (size_t)ctx->candidate_alloc * sizeof(*ctx->candidates));
  }
  ctx->candidates[ctx->candidate_count].x =
      (int32_t)(ctx->sieve_offset + pos)
      - (int32_t)ctx->params.half_interval;
  ctx->candidates[ctx->candidate_count].first_hit = SIQS_NO_INDEX;
  ctx->candidates[ctx->candidate_count].sieve_score =
      (uint8_t)logical_score;
  if (ctx->candidate_wide)
    ctx->candidate_at_wide[pos] = ctx->candidate_count + 1;
  else
    ctx->candidate_at[pos] = (uint16_t)(ctx->candidate_count + 1);
  ctx->candidate_count++;
}

static void siqs_find_candidates(siqs_ctx_t *ctx) {
  static const uint64_t high_bits = UINT64_C(0x8080808080808080);
  uint32_t pos, b;
  uint32_t bulk_length = ctx->active_sieve_length & ~31U;
  siqs_clear_candidate_map(ctx);
  for (pos = 0; pos < bulk_length; pos += 32) {
    uint64_t w0, w1, w2, w3;
    memcpy(&w0, ctx->sieve + pos,      sizeof(w0));
    memcpy(&w1, ctx->sieve + pos + 8,  sizeof(w1));
    memcpy(&w2, ctx->sieve + pos + 16, sizeof(w2));
    memcpy(&w3, ctx->sieve + pos + 24, sizeof(w3));
    if (!((w0 | w1 | w2 | w3) & high_bits))
      continue;
    for (b = 0; b < 32; b++)
      if (ctx->sieve[pos + b] & 0x80U)
        siqs_add_candidate(ctx, pos + b);
  }
  for (; pos < ctx->active_sieve_length; pos++)
    if (ctx->sieve[pos] & 0x80U)
      siqs_add_candidate(ctx, pos);
  ctx->total_candidates += ctx->candidate_count;
}

static void siqs_add_hit(siqs_ctx_t *ctx, uint32_t candidate,
                         uint32_t fb_index) {
  siqs_candidate_t *c = &ctx->candidates[candidate];
  if (ctx->hit_count == ctx->hit_alloc) {
    ctx->hit_alloc = ctx->hit_alloc ? ctx->hit_alloc * 2 : 8192;
    ctx->hits = (siqs_hit_t *)siqs_realloc(
        ctx->hits, (size_t)ctx->hit_alloc * sizeof(*ctx->hits));
  }
  ctx->hits[ctx->hit_count].fb_index = fb_index;
  ctx->hits[ctx->hit_count].next = c->first_hit;
  c->first_hit = ctx->hit_count++;
}

static INLINE void siqs_resieve_one_root16(siqs_ctx_t *ctx,
                                           uint32_t fb_index,
                                           uint32_t root, uint32_t p) {
  uint32_t pos;
  for (pos = root; pos < ctx->active_sieve_length; pos += p) {
    uint16_t candidate = ctx->candidate_at[pos];
    if (candidate != 0) {
      siqs_add_hit(ctx, candidate - 1, fb_index);
    }
  }
}

static INLINE void siqs_resieve_one_root32(siqs_ctx_t *ctx,
                                           uint32_t fb_index,
                                           uint32_t root, uint32_t p) {
  uint32_t pos;
  for (pos = root; pos < ctx->active_sieve_length; pos += p) {
    uint32_t candidate = ctx->candidate_at_wide[pos];
    if (candidate != 0) {
      siqs_add_hit(ctx, candidate - 1, fb_index);
    }
  }
}

static INLINE uint32_t siqs_reduce_u32(uint32_t n, uint32_t p,
                                       uint32_t reciprocal) {
  uint32_t q = (uint32_t)(((uint64_t)n * reciprocal) >> 32);
  uint32_t r = n - q * p;
  return r >= p ? r - p : r;
}

static void siqs_filter_q(mpz_t q, const siqs_poly_t *poly,
                          const mpz_t twice_b, int32_t x) {
  /* q(x) = d*A*x^2 + 2*B*x + C.  Horner form avoids the square and exact
   * division used later when constructing the relation's square root. */
  mpz_mul_si(q, poly->DA, x);
  mpz_add(q, q, twice_b);
  mpz_mul_si(q, q, x);
  mpz_add(q, q, poly->C);
}

/* Factor-base primes fit unsigned long, so avoid constructing an mpz divisor
 * and entering the general mpz_remove machinery for every resieve hit.  The
 * 2-adic valuation has a direct bit operation; odd primes normally occur to
 * the first power.  Both callers exclude zero, and returning zero for a
 * nonfactor preserves their invariant checks. */
static INLINE mp_bitcnt_t siqs_remove_ui(mpz_t value, unsigned long p) {
  mp_bitcnt_t exponent = 0;
  if (p == 2) {
    exponent = mpz_scan1(value, 0);
    mpz_tdiv_q_2exp(value, value, exponent);
    return exponent;
  }
  while (mpz_divisible_ui_p(value, p)) {
    mpz_divexact_ui(value, value, p);
    exponent++;
  }
  return exponent;
}

/* Replace the family initializer's expected omitted-prime contribution with
 * exact trial division, including powers, and compare the predicted residual
 * directly with the accepted smooth-product bound. */
static void siqs_filter_candidates(siqs_ctx_t *ctx,
                                   const siqs_poly_t *poly) {
  uint32_t original_count = ctx->candidate_count;
  uint32_t read, out = 0;
  uint32_t small_hits[SIQS_POSTFILTER_MAX_SMALL];
  mpz_t q, twice_b;
  double log_smooth = log((double)ctx->params.smooth_bound);

  if (ctx->params.sieve_start > SIQS_POSTFILTER_MAX_SMALL)
    croak("SIQS: too many omitted primes for the candidate postfilter");
#ifdef SIQS_DEBUG
  if (ctx->hit_count != 0)
    croak("SIQS: candidate postfilter started with resieve hits");
#endif
  mpz_init(q);
  mpz_init(twice_b);
  mpz_mul_2exp(twice_b, poly->B, 1);

  for (read = 0; read < original_count; read++) {
    siqs_candidate_t candidate = ctx->candidates[read];
    uint32_t absolute = (uint32_t)(candidate.x
                                 + (int32_t)ctx->params.half_interval);
    uint32_t local = absolute - ctx->sieve_offset;
    uint32_t hit_count = 0, i;
    double fine_margin;
    int would_pass;

#ifdef SIQS_DEBUG
    if (absolute < ctx->sieve_offset ||
        local >= ctx->active_sieve_length)
      croak("SIQS: postfilter candidate lies outside its sieve block");
#endif
    ctx->candidate_at[local] = 0;
    if (ctx->candidate_wide)
      ctx->candidate_at_wide[local] = 0;

    siqs_filter_q(q, poly, twice_b, candidate.x);
    if (mpz_sgn(q) == 0) {
      would_pass = 1;
    } else {
      mpz_abs(q, q);
      for (i = 0; i < ctx->params.sieve_start; i++) {
        uint32_t p = ctx->prime[i];
        uint32_t rem = siqs_reduce_u32(absolute, p,
                                      ctx->fb_reciprocal[i]);
        if (rem == ctx->root1[i] || rem == ctx->root2[i]) {
          mp_bitcnt_t exponent;
          exponent = siqs_remove_ui(q, p);
          if (exponent == 0)
            croak("SIQS: postfilter root is not a polynomial factor");
          small_hits[hit_count++] = i;
        }
      }
      fine_margin = (double)candidate.sieve_score
                  - (double)ctx->active_sieve_initial
                  - ctx->params.log_scale *
                    (siqs_log_mpz(q) - log_smooth);
      would_pass = fine_margin >= 0.0;
    }

    if (would_pass) {
      ctx->candidates[out] = candidate;
      ctx->candidates[out].first_hit = SIQS_NO_INDEX;
      if (ctx->candidate_wide)
        ctx->candidate_at_wide[local] = out + 1;
      else
        ctx->candidate_at[local] = (uint16_t)(out + 1);
      for (i = 0; i < hit_count; i++)
        siqs_add_hit(ctx, out, small_hits[i]);
      out++;
    }
  }
  ctx->candidate_count = out;
  mpz_clear(q);
  mpz_clear(twice_b);
}

/* Translate a whole-interval root to the current block without changing the
 * roots maintained by the polynomial code. */
static INLINE uint32_t siqs_local_root(const siqs_ctx_t *ctx,
                                       uint32_t fb_index, uint32_t root) {
  uint32_t p, offset_mod;
  if (ctx->sieve_offset == 0)
    return root;
  p = ctx->prime[fb_index];
  offset_mod = siqs_reduce_u32(ctx->sieve_offset, p,
                              ctx->fb_reciprocal[fb_index]);
  return root >= offset_mod ? root - offset_mod
                            : root + p - offset_mod;
}

static void siqs_bucket_count_root(siqs_ctx_t *ctx,
                                   uint32_t root, uint32_t p) {
  uint32_t pos;
  for (pos = root; pos < ctx->sieve_length; ) {
    uint32_t block = pos >> SIQS_SIEVE_BLOCK_BITS;
    if (ctx->bucket_bounds[block + 1] == UINT32_MAX)
      croak("SIQS: too many sieve bucket entries");
    ctx->bucket_bounds[block + 1]++;
    if (p > UINT32_MAX - pos)
      break;
    pos += p;
  }
}

static void siqs_bucket_fill_root(siqs_ctx_t *ctx, uint32_t fb_index,
                                  uint32_t root, uint32_t p) {
  uint32_t pos;
  for (pos = root; pos < ctx->sieve_length; ) {
    uint32_t block = pos >> SIQS_SIEVE_BLOCK_BITS;
    uint32_t slot = ctx->bucket_fill[block]++;
#ifdef SIQS_DEBUG
    if (slot >= ctx->bucket_event_count)
      croak("SIQS: sieve bucket write overflow");
#endif
    ctx->bucket_events[slot] =
        (fb_index << SIQS_SIEVE_BLOCK_BITS)
      | (pos & SIQS_SIEVE_BLOCK_MASK);
    if (p > UINT32_MAX - pos)
      break;
    pos += p;
  }
}

/* Large primes are visited once per actual sieve hit, rather than once per
 * block.  Keeping the factor-base index in each event lets resieving reuse
 * the same buckets. */
static void siqs_build_buckets(siqs_ctx_t *ctx) {
  uint32_t i, block, total, alloc;
  uint32_t bucket_start = ctx->block_large_index;
  /* The postfilter has already attached every hit below sieve_start. */
  if (bucket_start < ctx->params.sieve_start)
    bucket_start = ctx->params.sieve_start;
  memset(ctx->bucket_bounds, 0,
         (size_t)(ctx->block_count + 1) * sizeof(uint32_t));
  for (i = bucket_start; i < ctx->params.fb_size; i++) {
    uint32_t p = ctx->prime[i];
    siqs_bucket_count_root(ctx, ctx->root1[i], p);
    if (ctx->root2[i] != SIQS_NO_ROOT)
      siqs_bucket_count_root(ctx, ctx->root2[i], p);
  }
  for (block = 0; block < ctx->block_count; block++) {
    if (ctx->bucket_bounds[block + 1] >
        UINT32_MAX - ctx->bucket_bounds[block])
      croak("SIQS: sieve bucket index overflow");
    ctx->bucket_bounds[block + 1] += ctx->bucket_bounds[block];
  }
  total = ctx->bucket_bounds[ctx->block_count];
  if (total > ctx->bucket_event_alloc) {
    alloc = ctx->bucket_event_alloc ? ctx->bucket_event_alloc : 65536U;
    while (alloc < total) {
      if (alloc > UINT32_MAX / 2U) {
        alloc = total;
        break;
      }
      alloc *= 2U;
    }
    ctx->bucket_events = (uint32_t *)siqs_realloc(
        ctx->bucket_events, (size_t)alloc * sizeof(uint32_t));
    ctx->bucket_event_alloc = alloc;
  }
  ctx->bucket_event_count = total;
  memcpy(ctx->bucket_fill, ctx->bucket_bounds,
         (size_t)ctx->block_count * sizeof(uint32_t));
  for (i = bucket_start; i < ctx->params.fb_size; i++) {
    uint32_t p = ctx->prime[i];
    siqs_bucket_fill_root(ctx, i, ctx->root1[i], p);
    if (ctx->root2[i] != SIQS_NO_ROOT)
      siqs_bucket_fill_root(ctx, i, ctx->root2[i], p);
  }
#ifdef SIQS_DEBUG
  for (block = 0; block < ctx->block_count; block++)
    if (ctx->bucket_fill[block] != ctx->bucket_bounds[block + 1])
      croak("SIQS: incomplete sieve bucket");
#endif
}

static void siqs_run_sieve_block(siqs_ctx_t *ctx, uint32_t block) {
  uint32_t i, event, begin, end;
  uint32_t offset = block << SIQS_SIEVE_BLOCK_BITS;
  uint32_t remaining = ctx->sieve_length - offset;
  uint32_t length = remaining < SIQS_SIEVE_BLOCK_SIZE
                  ? remaining : SIQS_SIEVE_BLOCK_SIZE;
  ctx->sieve_offset = offset;
  ctx->active_sieve_length = length;
  memset(ctx->sieve, siqs_physical_sieve_initial(ctx), length);

  for (i = ctx->params.sieve_start; i < ctx->block_large_index; i++) {
    uint32_t p = ctx->prime[i];
    uint32_t root1 = siqs_local_root(ctx, i, ctx->root1[i]);
    if (ctx->root2[i] == SIQS_NO_ROOT) {
      siqs_sieve_one_root(ctx->sieve, length, root1, p, ctx->fb[i].logp);
    } else {
      uint32_t root2 = siqs_local_root(ctx, i, ctx->root2[i]);
      siqs_sieve_two_roots(ctx->sieve, length, root1, root2, p,
                           ctx->fb[i].logp);
    }
  }

  begin = ctx->bucket_bounds[block];
  end = ctx->bucket_bounds[block + 1];
  for (event = begin; event < end; event++) {
    uint32_t packed = ctx->bucket_events[event];
    uint32_t fb_index = packed >> SIQS_SIEVE_BLOCK_BITS;
    uint32_t pos = packed & SIQS_SIEVE_BLOCK_MASK;
#ifdef SIQS_DEBUG
    uint32_t absolute, rem;
    if (fb_index < ctx->block_large_index ||
        fb_index >= ctx->params.fb_size || pos >= length)
      croak("SIQS: invalid sieve bucket event");
    absolute = offset + pos;
    rem = absolute % ctx->prime[fb_index];
    if (rem != ctx->root1[fb_index] && rem != ctx->root2[fb_index])
      croak("SIQS: sieve bucket event is not a polynomial root");
#endif
    siqs_sieve_add(ctx->sieve + pos, ctx->fb[fb_index].logp);
  }
}

static void siqs_resieve_candidates(siqs_ctx_t *ctx,
                                    uint32_t factor_begin,
                                    uint32_t progression_end,
                                    uint32_t bucket_begin,
                                    uint32_t bucket_end) {
  uint32_t i, c, cutoff = factor_begin, cutoff_prime;
  if (ctx->candidate_count == 0) {
    return;
  }
  {
    uint64_t proposed =
        (uint64_t)SIQS_RESIEVE_CUTOFF_COEFFICIENT *
        ctx->active_sieve_length / ctx->candidate_count;
    cutoff_prime = proposed > UINT32_MAX
                 ? UINT32_MAX : (uint32_t)proposed;
  }
  if (cutoff_prime < SIQS_RESIEVE_CUTOFF_FLOOR)
    cutoff_prime = SIQS_RESIEVE_CUTOFF_FLOOR;
  while (cutoff < progression_end &&
         ctx->prime[cutoff] <= cutoff_prime)
    cutoff++;

  /* Rewalking the dense progressions of the smallest primes costs more than
   * testing the handful of candidates directly. */
  for (c = 0; c < ctx->candidate_count; c++) {
    uint32_t pos = (uint32_t)(ctx->candidates[c].x
                            + (int32_t)ctx->params.half_interval);
    for (i = factor_begin; i < cutoff; i++) {
      uint32_t rem = siqs_reduce_u32(pos, ctx->prime[i],
                                    ctx->fb_reciprocal[i]);
      if (rem == ctx->root1[i] || rem == ctx->root2[i]) {
        siqs_add_hit(ctx, c, i);
      }
    }
  }
  if (ctx->candidate_wide) {
    for (i = cutoff; i < progression_end; i++) {
      uint32_t root1 = siqs_local_root(ctx, i, ctx->root1[i]);
      siqs_resieve_one_root32(ctx, i, root1, ctx->prime[i]);
      if (ctx->root2[i] != SIQS_NO_ROOT) {
        siqs_resieve_one_root32(ctx, i,
            siqs_local_root(ctx, i, ctx->root2[i]), ctx->prime[i]);
      }
    }
  } else {
    for (i = cutoff; i < progression_end; i++) {
      uint32_t root1 = siqs_local_root(ctx, i, ctx->root1[i]);
      siqs_resieve_one_root16(ctx, i, root1, ctx->prime[i]);
      if (ctx->root2[i] != SIQS_NO_ROOT) {
        siqs_resieve_one_root16(ctx, i,
            siqs_local_root(ctx, i, ctx->root2[i]), ctx->prime[i]);
      }
    }
  }

  /* Every large-prime event is an actual progression hit, and carries the
   * factor-base index needed by the relation evaluator. */
  if (ctx->candidate_wide) {
    for (i = bucket_begin; i < bucket_end; i++) {
      uint32_t packed = ctx->bucket_events[i];
      uint32_t pos = packed & SIQS_SIEVE_BLOCK_MASK;
      uint32_t candidate = ctx->candidate_at_wide[pos];
      if (candidate != 0) {
        siqs_add_hit(ctx, candidate - 1,
                     packed >> SIQS_SIEVE_BLOCK_BITS);
      }
    }
  } else {
    for (i = bucket_begin; i < bucket_end; i++) {
      uint32_t packed = ctx->bucket_events[i];
      uint32_t pos = packed & SIQS_SIEVE_BLOCK_MASK;
      uint16_t candidate = ctx->candidate_at[pos];
      if (candidate != 0) {
        siqs_add_hit(ctx, candidate - 1,
                     packed >> SIQS_SIEVE_BLOCK_BITS);
      }
    }
  }
}

static int siqs_mpz_to_u64(const mpz_t n, uint64_t *value) {
  uint64_t out = 0;
  size_t count = 0;
  if (mpz_sgn(n) < 0 || mpz_sizeinbase(n, 2) > 64)
    return 0;
  mpz_export(&out, &count, -1, sizeof(out), 0, 0, n);
  *value = out;
  return 1;
}

static int siqs_u64_probable_prime(uint64_t n) {
  int result;
  mpz_t z;
  mpz_init(z);
  mpz_import(z, 1, -1, sizeof(n), 0, 0, &n);
  result = mpz_probab_prime_p(z, 3) != 0;
  mpz_clear(z);
  return result;
}

/* Resolve a cofactor into zero, one, or two graph large primes. */
static int siqs_resolve_cofactor(siqs_ctx_t *ctx, const mpz_t rest,
                                 uint64_t *lp1, uint64_t *lp2) {
  uint64_t n, pmax2;
  *lp1 = *lp2 = 1;
  if (mpz_cmp_ui(rest, 1) == 0)
    return 1;
  if (!siqs_mpz_to_u64(rest, &n))
    return 0;
  if (n > ctx->params.smooth_bound)
    return 0;

  /* Once all factor-base primes have been removed, a residual below pmax^2
   * cannot be composite.  This is the common one-large-prime path and avoids
   * a probable-prime test for every relation. */
  pmax2 = (uint64_t)ctx->largest_fb_prime * ctx->largest_fb_prime;
  if (n < pmax2) {
    if (n <= ctx->largest_fb_prime)
      return 0;
    if (n > ctx->params.large_prime_bound)
      return 0;
    *lp2 = n;
    return 1;
  }
  if (mpz_probab_prime_p(rest, 2)) {
    if (n <= ctx->largest_fb_prime)
      return 0;
    if (n > ctx->params.large_prime_bound)
      return 0;
    *lp2 = n;
    return 1;
  }
  if (ctx->params.max_large_primes < 2)
    return 0;
  ctx->split_attempts++;

#if BITS_PER_WORD == 64 && HAVE_STD_U64 && defined(__GNUC__) && defined(__x86_64__)
  {
    UV factors[2];
    int count = uvpbrent63((UV)n, factors, 20000,
                           (UV)(siqs_rand64(&ctx->cofactor_rng) | 1U));
    uint64_t a, b;
    if (count != 2) {
      ctx->split_failures++;
      return 0;
    }
    a = factors[0];
    b = factors[1];
    if (a == 0 || n % a != 0 || n / a != b) {
      ctx->split_failures++;
      return 0;
    }
    if (a <= ctx->largest_fb_prime || b <= ctx->largest_fb_prime) {
      ctx->split_failures++;
      return 0;
    }
    if (a > ctx->params.large_prime_bound ||
        b > ctx->params.large_prime_bound) {
      ctx->split_failures++;
      return 0;
    }
    if (!siqs_u64_probable_prime(a) || !siqs_u64_probable_prime(b)) {
      ctx->split_failures++;
      return 0;
    }
    if (a > b) {
      uint64_t t = a; a = b; b = t;
    }
    *lp1 = a;
    *lp2 = b;
    ctx->split_rho++;
    return 1;
  }
#else
  {
    mpz_t factor, quotient;
    uint32_t nbits = (uint32_t)mpz_sizeinbase(rest, 2);
    UV rounds = nbits <= 40 ? 20000
              : nbits <= 44 ? 50000
              : nbits <= 48 ? 100000
              : nbits <= 52 ? 200000 : 500000;
    uint64_t a, b;
    int success;
    mpz_init(factor);
    mpz_init(quotient);
    if (mpz_perfect_square_p(rest)) {
      mpz_sqrt(factor, rest);
      success = 1;
      ctx->split_squares++;
    } else {
      success = squfof126(rest, factor, rounds);
      if (success) {
        ctx->split_squfof++;
      } else {
        success = _GMP_pbrent_factor(
            rest, factor,
            (UV)(3 + (siqs_rand64(&ctx->cofactor_rng) & 0xffffU)),
            250000);
        if (success)
          ctx->split_rho++;
      }
    }
    if (!success) {
      ctx->split_failures++;
      mpz_clear(factor);
      mpz_clear(quotient);
      return 0;
    }
    if (!mpz_divisible_p(rest, factor)) {
      ctx->split_failures++;
      mpz_clear(factor);
      mpz_clear(quotient);
      return 0;
    }
    mpz_divexact(quotient, rest, factor);
    success = siqs_mpz_to_u64(factor, &a)
           && siqs_mpz_to_u64(quotient, &b);
    mpz_clear(factor);
    mpz_clear(quotient);
    if (!success) {
      ctx->split_failures++;
      return 0;
    }
    if (a <= ctx->largest_fb_prime || b <= ctx->largest_fb_prime) {
      ctx->split_failures++;
      return 0;
    }
    if (a > ctx->params.large_prime_bound ||
        b > ctx->params.large_prime_bound) {
      ctx->split_failures++;
      return 0;
    }
    if (ctx->params.smooth_bound / pmax2 >= ctx->largest_fb_prime &&
        (!siqs_u64_probable_prime(a) || !siqs_u64_probable_prime(b))) {
      ctx->split_failures++;
      return 0;
    }
    if (a > b) {
      uint64_t t = a; a = b; b = t;
    }
    *lp1 = a;
    *lp2 = b;
    return 1;
  }
#endif
}

static siqs_factor_t *siqs_eval_reserve_factors(siqs_ctx_t *ctx,
                                                uint32_t required) {
  siqs_eval_workspace_t *workspace = &ctx->eval;
  uint32_t alloc;
  if (required <= workspace->factor_alloc)
    return workspace->factors;
  alloc = workspace->factor_alloc;
  if (alloc == 0)
    alloc = SIQS_EVAL_INITIAL_FACTORS;
  else if (alloc <= UINT32_MAX / 2)
    alloc *= 2;
  if (alloc < required)
    alloc = required;
  workspace->factors = (siqs_factor_t *)siqs_realloc(
      workspace->factors, (size_t)alloc * sizeof(*workspace->factors));
  workspace->factor_alloc = alloc;
  return workspace->factors;
}

/* Resieving visits factor-base entries in increasing order and prepends each
 * hit, so a candidate's hit list is strictly descending.  The evaluator
 * stores that list backwards to make an ascending run, then merges the
 * sign/d/A run into it.  Keeping the hit run after its maximum-sized prefix
 * makes this in-place merge safe even when an extra factor sorts first. */
static uint32_t siqs_merge_candidate_factors(
    siqs_factor_t *factors, uint32_t hit_offset, uint32_t hit_count,
    const siqs_factor_t *extra, uint32_t extra_count) {
  uint32_t hit = 0, add = 0, out = 0;
  while (hit < hit_count || add < extra_count) {
    siqs_factor_t next;
    if (hit == hit_count ||
        (add < extra_count &&
         extra[add].row < factors[hit_offset + hit].row)) {
      next = extra[add++];
    } else if (add == extra_count ||
               factors[hit_offset + hit].row < extra[add].row) {
      next = factors[hit_offset + hit++];
    } else {
      next = extra[add++];
      next.exponent += factors[hit_offset + hit++].exponent;
    }
    factors[out++] = next;
  }
#ifdef SIQS_DEBUG
  for (hit = 1; hit < out; hit++)
    if (factors[hit - 1].row >= factors[hit].row)
      croak("SIQS: candidate factors are not strictly ordered");
#endif
  return out;
}

static void siqs_evaluate_candidate(siqs_ctx_t *ctx, const siqs_poly_t *poly,
                                    uint32_t candidate_index) {
  const siqs_candidate_t *candidate = &ctx->candidates[candidate_index];
  siqs_factor_t extra[SIQS_EVAL_MAX_EXTRA_FACTORS];
  siqs_factor_t *factors;
  uint32_t count, extra_count = 0, hit_count = 0, hit;
  uint32_t hit_offset = poly->q_count + 2;
  uint32_t required = hit_offset;
  uint64_t lp1, lp2;
  int32_t x = candidate->x;
  mpz_ptr y = ctx->eval.y;
  mpz_ptr q = ctx->eval.q;
  mpz_ptr rest = ctx->eval.rest;


  for (hit = candidate->first_hit; hit != SIQS_NO_INDEX;
       hit = ctx->hits[hit].next) {
    required++;
    hit_count++;
  }
  factors = siqs_eval_reserve_factors(ctx, required);

  mpz_mul_si(y, poly->DA, x);
  mpz_add(y, y, poly->B);
  /* q(x) = d*A*x^2 + 2*B*x + C.  Since y = d*A*x + B is
   * required by the relation, form the Horner coefficient as y + B and
   * avoid squaring y followed by an exact division by d*A. */
  mpz_add(q, y, poly->B);
  mpz_mul_si(q, q, x);
  mpz_add(q, q, poly->C);
#ifdef SIQS_DEBUG
  /* Keep the independent defining identity as a shadow check while the
   * evaluator rewrite is exercised by debug builds. */
  mpz_mul(rest, y, y);
  mpz_sub(rest, rest, ctx->kn);
  if (!mpz_divisible_p(rest, poly->DA))
    croak("SIQS: candidate polynomial division is not exact");
  mpz_divexact(rest, rest, poly->DA);
  if (mpz_cmp(q, rest) != 0)
    croak("SIQS: Horner candidate polynomial mismatch");
#endif
  if (mpz_sgn(q) == 0) {
    mpz_gcd(rest, y, ctx->n);
    if (mpz_cmp_ui(rest, 1) > 0 && mpz_cmp(rest, ctx->n) < 0 &&
        siqs_insert_divisor(ctx->result, rest))
      ctx->factor_found = 1;
    return;
  }
  if (mpz_sgn(q) < 0) {
    extra[extra_count].row = 0;
    extra[extra_count].exponent = 1;
    extra_count++;
    mpz_neg(rest, q);
  } else {
    mpz_set(rest, q);
  }

  {
    uint32_t write = hit_offset + hit_count;
#ifdef SIQS_DEBUG
    uint32_t previous_fb_index = SIQS_NO_INDEX;
#endif
    for (hit = candidate->first_hit; hit != SIQS_NO_INDEX;
         hit = ctx->hits[hit].next) {
      uint32_t fb_index = ctx->hits[hit].fb_index;
      uint32_t p = ctx->fb[fb_index].p;
      uint32_t exponent;
#ifdef SIQS_DEBUG
      if (previous_fb_index != SIQS_NO_INDEX &&
          fb_index >= previous_fb_index)
        croak("SIQS: candidate resieve hits are not strictly descending");
      previous_fb_index = fb_index;
#endif
      exponent = (uint32_t)siqs_remove_ui(rest, p);
      if (exponent == 0)
        croak("SIQS: resieve reported a nonfactor");
      write--;
      factors[write].row = fb_index + 1;
      factors[write].exponent = exponent;
    }
    if (write != hit_offset)
      croak("SIQS: candidate resieve hit count changed during evaluation");
  }

  {
    uint32_t i;
    if (ctx->params.poly_d == 2) {
      /* The numerator is d*A*q(x), so d contributes to every relation just
       * like the selected odd primes of A.  q(x) is always divisible by four;
       * its complete 2-adic exponent was removed through the root hit. */
      extra[extra_count].row = 1;
      extra[extra_count].exponent = 1;
      extra_count++;
    }
    for (i = 0; i < poly->q_count; i++) {
      extra[extra_count].row = poly->a_index[i] + 1;
      extra[extra_count].exponent = 1;
      extra_count++;
    }
  }
  count = siqs_merge_candidate_factors(
      factors, hit_offset, hit_count, extra, extra_count);
  if (siqs_resolve_cofactor(ctx, rest, &lp1, &lp2)) {
    siqs_raw_relation_t *relation;
    mpz_mod(y, y, ctx->n);
    relation = siqs_raw_relation_new(y, factors, count, lp1, lp2);
    siqs_accept_raw_relation(ctx, relation);
  }
}

static void siqs_sieve_polynomial(siqs_ctx_t *ctx, siqs_poly_t *poly) {
  uint32_t i, block;
  if (!ctx->use_block_sieve) {
    siqs_run_sieve(ctx);
    siqs_find_candidates(ctx);
    siqs_filter_candidates(ctx, poly);
    siqs_resieve_candidates(ctx,
        ctx->params.sieve_start,
        ctx->params.fb_size, 0, 0);
    for (i = 0; i < ctx->candidate_count && !ctx->factor_found; i++)
      siqs_evaluate_candidate(ctx, poly, i);
    siqs_clear_candidate_map(ctx);
    return;
  }

  siqs_build_buckets(ctx);
  for (block = 0; block < ctx->block_count && !ctx->factor_found; block++) {
    siqs_run_sieve_block(ctx, block);
    siqs_find_candidates(ctx);
    siqs_filter_candidates(ctx, poly);
    siqs_resieve_candidates(ctx,
                            ctx->params.sieve_start,
                            ctx->block_large_index,
                            ctx->bucket_bounds[block],
                            ctx->bucket_bounds[block + 1]);
    for (i = 0; i < ctx->candidate_count && !ctx->factor_found; i++)
      siqs_evaluate_candidate(ctx, poly, i);
    siqs_clear_candidate_map(ctx);
  }
}

/*----------------------------------------------------------------------------
 * Sparse matrix and square-root phase
 *----------------------------------------------------------------------------*/

static int siqs_column_cmp(const void *va, const void *vb) {
  const la_col_t *a = (const la_col_t *)va;
  const la_col_t *b = (const la_col_t *)vb;
  return a->weight < b->weight ? -1 : a->weight > b->weight ? 1 : 0;
}

#ifdef SIQS_DEBUG
static void siqs_verify_full_relation(const siqs_ctx_t *ctx,
                                      const siqs_full_relation_t *r) {
  uint32_t i;
  mpz_t lhs, rhs, t;
  mpz_init(lhs);
  mpz_init_set_ui(rhs, 1);
  mpz_init(t);
  mpz_mul(lhs, r->y, r->y);
  mpz_mod(lhs, lhs, ctx->n);
  for (i = 0; i < r->nfactors; i++) {
    uint32_t row = r->factors[i].row;
    if (row == 0) {
      if (r->factors[i].exponent & 1U)
        mpz_neg(rhs, rhs);
    } else {
      mpz_set_ui(t, ctx->fb[row - 1].p);
      mpz_powm_ui(t, t, r->factors[i].exponent, ctx->n);
      mpz_mul(rhs, rhs, t);
      mpz_mod(rhs, rhs, ctx->n);
    }
  }
  mpz_mod(rhs, rhs, ctx->n);
  if (mpz_cmp(lhs, rhs) != 0)
    croak("SIQS: full relation invariant failed");
  mpz_clear(lhs);
  mpz_clear(rhs);
  mpz_clear(t);
}
#endif

static la_col_t *siqs_build_matrix(siqs_ctx_t *ctx,
                                   unsigned long *nrows,
                                   unsigned long *ncols) {
  uint32_t i, j;
  la_col_t *columns = (la_col_t *)siqs_calloc(
      ctx->full_count, sizeof(*columns));
  *nrows = ctx->params.fb_size + 1;
  *ncols = ctx->full_count;
  for (i = 0; i < ctx->full_count; i++) {
    const siqs_full_relation_t *r = ctx->full[i];
    unsigned long weight = 0;
#ifdef SIQS_DEBUG
    siqs_verify_full_relation(ctx, r);
#endif
    for (j = 0; j < r->nfactors; j++)
      if (r->factors[j].exponent & 1U)
        weight++;
    columns[i].data = (unsigned long *)siqs_malloc(
        (size_t)weight * sizeof(unsigned long));
    columns[i].weight = 0;
    columns[i].orig = i;
    for (j = 0; j < r->nfactors; j++)
      if (r->factors[j].exponent & 1U)
        columns[i].data[columns[i].weight++] = r->factors[j].row;
  }
  qsort(columns, ctx->full_count, sizeof(*columns), siqs_column_cmp);
  return columns;
}

/* Determine whether singleton removal leaves enough columns for Lanczos.
 * This is the same 2-core that reduce_matrix() starts from, but is computed
 * without changing or copying the full relations.  Keeping an incidence list
 * from rows to columns makes the peeling pass linear in the number of odd
 * factor-base exponents. */
static int siqs_matrix_ready(const siqs_ctx_t *ctx,
                             uint32_t *surviving_rows,
                             uint32_t *surviving_cols) {
  const uint32_t nrows = ctx->params.fb_size + 1;
  const uint32_t ncols = ctx->full_count;
  uint32_t *counts, *incidence, *queue;
  size_t *offsets, *next;
  uint8_t *removed;
  size_t entries = 0;
  uint32_t row, col, head, tail, removed_count = 0;

  counts = (uint32_t *)siqs_calloc(nrows, sizeof(*counts));
  for (col = 0; col < ncols; col++) {
    const siqs_full_relation_t *r = ctx->full[col];
    uint32_t i;
    for (i = 0; i < r->nfactors; i++) {
      if (!(r->factors[i].exponent & 1U))
        continue;
      row = r->factors[i].row;
      if (row >= nrows)
        croak("SIQS: relation contains an invalid matrix row");
      counts[row]++;
      entries++;
    }
  }

  offsets = (size_t *)siqs_malloc((size_t)(nrows + 1) * sizeof(*offsets));
  next = (size_t *)siqs_malloc((size_t)nrows * sizeof(*next));
  offsets[0] = 0;
  for (row = 0; row < nrows; row++)
    offsets[row + 1] = offsets[row] + counts[row];
  memcpy(next, offsets, (size_t)nrows * sizeof(*next));
  incidence = (uint32_t *)siqs_malloc(entries * sizeof(*incidence));
  for (col = 0; col < ncols; col++) {
    const siqs_full_relation_t *r = ctx->full[col];
    uint32_t i;
    for (i = 0; i < r->nfactors; i++) {
      if (r->factors[i].exponent & 1U) {
        row = r->factors[i].row;
        incidence[next[row]++] = col;
      }
    }
  }
  free(next);

  removed = (uint8_t *)siqs_calloc(ncols, sizeof(*removed));
  queue = (uint32_t *)siqs_malloc((size_t)nrows * sizeof(*queue));
  head = tail = 0;
  for (row = 0; row < nrows; row++)
    if (counts[row] == 1)
      queue[tail++] = row;

  while (head < tail) {
    size_t j;
    row = queue[head++];
    if (counts[row] != 1)
      continue;
    for (j = offsets[row]; j < offsets[row + 1]; j++)
      if (!removed[incidence[j]])
        break;
    if (j == offsets[row + 1])
      continue;

    col = incidence[j];
    removed[col] = 1;
    removed_count++;
    {
      const siqs_full_relation_t *r = ctx->full[col];
      uint32_t i;
      for (i = 0; i < r->nfactors; i++) {
        if (!(r->factors[i].exponent & 1U))
          continue;
        row = r->factors[i].row;
        counts[row]--;
        if (counts[row] == 1)
          queue[tail++] = row;
      }
    }
  }

  *surviving_rows = 0;
  for (row = 0; row < nrows; row++)
    if (counts[row] != 0)
      (*surviving_rows)++;
  *surviving_cols = ncols - removed_count;

  free(queue);
  free(removed);
  free(incidence);
  free(offsets);
  free(counts);

  return (uint64_t)*surviving_cols >=
         (uint64_t)*surviving_rows + SIQS_MATRIX_EXTRA_RELS(ctx);
}

static int siqs_test_dependencies(siqs_ctx_t *ctx, const la_col_t *columns,
                                  unsigned long ncols,
                                  const uint64_t *nullrows, uint64_t mask) {
  uint32_t dependency;
  mpz_t lhs, rhs, power, delta, divisor;
  mpz_init(lhs);
  mpz_init(rhs);
  mpz_init(power);
  mpz_init(delta);
  mpz_init(divisor);

  for (dependency = 0; dependency < 64 && !ctx->factor_found; dependency++) {
    unsigned long i;
    uint32_t j;
    if (!(mask & (UINT64_C(1) << dependency)))
      continue;
    siqs_reset_touched_factors(ctx);
    mpz_set_ui(lhs, 1);
    mpz_set_ui(rhs, 1);
    for (i = 0; i < ncols; i++) {
      if (!getNullEntry(nullrows, i, dependency))
        continue;
      {
        const siqs_full_relation_t *r = ctx->full[columns[i].orig];
        mpz_mul(lhs, lhs, r->y);
        mpz_mod(lhs, lhs, ctx->n);
        for (j = 0; j < r->nfactors; j++)
          siqs_touch_factor(ctx, r->factors[j].row,
                           r->factors[j].exponent);
      }
    }
    for (j = 0; j < ctx->factor_touched_count; j++) {
      uint32_t row = ctx->factor_touched[j];
      uint32_t exponent = ctx->factor_counts[row];
      if (exponent & 1U) {
#ifdef SIQS_DEBUG
        uint32_t matrix_parity = 0;
        unsigned long c, k;
        for (c = 0; c < ncols; c++) {
          if (!getNullEntry(nullrows, c, dependency))
            continue;
          for (k = 0; k < columns[c].weight; k++)
            if (columns[c].data[k] == row)
              matrix_parity ^= 1U;
        }
        fprintf(stderr,
                "SIQS: dependency %u has odd exponent %u in row %u "
                "(matrix parity %u, %lu columns)\n",
                dependency, exponent, row, matrix_parity, ncols);
#endif
        croak("SIQS: linear solver returned a nonsquare dependency");
      }
      if (row == 0) {
        if ((exponent >> 1) & 1U)
          mpz_neg(rhs, rhs);
      } else if (exponent != 0) {
        mpz_set_ui(power, ctx->fb[row - 1].p);
        mpz_powm_ui(power, power, exponent >> 1, ctx->n);
        mpz_mul(rhs, rhs, power);
        mpz_mod(rhs, rhs, ctx->n);
      }
    }
    mpz_mod(rhs, rhs, ctx->n);

    mpz_sub(delta, lhs, rhs);
    mpz_gcd(divisor, delta, ctx->n);
    if (mpz_cmp_ui(divisor, 1) <= 0 || mpz_cmp(divisor, ctx->n) == 0) {
      mpz_add(delta, lhs, rhs);
      mpz_gcd(divisor, delta, ctx->n);
    }
    if (mpz_cmp_ui(divisor, 1) > 0 && mpz_cmp(divisor, ctx->n) < 0 &&
        siqs_insert_divisor(ctx->result, divisor)) {
      ctx->factor_found = 1;
      siqs_verify_partition(ctx->original_n, ctx->result);
    }
    siqs_reset_touched_factors(ctx);
  }
  mpz_clear(lhs);
  mpz_clear(rhs);
  mpz_clear(power);
  mpz_clear(delta);
  mpz_clear(divisor);
  return ctx->factor_found;
}

static int siqs_use_dense_solver(unsigned long ncols) {
  return SIQS_DENSE_SOLVER_MAX_COLS != 0 &&
         ncols <= SIQS_DENSE_SOLVER_MAX_COLS;
}

static int siqs_solve(siqs_ctx_t *ctx) {
  unsigned long original_cols, nrows, ncols, i;
  uint32_t seed1, seed2;
  uint64_t mask = 0;
  uint64_t *nullrows = NULL;
  la_col_t *columns;
  int dense_selected, dense_result = 0;
  if (ctx->full_count < SIQS_MATRIX_EXTRA_RELS(ctx))
    return 0;
  columns = siqs_build_matrix(ctx, &nrows, &ncols);
  original_cols = ncols;
  reduce_matrix(&nrows, &ncols, columns);
  if (ncols == 0) {
    for (i = 0; i < original_cols; i++)
      free(columns[i].data);
    free(columns);
    return 0;
  }
  dense_selected = siqs_use_dense_solver(ncols);
  if (dense_selected) {
    nullrows = dense_nullspace64(nrows, ncols, columns, &mask);
    dense_result = nullrows != NULL;
  }

  if (dense_result) {
    siqs_test_dependencies(ctx, columns, ncols, nullrows, mask);
    free(nullrows);
    nullrows = NULL;
  }

  /* A valid exact basis can, in principle, consist entirely of trivial
   * congruences.  Randomized Lanczos combinations give the same matrix one
   * more chance before relation collection is resumed; this also handles a
   * dense allocation or verification failure. */
  if (!dense_result || !ctx->factor_found) {
    seed1 = (uint32_t)siqs_rand64(&ctx->la_rng);
    seed2 = (uint32_t)siqs_rand64(&ctx->la_rng);
    nullrows = block_lanczos(nrows, 0, ncols, columns,
                             seed1, seed2, &mask);
  }
  if (nullrows != NULL) {
    siqs_test_dependencies(ctx, columns, ncols, nullrows, mask);
    free(nullrows);
  }
  for (i = 0; i < original_cols; i++)
    free(columns[i].data);
  free(columns);
  return ctx->factor_found;
}

/*----------------------------------------------------------------------------
 * Collection driver and context lifetime
 *----------------------------------------------------------------------------*/

static int siqs_collect_relations(siqs_ctx_t *ctx, siqs_poly_t *poly,
                                  uint32_t target,
                                  uint32_t *next_matrix_check,
                                  uint32_t *family_count,
                                  uint32_t *poly_count) {
  int verbose = get_verbose_level();
  int family_ok, have_next;
  uint32_t check_interval = ctx->params.fb_size / 128;
  uint32_t max_polynomials = ctx->params.fb_size > UINT32_MAX / 128U
                           ? UINT32_MAX
                           : ctx->params.fb_size * 128U;
  uint32_t next_report = ctx->full_count + (target - ctx->full_count) / 20 + 1;
  if (max_polynomials < 1000000U)
    max_polynomials = 1000000U;
  if (check_interval < SIQS_MATRIX_CHECK_MIN)
    check_interval = SIQS_MATRIX_CHECK_MIN;
  if (check_interval > SIQS_MATRIX_CHECK_MAX)
    check_interval = SIQS_MATRIX_CHECK_MAX;

  for (;;) {
    if (ctx->factor_found)
      break;
    if (ctx->full_count >= target)
      break;

    family_ok = siqs_new_family(ctx, poly);
    if (!family_ok) {
      return 0;
    }
    (*family_count)++;
    for (;;) {
      siqs_sieve_polynomial(ctx, poly);
      (*poly_count)++;

      if (ctx->full_count >= *next_matrix_check) {
        uint32_t core_rows, core_cols;
        int ready = siqs_matrix_ready(ctx, &core_rows, &core_cols);
        if (verbose > 3)
          printf("# siqs matrix core %u columns, %u rows%s\n",
                 core_cols, core_rows, ready ? ", ready" : "");
        *next_matrix_check = ctx->full_count + check_interval;
        if (ready) {
          return 1;
        }
      }
      if (verbose > 3 && ctx->full_count >= next_report) {
        printf("# siqs relations %u/%u, raw %llu, polys %u\n",
               ctx->full_count, target,
               (unsigned long long)(ctx->accepted_smooth
                                  + ctx->accepted_one_lp
                                  + ctx->accepted_two_lp),
               *poly_count);
        next_report += (target - next_report) / 16 + 1;
      }
      if (*poly_count >= max_polynomials) {
        return 0;
      }
      if (ctx->full_count >= target || ctx->factor_found)
        break;

      have_next = siqs_next_B(ctx, poly);
      if (!have_next)
        break;
    }

  }
  return ctx->full_count >= target || ctx->factor_found;
}

static void siqs_ctx_init(siqs_ctx_t *ctx, const mpz_t original,
                          const mpz_t n, siqs_factor_array_t *result) {
  uint64_t seed;
  memset(ctx, 0, sizeof(*ctx));
  ctx->original_n = original;
  ctx->result = result;
  mpz_init_set(ctx->n, n);
  mpz_init(ctx->kn);
  mpz_init(ctx->eval.y);
  mpz_init(ctx->eval.q);
  mpz_init(ctx->eval.rest);
  siqs_select_parameters(&ctx->params, n);
  seed = (uint64_t)mpz_fdiv_ui(n, 4294967291UL);
  seed ^= (uint64_t)mpz_fdiv_ui(n, 4294967279UL) << 32;
  seed ^= (uint64_t)ctx->params.bits << 17;
  /* A selection must not depend on how many candidate cofactors a parameter
   * variant happens to split.  Domain-separated streams make fixed-family
   * comparisons use the same A sequence and keep Lanczos independent too. */
  ctx->poly_rng.state = siqs_mix64(seed);
  ctx->cofactor_rng.state = siqs_mix64(
      seed ^ UINT64_C(0x6a09e667f3bcc909));
  ctx->la_rng.state = siqs_mix64(
      seed ^ UINT64_C(0xbb67ae8584caa73b));
  ctx->multiplier = siqs_choose_multiplier(n, ctx->params.fb_size);
  mpz_mul_ui(ctx->kn, n, ctx->multiplier);
  ctx->params.poly_d = mpz_fdiv_ui(ctx->kn, 8) == 1 ? 2U : 1U;
  if (ctx->params.poly_d == 2 && mpz_fdiv_ui(ctx->kn, 8) != 1)
    croak("SIQS: d=2 polynomial requires kN congruent to 1 modulo 8");
}

static uint64_t siqs_scaled_bound(uint64_t base, uint32_t multiplier,
                                  uint64_t maximum) {
  if (multiplier == 0)
    return 0;
  return base > maximum / multiplier
       ? maximum : base * (uint64_t)multiplier;
}

static int siqs_lp_floor_policy_active(const siqs_parameters_t *p) {
  return p->max_large_primes == 2 &&
         (p->lp_policy_multiplier_floor != 0 ||
          p->lp_policy_product_multiplier_floor != 0);
}

static void siqs_set_large_prime_bounds(siqs_ctx_t *ctx) {
  uint64_t pmax = ctx->largest_fb_prime;
  uint64_t pmax2 = pmax * pmax;
  uint64_t automatic_product_bound = ctx->params.smooth_bound;
  uint64_t automatic_large_prime_bound = automatic_product_bound;
  uint32_t one_lp_policy_multiplier =
      ctx->params.one_lp_policy_multiplier;
  int policy_active = siqs_lp_floor_policy_active(&ctx->params);
  uint64_t limit;

  if (ctx->params.max_large_primes == 2)
    automatic_large_prime_bound = automatic_product_bound / pmax;
  if (automatic_large_prime_bound > SIQS_LP_MAX)
    automatic_large_prime_bound = SIQS_LP_MAX;
  ctx->params.large_prime_bound = automatic_large_prime_bound;

  /* The full 250-bit corpus selected K=60/R=60 over the lower automatic
   * ratios.  Apply those floors whenever two-LP mode is active through the
   * measured 300-bit endpoint.  Automatic selection has already risen beyond
   * both floors near 270 bits, so the policy naturally rejoins the original
   * curve without an upper cap; above 300 bits it is wholly automatic. */
  if (policy_active) {
    limit = siqs_scaled_bound(pmax,
                              ctx->params.lp_policy_multiplier_floor,
                              SIQS_LP_MAX);
    if (ctx->params.large_prime_bound < limit)
      ctx->params.large_prime_bound = limit;
    limit = siqs_scaled_bound(
        pmax2, ctx->params.lp_policy_product_multiplier_floor,
        SIQS_RESIDUAL_PRODUCT_MAX);
    if (ctx->params.smooth_bound < limit)
      ctx->params.smooth_bound = limit;
  }

  if (one_lp_policy_multiplier != 0) {
    limit = siqs_scaled_bound(pmax, one_lp_policy_multiplier, SIQS_LP_MAX);
    ctx->params.large_prime_bound = limit;
    ctx->params.smooth_bound = limit;
  }

  /* The exact candidate postfilter must accept the complete residual range. */
  ctx->params.sieve_hit_bound = ctx->params.sieve_hit_bound_nominal;
  if (ctx->params.sieve_hit_bound <
      ctx->params.sieve_hit_residual_multiplier *
        (double)ctx->params.smooth_bound)
    ctx->params.sieve_hit_bound =
        ctx->params.sieve_hit_residual_multiplier *
        (double)ctx->params.smooth_bound;
}

static int siqs_ctx_allocate(siqs_ctx_t *ctx) {
  uint32_t rows;
  if (!siqs_build_factor_base(ctx))
    return 0;
  siqs_set_large_prime_bounds(ctx);
  ctx->prime = (uint32_t *)siqs_malloc(
      (size_t)ctx->params.fb_size * sizeof(uint32_t));
  ctx->root1 = (uint32_t *)siqs_malloc(
      (size_t)ctx->params.fb_size * sizeof(uint32_t));
  ctx->root2 = (uint32_t *)siqs_malloc(
      (size_t)ctx->params.fb_size * sizeof(uint32_t));
  siqs_set_log_weights(ctx);
  ctx->sieve_length = 2 * ctx->params.half_interval;
  ctx->active_sieve_length = ctx->sieve_length;
  ctx->sieve = (uint8_t *)siqs_malloc(ctx->sieve_length + 8);
  ctx->candidate_at = (uint16_t *)siqs_calloc(ctx->sieve_length,
                                              sizeof(uint16_t));
  ctx->fb_reciprocal = (uint32_t *)siqs_malloc(
      (size_t)ctx->params.fb_size * sizeof(uint32_t));
  {
    uint32_t i;
    for (i = 0; i < ctx->params.fb_size; i++) {
      ctx->prime[i] = ctx->fb[i].p;
      ctx->fb_reciprocal[i] =
          (uint32_t)(UINT64_C(0x100000000) / ctx->prime[i]);
    }
  }
  ctx->use_block_sieve =
      ctx->sieve_length >= SIQS_SIEVE_BLOCK_MIN &&
      ctx->params.fb_size < SIQS_BUCKET_FB_LIMIT;
  if (ctx->use_block_sieve) {
    ctx->block_count = (ctx->sieve_length + SIQS_SIEVE_BLOCK_SIZE - 1U)
                     >> SIQS_SIEVE_BLOCK_BITS;
    while (ctx->block_large_index < ctx->params.fb_size &&
           ctx->prime[ctx->block_large_index] <= SIQS_SIEVE_BLOCK_SIZE)
      ctx->block_large_index++;
    ctx->bucket_bounds = (uint32_t *)siqs_calloc(
        (size_t)ctx->block_count + 1, sizeof(uint32_t));
    ctx->bucket_fill = (uint32_t *)siqs_malloc(
        (size_t)ctx->block_count * sizeof(uint32_t));
  }
  rows = ctx->params.fb_size + 1;
  ctx->factor_counts = (uint32_t *)siqs_calloc(rows, sizeof(uint32_t));
  ctx->factor_touched = (uint32_t *)siqs_malloc(
      (size_t)rows * sizeof(uint32_t));
  if (siqs_use_one_lp_relation_path(ctx))
    siqs_one_lp_init(&ctx->one_lp);
  else
    siqs_graph_init(&ctx->graph);
  siqs_hashset_init(&ctx->relation_hashes, ctx->params.fb_size * 4);
  siqs_hashset_init(&ctx->a_hashes, 1024);
  return 1;
}

static void siqs_ctx_clear(siqs_ctx_t *ctx) {
  uint32_t i;
  for (i = 0; i < ctx->raw_count; i++)
    siqs_raw_relation_free(ctx->raw[i]);
  for (i = 0; i < ctx->full_count; i++)
    siqs_full_relation_free(ctx->full[i]);
  free(ctx->raw);
  free(ctx->full);
  free(ctx->fb);
  free(ctx->prime);
  free(ctx->root1);
  free(ctx->root2);
  free(ctx->fb_reciprocal);
  free(ctx->sieve);
  free(ctx->bucket_bounds);
  free(ctx->bucket_fill);
  free(ctx->bucket_events);
  free(ctx->candidate_at);
  free(ctx->candidate_at_wide);
  free(ctx->candidates);
  free(ctx->hits);
  free(ctx->eval.factors);
  free(ctx->factor_counts);
  free(ctx->factor_touched);
  siqs_one_lp_clear(&ctx->one_lp);
  siqs_graph_clear(&ctx->graph);
  siqs_hashset_clear(&ctx->relation_hashes);
  siqs_hashset_clear(&ctx->a_hashes);
  mpz_clear(ctx->eval.y);
  mpz_clear(ctx->eval.q);
  mpz_clear(ctx->eval.rest);
  mpz_clear(ctx->n);
  mpz_clear(ctx->kn);
  memset(ctx, 0, sizeof(*ctx));
}


static int siqs_run(siqs_ctx_t *ctx) {
  siqs_poly_t poly;
  uint32_t family_count = 0, poly_count = 0;
  uint32_t target = ctx->params.target_relations;
  /* Small exact matrices need only a modest supply of new dependencies after
   * a rare trivial first result.  FB/16 with an eight-relation floor is
   * cheaper than FB/4; only two of 90,000 audited low-bit inputs needed a
   * third small batch. */
  uint32_t retry_batch = ctx->params.fb_size / 16;
  uint32_t next_matrix_check = ctx->params.fb_size
                             - ctx->params.fb_size / 4;
  int verbose = get_verbose_level();
  if (retry_batch < 8)
    retry_batch = 8;
  if (retry_batch > SIQS_MATRIX_RETRY_BATCH_MAX)
    retry_batch = SIQS_MATRIX_RETRY_BATCH_MAX;
  if (next_matrix_check < 256)
    next_matrix_check = 256;
  if (!siqs_ctx_allocate(ctx)) {
    return ctx->factor_found;
  }
  siqs_poly_init(ctx, &poly);

  if (verbose > 2) {
    gmp_printf("# siqs trying %Zd (%u bits)\n", ctx->n, ctx->params.bits);
    printf("# siqs policy %s [%u-%u], %u-LP\n",
           ctx->params.policy_name, ctx->params.policy_first_bits,
           ctx->params.policy_last_bits, ctx->params.max_large_primes);
    printf("# siqs mult %lu, FB %u (pmax %u), M %u, q %u, d %u, "
           "LP %llu, residual %llu, "
           "pmin %u (init %u)\n",
           ctx->multiplier, ctx->params.fb_size, ctx->largest_fb_prime,
           ctx->params.half_interval, ctx->params.q_count,
           ctx->params.poly_d,
           (unsigned long long)ctx->params.large_prime_bound,
           (unsigned long long)ctx->params.smooth_bound,
           ctx->fb[ctx->params.sieve_start].p,
           (unsigned)ctx->params.sieve_initial);
    if (ctx->use_block_sieve)
      printf("# siqs block sieve %u-byte blocks, %u blocks, "
             "large FB index %u\n",
             SIQS_SIEVE_BLOCK_SIZE, ctx->block_count,
             ctx->block_large_index);
  }

  while (!ctx->factor_found &&
         target <= ctx->params.fb_size + 1 + SIQS_MAX_EXTRA_RELS) {
    if (!siqs_collect_relations(ctx, &poly, target, &next_matrix_check,
                                &family_count, &poly_count))
      break;
    if (ctx->factor_found)
      break;
    if (verbose > 2)
      printf("# siqs linear algebra with %u relations\n", ctx->full_count);
    if (siqs_solve(ctx))
      break;
    target += retry_batch;
    next_matrix_check = ctx->full_count + retry_batch;
  }
  if (verbose > 2)
    printf("# siqs used %u families, %u polynomials, %llu candidates, "
           "%u full relations\n",
           family_count, poly_count,
           (unsigned long long)ctx->total_candidates, ctx->full_count);
  if (verbose > 2)
    printf("# siqs accepted %llu smooth, %llu one-LP, %llu two-LP\n",
           (unsigned long long)ctx->accepted_smooth,
           (unsigned long long)ctx->accepted_one_lp,
           (unsigned long long)ctx->accepted_two_lp);
  if (verbose > 2 && ctx->split_attempts != 0)
    printf("# siqs split %llu composites: %llu squares, %llu SQUFOF, "
           "%llu rho, %llu rejected\n",
           (unsigned long long)ctx->split_attempts,
           (unsigned long long)ctx->split_squares,
           (unsigned long long)ctx->split_squfof,
           (unsigned long long)ctx->split_rho,
           (unsigned long long)ctx->split_failures);
  siqs_poly_clear(ctx, &poly);
  return ctx->factor_found;
}

/*----------------------------------------------------------------------------
 * Public entry point
 *----------------------------------------------------------------------------*/

mpz_t *_GMP_siqs(const mpz_t n, uint32_t *nfactors,
                  uint32_t trial_start) {
  siqs_factor_array_t result;
  siqs_ctx_t ctx;
  mpz_t work, divisor, root;
  uint32_t bits;

  if (nfactors == NULL)
    croak("SIQS: missing factor count output");
  siqs_factor_array_init(&result, n);
  if (mpz_cmp_ui(n, 1) <= 0) {
    *nfactors = result.count;
    return result.values;
  }
  mpz_init_set(work, n);
  mpz_init(divisor);
  mpz_init(root);

  if (trial_start < SIQS_TRIAL_LIMIT) {
    UV p, first = trial_start < 2 ? 2 : trial_start;
    PRIME_ITERATOR(iter);
    prime_iterator_setprime(&iter, first - 1);
    for (p = prime_iterator_next(&iter); p < SIQS_TRIAL_LIMIT;
         p = prime_iterator_next(&iter)) {
      if (mpz_cmp_ui(work, p * p) < 0)
        break;
      if (mpz_divisible_ui_p(work, p)) {
        mpz_set_ui(divisor, p);
        siqs_insert_divisor(&result, divisor);
        mpz_remove(work, work, divisor);
      }
    }
    prime_iterator_destroy(&iter);
    siqs_verify_partition(n, &result);
  }

  if (mpz_cmp_ui(work, 1) == 0 || siqs_all_factors_prime(&result))
    goto finish;
  if (_GMP_is_prob_prime(work))
    goto finish;

  bits = (uint32_t)mpz_sizeinbase(work, 2);
  if (bits < MPU_SIQS_MIN_BITS || bits > MPU_SIQS_MAX_BITS)
    goto finish;

  /* Prime powers yield only trivial QS dependencies.  Return a useful split. */
  if (mpz_perfect_power_p(work)) {
    uint32_t exponent;
    for (exponent = 2; exponent <= bits; exponent++)
      if (mpz_root(root, work, exponent)) {
        siqs_insert_divisor(&result, root);
        siqs_verify_partition(n, &result);
        goto finish;
      }
  }

  siqs_ctx_init(&ctx, n, work, &result);
  /* If the selected square-free multiplier completes a square, the zero of
   * the corresponding polynomial gives a factor directly and would receive
   * an unbounded number of byte-sieve hits. */
  if (mpz_perfect_square_p(ctx.kn)) {
    mpz_sqrt(root, ctx.kn);
    mpz_gcd(divisor, root, work);
    if (mpz_cmp_ui(divisor, 1) > 0 && mpz_cmp(divisor, work) < 0) {
      siqs_insert_divisor(&result, divisor);
      siqs_ctx_clear(&ctx);
      siqs_verify_partition(n, &result);
      goto finish;
    }
  }
  (void)siqs_run(&ctx);
  siqs_ctx_clear(&ctx);
  siqs_verify_partition(n, &result);

finish:
  mpz_clear(work);
  mpz_clear(divisor);
  mpz_clear(root);
  *nfactors = result.count;
  return result.values;
}
