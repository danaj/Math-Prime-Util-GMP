#include <gmp.h>
#include <stdlib.h>
#include "ptypes.h"

#include "factor.h"
#include "primality.h"
#include "prime_iterator.h"
#include "utility.h"
#include "tinyecm64.h"
#include "squfof126.h"
#include "ecm.h"
#include "siqs.h"
#include "lucas_seq.h"

static unsigned long _gcd_ui(unsigned long a, unsigned long b) {
  while (b != 0) {
    unsigned long t = a % b;
    a = b;
    b = t;
  }
  return a;
}

#define _GMP_ECM_FACTOR(n, f, b1, ncurves) \
   _GMP_ecm_factor_projective(n, f, b1, 0, ncurves, irand_u64(), get_verbose_level())

#define NPRIMES_SMALL 6500
static unsigned short primes_small[NPRIMES_SMALL];
static mpz_t _gcd_1k;
static mpz_t _gcd_4k;
static mpz_t _gcd_16k;
static mpz_t _gcd_32k;
static mpz_t _gcd_64k;
void _init_factor(void) {
  uint32_t pn;
  unsigned long nprimes, *pr;

  mpz_init_set_ui(_gcd_1k, 1);
  mpz_init_set_ui(_gcd_4k, 1);
  mpz_init_set_ui(_gcd_16k, 1);
  mpz_init_set_ui(_gcd_32k, 1);
  mpz_init_set_ui(_gcd_64k, 1);

  pr = sieve_to_n_ui(65063, &nprimes);

  /* fill in small primes static */
  MPUassert(nprimes >= NPRIMES_SMALL, "Not enough primes generated in GMP init_factor");
  primes_small[0] = 0;
  for (pn = 1; pn < NPRIMES_SMALL; pn++)
    primes_small[pn] = pr[pn-1];

  /* Construct products for use with gcd to do bulk divisibility testing. */
  mpz_product_ui( _gcd_1k, pr +    1,  167);
  mpz_product_ui( _gcd_4k, pr +  168,  550 -  168);
  mpz_product_ui(_gcd_16k, pr +  550, 1862 -  550);
  mpz_product_ui(_gcd_32k, pr + 1862, 3432 - 1862);
  mpz_product_ui(_gcd_64k, pr + 3432, 6413 - 3432);

#if 0
  for (pn = 0; pn < NPRIMES_SMALL; pn++) {
    UV p = pr[pn];
    if (p > 2 && p < 1000 && !mpz_divisible_ui_p(_gcd_1k,p)) croak("1k");
    if (p >= 1000 && p <= 4000 && !mpz_divisible_ui_p(_gcd_4k,p)) croak("4k");
    if (p >= 4000 && p <= 16000 && !mpz_divisible_ui_p(_gcd_16k,p)) croak("16k");
    if (p >= 16000 && p <= 32000 && !mpz_divisible_ui_p(_gcd_32k,p)) croak("32k");
    if (p >= 32000 && p <= 64000 && !mpz_divisible_ui_p(_gcd_64k,p)) croak("64k");
  }
#endif
  Safefree(pr);
}

void _destroy_factor(void) {
  mpz_clear(_gcd_1k);
  mpz_clear(_gcd_4k);
  mpz_clear(_gcd_16k);
  mpz_clear(_gcd_32k);
  mpz_clear(_gcd_64k);
}

/* Max number of factors on the unfactored stack, not the max total factors.
 * This is used when we split n into two or more composites.  Since we work
 * on the smaller of the composites first, this rarely goes above 10 even
 * with thousands of non-trivial factors. */

static int _mpz_cmp_asc(const void *av, const void *bv)
  { return mpz_cmp(*(const mpz_t*)av, *(const mpz_t*)bv); }

enum {
  FACTOR_METHOD_EARLY_PM1,
  FACTOR_METHOD_PM1,
  FACTOR_METHOD_ECM,
  FACTOR_METHOD_ECM_SCALED,
  FACTOR_METHOD_LARGE_ECM,
  FACTOR_METHOD_SIQS,
  FACTOR_METHOD_PBRENT,
  FACTOR_METHOD_HOLF,
  FACTOR_METHOD_CHEB
};

typedef struct {
  unsigned int bits;
  int method;
  UV arg1;
  UV arg2;
  const char *description;
} factor_step_t;

/* Pretests: BITS is a minimum; zero means unrestricted.  These are
 * cumulative, deliberately modest budgets before the first SIQS checkpoint.
 * Early p-1 continues while it splits, stopping for this factor() call
 * after its first miss. Other steps apply to each unresolved cofactor.
 * B1=600 uses four curves from 160 bits, plus twelve more from 176 bits.
 * Shorter stage-2 budgets are used through 224 bits; larger inputs retain
 * their previous bounds, curve counts and execution order. */
static const factor_step_t factor_pretests[] = {
  /* bits  method                      B1       B2/curves */
  {   0, FACTOR_METHOD_EARLY_PM1,       0,       0, NULL                    },
  { 128, FACTOR_METHOD_ECM,           200,       4, "ecm (200)"        },
  { 160, FACTOR_METHOD_ECM,           600,       4, "ecm (600)"        },
  { 176, FACTOR_METHOD_ECM,           600,      12, "ecm (600)"        },
  { 192, FACTOR_METHOD_ECM,          2000,      10, "ecm (2000)"       },
  { 208, FACTOR_METHOD_PM1,        200000, 3000000, "p-1 (200k)"            },
  { 208, FACTOR_METHOD_ECM,         20000,       2, "ecm (20k,2)"     },
  { 240, FACTOR_METHOD_ECM,         30000,      20, "ecm (30k,20)"    },
  { 256, FACTOR_METHOD_ECM,         40000,      40, "ecm (40k,40)"    },
  { 270, FACTOR_METHOD_ECM,         80000,      40, "ecm (80k,40)"    },
  { 288, FACTOR_METHOD_ECM,        160000,      80, "ecm (160k,80)"   }
};

/* Main ladder: BITS is a maximum; zero means unrestricted.  Method support
 * limits still apply.  ECM_SCALED/LARGE_ECM use (base-B1 multiplier, curves),
 * PM1 uses (B1, B2), PBRENT uses (a, rounds), HOLF uses (rounds, 0), and CHEB
 * uses (bound, initial x).  SIQS needs no additional arguments.
 *
 * Exactly how much ECM should be done before SIQS is a heuristic.
 * If more information is known about the input shape, directly calling
 * ecm_factor and qs_factor can be faster (e.g. if the input is known to be
 * a large balanced semiprime, then ECM wastes time, while if we know there
 * are smaller factors then deeper ECM is useful).
 * ECM_SCALED/LARGE_ECM multiply arg1 by factor_ecm_base_bound(bits).
 *
 * SIQS runs at most once for each cofactor; misses stay on the ladder. */
static const factor_step_t factor_ladder[] = {
  /* bits  method                     arg1      arg2 */
  { 320, FACTOR_METHOD_SIQS,              0,         0, "SIQS"             },
  {   0, FACTOR_METHOD_ECM_SCALED,        2,        20, NULL               },
  { 360, FACTOR_METHOD_SIQS,              0,         0, "SIQS"             },
  {   0, FACTOR_METHOD_PBRENT,            1,   1048576, "pbrent (1,1M)"    },
  {   0, FACTOR_METHOD_ECM_SCALED,        4,        20, NULL               },
  {   0, FACTOR_METHOD_ECM_SCALED,        8,        20, NULL               },
  { 400, FACTOR_METHOD_SIQS,              0,         0, "SIQS"             },
  {   0, FACTOR_METHOD_HOLF,        1048576,         0, "holf"             },
  {   0, FACTOR_METHOD_CHEB,          50000,         0, "cheb (50k)"       },
  {   0, FACTOR_METHOD_PM1,         5000000, 100000000, "p-1 (5M)"         },
  /* Above 400 bits SIQS is not end-to-end validated and can use substantial
   * memory/time.  Give ECM more opportunity before trying the supported range. */
  {   0, FACTOR_METHOD_ECM_SCALED,       32,        40, NULL               },
  {   0, FACTOR_METHOD_LARGE_ECM,         8,       100, NULL               },
  {   0, FACTOR_METHOD_ECM_SCALED,       16,       100, NULL               },
  {   0, FACTOR_METHOD_ECM_SCALED,       32,       100, NULL               },
  {   0, FACTOR_METHOD_SIQS,              0,         0, "SIQS"             },
  {   0, FACTOR_METHOD_ECM_SCALED,       64,       100, NULL               },
  {   0, FACTOR_METHOD_ECM_SCALED,      128,       100, NULL               },
  /* More ECM work could be done but at this point we are the wrong tool. */
};

static UV factor_early_pm1_bound(UV nbits)
{
  /* Scale or skip this fixed work where SIQS itself is very fast. */
  return nbits <  88 ?     0
       : nbits < 112 ?  1000
       : nbits < 176 ?  2000
                     : 20000;
}

static UV factor_ecm_base_bound(UV nbits)
{
  return nbits <  65 ?   5000
       : nbits < 160 ?  20000
       : nbits < 192 ?  30000
       : nbits < 224 ?  40000
       : nbits < 256 ?  80000
       : nbits < 512 ? 160000
                     : 320000;
}

/* arg1/arg2 are resolved call parameters, not the table's ECM multipliers.
 * Return a table label or format a dynamic label in the caller's buffer. */
static const char *factor_step_name(const factor_step_t *step, UV arg1, UV arg2,
                                    char *buffer, size_t size)
{
  if (step->method == FACTOR_METHOD_EARLY_PM1)
    gmp_snprintf(buffer, size, "p-1 (%luk)", (unsigned long)(arg1/1000));
  else if (step->method == FACTOR_METHOD_ECM_SCALED ||
           step->method == FACTOR_METHOD_LARGE_ECM)
    gmp_snprintf(buffer, size, "ecm (%luk,%lu)",
                 (unsigned long)(arg1/1000), (unsigned long)arg2);
  else
    return step->description;
  return buffer;
}

enum { FACTOR_UNKNOWN, FACTOR_COMPOSITE, FACTOR_PRIME };

typedef struct {
  mpz_t value;
  int exponent;
  unsigned char state;
} factor_piece_t;

/* A coprime multiplicative partition: product(value^exponent) is the input.
 * COMPOSITE means primality and power checks have already been done for
 * this value. Changing a value invalidates both cached classifications.
 * Entries may move on growth/removal; callers retain indices, not pointers. */
typedef struct {
  factor_piece_t *pieces;
  size_t count, alloc;
  mpz_t scratch;
  const char *error;
} factor_partition_t;

static void factor_partition_init(factor_partition_t *part)
{
  part->pieces = NULL;
  part->count = part->alloc = 0;
  part->error = NULL;
  mpz_init(part->scratch);
}

static void factor_partition_clear(factor_partition_t *part)
{
  while (part->count) mpz_clear(part->pieces[--part->count].value);
  Safefree(part->pieces);
  mpz_clear(part->scratch);
}

static int factor_partition_append(factor_partition_t *part, const mpz_t value,
                                   int exponent, unsigned char state)
{
  if (part->count == part->alloc) {
    size_t maximum = SIZE_MAX / sizeof(factor_piece_t), capacity;
    if (maximum > INT_MAX) maximum = INT_MAX;
    if (part->count == maximum) {
      part->error = "Too many factors";
      return 0;
    }
    capacity = part->alloc == 0 ? 4
             : part->alloc > maximum/2 ? maximum : 2*part->alloc;
    if (capacity > maximum) capacity = maximum;
    if (part->alloc == 0) New(0, part->pieces, capacity, factor_piece_t);
    else                 Renew(part->pieces, capacity, factor_piece_t);
    part->alloc = capacity;
  }
  mpz_init_set(part->pieces[part->count].value, value);
  part->pieces[part->count].exponent = exponent;
  part->pieces[part->count++].state = state;
  return 1;
}

static void factor_partition_remove(factor_partition_t *part, size_t index)
{
  size_t last = --part->count;
  if (index != last) {
    mpz_swap(part->pieces[index].value, part->pieces[last].value);
    part->pieces[index].exponent = part->pieces[last].exponent;
    part->pieces[index].state = part->pieces[last].state;
  }
  mpz_clear(part->pieces[last].value);
}

/* Merge equal pieces and split overlaps until all values are coprime.
 * A whole-piece divisor is removed to its full valuation at once: this
 * handles large powers without one entry or one gcd pass per copy. */
static int factor_partition_refine(factor_partition_t *part)
{
  size_t i, j;
  int a, b, ag, bg;
  mp_bitcnt_t copies;
again:
  for (i = 0; i < part->count; i++) {
    for (j = i+1; j < part->count; j++) {
      a = part->pieces[i].exponent;
      b = part->pieces[j].exponent;
      if (mpz_cmp(part->pieces[i].value, part->pieces[j].value) == 0) {
        if (a > INT_MAX-b) goto overflow;
        part->pieces[i].exponent = a+b;
        if (part->pieces[i].state == FACTOR_UNKNOWN)
          part->pieces[i].state = part->pieces[j].state;
        factor_partition_remove(part, j);
        goto again;
      }
      if (part->pieces[i].state == FACTOR_PRIME &&
          part->pieces[j].state == FACTOR_PRIME) continue;
      mpz_gcd(part->scratch, part->pieces[i].value, part->pieces[j].value);
      if (mpz_cmp_ui(part->scratch, 1) == 0) continue;
      ag = mpz_cmp(part->scratch, part->pieces[i].value) == 0;
      bg = mpz_cmp(part->scratch, part->pieces[j].value) == 0;
      if (ag) {
        copies = mpz_remove(part->pieces[j].value, part->pieces[j].value,
                            part->pieces[i].value);
        if (copies > (mp_bitcnt_t)((INT_MAX-a)/b)) goto overflow;
        part->pieces[i].exponent = a + b*(int)copies;
        part->pieces[j].state = FACTOR_UNKNOWN;
        if (mpz_cmp_ui(part->pieces[j].value, 1) == 0)
          factor_partition_remove(part, j);
      } else if (bg) {
        copies = mpz_remove(part->pieces[i].value, part->pieces[i].value,
                            part->pieces[j].value);
        if (copies > (mp_bitcnt_t)((INT_MAX-b)/a)) goto overflow;
        part->pieces[j].exponent = b + a*(int)copies;
        part->pieces[i].state = FACTOR_UNKNOWN;
        if (mpz_cmp_ui(part->pieces[i].value, 1) == 0)
          factor_partition_remove(part, i);
      } else {
        if (a > INT_MAX-b) goto overflow;
        mpz_divexact(part->pieces[i].value, part->pieces[i].value, part->scratch);
        mpz_divexact(part->pieces[j].value, part->pieces[j].value, part->scratch);
        part->pieces[i].state = part->pieces[j].state = FACTOR_UNKNOWN;
        if (!factor_partition_append(part, part->scratch, a+b, FACTOR_UNKNOWN))
          return 0;
      }
      goto again;
    }
  }
  return 1;
overflow:
  part->error = "factor() exponent overflow";
  return 0;
}

/* Classify each new value once and choose the smallest unresolved piece.
 * Power roots stay in the partition with a multiplied exponent. Any exact
 * prime-exponent root suffices; the next pass handles powers in that root.
 * This does not need power_factor's maximal-exponent search or recursion.
 * prime_cutoff is the square of a known lower bound on prime factors of
 * unresolved pieces (zero if none). Smaller pieces are therefore prime. */
static int factor_partition_prepare(factor_partition_t *part, mpz_t root,
                                    size_t *selected, unsigned long prime_cutoff)
{
  size_t i, bits;
  unsigned long degree, limit;
  *selected = part->count;
  for (i = 0; i < part->count; i++) {
    while (part->pieces[i].state == FACTOR_UNKNOWN) {
      if (prime_cutoff != 0 && mpz_cmp_ui(part->pieces[i].value, prime_cutoff) < 0) {
        part->pieces[i].state = FACTOR_PRIME;
      } else if (mpz_perfect_power_p(part->pieces[i].value)) {
        PRIME_ITERATOR(iter);
        bits = mpz_sizeinbase(part->pieces[i].value, 2);
        limit = (unsigned long)(INT_MAX / part->pieces[i].exponent);
        if (bits-1 < limit) limit = (unsigned long)(bits-1);
        for (degree = 2; degree <= limit; degree = prime_iterator_next(&iter))
          if (mpz_root(root, part->pieces[i].value, degree)) break;
        prime_iterator_destroy(&iter);
        if (degree > limit) {
          part->error = "factor() exponent overflow";
          return 0;
        }
        part->pieces[i].exponent *= (int)degree;
        mpz_set(part->pieces[i].value, root);
        /* A root divides the old value, so coprimality with other entries
         * is preserved. Its primality and remaining power are still unknown. */
      } else {
        part->pieces[i].state = _GMP_is_prob_prime(part->pieces[i].value)
                             ? FACTOR_PRIME : FACTOR_COMPOSITE;
      }
    }
    if (part->pieces[i].state == FACTOR_COMPOSITE &&
        (*selected == part->count ||
         mpz_cmp(part->pieces[i].value, part->pieces[*selected].value) < 0))
      *selected = i;
  }
  return 1;
}

/* d is external scratch, not an entry that may move during append. */
static int factor_partition_split(factor_partition_t *part, size_t index,
                                  const mpz_t d)
{
  int exponent = part->pieces[index].exponent;
  if (mpz_cmp_ui(d, 1) <= 0 || mpz_cmp(d, part->pieces[index].value) >= 0 ||
      !mpz_divisible_p(part->pieces[index].value, d)) {
    part->error = "Incorrect factoring";
    return 0;
  }
  mpz_divexact(part->scratch, part->pieces[index].value, d);
  mpz_set(part->pieces[index].value, d);
  part->pieces[index].state = FACTOR_UNKNOWN;
  return factor_partition_append(part, part->scratch, exponent, FACTOR_UNKNOWN)
      && factor_partition_refine(part);
}

/* Consume a verified SIQS partition through the same refinement path.
 * Return 1 on progress, 0 on a miss, -1 on invalid data or size overflow. */
static int factor_partition_siqs(factor_partition_t *part, size_t index,
                                 const mpz_t n, int verbose)
{
  mpz_t *values;
  uint32_t i, count;
  int exponent = part->pieces[index].exponent, ok = 1;
  values = gmp_siqs(n, &count, 64007, verbose - 2);
  if (count < 2) {
    gmp_siqs_free(values, count);
    return 0;
  }
  mpz_set_ui(part->scratch, 1);
  for (i = 0; i < count; i++) {
    if (mpz_cmp_ui(values[i], 1) <= 0) { ok = 0; break; }
    mpz_mul(part->scratch, part->scratch, values[i]);
  }
  if (!ok || mpz_cmp(part->scratch, n) != 0) {
    part->error = "SIQS: invalid result partition";
    gmp_siqs_free(values, count);
    return -1;
  }
  mpz_set(part->pieces[index].value, values[0]);
  part->pieces[index].state = FACTOR_UNKNOWN;
  for (i = 1; ok && i < count; i++)
    ok = factor_partition_append(part, values[i], exponent, FACTOR_UNKNOWN);
  gmp_siqs_free(values, count);
  if (!ok || !factor_partition_refine(part)) return -1;
  if (verbose) gmp_printf("SIQS found partition with %u factors\n", (unsigned)count);
  return 1;
}

/* Initial trial division proves these pieces prime and mutually coprime.
 * Keep the residual in slot zero, including a temporary unit until trial
 * division finishes, so the represented product is preserved throughout. */
static int factor_partition_trial(factor_partition_t *part, const mpz_t n,
                                  const mpz_t p, mp_bitcnt_t copies)
{
  if (copies > INT_MAX) {
    part->error = "factor() exponent overflow";
    return 0;
  }
  mpz_set(part->pieces[0].value, n);
  return factor_partition_append(part, p, (int)copies, FACTOR_PRIME);
}

static int factor_piece_cmp(const void *a, const void *b)
{
  return mpz_cmp(((const factor_piece_t *)a)->value,
                 ((const factor_piece_t *)b)->value);
}

#define FACTOR_TRIAL_SMALL(n, pn_lo, pn_hi) \
  do { unsigned long sp, p; mp_bitcnt_t copies; \
    for (sp = pn_lo, p = primes_small[sp]; \
         sp <= pn_hi && mpz_cmp_ui(n,p*p) >= 0; \
         p = primes_small[++sp]) { \
      if (mpz_divisible_ui_p(n, p)) { \
        mpz_set_ui(f, p); \
        copies = mpz_remove(n, n, f); \
        if (!factor_partition_trial(&part, n, f, copies)) goto FAIL; \
      } \
    } \
  } while (0)

int factor(const mpz_t input_n, mpz_t* pfactors[], int* pexponents[])
{
  factor_partition_t part;
  mpz_t n, f;
  mpz_t *factors = NULL;
  int *exponents = NULL, nfactors, early_pm1_done = 0;
  size_t index, i;
  mp_bitcnt_t copies;
  const char *error;

  *pfactors = NULL;
  *pexponents = NULL;
  mpz_init_set(n, input_n);
  mpz_abs(n, n);
  if (mpz_cmp_ui(n, 4) < 0) {
    nfactors = mpz_cmp_ui(n, 1) != 0;
    if (nfactors) {
      New(0, factors, 1, mpz_t);
      New(0, exponents, 1, int);
      mpz_init_set(factors[0], n);
      exponents[0] = 1;
    }
    mpz_clear(n);
    *pfactors = factors;
    *pexponents = exponents;
    return nfactors;
  }
  mpz_init(f);
  factor_partition_init(&part);
  if (!factor_partition_append(&part, n, 1, FACTOR_UNKNOWN)) goto FAIL;

  if (mpz_even_p(n)) {
    mpz_set_ui(f, 2);
    copies = mpz_remove(n, n, f);
    if (!factor_partition_trial(&part, n, f, copies)) goto FAIL;
  }
  if (mpz_gcd(f,n,_gcd_1k),mpz_cmp_ui(f,1)>0) FACTOR_TRIAL_SMALL(n,   2,  168);
  if (mpz_gcd(f,n,_gcd_4k),mpz_cmp_ui(f,1)>0) FACTOR_TRIAL_SMALL(n, 169,  550);
  if (mpz_gcd(f,n,_gcd_16k),mpz_cmp_ui(f,1)>0)FACTOR_TRIAL_SMALL(n, 551, 1862);
  if (mpz_gcd(f,n,_gcd_32k),mpz_cmp_ui(f,1)>0)FACTOR_TRIAL_SMALL(n,1863, 3432);
  if (mpz_gcd(f,n,_gcd_64k),mpz_cmp_ui(f,1)>0)FACTOR_TRIAL_SMALL(n,3433, 6413);
  if (mpz_cmp_ui(n, 1) == 0)
    factor_partition_remove(&part, 0);
  else if (mpz_cmp_ui(n, 64007UL*64007UL) < 0)
    part.pieces[0].state = FACTOR_PRIME;

  while (1) {
    int success = 0, progress = 0, siqs_attempted = 0;
    int o = get_verbose_level();
    UV B1 = 0, nbits;
    size_t stage, step_index;
    if (!factor_partition_prepare(&part, f, &index, 64007UL*64007UL)) goto FAIL;
    if (index == part.count) break;
    mpz_set(n, part.pieces[index].value);
    nbits = mpz_sizeinbase(n, 2);
    if (o > 1) { gmp_printf("Factoring %Zd...\n", n); fflush(stdout); }

    if (nbits <= 64) {
      success = _GMP_tinyecm64_factor(n, f, 0);
      if (success && o) gmp_printf("TINYECM64 found factor %Zd\n", f);
    }
    for (stage = 0; stage < 2 && !success && !progress; stage++) {
      const factor_step_t *steps = stage == 0 ? factor_pretests : factor_ladder;
      size_t nsteps = stage == 0 ? sizeof(factor_pretests)/sizeof(factor_pretests[0])
                               : sizeof(factor_ladder)/sizeof(factor_ladder[0]);
      for (step_index = 0; step_index < nsteps && !success && !progress; step_index++) {
        const factor_step_t *step = &steps[step_index];
        UV arg1 = step->arg1, arg2 = step->arg2;
        char method_buffer[64];
        const char *method_name = NULL;
        if (stage == 0 && nbits < step->bits) continue;
        if (stage != 0 && step->bits != 0 && nbits > step->bits) continue;
        if (step->method == FACTOR_METHOD_EARLY_PM1) {
          if (early_pm1_done) continue;
          arg1 = factor_early_pm1_bound(nbits);
          if (arg1 == 0) continue;
          arg2 = (nbits <= 224 ? 5 : 10) * arg1;
        } else if (step->method == FACTOR_METHOD_ECM_SCALED ||
                   step->method == FACTOR_METHOD_LARGE_ECM) {
          if (B1 == 0) B1 = factor_ecm_base_bound(nbits);
          arg1 *= B1;
        } else if (step->method == FACTOR_METHOD_SIQS &&
                   (siqs_attempted || nbits < MPU_SIQS_MIN_BITS || nbits > MPU_SIQS_MAX_BITS)) {
          continue;
        }
        if (o) method_name = factor_step_name(step, arg1, arg2, method_buffer, sizeof(method_buffer));
        if (o > 1) { gmp_printf("  Factor trying %s\n", method_name); fflush(stdout); }
        switch (step->method) {
          case FACTOR_METHOD_EARLY_PM1:
          case FACTOR_METHOD_PM1:
            success = _GMP_pminus1_factor(n, f, arg1, arg2);
            if (step->method == FACTOR_METHOD_EARLY_PM1 && !success)
              early_pm1_done = 1;
            break;
          case FACTOR_METHOD_ECM:
            success = _GMP_ecm_factor_projective(n, f, arg1,
                stage == 0 && nbits <= 224 ? 40*arg1 : 0,
                arg2, irand_u64(), get_verbose_level());
            break;
          case FACTOR_METHOD_LARGE_ECM:
            if (o) gmp_printf("starting large ECM on %Zd\n", n);
            /* FALLTHROUGH */
          case FACTOR_METHOD_ECM_SCALED:
            success = _GMP_ECM_FACTOR(n, f, arg1, arg2);
            break;
          case FACTOR_METHOD_PBRENT:
            success = _GMP_pbrent_factor(n, f, arg1, arg2);
            break;
          case FACTOR_METHOD_HOLF:
            success = _GMP_holf_factor(n, f, arg1);
            break;
          case FACTOR_METHOD_CHEB:
            success = _GMP_cheb_factor(n, f, arg1, arg2);
            break;
          case FACTOR_METHOD_SIQS:
            siqs_attempted = 1;
            progress = factor_partition_siqs(&part, index, n, o);
            if (progress < 0) goto FAIL;
            break;
          default:
            part.error = "Unknown factoring method";
            goto FAIL;
        }
        if (success && o) gmp_printf("%s found %Zd\n", method_name, f);
      }
    }
    if (success) {
      if (!factor_partition_split(&part, index, f)) goto FAIL;
    } else if (!progress) {
      gmp_printf("Failed to factor: %Zd\n", n);
      part.error = "factor() unable to factor";
      goto FAIL;
    }
  }

  nfactors = (int)part.count;
  qsort(part.pieces, part.count, sizeof(factor_piece_t), factor_piece_cmp);
  New(0, factors, part.count, mpz_t);
  New(0, exponents, part.count, int);
  /* Transfer limb storage instead of copying each finished prime. */
  for (i = 0; i < part.count; i++) {
    mpz_init(factors[i]);
    mpz_swap(factors[i], part.pieces[i].value);
    exponents[i] = part.pieces[i].exponent;
  }
  factor_partition_clear(&part);
  mpz_clear(f);
  mpz_clear(n);
  *pfactors = factors;
  *pexponents = exponents;
  return nfactors;

FAIL:
  error = part.error;
  factor_partition_clear(&part);
  mpz_clear(f);
  mpz_clear(n);
  croak("%s", error ? error : "factor() unable to factor");
  return 0;
}


void clear_factors(int nfactors, mpz_t* pfactors[], int* pexponents[])
{
  while (nfactors > 0)
    mpz_clear((*pfactors)[--nfactors]);
  Safefree(*pfactors);
  if (pexponents != 0)
    Safefree(*pexponents);
}


/*****************************************************************************/


/* TODO: Optimize to factor down to the last semiprime. */
static void _omega(const mpz_t n, uint32_t* omega, uint32_t* bigomega)
{
  mpz_t* factors;
  int i, nfactors, result, *exponents;

  nfactors = factor(n, &factors, &exponents);
  if (bigomega != 0) {
    for (i = 0, result = 0; i < nfactors; i++)
      result += exponents[i];
    *bigomega = result;
  }
  if (omega != 0)
    *omega = nfactors;
  clear_factors(nfactors, &factors, &exponents);
}

uint32_t omega(const mpz_t n)
{
  uint32_t o;
  _omega(n, &o, 0);
  return o;
}
uint32_t bigomega(const mpz_t n)
{
  uint32_t bo;
  _omega(n, 0, &bo);
  return bo;
}

void sopf(mpz_t res, const mpz_t n)
{
  mpz_t* factors;
  int i, nfactors, *exponents;

  if (mpz_cmp_ui(n, 1) <= 0) {
    mpz_set_ui(res, 0);
    return;
  }

  nfactors = factor(n, &factors, &exponents);
  mpz_set_ui(res, 0);
  for (i = 0; i < nfactors; i++)
    mpz_add(res, res, factors[i]);
  clear_factors(nfactors, &factors, &exponents);
}

void sopfr(mpz_t res, const mpz_t n)
{
  mpz_t* factors;
  int i, j, nfactors, *exponents;

  if (mpz_cmp_ui(n, 1) <= 0) {
    mpz_set_ui(res, 0);
    return;
  }

  nfactors = factor(n, &factors, &exponents);
  mpz_set_ui(res, 0);
  for (i = 0; i < nfactors; i++)
    for (j = 0; j < exponents[i]; j++)
      mpz_add(res, res, factors[i]);
  clear_factors(nfactors, &factors, &exponents);
}

static int _cmp_uint32_desc(const void *a, const void *b)
{
  uint32_t ai = *(const uint32_t*)a;
  uint32_t bi = *(const uint32_t*)b;
  return (ai < bi) ? 1 : (ai > bi) ? -1 : 0;
}

int prime_signature(mpz_t res, uint32_t **signature, const mpz_t n)
{
  mpz_t* factors;
  int i, nfactors, *exponents;
  uint32_t *sig = 0;

  if (mpz_sgn(n) == 0) {
    if (res != 0) mpz_set_ui(res, 2);
    if (signature != 0) {
      New(0, sig, 1, uint32_t);
      sig[0] = 1;
      *signature = sig;
    }
    return 1;
  }
  if (mpz_cmp_ui(n, 1) == 0) {
    if (res != 0) mpz_set_ui(res, 0);
    if (signature != 0) *signature = 0;
    return 0;
  }

  nfactors = factor(n, &factors, &exponents);
  New(0, sig, nfactors, uint32_t);
  for (i = 0; i < nfactors; i++)
    sig[i] = exponents[i];
  qsort(sig, nfactors, sizeof(*sig), _cmp_uint32_desc);

  if (res != 0) {
    mpz_t t;
    PRIME_ITERATOR(iter);

    mpz_init(t);
    mpz_set_ui(res, 1);
    prime_iterator_setprime(&iter, 1);
    for (i = 0; i < nfactors; i++) {
      unsigned long p = prime_iterator_next(&iter);
      mpz_ui_pow_ui(t, p, sig[i]);
      mpz_mul(res, res, t);
    }
    prime_iterator_destroy(&iter);
    mpz_clear(t);
  }

  clear_factors(nfactors, &factors, &exponents);
  if (signature != 0) {
    *signature = sig;
  } else {
    Safefree(sig);
  }
  return nfactors;
}

void sigma(mpz_t res, const mpz_t n, unsigned long k)
{
  mpz_t* factors;
  mpz_t pk, pke, fmult;
  int* exponents;
  int i, j, nfactors;

  if (mpz_cmp_ui(n, 1) <= 0) {
    mpz_set_ui(res, (mpz_sgn(n) > 0));
    return;
  }

  if (_GMP_is_prob_prime(n)) {
    mpz_pow_ui(res, n, k);
    mpz_add_ui(res, res, 1);
    return;
  }

  nfactors = factor(n, &factors, &exponents);
  if (k == 0) {
    for (i = 0; i < nfactors; i++) {
      mpz_set_ui(factors[i], exponents[i]+1);
    }
  } else if (k == 1) {
    mpz_init(pke);
    mpz_init(fmult);
    for (i = 0; i < nfactors; i++) {
      mpz_set(pke, factors[i]);
      mpz_add_ui(fmult, factors[i], 1);
      for (j = 1; j < exponents[i]; j++) {
        mpz_mul(pke, pke, factors[i]);
        mpz_add(fmult, fmult, pke);
      }
      mpz_set(factors[i], fmult);
    }
    mpz_clear(fmult);
    mpz_clear(pke);
  } else {
    mpz_init(pk);
    mpz_init(pke);
    mpz_init(fmult);
    for (i = 0; i < nfactors; i++) {
      mpz_pow_ui(pk, factors[i], k);
      mpz_add_ui(fmult, pk, 1);
      mpz_set(pke, pk);
      for (j = 1; j < exponents[i]; j++) {
        mpz_mul(pke, pke, pk);
        mpz_add(fmult, fmult, pke);
      }
      mpz_set(factors[i], fmult);
    }
    mpz_clear(fmult);
    mpz_clear(pke);
    mpz_clear(pk);
  }
  mpz_product(factors, 0, nfactors-1);
  mpz_set(res, factors[0]);
  clear_factors(nfactors, &factors, &exponents);
}

void dedekind_psi(mpz_t res, const mpz_t n)
{
  mpz_t* factors;
  mpz_t t, pke;
  int i, j, nfactors, *exponents;

  if (mpz_cmp_ui(n, 1) <= 0) {
    mpz_set_ui(res, (mpz_cmp_ui(n, 1) == 0) ? 1 : 0);
    return;
  }

  nfactors = factor(n, &factors, &exponents);
  mpz_init(t);
  mpz_init(pke);
  mpz_set_ui(res, 1);
  for (i = 0; i < nfactors; i++) {
    mpz_add_ui(t, factors[i], 1);
    if (exponents[i] > 1) {
      mpz_set(pke, factors[i]);
      for (j = 2; j < exponents[i]; j++)
        mpz_mul(pke, pke, factors[i]);
      mpz_mul(t, t, pke);
    }
    mpz_mul(res, res, t);
  }
  mpz_clear(pke);
  mpz_clear(t);
  clear_factors(nfactors, &factors, &exponents);
}

void aliquot_sum(mpz_t res, const mpz_t n)
{
  if (mpz_cmp_ui(n, 1) <= 0) {
    mpz_set_ui(res, 0);
  } else {
    mpz_t r;
    mpz_init(r);
    sigma(r, n, 1);
    mpz_sub(res, r, n);
    mpz_clear(r);
  }
}

void abundance(mpz_t res, const mpz_t n)
{
  if (mpz_sgn(n) <= 0) {
    mpz_set_ui(res, 0);
  } else {
    mpz_t r;
    mpz_init(r);
    sigma(r, n, 1);
    mpz_sub(r, r, n);
    mpz_sub(res, r, n);
    mpz_clear(r);
  }
}


static const signed char _smallmoebius[64] = {0,1,-1,-1,0,-1,1,-1,0,0,1,-1,0,-1,1,1,0,-1,0,-1,0,1,1,-1,0,0,1,0,0,-1,-1,-1,0,1,1,1,0,-1,1,1,0,-1,-1,-1,0,0,1,-1,0,0,0,1,0,-1,0,1,0,1,1,-1,0,-1,1,0};
static const unsigned long _smalldiv[] = {4, 9, 25, 49, 121, 169, 289};
int moebius(const mpz_t n)
{
  uint32_t i, o, bo;

  if (mpz_sgn(n) < 0) {
    mpz_t t;
    mpz_init(t);
    mpz_neg(t,n);
    o = moebius(t);
    mpz_clear(t);
    return o;
  }
  if (mpz_cmp_ui(n,64) < 0) return _smallmoebius[mpz_get_ui(n)];

  for (i = 0; i < 7; i++)
    if (mpz_divisible_ui_p(n, _smalldiv[i]))
      return 0;

  _omega(n, &o, &bo);
  return (o != bo) ?   0
      :  (o % 2)   ?  -1
                   :   1;
}

/* mpu 'for (0..256) { next if moebius($_) == 0; $b[$_ >> 5] |= (1 << ($_%32)); } say join ",",@b;' */
static const uint32_t _isf[8] = {3840601326U,1856556782U,3941394158U,2362371810U,3970362990U,3471729898U,4008603310U,3938642668U};
int is_square_free_ui(unsigned long n)
{
  uint32_t i, p, psq;

  if (n < 256)  return (_isf[n >> 5] & (1U << (n % 32))) != 0;

  for (i = 1; i < NPRIMES_SMALL; i++) {
    p = primes_small[i];
    if (p >= 65536) break;
    psq = p * p;
    if (psq > n)  return 1;
    if ((n % psq) == 0) return 0;
  }
  {
    mpz_t t;
    mpz_init_set_ui(t, n);
    i = (moebius(t) != 0);
    mpz_clear(t);
    return i;
  }
}
int is_square_free(const mpz_t n)
{
  if (mpz_cmp_ui(n, 1000000U) <= 0)
    return is_square_free_ui(mpz_get_ui(n));
  return moebius(n) != 0;
}

int liouville(const mpz_t n)
{
  uint32_t result = bigomega(n);
  return (result & 1)  ?  -1  : 1;
}

int is_safe_prime(const mpz_t n)
{
  mpz_t q;
  int ret;

  if (mpz_sgn(n) <= 0 || !_GMP_is_prime(n))
    return 0;

  mpz_init(q);
  mpz_sub_ui(q, n, 1);
  mpz_tdiv_q_2exp(q, q, 1);
  ret = _GMP_is_prime(q);
  mpz_clear(q);
  return ret ? 1 : 0;
}

void totient(mpz_t tot, const mpz_t n_input)
{
  mpz_t t, n;
  mpz_t* factors;
  int* exponents;
  int i, j, nfactors;

  if (mpz_sgn(n_input) <= 0)
    { mpz_set_ui(tot, 0); return; }
  if (mpz_cmp_ui(n_input, 1) == 0)
    { mpz_set_ui(tot, 1); return; }

  mpz_init_set(n, n_input);
  mpz_set_ui(tot, 1);
  /* Fast reduction of multiples of 2 */
  i = mpz_scan1(n, 0);
  if (i > 0) {
    if (i > 1)  mpz_mul_2exp(tot, tot, i-1);
    mpz_tdiv_q_2exp(n, n, i);
  }
  /* Now factor and calculate totient */
  nfactors = factor(n, &factors, &exponents);
  mpz_init(t);
  for (i = 0; i < nfactors; i++) {
    mpz_sub_ui(t, factors[i], 1);
    for (j = 1; j < exponents[i]; j++)
      mpz_mul(t, t, factors[i]);
    mpz_mul(tot, tot, t);
  }
  mpz_clear(t);
  clear_factors(nfactors, &factors, &exponents);
  mpz_clear(n);
}

void jordan_totient(mpz_t tot, const mpz_t n, unsigned long k)
{
  if (k == 0) {
    mpz_set_ui(tot, (mpz_cmp_ui(n, 1) == 0) ? 1 : 0);
  } else if (k == 1) {
    totient(tot, n);
  } else if (mpz_cmp_ui(n, 1) <= 0) {
    mpz_set_ui(tot, (mpz_cmp_ui(n, 1) == 0) ? 1 : 0);
  } else {
    mpz_t t;
    mpz_t* factors;
    int* exponents;
    int i, j, nfactors;

    nfactors = factor(n, &factors, &exponents);
    mpz_init(t);
    for (i = 0; i < nfactors; i++) {
      mpz_pow_ui(t, factors[i], k);
      mpz_sub_ui(factors[i], t, 1);
      for (j = 1; j < exponents[i]; j++)
        mpz_mul(factors[i], factors[i], t);
    }
    mpz_product(factors, 0, nfactors-1);
    mpz_set(tot, factors[0]);
    mpz_clear(t);
    clear_factors(nfactors, &factors, &exponents);
  }
}

void carmichael_lambda(mpz_t lambda, const mpz_t n)
{
  if (mpz_cmp_ui(n, 8) < 0) {
    totient(lambda, n);
  } else if (mpz_scan1(n, 0) == mpz_sizeinbase(n, 2)-1) {
    mpz_tdiv_q_2exp(lambda, n, 2);
  } else {
    mpz_t t;
    mpz_t* factors;
    int* exponents;
    int i, j, nfactors;

    nfactors = factor(n, &factors, &exponents);
    mpz_init(t);
    mpz_set_ui(lambda, 1);
    if (exponents[0] > 2 && mpz_cmp_ui(factors[0], 2) == 0) exponents[0]--;
    for (i = 0; i < nfactors; i++) {
      mpz_sub_ui(t, factors[i], 1);
      for (j = 1; j < exponents[i]; j++)
        mpz_mul(t, t, factors[i]);
      mpz_lcm(lambda, lambda, t);
    }
    mpz_clear(t);
    clear_factors(nfactors, &factors, &exponents);
  }
}

static const unsigned char _zntiny[7][7] = {
  {2},
  {0,2},
  {4,4,2},
  {0,0,0,2},
  {3,6,3,6,2},
  {0,2,0,2,0,2},
  {6,0,3,6,0,3,2} };
static const unsigned char _zn11[ 9] = {10,5,5,5,10,10,10,5,2};
static const unsigned char _zn13[11] = {12,3,6,4,12,12,4,3,6,12,2};
static const unsigned char _zn16[14] = {0,4,0,4,0,2,0,2,0,4,0,4,0,2};
static const unsigned char _zn17[15] = {8,16,4,16,16,16,8,8,16,16,16,4,16,8,2};

/* Do not alias any of the arguments */
static void _znorder1(mpz_t order, const mpz_t a, mpz_t p, int e, mpz_t t, mpz_t n)
{
  mpz_t phi, *factors;
  unsigned long s;
  int* exponents;
  int i, j, nfactors;

  mpz_pow_ui(n, p, e);

  /* Remove some simple cases */
  if (mpz_cmp_ui(n, 1) <= 0)  { mpz_set(order, n); return; }
  mpz_mod(t, a, n);
  if (mpz_cmp_ui(t, 1) <= 0)  { mpz_set(order, t); return; }
  /* These are purely for performance. */
  if (mpz_cmp_ui(n,17) <= 0) {
    unsigned long int nn = mpz_get_ui(n), aa = mpz_get_ui(t);
    if (nn <=  9)  { mpz_set_ui(order, _zntiny[nn-3][aa-2]); return; }
    if (nn == 11)  { mpz_set_ui(order, _zn11[aa-2]); return; }
    if (nn == 13)  { mpz_set_ui(order, _zn13[aa-2]); return; }
    if (nn == 16)  { mpz_set_ui(order, _zn16[aa-2]); return; }
    if (nn == 17)  { mpz_set_ui(order, _zn17[aa-2]); return; }
  }

  /* For 2^e, use the explicit structure of (Z/2^eZ)*. */
  if (mpz_cmp_ui(p, 2) == 0) {
    if (e <= 1) {
      mpz_set_ui(order, 1);
    } else if (e == 2) {
      mpz_set_ui(order, (mpz_fdiv_ui(t, 4) == 1) ? 1 : 2);
    } else if (mpz_fdiv_ui(t, 4) == 1) {
      mpz_sub_ui(t, t, 1);
      s = mpz_scan1(t, 0);
      mpz_set_ui(order, 1);
      if ((unsigned long)e > s) mpz_mul_2exp(order, order, (unsigned long)e - s);
    } else {
      mpz_add_ui(t, t, 1);
      s = mpz_scan1(t, 0);
      mpz_set_ui(order, 2);
      if ((unsigned long)e > s) mpz_mul_2exp(order, order, (unsigned long)e - s - 1);
    }
    return;
  }

  /* p is odd.  Compute and factor p-1. */
  mpz_init(phi);
  mpz_sub_ui(phi, p, 1);
  nfactors = factor(phi, &factors, &exponents);

  if (e == 1) {

    /* Abhijit Das, algorithm 1.7.  This is best for prime moduli. */
    mpz_set_ui(order, 1);
    for (i = 0; i < nfactors; i++) {
      mpz_divexact(t, phi, factors[i]);
      for (j = 1; j < exponents[i]; j++)
        mpz_divexact(t, t, factors[i]);
      mpz_powm(t, a, t, n);
      for (j = 0;  mpz_cmp_ui(t, 1) != 0;  mpz_powm(t, t, factors[i], n)) {
        if (j++ >= exponents[i]) {
          mpz_set_ui(order, 0);
          break;
        }
        mpz_mul(order, order, factors[i]);
      }
      if (j > exponents[i]) break;
    }

  } else {

    /* For odd p^e, compute ord_p(a), then lift using LTE:
     * ord_{p^e}(a) = ord_p(a) * p^max(0, e - v_p(a^ord_p(a)-1)).
     * This avoids factoring (p-1)*p^(e-1) and repeated powm modulo p^e. */
    mpz_mod(n, t, p);
    mpz_set(order, phi);
    for (i = 0; i < nfactors; i++) {
      for (j = 0; j < exponents[i]; j++) {
        mpz_divexact(phi, order, factors[i]);
        mpz_powm(t, n, phi, p);
        if (mpz_cmp_ui(t, 1) != 0)
          break;
        mpz_set(order, phi);
      }
    }
    s = 1;
    mpz_mul(n, p, p);
    while (s < (unsigned long)e) {
      mpz_powm(t, a, order, n);
      if (mpz_cmp_ui(t, 1) != 0)
        break;
      s++;
      if (s < (unsigned long)e)
        mpz_mul(n, n, p);
    }
    if ((unsigned long)e > s) {
      mpz_pow_ui(t, p, (unsigned long)e - s);
      mpz_mul(order, order, t);
    }

  }

  mpz_clear(phi);
  clear_factors(nfactors, &factors, &exponents);
}

int znorder(mpz_t res, const mpz_t ina, const mpz_t inn)
{
  mpz_t a, n, t, order;

  mpz_init(n);
  mpz_abs(n, inn);
  if (mpz_cmp_ui(n, 1) <= 0)
    { mpz_set(res, n); mpz_clear(n); return mpz_sgn(res); }

  mpz_init(a);
  mpz_mod(a, ina, n);
  if (mpz_cmp_ui(a, 1) <= 0)
    { mpz_set(res, a); mpz_clear(a); mpz_clear(n); return mpz_sgn(res); }

  mpz_init(t);
  mpz_gcd(t, a, n);
  if (mpz_cmp_ui(t, 1) != 0) {
    mpz_set_ui(res, 0);
    mpz_clear(t);
    mpz_clear(a);
    mpz_clear(n);
    return 0;
  }
  mpz_init_set_ui(order, 1);

  { /* Factor n, then lcm all the znorder(p^e). */
    mpz_t order1, t2, *factors;
    int* exponents;
    int i, nfactors;

    mpz_init(order1);
    mpz_init(t2);
    nfactors = factor(n, &factors, &exponents);
    for (i = 0; i < nfactors; i++) {
      _znorder1(order1, a, factors[i], exponents[i], t, t2);
      mpz_lcm(order, order, order1);
    }
    clear_factors(nfactors, &factors, &exponents);
    mpz_clear(t2);
    mpz_clear(order1);
  }
  mpz_set(res, order);
  mpz_clear(order);
  mpz_clear(t);
  mpz_clear(a);
  mpz_clear(n);
  return mpz_sgn(res);
}

static int _znprimroot_prime(mpz_t root, const mpz_t p,
                             const int only_odd_root, const int ispow)
{
  mpz_t t, a, phi, psquared, *factors;
  int i, nfactors, *exponents, found;

  /* It is assumed that p is an odd prime. */
  if (mpz_cmp_ui(p,3) == 0 && only_odd_root) { mpz_set_ui(root, 5); return 1; }

  mpz_init(t);
  mpz_init(a);
  mpz_init(phi);
  mpz_sub_ui(phi, p, 1);
  nfactors = factor(phi, &factors, &exponents);

  /* Replace each factor with phi/factor */
  for (i = 1; i < nfactors; i++)
    mpz_divexact(factors[i], phi, factors[i]);

  if (ispow) {
    mpz_init(psquared);
    mpz_pow_ui(psquared, p, 2);
  }

  for (mpz_set_ui(a,2);  mpz_cmp(a,p) < 0;  mpz_add_ui(a,a,1)) {
    if (only_odd_root && mpz_even_p(a)) continue;
    if (!mpz_cmp_ui(a,4) || !mpz_cmp_ui(a,8) || !mpz_cmp_ui(a,9)) continue;
    if (mpz_kronecker(a,p) != -1) continue;
    for (i = 1; i < nfactors; i++) {
      mpz_powm(t, a, factors[i], p);
      if (mpz_cmp_ui(t, 1) == 0)
        break;
    }
    if (i == nfactors) {
      if (ispow) { /* For p^k, k>1, g is not a root if g^(p-1) = 1 mod p^2. */
        mpz_powm(t, a, phi, psquared);
        if (mpz_cmp_ui(t, 1) == 0)
          continue;
      }
      break;
    }
  }
  found = (mpz_cmp(a,p) < 0);  /* We should always find one */
  if (found) mpz_set(root,a);
  clear_factors(nfactors, &factors, &exponents);
  if (ispow) mpz_clear(psquared);
  mpz_clear(phi);
  mpz_clear(a);
  mpz_clear(t);
  return found;
}

/* root is set to -1 if no root.  n should be positive. */
void znprimroot(mpz_t root, const mpz_t n)
{
  mpz_t pk, p;
  unsigned long k;
  int is_neven, kr;

  if (mpz_cmp_ui(n, 4) <= 0) {
    if (mpz_sgn(n) > 0) mpz_sub_ui(root, n, 1);
    else                mpz_set_si(root, -1);
    return;
  }
  if (mpz_divisible_ui_p(n, 4)) {
    mpz_set_si(root, -1);
    return;
  }

  mpz_init(pk);
  mpz_init(p);

  is_neven = mpz_even_p(n);
  if (is_neven)  mpz_tdiv_q_2exp(pk, n, 1);
  else           mpz_set(pk, n);

  k = prime_power(p, pk);
  mpz_clear(pk);
  if (k == 0) { mpz_set_si(root,-1);  mpz_clear(p);  return; }

  /* n is either p^k or 2p^k for p an odd prime.  A root exists. */
  kr = _znprimroot_prime(root, p, is_neven, k>1);
  if (kr == 0) {
    gmp_printf("  Failed to find primitive root for n %Zd\n", n);
    mpz_set_si(root, -1);
  }
  mpz_clear(p);
}

#if 0  /* Replaced by rtau(res,n) in gmp_main.c. */
static const int32_t tau_table[] = {
  0,1,-24,252,-1472,4830,-6048,-16744,84480,-113643,-115920,534612,-370944,-577738,401856,1217160,987136,-6905934,2727432,10661420,-7109760,-4219488,-12830688,18643272,21288960,-25499225,13865712,-73279080,24647168,128406630,-29211840,-52843168,-196706304,134722224,165742416,-80873520,167282496,-182213314,-255874080,-145589976,408038400,308120442,101267712,-17125708,-786948864,-548895690,-447438528
};
#define NTAU (sizeof(tau_table)/sizeof(tau_table[0]))
void ramanujan_tau(mpz_t res, const mpz_t n)
{
  mpz_t t, t1, t2, t3, t4, *factors;
  int i, nfactors, *exponents;
  unsigned long j, p2;

  if (mpz_cmp_ui(n, NTAU) < 0) {
    if (mpz_sgn(n) <= 0) mpz_set_si(res, 0);
    else                 mpz_set_si(res, tau_table[mpz_get_ui(n)]);
    return;
  }

  /* We are doing far too much work here for sigma5.  We could do it just
   * for primes then use the multiplicative property.  However that works
   * for prime *powers*, so it isn't quite so simple.  This solution also
   * gets to be high memory use. */

  /* Pari/GP does this using Hurwitz class numbers.  That is a more
   * complicated but far more efficient solution. */

  mpz_init(t); mpz_init(t1); mpz_init(t2); mpz_init(t3); mpz_init(t4);
  nfactors = factor(n, &factors, &exponents);
  for (i = 0; i < nfactors; i++) {
    /* t = tau(p) */
    if (mpz_cmp_ui(factors[i], NTAU) < 0) {
      mpz_set_si(t, tau_table[mpz_get_ui(factors[i])]);
    } else {
      mpz_pow_ui(t, factors[i], 11);   mpz_add_ui(t, t, 1); /* sigma(t,f,11) */
      mpz_mul_ui(t1, t, 65);
      mpz_pow_ui(t, factors[i],  5);   mpz_add_ui(t, t, 1); /* sigma(t,f, 5) */
      mpz_mul_ui(t2, t, 691);
      mpz_add(t1, t1, t2);

      /* t1 in use. t2 accumulate.  t3, t4, t free. */
      mpz_sub_ui(t, factors[i], 1);
      mpz_tdiv_q_2exp(t, t, 1);
      p2 = mpz_get_ui(t);
      mpz_set_ui(t2, 0);
      for (j = 1; j <= p2; j++) {
        mpz_set_uv(t, j);
        sigma(t3, t, 5);
        mpz_sub_ui(t, factors[i], j);
        sigma(t, t, 5);
        mpz_mul(t4, t3, t);
        mpz_add(t2, t2, t4);
      }
      mpz_mul_ui(t2, t2, 2*691*252);
      mpz_sub(t, t1, t2);
      mpz_tdiv_q_ui(t, t, 756);
    }
    /* t holds tau(p), all other temps are free. */

    if (exponents[i] > 1) {
      mpz_pow_ui(t1, t, exponents[i]);
      if (exponents[i] == 2) {
        mpz_pow_ui(t2, factors[i], 11);
      } else if (exponents[i] == 3) {
        mpz_pow_ui(t2, factors[i], 11);
        mpz_mul(t2, t2, t);
        mpz_mul_ui(t2, t2, 2);
      } else {
        /* t1 = t^e  t2 = neg sum,  t3 = prod,  t4 = temp */
        mpz_set_ui(t2, 0);
        for (j = 1; j <= ((unsigned long)exponents[i])>>1; j++) {
          mpz_set_si(t3, (j&1) ? -1 : 1);
          mpz_pow_ui(t4, factors[i], 11*j);
          mpz_mul(t3, t3, t4);
          mpz_bin_uiui(t4, exponents[i]-j, exponents[i]-2*j);
          mpz_mul(t3, t3, t4);
          mpz_pow_ui(t4, t, exponents[i]-2*j);
          mpz_mul(t3, t3, t4);
          mpz_sub(t2, t2, t3);
        }
      }
      mpz_sub(t, t1, t2);
    }
    mpz_set(factors[i], t);
  }
  mpz_product(factors, 0, nfactors-1);
  mpz_set(res, factors[0]);
  clear_factors(nfactors, &factors, &exponents);
  mpz_clear(t4); mpz_clear(t3); mpz_clear(t2); mpz_clear(t1); mpz_clear(t);
}
#endif

int is_semiprime(const mpz_t n)
{
  int ret;
  unsigned long k, div, lim = 6000;
  mpz_t t;

  if (mpz_cmp_ui(n,6) < 0)
    return (mpz_cmp_ui(n,4) == 0);

  mpz_init(t);

  div = _GMP_trial_factor(n, 2, lim);
  if (div > 0) {
    mpz_divexact_ui(t, n, div);
    ret = _GMP_is_prime(t);
    mpz_clear(t);
    return !!ret;
  }
  /* No small divisors */
  if (_GMP_BPSW(n))
    { mpz_clear(t); return 0; }
  if (mpz_ui_pow_ui(t,lim,3), mpz_cmp(n,t) < 0)
    { mpz_clear(t); return 1; }

  /* Quick check for power */
  k = power_factor(n, t);
  if (k >= 2) {
    ret = (k == 2 && _GMP_is_prime(t));
    mpz_clear(t);
    return ret;
  }

  /* Number is composite, isn't tiny, and has no small divisors */
  if (    0
       || _GMP_pbrent_factor(n, t, 1, 15000)
       || _GMP_pminus1_factor(n, t, 50000, 500000)
       || _GMP_ECM_FACTOR(n, t,     800, 10)
       || _GMP_ECM_FACTOR(n, t,    8000, 20)
       || _GMP_ECM_FACTOR(n, t,   80000, 40)
       || _GMP_ECM_FACTOR(n, t,  320000, 40)
       || _GMP_ECM_FACTOR(n, t, 1000000, 80)
     ) {
    ret = _GMP_BPSW(t);
    if (ret) {
      mpz_divexact(t, n, t);
      ret = _GMP_BPSW(t);
    }
    mpz_clear(t);
    return !!ret;
  }
  /* No luck finding a small factor.  Do it the hard way. */
  {
    mpz_t* factors;
    int* exponents;
    int nfactors, i, j;

    nfactors = factor(n, &factors, &exponents);
    for (i = 0, j = 0; i < nfactors; i++)
      j += exponents[i];
    clear_factors(nfactors, &factors, &exponents);
    mpz_clear(t);
    return (j == 2);
  }
}


/*****************************************************************************/

UV _GMP_trial_factor(const mpz_t n, UV from_n, UV to_n)
{
  size_t log2n = mpz_sizeinbase(n, 2);
  UV p = 0;
  PRIME_ITERATOR(iter);

  if (mpz_cmp_ui(n, 6) < 0) {
    unsigned long un = mpz_get_ui(n);
    if (un == 1) p = 1;
    else if (un == 4 && from_n <= 2 && to_n >= 2) p = 2;
    prime_iterator_destroy(&iter);
    return p;
  }
  if      (from_n <= 2 && to_n >= 2 && mpz_even_p(n))             p = 2;
  else if (from_n <= 3 && to_n >= 3 && mpz_divisible_ui_p(n, 3))  p = 3;
  else if (from_n <= 5 && to_n >= 5 && mpz_divisible_ui_p(n, 5))  p = 5;
  if (p != 0) {
    prime_iterator_destroy(&iter);
    return p;
  }

  if (from_n < 7)
    from_n = 7;
  if (from_n > to_n)
    { prime_iterator_destroy(&iter);  return 0; }
  /* p will be the next prime >= from_n */
  prime_iterator_setprime(&iter, from_n-1);
  p = prime_iterator_next(&iter);

  /* All native math if n fits in an unsigned long */
  if (log2n <= sizeof(unsigned long)*8) {
    unsigned long un = mpz_get_ui(n);
    unsigned long sqrtn = (unsigned long) sqrt((double)un);
    /* Be extra careful here, as we are using unsigned long, which may not
     * match a UV.  But GMP's ui is 'unsigned long' so that's what we have
     * to deal with.  We want to make sure we get the correct integer sqrt,
     * but also watch out for overflow. */
    while (sqrtn*sqrtn > un) sqrtn--;
    while ( (sqrtn+1)*(sqrtn+1) <= un
            && sqrtn < (1UL << 4*sizeof(unsigned long)) )
      sqrtn++;
    if (to_n > sqrtn)
      to_n = sqrtn;
    while (p <= to_n) {
      if ((un % p) == 0)
        break;
      p = prime_iterator_next(&iter);
    }
    prime_iterator_destroy(&iter);
    return (p <= to_n) ? p : 0;
  }

  /* For "small" numbers, this simple method is best. */
  {
    UV small_to = (log2n < 3000)  ?  to_n  :  30000;
    while (p <= small_to) {
      if (mpz_divisible_ui_p(n, p))
        break;
      p = prime_iterator_next(&iter);
    }
    if (p <= small_to || p > to_n) {
      prime_iterator_destroy(&iter);
      return (p <= small_to) ? p : 0;
    }
  }

  /* Simple treesieve.
   * This is much faster than simple divisibility for really big numbers.
   * Credit to Jens K Andersen for writing up the generic algorithm.
   *
   * This will search until the first group element is > to_n, which means
   * we will search a bit farther than to_n.
   */
  {
    UV found = 0;
    unsigned long* xn;        /* leaves */
    mpz_t* xtree[16+1];       /* the tree (maxdepth = 16) */
    mpz_t* xtemp;
    unsigned int i, j, d, depth, leafsize, nleaves;

    /* Decide on the tree depth (3-16) and number of leaves (10-31) */
    {
      unsigned int dp = log2n >> 10;
      depth = 0;
      while (dp >>= 1) depth++;
      if (depth < 3) depth = 3;
      if (depth > 16) depth = 16;
    }
    leafsize = log2n / (1U << depth) / 68;
    nleaves = 1 << depth;
    /* printf("log2n %lu  depth %u  leafsize %u  nleaves %u\n",log2n,depth,leafsize,nleaves); */

    New(0, xn, nleaves * leafsize, unsigned long);
    for (d = 0; d <= depth; d++) {
      unsigned int nodes = 1 << (depth - d);
      New(0, xtree[d], nodes, mpz_t);
      for (j = 0; j < nodes; j++)
        mpz_init(xtree[d][j]);
    }
    xtemp = xtree[1];   /* implies mindepth = 3 */

    while (!found && p <= to_n) {
      /* Create nleaves x[0] values, each the product of leafsize primes */
      for (i = 0; i < nleaves; i++) {
        for (j = 0; j < 4; j++)                  /* Create 4 sub-products */
          mpz_set_ui(xtemp[j], 1);
        for (j = 0; j < leafsize; j++) {
          xn[i*leafsize+j] = p;
          mpz_mul_ui(xtemp[j&3], xtemp[j&3], p);
          p = prime_iterator_next(&iter);
        }
        mpz_mul(xtemp[0], xtemp[0], xtemp[1]);   /* Combine for final product*/
        mpz_mul(xtemp[2], xtemp[2], xtemp[3]);
        mpz_mul(xtree[0][i], xtemp[0], xtemp[2]);
      }
      /* Multiply product tree, xtree[depth][0] has nleaves*leafsize product */
      for (d = 1; d <= depth; d++)
        for (i = 0; i < (1U << (depth-d)); i++)
          mpz_mul(xtree[d][i], xtree[d-1][2*i], xtree[d-1][2*i+1]);
      /* Go backwards replacing the products with remainders */
      mpz_tdiv_r(xtree[depth][0], n, xtree[depth][0]);
      for (d = 1; d <= depth; d++)
        for (i = 0; i < (1U << d); i++)
          mpz_tdiv_r(xtree[depth-d][i], xtree[depth-d+1][i>>1], xtree[depth-d][i]);
      /* Search each leaf for divisors */
      for (i = 0; !found && i < nleaves; i++)
        for (j = 0; j < leafsize; j++)
          if (mpz_divisible_ui_p(xtree[0][i], xn[i*leafsize+j]))
            { found = xn[i*leafsize+j]; break; }
    }
    p = found;
    for (d = 0; d <= depth; d++) {
      unsigned int nodes = 1U << (depth - d);
      for (j = 0; j < nodes; j++)
        mpz_clear(xtree[d][j]);
      Safefree(xtree[d]);
    }
    Safefree(xn);
    if (p > 0 && !mpz_divisible_ui_p(n, p))
      croak("incorrect trial factor\n");
  }
  prime_iterator_destroy(&iter);
  return p;
}


#define TEST_FOR_2357(n, f) \
  { \
    if (mpz_divisible_ui_p(n, 2)) { mpz_set_ui(f, 2); return 1; } \
    if (mpz_divisible_ui_p(n, 3)) { mpz_set_ui(f, 3); return 1; } \
    if (mpz_divisible_ui_p(n, 5)) { mpz_set_ui(f, 5); return 1; } \
    if (mpz_divisible_ui_p(n, 7)) { mpz_set_ui(f, 7); return 1; } \
    if (mpz_cmp_ui(n, 121) < 0) { return 0; } \
  }

int _GMP_prho_factor(const mpz_t n, mpz_t f, UV a, UV rounds)
{
  mpz_t U, V, oldU, oldV, m;
  int i;
  const UV inner = 256;

  TEST_FOR_2357(n, f);
  rounds = (rounds + inner - 1) / inner;
  mpz_init_set_ui(U, 7);
  mpz_init_set_ui(V, 7);
  mpz_init(m);
  mpz_init(oldU);
  mpz_init(oldV);
  while (rounds-- > 0) {
    mpz_set_ui(m, 1); mpz_set(oldU, U);  mpz_set(oldV, V);
    for (i = 0; i < (int)inner; i++) {
      mpz_mul(U, U, U);  mpz_add_ui(U, U, a);  mpz_tdiv_r(U, U, n);
      mpz_mul(V, V, V);  mpz_add_ui(V, V, a);  mpz_tdiv_r(V, V, n);
      mpz_mul(V, V, V);  mpz_add_ui(V, V, a);  mpz_tdiv_r(V, V, n);
      mpz_sub(f, U, V);
      mpz_mul(m, m, f);
      if ((i%4) == ((inner-1)%4)) mpz_tdiv_r(m, m, n);
    }
    mpz_gcd(f, m, n);
    if (!mpz_cmp_ui(f, 1))
      continue;
    if (!mpz_cmp(f, n)) {
      /* f == n, so we have to back up to see what factor got found */
      mpz_set(U, oldU); mpz_set(V, oldV);
      i = inner;
      do {
        mpz_mul(U, U, U);  mpz_add_ui(U, U, a);  mpz_tdiv_r(U, U, n);
        mpz_mul(V, V, V);  mpz_add_ui(V, V, a);  mpz_tdiv_r(V, V, n);
        mpz_mul(V, V, V);  mpz_add_ui(V, V, a);  mpz_tdiv_r(V, V, n);
        if (mpz_cmp(U, V) >= 0)  mpz_sub(f, U, V);
        else                     mpz_sub(f, V, U);
        mpz_gcd(f, f, n);
      } while (!mpz_cmp_ui(f, 1) && i-- != 0);
      if ( (!mpz_cmp_ui(f, 1)) || (!mpz_cmp(f, n)) )  break;
    }
    mpz_clear(U); mpz_clear(V); mpz_clear(m); mpz_clear(oldU); mpz_clear(oldV);
    return 1;
  }
  mpz_clear(U); mpz_clear(V); mpz_clear(m); mpz_clear(oldU); mpz_clear(oldV);
  mpz_set(f, n);
  return 0;
}

int _GMP_pbrent_factor(const mpz_t n, mpz_t f, UV a, UV rounds)
{
  mpz_t Xi, Xm, saveXi, m, t;
  UV i, r;
  const UV inner = 256;

  TEST_FOR_2357(n, f);
  mpz_init_set_ui(Xi, 2);
  mpz_init_set_ui(Xm, 2);
  mpz_init(m);
  mpz_init(t);
  mpz_init(saveXi);

  r = 1;
  mpz_set_ui(f,1);
  while (rounds > 0) {
    UV rleft = (r > rounds) ? rounds : r;
    while (rleft > 0) {   /* Do rleft rounds, inner at a time */
      UV dorounds = (rleft > inner) ? inner : rleft;
      mpz_set_ui(m, 1);
      mpz_set(saveXi, Xi);
      for (i = 0; i < dorounds; i++) {
        mpz_mul(t, Xi, Xi);  mpz_add_ui(t, t, a);  mpz_tdiv_r(Xi, t, n);
        mpz_sub(f, Xm, Xi);
        mpz_mul(m, m, f);
        if ((i%4) == ((dorounds-1)%4)) mpz_tdiv_r(m, m, n);
      }
      rleft -= dorounds;
      rounds -= dorounds;
      mpz_gcd(f, m, n);
      if (mpz_cmp_ui(f, 1) != 0)
        break;
    }
    if (!mpz_cmp_ui(f, 1)) {
      r *= 2;
      mpz_set(Xm, Xi);
      continue;
    }
    if (!mpz_cmp(f, n)) {
      /* f == n, so we have to back up to see what factor got found */
      mpz_set(Xi, saveXi);
      do {
        mpz_mul(t, Xi, Xi);  mpz_add_ui(t, t, a);  mpz_tdiv_r(Xi, t, n);
        mpz_sub(f, Xm, Xi); if (mpz_sgn(f) < 0) mpz_add(f,f,n);
        mpz_gcd(f, f, n);
      } while (!mpz_cmp_ui(f, 1) && r-- != 0);
    }
    break;
  }
  mpz_clear(Xi); mpz_clear(Xm); mpz_clear(m); mpz_clear(saveXi); mpz_clear(t);
  if (!mpz_cmp_ui(f, 1) || !mpz_cmp(f, n)) {
    mpz_set(f, n);
    return 0;
  }
  return 1;
}

/* References for P-1:
 *  Montgomery 1987:  https://cr.yp.to/bib/1987/montgomery.pdf
 *  Brent 1990:       http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.127.4316
 *
 * The main advantage of this over ECM is that it is *much* lower overhead,
 * so very cheap to run with relatively small B1,B2 values.  A disadvantage
 * is no continuation method, so subsequent calls with larger B1,B2 will
 * repeat all the previous work.  ECM is much better for harder factorisations,
 * so we typically want to try a little p-1 then move to ECM (or QS).
 */
int _GMP_pminus1_factor(const mpz_t n, mpz_t f, UV B1, UV B2)
{
  mpz_t a, savea, t, uvtmp;
  UV q, saveq, j, sqrtB1;
  int _verbose = get_verbose_level();
  PRIME_ITERATOR(iter);

  TEST_FOR_2357(n, f);
  if (B1 < 7) return 0;

  mpz_init(a);
  mpz_init(savea);
  mpz_init(t);
  /* GMP's _ui arguments can be narrower than Perl's UV. */
  if (sizeof(UV) > sizeof(unsigned long)) mpz_init(uvtmp);

  if (_verbose>2) gmp_printf("# p-1 trying %Zd (B1=%"UVuf" B2=%"UVuf")\n", n, B1, B2);

  /* STAGE 1
   * See Montgomery 1987 p249-250 or Brent 1990 p5.  We can take E to be the
   * lcm of integers to B1, then gcd(a^E-1,n) may be a factor of n.  While
   * we could actually calculate the LCM, it is quite inefficient to do so.
   * There are various ways to speed it up, but generally we prefer to do it
   * the way Brent indicates, which is one powmod for each prime p below B1,
   * a = a^(p^e) mod n, where e is the largest e such that p^e <= B1.
   *
   * One thing that can speed things up quite a bit is not running the GCD
   * on every step.  However with small factors this means we can easily end
   * up with multiple factors between GCDs, so we allow backtracking.  This
   * could also be added to stage 2, but it's far less likely to happen there.
   */
  j = 15;
  mpz_set_ui(a, 2);
  mpz_set_ui(savea, 2);
  saveq = 2;
  /* We could wrap this in a loop trying a few different a values, in case
   * the current one ended up going to 0. */
  q = 2;
  mpz_set_ui(t, 1);
  sqrtB1 = (UV) sqrt(B1);
  while (q <= B1) {
    UV k = q;
    if (q <= sqrtB1) {
      UV kmin = B1/q;
      while (k <= kmin)
        k *= q;
    }
    /* Accumulate powers for a, preserving the full UV value. */
    if (k <= ULONG_MAX) {
      mpz_mul_ui(t, t, k);
    } else {
      mpz_set_uv(uvtmp, k);
      mpz_mul(t, t, uvtmp);
    }
    if ( (j++ % 32) == 0) {
      mpz_powm(a, a, t, n);     /* a=a^(k1*k2*k3*...) mod n */
      if (mpz_sgn(a))  mpz_sub_ui(t, a, 1);
      else             mpz_sub_ui(t, n, 1);
      mpz_gcd(f, t, n);         /* f = gcd(a-1, n) */
      mpz_set_ui(t, 1);
      if (mpz_cmp(f, n) == 0)
        break;
      if (mpz_cmp_ui(f, 1) != 0)
        goto end_success;
      saveq = q;
      mpz_set(savea, a);
    }
    q = prime_iterator_next(&iter);
  }
  mpz_powm(a, a, t, n);
  if (mpz_sgn(a))  mpz_sub_ui(t, a, 1);
  else             mpz_sub_ui(t, n, 1);
  mpz_gcd(f, t, n);
  if (mpz_cmp(f, n) == 0) {
    /* We found multiple factors.  Loop one at a time. */
    prime_iterator_setprime(&iter, saveq);
    mpz_set(a, savea);
    for (q = saveq; q <= B1; q = prime_iterator_next(&iter)) {
      UV k = q;
      if (q <= sqrtB1) {
        UV kmin = B1/q;
        while (k <= kmin)
          k *= q;
      }
      if (k <= ULONG_MAX) {
        mpz_powm_ui(a, a, k, n);
      } else {
        mpz_set_uv(uvtmp, k);
        mpz_powm(a, a, uvtmp, n);
      }
      mpz_sub_ui(t, a, 1);
      mpz_gcd(f, t, n);
      if (mpz_cmp(f, n) == 0)
        goto end_fail;
      if (mpz_cmp_ui(f, 1) != 0)
        goto end_success;
    }
  }
  if ( (mpz_cmp_ui(f, 1) != 0) && (mpz_cmp(f, n) != 0) )
    goto end_success;

  /* STAGE 2
   * This is the standard continuation which replaces the powmods in stage 1
   * with two mulmods, with a GCD every 128 primes (no backtracking).
   * This is quite a bit faster than stage 1.
   * See Montgomery 1987, p250-253 for possible optimizations.
   * We quickly precalculate a few of the prime gaps, and lazily cache others
   * up to a gap of 222.  That's enough for a B2 value of 189 million.  We
   * still work above that, we just won't cache the value for big gaps.
   */
  if (B2 > B1 && q <= B2) {
    mpz_t b, bm, bmdiff;
    mpz_t precomp_bm[111];
    int   is_precomp[111] = {0};
    UV* primes = 0;
    UV sp = 1;

    mpz_init(bmdiff);
    mpz_init_set(bm, a);
    mpz_init_set_ui(b, 1);

    /* Set the first 22 differences */
    mpz_powm_ui(bmdiff, bm, 2, n);
    mpz_init_set(precomp_bm[0], bmdiff);
    is_precomp[0] = 1;
    for (j = 1; j < 22; j++) {
      mpz_mul(bmdiff, bmdiff, bm);
      mpz_mul(bmdiff, bmdiff, bm);
      mpz_tdiv_r(bmdiff, bmdiff, n);
      mpz_init_set(precomp_bm[j], bmdiff);
      is_precomp[j] = 1;
    }

    if (q <= ULONG_MAX) {
      mpz_powm_ui(a, a, q, n);
    } else {
      mpz_set_uv(uvtmp, q);
      mpz_powm(a, a, uvtmp, n);
    }
    if (B2 < 10000000) {
      /* grab all the primes at once.  Hack around non-perfect iterator. */
      primes = sieve_to_n(B2+300, 0);
      for (sp = B1>>4; primes[sp] <= q; sp++)  ;
      /* q is the first prime above B1; primes[sp] is its successor. */
    }

    j = 95;  /* First GCD after 34 primes, then every 128. */
    while (q <= B2) {
      UV lastq, qdiff;

      /* Include the first prime above B1 and any prime equal to B2. */
      if (mpz_sgn(a))  mpz_sub_ui(t, a, 1);
      else             mpz_sub_ui(t, n, 1);
      mpz_mul(b, b, t);
      if ((j % 2) == 0)           /* put off mods a little */
        mpz_tdiv_r(b, b, n);
      if ( (j++ % 128) == 0) {    /* GCD every so often */
        mpz_gcd(f, b, n);
        /* With no backtracking, a whole-input GCD cannot recover. */
        if (mpz_cmp_ui(f, 1) != 0)
          break;
      }

      lastq = q;
      q = primes ? primes[sp++] : prime_iterator_next(&iter);
      if (q > B2) break;
      qdiff = (q - lastq) / 2 - 1;

      if (qdiff < 111 && is_precomp[qdiff]) {
        mpz_mul(t, a, precomp_bm[qdiff]);
      } else if (qdiff < 111) {
        mpz_powm_ui(bmdiff, bm, q-lastq, n);
        mpz_init_set(precomp_bm[qdiff], bmdiff);
        is_precomp[qdiff] = 1;
        mpz_mul(t, a, bmdiff);
      } else {
        mpz_powm_ui(bmdiff, bm, q-lastq, n);  /* Big gap */
        mpz_mul(t, a, bmdiff);
      }
      mpz_tdiv_r(a, t, n);
    }
    mpz_gcd(f, b, n);
    mpz_clear(b);
    mpz_clear(bm);
    mpz_clear(bmdiff);
    for (j = 0; j < 111; j++) {
      if (is_precomp[j])
        mpz_clear(precomp_bm[j]);
    }
    if (primes != 0) Safefree(primes);
    if ( (mpz_cmp_ui(f, 1) != 0) && (mpz_cmp(f, n) != 0) )
      goto end_success;
  }

  end_fail:
    mpz_set(f,n);
  end_success:
    prime_iterator_destroy(&iter);
    mpz_clear(a);
    mpz_clear(savea);
    mpz_clear(t);
    if (sizeof(UV) > sizeof(unsigned long)) mpz_clear(uvtmp);
    if ( (mpz_cmp_ui(f, 1) != 0) && (mpz_cmp(f, n) != 0) ) {
      if (_verbose>2) gmp_printf("# p-1: %Zd\n", f);
      return 1;
    }
    if (_verbose>2) gmp_printf("# p-1: no factor\n");
    mpz_set(f, n);
    return 0;
}

static void pp1_pow(mpz_t X, mpz_t Y, unsigned long exp, const mpz_t n)
{
  mpz_t x0;
  unsigned long bit;
  {
    unsigned long v = exp;
    unsigned long b = 1;
    while (v >>= 1) b++;
    bit = 1UL << (b-2);
  }
  mpz_init_set(x0, X);
  mpz_mul(Y, X, X);
  mpz_sub_ui(Y, Y, 2);
  mpz_tdiv_r(Y, Y, n);
  while (bit) {
    if ( exp & bit ) {
      mpz_mul(X, X, Y);
      mpz_sub(X, X, x0);
      mpz_mul(Y, Y, Y);
      mpz_sub_ui(Y, Y, 2);
    } else {
      mpz_mul(Y, X, Y);
      mpz_sub(Y, Y, x0);
      mpz_mul(X, X, X);
      mpz_sub_ui(X, X, 2);
    }
    mpz_mod(X, X, n);
    mpz_mod(Y, Y, n);
    bit >>= 1;
  }
  mpz_clear(x0);
}

int _GMP_pplus1_factor(const mpz_t n, mpz_t f, UV P0, UV B1, UV B2)
{
  UV j, q, saveq, sqrtB1;
  mpz_t X, Y, saveX;
  PRIME_ITERATOR(iter);

  TEST_FOR_2357(n, f);
  if (B1 < 7) return 0;

  mpz_init_set_ui(X, P0);
  mpz_init(Y);
  mpz_init(saveX);

  /* Montgomery 1987 */
  if (P0 == 0) {
    mpz_set_ui(X, 7);
    if (mpz_invert(X, X, n)) {
      mpz_mul_ui(X, X, 2);
      mpz_mod(X, X, n);
    } else
      P0 = 1;
  }
  if (P0 == 1) {
    mpz_set_ui(X, 5);
    if (mpz_invert(X, X, n)) {
      mpz_mul_ui(X, X, 6);
      mpz_mod(X, X, n);
    } else
      P0 = 2;
  }
  if (P0 == 2) {
    mpz_set_ui(X, 11);
    if (mpz_invert(X, X, n)) {
      mpz_mul_ui(X, X, 23);
      mpz_mod(X, X, n);
    }
  }

  sqrtB1 = (UV) sqrt(B1);
  j = 8;
  q = 2;
  saveq = q;
  mpz_set(saveX, X);
  while (q <= B1) {
    UV k = q;
    if (q <= sqrtB1) {
      UV kmin = B1/q;
      while (k <= kmin)
        k *= q;
    }
    pp1_pow(X, Y, k, n);
    if ( (j++ % 16) == 0) {
      mpz_sub_ui(f, X, 2);
      if (mpz_sgn(f) == 0)        break;
      mpz_gcd(f, f, n);
      if (mpz_cmp(f, n) == 0)     break;
      if (mpz_cmp_ui(f, 1) > 0)   goto end_success;
      saveq = q;
      mpz_set(saveX, X);
    }
    q = prime_iterator_next(&iter);
  }
  mpz_sub_ui(f, X, 2);
  mpz_gcd(f, f, n);
  if (mpz_cmp_ui(X, 2) == 0 || mpz_cmp(f, n) == 0) {
    /* Backtrack */
    prime_iterator_setprime(&iter, saveq);
    mpz_set(X, saveX);
    for (q = saveq; q <= B1; q = prime_iterator_next(&iter)) {
      UV k = q;
      if (q <= sqrtB1) {
        UV kmin = B1/q;
        while (k <= kmin)
          k *= q;
      }
      pp1_pow(X, Y, k, n);
      mpz_sub_ui(f, X, 2);
      if (mpz_sgn(f) == 0)        goto end_fail;
      mpz_gcd(f, f, n);
      if (mpz_cmp(f, n) == 0)     break;
      if (mpz_cmp_ui(f, 1) > 0)   goto end_success;
    }
  }
  if ( (mpz_cmp_ui(f, 1) > 0) && (mpz_cmp(f, n) != 0) )
    goto end_success;
  /* TODO: stage 2 */
  end_fail:
    mpz_set(f,n);
  end_success:
    prime_iterator_destroy(&iter);
    mpz_clear(X);  mpz_clear(Y);  mpz_clear(saveX);
    return (mpz_cmp_ui(f, 1) != 0) && (mpz_cmp(f, n) != 0);
}

int _GMP_cheb_factor(const mpz_t n, mpz_t f, UV B, UV initx)
{
  unsigned long p;
  double logB;
  mpz_t x, inv, t, k, P, Q;
  PRIME_ITERATOR(iter);

  if (B == 0) { B = mpz_sizeinbase(n,2);  B = B*B*B; }

  TEST_FOR_2357(n, f);
  if (B < 7) return 0;

  logB = log(B);
  mpz_init_set_ui(inv, 2);
  mpz_invert(inv, inv, n);   /* multiplying by this will divide by two */
  mpz_init_set_ui(x, (initx == 0) ? 72 : initx);
  mpz_init(t);
  mpz_init(k);
  mpz_init(P);
  mpz_init_set_ui(Q, 1);

  mpz_set_ui(f, 1);
  for (p = 2; p <= B && mpz_cmp_ui(f,1) <= 0; p = prime_iterator_next(&iter)) {
    unsigned long lgbp = (unsigned long) (logB / log(p));   /* Alternately logint(mpzB,p) */
    if (lgbp > 1)  mpz_ui_pow_ui(k, p, lgbp);
    else           mpz_set_uv(k, p);
    mpz_mul_2exp(P, x, 1);
    lucasvmod(x, P, Q, k, n, t);
    mpz_mul(x, x, inv);
    mpz_mod(x, x, n);
    mpz_sub_ui(t, x, 1);
    mpz_gcd(f, t, n);
  }
  if (mpz_cmp_ui(f,1) <= 0)
    mpz_set(f,n);
  prime_iterator_destroy(&iter);
  mpz_clear(Q);  mpz_clear(P);  mpz_clear(k);  mpz_clear(t);
  mpz_clear(x);  mpz_clear(inv);
  return (mpz_cmp_ui(f, 1) > 0) && (mpz_cmp(f, n) < 0);
}

int _GMP_holf_factor(const mpz_t n, mpz_t f, UV rounds)
{
  mpz_t npre, s, m;
  UV i;

#define PREMULT 480   /* 1  2  6  12  480  151200 */

  TEST_FOR_2357(n, f);
  if (mpz_perfect_square_p(n)) {
    mpz_sqrt(f, n);
    return 1;
  }

  mpz_init(s);
  mpz_init(m);
  mpz_init(npre);
  mpz_mul_ui(npre, n, PREMULT);

  for (i = 1; i <= rounds; i++) {
    mpz_mul_ui(f, npre, i);    /* f = n*i */
    if (mpz_perfect_square_p(f)) {
      /* s^2 = n*i, so m = s^2 mod n = 0.  Hence f = GCD(n, s) = GCD(n, n*i) */
      mpz_gcd(f, f, n);
      mpz_clear(npre); mpz_clear(s); mpz_clear(m);
      if (mpz_cmp(f, n) == 0)  return 0;
      return 1;
    }
    mpz_sqrt(s, f);
    mpz_add_ui(s, s, 1);    /* s = ceil(sqrt(n*i)) */
    mpz_mul(m, s, s);
    mpz_sub(m, m, f);       /* m = s^2 mod n = s^2 - n*i */
    if (mpz_perfect_square_p(m)) {
      mpz_sqrt(f, m);
      mpz_sub(s, s, f);
      mpz_gcd(f, s, n);
      mpz_clear(npre); mpz_clear(s); mpz_clear(m);
      return (mpz_cmp_ui(f, 1) > 0);
    }
  }
  mpz_set(f, n);
  mpz_clear(npre); mpz_clear(s); mpz_clear(m);
  return 0;
}

static INLINE uint64_t mpz_get64(const mpz_t n) {
  uint64_t v = (uint64_t)mpz_getlimbn(n, 0);
#if GMP_NUMB_BITS < 64
  unsigned shift;
  mp_size_t i;
  for (i = 1, shift = GMP_NUMB_BITS; shift < 64; i++, shift += GMP_NUMB_BITS)
    v |= (uint64_t)mpz_getlimbn(n, i) << shift;
#endif
  return v;
}
static INLINE void mpz_set64(mpz_t n, uint64_t v) {
#if ULONG_MAX >= UINT64_MAX
  mpz_set_ui(n, (unsigned long)v);
#else
  if (v <= ULONG_MAX) {
    mpz_set_ui(n, (unsigned long)v);
  } else {
    mpz_set_ui(n, (unsigned long)(v >> 32));
    mpz_mul_2exp(n, n, 32);
    mpz_add_ui(n, n, (unsigned long)(uint32_t)v);
  }
#endif
}
int _GMP_tinyecm64_factor(const mpz_t n, mpz_t f, int try_small)
{
  uint64_t F64 = 0;
  if (mpz_cmp_ui(n, 4) >= 0 && mpz_sizeinbase(n, 2) <= 64) {
    F64 = tinyecm64(mpz_get64(n), try_small);
    if (F64 > 0)
      mpz_set64(f, F64);
  }
  return (F64 > 0);
}


#define POWER_FACTOR_SMALL_EXP_LIMIT 31
#define POWER_FACTOR_PRECHECK_TRIAL_LIMIT 127UL
#define POWER_FACTOR_EARLY_TRIAL_LIMIT 1000UL
#define POWER_FACTOR_TRIAL_LIMIT 1000000UL

static int _power_factor_trial_root(mpz_t root, const mpz_t cofactor,
                                    const unsigned long *pfactors,
                                    const unsigned long *exponents,
                                    unsigned long nfactors, unsigned long k)
{
  mpz_t pk;
  unsigned long i;

  for (i = 0; i < nfactors; i++)
    if ((exponents[i] % k) != 0)
      return 0;

  mpz_init(pk);
  if (mpz_cmp_ui(cofactor, 1) == 0) {
    mpz_set_ui(root, 1); /* We'll multiply by valuation factors below */
  } else {
    if (!mpz_root(pk, cofactor, k)) {
      mpz_clear(pk);
      return 0;
    }
    /* Only after confirmed will we write into 'root' */
    mpz_set(root, pk);
  }

  for (i = 0; i < nfactors; i++) {
    mpz_ui_pow_ui(pk, pfactors[i], exponents[i]/k);
    mpz_mul(root, root, pk);
  }
  mpz_clear(pk);
  return 1;
}

/* The three valuation helpers return 1 when no nontrivial exponent is found,
 * including when no usable valuation information is available.  This differs
 * from power_factor's return value of 0.  They set f only when returning k > 1.
 * Here g is the gcd of the recorded prime valuations (0 if none), and the
 * original input is cofactor * product(pfactors[i]^exponents[i]). */
static unsigned long _power_factor_from_valuations(const mpz_t cofactor, mpz_t f,
                                                   const unsigned long *pfactors,
                                                   const unsigned long *exponents,
                                                   unsigned long nfactors,
                                                   unsigned long g)
{
  unsigned long k = 1, q = 2, rem = g;

  if (g <= 1)
    return 1;
  if (mpz_cmp_ui(cofactor, 1) == 0) {
    _power_factor_trial_root(f, cofactor, pfactors, exponents, nfactors, g);
    return g;
  }

  {
    PRIME_ITERATOR(iter);
    while (rem > 1 && q <= rem/q) {
      if ((rem % q) == 0) {
        do { rem /= q; } while ((rem % q) == 0);
        while (k <= g/q && (g % (k*q)) == 0) {
          unsigned long cand = k*q;
          if (!_power_factor_trial_root(f, cofactor, pfactors, exponents, nfactors, cand))
            break;
          k = cand;
        }
      }
      q = prime_iterator_next(&iter);
    }
    if (rem > 1) {
      q = rem;
      while (k <= g/q && (g % (k*q)) == 0) {
        unsigned long cand = k*q;
        if (!_power_factor_trial_root(f, cofactor, pfactors, exponents, nfactors, cand))
          break;
        k = cand;
      }
    }
    prime_iterator_destroy(&iter);
  }
  return k;
}

static unsigned long _power_factor_from_large_valuation(const mpz_t n, mpz_t f,
                                                        unsigned long B,
                                                        unsigned long limit)
{
  mpz_t cofactor, zp;
  unsigned long pfactor = 0, exponent = 0;
  unsigned long p, ret = 1;
  PRIME_ITERATOR(iter);

  mpz_init(cofactor);
  mpz_init(zp);

  if (B >= 2 && mpz_even_p(n)) {
    exponent = mpz_scan1(n, 0);
    if (exponent > limit) {
      pfactor = 2;
      mpz_tdiv_q_2exp(cofactor, n, exponent);
      ret = _power_factor_from_valuations(cofactor, f, &pfactor, &exponent, 1, exponent);
      goto end;
    }
  }

  prime_iterator_setprime(&iter, 3);
  for (p = 3; p <= B; p = prime_iterator_next(&iter)) {
    if (mpz_divisible_ui_p(n, p)) {
      mpz_set_ui(zp, p);
      exponent = mpz_remove(cofactor, n, zp);
      if (exponent > limit) {
        pfactor = p;
        ret = _power_factor_from_valuations(cofactor, f, &pfactor, &exponent, 1, exponent);
        break;
      }
    }
  }

  end:
  prime_iterator_destroy(&iter);
  mpz_clear(zp);
  mpz_clear(cofactor);
  return ret;
}

static unsigned long _power_factor_from_trial_factors(const mpz_t n, mpz_t f, unsigned long B)
{
  mpz_t cofactor, zf;
  unsigned long *pfactors = 0, *exponents = 0;
  unsigned long alloc = 0, nfactors = 0;
  unsigned long p, v, g = 0;
  unsigned long k;
  PRIME_ITERATOR(iter);

  mpz_init_set(cofactor, n);
  mpz_init(zf);

  for (p = 2; p <= B && mpz_cmp_ui(cofactor, 1) > 0; p = prime_iterator_next(&iter)) {
    if (mpz_divisible_ui_p(cofactor, p)) {
      if (nfactors == alloc) {
        if (alloc == 0) {
          alloc = 8;
          New(0, pfactors, alloc, unsigned long);
          New(0, exponents, alloc, unsigned long);
        } else {
          alloc += 8;
          Renew(pfactors, alloc, unsigned long);
          Renew(exponents, alloc, unsigned long);
        }
      }
      mpz_set_ui(zf, p);
      v = mpz_remove(cofactor, cofactor, zf);
      pfactors[nfactors] = p;
      exponents[nfactors] = v;
      nfactors++;
      g = (g == 0) ? v : _gcd_ui(g, v);
      if (g == 1)
        break;
    }
  }
  prime_iterator_destroy(&iter);

  k = _power_factor_from_valuations(cofactor, f, pfactors, exponents, nfactors, g);
  if (pfactors != 0)  Safefree(pfactors);
  if (exponents != 0) Safefree(exponents);
  mpz_clear(zf);
  mpz_clear(cofactor);
  return k;
}

/* Return the maximal k > 1 such that n = f^k, or 0 if n <= 1 or n is not
 * a perfect power.  The output f is meaningful only on success and may
 * alias n. */
unsigned long power_factor(const mpz_t n, mpz_t f)
{
  unsigned long k = 1, b = 2;
  if (mpz_cmp_ui(n, 1) <= 0)
    return 0;

  k = _power_factor_from_large_valuation(n, f, POWER_FACTOR_PRECHECK_TRIAL_LIMIT,
                                               POWER_FACTOR_SMALL_EXP_LIMIT);
  if (k > 1)
    return k;

  if (mpz_perfect_power_p(n)) {
    int tried_factor = 0;
    mpz_t nf, tf;
    PRIME_ITERATOR(iter);

    mpz_init_set(nf, n);
    mpz_init(tf);
    k = _power_factor_from_trial_factors(nf, f, POWER_FACTOR_EARLY_TRIAL_LIMIT);
    if (k == 1) {
      while (1) {
        unsigned long ok = k;
        while (mpz_root(tf, nf, b)) {
          mpz_set(f, tf);
          mpz_set(nf, tf);
          k *= b;
        }
        if (ok != k && !mpz_perfect_power_p(nf)) break;
        if (mpz_cmp_ui(tf, 1) <= 0) break; /* Exit if we can't find the power */
        b = prime_iterator_next(&iter);
        if (!tried_factor && b > POWER_FACTOR_SMALL_EXP_LIMIT) {
          unsigned long vk;
          tried_factor = 1;
          /* Known small-factor valuations bound the remaining exponent. */
          vk = _power_factor_from_trial_factors(nf, f, POWER_FACTOR_TRIAL_LIMIT);
          if (vk > 1) {
            k *= vk;
            break;
          }
        }
      }
    }
    mpz_clear(tf);  mpz_clear(nf);
    prime_iterator_destroy(&iter);
  }
  return (k == 1) ? 0 : k;
}


mpz_t * divisor_list(int *num_divisors, const mpz_t n, const mpz_t maxd)
{
  mpz_t *factors, *divs, mult, t;
  int nfactors, ndivisors, i, j, k, count, *exponents;

  if (mpz_sgn(n) <= 0 || mpz_sgn(maxd) <= 0) {
    *num_divisors = 0;
    return 0;
  }

  nfactors = factor(n, &factors, &exponents);
  ndivisors = 1 + ((nfactors > 0) ? exponents[0] : 0);
  for (i = 1; i < nfactors; i++)
    ndivisors *= (exponents[i] + 1);

  mpz_init(mult);  mpz_init(t);
  New(0, divs, ndivisors, mpz_t);
  mpz_init_set_ui(divs[0], 1);
  for (count = 1, k = 0; k < nfactors; k++) {
    int scount = count;
    mpz_set_ui(mult, 1);
    for (j = 0; j < exponents[k]; j++) {
      mpz_mul(mult, mult, factors[k]);
      for (i = 0; i < scount; i++) {
        mpz_mul(t, divs[i], mult);
        if (mpz_cmp(t, maxd) <= 0) {
          mpz_init(divs[count]);
          mpz_mul(divs[count], divs[i], mult);
          count++;
        }
      }
    }
  }
  ndivisors = count;
  mpz_clear(t);  mpz_clear(mult);
  clear_factors(nfactors, &factors, &exponents);

  qsort(divs, ndivisors, sizeof(mpz_t), _mpz_cmp_asc);

  *num_divisors = ndivisors;
  return divs;
}

int is_smooth(const mpz_t n, const mpz_t k) {
  mpz_t *factors;
  mpz_t N;
  uint32_t khi, res;

  if (mpz_cmp_ui(n,1) <= 0) return 1;
  if (mpz_cmp_ui(k,1) <= 0) return 0;
  if (mpz_cmp(n, k) <= 0) return 1;

  mpz_init_set(N, n);
  khi = (mpz_cmp_ui(k, 10000000) >= 0) ? 10000000 : mpz_get_ui(k);
  res = 0;

  {
    void *iter = trial_factor_iterator_create(N, khi);
    unsigned long f;
    uint32_t e;

    while (!res && trial_factor_iterator_next(&f, &e, iter)) {
      while (e-- > 0)
        mpz_divexact_ui(N, N, f);
      if (mpz_cmp(N, k) <= 0)
        res = 1;
    }
    trial_factor_iterator_destroy(iter);
  }

  /* If we need to keep looking for factors larger than khi, full factor. */
  if (!res && mpz_cmp_ui(k, khi) > 0) {
    int i, nfactors, *exponents;
    nfactors = factor(N, &factors, &exponents);
    for (i = 0; !res && i < nfactors; i++) {
      if (mpz_cmp(factors[i], k) > 0)
        res = 1;
    }
    clear_factors(nfactors, &factors, &exponents);
  }

  mpz_clear(N);
  return res;
}

int is_rough(const mpz_t n, const mpz_t k) {
  mpz_t *factors;
  mpz_t N, f;
  int i, nfactors, *exponents;
  uint32_t khi, stage;

  if (mpz_sgn(n) == 0) return 0 + (mpz_sgn(k) == 0);
  if (mpz_cmp_ui(n,1) == 0 || mpz_cmp_ui(k,1) <= 0) return 1;
  if (mpz_cmp_ui(k,2) == 0) return 0 + (mpz_cmp_ui(n,1) >= 0);
  if (mpz_cmp_ui(k,3) == 0) return (mpz_sgn(n) > 0 && mpz_odd_p(n));
  if (mpz_cmp(n, k) < 0) return 0;  /* 1 < n < k all are zero */
  /* k > 3;  n >= k;  n >= 2 */

  if (mpz_even_p(n) || mpz_divisible_ui_p(n,3)) return 0;
  if (mpz_cmp_ui(k,5) <= 0) return 1;
  if (mpz_divisible_ui_p(n,5)) return 0;
  /* k > 5;  n >= k;  n not divisible by 2, 3, or 5 */

  /* 1. Trial division up to a limit */
  khi = (mpz_cmp_ui(k, 1000000) >= 0) ? 1000000 : mpz_get_ui(k);

  if (_GMP_trial_factor(n, 2, khi-1) > 0)
    return 0;
  if (mpz_cmp_ui(k, khi) <= 0)
    return 1;

  /* 2. Try to efficiently pull out small factors */
  mpz_init_set(N, n);
  mpz_init(f);
  for (stage = 1; stage <= 7; stage++) {
    int success = 0;
    if (stage == 5 && _GMP_BPSW(N)) { mpz_clear(f); mpz_clear(N); return 1; }
    switch (stage) {
      case 1: success = _GMP_pminus1_factor(N, f,   13,   130); break;
      case 2: success = _GMP_pminus1_factor(N, f,  150,  1500); break;
      case 3: success = _GMP_pminus1_factor(N, f, 1500, 15000); break;
      case 4: success = _GMP_ECM_FACTOR(N, f,  100, 1); break;
      case 5: success = _GMP_ECM_FACTOR(N, f,  400, 1); break;
      case 6: success = _GMP_ECM_FACTOR(N, f,  800, 1); break;
      case 7: success = _GMP_ECM_FACTOR(N, f, 1600, 1); break;
      default: break;
    }
    if (!success) continue;
    nfactors = factor(f, &factors, &exponents);
    for (i = 0; i < nfactors; i++) {
      if (mpz_cmp(factors[i], k) < 0)
        break;
    }
    clear_factors(nfactors, &factors, &exponents);
    if (i < nfactors) break;       /* Return 0: Found a small factor. */
    mpz_divexact(N, N, f);         /* Divide out the factors all checked */
    if (mpz_cmp(N, k) < 0) break;  /* Return 0: Reduced N to lower than k */
    if (stage >= 5 && _GMP_BPSW(N)) { mpz_clear(f); mpz_clear(N); return 1; }
  }
  mpz_clear(f);
  if (stage <= 7) {
    mpz_clear(N);
    return 0;
  }

  /* 3. Fully factor */
  nfactors = factor(N, &factors, &exponents);
  for (i = 0; i < nfactors; i++) {
    if (mpz_cmp(factors[i], k) < 0)
      break;
  }
  clear_factors(nfactors, &factors, &exponents);
  mpz_clear(N);
  return i >= nfactors;
}

int is_powerful(const mpz_t n, uint32_t k) {
  mpz_t *factors;
  mpz_t N, f;
  int i, nfactors, *exponents;
  uint32_t e, klo, khi, div;

  /* If n is not a positive integer, return 0 regardless of k.*/
  if (mpz_sgn(n) <= 0) return 0;

  if (k <= 1 || mpz_cmp_ui(n,1) == 0) return 1;

  mpz_init_set(N, n);
  e = mpz_scan1(N, 0);
  if (e > 0) {
    if (e < k) { mpz_clear(N); return 0; }
    mpz_tdiv_q_2exp(N, N, e);
    if (mpz_cmp_ui(N,1) == 0) { mpz_clear(N); return 1; }
  }

  /* TODO: https://github.com/trizen/sidef/commit/455eeedf5231e27a0a7a89b6f70055b80d185d17 */
  /* Can also use a factor iterator when we have one */

  for (i = 2; i <= 8; i++) {
    static const unsigned char pr[9] = {0,2,3,5,7,11,13,17,19};
    uint32_t p = pr[i];
    uint32_t pk = (k == 2) ? p*p : (k == 3) ? p*p*p : p*p*p*p;
    if (mpz_divisible_ui_p(N, p) && !mpz_divisible_ui_p(N, pk))
      { mpz_clear(N); return 0; }
  }

  mpz_init(f);
  mpz_root(f, N, k);
  klo = 3;
  khi = (mpz_cmp_ui(f, 1000000) >= 0) ? 1000000 : mpz_get_ui(f);

  while (klo <= khi && (div = _GMP_trial_factor(N, klo, khi)) > 0) {
    mpz_set_ui(f, div);
    if (mpz_remove(N, N, f) < k)
      { mpz_clear(N); mpz_clear(f); return 0; }
    if (mpz_cmp_ui(N,1) == 0 || power_factor(N, f) >= k)
      { mpz_clear(N); mpz_clear(f); return 1; }
    mpz_ui_pow_ui(f, div, 2*k);
    if (mpz_cmp(N,f) < 0)
      { mpz_clear(N); mpz_clear(f); return 0; }
    klo = div+1;
  }
  if (mpz_cmp_ui(N,1) == 0 || power_factor(N, f) >= k)
    { mpz_clear(N); mpz_clear(f); return 1; }
  mpz_ui_pow_ui(f, khi, 2*k);
  if (mpz_cmp(N,f) < 0)
    { mpz_clear(N); mpz_clear(f); return 0; }
  mpz_clear(f);

  /* 3. Fully factor */
  nfactors = factor(N, &factors, &exponents);
  for (i = 0; i < nfactors; i++) {
    if ((uint32_t)exponents[i] < k)
      break;
  }
  clear_factors(nfactors, &factors, &exponents);
  mpz_clear(N);
  return i >= nfactors;
}

int is_almost_prime(uint32_t k, const mpz_t n)
{
#if 0
  if (k == 0)  return (mpz_cmp_ui(n,1) == 0) ? 1 : 0;
  if (k == 1)  return _GMP_is_prime(n) ? 1 : 0;
  if (k == 2)  return is_semiprime(n);
#endif
  mpz_t nr;
  int res;

  if (mpz_sgn(n) <= 0) return 0;

  mpz_init_set(nr, n);
  while (k > 0 && mpz_even_p(nr))           { k--; mpz_divexact_ui(nr,nr,2); }
  while (k > 0 && mpz_divisible_ui_p(nr,3)) { k--; mpz_divexact_ui(nr,nr,3); }
  while (k > 0 && mpz_divisible_ui_p(nr,5)) { k--; mpz_divexact_ui(nr,nr,5); }
  while (k > 0 && mpz_divisible_ui_p(nr,7)) { k--; mpz_divexact_ui(nr,nr,7); }

  /* Optimally, we should be factoring one at a time.  This will let us
   * exit early if we find too many factors, and stop when we get a final
   * semiprime.  This is tedious without an efficient factor iterator. */

  if      (k == 0)  res = mpz_cmp_ui(nr,1) == 0;
  else if (k == 1)  res = _GMP_is_prime(nr);
  else if (k == 2)  res = is_semiprime(nr);
  else              res = bigomega(nr) == k;

  mpz_clear(nr);
  return res ? 1 : 0;
}

/******************************************************************************/

/* Simple treesieve, algorithm from Jens K Andersen. */
static unsigned long _treesieve(mpz_t n, prime_iterator *piter, unsigned long to_n, unsigned long *vextra, uint32_t *nextra)
{
  unsigned long p = 0,  found = 0,  log2n = mpz_sizeinbase(n,2);
  unsigned long* xn;        /* leaves */
  mpz_t* xtree[16+1];       /* the tree (maxdepth = 16) */
  mpz_t* xtemp;
  unsigned int i, j, d, depth, leafsize, nleaves;

  { /* Decide on the tree depth (3-16) and number of leaves (10-31) */
    unsigned int dp = log2n >> 10;
    depth = 0;
    while (dp >>= 1) depth++;
    if (depth < 3) depth = 3;
    if (depth > 16) depth = 16;
  }
  leafsize = log2n / (1U << depth) / 68;
  nleaves = 1 << depth;
  /* printf("log2n %lu  depth %u  leafsize %u  nleaves %u\n",log2n,depth,leafsize,nleaves); */
  //printf("log2n %lu  depth %u  leafsize %u  nleaves %u\n",log2n,depth,leafsize,nleaves);

  New(0, xn, nleaves * leafsize, unsigned long);
  for (d = 0; d <= depth; d++) {
    unsigned int nodes = 1 << (depth - d);
    New(0, xtree[d], nodes, mpz_t);
    for (j = 0; j < nodes; j++)
      mpz_init(xtree[d][j]);
  }
  xtemp = xtree[1];   /* implies mindepth = 3 */

  while (!found && p <= to_n) {
    /* Create nleaves x[0] values, each the product of leafsize primes */
    for (i = 0; i < nleaves && p <= to_n; i++) {
      for (j = 0; j < 4; j++)                  /* Create 4 sub-products */
        mpz_set_ui(xtemp[j], 1);
      for (j = 0; j < leafsize; j++) {
        p = prime_iterator_next(piter);
        if (p > to_n) break;
        xn[i*leafsize+j] = p;
        mpz_mul_ui(xtemp[j&3], xtemp[j&3], p);
      }
      mpz_mul(xtemp[0], xtemp[0], xtemp[1]);   /* Combine for final product*/
      mpz_mul(xtemp[2], xtemp[2], xtemp[3]);
      mpz_mul(xtree[0][i], xtemp[0], xtemp[2]);
    }
    /* Multiply product tree, xtree[depth][0] has nleaves*leafsize product */
    for (d = 1; d <= depth; d++)
      for (i = 0; i < (1U << (depth-d)); i++)
        mpz_mul(xtree[d][i], xtree[d-1][2*i], xtree[d-1][2*i+1]);
    /* Go backwards replacing the products with remainders */
    mpz_tdiv_r(xtree[depth][0], n, xtree[depth][0]);
    for (d = 1; d <= depth; d++)
      for (i = 0; i < (1U << d); i++)
        mpz_tdiv_r(xtree[depth-d][i], xtree[depth-d+1][i>>1], xtree[depth-d][i]);
    /* Search each leaf for divisors */
    for (i = 0; i < nleaves; i++)
      for (j = 0; j < leafsize; j++)
        if (mpz_divisible_ui_p(xtree[0][i], xn[i*leafsize+j])) {
          if (found) vextra[(*nextra)++] = xn[i*leafsize+j];
          else       found = xn[i*leafsize+j];
        }

    /* TODO:
     *   { nf++;  if (!found) found = xn[i*leafsize+j]; }
     * .. if (nf > 1) alloc and fill extra array
     * .. in main routine, store the extra array and counter
     */
  }
  for (d = 0; d <= depth; d++) {
    unsigned int nodes = 1U << (depth - d);
    for (j = 0; j < nodes; j++)
      mpz_clear(xtree[d][j]);
    Safefree(xtree[d]);
  }
  Safefree(xn);
  return found;
}

/******************************************************************************/

#define NFOUND 8

typedef struct {
  mpz_t n;           /* the number we're finding factors in */
  mpz_t gcd;         /* a variable we can use to store a gcd */
  mpz_t t;           /* temporary */
  unsigned long B;   /* the limit they set for trial factoring */
  int step;          /* where we are in the process */
  uint32_t pn;       /* prime number for small primes */
  prime_iterator piter;
  unsigned long found[NFOUND];
  uint32_t ifound;
  uint32_t nfound;
} tf_iterator_t;

void* trial_factor_iterator_create(const mpz_t n, unsigned long B)
{
  tf_iterator_t *iter;
  New(0, iter, 1, tf_iterator_t);
  mpz_init_set(iter->n, n);
  mpz_init(iter->gcd);
  mpz_init(iter->t);
  iter->B = B;
  iter->step = 0;
  prime_iterator_init(&(iter->piter));
  iter->ifound = iter->nfound = 0;
  return iter;
}


#define TFITER_FINAL_N(n) \
  do { \
    if (mpz_cmp_ui(n, B) > 0) return 0; \
    *f = mpz_get_ui(n); \
    if (e) *e = 1; \
    mpz_set_ui(n,1); \
    return 1; \
  } while (0)

#define TFITER_RETURN_UI(p) \
  do { \
    mpz_set_ui(iter->t, p); \
    v = mpz_remove(iter->n, iter->n, iter->t); \
    *f = p; \
    if (e) *e = v; \
    return 1; \
  } while (0)


#define TFITER_PRE_STEP(gcdv) \
  do { \
    unsigned long pstart = primes_small[iter->pn]; \
    if ((pstart) > B || mpz_cmp_ui(iter->n, (pstart)*(pstart)) < 0) return 0; \
    mpz_gcd(iter->gcd, iter->n, gcdv); \
    if (mpz_cmp_ui(iter->gcd, 1) != 0) iter->step += 1; \
    else                               iter->step += 2; \
  } while (0)

#define TFITER_CHK_STEP(pend) \
  do { \
    if (mpz_cmp_ui(iter->gcd,1)==0) croak("fail tf1 step %u gcd",iter->step); \
    if (mpz_sizeinbase(iter->n,2) <= sizeof(unsigned long)*8) { \
      unsigned long un = mpz_get_ui(iter->n); \
      while (iter->pn <= (pend)) { \
        unsigned long p = primes_small[iter->pn++]; \
        if (p > B || un < p*p) \
          TFITER_FINAL_N(iter->n); \
        if ((un % p) == 0) { \
          mpz_divexact_ui(iter->gcd, iter->gcd, p); \
          if (mpz_cmp_ui(iter->gcd, 1) == 0) \
            iter->step++; \
          TFITER_RETURN_UI(p); \
        } \
      } \
    } else { \
      while (iter->pn <= (pend)) { \
        unsigned long p = primes_small[iter->pn++]; \
        if (p > B || mpz_cmp_ui(iter->n, p*p) < 0) \
          TFITER_FINAL_N(iter->n); \
        if (mpz_divisible_ui_p(iter->gcd, p)) { \
          mpz_divexact_ui(iter->gcd, iter->gcd, p); \
          if (mpz_cmp_ui(iter->gcd, 1) == 0) \
            iter->step++; \
          TFITER_RETURN_UI(p); \
        } \
      } \
    } \
    croak("fail tf1 step %u gcd final",iter->step); \
    iter->step++; \
  } while (0)

int trial_factor_iterator_next(unsigned long *f, uint32_t *e, void* ctx)
{
  tf_iterator_t *iter = ctx;
  uint32_t v, log2n;
  unsigned long B = iter->B;
  prime_iterator *piter = &(iter->piter);

  if (iter->nfound > 0) {
    unsigned long p = iter->found[iter->ifound++];
    if (iter->ifound == iter->nfound)
      iter->ifound = iter->nfound = 0;
    TFITER_RETURN_UI(p);
  }

  if (mpz_cmp_ui(iter->n, 1) <= 0)
    return 0;

  if (iter->step == 101) {
    unsigned long un = mpz_get_ui(iter->n);
    unsigned long sqrtn = (unsigned long) sqrt((double)un);
    if (sqrtn*sqrtn > un) sqrtn--;
    if ( (sqrtn+1)*(sqrtn+1) <= un && sqrtn < (1UL << 4*sizeof(unsigned long)) )
      sqrtn++;
    if (sqrtn < B) B = sqrtn;

    while (1) {
      unsigned long p = prime_iterator_next(piter);
      if (p > B)
        break;
      if ((un % p) == 0)
        TFITER_RETURN_UI(p);
    }
    return 0;
  }

  /* Step 0: just starting.  Remove factors of 2. */
  if (iter->step == 0) {
    if (B < 2) return 0;
    iter->step = iter->pn = 1;
    v = mpz_scan1(iter->n, 0);
    if (v > 0) {
      mpz_tdiv_q_2exp(iter->n, iter->n, v);
      *f = 2;
      if (e) *e = v;
      return 1;
    }
  }
  if (iter->step ==  1) {
    if (B < 3) return 0;
    iter->step = iter->pn = 2;
    if (mpz_divisible_ui_p(iter->n, 3))  TFITER_RETURN_UI(3);
  }
  if (iter->step ==  2) {
    if (B < 5) return 0;
    iter->step = iter->pn = 3;
    if (mpz_divisible_ui_p(iter->n, 5))  TFITER_RETURN_UI(5);
  }

  if (iter->step ==  3)  TFITER_PRE_STEP(   _gcd_1k);    /* Prep */
  if (iter->step ==  4)  TFITER_CHK_STEP(   168);        /* Check */

  if (iter->step ==  5)  TFITER_PRE_STEP(   _gcd_4k);    /* Prep */
  if (iter->step ==  6)  TFITER_CHK_STEP(   550);        /* Check */

  if (iter->step ==  7)  TFITER_PRE_STEP(   _gcd_16k);   /* Prep */
  if (iter->step ==  8)  TFITER_CHK_STEP(  1862);        /* Check */

  if (iter->step ==  9)  TFITER_PRE_STEP(   _gcd_32k);   /* Prep */
  if (iter->step == 10)  TFITER_CHK_STEP(  3432);        /* Check */

  if (iter->step == 11)  TFITER_PRE_STEP(   _gcd_64k);   /* Prep */
  if (iter->step == 12)  TFITER_CHK_STEP(  6413);        /* Check */

  if (iter->step == 13) {
    mpz_set_ui(iter->gcd, 1);
    iter->pn = 6413;
    prime_iterator_setprime(piter, 64000);
    iter->step++;
  }

  if (B < 64007U || mpz_cmp_ui(iter->n, 64007U*64007U) < 0)
    return 0;

  if (_GMP_is_prime(iter->n))
    TFITER_FINAL_N(iter->n);

  log2n = mpz_sizeinbase(iter->n,2);

  if (log2n <= sizeof(unsigned long)*8) {
    iter->step = 101;
    return trial_factor_iterator_next(f, e, ctx);
  }

  /* Determine where to stop.  If B^2 <= n, B <= sqrt(n) so stop at B */
  mpz_ui_pow_ui(iter->t, B, 2);
  if (mpz_cmp(iter->t, iter->n) > 0) {
    /* B fits in an unsigned long.  sqrt(n) < B, so it must also fit. */
    mpz_sqrt(iter->t, iter->n);
    B = mpz_get_ui(iter->t);
  }

  /* TODO:
   *    separate this into:
   *      1) product tree (pass in primes?)
   *      2) remainder tree
   *      3) suck out factors
   */

  if (1 && log2n >= 3000) {
    unsigned long p = _treesieve(iter->n, piter, B, iter->found, &(iter->nfound));
    if (p != 0)
      TFITER_RETURN_UI(p);
    return 0;
  }

  while (1) {
    unsigned long p = prime_iterator_next(piter);
    if (p > B)
      break;
    if (mpz_divisible_ui_p(iter->n, p))
      TFITER_RETURN_UI(p);
  }

  return 0;
}

void trial_factor_iterator_n(mpz_t n, void* ctx)
{
  tf_iterator_t *iter = ctx;
  mpz_set(n, iter->n);
}
void trial_factor_iterator_destroy(void* ctx)
{
  tf_iterator_t *iter = ctx;
  mpz_clear(iter->n);
  mpz_clear(iter->gcd);
  mpz_clear(iter->t);
  prime_iterator_destroy(&(iter->piter));
  Safefree(ctx);
}
