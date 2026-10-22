/* Independent reference and allocation-failure checks; compile this file
 * alone. Example: cc -O2 -o /tmp/pr32-check xt/dev/siqs/pr32-iter-check.c
 * Add -DPR32_ITER_SEGMENT_BYTES=1 to stress short sieve tails.
 */
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int fail_after = -1;
static size_t live;
static int fail(void) {
  if (fail_after < 0) return 0;
  if (fail_after == 0) return 1;
  fail_after--; return 0;
}
static void *checked_malloc(size_t n) {
  void *p; if (fail()) return NULL;
  p = malloc(n); if (p != NULL) live++; return p;
}
static void *checked_realloc(void *p, size_t n) {
  void *q; int fresh = p == NULL; if (fail()) return NULL;
  q = realloc(p, n); if (q != NULL && fresh) live++; return q;
}
static void checked_free(void *p) {
  if (p != NULL) { assert(live > 0); live--; } free(p);
}
#define malloc checked_malloc
#define realloc checked_realloc
#define free checked_free
#include "../../../siqs/pr32_iter.c"
#undef malloc
#undef realloc
#undef free

static uint32_t reference_next(uint32_t n) {
  uint64_t p, d;
  if (n < 2) return 2;
  for (p = (uint64_t)n + 1; p <= UINT32_MAX; p++) {
    if ((p & 1) == 0) continue;
    for (d = 3; d*d <= p && p%d != 0; d += 2) {}
    if (d*d > p) return (uint32_t)p;
  }
  return 0;
}

int main(void) {
  enum { limit = 2000000 };
  unsigned char *composite = (unsigned char *)calloc(limit+1, 1);
  PR32_ITER(iter);
  uint32_t p, q, expected = 2, state = 123456789;
  size_t count = 0, i; int failure;
  assert(composite != NULL);
  for (p = 2; p*p <= limit; p++) if (!composite[p])
    for (q = p*p; q <= limit; q += p) composite[q] = 1;
  while ((p = pr32_iter_next(&iter)) <= limit) {
    assert(p != 0 && p == expected); count++;
    do { expected++; } while (expected <= limit && composite[expected]);
  }
  assert(count == 148933 && !pr32_iter_error(&iter));
  for (i = 0; i < 500; i++) {
    state ^= state << 13; state ^= state >> 17; state ^= state << 5;
    q = state % (limit - 1000); expected = q + 1;
    if (expected < 2) expected = 2;
    while (composite[expected]) expected++;
    pr32_iter_seek(&iter, q); assert(pr32_iter_next(&iter) == expected);
    /* Keep reading through cached words after repositioning the cursor. */
    {
      unsigned j;
      for (j = 0; j < 16; j++) {
        do { expected++; } while (expected <= limit && composite[expected]);
        assert(expected <= limit && pr32_iter_next(&iter) == expected);
      }
    }
  }
  for (q = 0; q < 512; q++) {
    pr32_iter_seek(&iter, q); assert(pr32_iter_next(&iter) == reference_next(q));
  }
  pr32_iter_seek(&iter, UINT32_MAX - 20000U);
  q = UINT32_MAX - 20000U;
  do { expected = reference_next(q); p = pr32_iter_next(&iter); assert(p == expected); q = p; } while (p != 0);
  assert(pr32_iter_next(&iter) == 0 && !pr32_iter_error(&iter));
  for (i = 0; i < 200; i++) {
    state ^= state << 13; state ^= state >> 17; state ^= state << 5;
    pr32_iter_seek(&iter, state); assert(pr32_iter_next(&iter) == reference_next(state));
  }
  pr32_iter_destroy(&iter); assert(live == 0);
  {
    static const uint32_t limits[] = {
      719, 720, 721, CACHE_LIMIT-1, CACHE_LIMIT, CACHE_LIMIT+1,
      CACHE_LIMIT+239, CACHE_LIMIT+240, 1000000, 2000000
    };
    size_t j;
    /* Exercise cutoffs inside words and across static/dynamic boundaries. */
    for (j = 0; j < 512 + sizeof(limits)/sizeof(*limits); j++) {
      uint32_t end = j < 512 ? (uint32_t)j : limits[j-512];
      uint32_t *primes = pr32_primes_upto(end, &count);
      size_t used = 0;
      for (p = 2; p <= end; p++) if (!composite[p]) {
        assert(primes != NULL && used < count && primes[used] == p);
        used++;
      }
      assert(used == count);
      checked_free(primes);
      assert(live == 0);
    }
  }
  for (failure = 0; failure <= 16; failure++) {
    uint32_t *primes; fail_after = failure; count = 99;
    primes = pr32_primes_upto(200000, &count);
    if (primes == NULL) assert(count == 0);
    else {
      assert(count == 17984); expected = 2;
      for (i = 0; i < count; i++) { assert(primes[i] == expected); do { expected++; } while (composite[expected]); }
      checked_free(primes);
    }
    assert(live == 0);
  }
  fail_after = 0; pr32_iter_seek(&iter, CACHE_LIMIT);
  assert(pr32_iter_next(&iter) == 0 && pr32_iter_error(&iter));
  fail_after = -1; pr32_iter_seek(&iter, CACHE_LIMIT);
  assert(pr32_iter_next(&iter) == reference_next(CACHE_LIMIT));
  pr32_iter_destroy(&iter); assert(live == 0);
  assert(pr32_primes_upto(0, &count) == NULL && count == 0);
  assert(pr32_primes_upto(1, &count) == NULL && count == 0);
  { uint32_t *primes = pr32_primes_upto(2, NULL); assert(primes && primes[0] == 2); checked_free(primes); }
  assert(live == 0); free(composite);
  puts("PASS pr32_iter: dense reference, random seeks, UINT32_MAX, bulk arrays, allocation failures and retry");
  return 0;
}
