/* Internal SIQS sieve template; not a separately compiled unit.
 * Included twice by siqs.c: ordinary byte addition, then saturation.
 * Keep the progression/root/padding logic identical in both variants. */

#ifdef SIQS_SATURATING_KERNEL
#define siqs_sieve_add siqs_sieve_add_sat
#define siqs_sieve_one_root siqs_sieve_one_root_sat
#define siqs_sieve_two_roots siqs_sieve_two_roots_sat
#define siqs_sieve_large siqs_sieve_large_sat
#define siqs_sieve_one_hit_local siqs_sieve_one_hit_local_sat
#define siqs_sieve_two_hit_local siqs_sieve_two_hit_local_sat
#define siqs_sieve_two_hit_range siqs_sieve_two_hit_range_sat
#define siqs_run_sieve_kernel siqs_run_sieve_kernel_sat
#endif

static INLINE void siqs_sieve_add(uint8_t *cell, uint8_t logp) {
#ifdef SIQS_SATURATING_KERNEL
  uint32_t value = (uint32_t)*cell + logp;
  /* Logical cells are <=255 and weights <=127: the carry is zero or one.
   * Expand it to an all-one mask, then clamp without a data-dependent branch. */
  *cell = (uint8_t)((value | (0U - (value >> SIQS_SIEVE_SCORE_BITS))) & SIQS_SIEVE_SCORE_MAX);
#elif defined(SIQS_DEBUG)
  uint32_t value = (uint32_t)*cell + logp;
  if (value > SIQS_SIEVE_SCORE_MAX)
    croak("SIQS: byte sieve score overflow");
  *cell = (uint8_t)value;
#else
  *cell += logp;
#endif
}

static INLINE void siqs_sieve_one_root(uint8_t *sieve, uint32_t length,
                                       uint32_t root, uint32_t p,
                                       uint8_t logp) {
  size_t pos = root, gap = p, p3 = (size_t)3 * p;
  while (pos + p3 < length) {
    siqs_sieve_add(sieve + pos, logp);  pos += gap;
    siqs_sieve_add(sieve + pos, logp);  pos += gap;
    siqs_sieve_add(sieve + pos, logp);  pos += gap;
    siqs_sieve_add(sieve + pos, logp);  pos += gap;
  }
  if (pos < length) { siqs_sieve_add(sieve + pos, logp);  pos += gap; }
  if (pos < length) { siqs_sieve_add(sieve + pos, logp);  pos += gap; }
  if (pos < length)   siqs_sieve_add(sieve + pos, logp);
}

/* Alternating-root-gap sieve structure inspired by PARI/GP's MPQS sieve. */
static INLINE void siqs_sieve_two_roots(uint8_t *sieve, uint32_t length,
                                        uint32_t root1, uint32_t root2,
                                        uint32_t p, uint8_t logp) {
  /* clang vectorizer likes size_t here */
  size_t final_offset;
  size_t pos = root1;
  size_t gap1 = root2 - root1;
  size_t gap2 = p - gap1;
  /* The eighth store is at pos + 3*p + gap1, not the next-cycle pos + 4*p. */
  final_offset = (size_t)3 * p + gap1;
  while (pos + final_offset < length) {
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

/* Fixed-hit large-prime tiers and their padded stores are inspired by JML's
 * sieving.  If p > length/count, a root has at most count locations in the
 * logical sieve.  The allocation below makes the final stores safe even when
 * they fall beyond that sieve, eliminating a position-range test per store.
 * Six locations over our combined +/- interval correspond to JML's three
 * locations in each separately processed half. */
static INLINE void siqs_sieve_large(uint8_t *sieve, uint32_t length,
                                     uint32_t root1, uint32_t root2,
                                     uint32_t p, uint8_t logp,
                                     uint32_t count) {
  uint32_t i;
  size_t pos = root1;
  for (i = 0; i < count; i++, pos += p) {
#ifdef SIQS_DEBUG
    if (pos < length)
      siqs_sieve_add(sieve + pos, logp);
#else
    siqs_sieve_add(sieve + pos, logp);
#endif
  }
  if (root2 != root1) {
    pos = root2;
    for (i = 0; i < count; i++, pos += p) {
#ifdef SIQS_DEBUG
      if (pos < length)
        siqs_sieve_add(sieve + pos, logp);
#else
      siqs_sieve_add(sieve + pos, logp);
#endif
    }
  }
#ifndef SIQS_DEBUG
  (void)length;
#endif
}

#define SIQS_ONE_HIT_LOCAL_MIN_PRIMES  5200U
#define SIQS_TWO_HIT_LOCAL_LENGTH_THRESHOLD (64U * 1024U)
#define SIQS_ONE_HIT_SINK_STRIPE_SIZE  4096U

/* Apply one root from the one-hit tier.  In release, redirect an
 * out-of-range store to the caller-supplied padding byte. */
static INLINE void siqs_sieve_one_hit_local(uint8_t *sieve,
                                            uint8_t logp, uint32_t root,
                                            uint32_t length, uint32_t sink) {
#ifdef SIQS_DEBUG
  if (root < length)
    siqs_sieve_add(sieve + root, logp);
  (void)sink;
#else
  uint32_t pos = root < length ? root : sink;
  siqs_sieve_add(sieve + pos, logp);
#endif
}

/* For length/2 < p <= length, the first hit is always logical.  Localize
 * only the optional second hit, using the same padding stripe as one-hit
 * primes.  Native-width addition keeps the address calculation unwrapped. */
static INLINE void siqs_sieve_two_hit_local(uint8_t *sieve,
                                            uint8_t logp, uint32_t root,
                                            uint32_t p, uint32_t length,
                                            uint32_t sink) {
  size_t pos = (size_t)root + p;
  siqs_sieve_add(sieve + root, logp);
#ifdef SIQS_DEBUG
  if (pos < length)
    siqs_sieve_add(sieve + pos, logp);
  (void)sink;
#else
  pos = pos < length ? pos : sink;
  siqs_sieve_add(sieve + pos, logp);
#endif
}

/* Keep this loop out of the generic kernel so small-interval loops retain
 * their own compiler unrolling and register allocation.  Called once per
 * eligible sieve, not once per prime. */
static NOINLINE size_t siqs_sieve_two_hit_range(
    uint8_t *RESTRICT sieve, uint32_t length,
    const uint32_t *RESTRICT prime, const uint32_t *RESTRICT root1,
    const uint32_t *RESTRICT root2, const uint8_t *RESTRICT logp,
    size_t first, size_t end) {
  size_t i;
  uint32_t sink = length;
  for (i = first; i < end && prime[i] <= length; i++) {
    siqs_sieve_two_hit_local(sieve, logp[i], root1[i], prime[i], length, sink);
    if (root1[i] != root2[i])
      siqs_sieve_two_hit_local(sieve, logp[i], root2[i], prime[i], length,
                               sink + 1U);
    sink += 2U;
    if (sink == length + SIQS_ONE_HIT_SINK_STRIPE_SIZE)
      sink = length;
  }
  return i;
}

/* Add prime contributions to already initialized logical sieve bytes.
 * Initialization stays with the caller: the sparse suffix must augment,
 * not erase, the completed dense blocks.  Expose invariant arrays directly
 * to the hot per-prime loops. */
static void siqs_run_sieve_kernel(
    uint8_t *RESTRICT sieve, uint32_t length,
    const uint32_t *RESTRICT prime, const uint32_t *RESTRICT root1,
    const uint32_t *RESTRICT root2, const uint8_t *RESTRICT logp,
    size_t first, size_t end) {
  size_t i;
  for (i = first; i < end && prime[i] <= length / 6U;
       i++) {
    if (root2[i] == root1[i]) {
      siqs_sieve_one_root(sieve, length, root1[i], prime[i], logp[i]);
    } else {
      siqs_sieve_two_roots(sieve, length, root1[i], root2[i], prime[i],
                           logp[i]);
    }
  }
# define SIQS_SIEVE_LARGE_RANGE(bound, count) do {                         \
    for (; i < end && prime[i] <= (bound); i++) {                           \
      siqs_sieve_large(sieve, length, root1[i], root2[i], prime[i],         \
                       logp[i], (count));                                   \
    }                                                                       \
  } while (0)
  SIQS_SIEVE_LARGE_RANGE(length / 5U, 6U);
  SIQS_SIEVE_LARGE_RANGE(length / 4U, 5U);
  SIQS_SIEVE_LARGE_RANGE(length / 3U, 4U);
  SIQS_SIEVE_LARGE_RANGE(length / 2U, 3U);
  /* Local padding is worthwhile only when ordinary dummy stores scatter
   * widely.  M1 Pro screens favored leaving intervals through 64 KiB alone;
   * this also keeps the default balanced dense blocks on their old path. */
  if (length <= SIQS_TWO_HIT_LOCAL_LENGTH_THRESHOLD) {
    SIQS_SIEVE_LARGE_RANGE(length, 2U);
  } else {
    i = siqs_sieve_two_hit_range(sieve, length, prime, root1, root2, logp,
                                i, end);
  }
  if (end - i < SIQS_ONE_HIT_LOCAL_MIN_PRIMES) {
    SIQS_SIEVE_LARGE_RANGE(UINT32_MAX, 1U);
  } else {
    /* Once the one-hit tier is large, keep its out-of-range stores in a
     * cache-local padding stripe.  Cycling through 4k avoids a
     * single-address read-modify-write dependency chain.
     * The minimum sieve allocation is twice an 8K logical interval,
     * so every sink position is inside the existing padding. */
    uint32_t sink = length;
    for (; i < end; i++) {
      siqs_sieve_one_hit_local(sieve, logp[i], root1[i], length, sink);
      if (root1[i] != root2[i])
        siqs_sieve_one_hit_local(sieve, logp[i], root2[i], length, sink + 1U);
      sink += 2U;
      if (sink == length + SIQS_ONE_HIT_SINK_STRIPE_SIZE)
        sink = length;
    }
  }
# undef SIQS_SIEVE_LARGE_RANGE
}


#ifdef SIQS_SATURATING_KERNEL
#undef siqs_sieve_add
#undef siqs_sieve_one_root
#undef siqs_sieve_two_roots
#undef siqs_sieve_large
#undef siqs_sieve_one_hit_local
#undef siqs_sieve_two_hit_local
#undef siqs_sieve_two_hit_range
#undef siqs_run_sieve_kernel
#endif
