/* Segmented prime iteration over the uint32_t range.
 * No initialization of shared state is required. Separate iterators may be
 * used concurrently; each owns its optional segment buffer. Do not copy a
 * live iterator or share it between threads without synchronization.
 * Copyright (c) 2026 Dana Jacobsen. See LICENSE for redistribution terms.
 */
#ifndef PR32_ITER_H
#define PR32_ITER_H

#include <stddef.h>
#include <stdint.h>

typedef struct {
  uint32_t p, segment_start, segment_bytes;
  uint8_t *segment;
  int error;
  /* Private word cursor; use seek rather than modifying these fields. */
  uint32_t next_byte;
  uint64_t pending;
  uint32_t word_base;
} pr32_iter_t;

#define PR32_ITER_INIT {1, 0, 0, NULL, 0, UINT32_MAX, 0, 0}
#define PR32_ITER(name) pr32_iter_t name = PR32_ITER_INIT

#ifdef __cplusplus
extern "C" {
#endif

/* Initialize before first use. init does not free an existing buffer;
 * destroy does, and resets the iterator to its initial state. */
void pr32_iter_init(pr32_iter_t *iter);
void pr32_iter_destroy(pr32_iter_t *iter);
/* The next call returns the smallest prime strictly greater than n.
 * seek retains allocated storage and clears a previous allocation error. */
void pr32_iter_seek(pr32_iter_t *iter, uint32_t n);
/* Return the next prime, or zero at exhaustion or allocation failure.
 * error is nonzero only for allocation failure, and stays set until seek
 * or destroy. Initial iteration starts with 2; exhaustion never wraps. */
uint32_t pr32_iter_next(pr32_iter_t *iter);
int pr32_iter_error(const pr32_iter_t *iter);

/* Return a malloc-allocated array of every prime <= limit. Free it with
 * free(). count may be NULL. For limit < 2 return NULL and count 0;
 * otherwise NULL means allocation failure, also with count 0. */
uint32_t *pr32_primes_upto(uint32_t limit, size_t *count);

#ifdef __cplusplus
}
#endif

#endif
