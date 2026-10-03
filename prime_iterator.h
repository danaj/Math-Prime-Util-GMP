#ifndef MPU_PITERATOR_H
#define MPU_PITERATOR_H

#include "ptypes.h"

typedef struct {
  UV p;
  UV segment_start;
  UV segment_bytes;
  const unsigned char* segment_mem;
} prime_iterator;

#define PRIME_ITERATOR(i)  prime_iterator i = {2, 0, 0, 0}

/* Shared-cache lifetime belongs to the host, not to each iterator/SIQS call.
 * Call startup once before any use (and before launching caller threads).
 * Do not repeat startup while live; it is not reference-counted. Shutdown
 * only after all users finish, iterators are destroyed, and threads joined.
 * Neither operation is safe to race with cache readers or with the other.
 * Separate iterators can use the stable cache concurrently; an individual
 * iterator's mutable state must not be shared unsynchronized. MPU-GMP's
 * normal host initialization/destruction already manages this cache. */
extern void prime_iterator_global_startup(void);
extern void prime_iterator_global_shutdown(void);

extern void prime_iterator_init(prime_iterator *iter);
extern void prime_iterator_destroy(prime_iterator *iter);
extern UV prime_iterator_next(prime_iterator *iter);
extern void prime_iterator_setprime(prime_iterator *iter, UV n);
extern int prime_iterator_isprime(prime_iterator *iter, UV n);

extern UV* sieve_to_n(UV n, UV* count);
extern unsigned long* sieve_to_n_ui(unsigned long n, unsigned long* count);

#endif
