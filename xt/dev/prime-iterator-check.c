/* Internal allocation-boundary checks, not part of the package or test suite.
 * Build from the repository root, then run each resulting executable:
 * cc -O1 -g -fsanitize=address,undefined xt/dev/prime-iterator-check.c -lm -o /tmp/prime-iterator-check
 * cc -O1 -g -fsanitize=address,undefined -DTEST_ALLOC32 xt/dev/prime-iterator-check.c -lm -o /tmp/prime-iterator-check32
 * TEST_ALLOC32 models a 32-bit allocation limit while retaining 64-bit UV.
 * Copyright (c) 2026 Dana Jacobsen. See LICENSE for redistribution terms. */
#include <setjmp.h>
#include <stdarg.h>
#include <string.h>
#define STANDALONE
#include "../../ptypes.h"

#ifdef TEST_ALLOC32
#if BITS_PER_WORD != 64
#error TEST_ALLOC32 requires a 64-bit UV build
#endif
#undef SIZE_MAX
#define SIZE_MAX UINT32_MAX
#endif

static jmp_buf error_trap;
static int expect_error;
static unsigned allocation_calls;
static char error_text[160];

static void test_croak(const char *format, ...)
{
  va_list args;
  va_start(args, format);
  vsnprintf(error_text, sizeof(error_text), format, args);
  va_end(args);
  if (expect_error) longjmp(error_trap, 1);
  fprintf(stderr, "Unexpected croak: %s\n", error_text);
  abort();
}

static void check(int condition, const char *message)
{
  if (condition) return;
  fprintf(stderr, "FAIL: %s\n", message);
  abort();
}

static void *test_allocate(UV count, size_t size)
{
  void *memory;
  allocation_calls++;
  check(count <= SIZE_MAX / size, "oversized count reached allocator");
  memory = malloc((size_t)count * size);
  check(memory != NULL || count == 0, "allocation failed");
  return memory;
}

#undef croak
#define croak(...) do { test_croak(__VA_ARGS__); abort(); } while (0)
#undef New
#define New(id, memory, count, type) ((memory) = (type *)test_allocate((count), sizeof(type)))
#include "../../prime_iterator.c"

#ifdef TEST_ALLOC32
static void check_overflow(int sieve)
{
  allocation_calls = 0;
  expect_error = 1;
  if (setjmp(error_trap) == 0) {
    if (sieve)
      (void)sieve_erat30((UV)(UINT64_C(30) << 32));
    else
      (void)sieve_to_n((UV)UINT64_C(17962970571202560), NULL);
    check(0, "oversized request did not croak");
  }
  expect_error = 0;
  check(allocation_calls == 0, "overflow must croak before allocating");
  check(strcmp(error_text, sieve ? "prime iterator: sieve allocation size overflow"
                                 : "prime iterator: prime-list allocation size overflow") == 0,
        "overflow diagnostic");
}
#endif

int main(void)
{
  UV count, i, *primes;
  unsigned long ui_count, *ui_primes;
  prime_iterator iter;
#ifdef TEST_ALLOC32
  check_overflow(0);
  check_overflow(1);
#endif
  prime_iterator_global_startup();
  primes = sieve_to_n(1000, &count);
  check(count == 168 && primes[167] == 997, "cached prime list");
  Safefree(primes);
  ui_primes = sieve_to_n_ui(1000, NULL);
  check(ui_primes[167] == 997, "unsigned-long cached list, null count");
  Safefree(ui_primes);
  primes = sieve_to_n(1000000, &count);
  check(count == 78498 && primes[78497] == 999983, "uncached prime list");
  ui_primes = sieve_to_n_ui(1000000, &ui_count);
  check(ui_count == count, "unsigned-long prime count");
  for (i = 0; i < count; i++)
    check(ui_primes[i] == primes[i], "unsigned-long prime-list conversion");
  Safefree(ui_primes);
  Safefree(primes);
  prime_iterator_init(&iter);
  check(prime_iterator_next(&iter) == 3, "first iterator step");
  for (i = 0; i < 2; i++) {
    prime_iterator_setprime(&iter, i);
    check(prime_iterator_next(&iter) == 2, "seek from 0 or 1 wraps to index zero");
  }
  prime_iterator_setprime(&iter, 1000);
  check(prime_iterator_next(&iter) == 1009, "iterator seek");
  check(!prime_iterator_isprime(&iter, 1001), "composite rejection");
  prime_iterator_destroy(&iter);
  prime_iterator_global_shutdown();
  puts("PASS prime iterator allocation and correctness checks");
  return 0;
}
