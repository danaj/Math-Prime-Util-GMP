#ifndef MPU_PTYPES_H
#define MPU_PTYPES_H

#ifndef _MSC_VER
#define __STDC_LIMIT_MACROS
#include <stdint.h>
#else
typedef unsigned __int8  uint8_t;
typedef unsigned __int16 uint16_t;
typedef unsigned __int32 uint32_t;
typedef unsigned __int64 uint64_t;
typedef __int64 int64_t;
typedef __int32 int32_t;
typedef __int16 int16_t;
typedef __int8 int8_t;
#endif

#ifndef HAVE_UINT128
#if defined(__SIZEOF_INT128__) && !defined(__CUDACC__)
  #define HAVE_UINT128 1
  typedef unsigned __int128 uint128_t;
  typedef   signed __int128  int128_t;
#elif (__GNUC__ >= 4) && (defined(__x86_64__) || defined(__powerpc64__))
  #if __clang__ && (__clang_major__ > 4 || (__clang_major__ == 4 && __clang_minor__ >= 2))
    #define HAVE_UINT128 1
    typedef unsigned __int128 uint128_t;
    typedef   signed __int128  int128_t;
  #elif __GNUC__ < 4 || (__GNUC__ == 4 && __GNUC_MINOR__ < 4)
    #define HAVE_UINT128 0
  #elif __GNUC__ == 4 && __GNUC_MINOR__ >= 4 && __GNUC_MINOR__ < 6
    #define HAVE_UINT128 1
    typedef unsigned int uint128_t __attribute__ ((__mode__ (TI)));
    typedef   signed int  int128_t __attribute__ ((__mode__ (TI)));
  #else
    #define HAVE_UINT128 1
    typedef unsigned __int128 uint128_t;
    typedef   signed __int128  int128_t;
  #endif
#elif defined(__BITINT_MAXWIDTH__) && __BITINT_MAXWIDTH__ >= 128
  #define HAVE_UINT128 1
  typedef unsigned _BitInt(128) uint128_t;
  typedef   signed _BitInt(128)  int128_t;
#else
  #define HAVE_UINT128 0
#endif
#endif

#ifndef MAYBE_UNUSED
# if defined(__GNUC__) || defined(__clang__)
#  define MAYBE_UNUSED __attribute__((unused))
# else
#  define MAYBE_UNUSED
# endif
#endif

#ifdef STANDALONE
  #include <limits.h>
  #include <stdio.h>
  #include <stdlib.h>
  typedef unsigned long UV;
  typedef   signed long IV;
  typedef        double NV;
  #define UV_MAX ULONG_MAX
  #define UVCONST(x) ((unsigned long)x##UL)
  #define UVuf "lu"
  #define IVdf "ld"
  /* Fatal diagnostics belong on stderr; some callers already supply '\n'. */
  #define croak(fmt,...)            do { \
    const char *const _mpu_croak_fmt = (fmt); \
    const char *_mpu_croak_end = _mpu_croak_fmt; \
    fprintf(stderr, _mpu_croak_fmt,##__VA_ARGS__); \
    while (*_mpu_croak_end != '\0') _mpu_croak_end++; \
    if (_mpu_croak_end == _mpu_croak_fmt || _mpu_croak_end[-1] != '\n') \
      fputc('\n', stderr); \
    exit(3); \
  } while(0)
  /* Empty allocations may return NULL; renewing to zero frees the allocation.
   * The typed macros pass sizeof(type), so size is always nonzero. */
  static MAYBE_UNUSED void *mpu_malloc(size_t count, size_t size)
  {
    void *mem;
    if (count > (size_t)-1 / size)
      croak("Allocation size overflow");
    mem = malloc(count * size);
    if (mem == NULL && count != 0)
      croak("Out of memory");
    return mem;
  }
  static MAYBE_UNUSED void *mpu_calloc(size_t count, size_t size)
  {
    void *mem;
    if (count > (size_t)-1 / size)
      croak("Allocation size overflow");
    mem = calloc(count, size);
    if (mem == NULL && count != 0)
      croak("Out of memory");
    return mem;
  }
  static MAYBE_UNUSED void *mpu_realloc(void *mem, size_t count, size_t size)
  {
    if (count > (size_t)-1 / size)
      croak("Allocation size overflow");
    if (count == 0) {
      free(mem);
      return NULL;
    }
    mem = realloc(mem, count * size);
    if (mem == NULL)
      croak("Out of memory");
    return mem;
  }
  #define New(id, mem, size, type)  ((mem) = (type*) mpu_malloc((size), sizeof(type)))
  #define Newz(id, mem, size, type) ((mem) = (type*) mpu_calloc((size), sizeof(type)))
  #define Renew(mem, size, type)    ((mem) = (type*) mpu_realloc((void*)(mem), (size), sizeof(type)))
  #define Safefree(mem) free((void*)(mem))
  /* iterator using mpz_nextprime, which is really slow
  #define PRIME_ITERATOR(i) mpz_t i; mpz_init_set_ui(i, 2)
  static UV prime_iterator_next(mpz_t *iter) { mpz_nextprime(*iter, *iter); return mpz_get_ui(*iter); }
  static void prime_iterator_destroy(mpz_t *iter) { mpz_clear(*iter); }
  static void prime_iterator_setprime(mpz_t *iter, UV n) {mpz_set_ui(*iter, n);}
  static int prime_iterator_isprime(mpz_t *iter, UV n) {int isp; mpz_t t; mpz_init_set_ui(t, n); isp = mpz_probab_prime_p(t, 10); mpz_clear(t); return isp;}
  */
#if ULONG_MAX >> 31 == 1
  #define BITS_PER_WORD  32
#elif ULONG_MAX >> 63 == 1
  #define BITS_PER_WORD  64
#else
  #error Unsupported bits per word (must be 32 or 64)
#endif

#else

#if defined(__clang__) && defined(__clang_major__) && __clang_major__ > 11
#pragma clang diagnostic ignored "-Wcompound-token-split-by-macro"
#endif

#include "EXTERN.h"
#include "perl.h"

/* From perl.h, wrapped in PERL_CORE */
#ifndef U32_CONST
# if INTSIZE >= 4
#  define U32_CONST(x) ((U32TYPE)x##U)
# else
#  define U32_CONST(x) ((U32TYPE)x##UL)
# endif
#endif

/* From perl.h, wrapped in PERL_CORE */
#ifndef U64_CONST
# ifdef HAS_QUAD
#  if INTSIZE >= 8
#   define U64_CONST(x) ((U64TYPE)x##U)
#  elif LONGSIZE >= 8
#   define U64_CONST(x) ((U64TYPE)x##UL)
#  elif QUADKIND == QUAD_IS_LONG_LONG
#   define U64_CONST(x) ((U64TYPE)x##ULL)
#  else /* best guess we can make */
#   define U64_CONST(x) ((U64TYPE)x##UL)
#  endif
# endif
#endif


#ifdef HAS_QUAD
  #define BITS_PER_WORD  64
  #define UVCONST(x)     U64_CONST(x)
#else
  #define BITS_PER_WORD  32
  #define UVCONST(x)     U32_CONST(x)
#endif

#endif

/* Try to determine if we have 64-bit available via uint64_t */
#if defined(UINT64_MAX) || defined(_UINT64_T) || defined(__UINT64_TYPE__)
  #define HAVE_STD_U64 1
#elif defined(_MSC_VER)   /* We set up the types earlier */
  #define HAVE_STD_U64 1
#else
  #define HAVE_STD_U64 0
#endif

#define MAXBIT        (BITS_PER_WORD-1)
#define NWORDS(bits)  ( ((bits)+BITS_PER_WORD-1) / BITS_PER_WORD )
#define NBYTES(bits)  ( ((bits)+8-1) / 8 )

#define MPUassert(c,text) if (!(c)) { croak("Math::Prime::Util internal error: " text); }

#undef INLINE
#undef RESTRICT
#undef NOINLINE
#undef ISCONSTFUNC

#if defined(__GNUC__)
  #define INLINE __inline__
#elif defined(_MSC_VER)
  #define INLINE __inline
#else
  #define INLINE
#endif

#if defined(_MSC_VER)
  #define RESTRICT __restrict
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
  #define RESTRICT restrict
#elif defined(__GNUC__) || defined(__clang__)
  #define RESTRICT __restrict__
#else
  #define RESTRICT
#endif

#if (defined(__GNUC__) || defined(__clang__)) && !defined(__INTEL_COMPILER)
  #define ISCONSTFUNC __attribute__((const))
  #define NOINLINE __attribute__((noinline))
#else
  #define ISCONSTFUNC
  #define NOINLINE
#endif

#endif
