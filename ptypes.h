#ifndef MPU_PTYPES_H
#define MPU_PTYPES_H

#ifndef _MSC_VER
#define __STDC_LIMIT_MACROS
#include <stdint.h>
#endif

#ifdef STANDALONE
  #include <limits.h>
  #include <stdio.h>
  #include <stdlib.h>
  typedef unsigned long UV;
  typedef   signed long IV;
  #define UV_MAX ULONG_MAX
  #define UVCONST(x) ((unsigned long)x##UL)
  #define UVuf "lu"
  #define croak(fmt,...)            { printf(fmt,##__VA_ARGS__); exit(1); }
  #define New(id, mem, size, type)  mem = (type*) malloc((size)*sizeof(type))
  #define Newz(id, mem, size, type) mem = (type*) calloc(size, sizeof(type))
  #define Renew(mem, size, type)    mem =(type*)realloc(mem,(size)*sizeof(type))
  #define Safefree(mem)             free((void*)mem)
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

/* Try to determine if we have 64-bit available via uint64_t */
#define HAVE_STD_U64 0
#if defined(UINT64_MAX) && defined(__UINT64_C)
  #if (UINT64_MAX >= __UINT64_C(18446744073709551615))
    #undef HAVE_STD_U64
    #define HAVE_STD_U64 1
  #endif
#endif

#endif

#define MAXBIT        (BITS_PER_WORD-1)
#define NWORDS(bits)  ( ((bits)+BITS_PER_WORD-1) / BITS_PER_WORD )
#define NBYTES(bits)  ( ((bits)+8-1) / 8 )

#define MPUassert(c,text) if (!(c)) { croak("Math::Prime::Util internal error: " text); }

#if defined(__GNUC__)
  #define INLINE inline
#elif defined(_MSC_VER)
  #define INLINE __inline
#else
  #define INLINE
#endif

#endif
