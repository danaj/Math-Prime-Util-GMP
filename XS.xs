
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
/* We're not using anything for which we need ppport.h */

#include <string.h>
#include <gmp.h>

#include "ptypes.h"

/* I think we're going to have to end up sucking in a lot of Math::BigInt::GMP
 * infrastructure.  For now manage everything through strings.
 */

static void validate_string_number(const char* f, const char* s)
{
  const char* p;
  if (s == 0)
    croak("%s: null string pointer as input", f);
  if (*s == 0)
    croak("%s: empty string as input", f);
  p = s;
  while (*p != 0) {
    if (!isdigit(*p))
      croak("%s: input '%s' must be a positive integer", f, s);
    p++;
  }
}

static const unsigned short primes_small[] =
  {0,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
   101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,
   193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,
   293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,
   409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,
   521,523,541,547,557,563,569,571,577,587,593,599,601,607,613,617,619,631,
   641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,
   757,761,769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,
   881,883,887,907,911,919,929,937,941,947,953,967,971,977,983,991,997,1009
  };
#define NPRIMES_SMALL (sizeof(primes_small)/sizeof(primes_small[0]))

#define TRIAL_FACTORS 400


MODULE = Math::Prime::Util::GMP		PACKAGE = Math::Prime::Util::GMP

PROTOTYPES: ENABLE

int
GMP_miller_rabin(IN char* strn, IN char* strbase)
  PREINIT:
    mpz_t n, a;
  CODE:
    validate_string_number("GMP_miller_rabin (n)", strn);
    validate_string_number("GMP_miller_rabin (base)", strbase);
    if (strn[1] == 0) {
      switch (strn[0]) {
        case '2': case '3': case '5': case '7': XSRETURN_IV(1); break;
        case '0': case '1': case '4': case '6': case '8': XSRETURN_IV(0); break;
        default:  break; /* let 9 fall through */
      }
    }
    mpz_init_set_str(n, strn, 10);
    mpz_init_set_str(a, strbase, 10);
    RETVAL = _GMP_miller_rabin(n, a);
    mpz_clear(n); mpz_clear(a);
  OUTPUT:
    RETVAL


int
is_strong_lucas_pseudoprime(IN char* strn)
  PREINIT:
    mpz_t n;
  CODE:
    validate_string_number("is_strong_lucas_pseudoprime (n)", strn);
    if (strn[1] == 0) {
      int q_is_prime = 0;
      switch (strn[0]) {
        case '2': case '3': case '5': case '7': q_is_prime = 1; break;
      }
      XSRETURN_IV(q_is_prime);
    }
    mpz_init_set_str(n, strn, 10);
    RETVAL = _GMP_is_strong_lucas_pseudoprime(n);
    mpz_clear(n);
  OUTPUT:
    RETVAL

int
GMP_trial_div(IN char* strn, IN UV to_n)
  PREINIT:
    mpz_t n;
  CODE:
    /* If this function returns 0, it is a composite.
     * If this function returns 1, it is a prime if n <= to_n*to_n.
     */
    validate_string_number("GMP_trial_div (n)", strn);
    if (strn[1] == 0) {
      int q_is_prime = 0;
      switch (strn[0]) {
        case '2': case '3': case '5': case '7': q_is_prime = 1; break;
      }
      XSRETURN_IV(q_is_prime);
    }
    /* n is >= 10 */
    mpz_init_set_str(n, strn, 10);
    RETVAL = _GMP_trial_div(n, to_n);
    mpz_clear(n);
  OUTPUT:
    RETVAL

int
is_prob_prime(IN char* strn)
  PREINIT:
    mpz_t n;
  CODE:
    validate_string_number("is_prob_prime (n)", strn);
    if (strn[1] == 0) {
      int q_is_prime = 0;
      switch (strn[0]) {
        case '2': case '3': case '5': case '7': q_is_prime = 2; break;
      }
      XSRETURN_IV(q_is_prime);
    }
    mpz_init_set_str(n, strn, 10);
    RETVAL = _GMP_is_prob_prime(n);
    mpz_clear(n);
  OUTPUT:
    RETVAL

int
is_prime(IN char* strn)
  PREINIT:
    mpz_t n;
  CODE:
    validate_string_number("is_prime (n)", strn);
    if (strn[1] == 0) {
      int q_is_prime = 0;
      switch (strn[0]) {
        case '2': case '3': case '5': case '7': q_is_prime = 2; break;
      }
      XSRETURN_IV(q_is_prime);
    }
    mpz_init_set_str(n, strn, 10);
    RETVAL = _GMP_is_prime(n);
    mpz_clear(n);
  OUTPUT:
    RETVAL

#define XPUSH_MPZ(n) \
  { \
    /* Push as a scalar if <= min(ULONG_MAX,UV_MAX), string otherwise */ \
    UV v = mpz_get_ui(n); \
    if (!mpz_cmp_ui(n, v)) { \
      XPUSHs(sv_2mortal(newSVuv( v ))); \
    } else { \
      char* str; \
      int nsize = mpz_sizeinbase(n, 10) + 2; \
      New(0, str, nsize, char); \
      mpz_get_str(str, 10, n); \
      XPUSHs(sv_2mortal(newSVpv(str, 0))); \
      Safefree(str); \
    } \
  }

SV *
next_prime(IN char* strn)
  PREINIT:
    mpz_t n;
    int nsize;
    char* str;
  PPCODE:
    validate_string_number("next_prime (n)", strn);
    mpz_init_set_str(n, strn, 10);

    _GMP_next_prime(n);

    XPUSH_MPZ(n);
    mpz_clear(n);

SV *
prev_prime(IN char* strn)
  PREINIT:
    mpz_t n;
    int nsize;
    char* str;
  PPCODE:
    validate_string_number("prev_prime (n)", strn);
    mpz_init_set_str(n, strn, 10);

    _GMP_prev_prime(n);

    XPUSH_MPZ(n);
    mpz_clear(n);

void
prho_factor(IN char* strn, IN UV maxrounds = 4*1024*1024)
  PREINIT:
    mpz_t n;
  PPCODE:
    validate_string_number("prho_factor (n)", strn);
    mpz_init_set_str(n, strn, 10);
    if (_GMP_is_prime(n)) {
      XPUSHs(sv_2mortal(newSVpv(strn, 0)));
    } else {
      mpz_t f;
      int success;

      mpz_init(f);
      success = _GMP_prho_factor(n, f, 3, maxrounds);
      if (!success) {
        XPUSHs(sv_2mortal(newSVpv(strn, 0)));
      } else {
        mpz_divexact(n, n, f);
        XPUSH_MPZ(n);
        XPUSH_MPZ(f);
      }
      mpz_clear(f);
    }
    mpz_clear(n);

void
pminus1_factor(IN char* strn, IN UV smoothness = 1000000)
  PREINIT:
    mpz_t n;
  PPCODE:
    validate_string_number("pminus1_factor (n)", strn);
    mpz_init_set_str(n, strn, 10);
    if (_GMP_is_prime(n)) {
      XPUSHs(sv_2mortal(newSVpv(strn, 0)));
    } else {
      mpz_t f;
      int success;

      mpz_init(f);
      success = _GMP_pminus1_factor(n, f, smoothness);
      //success = _GMP_pminus1_factor2(n, f, 1*1024*1024);
      if (!success) {
        XPUSHs(sv_2mortal(newSVpv(strn, 0)));
      } else {
        mpz_divexact(n, n, f);
        XPUSH_MPZ(n);
        XPUSH_MPZ(f);
      }
      mpz_clear(f);
    }
    mpz_clear(n);

void
holf_factor(IN char* strn, IN UV maxrounds = 256*1024*1024)
  PREINIT:
    mpz_t n;
  PPCODE:
    validate_string_number("holf_factor (n)", strn);
    mpz_init_set_str(n, strn, 10);
    if (_GMP_is_prime(n)) {
      XPUSHs(sv_2mortal(newSVpv(strn, 0)));
    } else {
      mpz_t f;
      int success;

      mpz_init(f);
      success = _GMP_holf_factor(n, f, maxrounds);
      if (!success) {
        XPUSHs(sv_2mortal(newSVpv(strn, 0)));
      } else {
        mpz_divexact(n, n, f);
        XPUSH_MPZ(n);
        XPUSH_MPZ(f);
      }
      mpz_clear(f);
    }
    mpz_clear(n);

void
squfof_factor(IN char* strn, IN UV rounds = 16*1024*1024)
  PREINIT:
    mpz_t n;
  PPCODE:
    validate_string_number("squfof_factor (n)", strn);
    mpz_init_set_str(n, strn, 10);
    if (_GMP_is_prime(n)) {
      XPUSHs(sv_2mortal(newSVpv(strn, 0)));
    } else {
      mpz_t f;
      int success;

      mpz_init(f);
      success = _GMP_squfof_factor(n, f, rounds);
      if (!success) {
        XPUSHs(sv_2mortal(newSVpv(strn, 0)));
      } else {
        mpz_divexact(n, n, f);
        XPUSH_MPZ(n);
        XPUSH_MPZ(f);
      }
      mpz_clear(f);
    }
    mpz_clear(n);

void
GMP_factor(IN char* strn)
  PREINIT:
    mpz_t n;
    UV f;
  PPCODE:
    validate_string_number("prho_factor (n)", strn);
    mpz_init_set_str(n, strn, 10);
    /* trial factor */
    f = 2;
    while (f < TRIAL_FACTORS) {
      f = _GMP_trial_factor(n, f, TRIAL_FACTORS);
      if (f == 0)  break;
      XPUSHs(sv_2mortal(newSVuv( f )));
      mpz_divexact_ui(n, n, f);
      if (mpz_cmp_ui(n,1) == 0)
        break;
    }
    if (mpz_cmp_ui(n, TRIAL_FACTORS*TRIAL_FACTORS) > 0) {
      mpz_t f;
      mpz_init(f);
      while (!_GMP_is_prime(n)) {
        int success = 0;
        int o=0;
        success =                _GMP_prho_factor(n, f, 3, 64*1024);
        if (!success)  success = _GMP_prho_factor(n, f, 5, 64*1024);
        if (!success)  success = _GMP_prho_factor(n, f, 7, 64*1024);
        if (!success)  success = _GMP_prho_factor(n, f,11, 64*1024);
        if (!success)  success = _GMP_prho_factor(n, f,13, 64*1024);
        if (success&&o) {gmp_printf("small prho found factor %Zd\n", f);o=0;}

        if (!success)  success = _GMP_pminus1_factor(n, f, 200000);
        if (success&&o) {gmp_printf("p-1 found factor %Zd\n", f);o=0;}

        if (!success)  success = _GMP_pbrent_factor(n, f, 1, 32*1024*1024);
        if (success&&o) {gmp_printf("pbrent found factor %Zd\n", f);o=0;}

        if (!success)  success = _GMP_prho_factor(n, f,17, 8*1024*1024);
        if (success&&o) {gmp_printf("prho found factor %Zd\n", f);o=0;}

        if (!success)  success = _GMP_squfof_factor(n, f, 256*1024*1024);
        if (success&&o) {gmp_printf("squfof found factor %Zd\n", f);o=0;}

        if (!success)  croak("Could not factor n: %s\n", strn);

        /* TODO: we have to break up f */
        XPUSH_MPZ(f);
        mpz_divexact(n, n, f);
      }
      mpz_clear(f);
    }
    if (mpz_cmp_ui(n, 1) > 0) {
      XPUSH_MPZ(n);
    }
    mpz_clear(n);
