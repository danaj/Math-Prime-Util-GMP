
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
/* We're not using anything for which we need ppport.h */

#include <string.h>
#include <gmp.h>

#include "ptypes.h"
#include "gmp_main.h"

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

#define TRIAL_LIM 400
#define MAX_FACTORS 128


MODULE = Math::Prime::Util::GMP		PACKAGE = Math::Prime::Util::GMP

PROTOTYPES: ENABLE

void
_GMP_set_verbose(IN int v)

int
_GMP_miller_rabin(IN char* strn, IN char* strbase)
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
is_prob_prime(IN char* strn)
  PREINIT:
    mpz_t n;
  CODE:
    if ((strn != 0) && (strn[0] == '-') )  /* Negative numbers return 0 */
      XSRETURN_IV(0);
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
    if ((strn != 0) && (strn[0] == '-') )  /* Negative numbers return 0 */
      XSRETURN_IV(0);
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


SV*
_GMP_trial_primes(IN char* strlow, IN char* strhigh)
  PREINIT:
    mpz_t low, high, curprime;
    AV* av = newAV();
  CODE:
    validate_string_number("trial_primes (low)", strlow);
    validate_string_number("trial_primes (high)", strhigh);
    mpz_init_set_str(low, strlow, 10);
    mpz_init_set_str(high, strhigh, 10);

    if (mpz_cmp(low, high) <= 0) {
      mpz_t curprime;
      char* str;
      int nsize = mpz_sizeinbase(high, 10) + 2;

      New(0, str, nsize, char);
      if (str == 0)  croak("Could not allocate space for return string");

      mpz_init_set(curprime, low);
      if (mpz_cmp_ui(curprime, 2) >= 0)
        mpz_sub_ui(curprime, curprime, 1);  /* Make sure low gets included */
      _GMP_next_prime(curprime);

      while (mpz_cmp(curprime, high) <= 0) {
        UV v = mpz_get_ui(curprime);     /* mpz_fits_ulong_p is nice, but if */
        if (!mpz_cmp_ui(curprime, v)) {  /* UV_MAX < ULONG_MAX then it fails */
          av_push(av,newSVuv(v));
        } else {
          mpz_get_str(str, 10, curprime);
          av_push(av,newSVpv(str, 0));
        }
        _GMP_next_prime(curprime);
      }
      Safefree(str);
      mpz_clear(curprime);
    }
    mpz_clear(low);
    mpz_clear(high);
    RETVAL = newRV_noinc( (SV*) av );
  OUTPUT:
    RETVAL



#define SIMPLE_FACTOR_START(name) \
    validate_string_number(name " (n)", strn); \
    mpz_init_set_str(n, strn, 10); \
    if (mpz_cmp_ui(n, 3) <= 0) { \
      XPUSH_MPZ(n); \
    } else { \
      /* Skip the trivial division tests */ \
      /* while (mpz_divisible_ui_p(n, 2)) { mpz_divexact_ui(n, n, 2); XPUSHs(sv_2mortal(newSVuv( 2 ))); } */ \
      /* while (mpz_divisible_ui_p(n, 3)) { mpz_divexact_ui(n, n, 3); XPUSHs(sv_2mortal(newSVuv( 3 ))); } */ \
      /* while (mpz_divisible_ui_p(n, 5)) { mpz_divexact_ui(n, n, 5); XPUSHs(sv_2mortal(newSVuv( 5 ))); } */ \
      if (mpz_cmp_ui(n, 1) == 0) { /* done */ } \
      else if (_GMP_is_prime(n)) { XPUSH_MPZ(n); } \
      else { \
        mpz_t f; \
        int success; \
        mpz_init(f);

#define SIMPLE_FACTOR_END \
        if (!success) { \
          XPUSHs(sv_2mortal(newSVpv(strn, 0))); \
        } else { \
          mpz_divexact(n, n, f); \
          XPUSH_MPZ(n); \
          XPUSH_MPZ(f); \
        } \
        mpz_clear(f); \
      } \
    } \
    mpz_clear(n);


void
prho_factor(IN char* strn, IN UV maxrounds = 64*1024*1024)
  PREINIT:
    mpz_t n;
  PPCODE:
    SIMPLE_FACTOR_START("prho_factor");
    success = _GMP_prho_factor(n, f, 3, maxrounds);
    SIMPLE_FACTOR_END;

void
pbrent_factor(IN char* strn, IN UV maxrounds = 64*1024*1024)
  PREINIT:
    mpz_t n;
  PPCODE:
    SIMPLE_FACTOR_START("pbrent_factor");
    success = _GMP_pbrent_factor(n, f, 3, maxrounds);
    SIMPLE_FACTOR_END;

void
pminus1_factor(IN char* strn, IN UV smoothness = 1000000, IN UV B2 = 0)
  PREINIT:
    mpz_t n;
  PPCODE:
    SIMPLE_FACTOR_START("pminus1_factor");
    if (B2 == 0) {
      /* Without a B2, increment up */
      UV B = 5;
      while (!success) {
        success = _GMP_pminus1_factor(n, f, B, 20*B);
        if (B == smoothness) break;
        B = 20*B;
        if (B > smoothness) B = smoothness;
      }
    } else {
      /* Given a B1 and B2, do just what they asked. */
      success = _GMP_pminus1_factor(n, f, smoothness, B2);
    }
    //success = _GMP_pminus1_factor2(n, f, smoothness);
    SIMPLE_FACTOR_END;

void
holf_factor(IN char* strn, IN UV maxrounds = 256*1024*1024)
  PREINIT:
    mpz_t n;
  PPCODE:
    SIMPLE_FACTOR_START("holf_factor");
    success = _GMP_holf_factor(n, f, maxrounds);
    SIMPLE_FACTOR_END;

void
squfof_factor(IN char* strn, IN UV maxrounds = 16*1024*1024)
  PREINIT:
    mpz_t n;
  PPCODE:
    SIMPLE_FACTOR_START("squfof_factor");
    success = _GMP_squfof_factor(n, f, maxrounds);
    SIMPLE_FACTOR_END;

void
_GMP_factor(IN char* strn)
  PREINIT:
    mpz_t n;
  PPCODE:
    validate_string_number("factor (n)", strn);
    mpz_init_set_str(n, strn, 10);
    if (mpz_cmp_ui(n, 4) < 0) {
      XPUSH_MPZ(n);
    } else {
      UV tf;
      mpz_t f;
      mpz_t tofac_stack[MAX_FACTORS];
      int ntofac = 0;

      { /* trial factor.  We could possibly use mpz_remove here. */
        tf = 2;
        while (tf < TRIAL_LIM) {
          tf = _GMP_trial_factor(n, tf, TRIAL_LIM);
          if (tf == 0)  break;
          XPUSHs(sv_2mortal(newSVuv( tf )));
          mpz_divexact_ui(n, n, tf);
          if (mpz_cmp_ui(n, 1) == 0)
            break;
        }
      }

      mpz_init(f);
      do { /* loop over each remaining factor */
        while ( mpz_cmp_ui(n, TRIAL_LIM*TRIAL_LIM) > 0 && !_GMP_is_prime(n) ) {
          int success = 0;
          int o = _GMP_get_verbose();

          /*
           * This set of operations is meant to provide good performance for
           * "random" numbers as input.  Hence we stack lots of effort up front
           * looking for small factors: prho and pbrent are ~ O(f^1/2) where
           * f is the smallest factor.  SQUFOF is O(N^1/4), so arguable not
           * any better.
           * 
           * On my small 32-bit workstation, these will factor a 778-digit
           * number consisting of 101 8-digit factors in under 10 seconds.
           * A 246-digit number with 21 12-digit factors took a little under
           * two minutes.  A 150-digit number consisting p12 * 10 p14's took
           * 4.5 minutes.  Times for all of these would be faster if the input
           * was a single semiprime.  For comparison, Pari took about 1 minute
           * on the same machine to factor the 150-digit number.
           *
           * A half-decent ECM, MPQS, or SIQS factoring method would be a
           * big help in factoring larger numbers.  This code will take far
           * too long to factor numbers where the smallest factor is 17+ digits,
           * while software like gmp-ecm, msieve, or yafu have no troubles.
           *
           * On the plus side, this GMP code is far, far faster than Perl
           * bigint code, and gives an easy way for Perl programs to factor
           * many large numbers easily and quickly.
           */

          success =                _GMP_prho_factor(n, f, 3, 64*1024);
          if (!success)  success = _GMP_prho_factor(n, f, 5, 64*1024);
          if (!success)  success = _GMP_prho_factor(n, f, 7, 64*1024);
          if (!success)  success = _GMP_prho_factor(n, f,11, 64*1024);
          if (!success)  success = _GMP_prho_factor(n, f,13, 64*1024);
          if (success&&o) {gmp_printf("small prho found factor %Zd\n", f);o=0;}

          /* Small p-1 */
          {
            UV B = 5;
            while (!success) {
              success = _GMP_pminus1_factor(n, f, B, 10*B);
              if (B == 5000) break;
              B = 10*B;
              if (B > 5000) B = 5000;
            }
          }
          if (success&&o) {gmp_printf("p-1 (50k) found factor %Zd\n", f);o=0;}

          if (!success)  success = _GMP_pbrent_factor(n, f, 1, 32*1024*1024);
          if (success&&o) {gmp_printf("pbrent (1,32M) found factor %Zd\n", f);o=0;}

          if (!success)  success = _GMP_prho_factor(n, f,17, 32*1024*1024);
          if (success&&o) {gmp_printf("prho (17,32M) found factor %Zd\n", f);o=0;}

          /* HOLF in case it's a near-ratio-of-perfect-square */
          if (!success)  success = _GMP_holf_factor(n, f, 1*1024*1024);
          if (success&&o) {gmp_printf("holf found factor %Zd\n", f);o=0;}

          /* Large p-1 with stage 2: B2 = 20*B1 */
          if (!success)  success = _GMP_pminus1_factor(n, f, 3000000, 3000000*20);
          if (success&&o) {gmp_printf("p-1 (3M) found factor %Zd\n", f);o=0;}

          if (!success)  success = _GMP_pbrent_factor(n, f, 3, 256*1024*1024);
          if (success&&o) {gmp_printf("pbrent (3,256M) found factor %Zd\n", f);o=0;}

          if (!success)  success = _GMP_pbrent_factor(n, f, 2, 512*1024*1024);
          if (success&&o) {gmp_printf("pbrent (2,512M) found factor %Zd\n", f);o=0;}

          if (!success) gmp_printf("starting squfof on %Zd\n", n);
          if (!success)  success = _GMP_squfof_factor(n, f, 256*1024*1024);
          if (success&&o) {gmp_printf("squfof found factor %Zd\n", f);o=0;}

          /* if (!success)  croak("Could not factor n: %s\n", strn); */

          if (success) {
            if (!mpz_divisible_p(n, f) || !mpz_cmp_ui(f, 1) || !mpz_cmp(f, n)) {
              gmp_printf("n = %Zd  f = %Zd\n", n, f);
              croak("Incorrect factoring");
            }
          }

          if (_GMP_is_prime(f)) {
            XPUSH_MPZ(f);
          } else {
            if (ntofac == MAX_FACTORS-1)
              croak("Too many factors\n");
            mpz_init_set(tofac_stack[ntofac], f);
            ntofac++;
          }
          mpz_divexact(n, n, f);
        }
        /* n is now prime or 1 */
        if (mpz_cmp_ui(n, 1) > 0) {
          XPUSH_MPZ(n);
          mpz_set_ui(n, 1);
        }
        if (ntofac-- > 0) {
          mpz_set(n, tofac_stack[ntofac]);
          mpz_clear(tofac_stack[ntofac]);
        }
      } while (mpz_cmp_ui(n, 1) > 0);
      mpz_clear(f);
    }
    mpz_clear(n);
