
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
/* We're not using anything for which we need ppport.h */

#include <string.h>
#include <gmp.h>

#include "ptypes.h"
#include "gmp_main.h"
#include "small_factor.h"
#include "ecm.h"
#include "simpqs.h"
#include "bls75.h"
#include "ecpp.h"
#include "utility.h"
#define _GMP_ECM_FACTOR _GMP_ecm_factor_projective

/* Instead of trying to suck in lots of Math::BigInt::GMP and be terribly
 * clever (and brittle), just do all C<->Perl bigints via strings.  It's
 * crude but seems to work pretty well.
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

#define VALIDATE_AND_SET(func, var, str) \
  validate_string_number(func " (" #var ")", str); \
  mpz_init_set_str(var, str, 10);

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
  PPCODE:
     set_verbose_level(v);

void
_GMP_init()

void
_GMP_destroy()


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

#define PRIMALITY_START(name, small_retval) \
    /* Negative numbers return 0 */ \
    if ((strn != 0) && (strn[0] == '-') ) \
      XSRETURN_IV(0); \
    validate_string_number(name " (n)", strn); \
    if (strn[1] == 0) { \
      int q_is_prime = 0; \
      switch (strn[0]) { \
        case '2': case '3': case '5': case '7': q_is_prime = small_retval; \
                                                break; \
      } \
      XSRETURN_IV(q_is_prime); \
    } \
    mpz_init_set_str(n, strn, 10);

int
is_lucas_pseudoprime(IN char* strn)
  ALIAS:
    is_strong_lucas_pseudoprime = 1
    is_extra_strong_lucas_pseudoprime = 2
  PREINIT:
    mpz_t n;
  CODE:
    if ((strn != 0) && (strn[0] == '-') )
      croak("Parameter '%s' must be a positive integer\n", strn);
    PRIMALITY_START("is_lucas_pseudoprime", 1);
    if (ix == 2)
      RETVAL = _GMP_is_extra_strong_lucas_pseudoprime(n);
    else
      RETVAL = _GMP_is_lucas_pseudoprime(n, ix);
    mpz_clear(n);
  OUTPUT:
    RETVAL

int
is_prime(IN char* strn)
  ALIAS:
    is_prob_prime = 1
    is_aks_prime = 2
    is_nminus1_prime = 3
    is_ecpp_prime = 4
  PREINIT:
    mpz_t n;
    int ret;
  CODE:
    PRIMALITY_START("is_prime", 2);
    switch (ix) {
      case 0: ret = _GMP_is_prime(n); break;
      case 1: ret = _GMP_is_prob_prime(n); break;
      case 2: ret = _GMP_is_aks_prime(n); break;
      case 3: ret = _GMP_primality_bls_nm1(n, 100, 0); break;
      case 4: ret = _GMP_ecpp(n, 0); break;
      default: croak("is_prime: Unknown function alias"); break;
    }
    RETVAL = ret;
    mpz_clear(n);
  OUTPUT:
    RETVAL


void
_is_provable_prime(IN char* strn, IN int wantproof = 0)
  PREINIT:
    int result;
    mpz_t n;
  PPCODE:
    PRIMALITY_START("is_provable_prime", 2);
    if (wantproof == 0) {
      result = _GMP_is_provable_prime(n, 0);
      XPUSHs(sv_2mortal(newSViv( result )));
    } else {
      char* prooftext = 0;
      result = _GMP_is_provable_prime(n, &prooftext);
      XPUSHs(sv_2mortal(newSViv( result )));
      if (prooftext) {
        XPUSHs(sv_2mortal(newSVpv(prooftext, 0)));
        Safefree(prooftext);
      } else {
        XPUSHs(sv_2mortal(newSVpv("", 0)));
      }
    }
    mpz_clear(n);

int
_validate_ecpp_curve(IN char* stra, IN char* strb, IN char* strn, IN char* strpx, IN char* strpy, IN char* strm, IN char* strq)
  PREINIT:
    mpz_t a, n, px, py, m, q, t1, t2;
  CODE:
    VALIDATE_AND_SET("_validate_ecpp_curve", a, stra);
    /* ignore b */
    VALIDATE_AND_SET("_validate_ecpp_curve", n, strn);
    VALIDATE_AND_SET("_validate_ecpp_curve", px, strpx);
    VALIDATE_AND_SET("_validate_ecpp_curve", py, strpy);
    VALIDATE_AND_SET("_validate_ecpp_curve", m, strm);
    VALIDATE_AND_SET("_validate_ecpp_curve", q, strq);
    mpz_init(t1);  mpz_init(t2);
    RETVAL = (ecpp_check_point(px, py, m, q, a, n, t1, t2) == 2) ? 1 : 0;
    mpz_clear(t1); mpz_clear(t2);
    mpz_clear(a); mpz_clear(n); mpz_clear(px); mpz_clear(py);
    mpz_clear(m); mpz_clear(q);
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

void
next_prime(IN char* strn)
  ALIAS:
    prev_prime = 1
  PREINIT:
    mpz_t n;
  PPCODE:
    VALIDATE_AND_SET("next_prime", n, strn);

    if (ix == 0) _GMP_next_prime(n);
    else         _GMP_prev_prime(n);

    XPUSH_MPZ(n);
    mpz_clear(n);


void
prime_count(IN char* strlow, IN char* strhigh)
  PREINIT:
    mpz_t low, high, count;
  PPCODE:
    VALIDATE_AND_SET("prime_count", low, strlow);
    VALIDATE_AND_SET("prime_count", high, strhigh);
    mpz_init_set_ui(count, 0);

    if (mpz_cmp(low, high) <= 0) {
      mpz_t curprime;
      mpz_init_set(curprime, low);
      if (mpz_cmp_ui(curprime, 2) >= 0)
        mpz_sub_ui(curprime, curprime, 1);  /* Make sure low gets included */
      _GMP_next_prime(curprime);
      while (mpz_cmp(curprime, high) <= 0) {
        mpz_add_ui(count, count, 1);
        _GMP_next_prime(curprime);
      }
      mpz_clear(curprime);
    }
    XPUSH_MPZ(count);
    mpz_clear(count);
    mpz_clear(high);
    mpz_clear(low);


void
consecutive_integer_lcm(IN UV B)
  PREINIT:
    mpz_t m;
  PPCODE:
    mpz_init(m);
    _GMP_lcm_of_consecutive_integers(B, m);
    XPUSH_MPZ(m);
    mpz_clear(m);

void
primorial(IN char* strn)
  PREINIT:
    mpz_t prim, n;
  PPCODE:
    VALIDATE_AND_SET("primorial", n, strn);
    mpz_init(prim);
    _GMP_primorial(prim, n);
    XPUSH_MPZ(prim);
    mpz_clear(n);
    mpz_clear(prim);

void
pn_primorial(IN UV n)
  PREINIT:
    mpz_t prim;
  PPCODE:
    mpz_init(prim);
    _GMP_pn_primorial(prim, n);
    XPUSH_MPZ(prim);
    mpz_clear(prim);

SV*
_GMP_trial_primes(IN char* strlow, IN char* strhigh)
  PREINIT:
    mpz_t low, high;
    AV* av = newAV();
  CODE:
    VALIDATE_AND_SET("trial_primes", low, strlow);
    VALIDATE_AND_SET("trial_primes", high, strhigh);

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
      else if (_GMP_is_prob_prime(n)) { XPUSH_MPZ(n); } \
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
trial_factor(IN char* strn, IN UV maxn = 0)
  PREINIT:
    mpz_t n;
    UV factor;
  PPCODE:
    SIMPLE_FACTOR_START("trial_factor");
    factor = _GMP_trial_factor(n, 2, (maxn == 0) ? 2147483647 : maxn);
    mpz_set_ui(f, factor);
    success = (factor != 0);
    SIMPLE_FACTOR_END;

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
pminus1_factor(IN char* strn, IN UV B1 = 5000000, IN UV B2 = 0)
  PREINIT:
    mpz_t n;
  PPCODE:
    SIMPLE_FACTOR_START("pminus1_factor");
    success = _GMP_pminus1_factor(n, f, B1, (B2 == 0) ? B1*10 : B2);
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
ecm_factor(IN char* strn, IN UV bmax = 0, IN UV ncurves = 100)
  PREINIT:
    mpz_t n;
  PPCODE:
    SIMPLE_FACTOR_START("ecm_factor");
    if (bmax == 0) {
                    success = _GMP_ecm_factor_projective(n, f,     1000, 40);
      if (!success) success = _GMP_ecm_factor_projective(n, f,    10000, 40);
      if (!success) success = _GMP_ecm_factor_projective(n, f,   100000, 40);
      if (!success) success = _GMP_ecm_factor_projective(n, f,  1000000, 40);
      if (!success) success = _GMP_ecm_factor_projective(n, f, 10000000, 100);
    } else {
      success = _GMP_ecm_factor_projective(n, f, bmax, ncurves);
    }
    SIMPLE_FACTOR_END;

void
qs_factor(IN char* strn)
  PREINIT:
    mpz_t n;
  PPCODE:
    /* Returns multiple factors, so do this separately */
    VALIDATE_AND_SET("qs_factor", n, strn);
    if (mpz_cmp_ui(n, 3) <= 0) {
      XPUSH_MPZ(n);
    } else {
      if (_GMP_is_prob_prime(n)) {
        XPUSH_MPZ(n);
      } else {
        mpz_t farray[66];
        int i, nfactors;
        for (i = 0; i < 66; i++)
          mpz_init(farray[i]);
        nfactors = _GMP_simpqs(n, farray);
        for (i = 0; i < nfactors; i++) {
          XPUSH_MPZ(farray[i]);
        }
        for (i = 0; i < 66; i++)
          mpz_clear(farray[i]);
      }
    }
    mpz_clear(n);

void
_GMP_factor(IN char* strn)
  PREINIT:
    mpz_t n;
  PPCODE:
    VALIDATE_AND_SET("factor", n, strn);
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
        while ( mpz_cmp_ui(n, TRIAL_LIM*TRIAL_LIM) > 0 && !_GMP_is_prob_prime(n) ) {
          int success = 0;
          int o = get_verbose_level();
          UV B1 = 5000;

          /*
           * This set of operations is meant to provide good performance for
           * "random" numbers as input.  Hence we stack lots of effort up front
           * looking for small factors: prho and pbrent are ~ O(f^1/2) where
           * f is the smallest factor.  SQUFOF is O(N^1/4), so arguable not
           * any better.  p-1 and ECM are quite useful for pulling out small
           * factors (6-20 digits).
           *
           * Factoring a 778-digit number consisting of 101 8-digit factors
           * should complete in under 3 seconds.  Factoring numbers consisting
           * of many 12-digit or 14-digit primes should take under 10 seconds.
           */

          if (mpz_cmp_ui(n, (unsigned long)(UV_MAX>>2)) < 0) {
            UV ui_n = mpz_get_ui(n);
            UV ui_factors[2];
            if (!mpz_cmp_ui(n, ui_n)) {
              success = racing_squfof_factor(ui_n, ui_factors, 200000)-1;
              if (success) {
                mpz_set_ui(f, ui_factors[0]);
              } else {
                if (o > 2) {gmp_printf("UV SQUFOF failed %Zd\n", n);}
              }
            }
            if (success&&o) {gmp_printf("UV SQUFOF found factor %Zd\n", f);o=0;}
          }

          /* Make sure it isn't a perfect power */
          if (!success)  success = _GMP_power_factor(n, f);
          if (success&&o) {gmp_printf("perfect power found factor %Zd\n", f);o=0;}

          if (!success)  success = _GMP_pminus1_factor(n, f, 10000, 200000);
          if (success&&o) {gmp_printf("p-1 (10k) found factor %Zd\n", f);o=0;}

          /* Really small ECM to find small factors */
          if (!success)  success = _GMP_ECM_FACTOR(n, f, 150, 50);
          if (success&&o) {gmp_printf("tiny ecm (150) found factor %Zd\n", f);o=0;}
          if (!success)  success = _GMP_ECM_FACTOR(n, f, 500, 30);
          if (success&&o) {gmp_printf("tiny ecm (500) found factor %Zd\n", f);o=0;}
          if (!success)  success = _GMP_ECM_FACTOR(n, f, 2000, 10);
          if (success&&o) {gmp_printf("tiny ecm (2000) found factor %Zd\n", f);o=0;}

          /* Small p-1 */
          if (!success)  success = _GMP_pminus1_factor(n, f, 200000, 4000000);
          if (success&&o) {gmp_printf("p-1 (200k) found factor %Zd\n", f);o=0;}

          /* ECM with a good chance of success */
          if (!success) {
            UV bits = mpz_sizeinbase(n, 2);
            UV curves = 20;
            if      (bits < 100)  B1 =   5000;     /* All need tuning */
            else if (bits < 128)  B1 =  10000;
            else if (bits < 160)  B1 =  20000;
            else if (bits < 192)  B1 =  30000;
            else if (bits < 224){ B1 =  40000; curves =  40; }
            else if (bits < 256){ B1 =  80000; curves =  40; }
            else if (bits < 512){ B1 = 160000; curves =  80; }
            else                { B1 = 320000; curves = 160; }
            success = _GMP_ECM_FACTOR(n, f, B1, curves);
            if (success&&o) {gmp_printf("small ecm (%luk,%lu) found factor %Zd\n", B1/1000, curves, f);o=0;}
          }

          /* QS (30+ digits).  Fantastic if it is a semiprime, but can be
           * slow and a memory hog if not (compared to ECM).  Restrict to
           * reasonable size numbers (< 91 digits).  Because of the way it
           * works, it will generate (possibly) multiple factors for the same
           * amount of work.  Go to some trouble to use them. */
          if (!success && mpz_sizeinbase(n,10)>=30 && mpz_sizeinbase(n,2)<300) {
            mpz_t farray[66];
            int i, nfactors;
            for (i = 0; i < 66; i++)
              mpz_init(farray[i]);
            nfactors = _GMP_simpqs(n, farray);
            mpz_set(f, farray[0]);
            if (nfactors > 2) {
              /* We found multiple factors */
              for (i = 2; i < nfactors; i++) {
                if (o){gmp_printf("SIMPQS found extra factor %Zd\n",farray[i]);}
                if (ntofac == MAX_FACTORS-1) croak("Too many factors\n");
                mpz_init_set(tofac_stack[ntofac], farray[i]);
                ntofac++;
                mpz_divexact(n, n, farray[i]);
              }
              /* f = farray[0], n = farray[1], farray[2..] pushed */
            }
            for (i = 0; i < 66; i++)
              mpz_clear(farray[i]);
            success = nfactors > 1;
            if (success&&o) {gmp_printf("SIMPQS found factor %Zd\n", f);o=0;}
          }

          if (!success)  success = _GMP_ECM_FACTOR(n, f, 2*B1, 20);
          if (success&&o) {gmp_printf("ecm (%luk,20) found factor %Zd\n",2*B1/1000,f);o=0;}

          if (!success)  success = _GMP_pbrent_factor(n, f, 1, 1*1024*1024);
          if (success&&o) {gmp_printf("pbrent (1,1M) found factor %Zd\n", f);o=0;}

          if (!success)  success = _GMP_ECM_FACTOR(n, f, 4*B1, 20);
          if (success&&o) {gmp_printf("ecm (%luk,20) ecm found factor %Zd\n", 4*B1,f);o=0;}

          if (!success)  success = _GMP_ECM_FACTOR(n, f, 8*B1, 20);
          if (success&&o) {gmp_printf("ecm (%luk,20) ecm found factor %Zd\n", 8*B1,f);o=0;}

          /* HOLF in case it's a near-ratio-of-perfect-square */
          if (!success)  success = _GMP_holf_factor(n, f, 1*1024*1024);
          if (success&&o) {gmp_printf("holf found factor %Zd\n", f);o=0;}

          /* Large p-1 with stage 2: B2 = 20*B1 */
          if (!success)  success = _GMP_pminus1_factor(n,f,5000000,5000000*20);
          if (success&&o) {gmp_printf("p-1 (5M) found factor %Zd\n", f);o=0;}

          if (!success)  success = _GMP_ECM_FACTOR(n, f, 32*B1, 40);
          if (success&&o) {gmp_printf("ecm (%luk,40) ecm found factor %Zd\n", 32*B1,f);o=0;}

          /*
          if (!success)  success = _GMP_pbrent_factor(n, f, 2, 512*1024*1024);
          if (success&&o) {gmp_printf("pbrent (2,512M) found factor %Zd\n", f);o=0;}

          if (!success)  success = _GMP_squfof_factor(n, f, 256*1024*1024);
          if (success&&o) {gmp_printf("squfof found factor %Zd\n", f);o=0;}
          */

          /* Our method of last resort: ECM with high bmax and many curves*/
          if (!success) {
            UV i;
            if (get_verbose_level()) gmp_printf("starting large ECM on %Zd\n",n);
            B1 *= 8;
            for (i = 0; i < 10; i++) {
              success = _GMP_ECM_FACTOR(n, f, B1, 100);
              if (success) break;
              B1 *= 2;
            }
            if (success&&o) {gmp_printf("ecm (%luk,100) ecm found factor %Zd\n", B1,f);o=0;}
          }

          if (success) {
            if (!mpz_divisible_p(n, f) || !mpz_cmp_ui(f, 1) || !mpz_cmp(f, n)) {
              gmp_printf("n = %Zd  f = %Zd\n", n, f);
              croak("Incorrect factoring");
            }
          }

          if (!success) {
            /* TODO: What to do with composites we can't factor?
             *       Push them as "C#####" ?
             *       For now, just push them as if we factored.
             */
            mpz_set(f, n);
            XPUSH_MPZ(f);
          } else if (_GMP_is_prob_prime(f)) {
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
