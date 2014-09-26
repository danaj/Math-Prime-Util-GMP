
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
#include "factor.h"
#define _GMP_ECM_FACTOR(n, f, b1, ncurves) \
   _GMP_ecm_factor_projective(n, f, b1, 0, ncurves)


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
_GMP_miller_rabin(IN char* strn, ...)
  ALIAS:
    is_pseudoprime = 1
  PREINIT:
    mpz_t n, a;
    char* strbase;
  CODE:
    validate_string_number("GMP_miller_rabin (n)", strn);
    strbase = (items == 1) ? "2" : SvPV_nolen(ST(1));  /* default base = 2 */
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
    if (ix == 0) {
      RETVAL = _GMP_miller_rabin(n, a);
    } else {
      mpz_t nm1; mpz_init(nm1); mpz_sub_ui(nm1, n, 1);
      mpz_powm(a, a, nm1, n);
      mpz_clear(nm1);
      RETVAL = !mpz_cmp_ui(a,1);
    }
    mpz_clear(a);
    mpz_clear(n);
  OUTPUT:
    RETVAL

int miller_rabin_random(IN char* strn, IN UV nbases, IN char* seedstr = 0)
  PREINIT:
    mpz_t n;
  CODE:
    VALIDATE_AND_SET("miller_rabin_random", n, strn);
    RETVAL = _GMP_miller_rabin_random(n, nbases, seedstr);
    mpz_clear(n);
  OUTPUT:
    RETVAL

#define PRIMALITY_START(name, small_retval, test_small_factors) \
    /* Negative numbers return 0 */ \
    if ((strn != 0) && (strn[0] == '-') ) \
      XSRETURN_IV(0); \
    validate_string_number(name " (n)", strn); \
    /* Handle single digit numbers */ \
    if (strn[1] == 0) { \
      int q_is_prime = 0; \
      switch (strn[0]) { \
        case '2': case '3': case '5': case '7': q_is_prime = small_retval; \
                                                break; \
      } \
      XSRETURN_IV(q_is_prime); \
    } \
    /* Test for small multiples while it is still a string */ \
    if (test_small_factors) { \
      UV digsum = 0; \
      int i, slen = strlen(strn); \
      /* Multiples of 2 and 5 return 0 */ \
      switch (strn[slen-1]-'0') { \
        case 0: case 2: case 4: case 5: case 6: case 8: \
           XSRETURN_IV(0); break; \
      } \
      /* Multiples of 3 return 0 */ \
      for (i = 0; i < slen; i++)  digsum += strn[i]-'0'; \
      if (digsum % 3 == 0)  XSRETURN_IV(0); \
    } \
    mpz_init_set_str(n, strn, 10);

int
is_lucas_pseudoprime(IN char* strn)
  ALIAS:
    is_strong_lucas_pseudoprime = 1
    is_extra_strong_lucas_pseudoprime = 2
    is_frobenius_underwood_pseudoprime = 3
    is_perrin_pseudoprime = 4
  PREINIT:
    mpz_t n;
  CODE:
    if ((strn != 0) && (strn[0] == '-') )
      croak("Parameter '%s' must be a positive integer\n", strn);
    PRIMALITY_START("is_lucas_pseudoprime", 1, 0);
    switch (ix) {
      case 0: RETVAL = _GMP_is_lucas_pseudoprime(n, 0); break;
      case 1: RETVAL = _GMP_is_lucas_pseudoprime(n, 1); break;
      case 2: RETVAL = _GMP_is_lucas_pseudoprime(n, 2); break;
      case 3: RETVAL = _GMP_is_frobenius_underwood_pseudoprime(n); break;
      case 4:
      default:RETVAL = is_perrin_pseudoprime(n); break;
    }
    mpz_clear(n);
  OUTPUT:
    RETVAL

int
is_almost_extra_strong_lucas_pseudoprime(IN char* strn, IN UV increment = 1)
  PREINIT:
    mpz_t n;
  CODE:
    if ((strn != 0) && (strn[0] == '-') )
      croak("Parameter '%s' must be a positive integer\n", strn);
    if (increment == 0 || increment > 65535)
      croak("Increment parameter must be >0 and < 65536");
    PRIMALITY_START("is_almost_extra_strong_lucas_pseudoprime", 1, 0);
    RETVAL = _GMP_is_almost_extra_strong_lucas_pseudoprime(n, increment);
    mpz_clear(n);
  OUTPUT:
    RETVAL

int
is_frobenius_pseudoprime(IN char* strn, IN IV P = 0, IN IV Q = 0)
  PREINIT:
    mpz_t n;
  CODE:
    if ((strn != 0) && (strn[0] == '-') )
      croak("Parameter '%s' must be a positive integer\n", strn);
    PRIMALITY_START("is_frobenius_pseudoprime", 1, 0);
    RETVAL = is_frobenius_pseudoprime(n, P, Q);
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
    is_bpsw_prime = 5
  PREINIT:
    mpz_t n;
    int ret;
  CODE:
    /* Returns arg for single-dig primes, 0 for multiples of 2, 3, 5, or neg */
    PRIMALITY_START("is_prime", 2, 1);
    switch (ix) {
      case 0: ret = _GMP_is_prime(n); break;
      case 1: ret = _GMP_is_prob_prime(n); break;
      case 2: ret = _GMP_is_aks_prime(n); break;
      case 3: ret = _GMP_primality_bls_nm1(n, 100, 0); break;
      case 4: ret = _GMP_ecpp(n, 0); break;
      case 5:
      default:ret = _GMP_BPSW(n); break;
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
    PRIMALITY_START("is_provable_prime", 2, 1);
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

UV
is_power(IN char* strn, IN UV a = 0)
  PREINIT:
    mpz_t n;
  CODE:
    validate_string_number("is_power (n)", strn);
    mpz_init_set_str(n, strn, 10);
    RETVAL = is_power(n, a);
    mpz_clear(n);
  OUTPUT:
    RETVAL


#define XPUSH_MPZ(n) \
  do { \
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
  } while (0)

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
primorial(IN char* strn)
  ALIAS:
    pn_primorial = 1
    consecutive_integer_lcm = 2
    exp_mangoldt = 3
    totient = 4
    carmichael_lambda = 5
    factorial = 6
    bernfrac = 7
    znprimroot = 8
  PREINIT:
    mpz_t res, n;
    UV un;
  PPCODE:
    if (strn != 0 && strn[0] == '-') { /* If input is negative... */
      if (ix == 3)  XSRETURN_IV(1);    /* exp_mangoldt return 1 */
      if (ix == 8)  strn++;            /* znprimroot flip sign */
    }
    VALIDATE_AND_SET("primorial", n, strn);
    un = mpz_get_ui(n);
    mpz_init(res);
    switch (ix) {
      case 0:  _GMP_primorial(res, n);  break;
      case 1:  _GMP_pn_primorial(res, un);  break;
      case 2:  _GMP_lcm_of_consecutive_integers(un, res);  break;
      case 3:  exp_mangoldt(res, n);  break;
      case 4:  totient(res, n);  break;
      case 5:  carmichael_lambda(res, n);  break;
      case 6:  mpz_fac_ui(res, un);  break;  /* Uses swing, so let's use it */
      case 7:  bernfrac(n, res, n);
               XPUSH_MPZ(n);
               break;
      case 8:
      default: znprimroot(res, n);  break;
    }
    if (ix == 8 && !mpz_sgn(res) && mpz_cmp_ui(n,1) != 0)
      {  mpz_clear(n);  mpz_clear(res);  XSRETURN_UNDEF;  }
    XPUSH_MPZ(res);
    mpz_clear(n);
    mpz_clear(res);


void
gcd(...)
  PROTOTYPE: @
  ALIAS:
    lcm = 1
    vecsum = 2
    vecprod = 3
  PREINIT:
    int i;
    mpz_t ret, n;
  PPCODE:
    if (items == 0) XSRETURN_IV( (ix == 3) ? 1 : 0);
    if (ix == 3) {
      mpz_t* list;
      New(0, list, items, mpz_t);
      for (i = 0; i < items; i++) {
        char* strn = SvPV_nolen(ST(i));
        validate_string_number("vecprod", (strn[0]=='-') ? strn+1 : strn);
        mpz_init_set_str(list[i], strn, 10);
      }
      mpz_product(list, 0, items-1);
      XPUSH_MPZ(list[0]);
      for (i = 0; i < items; i++)  mpz_clear(list[i]);
      Safefree(list);
      XSRETURN(1);
    }
    mpz_init(n);
    mpz_init_set_ui(ret, (ix == 1 || ix == 3) ? 1 : 0);
    for (i = 0; i < items; i++) {
      char* strn = SvPV_nolen(ST(i));
      validate_string_number("gcd/lcm", (strn[0]=='-') ? strn+1 : strn);
      if (ix <= 1 && strn != 0 && strn[0] == '-') strn++;
      mpz_set_str(n, strn, 10);
      switch (ix) {
        case 0:  mpz_gcd(ret, ret, n); break;
        case 1:  mpz_lcm(ret, ret, n); break;
        case 2:  mpz_add(ret, ret, n); break;
        case 3:
        default: mpz_mul(ret, ret, n); break;
      }
    }
    XPUSH_MPZ(ret);
    mpz_clear(n);
    mpz_clear(ret);

int
kronecker(IN char* stra, IN char* strb)
  ALIAS:
    valuation = 1
  PREINIT:
    mpz_t a, b;
  CODE:
    validate_string_number("kronecker", (stra[0]=='-') ? stra+1 : stra);
    validate_string_number("kronecker", (strb[0]=='-') ? strb+1 : strb);
    mpz_init_set_str(a, stra, 10);
    mpz_init_set_str(b, strb, 10);
    if (ix == 0) {
      RETVAL = mpz_kronecker(a, b);
    } else {
      mpz_abs(a,a);
      mpz_abs(b,b);
      if (mpz_cmp_ui(a,1) <= 0 || mpz_cmp_ui(b,1) <= 0)
        RETVAL = 0;
      else
        RETVAL = mpz_remove(a, a, b);
    }
    mpz_clear(b);
    mpz_clear(a);
  OUTPUT:
    RETVAL

void
moebius(IN char* strn, IN char* stro = 0)
  PREINIT:
    mpz_t n;
  PPCODE:
    VALIDATE_AND_SET("moebius", n, strn);
    if (stro == 0) {
      int result = moebius(n);
      mpz_clear(n);
      XSRETURN_IV(result);
    } else {   /* Ranged result */
      mpz_t nhi;
      VALIDATE_AND_SET("moebius high value", nhi, stro);
      while (mpz_cmp(n, nhi) <= 0) {
        XPUSHs(sv_2mortal(newSViv( moebius(n) )));
        mpz_add_ui(n, n, 1);
      }
      mpz_clear(n);
      mpz_clear(nhi);
    }

int
liouville(IN char* strn)
  PREINIT:
    mpz_t n;
  CODE:
    VALIDATE_AND_SET("liouville", n, strn);
    RETVAL = liouville(n);
    mpz_clear(n);
  OUTPUT:
    RETVAL

void
invmod(IN char* stra, IN char* strb)
  ALIAS:
    binomial = 1
    gcdext = 2
    jordan_totient = 3
    znorder = 4
  PREINIT:
    mpz_t a, b, t;
    int retundef;
  PPCODE:
    validate_string_number("invmod", (stra[0]=='-') ? stra+1 : stra);
    validate_string_number("invmod", (strb[0]=='-') ? strb+1 : strb);
    mpz_init_set_str(a, stra, 10);
    mpz_init_set_str(b, strb, 10);
    retundef = 0;
    if (ix == 0) {
      /* undef if a|b is zero, 0 if b is 1, otherwise result of mpz_invert */
      if (!mpz_sgn(b) || !mpz_sgn(a))  retundef = 1;
      else if (!mpz_cmp_ui(b,1))       mpz_set_ui(a,0);
      else                             retundef = !mpz_invert(a,a,b);
    } else if (ix == 1) {
      if (mpz_sgn(b) < 0) {   /* Handle negative k */
        if (mpz_sgn(a) >= 0 || mpz_cmp(b,a) > 0)  mpz_set_ui(a, 0);
        else                                      mpz_sub(b, a, b);
      }
      mpz_bin_ui(a, a, mpz_get_ui(b));
    } else if (ix == 2) {
      mpz_init(t);
      mpz_gcdext(a, t, b, a, b);
      XPUSH_MPZ(t);  XPUSH_MPZ(b);
      mpz_clear(t);
    } else if (ix == 3) {
      jordan_totient(a, b, mpz_get_ui(a));
    } else {
      znorder(a, a, b);
      if (!mpz_sgn(a)) retundef = 1;
    }
    if (!retundef) XPUSH_MPZ(a);
    mpz_clear(b); mpz_clear(a);
    if (retundef) XSRETURN_UNDEF;

void partitions(IN UV n)
  ALIAS:
    Pi = 1
  PREINIT:
    UV i, j, k;
  PPCODE:
    if (ix ==1) {
      if (n == 1)
        XSRETURN_IV(3);
      else if (n > 0) {
        char* pi = pidigits(n);
        XPUSHs(sv_2mortal(newSVpvn(pi, n+1)));
        Safefree(pi);
      }
    } else if (n == 0) {
      XPUSHs(sv_2mortal(newSVuv( 1 )));
    } else if (n <= 3) {
      XPUSHs(sv_2mortal(newSVuv( n )));
    } else {
      mpz_t psum;
      mpz_t* part;
      UV* pent;
      UV d = (UV) sqrt(n+1);

      New(0, pent, 2*d+2, UV);
      pent[0] = 0;
      pent[1] = 1;
      for (i = 1; i <= d; i++) {
        pent[2*i  ] = ( i   *(3*i+1)) / 2;
        pent[2*i+1] = ((i+1)*(3*i+2)) / 2;
      }
      New(0, part, n+1, mpz_t);
      mpz_init_set_ui(part[0], 1);
      mpz_init(psum);
      for (j = 1; j <= n; j++) {
        mpz_set_ui(psum, 0);
        for (k = 1; pent[k] <= j; k++) {
          if ((k+1) & 2) mpz_add(psum, psum, part[ j - pent[k] ]);
          else           mpz_sub(psum, psum, part[ j - pent[k] ]);
        }
        mpz_init_set(part[j], psum);
      }
      mpz_clear(psum);
      XPUSH_MPZ( part[n] );
      for (i = 0; i <= n; i++)
        mpz_clear(part[i]);
      Safefree(part);
      Safefree(pent);
    }

void
stirling(IN UV n, IN UV m, IN UV type = 1)
  PREINIT:
    mpz_t r;
  PPCODE:
    mpz_init(r);
    stirling(r, n, m, type);
    XPUSH_MPZ( r );
    mpz_clear(r);

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

#define TSTAVAL(arr, val)   (arr[(val) >> 6] & (1U << (((val)>>1) & 0x1F)))

void
sieve_primes(IN char* strlow, IN char* strhigh, IN UV k = 0)
  PREINIT:
    mpz_t low, high, t;
    UV i, length;
    int test_primality;
    uint32_t* comp;
  PPCODE:
    VALIDATE_AND_SET("sieve_primes", low, strlow);
    VALIDATE_AND_SET("sieve_primes", high, strhigh);
    test_primality = 0;
    if (k < 2) {
      test_primality = 1;
      k = 5000 * mpz_sizeinbase(high,2);
    }

    if (mpz_cmp_ui(low, k) < 0)    croak("TODO: small sieves");
    if (mpz_even_p(low))           mpz_add_ui(low, low, 1);
    if (mpz_even_p(high))          mpz_sub_ui(high, high, 1);

    if (mpz_cmp(low, high) <= 0) {
      mpz_init(t);

      mpz_sqrt(t, high);           /* No need for k to be > sqrt(high) */
      if (mpz_cmp_ui(t, k) < 0)
        k = mpz_get_ui(t);

      mpz_sub(t, high, low);
      length = mpz_get_ui(t) + 1;

      /* Get bit array of odds marked with composites(k) marked with 1 */
      comp = partial_sieve(low, length, k);
      /* Convert to corresponding mpz and send to output */
      for (i = 1; i <= length; i += 2) {
        if (!TSTAVAL(comp, i)) {
          mpz_add_ui(t, low, i);
          if (!test_primality || _GMP_BPSW(t))
            XPUSH_MPZ( t );
        }
      }
      mpz_clear(t);
      Safefree(comp);
    }
    mpz_clear(low);
    mpz_clear(high);

void
lucas_sequence(IN char* strn, IN IV P, IN IV Q, IN char* strk)
  PREINIT:
    mpz_t U, V, Qk, n, k, t;
  PPCODE:
    VALIDATE_AND_SET("lucas_sequence", n, strn);
    VALIDATE_AND_SET("lucas_sequence", k, strk);
    mpz_init(U);  mpz_init(V);  mpz_init(Qk);  mpz_init(t);

    _GMP_lucas_seq(U, V, n, P, Q, k, Qk, t);
    XPUSH_MPZ(U);
    XPUSH_MPZ(V);
    XPUSH_MPZ(Qk);

    mpz_clear(n);  mpz_clear(k);
    mpz_clear(U);  mpz_clear(V);  mpz_clear(Qk);  mpz_clear(t);


#define SET_UV_VIA_MPZ_STRING(uva, sva, name) \
  { \
      mpz_t t; \
      char* stra = SvPV_nolen(sva); \
      validate_string_number(name, stra); \
      mpz_init_set_str(t, stra, 10); \
      uva = mpz_get_ui(t); \
      mpz_clear(t); \
  }

void
trial_factor(IN char* strn, ...)
  ALIAS:
    prho_factor = 1
    pbrent_factor = 2
    pminus1_factor = 3
    pplus1_factor = 4
    holf_factor = 5
    squfof_factor = 6
    ecm_factor = 7
    qs_factor = 8
  PREINIT:
    mpz_t n;
    UV arg1, arg2, uf;
    static const UV default_arg1[] =
       {0,    64000000,64000000,5000000,5000000,256000000,16000000,0,  0  };
     /* Trial,Rho,     Brent,   P-1,    P+1,    HOLF,     SQUFOF,  ECM,QS */
  PPCODE:
    VALIDATE_AND_SET(" specific factor", n, strn);
    {
      int cmpr = mpz_cmp_ui(n,1);
      if (cmpr <= 0) {
        mpz_clear(n);
        XSRETURN_IV( (cmpr < 0)  ?  0  :  1 );
      }
    }
    arg1 = default_arg1[ix];
    arg2 = 0;
    if (items >= 2) SET_UV_VIA_MPZ_STRING(arg1, ST(1), "specific factor arg 1");
    if (items >= 3) SET_UV_VIA_MPZ_STRING(arg2, ST(2), "specific factor arg 2");
    while (mpz_even_p(n)) {
      XPUSHs(sv_2mortal(newSVuv(2)));
      mpz_divexact_ui(n, n, 2);
    }
    while (mpz_divisible_ui_p(n, 3)) {
      XPUSHs(sv_2mortal(newSVuv(3)));
      mpz_divexact_ui(n, n, 3);
    }
    while (mpz_divisible_ui_p(n, 5)) {
      XPUSHs(sv_2mortal(newSVuv(5)));
      mpz_divexact_ui(n, n, 5);
    }
    if (mpz_cmp_ui(n,1) > 0 && !_GMP_is_prob_prime(n)) {
      mpz_t f;
      int success = 0;

      mpz_init(f);
      switch (ix) {
        case 0: if (arg1 == 0) arg1 = 2147483647;
                uf = _GMP_trial_factor(n, 2, arg1);
                mpz_set_ui(f, uf);
                success = (uf > 0);
                break;
        case 1: success = _GMP_prho_factor(n, f, 3, arg1);        break;
        case 2: success = _GMP_pbrent_factor(n, f, 3, arg1);      break;
        case 3: if (arg2 == 0)  arg2 = arg1*10;
                success = _GMP_pminus1_factor(n, f, arg1,arg2);   break;
        case 4: if (arg2 == 0)  arg2 = arg1*10;
                success = _GMP_pplus1_factor(n, f, 0,arg1,arg2);  break;
        case 5: success = _GMP_holf_factor(n, f, arg1);           break;
        case 6: success = _GMP_squfof_factor(n, f, arg1);         break;
        case 7: if (arg2 == 0) arg2 = 100;
                if (arg1 == 0) {
                  success =    _GMP_ECM_FACTOR(n, f,     1000, 40)
                            || _GMP_ECM_FACTOR(n, f,    10000, 40)
                            || _GMP_ECM_FACTOR(n, f,   100000, 40)
                            || _GMP_ECM_FACTOR(n, f,  1000000, 40)
                            || _GMP_ECM_FACTOR(n, f, 10000000,100);
                } else {
                  success = _GMP_ECM_FACTOR(n, f, arg1, arg2);
                }
                break;
        case 8:
        default:{
                  mpz_t farray[66];
                  int i, nfactors;
                  for (i = 0; i < 66; i++)
                    mpz_init(farray[i]);
                  nfactors = _GMP_simpqs(n, farray);
                  for (i = 0; i < nfactors; i++)
                    XPUSH_MPZ(farray[i]);
                  for (i = 0; i < 66; i++)
                    mpz_clear(farray[i]);
                  if (nfactors > 0)
                    mpz_set_ui(n, 1);
                }
                break;
      }

      if (success) {
        XPUSH_MPZ(f);
        mpz_divexact(n, n, f);
      }
      mpz_clear(f);
    }
    if (mpz_cmp_ui(n,1) > 0)
      XPUSH_MPZ(n);
    mpz_clear(n);

void
_GMP_factor(IN char* strn)
  PREINIT:
    mpz_t n;
    mpz_t* factors;
    int* exponents;
    int nfactors, i, j;
  PPCODE:
    VALIDATE_AND_SET("factor", n, strn);
    nfactors = factor(n, &factors, &exponents);
    for (i = 0; i < nfactors; i++) {
      for (j = 0; j < exponents[i]; j++) {
        XPUSH_MPZ(factors[i]);
      }
    }
    clear_factors(nfactors, &factors, &exponents);
    mpz_clear(n);
