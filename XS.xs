#define PERL_NO_GET_CONTEXT 1 /* Define at top for more efficiency. */

#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
/* We're not using anything for which we need ppport.h */

#include <string.h>
#include <gmp.h>

#include "ptypes.h"
#include "gmp_main.h"
#include "factmod.h"
#include "primality.h"
#include "lucas_seq.h"
#include "squfof126.h"
#include "ecm.h"
#include "simpqs.h"
#include "bls75.h"
#include "ecpp.h"
#include "aks.h"
#include "rootmod.h"
#include "legendre_phi.h"
#include "znlog.h"
#include "utility.h"
#include "factor.h"
#include "isaac.h"
#include "random_prime.h"
#include "perfect_powers.h"
#include "powerfree.h"
#include "real.h"
#define _GMP_ECM_FACTOR(n, f, b1, ncurves) \
   _GMP_ecm_factor_projective(n, f, b1, 0, ncurves)


/* Instead of trying to suck in lots of Math::BigInt::GMP and be terribly
 * clever (and brittle), just do all C<->Perl bigints via strings.  It's
 * crude but seems to work pretty well.
 */

static const char* _subname(pTHX_ const CV *cv) { return GvNAME(CvGV(cv)); }
#define SUBNAME _subname(aTHX_ cv)

#define IFLAG_ANY    0x00
#define IFLAG_NONNEG 0x01
#define IFLAG_POS    0x02
#define IFLAG_ABS    0x04

static void _croak_bad_integer(pTHX_ CV* cv, const char* vname, const char* s, int iflag)
{
  const char* m = (iflag & IFLAG_POS)    ? "a positive integer"
                : (iflag & IFLAG_NONNEG) ? "a non-negative integer"
                                         : "an integer";
  croak("%s (%s): input '%s' must be %s", SUBNAME, vname, s, m);
}

static int _validate_integer_string(pTHX_ CV* cv, const char* vname,
                                    const char* s, int iflag)
{
  const char *p;
  int hasneg, nonzero;

  if (s == 0)  croak("%s (%s): null string pointer as input", SUBNAME, vname);
  if (*s == 0) croak("%s (%s): empty string as input", SUBNAME, vname);

  hasneg = (*s == '-');
  nonzero = 0;
  p = s + (*s == '-' || *s == '+');
  if (*p == 0)  /* Check for no digits, e.g. only a sign */
    _croak_bad_integer(aTHX_ cv, vname, s, iflag);

  for (; *p != 0; p++) {
    if (!isdigit((unsigned char)*p))
      _croak_bad_integer(aTHX_ cv, vname, s, iflag);
    nonzero |= (*p != '0');
  }

  if ((iflag & IFLAG_POS) && (hasneg || !nonzero))
    _croak_bad_integer(aTHX_ cv, vname, s, iflag);
  if ((iflag & IFLAG_NONNEG) && hasneg && nonzero)
    _croak_bad_integer(aTHX_ cv, vname, s, iflag);

  return hasneg && nonzero;
}
#define validate_integer_string(vname,s,iflag) \
  _validate_integer_string(aTHX_ cv, vname, s, iflag)

static NOINLINE int _set_integer_string(pTHX_ CV* cv, mpz_t v, const char* vname, const char* s, int iflag) {
  int wasneg;
  wasneg = validate_integer_string(vname, s, iflag);
  /* After validate_integer_string, s is non-zero and *s =~ /^[+-]?\d+$/ */
  mpz_init_set_str(v, s + (*s == '+'), 10);
  if ((iflag & IFLAG_ABS) && wasneg)
    mpz_abs(v, v);

  if ((iflag & IFLAG_POS) && mpz_sgn(v) <= 0) {
    mpz_clear(v);
    _croak_bad_integer(aTHX_ cv, vname, s, iflag);
  }

  return wasneg;
}
#define set_integer_string(v, vname, s, iflag) \
  _set_integer_string(aTHX_ cv, v, vname, s, iflag)

#define validate_and_set(var, iflag) \
  set_integer_string(var, #var, str ## var, iflag)

static SV* sv_return_for_mpz(pTHX_ const mpz_t n) {
  if (mpz_fits_uv_p(n)) {
    return newSVuv(mpz_get_uv(n));
  } else if (mpz_fits_iv_p(n)) {
    return newSViv(mpz_get_iv(n));
  } else {
    SV* sv;
    char* str;
    int nsize = mpz_sizeinbase(n, 10) + 2;
    New(0, str, nsize, char);
    mpz_get_str(str, 10, n);
    sv = newSVpv(str, 0);
    Safefree(str);
    return sv;
  }
}

#define XPUSH_MPZ(n) \
  XPUSHs(sv_2mortal( sv_return_for_mpz(aTHX_ n) ))
#define XPUSH_INT(n) \
  XPUSHs(sv_2mortal( newSViv(n) ))
#define XPUSH_UINT(n) \
  XPUSHs(sv_2mortal( newSVuv(n) ))

static int _sv_is_math_object(pTHX_ SV* sv)
{
  const char *stashname;

  if (!SvROK(sv) || !sv_isobject(sv))
    return 0;
  stashname = HvNAME(SvSTASH(SvRV(sv)));
  return stashname != 0 && strnEQ(stashname, "Math::", 6);
}

static void _croak_invalid_fromdigits_digit(pTHX_ const mpz_t base)
{
  SV *basesv = sv_2mortal(sv_return_for_mpz(aTHX_ base));
  croak("fromdigits: invalid digit for base %s", SvPV_nolen(basesv));
}

static char* cert_with_header(char* proof, mpz_t n) {
  char *str, *strptr;
  if (proof == 0 && mpz_sizeinbase(n,2) <= 64 && _GMP_is_prime(n)) {
    New(0, str, 50, char);
    gmp_sprintf(str, "Type Small\nN %Zd\n", n);
    return cert_with_header(str, n);
  } else if (proof == 0) {
    New(0, str, 1, char);
    str[0] = '\0';
  } else {
    New(0, str, strlen(proof) + 100 + mpz_sizeinbase(n,10), char);
    strptr = str;
    strptr += gmp_sprintf(strptr, "[MPU - Primality Certificate]\nVersion 1.0\n\nProof for:\nN %Zd\n\n", n);
    strcat(strptr, proof);
    Safefree(proof);
  }
  return str;
}


MODULE = Math::Prime::Util::GMP		PACKAGE = Math::Prime::Util::GMP

PROTOTYPES: ENABLE

void _GMP_init()

void _GMP_destroy()

void _GMP_memfree()

void _GMP_set_verbose(IN int v)
  PPCODE:
     set_verbose_level(v);

void seed_csprng(IN UV bytes, IN unsigned char* seed)
  PPCODE:
    isaac_init(bytes, seed);

UV irand()
  ALIAS:
    irand64 = 1
    is_csprng_well_seeded = 2
  CODE:
    switch (ix) {
#if BITS_PER_WORD >= 64
      case 0:  RETVAL = isaac_rand32(); break;
      case 1:  RETVAL = (((UV)isaac_rand32()) << 32) | isaac_rand32();  break;
#else
      case 0:
      case 1:  RETVAL = isaac_rand32(); break;
#endif
      case 2:
      default: RETVAL = isaac_seeded(); break;
    }
  OUTPUT:
    RETVAL

NV drand(NV m = 1.0)
  CODE:
    RETVAL = m * drand64();
  OUTPUT:
    RETVAL

int
is_pseudoprime(IN char* strn, ...)
  ALIAS:
    is_euler_pseudoprime = 1
    is_strong_pseudoprime = 2
  PREINIT:
    int i, nbases;
    mpz_t n, *A;
  CODE:
    validate_and_set(n, IFLAG_ANY);
    if (items == 1) {
      nbases = 1;
      New(0, A, nbases, mpz_t);
      mpz_init_set_ui(A[0], 2);
    } else {
      nbases = items-1;
      New(0, A, nbases, mpz_t);
      for (i = 1; i < items; i++) {
        set_integer_string(A[i-1], "base", SvPV_nolen(ST(i)), IFLAG_NONNEG);
        if (mpz_cmp_ui(A[i-1], 2) < 0)
          croak("Base %s is invalid", SvPV_nolen(ST(i)));
      }
    }
    RETVAL = (mpz_cmp_ui(n,2) >= 0);
    for (i = 0; RETVAL && i < nbases; i++) {
      switch (ix) {
        case 0:  RETVAL = is_pseudoprime(n, A[i]); break;
        case 1:  RETVAL = is_euler_pseudoprime(n, A[i]); break;
        case 2:
        default: RETVAL = miller_rabin(n, A[i]); break;
      }
    }
    for (i = 0; i < nbases; i++)
      mpz_clear(A[i]);
    Safefree(A);
    mpz_clear(n);
  OUTPUT:
    RETVAL

int miller_rabin_random(IN char* strn, IN IV nbases = 1, IN char* seedstr = 0)
  PREINIT:
    mpz_t n;
  CODE:
    if (nbases < 0)
      croak("Parameter '%"IVdf"' must be a positive integer\n", nbases);
    validate_and_set(n, IFLAG_ANY);
    RETVAL = miller_rabin_random(n, nbases, seedstr);
    mpz_clear(n);
  OUTPUT:
    RETVAL

#define PRIMALITY_START(name, small_retval, test_small_factors) \
    /* Negative numbers return 0 */ \
    if (validate_integer_string("n", strn, IFLAG_ANY)) \
      XSRETURN_IV(0); \
    if (*strn == '+') strn++; \
    /* Handle single digit numbers */ \
    if (strn[1] == 0) { \
      uint8_t v = strn[0] - '0'; \
      XSRETURN_IV((v == 2 || v == 3 || v == 5 || v == 7) ? (small_retval) : 0);\
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
    is_frobenius_khashin_pseudoprime = 4
    is_euler_plumb_pseudoprime = 5
  PREINIT:
    mpz_t n;
  CODE:
    PRIMALITY_START("is_lucas_pseudoprime", 1, 0);
    switch (ix) {
      case 0: RETVAL = _GMP_is_lucas_pseudoprime(n, 0); break;
      case 1: RETVAL = _GMP_is_lucas_pseudoprime(n, 1); break;
      case 2: RETVAL = _GMP_is_lucas_pseudoprime(n, 2); break;
      case 3: RETVAL = _GMP_is_frobenius_underwood_pseudoprime(n); break;
      case 4: RETVAL = _GMP_is_frobenius_khashin_pseudoprime(n); break;
      case 5:
      default:RETVAL = is_euler_plumb_pseudoprime(n); break;
    }
    mpz_clear(n);
  OUTPUT:
    RETVAL

int
is_almost_extra_strong_lucas_pseudoprime(IN char* strn, IN UV increment = 1)
  PREINIT:
    mpz_t n;
  CODE:
    if (increment == 0 || increment > 65535)
      croak("Increment parameter must be >0 and < 65536");
    PRIMALITY_START("is_almost_extra_strong_lucas_pseudoprime", 1, 0);
    RETVAL = _GMP_is_almost_extra_strong_lucas_pseudoprime(n, increment);
    mpz_clear(n);
  OUTPUT:
    RETVAL

int
is_frobenius_pseudoprime(IN char* strn, IN char* strP = 0, IN char* strQ = 0)
  PREINIT:
    mpz_t n, P, Q;
  CODE:
    if (items != 1 && items != 3)
      croak("is_frobenius_pseudoprime: expected 1 or 3 arguments");
    if (validate_and_set(n, IFLAG_ANY)) {
      RETVAL = 0;  /* Negative input */
    } else if (items == 1) {
      RETVAL = is_frobenius_pseudoprime(n);
    } else {
      validate_and_set(P, IFLAG_ANY);
      validate_and_set(Q, IFLAG_ANY);
      RETVAL = is_frobenius_pseudoprime_pq(n, P, Q);
      mpz_clear(Q); mpz_clear(P);
    }
    mpz_clear(n);
  OUTPUT:
    RETVAL

int
is_prime(IN char* strn)
  ALIAS:
    is_prob_prime = 1
    is_bpsw_prime = 2
    is_llr_prime = 3
    is_proth_prime = 4
    is_trial_prime = 5
    is_aks_prime = 6
    is_nminus1_prime = 7
    is_nplus1_prime = 8
    is_bls75_prime = 9
    is_ecpp_prime = 10
  PREINIT:
    mpz_t n;
    int ret;
  CODE:
    /* Returns arg for single-dig primes, 0 for multiples of 2, 3, 5, or neg */
    PRIMALITY_START("is_prime", (ix < 5) ? 2 : 1, 1);
    switch (ix) {
      case 0: ret = _GMP_is_prime(n); break;
      case 1: ret = _GMP_is_prob_prime(n); break;
      case 2: ret = _GMP_BPSW(n); break;
      case 3: ret = llr(n); break;
      case 4: ret = proth(n); break;
      case 5: ret = is_trial_prime(n); break;
      case 6: ret = is_aks_prime(n); break;
      case 7: ret = (BLS_primality_nm1(n, 100, 0) == 2) ? 1 : 0; break;
      case 8: ret = (BLS_primality_np1(n, 100, 0) == 2) ? 1 : 0; break;
      case 9: ret = (BLS_primality(n, 100, 0) == 2) ? 1 : 0; break;
      case 10:
      default:ret = (_GMP_ecpp(n, 0) == 2) ? 1 : 0; break;
    }
    RETVAL = ret;
    mpz_clear(n);
  OUTPUT:
    RETVAL

int is_perrin_pseudoprime(IN char* strn, IN int type = 0)
  PREINIT:
    mpz_t n;
  CODE:
    PRIMALITY_START("is_perrin_pseudoprime", 1, 0);
    RETVAL = is_perrin_pseudoprime(n, type);  /* Restricted or not */
    mpz_clear(n);
  OUTPUT:
    RETVAL

int is_miller_prime(IN char* strn, IN int assumegrh = 0)
  PREINIT:
    mpz_t n;
  CODE:
    PRIMALITY_START("is_miller_prime", 2, 1);
    RETVAL = is_miller_prime(n, assumegrh);
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
      XPUSH_INT(result);
    } else {
      char* prooftext = 0;
      result = _GMP_is_provable_prime(n, &prooftext);
      XPUSH_INT(result);
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
    mpz_t a, b, n, px, py, m, q, t1, t2;
  CODE:
    validate_and_set(a, IFLAG_NONNEG);
    validate_and_set(b, IFLAG_NONNEG);  /* Unused */
    validate_and_set(n, IFLAG_NONNEG);
    validate_and_set(px,IFLAG_NONNEG);
    validate_and_set(py,IFLAG_NONNEG);
    validate_and_set(m, IFLAG_NONNEG);
    validate_and_set(q, IFLAG_NONNEG);
    mpz_init(t1);  mpz_init(t2);
    RETVAL = (ecpp_check_point(px, py, m, q, a, n, t1, t2) == 2) ? 1 : 0;
    mpz_clear(t1); mpz_clear(t2);
    mpz_clear(a); mpz_clear(b); mpz_clear(n); mpz_clear(px); mpz_clear(py);
    mpz_clear(m); mpz_clear(q);
  OUTPUT:
    RETVAL

int
is_almost_prime(IN char* strk, IN char* strn)
  PREINIT:
    mpz_t k, n;
  CODE:
    validate_and_set(k, IFLAG_NONNEG);
    validate_and_set(n, IFLAG_ANY);
    if (mpz_sgn(n) <= 0)
      RETVAL = 0;
    else if (mpz_fits_ulong_p(k) && mpz_cmp_ui(k, UINT32_MAX) <= 0)
      RETVAL = is_almost_prime((uint32_t)mpz_get_ui(k), n);
    else if (mpz_sizeinbase(n,2) <= (size_t)UINT32_MAX)
      RETVAL = 0;
    else {
      mpz_clear(k); mpz_clear(n);
      croak("is_almost_prime: k too large");
    }
    mpz_clear(k);
    mpz_clear(n);
  OUTPUT:
    RETVAL

UV is_power(IN char* strn, IN SV* svk = 0)
  PREINIT:
    mpz_t n, k;
    int isneg, have_k = 0, auto_k = 0;
  CODE:
    RETVAL = 0;
    isneg = validate_and_set(n, IFLAG_ABS);
    if (items == 1 || !SvOK(svk)) {
      auto_k = 1;
    } else {
      set_integer_string(k, "k", SvPV_nolen(svk), IFLAG_NONNEG);
      have_k = 1;
      if (mpz_sgn(k) == 0) {
        auto_k = 1;
      } else if (mpz_cmp_ui(k, 1) == 0) {
        RETVAL = 1;
      } else if (isneg && mpz_even_p(k)) {
        RETVAL = 0;
      } else if (mpz_cmp_ui(n, 1) == 0) {
        RETVAL = 1;
      } else if (!mpz_fits_uv_p(k)) {
        RETVAL = 0;
      } else {
        RETVAL = is_power(n, mpz_get_uv(k));
      }
    }
    if (auto_k) {
      RETVAL = is_power(n, 0);
      if (isneg && RETVAL != 0) {
        UV r = RETVAL;
        while (!(r & 1)) r >>= 1;
        RETVAL = (r == 1) ? 0 : r;
      }
    }
    if (have_k) mpz_clear(k);
    mpz_clear(n);
  OUTPUT:
    RETVAL

int is_divisible(IN char* strn, IN char* strd, ...)
  PREINIT:
    mpz_t n, d;
    size_t i;
  CODE:
    validate_and_set(n, IFLAG_ANY);
    validate_and_set(d, IFLAG_ANY);
    RETVAL = !!mpz_divisible_p(n, d);
    for (i = 2; i < (size_t)items && RETVAL == 0; i++) {
      mpz_clear(d);
      set_integer_string(d, "d", SvPV_nolen(ST(i)), IFLAG_ANY);
      RETVAL = !!mpz_divisible_p(n, d);
    }
    mpz_clear(d);  mpz_clear(n);
  OUTPUT:
    RETVAL

int is_congruent(IN char* strn, IN char* strc, IN char* strd)
  PREINIT:
    mpz_t n, c, d;
  CODE:
    validate_and_set(n, IFLAG_ANY);
    validate_and_set(c, IFLAG_ANY);
    validate_and_set(d, IFLAG_ABS);
    RETVAL = !!mpz_congruent_p(n, c, d);
    mpz_clear(d);  mpz_clear(c);  mpz_clear(n);
  OUTPUT:
    RETVAL

void
next_prime(IN char* strn)
  ALIAS:
    prev_prime = 1
    next_twin_prime = 2
  PREINIT:
    mpz_t n;
  PPCODE:
    validate_and_set(n, IFLAG_NONNEG);
    if (ix == 1 && mpz_cmp_ui(n,3) < 0) { mpz_clear(n); XSRETURN_UNDEF; }
    if      (ix == 0) _GMP_next_prime(n);
    else if (ix == 1) _GMP_prev_prime(n);
    else              next_twin_prime(n, n);
    XPUSH_MPZ(n);
    mpz_clear(n);


void
random_prime(IN char* strlo, IN char* strhi = 0)
  PREINIT:
    mpz_t lo, hi, res;
    int retundef;
  PPCODE:
    if (items == 1) {
      validate_and_set(lo, IFLAG_NONNEG);
      mpz_init_set(hi, lo);
      mpz_set_ui(lo,0);
    } else {
      validate_and_set(lo, IFLAG_NONNEG);
      validate_and_set(hi, IFLAG_NONNEG);
    }
    mpz_init(res);
    retundef = !mpz_random_prime(res, lo, hi);
    if (!retundef) XPUSH_MPZ(res);
    mpz_clear(res); mpz_clear(hi); mpz_clear(lo);
    if (retundef) XSRETURN_UNDEF;

void
urandomr(IN char* strlo, IN char* strhi)
  PREINIT:
    mpz_t lo, hi, res;
  PPCODE:
    validate_and_set(lo, IFLAG_ANY);
    validate_and_set(hi, IFLAG_ANY);
    if (mpz_cmp(lo,hi) > 0) {
      mpz_clear(lo); mpz_clear(hi);
      XSRETURN_UNDEF;
    }
    if (mpz_sgn(lo) >= 0 && mpz_sgn(hi) >= 0 && mpz_sizeinbase(hi,2) <= 32) {
      uint32_t ulo = mpz_get_ui(lo),  uhi = mpz_get_ui(hi);
      if (uhi - ulo < UINT32_MAX) {
        mpz_clear(lo); mpz_clear(hi);
        XSRETURN_UV( ulo + isaac_rand(uhi-ulo+1) );
      }
    }
    mpz_init(res);
    mpz_sub(hi,hi,lo);
    mpz_add_ui(hi,hi,1);
    mpz_isaac_urandomm(res, hi);
    mpz_add(res,res,lo);
    XPUSH_MPZ(res);
    mpz_clear(res); mpz_clear(hi); mpz_clear(lo);

void prime_count(IN char* strlo, IN char* strhi = 0)
  ALIAS:
    prime_power_count = 1
    perfect_power_count = 2
  PREINIT:
    mpz_t lo, hi, res;
  PPCODE:
    mpz_init(res);
    if (items == 1) {
      validate_and_set(lo, IFLAG_NONNEG);
      switch (ix) {
        case  0: prime_count(res, lo); break;
        case  1: prime_power_count(res, lo); break;
        case  2: perfect_power_count(res, lo); break;
        default: break;
      }
      mpz_clear(lo);
    } else {
      validate_and_set(lo, IFLAG_NONNEG);
      validate_and_set(hi, IFLAG_NONNEG);
      switch (ix) {
        case  0: prime_count_range(res, lo, hi); break;
        case  1: prime_power_count_range(res, lo, hi); break;
        case  2: perfect_power_count_range(res, lo, hi); break;
        default: break;
      }
      mpz_clear(lo);
      mpz_clear(hi);
    }
    XPUSH_MPZ(res);
    mpz_clear(res);

void legendre_phi(IN char* strn, IN char* stra)
  PREINIT:
    mpz_t n, a, res;
  PPCODE:
    validate_and_set(n, IFLAG_NONNEG);
    validate_and_set(a, IFLAG_NONNEG);
    mpz_init(res);
    legendre_phi(res, n, a);
    XPUSH_MPZ(res);
    mpz_clear(res);
    mpz_clear(a);
    mpz_clear(n);

void nth_perfect_power(IN char* strn)
  ALIAS:
    nth_perfect_power_approx = 1
  PREINIT:
    mpz_t n, res;
  PPCODE:
    validate_and_set(n, IFLAG_NONNEG);
    mpz_init(res);
    if (ix == 0) nth_perfect_power(res, n);
    else         nth_perfect_power_approx(res, n);
    XPUSH_MPZ(res);
    mpz_clear(n);
    mpz_clear(res);

void next_perfect_power(IN char* strn)
  ALIAS:
    prev_perfect_power = 1
  PREINIT:
    mpz_t n, res;
  PPCODE:
    validate_and_set(n, IFLAG_ANY);
    mpz_init(res);
    if (ix == 0) next_perfect_power(res, n);
    else         prev_perfect_power(res, n);
    XPUSH_MPZ(res);
    mpz_clear(n);
    mpz_clear(res);

void totient(IN char* strn)
  ALIAS:
    carmichael_lambda = 1
    ramanujan_tau = 2
    sqrtint = 3
    prime_count_lower = 4
    prime_count_upper = 5
  PREINIT:
    mpz_t n;
  PPCODE:
    validate_and_set(n, IFLAG_NONNEG);
    switch (ix) {
      case 0:  totient(n, n);            break;
      case 1:  carmichael_lambda(n, n);  break;
      case 2:  rtau(n, n);               break;
      case 3:  mpz_sqrt(n, n);           break;
      case 4:  prime_count_lower(n, n);  break;
      case 5:  prime_count_upper(n, n);  break;
      default: break;
    }
    XPUSH_MPZ(n);
    mpz_clear(n);

void absint(IN char* strn)
  PREINIT:
    mpz_t n;
  PPCODE:
    validate_and_set(n, IFLAG_ABS);
    XPUSH_MPZ(n);
    mpz_clear(n);

void urandomm(IN char* strn)
  PREINIT:
    mpz_t n;
  PPCODE:
    validate_and_set(n, IFLAG_NONNEG);
    mpz_isaac_urandomm(n, n);
    XPUSH_MPZ(n);
    mpz_clear(n);

void is_prime_power(IN char* strn)
  PREINIT:
    mpz_t n;
    UV    R;
  PPCODE:
    validate_and_set(n, IFLAG_ANY);
    R = (mpz_sgn(n) <= 0) ? 0 : prime_power(n, n);
    mpz_clear(n);
    XSRETURN_UV(R);

void signint(IN char* strn)
  PREINIT:
    mpz_t n;
    int res;
  PPCODE:
    validate_and_set(n, IFLAG_ANY);
    res = mpz_sgn(n);
    mpz_clear(n);
    XSRETURN_IV(res);

void cmpint(IN char* stra, IN char* strb)
  ALIAS:
    cmpabsint = 1
  PREINIT:
    mpz_t a, b;
    int res;
  PPCODE:
    validate_and_set(a, IFLAG_ANY);
    validate_and_set(b, IFLAG_ANY);
    res = (ix == 0) ? mpz_cmp(a, b) : mpz_cmpabs(a, b);
    /* GMP 6.2 changed to only return -1,0,1 */
    /* Enforce -1, 0, 1 as our only return values. */
    if (res < 0) res = -1;
    if (res > 0) res = 1;
    mpz_clear(a);
    mpz_clear(b);
    XSRETURN_IV(res);

void setbit(IN char* strn, IN UV k)
  ALIAS:
    clrbit = 1
    notbit = 2
    tstbit = 3
  PREINIT:
    mpz_t n;
    int res;
  PPCODE:
    validate_and_set(n, IFLAG_ANY);
    switch (ix) {
      case 0:  mpz_setbit(n, k); break;
      case 1:  mpz_clrbit(n, k); break;
      case 2:  mpz_combit(n, k); break;
      case 3:  res = mpz_tstbit(n, k); break;
      default: break;
    }
    if (ix != 3) XPUSH_MPZ(n);
    mpz_clear(n);
    if (ix == 3) XSRETURN_IV(res);

void bitand(IN char* stra, IN char* strb)
  ALIAS:
    bitor = 1
    bitxor = 2
  PREINIT:
    mpz_t a, b;
  PPCODE:
    validate_and_set(a, IFLAG_ANY);
    validate_and_set(b, IFLAG_ANY);
    switch (ix) {
      case 0:  mpz_and(a, a, b); break;
      case 1:  mpz_ior(a, a, b); break;
      case 2:  mpz_xor(a, a, b); break;
      default: break;
    }
    XPUSH_MPZ(a);
    mpz_clear(a);
    mpz_clear(b);

void bitnot(IN char* strn)
  ALIAS:
    negint = 1
    add1int = 2
    sub1int = 3
    exp_mangoldt = 4
  PREINIT:
    mpz_t n;
  PPCODE:
    validate_and_set(n, IFLAG_ANY);
    switch (ix) {
      case 0:  mpz_com(n, n);       break;
      case 1:  mpz_neg(n, n);       break;
      case 2:  mpz_add_ui(n, n, 1); break;
      case 3:  mpz_sub_ui(n, n, 1); break;
      case 4:  exp_mangoldt(n, n);  break;
      default: break;
    }
    XPUSH_MPZ(n);
    mpz_clear(n);

void bernfrac(IN char* strn)
  ALIAS:
    harmfrac = 1
  PREINIT:
    mpz_t n, d;
  PPCODE:
    validate_and_set(n, IFLAG_NONNEG);
    mpz_init(d);
    if (ix == 0) bernfrac(n,d,n);
    else         harmfrac(n,d,n);
    XPUSH_MPZ(n);
    XPUSH_MPZ(d);
    mpz_clear(d);
    mpz_clear(n);

void harmreal(IN char* strn, IN UV prec = 40)
  ALIAS:
    bernreal = 1
    logreal = 2
    expreal = 3
    zeta = 4
    li = 5
    ei = 6
    riemannr = 7
    lambertw = 8
    surround_primes = 9
  PREINIT:
    mpz_t n;
    mpf_t f;
    char* res;
  PPCODE:
    if (ix == 9) {  /* surround_primes */
      UV prev, next;
      validate_and_set(n, IFLAG_NONNEG);
      next = 1 + (mpz_sgn(n)==0);
      if (mpz_cmp_ui(n,2) > 0) {
        surround_primes(n, &prev, &next, (items == 1) ? 0 : prec);
        XPUSH_UINT(prev);
      } else {
        XPUSHs(sv_2mortal(newSV(0)));
      }
      XPUSHs(sv_2mortal(newSVuv(next)));
      mpz_clear(n);
    } else if (ix <= 1) {
      validate_and_set(n, IFLAG_NONNEG);
      res = (ix == 0) ? harmreal(n, prec) : bernreal(n, prec);
      mpz_clear(n);
      XPUSHs(sv_2mortal(newSVpv(res, 0)));
      Safefree(res);
    } else {
      unsigned long bits  = 64 + (unsigned long)(3.32193 * prec);
      unsigned long bits2 = 64 + (unsigned long)(3.32193 * strlen(strn));
      if (bits2 > bits) bits = bits2;
      mpf_init2(f, bits);
      if (mpf_set_str(f, strn, 10) != 0)
        croak("Not valid base-10 floating point input: %s", strn);
      res = 0;
      switch (ix) {
        case 2:  res = logreal(f, prec); break;
        case 3:  res = expreal(f, prec); break;
        case 4:  res = zetareal(f, prec); break;
        case 5:  res = lireal(f, prec); break;
        case 6:  res = eireal(f, prec); break;
        case 7:  res = riemannrreal(f, prec); break;
        case 8:
        default: res = lambertwreal(f, prec); break;
      }
      mpf_clear(f);
      if (res == 0)
        XSRETURN_UNDEF;
      XPUSHs(sv_2mortal(newSVpv(res, 0)));
      Safefree(res);
    }

void powreal(IN char* strn, IN char* strx, IN UV prec = 40)
  ALIAS:
    rootreal = 1
    agmreal = 2
    addreal = 3
    subreal = 4
    mulreal = 5
    divreal = 6
  PREINIT:
    mpf_t n, x;
    char* res;
    unsigned long bits, bits2, bits3;
  PPCODE:
    bits  = 64 + (unsigned long)(3.32193 * prec);
    bits2 = 64 + (unsigned long)(3.32193 * strlen(strn));
    bits3 = 64 + (unsigned long)(3.32193 * strlen(strx));
    if (bits2 > bits) bits = bits2;
    if (bits3 > bits) bits = bits3;
    mpf_init2(n, bits);
    if (mpf_set_str(n, strn, 10) != 0)
      croak("Not valid base-10 floating point input: %s", strn);
    mpf_init2(x, bits);
    if (mpf_set_str(x, strx, 10) != 0)
      croak("Not valid base-10 floating point input: %s", strx);
    if ( (ix == 0 && mpf_sgn(n) < 0 && !mpf_integer_p(x)) ||
         (ix == 0 && mpf_sgn(n) == 0 && mpf_sgn(x) < 0) ||
         (ix == 1 && mpf_sgn(x) == 0) ||
         (ix == 1 && mpf_sgn(n) < 0 && mpf_cmp_ui(x,1)!=0 && mpf_cmp_si(x,-1)!=0) )
      XSRETURN_UNDEF;
    switch (ix) {
      case 0:  res = powreal(n, x, prec);  break;
      case 1:  res = rootreal(n, x, prec); break;
      case 2:  res = agmreal(n, x, prec); break;
      case 3:  res = addreal(n, x, prec); break;
      case 4:  res = subreal(n, x, prec); break;
      case 5:  res = mulreal(n, x, prec); break;
      case 6:
      default: res = divreal(n, x, prec);  break;
    }
    mpf_clear(n);
    mpf_clear(x);
    if (res == 0)
      XSRETURN_UNDEF;
    XPUSHs(sv_2mortal(newSVpv(res, 0)));
    Safefree(res);

void bernvec(IN UV n)
  PREINIT:
    const mpz_t *N, *D;
    UV i;
  PPCODE:
    bernvec(&N, &D, n);  /* Cached array, do not destroy */
    if (GIMME_V != G_VOID) {
      EXTEND(SP, (long)(n+1));
      for (i = 0; i <= n; i++) {
        AV* av = newAV();
        av_push(av, sv_return_for_mpz(aTHX_ N[i]));
        av_push(av, sv_return_for_mpz(aTHX_ D[i]));
        PUSHs( sv_2mortal(newRV_noinc( (SV*) av )) );
      }
    }

void
gcd(...)
  PROTOTYPE: @
  ALIAS:
    lcm = 1
    vecsum = 2
    vecprod = 3
  PREINIT:
    int i, negflag;
    mpz_t ret, n;
  PPCODE:
    if (items == 0) XSRETURN_IV( (ix == 1 || ix == 3) ? 1 : 0);
    negflag = (ix <= 1) ? IFLAG_ABS : IFLAG_ANY;
    if (ix == 1 || ix == 3) {
      mpz_t* list;
      New(0, list, items, mpz_t);
      for (i = 0; i < items; i++) {
        char* strn = SvPV_nolen(ST(i));
        set_integer_string(list[i], "arg", strn, negflag);
      }
      if (ix == 1) mpz_veclcm(list, 0, items-1);
      else         mpz_product(list, 0, items-1);
      XPUSH_MPZ(list[0]);
      for (i = 0; i < items; i++)  mpz_clear(list[i]);
      Safefree(list);
      XSRETURN(1);
    }
    mpz_init_set_ui(ret, (ix == 1 || ix == 3) ? 1 : 0);
    for (i = 0; i < items; i++) {
      char* strn = SvPV_nolen(ST(i));
      set_integer_string(n, "arg", strn, negflag);
      switch (ix) {
        case 0:  mpz_gcd(ret, ret, n); break;
        case 1:  mpz_lcm(ret, ret, n); break;
        case 2:  mpz_add(ret, ret, n); break;
        case 3:
        default: mpz_mul(ret, ret, n); break;
      }
      mpz_clear(n);
    }
    XPUSH_MPZ(ret);
    mpz_clear(ret);

void
vecprefixsum(...)
  PROTOTYPE: @
  PREINIT:
    AV *av;
    SV *sv;
    SV **svp;
    int plen;
    size_t i, len;
    mpz_t sum, n;
  PPCODE:
    av = 0;
    plen = -1;
    if (items > 0 && SvROK(ST(0)) && SvTYPE(SvRV(ST(0))) == SVt_PVAV) {
      if (items != 1)
        croak("vecprefixsum: expected integer list or single array reference");
      av = (AV*) SvRV(ST(0));
      plen = av_len(av);
      len = (plen < 0) ? 0 : (size_t)plen + 1;
    } else {
      len = (size_t)items;
    }

    if (GIMME_V == G_ARRAY) {
      EXTEND(SP, (long)len);
      mpz_init_set_ui(sum, 0);
    }
    for (i = 0; i < len; i++) {
      if (av != 0) {
        svp = av_fetch(av, (I32)i, 0);
        if (svp == 0 || !SvOK(*svp)) {
          if (GIMME_V == G_ARRAY) mpz_clear(sum);
          croak("Parameter must be defined");
        }
        sv = *svp;
      } else {
        sv = ST(i);
        if (!SvOK(sv)) {
          if (GIMME_V == G_ARRAY) mpz_clear(sum);
          croak("Parameter must be defined");
        }
      }
      set_integer_string(n, "arg", SvPV_nolen(sv), IFLAG_ANY);
      if (GIMME_V == G_ARRAY) {
        mpz_add(sum, sum, n);
        XPUSH_MPZ(sum);
      }
      mpz_clear(n);
    }
    if (GIMME_V == G_ARRAY)
      mpz_clear(sum);
    else if (GIMME_V == G_SCALAR)
      XPUSH_UINT(len);

int
kronecker(IN char* stra, IN char* strb)
  ALIAS:
    valuation = 1
    is_gaussian_prime = 2
    is_smooth = 3
    is_rough = 4
  PREINIT:
    mpz_t a, b;
  CODE:
    validate_and_set(a, IFLAG_ANY);
    validate_and_set(b, IFLAG_ANY);
    if (ix != 0) {
      mpz_abs(a,a);
      mpz_abs(b,b);
    }
    RETVAL = 0;
    switch (ix) {
      case 0: RETVAL = mpz_kronecker(a, b);
              break;
      case 1: if (mpz_cmp_ui(b,2) < 0)  croak("valuation: k must be > 1");
              if (mpz_cmp_ui(a,0) == 0) XSRETURN_UNDEF;
              RETVAL = (mpz_cmp_ui(b,2) == 0)  ?  mpz_scan1(a, 0)
                                               :  mpz_remove(a, a, b);
              break;
      case 2: if (mpz_sgn(a) == 0) {
                RETVAL = (mpz_fdiv_ui(b,4) == 3) ? _GMP_is_prime(b) : 0;
              } else if (mpz_sgn(b) == 0) {
                RETVAL = (mpz_fdiv_ui(a,4) == 3) ? _GMP_is_prime(a) : 0;
              } else {
                mpz_mul(a, a, a);
                mpz_mul(b, b, b);
                mpz_add(a, a, b);
                RETVAL = (!mpz_cmp_ui(a,2))      ?  2
                       : (mpz_fdiv_ui(a,4) == 1) ?  _GMP_is_prime(a)
                                                 :  0;
              }
              break;
      case 3: RETVAL = is_smooth(a, b);
              break;
      case 4: RETVAL = is_rough(a, b);
              break;
      default:break;
    }
    mpz_clear(b);
    mpz_clear(a);
  OUTPUT:
    RETVAL

void remove_factors(IN char* strn, IN char* strk)
  ALIAS:
    remove_factors_exp = 1
  PREINIT:
    mpz_t n, k;
  PPCODE:
    validate_and_set(n, IFLAG_ANY);
    validate_and_set(k, IFLAG_NONNEG);
    if (mpz_cmp_ui(k,2) < 0) {
      mpz_clear(k); mpz_clear(n);
      croak("%s: k must be > 1", SUBNAME);
    }
    if (mpz_sgn(n) == 0) {
      XPUSHs(&PL_sv_undef);
      if (ix == 1) XPUSHs(&PL_sv_undef);
    } else {
      UV e = mpz_remove(n, n, k);
      XPUSH_MPZ(n);
      if (ix == 1) XPUSH_UINT(e);
    }
    mpz_clear(k);
    mpz_clear(n);

void moebius(IN char* strn, IN char* strnhi = 0)
  PREINIT:
    mpz_t n, nhi;
  PPCODE:
    validate_and_set(n, IFLAG_ANY);
    if (items == 1) {
      XPUSH_INT(moebius(n));
    } else {
      validate_and_set(nhi, IFLAG_ANY);
      if (GIMME_V != G_ARRAY) {
        if (mpz_cmp(n,nhi) > 0) { mpz_set_ui(n,0); }
        else                    { mpz_sub(nhi,nhi,n); mpz_add_ui(n,nhi,1); }
        XPUSH_MPZ(n);
      } else {
        while (mpz_cmp(n, nhi) <= 0) {
          XPUSH_INT(moebius(n));
          mpz_add_ui(n, n, 1);
        }
      }
      mpz_clear(nhi);
    }
    mpz_clear(n);

void euler_phi(IN char* strn, IN char* strnhi = 0)
  PREINIT:
    mpz_t n, nhi, r;
  PPCODE:
    validate_and_set(n, IFLAG_ANY);
    if (items == 1) {
      totient(n, n);
      XPUSH_MPZ(n);
    } else {
      validate_and_set(nhi, IFLAG_ANY);
      if (GIMME_V != G_ARRAY) {
        if (mpz_cmp(n,nhi) > 0) { mpz_set_ui(n,0); }
        else                    { mpz_sub(nhi,nhi,n); mpz_add_ui(n,nhi,1); }
        XPUSH_MPZ(n);
      } else {
        mpz_init(r);
        while (mpz_cmp(n, nhi) <= 0) {
          totient(r, n);
          XPUSH_MPZ(r);
          mpz_add_ui(n, n, 1);
        }
        mpz_clear(r);
      }
      mpz_clear(nhi);
    }
    mpz_clear(n);

void lucasu(IN char* strp, IN char* strq, IN char* strk)
  ALIAS:
    lucasv = 1
    lucasuv = 2
  PREINIT:
    mpz_t u, v, p, q, k;
  PPCODE:
    validate_and_set(p, IFLAG_ANY);
    validate_and_set(q, IFLAG_ANY);
    validate_and_set(k, IFLAG_NONNEG);
    mpz_init(u);  mpz_init(v);
    lucasuv(u, v, p, q, k);
    switch (ix) {
      case 0:  XPUSH_MPZ(u);  break;
      case 1:  XPUSH_MPZ(v);  break;
      case 2:
      default: XPUSH_MPZ(u);  XPUSH_MPZ(v); break;
    }
    mpz_clear(v); mpz_clear(u);
    mpz_clear(k); mpz_clear(q); mpz_clear(p);

void lucasumod(IN char* strp, IN char* strq, IN char* strk, IN char* strn)
  ALIAS:
    lucasvmod = 1
    lucasuvmod = 2
  PREINIT:
    mpz_t u, v, t, p, q, k, n;
  PPCODE:
    validate_and_set(p, IFLAG_ANY);
    validate_and_set(q, IFLAG_ANY);
    validate_and_set(k, IFLAG_NONNEG);
    validate_and_set(n, IFLAG_ABS);
    if (mpz_cmpabs_ui(n,1) <= 0) {
      int retundef = (mpz_sgn(n) == 0);
      mpz_clear(n); mpz_clear(k); mpz_clear(q); mpz_clear(p);
      if (retundef)     XSRETURN_UNDEF;
      else if (ix != 2) XSRETURN_IV(0);
      else              { XPUSH_UINT(0); XPUSH_UINT(0); XSRETURN(2); }
    }
    mpz_init(t);
    if (ix == 0 || ix == 2) mpz_init(u);
    if (ix == 1 || ix == 2) mpz_init(v);
    switch (ix) {
      case 0:  lucasumod(u, p, q, k, n, t);  XPUSH_MPZ(u);  break;
      case 1:  lucasvmod(v, p, q, k, n, t);  XPUSH_MPZ(v);  break;
      case 2:
      default: lucasuvmod(u, v, p, q, k, n, t);
               XPUSH_MPZ(u);  XPUSH_MPZ(v);
               break;
    }
    if (ix == 0 || ix == 2) mpz_clear(u);
    if (ix == 1 || ix == 2) mpz_clear(v);
    mpz_clear(t);
    mpz_clear(n); mpz_clear(k); mpz_clear(q); mpz_clear(p);

void catalan_number(IN char* strn)
  PREINIT:
    mpz_t n;
    unsigned long un;
  PPCODE:
    validate_and_set(n, IFLAG_NONNEG);
    if (!mpz_fits_ulong_p(n))
      croak("catalan_number: argument too large");
    un = mpz_get_ui(n);
    if (un <= 1) {
      mpz_set_ui(n, 1);
    } else {
      mpz_mul_2exp(n, n, 1);        /* n = 2*un as mpz */
      mpz_bin_ui(n, n, un);         /* C(2un, un) */
      mpz_divexact_ui(n, n, un+1);  /* / (un+1) */
    }
    XPUSH_MPZ(n);
    mpz_clear(n);

void bell_number(IN char* strn)
  PREINIT:
    mpz_t n;
  PPCODE:
    validate_and_set(n, IFLAG_NONNEG);
    if (!mpz_fits_ulong_p(n))
      croak("bell_number: argument too large");
    bell_number(n, mpz_get_ui(n));
    XPUSH_MPZ(n);
    mpz_clear(n);

void fubini(IN char* strn)
  PREINIT:
    mpz_t n;
  PPCODE:
    validate_and_set(n, IFLAG_NONNEG);
    if (!mpz_fits_ulong_p(n))
      croak("fubini: argument too large");
    fubini(n, mpz_get_ui(n));
    XPUSH_MPZ(n);
    mpz_clear(n);

void fibonacci(IN char* strk)
  ALIAS:
    lucas_number = 1
  PREINIT:
    mpz_t k;
    unsigned long uk;
    int isneg;
  PPCODE:
    isneg = validate_and_set(k, IFLAG_ABS);
    if (!mpz_fits_ulong_p(k))
      croak("%s: argument too large", GvNAME(CvGV(cv)));
    uk = mpz_get_ui(k);
    if (ix == 0) {
      mpz_fib_ui(k, uk);
      if (isneg && (uk & 1) == 0) mpz_neg(k, k);  /* F(-n) = -F(n) if n even */
    } else {
      mpz_lucnum_ui(k, uk);
      if (isneg && (uk & 1) == 1) mpz_neg(k, k);  /* L(-n) = -L(n) if n odd */
    }
    XPUSH_MPZ(k);
    mpz_clear(k);

int
liouville(IN char* strn)
  PREINIT:
    mpz_t n;
  CODE:
    validate_and_set(n, IFLAG_NONNEG);
    RETVAL = liouville(n);
    mpz_clear(n);
  OUTPUT:
    RETVAL

int
is_square(IN char* strn)
  ALIAS:
    is_semiprime = 1
    is_totient = 2
    is_carmichael = 3
    is_practical = 4
    is_fundamental = 5
    is_perfect_power = 6
  PREINIT:
    mpz_t n;
    int isneg;
  CODE:
    isneg = validate_and_set(n, IFLAG_ANY);
    if (isneg && ix <= 4) {
      RETVAL = 0;
    } else {
      switch (ix) {
        case 0:  RETVAL = is_power(n,2);      break;
        case 1:  RETVAL = is_semiprime(n);    break;
        case 2:  RETVAL = is_totient(n);      break;
        case 3:  RETVAL = is_carmichael(n);   break;
        case 4:  RETVAL = is_practical(n);    break;
        case 5:  RETVAL = is_fundamental(n);  break;
        case 6:  RETVAL = is_perfect_power(n);break;
        default: RETVAL = 0; break;
      }
    }
    mpz_clear(n);
  OUTPUT:
    RETVAL

int
prime_omega(IN char* strn)
  ALIAS:
    prime_bigomega = 1
    hammingweight = 2
    is_square_free = 3
  PREINIT:
    mpz_t n;
  CODE:
    validate_and_set(n, IFLAG_ABS);
    switch (ix) {
      case 0:  RETVAL = omega(n);          break;
      case 1:  RETVAL = bigomega(n);       break;
      case 2:  RETVAL = mpz_popcount(n);   break;
      case 3:  RETVAL = is_square_free(n); break;
      default: RETVAL = 0; break;
    }
    mpz_clear(n);
  OUTPUT:
    RETVAL

void
sopf(IN char* strn)
  ALIAS:
    sopfr = 1
    dedekind_psi = 2
    aliquot_sum = 3
    abundance = 4
  PREINIT:
    mpz_t n;
  PPCODE:
    validate_and_set(n, ix == 2 ? IFLAG_ANY : IFLAG_NONNEG);
    switch (ix) {
      case 0:  sopf(n, n);          break;
      case 1:  sopfr(n, n);         break;
      case 2:  dedekind_psi(n, n);  break;
      case 3:  aliquot_sum(n, n);   break;
      case 4:  abundance(n, n);     break;
      default: break;
    }
    XPUSH_MPZ(n);
    mpz_clear(n);

void
prime_signature(IN char* strn)
  PREINIT:
    mpz_t n, r;
    uint32_t *signature;
    int i, nsig;
  PPCODE:
    validate_and_set(n, IFLAG_NONNEG);
    if (GIMME_V == G_ARRAY) {
      nsig = prime_signature(0, &signature, n);
      EXTEND(SP, nsig);
      for (i = 0; i < nsig; i++)
        XPUSH_UINT(signature[i]);
      if (signature != 0) Safefree(signature);
    } else {
      mpz_init(r);
      prime_signature(r, 0, n);
      XPUSH_MPZ(r);
      mpz_clear(r);
    }
    mpz_clear(n);

int
is_safe_prime(IN char* strn)
  PREINIT:
    mpz_t n;
  CODE:
    validate_and_set(n, IFLAG_ANY);
    RETVAL = is_safe_prime(n);
    mpz_clear(n);
  OUTPUT:
    RETVAL

int
is_powerful(IN char* strn, IN UV k = 2)
  ALIAS:
    is_powerfree = 1
  PREINIT:
    mpz_t n;
  CODE:
    validate_and_set(n, IFLAG_ANY);
    switch (ix) {
      case 0:  RETVAL = is_powerful(n,k);   break;
      case 1:  RETVAL = is_powerfree(n,k);  break;
      default: RETVAL = 0; break;
    }
    mpz_clear(n);
  OUTPUT:
    RETVAL

void
next_powerfree(IN char* strn, IN UV k = 2)
  ALIAS:
    prev_powerfree = 1
    nth_powerfree = 2
  PREINIT:
    mpz_t n;
    int retundef;
  PPCODE:
    validate_and_set(n, IFLAG_NONNEG);
    switch (ix) {
      case 0:  next_powerfree(n,n,k);  break;
      case 1:  prev_powerfree(n,n,k);  break;
      case 2:  nth_powerfree(n,n,k);  break;
      default: break;
    }
    retundef = (mpz_sgn(n) <= 0);
    if (!retundef) XPUSH_MPZ(n);
    mpz_clear(n);
    if (retundef) XSRETURN_UNDEF;

void addint(IN char* stra, IN char* strb)
  ALIAS:
    subint = 1
    mulint = 2
    powint = 3
    divint = 4
    modint = 5
    cdivint = 6
    divrem = 7
    tdivrem = 8
    fdivrem = 9
    cdivrem = 10
  PREINIT:
    mpz_t a, b;
  PPCODE:
    validate_and_set(a, IFLAG_ANY);
    validate_and_set(b, IFLAG_ANY);
    if (ix >= 4 && mpz_sgn(b) == 0)
      croak("%s: divide by zero", SUBNAME);
    switch (ix) {
      case 0: mpz_add(a, a, b); break;
      case 1: mpz_sub(a, a, b); break;
      case 2: mpz_mul(a, a, b); break;
      case 3: if (mpz_sgn(b) < 0)
                croak("powint: exponent must be non-negative");
              else if (mpz_sgn(a) == 0)
                mpz_set_si(a, mpz_sgn(b) == 0);
              else if (mpz_cmp_si(a,1) == 0)
                mpz_set_si(a, 1);
              else if (mpz_cmp_si(a,-1) == 0)
                mpz_set_si(a, mpz_odd_p(b) ? -1 : 1);
              else if (!mpz_fits_ulong_p(b))
                croak("powint: exponent must be <= ULONG_MAX");
              else
                mpz_pow_ui(a, a, mpz_get_ui(b));
              break;
      case 4: mpz_fdiv_q(a, a, b); break;
      case 5: mpz_fdiv_r(a, a, b); break;
      case 6: mpz_cdiv_q(a, a, b); break;
      case 7: mpz_ediv_qr(b, a, a, b); break;
      case 8: mpz_tdiv_qr(b, a, a, b); break;
      case 9: mpz_fdiv_qr(b, a, a, b); break;
      case 10:mpz_cdiv_qr(b, a, a, b); break;
      default:break;
    }
    if (ix >= 7)  XPUSH_MPZ(b);
    XPUSH_MPZ(a);
    mpz_clear(b); mpz_clear(a);

void rootint(IN char* strn, IN char* strk)
  PREINIT:
    mpz_t n, k;
  PPCODE:
    validate_and_set(n, IFLAG_NONNEG);
    validate_and_set(k, IFLAG_POS);
    mpz_rootint(n, n, k);
    XPUSH_MPZ(n);
    mpz_clear(n); mpz_clear(k);

void logint(IN char* strn, IN char* strb)
  PREINIT:
    mpz_t n, b;
  PPCODE:
    validate_and_set(n, IFLAG_POS);
    validate_and_set(b, IFLAG_POS);
    if (mpz_cmp_ui(b,2) < 0) croak("logint: base must be at least 2");
    mpz_logint(n, n, b);
    XPUSH_MPZ(n);
    mpz_clear(n); mpz_clear(b);

void invmod(IN char* stra, IN char* strn)
  ALIAS:
    negmod = 1
    sqrtmod = 2
    factorialmod = 3
    znorder = 4
    is_qr = 5
    is_primitive_root = 6
  PREINIT:
    mpz_t a, n;
    int retundef;
  PPCODE:
    validate_and_set(a, IFLAG_ANY);
    validate_and_set(n, IFLAG_ABS);
    if (mpz_sgn(n) == 0) {
      mpz_clear(n); mpz_clear(a);
      XSRETURN_UNDEF;
    }
    if (mpz_cmp_ui(n,1) == 0) {
      mpz_clear(n); mpz_clear(a);
      XSRETURN_UV(ix >= 4 ? 1 : 0);
    }
    retundef = 0;
    switch (ix) {
      case 0:  if (!mpz_sgn(a))  retundef = 1;
               else              retundef = !mpz_invert(a, a, n);  break;
      case 1:  mpz_neg(a,a);  mpz_mod(a,a,n);  break;
      case 2:  retundef = !sqrtmod(a, a, n);   break;
      case 3:  factorialmod(a, a, n);          break;
      case 4:  retundef = !znorder(a, a, n);   break;
      case 5:  mpz_set_ui(a, is_qr(a,n));      break;
      case 6:
      default: mpz_set_ui(a, is_primitive_root(a,n,0));  break;
    }
    if (retundef) {
      mpz_clear(n); mpz_clear(a);
      XSRETURN_UNDEF;
    }
    XPUSH_MPZ(a);
    mpz_clear(n); mpz_clear(a);

void znprimroot(IN char* strn)
  PREINIT:
    mpz_t n;
  PPCODE:
    validate_and_set(n, IFLAG_ABS);
    znprimroot(n, n);
    if (mpz_sgn(n) < 0)
      { mpz_clear(n);  XSRETURN_UNDEF; }
    XPUSH_MPZ(n);
    mpz_clear(n);

void znlog(IN char* stra, IN char* strg, IN char* strn)
  PREINIT:
    mpz_t a, g, n, r;
    int ok;
  PPCODE:
    validate_and_set(a, IFLAG_ANY);
    validate_and_set(g, IFLAG_ANY);
    validate_and_set(n, IFLAG_ABS);
    mpz_init(r);
    ok = _GMP_znlog(r, a, g, n);
    if (ok) XPUSH_MPZ(r);
    mpz_clear(r); mpz_clear(n); mpz_clear(g); mpz_clear(a);
    if (!ok) XSRETURN_UNDEF;

void multifactorial(IN char* strn, IN char* strm)
  PREINIT:
    mpz_t n, m;
  PPCODE:
    validate_and_set(n, IFLAG_NONNEG);
    validate_and_set(m, IFLAG_POS);
    mpz_mfac(n, n, m);
    XPUSH_MPZ(n);
    mpz_clear(m); mpz_clear(n);

void is_polygonal(IN char* stra, IN char* strb)
  ALIAS:
    polygonal_nth = 1
  PREINIT:
    mpz_t a, b;
  PPCODE:
    validate_and_set(a, IFLAG_ANY);
    validate_and_set(b, IFLAG_NONNEG);
    if (mpz_cmp_ui(b,3) < 0) croak("%s: k must be >= 3",SUBNAME);
    polygonal_nth(a, a, b);
    if (ix == 0) {
      int ret = mpz_sgn(a) > 0;
      mpz_clear(b); mpz_clear(a);
      XSRETURN_IV(ret);
    }
    XPUSH_MPZ(a);
    mpz_clear(b); mpz_clear(a);

void binomial(IN char* stra, IN char* strb)
  PREINIT:
    mpz_t a, b;
    unsigned long n, k;
  PPCODE:
    validate_and_set(a, IFLAG_ANY);
    validate_and_set(b, IFLAG_ANY);
    if (mpz_sgn(a) >= 0) {
      int reflect = 0;
      if (mpz_sgn(b) < 0 || mpz_cmp(b,a) > 0)
        { mpz_clear(a); mpz_clear(b); XSRETURN_IV(0); }
      /* Check reflection to possibly reduce b */
      mpz_mul_2exp(b,b,1);
      reflect = mpz_cmp(b,a) > 0;
      mpz_tdiv_q_2exp(b,b,1);
      if (reflect)
        mpz_sub(b, a, b);
    } else {
      if (mpz_sgn(b) < 0) {
        if (mpz_cmp(b,a) > 0) { mpz_clear(a); mpz_clear(b); XSRETURN_IV(0); }
        mpz_sub(b, a, b);
      }
    }
    if (mpz_sgn(b) < 0)
      croak("binomial internal error: k should be non-negative here");

    if (mpz_sgn(b) == 0) { mpz_clear(a); mpz_clear(b); XSRETURN_IV(1); }

    if (!mpz_fits_ulong_p(b))
      croak("binomial: k must fit in native unsigned long");
    k = mpz_get_ui(b);

    if (k == 1 || !mpz_cmp_ui(a,k-1)) {
      /* result = a */
    } else if (mpz_sgn(a) >= 0 && mpz_fits_ulong_p(a)) {
      n = mpz_get_ui(a);
      /* Note: mpz_bin_uiui is *much* faster than mpz_bin_ui.
       * It is a bit faster than our code for small values, and
       * a tiny bit slower for larger values. */
      if (n > 50000 && k > 1000)  binomial(a, n, k);
      else                        mpz_bin_uiui(a, n, k);
    } else {
      mpz_bin_ui(a, a, k);
    }
    XPUSH_MPZ(a);
    mpz_clear(b); mpz_clear(a);

void gcdext(IN char* stra, IN char* strb)
  PREINIT:
    mpz_t a, b, t;
  PPCODE:
    validate_and_set(a, IFLAG_ANY);
    validate_and_set(b, IFLAG_ANY);
    mpz_init(t);
    if (mpz_sgn(a) == 0 && mpz_sgn(b) == 0) {
      mpz_set_ui(t, 0);  /* This changed in GMP 5.1.2.  Enforce new result. */
    } else {
      gcdext(a, t, b, a, b);
    }
    XPUSH_MPZ(t);  XPUSH_MPZ(b);  XPUSH_MPZ(a);
    mpz_clear(t); mpz_clear(b); mpz_clear(a);

void muladdint(IN char* stra, IN char* strb, IN char* strc)
  ALIAS:
    mulsubint = 1
    addmulint = 2
    submulint = 3
  PREINIT:
    mpz_t a, b, c;
  PPCODE:
    validate_and_set(a, IFLAG_ANY);
    validate_and_set(b, IFLAG_ANY);
    validate_and_set(c, IFLAG_ANY);
    if      (ix == 0) { mpz_mul(a, a, b); mpz_add(a, a, c); } /* a * b + c */
    else if (ix == 1) { mpz_mul(a, a, b); mpz_sub(a, a, c); } /* a * b - c */
    else if (ix == 2) { mpz_addmul(a, b, c); }                /* a + b * c */
    else              { mpz_submul(a, b, c); }                /* a - b * c */
    XPUSH_MPZ(a);
    mpz_clear(c);  mpz_clear(b);  mpz_clear(a);

void
binomialmod(IN char* strn, IN char* strk, IN char* strm)
  PREINIT:
    mpz_t n, k, m;
  PPCODE:
    validate_and_set(n, IFLAG_ANY);
    validate_and_set(k, IFLAG_ANY);
    validate_and_set(m, IFLAG_ABS);
    if (mpz_sgn(m) == 0) {
      mpz_clear(n); mpz_clear(k); mpz_clear(m);
      XSRETURN_UNDEF;
    }
    binomialmod(n, n, k, m);
    XPUSH_MPZ(n);
    mpz_clear(m);
    mpz_clear(k);
    mpz_clear(n);

void
falling_factorial(IN char* strx, IN char* strn)
  ALIAS:
    rising_factorial = 1
  PREINIT:
    mpz_t x, n, r;
  PPCODE:
    validate_and_set(x, IFLAG_ANY);
    validate_and_set(n, IFLAG_NONNEG);
    mpz_init(r);
    if (ix == 0)  falling_factorial(r, x, n);
    else          rising_factorial(r, x, n);
    mpz_clear(x);
    mpz_clear(n);
    XPUSH_MPZ(r);
    mpz_clear(r);

void powersum(IN char* stra, IN char* strb)
  ALIAS:
    faulhaber_sum = 1
    jordan_totient = 2
  PREINIT:
    mpz_t a, b;
  PPCODE:
    validate_and_set(a, IFLAG_NONNEG);
    validate_and_set(b, IFLAG_NONNEG);
    if (ix == 0 || ix == 1) {
      if (!mpz_fits_ulong_p(b)) croak("%s: power too large", SUBNAME);
      faulhaber_sum(a, a, mpz_get_ui(b));
    } else {
      if (!mpz_fits_ulong_p(a)) croak("%s: power too large", SUBNAME);
      jordan_totient(a, b, mpz_get_ui(a));
    }
    XPUSH_MPZ(a);
    mpz_clear(b); mpz_clear(a);

void
lshiftint(IN char* strn, IN long k = 1)
  ALIAS:
    rshiftint = 1
    rashiftint = 2
  PREINIT:
    mpz_t n;
    int nix;
  PPCODE:
    validate_and_set(n, IFLAG_ANY);
    nix = ix;
    if (k < 0) {
      k = -k;
      nix = !nix;  /* left => right, right or arith_right => left */
    }
    switch (nix) {
      case 0:  mpz_mul_2exp(n, n, k);      break;
      case 1:  mpz_tdiv_q_2exp(n, n, k);   break;
      case 2:
      default: mpz_fdiv_q_2exp(n, n, k);   break;
    }
    XPUSH_MPZ(n);
    mpz_clear(n);

void
powerful_count(IN char* strn, IN int k = 2)
  ALIAS:
    powerfree_count = 1
  PREINIT:
    mpz_t n, r;
  PPCODE:
    validate_and_set(n, IFLAG_ANY);
    mpz_init(r);
    switch (ix) {
      case 0:  powerful_count(r, n, (unsigned long) k);  break;
      case 1:  powerfree_count(r, n, (uint32_t) k);  break;
      default: break;
    }
    XPUSH_MPZ(r);
    mpz_clear(r);
    mpz_clear(n);

void
addmod(IN char* stra, IN char* strb, IN char* strn)
  ALIAS:
    submod = 1
    mulmod = 2
    powmod = 3
    divmod = 4
    rootmod = 5
  PREINIT:
    mpz_t a, b, n;
    int retundef;
  PPCODE:
    validate_and_set(a, IFLAG_ANY);
    validate_and_set(b, IFLAG_ANY);
    validate_and_set(n, IFLAG_ABS);
    retundef = (mpz_sgn(n) <= 0);
    if (!retundef && ix == 4) {
      if (mpz_cmp_ui(n,1) > 0) {  /* if n is 1, let the mod turn it into zero */
        mpz_mod(b, b, n);         /* Get b between 0 and n-1. */
        if (mpz_sgn(b) == 0)           retundef = 1;
        else if (mpz_cmp_ui(b,1) > 0)  retundef = !mpz_invert(b,b,n);
      }
    }
    if (!retundef && (ix == 3 || ix == 5) && mpz_sgn(b) < 0) {
      if (!mpz_cmp_ui(n,1))       mpz_set_ui(b,0);
      else                        retundef = !mpz_invert(a,a,n);
      mpz_abs(b,b);
    }
    if (retundef) {
      mpz_clear(n); mpz_clear(b); mpz_clear(a);
      XSRETURN_UNDEF;
    }
    if (ix == 0) {
      mpz_add(a,a,b);
      mpz_mod(a,a,n);
    } else if (ix == 1) {
      mpz_sub(a,a,b);
      mpz_mod(a,a,n);
    } else if (ix == 2 || ix == 4) {
      mpz_mul(a,a,b);
      mpz_mod(a,a,n);
    } else if (ix == 3) {
      mpz_powm(a, a, b, n);
    } else if (ix == 5) {
      retundef = !rootmod(a, a, b, n);
    }
    if (retundef) {
      mpz_clear(n); mpz_clear(b); mpz_clear(a);
      XSRETURN_UNDEF;
    }
    XPUSH_MPZ(a);
    mpz_clear(n); mpz_clear(b); mpz_clear(a);

void allsqrtmod(IN char* stra, IN char* strn)
  PREINIT:
    mpz_t a, n;
    mpz_t *roots;
    UV i, nroots;
  PPCODE:
    validate_and_set(a, IFLAG_ANY);
    validate_and_set(n, IFLAG_ABS);
    if (mpz_sgn(n) == 0) {
      mpz_clear(n); mpz_clear(a);
      if (GIMME_V == G_ARRAY) XSRETURN_EMPTY;
      XSRETURN_UV(0);
    }
    if (GIMME_V != G_ARRAY) {
      nroots = allsqrtmod_count(a, n);
      mpz_clear(n); mpz_clear(a);
      XSRETURN_UV(nroots);
    }
    roots = allsqrtmod(&nroots, a, n);
    EXTEND(SP, (long)nroots);
    for (i = 0; i < nroots; i++)
      XPUSH_MPZ(roots[i]);
    clear_rootmod_list(roots, nroots);
    mpz_clear(n); mpz_clear(a);

void allrootmod(IN char* stra, IN char* strk, IN char* strn)
  PREINIT:
    mpz_t a, k, n;
    mpz_t *roots;
    UV i, nroots;
  PPCODE:
    validate_and_set(a, IFLAG_ANY);
    validate_and_set(k, IFLAG_ANY);
    validate_and_set(n, IFLAG_ABS);
    if (mpz_sgn(n) == 0) {
      mpz_clear(n); mpz_clear(k); mpz_clear(a);
      if (GIMME_V == G_ARRAY) XSRETURN_EMPTY;
      XSRETURN_UV(0);
    }
    if (GIMME_V != G_ARRAY) {
      nroots = allrootmod_count(a, k, n);
      mpz_clear(n); mpz_clear(k); mpz_clear(a);
      XSRETURN_UV(nroots);
    }
    roots = allrootmod(&nroots, a, k, n);
    EXTEND(SP, (long)nroots);
    for (i = 0; i < nroots; i++)
      XPUSH_MPZ(roots[i]);
    clear_rootmod_list(roots, nroots);
    mpz_clear(n); mpz_clear(k); mpz_clear(a);

void muladdmod(IN char* stra, IN char* strb, IN char* strc, IN char* strn)
  ALIAS:
    mulsubmod = 1
  PREINIT:
    mpz_t a, b, c, n;
  PPCODE:
    validate_and_set(a, IFLAG_ANY);
    validate_and_set(b, IFLAG_ANY);
    validate_and_set(c, IFLAG_ANY);
    validate_and_set(n, IFLAG_ABS);
    if (mpz_sgn(n) <= 0) {
      mpz_clear(n); mpz_clear(c); mpz_clear(b); mpz_clear(a);
      XSRETURN_UNDEF;
    }
    mpz_mul(a,a,b);
    if (ix == 0)  mpz_add(a, a, c);
    else          mpz_sub(a, a, c);
    mpz_mod(a,a,n);
    XPUSH_MPZ(a);
    mpz_clear(n); mpz_clear(c); mpz_clear(b); mpz_clear(a);


int is_mersenne_prime(IN UV n)
  CODE:
    RETVAL = lucas_lehmer(n);
  OUTPUT:
    RETVAL

void Pi(IN UV n)
  ALIAS:
    Euler = 1
    random_bytes = 2
  PREINIT:
    UV prec;
  PPCODE:
    if (ix == 2) {  /* random_bytes */
      char* sptr;
      SV* sv = newSV(n == 0 ? 1 : n);
      SvPOK_only(sv);
      SvCUR_set(sv, n);
      sptr = SvPVX(sv);
      isaac_rand_bytes(n, (unsigned char*)sptr);
      sptr[n] = '\0';
      PUSHs(sv_2mortal(sv));
      XSRETURN(1);
    }
    if (ix == 0 && n == 0) XSRETURN(0);
    if (ix == 0 && n == 1) XSRETURN_IV(3);
    if (ix == 1 && n == 0) XSRETURN_IV(1);
    prec = (ix == 0) ? n+1 : n+2;
    if (GIMME_V == G_VOID) {
      mpf_t c;
      mpf_init2(c, 7+prec*3.32193);
      if (ix == 0)  const_pi(c, prec);
      else          const_euler(c, prec);
      mpf_clear(c);
    } else {
      char* cstr = (ix == 0) ? piconst(n) : eulerconst(n);
      XPUSHs(sv_2mortal(newSVpvn(cstr, prec)));
      Safefree(cstr);
    }

void random_nbit_prime(IN char* strN)
  ALIAS:
    random_safe_prime = 1
    random_strong_prime = 2
    random_maurer_prime = 3
    random_maurer_prime_with_cert = 4
    random_shawe_taylor_prime = 5
    random_shawe_taylor_prime_with_cert = 6
    random_ndigit_prime = 7
    urandomb = 8
    factorial = 9
    factorial_sum = 10
    subfactorial = 11
    partitions = 12
    partitionsq = 13
    primorial = 14
    pn_primorial = 15
    consecutive_integer_lcm = 16
  PREINIT:
    mpz_t p, N;
    UV n;
    char* proof;
  PPCODE:
    validate_and_set(N, IFLAG_NONNEG);
    if (!mpz_fits_uv_p(N))
      { mpz_clear(N); croak("%s: argument too large",SUBNAME); }
    n = mpz_get_uv(N);
    mpz_clear(N);
    if (ix == 8 && n <= BITS_PER_WORD) {
      UV v = irand64(n);
      ST(0) = sv_2mortal(newSVuv(v));
      XSRETURN(1);
    }
    mpz_init(p);
    proof = 0;
    switch (ix) {
      case 0:  mpz_random_nbit_prime(p, n); break;
      case 1:  mpz_random_safe_prime(p, n); break;
      case 2:  mpz_random_strong_prime(p, n); break;
      case 3:  mpz_random_maurer_prime(p, n, 0); break;
      case 4:  mpz_random_maurer_prime(p, n, &proof);
               proof = cert_with_header(proof, p);
               break;
      case 5:  mpz_random_shawe_taylor_prime(p, n, 0); break;
      case 6:  mpz_random_shawe_taylor_prime(p, n, &proof);
               proof = cert_with_header(proof, p);
               break;
      case 7:  mpz_random_ndigit_prime(p, n); break;
      case 8:  mpz_isaac_urandomb(p, n); break;
      case 9:  mpz_fac_ui(p, n); break;   /* swing impl in 5.1+, so fast */
      case 10: factorial_sum(p, n); break;
      case 11: subfactorial(p, n); break;
      case 12: partitions(p, n); break;
      case 13: partitionsq(p, n); break;
      case 14: _GMP_primorial(p, n);  break;
      case 15: _GMP_pn_primorial(p, n);  break;
      case 16:
      default: consecutive_integer_lcm(p, n);  break;
    }
    XPUSH_MPZ(p);
    mpz_clear(p);
    if (proof) {
      XPUSHs(sv_2mortal(newSVpv(proof, 0)));
      Safefree(proof);
    }

void
stirling(IN char* strn, IN char* strm, IN char* strtype = 0)
  PREINIT:
    mpz_t n, m, type;
    int stype = 1;
  PPCODE:
    validate_and_set(n, IFLAG_NONNEG);
    validate_and_set(m, IFLAG_NONNEG);
    if (items > 2) {
      validate_and_set(type, IFLAG_ANY);
      stype = mpz_fits_sint_p(type) ? (int)mpz_get_si(type) : 0;
      mpz_clear(type);
    }
    mpz_stirling(n, n, m, stype);
    XPUSH_MPZ(n);
    mpz_clear(m);
    mpz_clear(n);

void chinese(...)
  ALIAS:
    chinese2 = 1
  PROTOTYPE: @
  PREINIT:
    int i, doretval;
    mpz_t* an;
    mpz_t ret, lcm;
  PPCODE:
    if (items == 0) {
      if (ix == 0)  XSRETURN_IV(0);
      XPUSH_UINT(0);
      XPUSH_UINT(0);
      XSRETURN(2);
    }
    mpz_init_set_ui(ret, 0);
    New(0, an, 2*items, mpz_t);
    for (i = 0; i < items; i++) {
      AV* av;
      SV** psva;
      SV** psvn;
      char* strn;
      if (!SvROK(ST(i)) || SvTYPE(SvRV(ST(i))) != SVt_PVAV || av_len((AV*)SvRV(ST(i))) != 1)
        croak("chinese arguments are two-element array references");
      av = (AV*) SvRV(ST(i));
      psva = av_fetch(av, 0, 0);
      psvn = av_fetch(av, 1, 0);

      strn = SvPV_nolen(*psva);
      set_integer_string(an[i+0], "a", strn, IFLAG_ANY);

      strn = SvPV_nolen(*psvn);
      set_integer_string(an[i+items], "b", strn, IFLAG_ANY);
    }
    mpz_init(lcm);
    doretval = chinese(ret, lcm, an, an+items, items);
    if (doretval) {
      XPUSH_MPZ(ret);
      if (ix == 1) XPUSH_MPZ(lcm);
    }
    for (i = 0; i < items; i++) {
      mpz_clear(an[i+0]);
      mpz_clear(an[i+items]);
    }
    Safefree(an);
    mpz_clear(lcm);
    mpz_clear(ret);
    if (!doretval) {
      if (ix == 0) {
        XSRETURN_UNDEF;
      } else {
        XPUSHs(&PL_sv_undef);
        XPUSHs(&PL_sv_undef);
        XSRETURN(2);
      }
    }

void
permtonum(SV* svp)
  PREINIT:
    AV *av;
    char* seen;
    UV val, *V;
    int plen, n, i, j, k;
    mpz_t f, t, num;
  PPCODE:
    if ((!SvROK(svp)) || (SvTYPE(SvRV(svp)) != SVt_PVAV))
      croak("permtonum argument must be an array reference");
    av = (AV*) SvRV(svp);
    plen = av_len(av);
    if (plen < 0) XSRETURN_IV(0);
    Newz(0, seen, plen+1, char);
    New(0, V, plen+1, UV);
    for (i = 0; i <= plen; i++) {
      SV **iv = av_fetch(av, i, 0);
      if (iv == 0) break;
      val = SvUV(*iv);
      if (val > (UV)plen || seen[val] != 0) break;
      seen[val] = 1;
      V[i] = val;
    }
    Safefree(seen);
    if (i <= plen)
      croak("permtonum invalid permutation array");

    mpz_init(f);  mpz_init(t);
    mpz_init_set_ui(num, 0);
    n = plen+1;
    mpz_fac_ui(f, n-1);
    for (i = 0; i < n-1; i++) {
      for (j = i+1, k = 0; j < n; j++)
        if (V[j] < V[i])
          k++;
      mpz_mul_ui(t, f, k);
      mpz_add(num, num, t);
      mpz_divexact_ui(f, f, n-i-1);
    }
    Safefree(V);
    XPUSH_MPZ(num);
    mpz_clear(num);  mpz_clear(t);  mpz_clear(f);

void numtoperm(IN UV n, IN char* strk)
  PREINIT:
    mpz_t k, f, p;
    UV i, j, tv, *perm;
  PPCODE:
    if (n == 0)
      XSRETURN_EMPTY;
    validate_and_set(k, IFLAG_ANY);
    mpz_init(f);  mpz_init(p);
    New(0, perm, n, UV);
    for (i = 0; i < n; i++)
      perm[i] = i;
    mpz_fac_ui(f, n);
    mpz_mod(k,k,f);
    for (i = 0; i < n-1; i++) {
      mpz_divexact_ui(f, f, n-i);
      mpz_tdiv_qr(p, k, k, f);
      if (mpz_sgn(p)) {
        for (j = i + mpz_get_ui(p), tv = perm[j]; j > i; j--)
          perm[j] = perm[j-1];
        perm[i] = tv;
      }
    }
    EXTEND(SP, (IV)n);
    for (i = 0; i < n; i++)
      PUSHs(sv_2mortal(newSVuv( perm[i] )));
    Safefree(perm);
    mpz_clear(p);  mpz_clear(f); mpz_clear(k);

void
sieve_prime_cluster(IN char* strlow, IN char* strhigh, ...)
  ALIAS:
    sieve_primes = 1
    sieve_twin_primes = 2
  PREINIT:
    mpz_t low, seghigh, high, t;
    UV i, nc, nprimes, maxseg, *list;
  PPCODE:
    validate_and_set(low, IFLAG_NONNEG);
    validate_and_set(high, IFLAG_NONNEG);
    mpz_init(seghigh);
    mpz_init_set_ui(t,0);

    nc = items-1;
    maxseg = ((UV_MAX > ULONG_MAX) ? ULONG_MAX : UV_MAX);

    /* Loop as needed */
    while (mpz_cmp(low, high) <= 0) {
      mpz_add_ui(seghigh, low, maxseg - 1);
      if (mpz_cmp(seghigh, high) > 0)
        mpz_set(seghigh, high);
      if (ix == 1) {
        UV k = (nc <= 1) ? 0 : SvUV(ST(2));
        list = sieve_primes(low, seghigh, k, &nprimes);
      } else if (ix == 2) {
        list = sieve_twin_primes(low, seghigh, 2, &nprimes);
      } else {
        uint32_t *cl;
        New(0, cl, nc, uint32_t);
        cl[0] = 0;
        for (i = 1; i < nc; i++) {
          UV cval = SvUV(ST(1+i));
          if (cval & 1) croak("sieve_prime_cluster: values must be even");
          if (cval > 2147483647UL) croak("sieve_prime_cluster: values must be 31-bit");
          if (cval <= cl[i-1]) croak("sieve_prime_cluster: values must be increasing");
          cl[i] = cval;
        }
        list = sieve_cluster(low, seghigh, cl, nc, &nprimes);
        Safefree(cl);
      }

      if (list != 0) {
        if (GIMME_V != G_ARRAY) {
          mpz_add_uv(t, t, nprimes);
        } else {
          for (i = 0; i < nprimes; i++) {
            mpz_add_ui(t, low, list[i]);
            XPUSH_MPZ( t );
          }
        }
        Safefree(list);
      }
      mpz_add_ui(low, seghigh, 1);
    }
    if (GIMME_V != G_ARRAY)
      XPUSH_MPZ(t);
    mpz_clear(t);
    mpz_clear(seghigh);
    mpz_clear(high);
    mpz_clear(low);

void
primes(IN SV* svlo, IN SV* svhi = 0)
  ALIAS:
    twin_primes = 1
  PREINIT:
    AV* av;
    mpz_t low, seghigh, high, t;
    UV i, nprimes, maxseg, *list;
  PPCODE:
    if (!SvOK(svlo) || (items >= 2 && !SvOK(svhi)))
      croak("Parameter must be defined");

    if (items == 1) {
      set_integer_string(high, "hi", SvPV_nolen(svlo), IFLAG_NONNEG);
      mpz_init_set_ui(low, 2);
    } else {
      set_integer_string(low, "lo", SvPV_nolen(svlo), IFLAG_NONNEG);
      set_integer_string(high, "hi", SvPV_nolen(svhi), IFLAG_NONNEG);
    }

    av = newAV();
    mpz_init(seghigh);
    mpz_init(t);
    maxseg = ((UV_MAX > ULONG_MAX) ? ULONG_MAX : UV_MAX);
    if (ix == 1) maxseg -= 2;  /* sieve_twin_primes uses length + 2 */

    while (mpz_cmp(low, high) <= 0) {
      mpz_add_ui(seghigh, low, maxseg - 1);
      if (mpz_cmp(seghigh, high) > 0)
        mpz_set(seghigh, high);

      if (ix == 0) list = sieve_primes(low, seghigh, 0, &nprimes);
      else         list = sieve_twin_primes(low, seghigh, 2, &nprimes);
      if (list != 0) {
        for (i = 0; i < nprimes; i++) {
          mpz_add_ui(t, low, list[i]);
          av_push(av, sv_return_for_mpz(aTHX_ t));
        }
        Safefree(list);
      }
      mpz_add_ui(low, seghigh, 1);
    }

    mpz_clear(t);
    mpz_clear(seghigh);
    mpz_clear(high);
    mpz_clear(low);
    XPUSHs(sv_2mortal(newRV_noinc((SV*)av)));

void
sieve_range(IN char* strlow, IN char* strwidth, IN char* strdepth)
  PREINIT:
    mpz_t low, width, depth, high, seghigh;
    UV udepth, uwidth, i, nprimes, *list, offset = 0;
  PPCODE:
    validate_and_set(low, IFLAG_NONNEG);
    validate_and_set(width, IFLAG_NONNEG);
    validate_and_set(depth, IFLAG_NONNEG);

    if (mpz_sgn(width) == 0) {
      mpz_clear(low); mpz_clear(width); mpz_clear(depth);
      XSRETURN_EMPTY;
    }
    if (!mpz_fits_uv_p(depth)) {
      mpz_clear(low); mpz_clear(width); mpz_clear(depth);
      croak("sieve_range: depth must fit in UV");
    }
    if (!mpz_fits_uv_p(width)) {
      mpz_clear(low); mpz_clear(width); mpz_clear(depth);
      croak("sieve_range: width must fit in UV");
    }
    udepth = mpz_get_uv(depth);    mpz_clear(depth);
    uwidth = mpz_get_uv(width);    mpz_clear(width);

    /* sieve_range returns prime candidates, so values below 2 are skipped. */
    if (mpz_cmp_ui(low,2) < 0)
      offset = 2 - mpz_get_ui(low);

    if (udepth < 2) {  /* Depth 0 or 1 does no divisibility sieving. */
      for (i = offset; i < uwidth; i++)
        XPUSH_UINT(i);
      mpz_clear(low);
      XSRETURN(i - offset);
    }

    mpz_init(high);
    mpz_init(seghigh);

    mpz_add_uv(high, low, uwidth-1);
    if (offset != 0)
      mpz_set_ui(low,2);
    /* Loop processing size-2^32-1 segments from low to high */
    while (mpz_cmp(low, high) <= 0) {
      mpz_add_ui(seghigh, low, 4294967295U - 1);
      if (mpz_cmp(seghigh, high) > 0)
        mpz_set(seghigh, high);
      list = sieve_primes(low, seghigh, udepth, &nprimes);

      if (list != 0) {
        for (i = 0; i < nprimes; i++)
          XPUSH_UINT(offset + list[i]);
        Safefree(list);
      }
      mpz_add_ui(low, seghigh, 1);
      offset += 4294967295U;
    }
    mpz_clear(seghigh);
    mpz_clear(high);
    mpz_clear(low);

void
lucas_sequence(IN char* strn, IN char* strP, IN char* strQ, IN char* strk)
  PREINIT:
    mpz_t U, V, Qk, n, P, Q, k, t;
  PPCODE:
    validate_and_set(n, IFLAG_POS);
    validate_and_set(P, IFLAG_ANY);
    validate_and_set(Q, IFLAG_ANY);
    validate_and_set(k, IFLAG_NONNEG);
    mpz_init(U);  mpz_init(V);  mpz_init(Qk);  mpz_init(t);
    lucasuvmod(U, V, P, Q, k, n, t);
    mpz_powm(Qk, Q, k, n);
    XPUSH_MPZ(U);
    XPUSH_MPZ(V);
    XPUSH_MPZ(Qk);
    mpz_clear(n);  mpz_clear(P);  mpz_clear(Q);  mpz_clear(k);
    mpz_clear(U);  mpz_clear(V);  mpz_clear(Qk);  mpz_clear(t);

void
trial_factor(IN char* strn, ...)
  ALIAS:
    prho_factor = 1
    pbrent_factor = 2
    pminus1_factor = 3
    pplus1_factor = 4
    cheb_factor = 5
    holf_factor = 6
    squfof_factor = 7
    ecm_factor = 8
    qs_factor = 9
  PREINIT:
    mpz_t n, t;
    UV arg1, arg2, uf;
    static const UV default_arg1[] =
      {0,   64000000,64000000,5000000,5000000,0,256000000,100000000,0,  0  };
    /*Trial,Rho,     Brent,   P-1,    P+1,    Cheb, HOLF, SQUFOF,   ECM,QS */
  PPCODE:
    validate_and_set(n, IFLAG_NONNEG);
    {
      int cmpr = mpz_cmp_ui(n,1);
      if (cmpr <= 0) {
        mpz_clear(n);
        if (cmpr < 0) XSRETURN_IV(0);
        XSRETURN_EMPTY;
      }
    }
    arg1 = default_arg1[ix];
    arg2 = 0;
    if (items >= 2) {
      set_integer_string(t, "arg1", SvPV_nolen(ST(1)), IFLAG_NONNEG);
      if (!mpz_fits_uv_p(t))
        { mpz_clear(t); croak("%s: arg1 must fit in UV",SUBNAME); }
      arg1 = mpz_get_uv(t);
      mpz_clear(t);
    }
    if (items >= 3) {
      set_integer_string(t, "arg2", SvPV_nolen(ST(2)), IFLAG_NONNEG);
      if (!mpz_fits_uv_p(t))
        { mpz_clear(t); croak("%s: arg2 must fit in UV",SUBNAME); }
      arg2 = mpz_get_uv(t);
      mpz_clear(t);
    }
    while (mpz_even_p(n)) {
      XPUSH_UINT(2);
      mpz_divexact_ui(n, n, 2);
    }
    while (mpz_divisible_ui_p(n, 3)) {
      XPUSH_UINT(3);
      mpz_divexact_ui(n, n, 3);
    }
    while (mpz_divisible_ui_p(n, 5)) {
      XPUSH_UINT(5);
      mpz_divexact_ui(n, n, 5);
    }
    if (mpz_cmp_ui(n,1) > 0 && !_GMP_is_prob_prime(n)) {
      mpz_t f;
      int success = 0;

      mpz_init(f);
      switch (ix) {
        case 0: if (arg1 == 0) arg1 = 2147483647;
                uf = _GMP_trial_factor(n, 2, arg1);
                mpz_set_uv(f, uf);
                success = (uf > 0);
                break;
        case 1: success = _GMP_prho_factor(n, f, 3, arg1);        break;
        case 2: success = _GMP_pbrent_factor(n, f, 3, arg1);      break;
        case 3: if (arg2 == 0)  arg2 = arg1*10;
                success = _GMP_pminus1_factor(n, f, arg1,arg2);   break;
        case 4: if (arg2 == 0)  arg2 = arg1*10;
                success = _GMP_pplus1_factor(n, f, 0,arg1,arg2);  break;
        case 5: success = _GMP_cheb_factor(n, f, arg1,arg2);      break;
        case 6: success = _GMP_holf_factor(n, f, arg1);           break;
        case 7: success = squfof126(n, f, arg1);                  break;
        case 8: if (arg2 == 0) arg2 = 100;
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
        case 9:
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
        mpz_divexact(n, n, f);
        if (mpz_cmp(f,n) > 0)  /* print smallest first */
          mpz_swap(n, f);
        XPUSH_MPZ(f);
      }
      mpz_clear(f);
    }
    if (mpz_cmp_ui(n,1) > 0)
      XPUSH_MPZ(n);
    mpz_clear(n);

void
factor(IN char* strn)
  PREINIT:
    mpz_t n;
    mpz_t* factors;
    int* exponents;
    int nfactors, i, j, isneg;
  PPCODE:
    isneg = validate_and_set(n, IFLAG_ABS);
    if (GIMME_V != G_VOID) {
      nfactors = factor(n, &factors, &exponents);
      if (GIMME_V == G_SCALAR) {
        for (i = 0, j = 0; i < nfactors; i++)
          j += exponents[i];
        PUSHs(sv_2mortal(newSVuv(j)));
      } else {
        if (isneg)
          XPUSH_INT(-1);
        for (i = 0; i < nfactors; i++) {
          for (j = 0; j < exponents[i]; j++) {
            XPUSH_MPZ(factors[i]);
          }
        }
      }
      clear_factors(nfactors, &factors, &exponents);
    }
    mpz_clear(n);

void divisors(IN char* strn, IN char* strk = 0)
  PREINIT:
    mpz_t n, k;
    mpz_t* divs;
    int ndivisors, i;
  PPCODE:
    validate_and_set(n, IFLAG_NONNEG);
    if (strk == 0) {
      mpz_init_set(k, n);
    } else {
      validate_and_set(k, IFLAG_NONNEG);
    }
    if (GIMME_V == G_VOID) {
      /* Nothing */
    } else if (GIMME_V == G_SCALAR && mpz_cmp(k, n) >= 0) {
      sigma(n, n, 0);
      XPUSH_MPZ(n);
    } else {
      divs = divisor_list(&ndivisors, n, k);
      if (GIMME_V == G_SCALAR) {
        XPUSH_INT( ndivisors );
      } else {
        EXTEND(SP, ndivisors);
        for (i = 0; i < ndivisors; i++)
          XPUSH_MPZ(divs[i]);
      }
      for (i = 0; i < ndivisors; i++)
        mpz_clear(divs[i]);
      Safefree(divs);
    }
    mpz_clear(k);
    mpz_clear(n);

void
sigma(IN char* strn, IN UV k = 1)
  PREINIT:
    mpz_t n;
  PPCODE:
    validate_and_set(n, IFLAG_NONNEG);
    sigma(n, n, k);
    XPUSH_MPZ(n);
    mpz_clear(n);

void
todigits(IN char* strn, unsigned int base=10, int length=-1)
  PREINIT:
    mpz_t n;
    uint32_t d, *digits;
  PPCODE:
    if (base < 2 || base > 0xFFFFFFFFU) croak("invalid base: %u\n", base);
    validate_integer_string("n", strn, IFLAG_ANY);
    if (strn[0] == '-' || strn[0] == '+')  strn++;
    if (base == 10) {
      uint32_t l = strlen(strn);
      New(0, digits, l, uint32_t);
      for (d = 0; d < l; d++)
        digits[d] = strn[d]-'0';
    } else {
      mpz_init_set_str(n, strn, 10);
      digits = todigits(&d, n, base);
      mpz_clear(n);
    }
    if (length > 0 || d > 1 || digits[0] != 0) {
      if (length < 0) length = d;
      EXTEND(SP, length);
      for (; length > (int)d; length--)
        PUSHs(sv_2mortal(newSVuv( 0 )));
      for (; length > 0; length--)
        PUSHs(sv_2mortal(newSVuv( digits[d-length] )));
    }
    Safefree(digits);

void
fromdigits(IN SV* svp, IN SV* svbase = 0)
  PREINIT:
    AV *av;
    const char *ds;
    int i, plen;
    size_t j, len;
    mpz_t n, base, *digits;
  PPCODE:
    if (!SvOK(svp)) croak("Parameter must be defined");
    if (items > 1 && SvOK(svbase)) {
      set_integer_string(base, "base", SvPV_nolen(svbase), IFLAG_NONNEG);
    } else {
      mpz_init_set_ui(base, 10);
    }
    if (mpz_cmp_ui(base, 2) < 0) {
      SV *basesv = sv_2mortal(sv_return_for_mpz(aTHX_ base));
      mpz_clear(base);
      croak("fromdigits: invalid base: %s", SvPV_nolen(basesv));
    }
    mpz_init(n);
    if (SvROK(svp) && SvTYPE(SvRV(svp)) == SVt_PVAV) {
      av = (AV*) SvRV(svp);
      plen = av_len(av);
      if (plen < 0) {
        mpz_set_ui(n, 0);
      } else {
        len = (size_t) plen + 1;
        New(0, digits, len, mpz_t);
        for (j = 0; j < len; j++)
          mpz_init_set_ui(digits[j], 0);
        for (i = 0; i <= plen; i++) {
          SV **iv = av_fetch(av, i, 0);
          if (iv == 0 || !SvOK(*iv)) {
            for (j = 0; j < len; j++)
              mpz_clear(digits[j]);
            Safefree(digits);
            mpz_clear(n);
            mpz_clear(base);
            croak("Parameter must be defined");
          }
          ds = SvPV_nolen(*iv);
          validate_integer_string("digit", ds, IFLAG_ANY);
          mpz_set_str(digits[plen-i], ds + (*ds == '+'), 10);
        }
        mpz_fromdigits(n, digits, len, base);
        for (j = 0; j < len; j++)
          mpz_clear(digits[j]);
        Safefree(digits);
      }
    } else if (!SvROK(svp) || _sv_is_math_object(aTHX_ svp)) { /* string */
      if (!mpz_fromdigits_str(n, SvPV_nolen(svp), base))
        _croak_invalid_fromdigits_digit(aTHX_ base);
    } else {
      mpz_clear(n);
      mpz_clear(base);
      croak("fromdigits: first argument must be a string or array reference");
    }
    XPUSH_MPZ(n);
    mpz_clear(n);
    mpz_clear(base);
