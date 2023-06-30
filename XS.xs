
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
/* We're not using anything for which we need ppport.h */

#include <string.h>
#include <gmp.h>

#include "ptypes.h"
#include "gmp_main.h"
#include "primality.h"
#include "lucas_seq.h"
#include "squfof126.h"
#include "ecm.h"
#include "simpqs.h"
#include "bls75.h"
#include "ecpp.h"
#include "aks.h"
#include "rootmod.h"
#include "utility.h"
#include "factor.h"
#include "isaac.h"
#include "random_prime.h"
#include "real.h"
#define _GMP_ECM_FACTOR(n, f, b1, ncurves) \
   _GMP_ecm_factor_projective(n, f, b1, 0, ncurves)


/* Instead of trying to suck in lots of Math::BigInt::GMP and be terribly
 * clever (and brittle), just do all C<->Perl bigints via strings.  It's
 * crude but seems to work pretty well.
 */

static void validate_string_number(CV* cv, const char* var, const char* s)
{
  const char* p;
  if (s == 0)
    croak("%s (%s): null string pointer as input", GvNAME(CvGV(cv)),var);
  if (*s == 0)
    croak("%s (%s): empty string as input", GvNAME(CvGV(cv)),var);
  p = s;
  while (*p != 0) {
    if (!isdigit(*p))
      croak("%s (%s): input '%s' must be a positive integer", GvNAME(CvGV(cv)), var, s);
    p++;
  }
}

#define VALIDATE_AND_SET(var, str) \
  do { \
    const char* s = str; \
    if (*s == '+') s++; \
    validate_string_number(cv,#var,s); \
    mpz_init_set_str(var, s, 10); \
  } while (0)

#define VSETNEG_ERR 0
#define VSETNEG_POS 1
#define VSETNEG_OK  2
static int validate_and_set_signed(CV* cv, mpz_t v, const char* vname, const char* s, int negflag) {
  int neg = (s && *s == '-');
  if (s && *s == '+') s++;
  validate_string_number(cv, vname, (neg && negflag != VSETNEG_ERR) ? s+1 : s);
  mpz_init_set_str(v, (neg && negflag == VSETNEG_POS) ? s+1 : s, 10);
  return neg;
}

static char* cert_with_header(char* proof, mpz_t n) {
  char *str, *strptr;
  if (proof == 0) {
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

static SV* sv_return_for_mpz(const mpz_t n) {
  SV* sv = 0;
  UV v = mpz_get_ui(n);
  if (!mpz_cmp_ui(n, v)) {  /* Try to use a scalar */
    sv = newSVuv(v);
  } else {
    char* str;
    int nsize = mpz_sizeinbase(n, 10) + 2;
    New(0, str, nsize, char);
    mpz_get_str(str, 10, n);
    sv = newSVpv(str, 0);
    Safefree(str);
  }
  return sv;
}

#define XPUSH_MPZ(n) \
  XPUSHs(sv_2mortal( sv_return_for_mpz(n) ))


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
    int i;
    mpz_t n, a;
  CODE:
    if (items < 2) croak("%s: no bases", GvNAME(CvGV(cv)));
    validate_string_number(cv,"n",strn);
    for (i = 1; i < items; i++) {
      const char* strbase = SvPV_nolen(ST(i));
      validate_string_number(cv, "base", strbase);
      if (strbase[1] == '\0' && (strbase[0] == '0' || strbase[0] == '1'))
        croak("Base %s is invalid", strbase);
    }
    mpz_init_set_str(n, strn, 10);
    for (i = 1; i < items; i++) {
      mpz_init_set_str(a, SvPV_nolen(ST(i)), 10);
      switch (ix) {
        case 0:  RETVAL = is_pseudoprime(n, a); break;
        case 1:  RETVAL = is_euler_pseudoprime(n, a); break;
        case 2:
        default: RETVAL = miller_rabin(n, a); break;
      }
      mpz_clear(a);
      if (!RETVAL) break;
    }
    mpz_clear(n);
  OUTPUT:
    RETVAL

int miller_rabin_random(IN char* strn, IN IV nbases = 1, IN char* seedstr = 0)
  PREINIT:
    mpz_t n;
  CODE:
    if (nbases < 0)
      croak("Parameter '%"IVdf"' must be a positive integer\n", nbases);
    VALIDATE_AND_SET(n, strn);
    RETVAL = miller_rabin_random(n, nbases, seedstr);
    mpz_clear(n);
  OUTPUT:
    RETVAL

#define PRIMALITY_START(name, small_retval, test_small_factors) \
    /* Negative numbers return 0 */ \
    if ((strn != 0) && (strn[0] == '-') ) \
      XSRETURN_IV(0); \
    validate_string_number(cv, "n", strn); \
    /* Handle single digit numbers */ \
    if (strn[1] == 0) { \
      int q_is_prime = 0; \
      switch (strn[0]) { \
        case '2': case '3': case '5': case '7': q_is_prime = (small_retval); \
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
    is_frobenius_khashin_pseudoprime = 4
    is_euler_plumb_pseudoprime = 5
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


void
_is_provable_prime(IN char* strn, IN int wantproof = 0)
  ALIAS:
    is_miller_prime = 1
    is_perrin_pseudoprime = 2
  PREINIT:
    int result;
    mpz_t n;
  PPCODE:
    PRIMALITY_START("is_provable_prime", 2, ix != 2);
    if (ix == 1) {
      result = is_miller_prime(n, wantproof);  /* Assume GRH or not */
      mpz_clear(n);
      XSRETURN_IV(result);
    }
    if (ix == 2) {
      result = is_perrin_pseudoprime(n, wantproof);  /* Restricted or not */
      mpz_clear(n);
      XSRETURN_IV(result);
    }
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
    mpz_t a, b, n, px, py, m, q, t1, t2;
  CODE:
    VALIDATE_AND_SET(a, stra);
    VALIDATE_AND_SET(b, strb);  /* Unused */
    VALIDATE_AND_SET(n, strn);
    VALIDATE_AND_SET(px, strpx);
    VALIDATE_AND_SET(py, strpy);
    VALIDATE_AND_SET(m, strm);
    VALIDATE_AND_SET(q, strq);
    mpz_init(t1);  mpz_init(t2);
    RETVAL = (ecpp_check_point(px, py, m, q, a, n, t1, t2) == 2) ? 1 : 0;
    mpz_clear(t1); mpz_clear(t2);
    mpz_clear(a); mpz_clear(b); mpz_clear(n); mpz_clear(px); mpz_clear(py);
    mpz_clear(m); mpz_clear(q);
  OUTPUT:
    RETVAL

int
is_almost_prime(IN unsigned int k, IN char* strn)
  PREINIT:
    mpz_t n;
  CODE:
    VALIDATE_AND_SET(n, strn);
    RETVAL = is_almost_prime(k, n);
    mpz_clear(n);
  OUTPUT:
    RETVAL

UV is_power(IN char* strn, IN UV a = 0)
  PREINIT:
    mpz_t n;
    int isneg;
  CODE:
    isneg = validate_and_set_signed(cv, n, "n", strn, VSETNEG_POS);
    RETVAL = 0;
    if (!isneg || (a == 0 || a & 1)) {
      RETVAL = is_power(n, a);
    }
    if (isneg && a == 0 && RETVAL != 0) {
      UV r = RETVAL;
      while (!(r & 1)) r >>= 1;
      RETVAL = (r == 1) ? 0 : r;
    }
    mpz_clear(n);
  OUTPUT:
    RETVAL

int is_divisible(IN char* strn, IN char* strd)
  PREINIT:
    mpz_t n, d;
  CODE:
    validate_and_set_signed(cv, n, "n", strn, VSETNEG_OK);
    validate_and_set_signed(cv, d, "d", strd, VSETNEG_OK);
    RETVAL = !!mpz_divisible_p(n, d);
    mpz_clear(d);  mpz_clear(n);
  OUTPUT:
    RETVAL

int is_congruent(IN char* strn, IN char* strc, IN char* strd)
  PREINIT:
    mpz_t n, c, d;
  CODE:
    validate_and_set_signed(cv, n, "n", strn, VSETNEG_OK);
    validate_and_set_signed(cv, c, "c", strc, VSETNEG_OK);
    validate_and_set_signed(cv, d, "d", strd, VSETNEG_OK);
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
    VALIDATE_AND_SET(n, strn);
    if (ix == 1 && mpz_cmp_ui(n,3) < 0) { mpz_clear(n); XSRETURN_UNDEF; }
    if      (ix == 0) _GMP_next_prime(n);
    else if (ix == 1) _GMP_prev_prime(n);
    else              next_twin_prime(n, n);
    XPUSH_MPZ(n);
    mpz_clear(n);


void
random_prime(IN char* strlo, IN char* strhi = 0)
  ALIAS:
    urandomr = 1
  PREINIT:
    mpz_t lo, hi, res;
    int retundef;
  PPCODE:
    if (strhi == 0) {
      mpz_init_set_ui(lo, 0);
      VALIDATE_AND_SET(hi, strlo);
    } else {
      if (*strlo == '+')  strlo++;
      if (*strhi == '+')  strhi++;
      validate_string_number(cv, "lo", strlo);
      validate_string_number(cv, "hi", strhi);
      mpz_init_set_str(lo, strlo, 10);
      mpz_init_set_str(hi, strhi, 10);
    }
    if (ix == 1 && mpz_sizeinbase(hi,2) <= 32) {
      uint32_t ulo = mpz_get_ui(lo),  uhi = mpz_get_ui(hi);
      if (ulo <= uhi) {
        mpz_clear(lo); mpz_clear(hi);
        XSRETURN_IV( ulo + isaac_rand(uhi-ulo+1) );
      }
    }
    retundef = 0;
    mpz_init(res);
    if (ix == 0) {
      retundef = !mpz_random_prime(res, lo, hi);
    } else {
      if (mpz_cmp(lo,hi) > 0) {
        retundef = 1;
      } else {
        mpz_sub(hi,hi,lo);
        mpz_add_ui(hi,hi,1);
        mpz_isaac_urandomm(res, hi);
        mpz_add(res,res,lo);
      }
    }
    if (!retundef) XPUSH_MPZ(res);
    mpz_clear(res);
    mpz_clear(hi);
    mpz_clear(lo);
    if (retundef) XSRETURN_UNDEF;

void prime_count(IN char* strlo, IN char* strhi = 0)
  ALIAS:
    prime_power_count = 1
    perfect_power_count = 2
  PREINIT:
    mpz_t n, lo, hi, res;
  PPCODE:
    mpz_init(res);
    if (strhi == 0) {
      VALIDATE_AND_SET(n, strlo);
      switch (ix) {
        case  0: prime_count(res, n); break;
        case  1: prime_power_count(res, n); break;
        case  2: perfect_power_count(res, n); break;
        default: break;
      }
      mpz_clear(n);
    } else {
      VALIDATE_AND_SET(lo, strlo);
      VALIDATE_AND_SET(hi, strhi);
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

void
totient(IN char* strn)
  ALIAS:
    carmichael_lambda = 1
    exp_mangoldt = 2
    bernfrac = 3
    harmfrac = 4
    znprimroot = 5
    ramanujan_tau = 6
    sqrtint = 7
    is_prime_power = 8
    prime_count_lower = 9
    prime_count_upper = 10
    urandomm = 13
    add1int = 14
    sub1int = 15
  PREINIT:
    mpz_t res, n;
  PPCODE:
    if (strn != 0 && strn[0] == '-') { /* If input is negative... */
      if (ix == 2)  XSRETURN_IV(1);    /* exp_mangoldt return 1 */
      if (ix == 5)  strn++;            /* znprimroot flip sign */
      if (ix == 8)  XSRETURN_IV(0);    /* is_prime_power return 0 */
    }
    validate_and_set_signed(cv, n, "n", strn,
                            (ix==14 || ix==15) ? VSETNEG_OK : VSETNEG_ERR);
    mpz_init(res);
    switch (ix) {
      case 0:  totient(res, n);  break;
      case 1:  carmichael_lambda(res, n);  break;
      case 2:  exp_mangoldt(res, n);  break;
      case 3:  bernfrac(n, res, n);
               XPUSH_MPZ(n);
               break;
      case 4:  harmfrac(n, res, n);
               XPUSH_MPZ(n);
               break;
      case 5:  znprimroot(res, n);
               if (!mpz_sgn(res) && mpz_cmp_ui(n,1) != 0) {
                 mpz_clear(n);  mpz_clear(res);
                 XSRETURN_UNDEF;
               }
               break;
      case 6:  ramanujan_tau(res, n);  break;
      case 7:  mpz_sqrt(res, n);  break;
      case 8:  mpz_set_uv(res, prime_power(res, n)); break;
      case 9:  prime_count_lower(res, n); break;
      case 10: prime_count_upper(res, n); break;
      case 13: mpz_isaac_urandomm(res, n); break;
      case 14: mpz_add_ui(res, n, 1); break;
      case 15: mpz_sub_ui(res, n, 1); break;
      default: break;
    }
    XPUSH_MPZ(res);
    mpz_clear(n);
    mpz_clear(res);

void absint(IN char* strn)
  ALIAS:
    negint = 1
  PREINIT:
    mpz_t res;
  PPCODE:
    if (strn != 0 && strn[0] == '-') {
      VALIDATE_AND_SET(res, strn+1);
    } else {
      VALIDATE_AND_SET(res, strn);
      if (ix == 1) mpz_neg(res, res);
    }
    XPUSH_MPZ(res);
    mpz_clear(res);

void signint(IN char* strn)
  PREINIT:
    mpz_t n;
    int res;
  PPCODE:
    VALIDATE_AND_SET(n, strn);
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
    VALIDATE_AND_SET(a, stra);
    VALIDATE_AND_SET(b, strb);
    res = (ix == 0) ? mpz_cmp(a, b) : mpz_cmpabs(a, b);
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
    VALIDATE_AND_SET(n, strn);
    switch (ix) {
      case 0:  mpz_setbit(n, k); break;
      case 1:  mpz_clrbit(n, k); break;
      case 2:  mpz_combit(n, k); break;
      case 3:  res = mpz_tstbit(n, k); break;
      default: break;
    }
    if (ix != 3) XPUSH_MPZ(n);
    mpz_clear(n);
    if (ix == 3) XSRETURN_UV(res);

void bitand(IN char* stra, IN char* strb)
  ALIAS:
    bitor = 1
    bitxor = 2
  PREINIT:
    mpz_t a, b;
  PPCODE:
    VALIDATE_AND_SET(a, stra);
    VALIDATE_AND_SET(b, strb);
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
  PREINIT:
    mpz_t n;
  PPCODE:
    VALIDATE_AND_SET(n, strn);
    mpz_com(n, n);
    XPUSH_MPZ(n);
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
      VALIDATE_AND_SET(n, strn);
      next = 1 + (mpz_sgn(n)==0);
      if (mpz_cmp_ui(n,2) > 0) {
        surround_primes(n, &prev, &next, (items == 1) ? 0 : prec);
        XPUSHs(sv_2mortal(newSVuv(prev)));
      } else {
        XPUSHs(sv_2mortal(newSV(0)));
      }
      XPUSHs(sv_2mortal(newSVuv(next)));
      mpz_clear(n);
    } else if (ix <= 1) {
      VALIDATE_AND_SET(n, strn);
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
      EXTEND(SP, n+1);
      for (i = 0; i <= n; i++) {
        AV* av = newAV();
        av_push(av, sv_return_for_mpz(N[i]));
        av_push(av, sv_return_for_mpz(D[i]));
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
    negflag = (ix <= 1) ? VSETNEG_POS : VSETNEG_OK;
    if (ix == 1 || ix == 3) {
      mpz_t* list;
      New(0, list, items, mpz_t);
      for (i = 0; i < items; i++) {
        char* strn = SvPV_nolen(ST(i));
        validate_and_set_signed(cv, list[i], "arg", strn, negflag);
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
      validate_and_set_signed(cv, n, "arg", strn, negflag);
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
    validate_and_set_signed(cv, a, "a", stra, VSETNEG_OK);
    validate_and_set_signed(cv, b, "b", strb, VSETNEG_OK);
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

void
moebius(IN char* strn, IN char* stro = 0)
  PREINIT:
    mpz_t n;
  PPCODE:
    validate_and_set_signed(cv, n, "n", strn, VSETNEG_OK);
    if (stro == 0) {
      int result = moebius(n);
      mpz_clear(n);
      XSRETURN_IV(result);
    } else {   /* Ranged result */
      mpz_t nhi;
      validate_and_set_signed(cv, nhi, "nhi", stro, VSETNEG_OK);
      while (mpz_cmp(n, nhi) <= 0) {
        XPUSHs(sv_2mortal(newSViv( moebius(n) )));
        mpz_add_ui(n, n, 1);
      }
      mpz_clear(n);
      mpz_clear(nhi);
    }

void lucasu(IN char* strp, IN char* strq, IN char* strk)
  ALIAS:
    lucasv = 1
    lucasuv = 2
  PREINIT:
    mpz_t u, v, p, q, k;
  PPCODE:
    validate_and_set_signed(cv, p, "P", strp, VSETNEG_OK);
    validate_and_set_signed(cv, q, "Q", strq, VSETNEG_OK);
    VALIDATE_AND_SET(k, strk);
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
    validate_and_set_signed(cv, p, "P", strp, VSETNEG_OK);
    validate_and_set_signed(cv, q, "Q", strq, VSETNEG_OK);
    VALIDATE_AND_SET(k, strk);
    VALIDATE_AND_SET(n, strn);
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

int
liouville(IN char* strn)
  ALIAS:
    is_square = 1
    is_semiprime = 2
    is_totient = 3
    is_carmichael = 4
    is_practical = 5
    is_fundamental = 6
    hammingweight = 7
    prime_omega = 8
    prime_bigomega = 9
  PREINIT:
    mpz_t n;
    int isneg;
  CODE:
    isneg = validate_and_set_signed( cv, n, "n", strn,
                                     (ix == 0) ? VSETNEG_ERR
                                   : (ix  < 7) ? VSETNEG_OK
                                   :             VSETNEG_POS );
    if (isneg && (ix >= 1 && ix <= 5)) {
      RETVAL = 0;
    } else {
      switch (ix) {
        case 0:  RETVAL = liouville(n);      break;
        case 1:  RETVAL = is_power(n,2);     break;
        case 2:  RETVAL = is_semiprime(n);   break;
        case 3:  RETVAL = is_totient(n);     break;
        case 4:  RETVAL = is_carmichael(n);  break;
        case 5:  RETVAL = is_practical(n);   break;
        case 6:  RETVAL = is_fundamental(n); break;
        case 7:  RETVAL = mpz_popcount(n);   break;
        case 8:  RETVAL = omega(n);          break;
        case 9:
        default: RETVAL = bigomega(n);       break;
      }
    }
    mpz_clear(n);
  OUTPUT:
    RETVAL

int
is_powerful(IN char* strn, IN UV k = 0)
  PREINIT:
    mpz_t n;
  CODE:
    VALIDATE_AND_SET(n, strn);
    RETVAL = is_powerful(n, (k == 0) ? 2 : k);
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
    sqrtmod = 5
    is_primitive_root = 6
    is_polygonal = 7
    polygonal_nth = 8
    rootint = 9
    logint = 10
    powint = 11
    mulint = 12
    addint = 13
    subint = 14
    divint = 15
    modint = 16
    cdivint = 17
    tdivrem = 18
    fdivrem = 19
    cdivrem = 20
    divrem = 21
    factorialmod = 22
    multifactorial = 23
  PREINIT:
    mpz_t a, b, t;
    int retundef;
  PPCODE:
    validate_and_set_signed(cv, a, "a", stra, VSETNEG_OK);
    validate_and_set_signed(cv, b, "b", strb, VSETNEG_OK);
    retundef = 0;
    switch (ix) {
               /* undef if a|b = 0, 0 if b is 1, else result of mpz_invert */
      case 0:{ if (!mpz_sgn(b) || !mpz_sgn(a))  retundef = 1;
               else if (!mpz_cmp_ui(b,1))       mpz_set_ui(a,0);
               else                             retundef = !mpz_invert(a,a,b);
             } break;
      case 1:{ unsigned long n, k;
               if (mpz_sgn(b) < 0) {   /* Handle negative k */
                 if (mpz_sgn(a) >= 0 || mpz_cmp(b,a) > 0)  mpz_set_ui(a, 0);
                 else                                      mpz_sub(b, a, b);
               }
               n = mpz_get_ui(a);
               k = mpz_get_ui(b);
               if (k > n || k == 0 || k == n || mpz_sgn(a) < 0) {
                 mpz_bin_ui(a, a, k);
               } else {
                 if (k > n/2) k = n-k;
                 /* Note: mpz_bin_uiui is *much* faster than mpz_bin_ui.
                  * It is a bit faster than our code for small values, and
                  * a tiny bit slower for larger values. */
                 if (n > 50000 && k > 1000)  binomial(a, n, k);
                 else                        mpz_bin_uiui(a, n, k);
               }
             } break;
      case 2:{ mpz_init(t);
               if (mpz_sgn(a) == 0 && mpz_sgn(b) == 0) {
                 mpz_set_ui(t, 0);  /* This changed in GMP 5.1.2.  Enforce new result. */
               } else {
                 gcdext(a, t, b, a, b);
               }
               XPUSH_MPZ(t);  XPUSH_MPZ(b);
               mpz_clear(t);
             } break;
      case 3: jordan_totient(a, b, mpz_get_ui(a));
              break;
      case 4: znorder(a, a, b);
              if (!mpz_sgn(a)) retundef = 1;
              break;
      case 5: mpz_abs(b, b);
              retundef = !sqrtmod(a, a, b);
              break;
      case 6: mpz_set_si(a, is_primitive_root(a, b, 0) );
              break;
      case 7: if (mpz_cmp_ui(b,3) < 0) croak("is_polygonal: k must be >= 3");
              polygonal_nth(a, a, b);
              mpz_set_si(a, mpz_sgn(a));
              break;
      case 8: if (mpz_cmp_ui(b,3) < 0) croak("polygonal_nth: k must be >= 3");
              polygonal_nth(a, a, b);
              break;
      case 9: if (mpz_sgn(b) <= 0) croak("rootint: k must be > 0");
              if (mpz_sgn(a) <  0) croak("rootint: n must be >= 0");
              mpz_root(a, a, mpz_get_ui(b));
              break;
      case 10:if (mpz_cmp_ui(b,2) < 0) croak("rootint: base must be > 1");
              if (mpz_sgn(a) <=  0) croak("rootint: n must be > 0");
              mpz_set_uv(a, logint(a, mpz_get_uv(b)));
              break;
      case 11:if (mpz_sgn(b) < 0) croak("powint: exponent must be >= 0");
              mpz_pow_ui(a, a, mpz_get_ui(b));
              break;
      case 12:mpz_mul(a, a, b);
              break;
      case 13:mpz_add(a, a, b);
              break;
      case 14:mpz_sub(a, a, b);
              break;
      case 15:if (mpz_sgn(b) == 0) croak("divint: divide by zero");
              mpz_fdiv_q(a, a, b);
              break;
      case 16:if (mpz_sgn(b) == 0) croak("modint: divide by zero");
              mpz_fdiv_r(a, a, b);
              break;
      case 17:if (mpz_sgn(b) == 0) croak("cdivint: divide by zero");
              mpz_cdiv_q(a, a, b);
              break;
      case 18:if (mpz_sgn(b) == 0) croak("tdivrem: divide by zero");
              mpz_tdiv_qr(b, a, a, b);  /* t is t-quotient, a is t-remainder */
              XPUSH_MPZ(b);
              break;
      case 19:if (mpz_sgn(b) == 0) croak("fdivrem: divide by zero");
              mpz_fdiv_qr(b, a, a, b);  /* t is f-quotient, a is f-remainder */
              XPUSH_MPZ(b);
              break;
      case 20:if (mpz_sgn(b) == 0) croak("cdivrem: divide by zero");
              mpz_cdiv_qr(b, a, a, b);  /* t is c-quotient, a is c-remainder */
              XPUSH_MPZ(b);
              break;
      case 21:mpz_init_set(t, b);
              if (mpz_sgn(b) == 0) croak("divrem: divide by zero");
              mpz_tdiv_qr(t, a, a, b);  /* t is t-quotient, a is t-remainder */
              if (mpz_sgn(a) < 0) {  /* Change from trunc to Euclidean */
                if (mpz_sgn(b) > 0) {
                  mpz_sub_ui(t, t, 1);
                  mpz_add(a, a, b);
                } else {
                  mpz_add_ui(t, t, 1);
                  mpz_sub(a, a, b);
                }
              }
              XPUSH_MPZ(t);
              mpz_clear(t);
              break;
      case 22:if (mpz_sgn(b) < 0) retundef = 1;
              else                factorialmod(a, mpz_get_ui(a), b);
              break;
      case 23:
      default:if (mpz_sgn(a) < 0 || mpz_sgn(b) < 0) retundef = 1;
              else                multifactorial(a, mpz_get_ui(a), mpz_get_ui(b));
              break;

    }
    if (!retundef) XPUSH_MPZ(a);
    mpz_clear(b); mpz_clear(a);
    if (retundef) XSRETURN_UNDEF;

int is_qr(IN char* stra, IN char* strn)
  PREINIT:
    mpz_t a, n;
    int retval;
  PPCODE:
    validate_and_set_signed(cv, a, "a", stra, VSETNEG_OK);
    validate_and_set_signed(cv, n, "n", strn, VSETNEG_OK);
    retval = -1;
    if (mpz_sgn(n) != 0) {
      mpz_abs(n,n);
      retval = sqrtmod(a,a,n);
    }
    mpz_clear(n); mpz_clear(a);
    if (retval == -1)
      XSRETURN_UNDEF;
    else
      XSRETURN_IV(retval);

void powersum(IN char* stra, IN char* strb)
  ALIAS:
    faulhaber_sum = 1
  PREINIT:
    mpz_t a, b;
  PPCODE:
    validate_and_set_signed(cv, a, "a", stra, VSETNEG_OK);
    validate_and_set_signed(cv, b, "b", strb, VSETNEG_OK);
    if (mpz_sgn(a) < 0 || mpz_sgn(b) < 0) croak("powersum: negative argument");
    if (ix == 0 || ix == 1)
      faulhaber_sum(a, a, mpz_get_ui(b));
    XPUSH_MPZ(a);
    mpz_clear(b); mpz_clear(a);

void
lshiftint(IN char* strn, IN unsigned long k = 1)
  ALIAS:
    rshiftint = 1
    rashiftint = 2
  PREINIT:
    mpz_t n;
  PPCODE:
    validate_and_set_signed(cv, n, "n", strn, VSETNEG_OK);
    switch (ix) {
      case 0:  mpz_mul_2exp(n, n, k);      break;
      case 1:  mpz_tdiv_q_2exp(n, n, k);   break;
      case 2:
      default: mpz_fdiv_q_2exp(n, n, k);   break;
    }
    XPUSH_MPZ(n);
    mpz_clear(n);

void
powerful_count(IN char* strn, IN int k = 2)
  PREINIT:
    mpz_t n, r;
  PPCODE:
    validate_and_set_signed(cv, n, "n", strn, VSETNEG_ERR);
    mpz_init(r);
    powerful_count(r, n, (unsigned long) k);
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
  PREINIT:
    mpz_t a, b, n;
    int retundef;
  PPCODE:
    validate_and_set_signed(cv, a, "a", stra, VSETNEG_OK);
    validate_and_set_signed(cv, b, "b", strb, VSETNEG_OK);
    validate_and_set_signed(cv, n, "n", strn, VSETNEG_ERR);
    retundef = (mpz_sgn(n) <= 0);
    if (!retundef && ix == 4) {
      if (mpz_cmp_ui(n,1) > 0) {  /* if n is 1, let the mod turn it into zero */
        mpz_mod(b, b, n);         /* Get b between 0 and n-1. */
        if (mpz_sgn(b) == 0)           retundef = 1;
        else if (mpz_cmp_ui(b,1) > 0)  retundef = !mpz_invert(b,b,n);
      }
    }
    if (!retundef && ix == 3 && mpz_sgn(b) < 0) {
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
    }
    XPUSH_MPZ(a);
    mpz_clear(n); mpz_clear(b); mpz_clear(a);

void muladdmod(IN char* stra, IN char* strb, IN char* strc, IN char* strn)
  ALIAS:
    mulsubmod = 1
  PREINIT:
    mpz_t a, b, c, n;
  PPCODE:
    validate_and_set_signed(cv, a, "a", stra, VSETNEG_OK);
    validate_and_set_signed(cv, b, "b", strb, VSETNEG_OK);
    validate_and_set_signed(cv, c, "c", strc, VSETNEG_OK);
    validate_and_set_signed(cv, n, "n", strn, VSETNEG_ERR);
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

void random_nbit_prime(IN UV n)
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
    primorial = 13
    pn_primorial = 14
    consecutive_integer_lcm = 15
  PREINIT:
    mpz_t p;
    char* proof;
  PPCODE:
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
      case 13: _GMP_primorial(p, n);  break;
      case 14: _GMP_pn_primorial(p, n);  break;
      case 15:
      default: consecutive_integer_lcm(p, n);  break;
    }
    XPUSH_MPZ(p);
    mpz_clear(p);
    if (proof) {
      XPUSHs(sv_2mortal(newSVpv(proof, 0)));
      Safefree(proof);
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
      if (ix == 0)  XSRETURN_UV(0);
      XPUSHs(sv_2mortal(newSVuv( 0 )));
      XPUSHs(sv_2mortal(newSVuv( 0 )));
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
      validate_and_set_signed(cv, an[i+0], "a", strn, VSETNEG_OK);

      strn = SvPV_nolen(*psvn);
      validate_and_set_signed(cv, an[i+items], "b", strn, VSETNEG_OK);
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
    validate_and_set_signed(cv, k, "k", strk, VSETNEG_OK);
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
    VALIDATE_AND_SET(low, strlow);
    VALIDATE_AND_SET(high, strhigh);
    mpz_init(seghigh);
    mpz_init(t);

    nc = items-1;
    maxseg = ((UV_MAX > ULONG_MAX) ? ULONG_MAX : UV_MAX);

    /* Loop as needed */
    while (mpz_cmp(low, high) <= 0) {
      mpz_add_ui(seghigh, low, maxseg - 1);
      if (mpz_cmp(seghigh, high) > 0)
        mpz_set(seghigh, high);
      mpz_set(t, seghigh);  /* Save in case it is modified */
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
      mpz_set(seghigh, t);  /* Restore the value we used */

      if (list != 0) {
        for (i = 0; i < nprimes; i++) {
          mpz_add_ui(t, low, list[i]);
          XPUSH_MPZ( t );
        }
        Safefree(list);
      }
      mpz_add_ui(low, seghigh, 1);
    }
    mpz_clear(t);
    mpz_clear(seghigh);
    mpz_clear(high);
    mpz_clear(low);

void
sieve_range(IN char* strn, IN UV width, IN UV depth)
  PREINIT:
    mpz_t low, seghigh, high, t;
    UV i, nprimes, maxseg, offset, *list;
  PPCODE:
    if (width == 0) XSRETURN(0);
    if (depth == 0) depth = 1;

    VALIDATE_AND_SET(low, strn);
    mpz_init(high);
    mpz_add_ui(high, low, width-1);
    mpz_init(seghigh);
    mpz_init(t);
    maxseg = ((UV_MAX > ULONG_MAX) ? ULONG_MAX : UV_MAX);
    offset = 0;

    /* Deal with 0 and 1 inside range */
    if (mpz_cmp_ui(low,2) < 0) {
      offset = 2 - mpz_get_ui(low);
      width = (width < offset) ? 0 : width - offset;
      mpz_set_ui(low,2);
    }
    /* Deal with depth < 2 (no sieving) */
    if (depth < 2) {
      for (i = 0; i < width; i++)
        XPUSHs(sv_2mortal(newSVuv(offset + i)));
      mpz_add_ui(low, high, 1);
    }
    /* Loop as needed */
    while (mpz_cmp(low, high) <= 0) {
      mpz_add_ui(seghigh, low, maxseg - 1);
      if (mpz_cmp(seghigh, high) > 0)
        mpz_set(seghigh, high);
      mpz_set(t, seghigh);  /* Save in case it is modified */
      list = sieve_primes(low, seghigh, depth, &nprimes);
      mpz_set(seghigh, t);  /* Restore the value we used */

      if (list != 0) {
        for (i = 0; i < nprimes; i++) {
          XPUSHs(sv_2mortal(newSVuv( offset + list[i] )));
        }
        Safefree(list);
      }
      mpz_add_ui(low, seghigh, 1);
      offset += maxseg;
    }
    mpz_clear(t);
    mpz_clear(seghigh);
    mpz_clear(high);
    mpz_clear(low);

void
lucas_sequence(IN char* strn, IN IV P, IN IV Q, IN char* strk)
  PREINIT:
    mpz_t U, V, Qk, n, k, t;
  PPCODE:
    VALIDATE_AND_SET(n, strn);
    VALIDATE_AND_SET(k, strk);
    mpz_init(U);  mpz_init(V);  mpz_init(Qk);  mpz_init(t);
    lucas_seq(U, V, n, P, Q, k, Qk, t);
    XPUSH_MPZ(U);
    XPUSH_MPZ(V);
    XPUSH_MPZ(Qk);

    mpz_clear(n);  mpz_clear(k);
    mpz_clear(U);  mpz_clear(V);  mpz_clear(Qk);  mpz_clear(t);


#define SET_UV_VIA_MPZ_STRING(uva, sva, name) \
  { \
      mpz_t t; \
      char* stra = SvPV_nolen(sva); \
      validate_string_number(cv, name, stra); \
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
    cheb_factor = 5
    holf_factor = 6
    squfof_factor = 7
    ecm_factor = 8
    qs_factor = 9
  PREINIT:
    mpz_t n;
    UV arg1, arg2, uf;
    static const UV default_arg1[] =
      {0,   64000000,64000000,5000000,5000000,0,256000000,100000000,0,  0  };
    /*Trial,Rho,     Brent,   P-1,    P+1,    Cheb, HOLF, SQUFOF,   ECM,QS */
  PPCODE:
    VALIDATE_AND_SET(n, strn);
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
  ALIAS:
    divisors = 1
  PREINIT:
    mpz_t n;
    mpz_t* factors;
    int* exponents;
    int nfactors, i, j;
  PPCODE:
    VALIDATE_AND_SET(n, strn);
    if (ix == 0) {
      nfactors = factor(n, &factors, &exponents);
      if (GIMME_V == G_SCALAR) {
        for (i = 0, j = 0; i < nfactors; i++)
          j += exponents[i];
        PUSHs(sv_2mortal(newSVuv(j)));
      } else {
        for (i = 0; i < nfactors; i++) {
          for (j = 0; j < exponents[i]; j++) {
            XPUSH_MPZ(factors[i]);
          }
        }
      }
      clear_factors(nfactors, &factors, &exponents);
    } else {
      if (GIMME_V == G_SCALAR) {
        sigma(n, n, 0);
        XPUSH_MPZ(n);
      } else {
        factors = divisor_list(&nfactors, n);
        EXTEND(SP, nfactors);
        for (i = 0; i < nfactors; i++) {
          XPUSH_MPZ(factors[i]);
          mpz_clear(factors[i]);
        }
        Safefree(factors);
      }
    }
    mpz_clear(n);

void
sigma(IN char* strn, IN UV k = 1)
  PREINIT:
    mpz_t n;
  PPCODE:
    VALIDATE_AND_SET(n, strn);
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
    if (strn[0] == '-' || strn[0] == '+')  strn++;
    validate_string_number(cv, "n", strn);
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
fromdigits(IN SV* svp, unsigned int base=10)
  PREINIT:
    AV *av;
    int i, plen;
    uint32_t *digits;
    mpz_t n;
  PPCODE:
    if (base < 2 || base > 0xFFFFFFFFU) croak("invalid base: %u\n", base);
    mpz_init(n);
    if (!SvROK(svp)) { /* string */
      fromdigits_str(n, SvPV_nolen(svp), base);
    } else {
      if (SvTYPE(SvRV(svp)) != SVt_PVAV)
        croak("fromdigits argument must be a string or array reference");
      av = (AV*) SvRV(svp);
      plen = av_len(av);
      if (plen < 0) XSRETURN_IV(0);
      New(0, digits, plen+1, uint32_t);
      for (i = 0; i <= plen; i++) {
        SV **iv = av_fetch(av, i, 0);
        if (iv == 0) break;
        digits[plen-i] = SvUV(*iv); /* TODO: anything other than 32-bit */
      }
      if (i >= plen)
        fromdigits(n, digits, plen+1, base);
      Safefree(digits);
    }
    XPUSH_MPZ(n);
    mpz_clear(n);
