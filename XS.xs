
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



MODULE = Math::Prime::Util::GMP		PACKAGE = Math::Prime::Util::GMP

PROTOTYPES: ENABLE

int
miller_rabin(IN char* strn, IN char* strbase)
  PREINIT:
    mpz_t n, a, nminus1, d, x;
    UV s, r;
  CODE:
    validate_string_number("miller_rabin (n)", strn);
    validate_string_number("miller_rabin (base)", strbase);
    if (strn[1] == 0) {
      switch (strn[0]) {
        case '2': case '3': case '5': case '7': XSRETURN_IV(1); break;
        case '0': case '1': case '4': case '6': case '8': XSRETURN_IV(0); break;
        default:  break; /* let 9 fall through */
      }
    }
    mpz_init_set_str(n, strn, 10);
    /* gmp_printf("Computing MR base %Zd on %Zd\n", a, n); */
    if (mpz_even_p(n)) {                     /* multiple of 2 */
      mpz_clear(n); XSRETURN_IV(0);
    }
    mpz_init_set_str(a, strbase, 10);
    mpz_init_set(nminus1, n);
    mpz_sub_ui(nminus1, nminus1, 1);
    mpz_init_set(d, nminus1);
    if (1) {
      s = 0;
      while (mpz_even_p(d)) {
        s++;
        mpz_divexact_ui(d, d, 2);
      }
    } else {
      /* faster way, verify s is identical */
      s = mpz_scan1(d, 0);
      mpz_tdiv_q_2exp(d, d, s);
    }
    mpz_init(x);
    mpz_powm(x, a, d, n);
    mpz_clear(a);  mpz_clear(d); /* done with a and d */
    RETVAL = 0;
    if (!mpz_cmp_ui(x, 1) || !mpz_cmp(x, nminus1)) {
      RETVAL = 1;
    } else {
      for (r = 0; r < s; r++) {
        mpz_powm_ui(x, x, 2, n);
        if (!mpz_cmp_ui(x, 1)) {
          break;
        }
        if (!mpz_cmp(x, nminus1)) {
          RETVAL = 1;
          break;
        }
      }
    }
    mpz_clear(n); mpz_clear(nminus1); mpz_clear(x);
  OUTPUT:
    RETVAL

int
is_strong_lucas_pseudoprime(IN char* strn)
  PREINIT:
    mpz_t n, d, U, V, Qkd;
    IV D;
    UV P = 1;
    IV Q;
    UV s;
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
    if (mpz_even_p(n) || mpz_perfect_square_p(n)) { /* n even or perfect sq */
      mpz_clear(n);
      XSRETURN_IV(0);
    }
    /* Determine Selfridge D, P, Q parameters */
    {
      mpz_t t;
      UV D_ui = 5;
      IV sign = 1;
      mpz_init(t);
      while (1) {
        UV gcd, j;
        gcd = mpz_gcd_ui(NULL, n, D_ui);
        if ((gcd > 1) && mpz_cmp_ui(n, gcd) != 0) {
          D_ui = 0;
          break;
        }
        mpz_set_si(t, (IV)D_ui * sign);
        j = mpz_jacobi(t, n);
        if (j == -1)  break;
        D_ui += 2;
        sign = -sign;
      }
      mpz_clear(t);
      if (D_ui == 0) {
        mpz_clear(n);
        XSRETURN_IV(0);
      }
      D = (IV)D_ui * sign;
    }
    Q = (1 - D) / 4;
    //gmp_printf("N: %Zd  D: %ld  P: %lu  Q: %ld\n", n, D, P, Q);
    if (D != P*P - 4*Q)  croak("incorrect DPQ\n");
    /* Now start on the lucas sequence */
    mpz_init_set(d, n);
    mpz_add_ui(d, d, 1);
    s = 0;
    while (mpz_even_p(d)) {
      s++;
      mpz_divexact_ui(d, d, 2);
    }
    mpz_init_set_ui(U, 1);
    mpz_init_set_ui(V, P);
    {
      mpz_t U2m, V2m, Qm, T1, T2;
      mpz_init_set(U2m, U);
      mpz_init_set(V2m, V);
      mpz_init_set_si(Qm, Q);
      mpz_init_set(Qkd, Qm);
      mpz_tdiv_q_ui(d, d, 2);
      mpz_init(T1);
      mpz_init(T2);
      while (mpz_sgn(d) > 0) {
        //gmp_printf("U=%Zd  V=%Zd  Qm=%Zd\n", U, V, Qm);
        mpz_mul(U2m, U2m, V2m);
        mpz_mod(U2m, U2m, n);
        mpz_powm_ui(V2m, V2m, 2, n);
        mpz_submul_ui(V2m, Qm, 2);
        mpz_mod(V2m, V2m, n);
        //gmp_printf("  l  U2m=%Zd  V2m=%Zd\n", U2m, V2m);
        mpz_powm_ui(Qm, Qm, 2, n);
        if (mpz_odd_p(d)) {
          mpz_mul(T1, U2m, V);
          mpz_mul(T2, U2m, U);
          mpz_mul_si(T2, T2, D);
          //gmp_printf("      T1 %Zd  T2 %Zd\n", T1, T2);
          /* U */
          mpz_mul(U, U, V2m);
          mpz_add(U, U, T1);
          if (mpz_odd_p(U)) mpz_add(U, U, n);
          mpz_fdiv_q_ui(U, U, 2);
          mpz_mod(U, U, n);
          /* V */
          mpz_mul(V, V, V2m);
          mpz_add(V, V, T2);
          if (mpz_odd_p(V)) mpz_add(V, V, n);
          mpz_fdiv_q_ui(V, V, 2);
          mpz_mod(V, V, n);
          /* Qkd */
          mpz_mul(Qkd, Qkd, Qm);
          mpz_mod(Qkd, Qkd, n);
        }
        mpz_tdiv_q_ui(d, d, 2);
      }
      mpz_clear(U2m); mpz_clear(V2m); mpz_clear(Qm); mpz_clear(T1); mpz_clear(T2);
    }
    //gmp_printf("l0 U=%Zd  V=%Zd\n", U, V);
    RETVAL = 0;
    if ( (mpz_sgn(U) == 0) || (mpz_sgn(V) == 0) ) {
      RETVAL = 1;
      s = 0;
    }
    /* Powers of V */
    while (s--) {
      mpz_mul(V, V, V);
      mpz_submul_ui(V, Qkd, 2);
      mpz_mod(V, V, n);
      if (mpz_sgn(V) == 0) {
        RETVAL = 1;
        break;
      }
      if (s) {
        mpz_powm_ui(Qkd, Qkd, 2, n);
      }
    }
    mpz_clear(n); mpz_clear(d); mpz_clear(U); mpz_clear(V); mpz_clear(Qkd);
  OUTPUT:
    RETVAL

int
_GMP_trial_div(IN char* strn, UV to_n)
  PREINIT:
    mpz_t n;
  CODE:
    /* If this function returns 0, it is a composite.
     * If this function returns 1, it is a prime if n <= to_n*to_n.
     */
    validate_string_number("_GMP_trial_div (n)", strn);
    if (strn[1] == 0) {
      int q_is_prime = 0;
      switch (strn[0]) {
        case '2': case '3': case '5': case '7': q_is_prime = 1; break;
      }
      XSRETURN_IV(q_is_prime);
    }
    /* n is >= 10 */
    mpz_init_set_str(n, strn, 10);
    if (mpz_even_p(n)) {
      mpz_clear(n); XSRETURN_IV(0);
    }
    {
      int small_n = 0;
      int primei = 2;
      UV f = primes_small[primei];
      if (mpz_cmp_ui(n, to_n*to_n) < 0)
        small_n = 1;
      while (f <= to_n) {
        if (f > to_n) {
          mpz_clear(n); XSRETURN_IV(1);
        }
        if (small_n && mpz_cmp_ui(n, f*f) < 0) {
          mpz_clear(n); XSRETURN_IV(1);
        }
        //gmp_printf("  is %Zd divisible by %lu?\n", n, primes_small[primei]);
        if (mpz_divisible_ui_p(n, primes_small[primei])) {
          mpz_clear(n); XSRETURN_IV(0);
        }
        if (++primei < NPRIMES_SMALL) {
          f = primes_small[primei];
        } else {
          do {
            f += 2;
          } while ( !(f % 3) || !(f % 5) || !(f % 7) );
        }
      }
    }
    mpz_clear(n);
    RETVAL = 1;
  OUTPUT:
    RETVAL
