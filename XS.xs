
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
/* We're not using anything for which we need ppport.h */

#include <gmp.h>

#include "ptypes.h"

/* I think we're going to have to end up sucking in a lot of Math::BigInt::GMP
 * infrastructure. */

MODULE = Math::Prime::Util::GMP		PACKAGE = Math::Prime::Util::GMP

PROTOTYPES: ENABLE

int
miller_rabin(IN char* strn, IN char* strbase)
  PREINIT:
    mpz_t n, a, nminus1, d, x;
    UV s, r;
  CODE:
    mpz_init_set_str(n, strn, 10);
    mpz_init_set_str(a, strbase, 10);
    /* gmp_printf("Computing MR base %Zd on %Zd\n", a, n); */
    if (mpz_cmp_ui(n, 1) <= 0)    XSRETURN_IV(0); /* n <= 1 is composite */
    if (mpz_cmp_ui(n, 3) <= 0)    XSRETURN_IV(1); /* 2 and 3 are prime */
    if (mpz_even_p(n))            XSRETURN_IV(0); /* multiple of 2 */
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
    if (!mpz_cmp_ui(x, 1) || !mpz_cmp(x, nminus1))
      XSRETURN_IV(1);
    for (r = 0; r < s; r++) {
      mpz_powm_ui(x, x, 2, n);
      if (!mpz_cmp_ui(x, 1))
        XSRETURN_IV(0);
      if (!mpz_cmp(x, nminus1))
        XSRETURN_IV(1);
    }
    RETVAL = 0;
  OUTPUT:
    RETVAL
