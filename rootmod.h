#ifndef MPU_GMPROOTMOD_H
#define MPU_GMPROOTMOD_H

#include <gmp.h>
#include "ptypes.h"

extern int sqrtmod(mpz_t r, const mpz_t a, const mpz_t n);  /* sqrt(a) mod n */
extern int rootmod(mpz_t r, const mpz_t a, const mpz_t k, const mpz_t n);  /* a^(1/k) mod n */

extern mpz_t* allsqrtmod(UV* nroots, const mpz_t a, const mpz_t n);
extern mpz_t* allrootmod(UV* nroots, const mpz_t a, const mpz_t k, const mpz_t n); /* all results */

extern UV allsqrtmod_count(const mpz_t a, const mpz_t n);
extern UV allrootmod_count(const mpz_t a, const mpz_t k, const mpz_t n);

extern void clear_rootmod_list(mpz_t* roots, UV nroots);


/* Special case when we know the modulus is prime */
extern int sqrtmodp(mpz_t r, const mpz_t a, const mpz_t p);  /* sqrt(a) mod p */
/* No aliasing allowed and pass in 4 temps. */
extern int sqrtmodp_t(mpz_t r,  const mpz_t a,  const mpz_t p,
                      mpz_t t1, mpz_t t2, mpz_t t3, mpz_t t4);


#endif
