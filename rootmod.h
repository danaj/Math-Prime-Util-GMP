#ifndef MPU_GMPROOTMOD_H
#define MPU_GMPROOTMOD_H

#include <gmp.h>
#include "ptypes.h"

extern int sqrtmodp(mpz_t r, mpz_t a, mpz_t p);   /* sqrt(a) mod p */
extern int sqrtmod (mpz_t r, mpz_t a, mpz_t n);   /* sqrt(a) mod n */

/* No aliasing allowed and pass in 4 temps. */
extern int sqrtmodp_t(mpz_t r,  mpz_t a,  mpz_t p,
                      mpz_t t1, mpz_t t2, mpz_t t3, mpz_t t4);


#if 0
extern int rootmodp(mpz_t r, mpz_t a, mpz_t k, mpz_t p);  /* a^(1/k) mod p */
extern int rootmod(mpz_t r, mpz_t a, mpz_t k, mpz_t n);   /* a^(1/k) mod n */


extern UV allsqrtmod(mpz_t* roots, mpz_t a, mpz_t n);          /* all results */
extern UV allrootmod(mpz_t* roots, mpz_t a, mpz_t k, mpz_t n); /* all results */
#endif

#endif
