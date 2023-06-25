#ifndef MPU_BLS75_H
#define MPU_BLS75_H

#include <gmp.h>
#include "ptypes.h"

/* Return values of 0 do not indicate composite with certainty.
 * A zero result means that we were unable to construct a proof.
 */

/* These check the theorem conditions for given n and a factor. */
/* This is NOT a full proof as the factor isn't verified. */
/* Check BLS75 theorem  3 conditions */
extern int BLS_check_T3(mpz_t n, mpz_t p, UV* a);
/* Check BLS75 theorem 15 conditions */
extern int BLS_check_T15(mpz_t n, mpz_t q, IV* lp, IV* lq);


/* These construct a complete recursive proof. */

/* BLS75 theorem 5/7 complete proof */
extern int BLS_primality_nm1(mpz_t n, int effort, char ** prooftextptr);

/* BLS75 theorem 17 complete proof (N+1) */
extern int BLS_primality_np1(mpz_t n, int effort, char** prooftextptr);

/* BLS75 theorem 20 complete proof (N-1 and N+1) */
extern int BLS_primality(mpz_t n, int effort, char** prooftextptr);

#endif
