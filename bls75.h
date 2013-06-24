#ifndef MPU_BLS75_H
#define MPU_BLS75_H

#include <gmp.h>
#include "ptypes.h"

/* extern int _GMP_primality_pocklington(mpz_t n, int do_quick); */

/* Note that the theorem 3 and 15 checks, as well as the splitters:
 *  1) do not do a full proof.  You must verify q.
 *  2) do not indicate compositeness on failure.
 */

/* These will check the theorem conditions for given n and factor. */
/* Check BLS75 theorem  3 conditions */
extern int _GMP_primality_bls_3(mpz_t n, mpz_t p, UV* a);
/* Check BLS75 theorem 15 conditions */
extern int _GMP_primality_bls_15(mpz_t n, mpz_t q, IV* lp, IV* lq);

/* These will try to factor and check the theorem conditions. */
/* BLS75 theorem  3, you must verify q for a proof */
extern int _GMP_primality_bls_nm1_split(mpz_t n, int effort, mpz_t q, UV* a);
/* BLS75 theorem 15, you must verify q for a proof */
extern int _GMP_primality_bls_np1_split(mpz_t n, int effort, mpz_t q, IV* lp, IV* lq);

/* This does a complete recursive proof */
/* BLS75 theorem 5/7 complete proof */
extern int _GMP_primality_bls_nm1(mpz_t n, int effort, char ** prooftextptr);

#endif
