#ifndef MPU_BLS75_H
#define MPU_BLS75_H

#include <gmp.h>
#include "ptypes.h"

/* extern int _GMP_primality_pocklington(mpz_t n, int do_quick); */
extern int _GMP_primality_bls_nm1(mpz_t n, int effort, char ** prooftextptr);

#endif
