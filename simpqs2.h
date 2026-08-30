#ifndef MPU_SIMPQS2_H
#define MPU_SIMPQS2_H

#include <gmp.h>
#include "ptypes.h"

/* trial_start is the first candidate not already checked for small factors. */
extern mpz_t *_GMP_simpqs2(const mpz_t n, uint32_t *nfactors,
                            uint32_t trial_start);
extern void _GMP_simpqs2_free(mpz_t *factors, uint32_t nfactors);

#endif
