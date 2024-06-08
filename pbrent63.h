#ifndef MPU_PBRENT63_H
#define MPU_PBRENT63_H

#include <gmp.h>
#include "ptypes.h"

extern int pbrent63(const mpz_t n, mpz_t f, UV rounds);

extern int uvpbrent63(UV n, UV *factors, UV rounds, UV a);

#endif
