#ifndef MPU_POWERFREE_H
#define MPU_POWERFREE_H

#include <gmp.h>
#include "ptypes.h"

extern int is_powerfree(const mpz_t n, uint32_t k);

extern void next_powerfree(mpz_t next, const mpz_t n, uint32_t k);
extern void prev_powerfree(mpz_t prev, const mpz_t n, uint32_t k);

extern void powerfree_count(mpz_t r, const mpz_t n, uint32_t k);
extern void powerfree_count_range(mpz_t r, const mpz_t lo, const mpz_t hi, uint32_t k);

extern void nth_powerfree(mpz_t nth, const mpz_t n, uint32_t k);

#endif
