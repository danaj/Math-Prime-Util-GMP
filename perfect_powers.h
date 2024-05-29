#ifndef MPU_PERFECT_POWERS_H
#define MPU_PERFECT_POWERS_H

#include <gmp.h>
#include "ptypes.h"

extern int is_perfect_power(mpz_t n);

extern void next_perfect_power(mpz_t next, mpz_t n);
extern void prev_perfect_power(mpz_t prev, mpz_t n);

extern void perfect_power_count(mpz_t r, mpz_t n);
extern void perfect_power_count_range(mpz_t r, mpz_t lo, mpz_t hi);

extern void nth_perfect_power(mpz_t nth, mpz_t n);
extern void nth_perfect_power_approx(mpz_t nth, mpz_t n);

#endif
