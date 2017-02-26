#ifndef MPU_RANDOM_PRIME_H
#define MPU_RANDOM_PRIME_H

#include "ptypes.h"

extern int mpz_random_prime(mpz_t p, mpz_t lo, mpz_t hi);
extern void mpz_random_nbit_prime(mpz_t p, UV n);
extern void mpz_random_ndigit_prime(mpz_t p, UV n);

#endif
