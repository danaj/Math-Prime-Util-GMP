#ifndef MPU_SIQS_DEP_H
#define MPU_SIQS_DEP_H

#include <gmp.h>
#include "ptypes.h"

#ifdef STANDALONE

extern int siqs_verbose_level(void);
extern int siqs_is_prob_prime(const mpz_t n);
extern int siqs_pbrent_factor(const mpz_t n, mpz_t f, UV a, UV rounds);

#else

#include "utility.h"
#include "primality.h"
#include "factor.h"

#define siqs_verbose_level()         get_verbose_level()
#define siqs_is_prob_prime(n)        _GMP_is_prob_prime(n)
#define siqs_pbrent_factor(n,f,a,r)  _GMP_pbrent_factor(n,f,a,r)

#endif

#endif
