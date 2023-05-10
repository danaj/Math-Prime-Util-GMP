#ifndef MPU_GMPLUCASSEQ_H
#define MPU_GMPLUCASSEQ_H

#include <gmp.h>
#include "ptypes.h"

extern void lucasuv(mpz_t U, mpz_t V, mpz_t P, mpz_t Q, mpz_t k);

extern void internal_lucas_vmod_q1(mpz_t V, mpz_t W, mpz_t P, mpz_t k, mpz_t n);
extern void lucasuvmod(mpz_t U, mpz_t V, mpz_t P, mpz_t Q, mpz_t k, mpz_t n, mpz_t t);
extern void lucas_seq(mpz_t U, mpz_t V, mpz_t n, IV P, IV Q, mpz_t k,
                      mpz_t Qk, mpz_t t);

extern void lucasumod(mpz_t U, mpz_t P, mpz_t Q, mpz_t k, mpz_t n, mpz_t t);
extern void lucasvmod(mpz_t V, mpz_t P, mpz_t Q, mpz_t k, mpz_t n, mpz_t t);

#endif
