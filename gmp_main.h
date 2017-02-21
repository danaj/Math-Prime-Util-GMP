#ifndef MPU_GMPMAIN_H
#define MPU_GMPMAIN_H

#include <gmp.h>
#include "ptypes.h"

extern void _GMP_init(void);
extern void _GMP_destroy(void);

extern int  primality_pretest(mpz_t n);

extern void _GMP_next_prime(mpz_t n);
extern void _GMP_prev_prime(mpz_t n);
extern void surround_primes(mpz_t n, UV* prev, UV* next, UV skip_width);

extern void _GMP_pn_primorial(mpz_t prim, UV n);
extern void _GMP_primorial(mpz_t prim, UV n);
extern void _GMP_lcm_of_consecutive_integers(UV B, mpz_t m);
extern void bernfrac(mpz_t num, mpz_t den, mpz_t n);
extern void harmfrac(mpz_t num, mpz_t den, mpz_t n);
extern void stirling(mpz_t r, unsigned long n, unsigned long m, UV type);
extern void binomial(mpz_t r, UV n, UV k);
extern void partitions(mpz_t npart, UV n);

extern UV   is_power(mpz_t n, UV a);

extern void exp_mangoldt(mpz_t res, mpz_t n);

extern uint32_t* partial_sieve(mpz_t start, UV length, UV maxprime);
extern char* pidigits(UV n);
extern char* bernreal(mpz_t zn, unsigned long prec);
extern char* harmreal(mpz_t zn, unsigned long prec);
extern char* zetareal(mpf_t r, unsigned long prec);
extern char* riemannrreal(mpf_t r, unsigned long prec);
extern char* lambertwreal(mpf_t r, unsigned long prec);

extern UV* sieve_primes(mpz_t low, mpz_t high, UV k, UV *rn);
extern UV* sieve_twin_primes(mpz_t low, mpz_t high, UV twin, UV *rn);
extern UV* sieve_cluster(mpz_t low, mpz_t high, uint32_t* cl, UV nc, UV *rn);

#endif
