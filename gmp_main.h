#ifndef MPU_GMPMAIN_H
#define MPU_GMPMAIN_H

#include <gmp.h>
#include "ptypes.h"

extern void _GMP_init(void);
extern void _GMP_destroy(void);
extern void _GMP_memfree(void);

extern int  primality_pretest(mpz_t n);
extern int  is_trial_prime(mpz_t n);

extern void _GMP_next_prime(mpz_t n);
extern void _GMP_prev_prime(mpz_t n);
extern void surround_primes(mpz_t n, UV* prev, UV* next, UV skip_width);

extern void _GMP_pn_primorial(mpz_t prim, UV n);
extern void _GMP_primorial(mpz_t prim, UV n);
extern void consecutive_integer_lcm(mpz_t m, unsigned long B);
extern void stirling(mpz_t r, unsigned long n, unsigned long m, UV type);
extern void binomial(mpz_t r, UV n, UV k);
extern void partitions(mpz_t npart, UV n);
extern void factorialmod(mpz_t r, UV n, mpz_t m);
extern void multifactorial(mpz_t r, unsigned long n, unsigned long k);
extern void factorial_sum(mpz_t r, unsigned long n);
extern void subfactorial(mpz_t r, unsigned long n);
extern void rising_factorial(mpz_t r, unsigned long x, unsigned long n);
extern void falling_factorial(mpz_t r, unsigned long x, unsigned long n);

extern void faulhaber_sum(mpz_t sum, mpz_t zn, unsigned long p);

extern void powerful_count(mpz_t r, mpz_t n, unsigned long k);

extern int  is_carmichael(mpz_t n);
extern int  is_fundamental(mpz_t n);
extern int  is_practical(mpz_t n);
extern int  is_totient(mpz_t n);
extern void polygonal_nth(mpz_t r, mpz_t n, mpz_t k);

extern void exp_mangoldt(mpz_t res, mpz_t n);

extern uint32_t* partial_sieve(mpz_t start, UV length, UV maxprime);

extern void prime_count_lower(mpz_t pc, mpz_t n);
extern void prime_count_upper(mpz_t pc, mpz_t n);
extern UV* sieve_primes(mpz_t low, mpz_t high, UV k, UV *rn);
extern UV* sieve_twin_primes(mpz_t low, mpz_t high, UV twin, UV *rn);
extern UV* sieve_cluster(mpz_t low, mpz_t high, uint32_t* cl, UV nc, UV *rn);

extern void prime_count(mpz_t count, mpz_t hi);
extern void prime_count_range(mpz_t count, mpz_t lo, mpz_t hi);

extern void prime_power_count(mpz_t r, mpz_t n);
extern void prime_power_count_range(mpz_t r, mpz_t lo, mpz_t hi);

extern void perfect_power_count(mpz_t r, mpz_t n);
extern void perfect_power_count_range(mpz_t r, mpz_t lo, mpz_t hi);

extern void next_twin_prime(mpz_t res, mpz_t n);

extern uint32_t* todigits(uint32_t *ndigits, mpz_t n, uint32_t base);
extern void fromdigits(mpz_t n, uint32_t *d, uint32_t len, uint32_t base);
extern void fromdigits_str(mpz_t n, const char* s, uint32_t base);

#endif
