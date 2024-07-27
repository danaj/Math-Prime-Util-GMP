#ifndef MPU_GMPMAIN_H
#define MPU_GMPMAIN_H

#include <gmp.h>
#include "ptypes.h"

extern void _GMP_init(void);
extern void _GMP_destroy(void);
extern void _GMP_memfree(void);

extern int  primality_pretest(const mpz_t n);
extern int  is_trial_prime(const mpz_t n);

extern void _GMP_next_prime(mpz_t n);
extern void _GMP_prev_prime(mpz_t n);
extern void surround_primes(const mpz_t n, UV* prev, UV* next, UV skip_width);

extern void _GMP_pn_primorial(mpz_t prim, UV n);
extern void _GMP_primorial(mpz_t prim, UV n);
extern void consecutive_integer_lcm(mpz_t m, unsigned long B);
extern void stirling(mpz_t r, unsigned long n, unsigned long m, UV type);
extern void binomial(mpz_t r, UV n, UV k);
extern void partitions(mpz_t npart, UV n);
extern void factorialmod(mpz_t r, UV n, const mpz_t m);
extern void multifactorial(mpz_t r, unsigned long n, unsigned long k);
extern void factorial_sum(mpz_t r, unsigned long n);
extern void subfactorial(mpz_t r, unsigned long n);
extern void rising_factorial(mpz_t r, mpz_t x, mpz_t n);
extern void falling_factorial(mpz_t r, mpz_t x, mpz_t n);

extern void faulhaber_sum(mpz_t sum, const mpz_t zn, unsigned long p);

extern void hclassno(mpz_t res, const mpz_t n);
extern void rtau(mpz_t res, const mpz_t n);

extern void powerful_count(mpz_t r, const mpz_t n, unsigned long k);

extern int  is_carmichael(const mpz_t n);
extern int  is_fundamental(const mpz_t n);
extern int  is_practical(const mpz_t n);
extern int  is_totient(const mpz_t n);
extern void polygonal_nth(mpz_t r, const mpz_t n, const mpz_t k);

extern void exp_mangoldt(mpz_t res, const mpz_t n);

extern void prime_count_lower(mpz_t pc, const mpz_t n);
extern void prime_count_upper(mpz_t pc, const mpz_t n);

extern void prime_count(mpz_t count, const mpz_t hi);
extern void prime_count_range(mpz_t count, const mpz_t lo, const mpz_t hi);

extern void prime_power_count(mpz_t r, const mpz_t n);
extern void prime_power_count_range(mpz_t r, const mpz_t lo, const mpz_t hi);

extern void next_twin_prime(mpz_t res, const mpz_t n);

extern uint32_t* todigits(uint32_t *ndigits, const mpz_t n, uint32_t base);
extern void fromdigits(mpz_t n, uint32_t *d, uint32_t len, uint32_t base);
extern void fromdigits_str(mpz_t n, const char* s, uint32_t base);

/* Partial sieve, used by many functions in this file.
 *
 * start must be odd.  It is a mpz_t type and can be very large.
 * length so start+length-1 is the last value checked.  Must be > 0.
 *        length = hi-lo+1.  hi = lo+length-1.
 * maxprime is the maximum prime for sieving.
 *          Reduced to sqrt(start+length) if larger.
 * Returns an array of odd values, where 1 bits indicate composite.
 * Since the array of odds, 2 is always implicitly sieved.
 */

extern uint32_t* partial_sieve(const mpz_t start, UV length, UV maxprime);

/* Sieving for primes.
 * low and high can be any values.
 * k indicates how much sieving should be done before primality testing.
 *   set k=0 to let the function figure something out.
 * rn is the number of primes returned.
 * Return value are the primes as offsets from low.
 */
extern UV* sieve_primes(const mpz_t low, const mpz_t high, UV k, UV *rn);
extern UV* sieve_twin_primes(const mpz_t low, const mpz_t high, UV twin, UV *rn);
extern UV* sieve_cluster(const mpz_t low, const mpz_t high, uint32_t* cl, UV nc, UV *rn);

#endif
