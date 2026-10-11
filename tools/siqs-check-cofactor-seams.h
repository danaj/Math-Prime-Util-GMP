/* Test-only bounded-splitter misses. Include after the host adapters and
 * before siqs.c; undefine the interposition immediately afterward. Flags
 * stay immutable while any concurrent cofactor checks are running. */
#include "../squfof126.h"
static int cofactor_miss_squfof;
static int cofactor_miss_siqs, cofactor_trace;
static unsigned int cofactor_squfof_calls, cofactor_prime_calls;
static int check_cofactor_squfof(const mpz_t n, mpz_t f, UV rounds) {
  if (cofactor_trace) cofactor_squfof_calls++;
  return cofactor_miss_squfof ? 0 : squfof126(n, f, rounds);
}
/* A false "prime" answer makes the inner SIQS return its unsplit partition.
 * The resolver's real outer primality/acceptance checks are not interposed. */
static int check_cofactor_prime(const mpz_t n) {
  if (cofactor_trace) cofactor_prime_calls++;
  return cofactor_miss_siqs ? 1 : siqs_is_prob_prime(n);
}
#define squfof126 check_cofactor_squfof
#define siqs_is_prob_prime check_cofactor_prime
