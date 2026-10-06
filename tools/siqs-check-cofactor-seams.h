/* Test-only bounded-splitter misses. Include after the host adapters and
 * before siqs.c; undefine the interposition immediately afterward. Flags
 * stay immutable while any concurrent cofactor checks are running. */
#include "../squfof126.h"
#include "../pbrent63.h"
static int cofactor_miss_native, cofactor_miss_squfof, cofactor_miss_rho;
static int cofactor_miss_siqs, cofactor_trace;
static unsigned int cofactor_squfof_calls, cofactor_prime_calls, cofactor_rho_calls;
static int check_cofactor_squfof(const mpz_t n, mpz_t f, UV rounds) {
  if (cofactor_trace) cofactor_squfof_calls++;
  return cofactor_miss_squfof ? 0 : squfof126(n, f, rounds);
}
static int check_cofactor_rho(const mpz_t n, mpz_t f, UV a, UV rounds) {
  if (cofactor_trace) cofactor_rho_calls++;
  return cofactor_miss_rho ? 0 : siqs_pbrent_factor(n, f, a, rounds);
}
/* A false "prime" answer makes the inner SIQS return its unsplit partition.
 * The resolver's real GMP pretest/acceptance checks are not interposed. */
static int check_cofactor_prime(const mpz_t n) {
  if (cofactor_trace) cofactor_prime_calls++;
  return cofactor_miss_siqs ? 1 : siqs_is_prob_prime(n);
}
#if BITS_PER_WORD == 64 && HAVE_STD_U64 && defined(__GNUC__) && defined(__x86_64__)
static int check_cofactor_native(UV n, UV *factors, UV rounds, UV a) {
  return cofactor_miss_native ? 0 : uvpbrent63(n, factors, rounds, a);
}
# define uvpbrent63 check_cofactor_native
#endif
#define squfof126 check_cofactor_squfof
#define siqs_pbrent_factor check_cofactor_rho
#define siqs_is_prob_prime check_cofactor_prime
