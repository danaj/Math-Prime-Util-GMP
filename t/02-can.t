#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util::GMP;

use Test::More  tests => 1;

my @functions = qw(
                     is_prime
                     is_prob_prime
                     is_bpsw_prime
                     is_provable_prime
                     is_provable_prime_with_cert
                     is_aks_prime
                     is_nminus1_prime
                     is_nplus1_prime
                     is_bls75_prime
                     is_ecpp_prime
                     is_pseudoprime
                     is_euler_pseudoprime
                     is_euler_plumb_pseudoprime
                     is_strong_pseudoprime
                     is_lucas_pseudoprime
                     is_strong_lucas_pseudoprime
                     is_extra_strong_lucas_pseudoprime
                     is_almost_extra_strong_lucas_pseudoprime
                     is_perrin_pseudoprime
                     is_frobenius_pseudoprime
                     is_frobenius_underwood_pseudoprime
                     is_frobenius_khashin_pseudoprime
                     is_mersenne_prime
                     is_llr_prime
                     is_proth_prime
                     is_miller_prime
                     miller_rabin_random
                     is_gaussian_prime
                     lucas_sequence  lucasuv  lucasu  lucasv
                     lucasuvmod  lucasumod  lucasvmod
                     primes
                     sieve_primes
                     sieve_twin_primes
                     sieve_prime_cluster
                     sieve_range
                     next_prime
                     prev_prime
                     surround_primes
                     trial_factor
                     prho_factor
                     pbrent_factor
                     pminus1_factor
                     pplus1_factor
                     holf_factor
                     squfof_factor
                     ecm_factor
                     qs_factor
                     factor
                     divisors
                     sigma
                     chinese
                     moebius
                     prime_count prime_count_lower prime_count_upper
                     primorial
                     pn_primorial
                     factorial subfactorial multifactorial factorial_sum
                     factorialmod
                     consecutive_integer_lcm
                     partitions bernfrac bernreal harmfrac harmreal stirling
                     zeta li ei riemannr lambertw
                     addreal subreal mulreal divreal
                     logreal expreal powreal rootreal agmreal
                     gcd lcm kronecker valuation binomial gcdext hammingweight
                     is_qr
                     invmod sqrtmod addmod submod mulmod divmod powmod
                     muladdmod mulsubmod
                     vecsum vecprod
                     exp_mangoldt
                     liouville
                     totient
                     jordan_totient
                     carmichael_lambda
                     prime_omega prime_bigomega
                     sqrtint rootint logint powint mulint addint subint
                     divint modint cdivint divrem tdivrem fdivrem cdivrem
                     add1int sub1int
                     negint absint signint cmpint cmpabsint
                     lshiftint rshiftint rashiftint
                     setbit clrbit tstbit
                     bitand bitor bitxor
                     is_power is_prime_power is_semiprime is_almost_prime
                     is_square is_smooth is_rough is_powerful is_practical
                     is_carmichael is_fundamental is_totient
                     is_primitive_root
                     is_polygonal polygonal_nth
                     powerful_count
                     perfect_power_count prime_power_count
                     znorder
                     znprimroot
                     ramanujan_tau
                     Pi Euler
                     todigits fromdigits
                     random_prime random_nbit_prime random_ndigit_prime
                     random_safe_prime random_strong_prime
                     random_maurer_prime random_shawe_taylor_prime
                     random_maurer_prime_with_cert
                     random_shawe_taylor_prime_with_cert
                     seed_csprng is_csprng_well_seeded
                     irand irand64 drand urandomb urandomm urandomr random_bytes
                     permtonum numtoperm
);
can_ok( 'Math::Prime::Util::GMP', @functions);
