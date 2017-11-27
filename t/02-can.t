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
                     is_ecpp_prime
                     is_pseudoprime
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
                     miller_rabin_random
                     lucas_sequence  lucasu  lucasv
                     primes
                     sieve_primes
                     sieve_twin_primes
                     sieve_prime_cluster
                     sieve_range
                     next_prime
                     prev_prime
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
                     sigma
                     chinese
                     moebius
                     prime_count prime_count_lower prime_count_upper
                     primorial
                     pn_primorial
                     factorial factorialmod
                     consecutive_integer_lcm
                     partitions bernfrac bernreal harmfrac harmreal stirling
                     zeta li ei riemannr lambertw
                     logreal expreal powreal agmreal
                     gcd lcm kronecker valuation binomial gcdext hammingweight
                     invmod sqrtmod addmod mulmod divmod powmod
                     vecsum vecprod
                     exp_mangoldt
                     liouville
                     totient
                     jordan_totient
                     carmichael_lambda
                     sqrtint rootint logint
                     is_power is_prime_power is_semiprime is_square
                     is_carmichael is_fundamental is_totient
                     is_primitive_root
                     is_polygonal polygonal_nth
                     znorder
                     znprimroot
                     ramanujan_tau
                     Pi Euler
                     todigits
                     random_prime random_nbit_prime random_ndigit_prime
                     random_strong_prime
                     random_maurer_prime random_shawe_taylor_prime
                     random_maurer_prime_with_cert
                     random_shawe_taylor_prime_with_cert
                     seed_csprng is_csprng_well_seeded
                     irand irand64 drand urandomb urandomm urandomr random_bytes
                     permtonum numtoperm
);
can_ok( 'Math::Prime::Util::GMP', @functions);
