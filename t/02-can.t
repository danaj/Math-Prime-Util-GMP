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
                     prime_count
                     primorial
                     pn_primorial
                     factorial
                     consecutive_integer_lcm
                     partitions bernfrac harmfrac harmreal stirling
                     gcd lcm kronecker valuation invmod binomial gcdext
                     invmod sqrtmod addmod mulmod divmod powmod
                     vecsum vecprod
                     exp_mangoldt
                     liouville
                     totient
                     jordan_totient
                     carmichael_lambda
                     is_power
                     is_primitive_root
                     znorder
                     znprimroot
                     ramanujan_tau
                     Pi
);
can_ok( 'Math::Prime::Util::GMP', @functions);
