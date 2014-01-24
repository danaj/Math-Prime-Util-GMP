#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util::GMP;

use Test::More  tests => 1;

my @functions = qw(
  is_prime is_prob_prime is_provable_prime is_provable_prime_with_cert
  is_aks_prime is_nminus1_prime is_ecpp_prime
  is_strong_pseudoprime is_lucas_pseudoprime is_strong_lucas_pseudoprime
  is_extra_strong_lucas_pseudoprime is_almost_extra_strong_lucas_pseudoprime
  is_frobenius_underwood_pseudoprime miller_rabin_random lucas_sequence
  primes next_prime prev_prime
  trial_factor prho_factor pbrent_factor pminus1_factor pplus1_factor
  holf_factor squfof_factor ecm_factor factor
  prime_count
  primorial pn_primorial
  consecutive_integer_lcm partitions gcd lcm kronecker
);
can_ok( 'Math::Prime::Util::GMP', @functions);
