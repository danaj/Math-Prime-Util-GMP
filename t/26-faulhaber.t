#!/usr/bin/env perl
use strict;
use warnings;

# Note: we are now using the name powersum for this function.

use Test::More;
use Math::Prime::Util::GMP qw/faulhaber_sum/;

plan tests => 0
            + 4
            + 5
            + 0;

is_deeply( [map { faulhaber_sum(0,$_) } 0..6], [map { 0 } 0..6], "faulhaber_sum(0,n) = 0" );
is_deeply( [map { faulhaber_sum(1,$_) } 0..6], [map { 1 } 0..6], "faulhaber_sum(1,n) = 1" );
is_deeply( [map { faulhaber_sum($_,0) } 0..6], [map { $_ } 0..6], "faulhaber_sum(n,0) = n" );
is_deeply( [map { faulhaber_sum($_,1) } 0..6], [map { ($_*($_+1))>>1 } 0..6], "faulhaber_sum(n,1) = n*(n+1)/2" );

is(faulhaber_sum(27,2), 6930, "faulhaber(27,2) = 1^2 + 2^2 + 3^2 + ... + 27^2");
is(faulhaber_sum(24,6), 754740700, "faulhaber(24,6) = 1^6 + 2^6 + 3^6 + ... + 24^6");
is(faulhaber_sum(15,13), "3197503726489920", "faulhaber(15,13) = 1^13 + 2^13 + ... + 15^13");
is(faulhaber_sum(1000000,3), "250000500000250000000000", "faulhaber(1000000,3) = 1^3 + 2^3 + ... + 1000000^3");
is(faulhaber_sum(3,75), "608266787713395488051546949780570876", "faulhaber(3,75)");
