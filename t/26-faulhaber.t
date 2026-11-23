#!/usr/bin/env perl
use strict;
use warnings;

# Note: we are now using the name powersum for this function.

use Test::More;
use Math::Prime::Util::GMP qw/faulhaber_sum primepowersum prime_count primes powint addint/;

plan tests => 0
            + 4
            + 5
            + 7
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

is_deeply([map { primepowersum(0,$_), primepowersum(1,$_) } 0..6],
          [(0) x 14], "primepowersum is zero below 2");
is(primepowersum(100,0), prime_count(100), "primepowersum with power zero counts primes");
is(primepowersum(100000000,1), "279209790387276", "prime sum through 10^8");
is(primepowersum(100000000,2), "18433608754948081174274", "sum of prime squares through 10^8");
is(primepowersum(2,160), "1461501637330902918203684832716283019655932542976",
   "primepowersum supports results beyond 128 bits");

my @bounds = (2, 3, 4, 8, 9, 11, 12, 24, 25, 29, 30, 48, 49, 55, 56, 97, 100);
my (@got, @expected);
for my $k (1..8) {
  for my $n (@bounds) {
    my $sum = 0;
    $sum = addint($sum, powint($_,$k)) for @{primes($n)};
    push @got, primepowersum($n,$k);
    push @expected, $sum;
  }
}
is_deeply(\@got, \@expected, "primepowersum matches direct sums around quotient table boundaries");
ok(!eval { primepowersum("18446744073709551616",1); 1 } && $@ =~ /must fit in 64 bits/,
   "primepowersum rejects bounds wider than 64 bits for positive powers");
