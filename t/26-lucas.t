#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/fibonacci lucas_number/;

my @fibs = (0,1,1,2,3,5,8,13,21,34,55,89,144,233,377,610,987,1597,2584,4181,6765);
my @lucas = (2,1,3,4,7,11,18,29,47,76,123,199,322,521,843,1364,2207,3571,5778,9349,15127);

# F(-n) = (-1)^(n+1) * F(n): negate F(n) when n is even
my @neg_fibs  = map { ($_ & 1) ? $fibs[$_] : -$fibs[$_] } 1 .. $#fibs;

# L(-n) = (-1)^n * L(n): negate L(n) when n is odd
my @neg_lucas = map { ($_ & 1) ? -$lucas[$_] : $lucas[$_] } 1 .. $#lucas;

plan tests => 0
            + 1 + scalar(@neg_fibs)  + 2   # fibonacci
            + 1 + scalar(@neg_lucas) + 2   # lucas_number
            + 0;

###### fibonacci
is_deeply( [map { fibonacci($_) } 0 .. $#fibs], \@fibs, "fibonacci(0..$#fibs)" );
for my $k (1 .. $#fibs) {
  is( fibonacci(-$k), $neg_fibs[$k-1], "fibonacci(-$k)" );
}
is( fibonacci(100),  "354224848179261915075", "fibonacci(100)" );
is( fibonacci(-100), "-354224848179261915075", "fibonacci(-100)" );

###### lucas_number
is_deeply( [map { lucas_number($_) } 0 .. $#lucas], \@lucas, "lucas_number(0..$#lucas)" );
for my $k (1 .. $#lucas) {
  is( lucas_number(-$k), $neg_lucas[$k-1], "lucas_number(-$k)" );
}
is( lucas_number(100),  "792070839848372253127", "lucas_number(100)" );
is( lucas_number(-100), "792070839848372253127", "lucas_number(-100)" );
