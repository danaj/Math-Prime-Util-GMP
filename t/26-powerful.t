#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/is_powerful powerful_count
                              factor/;

plan tests => 2+1+10+5+4+2+1
            + 8 + 2;

{
  my @exp = map { fac_is_powerful($_, 2) } 0 .. 258;
  is_deeply( [map { is_powerful($_,2) } 0..258], \@exp, "is_powerful(0..258,2)");
  is_deeply( [map { is_powerful($_) } 0..258], \@exp, "is_powerful(0..258)");
}

is( scalar(grep { is_powerful($_,1) } 0..32), 32, "is_powerful(n,1) = 1 for positive n" );

for my $k (3 .. 12) {
  my @nums = (227411960,105218838,79368063,58308379,210322300,44982156,67831696,165946352,243118692,128757041,150085583);
  my @exp = map { fac_is_powerful($_, $k) } 0 .. 32, @nums;
  my @got = map {     is_powerful($_, $k) } 0 .. 32, @nums;
  is_deeply(\@got, \@exp, "is_powerful(n,$k) for 0..32 and 11 larger nums");
}

is( is_powerful(0), 0, "is_powerful(0) returns false");
is( is_powerful(-16,0), 0, "is_powerful(-16,0) returns false");
is( is_powerful(-32,1), 0, "is_powerful(-32,1) returns false");
is( is_powerful(-64,2), 0, "is_powerful(-64,2) returns false");
is( is_powerful(-128,3), 0, "is_powerful(-128,3) returns false");

{
  my @pow2 = map { 5*5 * $_*$_ } 1..50;
  my @npow2 = map { 149 * $_*$_ } 1..50;
  my @pow3 = map { 7*7*7 * $_*$_*$_ } 1..50;
  my @npow3 = map { 4489 * $_*$_*$_ } 1..50;

  is( scalar(grep { is_powerful($_,2) } @pow2), scalar(@pow2), "small is_powerful(n,2), n powerful" );
  is( scalar(grep { is_powerful($_,3) } @pow3), scalar(@pow3), "small is_powerful(n,3), n powerful" );
  is( scalar(grep { is_powerful($_,2) } @npow2), 0, "small is_powerful(n,2), n not powerful" );
  is( scalar(grep { is_powerful($_,3) } @npow3), 0, "small is_powerful(n,3), n not powerful" );
}

is( is_powerful("1377276413364943226363244108454842276965894752197358387200000"), 0, "large easy non-powerful number" );
is( is_powerful("2346889178458529643625998598305409091755415961600000"), 1, "large easy powerful number" );

is( is_powerful("56648008573112538662596929676588737208124071038924666321487873929306609840197", 30), 0, "256-bit semiprime is not 30-powerful, without factoring" );


###### powerful_count
is_deeply( [map { powerful_count($_,1) } 1..20], [1..20], "powerful_count(n,1)=n" );

is_deeply( [map { powerful_count($_) } 1..20],
           [1,1,1,2,2,2,2,3,4,4,4,4,4,4,4,5,5,5,5,5],
           "powerful_count(1..20)" );
is_deeply( [map { powerful_count($_,3) } 1..20],
           [1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3],
           "powerful_count(1..20,3)" );

# test with n <= 0 and k = 0,1,2
is_deeply( [map { powerful_count($_,0) } -100, -10, -4, -1, 0, 1, 4, 10, 100],
           [0, 0, 0, 0, 0, 1, 4, 10, 100], "powerful_count(+/- n, 0)" );
is_deeply( [map { powerful_count($_,0) } -100, -10, -4, -1, 0, 1, 4, 10, 100],
           [0, 0, 0, 0, 0, 1, 4, 10, 100], "powerful_count(+/- n, 1)" );
is_deeply( [map { powerful_count($_,2) } -100, -10, -4, -1, 0, 1, 4, 10, 100],
           [0, 0, 0, 0, 0, 1, 2,  4 , 14], "powerful_count(+/- n, 2)" );
{
  #my $maxk = ($extra) ? 30 : 15;
  my $maxk = 30;
  my @pow14=(0,14,100,432, 2048, 6561, 16384, 59049, 131072, 262144, 531441, 1594323, 4194304, 8388608, 16777216, 43046721, 129140163, 268435456, 536870912, 1162261467, 3486784401, 8589934592, 17179869184, 34359738368, 68719476736, 137438953472, 274877906944, 549755813888, 1099511627776, 2199023255552, 4398046511104);
  is_deeply( [map { powerful_count($pow14[$_],$_) } 1..$maxk],
             [map { 14 } 1..$maxk],
             "powerful_count(x,1..$maxk) = 14" );
  is_deeply( [map { powerful_count($pow14[$_]-1,$_) } 1..$maxk],
             [map { 13 } 1..$maxk],
             "powerful_count(x-1,1..$maxk) = 13" );
}

if (1) {
  my @exp = (4,14,54,185,619,2027,6553,21044,67231,214122,680330,2158391,6840384,21663503);
  my @got = map { powerful_count(10**$_) } 1..14;
  is_deeply(\@got, \@exp, "2-powerful_count 10^1, 10^2, ..., 10^14");
}
if (1) {
  my @exp = (1, 1, 4, 10, 16, 26, 46, 77, 129, 204, 318, 495, 761, 1172, 1799, 2740, 4128, 6200, 9224, 13671, 20205, 29764);
  my @got = map { powerful_count("1".("0"x$_),7) } 1..22;
  is_deeply(\@got, \@exp, "7-powerful_count 10^1, 10^2, ..., 10^22");
}


############################################

sub fac_is_powerful {
  my($n, $k) = @_;
  $k = 2 if !defined $k;
  return 0 if $n < 1;
  return 1 if $n == 1 || $k <= 1;
  return 0 if $n < (1<<$k);
  return 0 if (!($n%2)) && ($n%4);
  my %exponents;
  my @factors = grep { !$exponents{$_}++ } factor($n);
  for my $f (keys %exponents) {
    return 0 if $exponents{$f} < $k;
  }
  1;
}
