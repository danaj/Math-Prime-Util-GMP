#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/is_powerful
                              factor/;

plan tests => 3+1+10+4+2+1;

{
  my @exp = map { fac_is_powerful($_, 2) } 0 .. 258;
  is_deeply( [map { is_powerful($_,2) } 0..258], \@exp, "is_powerful(0..258,2)");
  is_deeply( [map { is_powerful($_) } 0..258], \@exp, "is_powerful(0..258)");
  is_deeply( [map { is_powerful($_,0) } 0..258], \@exp, "is_powerful(0..258,0)");
}

is( scalar(grep { is_powerful($_,1) } 0..32), 33, "is_powerful(n,1) = 1" );

for my $k (3 .. 12) {
  my @nums = (227411960,105218838,79368063,58308379,210322300,44982156,67831696,165946352,243118692,128757041,150085583);
  my @exp = map { fac_is_powerful($_, $k) } 0 .. 32, @nums;
  my @got = map {     is_powerful($_, $k) } 0 .. 32, @nums;
  is_deeply(\@got, \@exp, "is_powerful(n,$k) for 0..32 and 11 larger nums");
}

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


sub fac_is_powerful {
  my($n, $k) = @_;
  $k = 2 if !defined $k || $k == 0;
  return 1 if $n <= 1 || $k <= 1;
  return 0 if $n < (1<<$k);
  return 0 if (!($n%2)) && ($n%4);
  my %exponents;
  my @factors = grep { !$exponents{$_}++ } factor($n);
  for my $f (keys %exponents) {
    return 0 if $exponents{$f} < $k;
  }
  1;
}
