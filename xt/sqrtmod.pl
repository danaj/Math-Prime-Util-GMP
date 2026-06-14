#!/usr/bin/perl
use warnings;
use strict;
use v5.32;
use Math::Prime::Util ":all";

my @fails;

for my $n (0 .. 1000) {
  my $ok = 1;
  for my $a (0 .. $n-1) { $ok = compare($a,$n); }
  for my $a ($n .. 2*$n) { $ok = compare($a,$n); }
  for my $a (-2*$n .. -1) { $ok = compare($a,$n); }
}

say "Total ",scalar(@fails)," failures";
$#fails = 4 if @fails > 5;
print join "\n", @fails;

sub compare {
  my($a, $n) = @_;
  my $s1 = sqrtmod($a,$n);
  my $s2 = Math::Prime::Util::GMP::sqrtmod($a,$n);
  return 1 if !defined $s1 && !defined $s2;
  if (defined $s1 && defined $s2) {
    return 1 if $s1 == $s2;
    return 1 if mulmod($s1,$s1,$n) == modint($a,$n) && mulmod($s2,$s2,$n) == modint($a,$n);
  }
  $s1 = '<undef>' unless defined $s1;
  $s2 = '<undef>' unless defined $s2;
  push @fails, "sqrtmod($a,$n) = $s1  (GMP gave $s2)";
  return 0;
}
