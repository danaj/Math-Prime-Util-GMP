#!/usr/bin/perl
use warnings;
use strict;
use v5.32;
use Math::Prime::Util ":all";
use Math::Prime::Util::GMP;
use Math::Prime::Util::PP;

#say Math::Prime::Util::GMP::is_pseudoprime(3,3);  exit(1);
#say _is_pseudoprime(3,3);  exit(1);

#say Math::Prime::Util::is_strong_pseudoprime(367,1101);  exit(1);
#say Math::Prime::Util::GMP::is_strong_pseudoprime(367,1101);  exit(1);
#say _is_strong_pseudoprime(367,1101);  exit(1);

{
  my(@IM, @IG, @IP);
  for my $n (2..500) {
    for my $base (2 .. 1500) {
      my $p1 = _is_pseudoprime($n,$base);
      my $p2 = is_pseudoprime($n,$base);
      my $p3 = Math::Prime::Util::GMP::is_pseudoprime($n,$base);
      my $p4 = Math::Prime::Util::PP::is_pseudoprime($n,$base);
      push @IM, "$n $base     $p1  $p2" unless $p1 == $p2;
      push @IG, "$n $base     $p1  $p3" unless $p1 == $p3;
      push @IP, "$n $base     $p1  $p4" unless $p1 == $p4;
    }
  }
  showerr("is_pseudoprime", "MPU", @IM);
  showerr("is_pseudoprime", "GMP", @IG);
  showerr("is_pseudoprime", "PP",  @IP);
}
print "\n";
{
  my(@IM, @IG, @IP);
  for my $n (2..500) {
    for my $base (2 .. 1500) {
      my $p1 = _is_euler_pseudoprime($n,$base);
      my $p2 = is_euler_pseudoprime($n,$base);
      my $p3 = Math::Prime::Util::GMP::is_euler_pseudoprime($n,$base);
      my $p4 = Math::Prime::Util::PP::is_euler_pseudoprime($n,$base);
      push @IM, "$n $base     $p1  $p2" unless $p1 == $p2;
      push @IG, "$n $base     $p1  $p3" unless $p1 == $p3;
      push @IP, "$n $base     $p1  $p4" unless $p1 == $p4;
    }
  }
  showerr("is_euler_pseudoprime", "MPU", @IM);
  showerr("is_euler_pseudoprime", "GMP", @IG);
  showerr("is_euler_pseudoprime", "PP",  @IP);
}
print "\n";
{
  my(@IM, @IG, @IP);
  for my $n (2..500) {
    for my $base (2 .. 1500) {
      my $p1 = _is_strong_pseudoprime($n,$base);
      my $p2 = is_strong_pseudoprime($n,$base);
      my $p3 = Math::Prime::Util::GMP::is_strong_pseudoprime($n,$base);
      my $p4 = Math::Prime::Util::PP::is_strong_pseudoprime($n,$base);
      push @IM, "$n $base     $p1  $p2" unless $p1 == $p2;
      push @IG, "$n $base     $p1  $p3" unless $p1 == $p3;
      push @IP, "$n $base     $p1  $p4" unless $p1 == $p4;
    }
  }
  showerr("is_strong_pseudoprime", "MPU", @IM);
  showerr("is_strong_pseudoprime", "GMP", @IG);
  showerr("is_strong_pseudoprime", "PP",  @IP);
}


sub showerr {
  my($test, $code, @arr) = @_;
  my $nerr = scalar(@arr);
  if ($nerr == 0) {
    print "$test $code  OK\n";
  } else {
    print "$test $code has $nerr mismatches.\n     $arr[0]\n     $arr[1]\n"
  }
}
sub _is_pseudoprime {
  my($n,$a)=@_;
  return 1 if $n == 2;
  #return 0 if ($n & 1) == 0;
  #$a = modint($a,$n);
  #return $a if $a <= 1;
  (powmod($a,$n-1,$n) == 1) ? 1 : 0;
}
sub _is_euler_pseudoprime {
  my($n,$a)=@_;
  return 1 if $n == 2;
  return 0 if ($n & 1) == 0;
  #$a = modint($a,$n);
  #return $a if $a <= 1;
  return 0 if gcd($a,$n) != 1;
  (modint(kronecker($a,$n),$n) == powmod($a, ($n-1)>>1, $n)) ? 1 : 0;
}
sub _is_strong_pseudoprime {
  my($n,$a)=@_;
  return 1 if $n == 2; 
  return 0 if ($n & 1) == 0;
  my $nm1 = $n-1;

  # Must make a decision here about what to do if a = 0 mod n.
  # If we choose to return 0, then actual primes will return 0.
  # If we choose to return 1, then we aren't following the strict and simple
  #    definition, though we can argue this is outside the scope.
  #
  # Crandall and Pomerance, page 135, theorem "If a is not divisible by n..."
  # and later define the test only for 1 < a < n-1  (i.e. a between 2 and n-2)..
  #
  # Therefore I suggest it makes sense to return 1 when a is divisible by n.

  # The best way to do it:
  #$a = modint($a,$n);
  #return 1 if $a == $nm1;
  #return 1 if $a <= 1;   # compare with:   return $a if $a <= 1;

  # The most basic change:
  return 1 if $a % $n == 0;    # Do as little as possible

  my $d = $nm1;
  my $s = 0;
  do {
    $s++;
    $d >>= 1; 
  } until ($d & 1);
  my $t = powmod($a,$d,$n);
  return 1 if $t == 1 || $t == $nm1;
  for my $r (1 .. $s-1) {
    return 1 if powmod($a,$d << $r,$n) == $nm1;
  }
  0;
}
