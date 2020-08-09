#!/usr/bin/perl
use warnings;
use strict;
use v5.32;

use Math::BigInt try=>"GMP";
#use Math::Prime::Util qw/addmod mulmod divmod powmod invmod sqrtmod/;
use Math::Prime::Util::GMP qw/addmod mulmod divmod powmod invmod sqrtmod/;

for my $p (1 .. 20) {       my $P = Math::BigInt->new($p);
  for my $a (-20..20) {     my $A = Math::BigInt->new($a);
    for my $b (-20..20) {   my $B = Math::BigInt->new($b);
      if (1) {
        die "addmod($a,$b,$p)" unless addmod($a,$b,$p) == ($A + $B) % $P;
      }
      if (1) {
        die "mulmod($a,$b,$p)" unless mulmod($a,$b,$p) == ($A * $B) % $P;
      }
      if (1) {
        my $pm1 = powmod($a,$b,$p);
        my $pm2 = $A->copy->bmodpow($B,$P);
        die "powmod($a,$b,$p) unequal" if defined $pm1 && $pm1 != $pm2;
        die "powmod($a,$b,$p) missing" if !defined $pm1 && $pm2->is_int();
      }
      if (1) {
        # Both Pari and Math::BigInt say:
        #   divmod(a,0,p) = 0 if p=1, range error otherwise
        my $dm1 = divmod($a,$b,$p);
        my $dm2 = $B->copy->bmodinv($P)->bmul($A)->bmod($P);
        die "divmod($a,$b,$p) unequal" if defined $dm1 && $dm1 != $dm2;
        die "divmod($a,$b,$p) missing" if !defined $dm1 && $dm2->is_int();
      }
      # invmod
      # sqrtmod
    }
  }
}

die "Bad divmod 1,1,1000000001" unless divmod(1,1,1000000001) == 1;
die "Bad divmod 43,1,1000000001" unless divmod(43, 1, 1000000001) == 43;
die "Bad divmod 127,1,97" unless divmod(127, 1, 97) == 30;

# TODO test handle zero P
