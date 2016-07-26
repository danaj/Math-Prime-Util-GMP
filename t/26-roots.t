#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/sqrtint rootint/;
use Math::BigInt;  # Don't use GMP so we don't have to work around bug

my $use64 = (~0 > 4294967296 && 18446744073709550592 != ~0);

plan tests => 0
            + 3 #sqrtint
            + 7 #rootint
            + 0;

###### sqrtint
is_deeply( [map { sqrtint($_**2) } 0..10], [0..10], "sqrtint 0-10" );
is_deeply( [map { sqrtint($_**2-1) } 2..20], [2-1..20-1], "sqrtint 2-20 -1" );
is( sqrtint("647307989872865201422284359961937038113215496061434545237"), "25442248129299918480155830589", "sqrtint(13^51)" );

###### rootint
is( rootint("43091031920942300256108314560009772304748698124094750326895058640841523270081624169128280918534127523222564290447104831706207227117677890695945149868732770531628297914633063561406978145215542597509491443634033203125", 23), 2147483645, "rootint( (2^31-3)^23, 23) = 2^31-3" );
is( rootint("43091031920942300256108314560009772304748698124094750326895058640841523270081624169128280918534127523222564290447104831706207227117677890695945149868732770531628297914633063561406978145215542597509491443634033203124", 23), 2147483644, "rootint( (2^31-3)^23-1, 23) = 2^31-3-1" );
is( rootint(Math::BigInt->new(10)->bpow(1000), 1001), 9, "rootint(10^1000,1001) = 9");
is( rootint(Math::BigInt->new(2)->bpow(240), 9), 106528681, "rootint(2^240,9) = 106528681");
is( root_test(2,100), 0, "roots of powers of 2" );
is( root_test(Math::BigInt->new("4294967297"),10), 0, "roots of powers of 2^32+1" );
is( root_test(Math::BigInt->new("18446744073709551617"),10), 0, "roots of powers of 2^64+1" );

# Similar to Pari/GP's test
sub root_test {
  my($base, $max) = @_;
  my $B = Math::BigInt->new($base);
  for my $i (2 .. $max) {
    $B *= $base;
    my $root = rootint($B+1, $i);
    return "$i $B: $root should be $base" if $root != $base;
  }
  0;
}
