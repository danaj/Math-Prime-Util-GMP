#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/is_perfect_power
                              next_perfect_power prev_perfect_power
                              perfect_power_count
                              nth_perfect_power nth_perfect_power_approx/;

my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};


my @A069623 = (1, 1, 1, 2, 2, 2, 2, 3, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12);
my @A070428 = (1, 4, 13, 41, 125, 367, 1111, 3395, 10491, 32670, 102231, 320990, 1010196, 3184138, 10046921, 31723592, 100216745, 316694005, 1001003332, 3164437425, 10004650118, 31632790244, 100021566157, 316274216762, 1000100055684);
my @A001597 = (1, 4, 8, 9, 16, 25, 27, 32, 36, 49, 64, 81, 100, 121, 125, 128, 144, 169, 196, 216, 225, 243, 256, 289, 324, 343, 361, 400, 441, 484, 512, 529, 576, 625, 676, 729, 784, 841, 900, 961, 1000, 1024, 1089, 1156, 1225, 1296, 1331, 1369, 1444, 1521, 1600, 1681, 1728, 1764);

$#A069623 = 40;
$#A070428 = 10;

plan tests => 0
            + 4    # is_perfect_power
            + 9    # next / prev
            + 4    # count  basic tests
            + 1    # count  large value
            + 2    # count  ranges
            + 3    # nth
            + 3    # approx
            + 0;

######  is_perfect_power

is_deeply( [map { is_perfect_power($_) } 0..10], [0,1,0,0,1,0,0,0,1,1,0], "is_perfect_power(0 .. 10)" );
is_deeply( [grep { is_perfect_power($_) } -100..100], [qw/-64 -32 -27 -8 -1 1 4 8 9 16 25 27 32 36 49 64 81 100/] , "is_perfect_power(-100 .. 100)" );
is( is_perfect_power("18446744065119617025"), 1, "is_perfect_power(18446744065119617025)" );
is( is_perfect_power("18446744073709551616"), 1, "is_perfect_power(18446744073709551616)" );


######  next / prev

is_deeply( [map { next_perfect_power($_) } 0..20],
           [1,4,4,4,8,8,8,8,9,16,16,16,16,16,16,16,25,25,25,25,25],
           "next perfect power with small inputs" );
is_deeply( [map { prev_perfect_power($_) } 0..20],
           [-1,-1,1,1,1,4,4,4,4,8,9,9,9,9,9,9,9,16,16,16,16],
           "prev perfect power with small inputs" );
is( next_perfect_power("18446744065119617025"), "18446744073709551616", "next_perfect_power(18446744065119617025)" );
is( prev_perfect_power("18446744073709551616"), "18446744065119617025", "prev_perfect_power(18446744073709551616)" );

is_deeply( [map { next_perfect_power($_) } -9 .. 9],
           [-8,-1,-1,-1,-1,-1,-1,-1,1,1,4,4,4,8,8,8,8,9,16],
           "next perfect power with small inputs around zero" );
is_deeply( [map { prev_perfect_power($_) } -9 .. 9],
           [-27,-27,-8,-8,-8,-8,-8,-8,-8,-1,-1,1,1,1,4,4,4,4,8],
           "prev perfect power with small inputs around zero" );
is( prev_perfect_power(-8), -27, "prev_perfect_power(-8) = -27" );
is( next_perfect_power(-64), -32, "prev_perfect_power(-64) = -32" );
is( prev_perfect_power(-64), -125, "prev_perfect_power(-64) = -125" );


######  perfect_power_count
is(perfect_power_count(0), 0, "perfect_power_count(0) = 0");
is(perfect_power_count(1), 1, "perfect_power_count(1) = 1");
is_deeply( [map { perfect_power_count(1+$_) } 0..$#A069623], \@A069623,  "perfect_power_count(n) for 1..".scalar(@A069623) );
is_deeply( [map { perfect_power_count(10**$_) } 0..$#A070428], \@A070428,  "perfect_power_count(10^n) for 0..$#A070428" );

# mpu 'say 1+vecsum(map{!!is_power($_)}1..12345678)'
is(perfect_power_count(12345678), 3762, "perfect_power_count(12345678) = 3762");


is( perfect_power_count(123456, 133332), 17, "perfect_power_count(123456,133332) = 17" );
is_deeply( [map { perfect_power_count($_,16) } 8,9,10],
           [3,2,1],
           "perfect_power_count(8..10,16) = 3,2,1" );

######  nth_perfect_power

is_deeply( [map { nth_perfect_power($_) } 1 .. scalar(@A001597)],
           \@A001597,
           "nth perfect_powers creates A001597" );

is_deeply( [map { nth_perfect_power($_) } 67224..67229],
           [qw/4294574089 4294705156 4294836225 4294967296 4295098369 4295229444/],
           "nth perfect powers with results around 2^32" );
is_deeply( [map { nth_perfect_power($_) } 4297615579..4297615582],
           [qw/18446744047939747849 18446744056529682436 18446744065119617025 18446744073709551616/],
           "nth perfect powers with results around 2^64" );

######  approx for count and nth

is( approx_in_range(1571,2048383), 2048383, "perfect power approx for 1571" );
is( approx_in_range(59643,3373286400), 3373286400, "perfect power approx for 59643" );
is( approx_in_range(15964377,"252826822479841"), "252826822479841", "perfect power approx for 15964377" );



sub approx_in_range {
  my($n,$rn) = @_;
  my $arn = nth_perfect_power_approx($n);
  #my $an  = perfect_power_count_approx($rn);
  return 'nth approx too low' if $arn < ($rn-$rn/100);
  return 'nth approx too high' if $arn > ($rn+$rn/100);
  #return 'count approx too low' if $an < ($n-$n/100);
  #return 'count approx too high' if $an > ($n+$n/100);
  $rn;
}
