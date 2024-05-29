#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/prime_power_count/;
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};


my @A025528 = (0, 1, 2, 3, 4, 4, 5, 6, 7, 7, 8, 8, 9, 9, 9, 10, 11, 11, 12, 12, 12, 12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 18, 18, 18, 18, 18, 19, 19, 19, 19, 20, 20, 21, 21, 21, 21, 22, 22, 23, 23, 23, 23, 24, 24, 24, 24, 24, 24, 25, 25, 26, 26, 26, 27, 27, 27, 28, 28, 28, 28, 29, 29, 30, 30);
my @A267712 = (7, 35, 193, 1280, 9700, 78734, 665134, 5762859, 50851223, 455062595, 4118082969, 37607992088, 346065767406, 3204942420923, 29844572385358, 279238346816392, 2623557174778438, 24739954338671299, 234057667428388198, 2220819603016308079);

$#A025528 = 40;
$#A267712 = ($extra) ? 7 : 5;

plan tests => 0
            + 4  # prime_power_count  simple
            + 1  # large values
            + 0;

is(prime_power_count(0), 0, "prime_power_count(0) = 0");
is(prime_power_count(1), 0, "prime_power_count(1) = 0");
is_deeply( [map { prime_power_count(1+$_) } 0..$#A025528], \@A025528,  "prime_power_count(n) for 1..".scalar(@A025528) );
is_deeply( [map { prime_power_count(10**(1+$_)) } 0..$#A267712], \@A267712,  "prime_power_count(10^n) for 1..".scalar(@A267712) );

# mpu 'say vecsum(map{!!is_prime_power($_)}1..12345678)'
#is(prime_power_count(12345678), 809830, "prime_power_count(12345678) = 809830");
is(prime_power_count(1234567), 95618, "prime_power_count(1234567) = 95618");
