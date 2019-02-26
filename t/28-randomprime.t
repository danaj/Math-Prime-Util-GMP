#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::BigInt try=>"GMP,Pari";
use Math::Prime::Util::GMP
    qw/random_prime random_ndigit_prime random_nbit_prime
       random_safe_prime random_strong_prime
       random_maurer_prime random_shawe_taylor_prime
       random_maurer_prime_with_cert
       random_shawe_taylor_prime_with_cert
       is_prime
       seed_csprng/;

my $use64 = (~0 > 4294967295);
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $maxbits = $use64 ? 64 : 32;

my @random_to = (2, 3, 4, 5, 6, 7, 8, 9, 100, 1000, 1000000, 4294967295);

my @random_nbit_tests = (2 .. 6, 10, 30 .. 34, 62 .. 66, 126 .. 130);
push @random_nbit_tests, (256, 512, 1024, 2048) if $extra;

my @random_ndigit_tests = (1 .. 25);

my @random_safe_tests   = (16, 32, 33, 34, 64, 128, 255, 256, 512);
my @random_strong_tests = (128, 255, 256, 512);

my %ranges = (
  "2 to 20" => [2,19],
  "3 to 7" => [3,7],
  "20 to 100" => [23,97],
  "5678 to 9876" => [5683,9871],
  "27767 to 88493" => [27767,88493],
  "27764 to 88498" => [27767,88493],
  "27764 to 88493" => [27767,88493],
  "27767 to 88498" => [27767,88493],
  "17051687 to 17051899" => [17051687,17051899],
  "17051688 to 17051898" => [17051707,17051887],
);

my %range_edge = (
  "0 to 2" => [2,2],
  "2 to 2" => [2,2],
  "2 to 3" => [2,3],
  "3 to 5" => [3,5],
  "10 to 20" => [11,19],
  "8 to 12" => [11,11],
  "10 to 12" => [11,11],
  "16706143 to 16706143" => [16706143,16706143],
  "16706142 to 16706144" => [16706143,16706143],
  "3842610773 to 3842611109" => [3842610773,3842611109],
  "3842610772 to 3842611110" => [3842610773,3842611109],
);
my %range_edge_empty = (
  "0 to 0" => [],
  "0 to 1" => [],
  "2 to 1" => [],
  "3 to 2" => [],
  "1294268492 to 1294268778" => [],
  "3842610774 to 3842611108" => [],
);

plan tests => 0
              + (1 * scalar (keys %range_edge_empty))
              + (3 * scalar (keys %range_edge))
              + (2 * scalar (keys %ranges))
              + (2 * scalar @random_to)
              + (1 * scalar @random_ndigit_tests)
              + (1 * scalar @random_nbit_tests)
              + (2 * scalar @random_safe_tests)
              + (1 * scalar @random_strong_tests)
              + (2 * scalar @random_nbit_tests)   # proven primes
              + 6
              + 0;

my $infinity = 20**20**20;
my $nrandom_range_samples = $extra ? 1000 : 50;

while (my($range, $expect) = each (%range_edge_empty)) {
  my($low,$high) = $range =~ /(\d+) to (\d+)/;
  is( random_prime($low,$high), undef, "primes($low,$high) should return undef" );
}

while (my($range, $expect) = each (%range_edge)) {
  my($low,$high) = $range =~ /(\d+) to (\d+)/;
  my $got = random_prime($low,$high);
  ok( is_prime($got), "Prime in range $low-$high is indeed prime" );
  cmp_ok( $got, '>=', $expect->[0], "random_prime($low,$high) >= $expect->[0]");
  cmp_ok( $got, '<=', $expect->[1], "random_prime($low,$high) <= $expect->[1]");
}

while (my($range, $expect) = each (%ranges)) {
  my($low,$high) = $range =~ /(\d+) to (\d+)/;
  my $isprime = 1;
  my $inrange = 1;
  for (1 .. $nrandom_range_samples) {
    my $got = random_prime($low,$high);
    $isprime *= is_prime($got) ? 1 : 0;
    $inrange *= (($got >= $expect->[0]) && ($got <= $expect->[1])) ? 1 : 0;
  }
  ok($isprime, "All returned values for $low-$high were prime" );
  ok($inrange, "All returned values for $low-$high were in the range" );
}

foreach my $high (@random_to) {
  my $isprime = 1;
  my $inrange = 1;
  for (1 .. $nrandom_range_samples) {
    my $got = random_prime(0,$high);
    $isprime *= is_prime($got) ? 1 : 0;
    $inrange *= (($got >= 2) && ($got <= $high)) ? 1 : 0;
  }
  ok($isprime, "All returned values for $high were prime" );
  ok($inrange, "All returned values for $high were in the range" );
}

foreach my $digits ( @random_ndigit_tests ) {
  my $n = random_ndigit_prime($digits);
  ok ( length($n) == $digits && is_prime($n),
       "$digits-digit random prime '$n' is in range and prime");
}

foreach my $bits ( @random_nbit_tests ) {
  check_bits( random_nbit_prime($bits), $bits, "random $bits-bit" );
}

foreach my $bits ( @random_safe_tests ) {
  my $p = random_safe_prime($bits);
  my $q = (Math::BigInt->new("$p") - 1) >> 1;
  check_bits( $p, $bits, "random $bits-bit safe (p)" );
  check_bits( $q, $bits-1, "random $bits-bit safe (q)" );
}

foreach my $bits ( @random_strong_tests ) {
  # TODO: Check that p-1 and p+1 have a big factor.  Might be really slow.
  check_bits( random_strong_prime($bits), $bits, "random $bits-bit strong" );
}

foreach my $bits ( @random_nbit_tests ) {
  check_bits( random_maurer_prime($bits), $bits, "random $bits-bit proven (Maurer)" );
  check_bits( random_shawe_taylor_prime($bits), $bits, "random $bits-bit proven (Shawe-Taylor)" );
}

sub check_bits {
  my($n, $bits, $what) = @_;
  $n = Math::BigInt->new("$n");
  my $min = Math::BigInt->new(1)->blsft($bits-1);
  my $max = Math::BigInt->new(1)->blsft($bits)->bdec;
  ok ( $n >= $min && $n <= $max && is_prime($n),
       "$bits-bit random $what prime '$n' is in range and prime");
}

# Now check with seed
seed_csprng(3,"xyz");
is( random_nbit_prime(24), 10207999, "random 20-bit prime with seeded rng" );
is( random_ndigit_prime(9), 842208331, "random 9-digit with seeded rng" );

{
  my($n,$cert) = random_maurer_prime_with_cert(139);
  is( $n, "461450530708271369919467043900581427327289", "random Maurer prime" );
  is( $cert, "[MPU - Primality Certificate]
Version 1.0

Proof for:
N 461450530708271369919467043900581427327289

Type BLS3
N 461450530708271369919467043900581427327289
Q 399271657906884417934953679
A 23

Type BLS3
N 399271657906884417934953679
Q 125756461377167
A 3

Type BLS3
N 125756461377167
Q 18944413
A 5
", "random Maurer prime certificate" );
}

{
  my($n,$cert) = random_shawe_taylor_prime_with_cert(147);
  is( $n, "175813161209418279279866476297230082927452359", "random Shawe-Taylor prime" );
  is( $cert, "[MPU - Primality Certificate]
Version 1.0

Proof for:
N 175813161209418279279866476297230082927452359

Type Pocklington
N 175813161209418279279866476297230082927452359
Q 34388730545435187832661
A 67309515455075813303416848806399820637730471

Type Pocklington
N 34388730545435187832661
Q 482980495961
A 23452682901285108721511

Type Pocklington
N 482980495961
Q 1061869
A 178206865367
", "random Shawe-Taylor prime certificate" );
}
