#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/seed_csprng irand irand64 drand/;

my $use64 = (~0 > 4294967295);
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $maxbits = $use64 ? 64 : 32;

my $samples = $extra ? 100000 :  10000;

plan tests => 0
            + 2
            + ($use64 ? 2 : 0)
            + 2;

########

# seed_csprng(55,"BLAKEGrostlJHKeccakSkein--RijndaelSerpentTwofishRC6MARS");

########

{
  my @s = map { irand } 1 .. $samples;
  is( scalar(grep { $_ > 4294967295 } @s), 0, "irand values are 32-bit" );
  is( scalar(grep { $_ != int($_) } @s), 0, "irand values are integers" );
}

########

if ($use64) {
  my $bits_on  = 0;
  my $bits_off = 0;
  my $iter = 0;
  for (1 .. 6400) {
    $iter++;
    my $v = irand64;
    $bits_on |= $v;
    $bits_off |= (~$v);
    last if ~$bits_on == 0 && ~$bits_off == 0;
  }
  is( ~$bits_on,  0, "irand64 all bits on in $iter iterations" );
  is( ~$bits_off, 0, "irand64 all bits off in $iter iterations" );
}

########

{
  my $mask = 0;
  my $v;
  for (1..1024) {
    $v = drand;
    last if $v >= 1;
    next if $v < .5;
    for my $b (0..127) {
      last unless $v;
      $v *= 2;
      if ($v >= 1) {
        $mask |= (1 << $b);
        $v -= 1;
      }
    }
  }
  ok($v < 1, "drand values between 0 and 1-eps");
  my $k = 0; while ($mask) { $k++; $mask >>= 1; }
  # Assuming drand is working properly:
  #   k = 24   NV is float
  #   k = 53   NV is double
  #   k = 64   NV is long double
  # If we used drand48 we'd get 48 with double or long double.
  ok($k >= 21, "drand supplies at least 21 bits (got $k)");
}

