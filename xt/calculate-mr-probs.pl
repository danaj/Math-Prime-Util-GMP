#!/usr/bin/env perl
use warnings;
use strict;
use List::Util qw/min max/;

# From Damgård, Landrock, and Pomerance, 1993.  Page 178.

# k = bits in n
my $k = shift;
if (!defined $k || $k < 2) {
  die "Usage: $0 <k> [<maxtests>]\n k = the number of bits in n\n";
}
my $maxtests = shift || 30;

foreach my $t (1 .. $maxtests) {
  printf("k: %d t:%4d  p < %lg\n", $k, $t, min(mr_prob($k, $t)));
  #my @p = mr_prob($k, $t); print "k: $k t: $t p: @p\n";
}


sub mr_prob {
  my($k, $t) = @_;
  die "k must be >= 2" unless $k >= 2;
  die "t must be >= 1" unless $t >= 1;
  my @probs;

  push @probs, 4**-$t;

  if ($t == 1) {
    push @probs, $k*$k * 4**(2.0-sqrt($k));
  }

  if (($t == 2 && $k >= 88) || ($t >= 3 && $k >= 21 && $t <= $k/9)) {
    push @probs, $k**1.5 * 2**$t * $t**-.5 * 4**(2.0-sqrt($t*$k));
  }

  if ($t >= $k/9 && $k >= 21) {
    push @probs, (   (7.0/20.0) * $k * 2**(-5*$t)
                   + (1.0/7.0) * $k**(15.0/4.0) * 2**(-$k/2 - 2*$t)
                   + 12 * $k * 2**(-$k/4 - 3*$t) );
  }

  if ($t >= $k/4 && $k >= 21) {
    push @probs, (1.0/7.0) * $k**(15.0/4.0) * 2**(-$k/2 - 2*$t);
  }

  return @probs;
}
