#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/znlog/;

my @znlogs = (
  [ [5,2,1019], 10 ],
  [ [2,4,17], undef ],
  [ [7,17,36], undef ],
  [ [1,8,9], [0,2,4,6,8] ],
  [ [3,3,8], [1,3,5,7] ],
  [ [10,2,101], 25 ],
  [ [5,2,401], [48,248] ],
  [ [5678,5,10007], 8620 ],
  [ [0,30,100], 2 ],
  [ [0,2,8], 3 ],
  [ [8,2,102], 3 ],
  [ [18,18,102], 1 ],
  [ [130,85,177], 15 ],
  [ [79,92,129], 2 ],
  [ [12,42,122], 13 ],
  [ [36,44,50], 2 ],
  [ [34,170,187], 5 ],
  [ [3,4,7], undef ],
  [ [3,2,4], undef ],
  [ [15,2,"1000000000000000000117"], "77714890122519843915" ],
);

plan tests => scalar @znlogs;

foreach my $arg (@znlogs) {
  my($aref, $exp) = @$arg;
  my($a, $g, $n) = @$aref;
  my $k = znlog($a, $g, $n);

  if (defined $exp && ref($exp)) {
    ok(is_one_of($k, @$exp), "znlog($a,$g,$n)");
  } else {
    is($k, $exp, "znlog($a,$g,$n)");
  }
}

sub is_one_of {
  my($n, @list) = @_;
  if (defined $n) {
    for (@list) {
      return 1 if defined $_ && "$n" eq "$_";
    }
  } else {
    for (@list) {
      return 1 if !defined $_;
    }
  }
  0;
}
