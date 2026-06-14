#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/sqrtmod allsqrtmod rootmod allrootmod powmod invmod/;

sub is_one_of {
  my($n, @list) = @_;
  return 0 unless defined $n;
  for my $v (@list) {
    return 1 if "$n" eq "$v";
  }
  return 0;
}

my @allsqrt = (
  [ 0, 0, [] ],
  [ 2, 1, [0] ],
  [ 3, 9, [] ],
  [ 1, 8, [1,3,5,7] ],
  [ 5, 404, [45,157,247,359] ],
  [ 9, 400, [3,53,147,197,203,253,347,397] ],
  [ 0, 36, [0,6,12,18,24,30] ],
  [ 4, 16, [2,6,10,14] ],
  [ "9223372036854775808", "5675921253449092823",
    ["22172359690642254","5653748893758450569"] ],
);

my @allroot = (
  [ [2,1,0], [] ],
  [ [2,-1,1], [0] ],
  [ [2,0,1], [0] ],
  [ [1,0,7], [0,1,2,3,4,5,6] ],
  [ [2,0,7], [] ],
  [ [0,-2,7], [] ],
  [ [2,-1,4], [] ],
  [ [0,1,7], [0] ],
  [ [2,1,7], [2] ],
  [ [2,3,7], [] ],
  [ [1,2,8], [1,3,5,7] ],
  [ [13,6,107], [24,83] ],
  [ [13,-6,107], [49,58] ],
  [ [5,3,13], [7,8,11] ],
  [ [53,3,151], [15,27,109] ],
  [ [581,5,151], [34,42,43,62,121] ],
  [ [8,3,27], [2,11,20] ],
  [ [22,3,1505], [148,578,673,793,813,1103,1243,1318,1458] ],
  [ [1,4,20], [1,3,7,9,11,13,17,19] ],
  [ [96,5,128], [6,14,22,30,38,46,54,62,70,78,86,94,102,110,118,126] ],
  [ [26,5,625], [81,206,331,456,581] ],
  [ ["9833625071",3,"10000000071"], [qw/3333332807 6666666164 9999999521/] ],
);

plan tests => 0
            + 2 * @allsqrt
            + 2 * @allroot
            + 7;

for my $t (@allsqrt) {
  my($a, $n, $exp) = @$t;
  is_deeply([map {"$_"} allsqrtmod($a,$n)], [map {"$_"} @$exp],
            "allsqrtmod($a,$n)");
  is(scalar allsqrtmod($a,$n), scalar(@$exp),
     "scalar allsqrtmod($a,$n)");
}

for my $t (@allroot) {
  my($args, $exp) = @$t;
  my($a, $k, $n) = @$args;
  is_deeply([map {"$_"} allrootmod($a,$k,$n)], [map {"$_"} @$exp],
            "allrootmod($a,$k,$n)");
  is(scalar allrootmod($a,$k,$n), scalar(@$exp),
     "scalar allrootmod($a,$k,$n)");
}

is(rootmod(0,5,0), undef, "rootmod(a,k,0) is undef");
is(rootmod(3,5,1), 0, "rootmod(a,k,1) is 0");
is(rootmod(1,0,17), 1, "rootmod(1,0,17)");
is(rootmod(2,0,17), undef, "rootmod(2,0,17)");
is(rootmod(2,-11,4725), 4412, "rootmod negative k");
is(rootmod(12,41,1147), 1106, "rootmod one-root composite case");

{
  my $r = rootmod(13,6,107);
  ok(is_one_of($r, 24, 83) && powmod($r, 6, 107) == 13,
     "rootmod returns a valid root");
}
