#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/todigits fromdigits/;
use Math::BigInt;  # Don't use GMP so we don't have to work around bug

plan tests =>  0
            + 13 + 1      # todigits
            +  6 + 4 + 22 # fromdigits
            +  1          # combined
            +  0;

###### todigits
is_deeply([todigits(0)], [], "todigits 0");
is_deeply([todigits(1)], [1], "todigits 1");
is_deeply([todigits(77)], [7,7], "todigits 77");
is_deeply([todigits(77,2)], [1,0,0,1,1,0,1], "todigits 77 base 2");
is_deeply([todigits(77,3)], [2,2,1,2], "todigits 77 base 3");
is_deeply([todigits(77,21)], [3,14], "todigits 77 base 21");
is_deeply([todigits(900,2)], [1,1,1,0,0,0,0,1,0,0], "todigits 900 base 2");
is_deeply([todigits(900,2,0)], [], "todigits 900 base 2 len 0");
is_deeply([todigits(900,2,3)], [1,0,0], "todigits 900 base 2 len 3");
is_deeply([todigits(900,2,32)], [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,1,0,0], "todigits 900 base 2 len 32");
is_deeply([todigits(58127,16)], [14,3,0,15], "todigits 58127 base 16");
is_deeply([todigits(6345354,10,4)], [5,3,5,4], "todigits 6345354 base 10 len 4");
is_deeply([todigits("188661215071572375748916455613",503)], [181, 488, 270, 406, 138, 112, 263, 156, 399, 236, 416], "todigits 30-digit base 503");

is_deeply([todigits(-24)], [2,4], "todigits ignores sign");

###### fromdigits
is(fromdigits([]), 0, "fromdigits([]) = 0");
is(fromdigits([1]), 1, "fromdigits([1]) = 1");
is(fromdigits([1,0,1],2), 5, "101 base 2 = 5");
is(fromdigits([1,1,2,1,2,0,2,0,1,0,1,1,1,2,0],3), 7749393, "fromdigits of 7749393 in base 3");
is(fromdigits([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,1,0,0],2), 900, "handle leading zeros");
is(fromdigits([14,3,0,15],16), 58127, "fromdigits of 58127 base 16");

is(fromdigits(""), 0, "fromdigits empty string returns 0");
is(fromdigits("1f",16), 31, "fromdigits hex string");
is(fromdigits("24"), 24, "fromdigits decimal");
is(fromdigits("zzzyzzzyzzzyzzzy",36), "7958656371562241451187966", "fromdigits with Large base 36 number");
is(fromdigits("101", undef), 101, "fromdigits with undef base defaults to 10");
is(fromdigits([7,999,44],5), 5214, "fromdigits array entries can exceed base");
is(fromdigits([1,-2,3],10), 83, "fromdigits array entries can be negative");
is(fromdigits([-1,2],10), -8, "fromdigits array result can be negative");
is("".fromdigits(["18446744073709551616",3],10), "184467440737095516163", "fromdigits bigint array entry");
{
  my $bigbase = Math::BigInt->new("18446744073709551616");
  my $expect = $bigbase->copy->bmul($bigbase)
               ->badd($bigbase->copy->bmul(2))
               ->badd(3);
  is("".fromdigits([1,2,3], "$bigbase"), "$expect", "fromdigits bigint base string");
  is("".fromdigits([1,2,3], $bigbase), "$expect", "fromdigits bigint base object");
}
is(fromdigits(Math::BigInt->new(123),10), 123, "fromdigits bigint object as digit string");
my @bad_fromdigits = (
  [ sub { fromdigits(undef) },          qr/Parameter must be defined/, "undef input" ],
  [ sub { fromdigits({}) },             qr/string or array reference/, "hash ref input" ],
  [ sub { fromdigits([1,undef,2],10) }, qr/Parameter must be defined/, "undef array entry" ],
  [ sub { my @d = (1,2,3); delete $d[1]; fromdigits(\@d,10) },
                                            qr/Parameter must be defined/, "sparse array entry" ],
  [ sub { fromdigits([1,"foo",2],10) }, qr/must be an integer/,       "non-integer array entry" ],
  [ sub { fromdigits("-101",2) },       qr/invalid digit/,            "negative string" ],
  [ sub { fromdigits("+101",2) },       qr/invalid digit/,            "positive string" ],
  [ sub { fromdigits("1!",1000) },      qr/invalid digit/,            "punctuation in high base" ],
  [ sub { fromdigits("2",2) },          qr/invalid digit/,            "digit equal to base" ],
  [ sub { fromdigits("101",1) },        qr/invalid base/,             "base one" ],
  [ sub { fromdigits("101","-0") },     qr/invalid base/,             "signed zero base" ],
  [ sub { fromdigits("101","-00") },    qr/invalid base/,             "signed zero base with leading zeros" ],
  [ sub { fromdigits("101",-2) },       qr/non-negative/,             "negative base" ],
  [ sub { fromdigits("101","-01") },    qr/non-negative/,             "negative base with leading zero" ],
);
foreach my $t (@bad_fromdigits) {
  my ($sub, $err, $name) = @$t;
  my $ok = eval { $sub->(); 1 };
  like($@, $err, "fromdigits croaks on $name") if !$ok;
  fail("fromdigits croaks on $name") if $ok;
}

###### more from/to
is(fromdigits([todigits(56,2,8)],2), 56, "fromdigits of previous");
