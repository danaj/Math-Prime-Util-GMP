#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/todigits fromdigits/;
use Math::BigInt;  # Don't use GMP so we don't have to work around bug

plan tests =>  0
            + 12 + 1      # todigits
            +  6 + 4      # fromdigits
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

###### more from/to
is(fromdigits([todigits(56,2,8)],2), 56, "fromdigits of previous");
