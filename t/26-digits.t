#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/todigits/;
use Math::BigInt;  # Don't use GMP so we don't have to work around bug

plan tests =>  0
            + 12
            +  1
            +  0;

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

