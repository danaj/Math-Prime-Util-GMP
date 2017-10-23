#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/is_carmichael is_fundamental is_totient is_polygonal
                              polygonal_nth/;

plan tests => 0
            + 3 # is_carmichael
            + 4 # is_fundamental
            + 8 # is_totient
            + 4 # is_polygonal
            + 0;

###### is_carmichael
is_deeply( [grep { is_carmichael($_) } 1 .. 20000],
           [561,1105,1729,2465,2821,6601,8911,10585,15841],
           "Carmichael numbers to 20000" );
ok( is_carmichael("1298392318741906953539071949881"), "Large Carmichael" );
ok( is_carmichael("341627175004511735787409078802107169251"), "Larger Carmichael" );

###### is_fundamental
is_deeply( [grep { is_fundamental($_) } -50 .. 0],
           [-47,-43,-40,-39,-35,-31,-24,-23,-20,-19,-15,-11,-8,-7,-4,-3],
           "is_fundamental(-50 .. 0)" );
is_deeply( [grep { is_fundamental($_) } 0 .. 50],
           [1,5,8,12,13,17,21,24,28,29,33,37,40,41,44],
           "is_fundamental(0 .. 50)" );
is( is_fundamental("147573952589676412937"), 1, "is_fundamental(2^67+9)" );
is( is_fundamental("-147573952589676412884"), 1, "is_fundamental(-2^67+44)" );

###### is_totient
is_deeply( [map { is_totient($_) } 0..40],
           [0,1,1,0,1,0,1,0,1,0,1,0,1,0,0,0,1,0,1,0,1,0,1,0,1,0,0,0,1,0,1,0,1,0,0,0,1,0,0,0,1],
           "is_totient 0 .. 40" );
is_deeply( [grep { is_totient( 2**29 + $_ ) } 1 .. 80],
           [4,10,12,16,32,38,48,64,68,72],
           "is_fundamental(2^29_1 .. 2^29+80)" );
is( is_totient("9223372036854775836"), 1, "is_totient(2^63+28)" );
is( is_totient("9223372036854775828"), 1, "is_totient(2^63+20)" );
is( is_totient("9223372036854775832"), 0, "is_totient(2^63+24)" );

is( is_totient("9671406556917033397649496"), 1, "is_totient(2^83+88)" );
is( is_totient("9671406556917033397649458"), 0, "is_totient(2^83+50)" );
is( is_totient("9671406556917033397649472"), 1, "is_totient(2^83+64)" );

###### is_polygonal
is_deeply( [grep { is_polygonal($_,3) } 1..55], [1,3,6,10,15,21,28,36,45,55], "first 10 triangular numbers" );
is_deeply( [grep { is_polygonal($_,23) } 1..955], [1,23,66,130,215,321,448,596,765,955], "first 10 23-gonal numbers" );
is( polygonal_nth("140737496743936",3), 16777216, "140737496743936 is the 16777216-th triangular number");
is( polygonal_nth("228623681298582551246684960361911294205",5), "12345678901234567890", "identified the 12345678901234567890-th pentagonal number");
