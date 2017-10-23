#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/permtonum/;

plan tests => 0
            + 7 # permtonum
            + 0;

###### permtonum
is(permtonum([]),0,"permtonum([])");
is(permtonum([0]),0,"permtonum([0])");
is(permtonum([1,0]),1,"permtonum([1,0])");
is(permtonum([6,3,4,2,5,0,1]),4768,"permtonum([6,3,4,2,5,0,1])");
is(permtonum([reverse(0..14),15..19]),"1790774578500738480","permtonum( 20 )");
is(permtonum([reverse(0..12),reverse(13..25)]),"193228515634198442606207999","permtonum( 26 )");
is(permtonum([16,10,18,0,33,25,34,30,8,13,11,27,22,26,21,17,4,12,29,5,39,31,1,6,7,14,32,28,19,15,23,24,9,20,38,35,3,36,37,2]),"331816865634235984753115643542286324917379756177","permtonum( 40 )");
