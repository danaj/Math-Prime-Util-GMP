#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/permtonum numtoperm/;

plan tests => 0
            + 7 # permtonum
            + 7 # numtoperm
            + 0;

###### permtonum
is(permtonum([]),0,"permtonum([])");
is(permtonum([0]),0,"permtonum([0])");
is(permtonum([1,0]),1,"permtonum([1,0])");
is(permtonum([6,3,4,2,5,0,1]),4768,"permtonum([6,3,4,2,5,0,1])");
is(permtonum([reverse(0..14),15..19]),"1790774578500738480","permtonum( 20 )");
is(permtonum([reverse(0..12),reverse(13..25)]),"193228515634198442606207999","permtonum( 26 )");
is(permtonum([16,10,18,0,33,25,34,30,8,13,11,27,22,26,21,17,4,12,29,5,39,31,1,6,7,14,32,28,19,15,23,24,9,20,38,35,3,36,37,2]),"331816865634235984753115643542286324917379756177","permtonum( 40 )");

is_deeply([numtoperm(0,0)],[],"numtoperm(0,0)");
is_deeply([numtoperm(1,0)],[0],"numtoperm(1,0)");
is_deeply([numtoperm(1,1)],[0],"numtoperm(1,1)");
is_deeply([numtoperm(5,15)],[0,3,2,4,1],"numtoperm(5,15)");
is_deeply([numtoperm(5,-2)],[4,3,2,0,1],"numtoperm(5,-2)");
is_deeply([numtoperm(24,987654321)],[0,1,2,3,4,5,6,7,8,9,10,13,11,21,14,20,17,15,12,22,18,19,23,16],"numtoperm(24,987654321)");
is_deeply([numtoperm(40,"168413612234311889416542084240799668239048237663")],[8,11,0,1,35,34,29,38,39,3,32,13,33,25,4,10,17,16,6,2,30,37,12,24,28,36,22,26,18,19,20,23,31,9,27,5,21,15,14,7],"numtoperm(40,...)");
