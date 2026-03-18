#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/catalan_number/;

my @catalan = (1,1,2,5,14,42,132,429,1430,4862,16796);

plan tests => 0
            + 1      # small values
            + 2      # larger values
            + 0;

###### catalan_number
is_deeply( [map { catalan_number($_) } 0 .. $#catalan], \@catalan, "catalan_number(0..$#catalan)" );
is( catalan_number(20), "6564120420",      "catalan_number(20)" );
is( catalan_number(30), "3814986502092304","catalan_number(30)" );
