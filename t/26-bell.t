#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/bell_number/;

my @bell = (1,1,2,5,15,52,203,877,4140,21147,115975);

plan tests => 0
            + 1   # small values
            + 2   # larger values
            + 0;

###### bell_number
is_deeply( [map { bell_number($_) } 0 .. $#bell], \@bell, "bell_number(0..$#bell)" );
is( bell_number(20), "51724158235372",             "bell_number(20)" );
is( bell_number(30), "846749014511809332450147",   "bell_number(30)" );
