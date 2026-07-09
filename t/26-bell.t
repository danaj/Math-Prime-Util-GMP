#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/bell_number fubini/;

my @bell = (1,1,2,5,15,52,203,877,4140,21147,115975);
my @fubini = (1,1,3,13,75,541,4683,47293,545835,7087261,102247563);

plan tests => 0
            + 1   # small values
            + 2   # larger values
            + 1   # small fubini values
            + 2   # larger fubini values
            + 0;

###### bell_number
is_deeply( [map { bell_number($_) } 0 .. $#bell], \@bell, "bell_number(0..$#bell)" );
is( bell_number(20), "51724158235372",             "bell_number(20)" );
is( bell_number(30), "846749014511809332450147",   "bell_number(30)" );

###### fubini
is_deeply( [map { fubini($_) } 0 .. $#fubini], \@fubini, "fubini(0..$#fubini)" );
is( fubini(20), "2677687796244384203115",        "fubini(20)" );
is( fubini(30), "11403568794011880483742464196184901963", "fubini(30)" );
