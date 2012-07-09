#!/usr/bin/env perl
use strict;
use warnings;
use Math::Prime::Util::GMP;

use Test::More  tests => 1;

my @functions =  qw(
                     is_strong_pseudoprime
                   );
can_ok( 'Math::Prime::Util::GMP', @functions);
