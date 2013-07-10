#!/usr/bin/env perl
use strict;
use warnings;

BEGIN {
  unless ($ENV{RELEASE_TESTING}) {
    require Test::More;
    Test::More::plan(skip_all => 'these tests are for release candidate testing');
  }
}

#---------------------------------------------------------------------


use Test::More;
eval "use Test::Spellunker";
plan skip_all => "Test::Spellunker required for testing POD spelling" if $@;

add_stopwords(qw/bigint bigints bignum bignums primorial
                 pseudoprime pseudoprimes
                 semiprime semiprimes
                 precalculated premultiplier
                 PSP-2
                 pp/);

all_pod_files_spelling_ok();
