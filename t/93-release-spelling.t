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
                 subfactorial multifactorials
                 gcd lcm kronecker invmod exp
                 irand irand64 drand drand64 urandomm urandomb
                 factorialmod hammingweight numtoperm permtonum
                 semiprime semiprimes coprime k-tuples
                 precalculated premultiplier multipoint
                 pseudoprime pseudoprimes compositeness
                 powersum
                 signint cmpint cmpabsint
                 setbit clrbit notbit tstbit
                 bitand bitor bitxor bitnot
                 p-adic bitwise
                 0-th -th
                 -inf
                 von
                 PSP-2
                 pp/);

all_pod_files_spelling_ok();
