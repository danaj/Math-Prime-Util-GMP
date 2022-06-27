#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/is_tau/;

my @ok = (
	[ 1, 1 ], [ 2, 2 ], [ 4, 3 ], [ 6, 4 ], [ 8, 4 ], [ 18, 6 ], [ 12, 6 ],
	[ 510510, 128 ],
	[ '10868740069638250502059754282498', 30 ], # 2.17^4.8066344400287^2
);
my @nok = (
	[ 0, 1 ], [ 0, 2 ], [ 1, 2 ], [ 256, 3 ], [ 121, 4 ], [ 121, 9 ],
	[ 125, 16 ],
);

plan tests => @ok + @nok;

for (@ok) {
	my($n, $k) = @$_;
	is(is_tau($n, $k), 1, "tau($n) == $k");
}
for (@nok) {
	my($n, $k) = @$_;
    is(is_tau($n, $k), 0, "tau($n) != $k");
}

done_testing();
