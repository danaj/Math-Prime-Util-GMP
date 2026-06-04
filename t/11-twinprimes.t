#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/twin_primes/;

my @small_twins = qw/
3 5 11 17 29 41 59 71 101 107 137 149 179 191 197 227 239 269 281 311
347 419 431 461 521 569 599 617 641 659 809 821 827 857 881 1019 1031
1049 1061 1091 1151 1229 1277 1289 1301 1319 1427 1451 1481 1487 1607
/;

my %small_range = (
  "0 to 0" => [],
  "0 to 2" => [],
  "0 to 3" => [3],
  "0 to 5" => [3,5],
  "1 to 11" => [3,5,11],
  "2 to 11" => [3,5,11],
  "3 to 11" => [3,5,11],
  "4 to 11" => [5,11],
  "5 to 10" => [5],
  "5 to 11" => [5,11],
  "5 to 12" => [5,11],
  "6 to 10" => [],
  "11 to 11" => [11],
  "12 to 17" => [17],
  "29 to 29" => [29],
  "29 to 31" => [29],
  "30 to 31" => [],
  "31 to 29" => [],
  "213897 to 213997" => [213947],
  "134217228 to 134217728" => [134217401,134217437],
  "4294957296 to 4294957796" => [4294957307,4294957397,4294957697],
);

plan tests => 12 + 3 + 1 + scalar(keys %small_range) + 2;

ok(!eval { twin_primes(undef); },   "twin_primes(undef)");
ok(!eval { twin_primes("a"); },     "twin_primes(a)");
ok(!eval { twin_primes(-4); },      "twin_primes(-4)");
ok(!eval { twin_primes(2,undef); }, "twin_primes(2,undef)");
ok(!eval { twin_primes(2,'x'); },   "twin_primes(2,x)");
ok(!eval { twin_primes(2,-4); },    "twin_primes(2,-4)");
ok(!eval { twin_primes(undef,7); }, "twin_primes(undef,7)");
ok(!eval { twin_primes('x',7); },   "twin_primes(x,7)");
ok(!eval { twin_primes(-10,7); },   "twin_primes(-10,7)");
ok(!eval { twin_primes(undef,undef); },  "twin_primes(undef,undef)");
ok(!eval { twin_primes('x','x'); }, "twin_primes(x,x)");
ok(!eval { twin_primes(-10,-4); },  "twin_primes(-10,-4)");
is_deeply( twin_primes("-0", 5), [3,5], 'twin_primes("-0",5)' );
is_deeply( twin_primes("-00", 5), [3,5], 'twin_primes("-00",5)' );
ok(!eval { twin_primes("-01", 5) }, 'twin_primes("-01",5)');

is_deeply( twin_primes($small_twins[-1]), \@small_twins,
           "twin_primes($small_twins[-1])" );

while (my($range, $expect) = each (%small_range)) {
  my($low,$high) = $range =~ /(\d+) to (\d+)/;
  is_deeply( twin_primes($low, $high), $expect,
             "twin_primes($low,$high) should return [@{$expect}]");
}

is_deeply( twin_primes("1000000000000000000000000000000",
                       "1000000000000000000000000020000"),
           [qw/1000000000000000000000000001681
               1000000000000000000000000004831
               1000000000000000000000000018739
               1000000000000000000000000019171/],
           "twin_primes 10^30 10^30+20000" );
is_deeply( twin_primes("1000000000000000000000000004832",
                       "1000000000000000000000000018738"),
           [],
           "twin_primes empty large range" );
