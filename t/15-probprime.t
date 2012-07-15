#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/is_prime is_prob_prime/;

plan tests => 0 + 6
                + 5
                + 5
                + 29
                + 22
                + 23
                + 16
                + 15
                + 28
                + 0;

# Some of these tests were inspired by Math::Primality's tests

ok(is_prob_prime(2)  == 2,  '2 is prime');
ok(is_prob_prime(1)  == 0,  '1 is not prime');
ok(is_prob_prime(0)  == 0,  '0 is not prime');
ok(is_prob_prime(-1) == 0, '-1 is not prime');
ok(is_prob_prime(-2) == 0, '-2 is not prime');
ok(is_prob_prime(20) == 0, '20 is not prime');

map { ok(!is_prob_prime($_), "A006945 number $_ is not prime") }
  qw/9 2047 1373653 25326001 3215031751/;

map { ok(!is_prob_prime($_), "A006945 number $_ is not prime") }
  qw/2152302898747 3474749660383 341550071728321 341550071728321 3825123056546413051/;

map { ok(!is_prob_prime($_), "Carmichael Number $_ is not prime") }
  qw/561 1105 1729 2465 2821 6601 8911 10585 15841 29341 41041 46657 52633
     62745 63973 75361 101101 340561 488881 852841 1857241 6733693
     9439201 17236801 23382529 34657141 56052361 146843929 216821881/;

map { ok(!is_prob_prime($_), "Pseudoprime (base 2) $_ is not prime" ) }
  qw/341 561 645 1105 1387 1729 1905 2047 2465 2701 2821 3277 4033 4369 4371
     4681 5461 6601 7957 8321 52633 88357/;

map { ok(!is_prob_prime($_), "Pseudoprime (base 3) $_ is not prime" ) }
  qw/121 703 1891 3281 8401 8911 10585 12403 16531 18721 19345 23521 31621
     44287 47197 55969 63139 74593 79003 82513 87913 88573 97567/;

map { ok(!is_prob_prime($_), "Pseudoprime (base 5) $_ is not prime" ) }
  qw/781 1541 5461 5611 7813 13021 14981 15751 24211 25351 29539 38081
     40501 44801 53971 79381/;

map { ok(is_prob_prime($_)==2, "Primegap start $_ is prime" ) }
  qw/2 3 7 23 89 113 523 887 1129 1327 9551 15683 19609 31397 155921/;

map { ok(is_prob_prime($_)==2, "Primegap end $_ is prime" ) }
  qw/5 11 29 97 127 541 907 1151 1361 9587 15727 19661 31469 156007 360749
     370373 492227 1349651 1357333 2010881 4652507 17051887 20831533 47326913
     122164969 189695893 191913031 10726905041/;
