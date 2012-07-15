#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/factor is_prime/;

my $extra = defined $ENV{RELEASE_TESTING} && $ENV{RELEASE_TESTING};


my @testn = qw/0 1 2 3 4 5 6 7 8 16 57 64 377 9592 30107 78498 664579 5761455
               114256942 2214143 999999929 50847534 455052511 2147483647
               4118054813
               30 210 2310 30030 510510 9699690 223092870
               1363 989 779 629 403
               547308031
               808 2727 12625 34643 134431 221897 496213 692759 1228867
               2463289 3008891 5115953 6961021 8030207 10486123
               10893343 12327779 701737021
              /;

my @testn64 = qw/37607912018 346065536839 600851475143
                 3204941750802 29844570422669
                 279238341033925 2623557157654233 24739954287740860
                 3369738766071892021 10023859281455311421
                 9007199254740991 9007199254740992 9007199254740993
                 6469693230 200560490130 7420738134810 304250263527210
                 13082761331670030 614889782588491410
                /;

push @testn, @testn64;

push @testn, qw/9999986200004761 99999989237606677 999999866000004473/;

plan tests =>  (2 * scalar @testn) + 2 + 6*5;

foreach my $n (@testn) {
  my @f = factor($n);
  my $facstring = join(' * ',@f);

  # Do they multiply to the number?
  my $product = 1;  $product *= $_ for @f;
  is( $product, $n, "$n = [ $facstring ]" );

  # Are they all prime?
  my $isprime = 1; $isprime *= is_prime($_) for @f;
  if ($n < 2) {
    ok( !$isprime, "All factors [ $facstring ] of $n are not prime" );
  } else {
    ok( $isprime, "All factors [ $facstring ] of $n are prime" );
  }
};

is_deeply( [ sort {$a<=>$b} Math::Prime::Util::GMP::pminus1_factor('22095311209999409685885162322219') ], ['3916587618943361', '5641469912004779'], "p-1 factors 22095311209999409685885162322219" );

is_deeply( [ sort {$a<=>$b} Math::Prime::Util::GMP::holf_factor('185486767418172501041516225455805768237366368964328490571098416064672288855543059138404131637447372942151236559829709849969346650897776687202384767704706338162219624578777915220190863619885201763980069247978050169295918863') ], ['192606732705880508138303165129171270891951231683030125996296974238495711578947569589234612013165893468683239489', '963033663529402540691515825645856354459756158415150629981484871192478557894737847946173060065829467343416197967'], "HOLF factors poorly formed 222-digit semiprime" );


extra_factor_test("prho_factor",   sub {Math::Prime::Util::GMP::prho_factor(shift)});
extra_factor_test("pbrent_factor", sub {Math::Prime::Util::GMP::pbrent_factor(shift)});
extra_factor_test("pminus1_factor",sub {Math::Prime::Util::GMP::pminus1_factor(shift)});
extra_factor_test("holf_factor",   sub {Math::Prime::Util::GMP::holf_factor(shift)});
extra_factor_test("squfof_factor", sub {Math::Prime::Util::GMP::squfof_factor(shift)});

sub extra_factor_test {
  my $fname = shift;
  my $fsub = shift;

  is_deeply( [ sort {$a<=>$b} $fsub->(0)   ], [0],       "$fname(0)" );
  is_deeply( [ sort {$a<=>$b} $fsub->(1)   ], [1],       "$fname(1)" );
  is_deeply( [ sort {$a<=>$b} $fsub->(2)   ], [2],       "$fname(2)" );
  is_deeply( [ sort {$a<=>$b} $fsub->(403) ], [13, 31],  "$fname(403)" );
  is_deeply( [ sort {$a<=>$b} $fsub->(53936983) ], [7013, 7691],  "$fname(53936983)" );
  is_deeply( [ sort {$a<=>$b} $fsub->('1754012594703269855671') ], ['41110234981', '42666080491'],  "$fname(1754012594703269855671)" );
}

