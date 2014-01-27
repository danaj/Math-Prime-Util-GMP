#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/factor is_prime/;

plan tests => 0 + 57
                + 24
                + 2
                + 6    # individual tets for factoring methods
                + 7*7  # factor extra tests
                + 8    # factor in scalar context
                + 0;

# On a 64-bit machine, put all 32-bit nums in /tmp/foo, 64-bit in /tmp/foo2
#   for i in `sort -n /tmp/foo`; do perl -MMath::Factor::XS=:all -E "say 'is_deeply( [ factor(', $i, ') ], [', join(',', prime_factors($i)), '], \"factor($i)\" );';"; done
#   for i in `sort -n /tmp/foo2`; do perl -MMath::Factor::XS=:all -E "say 'is_deeply( [ factor(\'', $i, '\') ], [', join(',', prime_factors($i)), '], \"factor($i)\" );';"; done
#
# For the latter, check every factor to make sure it fits in 32-bit, quote if
# not.  Run anything larger than 64-bit through Yafu or Pari.
#
# The obvious point here is that we shouldn't generate tests using our own code,
# unless we want to hand verify each case (admittedly not that hard).
#
#diag "factoring 32-bit numbers";
is_deeply( [ factor(0) ], [0], "factor(0)" );
is_deeply( [ factor(1) ], [1], "factor(1)" );
is_deeply( [ factor(2) ], [2], "factor(2)" );
is_deeply( [ factor(3) ], [3], "factor(3)" );
is_deeply( [ factor(4) ], [2,2], "factor(4)" );
is_deeply( [ factor(5) ], [5], "factor(5)" );
is_deeply( [ factor(6) ], [2,3], "factor(6)" );
is_deeply( [ factor(7) ], [7], "factor(7)" );
is_deeply( [ factor(8) ], [2,2,2], "factor(8)" );
is_deeply( [ factor(16) ], [2,2,2,2], "factor(16)" );
is_deeply( [ factor(30) ], [2,3,5], "factor(30)" );
is_deeply( [ factor(57) ], [3,19], "factor(57)" );
is_deeply( [ factor(64) ], [2,2,2,2,2,2], "factor(64)" );
is_deeply( [ factor(210) ], [2,3,5,7], "factor(210)" );
is_deeply( [ factor(377) ], [13,29], "factor(377)" );
is_deeply( [ factor(403) ], [13,31], "factor(403)" );
is_deeply( [ factor(629) ], [17,37], "factor(629)" );
is_deeply( [ factor(779) ], [19,41], "factor(779)" );
is_deeply( [ factor(808) ], [2,2,2,101], "factor(808)" );
is_deeply( [ factor(989) ], [23,43], "factor(989)" );
is_deeply( [ factor(1363) ], [29,47], "factor(1363)" );
is_deeply( [ factor(2310) ], [2,3,5,7,11], "factor(2310)" );
is_deeply( [ factor(2727) ], [3,3,3,101], "factor(2727)" );
is_deeply( [ factor(9592) ], [2,2,2,11,109], "factor(9592)" );
is_deeply( [ factor(12625) ], [5,5,5,101], "factor(12625)" );
is_deeply( [ factor(30030) ], [2,3,5,7,11,13], "factor(30030)" );
is_deeply( [ factor(30107) ], [7,11,17,23], "factor(30107)" );
is_deeply( [ factor(34643) ], [7,7,7,101], "factor(34643)" );
is_deeply( [ factor(78498) ], [2,3,3,7,7,89], "factor(78498)" );
is_deeply( [ factor(134431) ], [11,11,11,101], "factor(134431)" );
is_deeply( [ factor(221897) ], [13,13,13,101], "factor(221897)" );
is_deeply( [ factor(496213) ], [17,17,17,101], "factor(496213)" );
is_deeply( [ factor(510510) ], [2,3,5,7,11,13,17], "factor(510510)" );
is_deeply( [ factor(664579) ], [664579], "factor(664579)" );
is_deeply( [ factor(692759) ], [19,19,19,101], "factor(692759)" );
is_deeply( [ factor(1228867) ], [23,23,23,101], "factor(1228867)" );
is_deeply( [ factor(2214143) ], [1487,1489], "factor(2214143)" );
is_deeply( [ factor(2463289) ], [29,29,29,101], "factor(2463289)" );
is_deeply( [ factor(3008891) ], [31,31,31,101], "factor(3008891)" );
is_deeply( [ factor(5115953) ], [37,37,37,101], "factor(5115953)" );
is_deeply( [ factor(5761455) ], [3,5,7,37,1483], "factor(5761455)" );
is_deeply( [ factor(6961021) ], [41,41,41,101], "factor(6961021)" );
is_deeply( [ factor(8030207) ], [43,43,43,101], "factor(8030207)" );
is_deeply( [ factor(9699690) ], [2,3,5,7,11,13,17,19], "factor(9699690)" );
is_deeply( [ factor(10486123) ], [47,47,47,101], "factor(10486123)" );
is_deeply( [ factor(10893343) ], [1327,8209], "factor(10893343)" );
is_deeply( [ factor(12327779) ], [1627,7577], "factor(12327779)" );
is_deeply( [ factor(50847534) ], [2,3,3,3,19,49559], "factor(50847534)" );
is_deeply( [ factor(114256942) ], [2,57128471], "factor(114256942)" );
is_deeply( [ factor(223092870) ], [2,3,5,7,11,13,17,19,23], "factor(223092870)" );
is_deeply( [ factor(455052511) ], [97,331,14173], "factor(455052511)" );
is_deeply( [ factor(547308031) ], [547308031], "factor(547308031)" );
is_deeply( [ factor(701737021) ], [25997,26993], "factor(701737021)" );
is_deeply( [ factor(999999929) ], [999999929], "factor(999999929)" );
is_deeply( [ factor(2147483647) ], [2147483647], "factor(2147483647)" );
is_deeply( [ factor(4118054813) ], [19,216739727], "factor(4118054813)" );
is_deeply( [ factor(4294967293) ], [9241,464773], "factor(4294967293)" );

#diag "factoring 64-bit numbers";
is_deeply( [ factor('6469693230') ], [2,3,5,7,11,13,17,19,23,29], "factor(6469693230)" );
is_deeply( [ factor('17179869172') ], [2,2,9241,464773], "factor(17179869172)" );
is_deeply( [ factor('37607912018') ], [2,'18803956009'], "factor(37607912018)" );
is_deeply( [ factor('200560490130') ], [2,3,5,7,11,13,17,19,23,29,31], "factor(200560490130)" );
is_deeply( [ factor('346065536839') ], [11,11,163,373,47041], "factor(346065536839)" );
is_deeply( [ factor('600851475143') ], [71,839,1471,6857], "factor(600851475143)" );
is_deeply( [ factor('3204941750802') ], [2,3,3,3,11,277,719,27091], "factor(3204941750802)" );
is_deeply( [ factor('7420738134810') ], [2,3,5,7,11,13,17,19,23,29,31,37], "factor(7420738134810)" );
is_deeply( [ factor('29844570422669') ], [19,19,27259,3032831], "factor(29844570422669)" );
is_deeply( [ factor('279238341033925') ], [5,5,7,13,194899,629773], "factor(279238341033925)" );
is_deeply( [ factor('304250263527210') ], [2,3,5,7,11,13,17,19,23,29,31,37,41], "factor(304250263527210)" );
is_deeply( [ factor('2623557157654233') ], [3,113,136841,56555467], "factor(2623557157654233)" );
is_deeply( [ factor('9007199254740991') ], [6361,69431,20394401], "factor(9007199254740991)" );
is_deeply( [ factor('9007199254740992') ], [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2], "factor(9007199254740992)" );
is_deeply( [ factor('9007199254740993') ], [3,107,'28059810762433'], "factor(9007199254740993)" );
is_deeply( [ factor('9999986200004761') ], [99999931,99999931], "factor(9999986200004761)" );
is_deeply( [ factor('13082761331670030') ], [2,3,5,7,11,13,17,19,23,29,31,37,41,43], "factor(13082761331670030)" );
is_deeply( [ factor('24739954287740860') ], [2,2,5,7,1123,'157358823863'], "factor(24739954287740860)" );
is_deeply( [ factor('99999989237606677') ], [316227731,316227767], "factor(99999989237606677)" );
is_deeply( [ factor('614889782588491410') ], [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47], "factor(614889782588491410)" );
is_deeply( [ factor('999999866000004473') ], [999999929,999999937], "factor(999999866000004473)" );
is_deeply( [ factor('3369738766071892021') ], [204518747,'16476429743'], "factor(3369738766071892021)" );
is_deeply( [ factor('10023859281455311421') ], ['1308520867','7660450463'], "factor(10023859281455311421)" );
is_deeply( [ factor('18446744073709551611') ], [11,59,'98818999','287630261'], "factor(18446744073709551611)" );

# Check perfect squares that make it past early testing
is_deeply( [ factor('1524157875323973084894790521049') ], ['1234567890123493','1234567890123493'], "factor(1234567890123493^2)" );
is_deeply( [ factor('823543') ], [qw/7 7 7 7 7 7 7/], "factor 7^7" );
# A good test, but it's slow:
# is_deeply( [ factor('148675665359980297048795508874724049089546782584077934753925649') ], ['1234567890123493', '1234567890123493', '9876543210987701', '9876543210987701'], "factor(1234567890123493^2 * 9876543210987701^2)" );


#diag "factor 105-bit number with p-1";
Math::Prime::Util::GMP::_GMP_set_verbose(4);
is_deeply( [ sort {$a<=>$b} Math::Prime::Util::GMP::pminus1_factor('22095311209999409685885162322219') ], ['3916587618943361', '5641469912004779'], "p-1 factors 22095311209999409685885162322219" );
Math::Prime::Util::GMP::_GMP_set_verbose(0);

is_deeply( [ sort {$a<=>$b} Math::Prime::Util::GMP::pplus1_factor('22095311209999409685885162322219') ], ['3916587618943361', '5641469912004779'], "p+1 factors 22095311209999409685885162322219" );

is_deeply( [ sort {$a<=>$b} Math::Prime::Util::GMP::ecm_factor('16049407357301026788959025956634678743968244330856613525782006075043') ], [qw/99151111 161868154531329727500068314480456792299263740280798402004613/], "ECM factors p8*p60" );

is_deeply( [ sort {$a<=>$b} Math::Prime::Util::GMP::qs_factor('22095311209999409685885162322219') ], ['3916587618943361', '5641469912004779'], "QS factors 22095311209999409685885162322219" );

#diag "factor 736-bit number with HOLF";
is_deeply( [ sort {$a<=>$b} Math::Prime::Util::GMP::holf_factor('185486767418172501041516225455805768237366368964328490571098416064672288855543059138404131637447372942151236559829709849969346650897776687202384767704706338162219624578777915220190863619885201763980069247978050169295918863') ], ['192606732705880508138303165129171270891951231683030125996296974238495711578947569589234612013165893468683239489', '963033663529402540691515825645856354459756158415150629981484871192478557894737847946173060065829467343416197967'], "HOLF factors poorly formed 222-digit semiprime" );

# Test stage 2 of pminus1
is_deeply( [ sort {$a<=>$b} Math::Prime::Util::GMP::pminus1_factor('23113042053749572861737011', 100, 100000) ], ['694059980329', '33301217054459'], "p-1 factors 23113042053749572861737011 in stage 2");

#diag "extra tests for each method";
extra_factor_test("prho_factor",   sub {Math::Prime::Util::GMP::prho_factor(shift)});
extra_factor_test("pbrent_factor", sub {Math::Prime::Util::GMP::pbrent_factor(shift)});
extra_factor_test("pminus1_factor",sub {Math::Prime::Util::GMP::pminus1_factor(shift)});
extra_factor_test("pplus1_factor", sub {Math::Prime::Util::GMP::pplus1_factor(shift)});
extra_factor_test("holf_factor",   sub {Math::Prime::Util::GMP::holf_factor(shift)});
extra_factor_test("squfof_factor", sub {Math::Prime::Util::GMP::squfof_factor(shift)});
extra_factor_test("ecm_factor",    sub {Math::Prime::Util::GMP::ecm_factor(shift)});

sub extra_factor_test {
  my $fname = shift;
  my $fsub = shift;

  is_deeply( [ sort {$a<=>$b} $fsub->(0)   ], [0],       "$fname(0)" );
  is_deeply( [ sort {$a<=>$b} $fsub->(1)   ], [1],       "$fname(1)" );
  is_deeply( [ sort {$a<=>$b} $fsub->(2)   ], [2],       "$fname(2)" );
  is_deeply( [ sort {$a<=>$b} $fsub->(13)  ], [13],      "$fname(13)" );
  is_deeply( [ sort {$a<=>$b} $fsub->(403) ], [13, 31],  "$fname(403)" );
  is_deeply( [ sort {$a<=>$b} $fsub->(53936983) ], [7013, 7691],  "$fname(53936983)" );
  is_deeply( [ sort {$a<=>$b} $fsub->('1754012594703269855671') ], ['41110234981', '42666080491'],  "$fname(1754012594703269855671)" );
}

# Factor in scalar context
is( scalar factor(0), 1, "scalar factor(0) should be 1" );
is( scalar factor(1), 1, "scalar factor(1) should be 1" );
is( scalar factor(3), 1, "scalar factor(3) should be 1" );
is( scalar factor(4), 2, "scalar factor(4) should be 2" );
is( scalar factor(5), 1, "scalar factor(5) should be 1" );
is( scalar factor(6), 2, "scalar factor(6) should be 2" );
is( scalar factor(30107), 4, "scalar factor(30107) should be 4" );
is( scalar factor(174636000), 15, "scalar factor(174636000) should be 15" );
