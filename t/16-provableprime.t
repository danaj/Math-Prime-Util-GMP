#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/is_provable_prime is_provable_prime_with_cert
                              is_trial_prime
                              is_llr_prime is_proth_prime
                              is_aks_prime is_miller_prime is_ecpp_prime
                              is_nminus1_prime is_nplus1_prime is_bls75_prime/;

my @llrs = (
  [202, 0],
  [159807057, 0],
  [10000000019,-1],
  [10051583, 2],
  [10072063, 2],
  [2097151, 0],
  [2147483647, 2],
  [805306367, 0],
  ["26388279066623", 2],
  [1064959, 0],
  [1114111, 2],
  [4349951, 0],
  [4374527, 2],
);
my @prs = (
  [10072063, -1],
  [8642561, 0],
  [8650753, 2],
  [16785409, 0],  # Perfect square
  [22560769, 0],
  [56770561, 2],
  ["38335150030849", 0],    # Harder to find
  ["47270099681281", 0],
  ["302111489261569", 0],
  ["2372913730682881", 2],
  ["208987568115548161", 2],
  ["19578524666953729", 2],
);
my @ecpps = (
  ["340282366920938463463374607431768211507", 1],
  ["4546500098776576231268807308545439", 1],
  ["12985198116842947666516311049464592230676113912623", 1],
  ["52156071497798034055409940782395501364357", 1],
  ["50117127692312893981391715478615446797663", 1],
  ["24665048762541973552613860190140203906293", 1],
  ["43891111165377552467052180838054904286263", 1],
);
my @bls75s = (
  ["1000000000177", 1],
  ["1000000000000045819", 1],
  ["57850216533360484368293", 1],
  ["19568952034128395861091890269105913923337787205640409156470109155604436042237347889151", 1],
  ["2389755648366934394192070365850201237857", 0],
  ["32344792936896827502551761860817", 0],
);
my @np1s = (
  [391, 0],
  ["63699643930293116661668059033734770664712983894089510286262271", 1],
  ["17113454194771827263776721", 1],
);
my @akss = (  # Cover a number of the various tests before the big loop
  [1, 0],
  [15, 0],
  [44165497, 0],
  [136804519, 0],
  [31*31*31, 0],
  [51019, 0],
  [3, 1],
  [40841, 1],
  [74903, 1],   # Small input that actually runs the full test
);

my @composites = (35, 247, 377, 391, 527, 567, 2627, 5543, 13919, 14299, 23939, 47627, 86519, 92819);


plan tests => 0 + 6
                + 38
                + 43
                + 34
                + 2
                + 7   # _with_cert
                + 4   # Trial, Miller, N-1
                + scalar(@np1s)   # BLS75 N+1
                + scalar(@bls75s) # BLS75 hybrid
                + scalar(@ecpps)  # ecpp
                + scalar(@llrs)   # llr
                + scalar(@prs)    # proth
                + scalar(@akss)   # AKS
                + scalar(@composites)  # various composites
                + 2   # _validate_ecpp_curve
                + 0;

is(is_provable_prime(2) , 2,  '2 is prime');
is(is_provable_prime(1) , 0,  '1 is not prime');
is(is_provable_prime(0) , 0,  '0 is not prime');
is(is_provable_prime(-1), 0, '-1 is not prime');
is(is_provable_prime(-2), 0, '-2 is not prime');
is(is_provable_prime(20), 0, '20 is not prime');

foreach my $n (
  qw/2152302898747 3474749660383 341550071728321 341550071728321
     3825123056546413051 561 1105 1729 2465 2821 6601 8911 10585 15841 29341
     41041 46657 52633 4681 5461 6601 7957 8321 52633 88357 44287 47197 55969
     63139 74593 79003 82513 87913 88573 97567 44801 53971 79381
    /) {
  is(is_provable_prime($n), 0, "$n is not prime");
}

foreach my $n (
  qw/2 3 7 23 89 113 523 887 1129 1327 9551 15683 19609 31397 155921
     5 11 29 97 127 541 907 1151 1361 9587 15727 19661 31469 156007 360749
     370373 492227 1349651 1357333 2010881 4652507 17051887 20831533 47326913
     122164969 189695893 191913031 10726905041
    /) {
  is(is_provable_prime($n), 2, "$n is prime");
}

# Generated using:
#  perl -Iblib/lib -Iblib/arch -Mbignum -MMath::Prime::Util::GMP=:all -E '$|=1; use Time::HiRes qw(gettimeofday tv_interval); foreach my $e (62 .. 1000) { my $n = Math::BigInt->new(2)->bpow($e); $n = next_prime($n); my $s = [gettimeofday]; my $isp = is_prime($n); my $sec = tv_interval($s); printf("%d  %s  %d  %6.2fms\n", $e, "$n", $isp, $sec*1000.0); }' | grep ' 2 '
# after having is_prime do a quick test on all inputs.
# Remove all the inputs that take too long.
# Then run through GMP-ECPP to make sure.

foreach my $n (
  qw/9223372036854775837
     18446744073709551629
     73786976294838206473
     147573952589676412931
     295147905179352825889
     590295810358705651741
     1180591620717411303449
     2361183241434822606859
     18889465931478580854821
     37778931862957161709601
     75557863725914323419151
     302231454903657293676551
     604462909807314587353111
     38685626227668133590597803
     1237940039285380274899124357
     9903520314283042199192993897
     316912650057057350374175801351
     2535301200456458802993406410833
     162259276829213363391578010288167
     1298074214633706907132624082305051
     10384593717069655257060992658440473
     1329227995784915872903807060280345027
     680564733841876926926749214863536422929
     43556142965880123323311949751266331066401
     87112285931760246646623899502532662132821
     713623846352979940529142984724747568191373381
     2854495385411919762116571938898990272765493293
     196159429230833773869868419475239575503198607639501078831
     3138550867693340381917894711603833208051177722232017256453
     12554203470773361527671578846415332832204710888928069025857
     102844034832575377634685573909834406561420991602098741459288097
     210624583337114373395836055367340864637790190801098222508621955201
     14821387422376473014217086081112052205218558037201992197050570753012880593911817
     3351951982485649274893506249551461531869841455148098344430890360930441007518386744200468574541725856922507964546621512713438470702986642486608412251521039
    /) {
  is(is_provable_prime($n), 2, "$n is prime");
}

is(is_provable_prime('340282366920938463463374607431768211507'), 2, "is_prime(2**128+51) = 2");
is(is_provable_prime('340282366920938463463374607431768211621'), 2, "is_provable_prime(2**128+165) == 2");

is_deeply( [is_provable_prime_with_cert(0)],
           [0, ''],
           "is_provable_prime_with_cert(0)");
{ my @scert = is_provable_prime_with_cert(2);
  ok($scert[0] == 2 && $scert[1] =~ /\bType Small\nN 2\b/,
     "is_provable_prime_with_cert(2)"); }
{ my @scert = is_provable_prime_with_cert(96953);
  ok($scert[0] == 2 && $scert[1] =~ /\bType Small\nN 96953\b/,
     "is_provable_prime_with_cert(96953)"); }

my $proof;
$proof = <<EOPROOF;
[MPU - Primality Certificate]
Version 1.0

Proof for:
N 848301847013166693538593241183

Type BLS5
N  848301847013166693538593241183
Q[1]  77868535196553293845507
A[0]  3
----

Type BLS5
N  77868535196553293845507
Q[1]  1854012742775078424893
----

Type BLS5
N  1854012742775078424893
Q[1]  2373421880872807
----
EOPROOF
$proof =~ s/\n$//;
is_deeply( [is_provable_prime_with_cert("848301847013166693538593241183")],
           [2, $proof],
           "is_provable_prime_with_cert(848301847013166693538593241183)");

$proof = <<EOPROOF;
[MPU - Primality Certificate]
Version 1.0

Proof for:
N 316912650057057350374175801351

Type BLS5
N  316912650057057350374175801351
Q[1]  1451
Q[2]  487
Q[3]  443
Q[4]  5
A[0]  7
----
EOPROOF
$proof =~ s/\n$//;
is_deeply( [is_provable_prime_with_cert("316912650057057350374175801351")],
           [2, $proof],
           "is_provable_prime_with_cert(316912650057057350374175801351)");

$proof = <<EOPROOF;
[MPU - Primality Certificate]
Version 1.0

Proof for:
N 3138550867693340381917894711603833208051177722232017256453

Type BLS5
N  3138550867693340381917894711603833208051177722232017256453
Q[1]  119827
Q[2]  87211
Q[3]  5419
Q[4]  379
Q[5]  43
Q[6]  19
----
EOPROOF
$proof =~ s/\n$//;
my($isp, $cert) = is_provable_prime_with_cert("3138550867693340381917894711603833208051177722232017256453");
is($isp, 2, "is_provable_prime_with_cert(3138550867693340381917894711603833208051177722232017256453) is prime");
if ($cert =~ /\bType BLS5\b/) {
  is($cert, $proof, "is_provable_prime_with_cert(3138550867693340381917894711603833208051177722232017256453)");
} else {
  like($cert, qr/\nType BLS15\nN  3138550867693340381917894711603833208051177722232017256453\nQ  120713494911282322381457488907839738771199143162769894479\nLP 1\nLQ 6\n\n/, "is_provable_prime_with_cert(3138550867693340381917894711603833208051177722232017256453)");
}


#####################
# Individual routines

# Trial
ok( is_trial_prime(4539892831), "is_trial_prime(4539892831)" );

# Unconditional and conditional Miller test
ok( is_miller_prime("4835703278458516698824747"), "is_miller_prime(4835703278458516698824747)" );
ok( is_miller_prime("4835703278458516698824747",1), "is_miller_prime(4835703278458516698824747,1)" );

# BLS75 n-1
ok( is_nminus1_prime("340282366920938463463374607431768211507"), "is_nminus1_prime(340282366920938463463374607431768211507)" );

###### BLS75 (N+1)
for my $d (@np1s) {
  my($n,$exp) = @$d;
  is(is_nplus1_prime($n), $exp, "is_nplus1_prime($n) = $exp");
}
###### BLS75 (hybrid methods)
for my $d (@bls75s) {
  my($n,$exp) = @$d;
  is(is_bls75_prime($n), $exp, "is_bls75_prime($n) = $exp");
}
###### ECPP
for my $d (@ecpps) {
  my($n,$exp) = @$d;
  is(is_ecpp_prime($n), $exp, "is_ecpp_prime($n) = $exp");
}
###### llr
for my $d (@llrs) {
  my($n,$exp) = @$d;
  is(is_llr_prime($n), $exp, "is_llr_prime($n) = $exp");
}
###### proth
for my $d (@prs) {
  my($n,$exp) = @$d;
  is(is_proth_prime($n), $exp, "is_proth_prime($n) = $exp");
}
###### AKS
for my $d (@akss) {
  my($n,$exp) = @$d;
  is(is_aks_prime($n), $exp, "is_aks_prime($n) = $exp");
}

###### Test small composites
for my $n (@composites) {
  is_deeply(
    [is_provable_prime($n), is_trial_prime($n),
     is_miller_prime($n), is_ecpp_prime($n),
     is_nminus1_prime($n), is_nplus1_prime($n), is_bls75_prime($n)],
    [0,0,0,0,0,0,0],
    "No method says $n is prime" );
}

###### _validate_ecpp_curve (used by validator)
{
  my($a, $b, $n, $px, $py, $m, $q) = (qw/11974582979013017040544030800350811348567083551553 3745895188306636099862491686263587139808164603839 31502364452480398675892757951692134470845712114091 1130825859 12488879584249846189936346549055924440813808564924 31502364452480398675892768464573175914770636828425 166158920647084013166402856268513699557/);
  ok( Math::Prime::Util::GMP::_validate_ecpp_curve($a,$b,$n,$px,$py,$m,$q), "_validate_ecpp_curve ok" );
  $px = "1130825861";
  ok( !Math::Prime::Util::GMP::_validate_ecpp_curve($a,$b,$n,$px,$py,$m,$q), "_validate_ecpp_curve with different values fails" );
}
