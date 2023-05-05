#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/powint mulint addint subint divint modint cdivint divrem tdivrem fdivrem cdivrem absint negint lshiftint rshiftint rashiftint/;
use Math::BigInt;  # Don't use GMP so we don't have to work around bug

my $use64 = (~0 > 4294967296 && 18446744073709550592 != ~0);
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};

my @powints = (
 [5, 6, 15625],
 [2, 16, 65536],
);
my @mulints = (
  ["13282407956253574712","14991082624209354397","199117675120653046511338473800925208664"],
);
my @addints = (
  ["1178630961471601951655862","827639478068904540012","1179458600949670856195874"],
  ["-2555488174170453670799","1726145541361106236340","-829342632809347434459"],
);
my @subints = (
  ["68719214592","281474976448512","-281406257233920"],
  ["38631281077","12191281349924010278","-12191281311292729201"],
  ["-38631281077","12191281349924010278","-12191281388555291355"],
  ["-38631281077","-12191281349924010278","12191281311292729201"],
  ["9686117847286759","419039659798583","9267078187488176"],
  ["14888606332669627740937300680965976203","14888605897080617527808122501731945103","435589010213129178179234031100"],
);

my @quotients = (  # trunc, floor, ceil, euclidian
  ["S + +", "9949242744253247", "64", "155456917878956", "155456917878956", "155456917878957", "155456917878956"],
  ["S - +", "-9949242744253247", "64", "-155456917878956", "-155456917878957", "-155456917878956", "-155456917878957"],
  ["L + +",  "39458349850349850394853049583049",  "85889",  "459410982202026457344398579",  "459410982202026457344398579",  "459410982202026457344398580",  "459410982202026457344398579"],
  ["L + -",  "39458349850349850394853049583049", "-85889", "-459410982202026457344398579", "-459410982202026457344398580", "-459410982202026457344398579", "-459410982202026457344398579"],
  ["L - +", "-39458349850349850394853049583049",  "85889", "-459410982202026457344398579", "-459410982202026457344398580", "-459410982202026457344398579", "-459410982202026457344398580"],
  ["L - -", "-39458349850349850394853049583049", "-85889",  "459410982202026457344398579",  "459410982202026457344398579",  "459410982202026457344398580",  "459410982202026457344398580"],
);

my @negshifts = (
  # n, k,  >>, >>arith
  [ 0, 1,  0, 0],
  [-1, 1,  0, -1],
  [-5, 1,  -2, -3],
  [-8, 2,  -2, -2],
  ["-307385513", 6, -4802898, -4802899],
  ["-637526413", 6, -9961350, -9961351],
  ["-2045651239", 6, -31963300, -31963301],
  ["-3675663743", 6, -57432245, -57432246],
  ["-2332267979728172537", 6, "-36441687183252695", "-36441687183252696"],
  ["-8408654401686460807", 6, "-131385225026350950", "-131385225026350951"],
  ["-17640827963513397449", 6, "-275637936929896835", "-275637936929896836"],
  ["-32659506018295865747", 6, "-510304781535872902", "-510304781535872903"],
  ["-79231600218559026832557301750107210001", 6, "-1237993753414984794258707839845425156", "-1237993753414984794258707839845425157"],
  ["-131954888069700539887213633881194728277", 6, "-2061795126089070935737713029393667629", "-2061795126089070935737713029393667630"],
  ["-254262665582332530470619504253273698569", 6, "-3972854149723945788603429753957401540", "-3972854149723945788603429753957401541"],
  ["-416649423645764932216789232242651032187", 6, "-6510147244465077065887331753791422377", "-6510147244465077065887331753791422378"],
);

plan tests => 0
            + 7*4 + scalar(@powints) + 2     # powint
            + 1 + scalar(@mulints)           # mulint
            + 1 + scalar(@addints)           # addint
            + 1 + scalar(@subints)           # subint
            + 2 + 2                          # divint
            + 2 + 2                          # modint
            + 2 + 1                          # cdivint
            + 2                              # divrem
            + 2                              # tdivrem
            + 2                              # fdivrem
            + 2                              # cdivrem
            + 7 * scalar(@quotients)         # signed bigint division
            + 4 + 3*scalar(@negshifts)       # shiftint
            + 1                              # absint
            + 1                              # negint
            + 0;

###### powint
for my $a (-3 .. 3) {
  is(powint($a, 0), 1, "powint($a,0) = 1");
  is(powint($a, 1), $a, "powint($a,1) = $a");
  is(powint($a, 2), $a*$a, "powint($a,2) = " . $a*$a);
  is(powint($a, 3), $a*$a*$a, "powint($a,3) = " . $a*$a*$a);
}
foreach my $r (@powints) {
  my($a, $b, $exp) = @$r;
  is( powint($a,$b), $exp, "powint($a,$b) = ".((defined $exp)?$exp:"<undef>") );
}
is(powint(powint(2,32),3),"79228162514264337593543950336","(2^32)^3");
is(powint(3,powint(2,7)),"11790184577738583171520872861412518665678211592275841109096961","3^(2^7)");

###### mulint
{ my(@got,@exp);
  for my $a (-3 .. 3) {
    for my $b (-3 .. 3) {
      push @got, mulint($a,$b);
      push @exp, $a*$b;
    }
  }
  is_deeply( \@got, \@exp, "mulint( -3 .. 3, -3 .. 3)" );
}
foreach my $r (@mulints) {
  my($a, $b, $exp) = @$r;
  is( mulint($a,$b), $exp, "mulint($a,$b) = ".((defined $exp)?$exp:"<undef>") );
}

###### addint
{ my(@got,@exp);
  for my $a (-3 .. 3) {
    for my $b (-3 .. 3) {
      push @got, addint($a,$b);
      push @exp, $a+$b;
    }
  }
  is_deeply( \@got, \@exp, "addint( -3 .. 3, -3 .. 3)" );
}
foreach my $r (@addints) {
  my($a, $b, $exp) = @$r;
  is( addint($a,$b), $exp, "addint($a,$b) = ".((defined $exp)?$exp:"<undef>") );
}

###### subint
{ my(@got,@exp);
  for my $a (-3 .. 3) {
    for my $b (-3 .. 3) {
      push @got, subint($a,$b);
      push @exp, $a-$b;
    }
  }
  is_deeply( \@got, \@exp, "subint( -3 .. 3, -3 .. 3)" );
}
foreach my $r (@subints) {
  my($a, $b, $exp) = @$r;
  is( subint($a,$b), $exp, "subint($a,$b) = ".((defined $exp)?$exp:"<undef>") );
}

###### divint
ok(!eval { divint(0,0); }, "divint(0,0)");
ok(!eval { divint(1,0); }, "divint(1,0)");

# For negative inputs, the div and mod operations might be different than Perl's builtins.
# It matches Math::BigInt bdiv / bmod (post 1.997 Sep 2015).

my @qpos1024 = map { int(1024/$_) } 1 .. 1025;
my @qneg1024 = map { my $d=-1024/$_; my $i = int($d);  ($d==$i) ? $i : $i-1; } 1 .. 1025;

my @rpos1024 = map {  1024 - $_ * $qpos1024[$_-1] } 1 .. 1025;
my @rneg1024 = map { -1024 - $_ * $qneg1024[$_-1] } 1 .. 1025;

is_deeply( [map { divint(1024,$_) } 1..1025], \@qpos1024, "divint(1024,x) for 1 .. 1025" );
is_deeply( [map { divint(-1024,$_) } 1..1025], \@qneg1024, "divint(-1024,x) for 1 .. 1025" );

###### modint
ok(!eval { modint(0,0); }, "modint(0,0)");
ok(!eval { modint(1,0); }, "modint(1,0)");

is_deeply( [map { modint(1024,$_) } 1..1025], \@rpos1024, "modint(1024,x) for 1 .. 1025" );
is_deeply( [map { modint(-1024,$_) } 1..1025], \@rneg1024, "modint(-1024,x) for 1 .. 1025" );

###### cdivint
ok(!eval { cdivint(0,0); }, "cdivint(0,0)");
ok(!eval { cdivint(1,0); }, "cdivint(1,0)");

is_deeply([cdivint(7,3),cdivint(7,-3),cdivint(-7,3),cdivint(-7,-3)],
          [3,-2,-2,3],
          "cdivint with all signs of 7,3");

###### divrem
ok(!eval { divrem(0,0); }, "divrem(0,0)");
ok(!eval { divrem(1,0); }, "divrem(1,0)");

###### tdivrem
ok(!eval { tdivrem(0,0); }, "tdivrem(0,0)");
ok(!eval { tdivrem(1,0); }, "tdivrem(1,0)");

###### fdivrem
ok(!eval { fdivrem(0,0); }, "fdivrem(0,0)");
ok(!eval { fdivrem(1,0); }, "fdivrem(1,0)");

###### cdivrem
ok(!eval { cdivrem(0,0); }, "cdivrem(0,0)");
ok(!eval { cdivrem(1,0); }, "cdivrem(1,0)");

###### large values through divint, cdivint, modint,
######                      divrem, tdivrem, fdivrem, cdivrem
for my $s (@quotients) {
  my($signs, $n, $m, $qt, $qf, $qc, $qe) = @$s;
  my($bn,$bm) = map { Math::BigInt->new($_) } ($n,$m);
  my($rt, $rf, $rc, $re) = map { $bn - $bm * $_ } ($qt, $qf, $qc, $qe);
  is( divint($n, $m), $qf, "large divint  $signs" );
  is( modint($n, $m), $rf, "large modint  $signs" );
  is( cdivint($n, $m), $qc, "large divint  $signs" );
  is_deeply( [divrem($n, $m)], [$qe, $re], "large divrem  $signs" );
  is_deeply( [tdivrem($n, $m)], [$qt, $rt], "large tdivrem $signs" );
  is_deeply( [fdivrem($n, $m)], [$qf, $rf], "large fdivrem $signs" );
  is_deeply( [cdivrem($n, $m)], [$qc, $rc], "large cdivrem $signs" );
}

###### lshiftint
is_deeply([map { lshiftint($_) } 0..50], [map { $_ << 1 } 0..50], "lshiftint(0..50)");
is_deeply([map { rshiftint($_) } 0..50], [map { $_ >> 1 } 0..50], "rshiftint(0..50)");
is_deeply([map { rashiftint($_) } 0..50], [map { $_ >> 1 } 0..50], "rashiftint(0..50)");
is_deeply([map { lshiftint($_,5) } -65 .. 65], [map { $_ * 32 } -65 .. 65], "lshiftint(-65 .. 65, 5)");

for my $d (@negshifts) {
  my($n, $k, $rs, $ras) = @$d;
  my $ls = mulint($n, powint(2,$k));
  is( lshiftint($n,$k), $ls, "lshiftint($n,$k) = $ls" );
  is( rshiftint($n,$k), $rs, "rshiftint($n,$k) = $rs" );
  is( rashiftint($n,$k), $ras, "rashiftint($n,$k) = $ras" );
}

###### absint
is_deeply([map { absint($_) } -9..9], [map { abs($_) } -9..9], "absint(-9..9)");
###### negint
is_deeply([map { negint($_) } -9..9], [map { -$_ } -9..9], "negint(-9..9)");
