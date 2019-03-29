#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/powint mulint addint subint divint modint divrem tdivrem absint negint/;
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

my @quotients = (  # trunc, floor, euclidian
  ["+ +",  "39458349850349850394853049583049",  "85889",  "459410982202026457344398579",  "459410982202026457344398579",  "459410982202026457344398579"],
  ["+ -",  "39458349850349850394853049583049", "-85889", "-459410982202026457344398579", "-459410982202026457344398580", "-459410982202026457344398579"],
  ["- +", "-39458349850349850394853049583049",  "85889", "-459410982202026457344398579", "-459410982202026457344398580", "-459410982202026457344398580"],
  ["- -", "-39458349850349850394853049583049", "-85889",  "459410982202026457344398579",  "459410982202026457344398579",  "459410982202026457344398580"],
);

plan tests => 0
            + 7*4 + scalar(@powints) + 2     # powint
            + 1 + scalar(@mulints)           # mulint
            + 1 + scalar(@addints)           # addint
            + 1 + scalar(@subints)           # subint
            + 2 + 2                          # divint
            + 2 + 2                          # modint
            + 2                              # divrem
            + 2                              # tdivrem
            + 4 * scalar(@quotients)         # signed bigint division
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
ok(!eval { divint(0,0); }, "divint(1,0)");
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
ok(!eval { modint(0,0); }, "modint(1,0)");
ok(!eval { modint(1,0); }, "modint(1,0)");

is_deeply( [map { modint(1024,$_) } 1..1025], \@rpos1024, "modint(1024,x) for 1 .. 1025" );
is_deeply( [map { modint(-1024,$_) } 1..1025], \@rneg1024, "modint(-1024,x) for 1 .. 1025" );

###### divrem
ok(!eval { divrem(0,0); }, "divrem(1,0)");
ok(!eval { divrem(1,0); }, "divrem(1,0)");

###### tdivrem
ok(!eval { tdivrem(0,0); }, "tdivrem(1,0)");
ok(!eval { tdivrem(1,0); }, "tdivrem(1,0)");


###### large values through divint, modint, divrem, tdivrem
for my $s (@quotients) {
  my($signs, $n, $m, $qt, $qf, $qe) = @$s;
  my($bn,$bm) = map { Math::BigInt->new($_) } ($n,$m);
  my($rt, $rf, $re) = map { $bn - $bm * $_ } ($qt, $qf, $qe);
  is( divint($n, $m), $qf, "large divint  $signs" );
  is( modint($n, $m), $rf, "large modint  $signs" );
  is_deeply( [divrem($n, $m)], [$qe, $re], "large divrem  $signs" );
  is_deeply( [tdivrem($n, $m)], [$qt, $rt], "large tdivrem $signs" );
}

###### absint
is_deeply([map { absint($_) } -9..9], [map { abs($_) } -9..9], "absint(-9..9)");
###### negint
is_deeply([map { negint($_) } -9..9], [map { -$_ } -9..9], "negint(-9..9)");
