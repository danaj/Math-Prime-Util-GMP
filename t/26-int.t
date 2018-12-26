#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::Prime::Util::GMP qw/powint mulint addint divint remint divrem tdivrem/;
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
            + 2 + 2 + scalar(@quotients)     # divint
            + 2 + 2 + scalar(@quotients)     # remint
            + 2 + scalar(@quotients)         # divrem
            + 2 + scalar(@quotients)         # tdivrem
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

###### divint
ok(!eval { divint(0,0); }, "divint(1,0)");
ok(!eval { divint(1,0); }, "divint(1,0)");

is_deeply( [map { int(1024/$_) } 1..1025], [map { divint(1024,$_) } 1..1025], "divint(1024,x) for 1 .. 1025" );
# Note this will be DIFFERENT than int(-1024/$_)
{ my $num = Math::BigInt->new(-1024);
  is_deeply( [map { $num/$_ } 1..1025], [map { divint(-1024,$_) } 1..1025], "divint(-1024,x) for 1 .. 1025" );
}

for my $s (@quotients) {
  my($signs, $n, $m, $qt, $qf, $qe) = @$s;
  is( divint($n, $m), $qf, "large divint: $signs" );
}

###### remint
ok(!eval { remint(0,0); }, "remint(1,0)");
ok(!eval { remint(1,0); }, "remint(1,0)");

is_deeply( [map { int(1024%$_) } 1..1025], [map { remint(1024,$_) } 1..1025], "remint(1024,x) for 1 .. 1025" );
# Note this will be DIFFERENT than int(-1024 % $_)
{ my $num = Math::BigInt->new(-1024);
  is_deeply( [map { $num%$_ } 1..1025], [map { remint(-1024,$_) } 1..1025], "remint(-1024,x) for 1 .. 1025" );
}

for my $s (@quotients) {
  my($signs, $n, $m, $qt, $qf, $qe) = @$s;
  my $rf = $n - $m * $qf;
  is( remint($n, $m), $rf, "large remint: $signs" );
}

###### divrem
ok(!eval { divrem(0,0); }, "divrem(1,0)");
ok(!eval { divrem(1,0); }, "divrem(1,0)");

for my $s (@quotients) {
  my($signs, $n, $m, $qt, $qf, $qe) = @$s;
  my $re = $n - $m * $qe;
  is_deeply( [divrem($n, $m)], [$qe, $re], "large divrem: $signs" );
}

###### tdivrem
ok(!eval { tdivrem(0,0); }, "tdivrem(1,0)");
ok(!eval { tdivrem(1,0); }, "tdivrem(1,0)");

for my $s (@quotients) {
  my($signs, $n, $m, $qt, $qf, $qe) = @$s;
  my $rt = $n - $m * $qt;
  is_deeply( [tdivrem($n, $m)], [$qt, $rt], "large tdivrem: $signs" );
}
