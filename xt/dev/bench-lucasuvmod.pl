#!/usr/bin/perl
use warnings;
use strict;
use v5.32;
use Math::Prime::Util qw/vecsum cmpint prime_set_config/;
use Math::Prime::Util::PP;
use Math::GMP;
use Time::HiRes qw(gettimeofday tv_interval);

# Arbitrary 512-bit prime
my $n512 = Math::GMP->new('12700353393137284892745741745829061113670647614984104321625902290988934519411018386240707141830889777521351373268205080514233082985759465255049533377038287');
my $n64 = 13872796489734606043;
my $p128 = Math::GMP->new('270586991526651512219299876930645792417');
my $q128 = Math::GMP->new('193706242648041580937439909991206419311');

#prime_set_config( bigint => "Math::GMP" );
*lucasuvmod = \&Math::Prime::Util::lucasuvmod;
#*lucasuvmod = \&Math::Prime::Util::GMP::lucasuvmod;
#*lucasuvmod = \&Math::Prime::Util::PP::lucasuvmod;
my ($s1,$e1) = (1<<12, (1<<13)-1);
my ($s2,$e2) = (1<<12, (1<<13)-1);

my $donative = 1;
my $dobigint = 1;

if ($donative) {
  my @R;
  my $t0 = [gettimeofday];
  for my $k ($s1 .. $e1) {
    push @R, lucasuvmod($k, 1, $k, $n64);
  }
  print "====== native (x,1)   ";
  say "ns/call: ", 1e9 * tv_interval($t0, [gettimeofday])/($e1-$s1+1);
  say "    !!!!! Failed ",vecsum(@R),"\n\n" if cmpint(vecsum(@R),'56220166018828615152921');
}

if ($donative) {
  my @R;
  my $t0 = [gettimeofday];
  for my $k ($s1 .. $e1) {
    push @R, lucasuvmod($k, -1, $k, $n64);
  }
  print "====== native (x,-1)  ";
  say "ns/call: ", 1e9 * tv_interval($t0, [gettimeofday])/($e1-$s1+1);
  say "    !!!!! Failed ",vecsum(@R),"\n\n" if cmpint(vecsum(@R),'56421815297373700477619');
}

if ($donative) {
  my @R;
  my $t0 = [gettimeofday];
  for my $k ($s1 .. $e1) {
    push @R, lucasuvmod($k, $e1-$k, $k, $n64);
  }
  print "====== native (x,y)   ";
  say "ns/call: ", 1e9 * tv_interval($t0, [gettimeofday])/($e1-$s1+1);
  say "    !!!!! Failed ",vecsum(@R),"\n\n" if cmpint(vecsum(@R),'56147979494347120335058');
}

if ($donative) {
  my @R;
  my $t0 = [gettimeofday];
  for my $k ($s1 .. $e1) {
    push @R, lucasuvmod($k, 1, $k, $n64-1);
  }
  print "====== neven (x,1)    ";
  say "ns/call: ", 1e9 * tv_interval($t0, [gettimeofday])/($e1-$s1+1);
  say "    !!!!! Failed ",vecsum(@R),"\n\n" if cmpint(vecsum(@R),'56442715926104594453086');
}

if ($donative) {
  my @R;
  my $t0 = [gettimeofday];
  for my $k ($s1 .. $e1) {
    push @R, lucasuvmod($k, -1, $k, $n64-1);
  }
  print "====== neven (x,-1)   ";
  say "ns/call: ", 1e9 * tv_interval($t0, [gettimeofday])/($e1-$s1+1);
  say "    !!!!! Failed ",vecsum(@R),"\n\n" if cmpint(vecsum(@R),'56787171408989331382208');
}

if ($donative) {
  my @R;
  my $t0 = [gettimeofday];
  for my $k ($s1 .. $e1) {
    push @R, lucasuvmod($k, $e1-$k, $k, $n64-1);
  }
  print "====== neven (x,y)    ";
  say "ns/call: ", 1e9 * tv_interval($t0, [gettimeofday])/($e1-$s1+1);
  say "    !!!!! Failed ",vecsum(@R),"\n\n" if cmpint(vecsum(@R),'55983778281973304712626');
}



if ($dobigint) {
  my @R;
  my $t0 = [gettimeofday];
  for my $k ($s2 .. $e2) {
    push @R, lucasuvmod($p128+$k, 1, $k, $n512);
  }
  print "====== bigint (x,1)   ";
  say "ns/call: ", 1e9 * tv_interval($t0, [gettimeofday])/($e2-$s2+1);
  say "    !!!!! Failed ",vecsum(@R),"\n\n" if cmpint(vecsum(@R),'51955411897903644162368223173279166218537966348158568194465494079323596381651343897362523073645943266210069212540811101144445670324100791137644806738785766024');
}

if ($dobigint) {
  my @R;
  my $t0 = [gettimeofday];
  for my $k ($s2 .. $e2) {
    push @R, lucasuvmod($p128+$k, -1, $k, $n512);
  }
  print "====== bigint (x,-1)  ";
  say "ns/call: ", 1e9 * tv_interval($t0, [gettimeofday])/($e2-$s2+1);
  say "    !!!!! Failed ",vecsum(@R),"\n\n" if cmpint(vecsum(@R),'52273544602498573975470709515421047559941670767294728202447967843063610032765457849448319782712107103134313065445733589688268308946471399284566976410893409342');
}

if ($dobigint) {
  my @R;
  my $t0 = [gettimeofday];
  for my $k ($s2 .. $e2) {
    push @R, lucasuvmod($p128+$k, $q128+$k, $k, $n512);
  }
  print "====== bigint (x,y)   ";
  say "ns/call: ", 1e9 * tv_interval($t0, [gettimeofday])/($e2-$s2+1);
  say "    !!!!! Failed ",vecsum(@R),"\n\n" if cmpint(vecsum(@R),'51678854322814467640897779087169852143398839526233725755011603994407883436087766527182076641509939909200229468050193338674621439729380438506212357429283410740');
}



# Using definition:
#
# mmpu 'use bigint lib=>"GMP"; ($P,$Q,$k,$n)=(1,-1,27,7777); ($P,$Q)=($P % $n,$Q % $n); @U=(0,1); @V=(2,$P); for (2..$k) { $U[$_]=($P*$U[$_-1]-$Q*$U[$_-2]) % $n; $V[$_]=($P*$V[$_-1]-$Q*$V[$_-2]) % $n; } say $U[$k], " ", $V[$k]'
