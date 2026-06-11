#!/usr/bin/env perl
use strict;
use warnings;

use Benchmark qw/cmpthese/;
use Getopt::Long qw/GetOptions/;
#use FindBin;
#use lib "$FindBin::Bin/../blib/lib", "$FindBin::Bin/../blib/arch";

use Math::Prime::Util::GMP qw/is_power powint addint negint/;

my $seconds = 3;
my $count;
my $heavy = 0;
my $help = 0;

GetOptions(
  "seconds|s=f" => \$seconds,
  "count|n=i"   => \$count,
  "heavy!"      => \$heavy,
  "help|h"      => \$help,
) or usage();
usage() if $help;

sub usage {
  print <<"END_USAGE";
usage: perl -Mblib xt/bench-ispower.pl [options]

Options:
  --seconds N   Run each benchmark group for about N seconds (default: 3)
  --count N     Run each benchmark case N iterations instead of timed mode
  --heavy       Include very large cases
  --help        Show this message
END_USAGE
  exit 0;
}

sub bench_count {
  return defined $count ? $count : -$seconds;
}

sub digits {
  my ($n) = @_;
  return length("$n");
}

sub add_case {
  my ($cases, $name, $n, $expect) = @_;
  push @$cases, {
    name   => $name,
    n      => $n,
    digits => digits($n),
    expect => $expect,
    code   => sub { is_power($n) },
  };
}

sub add_case_k {
  my ($cases, $name, $n, $k, $expect) = @_;
  push @$cases, {
    name   => $name,
    n      => $n,
    digits => digits($n),
    expect => $expect,
    code   => sub { is_power($n, $k) },
  };
}

sub run_group {
  my ($title, $cases) = @_;
  my %bench;

  print "\n# $title\n";
  for my $case (@$cases) {
    my $got = $case->{code}->();
    die "$case->{name}: got $got, expected $case->{expect}\n"
      unless "$got" eq "$case->{expect}";
    printf "# %-28s %8d digits  => %s\n",
      $case->{name}, $case->{digits}, $case->{expect};
    $bench{$case->{name}} = $case->{code};
  }
  cmpthese(bench_count(), \%bench);
}

my $p101_e101       = powint(101, 101);
my $c202_e101       = powint(202, 101);
my $c202_e10001     = powint(202, 10001);
my $p101_e10001     = powint(101, 10001);
my $two_e10001      = powint(2, 10001);
my $near_202_e10001 = addint($c202_e10001, 1);
my $neg_202_e10001  = negint($c202_e10001);
my $square_large    = powint("29905047121918201644964877983907", 2);
my $e16_composite   = powint(903111, 16);

my @auto;
add_case(\@auto, "non-power-64bit", "18475335773296164196", 0);
add_case(\@auto, "large-square", $square_large, 2);
add_case(\@auto, "903111^16", $e16_composite, 16);
add_case(\@auto, "101^101", $p101_e101, 101);
add_case(\@auto, "202^101", $c202_e101, 101);
add_case(\@auto, "2^10001", $two_e10001, 10001);
add_case(\@auto, "101^10001", $p101_e10001, 10001);
add_case(\@auto, "202^10001", $c202_e10001, 10001);
add_case(\@auto, "202^10001+1", $near_202_e10001, 0);
add_case(\@auto, "-202^10001", $neg_202_e10001, 10001);

run_group("is_power(n), auto exponent", \@auto);

my @fixed;
add_case_k(\@fixed, "202^10001,k=2",     $c202_e10001, 2,     0);
add_case_k(\@fixed, "202^10001,k=71",    $c202_e10001, 71,    0);
add_case_k(\@fixed, "202^10001,k=73",    $c202_e10001, 73,    1);
add_case_k(\@fixed, "202^10001,k=10001", $c202_e10001, 10001, 1);
add_case_k(\@fixed, "202^10001,k=10002", $c202_e10001, 10002, 0);
add_case_k(\@fixed, "-202^10001,k=10001", $neg_202_e10001, 10001, 1);

run_group("is_power(n,k), fixed exponent", \@fixed);

if ($heavy) {
  my @heavy;
  my $p12347_e12347   = powint(12347, 12347);
  my $p123457_e123457 = powint(123457, 123457);
  my $near_heavy      = addint($p12347_e12347, 1);

  add_case(\@heavy, "12347^12347", $p12347_e12347, 12347);
  add_case(\@heavy, "12347^12347+1", $near_heavy, 0);
  add_case(\@heavy, "123457^123457", $p123457_e123457, 123457);

  run_group("is_power(n), heavy", \@heavy);
}
