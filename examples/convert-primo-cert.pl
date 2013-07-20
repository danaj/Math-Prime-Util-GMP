#!/usr/bin/env perl
use warnings;
use strict;
use Math::BigInt try=>"GMP,Pari";

# Convert a Primo certificate to a MPU certificate.
# Written by Dana Jacobsen, 2013.

# This uses type ECPP3 and ECPP4 to make smaller files, letting the
# MPU certificate prover handle the check and convert to ECPP type.

print "[MPU - Primality Certificate]\n";
print "Version 1.0\n";
print "\n";
print "# Converted from Primo by convert-primo version 1.0\n";
print "\n";

my $n;
while (<>) {
  if (/^N\$=(\S+)/) {
    $n = Math::BigInt->new('0x' . $1);
    print "Proof for:\n";
    print "N $n\n";
    print "\n";
  } elsif (/^Type=(\d)/) {
    my $type = $1;
    if ($type == 4) {
      my ($s, $r, $j, $t) = read_vars('4', qw/S R J T/);
      print "Type ECPP4\n";
      print "N $n\nS $s\nR $r\nJ $j\nT $t\n\n";
      $n = $r;
    } elsif ($type == 3) {
      my ($s, $r, $a, $b, $t) = read_vars('3', qw/S R A B T/);
      print "Type ECPP3\n";
      print "N $n\nS $s\nR $r\nA $a\nB $b\nT $t\n\n";
      $n = $r;
    } elsif ($type == 2) {
      my ($s, $r, $q) = read_vars('2', qw/S R Q/);
      print "Type BLS15\n";
      my $p = ($q->is_odd()) ? 2 : 1;
      print "N $n\nQ $r\nLP $p\nLQ $q\n\n";
      $n = $r;
    } elsif ($type == 1) {
      my ($s, $r, $b) = read_vars('1', qw/S R B/);
      print "Type Pocklington\n";
      print "N $n\nQ $r\nA $b\n\n";
      $n = $r;
    } elsif ($type == 0) {
      # End
    } else {
      die "Unkown type: $type\n";
    }
  }
}

sub read_vars {
  my $type = shift;
  my %vars = map { $_ => 1 } @_;
  my %return;
  while (scalar keys %vars) {
    my $line = <>;
    die "end of file during type $type\n" unless defined $line;
    # Skip comments and blank lines
    next if $line =~ /^\s*#/ or $line =~ /^\s*$/;
    chomp($line);
    die "Still missing values in type $type\n" if $line =~ /^Type/;
    if ($line =~ /^(\S+)\$\s*=\s*(\S+)/) {
      my ($var, $val) = ($1, $2);
      $var =~ tr/a-z/A-Z/;
      die "Type $type: repeated or inappropriate var: $line\n" unless defined $vars{$var};
      $return{$var} = $val;
      delete $vars{$var};
    } else {
      die "Unrecognized line: $line\n";
    }
  }
  # Now return them in the order given, turned into bigints.
  my @ret;
  foreach my $var (@_) {
    my $sign = 1;
    $sign = -1 if $return{$var} =~ s/^-//;
    push @ret, Math::BigInt->new('0x' . $return{$var}) * $sign;
  }
  return @ret;
}
