#!/usr/bin/env perl
use warnings;
use strict;
use Math::BigInt try=>"GMP,Pari";

# Convert the output of GMP-ECPP to a MPU certificate.
# Written by Dana Jacobsen, 2013.

# This does no error checking and is very simple.  It's up to verify-cert.pl
# to make sure the resulting MPU certificate is parsable and valid.  It might
# be nice to have this do a little more error checking just so any format
# errors could get caught right away.
# As an example, MPU certificates really don't care what order things arrive,
# as it just follows the tree of "N is prime if [Q1,Q2,Q3,...] is prime" down.
# PRIMO and GMP-ECPP both restrict the tree to a single child and in-order.
# This converter ought to verify that each N is the previous step's Q.

print "[MPU - Primality Certificate]\n";
print "Version 1.0\n";
print "\n";
print "# Converted from GMP-ECPP by convert-gmpecpp version 1.0\n";
print "\n";

while (<>) {
  if (/^N\[(\d+)\]\s*=\s*(\d+)/) {
    my($step, $n) = ($1, $2);
    if ($step == 0) {
      print "Proof for:\n";
      print "N $n\n";
      print "\n";
    }
    print "Type ECPP\n";
    print "N $n\n";
  }
  elsif (/^a\s*=\s*(\d+)/)       { print "A $1\n"; }
  elsif (/^b\s*=\s*(\d+)/)       { print "B $1\n"; }
  elsif (/^m\s*=\s*(\d+)/)       { print "M $1\n"; }
  elsif (/^q\s*=\s*(\d+)/)       { print "Q $1\n"; }
  elsif (/^P\s*=\s*\(\s*(\d+)\s*,\s*(\d+)\s*\)/) {
    print "X $1\n";
    print "Y $2\n";
    print "\n";
  }
}
