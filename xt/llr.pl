#!/usr/bin/env perl
use warnings;
use strict;
use Math::Prime::Util::GMP ":all";
use Math::GMPz;
$|=1;

llr_test(10,10+8*80-1,1,100);
llr_test(10,10+4*80-1,101,1000);
llr_test(14,14+2*80-1,1001,10000);
llr_test(17,17+1*80-1,10001,100000);


sub llr_test {
  my($b1, $b2, $k1, $k2) = @_;

  print "Testing LLR with {$k1..$k2} * 2^{$b1..$b2} - 1\n";
  for my $b ($b1..$b2) {
    print ".";
    my $b2 = Math::GMPz->new(2)**$b;
    for my $k ($k1..$k2) {
      my $n = $k * $b2 - 1;
      my $t1 = is_bpsw_prime($n) ? 2 : 0;
      my $t2 = is_llr_prime($n);
      die "$t1 $t2 ${k} * 2^$b -1" unless $t1==$t2;
    }
  }
  print "\n";
}
