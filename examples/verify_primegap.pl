#!/usr/bin/perl

# Verify a prime gap given the start (as number or expression)

# If start is less than 200 digits, then endpoints will be proven.
#
# Otherwise, all values will use BPSW + 1 M-R:
#    a) divisibility by small factors
#    b) SPSP-2
#    c) extra strong Lucas test
#    d) 1 Miller-Rabin test with a random base 3 .. N-2
#
# No composite has been found that passes (b) + (c), so the combination
# plus the extra M-R should give quite a bit of certainty that the
# endpoints really are prime.
#
# Runs in parallel for number larger than 200 digits.

# Examples:
#
#   verify_primegap.pl "9169*439#/55230 - 6926"
#   verify_primegap.pl "1931*1933#/7230 - 30244"  2
#   verify_primegap.pl "5557*4973#/2310 - 83542" 12

use warnings;
use strict;
use threads;
use threads::shared;
use feature 'say';
use Math::Prime::Util ':all';
use Math::Prime::Util::GMP;
use Math::BigInt try=>"GMP";
use Math::BigFloat;
$|=1;
my @mpu_funcs = (qw/next_prime prev_prime prime_count nth_prime random_prime
                    random_ndigit_prime random_nbit_prime random_strong_prime
                    random_maurer_prime primorial pn_primorial moebius mertens
                    euler_phi jordan_totient exp_mangoldt divisor_sum
                    consecutive_integer_lcm/);
my %mpu_func_map;

my $N = shift || dieusage();
my $nthreads = shift || 8;

$N =~ s/\s*$//;  $N =~ s/^\s*//;
$N = eval_expr($N) unless $N =~ /^\d+$/;
$N = Math::BigInt->new("$N") unless ref($N) eq 'Math::BigInt';
$nthreads = 1 if length($N) <= 200;
say_primality("start (" . length($N) . " digits)", $N);

if ($nthreads <= 1) {
  my $end = next_prime($N);
  my $gap = $end-$N;
  my $merit = merit($N, $gap);
  say_primality("Endpoint (merit $merit) n+$gap", $end);
  exit(0);
}

my $searchto :shared = 1;
my $found    :shared = 0;

my @threads;
push @threads, threads->create('search', $_) for 1 .. $nthreads;
$_->join() for (@threads);
{
  my $gap = $found;
  my $end = $N+$found;
  my $merit = merit($N, $gap);
  say_primality("Endpoint (merit $merit) n+$gap", $end);
  exit(0);
}

sub search {
  my $tnum = shift;
  while (!$found) {
    my ($j, $n);
    {
      lock($searchto);
      $j = $searchto++;
    }
    #print "  +$j($tnum)  ";
    if ( Math::Prime::Util::GMP::is_prime($N + $j) ) {
      lock($found);
      if ($found == 0 || $found > $j) {
        $found = $j;
      }
      last;
    }
  }
}


sub eval_expr {
  my $expr = shift;
  die "$expr cannot be evaluated" if $expr =~ /:/;  # Use : for escape
  if (scalar(keys %mpu_func_map) == 0) {
    my $n = 10;
    foreach my $func (@mpu_funcs) {
      $mpu_func_map{$func} = sprintf("%03d", $n++);
    }
  }
  $expr =~ s/\^/**/g;
  $expr =~ s/\b(\d+)#/primorial($1)/g;
  $expr =~ s/\blog\(/:001(/g;
  foreach my $func (@mpu_funcs) {
    $expr =~ s/\b$func\(/:$mpu_func_map{$func}(/g;
  }
  die "$expr cannot be evaluated" if $expr =~ tr|-0123456789+*/() :||c;
  $expr =~ s/:001/log/g;
  foreach my $func (@mpu_funcs) {
    $expr =~ s/:$mpu_func_map{$func}\(/ Math::BigInt->bone*Math::Prime::Util::$func(/g;
  }
  $expr =~ s/(\d+)/ Math::BigInt->new("$1") /g;
  my $res = eval $expr; ## no critic
  die "Cannot eval: $expr\n" if !defined $res;
  $res = int($res->bstr) if ref($res) eq 'Math::BigInt' && $res <= ~0;
  $res;
}

sub say_primality {
  my($text, $n) = @_;
  print "$text is ";
  my $is_prime = (length($n) <= 200) ? is_provable_prime($n) : is_prime($n);
  print "", ("composite",
             "probably prime (BPSW + 1 M-R)",
             "proven prime")[$is_prime];
  print "\n";
}

sub merit {
  my($n, $gap) = @_;
  my $fgap = Math::BigFloat->new("$gap");
  my $fn   = Math::BigFloat->new("$n");
  return sprintf("%7.4f", $fgap / log($fn));
}

sub dieusage {
  die "Usage: $0  \"expression\"  [number-of-threads]\n";
}
