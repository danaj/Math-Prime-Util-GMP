#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Math::BigInt try=>"GMP,Pari";
use Math::Prime::Util::GMP qw/urandomb urandomr   seed_csprng random_bytes/;

my $use64 = (~0 > 4294967295);
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $maxbits = $use64 ? 64 : 32;

my @nbit_check = (0..5,8,20,31,32,33,40);
my %nbit_range = map { $_ => 1 } grep { $_ <= 8 } @nbit_check;
my @large_nbit = (64,128,255,256,257,512,1024,2048,4096,8192,73100);

my $nsamples = $extra ? 30000 : 1000;
my $rsamples = $extra ? 10000 :  200;

plan tests => 0
            + scalar(@nbit_check)
            + scalar(@large_nbit)
            + scalar(keys %nbit_range)
            + 4
            + 5
            + 3;

########

check_nbit_range($_) for @nbit_check;

########

for my $b (@large_nbit) {
  my $n = Math::BigInt->new(urandomb($b));
  my $cmp = Math::BigInt->new(1)->blsft($b);
  cmp_ok($n, '<', $cmp, "Random $b-bit in range");
}

########

check_range(100,110);
check_range(128,255);
check_range(2**24, 2**25-1);
check_range(Math::BigInt->new(10)**24, Math::BigInt->new(10)**25-1);

########

ok(!eval { urandomr(-10,100); }, "urandomr(-10,x)");
ok(!eval { urandomr(100,-10); }, "urandomr(x,-10)");
ok(!eval { urandomr(-1,-1); }, "urandomr(-1,-1)");
is(urandomr(123456,123456), 123456, "urandomr(x,x)=x");
is(urandomr(123457,123456), undef, "urandomr(x,y)=undef if x > y");

########

seed_csprng(55,"BLAKEGrostlJHKeccakSkein--RijndaelSerpentTwofishRC6MARS");
is(unpack("h*",random_bytes( 4)),"538e1f65","random_bytes(4)");
is(unpack("h*",random_bytes(11)),"cbbac4ba12e6aa77bcfe6f","random_bytes(11)");
is(unpack("h*",random_bytes( 0)),"","random_bytes(0)");


########

sub check_nbit_range {
  my $b = shift;
  my $over = ($b < $maxbits) ? (1 << $b) : (Math::BigInt->new(1) << $b);
  if (!$nbit_range{$b}) {
    my @s = map { urandomb($b) } 1 .. $nsamples;
    is(scalar(grep { $_ >= $over } @s), 0, "Random $b-bit values are in range");
  } else {
    # For $b=8 and a uniform random generator, the probability of a given
    # 8-bit value not being selected is 1-1/256 = 0.99609375.  After 1000
    # tries, this is 0.01996 or about 2%.  After 1000+2000 tries it's ~8e-6.
    # The proper calculation is more tedious but for $b=8 and $nsamples=1000
    # we will find it likely we need a second pass but no more.
    my %t;
    for my $try (1..10) {
      $t{ urandomb($b) }++ for 1 .. $try * $nsamples;
      last if scalar(keys %t) >= $over;
    }
    my @keys = keys %t;
    is(scalar(grep { $_ >= $over } @keys), 0, "Random $b-bit values are in range");
    is(scalar(@keys), $over, "Random $b-bit produces all values in range");
  }
}

sub check_range {
  my($lo,$hi) = @_;
  my @s = map { urandomr($lo,$hi) } 1 .. $rsamples;
  @s = map { ref($hi)->new("$_") } @s if ref($hi);
  is( scalar(grep { $_ < $lo || $_ > $hi } @s), 0, "All random values from $lo to $hi in range" );
}
