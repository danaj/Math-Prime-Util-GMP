#!/usr/bin/env perl
use strict;
use warnings;

use Test::More;
use Config;
use Math::Prime::Util::GMP qw/seed_csprng is_csprng_well_seeded
                              irand irand64 drand random_bytes/;

my $use64 = (~0 > 4294967295);
my $extra = defined $ENV{EXTENDED_TESTING} && $ENV{EXTENDED_TESTING};
my $maxbits = $use64 ? 64 : 32;

my $samples = $extra ? 100000 :  10000;

plan tests => 0
            + 2
            + ($use64 ? 2 : 0)
            + 4
            + 4;

########

# seed_csprng(55,"BLAKEGrostlJHKeccakSkein--RijndaelSerpentTwofishRC6MARS");

########

diag "The library thinks the CSPRNG is" . (is_csprng_well_seeded() ? " " : " not ") . "well seeded.";

{
  my @s = map { irand } 1 .. $samples;
  is( scalar(grep { $_ > 4294967295 } @s), 0, "irand values are 32-bit" );
  is( scalar(grep { $_ != int($_) } @s), 0, "irand values are integers" );
}

########

if ($use64) {
  my $bits_on  = 0;
  my $bits_off = 0;
  my $iter = 0;
  for (1 .. 6400) {
    $iter++;
    my $v = irand64;
    $bits_on |= $v;
    $bits_off |= (~$v);
    last if ~$bits_on == 0 && ~$bits_off == 0;
  }
  is( ~$bits_on,  0, "irand64 all bits on in $iter iterations" );
  is( ~$bits_off, 0, "irand64 all bits off in $iter iterations" );
}

########

{
  my $mask = 0;
  my $v;
  for (1..1024) {
    $v = drand;
    last if $v >= 1;
    next if $v < .5;
    for my $b (0..127) {
      last unless $v;
      $v *= 2;
      if ($v >= 1) {
        $mask |= (1 << $b);
        $v -= 1;
      }
    }
  }
  ok($v < 1, "drand values between 0 and 1-eps");
  my $k = 0; while ($mask) { $k++; $mask >>= 1; }
  # Assuming drand is working properly:
  #   k = 24   NV is float
  #   k = 53   NV is double
  #   k = 64   NV is long double
  # If we used drand48 we'd get 48 with double or long double.
  ok($k >= 21, "drand supplies at least 21 bits (got $k)");
}

{
  my $seed = pack("V", 0x521974A3);
  seed_csprng(length($seed), $seed);
  my $expected = drand();
  my @got;
  for my $zero (0, "0.0", "0E0") {
    seed_csprng(length($seed), $seed);
    push @got, drand($zero);
  }
  is_deeply(\@got, [($expected) x 3],
            "drand treats every numeric zero limit as an omitted limit");

  my @negative = map { drand(-10) } 1 .. 100;
  ok(!(grep { $_ > 0 || $_ <= -10 } @negative),
     "drand with a negative limit returns values in (limit,0]");
}

SKIP: {
  my $have_fork = ($Config{d_fork} || '') eq 'define';
  my $have_atfork = ($Config{d_pthread_atfork} || '') eq 'define';
  my $have_pseudofork = ($Config{d_pseudofork} || '') eq 'define';
  skip "native fork CSPRNG reseeding is not available", 4
    unless $have_fork && !$have_pseudofork && $have_atfork;

  require POSIX;
  my $seed = "fork stream regression";
  seed_csprng(length($seed), $seed);
  irand();
  my $expected = unpack("H*", random_bytes(32));
  seed_csprng(length($seed), $seed);
  irand();

  pipe(my $reader, my $writer) or die "pipe failed: $!";
  my $pid = fork();
  die "fork failed: $!" unless defined $pid;
  if ($pid == 0) {
    close $reader;
    my $child = eval { unpack("H*", random_bytes(32)) };
    $child = "ERROR: $@" unless defined $child;
    print {$writer} "$child\n";
    close $writer;
    POSIX::_exit(0);
  }

  close $writer;
  my $parent = unpack("H*", random_bytes(32));
  my $child = <$reader>;
  close $reader;
  waitpid($pid, 0);
  my $status = $?;
  chomp $child if defined $child;

  is($parent, $expected, "fork does not alter the parent CSPRNG stream");
  like($child, qr/\A[0-9a-f]{64}\z/, "child produced a valid CSPRNG stream");
  isnt($child, $expected, "child reseeds the inherited CSPRNG stream");
  is($status, 0, "child exited successfully after CSPRNG reseed");
}
