#!/usr/bin/env perl
use 5.010;
use strict;
use warnings;
use Config;
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Getopt::Long qw(GetOptions);
use Text::ParseWords qw(shellwords);

my $cc = $ENV{CC} // $Config{cc};
my ($cflags, $ldflags, $sanitize) = ('-O2 -g', '-lgmp -lm', 0);
GetOptions('cc=s' => \$cc, 'cflags=s' => \$cflags, 'ldflags=s' => \$ldflags,
           'sanitize' => \$sanitize) && !@ARGV
  or die "usage: ecm-check.pl [--sanitize] [--cc CC] [--cflags FLAGS] [--ldflags FLAGS]\n";
my $root = File::Spec->rel2abs(File::Spec->catdir($Bin, '..', '..'));
my $temporary = tempdir('ecm-check-XXXXXX', TMPDIR => 1, CLEANUP => 1);
my $checks = File::Spec->catfile($Bin, 'ecm-check.c');
my $iterator = File::Spec->catfile($root, 'prime_iterator.c');
my @flags = ('-DSTANDALONE', "-I$root", shellwords($cflags));
push @flags, '-fsanitize=address,undefined,float-cast-overflow',
             '-fno-sanitize-recover=all' if $sanitize;

sub build_run {
  my ($name, @extra) = @_;
  my $binary = File::Spec->catfile($temporary, $name);
  system(shellwords($cc), @flags, @extra, $checks, $iterator,
         shellwords($ldflags), '-o', $binary) == 0
    or die "$name compilation failed\n";
  system($binary) == 0 or die "$name failed\n";
}

build_run('ecm-check');

# Exercise 32-bit UV arithmetic on a 64-bit host, without claiming to model
# a native 32-bit ABI or allocator.  Both translation units use the same UV.
my $uv32 = File::Spec->catfile($temporary, 'uv32.h');
open my $header, '>', $uv32 or die "$uv32: $!\n";
print {$header} <<'UV32';
#include <inttypes.h>
#include "ptypes.h"
/* Use the typed-copy wrapper, not its true-standalone UV=ulong shortcut. */
#undef STANDALONE
#define UV uint32_t
#undef UV_MAX
#define UV_MAX UINT32_MAX
#undef UVCONST
#define UVCONST(n) UINT32_C(n)
#undef UVuf
#define UVuf PRIu32
#undef BITS_PER_WORD
#define BITS_PER_WORD 32
UV32
close $header or die "$uv32: $!\n";
build_run('ecm-check-uv32', '-include', $uv32);

# Model a 32-bit allocation ceiling while preserving full 64-bit UV counts.
# Allocation probes stop before Newz allocates; this is not a native ABI test.
my $alloc32 = File::Spec->catfile($temporary, 'alloc32.h');
open my $alloc_header, '>', $alloc32 or die "$alloc32: $!\n";
print {$alloc_header} <<'ALLOC32';
#include "ptypes.h"
#undef SIZE_MAX
#define SIZE_MAX UINT32_MAX
ALLOC32
close $alloc_header or die "$alloc32: $!\n";
build_run('ecm-check-alloc32', '-DECM_CHECK_ALLOC32', '-include', $alloc32);

# Exercise the independent, standalone certificate verifier's actual affine
# code, not a second copy supplied by this test.  It has no projective ECM.
my $verifier = File::Spec->catfile($root, 'examples', 'vcert.c');
open my $input, '<', $verifier or die "$verifier: $!\n";
my $source = do { local $/; <$input> };
close $input;
my ($affine) = $source =~ /(^struct ec_affine_point\b.*?)(?=^#else)/ms;
defined $affine or die "Cannot locate certificate-verifier affine code\n";
my $affine_path = File::Spec->catfile($temporary, 'verifier-affine.h');
open my $output, '>', $affine_path or die "$affine_path: $!\n";
print {$output} $affine;
close $output or die "$affine_path: $!\n";
build_run('verifier-affine-check', '-DECM_CHECK_AFFINE_ONLY', "-I$temporary");
print "All 4 ECM arithmetic checks passed.\n";
