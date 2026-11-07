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

# Compile the factor-loop prefix unchanged with scripted splitters, plus
# the actual p-1 code for boundary checks.  This exercises scheduling without
# expensive factorizations or production test seams. Generated files are disposable.
my $cc = $ENV{CC} // $Config{cc};
my ($cflags, $ldflags, $sanitize) = ('-O2 -g', '-lgmp -lm', 0);
GetOptions('cc=s' => \$cc, 'cflags=s' => \$cflags, 'ldflags=s' => \$ldflags,
           'sanitize' => \$sanitize) && !@ARGV
  or die "usage: factor-check.pl [--sanitize] [--cc CC] [--cflags FLAGS] [--ldflags FLAGS]\n";
my $root = File::Spec->rel2abs(File::Spec->catdir($Bin, '..', '..', '..'));
my $temporary = tempdir('factor-check-XXXXXX', TMPDIR => 1, CLEANUP => 1);
my $source_path = File::Spec->catfile($root, 'factor.c');
open my $input, '<', $source_path or die "$source_path: $!\n";
my $source = do { local $/; <$input> };
close $input;
my ($pm1_small) = $source =~ /(^#define TEST_FOR_2357\b.*?)(?=^int )/ms;
my ($pm1_source) = $source =~ /(^int _GMP_pminus1_factor\b.*?^\})/ms;
defined $pm1_small && defined $pm1_source
  or die "Cannot locate the p-1 implementation\n";
$source =~ s/\nstatic void _omega\b.*\z/\n/s
  or die "Cannot locate the end of the factor-loop prefix\n";
my $ecm_path = File::Spec->catfile($root, 'ecm.c');
open my $ecm_input, '<', $ecm_path or die "$ecm_path: $!\n";
my $ecm_source = do { local $/; <$ecm_input> };
close $ecm_input;
my ($bounds) = $ecm_source =~ /(^static UV ecm_resolve_b2\b.*?^\})/ms;
defined $bounds or die "Cannot locate the ECM bound resolver\n";
my $generated = File::Spec->catfile($temporary, 'factor-check.c');
my $checks = File::Spec->catfile($Bin, 'factor-check.inc.c');
open my $output, '>', $generated or die "$generated: $!\n";
# Test-only exception/allocation hooks; production factor.c stays unchanged.
print {$output} <<'CHECK_HOOKS';
#include <stdarg.h>
#include <stdio.h>
#include <gmp.h>
#include "ptypes.h"
#include <setjmp.h>
static void check_croak(const char *format, ...);
static int check_printf(const char *format, ...);
#undef croak
#define croak(...) check_croak(__VA_ARGS__)
#undef gmp_printf
#define gmp_printf check_printf
CHECK_HOOKS
print {$output} $source, "\n#include <stdint.h>\n";
print {$output} $pm1_small;
print {$output} "\nstatic unsigned check_pm1_stage2_terms;\n";
for my $mode ('list', 'iterator') {
  my $pm1 = $pm1_source;
  $pm1 =~ s/\b_GMP_pminus1_factor\b/check_pm1_$mode/g;
  $pm1 =~ s/B2 < 10000000/0/ if $mode eq 'iterator';
  $pm1 =~ s/UV lastq, qdiff;/UV lastq, qdiff; check_pm1_stage2_terms++;/;
  print {$output} "\n", $pm1, "\n";
}
# Model GMP's narrower unsigned-long arguments on a UV64 host.
print {$output} <<'CHECK_UI32';
#if BITS_PER_WORD == 64
static void check_mul_ui32(mpz_ptr r, mpz_srcptr a, unsigned long k) {
  mpz_mul_ui(r, a, (uint32_t)k);
}
static void check_powm_ui32(mpz_ptr r, mpz_srcptr a, unsigned long k, mpz_srcptr n) {
  mpz_powm_ui(r, a, (uint32_t)k, n);
}
CHECK_UI32
my $pm1_ui32 = $pm1_source;
$pm1_ui32 =~ s/\b_GMP_pminus1_factor\b/check_pm1_ui32/g;
$pm1_ui32 =~ s/\bULONG_MAX\b/UINT32_MAX/g;
$pm1_ui32 =~ s/sizeof\(unsigned long\)/sizeof(uint32_t)/g;
$pm1_ui32 =~ s/\bmpz_mul_ui\b/check_mul_ui32/g;
$pm1_ui32 =~ s/\bmpz_powm_ui\b/check_powm_ui32/g;
print {$output} "\n", $pm1_ui32, "\n#endif\n";
# Extract the actual XS default expression, rather than duplicating its policy.
open my $xs_input, '<', File::Spec->catfile($root, 'XS.xs') or die "$!\n";
my $xs_source = do { local $/; <$xs_input> };
close $xs_input;
my ($pm1_b2) = $xs_source =~ /case 3:\s*if \(arg2 == 0\)\s*arg2 = ([^;]+);\s*success = _GMP_pminus1_factor/;
defined $pm1_b2 or die "Cannot locate the p-1 default B2 expression\n";
for my $bits (32, 64) {
  my $expression = $pm1_b2;
  $expression =~ s/\bUV_MAX\b/UINT${bits}_MAX/g;
  print {$output} "static uint${bits}_t check_pm1_b2_$bits(uint${bits}_t arg1) { return $expression; }\n";
}
# Exercise the actual resolver with both UV widths, even on a 64-bit host.
for my $bits (32, 64) {
  my $typed_bounds = $bounds;
  $typed_bounds =~ s/\becm_resolve_b2\b/ecm_resolve_b2_$bits/g;
  $typed_bounds =~ s/\bUV_MAX\b/UINT${bits}_MAX/g;
  $typed_bounds =~ s/\bUV\b/uint${bits}_t/g;
  print {$output} "\n", $typed_bounds, "\n";
}
print {$output} "\n#include \"", escaped($checks), "\"\n";
close $output or die "$generated: $!\n";
my $binary = File::Spec->catfile($temporary, 'factor-check');
my @sanitizers = $sanitize ? ('-fsanitize=address,undefined', '-fno-omit-frame-pointer') : ();
my @command = (shellwords($cc), shellwords($cflags), @sanitizers,
               '-DSTANDALONE', '-I'.$root, '-o', $binary, $generated,
               File::Spec->catfile($root, 'prime_iterator.c'), shellwords($ldflags));
my $log = File::Spec->catfile($temporary, 'build.log');
my $status;
{
  local *STDOUT;
  local *STDERR;
  open STDOUT, '>', $log or die "$log: $!\n";
  open STDERR, '>&', \*STDOUT or die "redirect compiler output: $!\n";
  $status = system(@command);
}
if ($status != 0) {
  open my $errors, '<', $log or die "$log: $!\n";
  my $lines = 0;
  while (my $line = <$errors>) { print STDERR $line; last if ++$lines == 35; }
  die "factor-check compilation failed\n";
}
$status = system { $binary } $binary;
die "factor-check failed (status $status)\n" if $status != 0;

sub escaped {
  my ($path) = @_;
  $path =~ s/\\/\\\\/g;
  $path =~ s/"/\\"/g;
  return $path;
}
