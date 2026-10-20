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

# Compile the factor-loop prefix, unchanged, with scripted splitters.  This
# exercises high-bit scheduling without doing expensive factorizations or
# adding test seams to production code.  Everything generated is disposable.
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
