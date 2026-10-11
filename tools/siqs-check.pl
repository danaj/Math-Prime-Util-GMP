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

# Build a standalone checker in a disposable directory.  No blib, XS build,
# MPU modules, or source rewriting.  The C executable owns suite selection.
my ($extended, $verbose, $sanitize, $debug, $keep, $list, $help, $block, $threaded);
my $suite = 'all';
my $cc = $ENV{CC} // $Config{cc};
my $cflags = '-O3';
my $ldflags = '-lgmp -lm';
Getopt::Long::Configure(qw(no_auto_abbrev no_ignore_case));
GetOptions('suite=s' => \$suite, 'extended' => \$extended, 'verbose' => \$verbose,
  'sanitize' => \$sanitize, 'debug' => \$debug, 'keep' => \$keep,
  'list' => \$list, 'help' => \$help, 'block-size=i' => \$block,
  'threaded' => \$threaded,
  'cc=s' => \$cc, 'cflags=s' => \$cflags, 'ldflags=s' => \$ldflags)
  or usage(2);
usage(0) if $help;
usage(2) if @ARGV;
die "block size must be zero or 4096..1048576 bytes\n"
  if defined($block) && $block != 0 && ($block < 4096 || $block > 1048576);
my $temporary = tempdir('siqs-check-XXXXXX', TMPDIR => 1, CLEANUP => !$keep);
my $root = File::Spec->rel2abs(File::Spec->catdir($Bin, '..'));
my $binary = File::Spec->catfile($temporary, 'siqs-check');
my @command = (shellwords($cc), shellwords($cflags));
push @command, '-O1', '-g', '-fsanitize=address,undefined', '-fno-sanitize-recover=all'
  if $sanitize;
push @command, '-DSTANDALONE', "-I$root";
push @command, '-DSIQS_DEBUG' if $debug;
push @command, '-DPSIQS', '-pthread' if $threaded;
push @command, "-DSIQS_SIEVE_BLOCK_SIZE=$block" if defined $block;
push @command, '-o', $binary, map {File::Spec->catfile($root, $_)}
  qw(tools/siqs-check.c prime_iterator.c squfof126.c);
push @command, shellwords($ldflags);
print "Building standalone SIQS checker", ($sanitize ? ' (ASan/UBSan)' : ''), "\n";
$| = 1;
system(@command) == 0 or die "SIQS checker build failed\n";
my @options = ('--suite', $suite);
push @options, '--extended' if $extended;
push @options, '--verbose' if $verbose;
push @options, '--list' if $list;
system($binary, @options);
my $status = $?;
print "Checker retained in $temporary\n" if $keep;
die "could not execute SIQS checker: $!\n" if $status == -1;
die "SIQS checker terminated by signal " . ($status & 127) . "\n" if $status & 127;
exit($status >> 8);

sub usage {
  my ($status) = @_;
  print <<'USAGE';
usage: perl tools/siqs-check.pl [options]
  --suite NAME          Run policies, sieve, relations, matrix, workers, cofactors, or all.
  --extended           Broader fixtures and more polynomials, not full factors.
  --verbose            Report individual polynomial and matrix fixtures.
  --threaded           Also check pthread Lanczos and worker pools.
  --block-size N       Build with another block maximum (0 disables blocking).
  --sanitize           Enable ASan/UBSan; compiler/runtime support required.
  --debug              Also enable production SIQS_DEBUG assertions.
  --cc COMMAND         Compiler command (default CC or Perl's configured cc).
  --cflags FLAGS       Compiler flags (default -O3).
  --ldflags FLAGS      Linker flags (default -lgmp -lm).
  --keep               Keep the temporary binary for further checks.
  --list               List the suites in the compiled checker.
USAGE
  exit $status;
}
