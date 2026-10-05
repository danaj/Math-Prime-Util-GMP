#!/usr/bin/env perl
use 5.006002;
use strict;
use warnings;
use File::Spec;
use FindBin qw($Bin);
use Getopt::Long qw(GetOptions);

# Compilation deliberately lives in the Makefile, not in this wrapper.
# Only core Perl modules are needed. Keep measured output off stdout in the
# flags-only mode so a caller can safely capture the complete result.
my $root = File::Spec->rel2abs(File::Spec->catdir($Bin, '..'));
my $bench = File::Spec->catfile($root, 'siqs-sieve-bench');
my ($threads, $sizes, $seconds, $output, $flags_only, $flags_progress, $help) =
   (1, '0,16,32,48,64,96,128', 2, undef, 0, 0, 0);
Getopt::Long::Configure(qw(no_auto_abbrev no_ignore_case));
GetOptions('threads=i' => \$threads, 'sizes=s' => \$sizes,
           'work-seconds=f' => \$seconds, 'output=s' => \$output,
           'bench=s' => \$bench,
           'flags-only|onlyflags' => \$flags_only,
           'flags' => \$flags_progress, 'help' => \$help)
  or die "Try --help for usage.\n";
if ($help) {
  print <<'HELP';
usage: tools/siqs-autotune.pl [options] [INTEGER ...]
  --threads N       Intended sieve concurrency (default 1).
  --sizes LIST      Block maxima in KiB (default 0,16,32,48,64,96,128).
  --work-seconds S  Serial collection probe (default 2; 0 skips it).
  --output FILE     Raw block-screen TSV; refuses existing files.
  --flags-only      Print only suggested -D flags (alias --onlyflags).
  --flags           Flags on stdout, concise progress/decision on stderr.
  --bench FILE      Alternative siqs-sieve-bench executable.

Build first with 'make siqs-sieve-bench'. Defaults are RSA-100/RSA-110
(330/364 bits); INTEGER arguments replace them (up to 16 distinct inputs).
Combine per-input relative timings with an equal-weight geometric mean.
Screen for 200 ms serial, 400 ms at --threads, then two 400 ms final sweeps.
Close/noisy aggregate results receive at most two additional 1-second sweeps.
The work-rate probe uses the first input's serial collection pipeline, with
no solving. Run on an idle machine.
No files or build settings are changed, except an explicitly requested TSV.

'make siqs-tuned' measures and rebuilds msiqs with the suggested flags.
Use SIQS_TUNE_ARGS='--threads 8' to screen at another concurrency.
See tools/README-siqs-autotune.txt for caveats and manual build examples.
HELP
  exit 0;
}
die "Choose --flags or --onlyflags, not both.\n" if $flags_only && $flags_progress;
die "Supply at most 16 positive decimal INTEGER arguments.\n"
  if @ARGV > 16 || grep { !/\A[0-9]+\z/ } @ARGV;
die "--threads must be positive.\n" unless $threads > 0;
die "--work-seconds must be in [0,60].\n"
  unless $seconds >= 0 && $seconds <= 60;
$bench = File::Spec->rel2abs($bench);
die "Cannot execute $bench; build with 'make siqs-sieve-bench' first.\n"
  unless -x $bench;
my @command = ($bench, '--autotune', '--threads', $threads,
               '--sizes', $sizes, '--work-seconds', $seconds);
push @command, '--output', $output if defined $output;
push @command, @ARGV;
$| = 1;
open my $pipe, '-|', @command or die "Cannot start $bench: $!\n";
my $flags;
while (my $line = <$pipe>) {
  if ($line =~ /\ASIQS_AUTOTUNE_FLAGS\t([^\r\n]*)\r?\n\z/) {
    die "Duplicate autotune result.\n" if defined $flags;
    $flags = $1;
  } elsif ($line =~ /\ASIQS_AUTOTUNE_PROGRESS\t(.*)\r?\n\z/) {
    if ($flags_progress) { print STDERR "$1\n"; }
    elsif (!$flags_only) { print "$1\n"; }
  } else {
    print $line unless $flags_only || $flags_progress;
  }
}
close $pipe;
my $status = $?;
die "Autotune terminated by signal " . ($status & 127) . ".\n" if $status & 127;
die "Autotune failed (exit " . ($status >> 8) . "); no flags emitted.\n"
  if $status != 0;
die "Autotune produced no usable recommendation; use a larger input.\n"
  unless defined($flags) && length($flags);
# Treat the helper's output as data, not a shell command. In particular, failed
# or unexpected output must not get silently incorporated into a tuned build.
my %seen;
for my $flag (split / /, $flags) {
  die "Unexpected autotune flag: $flag\n"
    unless $flag =~ /\A-D(SIQS_SIEVE_BLOCK_SIZE)=([0-9]+)U\z/ ||
           $flag =~ /\A-D(SIQS_PROGRESS_WORK_PER_SEC)=([0-9]+)ULL\z/;
  my ($name, $value) = ($1, $2);
  die "Duplicate autotune flag: $name\n" if $seen{$name}++;
  die "Invalid autotune block size: $value\n"
    if $name eq 'SIQS_SIEVE_BLOCK_SIZE' &&
       $value != 0 && ($value < 4096 || $value > 1048576);
  die "Invalid autotune work rate: $value\n"
    if $name eq 'SIQS_PROGRESS_WORK_PER_SEC' && ($value < 1 || $value > 1e15);
}
print "Suggested compiler flags:\n  " unless $flags_only || $flags_progress;
print "$flags\n";
