#!/usr/bin/env perl
use 5.006002;
use strict;
use warnings;
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use Getopt::Long qw(GetOptions);
use Text::ParseWords qw(shellwords);

# Called only when make needs a build. Probe by linking, never by executing,
# so this also works with cross compilers. No build-mode stamp is retained.
my ($cc, $cflags, $output, $serial) = ('cc', '-O3', 'msiqs', 0);
Getopt::Long::Configure(qw(no_auto_abbrev no_ignore_case));
GetOptions('cc=s' => \$cc, 'cflags=s' => \$cflags,
           'output=s' => \$output, 'serial' => \$serial)
  && @ARGV && length($output)
  or die "usage: build-siqs.pl [--cc COMMAND] [--cflags FLAGS] " .
         "[--output FILE] [--serial] SOURCE ...\n";
my @compiler = shellwords($cc);
die "empty compiler command\n" unless @compiler;
my @flags = shellwords($cflags);
my $threaded = 0;
unless ($serial) {
  my $temporary = tempdir('siqs-pthread-XXXXXX', TMPDIR => 1, CLEANUP => 1);
  my $probe = File::Spec->catfile($temporary, 'pthread-probe');
  my $log = File::Spec->catfile($temporary, 'probe.log');
  {
    local *STDOUT;
    local *STDERR;
    open STDOUT, '>', $log or die "cannot open $log: $!\n";
    open STDERR, '>&', \*STDOUT or die "cannot redirect probe errors: $!\n";
    $threaded = system(@compiler, @flags, '-pthread', '-o', $probe,
                      File::Spec->catfile($Bin, 'pthread-probe.c')) == 0;
  }
}

print "Building $output ", ($threaded ? 'with pthreads' : 'serial-only'), "\n";
$| = 1;
my @command = (@compiler, @flags, '-DSTANDALONE', '-UPSIQS');
push @command, '-DPSIQS', '-pthread' if $threaded;
push @command, '-o', $output, @ARGV, '-lgmp', '-lm';
print join(' ', map {shell_quote($_)} @command), "\n";
system(@command);
my $status = $?;
die "could not execute compiler: $!\n" if $status == -1;
die "SIQS build terminated by signal " . ($status & 127) . "\n"
  if $status & 127;
# A real build failure must not silently fall back to the other mode.
exit($status >> 8);

# Display a copyable shell command; execution still uses the argument list.
sub shell_quote {
  my ($argument) = @_;
  return $argument if $argument =~ m{\A[-A-Za-z0-9_./:=+,@%]+\z};
  $argument =~ s/'/'\\''/g;
  return "'$argument'";
}
