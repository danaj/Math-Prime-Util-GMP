use strict;
use warnings;
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($Bin);
use IPC::Open3;
use Symbol qw(gensym);
use Test::More;

# Optional wrapper protocol checks: no compiler, GMP, timing or factoring.
# The fixture output is small, so sequential pipe reads cannot fill a pipe.
my $dir = tempdir('siqs-autotune-test-XXXXXX', TMPDIR => 1, CLEANUP => 1);
my $fixture = File::Spec->catfile($dir, 'bench');
my $tool = File::Spec->catfile($Bin, '..', 'tools', 'siqs-autotune.pl');
my $flags = '-DSIQS_SIEVE_BLOCK_SIZE=65536U -DSIQS_PROGRESS_WORK_PER_SEC=600000000ULL';
open my $f, '>', $fixture or die "$fixture: $!";
print {$f} '#!', $^X, "\n", <<'FIXTURE';
use strict;
use warnings;
my $case = $ENV{MPU_SIQS_AUTOTUNE_TEST_CASE};
print "measured stages\n";
print "SIQS_AUTOTUNE_PROGRESS\tTesting block stage 1: 0 16 32 KiB\n";
print "SIQS_AUTOTUNE_PROGRESS\tTied: choosing 96k over 128k\n";
exit 2 if $case eq 'missing';
my $flags = '-DSIQS_SIEVE_BLOCK_SIZE=65536U -DSIQS_PROGRESS_WORK_PER_SEC=600000000ULL';
$flags = '-DUNEXPECTED=1' if $case eq 'unknown';
$flags = '-DSIQS_SIEVE_BLOCK_SIZE=1U' if $case eq 'bad-size';
$flags = '-DSIQS_PROGRESS_WORK_PER_SEC=0ULL' if $case eq 'bad-rate';
$flags = '-DSIQS_SIEVE_BLOCK_SIZE=65536U -DSIQS_SIEVE_BLOCK_SIZE=32768U'
  if $case eq 'duplicate-flag';
$flags = '' if $case eq 'empty';
print "SIQS_AUTOTUNE_FLAGS\t$flags\n";
print "SIQS_AUTOTUNE_FLAGS\t$flags\n" if $case eq 'duplicate-result';
exit($case eq 'failed-after-flags' ? 1 : 0);
FIXTURE
close $f or die "$fixture: $!";
chmod 0755, $fixture or die "chmod $fixture: $!";

sub run_tool {
  my ($case, @options) = @_;
  local $ENV{MPU_SIQS_AUTOTUNE_TEST_CASE} = $case;
  my $err = gensym;
  my $pid = open3(undef, my $out, $err, $^X, $tool, '--bench', $fixture, @options);
  local $/;
  my $stdout = <$out>;
  my $stderr = <$err>;
  waitpid $pid, 0;
  return ($?, defined($stdout) ? $stdout : '', defined($stderr) ? $stderr : '');
}

for my $option ('--flags-only', '--onlyflags') {
  my ($status, $out, $err) = run_tool('success', $option);
  is($status, 0, "$option succeeds");
  is($out, "$flags\n", "$option prints exactly one flags line");
  is($err, '', "$option has no extra text on stderr");
}
my ($status, $out, $err) = run_tool('success');
is($status, 0, 'normal output succeeds');
like($out, qr/measured stages/, 'normal output retains diagnostics');
like($out, qr/Suggested compiler flags:\n  \Q$flags\E\n\z/, 'normal output shows flags');
unlike($out, qr/SIQS_AUTOTUNE_(?:FLAGS|PROGRESS)/, 'internal protocol tags are hidden');
{
  my ($success, $stdout, $stderr) = run_tool('success', '--flags');
  is($success, 0, '--flags succeeds');
  is($stdout, "$flags\n", '--flags stdout contains only compiler flags');
  is($stderr, "Testing block stage 1: 0 16 32 KiB\nTied: choosing 96k over 128k\n",
    '--flags stderr contains concise progress and decision only');
}
{
  my ($failure, $stdout, $stderr) = run_tool('failed-after-flags', '--flags');
  isnt($failure, 0, '--flags rejects a failure after the result tag');
  is($stdout, '', '--flags emits no partial flags on failure');
  like($stderr, qr/Testing block stage 1/, 'progress remains visible on failure');
  like($stderr, qr/Autotune failed/, 'failure is explained after progress');
}
{
  my ($failure, $stdout, $stderr) = run_tool('success', '--flags', '--onlyflags');
  isnt($failure, 0, 'conflicting output modes fail');
  is($stdout, '', 'conflicting modes emit no compiler flags');
}
for my $case (qw(missing failed-after-flags unknown bad-size bad-rate
                duplicate-flag duplicate-result empty)) {
  my ($failure, $stdout, $stderr) = run_tool($case, '--flags-only');
  isnt($failure, 0, "$case fails");
  is($stdout, '', "$case emits no partial flags");
  isnt($stderr, '', "$case explains its failure");
}
done_testing();
