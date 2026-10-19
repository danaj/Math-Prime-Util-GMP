#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 1;
use FindBin qw($Bin);
use File::Spec;

my $checker = File::Spec->catfile($Bin, '..', 'tools', 'factor-check.pl');
open my $run, '-|', $^X, $checker or die "Cannot start factor checker: $!";
my $output = do { local $/; <$run> };
my $ok = close $run;
diag($output) if defined $output && length $output;
ok($ok, 'scripted factor ladder, SIQS helper and cofactor bookkeeping');
