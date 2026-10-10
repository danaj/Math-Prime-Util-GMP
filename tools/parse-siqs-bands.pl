#!/usr/bin/env perl
use 5.010;
use strict;
use warnings;
use File::Basename qw(dirname);
use File::Spec;
use FindBin qw($Bin);
use Getopt::Long qw(GetOptions);

# Read the primary policy, not a saved copy of its values.  The size formulas
# below mirror siqs_select_parameters; no factoring or module build is needed.
my ($details, $describe, $help, $header_path);
Getopt::Long::Configure(qw(no_auto_abbrev no_ignore_case));
GetOptions('details' => \$details, 'describe|verbose|v' => \$describe, 'help' => \$help,
           'header=s' => \$header_path) or usage(2);
usage(0) if $help;
usage(2) if @ARGV > 1;
my $source_path = shift @ARGV // File::Spec->catfile($Bin, '..', 'siqs.c');
$header_path //= File::Spec->catfile(dirname($source_path), 'siqs.h');
my $source = read_source($source_path);
my $header = read_source($header_path);

my %define;
my $definitions = "$header\n$source";
while ($definitions =~ /^\s*#\s*define\s+(\w+)[ \t]+([^\n]+)/mg) {
  $define{$1} = $2;
}
my $number = qr/(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?/;
my $align = numeric('SIQS_SIEVE_ALIGN');
my $srange = '-';  # text for ranges, e.g. ' -> ' or '-' or '..' or '->'
die "invalid SIQS_SIEVE_ALIGN\n" unless $align > 0 && int($align) == $align;

# Extract these constants from the implementation rather than duplicating
# them.  Reject an unfamiliar formula instead of silently reporting old sizes.
my ($interval_offset, $interval_exponent) = $source =~
  /interval\s*=\s*($number)\s*\+\s*exp\(\s*($number)\s*\*\s*ln_term\s*\)/;
my ($interval_floor) = $source =~
  /if\s*\(p->half_interval\s*<\s*(\d+)\)\s*p->half_interval\s*=\s*\1\s*;/;
die "unrecognized SIQS size formulas in $source_path\n"
  unless defined($interval_offset) && defined($interval_floor) &&
         $source =~ /ln_n\s*=\s*p->bits\s*\*\s*M_LN2\s*;/ &&
         $source =~ /ln_term\s*=\s*sqrt\(ln_n\s*\*\s*log\(ln_n\)\)\s*;/ &&
         $source =~ /fb\s*=\s*exp\(p->fb_coefficient\s*\*\s*ln_term\)\s*;/;

my ($struct) = $source =~
  /typedef\s+struct\s*\{([^{}]*)\}\s*siqs_policy_band_t\s*;/s;
my ($array) = $source =~
  /\bsiqs_policy_bands\s*\[\s*\]\s*=\s*\{(.*?)\}\s*;/s;
die "cannot find siqs_policy_band_t and siqs_policy_bands in $source_path\n"
  unless defined($struct) && defined($array);
my @fields;
for my $field (split /;/, $struct) {
  next unless $field =~ /\S/;
  $field =~ /(\w+)\s*$/ or die "unrecognized band field: $field\n";
  push @fields, $1;
}
my @bands;
for my $row (split_top_level($array)) {
  $row =~ /^\s*\{(.*)\}\s*$/s or die "unrecognized band row: $row\n";
  my @values = split_top_level($1);
  die "band row has " . @values . " values, expected " . @fields . "\n"
    unless @values == @fields;
  my %band;
  for my $index (0 .. $#fields) {
    my $value = $values[$index];
    if ($fields[$index] eq 'name') {
      $value =~ /^"([^"\\]*)"$/ or die "invalid band name: $value\n";
      $value = $1;
    } elsif ($value eq 'SIQS_R_DERIVED_FROM_LP') {
      $value = 'derived';
    } elsif ($value =~ /^SIQS_POLICY_(LINEAR|RAMP)\((.*)\)$/s) {
      my ($kind, $arguments) = ($1, $2);
      my @args = split_top_level($arguments);
      die "invalid $kind policy: $value\n" unless @args == ($kind eq 'LINEAR' ? 3 : 4);
      $value = {kind => $kind, args => [map { numeric($_) } @args]};
    } elsif ($value =~ /^SIQS_POLICY_RATIO\((.*)\)$/s) {
      my @args = split_top_level($1);
      die "invalid ratio policy: $value\n" unless @args == 4;
      $value = [map { numeric($_) } @args];
    } else {
      $value = numeric($value);
    }
    $band{$fields[$index]} = $value;
  }
  die "invalid bit range in $band{name}\n"
    unless $band{first_bits} >= 1 && $band{last_bits} >= $band{first_bits};
  die "unsupported LP count in $band{name}\n"
    unless $band{production_large_primes} == 1 ||
           $band{production_large_primes} == 2;
  push @bands, \%band;
}
die "no production bands found\n" unless @bands;

my @headings = ('bits', 'type', 'q', 'LP', 'M', 'FB', 'res');
push @headings, ('bias', 'headroom', 'pmin floor', 'rel +') if $details;
my @rows;
for my $band (@bands) {
  my ($first, $last) = @{$band}{qw(first_bits last_bits)};
  my @sizes = map { sizes_at($band, $_) } ($first, $last);
  my ($lp, $residual, $type);
  my @k = map { curve_at($band->{k_l}, $_) } ($first, $last);
  die "invalid K_L in $band->{name}\n" if grep { $_ < 1 } @k;
  $type = $band->{production_large_primes} == 2 ? '2LP'
        : $k[0] == 1 && $k[1] == 1 ? 'smooth' : '1LP';
  $lp = coefficient_range('P', @k);
  if (!ref($band->{k_r}) && $band->{k_r} eq 'derived') {
    $residual = $type eq '2LP' ? 'LP*P' : 'LP';
  } else {
    my @r = map { curve_at($band->{k_r}, $_) } ($first, $last);
    die "invalid K_R in $band->{name}\n" if grep { $_ <= 0 } @r;
    $residual = coefficient_range($type eq '2LP' ? 'P^2' : 'P', @r);
  }
  my @row = (sprintf('%3d-%-3d', $first, $last), $type, $band->{q_count},
             $lp, range_text(map { $_->[0] } @sizes),
             range_text(map { $_->[1] } @sizes), $residual);
  if ($details) {
    my @bias = map {
      $band->{stage1_bias_base} + ($band->{stage1_bias_step_bits}
        ? int(($_ - $band->{stage1_bias_origin_bits}) /
              $band->{stage1_bias_step_bits}) : 0)
    } ($first, $last);
    push @row, range_text(@bias), $band->{sieve_free_units},
      $band->{sieve_start_prime_floor} || '-', $band->{relation_extra};
  }
  push @rows, \@row;
}
my @width = map { length($_) } @headings;
for my $row (@rows) {
  for my $column (0 .. $#width) {
    $width[$column] = length($row->[$column])
      if length($row->[$column]) > $width[$column];
  }
}
my $format = join(' | ', map { "%-${_}s" } @width) . "\n";
if ($describe) {
  say 'bits:   Bit range of the post-trial cofactor N, before multiplying by k.';
  say 'type:   Smooth, or up to one/two large primes per raw relation (1LP/2LP).';
  say 'q:      Number of distinct factor-base primes used in polynomial A.';
  say 'LP:     Maximum individual large prime, from K_L times the actual prime P.';
  say 'M:      Half the sieve interval; the sieve covers 2M positions.';
  say 'FB:     Number of primes used in the factor base.';
  say 'res:    Maximum cofactor after factor-base division (a product in 2LP).';
  if ($details) {
    say 'bias:           Additional coarse candidate-filter headroom (sieve-byte units).';
    say 'headroom:       Sieve-byte initialization headroom (sieve_free_units).';
    say 'pmin floor:     Lower bound for the first sieved prime, not its actual value.';
    say 'rel +:          Initial relation surplus above FB + 1; readiness may stop earlier.';
  }
  say '';
}
printf $format, @headings;
say join('-+-', map { '-' x $_ } @width);
printf $format, @$_ for @rows;
if ($details) {
  say "\nM and FB show first -> last endpoint values, not sorted extrema; sieve length is 2M.";
  say 'P is the largest factor-base prime; smooth rows allow no primes above P.';
  say 'LP coefficients show first -> last values of K_L; products are rounded down.';
  say 'LP is capped at SIQS_LP_MAX before coupled R is derived: LP in 1LP, LP*P in 2LP.';
  say 'Explicit K_R uses P in 1LP or P^2 in 2LP; R is capped at SIQS_RESIDUAL_PRODUCT_MAX.';
  say 'Primary production bands only; low-end recovery profiles are not shown.';
  say 'At the unused 1-bit SIQS endpoint, the FB floor is shown (pretests return first).'
  if $bands[0]{first_bits} == 1;
    say 'pmin floor is a lower bound, not the actual first sieved prime; rel + is the initial surplus.'
}

sub usage {
  my ($status) = @_;
  my $out = $status ? *STDERR : *STDOUT;
  print {$out} "Usage: perl $0 [--describe] [--details] [--header siqs.h] [siqs.c]\n",
    "  Default source: siqs.c beside the repository's misc directory.\n",
    "  --describe prints column descriptions above the table (aliases: --verbose, -v).\n",
    "  --details adds byte-sieve bias/headroom, prime floor, and relation surplus.\n",
    "  --header supplies the header defining the bit limits (default: beside siqs.c).\n";
  exit $status;
}

sub read_source {
  my ($path) = @_;
  open my $in, '<', $path or die "$path: $!\n";
  local $/;
  my $text = <$in>;
  close $in or die "$path: $!\n";
  # Preserve quoted strings while removing C comments and line continuations.
  $text =~ s{("(?:\\.|[^"\\])*"|'(?:\\.|[^'\\])*')|/\*.*?\*/|//[^\n]*}
            {defined($1) ? $1 : ' '}gse;
  $text =~ s/\\\r?\n/ /g;
  return $text;
}

sub numeric {
  my ($expression, $depth) = @_;
  die "recursive numeric macro: $expression\n" if ($depth // 0) > 20;
  # Only numeric tokens/operators reach eval; source text cannot execute Perl.
  $expression =~ s/\b(0[xX][0-9a-fA-F]+|\d+)[uUlL]+\b/$1/g;
  $expression =~ s{\b([A-Za-z_]\w*)\b}{
    exists($define{$1}) ? '(' . numeric($define{$1}, ($depth // 0) + 1) . ')'
                       : die("unknown numeric macro: $1\n")}ge;
  $expression =~ s/\b0[xX]([0-9a-fA-F]+)\b/hex($1)/ge;
  die "unsupported numeric expression: $expression\n"
    unless $expression =~ /\A(?:$number|[\s()+*\/<>-])+\z/;
  my $value = eval $expression;
  die "invalid numeric expression: $expression ($@)\n" if $@ || !defined($value);
  die "nonfinite numeric expression: $expression\n"
    unless $value == $value && abs($value) <= 1.7976931348623157e308;
  return $value;
}

sub split_top_level {
  my ($text) = @_;
  my (@parts, @close);
  my %matching = ('(' => ')', '{' => '}', '[' => ']');
  my $start = 0;
  while ($text =~ /"(?:\\.|[^"\\])*"|[(){}\[\],]/g) {
    my $token = $&;
    if (exists $matching{$token}) {
      push @close, $matching{$token};
    } elsif ($token =~ /^[)}\]]$/) {
      die "unbalanced policy expression\n" unless @close && pop(@close) eq $token;
    } elsif ($token eq ',' && !@close) {
      push @parts, substr($text, $start, $-[0] - $start);
      $start = $+[0];
    }
  }
  die "unbalanced policy expression\n" if @close;
  push @parts, substr($text, $start);
  s/^\s+|\s+$//g for @parts;
  pop @parts if @parts && $parts[-1] eq '';
  return @parts;
}

sub curve_at {
  my ($curve, $bits) = @_;
  die "derived R is not a numeric curve\n" unless ref($curve) eq 'HASH';
  my ($begin, $step_or_end, $first, $last) = @{$curve->{args}};
  return $begin + $step_or_end * ($bits - $first) if $curve->{kind} eq 'LINEAR';
  die "invalid endpoint policy ramp\n"
    if $first > $last || ($first == $last && $begin != $step_or_end);
  die "cannot extrapolate policy ramp to $bits bits\n" if $bits < $first || $bits > $last;
  return $begin if $bits == $first;
  return $step_or_end if $bits == $last;
  return $begin + ($step_or_end - $begin) * (($bits - $first) / ($last - $first));
}

sub sizes_at {
  my ($band, $bits) = @_;
  # ln(N)*ln(ln(N)) is negative at one bit, but SIQS never runs there.
  my $ln_n = $bits * numeric('M_LN2');
  my $ln_term = $ln_n >= 1 ? sqrt($ln_n * log($ln_n)) : 0;
  my $fb = exp(curve_at($band->{fb_coefficient}, $bits) * $ln_term);
  $fb = $band->{fb_floor} if $fb < $band->{fb_floor};
  my $m = $band->{fixed_half_interval};
  if (!$m) {
    my ($amount, $numerator, $step, $denominator) = @{$band->{interval_lp_scale}};
    my $ratio = $denominator
      ? 1 + $amount * ($numerator + $step * ($bits - $band->{first_bits})) / $denominator
      : 1;
    my $scale = curve_at($band->{interval_base_scale}, $bits) * $ratio;
    my $interval = ($interval_offset + exp($interval_exponent * $ln_term)) * $scale;
    die "negative interval in $band->{name} at $bits bits\n" if $interval < 0;
    $m = int((int($interval) + $align - 1) / $align) * $align;
    $m = $interval_floor if $m < $interval_floor;
  }
  return [$m, int($fb)];
}

sub coefficient_range {
  my ($unit, $first, $last) = @_;
  my @text = map { $_ == 1 ? "" : sprintf('%.6g', $_) } ($first, $last);
  my $out = $text[0] eq $text[1] ? $text[0] : "$text[0]${srange}$text[1]";
  return $out.$unit;
}

sub integer_text {
  my ($value) = @_;
  my $text = sprintf('%.0f', $value);
  $text =~ s/(?<=\d)(?=(?:\d{3})+$)/,/g;
  return $text;
}

sub range_text {
  my ($first, $last) = map { integer_text($_) } @_;
  return $first eq $last ? $first : "$first${srange}$last";
}
