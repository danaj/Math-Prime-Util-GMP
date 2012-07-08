package Math::Prime::Util::GMP;
use strict;
use warnings;
use Carp qw/croak confess carp/;

BEGIN {
  $Math::Prime::Util::GMP::AUTHORITY = 'cpan:DANAJ';
  $Math::Prime::Util::GMP::VERSION = '0.01';
}

# parent is cleaner, and in the Perl 5.10.1 / 5.12.0 core, but not earlier.
# use parent qw( Exporter );
use base qw( Exporter );
our @EXPORT_OK = qw(
                     is_strong_pseudoprime
                   );
our %EXPORT_TAGS = (all => [ @EXPORT_OK ]);

BEGIN {

  eval {
    require XSLoader;
    XSLoader::load(__PACKAGE__, $Math::Prime::Util::GMP::VERSION);
    1;
  } or do {
    die $@;
  }
}

sub is_strong_pseudoprime {
  my($n) = shift;
  foreach my $base (@_) {
    return 0 unless miller_rabin("$n", "$base");
  }
  1;
}
1;

__END__


# ABSTRACT: Utilities related to prime numbers, using GMP

=pod

=encoding utf8


=head1 NAME

Math::Prime::Util::GMP - Utilities related to prime numbers, using GMP


=head1 VERSION

Version 0.01


=head1 SYNOPSIS

  use Math::Prime::Util::GMP ':all';

  my $n = "115792089237316195423570985008687907853269984665640564039457584007913129639937";
  say "$n is a prime or spsp-2" if is_strong_pseudoprime($n, 2);


=head1 DESCRIPTION

A set of utilities related to prime numbers, using GMP.


=head1 FUNCTIONS

=head2 is_strong_pseudoprime

  my $maybe_prime = is_strong_pseudoprime($n, 2);
  my $probably_prime = is_strong_pseudoprime($n, 2, 3, 5, 7, 11, 13, 17);

Takes a positive number as input and one or more bases.  The bases must be
between C<2> and C<n - 2>.  Returns 1 is C<n> is a prime or a strong
pseudoprime to all of the bases, and 0 if not.

If 0 is returned, then the number really is a composite.  If 1 is returned,
then it is either a prime or a strong pseudoprime to all the given bases.
Given enough distinct bases, the chances become very, very strong that the
number is actually prime.

This is usually used in combination with other tests to make either stronger
tests (e.g. the strong BPSW test) or deterministic results for numbers less
than some verified limit (e.g. it has long been known that no more than three
selected bases are required to give correct primality test results for any
32-bit number).  Given the small chances of passing multiple bases, there
are some math packages that just use multiple MR tests for primality testing.

Even numbers other than 2 will always return 0 (composite).  While the
algorithm does run with even input, most sources define it only on odd input.
Returning composite for all non-2 even input makes the function match most
other implementations including L<Math::Primality>'s C<is_strong_pseudoprime>
function.

=head1 SEE ALSO

  L<Math::Prime::Util>

=head1 AUTHORS

Dana Jacobsen E<lt>dana@acm.orgE<gt>


=head1 COPYRIGHT

Copyright 2011-2012 by Dana Jacobsen E<lt>dana@acm.orgE<gt>

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

=cut
