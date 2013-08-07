package Math::Prime::Util::GMP;
use strict;
use warnings;
use Carp qw/croak confess carp/;

BEGIN {
  $Math::Prime::Util::GMP::AUTHORITY = 'cpan:DANAJ';
  $Math::Prime::Util::GMP::VERSION = '0.14';
}

# parent is cleaner, and in the Perl 5.10.1 / 5.12.0 core, but not earlier.
# use parent qw( Exporter );
use base qw( Exporter );
our @EXPORT_OK = qw(
                     is_prime
                     is_prob_prime
                     is_provable_prime
                     is_provable_prime_with_cert
                     is_aks_prime
                     is_nminus1_prime
                     is_ecpp_prime
                     is_strong_pseudoprime
                     is_lucas_pseudoprime
                     is_strong_lucas_pseudoprime
                     is_extra_strong_lucas_pseudoprime
                     is_almost_extra_strong_lucas_pseudoprime
                     is_frobenius_underwood_pseudoprime
                     lucas_sequence
                     primes
                     next_prime
                     prev_prime
                     trial_factor
                     prho_factor
                     pbrent_factor
                     pminus1_factor
                     pplus1_factor
                     holf_factor
                     squfof_factor
                     ecm_factor
                     factor
                     prime_count
                     primorial
                     pn_primorial
                     consecutive_integer_lcm
                   );
                   # Should add:
                   # nth_prime
our %EXPORT_TAGS = (all => [ @EXPORT_OK ]);

BEGIN {
  eval {
    require XSLoader;
    XSLoader::load(__PACKAGE__, $Math::Prime::Util::GMP::VERSION);
    _GMP_init();
    1;
  } or do {
    die $@;
  }
}
END {
  _GMP_destroy();
}

sub _validate_positive_integer {
  my($n, $min, $max) = @_;
  croak "Parameter must be defined" if !defined $n;
  if (ref($n) eq 'Math::BigInt' && $n->can("sign")) {
    croak "Parameter '$n' must be a positive integer" unless $n->sign() eq '+';
  } else {
    croak "Parameter '$n' must be a positive integer"
          if $n eq '' || $n =~ tr/0123456789//c;
  }
  croak "Parameter '$n' must be >= $min" if defined $min && $n < $min;
  croak "Parameter '$n' must be <= $max" if defined $max && $n > $max;
  1;
}


sub is_strong_pseudoprime {
  my($n, @bases) = @_;
  _validate_positive_integer($n);
  croak "No bases given to is_strong_pseudoprime" unless @bases;
  foreach my $base (@bases) {
    _validate_positive_integer($base);
    return 0 unless _GMP_miller_rabin("$n", "$base");
  }
  1;
}

sub is_provable_prime {
  my ($n) = @_;
  return 0 if $n < 2;
  return _is_provable_prime($n);
}

sub is_provable_prime_with_cert {
  my ($n) = @_;
  my @composite = (0, '');
  return @composite if $n < 2;

  my ($result, $text) = _is_provable_prime($n, 1);
  return @composite if $result == 0;
  return ($result, '') if $result != 2;
  $text = "Type Small\nN $n\n" if !defined $text || $text eq '';
  $text =~ s/\n$//;
  $text = "[MPU - Primality Certificate]\nVersion 1.0\n\nProof for:\nN $n\n\n$text";
  return ($result, $text);
}

sub factor {
  my ($n) = @_;
  return ($n) if $n < 4;

  my @factors = sort {$a<=>$b} _GMP_factor($n);
  return @factors;
}

sub primes {
  my $optref = (ref $_[0] eq 'HASH')  ?  shift  :  {};
  croak "no parameters to primes" unless scalar @_ > 0;
  croak "too many parameters to primes" unless scalar @_ <= 2;
  my $low = (@_ == 2)  ?  shift  :  2;
  my $high = shift;
  my $sref = [];

  _validate_positive_integer($low);
  _validate_positive_integer($high);

  return $sref if ($low > $high) || ($high < 2);

  # Simple trial method for now.
  return _GMP_trial_primes($low, $high);

  # Trial primes without the XS code.  Works fine and is a lot easier than the
  # XS code (duh -- it's Perl).  But 30-40% slower, mostly due to lots of
  # string -> mpz -> string conversions and little memory allocations.
  #
  #my @primes;
  #my $curprime = is_prime($low)  ?  $low  :  next_prime($low);
  #while ($curprime <= $high) {
  #  push @primes, $curprime;
  #  $curprime = next_prime($curprime);
  #}
  #return \@primes;
}

1;

__END__


# ABSTRACT: Utilities related to prime numbers and factoring, using GMP

=pod

=encoding utf8

=head1 NAME

Math::Prime::Util::GMP - Utilities related to prime numbers and factoring, using GMP


=head1 VERSION

Version 0.14


=head1 SYNOPSIS

  use Math::Prime::Util::GMP ':all';
  my $n = "115792089237316195423570985008687907853269984665640564039457584007913129639937";

  # This doesn't impact the operation of the module at all, but does let you
  # enter big number arguments directly as well as enter (e.g.): 2**2048 + 1.
  use bigint;

  # These return 0 for composite, 2 for prime, and 1 for probably prime
  # Numbers under 2^64 will return 0 or 2.
  # is_prob_prime does a BPSW primality test for numbers > 2^64
  # is_prime adds some MR tests and a quick test to try to prove the result
  # is_provable_prime will spend a lot of effort on proving primality

  say "$n is probably prime"    if is_prob_prime($n);
  say "$n is ", qw(composite prob_prime def_prime)[is_prime($n)];
  say "$n is definitely prime"  if is_provable_prime($n) == 2;

  # Miller-Rabin and strong Lucas-Selfridge pseudoprime tests
  say "$n is a prime or spsp-2/7/61" if is_strong_pseudoprime($n, 2, 7, 61);
  say "$n is a prime or slpsp"       if is_strong_lucas_pseudoprime($n);
  say "$n is a prime or eslpsp"      if is_extra_strong_lucas_pseudoprime($n);

  # Return array reference to primes in a range.
  my $aref = primes( 10 ** 200, 10 ** 200 + 10000 );

  $next = next_prime($n);    # next prime > n
  $prev = prev_prime($n);    # previous prime < n

  # Primorials and lcm
  say "23# is ", primorial(23);
  say "The product of the first 47 primes is ", pn_primorial(47);
  say "lcm(1..1000) is ", consecutive_integer_lcm(1000);


  # Find prime factors of big numbers
  @factors = factor(5465610891074107968111136514192945634873647594456118359804135903459867604844945580205745718497);

  # Finer control over factoring.
  # These stop after finding one factor or exceeding their limit.
  #                               # optional arguments o1, o2, ...
  @factors = trial_factor($n);    # test up to o1
  @factors = prho_factor($n);     # no more than o1 rounds
  @factors = pbrent_factor($n);   # no more than o1 rounds
  @factors = holf_factor($n);     # no more than o1 rounds
  @factors = squfof_factor($n);   # no more than o1 rounds
  @factors = pminus1_factor($n);  # o1 = smoothness limit, o2 = stage 2 limit
  @factors = ecm_factor($n);      # o1 = B1, o2 = # of curves
  @factors = qs_factor($n);       # (no arguments)

=head1 DESCRIPTION

A set of utilities related to prime numbers, using GMP.  This includes
primality tests, getting primes in a range, and factoring.

While it certainly can be used directly, the main purpose of this module is
for L<Math::Prime::Util>.  That module will automatically load this if it is
installed, greatly speeding up many of its operations on big numbers.

Inputs and outputs for big numbers are via strings, so you do not need to
use a bigint package in your program.  However if you do use bigints, inputs
will be converted internally so there is no need to convert before a call.
Output results are returned as either Perl scalars (for native-size) or
strings (for bigints).  L<Math::Prime::Util> tries to reconvert
all strings back into the callers bigint type if possible, which makes it more
convenient for calculations.

The various C<is_*_pseudoprime> tests are more appropriately called
C<is_*_probable_prime> or C<is_*_prp>.  They return 1 if the input is a
probable prime based on their test.  The naming convention is historical
and follows Pari and some other math packages.  The modern definition of
pseudoprime is a I<composite> that passes the test, rather than any number.


=head1 FUNCTIONS

=head2 is_prob_prime

  my $prob_prime = is_prob_prime($n);
  # Returns 0 (composite), 2 (prime), or 1 (probably prime)

Takes a positive number as input and returns back either 0 (composite),
2 (definitely prime), or 1 (probably prime).

For inputs below C<2^64> the test is deterministic, so the possible return
values are 0 (composite) or 2 (definitely prime).

For inputs above C<2^64>, a probabilistic test is performed.  Only 0 (composite)
and 1 (probably prime) are returned.  The current implementation uses the
Baillie-PSW (BPSW) test.  There is a possibility that composites may be
returned marked prime, but since the test was published in 1980, not a
single BPSW pseudoprime has been found, so it is extremely likely to be prime.
While we believe (Pomerance 1984) that an infinite number of counterexamples
exist, there is a weak conjecture (Martin) that none exist under 10000 digits.

In more detail, we are using the extra-strong Lucas test (Grantham 2000)
using the Baillie parameter selection method (see OEIS A217719).  Previous
versions of this module used the strong Lucas test with Selfridge parameters,
but the extra-strong version produces fewer pseudoprimes while running
1.2 - 1.5x faster.  It is slightly stronger than the test used in
L<Pari|http://pari.math.u-bordeaux.fr/faq.html#primetest>.


=head2 is_prime

  say "$n is prime!" if is_prime($n);

Takes a positive number as input and returns back either 0 (composite),
2 (definitely prime), or 1 (probably prime).  Composites will act exactly
like C<is_prob_prime>, as will numbers less than C<2^64>.  For numbers
larger than C<2^64>, some additional tests are performed on probable primes
to see if they can be proven by another means.

As with L</is_prob_prime>, a BPSW test is first performed.  If this
indicates "probably prime" then a small number of Miller-Rabin tests
with random bases are performed.  For numbers under 200 bits, a quick
BLS75 C<n-1> primality proof is attempted.  This is tuned to give up
if the result cannot be quickly determined, and results in approximately
30% success rate at 128-bits.

The result is that many numbers will return 2 (definitely prime),
and the numbers that return 1 (probably prime) have gone through more
tests than L</is_prob_prime> while not taking too long.


=head2 is_provable_prime

  say "$n is definitely prime!" if is_provable_prime($n) == 2;

Takes a positive number as input and returns back either 0 (composite),
2 (definitely prime), or 1 (probably prime).  A great deal of effort is
taken to return either 0 or 2 for all numbers.

The current method first uses BPSW and a small number of Miller-Rabin
tests with random bases to weed out composites and provide a deterministic
answer for tiny numbers (under C<2^64>).  A quick BLS75 C<n-1> test is
attempted, followed by ECPP.

The time required for primes of different input sizes on a circa-2009
workstation averages about C<3ms> for 30-digits, C<5ms> for 40-digit,
C<20ms> for 60-digit, C<50ms> for 80-digit, C<100ms> for 100-digit,
C<2s> for 200-digit, and 400-digit inputs about a minute.
Expect a lot of time variation for larger inputs.  You can see progress
indication if verbose is turned on (some at level 1, and a lot at level 2).

A certificate can be obtained along with the result using the
L</is_provable_prime_with_cert> method.


=head2 is_provable_prime_with_cert

Takes a positive number as input and returns back an array with two elements.
The result will be one of:

  (0, '')      The input is composite.

  (1, '')      The input is probably prime but we could not prove it.
               This is a failure in our ability to factor some necessary
               element in a reasonable time, not a significant proof
               failure (in other words, it remains a probable prime).

  (2, '...')   The input is prime, and the certificate contains all the
               information necessary to verify this.

The certificate is a text representation containing all the necessary
information to verify the primality of the input in a reasonable time.
The result can be used with L<Math::Prime::Util/verify_prime> for
verification.  Proof types used include:

  ECPP
  BLS3
  BLS15
  BLS5
  Small

=head2 is_strong_pseudoprime

  my $maybe_prime = is_strong_pseudoprime($n, 2);
  my $probably_prime = is_strong_pseudoprime($n, 2, 3, 5, 7, 11, 13, 17);

Takes a positive number as input and one or more bases.  Returns 1 if the
input is a prime or a strong pseudoprime to all of the bases, and 0 if not.
The base must be a positive integer.  This is often called the Miller-Rabin
test.

If 0 is returned, then the number really is a composite.  If 1 is returned,
then it is either a prime or a strong pseudoprime to all the given bases.
Given enough distinct bases, the chances become very strong that the number
is actually prime.

Both the input number and the bases may be big integers.  If base modulo n
E<lt>= 1 or base modulo n = n-1, then the result will be 1.  This allows the
bases to be larger than n if desired, while still returning meaningful results.
For example,

  is_strong_pseudoprime(367, 1101)

would incorrectly return 0 if this was not done properly.  A 0 result should
be returned only if n is composite, regardless of the base.

This is usually used in combination with other tests to make either stronger
tests (e.g. the strong BPSW test) or deterministic results for numbers less
than some verified limit (e.g. Jaeschke showed in 1993 that no more than three
selected bases are required to give correct primality test results for any
32-bit number).  Given the small chances of passing multiple bases, there
are some math packages that just use multiple MR tests for primality testing,
though in the early 1990s almost all serious software switched to the
BPSW test.

Even numbers other than 2 will always return 0 (composite).  While the
algorithm works with even input, most sources define it only on odd input.
Returning composite for all non-2 even input makes the function match most
other implementations including L<Math::Primality>'s C<is_strong_pseudoprime>
function.


=head2 is_lucas_pseudoprime

=head2 is_strong_lucas_pseudoprime

Takes a positive number as input, and returns 1 if the input is a standard
or strong Lucas probable prime.  The Selfridge method of choosing D, P, and
Q are used (some sources call this a Lucas-Selfridge test).  This is one
half of the BPSW primality test (the Miller-Rabin strong probable prime test
with base 2 being the other half).  The canonical BPSW test (page 1401 of
Baillie and Wagstaff (1980)) uses the strong Lucas test with Selfridge
parameters, but in practice a variety of Lucas tests with different
parameters are used by tests calling themselves BPSW.


=head2 is_extra_strong_lucas_pseudoprime

Takes a positive number as input, and returns 1 if the input is an
extra-strong Lucas probable prime.  This is defined in Grantham (2000),
and is a slightly more stringent test than the strong Lucas test, though
because different parameters are used the pseudoprimes are not a subset.
As expected by the extra conditions, the number of pseudoprimes is less
than 2/3 that of the strong Lucas-Selfridge test.
Runtime performance is 1.2 to 1.5x faster than the strong Lucas test.

The parameters are selected using the Baillie-OEIS method:

  P = 3;
  Q = 1;
  while ( jacobi( P*P-4, n ) != -1 )
    P += 1;


=head2 is_almost_extra_strong_lucas_pseudoprime

Takes a positive number as input and returns 1 if the input is an "almost"
extra-strong Lucas probable prime.  This is the classic extra-strong Lucas
test but without calculating the U sequence.  This makes it very fast,
although as the input increases in size the time converges to the conventional
extra-strong implementation:  at 30 digits this routine is about 15% faster,
at 300 digits it is only 2% faster.

With the current implementations, there is little reason to prefer this unless
trying to reproduce specific results.  The extra-strong implementation has been
optimized to use similar features, removing most of the performance advantage.

An optional second argument (must be between 1 and 256) indicates the
increment amount for P parameter selection.  The default value of one yields
the method described in L</is_extra_strong_lucas_pseudoprime>.  A value of
2 yields the method used in
L<Pari|http://pari.math.u-bordeaux.fr/faq.html#primetest>.

Because the C<U = 0> condition is ignored, this produces about 5% more
pseudoprimes than the extra-strong Lucas test.  However this is still only
66% of the number produced by the strong Lucas-Selfridge test.  No BPSW
counterexamples have been found with any of the Lucas tests described.


=head2 is_frobenius_underwood_pseudoprime

Takes a positive number as input, and returns 1 if the input passes the minimal
lambda+2 test (see Underwood 2012 "Quadratic Compositeness Tests"), where
C<(L+2)^(n-1) = 5 + 2x mod (n, L^2 - Lx + 1)>.  There are no known
counterexamples, but this is not a well studied test.

The computational cost is about 2.5x the cost of a strong pseudoprime test
(this will vary somewhat with platform and input size).  It is typically a
little slower than an extra-strong Lucas test, and faster than a strong
Lucas test.


=head2 is_aks_prime

  say "$n is definitely prime" if is_aks_prime($n);

Takes a positive number as input, and returns 1 if the input passes the
Agrawal-Kayal-Saxena (AKS) primality test.  This is a deterministic
unconditional primality test which runs in polynomial time for general input.
In practice, L</is_nminus1_prime> and L</is_ecpp_prime> are B<much> faster.

Typically you should use L</is_provable_prime> and let it decide the method.

=head2 is_nminus1_prime

  say "$n is definitely prime" if is_nminus1_prime($n);

Takes a positive number as input, and returns 1 if the input passes either
theorem 5 or theorem 7 of the Brillhart-Lehmer-Selfridge primality test.
This is a deterministic unconditional primality test which requires factoring
C<n-1> to a linear factor less than the cube root of the input.  For small
inputs (under 40 digits) this is typically very easy, and some numbers will
naturally lead to this being very fast.  As the input grows, this method
slows down rapidly.

Typically you should use L</is_provable_prime> and let it decide the method.

=head2 is_ecpp_prime

  say "$n is definitely prime" if is_ecpp_prime($n);

Takes a positive number as input, and returns 1 if the input passes the
ECPP primality test.  This is the Atkin-Morain Elliptic Curve Primality
Proving algorithm.  It is the fastest primality proving method in
Math::Prime::Util.

This implementation uses a "factor all strategy" (FAS) with backtracking.
A limited set of about 500 precalculated discriminants are used, which works
well for inputs up to 300 digits, and for many inputs up to one thousand
digits.  Having a larger set will help with large numbers (a set of 2650
is available on github in the C<xt/> directory).  A future implementation
may include code to generate class polynomials as needed.

Typically you should use L</is_provable_prime> and let it decide the method.


=head2 primes

  my $aref1 = primes( 1_000_000 );
  my $aref2 = primes( 2 ** 448, 2 ** 448 + 10000 );
  say join ",", @{primes( 2**2048, 2**2048 + 10000 )};

Returns all the primes between the lower and upper limits (inclusive), with
a lower limit of C<2> if none is given.

An array reference is returned (with large lists this is much faster and uses
less memory than returning an array directly).

The current implementation uses repeated calls to C<next_prime>, which is
good for very small ranges, but not good for large ranges.  A future release
may use a multi-segmented sieve when appropriate.


=head2 next_prime

  $n = next_prime($n);

Returns the next prime greater than the input number.  The function
L</is_prob_prime> is used to determine when a prime is found.


=head2 prev_prime

  $n = prev_prime($n);

Returns the prime smaller than the input number.  0 is returned if the
input is C<2> or lower.  The function L</is_prob_prime> is used to
determine when a prime is found.


=head2 lucas_sequence

  my($U, $V, $Qk) = lucas_sequence($n, $P, $Q, $k)

Computes C<U_k>, C<V_k>, and C<Q_k> for the Lucas sequence defined by
C<P>,C<Q>, modulo C<n>.  The modular Lucas sequence is used in a
number of primality tests and proofs.

The following conditions must hold:
  - C<< D = P*P - 4*Q  !=  0 >>
  - C<< P > 0 >>
  - C<< P < n >>
  - C<< Q < n >>
  - C<< k >= 0 >>
  - C<< n >= 2 >>


=head2 primorial

  $p = primorial($n);

Given an unsigned integer argument, returns the product of the prime numbers
which are less than or equal to C<n>.  This definition of C<n#> follows
L<OEIS series A034386|http://oeis.org/A034386> and
L<Wikipedia: Primorial definition for natural numbers|http://en.wikipedia.org/wiki/Primorial#Definition_for_natural_numbers>.

=head2 pn_primorial

  $p = pn_primorial($n)

Given an unsigned integer argument, returns the product of the first C<n>
prime numbers.  This definition of C<p_n#> follows
L<OEIS series A002110|http://oeis.org/A002110> and
L<Wikipedia: Primorial definition for prime numbers|http://en.wikipedia.org/wiki/Primorial#Definition_for_prime_numbers>.

The two are related with the relationships:

  pn_primorial($n)  ==   primorial( nth_prime($n) )
  primorial($n)     ==   pn_primorial( prime_count($n) )


=head2 consecutive_integer_lcm

  $lcm = consecutive_integer_lcm($n);

Given an unsigned integer argument, returns the least common multiple of all
integers from 1 to C<n>.  This can be done by manipulation of the primes up
to C<n>, resulting in much faster and memory-friendly results than using
n factorial.


=head2 factor

  @factors = factor(640552686568398413516426919223357728279912327120302109778516984973296910867431808451611740398561987580967216226094312377767778241368426651540749005659);
  # Returns an array of 11 factors

Returns a list of prime factors of a positive number, in numerical order.  The
special cases of C<n = 0> and C<n = 1> will return C<n>.

Like most advanced factoring programs, a mix of methods is used.  This
includes trial division for small factors, perfect power detection,
Pollard's Rho, Pollard's P-1 with various smoothness and stage settings,
Hart's OLF (a Fermat variant), ECM (elliptic curve method), and
QS (quadratic sieve).
Certainly improvements could be designed for this algorithm
(suggestions are welcome).

In practice, this factors 26-digit semiprimes in under C<100ms>, 36-digit
semiprimes in under one second.  Arbitrary integers are factored faster.
It is many orders of magnitude faster than any other factoring module on
CPAN circa early 2013.  It is comparable in speed to Math::Pari's C<factorint>
for most inputs.

If you want better factoring in general, I recommend looking at the
standalone programs
L<yafu|http://sourceforge.net/projects/yafu/>,
L<msieve|http://sourceforge.net/projects/msieve/>,
L<gmp-ecm|http://ecm.gforge.inria.fr/>, and
L<GGNFS|http://sourceforge.net/projects/ggnfs/>.


=head2 trial_factor

  my @factors = trial_factor($n);
  my @factors = trial_factor($n, 1000);

Given a positive number input, tries to discover a factor using trial division.
The resulting array will contain either two factors (it succeeded) or the
original number (no factor was found).  In either case, multiplying @factors
yields the original input.  An optional divisor limit may be given as the
second parameter.  Factoring will stop when the input is a prime, one factor
is found, or the input has been tested for divisibility with all primes less
than or equal to the limit.  If no limit is given, then C<2**31-1> will be used.

This is a good and fast initial test, and will be very fast for small numbers
(e.g. under 1 million).  It becomes unreasonably slow in the general case as
the input size increases.


=head2 prho_factor

  my @factors = prho_factor($n);
  my @factors = prho_factor($n, 100_000_000);

Given a positive number input, tries to discover a factor using Pollard's Rho
method.  The resulting array will contain either two factors (it succeeded)
or the original number (no factor was found).  In either case, multiplying
@factors yields the original input.  An optional number of rounds may be
given as the second parameter.  Factoring will stop when the input is a prime,
one factor has been found, or the number of rounds has been exceeded.

This is the Pollard Rho method with C<f = x^2 + 3> and default rounds 64M.  It
is very good at finding small factors.


=head2 pbrent_factor

  my @factors = pbrent_factor($n);
  my @factors = pbrent_factor($n, 100_000_000);

Given a positive number input, tries to discover a factor using Pollard's Rho
method with Brent's algorithm.  The resulting array will contain either two
factors (it succeeded) or the original number (no factor was found).  In
either case, multiplying @factors yields the original input.  An optional
number of rounds may be given as the second parameter.  Factoring will stop
when the input is a prime, one factor has been found, or the number of
rounds has been exceeded.

This is the Pollard Rho method using Brent's modified cycle detection and
backtracking.  It is essentially Algorithm P''2 from Brent (1980).  Parameters
used are C<f = x^2 + 3> and default rounds 64M.  It is very good at finding
small factors.


=head2 pminus1_factor

  my @factors = pminus1_factor($n);

  # Set B1 smoothness to 10M, second stage automatically set.
  my @factors = pminus1_factor($n, 10_000_000);

  # Run p-1 with B1 = 10M, B2 = 100M.
  my @factors = pminus1_factor($n, 10_000_000, 100_000_000);

Given a positive number input, tries to discover a factor using Pollard's
C<p-1> method.  The resulting array will contain either two factors (it
succeeded) or the original number (no factor was found).  In either case,
multiplying @factors yields the original input.  An optional first stage
smoothness factor (B1) may be given as the second parameter.  This will be
the smoothness limit B1 for the first stage, and will use C<10*B1> for
the second stage limit B2.  If a third parameter is given, it will be used
as the second stage limit B2.
Factoring will stop when the input is a prime, one factor has been found, or
the algorithm fails to find a factor with the given smoothness.

This is Pollard's C<p-1> method using a default smoothness of 5M and a
second stage of C<B2 = 10 * B1>.  It can quickly find a factor C<p> of the input
C<n> if the number C<p-1> factors into small primes.  For example
C<n = 22095311209999409685885162322219> has the factor C<p = 3916587618943361>,
where C<p-1 = 2^7 * 5 * 47 * 59 * 3137 * 703499>, so this method will find
a factor in the first stage if C<B1 E<gt>= 703499> or in the second stage if
C<B1 E<gt>= 3137> and C<B2 E<gt>= 703499>.

The implementation is written from scratch using the basic algorithm including
a second stage as described in Montgomery 1987.  It is faster than most simple
implementations I have seen (many of which are written assuming native
precision inputs), but slower than Ben Buhrow's code used in earlier
versions of L<yafu|http://sourceforge.net/projects/yafu/>, and nowhere close
to the speed of the version included with modern GMP-ECM.


=head2 pplus1_factor

  my @factors = pplus1_factor($n);

Given a positive number input, tries to discover a factor using Williams'
C<p+1> method.  The resulting array will contain either two factors (it
succeeded) or the original number (no factor was found).  In either case,
multiplying @factors yields the original input.  An optional first stage
smoothness factor (B1) may be given as the second parameter.  This will be
the smoothness limit B1 for the first stage.
Factoring will stop when the input is a prime, one factor has been found, or
the algorithm fails to find a factor with the given smoothness.



=head2 holf_factor

  my @factors = holf_factor($n);
  my @factors = holf_factor($n, 100_000_000);

Given a positive number input, tries to discover a factor using Hart's OLF
method.  The resulting array will contain either two factors (it succeeded)
or the original number (no factor was found).  In either case, multiplying
@factors yields the original input.  An optional number of rounds may be
given as the second parameter.  Factoring will stop when the input is a
prime, one factor has been found, or the number of rounds has been exceeded.

This is Hart's One Line Factorization method, which is a variant of Fermat's
algorithm.  A premultiplier of 480 is used.  It is very good at factoring
numbers that are close to perfect squares, or small numbers.  Very naive
methods of picking RSA parameters sometimes yield numbers in this form, so
it can be useful to run this a few rounds to check.  For example, the number:

  18548676741817250104151622545580576823736636896432849057 \
  10984160646722888555430591384041316374473729421512365598 \
  29709849969346650897776687202384767704706338162219624578 \
  777915220190863619885201763980069247978050169295918863

was proposed by someone as an RSA key.  It is indeed composed of two distinct
prime numbers of similar bit length.  Most factoring methods will take a
B<very> long time to break this.  However one factor is almost exactly 5x
larger than the other, allowing HOLF to factor this 222-digit semiprime in
only a few milliseconds.


=head2 squfof_factor

  my @factors = squfof_factor($n);
  my @factors = squfof_factor($n, 100_000_000);

Given a positive number input, tries to discover a factor using Shanks'
square forms factorization method (usually known as SQUFOF).  The resulting
array will contain either two factors (it succeeded) or the original number
(no factor was found).  In either case, multiplying @factors yields the
original input.  An optional number of rounds may be given as the second
parameter.  Factoring will stop when the input is a prime, one factor has
been found, or the number of rounds has been exceeded.

This is Daniel Shanks' SQUFOF (square forms factorization) algorithm.  The
particular implementation is a non-racing multiple-multiplier version, based
on code ideas of Ben Buhrow and Jason Papadopoulos as well as many others.
SQUFOF is often the preferred method for small numbers, and L<Math::Prime::Util>
as well as many other packages use it was the default method for native size
(e.g. 32-bit or 64-bit) numbers after trial division.  The GMP version used
in this module will work for larger values, but my testing is showing that it
is not faster than the C<prho> and C<pbrent> methods in general.


=head2 ecm_factor

  my @factors = ecm_factor($n);
  my @factors = ecm_factor($n, 12500);      # B1 = 12500
  my @factors = ecm_factor($n, 12500, 10);  # B1 = 12500, curves = 10

Given a positive number input, tries to discover a factor using ECM.  The
resulting array will contain either two factors (it succeeded) or the original
number (no factor was found).  In either case, multiplying @factors yields the
original input.  An optional maximum smoothness may be given as the second
parameter, which relates to the size of factor to search for.  An optional
third parameter indicates the number of random curves to use at each
smoothness value being searched.

This is an implementation of Hendrik Lenstra's elliptic curve factoring
method, usually referred to as ECM.  The implementation is reasonable,
using projective coordinates, Montgomery's PRAC heuristic for EC
multiplication, and two stages.
It is much slower than the latest GMP-ECM, but still quite useful for
factoring reasonably sized inputs.


=head2 qs_factor

  my @factors = qs_factor($n);

Given a positive number input, tries to discover factors using QS (the
quadratic sieve).  The resulting array will contain one or more numbers such
that multiplying @factors yields the original input.  Typically multiple
factors will be produced, unlike the other C<..._factor> routines.

The current implementation is a modified version of SIMPQS, a predecessor to
the QS in FLINT, and was written by William Hart in 2006.  It will not operate
on input less than 30 digits.  The memory use for large inputs is more than
desired, so other methods such as L</pbrent_factor>, L</pminus1_factor>, and
L</ecm_factor> are recommended to begin with to filter out small factors.
However, it is substantially faster than the other methods on large inputs
having large factors, and is the method of choice for 35+ digit semiprimes.


=head1 SEE ALSO

=over 4

=item L<Math::Prime::Util>.
Has many more functions, lots of fast code for dealing with native-precision
arguments (including much faster primes using sieves), and will use this
module when needed for big numbers.  Using L<Math::Prime::Util> rather than
this module directly is recommended.

=item L<Math::Primality> (version 0.08)
A Perl module with support for the strong Miller-Rabin test, strong
Lucas-Selfridge test, the BPSW probable prime test, next_prime / prev_prime,
the AKS primality test, and prime_count.  It uses L<Math::GMPz> to do all
the calculations, so is faster than pure Perl bignums, but a little slower
than XS+GMP.  The prime_count function is only usable for very small inputs,
but the other functions are quite good for big numbers.  Make sure to use
version 0.05 or newer.

=item L<Math::Pari>
Supports quite a bit of the same functionality (and much more).  See
L<Math::Prime::Util/"SEE ALSO"> for more detailed information on how the
modules compare.

=item L<yafu|http://sourceforge.net/projects/yafu/>,
L<msieve|http://sourceforge.net/projects/msieve/>,
L<gmp-ecm|http://ecm.gforge.inria.fr/>,
L<GGNFS|http://sourceforge.net/projects/ggnfs/>
Good general purpose factoring utilities.  These will be faster than this
module, and B<much> better as the factor increases in size.

=item L<Primo|http://www.ellipsa.eu/public/primo/primo.html> is the state
of the art in freely available (though not open source!) primality proving
programs.  If you have 1000+ digit numbers to prove, you want to use this.

=back


=head1 REFERENCES

=over 4

=item Robert Baillie and Samuel S. Wagstaff, Jr., "Lucas Pseudoprimes", Mathematics of Computation, v35 n152, October 1980, pp 1391-1417.  L<http://mpqs.free.fr/LucasPseudoprimes.pdf>

=item Jon Grantham, "Frobenius Pseudoprimes", Mathematics of Computation, v70 n234, March 2000, pp 873-891.  L<http://www.ams.org/journals/mcom/2001-70-234/S0025-5718-00-01197-2/>

=item John Brillhart, D. H. Lehmer, and J. L. Selfridge, "New Primality Criteria and Factorizations of 2^m +/- 1", Mathematics of Computation, v29, n130, Apr 1975, pp 620-647.  L<http://www.ams.org/journals/mcom/1975-29-130/S0025-5718-1975-0384673-1/S0025-5718-1975-0384673-1.pdf>

=item Richard P. Brent, "An improved Monte Carlo factorization algorithm", BIT 20, 1980, pp. 176-184.  L<http://www.cs.ox.ac.uk/people/richard.brent/pd/rpb051i.pdf>

=item Peter L. Montgomery, "Speeding the Pollard and Elliptic Curve Methods of Factorization", Mathematics of Computation, v48, n177, Jan 1987, pp 243-264.  L<http://www.ams.org/journals/mcom/1987-48-177/S0025-5718-1987-0866113-7/>

=item Richard P. Brent, "Parallel Algorithms for Integer Factorisation", in Number Theory and Cryptography, Cambridge University Press, 1990, pp 26-37.  L<http://www.cs.ox.ac.uk/people/richard.brent/pd/rpb115.pdf>

=item Richard P. Brent, "Some Parallel Algorithms for Integer Factorisation", in Proc. Third Australian Supercomputer Conference, 1999. (Note: there are multiple versions of this paper)  L<http://www.cs.ox.ac.uk/people/richard.brent/pd/rpb193.pdf>

=item William B. Hart, "A One Line Factoring Algorithm", preprint.  L<http://wstein.org/home/wstein/www/home/wbhart/onelinefactor.pdf>

=item Daniel Shanks, "SQUFOF notes", unpublished notes, transcribed by Stephen McMath.  L<http://www.usna.edu/Users/math/wdj/mcmath/shanks_squfof.pdf>

=item Jason E. Gower and Samuel S. Wagstaff, Jr, "Square Form Factorization", Mathematics of Computation, v77, 2008, pages 551-588.  L<http://homes.cerias.purdue.edu/~ssw/squfof.pdf>

=item A.O.L. Atkin and F. Morain, "Elliptic Curves and primality proving", Mathematics of Computation, v61, 1993, pages 29-68.  L<http://www.ams.org/journals/mcom/1993-61-203/S0025-5718-1993-1199989-X/>

=item R.G.E. Pinch, "Some Primality Testing Algorithms", June 1993.  Describes the primality testing methods used by many CAS systems and how most were compromised.  Gives recommendations for primality testing APIs.  L<http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.33.4409>

=back


=head1 AUTHORS

Dana Jacobsen E<lt>dana@acm.orgE<gt>

William Hart wrote the SIMPQS code which is the basis for the QS code.


=head1 ACKNOWLEDGEMENTS

Obviously none of this would be possible without the mathematicians who
created and published their work.  Eratosthenes, Gauss, Euler, Riemann,
Fermat, Lucas, Baillie, Pollard, Brent, Montgomery, Shanks, Hart, Wagstaff,
Dixon, Pomerance, A.K. Lenstra, H. W. Lenstra Jr., Atkin, Knuth, etc.

The GNU GMP team, whose product allows me to concentrate on coding high-level
algorithms and not worry about any of the details of how modular exponentiation
and the like happen, and still get decent performance for my purposes.

Ben Buhrow and Jason Papadopoulos deserve special mention for their open
source factoring tools, which are both readable and fast.  In particular I am
leveraging their SQUFOF work in the current implementation.  They are a huge
resource to the community.

Jonathan Leto and Bob Kuo, who wrote and distributed the L<Math::Primality>
module on CPAN.  Their implementation of BPSW provided the motivation I needed
to do it in this module and L<Math::Prime::Util>.  I also used their
module quite a bit for testing against.

Paul Zimmermann's papers and GMP-ECM code were of great value for my projective
ECM implementation, as well as the papers by Brent and Montgomery.


=head1 COPYRIGHT

Copyright 2011-2013 by Dana Jacobsen E<lt>dana@acm.orgE<gt>

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

SIMPQS Copyright 2006, William Hart.  SIMPQS is distributed under GPL v2+.

=cut
