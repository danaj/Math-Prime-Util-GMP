use ExtUtils::MakeMaker;
use lib 'inc'; # load our bundled version of Devel::CheckLib
use Devel::CheckLib;

# We used to use Devel::CheckLib's "check_lib_or_exit" here.
#
# If both CC and D::CL find the GMP header/library, then everything works,
# but using D::CL didn't help us in any way.
#
# If both CC and D::CL can't find GMP, then we get a slightly nicer message.
#
# If CC can compile GMP but Devel::CheckLib isn't finding either the header
# or library then we have pointlessly obstructed the user.

if (!check_lib(lib => 'gmp', header => 'gmp.h')) {
  warn "\n  It looks like you don't have the GMP library.\n\n";
  warn "  Going to try to compile anyway.  Expect failure.\n";
}

WriteMakefile1(
    NAME         => 'Math::Prime::Util::GMP',
    ABSTRACT     => 'Utilities related to prime numbers, using GMP',
    VERSION_FROM => 'lib/Math/Prime/Util/GMP.pm',
    LICENSE      => 'perl_5',
    AUTHOR       => 'Dana A Jacobsen <dana@acm.org>',

    OBJECT       => 'XS.o '             .
                    'factor.o '         .
                    'gmp_main.o '       .
                    'tinyqs.o '         .
                    'simpqs.o '         .
                    'lanczos.o '        .
                    'simpqs2.o '        .
                    'siqs.o '           .
                    'primality.o '      .
                    'bls75.o '          .
                    'ecpp.o '           .
                    'aks.o '            .
                    'utility.o '        .
                    'misc_ui.o '        .
                    'factmod.o '        .
                    'znlog.o '          .
                    'rootmod.o '        .
                    'legendre_phi.o '   .
                    'prime_iterator.o ' .
                    'poly.o '           .
                    'pbrent63.o '       .
                    'squfof126.o '      .
                    'ecm.o '            .
                    'lucas_seq.o '      .
                    'real.o '           .
                    'isaac.o '          .
                    'random_prime.o '   .
                    'perfect_powers.o ' .
                    'powerfree.o',
    LIBS         => ['-lgmp -lm'],
    clean        => { FILES => 'mpu-siqs mpu-psiqs' },

    TEST_REQUIRES=> {
                      'Math::BigInt'     => '1.88',  # try && bug fixes
                      'Test::More'       => '0.45',
                    },
    PREREQ_PM    => {
                      'Exporter'         => '5.57',
                      'XSLoader'         => '0.01',
                      'Fcntl'            => 0,
                      'Carp'             => 0,
                      'base'             => 0,
                    },

    META_MERGE   => {
               'meta-spec'  => {
                  version => '2',
                  url     => 'http://search.cpan.org/perldoc?CPAN::Meta::Spec',
               },
               dynamic_config => 1,
               resources    => {
                  license    => [ 'http://dev.perl.org/licenses/' ],
                  homepage   => 'https://github.com/danaj/Math-Prime-Util-GMP',
                  repository => {
                    url  => 'https://github.com/danaj/Math-Prime-Util-GMP',
                  },
               },
               prereqs => {
                 runtime => {
                   recommends => {
                     'Math::Prime::Util' => 0.68,
                   },
                 },
               },
               provides => {
                 'Math::Prime::Util::GMP' => {
                    version => '0.55',
                    file => 'lib/Math/Prime/Util/GMP.pm',
                 },
               },
             },

    MIN_PERL_VERSION => 5.006002,
);

sub WriteMakefile1 {   # Cribbed from eumm-upgrade by Alexandr Ciornii
  my %params = @_;
  my $eumm_version = $ExtUtils::MakeMaker::VERSION;
  $eumm_version = eval $eumm_version;

  if ($params{TEST_REQUIRES} and $eumm_version < 6.6303) {
      $params{BUILD_REQUIRES}={ %{$params{BUILD_REQUIRES} || {}} , %{$params{TEST_REQUIRES}} };
      delete $params{TEST_REQUIRES};
  }
  if ($params{BUILD_REQUIRES} and $eumm_version < 6.5503) {
      #EUMM 6.5502 has problems with BUILD_REQUIRES
      $params{PREREQ_PM}={ %{$params{PREREQ_PM} || {}} , %{$params{BUILD_REQUIRES}} };
      delete $params{BUILD_REQUIRES};
  }
  delete $params{CONFIGURE_REQUIRES} if $eumm_version < 6.52;
  delete $params{MIN_PERL_VERSION} if $eumm_version < 6.48;
  delete $params{META_MERGE} if $eumm_version < 6.46;
  delete $params{META_ADD} if $eumm_version < 6.46;
  delete $params{LICENSE} if $eumm_version < 6.31;
  delete $params{AUTHOR} if $] < 5.005;
  delete $params{ABSTRACT_FROM} if $] < 5.005;
  delete $params{BINARY_LOCATION} if $] < 5.005;

  WriteMakefile(%params);
}

sub MY::postamble {
  return <<'MAKE_FRAG';

SIQS_STANDALONE_SOURCES = mpu-siqs.c siqs.c lanczos.c prime_iterator.c \
	squfof126.c pbrent63.c
SIQS_STANDALONE_HEADERS = ptypes.h siqs.h siqs_dep.h lanczos.h \
	prime_iterator.h squfof126.h pbrent63.h utility.h misc_ui.h

SIQS_ARCH_FLAGS = -march=native
SIQS_CFLAGS = -O3 $(SIQS_ARCH_FLAGS)
SIQS_CHECK_ARGS =

.PHONY: siqs psiqs check-siqs check-psiqs test-siqs
siqs: mpu-siqs
psiqs: mpu-psiqs

# Optional standalone health checks; not dependencies of Perl's test target.
check-siqs:
	$(PERLRUN) tools/siqs-check.pl --cc "$(CC)" --cflags "$(SIQS_CFLAGS)" $(SIQS_CHECK_ARGS)

check-psiqs:
	$(PERLRUN) tools/siqs-check.pl --cc "$(CC)" --cflags "$(SIQS_CFLAGS)" --threaded $(SIQS_CHECK_ARGS)

test-siqs: check-siqs

mpu-siqs: $(SIQS_STANDALONE_SOURCES) $(SIQS_STANDALONE_HEADERS)
	$(CC) $(SIQS_CFLAGS) -DSTANDALONE -o $@ $(SIQS_STANDALONE_SOURCES) -lgmp -lm

mpu-psiqs: $(SIQS_STANDALONE_SOURCES) $(SIQS_STANDALONE_HEADERS) psiqs_inc.c planczos_inc.c
	$(CC) $(SIQS_CFLAGS) -DSTANDALONE -DPSIQS -pthread -o $@ $(SIQS_STANDALONE_SOURCES) -lgmp -lm

# Build only: an optional synthetic matrix scaling tool, not a Perl test.
lanczos-bench: tools/lanczos-bench.c lanczos.c lanczos.h planczos_inc.c siqs_dep.h ptypes.h
	$(CC) $(SIQS_CFLAGS) -DSTANDALONE -DPSIQS -pthread -o $@ tools/lanczos-bench.c -lm

# Build only: a fixed-work byte-sieve block-size/concurrency screen.
siqs-sieve-bench: tools/siqs-sieve-bench.c $(SIQS_STANDALONE_SOURCES) $(SIQS_STANDALONE_HEADERS) psiqs_inc.c planczos_inc.c
	$(CC) $(SIQS_CFLAGS) -DSTANDALONE -DPSIQS -pthread -o $@ tools/siqs-sieve-bench.c lanczos.c prime_iterator.c squfof126.c pbrent63.c -lgmp -lm

MAKE_FRAG
}
