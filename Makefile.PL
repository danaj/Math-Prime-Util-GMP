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
                    'pr32_iter.o '      .
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
    clean        => { FILES => 'msiqs' },

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

# The module uses Perl's compiler/flags and root-level objects. Standalone
# programs and development targets have their own Makefile under siqs/.
siqs$(OBJ_EXT): siqs/siqs.c siqs/siqs.h siqs/siqs_dep.h siqs/lanczos.h siqs/siqs_choose_inc.c siqs/siqs_sieve_kernel_inc.c siqs/psiqs_inc.c siqs/pr32_iter.h ptypes.h mont64.h squfof126.h misc_ui.h
	$(CCCMD) $(CCCDLFLAGS) "-I$(PERL_INC)" $(PASTHRU_DEFINE) $(DEFINE) siqs/siqs.c

lanczos$(OBJ_EXT): siqs/lanczos.c siqs/lanczos.h siqs/planczos_inc.c ptypes.h
	$(CCCMD) $(CCCDLFLAGS) "-I$(PERL_INC)" $(PASTHRU_DEFINE) $(DEFINE) siqs/lanczos.c

pr32_iter$(OBJ_EXT): siqs/pr32_iter.c siqs/pr32_iter.h
	$(CCCMD) $(CCCDLFLAGS) "-I$(PERL_INC)" $(PASTHRU_DEFINE) $(DEFINE) siqs/pr32_iter.c

XS$(OBJ_EXT) factor$(OBJ_EXT) bls75$(OBJ_EXT): siqs/siqs.h
simpqs2$(OBJ_EXT): siqs/lanczos.h

# Convenience only: standalone flags and all other targets belong to siqs/.
.PHONY: siqs
siqs:
	$(MAKE) -C siqs siqs
	$(CP) siqs/msiqs msiqs

MAKE_FRAG
}
