use ExtUtils::MakeMaker;
use lib 'inc'; # load our bundled version of Devel::CheckLib
use Devel::CheckLib;

check_lib_or_exit(lib => 'gmp', header => 'gmp.h');

WriteMakefile1(
    NAME         => 'Math::Prime::Util::GMP',
    ABSTRACT     => 'Utilities related to prime numbers, using GMP',
    VERSION_FROM => 'lib/Math/Prime/Util/GMP.pm',
    LICENSE      => 'perl_5',
    AUTHOR       => 'Dana A Jacobsen <dana@acm.org>',

    OBJECT       => 'prime_iterator.o ' .
                    'utility.o '        .
                    'primality.o '      .
                    'lucas_seq.o '      .
                    'rootmod.o '        .
                    'factor.o '         .
                    'pbrent63.o '       .
                    'squfof126.o '      .
                    'ecm.o '            .
                    'tinyqs.o '         .
                    'simpqs.o '         .
                    'bls75.o '          .
                    'ecpp.o '           .
                    'aks.o '            .
                    'gmp_main.o '       .
                    'real.o '           .
                    'isaac.o '          .
                    'random_prime.o '   .
                    'XS.o',
    LIBS         => ['-lgmp -lm'],

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
                    version => '0.53',
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
