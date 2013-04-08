use ExtUtils::MakeMaker;
use lib 'inc'; # load our bundled version of Devel::CheckLib
use Devel::CheckLib;

check_lib_or_exit(lib => 'gmp', header => 'gmp.h');

WriteMakefile1(
    NAME         => 'Math::Prime::Util::GMP',
    ABSTRACT     => 'Utilities related to prime numbers, using GMP',
    VERSION_FROM => 'lib/Math/Prime/Util/GMP.pm',
    LICENSE      => 'perl',
    AUTHOR       => 'Dana A Jacobsen <dana@acm.org>',

    OBJECT       => 'prime_iterator.o ' .
                    'small_factor.o '   .
                    'utility.o '        .
                    'ecm.o '            .
                    'simpqs.o '         .
                    'gmp_main.o '       .
                    'XS.o',
    LIBS         => ['-lgmp -lm'],

    BUILD_REQUIRES=>{
                      'Test::More'       => '0.45',
                    },
    PREREQ_PM    => {
                      'Exporter'         => '5.562',
                      'XSLoader'         => '0.01',
                      'Carp'             => 0,
                      'base'             => 0,
                    },

    META_MERGE   => {
               resources => {
                  homepage   => 'https://github.com/danaj/Math-Prime-Util-GMP',
                  repository => 'https://github.com/danaj/Math-Prime-Util-GMP',
               },
               recommends => {
                  'Math::Prime::Util' => 0.024,
               },
             },

    MIN_PERL_VERSION => 5.006002,
);

sub WriteMakefile1 {   # Cribbed from eumm-upgrade by Alexandr Ciornii
  my %params = @_;
  my $eumm_version = $ExtUtils::MakeMaker::VERSION;
  $eumm_version = eval $eumm_version;

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
