#!/bin/bash

if [ ! -d standalone ]
then
  mkdir standalone
fi

cp -p ptypes.h standalone/
cp -p ecpp.[ch] bls75.[ch] ecm.[ch] prime_iterator.[ch] standalone/
cp -p gmp_main.[ch] small_factor.[ch] utility.[ch] standalone/
cp -p xt/proof-text-format.txt standalone/
cp -p examples/verify-cert.pl standalone/
cp -p examples/vcert.c standalone/

if [ -f xt/class_poly_data_big.h ]
then
  echo Using large poly set
  cp -p xt/class_poly_data_big.h standalone/class_poly_data.h
else
  echo Using small poly set
  cp -p class_poly_data.h standalone/
fi

# Standalone ECPP doesn't need SIMPQS, so let's not include it.
cat << 'EOSIMPQSH' > standalone/simpqs.h
#ifndef MPU_SIMPQS_H
#define MPU_SIMPQS_H
#include <gmp.h>
static int _GMP_simpqs(mpz_t n, mpz_t* farray) { return 0; }
#endif
EOSIMPQSH

# gcc -O3 -fomit-frame-pointer -DSTANDALONE -DSTANDALONE_ECPP ecpp.c bls75.c ecm.c prime_iterator.c gmp_main.c small_factor.c utility.c -o ecpp-dj -lgmp -lm

cat << 'EOM' > standalone/Makefile
TARGET = ecpp-dj
CC = gcc
DEFINES = -DSTANDALONE -DSTANDALONE_ECPP
CFLAGS = -O3 -g -Wall $(DEFINES)
LIBS = -lgmp -lm

OBJ = ecpp.o bls75.o ecm.o prime_iterator.o gmp_main.o small_factor.o utility.o
HEADERS = ptypes.h class_poly_data.h

.PHONY: default all clean

default: $(TARGET) vcert
all: default

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJ)
	$(CC) $^ $(LIBS) -o $@

vcert: vcert.o
	$(CC) $^ $(LIBS) -o $@

clean:
	-rm -f *.o

realclean distclean: clean
	-rm -f $(TARGET) vcert

TEST1 = 11739771271677308623
TEST2 = 4101186565771483058393796013015990306873
TEST3 = 35393208195490503043982975376709535956301978008095654379254570182486977128836509
TEST4 = 935006864223353238737221824562525122013245238547726726389470209512876196672409305516259055214126542251267716102384835899
TEST5 = 25433248801892216938808041984973587415848723989586550342338936524369752729114437550630504573611493987544033285909689162724975971196674890122730486522513928304389686931651745712950504486324914805161849
TEST6 = 9805395078382368090505040572030888128042903590778595659611793383233089744266057832338719254777179628060930905633717264804767081354556273122188172463117096876389627469817121375789418174467686118164854924603273433083668455457307317475867158494243611475656192285904836693990000951833
TEST7 = 14783869341310731862535181711924146413209106702463668486125206627959998111552028525869806779180480078890581625378293879184823213870417286554264959032457904699577075117266041686053691857746703656877556879868345571593862758480658145938864170620882703502631810038072518245192559323694333232320523172244347200606659255752705575709361904749964234503630803
TEST8 = 516069211402263086362802849056712149142247708608707602810064116323740967320507238792597022510888972176570219148745538673939283056074749627707686516295450338874349093022059487090007604327705115534233283204877186839852673049019899452650976370088509524263064316558360307931746214148027462933967332118502005274436360122078931230288854613802443446620402918656243749843390157467993315073584216791188709452340255395988925217870997387615378161410168698419453325449311

test: $(TARGET) vcert
	@for n in $(TEST1) $(TEST2) $(TEST3) $(TEST4) $(TEST5) $(TEST6) $(TEST7) $(TEST8) $(TEST9); do echo -n "$${#n} digits ..."; ./ecpp-dj -c $$n | ./vcert -q -; if [ $$? -eq 0 ]; then echo "Pass"; else echo "FAIL"; fi; done
EOM

cat << 'EOREADME' > standalone/README

ECPP-DJ:  Elliptic Curve Primality Proof.

Copyright 2012-2014, Dana Jacobsen (dana@acm.org).

Let me know if you find this software useful, and suggestions, comments, and
patches are welcome.

INSTALLATION:
     # See note 3 at the end of this file if you want to enable APRCL.
     make
     make test           (optional)
     ./ecpp-dj -help     (shows usage)

     # If you plan on doing proofs with numbers over 800 digits, consider:
     #   wget http://probableprime.org/ecpp/cpd/huge/class_poly_data.h.gz
     #   gunzip class_poly_data.h.gz

This is a standalone version of the ECPP implemention written for the Perl
module Math::Prime::Util::GMP in 2013.  This uses a "Factor All" strategy, and
closely follows the papers by Atkin and Morain.  Most of the utility functions
closely follow the algorithms presented in Cohen's book "A Course in
Computational Algebraic Number Theory".  Almost all the factoring is done
with my p-1 implementation.  The ECM factoring and manipulation was heavily
insipired by GMP-ECM by Phil Zimmerman and many others.

This includes a BPSW test (strong PRP-2 test followed by extra strong
Lucas-Selfridge test).  We use this to (1) detect composites, and (2) prove
small numbers.  BPSW has no counterexamples below 2^64, and this implementation
has been verified against Feitsma's database.

Using the -c option enables certificate generation.  The format is described
in the file proof-text-format.txt.  Two verification programs are supplied:
  verify-cert.pl  (Perl)
  vert.c          (C+GMP)
The Makefile should create the 'vcert' executable for you.  Both programs
will read both the MPU format and the PRIMO format.

Performance is quite good for most sub-1000 digit numbers, but starts getting
uneven over 500 digits.  Much more work is possible here.

For production proving of multi-thousand digit numbers, I recommend:
   Primo   http://www.ellipsa.eu/public/primo/primo.html
because of its large speed advantage for 1000+ digit numbers, especially on
multi-processor machines.

Quick performance comparisons:
 - Primo is slower under ~300 digits, *much* faster as input grows.
 - GMP-ECPP 2.49 is very, very slow.  Nearly unusable once over 500 digits.
 - Morain's ancient 6.4.5a ECPP is similar under 300 digits, but gets slower.
 - David Cleaver's mpz_aprcl is a tiny bit slower under ~700 digits, but gets
   faster for larger inputs (2-3x faster at 2000 digits).  Note that APR-CL
   does not produce a certificate.
 - AKS is not currently a viable method, with all known implementations being
   millions of times slower than alternative methods (ECPP or APR-CL).

Some areas to concentrate on:

 1. The polynomials.  We ought to generate Weber polynomials on the fly.  I
    think it still makes a lot of sense to include a fixed set (e.g. all polys
    of degree 6 or smaller) for speed.  However the lack of polynomials is a
    big issue with titanic numbers, as we run a good chance of not finding an
    easily splitable 'm' value and then get bogged down in factoring.  The CM
    package from http://cm.multiprecision.org/ would be an excellent choice,
    with my only concern being the large dependency chain.

 2. The factoring.  In most cases this will stay in factoring stage 1 the
    entire time, meaning we are running my _GMP_pminus1_factor code with small
    parameters.  Optimizing this would help (the stage 2 code certainly could
    be faster).  I have tried GMP-ECM's n-1 and it is quite a bit slower for
    this purpose (let me be clear: for large B1/B2 values, GMP-ECM rocks, but
    in this application with small B1/B2, it ran slower for me).  If you add
    the define USE_LIBECM, then GMP-ECM's ECM is used.  This will probably be
    slower.

    Where using GMP-ECM would really help (I think) is in later factoring
    stages where we're in trouble and need to work hard to find a factor since
    there just aren't any easy ones to be found.  At this point we want to
    unleash the full power of GMP-ECM.  I have not tested this in this
    application, but for general factoring, GMP-ECM is much faster than my
    ECM code so I would expect something similar here.

    Note that this interacts with #1, as if we can efficiently generate polys
    or have a giant set, then we must trade off doing more factoring vs. more
    polys.  Morain's fastECPP, for example, uses stupendously more
    discriminants than we do, so factoring isn't a big issue for them.  They
    have very different polynomials and root finding algorithms however.

 3. Parallelism.  Currently the entire code is single threaded.  There are
    many opportunities here.
    - Something I think would be useful and not too much work is parallelizing
    the dlist loop in ecpp_down, so all the discriminants can be searched at
    once.  A more complicated solution would be a work queue including pruning
    so we could recurse down many trees at once.
    - If we have to run ECM, then clearly we can run multiple curves at once.
    - Finally, once we've hit STEP 2 of ecpp_down, we could do curve finding
    for the entire chain in parallel.  This would be especially useful for
    titanic numbers where this is a large portion of the total time.

 4. ecpp_down.  There are a lot of little things here that can have big
    impacts on performance.  For instance the decisions on when to keep
    searching polys vs. backtracking.

    The current code, for smallish numbers, will use a cheap factoring stage
    zero for a while before switching to stage 1.  There is a lot of repeated
    work here that a rewrite could avoid.

 5. The poly root finding takes a long time for large degree polys, and perhaps
    we could speed it up.


Note 1: AKS is also included and you can use it via the '-aks' option.  The
        default implentation includes improvements from Bernstein and Voloch,
        with a tuning heuristic from Bornemann.  It is much faster than the
        version used by Brent (2006) for instance, but it is still very slow.

Note 2: You can also force use of the BLS75 theorem 5/7 n-1 proof using the
        '-nm1' option.  This is good for up to 70-80 digits or so.  It performs
        similarly to Pari's n-1 implementation (though is presumably very
        different internally).

Note 3: You can also run APRCL.  Get WraithX's code from:
               http://sourceforge.net/projects/mpzaprcl/
        , put the files in this directory, then add -DUSE_APRCL to the
        DEFINES in the Makefile.  The '-aprcl' option will now be available.

EOREADME
