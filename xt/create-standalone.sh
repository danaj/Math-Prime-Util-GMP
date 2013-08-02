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
CFLAGS = -O3 -fomit-frame-pointer -Wall $(DEFINES)
LIBS = -lgmp -lm

OBJ = ecpp.o bls75.o ecm.o prime_iterator.o gmp_main.o small_factor.o utility.o
HEADERS = ptypes.h class_poly_data.h

.PHONY: default all clean

default: $(TARGET)
all: default

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJ)
	$(CC) $(OBJ) $(LIBS) -o $@

clean:
	-rm -f *.o

realclean distclean: clean
	-rm -f $(TARGET)
EOM

cat << 'EOREADME' > standalone/README

ECPP-DJ:  Elliptic Curve Primality Proof.

Dana Jacobsen (dana@acm.org), 2012-2013

Let me know if you find this software useful, and suggestions, comments, and
patches are welcome.

This is a standalone version of the ECPP implemention written for the Perl
module Math::Prime::Util::GMP in 2013.  This uses a "Factor All" strategy, and
closely follows the papers by Atkin and Morain.  Most of the utility functions
closely follow the algorithms presented in Cohen's book "A Course in
Computational Algebraic Number Theory".  Almost all the factoring is done
with my p-1 implementation.  The ECM factoring and manipulation was heavily
insipired by GMP-ECM by Phil Zimmerman and many others.

This includes a strong BPSW test (e.g. strong PRP-2 test followed by strong
Lucas-Selfridge test).  We use this to (1) detect composites, and (2) prove
small numbers.  BPSW has no counterexamples below 2^64, and this implementation
has been verified against Feitsma's database.

Using the -c option enables certificate generation.  The format is described
in the file proof-text-format.txt, and a verification program for both this
format and PRIMO is supplied in verify-cert.pl, though it needs a very new
version of the Math::Prime::Util module.

Performance is quite good for most sub-1000 digit numbers, but starts getting
uneven over 500 digits.  Much more work is possible here.

In my testing, it is much, much faster than GMP-ECPP 2.46.  It is fairly
similar in speed under ~300 digits to Morain's old 6.4.5a ECPP, but is
substantially faster for larger numbers.  It is faster than WraithX's APRCL
to about 1000 digits, then starts getting slower.  Note that APRCL does not
produce a certificate.  AKS is currently not a viable proof method, being
millions of times slower than these methods.

Some areas to concentrate on:

 1. The polynomials.  We ought to generate Weber polynomials on the fly.  I
    this it still makes a lot of sense to include a fixed set (e.g. all polys
    of degree 6 or smaller) for speed.  However the lack of polynomials is a
    big issue with titanic numbers, as we run a good chance of not finding an
    easily splitable 'm' value and then get bogged down in factoring.

 2. The factoring.  In most cases this will stay in factoring stage 1 the
    entire time, meaning we are running my _GMP_pminus1_factor code with small
    parameters.  Optimizing this would help (the stage 2 code certainly could
    be faster).  I have tried GMP-ECM's n-1 and it is quite a bit slower for
    this purpose (let me be clear: for large B1/B2 values, GMP-ECM rocks, but
    in this application with small B1/B2, it ran slower for me).  If you add
    the define USE_LIBECM, then GMP-ECM's ECM is used.  This may or may not
    be faster, and needs tuning.

    Where using GMP-ECM would really help (I think) is in later factoring
    stages where we're in trouble and need to work hard to find a factor since
    there just aren't any easy ones to be found.  At this point we want to
    unleash the full power of GMP-ECM.  I have not tested this in this
    application, but for general factoring, GMP-ECM is much faster than my
    ECM code so I would expect something similar here.

    Note that this interacts with #1, as if we can efficiently generate polys
    or have a giant set, then we must trade off doing more factoring vs. more
    polys.

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
    searching polys vs. backtracking.  It may be worthwhile for largish
    numbers (e.g. 300+ digits) to find all the 'q' values at this stage, then
    select the one with the smallest 'q' (this is an idea from Morain, I
    believe).  The result would be more time spent searching at a given level,
    but we'd get a shallower tree.  This is not hard with a structure like my
    FPS code uses, but would take some jiggering in the simple FAS loop (since
    we'd have to be prepared for backtracking).

 5. The poly root finding takes a long time for large degree polys, and perhaps
    we could speed it up.  There is a little speedup applied, where we exit
    early after finding 4 roots, since we really only need 1.  If we could get
    polyz_pow_polymod to run faster this would help here in general.  For
    titanic numbers this becomes a big bottleneck.


Note 1: AKS is also included and you can use it via the '-aks' option.
        This runs about the same speed indicated by Brent 2010, which means
        absurdly slow.  Let me know if you find anything faster (the
        Berstein-Lenstra r selection would help).

Note 2: You can also force use of the BLS75 theorem 5/7 n-1 proof using the
        '-nm1' option.  This is good for up to 70-80 digits or so.  It performs
        similarly to Pari's n-1 implementation (though is presumably very
        different internally).

EOREADME
