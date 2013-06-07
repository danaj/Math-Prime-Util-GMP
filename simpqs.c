/*============================================================================

   Quadratic Sieve

   This is derived from SIMPQS, copyright 2006 William Hart.

   Modifications made in 2013 by Dana Jacobsen:
     - returns all coprime factors found
     - put it in one file
     - merge some of the 2.0 changes
     - make it work with smaller values
     - fix some memory errors
     - free memory all over
     - fewer globals
     - Use prime_iterator -- much faster than mpz_nextprime
     - Alternate multiplier selection routine.
     - lots of little changes / optimizations

   Version 2.0 scatters temp files everywhere, but that could be solved.
   The main benefits left in 2.0 are:
      (1) combining partial relations (this is huge for large inputs)
      (2) much less memory use, though partly due to using temp files
      (3) jasonp's block Lanczos routine.
   This code goes through curves slightly faster than v2.0, but with big
   inputs it ends up needing 2x the time because of not combining partials
   as well as the final linear algebra time.
============================================================================*/

/*============================================================================

    SIMPQS is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    SIMPQS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SIMPQS; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

============================================================================*/


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <limits.h>
#include <math.h>
#include <gmp.h>

#ifdef STANDALONE_SIMPQS
  typedef unsigned long UV;
  typedef   signed long IV;
  #define INLINE
  #define UV_MAX ULONG_MAX
  #define UVCONST(x) ((unsigned long)x##UL)
  #define croak(fmt,...)            { printf(fmt,##__VA_ARGS__); exit(1); }
  #define New(id, mem, size, type)  mem = (type*) malloc((size)*sizeof(type))
  #define Newz(id, mem, size, type) mem = (type*) calloc(size, sizeof(type))
  #define Safefree(mem)             free((void*)mem)
  #define PRIME_ITERATOR(i) mpz_t i; mpz_init_set_ui(i, 2)
  static UV prime_iterator_next(mpz_t *iter) { mpz_nextprime(*iter, *iter); return mpz_get_ui(*iter); }
  static void prime_iterator_destroy(mpz_t *iter) { mpz_clear(*iter); }
  static void prime_iterator_setprime(mpz_t *iter, UV n) {mpz_set_ui(*iter, n);}
  /* static int prime_iterator_isprime(mpz_t *iter, UV n) {int isp; mpz_t t; mpz_init_set_ui(t, n); isp = mpz_probab_prime_p(t, 10); mpz_clear(t); return isp;} */
  static int _verbose = 0;
  static int get_verbose_level(void) { return _verbose; }
#else
  #include "ptypes.h"
  #include "simpqs.h"
  #include "prime_iterator.h"
#endif

#include "utility.h"

/* DANAJ: Modify matrix code to do 64-bit-padded character arrays */
typedef unsigned char* row_t;  /* row of an F2 matrix */
typedef row_t* matrix_t;       /* matrix as a list of pointers to rows */

#define insertEntry(m, i, j)   m[i][(j)/8] |= (1U << ((j)%8))
#define xorEntry(m, i, j)      m[i][(j)/8] ^= (1U << ((j)%8))
#define getEntry(m, i, j)     (m[i][(j)/8] &  (1U << ((j)%8)))
#define swapRows(m, x, y) \
  do { row_t temp = m[x];   m[x] = m[y];  m[y] = temp; } while (0)

#define matBytes(numcols) (((numcols+63)/64) * 8)
#define rightMatrixOffset(numcols)  (8 * matBytes(numcols))

/* Clear just the left side */
static INLINE void clearRow(matrix_t m, unsigned int numcols, unsigned int row)
{
  memset( m[row], 0, matBytes(numcols) );
}

/* bitwise xor of two rows, both left and right matrices */
static void xorRows(matrix_t m, unsigned int numcols, unsigned int source, unsigned int dest)
{
  unsigned int i, q;
  UV* x = (UV*) m[dest];
  UV* y = (UV*) m[source];
  size_t nwords = (2 * matBytes(numcols)) / sizeof(UV);

  q = 8 * (nwords / 8);
  for (i = 0; i < q; i += 8) {
    x[i+0] ^= y[i+0];  x[i+1] ^= y[i+1];  x[i+2] ^= y[i+2];  x[i+3] ^= y[i+3];
    x[i+4] ^= y[i+4];  x[i+5] ^= y[i+5];  x[i+6] ^= y[i+6];  x[i+7] ^= y[i+7];
  }
  for ( ; i < nwords; i++)
    x[i] ^= y[i];
}

static matrix_t constructMat(unsigned int cols, unsigned int rows)
{
  unsigned int i;
  matrix_t m;
  size_t nbytes = matBytes(cols);
  unsigned int mat2offset = rightMatrixOffset(cols);

  /* printf("construct mat %u %u (%lu bytes)\n", cols, rows, rows*sizeof(row) + rows*(2*nbytes)); */
  /* If cols > rows, we write off the array */
  if (cols < rows) croak("SIMPQS:  cols %u > rows %u\n", cols, rows);
  New(0, m, rows, row_t);
  if (m == 0) croak("SIMPQS: Unable to allocate memory for matrix!\n");

  for (i = 0; i < rows; i++) { /* two matrices, side by side */
    Newz(0, m[i], 2*nbytes, unsigned char);
    if (m[i] == 0) croak("SIMPQS: Unable to allocate memory for matrix!\n");
  }

  /* make second matrix identity, i.e. 1's along diagonal */
  for (i = 0; i < rows; i++)
    insertEntry(m, i, mat2offset + i);

  return m;
}

static void destroyMat(matrix_t m, unsigned int rows)
{
  unsigned int i;
  for (i = 0; i < rows; i++)
    Safefree(m[i]);
  Safefree(m);
}

#if 0
static void displayRow(matrix_t m, unsigned int row, unsigned int numcols)
{
  int j;
  unsigned int mat2offset = rightMatrixOffset(numcols);

  printf("[");
  for (j = 0; j < numcols; j++)
    printf("%c", getEntry(m,row,j) ? '1' : '0');
  printf("  ");
  for (j = 0; j < numcols; j++)
    printf("%c", getEntry(m,row,mat2offset+j) ? '1' : '0');
  printf("]\n");
}
#endif

/* gaussReduce:  Apply Gaussian elimination to a matrix. */
static unsigned int gaussReduce(matrix_t m, unsigned int cols, unsigned int rows)
{
  unsigned int rowUpto = 0;
  unsigned int irow, checkRow;
  int icol;

  for (icol = cols-1; icol >= 0; icol--) {
    irow = rowUpto;

    while ( (irow < rows) && (getEntry(m,irow,icol) == 0) )
      irow++;

    if (irow < rows) {
      swapRows(m,rowUpto,irow);
      for (checkRow = rowUpto+1; checkRow < rows; checkRow++) {
        if (getEntry(m,checkRow,icol) != 0)
          xorRows(m, cols, rowUpto, checkRow);
      }
      rowUpto++;
    }
  }
  return rowUpto;
}

/*===========================================================================*/
 /* Uncomment these for various pieces of debugging information */

 /* Shows the number of relations generated and curves used during sieving */
/* #define COUNT */
 /* Shows the actual factorizations of the relations */
/* #define RELPRINT */
 /* Error if relation should be divisible by a prime but isn't */
/* #define ERRORS */
 /* Shows the polynomials being used by the sieve */
/* #define POLS */
 /* Prints some details about the factors of the A coefficients of the polys */
/* #define ADETAILS */
 /* Prints the size of the largest factorbase prime */
/* #define LARGESTP */
 /* Prints the number of curves used and number of partial relations */
/* #define CURPARTS */
 /* Report sieve size, multiplier and number of primes used */
/* #define REPORT */

#ifdef ERRORS
  #define CHECK_EXPONENT(exponent,k) \
     if (exponent==0) printf("Error with prime %u!\n", factorBase[k]);
#else
  #define CHECK_EXPONENT(exponent,k)
#endif
#ifdef RELPRINT
  #define PRINT_FB(exponent, k) \
     do { if (exponent > 0) printf(" %u", factorBase[k]); \
          if (exponent > 1) printf("^%d", exponent); } while (0)
#else
  #define PRINT_FB(exponent, k)
#endif

/*===========================================================================*/
/* Architecture dependent fudge factors */

#if ULONG_MAX == 4294967295UL
#define SIEVEMASK 0xC0C0C0C0UL
#define SIEVEDIV 1
#elif ULONG_MAX == 18446744073709551615UL
#define SIEVEMASK 0xC0C0C0C0C0C0C0C0UL
#define SIEVEDIV 1
#else
 #error Cannot determine ulong size
#endif

/* Should be a little less than the L1/L2 cache size and a multiple of 64000 */
#define CACHEBLOCKSIZE 64000
/* Make lower for slower machines */
#define SECONDPRIME    6000
/* Used for tweaking the bit size calculation for FB primes */
#define SIZE_FUDGE     0.15

/* Will not factor numbers with less than this number of decimal digits */
#define MINDIG 30

/*===========================================================================*/
/*  Large prime cutoffs, in thousands */
static const unsigned int largeprimes[] =
{
   100,  100,  125,   125,   150,   150,   175,   175,   200,   200, /* 30-39 */
   250,  300,  370,   440,   510,   580,   650,   720,   790,  8600, /* 40-49 */
   930, 1000, 1700,  2400,  3100,  3800,  4500,  5200,  5900,  6600, /* 50-59 */
  7300, 8000, 8900, 10000, 11300, 12800, 14500, 16300, 18100, 20000, /* 60-69 */
   22000,  24000,  27000,  32000,  39000, /* 70-74 */
   53000,  65000,  75000,  87000, 100000, /* 75-79 */
  114000, 130000, 150000, 172000, 195000, /* 80-84 */
  220000, 250000, 300000, 350000, 400000, /* 85-89 */
  450000, 500000 /* 90-91 */
};

/*===========================================================================*/
/* Number of primes to use in factor base */
static const unsigned int primesNo[] =
{
     1500, 1600, 1600, 1600, 1600, 1600, 1600, 1600, 1600, 1600, /* 30-39 */
     1600, 1600, 1600, 1700, 1750, 1800, 1900, 2000, 2050, 2100, /* 40-49 */
     2150, 2200, 2250, 2300, 2400, 2500, 2600, 2700, 2800, 2900, /* 50-59 */
     3000, 3150, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, /* 60-69 */
      9500, 10000, 11500, 13000, 15000, /* 70-74 */
     17000, 24000, 27000, 30000, 37000, /* 75-79 */
     45000, 47000, 53000, 57000, 58000, /* 80-84 */
     59000, 60000, 64000, 68000, 72000, /* 85-89 */
     76000, 80000 /* 90-91 */
};

/*===========================================================================*/
/* First prime actually sieved for */
static const unsigned int firstPrimes[] =
{
      3,  3,  3,  3,  3,  3,  3,  3,  3,  3, /* 30-39 */
      3,  3,  3,  4,  6,  6,  7,  8,  9, 10, /* 40-49 */
     11, 11, 11, 11, 11, 12, 12, 12, 12, 12, /* 50-59 */
     14, 14, 14, 14, 14, 14, 14, 14, 15, 17, /* 60-69 */
     19, 21, 22, 22, 23, /* 70-74 */
     24, 25, 25, 26, 26, /* 75-79 */
     27, 27, 27, 27, 28, /* 80-84 */
     28, 28, 28, 29, 29, /* 85-89 */
     29, 29 /* 90-91 */
};

/*===========================================================================*/
/* Logs of primes are rounded and errors accumulate
 * This specifies how great an error to allow */
static const unsigned int errorAmounts[] =
{
     10, 10, 10, 11, 13, 14, 14, 15, 15, 16, /* 30-39 */
     16, 17, 17, 18, 18, 19, 19, 19, 20, 20, /* 40-49 */
     21, 21, 21, 22, 22, 22, 23, 23, 23, 24, /* 50-59 */
     24, 24, 25, 25, 25, 25, 26, 26, 26, 26, /* 60-69 */
     27, 27, 28, 28, 29, /* 70-74 */
     29, 30, 30, 30, 31, /* 75-79 */
     31, 31, 31, 32, 32, /* 80-84 */
     32, 32, 32, 33, 33, /* 85-89 */
     33, 33 /* 90-91 */
};

/*===========================================================================*/
/* Threshold the sieve value must exceed to be considered for smoothness */
static const unsigned int thresholds[] =
{
     63, 63, 63, 64, 64, 64, 65, 65, 65, 66, /* 30-39 */
     66, 67, 67, 68, 68, 68, 69, 69, 69, 69, /* 40-49 */
     70, 70, 70, 71, 71, 71, 72, 72, 73, 73, /* 50-59 */
     74, 74, 75, 75, 76, 76, 77, 77, 78, 79, /* 60-69 */
     80, 81, 82, 83, 84, /* 70-74 */
     85, 86, 87, 88, 89, /* 75-79 */
     91, 92, 93, 93, 94, /* 80-84 */
     95, 96, 97, 98,100, /* 85-89 */
     101, 102 /* 90-91 */
};

/*===========================================================================*/
/* Size of sieve to use divided by 2
 * Optimal if chosen to be a multiple of 32000 */
static const unsigned int sieveSize[] =
{
     64000,64000,64000,64000,64000,64000,64000,64000,64000,64000, /* 30-39 */
     64000,64000,64000,64000,64000,64000,64000,64000,64000,64000, /* 40-49 */
     64000,64000,64000,64000,64000,64000,64000,64000,64000,64000, /* 50-59 */
     64000,64000,64000,64000,64000,64000,64000,64000,64000,64000, /* 60-69 */
      64000,  64000,  64000,  64000,  64000, /* 70-74 */
      96000,  96000,  96000, 128000, 128000, /* 75-79 */
     160000, 160000, 160000, 160000, 160000, /* 80-84 */
     192000, 192000, 192000, 192000, 192000, /* 85-89 */
     192000, 192000 /* 90-91 */
};

/*===========================================================================*/
static unsigned int secondprime; /* cutoff for using flags when sieving */
static unsigned int firstprime;  /* first prime actually sieved with */
static unsigned char errorbits;  /* first prime actually sieved with */
static unsigned char threshold;  /* sieve threshold cutoff for smth relations */
static unsigned int largeprime;

static unsigned int *factorBase; /* array of factor base primes */
static unsigned char * primeSizes; /* array of sizes in bits of fb primes */

#define RELATIONS_PER_PRIME 100
static INLINE void set_relation(unsigned long* rel, unsigned int prime, unsigned int nrel, unsigned long val)
{
  if (nrel < RELATIONS_PER_PRIME)
    rel[ prime*RELATIONS_PER_PRIME + nrel ] = val;
}
static INLINE unsigned long get_relation(unsigned long* rel, unsigned int prime, unsigned int nrel)
{
  return rel[ prime*RELATIONS_PER_PRIME + nrel ];
}


/*=========================================================================
   Knuth_Schroeppel Multiplier:

   This is derived from Jason Papadopoulos's mpqs K-S method.  I believe it
   does a slightly better job than the K-S in FLINT 2.3, but that's debatable.
   An alternative would be to implement the method directly from Silverman 1987.

==========================================================================*/
/* Multiplers should be small square-free numbers, i.e.
 *    do { say $_ if moebius($_) != 0 } for 1..100
 * but SIMPQS doesn't deal well with composite multipliers.  So, just primes.
 */
static const unsigned long multipliers[] = {
  1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59,
  61, 67, 71, 73, 79, 83, 89, 97};
#define NUMMULTS (sizeof(multipliers)/sizeof(unsigned long))

#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif

static unsigned long knuthSchroeppel(mpz_t n, unsigned long numPrimes)
{
  unsigned int i, j, best_mult, knmod8;
  unsigned int maxprimes = (2*numPrimes <= 1000) ? 2*numPrimes : 1000;
  float best_score, contrib;
  float scores[NUMMULTS];
  mpz_t temp;

  mpz_init(temp);

  for (i = 0; i < NUMMULTS; i++) {
    scores[i] = 0.5 * logf((float)multipliers[i]);
    mpz_mul_ui(temp, n, multipliers[i]);
    knmod8 = mpz_mod_ui(temp, temp, 8);
    switch (knmod8) {
      case 1:  scores[i] -= 2 * M_LN2;  break;
      case 5:  scores[i] -= M_LN2;      break;
      case 3:
      case 7:  scores[i] -= 0.5 * M_LN2; break;
      default: break;
    }
  }

  {
    unsigned long prime, modp, knmodp;
    PRIME_ITERATOR(iter);
    for (i = 1; i < maxprimes; i++) {
      prime = prime_iterator_next(&iter);
      modp = mpz_mod_ui(temp, n, prime);
      contrib = logf((float)prime) / (float)(prime-1);

      for (j = 0; j < NUMMULTS; j++) {
        knmodp = (modp * multipliers[j]) % prime;
        if (knmodp == 0) {
          scores[j] -= contrib;
        } else {
          mpz_set_ui(temp, knmodp);
          if (mpz_kronecker_ui(temp, prime) == 1)
            scores[j] -= 2*contrib;
        }
      }
    }
    prime_iterator_destroy(&iter);
  }
  mpz_clear(temp);

  best_score = 1000.0;
  best_mult = 1;
  for (i = 0; i < NUMMULTS; i++) {
    float score = scores[i];
    if (score < best_score) {
      best_score = score;
      best_mult = multipliers[i];
    }
  }
  /* gmp_printf("%Zd mult %lu\n", n, best_mult); */
  return best_mult;
}


/*========================================================================
   Initialize Quadratic Sieve:

   Function: Initialises the global gmp variables.

========================================================================*/
static void initFactorBase(void)
{
    factorBase = 0;
    primeSizes = 0;
}
static void clearFactorBase(void)
{
    if (factorBase) { Safefree(factorBase);  factorBase = 0; }
    if (primeSizes) { Safefree(primeSizes);  primeSizes = 0; }
}

/*========================================================================
   Compute Factor Base:

   Function: Computes primes p up to B for which n is a square mod p,
   allocates memory and stores them in an array pointed to by factorBase.
   Additionally allocates and computes the primeSizes array.
   Returns: number of primes actually in the factor base

========================================================================*/
static void computeFactorBase(mpz_t n, unsigned long B,unsigned long multiplier)
{
  UV p;
  UV primesinbase = 0;
  PRIME_ITERATOR(iter);

  if (factorBase) { Safefree(factorBase);  factorBase = 0; }
  New(0, factorBase, B, unsigned int);

  factorBase[primesinbase++] = multiplier;
  if (multiplier != 2)
    factorBase[primesinbase++] = 2;
  prime_iterator_setprime(&iter, 3);
  for (p = 3; primesinbase < B; p = prime_iterator_next(&iter)) {
    if (mpz_kronecker_ui(n, p) == 1)
      factorBase[primesinbase++] = p;
  }
  prime_iterator_destroy(&iter);
#ifdef LARGESTP
  gmp_printf("Largest prime less than %Zd\n",p);
#endif

  /* Allocate and compute the number of bits required to store each prime */
  New(0, primeSizes, B, unsigned char);
  for (p = 0; p < B; p++)
    primeSizes[p] =
      (unsigned char) floor( log(factorBase[p]) / log(2.0) - SIZE_FUDGE + 0.5 );
}

/*===========================================================================
   Tonelli-Shanks:

   Function: Performs Tonelli-Shanks on n mod every prime in the factor base

===========================================================================*/
static void tonelliShanks(unsigned long numPrimes, mpz_t n, mpz_t * sqrts)
{
  unsigned long i;
  mpz_t fbprime, t1, t2, t3, t4;

  mpz_init(fbprime);
  mpz_init(t1); mpz_init(t2); mpz_init(t3); mpz_init(t4);

  mpz_set_ui(sqrts[0], 0);
  for (i = 1; i < numPrimes; i++) {
    mpz_set_ui(fbprime, factorBase[i]);
    sqrtmod(sqrts[i], n, fbprime, t1, t2, t3, t4);
  }
  mpz_clear(t1); mpz_clear(t2); mpz_clear(t3); mpz_clear(t4);
  mpz_clear(fbprime);
}

/*==========================================================================
   evaluateSieve:

   Function: searches sieve for relations and sticks them into a matrix, then
             sticks their X and Y values into two arrays XArr and YArr

===========================================================================*/
static void evaluateSieve(
    unsigned long numPrimes,
    unsigned long Mdiv2,
    unsigned long * relations,
    unsigned long ctimesreps,
    unsigned long M,
    unsigned char * sieve,
    mpz_t A,
    mpz_t B,
    mpz_t C,
    unsigned long * soln1,
    unsigned long * soln2,
    unsigned char * flags,
    matrix_t m,
    mpz_t * XArr,
    unsigned long * aind,
    int min,
    int s,
    int * exponents,
    unsigned long * npartials,
    unsigned long * nrelsfound,
    unsigned long * nrelssought,
    mpz_t temp,
    mpz_t temp2,
    mpz_t temp3,
    mpz_t res)
{
     long i,j,ii;
     unsigned int k;
     unsigned int exponent, vv;
     unsigned char extra;
     unsigned int modp;
     unsigned long * sieve2;
     unsigned char bits;
     int numfactors;
     unsigned long relsFound = *nrelsfound;
     unsigned long relSought = *nrelssought;

     mpz_set_ui(temp, 0);
     mpz_set_ui(temp2, 0);
     mpz_set_ui(temp3, 0);
     mpz_set_ui(res, 0);
     i = 0;
     j = 0;
     sieve2 = (unsigned long *) sieve;
#ifdef POLS
     gmp_printf("%Zdx^2%+Zdx\n%+Zd\n",A,B,C);
#endif

     while ( (unsigned long)j < M/sizeof(unsigned long))
     {
        do
        {
           while (!(sieve2[j] & SIEVEMASK)) j++;
           i = j * sizeof(unsigned long);
           j++;
           while (((unsigned long)i < j*sizeof(unsigned long)) && (sieve[i] < threshold)) i++;
        } while (sieve[i] < threshold);

        if (((unsigned long)i<M) && (relsFound < relSought))
        {
           mpz_set_ui(temp,i+ctimesreps);
           mpz_sub_ui(temp, temp, Mdiv2); /* X         */

           mpz_set(temp3, B);             /* B          */
           mpz_addmul(temp3, A, temp);    /* AX+B       */
           mpz_add(temp2, temp3, B);      /* AX+2B      */
           mpz_mul(temp2, temp2, temp);   /* AX^2+2BX   */
           mpz_add(res, temp2, C);        /* AX^2+2BX+C */

           bits = mpz_sizeinbase(res,2) - errorbits;

           numfactors=0;
           extra = 0;
           memset(exponents, 0, firstprime * sizeof(int));

           if (factorBase[0] != 1 && mpz_divisible_ui_p(res, factorBase[0]))
           {
             extra += primeSizes[0];
             if (factorBase[0] == 2) {
                exponent = mpz_scan1(res, 0);
                mpz_tdiv_q_2exp(res, res, exponent);
             } else {
               mpz_set_ui(temp,factorBase[0]);
               exponent = mpz_remove(res,res,temp);
             }
             exponents[0] = exponent;
           }

           exponents[1] = 0;
           if (mpz_divisible_ui_p(res, factorBase[1]))
           {
             extra += primeSizes[1];
             if (factorBase[1] == 2) {
                exponent = mpz_scan1(res, 0);
                mpz_tdiv_q_2exp(res, res, exponent);
             } else {
               mpz_set_ui(temp,factorBase[1]);
               exponent = mpz_remove(res,res,temp);
             }
             exponents[1] = exponent;
           }

           for (k = 2; k < firstprime; k++)
           {
              modp=(i+ctimesreps)%factorBase[k];

              exponents[k] = 0;
              if (soln2[k] != (unsigned long)-1)
              {
                 if ((modp==soln1[k]) || (modp==soln2[k]))
                 {
                    extra+=primeSizes[k];
                    mpz_set_ui(temp,factorBase[k]);
                    exponent = mpz_remove(res,res,temp);
                    CHECK_EXPONENT(exponent, k);
                    PRINT_FB(exponent, k);
                    exponents[k] = exponent;
                 }
              } else if (mpz_divisible_ui_p(res, factorBase[k]))
              {
                 extra += primeSizes[k];
                 mpz_set_ui(temp,factorBase[k]);
                 exponent = mpz_remove(res,res,temp);
                 PRINT_FB(exponent, k);
                 exponents[k] = exponent;
              }
           }
           sieve[i]+=extra;
           if (sieve[i] >= bits)
           {
              vv=((unsigned char)1<<(i&7));
              for (k = firstprime; (k<secondprime)&&(extra<sieve[i]); k++)
              {
                 modp=(i+ctimesreps)%factorBase[k];
                 if (soln2[k] != (unsigned long)-1)
                 {
                    if ((modp==soln1[k]) || (modp==soln2[k]))
                    {
                       extra+=primeSizes[k];
                       mpz_set_ui(temp,factorBase[k]);
                       exponent = mpz_remove(res,res,temp);
                       CHECK_EXPONENT(exponent, k);
                       PRINT_FB(exponent, k);
                       if (exponent)
                         for (ii = 0; ii < (long)exponent; ii++)
                           set_relation(relations, relsFound, ++numfactors, k);
                       if (exponent & 1)
                         insertEntry(m,relsFound,k);
                    }
                 } else if (mpz_divisible_ui_p(res, factorBase[k]))
                 {
                    extra += primeSizes[k];
                    mpz_set_ui(temp,factorBase[k]);
                    exponent = mpz_remove(res,res,temp);
                    PRINT_FB(exponent, k);
                    for (ii = 0; ii < (long)exponent; ii++)
                      set_relation(relations, relsFound, ++numfactors, k);
                    if (exponent & 1)
                      insertEntry(m,relsFound,k);
                 }
              }


              for (k = secondprime; (k<numPrimes)&&(extra<sieve[i]); k++)
              {
                 if (flags[k]&vv)
                 {
                    modp=(i+ctimesreps)%factorBase[k];
                    if ((modp==soln1[k]) || (modp==soln2[k]))
                    {
                       extra+=primeSizes[k];
                       mpz_set_ui(temp,factorBase[k]);
                       exponent = mpz_remove(res,res,temp);
                       CHECK_EXPONENT(exponent, k);
                       PRINT_FB(exponent, k);
                       if (exponent)
                         for (ii = 0; ii < (long)exponent; ii++)
                           set_relation(relations, relsFound, ++numfactors, k);
                       if (exponent & 1)
                         insertEntry(m,relsFound,k);
                    }
                 }
              }

              for (ii =0; ii<s; ii++)
              {
                 xorEntry(m,relsFound,aind[ii]+min);
                 set_relation(relations, relsFound, ++numfactors, aind[ii]+min);
              }

              if (mpz_cmp_ui(res,1000)>0)
              {
                 if (mpz_cmp_ui(res,largeprime)<0)
                 {
                    (*npartials)++;
                 }
                 clearRow(m,numPrimes,relsFound);
#ifdef RELPRINT
                 gmp_printf(" %Zd\n",res);
#endif
              } else
              {
                 mpz_neg(res,res);
                 if (mpz_cmp_ui(res,1000)>0)
                 {
                    if (mpz_cmp_ui(res,largeprime)<0)
                    {
                       (*npartials)++;
                    }
                    clearRow(m,numPrimes,relsFound);
#ifdef RELPRINT
                    gmp_printf(" %Zd\n",res);
#endif
                 } else
                 {
#ifdef RELPRINT
                    printf("....R\n");
#endif
                    for (ii = 0; ii < (long)firstprime; ii++)
                    {
                       int jj;
                       for (jj = 0; jj < exponents[ii]; jj++)
                         set_relation(relations, relsFound, ++numfactors, ii);
                       if (exponents[ii] & 1)
                         insertEntry(m,relsFound,ii);
                    }
                    set_relation(relations, relsFound, 0, numfactors);

                    mpz_init_set(XArr[relsFound], temp3);  /* (AX+B) */

                    relsFound++;
#ifdef COUNT
                    if (relsFound%20==0) fprintf(stderr,"%lu relations, %lu partials.\n", relsFound, *npartials);
#endif
                 }
              }
           } else
           {
              clearRow(m,numPrimes,relsFound);
#ifdef RELPRINT
              printf("\r                                                                    \r");
#endif

           }
           i++;

        } else if (relsFound >= relSought) i++;
     }
     /* Update caller */
     *nrelsfound = relsFound;
     *nrelssought = relSought;
}


static void update_solns(unsigned long first, unsigned long limit, unsigned long * soln1, unsigned long * soln2, int polyadd, const unsigned long * polycorr)
{
  unsigned int prime;
  unsigned long p, correction;

  for (prime = first; prime < limit; prime++) {
    if (soln2[prime] == (unsigned long) -1) continue;
    p = factorBase[prime];
    correction = (polyadd) ? p - polycorr[prime] : polycorr[prime];
    soln1[prime] += correction;
    while (soln1[prime] >= p)  soln1[prime] -= p;
    soln2[prime] += correction;
    while (soln2[prime] >= p)  soln2[prime] -= p;
  }
}

static void set_offsets(unsigned char * const sieve, const unsigned long * const soln1, const unsigned long * const soln2, unsigned char * * offsets1, unsigned char * * offsets2)
{
  unsigned int prime;
  for (prime = firstprime; prime < secondprime; prime++) {
    if (soln2[prime] == (unsigned long) -1) {
      offsets1[prime] = 0;
      offsets2[prime] = 0;
    } else {
      offsets1[prime] = sieve+soln1[prime];
      offsets2[prime] = sieve+soln2[prime];
    }
  }
}

/*=============================================================================
   Sieve:

   Function: Allocates space for a sieve of M integers and sieves the interval
             starting at start

=============================================================================*/
static void sieveInterval(unsigned long M, unsigned char * sieve, int more, unsigned char * * offsets1, unsigned char * * offsets2)
{
  unsigned int prime, p;
  unsigned char size;
  unsigned char * pos1;
  unsigned char * pos2;
  unsigned char * end = sieve + M;
  unsigned char * bound;
  ptrdiff_t diff;

  for (prime = firstprime; prime < secondprime; prime++)
  {
    if (offsets1[prime] == 0) continue;
    p    = factorBase[prime];
    size = primeSizes[prime];
    pos1 = offsets1[prime];
    pos2 = offsets2[prime];
    diff = pos2 - pos1;
    /* if pos1 < bound, then both *pos1 and *pos2 can be written to. */
    bound = (diff >= 0) ? end-diff : end;

    /* Write both values, unrolled 4 times. */
    bound -= (4-1)*p;
    while (pos1 < bound) {
      pos1[0  ] += size;  pos1[    diff] += size;
      pos1[1*p] += size;  pos1[1*p+diff] += size;
      pos1[2*p] += size;  pos1[2*p+diff] += size;
      pos1[3*p] += size;  pos1[3*p+diff] += size;
      pos1 += 4*p;
    }
    bound += (4-1)*p;
    /* Write both values */
    while (pos1 < bound) {
      pos1[0] += size;  pos1[diff] += size;  pos1 += p;
    }
    pos2 = pos1 + diff;    /* Restore pos2 */

    /* Finish writing to pos1 and pos2 */
    while (pos1 < end) {
      *pos1 += size; pos1 += p;
    }
    while (pos2 < end) {
      *pos2 += size; pos2 += p;
    }
    if (more) {
      offsets1[prime] = pos1;
      offsets2[prime] = pos2;
    }
  }
}

/*===========================================================================
   Sieve 2:

   Function: Second sieve for larger primes

=========================================================================== */
static void sieve2(unsigned long M, unsigned long numPrimes, unsigned char * sieve, const unsigned long * soln1, const unsigned long * soln2, unsigned char * flags)
{
     unsigned int prime;
     unsigned char *end = sieve + M;

     memset(flags, 0, numPrimes*sizeof(unsigned char));

     for (prime = secondprime; prime < numPrimes; prime++)
     {
        unsigned int  p    = factorBase[prime];
        unsigned char size = primeSizes[prime];
        unsigned char* pos1 = sieve + soln1[prime];
        unsigned char* pos2 = sieve + soln2[prime];

        if (soln2[prime] == (unsigned long)-1 ) continue;
        while (end - pos1 > 0)
        {
              flags[prime] |= ((unsigned char)1<<((pos1-sieve)&7));
              *pos1 += size;  pos1 += p;
        }

        while (end - pos2 > 0)
        {
              flags[prime] |= ((unsigned char)1<<((pos2-sieve)&7));
              *pos2 += size;  pos2 += p;
        }
     }
}

/*============================================================================

   random:

   Function: Generates a pseudo-random integer between 0 and n-1 inclusive

============================================================================*/
static unsigned long randval = 2994439072U;
static unsigned long silly_random(unsigned long upto)
{
   randval = ((unsigned long)randval*1025416097U+286824428U)%(unsigned long)4294967291U;
   return randval%upto;
}

/*============================================================================

  danaj: added these routines to reduce the set of factors to co-primes.
  It's not the most efficient solution, but it's trivial in time compared
  to the loop it's in, much less the rest of the QS.  It gives us a nice
  set of factors back, which is much more useful than the essentially
  random combinations we discover.

============================================================================*/

/* Verify that the factor reduction hasn't broken anything */
static void verify_factor_array(mpz_t n, mpz_t* farray, int nfacs)
{
  int i, j;
  mpz_t t;
  mpz_init_set_ui(t, 1);
  /* Assert we don't have duplicates */
  for (i = 0; i < nfacs; i++) {
    for (j = i+1; j < nfacs; j++) {
      if (mpz_cmp(farray[i],farray[j]) == 0) { gmp_printf("duplicate: F[%d] = F[%d] = %Zd\n", i, j, farray[i]); croak("assert"); }
    }
  }
  /* Assert that all factors multiply to n */
  for (i = 0; i < nfacs; i++)
    mpz_mul(t, t, farray[i]);
  if (mpz_cmp(t, n) != 0) { gmp_printf("farray doesn't multiply: n=%Zd t=%Zd\n", n, t); croak("assert"); }
  /* Assert that gcd of each non-identical factor is 1 */
  for (i = 0; i < nfacs; i++) {
    for (j = i+1; j < nfacs; j++) {
      if (mpz_cmp(farray[i],farray[j]) != 0) {
        mpz_gcd(t, farray[i], farray[j]);
        if (mpz_cmp_ui(t, 1) != 0) { gmp_printf("gcd: farray[%d] = %Zd  farray[%d] = %Zd\n", i, farray[i], j, farray[j]); croak("assert"); }
      }
    }
  }
  mpz_clear(t);
}
static int allprime_factor_array(mpz_t* farray, int nfacs)
{
  int i;
  for (i = 0; i < nfacs; i++) {
    if (!mpz_probab_prime_p(farray[i], 5))   /* Be lazy */
      return 0;
  }
  return 1;
}

static int insert_factor(mpz_t n, mpz_t* farray, int nfacs, mpz_t f)
{
  int i, j;
  mpz_t t, t2;

  if (mpz_cmp_ui(f, 1) <= 0)
    return nfacs;

  /* skip duplicates */
  for (i = 0; i < nfacs; i++)
    if (mpz_cmp(farray[i], f) == 0)
      break;
  if (i != nfacs) { return nfacs; }

  /* Look for common factors in all the existing set */
  /* for (i = 0; i < nfacs; i++) gmp_printf("  F[%d] = %Zd\n", i, farray[i]); */
  mpz_init(t);  mpz_init(t2);
  for (i = 0; i < nfacs; i++) {
    mpz_gcd(t, farray[i], f);
    if (mpz_cmp_ui(t, 1) == 0) /* t=1:   F and f unchanged */
      continue;
    mpz_divexact(t2, farray[i], t);    /* t2 = F/t */
    mpz_divexact(f, f, t);             /* f  = f/t */
    /* Remove the old farray[i] */
    for (j = i+1; j < nfacs; j++)
      mpz_set(farray[j-1], farray[j]);
    mpz_set_ui(farray[nfacs--], 0);
    /* Insert F/t, t, f/t */
    nfacs = insert_factor(n, farray, nfacs, t2);
    nfacs = insert_factor(n, farray, nfacs, t);
    nfacs = insert_factor(n, farray, nfacs, f);
    i=0;
    break;
  }
  /* If nothing common, insert it. */
  if (i == nfacs)
    mpz_set(farray[nfacs++], f);
  mpz_clear(t);  mpz_clear(t2);
  return nfacs;
}


/*============================================================================
   mainRoutine:

   Function: Generates the polynomials, initialises and calls the sieve,
             implementing cache blocking (breaking the sieve interval into
             small blocks for the small primes.

============================================================================*/
static int mainRoutine(
  unsigned long numPrimes,
  unsigned long Mdiv2,
  unsigned long relSought,
  mpz_t n,
  mpz_t* farray,
  unsigned long multiplier)
{
    mpz_t A, B, C, D, Bdivp2, q, r, nsqrtdiv, temp, temp2, temp3, temp4;
    int i, j, l, s, fact, span, min, nfactors, verbose;
    unsigned long u1, p, reps, numRelations, M;
    unsigned long curves = 0;
    unsigned long npartials = 0;
    unsigned long relsFound = 0;
    unsigned long  * relations;
    unsigned short * primecount;
    unsigned char  * sieve;
    int            * exponents;
    unsigned long  * aind;
    unsigned long  * amodp;
    unsigned long  * Ainv;
    unsigned long  * soln1;
    unsigned long  * soln2;
    unsigned char  * flags;
    unsigned long ** Ainv2B;
    unsigned char ** offsets;
    unsigned char ** offsets2;
    mpz_t          * XArr;
    mpz_t          * Bterms;
    mpz_t          * sqrts;
    matrix_t m;

    verbose = get_verbose_level();
    s = mpz_sizeinbase(n,2)/28+1;

    New(  0, exponents, firstprime, int );
    Newz( 0, aind,          s, unsigned long );
    Newz( 0, amodp,         s, unsigned long );
    Newz( 0, Ainv,  numPrimes, unsigned long );
    Newz( 0, soln1, numPrimes, unsigned long );
    Newz( 0, soln2, numPrimes, unsigned long );
    Newz( 0, Ainv2B,        s, unsigned long*);
    Newz( 0, XArr,  relSought, mpz_t );
    New(  0, Bterms,        s, mpz_t );
    if (exponents == 0 || aind == 0 || amodp == 0 || Ainv == 0 ||
        soln1 == 0 || soln2 == 0 || Ainv2B == 0 || Bterms == 0 ||
        XArr == 0)
      croak("SIMPQS: Unable to allocate memory!\n");

    flags = 0;
    if (secondprime < numPrimes) {
      New(0, flags, numPrimes, unsigned char);
      if (flags == 0) croak("SIMPQS: Unable to allocate memory!\n");
    }

    for (i=0; i<s; i++)
    {
       New(0, Ainv2B[i], numPrimes, unsigned long);
       if (Ainv2B[i] == 0) croak("SIMPQS: Unable to allocate memory!\n");
       mpz_init(Bterms[i]);
    }

    m = constructMat(numPrimes, relSought);

    /* One extra word for sentinel */
    Newz(0, sieve,     Mdiv2*2 + sizeof(unsigned long), unsigned char);
    New( 0, offsets,   secondprime, unsigned char*);
    New( 0, offsets2,  secondprime, unsigned char*);
    Newz(0, relations, relSought * RELATIONS_PER_PRIME, unsigned long);

    if (sieve == 0 || offsets == 0 || offsets2 == 0 || relations == 0)
      croak("SIMPQS: Unable to allocate memory!\n");

    mpz_init(A); mpz_init(B); mpz_init(C); mpz_init(D);
    mpz_init(Bdivp2); mpz_init(q); mpz_init(r); mpz_init(nsqrtdiv);
    mpz_init(temp); mpz_init(temp2); mpz_init(temp3); mpz_init(temp4);

    /* Compute sqrt(n) mod factorbase[i] */
    New(0, sqrts, numPrimes, mpz_t);
    if (sqrts == 0) croak("SIMPQS: Unable to allocate memory!\n");
    for (p = 0; p < numPrimes; p++)
      mpz_init(sqrts[p]);
    tonelliShanks(numPrimes, n, sqrts);

    /* Compute min A_prime and A_span */

    mpz_mul_ui(temp,n,2);
    mpz_sqrt(temp,temp);
    mpz_div_ui(nsqrtdiv,temp,Mdiv2);
    mpz_root(temp,nsqrtdiv,s);
    for (fact = 0; mpz_cmp_ui(temp,factorBase[fact])>=0; fact++);
    span = numPrimes/s/s/2;
    min=fact-span/2;
    while ( min > 0 && (fact*fact)/min - min < span )
      min--;

#ifdef ADETAILS
    printf("s = %d, fact = %d, min = %d, span = %d\n",s,fact,min,span);
#endif

    /* Compute first polynomial and adjustments */

    while (relsFound < relSought)
    {
        int polyindex;
        mpz_set_ui(A,1);
        for (i = 0; i < s-1; )
        {
           unsigned long ran = span/2+silly_random(span/2);
           j=-1L;
           while (j!=i)
           {
              ran++;
              for (j=0;((j<i)&&(aind[j]!=ran));j++);
           }
           aind[i] = ran;
           mpz_mul_ui(A,A,factorBase[ran+min]);
           i++;
           if (i < s-1)
           {
              j=-1L;
              ran = ((min+span/2)*(min+span/2))/(ran+min) - silly_random(10)-min;
              while (j!=i)
              {
                 ran++;
                 for (j=0;((j<i)&&(aind[j]!=ran));j++);
              }
              aind[i] = ran;
              mpz_mul_ui(A,A,factorBase[ran+min]);
              i++;
           }
        }
        mpz_div(temp,nsqrtdiv,A);
        for (fact = 1; mpz_cmp_ui(temp,factorBase[fact])>=0; fact++);
        fact-=min;
        do
        {
           for (j=0;((j<i)&&(aind[j]!=(unsigned long)fact));j++);
           fact++;
        } while (j!=i);
        fact--;
        aind[i] = fact;
        mpz_mul_ui(A,A,factorBase[fact+min]);

        for (i=0; i<s; i++)
        {
           p = factorBase[aind[i]+min];
           mpz_div_ui(temp,A,p);
           amodp[i] = mpz_fdiv_r_ui(temp,temp,p);

           mpz_set_ui(temp,modinverse(mpz_get_ui(temp),p));
           mpz_mul(temp, temp, sqrts[aind[i]+min]);
           mpz_fdiv_r_ui(temp, temp, p);
           if (mpz_cmp_ui(temp,p/2)>0)
           {
              mpz_sub_ui(temp,temp,p);
              mpz_neg(temp,temp);
           }
           mpz_mul(temp,temp,A);
           mpz_div_ui(Bterms[i],temp,p);
        }

        mpz_set(B,Bterms[0]);
        for (i = 1; i < s; i++)
        {
           mpz_add(B,B,Bterms[i]);
        }

        for (i = 0; i < (int)numPrimes; i++)
        {
           p = factorBase[i];
           Ainv[i] = modinverse(mpz_fdiv_r_ui(temp,A,p),p);

           for (j=0; j<s; j++)
           {
              mpz_fdiv_r_ui(temp,Bterms[j],p);
              mpz_mul_ui(temp,temp,2*Ainv[i]);
              Ainv2B[j][i] = mpz_fdiv_r_ui(temp,temp,p);
           }

           mpz_fdiv_r_ui(temp,B,p);
           mpz_sub(temp,sqrts[i],temp);
           mpz_add_ui(temp,temp,p);
           mpz_mul_ui(temp,temp,Ainv[i]);
           mpz_add_ui(temp,temp,Mdiv2);
           soln1[i] = mpz_fdiv_r_ui(temp,temp,p);
           mpz_sub_ui(temp,sqrts[i],p);
           mpz_neg(temp,temp);
           mpz_mul_ui(temp,temp,2*Ainv[i]);
           soln2[i] = mpz_fdiv_r_ui(temp,temp,p)+soln1[i];
        }

        for (polyindex=1; polyindex<(1<<(s-1))-1; polyindex++)
        {
           int polyadd;
           unsigned long * polycorr;
           for (j=0; j<s; j++)
           {
              if (((polyindex>>j)&1)!=0) break;
           }
           if ((polyadd = (((polyindex>>j)&2)!=0)))
           {
              mpz_add(B,B,Bterms[j]);
              mpz_add(B,B,Bterms[j]);
           } else
           {
              mpz_sub(B,B,Bterms[j]);
              mpz_sub(B,B,Bterms[j]);
           }
           polycorr = Ainv2B[j];

           for (j=0; j<s; j++)
           {
              int findex = aind[j]+min;
              p = factorBase[findex];
              mpz_fdiv_r_ui(D,n,p*p);
              mpz_fdiv_r_ui(Bdivp2,B,p*p);
              mpz_mul_ui(temp,Bdivp2,amodp[j]);
              mpz_fdiv_r_ui(temp,temp,p);
              u1 = modinverse(mpz_fdiv_r_ui(temp,temp,p),p);
              mpz_mul(temp,Bdivp2,Bdivp2);
              mpz_sub(temp,temp,D);
              mpz_neg(temp,temp);
              mpz_div_ui(temp,temp,p);
              mpz_mul_ui(temp,temp,u1);
              mpz_add_ui(temp,temp,Mdiv2);
              mpz_add_ui(temp,temp,p);
              soln1[findex]=mpz_fdiv_r_ui(temp,temp,p);
              soln2[findex] = (unsigned long) -1;
           }

           /* Count the number of polynomial curves used so far and compute
            * the C coefficient of our polynomial */

           curves++;

           mpz_mul(C,B,B);
           mpz_sub(C,C,n);
           mpz_divexact(C,C,A);

           /* Do the sieving and relation collection */

           mpz_set_ui(temp,Mdiv2*2);
           mpz_fdiv_qr_ui(q,r,temp,CACHEBLOCKSIZE);
           M = mpz_get_ui(temp);

           /* set the solns1 and solns2 arrays */
           update_solns(1, numPrimes, soln1, soln2, polyadd, polycorr);
           /* Clear sieve and insert sentinel at end (used in evaluateSieve) */
           memset(sieve, 0, M*sizeof(unsigned char));
           sieve[M] = 255;
           /* Sieve [secondprime , numPrimes) */
           if (secondprime < numPrimes)
             sieve2(M, numPrimes, sieve, soln1, soln2, flags);
           /* Set the offsets and offsets2 arrays used for small sieve */
           set_offsets(sieve, soln1, soln2, offsets, offsets2);
           /* Sieve [firstprime , secondprime) */
           sieveInterval(CACHEBLOCKSIZE,sieve,1,offsets,offsets2);
           if (mpz_cmp_ui(q,1)>0)
           {
              unsigned long maxreps = mpz_get_ui(q)-1;
              for (reps = 1; reps < maxreps; reps++)
              {
                 sieveInterval(CACHEBLOCKSIZE,sieve+CACHEBLOCKSIZE*reps,1,offsets,offsets2);
              }
              if (mpz_cmp_ui(r,0)==0)
              {
                 sieveInterval(CACHEBLOCKSIZE,sieve+CACHEBLOCKSIZE*reps,0,offsets,offsets2);
              } else
              {
                 sieveInterval(CACHEBLOCKSIZE,sieve+CACHEBLOCKSIZE*reps,1,offsets,offsets2);
                 reps++;
                 sieveInterval(mpz_get_ui(r),sieve+CACHEBLOCKSIZE*reps,0,offsets,offsets2);
              }
           }

           evaluateSieve(
              numPrimes, Mdiv2,
              relations, 0, M, sieve, A, B, C,
              soln1, soln2, flags, m, XArr, aind,
              min, s, exponents,
              &npartials, &relsFound, &relSought,
              temp, temp2, temp3, temp4
           );
        }

#ifdef COUNT
        if (curves%20==0) printf("%ld curves.\n",(long)curves);
#endif
    }

#ifdef CURPARTS
    printf("%lu curves, %lu partials.\n", curves, npartials);
#endif

#ifdef REPORT
    printf("Done with sieving!\n");
#endif
    if (verbose>3) printf("# qs done sieving\n");

    /* Free everything we don't need for the linear algebra */

    for (p = 0; p < numPrimes; p++)
      mpz_clear(sqrts[p]);
    Safefree(sqrts);
    for (i = 0; i < s; i++) {
      Safefree(Ainv2B[i]);
      mpz_clear(Bterms[i]);
    }
    Safefree(exponents);
    Safefree(aind);
    Safefree(amodp);
    Safefree(Ainv);
    Safefree(soln1);
    Safefree(soln2);
    Safefree(Ainv2B);
    Safefree(Bterms);
    if (flags) Safefree(flags);

    Safefree(sieve);    sieve = 0;
    Safefree(offsets);  offsets = 0;
    Safefree(offsets2); offsets2 = 0;

    mpz_clear(A);  mpz_clear(B);  mpz_clear(C);  mpz_clear(D);
    mpz_clear(q);  mpz_clear(r);
    mpz_clear(Bdivp2); mpz_clear(nsqrtdiv);

    /* Do the matrix algebra step */

    numRelations = gaussReduce(m, numPrimes, relSought);
#ifdef REPORT
    printf("%ld relations in kernel.\n", numRelations);
#endif
    if (verbose>3) printf("# qs found %lu relations in kernel\n", numRelations);

    /* We want factors of n, not kn, so divide out by the multiplier */

    mpz_div_ui(n,n,multiplier);

    /* Now do the "sqrt" and GCD steps hopefully obtaining factors of n */
    mpz_set(farray[0], n);
    nfactors = 1;  /* We have one result -- n */
    New( 0, primecount, numPrimes, unsigned short);
    if (primecount == 0) croak("SIMPQS: Unable to allocate memory!\n");
    for (l = (int)relSought-64; l < (int)relSought; l++)
    {
        unsigned int mat2offset = rightMatrixOffset(numPrimes);
        mpz_set_ui(temp,1);
        mpz_set_ui(temp2,1);
        memset(primecount,0,numPrimes*sizeof(unsigned short));
        for (i = 0; i< (int)numPrimes; i++)
        {
           if (getEntry(m,l,mat2offset+i))
           {
              int nrelations = get_relation(relations, i, 0);
              if (nrelations >= RELATIONS_PER_PRIME)
                nrelations = RELATIONS_PER_PRIME-1;
              mpz_mul(temp2,temp2,XArr[i]);
              for (j = 1; j <= nrelations; j++)
                primecount[ get_relation(relations, i, j) ]++;
           }
           if (i%16==0) mpz_mod(temp2,temp2,n);
        }
        for (j = 0; j < (int)numPrimes; j++)
        {
           mpz_set_ui(temp3,factorBase[j]);
           mpz_pow_ui(temp3,temp3,primecount[j]/2);
           mpz_mul(temp,temp,temp3);
           if (j%16==0) mpz_mod(temp,temp,n);
        }
        mpz_sub(temp,temp2,temp);
        mpz_gcd(temp,temp,n);
        /* Only non-trivial factors */
        if (mpz_cmp_ui(temp,1) && mpz_cmp(temp,n) && mpz_divisible_p(n,temp) ) {
          if (verbose>4) gmp_printf("# qs factor %Zd\n", temp);
          nfactors = insert_factor(n, farray, nfactors, temp);
          verify_factor_array(n, farray, nfactors);
          if (allprime_factor_array(farray, nfactors))
            break;
        }
    }

    /* Free everything remaining */
    Safefree(primecount);

    destroyMat(m, relSought);
    Safefree(relations);

    for (i = 0; i < (int)relSought; i++) {
      mpz_clear(XArr[i]);
    }
    Safefree(XArr);

    mpz_clear(temp);  mpz_clear(temp2);  mpz_clear(temp3);  mpz_clear(temp4);

    return nfactors;
}

int _GMP_simpqs(mpz_t n, mpz_t* farray)
{
  unsigned long numPrimes, Mdiv2, multiplier, decdigits, relSought;
  int result = 0;
  int verbose = get_verbose_level();

  mpz_set(farray[0], n);
  decdigits = mpz_sizeinbase(n,10); /* often 1 too big */
  if (decdigits < MINDIG)
    return 0;

  if (verbose>2) gmp_printf("# qs trying %Zd (%lu digits)\n", n, decdigits);
#ifdef REPORT
  gmp_printf("%Zd (%ld decimal digits)\n", n, decdigits);
#endif

  /* It's important to remove small factors. */
  {
    UV p;
    PRIME_ITERATOR(iter);
    for (p = 2; p < 1000; p = prime_iterator_next(&iter)) {
      if (mpz_cmp_ui(n, p*p) < 0) break;
      while (mpz_divisible_ui_p(n, p)) {
        mpz_set_ui(farray[result++], p);
        mpz_divexact_ui(n, n, p);
      }
    }
    decdigits = mpz_sizeinbase(n,10);
    if (decdigits < MINDIG)
      return result;
    mpz_set(farray[result], n);
  }

  /* Get a preliminary number of primes, pick a multiplier, apply it */
  numPrimes = (decdigits <= 91) ? primesNo[decdigits-MINDIG] : 64000;
  multiplier = knuthSchroeppel(n, numPrimes);
  mpz_mul_ui(n, n, multiplier);
  decdigits = mpz_sizeinbase(n, 10);

  if (decdigits<=91) {
    numPrimes=primesNo[decdigits-MINDIG];

    Mdiv2 = sieveSize[decdigits-MINDIG]/SIEVEDIV;
    if (Mdiv2*2 < CACHEBLOCKSIZE) Mdiv2 = CACHEBLOCKSIZE/2;
    largeprime = 1000 * largeprimes[decdigits-MINDIG];

    secondprime = (numPrimes < SECONDPRIME) ? numPrimes : SECONDPRIME;

    firstprime = firstPrimes[decdigits-MINDIG];
    errorbits = errorAmounts[decdigits-MINDIG];
    threshold = thresholds[decdigits-MINDIG];
  } else {
    numPrimes = 64000;
    Mdiv2 = 192000/SIEVEDIV;
    largeprime = numPrimes*10*decdigits;

    secondprime = SECONDPRIME;
    firstprime = 30;
    errorbits = decdigits/4 + 2;
    threshold = 43+(7*decdigits)/10;
  }

#ifdef REPORT
  printf("Using multiplier: %lu\n",multiplier);
  printf("%lu primes in factor base.\n",numPrimes);
  printf("Sieving interval M = %lu\n",Mdiv2*2);
  printf("Large prime cutoff = factorBase[%u]\n",largeprime);
#endif
  if (verbose>2) gmp_printf("# qs    mult %lu, digits %lu, sieving %lu, primes %lu\n", multiplier, decdigits, Mdiv2*2, numPrimes);

  /* We probably need fewer than this */
  relSought = numPrimes;
  initFactorBase();
  computeFactorBase(n, numPrimes, multiplier);

  result += mainRoutine(numPrimes, Mdiv2, relSought, n, farray+result, multiplier);

  clearFactorBase();
  if (verbose>2) {
    int i;
    gmp_printf("# qs:");
    for (i = 0; i < result; i++)
      gmp_printf(" %Zd", farray[i]);
    gmp_printf("%s\n", (result) ? "" : " no factors");
  }
  /* if (!result) gmp_printf("QS Fail: %Zd (%ld digits)\n", n, decdigits); */
  return result;
}

#ifdef STANDALONE_SIMPQS
/*===========================================================================
   Main Program:

   Function: Factors a user specified number using a quadratic sieve

===========================================================================*/
int main(int argc, char **argv)
{
  int i, nfactors;
  mpz_t n;
  mpz_t* farray;

  mpz_init(n);
  New(0, farray, 64, mpz_t);
  for (i = 0; i < 64; i++)
    mpz_init_set_ui(farray[i], 0);

  printf("Input number to factor [ >=%d decimal digits]: ", MINDIG);
  gmp_scanf("%Zd",n);getchar();

  if (mpz_sizeinbase(n,10) < MINDIG)
    croak("SIMPQS: Error in input or number has too few digits.\n");

  nfactors = _GMP_simpqs(n, farray);

  for (i = 0; i < nfactors; i++)
    gmp_printf("  %Zd\n", farray[i]);

  for (i = 0; i < 64; i++)
    mpz_clear(farray[i]);
  Safefree(farray);

  return 0;
}
#endif
