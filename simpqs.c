/*============================================================================
    Copyright 2006 William Hart

    This file is part of SIMPQS.

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

/*============================================================================
   Some modifications made in 2013 by Dana Jacobsen:
     - put it in one file
     - merge some of the 2.0 changes
     - make it work with smaller values
     - fix some memory errors
     - free memory all over
     - fewer globals
     - mpz_nextprime is slow, slow, slow.  Use prime_iterator.
     - Alternate multiplier selection routine.
     - lots of little changes / optimizations

   Version 2.0 scatters temp files everywhere, but that could be solved.
   The main benefits left in 2.0 are:
      (1) combining partial relations (this is huge for large inputs)
      (2) much less memory use, though partly due to using temp files
      (3) jasonp's block Lanczos routine.
   This code goes through curves faster than v2.0, but with big inputs it
   ends up needing 2x the time because of not combining partials as well as
   the final linear algebra time.
   TODO: Tune 25-40 digit parameters
============================================================================*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <gmp.h>

#include "ptypes.h"
#include "simpqs.h"
#include "prime_iterator.h"

#ifdef STANDALONE_SIMPQS
  typedef unsigned long UV;
  typedef   signed long IV;
  #define UV_MAX ULONG_MAX
  #define UVCONST(x) ((unsigned long)x##UL)
  #define croak(fmt,...)            { printf(fmt,##__VA_ARGS__); exit(1); }
#else
  #include "EXTERN.h"
  #include "perl.h"
  #include "XSUB.h"
#endif


static void modmul(mpz_t ab, mpz_t a, mpz_t b, mpz_t p)
{
     mpz_mul(ab,a,b);
     mpz_fdiv_r(ab,ab,p);
}

/* variables for sqrtmod */
static mpz_t two;
static mpz_t p1;
static mpz_t b;
static mpz_t g;
static mpz_t mk;
static mpz_t bpow;
static mpz_t gpow;

/* - Initialises variables */
static void TonelliInit(void)
{
    mpz_init_set_ui(two, 2);
    mpz_init(p1);
    mpz_init(b);
    mpz_init(g);
    mpz_init(mk);
    mpz_init(bpow);
    mpz_init(gpow);
}
static void TonelliDestroy(void)
{
    mpz_clear(two);
    mpz_clear(p1);
    mpz_clear(b);
    mpz_clear(g);
    mpz_clear(mk);
    mpz_clear(bpow);
    mpz_clear(gpow);
}

/* - Tonelli-Shanks: sets asqrt to a square root of a modulo p */
/*- Return: 0 if a is not a square mod p, 1 otherwise. */
static int sqrtmod(mpz_t asqrt, mpz_t a, mpz_t p)
{
     int r,k,m,i;

     if (mpz_kronecker(a,p)!=1)
     {
         mpz_set_ui(asqrt,0);
         return 0;   /* return 0 if a is not a square mod p */
     }

     mpz_sub_ui(p1,p,1);
     r = mpz_remove(p1,p1,two);
     mpz_powm(b,a,p1,p);
     for (k=2; ;k++)
     {
         if (mpz_ui_kronecker(k,p) == -1) break;
     }
     mpz_set_ui(mk,k);
     mpz_powm(g,mk,p1,p);
     mpz_add_ui(p1,p1,1);
     mpz_divexact_ui(p1,p1,2);
     mpz_powm(asqrt,a,p1,p);
     if (!mpz_cmp_ui(b,1))
     {
          return 1;
     }

     while (mpz_cmp_ui(b,1))
     {
           mpz_set(bpow,b);
           for (m=1; (m<=r-1) && (mpz_cmp_ui(bpow,1));m++)
           {
               mpz_powm_ui(bpow,bpow,2,p);
           }
           mpz_set(gpow,g);
           for (i = 1;i<= r-m-1;i++)
           {
               mpz_powm_ui(gpow,gpow,2,p);
           };
           modmul(asqrt,asqrt,gpow,p);
           mpz_powm_ui(gpow,gpow,2,p);
           modmul(b,b,gpow,p);
           mpz_set(gpow,g);
           r = m;
     }

     return 1;
}

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
static void clearRow(matrix_t m, unsigned int numcols, unsigned int row)
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
  m = (row_t *) malloc( rows * sizeof(row_t) );
  if (m == 0) croak("SIMPQS: Unable to allocate memory for matrix!\n");

  for (i = 0; i < rows; i++) { /* two matrices, side by side */
    m[i] = (row_t) calloc( 2*nbytes, sizeof(unsigned char) );
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
    free(m[i]);
  free(m);
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
static unsigned int gaussReduce(matrix_t m, unsigned int numPrimes, unsigned int relSought)
{
  unsigned int rowUpto = 0;
  unsigned int irow, checkRow;
  int icol;

  for (icol = numPrimes-1; icol >= 0; icol--) {
    irow = rowUpto;

    while ( (irow < relSought) && (getEntry(m,irow,icol) == 0) )
      irow++;

    if (irow < relSought) {
      swapRows(m,rowUpto,irow);
      for (checkRow = rowUpto+1; checkRow < relSought; checkRow++) {
        if (getEntry(m,checkRow,icol) != 0)
          xorRows(m, numPrimes, rowUpto, checkRow);
      }
      rowUpto++;
    }
  }
  return rowUpto;
}

//===========================================================================
//Uncomment these for various pieces of debugging information

//#define COUNT    // Shows the number of relations generated and curves used during sieving
//#define RELPRINT // Shows the actual factorizations of the relations
//#define ERRORS   // Error if relation should be divisible by a prime but isn't
//#define POLS     // Shows the polynomials being used by the sieve
//#define ADETAILS // Prints some details about the factors of the A coefficients of the polys
//#define LARGESTP // Prints the size of the largest factorbase prime
//#define CURPARTS // Prints the number of curves used and number of partial relations
//#define REPORT //report sieve size, multiplier and number of primes used

//===========================================================================
//Architecture dependent fudge factors

#if ULONG_MAX == 4294967295UL
#define SIEVEMASK 0xC0C0C0C0UL
#define SIEVEDIV 1
#elif ULONG_MAX == 18446744073709551615UL
#define SIEVEMASK 0xC0C0C0C0C0C0C0C0UL
#define SIEVEDIV 1
#else
 #error Cannot determine ulong size
#endif

#define CACHEBLOCKSIZE 64000 //Should be a little less than the L1/L2 cache size
                             //and a multiple of 64000
#define SECONDPRIME    6000 //This should be lower for slower machines
#define FUDGE          0.15 //Every program needs a mysterious fudge factor

#define MINDIG 30 //Will not factor numbers with less than this number of decimal digits

//===========================================================================
// Large prime cutoffs

static unsigned int largeprimes[] =
{
     100000, 100000, 125000, 125000, 150000, 150000, 175000, 175000, 200000, 200000, //30-39
     250000, 300000, 370000, 440000, 510000, 580000, 650000, 720000, 790000, 8600000, //40-49
     930000, 1000000, 1700000, 2400000, 3100000, 3800000, 4500000, 5200000, 5900000, 6600000, //50-59
     7300000, 8000000, 8900000, 10000000, 11300000, 12800000, 14500000, 16300000, 18100000, 20000000, //60-69
     22000000, 24000000, 27000000, 32000000, 39000000,  //70-74
     53000000, 65000000, 75000000, 87000000, 100000000, //75-79
     114000000, 130000000, 150000000, 172000000, 195000000, //80-84
     220000000, 250000000, 300000000, 350000000, 400000000, //85-89
     450000000, 500000000 //90-91
};

//============================================================================
// Number of primes to use in factor base, given the number of decimal digits specified
static unsigned int primesNo[] =
{
     1500, 1600, 1600, 1600, 1600, 1600, 1600, 1600, 1600, 1600, //30-39
     1600, 1600, 1600, 1700, 1750, 1800, 1900, 2000, 2050, 2100, //40-49
     2150, 2200, 2250, 2300, 2400, 2500, 2600, 2700, 2800, 2900, //50-59
     3000, 3150, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, //60-69
     9500, 10000, 11500, 13000, 15000, //70-74
     17000, 24000, 27000, 30000, 37000, //75-79
     45000, 47000, 53000, 57000, 58000,  //80-84
     59000, 60000, 64000, 68000, 72000,  //85-89
     76000, 80000 //90-91
};

//============================================================================
// First prime actually sieved for
static unsigned int firstPrimes[] =
{
     5, 5, 5, 6, 6, 6, 6, 7, 7, 7, //30-39
     8, 8, 8, 8, 8, 8, 8, 8, 8, 8, //40-49
     9, 8, 9, 9, 9, 9, 10, 10, 10, 10, //50-59
     10, 10, 11, 11, 12, 12, 13, 14, 15, 17, //60-69  //10
     19, 21, 22, 22, 23, //70-74
     24, 25, 25, 26, 26, //75-79
     27, 27, 27, 27, 28, //80-84
     28, 28, 28, 29, 29, //85-89
     29, 29 //90-91
};

//============================================================================
// Logs of primes are rounded and errors accumulate; this specifies how great an error to allow
static unsigned int errorAmounts[] =
{
     10, 10, 10, 11, 13, 14, 14, 15, 15, 16, //30-39
     16, 17, 17, 18, 18, 19, 19, 19, 20, 20, //40-49
     21, 21, 21, 22, 22, 22, 23, 23, 23, 24, //50-59
     24, 24, 25, 25, 25, 25, 26, 26, 26, 26, //60-69 //24
     27, 27, 28, 28, 29, //70-74
     29, 30, 30, 30, 31, //75-79
     31, 31, 31, 32, 32, //80-84
     32, 32, 32, 33, 33, //85-89
     33, 33 //90-91
};

//============================================================================
// This is the threshold the sieve value must exceed in order to be considered for smoothness
static unsigned int thresholds[] =
{
     63, 63, 63, 64, 64, 64, 65, 65, 65, 66, //30-39
     66, 67, 67, 68, 68, 68, 69, 69, 69, 69, //40-49
     70, 70, 70, 71, 71, 71, 72, 72, 73, 73, //50-59
     74, 74, 75, 75, 76, 76, 77, 77, 78, 79, //60-69 //74
     80, 81, 82, 83, 84, //70-74
     85, 86, 87, 88, 89, //75-79
     91, 92, 93, 93, 94, //80-84
     95, 96, 97, 98, 100, //85-89
     101, 102 //90-91
};

//============================================================================
// Size of sieve to use divided by 2, given the number of decimal digits specified
//N.B: probably optimal if chosen to be a multiple of 32000, though other sizes are supported
static unsigned int sieveSize[] =
{
     64000, 64000, 64000, 64000, 64000, 64000, 64000, 64000, 64000, 64000, //30-39
     64000, 64000, 64000, 64000, 64000, 64000, 64000, 64000, 64000, 64000, //40-49
     64000, 64000, 64000, 64000, 64000, 64000, 64000, 64000, 64000, 64000, //50-59
     64000, 64000, 64000, 64000, 64000, 64000, 64000, 64000, 64000, 64000, //60-69
     64000, 64000, 64000, 64000, 64000, //70-74
     96000, 96000, 96000, 128000, 128000, //75-79
     160000, 160000, 160000, 160000, 160000, //80-84
     192000, 192000, 192000, 192000, 192000, //85-89
     192000, 192000 //90-91
};

//============================================================================
static unsigned int decdigits;   //number of decimal digits of n
static unsigned int secondprime; //min(numprimes, SECONDPRIME) = cutoff for using flags when sieving
static unsigned int firstprime;  //first prime actually sieved with
static unsigned char errorbits;  //first prime actually sieved with
static unsigned char threshold;  //sieve threshold cutoff for smooth relations
static unsigned int largeprime;

static unsigned int *factorBase; //array of factor base primes
//static unsigned int numPrimes; //number of primes in factor base
static unsigned int relSought; //number of relations sought, i.e. a "few" more than numPrimes
static unsigned char * primeSizes; //array of sizes in bits, of the factor base primes
static unsigned int relsFound =0; //number of relations found so far
static unsigned char * flags; //flags used for speeding up sieving for large primes
static unsigned int partials = 0; //number of partial relations

static mpz_t * sqrts; //square roots of n modulo each prime in the factor base

#define RELATIONS_PER_PRIME 100
static void set_relation(unsigned long* rel, unsigned int prime, unsigned int nrel, unsigned long val)
{
  if (nrel < RELATIONS_PER_PRIME)
    rel[ prime*RELATIONS_PER_PRIME + nrel ] = val;
}
static unsigned long get_relation(unsigned long* rel, unsigned int prime, unsigned int nrel)
{
  return rel[ prime*RELATIONS_PER_PRIME + nrel ];
}


/*========================================================================
   Modular Inversion:

   Function: GMP has a modular inverse function, but believe it or not,
             this clumsy implementation is apparently quite a bit faster.
             It inverts the value a, modulo the prime p, using the extended
             gcd algorithm.

========================================================================*/

static unsigned long modinverse(unsigned long a, unsigned long p)
{
  long u1, u3;
  long v1, v3;
  long t1, t3, quot;
   u1=1; u3=a;
   v1=0; v3=p;
   t1=0; t3=0;
   while (v3)
   {
      quot=u3-v3;
      if (u3 < (v3<<2))
      {
         if (quot < v3)
         {
            if (quot < 0)
            {
               t1 = u1; u1 = v1; v1 = t1;
               t3 = u3; u3 = v3; v3 = t3;
            } else
            {
               t1 = u1 - v1; u1 = v1; v1 = t1;
               t3 = u3 - v3; u3 = v3; v3 = t3;
            }
         } else if (quot < (v3<<1))
         {
            t1 = u1 - (v1<<1); u1 = v1; v1 = t1;
            t3 = u3 - (v3<<1); u3 = v3; v3 = t3;
         } else
         {
            t1 = u1 - v1*3; u1 = v1; v1 = t1;
            t3 = u3 - v3*3; u3 = v3; v3 = t3;
         }
      } else
      {
         quot=u3/v3;
         t1 = u1 - v1*quot; u1 = v1; v1 = t1;
         t3 = u3 - v3*quot; u3 = v3; v3 = t3;
      }
   }

   if (u1<0) u1+=p;

   return u1;
}

/*=========================================================================
   Knuth_Schroeppel Multiplier:

   This is derived from Jason Papadopoulos's mpqs K-S method.  I believe it
   does a slightly better job than the K-S in FLINT 2.3, but that's debatable.
   An alternative would be to implement the method directly from Silverman 1987.

==========================================================================*/
/* small square-free numbers:  do { say $_ if moebius($_) != 0 } for 1..101 */
static const unsigned long multipliers[] = {
  1, 2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, 19, 21, 22, 23, 26, 29, 30, 31,
  33, 34, 35, 37, 38, 39, 41, 42, 43, 46, 47, 51, 53, 55, 57, 58, 59, 61,
  62, 65, 66, 67, 69, 70, 71, 73, 74, 77, 78, 79, 82, 83, 85, 86, 87, 89,
  91, 93, 94, 95, 97, 101};
#define NUMMULTS (sizeof(multipliers)/sizeof(unsigned long))
#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif

static unsigned long knuthSchroeppel(mpz_t n, unsigned long numPrimes)
{
  unsigned int i, j, best_mult, num_mults, knmod8;
  unsigned int maxprimes = (2*numPrimes <= 300) ? 2*numPrimes : 300;
  float best_score;
  float scores[NUMMULTS];
  mpz_t temp;

  mpz_init(temp);

  for (i = 0; i < NUMMULTS; i++) {
    unsigned int curr_mult = multipliers[i];
    scores[i] = 0.5 * logf((float)curr_mult);
    mpz_mul_ui(temp, n, curr_mult);
    knmod8 = mpz_mod_ui(temp, temp, 8);
    switch (knmod8) {
      case 1:  scores[i] -= 2 * M_LN2;  break;
      case 5:  scores[i] -= M_LN2;      break;
      case 3:
      case 7:  scores[i] -= 0.5 * M_LN2; break;
      default: break;
    }
  }
  num_mults = i;

  {
    PRIME_ITERATOR(iter);
    for (i = 1; i < maxprimes; i++) {
      unsigned int prime = prime_iterator_next(&iter);
      float contrib = logf((float)prime) / (float)(prime-1);
      unsigned int modp = mpz_mod_ui(temp, n, prime);

      for (j = 0; j < num_mults; j++) {
        unsigned int curr_mult = multipliers[j];
        mpz_set_ui(temp, modp);
        mpz_mul_ui(temp, temp, curr_mult);
        mpz_mod_ui(temp, temp, prime);
        if (mpz_sgn(temp) == 0) {
          scores[j] -= contrib;
        } else if (mpz_kronecker_ui(temp, prime) == 1) {
          scores[j] -= 2*contrib;
        }
      }
    }
    prime_iterator_destroy(&iter);
  }
  mpz_clear(temp);

  best_score = 1000.0;
  best_mult = 1;
  for (i = 0; i < num_mults; i++) {
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
static void initSieve(void)
{
    factorBase = 0;
    primeSizes = 0;
    flags = 0;
    sqrts = 0;
}
static void clearSieve(unsigned long numPrimes)
{
    if (factorBase) { free(factorBase);  factorBase = 0; }
    if (primeSizes) { free(primeSizes);  primeSizes = 0; }
    if (flags) { free(flags);  flags = 0; }
    if (sqrts) {
      unsigned int i;
      for (i = 0; i < numPrimes; i++) {
        mpz_clear(sqrts[i]);
      }
      free(sqrts);  sqrts = 0;
    }
}

/*========================================================================
   Compute Factor Base:

   Function: Computes primes p up to B for which n is a square mod p,
   allocates memory and stores them in an array pointed to by factorBase
   Returns: number of primes actually in the factor base

========================================================================*/
static void computeFactorBase(mpz_t n, unsigned long B,unsigned long multiplier)
{
  UV p;
  UV primesinbase = 0;
  PRIME_ITERATOR(iter);

  if (factorBase) { free(factorBase);  factorBase = 0; }
  factorBase = (unsigned int *) malloc( B * sizeof(unsigned int));

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
}

/*===========================================================================
   Compute Prime Sizes:

   Function: Computes the size in bits of each prime in the factor base
     allocates memory for an array, primeSizes, to store the sizes
     stores the size for each of the numPrimes primes in the array

===========================================================================*/
static void computeSizes(unsigned long numPrimes)
{
  unsigned long i;
  primeSizes = (unsigned char *) malloc( numPrimes * sizeof(unsigned char));
  for (i = 0; i < numPrimes; i++)
    primeSizes[i]=(unsigned char)floor(log(factorBase[i])/log(2.0)-FUDGE+0.5);
}

/*===========================================================================
   Tonelli-Shanks:

   Function: Performs Tonelli-Shanks on n mod every prime in the factor base
      allocates memory for the results to be stored in the array sqrts

===========================================================================*/
static void tonelliShanks(unsigned long numPrimes, mpz_t n)
{
  unsigned long i;
  mpz_t temp;
  sqrts = (mpz_t *) malloc( numPrimes * sizeof(mpz_t) );
  mpz_init_set_ui(sqrts[0], 0);
  mpz_init(temp);

  for (i = 1; i<numPrimes; i++) {
    mpz_init(sqrts[i]);
    mpz_set_ui(temp,factorBase[i]);
    sqrtmod(sqrts[i],n,temp);
  }
  mpz_clear(temp);
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
    matrix_t m,
    mpz_t * XArr,
    unsigned long * aind,
    int min,
    int s,
    int * exponents,
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
           while ((i < j*sizeof(unsigned long)) && (sieve[i] < threshold)) i++;
        } while (sieve[i] < threshold);

        if ((i<M) && (relsFound < relSought))
        {
           mpz_set_ui(temp,i+ctimesreps);
           mpz_sub_ui(temp,temp,Mdiv2);

           mpz_set(temp3,B);  //B
           mpz_addmul(temp3,A,temp);  //AX+B
           mpz_add(temp2,temp3,B);  //AX+2B
           mpz_mul(temp2,temp2,temp);  //AX^2+2BX
           mpz_add(res,temp2,C);  //AX^2+2BX+C

           bits=mpz_sizeinbase(res,2);
           bits-=errorbits;

           numfactors=0;
           extra = 0;

           memset(exponents, 0, firstprime * sizeof(int));

           if (factorBase[0] != 1 && mpz_divisible_ui_p(res, factorBase[0]))
           {
             extra += primeSizes[0];
             mpz_set_ui(temp,factorBase[0]);
             exponent = mpz_remove(res,res,temp);
             exponents[0] = exponent;
           }

           exponents[1] = 0;
           if (mpz_divisible_ui_p(res, factorBase[1]))
           {
             extra += primeSizes[1];
             mpz_set_ui(temp,factorBase[1]);
             exponent = mpz_remove(res,res,temp);
             exponents[1] = exponent;
           }

           for (k = 2; k < firstprime; k++)
           {
              modp=(i+ctimesreps)%factorBase[k];

              exponents[k] = 0;
              if (soln2[k] != (unsigned long)-1 )
              {
                 if ((modp==soln1[k]) || (modp==soln2[k]))
                 {
                    extra+=primeSizes[k];
                    mpz_set_ui(temp,factorBase[k]);
                    exponent = mpz_remove(res,res,temp);

#ifdef ERRORS
                    if (exponent==0) printf("Error!\n");
#endif
#ifdef RELS
                    if (exponent > 0) printf(" %ld",(long)factorBase[k]);
                    if (exponent > 1) printf("^%d",exponent);
#endif
                    exponents[k] = exponent;
                 }
              } else if (mpz_divisible_ui_p(res, factorBase[k]))
              {
                 extra += primeSizes[k];
                 mpz_set_ui(temp,factorBase[k]);
                 exponent = mpz_remove(res,res,temp);
#ifdef RELS
                 if (exponent > 0) printf(" %ld",(long)factorBase[k]);
                 if (exponent > 1) printf("^%d",exponent);
#endif
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

#ifdef ERRORS
                       if (exponent==0) printf("Error!\n");
#endif

#ifdef RELS
                       if (exponent > 0) printf(" %ld",(long)factorBase[k]);
                       if (exponent > 1) printf("^%d",exponent);
#endif
                       if (exponent)
                         for (ii = 0; ii < exponent; ii++)
                           set_relation(relations, relsFound, ++numfactors, k);
                       if (exponent & 1)
                         insertEntry(m,relsFound,k);
                    }
                 } else if (mpz_divisible_ui_p(res, factorBase[k]))
                 {
                    extra += primeSizes[k];
                    mpz_set_ui(temp,factorBase[k]);
                    exponent = mpz_remove(res,res,temp);
#ifdef RELS
                    if (exponent > 0) printf(" %ld",(long)factorBase[k]);
                    if (exponent > 1) printf("^%d",exponent);
#endif
                    for (ii = 0; ii < exponent; ii++)
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
#ifdef ERRORS
                       if (exponent==0) printf("Error!\n");
#endif
#ifdef RELS
                       if (exponent > 0) printf(" %ld",(long)factorBase[k]);
                       if (exponent > 1) printf("^%d",exponent);
#endif
                       if (exponent)
                         for (ii = 0; ii < exponent; ii++)
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
                    partials++;
                 }
                 clearRow(m,numPrimes,relsFound);
#ifdef RELS
                 gmp_printf(" %Zd\n",res);
#endif
              } else
              {
                 mpz_neg(res,res);
                 if (mpz_cmp_ui(res,1000)>0)
                 {
                    if (mpz_cmp_ui(res,largeprime)<0)
                    {
                       partials++;
                    }
                    clearRow(m,numPrimes,relsFound);
#ifdef RELS
                    gmp_printf(" %Zd\n",res);
#endif
                 } else
                 {
#ifdef RELS
                    printf("....R\n");
#endif
                    for (ii = 0; ii<firstprime; ii++)
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
                    if (relsFound%20==0) fprintf(stderr,"%ld relations, %ld partials.\n",(long)relsFound,(long)partials);
#endif
                 }
              }
           } else
           {
              clearRow(m,numPrimes,relsFound);
#ifdef RELS
              printf("\r                                                                    \r");
#endif

           }
           i++;

        } else if (relsFound >= relSought) i++;
     }
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
static void sieve2(unsigned long M, unsigned long numPrimes, unsigned char * sieve, const unsigned long * soln1, const unsigned long * soln2)
{
     unsigned int prime;
     unsigned char *end = sieve + M;

     memset(sieve, 0, M*sizeof(unsigned char));
     memset(flags, 0, numPrimes*sizeof(unsigned char));
     *end = 255; //sentinel to speed up sieve evaluators inner loop

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
   mainRoutine:

   Function: Generates the polynomials, initialises and calls the sieve,
             implementing cache blocking (breaking the sieve interval into
             small blocks for the small primes.

============================================================================*/
static int mainRoutine(
  unsigned long numPrimes,
  unsigned long Mdiv2,
  mpz_t n,
  mpz_t f,
  unsigned long multiplier)
{
    mpz_t A, B, C, D, Bdivp2, q, r, nsqrtdiv, temp, temp2, temp3, temp4;
    int result = 0;
    int i, j, l, s, fact, span, min;
    unsigned long u1, p, reps, numRelations, M;
    unsigned long curves = 0;
    unsigned long  * relations;
    unsigned short * primecount;
    unsigned char  * sieve;
    int            * exponents;
    unsigned long  * aind;
    unsigned long  * amodp;
    unsigned long  * Ainv;
    unsigned long  * soln1;
    unsigned long  * soln2;
    unsigned long ** Ainv2B;
    unsigned char ** offsets;
    unsigned char ** offsets2;
    mpz_t          * XArr;
    mpz_t          * Bterms;
    matrix_t m;


    exponents = (int *) malloc(firstprime * sizeof(int));
    if (exponents==NULL) croak("SIMPQS: Unable to allocate memory!\n");

    s = mpz_sizeinbase(n,2)/28+1;

    aind = (unsigned long*) calloc(sizeof(unsigned long),s);
    amodp = (unsigned long*) calloc(sizeof(unsigned long),s);
    Ainv = (unsigned long*) calloc(sizeof(unsigned long),numPrimes);
    soln1 = (unsigned long*) calloc(sizeof(unsigned long),numPrimes);
    soln2 = (unsigned long*) calloc(sizeof(unsigned long),numPrimes);
    Ainv2B = (unsigned long**) calloc(sizeof(unsigned long*),s);
    if (Ainv2B == 0) croak("SIMPQS: Unable to allocate memory!\n");
    Bterms = (mpz_t*) malloc( s * sizeof(mpz_t));
    if (Bterms == 0) croak("SIMPQS: Unable to allocate memory!\n");

    for (i=0; i<s; i++)
    {
       Ainv2B[i] = (unsigned long *) calloc(sizeof(unsigned long),numPrimes);
       if (Ainv2B[i] == 0) croak("SIMPQS: Unable to allocate memory!\n");
       mpz_init(Bterms[i]);
    }

    XArr = (mpz_t*) calloc( sizeof(mpz_t), relSought );

    m = constructMat(numPrimes, relSought);
    relsFound = 0;

    //one dword extra for sentinel to speed up sieve evaluation loop
    sieve = (unsigned char *) calloc(sizeof(unsigned char),Mdiv2*2 + sizeof(unsigned long));
    if (sieve==NULL)
      croak("SIMPQS: Unable to allocate memory for sieve!\n");

    flags = (unsigned char*) malloc(numPrimes * sizeof(unsigned char));

    offsets = (unsigned char **) malloc(secondprime * sizeof(unsigned char *));
    offsets2 = (unsigned char **)malloc(secondprime * sizeof(unsigned char *));

    relations = (unsigned long *) calloc(relSought * RELATIONS_PER_PRIME, sizeof(unsigned long));

    primecount = (unsigned short *) malloc(numPrimes * sizeof(unsigned short));

    mpz_init(A); mpz_init(B); mpz_init(C); mpz_init(D);
    mpz_init(Bdivp2); mpz_init(q); mpz_init(r); mpz_init(nsqrtdiv);
    mpz_init(temp); mpz_init(temp2); mpz_init(temp3); mpz_init(temp4);

//Compute min A_prime and A_span

    mpz_mul_ui(temp,n,2);
    mpz_sqrt(temp,temp);
    mpz_div_ui(nsqrtdiv,temp,Mdiv2);
    mpz_root(temp,nsqrtdiv,s);
    for (fact = 0; mpz_cmp_ui(temp,factorBase[fact])>=0; fact++);
    span = numPrimes/s/s/2;
    min=fact-span/2;
    while ((fact*fact)/min - min < span)
      min--;

#ifdef ADETAILS
    printf("s = %d, fact = %d, min = %d, span = %d\n",s,fact,min,span);
#endif

//Compute first polynomial and adjustments

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
           for (j=0;((j<i)&&(aind[j]!=fact));j++);
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
           mpz_mul(temp,temp,sqrts[aind[i]+min]);
           mpz_fdiv_r_ui(temp,temp,p);
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

// Count the number of polynomial curves used so far and compute the C coefficient of our polynomial

           curves++;

           mpz_mul(C,B,B);
           mpz_sub(C,C,n);
           mpz_divexact(C,C,A);

// Do the sieving and relation collection

           mpz_set_ui(temp,Mdiv2*2);
           mpz_fdiv_qr_ui(q,r,temp,CACHEBLOCKSIZE);
           M = mpz_get_ui(temp);

           /* set the solns1 and solns2 arrays */
           update_solns(1, numPrimes, soln1, soln2, polyadd, polycorr);
           /* Sieve [secondprime , numPrimes) */
           sieve2(M, numPrimes, sieve, soln1, soln2);
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
              soln1, soln2, m, XArr, aind,
              min, s, exponents,
              temp, temp2, temp3, temp4
           );
        }

#ifdef COUNT
        if (curves%20==0) printf("%ld curves.\n",(long)curves);
#endif
    }

#ifdef CURPARTS
    printf("%ld curves, %ld partials.\n",(long)curves,(long)partials);
#endif

#ifdef REPORT
    printf("Done with sieving!\n");
#endif

// Do the matrix algebra step

    numRelations = gaussReduce(m, numPrimes, relSought);
#ifdef REPORT
    printf("%ld relations in kernel.\n", numRelations);
#endif
    numRelations += 0;

// We want factors of n, not kn, so divide out by the multiplier

    mpz_div_ui(n,n,multiplier);

// Now do the "square root" and GCD steps hopefully obtaining factors of n
    for (l = (int)relSought-40; l < (int)relSought; l++)
    {
        unsigned int mat2offset = rightMatrixOffset(numPrimes);
        mpz_set_ui(temp,1);
        mpz_set_ui(temp2,1);
        memset(primecount,0,numPrimes*sizeof(unsigned short));
        for (i = 0; i< (int)numPrimes; i++)
        {
           if (getEntry(m,l,mat2offset+i))
           {
              mpz_mul(temp2,temp2,XArr[i]);
              for (j = 1; j <= (int)get_relation(relations, i, 0); j++)
              {
                 primecount[ get_relation(relations, i, j) ]++;
              }
           }
           if (i%30==0) mpz_mod(temp2,temp2,n);
        }
        for (j = 0; j < (int)numPrimes; j++)
        {
           mpz_set_ui(temp3,factorBase[j]);
           mpz_pow_ui(temp3,temp3,primecount[j]/2);
           mpz_mul(temp,temp,temp3);
           if (j%30==0) mpz_mod(temp,temp,n);
        }
        mpz_sub(temp,temp2,temp);
        mpz_gcd(temp,temp,n);
        /* Only non-trivial factors */
        if (mpz_cmp_ui(temp,1) && mpz_cmp(temp,n) && mpz_divisible_p(n,temp) ) {
          result = 1;
          mpz_set(f, temp);
#ifdef STANDALONE_SIMPQS
          gmp_printf("%Zd\n",temp);
#else
          break;
#endif
          }
    }

    destroyMat(m, relSought);
    free(primecount);
    free(relations);

    for (i = 0; i < s; i++) {
      free(Ainv2B[i]);
      mpz_clear(Bterms[i]);
    }
    free(exponents);  
    free(aind);
    free(amodp);
    free(Ainv);
    free(soln1);
    free(soln2);
    free(Ainv2B);
    free(Bterms);

    free(sieve);    sieve = 0;
    free(flags);    flags = 0;
    free(offsets);  offsets = 0;
    free(offsets2); offsets2 = 0;
    for (i = 0; i < (int)relSought; i++) {
      mpz_clear(XArr[i]);
    }
    free(XArr);

    mpz_clear(A);  mpz_clear(B);  mpz_clear(C);  mpz_clear(D);
    mpz_clear(q);  mpz_clear(r);
    mpz_clear(temp);  mpz_clear(temp2);  mpz_clear(temp3);  mpz_clear(temp4);
    mpz_clear(Bdivp2); mpz_clear(nsqrtdiv);

    return result;
}

int _GMP_simpqs(mpz_t n, mpz_t f)
{
  unsigned long numPrimes;
  unsigned long Mdiv2;
  unsigned long multiplier;
  int result;

  initSieve();
  decdigits = mpz_sizeinbase(n,10); /* often 1 too big */
  if (decdigits < MINDIG) {
    mpz_set(f, n);
    return 0;
  }

#ifdef REPORT
  gmp_printf("%Zd (%ld decimal digits)\n", n, decdigits);
#endif

  /* Get a preliminary number of primes, pick a multiplier, apply it */
  numPrimes = (decdigits <= 91) ? primesNo[decdigits-MINDIG] : 64000;
  multiplier = knuthSchroeppel(n, numPrimes);
  mpz_mul_ui(n, n, multiplier);
  decdigits = mpz_sizeinbase(n, 10);

  if (decdigits<=91) {
    numPrimes=primesNo[decdigits-MINDIG];

    Mdiv2 = sieveSize[decdigits-MINDIG]/SIEVEDIV;
    if (Mdiv2*2 < CACHEBLOCKSIZE) Mdiv2 = CACHEBLOCKSIZE/2;
    largeprime = largeprimes[decdigits-MINDIG];

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

  /* We probably need fewer than this */
  relSought = numPrimes;
  computeFactorBase(n, numPrimes, multiplier);

  computeSizes(numPrimes);

  TonelliInit();
  tonelliShanks(numPrimes,n);
  TonelliDestroy();

  result = mainRoutine(numPrimes, Mdiv2, n, f, multiplier);
  if (!result)
    mpz_set(f, n);

  clearSieve(numPrimes);
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
  int result;
  mpz_t n, f;

  mpz_init(n);
  mpz_init(f);

  printf("Input number to factor [ >=%d decimal digits]: ", MINDIG);
  gmp_scanf("%Zd",n);getchar();

  decdigits = mpz_sizeinbase(n,10);
  if (decdigits < MINDIG)
    croak("SIMPQS: Error in input or number has too few digits.\n");

  result = _GMP_simpqs(n, f);

  if (result) {
    gmp_printf("SUCCESS: %Zd\n", f);
  } else {
    gmp_printf("FAILURE\n");
  }

  return 0;
}
#endif
