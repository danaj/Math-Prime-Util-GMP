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
     - merge some of the 1.0 and 2.0 changes
     - make it work with smaller values
     - fix some memory errors
     - free memory all over
     - fewer globals
     - mpz_nextprime is slow, slow, slow.  Use prime_iterator.
   This does not use jasonp's block Lanczos code that v2.0 uses.  That code
   litters temporary files everywhere, but that's a solvable problem.
   TODO: Tune 25-40 digit parameters
   TODO: Portability
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


/* DANAJ: TODO: for this matrix code I'm doing a naughty thing:
 * assuming sizeof(int)==4.  The original code from Hart uses unsigned long,
 * which is a BSDism.  uint should be reasonably portable, though
 * some testers still have old non-C99 compilers.
 */
typedef unsigned int * row;  /* row of an F2 matrix */
typedef row * matrix;        /* matrix as a list of pointers to rows */


static unsigned int bitPattern[]  =
{
  0x80000000, 0x40000000, 0x20000000, 0x10000000,
  0x08000000, 0x04000000, 0x02000000, 0x01000000,
  0x00800000, 0x00400000, 0x00200000, 0x00100000,
  0x00080000, 0x00040000, 0x00020000, 0x00010000,
  0x00008000, 0x00004000, 0x00002000, 0x00001000,
  0x00000800, 0x00000400, 0x00000200, 0x00000100,
  0x00000080, 0x00000040, 0x00000020, 0x00000010,
  0x00000008, 0x00000004, 0x00000002, 0x00000001
};

static void insertEntry(matrix m, unsigned int i, unsigned int j)
{
     m[i][j / 32] |= bitPattern[j % 32];
}

static void xorEntry(matrix m, unsigned int i, unsigned int j)
{
     m[i][j / 32] ^= bitPattern[j % 32];
}

static unsigned int getEntry(matrix m, unsigned int i, unsigned int j)
{
     return m[i][j / 32] & bitPattern[j % 32];
}

static void swapRows(matrix m, unsigned int x, unsigned int y)
{
     row temp = m[x];   m[x] = m[y];  m[y] = temp;
}


static void clearRow(matrix m, unsigned int numcols, unsigned int row)
{
    int dwords = numcols/32;

    if (numcols%32) dwords++;
    memset( m[row], 0, dwords*4);
}

#if 0
static void displayRow(matrix m, unsigned int row, unsigned int numPrimes)
{
     int j;
     int length = numPrimes/32;
     if (numPrimes%32) length++;
     length*=64;

     printf("[");
     for (j = 0; j < length/2; j++)
       printf("%c", getEntry(m,row,j) ? '1' : '0');
     printf("  ");
     for (j = length/2; j < length; j++)
       printf("%c", getEntry(m,row,j) ? '1' : '0');
     printf("]\n");
}
#endif

static void xorRows(matrix m, unsigned int source, unsigned int dest, unsigned int length)
{
  unsigned int i, q, r;
  row x = m[dest];
  row y = m[source];

  r = length % 8; q = length - r;
  for (i=0; i < q; i += 8)
  {
    x[i] ^= y[i]; x[1+i] ^= y[1+i]; x[2+i] ^= y[2+i]; x[3+i] ^= y[3+i];
    x[4+i] ^= y[4+i]; x[5+i] ^= y[5+i]; x[6+i] ^= y[6+i]; x[7+i] ^= y[7+i];
  }
  switch (r)
  {
    case 7: x[i] ^= y[i]; i++;
    case 6: x[i] ^= y[i]; i++;
    case 5: x[i] ^= y[i]; i++;
    case 4: x[i] ^= y[i]; i++;
    case 3: x[i] ^= y[i]; i++;
    case 2: x[i] ^= y[i]; i++;
    case 1: x[i] ^= y[i]; i++;
  }
}

static matrix constructMat(unsigned int cols, unsigned int rows)
{
     unsigned int i;
     matrix m;
     unsigned int dwords = cols/32 + ((cols%32) ? 1 : 0);

     /* printf("construct mat %u %u with %u words\n", cols, rows, dwords); */
     /* If cols > rows, we write off the array */
     if (cols < rows) croak("SIMPQS:  cols %u > rows %u\n", cols, rows);
     m = (row *) calloc(sizeof(row),rows);
     if (m==NULL)
       croak("SIMPQS: Unable to allocate memory for matrix!\n");

     for (i = 0; i < rows; i++)
     {
         /* two matrices, side by side */
         m[i] = (row) calloc(2*dwords,sizeof(unsigned int));
     }
     if (m[rows-1]==NULL)
       croak("SIMPQS: Unable to allocate memory for matrix!\n");

     /* make second matrix identity, i.e. 1's along diagonal */
     for (i = 0; i < rows; i++)
     {
        insertEntry(m,i,i+32*dwords);
     }

     return m;
}

static void destroyMat(matrix m, unsigned int cols, unsigned int rows)
{
  unsigned int i;
  for (i = 0; i < rows; i++)
    free(m[i]);
  free(m);
}

/* gaussReduce:  Apply Gaussian elimination to a matrix. */
static unsigned int gaussReduce(matrix m, unsigned int numPrimes, unsigned int relSought, int extras)
{
     unsigned int rowUpto = 0;
     unsigned int irow;
     unsigned int length = (numPrimes+extras)/32;
     int icol;
     unsigned int checkRow;

     if (numPrimes%32) length++;
     length*=2;

     for (icol = numPrimes-1; icol >= 0; icol--)
     {
         irow = rowUpto;
         while ((irow < relSought)&&(getEntry(m,irow,icol)==0UL)) irow++;
         if (irow < relSought)
         {

             swapRows(m,rowUpto,irow);

             for (checkRow = rowUpto+1; checkRow<relSought; checkRow++)
             {
                 if (getEntry(m,checkRow,icol)!=0UL)
                 {
                    xorRows(m,rowUpto,checkRow,length);
                 }
             }

             rowUpto++;
         }
     }

     return rowUpto;
}

//===========================================================================
//Uncomment these for various pieces of debugging information

//#define COUNT    // Shows the number of relations generated and curves used during sieving
//#define RELPRINT     // Shows the actual factorizations of the relations
//#define ERRORS   // Error if relation should be divisible by a prime but isn't
//#define POLS     // Shows the polynomials being used by the sieve
//#define ADETAILS // Prints some details about the factors of the A coefficients of the polys
//#define LARGESTP // Prints the size of the largest factorbase prime
//#define CURPARTS // Prints the number of curves used and number of partial relations
//#define REPORT //report sieve size, multiplier and number of primes used

//===========================================================================
//Architecture dependent fudge factors

#if ULONG_MAX == 4294967295U
#define SIEVEMASK 0xC0C0C0C0U
#define MIDPRIME 1500
#define SIEVEDIV 1
#elif ULONG_MAX == 18446744073709551615U
#define SIEVEMASK 0xC0C0C0C0C0C0C0C0U
#define MIDPRIME       1500
#define SIEVEDIV 1
#else
 #error Cannot determine ulong size
#endif

#define CACHEBLOCKSIZE 64000 //Should be a little less than the L1/L2 cache size
                             //and a multiple of 64000
#define LARGEPRIME 6000000
#define MEDIUMPRIME    900
#define SECONDPRIME    6000 //This should be lower for slower machines
#define FUDGE          0.15 //Every program needs a mysterious fudge factor

#define MINDIG 30 //Will not factor numbers with less than this number of decimal digits

//===========================================================================
//Knuth-Schroeppel multipliers and a macro to count them

static const unsigned long multipliers[] = {1, 2, 3, 5, 7, 11, 13, 17, 19,
                                                23, 29, 31, 37, 41, 43};

#define NUMMULTS (sizeof(multipliers)/sizeof(unsigned long))

//===========================================================================
// Large prime cutoffs

static unsigned long largeprimes[] =
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
static unsigned long primesNo[] =
{
     1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, //30-39
     1500, 1500, 1600, 1700, 1750, 1800, 1900, 2000, 2050, 2100, //40-49
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
static unsigned long firstPrimes[] =
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
static unsigned long errorAmounts[] =
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
static unsigned long thresholds[] =
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
static unsigned long sieveSize[] =
{
     64000, 64000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, //30-39
     32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, //40-49
     32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, //50-59
     32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, 32000, //60-69
     32000, 32000, 64000, 64000, 64000, //70-74
     96000, 96000, 96000, 128000, 128000, //75-79
     160000, 160000, 160000, 160000, 160000, //80-84
     192000, 192000, 192000, 192000, 192000, //85-89
     192000, 192000 //90-91
};

//============================================================================
static long decdigits; //number of decimal digits of n
static unsigned long secondprime; //min(numprimes, SECONDPRIME) = cutoff for using flags when sieving
static unsigned long firstprime;  //first prime actually sieved with
static unsigned char errorbits;  //first prime actually sieved with
static unsigned char threshold;  //sieve threshold cutoff for smooth relations
static unsigned long midprime;
unsigned long largeprime;

static unsigned long * factorBase; //array of factor base primes
static unsigned long numPrimes; //number of primes in factor base
static unsigned long relSought; //number of relations sought, i.e. a "few" more than numPrimes
static unsigned char * primeSizes; //array of sizes in bits, of the factor base primes
static unsigned char * sieve; //actual array where sieving takes place
static unsigned char * * offsets; //offsets for each prime to use in sieve
static unsigned char * * offsets2; //offsets for each prime to use in sieve (we switch between these)
static unsigned long relsFound =0; //number of relations found so far
static unsigned long potrels = 0; //potential relations (including duplicates)
static unsigned char * flags; //flags used for speeding up sieving for large primes
static unsigned long partials = 0; //number of partial relations
static unsigned long Mdiv2; //size of sieving interval divide 2
static unsigned long mat2off; //offset of second square block in matrix

static mpz_t * sqrts; //square roots of n modulo each prime in the factor base

//static mpz_t n; //number to be factored



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

   Function: Find the best multiplier to use (allows 2 as a multiplier).
             The general idea is to find a multiplier k such that kn will
             be faster to factor. This is achieved by making kn a square
             modulo lots of small primes. These primes will then be factor
             base primes, and the more small factor base primes, the faster
             relations will accumulate, since they hit the sieving interval
             more often.

==========================================================================*/
static unsigned long knuthSchroeppel(mpz_t n)
{
    float bestFactor = -10.0f;
    unsigned long multiplier = 1;
    unsigned long nmod8;
    unsigned long multindex;
    float factors[NUMMULTS];
    float logpdivp;
    mpz_t r, mult;
    long kron;
    UV p;
    PRIME_ITERATOR(iter);

    mpz_init(r);
    mpz_init(mult);

    nmod8 = mpz_fdiv_r_ui(r,n,8);

    for (multindex = 0; multindex < NUMMULTS; multindex++)
    {
       long mod = nmod8*multipliers[multindex]%8;
       factors[multindex] = 0.34657359; // ln2/2
       if (mod == 1) factors[multindex] *= 4.0;
       if (mod == 5) factors[multindex] *= 2.0;
       factors[multindex] -= (log((float) multipliers[multindex]) / 2.0);
    }

    prime_iterator_setprime(&iter, 3);
    for (p = 3; p < 10000; p = prime_iterator_next(&iter)) {
          logpdivp = log((float)p) / p;
          kron = mpz_kronecker_ui(n,p);
          for (multindex = 0; multindex < NUMMULTS; multindex++)
          {
              mpz_set_ui(mult,multipliers[multindex]);
              switch (kron*mpz_kronecker_ui(mult,p))
              {
                 case 0:
                 {
                      factors[multindex] += logpdivp;
                 } break;
                 case 1:
                 {
                      factors[multindex] += 2.0*logpdivp;
                 } break;
                 default: break;
              }
          }
    }
    prime_iterator_destroy(&iter);

    for (multindex=0; multindex<NUMMULTS; multindex++)
    {
      if (factors[multindex] > bestFactor)
      {
        bestFactor = factors[multindex];
        multiplier = multipliers[multindex];
      }
    }

    mpz_clear(r);
    mpz_clear(mult);

    return multiplier;
}



/*========================================================================
   Initialize Quadratic Sieve:

   Function: Initialises the global gmp variables.

========================================================================*/
static void initSieve(void)
{
    factorBase = 0;
    primeSizes = 0;
    sieve = 0;
    offsets = 0;
    offsets2 = 0;
    flags = 0;
    sqrts = 0;
}
static void clearSieve(void)
{
    if (factorBase) { free(factorBase);  factorBase = 0; }
    if (primeSizes) { free(primeSizes);  primeSizes = 0; }
    if (sieve) { free(sieve);  sieve = 0; }
    if (offsets) { free(offsets);  offsets = 0; }
    if (offsets2) { free(offsets2);  offsets2 = 0; }
    if (flags) { free(flags);  flags = 0; }
    if (sqrts) {
      int i;
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
  factorBase = (unsigned long *) malloc( B * sizeof(unsigned long));

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
  primeSizes = (unsigned char *) calloc(sizeof(unsigned char),numPrimes);
  for (i = 0; i < numPrimes; i++)
    primeSizes[i]=(unsigned char)floor(log(factorBase[i])/log(2.0)-FUDGE+0.5);
}

/*===========================================================================
   Tonelli-Shanks:

   Function: Performs Tonelli-Shanks on n mod every prime in the factor base
      allocates memory for the results to be stored in the array sqrts

===========================================================================*/
static void tonelliShanks(unsigned long numPrimes,mpz_t n)
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
    unsigned long ** relations,
    unsigned long ctimesreps,
    unsigned long M,
    unsigned char * sieve,
    mpz_t A,
    mpz_t B,
    mpz_t C,
    unsigned long * soln1,
    unsigned long * soln2,
    int polyadd,
    unsigned long * polycorr,
    matrix m,
    mpz_t * XArr,
    unsigned long * aind,
    int min,
    int s,
    unsigned long multiplier,
    int * exponents)
{
     long i,j,ii;
     unsigned int k;
     unsigned int exponent, vv;
     unsigned char extra;
     unsigned int modp;
     unsigned long * sieve2;
     unsigned char bits;
     int numfactors;
     mpz_t temp, temp2, temp3, res;

     mpz_init(temp);  mpz_init(temp2);  mpz_init(temp3);  mpz_init(res);
     i = 0;
     j=0;
     sieve2 = (unsigned long *) sieve;
#ifdef POLS
     gmp_printf("%Zdx^2%+Zdx\n%+Zd\n",A,B,C);
#endif

     while ( (unsigned long)j < M/sizeof(unsigned long))
     {
        do
        {
           while (!(sieve2[j] & SIEVEMASK)) j++;
           i=j*sizeof(unsigned long);
           j++;
           while ((i<j*sizeof(unsigned long))&&(sieve[i] < threshold)) i++;
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
           if (factorBase[0]!=1)
           {
              mpz_set_ui(temp,factorBase[0]);
              exponent = mpz_remove(res,res,temp);
              if (exponent)
              {
                 extra+=primeSizes[0];
                 for (ii = 0; ii<exponent; ii++)
                   relations[relsFound][++numfactors] = 0;
              }
              if (exponent & 1) insertEntry(m,relsFound,0);
           }

           mpz_set_ui(temp,factorBase[1]);
           exponent = mpz_remove(res,res,temp);
           if (exponent)
           {
              extra+=primeSizes[1];
              for (ii = 0; ii<exponent; ii++)
                relations[relsFound][++numfactors] = 1;
           }
           if (exponent & 1) insertEntry(m,relsFound,1);

           for (k = 2; k<firstprime; k++)
           {
              modp=(i+ctimesreps)%factorBase[k];

              if (soln2[k]!=0xFFFFFFFFl)
              {
                 if ((modp==soln1[k]) || (modp==soln2[k]))
                 {
                    mpz_set_ui(temp,factorBase[k]);
                    exponent = mpz_remove(res,res,temp);

#ifdef ERRORS
                    if (exponent==0) printf("Error!\n");
#endif
                    extra+=primeSizes[k];
#ifdef RELS
                    if (exponent > 0) printf(" %ld",(long)factorBase[k]);
                    if (exponent > 1) printf("^%d",exponent);
#endif
                    exponents[k] = exponent;
                    if (exponent&1) insertEntry(m,relsFound,k);
                 } else exponents[k] = 0;
              } else
              {
                 mpz_set_ui(temp,factorBase[k]);
                 exponent = mpz_remove(res,res,temp);
                 if (exponent) extra+=primeSizes[k];
#ifdef RELS
                 if (exponent > 0) gmp_printf(" %Zd",factorBase[k]);
                 if (exponent > 1) printf("^%d",exponent);
#endif
                 exponents[k] = exponent;
                 if (exponent &1) insertEntry(m,relsFound,k);
              }
           }
           sieve[i]+=extra;
           if (sieve[i] >= bits)
           {
              vv=((unsigned char)1<<(i&7));
              for (k = firstprime; (k<secondprime)&&(extra<sieve[i]); k++)
              {
                 modp=(i+ctimesreps)%factorBase[k];
                 if (soln2[k]!=0xFFFFFFFFl)
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
                       {
                          for (ii = 0; ii<exponent; ii++)
                            relations[relsFound][++numfactors] = k;
                       }
                       if (exponent&1) insertEntry(m,relsFound,k);
                    }
                 } else
                 {
                    mpz_set_ui(temp,factorBase[k]);
                    exponent = mpz_remove(res,res,temp);
                    if (exponent) extra+=primeSizes[k];

#ifdef RELS
                    if (exponent > 0) printf(" %ld",(long)factorBase[k]);
                    if (exponent > 1) printf("^%d",exponent);
#endif
                    if (exponent)
                    {
                       for (ii = 0; ii<exponent; ii++)
                         relations[relsFound][++numfactors] = k;
                    }
                    if (exponent &1) insertEntry(m,relsFound,k);
                 }
              }


              for (k = secondprime; (k<numPrimes)&&(extra<sieve[i]); k++)
              {
                 if (flags[k]&vv)
                 {
                    modp=(i+ctimesreps)%factorBase[k];
                    if ((modp==soln1[k]) || (modp==soln2[k]))
                    {
                       mpz_set_ui(temp,factorBase[k]);
                       exponent = mpz_remove(res,res,temp);
#ifdef ERRORS
                       if (exponent==0) printf("Error!\n");
#endif
                       extra+=primeSizes[k];
#ifdef RELS
                       if (exponent > 0) printf(" %ld",(long)factorBase[k]);
                       if (exponent > 1) printf("^%d",exponent);
#endif
                       if (exponent)
                       {
                          for (ii = 0; ii<exponent; ii++)
                            relations[relsFound][++numfactors] = k;
                       }
                       if (exponent&1) insertEntry(m,relsFound,k);
                    }
                 }
              }

              for (ii =0; ii<s; ii++)
              {
                 xorEntry(m,relsFound,aind[ii]+min);
                 relations[relsFound][++numfactors] = aind[ii]+min;
              }

              if (mpz_cmp_ui(res,1000)>0)
              {
                 if (mpz_cmp_ui(res,LARGEPRIME)<0)
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
                    if (mpz_cmp_ui(res,LARGEPRIME)<0)
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
                    for (ii = 1; ii<firstprime; ii++)
                    {
                       int jj;
                       for (jj = 0; jj < exponents[ii]; jj++)
                         relations[relsFound][++numfactors] = ii;
                    }
                    relations[relsFound][0] = numfactors;

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
     mpz_clear(temp);  mpz_clear(temp2);  mpz_clear(temp3);  mpz_clear(res);
}


/*=============================================================================
   Sieve:

   Function: Allocates space for a sieve of M integers and sieves the interval
             starting at start

=============================================================================*/
static void sieveInterval(unsigned long M, unsigned long numPrimes, unsigned char * sieve, long last, long first, long polyadd, unsigned long * soln1, unsigned long * soln2, unsigned long * polycorr, unsigned char * * offsets, unsigned char * * offsets2)
{
     unsigned char currentprimesize;
     unsigned long currentprime;
     unsigned long prime;
     unsigned char * position2;
     unsigned char * position;
     long diff;
     unsigned char * end;
     unsigned long ptimes4;
     long correction;

     end = sieve+M;

     if (first)
     {
        for (prime=1; prime<firstprime; prime++)
        {
            if (soln2[prime] == 0xFFFFFFFF) continue;
            currentprime = factorBase[prime];
            correction = polyadd ? -polycorr[prime]+currentprime : polycorr[prime];
            soln1[prime]+=correction;
            while (soln1[prime]>=currentprime) soln1[prime]-=currentprime;
            soln2[prime]+=correction;
            while (soln2[prime]>=currentprime) soln2[prime]-=currentprime;
        }
     }

     for (prime=firstprime; prime<MEDIUMPRIME; prime++)
     {
        if (soln2[prime] == 0xFFFFFFFF) continue;
        currentprime = factorBase[prime];
        currentprimesize = primeSizes[prime];

        if (first)
        {
           correction = polyadd ? -polycorr[prime]+currentprime : polycorr[prime];
           soln1[prime]+=correction;
           while (soln1[prime]>=currentprime) soln1[prime]-=currentprime;
           soln2[prime]+=correction;
           while (soln2[prime]>=currentprime) soln2[prime]-=currentprime;

           position = sieve+soln1[prime];
           position2 = sieve+soln2[prime];
        } else
        {
           position = offsets[prime];
           position2 = offsets2[prime];
        }
        diff=position2-position;

        ptimes4 = currentprime*4;
        unsigned char * bound=end-ptimes4;
        while (bound - position > 0)
        {
	      (* position)+=currentprimesize,(* (position+diff))+=currentprimesize, position+=currentprime;
	      (* position)+=currentprimesize,(* (position+diff))+=currentprimesize, position+=currentprime;
	      (* position)+=currentprimesize,(* (position+diff))+=currentprimesize, position+=currentprime;
	      (* position)+=currentprimesize,(* (position+diff))+=currentprimesize, position+=currentprime;
        }
        while ((end - position > 0)&&(end - position - diff > 0))
        {
              (* position)+=currentprimesize,(* (position+diff))+=currentprimesize, position+=currentprime;

        }
        position2 = position+diff;
        if (end - position2 > 0)
        {
              (* position2)+=currentprimesize, position2+=currentprime;
        }
        if (end - position > 0)
        {
              (* position)+=currentprimesize, position+=currentprime;
        }

        if (!last)
        {
           offsets[prime] = position;
           offsets2[prime] = position2;
        }
     }

     for (prime=MEDIUMPRIME; prime<midprime; prime++)
     {
        currentprime = factorBase[prime];
        currentprimesize = primeSizes[prime];

        if (first)
        {
           correction = polyadd ? -polycorr[prime]+currentprime : polycorr[prime];
           soln1[prime]+=correction;
           while (soln1[prime]>=currentprime) soln1[prime]-=currentprime;
           soln2[prime]+=correction;
           while (soln2[prime]>=currentprime) soln2[prime]-=currentprime;

           position = sieve+soln1[prime];
           position2 = sieve+soln2[prime];
        } else
        {
           position = offsets[prime];
           position2 = offsets2[prime];
        }
        diff=position2-position;

        ptimes4 = 2*currentprime;
        unsigned char * bound=end-ptimes4;
        while (bound - position > 0)
        {
              (* position)+=currentprimesize,(* (position+diff))+=currentprimesize, position+=currentprime;
              (* position)+=currentprimesize,(* (position+diff))+=currentprimesize, position+=currentprime;
        }
        position2 = position+diff;
        while ((end - position > 0)&&(end - position2 > 0))
        {
              (* position)+=currentprimesize, position+=currentprime, (* position2)+=currentprimesize, position2+=currentprime;

        }

        if (end - position2 > 0)
        {
              (* position2)+=currentprimesize, position2+=currentprime;
        }
        if (end - position > 0)
        {
              (* position)+=currentprimesize, position+=currentprime;
        }
        if (!last)
        {
           offsets[prime] = position;
           offsets2[prime] = position2;
        }
     }

     return;
}

/*===========================================================================
   Sieve 2:

   Function: Second sieve for larger primes

=========================================================================== */
static void sieve2(unsigned long M, unsigned long numPrimes, unsigned char * sieve, long last, long first, long polyadd, unsigned long * soln1, unsigned long * soln2, unsigned long * polycorr, unsigned char * * offsets, unsigned char * * offsets2)
{
     unsigned char currentprimesize;
     unsigned long currentprime, prime;
     unsigned char * position2;
     unsigned char * position;
     unsigned char * end;
     long correction;

     memset(sieve,0,M*sizeof(unsigned char));
     memset(flags,0,numPrimes*sizeof(unsigned char));
     end = sieve+M;
     *end = 255; //sentinel to speed up sieve evaluators inner loop

     for (prime=midprime; prime<secondprime; prime++)
     {
        currentprime = factorBase[prime];
        currentprimesize = primeSizes[prime];

        correction = polyadd ? -polycorr[prime]+currentprime : polycorr[prime];
        soln1[prime]+=correction;
        while (soln1[prime]>=currentprime) soln1[prime]-=currentprime;
        soln2[prime]+=correction;
        while (soln2[prime]>=currentprime) soln2[prime]-=currentprime;

        position = sieve+soln1[prime];
        position2 = sieve+soln2[prime];

        while ((end - position > 0)&&(end - position2 > 0))
        {
             (* position)+=currentprimesize, position+=currentprime, (* position2)+=currentprimesize, position2+=currentprime;
        }

        if (end - position2 > 0)
        {
              (* position2)+=currentprimesize;
        }
        if (end - position > 0)
        {
              (* position)+=currentprimesize;
        }
     }

     for (prime=secondprime; prime<numPrimes; prime++)
     {
        currentprime = factorBase[prime];
        currentprimesize = primeSizes[prime];

        correction = polyadd ? -polycorr[prime]+currentprime : polycorr[prime];
        soln1[prime]+=correction;
        while (soln1[prime]>=currentprime) soln1[prime]-=currentprime;
        soln2[prime]+=correction;
        while (soln2[prime]>=currentprime) soln2[prime]-=currentprime;

        position = sieve+soln1[prime];
        position2 = sieve+soln2[prime];

        while (end - position > 0)
        {
              flags[prime]|=((unsigned char)1<<((position-sieve)&7)), (* position)+=currentprimesize, position+=currentprime;
        }

        while (end - position2 > 0)
        {
              flags[prime]|=((unsigned char)1<<((position2-sieve)&7)), (* position2)+=currentprimesize, position2+=currentprime;
        }
     }

     return;
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
static int mainRoutine(unsigned long Mdiv2, mpz_t n, mpz_t f, unsigned long multiplier)
{
    mpz_t A; mpz_init(A);
    mpz_t B; mpz_init(B);
    mpz_t C; mpz_init(C);
    mpz_t D; mpz_init(D);
    mpz_t Bdivp2; mpz_init(Bdivp2);
    mpz_t q; mpz_init(q);
    mpz_t r; mpz_init(r);
    mpz_t temp, temp2, temp3;

    int result = 0;

    /* gmp_randstate_t state;
    gmp_randinit_default(state); */

    unsigned long u1;

    int i, j, l, s, fact, span, min;
    unsigned long p;
    unsigned long reps;

    unsigned long curves = 0;

    unsigned long ** relations;
    int * primecount;

    int * exponents = (int *) calloc(firstprime,sizeof(int));
    if (exponents==NULL)
      croak("SIMPQS: Unable to allocate memory!\n");

    s = mpz_sizeinbase(n,2)/28+1;

    unsigned long * aind = (unsigned long*) calloc(sizeof(unsigned long),s);
    unsigned long * amodp = (unsigned long*) calloc(sizeof(unsigned long),s);
    unsigned long * Ainv = (unsigned long*) calloc(sizeof(unsigned long),numPrimes);
    unsigned long * soln1 = (unsigned long*) calloc(sizeof(unsigned long),numPrimes);
    unsigned long * soln2 = (unsigned long*) calloc(sizeof(unsigned long),numPrimes);
    unsigned long ** Ainv2B = (unsigned long**) calloc(sizeof(unsigned long*),s);
    if (Ainv2B == 0) croak("SIMPQS: Unable to allocate memory!\n");
    mpz_t* Bterms = (mpz_t*) malloc( s * sizeof(mpz_t));
    if (Bterms == 0) croak("SIMPQS: Unable to allocate memory!\n");

    for (i=0; i<s; i++)
    {
       Ainv2B[i] = (unsigned long *) calloc(sizeof(unsigned long),numPrimes);
       if (Ainv2B[i] == 0) croak("SIMPQS: Unable to allocate memory!\n");
       mpz_init(Bterms[i]);
    }

    matrix m;
    m = constructMat(numPrimes, relSought);

    relsFound = 0;

    mpz_t XArr[relSought];

    //one dword extra for sentinel to speed up sieve evaluation loop
    sieve = (unsigned char *) calloc(sizeof(unsigned char),Mdiv2*2 + sizeof(unsigned long));
    if (sieve==NULL)
      croak("SIMPQS: Unable to allocate memory for sieve!\n");

    flags = (unsigned char*) calloc(sizeof(unsigned char),numPrimes);

    offsets = (unsigned char * *)calloc(numPrimes,sizeof(unsigned char *));
    offsets2 = (unsigned char * *)calloc(numPrimes,sizeof(unsigned char *));

    relations = (unsigned long * *) calloc(numPrimes,sizeof(unsigned long *));
    for (i = 0; i < numPrimes; i++)
    {
       relations[i] = (unsigned long *) calloc(50, sizeof(unsigned long));
    }

    primecount = (int *) calloc(numPrimes, sizeof(int));

    mpz_init(temp);  mpz_init(temp2);  mpz_init(temp3);

//Compute min A_prime and A_span

    mpz_mul_ui(temp,n,2);
    mpz_sqrt(temp,temp);
    mpz_div_ui(temp2,temp,Mdiv2);
    mpz_root(temp,temp2,s);
    for (fact = 0; mpz_cmp_ui(temp,factorBase[fact])>=0; fact++);
    span = numPrimes/s/s/2;
    min=fact-span/2;

#ifdef ADETAILS
    printf("s = %d, fact = %d, min = %d, span = %d\n",s,fact,min,span);
#endif

//Compute first polynomial and adjustments

    while (relsFound < relSought)
    {
        int polyindex;
        mpz_set_ui(A,1);
#if 0
        for (i=0; i<s-1; )
        {
           int j,ran;
           mpz_set_ui(temp,span);
           mpz_urandomm(temp,state,temp);
           ran = mpz_get_ui(temp);
           for (j=0;((j<i)&&(aind[j]!=ran));j++);
           if (j==i)
           {
              aind[i] = ran;
              mpz_mul_ui(A,A,factorBase[ran+min]);
              i++;
           }
        }
#else
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
#endif
        mpz_div(temp,temp2,A);
        for (fact = 0; mpz_cmp_ui(temp,factorBase[fact])>=0; fact++);
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
           //mpz_fdiv_r_ui(temp,temp,p*p);
	   amodp[i] = mpz_fdiv_r_ui(temp,temp,p);

           //mpz_set_ui(temp3,p);
	   //mpz_invert(temp,temp,temp3);
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
        for (i=1; i<s; i++)
        {
           mpz_add(B,B,Bterms[i]);
        }

        for (i=0; i<numPrimes; i++)
        {
           p = factorBase[i];
           //mpz_set_ui(temp3,p);
	   //mpz_fdiv_r_ui(temp,A,(unsigned long64_t)p*(unsigned long64_t)p);
	   //mpz_fdiv_r_ui(temp,temp,p);
	   //mpz_fdiv_r(temp,A,temp3);
	   //mpz_invert(temp3,temp,temp3);
	   //Ainv[i] = mpz_get_ui(temp3);
	   Ainv[i] = modinverse(mpz_fdiv_r_ui(temp,A,p),p);

           for (j=0; j<s; j++)
           {
              mpz_fdiv_r_ui(temp,Bterms[j],p);
	      //mpz_fdiv_r_ui(temp,temp,p);
              mpz_mul_ui(temp,temp,2*Ainv[i]);
              Ainv2B[j][i] = mpz_fdiv_r_ui(temp,temp,p);
           }

           mpz_fdiv_r_ui(temp,B,p);
	   //mpz_fdiv_r_ui(temp,temp,p);
           mpz_sub(temp,sqrts[i],temp);
           mpz_add_ui(temp,temp,p);
           mpz_mul_ui(temp,temp,Ainv[i]);
           mpz_add_ui(temp,temp,Mdiv2);
           soln1[i] = mpz_fdiv_r_ui(temp,temp,p);
           //mpz_set(temp,sqrts[i]);
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

           int index;
           for (j=0; j<s; j++)
           {
              index = aind[j]+min;
              p = factorBase[index];
              //mpz_set_ui(temp,p*p);
              //mpz_mul(temp,temp,temp);
              mpz_fdiv_r_ui(D,n,p*p);
              mpz_fdiv_r_ui(Bdivp2,B,p*p);
              mpz_mul_ui(temp,Bdivp2,amodp[j]);
              mpz_realloc2(temp3,64);
	      //mpz_set_ui(temp3,p);
	      //mpz_fdiv_r_ui(temp,temp,p*p);
	      mpz_fdiv_r_ui(temp,temp,p);
	      //mpz_invert(temp3,temp,temp3);
	      //u1=mpz_get_ui(temp3);
	      u1 = modinverse(mpz_fdiv_r_ui(temp,temp,p),p);
              mpz_mul(temp,Bdivp2,Bdivp2);
              mpz_sub(temp,temp,D);
              mpz_neg(temp,temp);
              mpz_div_ui(temp,temp,p);
              mpz_mul_ui(temp,temp,u1);
              mpz_add_ui(temp,temp,Mdiv2);
              mpz_add_ui(temp,temp,p);
              //mpz_fdiv_r_ui(temp,temp,p*p);
	      soln1[index]=mpz_fdiv_r_ui(temp,temp,p);
              //soln1[index] = mpz_get_ui(temp);
              soln2[index]=0xFFFFFFFFl;
           }

// Count the number of polynomial curves used so far and compute the C coefficient of our polynomial

           curves++;

           mpz_mul(C,B,B);
           mpz_sub(C,C,n);
           mpz_divexact(C,C,A);

// Do the sieving and relation collection

           mpz_set_ui(temp,Mdiv2*2);
           mpz_fdiv_qr_ui(q,r,temp,CACHEBLOCKSIZE);
           sieve2(mpz_get_ui(temp),numPrimes,sieve,1,1,polyadd,soln1,soln2,polycorr,offsets,offsets2);
           sieveInterval(CACHEBLOCKSIZE,numPrimes,sieve,0,1,polyadd,soln1,soln2,polycorr,offsets,offsets2);
           if (mpz_cmp_ui(q,1)>0)
           {
              unsigned long maxreps = mpz_get_ui(q)-1;
              for (reps = 1; reps < maxreps; reps++)
              {
                 sieveInterval(CACHEBLOCKSIZE,numPrimes,sieve+CACHEBLOCKSIZE*reps,0,0,polyadd,soln1,soln2,polycorr,offsets,offsets2);
              }
              if (mpz_cmp_ui(r,0)==0)
              {
                 sieveInterval(CACHEBLOCKSIZE,numPrimes,sieve+CACHEBLOCKSIZE*reps,1,0,polyadd,soln1,soln2,polycorr,offsets,offsets2);
              } else
              {
                 sieveInterval(CACHEBLOCKSIZE,numPrimes,sieve+CACHEBLOCKSIZE*reps,0,0,polyadd,soln1,soln2,polycorr,offsets,offsets2);
                 reps++;
                 sieveInterval(mpz_get_ui(r),numPrimes,sieve+CACHEBLOCKSIZE*reps,1,0,polyadd,soln1,soln2,polycorr,offsets,offsets2);
              }
           }

           evaluateSieve(relations,0,mpz_get_ui(temp),sieve,A,B,C,soln1, soln2, polyadd, polycorr, m,XArr,aind,min,s,multiplier,exponents);
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

    unsigned long kernel = gaussReduce(m, numPrimes, relSought, 0);
#ifdef REPORT
    printf("%ld relations in kernel.\n",kernel);
#endif


// We want factors of n, not kn, so divide out by the multiplier

    mpz_div_ui(n,n,multiplier);

// Now do the "square root" and GCD steps hopefully obtaining factors of n
    for (l = relSought-40;l<relSought;l++)
    {
        mpz_set_ui(temp,1);
        mpz_set_ui(temp2,1);
        mat2off = numPrimes/32;
        if (numPrimes%32) mat2off++;
        mat2off*=32;
        memset(primecount,0,numPrimes*sizeof(int));
        for (i = 0; i<(numPrimes); i++)
        {
           if (getEntry(m,l,i+mat2off))
           {
              mpz_mul(temp2,temp2,XArr[i]);
              for (j=1; j<=relations[i][0]; j++)
              {
                 primecount[relations[i][j]]++;
              }
           }
           if (i%30==0) mpz_mod(temp2,temp2,n);
        }
        for (j = 0; j<numPrimes; j++)
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

    destroyMat(m, numPrimes, relSought);
    free(primecount);
    for (i = 0; i < numPrimes; i++)
      free(relations[i]);
    free(relations);

    for (i=0; i<s; i++) {
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
    for (i = 0; i < relSought; i++) {
      mpz_clear(XArr[i]);
    }

    mpz_clear(A);  mpz_clear(B);  mpz_clear(C);  mpz_clear(D);
    mpz_clear(q);  mpz_clear(r);
    mpz_clear(temp);  mpz_clear(temp2);  mpz_clear(temp3);
    mpz_clear(Bdivp2);

    return result;
}

int _GMP_simpqs(mpz_t n, mpz_t f)
{
  unsigned long multiplier;
  int result;

  initSieve();
  decdigits = mpz_sizeinbase(n,10); /* often 1 too big */
  if (decdigits < MINDIG) {
    mpz_set(f, n);
    return 0;
  }

  multiplier = knuthSchroeppel(n);
  mpz_mul_ui(n, n, multiplier);

  if (decdigits<=91) {
    numPrimes=primesNo[decdigits-MINDIG];

    Mdiv2 = sieveSize[decdigits-MINDIG]/SIEVEDIV;
    if (Mdiv2*2 < CACHEBLOCKSIZE) Mdiv2 = CACHEBLOCKSIZE/2;
    largeprime = largeprimes[decdigits-MINDIG];

    secondprime = (numPrimes < SECONDPRIME) ? numPrimes : SECONDPRIME;
    midprime    = (numPrimes < MIDPRIME)    ? numPrimes : MIDPRIME;

    firstprime = firstPrimes[decdigits-MINDIG];
    errorbits = errorAmounts[decdigits-MINDIG];
    threshold = thresholds[decdigits-MINDIG];
  } else {
    numPrimes = 64000;
    Mdiv2 = 192000/SIEVEDIV;
    largeprime = numPrimes*10*decdigits;

    secondprime = SECONDPRIME;
    midprime = MIDPRIME;
    firstprime = 30;
    errorbits = decdigits/4 + 2;
    threshold = 43+(7*decdigits)/10;
  }

#ifdef REPORT
  gmp_printf("%Zd (%ld decimal digits)\n", n, decdigits);
  printf("Using multiplier: %ld\n",(long)multiplier);
  printf("%ld primes in factor base.\n",(long)numPrimes);
  printf("Sieving interval M = %ld\n",(long)Mdiv2*2);
  printf("Large prime cutoff = factorBase[%ld]\n",largeprime);
#endif

  /* We probably need fewer than this */
  relSought = numPrimes;
  computeFactorBase(n, numPrimes, multiplier);

  computeSizes(numPrimes);

  TonelliInit();
  tonelliShanks(numPrimes,n);
  TonelliDestroy();

  result = mainRoutine(Mdiv2, n, f, multiplier);
  if (!result)
    mpz_set(f, n);

  clearSieve();
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
