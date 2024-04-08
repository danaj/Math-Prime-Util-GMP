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

   Modifications made in 2024 by Hugo van der Sanden:
     - combining partial relations (this is huge for large inputs)
     - jasonp's block Lanczos routine
     - much less memory use for large inputs
     - further little changes / optimizations

   There may be further improvements to find in msieve, particularly to
   the block Lanczos code: https://sourceforge.net/projects/msieve/

   To compile standalone:
   gcc -O2 -DSTANDALONE_SIMPQS -DSTANDALONE simpqs.c utility.c \
     bls75.c ecm.c ecpp.c factor.c gmp_main.c isaac.c lucas_seq.c pbrent63.c \
     primality.c prime_iterator.c random_prime.c real.c rootmod.c squfof126.c \      tinyqs.c -lgmp -lm

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

/*============================================================================
    Block Lanczos code Copyright 2006 Jason Papadopoulos

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

===============================================================================

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may
benefit from your work.
                       --jasonp@boo.net 9/8/06
--------------------------------------------------------------------*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <gmp.h>

/* random failure rate of block_lanczos() appears to be about 4%, but
 * a matrix that fails once appears to have a higher probability of failing
 * more times - so there's a risk that some matrix will fail every time. */
#define BL_MAX_FAIL 100

#include "ptypes.h"
#ifdef STANDALONE_SIMPQS
# include "gmp_main.h"
# define UV_MAX ULONG_MAX
# define UVCONST(x) ((unsigned long)x##UL)
# define New(id, mem, size, type)  mem = (type*) malloc((size)*sizeof(type))
# define Newz(id, mem, size, type) mem = (type*) calloc(size, sizeof(type))
# define Safefree(mem)             free((void*)mem)
#else
# include "simpqs.h"
#endif

#include "prime_iterator.h"
#include "utility.h"
#include "rootmod.h"

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

#ifdef ERRORS
# define CHECK_EXPONENT(exponent, k) \
    if (exponent == 0)               \
      printf("Error with prime %u!\n", factorBase[k]);
#else
# define CHECK_EXPONENT(exponent, k)
#endif
#ifdef RELPRINT
# define PRINT_FB(exponent, k) do { \
    if (exponent > 0)               \
      printf(" %u", factorBase[k]); \
    if (exponent > 1)               \
      printf("^%u", exponent);      \
  } while (0)
#else
# define PRINT_FB(exponent, k)
#endif

/* Architecture dependent fudge factors
 */
#if ULONG_MAX == 4294967295UL
# define SIEVEMASK 0xC0C0C0C0UL
# define SIEVEDIV 1
#elif ULONG_MAX == 18446744073709551615UL
# define SIEVEMASK 0xC0C0C0C0C0C0C0C0UL
# define SIEVEDIV 1
#else
# error Cannot determine ulong size
#endif

/* Should be a little less than the L1/L2 cache size and a multiple of 64000 */
#define CACHEBLOCKSIZE 64000
#define MIDPRIME 1500
/* Make lower for slower machines */
#define SECONDPRIME 6000
/* Used for tweaking the bit size calculation for factorBase primes */
#define SIZE_FUDGE 0.15

/* Will not factor numbers with less than this number of decimal digits */
#define MINDIG 30

/* Large prime cutoffs, in thousands */
static const unsigned int largeprimes[] = {
   100,  100,  125,   125,   150,   150,   175,   175,   200,   200, /* 30-39 */
   250,  300,  370,   440,   510,   580,   650,   720,   790,   860, /* 40-49 */
   930, 1000, 1700,  2400,  3100,  3800,  4500,  5200,  5900,  6600, /* 50-59 */
  7300, 8000, 8900, 10000, 11300, 12800, 14500, 16300, 18100, 20000, /* 60-69 */
   22000,  24000,  27000,  32000,  39000, /* 70-74 */
   53000,  65000,  75000,  87000, 100000, /* 75-79 */
  114000, 130000, 150000, 172000, 195000, /* 80-84 */
  220000, 250000, 300000, 350000, 400000, /* 85-89 */
  450000, 500000 /* 90-91 */
};

/* Number of primes to use in factor base */
static const unsigned int primesNo[] = {
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

/* First prime actually sieved for */
static const unsigned int firstPrimes[] = {
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

/* Logs of primes are rounded and errors accumulate
 * This specifies how great an error to allow */
static const unsigned int errorAmounts[] = {
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

/* Threshold the sieve value must exceed to be considered for smoothness.
 * SIEVEMASK implies at least 64. */
static const unsigned int thresholds[] = {
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

/* Size of sieve to use divided by 2
 * Probably optimal if chosen to be a multiple of 32000, though other sizes
 * are supported */
static const unsigned int sieveSize[] = {
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

static unsigned int secondprime;  /* cutoff for using flags when sieving */
static unsigned int firstprime;   /* first prime actually sieved with */
static unsigned char errorbits;
static unsigned char threshold;   /* sieve threshold cutoff for smooth relations */
static unsigned int midprime;
static unsigned int largeprime;
static unsigned int *factorBase;  /* array of factor base primes */
static unsigned char *primeSizes; /* array of sizes in bits of fb primes */

/* lanczos.c */

typedef struct {
  unsigned long *data;        /* The list of occupied rows in this column */
  unsigned long weight;       /* Number of nonzero entries in this column */
  unsigned long orig;         /* Original relation number */
} la_col_t;

/* insertColEntry: insert an entry into a column of the matrix, reallocating
 * the space for the column if necessary.
 * Note: we keep the data sized to a multiple of 16 entries, and assume it
 * must be reallocated any time we find there is exactly a multiple of 16
 * present before insert. This means we may undergo unnecessary churn if
 * xorColEntry() frequently removes entries across a 16-boundary.
 * Probably better to track size and grow exponentially in any case.
*/
static inline void insertColEntry(
  la_col_t* colarray, unsigned long colNum, unsigned long entry
) {
  unsigned long* temp;

  if ((colarray[colNum].weight & 0x0f) == 0) {
    /* need more space */
    temp = colarray[colNum].data;
    colarray[colNum].data = (unsigned long*)malloc(
        (colarray[colNum].weight + 16) * sizeof(unsigned long));
    for (long i = 0; i < colarray[colNum].weight; ++i)
      colarray[colNum].data[i] = temp[i];
    free(temp);
  }

  colarray[colNum].data[colarray[colNum].weight] = entry;
  ++colarray[colNum].weight;
  colarray[colNum].orig = colNum;
}

/* xorColEntry: xor entry corresponding to a prime dividing A */
static inline void xorColEntry(
  la_col_t* colarray, unsigned long colNum, unsigned long entry
) {
  for (long i = 0; i < colarray[colNum].weight; ++i)
    if (colarray[colNum].data[i] == entry) {
      /* found, so remove it */
      for (unsigned long j = i; j < colarray[colNum].weight - 1; ++j)
        colarray[colNum].data[j] = colarray[colNum].data[j + 1];
      --colarray[colNum].weight;
      return;
    }
  /* not found, so insert it */
  insertColEntry(colarray, colNum, entry);
}

/* clearCol: Function: clear a column */
static inline void clearCol(la_col_t* colarray, unsigned long colNum) {
   colarray[colNum].weight = 0;
}

#define NUM_EXTRA_RELATIONS 64

#define BIT(x) (((u_int64_t)1) << (x))

static const u_int64_t bitmask[64] = {
  BIT( 0), BIT( 1), BIT( 2), BIT( 3), BIT( 4), BIT( 5), BIT( 6), BIT( 7),
  BIT( 8), BIT( 9), BIT(10), BIT(11), BIT(12), BIT(13), BIT(14), BIT(15),
  BIT(16), BIT(17), BIT(18), BIT(19), BIT(20), BIT(21), BIT(22), BIT(23),
  BIT(24), BIT(25), BIT(26), BIT(27), BIT(28), BIT(29), BIT(30), BIT(31),
  BIT(32), BIT(33), BIT(34), BIT(35), BIT(36), BIT(37), BIT(38), BIT(39),
  BIT(40), BIT(41), BIT(42), BIT(43), BIT(44), BIT(45), BIT(46), BIT(47),
  BIT(48), BIT(49), BIT(50), BIT(51), BIT(52), BIT(53), BIT(54), BIT(55),
  BIT(56), BIT(57), BIT(58), BIT(59), BIT(60), BIT(61), BIT(62), BIT(63),
};

/* Returns true if the entry with indices i,l is 1 in the
 * supplied 64xN matrix. This is used to read the nullspace
 * vectors which are output by the Lanczos routine
 */
u_int64_t getNullEntry(u_int64_t *nullrows, long i, long l) {
  return nullrows[i] & bitmask[l];
}

/* Poor man's random number generator. It satisfies no
 * particularly good randomness properties, but is good
 * enough for this application
 */
u_long random32(void) {
  static unsigned long randval = 4035456057U;
  randval = ((u_int64_t)randval * 1025416097U + 286824428U)
          % (u_int64_t)4294967291U;
  return randval;
}

/* Returns the maximum of two unsigned long's */
unsigned long max(unsigned long a, unsigned long b) {
   return (a < b) ? b : a;
}

/* Perform light filtering on the nrows x ncols matrix specified by cols[].
 * The processing here is limited to deleting columns that contain a singleton
 * row, then resizing the matrix to have a few more columns than rows.
 * Because deleting a column reduces the counts in several different rows,
 * the process must iterate to convergence.
 *
 * Note that this step is not intended to make the Lanczos iteration run
 * any faster (though it will); it's just that if we don't go to this trouble
 * then there are factorizations for which the matrix step will fail outright.
 */
void reduce_matrix(unsigned long *nrows, unsigned long *ncols, la_col_t *cols) {
  unsigned long r, c, i, j, k;
  unsigned long passes;
  unsigned long *counts;
  unsigned long reduced_rows;
  unsigned long reduced_cols;

  /* count the number of nonzero entries in each row */
  counts = (unsigned long *)calloc((size_t)*nrows, sizeof(unsigned long));
  for (i = 0; i < *ncols; ++i)
    for (j = 0; j < cols[i].weight; ++j)
      ++counts[cols[i].data[j]];

  reduced_rows = *nrows;
  reduced_cols = *ncols;
  passes = 0;

  do {
    r = reduced_rows;

    /* remove any columns that contain the only entry in one or more rows,
     * then update the row counts to reflect the missing column.
     * Iterate until no more columns can be deleted */
    do {
      c = reduced_cols;
      for (i = j = 0; i < reduced_cols; ++i) {
        la_col_t *col = cols + i;
        for (k = 0; k < col->weight; ++k)
          if (counts[col->data[k]] < 2)
            break;

        if (k < col->weight) {
          for (k = 0; k < col->weight; ++k)
            --counts[col->data[k]];
          free(col->data);
          col->data = NULL;
        } else {
          if (j != i) {
            /* j lags i, will never have data to free */
            cols[j] = cols[i];
            cols[i].data = NULL;
          }
          ++j;
        }
      }
      reduced_cols = j;
    } while (c != reduced_cols);

    /* count the number of rows that contain a nonzero entry */
    for (i = reduced_rows = 0; i < *nrows; ++i)
      if (counts[i])
        ++reduced_rows;

    /* Because deleting a column reduces the weight of many rows, the
     * number of nonzero rows may be much less than the number of columns.
     * Delete more columns until the matrix has the correct aspect ratio.
     * Columns at the end of cols[] are the heaviest, so delete those
     * (and update the row counts again) */
    if (reduced_cols > reduced_rows + NUM_EXTRA_RELATIONS) {
      for (i = reduced_rows + NUM_EXTRA_RELATIONS;
        i < reduced_cols; ++i
      ) {
        la_col_t *col = &cols[i];
        for (j = 0; j < col->weight; ++j)
          --counts[col->data[j]];
        free(col->data);
        col->data = NULL;
      }
      reduced_cols = reduced_rows + NUM_EXTRA_RELATIONS;
    }

    /* if any columns were deleted in the previous step, then the matrix
     * is less dense and more columns can be deleted; iterate until no
     * further deletions are possible */
    ++passes;
  } while (r != reduced_rows);

  if (get_verbose_level() > 3)
    printf("reduced to %lu x %lu in %lu passes\n",
        reduced_rows, reduced_cols, passes);

  free(counts);

  /* Record the final matrix size. Note that we can't touch nrows because
   * all the column data (and the sieving relations that produced it) would
   * have to be updated */
  *ncols = reduced_cols;
}

/* c[][] := x[][] * y[][], where all operands are 64 x 64 (i.e. contain 64
 * words of 64 bits each). The result may overwrite a or b.
 */
static void mul_64x64_64x64(u_int64_t *a, u_int64_t *b, u_int64_t *c) {
  u_int64_t ai, accum;
  u_int64_t tmp[64];
  unsigned long i, j;

  for (i = 0; i < 64; ++i) {
    j = 0;
    accum = 0;
    ai = a[i];
    while (ai) {
      if (ai & 1)
        accum ^= b[j];
      ai >>= 1;
      ++j;
    }
    tmp[i] = accum;
  }
  memcpy(c, tmp, sizeof(tmp));
}

/* Let x[][] be a 64 x 64 matrix in GF(2), represented as 64 words of 64 bits
 * each. Let c[][] be an 8 x 256 matrix of 64-bit words. This code fills c[][]
 * with a bunch of "partial matrix multiplies". For 0 <= i < 256, the j_th row
 * of c[][] contains the matrix product:
 *
 *   (i << (8 * j)) * x[][]
 *
 * where the quantity in parentheses is considered a 1 x 64 vector of elements
 * in GF(2). The resulting table can dramatically speed up matrix multiplies
 * by x[][].
 */
static void precompute_Nx64_64x64(u_int64_t *x, u_int64_t *c) {
  u_int64_t accum;
  unsigned long i, j, k, index;

  for (j = 0; j < 8; ++j) {
    for (i = 0; i < 256; ++i) {
      k = 0;
      index = i;
      accum = 0;
      while (index) {
        if (index & 1)
          accum ^= x[k];
        index >>= 1;
        ++k;
      }
      c[i] = accum;
    }
    x += 8;
    c += 256;
  }
}

/* Let v[][] be an n x 64 matrix with elements in GF(2), represented as an
 * array of n 64-bit words. Let c[][] be an 8 x 256 scratch matrix of 64-bit
 * words. This code multiplies v[][] by the 64x64 matrix x[][], then XORs
 * the n x 64 result into y[][].
 */
static void mul_Nx64_64x64_acc(
  u_int64_t *v, u_int64_t *x, u_int64_t *c, u_int64_t *y, unsigned long n
) {
  unsigned long i;
  u_int64_t word;

  precompute_Nx64_64x64(x, c);
  for (i = 0; i < n; ++i) {
    word = v[i];
    y[i] ^=  c[0 * 256 + ((word >>  0) & 0xff)]
           ^ c[1 * 256 + ((word >>  8) & 0xff)]
           ^ c[2 * 256 + ((word >> 16) & 0xff)]
           ^ c[3 * 256 + ((word >> 24) & 0xff)]
           ^ c[4 * 256 + ((word >> 32) & 0xff)]
           ^ c[5 * 256 + ((word >> 40) & 0xff)]
           ^ c[6 * 256 + ((word >> 48) & 0xff)]
           ^ c[7 * 256 + ((word >> 56)       )];
  }
}

/* Let x and y be n x 64 matrices. This routine computes the 64 x 64 matrix
 * xy[][] given by transpose(x) * y. c[][] is a 256 x 8 scratch matrix of
 * 64-bit words.
*/
static void mul_64xN_Nx64(
  u_int64_t *x, u_int64_t *y, u_int64_t *c, u_int64_t *xy, unsigned long n
) {
  unsigned long i;

  memset(c, 0, 256 * 8 * sizeof(u_int64_t));
  memset(xy, 0, 64 * sizeof(u_int64_t));

  for (i = 0; i < n; ++i) {
    u_int64_t xi = x[i];
    u_int64_t yi = y[i];
    c[0 * 256 + ( xi        & 0xff)] ^= yi;
    c[1 * 256 + ((xi >>  8) & 0xff)] ^= yi;
    c[2 * 256 + ((xi >> 16) & 0xff)] ^= yi;
    c[3 * 256 + ((xi >> 24) & 0xff)] ^= yi;
    c[4 * 256 + ((xi >> 32) & 0xff)] ^= yi;
    c[5 * 256 + ((xi >> 40) & 0xff)] ^= yi;
    c[6 * 256 + ((xi >> 48) & 0xff)] ^= yi;
    c[7 * 256 + ((xi >> 56)       )] ^= yi;
  }

  for (i = 0; i < 8; ++i) {
    unsigned long j;
    u_int64_t a0, a1, a2, a3, a4, a5, a6, a7;

    a0 = a1 = a2 = a3 = 0;
    a4 = a5 = a6 = a7 = 0;
    for (j = 0; j < 256; ++j) {
      if ((j >> i) & 1) {
        a0 ^= c[0 * 256 + j];
        a1 ^= c[1 * 256 + j];
        a2 ^= c[2 * 256 + j];
        a3 ^= c[3 * 256 + j];
        a4 ^= c[4 * 256 + j];
        a5 ^= c[5 * 256 + j];
        a6 ^= c[6 * 256 + j];
        a7 ^= c[7 * 256 + j];
      }
    }

    xy[ 0] = a0; xy[ 8] = a1; xy[16] = a2; xy[24] = a3;
    xy[32] = a4; xy[40] = a5; xy[48] = a6; xy[56] = a7;
    ++xy;
  }
}

/* Given a 64x64 matrix t[][] (i.e. sixty-four 64-bit words) and a list of
 * 'last_dim' column indices enumerated in last_s[]:
 *  - find a submatrix of t that is invertible
 *  - invert it and copy to w[][]
 *  - enumerate in s[] the columns represented in w[][]
 */
static unsigned long find_nonsingular_sub(
  u_int64_t *t, unsigned long *s, unsigned long *last_s,
  unsigned long last_dim, u_int64_t *w
) {
  unsigned long i, j;
  unsigned long dim;
  unsigned long cols[64];
  u_int64_t M[64][2];
  u_int64_t mask, *row_i, *row_j;
  u_int64_t m0, m1;

  /* M = [t | I] for I the 64x64 identity matrix */
  for (i = 0; i < 64; ++i) {
    M[i][0] = t[i];
    M[i][1] = bitmask[i];
  }

  /* put the column indices from last_s[] into the back of cols[], and copy
   * to the beginning of cols[] any column indices not in last_s[] */
  mask = 0;
  for (i = 0; i < last_dim; ++i) {
    cols[63 - i] = last_s[i];
    mask |= bitmask[last_s[i]];
  }
  for (i = j = 0; i < 64; ++i)
    if (!(mask & bitmask[i]))
      cols[j++] = i;

  /* compute the inverse of t[][] */
  for (i = dim = 0; i < 64; ++i) {
    /* find the next pivot row and put in row i */
    mask = bitmask[cols[i]];
    row_i = M[cols[i]];

    for (j = i; j < 64; ++j) {
      row_j = M[cols[j]];
      if (row_j[0] & mask) {
        m0 = row_j[0];
        m1 = row_j[1];
        row_j[0] = row_i[0];
        row_j[1] = row_i[1];
        row_i[0] = m0;
        row_i[1] = m1;
        break;
      }
    }

    /* if a pivot row was found, eliminate the pivot column from all
     * other rows */
    if (j < 64) {
      for (j = 0; j < 64; ++j) {
        row_j = M[cols[j]];
        if ((row_i != row_j) && (row_j[0] & mask)) {
          row_j[0] ^= row_i[0];
          row_j[1] ^= row_i[1];
        }
      }

      /* add the pivot column to the list of accepted columns */
      s[dim++] = cols[i];
      continue;
    }

    /* otherwise, use the right-hand half of M[] to compensate for
     * the absence of a pivot column */
    for (j = i; j < 64; ++j) {
      row_j = M[cols[j]];
      if (row_j[1] & mask) {
        m0 = row_j[0];
        m1 = row_j[1];
        row_j[0] = row_i[0];
        row_j[1] = row_i[1];
        row_i[0] = m0;
        row_i[1] = m1;
        break;
      }
    }

    if (j == 64) {
      printf("lanczos error: submatrix is not invertible\n");
      return 0;
    }

    /* eliminate the pivot column from the other rows of the inverse */
    for (j = 0; j < 64; ++j) {
      row_j = M[cols[j]];
      if (row_i != row_j && (row_j[1] & mask)) {
        row_j[0] ^= row_i[0];
        row_j[1] ^= row_i[1];
      }
    }

    /* wipe out the pivot row */
    row_i[0] = row_i[1] = 0;
  }

  /* the right-hand half of M[] is the desired inverse */
  for (i = 0; i < 64; ++i)
    w[i] = M[i][1];

  /* The block Lanczos recurrence depends on all columns of t[][] appearing
   * in s[] and/or last_s[]. Verify that condition here */
  mask = 0;
  for (i = 0; i < dim; ++i)
    mask |= bitmask[s[i]];
  for (i = 0; i < last_dim; ++i)
    mask |= bitmask[last_s[i]];
  if (mask != (u_int64_t)(-1)) {
    if (get_verbose_level() > 3)
      printf("lanczos error: not all columns used\n");
    return 0;
  }
  return dim;
}

/* Multiply the vector x[] by the matrix A (stored columnwise) and put
 * the result in b[]. vsize refers to the number of u_int64_t's allocated for
 * x[] and b[]; vsize is probably different from ncols.
 */
void mul_MxN_Nx64(
  unsigned long vsize, unsigned long dense_rows, unsigned long ncols,
  la_col_t *A, u_int64_t *x, u_int64_t *b
) {
  unsigned long i, j;

  memset(b, 0, vsize * sizeof(u_int64_t));

  for (i = 0; i < ncols; ++i) {
    la_col_t *col = &A[i];
    unsigned long *row_entries = col->data;
    u_int64_t tmp = x[i];

    for (j = 0; j < col->weight; ++j)
      b[row_entries[j]] ^= tmp;
  }

  if (dense_rows)
    for (i = 0; i < ncols; ++i) {
      la_col_t *col = &A[i];
      unsigned long *row_entries = col->data + col->weight;
      u_int64_t tmp = x[i];

      for (j = 0; j < dense_rows; ++j)
        if (row_entries[j / 32] & ((unsigned long)1 << (j % 32)))
          b[j] ^= tmp;
    }
}

/* Multiply the vector x[] by the transpose of the matrix A and put the result
 * in b[]. Since A is stored by columns, this is just a matrix-vector product.
 */
void mul_trans_MxN_Nx64(
  unsigned long dense_rows, unsigned long ncols,
  la_col_t *A, u_int64_t *x, u_int64_t *b
) {
  unsigned long i, j;

  for (i = 0; i < ncols; ++i) {
    la_col_t *col = &A[i];
    unsigned long *row_entries = col->data;
    u_int64_t accum = 0;

    for (j = 0; j < col->weight; ++j)
      accum ^= x[row_entries[j]];
    b[i] = accum;
  }

  if (dense_rows)
    for (i = 0; i < ncols; ++i) {
      la_col_t *col = &A[i];
      unsigned long *row_entries = col->data + col->weight;
      u_int64_t accum = b[i];

      for (j = 0; j < dense_rows; ++j)
        if (row_entries[j / 32] & ((unsigned long)1 << (j % 32)))
          accum ^= x[j];
      b[i] = accum;
    }
}

/* Hideously inefficent routine to transpose a vector v[] of 64-bit words
 * into a 2-D array trans[][] of 64-bit words.
 */
static void transpose_vector(
  unsigned long ncols, u_int64_t *v, u_int64_t **trans
) {
  unsigned long i, j;
  unsigned long col;
  u_int64_t mask, word;

  for (i = 0; i < ncols; ++i) {
    col = i / 64;
    mask = bitmask[i % 64];
    word = v[i];
    j = 0;
    while (word) {
      if (word & 1)
        trans[j][col] |= mask;
      word = word >> 1;
      ++j;
    }
  }
}

/* Once the block Lanczos iteration has finished, x[] and v[] will contain
 * mostly nullspace vectors between them, as well as possibly some columns
 * that are linear combinations of nullspace vectors. Given vectors ax[]
 * and av[] that are the result of multiplying x[] and v[] by the matrix,
 * this routine will use Gauss elimination on the columns of [ax | av]
 * to find all of the linearly dependent columns. The column operations
 * needed to accomplish this are mirrored in [x | v] and the columns that
 * are independent are skipped. Finally, the dependent columns are copied
 * back into x[] and represent the nullspace vector output of the block
 * Lanczos code.
 *
 * v[] and av[] can be NULL, in which case the elimination process assumes
 * 64 dependencies instead of 128.
 */
void combine_cols(
  unsigned long ncols, u_int64_t *x, u_int64_t *v,
  u_int64_t *ax, u_int64_t *av
) {
  unsigned long i, j, k, bitpos, col, col_words, num_deps;
  u_int64_t mask;
  u_int64_t *matrix[128], *amatrix[128], *tmp;

  num_deps = 128;
  if (v == NULL || av == NULL)
    num_deps = 64;
  col_words = (ncols + 63) / 64;

  for (i = 0; i < num_deps; ++i) {
    matrix[i] = (u_int64_t *)calloc((size_t)col_words, sizeof(u_int64_t));
    amatrix[i] = (u_int64_t *)calloc((size_t)col_words, sizeof(u_int64_t));
  }

  /* operations on columns can more conveniently become operations on rows
   * if all the vectors are first transposed */
  transpose_vector(ncols, x, matrix);
  transpose_vector(ncols, ax, amatrix);
  if (num_deps == 128) {
    transpose_vector(ncols, v, matrix + 64);
    transpose_vector(ncols, av, amatrix + 64);
  }

  /* Keep eliminating rows until the unprocessed part of amatrix[][] is
   * all zero. The rows where this happens correspond to linearly dependent
   * vectors in the nullspace */
  for (i = bitpos = 0; i < num_deps && bitpos < ncols; ++bitpos) {
    /* find the next pivot row */
    mask = bitmask[bitpos % 64];
    col = bitpos / 64;
    for (j = i; j < num_deps; ++j)
      if (amatrix[j][col] & mask) {
        tmp = matrix[i];
        matrix[i] = matrix[j];
        matrix[j] = tmp;
        tmp = amatrix[i];
        amatrix[i] = amatrix[j];
        amatrix[j] = tmp;
        break;
      }
    if (j == num_deps)
      continue;

    /* a pivot was found; eliminate it from the remaining rows */
    for (++j; j < num_deps; ++j)
      if (amatrix[j][col] & mask) {
        /* Note that the entire row, *not* just the nonzero part of it,
         * must be eliminated; this is because the corresponding (dense)
         * row of matrix[][] must have the same operation applied */
        for (k = 0; k < col_words; ++k) {
          amatrix[j][k] ^= amatrix[i][k];
          matrix[j][k] ^= matrix[i][k];
        }
      }
    ++i;
  }

  /* transpose rows i to 64 back into x[] */
  for (j = 0; j < ncols; ++j) {
    u_int64_t word = 0;

    col = j / 64;
    mask = bitmask[j % 64];

    for (k = i; k < 64; ++k)
      if (matrix[k][col] & mask)
        word |= bitmask[k];
    x[j] = word;
  }

  for (i = 0; i < num_deps; ++i) {
    free(matrix[i]);
    free(amatrix[i]);
  }
}

/* Solve Bx = 0 for some nonzero x; the computed solution, containing up to
 * 64 of these nullspace vectors, is returned.
 */
u_int64_t * block_lanczos(
  unsigned long nrows, unsigned long dense_rows, unsigned long ncols,
  la_col_t *B
) {
  u_int64_t *vnext, *v[3], *x, *v0;
  u_int64_t *winv[3];
  u_int64_t *vt_a_v[2], *vt_a2_v[2];
  u_int64_t *scratch;
  u_int64_t *d, *e, *f, *f2;
  u_int64_t *tmp;
  unsigned long s[2][64];
  unsigned long i, iter;
  unsigned long n = ncols;
  unsigned long dim0, dim1;
  u_int64_t mask0, mask1;
  unsigned long vsize;

  /* allocate all of the size-n variables. Note that because B has been
   * preprocessed to ignore singleton rows, the number of rows may really
   * be less than nrows and may be greater than ncols. vsize is the maximum
   * of these two numbers. */
  vsize = max(nrows, ncols);
  v[0] = (u_int64_t *)malloc(vsize * sizeof(u_int64_t));
  v[1] = (u_int64_t *)malloc(vsize * sizeof(u_int64_t));
  v[2] = (u_int64_t *)malloc(vsize * sizeof(u_int64_t));
  vnext = (u_int64_t *)malloc(vsize * sizeof(u_int64_t));
  x = (u_int64_t *)malloc(vsize * sizeof(u_int64_t));
  v0 = (u_int64_t *)malloc(vsize * sizeof(u_int64_t));
  scratch = (u_int64_t *)malloc(max(vsize, 256 * 8) * sizeof(u_int64_t));

  /* allocate all the 64x64 variables */
  winv[0] = (u_int64_t *)malloc(64 * sizeof(u_int64_t));
  winv[1] = (u_int64_t *)malloc(64 * sizeof(u_int64_t));
  winv[2] = (u_int64_t *)malloc(64 * sizeof(u_int64_t));
  vt_a_v[0] = (u_int64_t *)malloc(64 * sizeof(u_int64_t));
  vt_a_v[1] = (u_int64_t *)malloc(64 * sizeof(u_int64_t));
  vt_a2_v[0] = (u_int64_t *)malloc(64 * sizeof(u_int64_t));
  vt_a2_v[1] = (u_int64_t *)malloc(64 * sizeof(u_int64_t));
  d = (u_int64_t *)malloc(64 * sizeof(u_int64_t));
  e = (u_int64_t *)malloc(64 * sizeof(u_int64_t));
  f = (u_int64_t *)malloc(64 * sizeof(u_int64_t));
  f2 = (u_int64_t *)malloc(64 * sizeof(u_int64_t));

  /* The iterations computes v[0], vt_a_v[0], vt_a2_v[0], s[0] and winv[0].
   * Subscripts larger than zero represent past versions of these
   * quantities, which start off empty (except for the past version of s[],
   * which contains all the column indices */
  memset(v[1], 0, vsize * sizeof(u_int64_t));
  memset(v[2], 0, vsize * sizeof(u_int64_t));
  for (i = 0; i < 64; ++i) {
    s[1][i] = i;
    vt_a_v[1][i] = 0;
    vt_a2_v[1][i] = 0;
    winv[1][i] = 0;
    winv[2][i] = 0;
  }
  dim0 = 0;
  dim1 = 64;
  mask1 = (u_int64_t)-1;
  iter = 0;

  /* The computed solution 'x' starts off random, and v[0] starts off
   * as B*x. This initial copy of v[0] must be saved off separately */
  for (i = 0; i < n; ++i)
    v[0][i] = (u_int64_t)(random32()) << 32
        | (u_int64_t)(random32());

  memcpy(x, v[0], vsize * sizeof(u_int64_t));
  mul_MxN_Nx64(vsize, dense_rows, ncols, B, v[0], scratch);
  mul_trans_MxN_Nx64(dense_rows, ncols, B, scratch, v[0]);
  memcpy(v0, v[0], vsize * sizeof(u_int64_t));

  /* perform the iteration */
  while (1) {
    ++iter;

    /* multiply the current v[0] by a symmetrized version of B,
     * or B'B (apostrophe means transpose). Use "A" to refer to B'B  */
    mul_MxN_Nx64(vsize, dense_rows, ncols, B, v[0], scratch);
    mul_trans_MxN_Nx64(dense_rows, ncols, B, scratch, vnext);

    /* compute v0'*A*v0 and (A*v0)'(A*v0) */
    mul_64xN_Nx64(v[0], vnext, scratch, vt_a_v[0], n);
    mul_64xN_Nx64(vnext, vnext, scratch, vt_a2_v[0], n);

    /* if the former is orthogonal to itself, then the iteration has
     * finished */
    for (i = 0; i < 64; ++i)
      if (vt_a_v[0][i] != 0)
        break;
    if (i == 64)
      break;

    /* Find the size-'dim0' nonsingular submatrix of v0'*A*v0, invert it,
     * and list the column indices present in the submatrix */
    dim0 = find_nonsingular_sub(vt_a_v[0], s[0], s[1], dim1, winv[0]);
    if (dim0 == 0)
      break;

    /* mask0 contains one set bit for every column that participates
     * in the inverted submatrix computed above */
    mask0 = 0;
    for (i = 0; i < dim0; ++i)
      mask0 |= bitmask[s[0][i]];

    /* compute d */
    for (i = 0; i < 64; ++i)
      d[i] = (vt_a2_v[0][i] & mask0) ^ vt_a_v[0][i];
    mul_64x64_64x64(winv[0], d, d);
    for (i = 0; i < 64; ++i)
      d[i] = d[i] ^ bitmask[i];

    /* compute e */
    mul_64x64_64x64(winv[1], vt_a_v[0], e);
    for (i = 0; i < 64; ++i)
      e[i] = e[i] & mask0;

    /* compute f */
    mul_64x64_64x64(vt_a_v[1], winv[1], f);
    for (i = 0; i < 64; ++i)
      f[i] = f[i] ^ bitmask[i];
    mul_64x64_64x64(winv[2], f, f);
    for (i = 0; i < 64; ++i)
      f2[i] = ((vt_a2_v[1][i] & mask1) ^ vt_a_v[1][i]) & mask0;
    mul_64x64_64x64(f, f2, f);

    /* compute the next v */
    for (i = 0; i < n; ++i)
      vnext[i] = vnext[i] & mask0;
    mul_Nx64_64x64_acc(v[0], d, scratch, vnext, n);
    mul_Nx64_64x64_acc(v[1], e, scratch, vnext, n);
    mul_Nx64_64x64_acc(v[2], f, scratch, vnext, n);

    /* update the computed solution 'x' */
    mul_64xN_Nx64(v[0], v0, scratch, d, n);
    mul_64x64_64x64(winv[0], d, d);
    mul_Nx64_64x64_acc(v[0], d, scratch, x, n);

    /* rotate all the variables */
    tmp = v[2];
    v[2] = v[1];
    v[1] = v[0];
    v[0] = vnext;
    vnext = tmp;
    tmp = winv[2];
    winv[2] = winv[1];
    winv[1] = winv[0];
    winv[0] = tmp;
    tmp = vt_a_v[1]; vt_a_v[1] = vt_a_v[0]; vt_a_v[0] = tmp;
    tmp = vt_a2_v[1]; vt_a2_v[1] = vt_a2_v[0]; vt_a2_v[0] = tmp;
    memcpy(s[1], s[0], 64 * sizeof(unsigned long));
    mask1 = mask0;
    dim1 = dim0;
  }

  if (get_verbose_level() > 3)
    printf("lanczos halted after %lu iterations\n", iter);

  /* free unneeded storage */
  free(vnext);
  free(scratch);
  free(v0);
  free(vt_a_v[0]);
  free(vt_a_v[1]);
  free(vt_a2_v[0]);
  free(vt_a2_v[1]);
  free(winv[0]);
  free(winv[1]);
  free(winv[2]);
  free(d);
  free(e);
  free(f);
  free(f2);

  /* if a recoverable failure occurred, start everything over again */
  if (dim0 == 0) {
    if (get_verbose_level() > 3)
      printf("linear algebra failed; retrying...\n");
    free(x);
    free(v[0]);
    free(v[1]);
    free(v[2]);
    return NULL;
  }

  /* convert the output of the iteration to an actual collection of
   * nullspace vectors */
  mul_MxN_Nx64(vsize, dense_rows, ncols, B, x, v[1]);
  mul_MxN_Nx64(vsize, dense_rows, ncols, B, v[0], v[2]);
  combine_cols(ncols, x, v[0], v[1], v[2]);

  /* verify that these really are linear dependencies of B */
  mul_MxN_Nx64(vsize, dense_rows, ncols, B, x, v[0]);
  for (i = 0; i < ncols; ++i)
    if (v[0][i] != 0)
      break;
  if (i < ncols) {
    printf("lanczos error: dependencies don't work %lu\n",i);
    abort();
  }

  free(v[0]);
  free(v[1]);
  free(v[2]);
  return x;
}

/* small prime power */
typedef struct {
  unsigned int p;
  unsigned int e;
} spp_t;

/* dynamic array of prime powers */
typedef struct {
  unsigned int size;
  unsigned int count;
  spp_t *fact;
} fact_t;

/* relation */
typedef struct {
  fact_t f;
  mpz_t X;
  /* full relation if Q = 0, else Q is the partial */
  unsigned long Q;
} rel_t;

/* dynamic array of relations */
typedef struct {
  unsigned int size;
  unsigned int count;
  rel_t **r;
} arel_t;

rel_t *new_rel(void) {
  rel_t *r = calloc(1, sizeof(rel_t));
  mpz_init(r->X);
  return r;
}

void resize_fact(fact_t *f, unsigned int size) {
  if (size > f->size) {
    f->fact = realloc(f->fact, size * sizeof(spp_t));
    f->size = size;
  }
}

void add_factor(rel_t *r, unsigned int p, unsigned int e) {
  unsigned int i = r->f.count;
  if (i >= r->f.size)
    resize_fact(&r->f, r->f.size ? (r->f.size * 3 / 2) : 64);
  r->f.fact[i].p = p;
  r->f.fact[i].e = e;
  ++r->f.count;
}

void resize_arel(arel_t *ar, unsigned int size) {
  if (size > ar->size) {
    ar->r = realloc(ar->r, size * sizeof(rel_t *));
    ar->size = size;
  }
}
void reset_factors(rel_t *r) {
  r->f.count = 0;
}
void reset_arel(arel_t *ar) {
  ar->count = 0;
}
arel_t *new_arel(void) {
  return (arel_t *)calloc(1, sizeof(arel_t));
}
void free_rel(rel_t *r) {
  mpz_clear(r->X);
  free(r->f.fact);
  free(r);
}
/* note, does not free rel_t relations referenced */
void steal_arel(arel_t *ar) {
  free(ar->r);
  free(ar);
}
void free_arel(arel_t *ar) {
  for (unsigned int i = 0; i < ar->count; ++i)
    free_rel(ar->r[i]);
  steal_arel(ar);
}

int fact_cmp(const void *va, const void *vb) {
  spp_t *a = (spp_t *)va;
  spp_t *b = (spp_t *)vb;
  return (a->p < b->p) ? -1 : (a->p == b->p) ? 0 : 1;
}
int arel_cmp(const void *va, const void *vb) {
  rel_t **a = (rel_t **)va;
  rel_t **b = (rel_t **)vb;
  unsigned long aQ = (*a)->Q;
  unsigned long bQ = (*b)->Q;
  return (aQ == bQ) ? mpz_cmp((*a)->X, (*b)->X) : (aQ < bQ) ? -1 : 1;
}

void sort_fact(fact_t *f) {
  qsort(f->fact, f->count, sizeof(spp_t), fact_cmp);
}
void sort_rel(arel_t *ar) {
  qsort(ar->r, ar->count, sizeof(rel_t *), arel_cmp);
}

void save_rel(arel_t *ar, rel_t *r) {
  unsigned int i = ar->count;
#if RELPRINT
  if (r->Q)
    gmp_printf("save partial Q=%lu X=%Zd\n", r->Q, r->X);
  else
    gmp_printf("save full X=%Zd\n", r->X);
#endif
  if (i >= ar->size)
    resize_arel(ar, ar->size ? (ar->size * 3 / 2) : 64);
  ar->r[ar->count++] = r;
  return;
}

unsigned int merge_full(arel_t *ara, arel_t *arb) {
  arel_t *arc = new_arel();
  resize_arel(arc, ara->count + arb->count);
  unsigned int ia = 0, ib = 0, ic = 0;
  rel_t *next_ra = ara->count ? ara->r[ia++] : NULL;
  rel_t *next_rb = arb->count ? arb->r[ib++] : NULL;
  mpz_t *lastp = NULL;

  while (next_ra && next_rb) {
    if (mpz_cmp(next_ra->X, next_rb->X) <= 0) {
      if (lastp && mpz_cmp(next_ra->X, *lastp) == 0) {
        free_rel(next_ra);
      } else {
        arc->r[ic++] = next_ra;
        lastp = &next_ra->X;
      }
      next_ra = (ia < ara->count) ? ara->r[ia++] : NULL;
    } else {
      if (lastp && mpz_cmp(next_rb->X, *lastp) == 0) {
        free_rel(next_rb);
      } else {
        arc->r[ic++] = next_rb;
        lastp = &next_rb->X;
      }
      next_rb = (ib < arb->count) ? arb->r[ib++] : NULL;
    }
  }
  while (next_ra) {
    arc->r[ic++] = next_ra;
    lastp = &next_ra->X;
    next_ra = (ia < ara->count) ? ara->r[ia++] : NULL;
  }
  while (next_rb) {
    if (lastp && mpz_cmp(next_rb->X, *lastp) == 0) {
      free_rel(next_rb);
    } else {
      arc->r[ic++] = next_rb;
      lastp = &next_rb->X;
    }
    next_rb = (ib < arb->count) ? arb->r[ib++] : NULL;
  }
  arc->count = ic;
  /* now swap contents of ara and arc, and free the old ara */
  arel_t tmp = *ara;
  *ara = *arc;
  *arc = tmp;
  steal_arel(arc);
  return ara->count;
}

unsigned int combine_partial(arel_t *comb, rel_t *ra, rel_t *rb,
  unsigned long numPrimes, mpz_t n
) {
  /* same partial, nothing to do */
  if (mpz_cmpabs(ra->X, rb->X) == 0)
    return 0;

  mpz_t X;
  mpz_init_set_ui(X, ra->Q);
  if (!mpz_invert(X, X, n)) {
    /* We have found a factor. It could be N when N is quite small;
     * or we might just have found a divisor by sheer luck. */
    mpz_gcd_ui(X, n, ra->Q);
    if (mpz_cmp(X, n) == 0) {
      /* it was N, nothing to see here */
      mpz_clear(X);
      return 0;
    }
    /* FIXME: record this factor and let the world know */
    gmp_fprintf(stderr, "Early factor found %Zd (not necessarily prime)\n", X);
    mpz_clear(X);
    return 0;
  }
  mpz_mul(X, X, ra->X);
  mpz_mul(X, X, rb->X);
  mpz_mod(X, X, n);

  /* prefer the smaller of (X, n - X) */
  mpz_t X2;
  mpz_init(X2);
  mpz_sub(X2, n, X);
  if (mpz_cmp(X2, X) < 0)
    mpz_set(X, X2);
  mpz_clear(X2);

  rel_t *r = new_rel();
  mpz_set(r->X, X);
  fact_t *fa = &ra->f;
  fact_t *fb = &rb->f;
  resize_fact(&r->f, fa->count + fb->count);

  unsigned int ia = 0, ib = 0;
  while (ia < fa->count && ib < fb->count) {
    if (fa->fact[ia].p == fb->fact[ib].p) {
      add_factor(r, fa->fact[ia].p, fa->fact[ia].e + fb->fact[ib].e);
      ++ia;
      ++ib;
    } else if (fa->fact[ia].p < fb->fact[ib].p) {
      add_factor(r, fa->fact[ia].p, fa->fact[ia].e);
      ++ia;
    } else {
      add_factor(r, fb->fact[ib].p, fb->fact[ib].e);
      ++ib;
    }
  }
  while (ia < ra->f.count) {
    add_factor(r, fa->fact[ia].p, fa->fact[ia].e);
    ++ia;
  }
  while (ib < rb->f.count) {
    add_factor(r, fb->fact[ib].p, fb->fact[ib].e);
    ++ib;
  }
  save_rel(comb, r);
  mpz_clear(X);
  return 1;
}

unsigned int merge_partial(arel_t *ara, arel_t *arb, arel_t *comb,
  unsigned long numPrimes, mpz_t n
) {
  arel_t *arc = new_arel();
  /* make room for worst case, if no full relations are extracted */
  resize_arel(arc, ara->count + arb->count);
  unsigned int i, ia = 0, ib = 0, ic = 0;
  rel_t *ra = ara->count ? ara->r[ia++] : NULL;
  rel_t *rb = arb->count ? arb->r[ib++] : NULL;

  while (rb) {
    if (ra) {
      if (ra->Q < rb->Q) {
        arc->r[ic++] = ra;
        ra = (ia < ara->count) ? ara->r[ia++] : NULL;
        continue;
      }
      if (ra->Q == rb->Q) {
        while (rb && ra->Q == rb->Q) {
          combine_partial(comb, ra, rb, numPrimes, n);
          free_rel(rb);
          rb = (ib < arb->count) ? arb->r[ib++] : NULL;
        }
        arc->r[ic++] = ra;
        ra = (ia < ara->count) ? ara->r[ia++] : NULL;
        continue;
      }
    }
    /* rb is next, but may combine with rb' */
    if (ib < arb->count && rb->Q == arb->r[ib]->Q) {
      arc->r[ic++] = rb;
      rel_t *rb2 = arb->r[ib++];
      while (rb2 && rb->Q == rb2->Q) {
        combine_partial(comb, rb, rb2, numPrimes, n);
        free_rel(rb2);
        rb2 = (ib < arb->count) ? arb->r[ib++] : NULL;
      }
      rb = rb2;
      continue;
    }
    /* rb is next, and does not combine */
    arc->r[ic++] = rb;
    rb = (ib < arb->count) ? arb->r[ib++] : NULL;
  }
  if (ra) {
    for (i = ia - 1; i < ara->count; ++i)
      arc->r[ic++] = ara->r[i];
  }
  arc->count = ic;
  /* now swap contents of ara and arc, and free the old ara */
  arel_t tmp = *ara;
  *ara = *arc;
  *arc = tmp;
  steal_arel(arc);
  return ara->count;
}

unsigned long read_matrix(
  arel_t *full,
  la_col_t *colarray,
  unsigned long relsFound,
  unsigned long relSought,
  mpz_t n
) {
#ifdef ERROR
  mpz_t test1, test2;
  mpz_init(test1);
  mpz_init(test2);
#endif

  for (unsigned int i = 0; i < full->count && relsFound < relSought; ++i) {
    rel_t *r = full->r[i];
    for (unsigned int j = 0; j < r->f.count; ++j) {
      spp_t *f = &r->f.fact[j];
      if (f->e & 1)
        xorColEntry(colarray, relsFound, (unsigned long)f->p);
    }

#ifdef ERROR
    mpz_set_ui(test1,1);
    for (unsigned long j = 0; j < r->f.count; ++j) {
      mpz_set_ui(test2, r->f.fact[j].p);
      mpz_powm_ui(test2, test2, r->f.fact[j].e, n);
      mpz_mul(test1, test1, test2);
      if ((j % 30) == 0)
        mpz_mod(test1, test1, n);
    }
    mpz_mod(test1, test1, n);
    mpz_mul(test2, r->X, r->X);
    mpz_mod(test2, test2, n);
    if (mpz_cmp(test1, test2) != 0) {
      mpz_add(test1, test1, test2);
      if (mpz_cmp(test1, n) != 0) {
        clearCol(colarray, relsFound);
        mpz_sub(test1, test1, test2);
        gmp_printf("Product mismatch: for %Zu got %Zu with factors",
            test2, test1);
        for (unsigned int j = 0; j < r->f.count; ++j) {
          spp_t *f = &r->f.fact[j];
          printf(" %u", factorBase[f->p]);
          if (f->e > 1)
            printf("^%u", f->e);
        }
        printf("\n");
      } else
        ++relsFound;
    } else
      ++relsFound;
#else
    ++relsFound;
#endif
  }
#ifdef ERROR
  mpz_clear(test1);
  mpz_clear(test2);
#endif
  return relsFound;
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
  61, 67, 71, 73, 79, 83, 89, 97
};
#define NUMMULTS (sizeof(multipliers) / sizeof(unsigned long))

#ifndef M_LN2
# define M_LN2 0.69314718055994530942
#endif

static unsigned long knuthSchroeppel(mpz_t n, unsigned long numPrimes) {
  unsigned int i, j, best_mult, knmod8;
  unsigned int maxprimes = (2 * numPrimes <= 1000) ? 2 * numPrimes : 1000;
  float best_score, contrib;
  float scores[NUMMULTS];
  mpz_t temp;

  mpz_init(temp);

  for (i = 0; i < NUMMULTS; ++i) {
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
    for (i = 1; i < maxprimes; ++i) {
      prime = prime_iterator_next(&iter);
      modp = mpz_mod_ui(temp, n, prime);
      contrib = logf((float)prime) / (float)(prime - 1);

      for (j = 0; j < NUMMULTS; ++j) {
        knmodp = (modp * multipliers[j]) % prime;
        if (knmodp == 0) {
          scores[j] -= contrib;
        } else {
          mpz_set_ui(temp, knmodp);
          if (mpz_kronecker_ui(temp, prime) == 1)
            scores[j] -= 2 * contrib;
        }
      }
    }
    prime_iterator_destroy(&iter);
  }
  mpz_clear(temp);

  best_score = 1000.0;
  best_mult = 1;
  for (i = 0; i < NUMMULTS; ++i) {
    float score = scores[i];
    if (score < best_score) {
      best_score = score;
      best_mult = multipliers[i];
    }
  }
  /* gmp_printf("%Zd mult %u\n", n, best_mult); */
  return best_mult;
}

static void initFactorBase(void) {
  factorBase = 0;
  primeSizes = 0;
}
static void clearFactorBase(void) {
  if (factorBase) {
    Safefree(factorBase);
    factorBase = 0;
  }
  if (primeSizes) {
    Safefree(primeSizes);
    primeSizes = 0;
  }
}

/*========================================================================
   Compute Factor Base:

   Function: Computes B primes p for which n is a square mod p,
   allocates memory and stores them in an array pointed to by factorBase.
   Additionally allocates and computes the primeSizes array.

========================================================================*/
static void computeFactorBase(
  mpz_t n, unsigned long B, unsigned long multiplier
) {
  UV p;
  UV primesinbase = 0;
  PRIME_ITERATOR(iter);

  if (factorBase) {
    Safefree(factorBase);
    factorBase = 0;
  }
  New(0, factorBase, B, unsigned int);

  factorBase[primesinbase++] = multiplier;
  if (multiplier != 2)
    factorBase[primesinbase++] = 2;
  prime_iterator_setprime(&iter, 3);
  for (p = 3; primesinbase < B; p = prime_iterator_next(&iter))
    if (mpz_kronecker_ui(n, p) == 1)
      factorBase[primesinbase++] = p;
  prime_iterator_destroy(&iter);
#ifdef LARGESTP
  gmp_printf("Largest prime less than %lu\n", (unsigned long)p);
#endif

  /* Allocate and compute the number of bits required to store each prime */
  New(0, primeSizes, B, unsigned char);
  for (p = 0; p < B; ++p)
    primeSizes[p] = (unsigned char)floor(
      log(factorBase[p]) / log(2.0) - SIZE_FUDGE + 0.5
    );
}

/*===========================================================================
   Tonelli-Shanks:

   Function: Performs Tonelli-Shanks on n mod every prime in the factor base

===========================================================================*/
static void tonelliShanks(unsigned long numPrimes, mpz_t n, mpz_t *sqrts) {
  unsigned long i;
  mpz_t fbprime, t1, t2, t3, t4;

  mpz_init(fbprime);
  mpz_init(t1); mpz_init(t2); mpz_init(t3); mpz_init(t4);

  mpz_set_ui(sqrts[0], 0);
  for (i = 1; i < numPrimes; ++i) {
    mpz_set_ui(fbprime, factorBase[i]);
    sqrtmodp_t(sqrts[i], n, fbprime, t1, t2, t3, t4);
  }
  mpz_clear(t1); mpz_clear(t2); mpz_clear(t3); mpz_clear(t4);
  mpz_clear(fbprime);
}

/*==========================================================================
   evaluateSieve:

   Function: searches sieve for relations and sticks them into a list.

===========================================================================*/
static void evaluateSieve(
  unsigned long numPrimes,
  unsigned long Mdiv2,
  unsigned long ctimesreps,
  unsigned long M,
  unsigned char *sieve,
  mpz_t A,
  mpz_t B,
  mpz_t C,
  unsigned long *soln1,
  unsigned long *soln2,
  unsigned char *flags,
  la_col_t *colarray,
  unsigned long *aind,
  int min,
  int s,
  arel_t *rels,   /* to add full relations */
  arel_t *lpnew,  /* to add partial relations */
  mpz_t temp,
  mpz_t temp2,
  mpz_t temp3,
  mpz_t res
) {
  long i, j;
  unsigned int k;
  unsigned int exponent;
  unsigned char vv;
  unsigned char extra;
  unsigned int modp;
  unsigned long *sieve2;
  unsigned char bits;
  rel_t *rel = NULL;  /* new relation to add */

  i = 0;
  j = 0;
  sieve2 = (unsigned long *)sieve;
#ifdef POLS
  gmp_printf("%Zdx^2%+Zdx\n%+Zd\n", A, B, C);
#endif

  while ((unsigned long)j < M / sizeof(unsigned long)) {
    do {
      while (!(sieve2[j] & SIEVEMASK))
        ++j;
      i = j * sizeof(unsigned long);
      ++j;
      while ((unsigned long)i < j * sizeof(unsigned long)
        && sieve[i] < threshold
      )
        ++i;
    } while (sieve[i] < threshold);

    if ((unsigned long)i < M) {
      mpz_set_ui(temp, i + ctimesreps);
      mpz_sub_ui(temp, temp, Mdiv2); /* X              */
      mpz_set(temp3, B);             /* B              */
      mpz_addmul(temp3, A, temp);    /* AX + B         */
      mpz_add(temp2, temp3, B);      /* AX + 2B        */
      mpz_mul(temp2, temp2, temp);   /* AX^2 + 2BX     */
      mpz_add(res, temp2, C);        /* AX^2 + 2BX + C */

      bits = mpz_sizeinbase(res, 2) - errorbits;
      extra = 0;

      if (rel == NULL)
        rel = new_rel();
      else
        reset_factors(rel);
      if (factorBase[0] != 1 && mpz_divisible_ui_p(res, factorBase[0])) {
        extra += primeSizes[0];
        if (factorBase[0] == 2) {
          exponent = mpz_scan1(res, 0);
          mpz_tdiv_q_2exp(res, res, exponent);
        } else {
          mpz_set_ui(temp, factorBase[0]);
          exponent = mpz_remove(res, res, temp);
        }
        add_factor(rel, 0, exponent);
      }

      if (mpz_divisible_ui_p(res, factorBase[1])) {
        if (factorBase[1] == 2) {
          exponent = mpz_scan1(res, 0);
          mpz_tdiv_q_2exp(res, res, exponent);
        } else {
          mpz_set_ui(temp, factorBase[1]);
          exponent = mpz_remove(res, res, temp);
        }
        add_factor(rel, 1, exponent);
        /* (hv) simpqs-2.0 adds exponent rather than primeSizes[1]
         * here, no idea why */
        extra += exponent;
      }

      for (k = 2; k < firstprime; ++k) {
        modp = (i + ctimesreps) % factorBase[k];
        if (soln2[k] != (unsigned long)-1) {
          if (modp == soln1[k] || modp == soln2[k]) {
            mpz_set_ui(temp, factorBase[k]);
            exponent = mpz_remove(res, res, temp);
            CHECK_EXPONENT(exponent, k);
            PRINT_FB(exponent, k);
            extra += primeSizes[k];
            add_factor(rel, k, exponent);
          }
        } else {
          mpz_set_ui(temp, factorBase[k]);
          exponent = mpz_remove(res, res, temp);
          if (exponent) {
            PRINT_FB(exponent, k);
            extra += primeSizes[k];
            add_factor(rel, k, exponent);
          }
        }
      }
      sieve[i] += extra;
      if (sieve[i] >= bits) {
        vv = (unsigned char)1 << (i & 7);
        for (k = firstprime; k < secondprime && extra < sieve[i]; ++k) {
          modp = (i + ctimesreps) % factorBase[k];
          if (soln2[k] != (unsigned long)-1) {
            if (modp == soln1[k] || modp == soln2[k]) {
              mpz_set_ui(temp, factorBase[k]);
              exponent = mpz_remove(res, res, temp);
              CHECK_EXPONENT(exponent, k);
              PRINT_FB(exponent, k);
              extra += primeSizes[k];
              add_factor(rel, k, exponent);
            }
          } else {
            mpz_set_ui(temp, factorBase[k]);
            exponent = mpz_remove(res, res, temp);
            if (exponent) {
              PRINT_FB(exponent, k);
              extra += primeSizes[k];
              add_factor(rel, k, exponent);
            }
          }
        }

        for (k = secondprime; k < numPrimes && extra < sieve[i]; ++k) {
          if (flags[k] & vv) {
            modp = (i + ctimesreps) % factorBase[k];
            if (modp == soln1[k] || modp == soln2[k]) {
              mpz_set_ui(temp, factorBase[k]);
              exponent = mpz_remove(res, res, temp);
              CHECK_EXPONENT(exponent, k);
              PRINT_FB(exponent, k);
              extra += primeSizes[k];
              add_factor(rel, k, exponent);
            }
          }
        }

        if (mpz_cmp_ui(res, 1000) > 0) {
          if (mpz_cmp_ui(res, largeprime) < 0) {
#ifdef RELPRINT
            gmp_printf(" %Zd\n", res);
#endif
            for (int i = 0; i < s; ++i)
              add_factor(rel, aind[i] + min, 1);
            sort_fact(&rel->f);
            mpz_set(rel->X, temp3);
            rel->Q = mpz_get_ui(res);
            save_rel(lpnew, rel);
            rel = NULL;
          }
        } else {
          mpz_neg(res, res);
          if (mpz_cmp_ui(res, 1000) > 0) {
            if (mpz_cmp_ui(res, largeprime) < 0) {
#ifdef RELPRINT
              gmp_printf(" %Zd\n", res);
#endif
              for (int i = 0; i < s; ++i)
                add_factor(rel, aind[i] + min, 1);
              sort_fact(&rel->f);
              mpz_set(rel->X, temp3);
              rel->Q = mpz_get_ui(res);
              save_rel(lpnew, rel);
              rel = NULL;
            }
          } else {
#ifdef RELPRINT
            printf("....R\n");
#endif
            for (int i = 0; i < s; ++i)
              add_factor(rel, aind[i] + min, 1);
            sort_fact(&rel->f);
            mpz_set(rel->X, temp3);
            save_rel(rels, rel);
            rel = NULL;
          }
        }
      } else {
#ifdef RELPRINT
        printf("\n");
#endif
      }
      ++i;
    }
  }
  if (rel != NULL)
    free_rel(rel);
}

static void update_solns(
  unsigned long first, unsigned long limit,
  unsigned long *soln1, unsigned long *soln2,
  int polyadd, const unsigned long *polycorr
) {
  unsigned int prime;
  unsigned long p, correction;

  for (prime = first; prime < limit; ++prime) {
    if (soln2[prime] == (unsigned long)-1)
      continue;
    p = factorBase[prime];
    correction = (polyadd) ? p - polycorr[prime] : polycorr[prime];
    soln1[prime] += correction;
    while (soln1[prime] >= p)
      soln1[prime] -= p;
    soln2[prime] += correction;
    while (soln2[prime] >= p)
      soln2[prime] -= p;
  }
}

static void set_offsets(
  unsigned char *const sieve,
  const unsigned long *const soln1, const unsigned long *const soln2,
  unsigned char **offsets1, unsigned char **offsets2
) {
  unsigned int prime;
  for (prime = firstprime; prime < midprime; ++prime)
    if (soln2[prime] == (unsigned long)-1) {
      offsets1[prime] = 0;
      offsets2[prime] = 0;
    } else {
      offsets1[prime] = sieve + soln1[prime];
      offsets2[prime] = sieve + soln2[prime];
    }
}

/*=============================================================================
   Sieve:

   Function: Allocates space for a sieve of M integers and sieves the interval
             starting at start

=============================================================================*/
static void sieveInterval(
  unsigned long M, unsigned char *sieve, int more,
  unsigned char **offsets1, unsigned char **offsets2
) {
  unsigned int prime, p;
  unsigned char size;
  unsigned char *pos1;
  unsigned char *pos2;
  unsigned char *end = sieve + M;
  unsigned char *bound;
  ptrdiff_t diff;

  for (prime = firstprime; prime < midprime; ++prime) {
    if (offsets1[prime] == 0)
      continue;
    p = factorBase[prime];
    size = primeSizes[prime];
    pos1 = offsets1[prime];
    pos2 = offsets2[prime];
    diff = pos2 - pos1;
    /* if pos1 < bound, then both *pos1 and *pos2 can be written to. */
    bound = end - 4 * p;

    /* Write both values, unrolled 4 times. */
    while (pos1 < bound) {
      pos1[0    ] += size;  pos1[        diff] += size;
      pos1[1 * p] += size;  pos1[1 * p + diff] += size;
      pos1[2 * p] += size;  pos1[2 * p + diff] += size;
      pos1[3 * p] += size;  pos1[3 * p + diff] += size;
      pos1 += 4 * p;
    }

    /* Write both values */
    while (pos1 < end && pos1 + diff < end) {
      pos1[0] += size;
      pos1[diff] += size;
      pos1 += p;
    }
    pos2 = pos1 + diff;    /* Restore pos2 */

    /* Finish writing to pos1 and pos2 */
    if (pos1 < end) {
      *pos1 += size;
      pos1 += p;
    }
    if (pos2 < end) {
      *pos2 += size;
      pos2 += p;
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
static void sieve2(
  unsigned long M, unsigned long numPrimes, unsigned char *sieve,
  const unsigned long *soln1, const unsigned long *soln2,
  unsigned char *flags
) {
  unsigned int prime, p;
  unsigned char size;
  unsigned char *pos1;
  unsigned char *pos2;
  unsigned char *end = sieve + M;

  for (prime = midprime; prime < secondprime; ++prime) {
    /* (hv) is this possible? needs aind[0..s-1] + min >= midprime */
    if (soln2[prime] == (unsigned long)-1)
      continue;

    p = factorBase[prime];
    size = primeSizes[prime];
    pos1 = sieve + soln1[prime];
    pos2 = sieve + soln2[prime];
    while (end - pos1 > 0 && end - pos2 > 0) {
      *pos1 += size;  pos1 += p;
      *pos2 += size;  pos2 += p;
    }
    if (end - pos2 > 0)
      *pos2 += size;
    if (end - pos1 > 0)
      *pos1 += size;
  }

  for (prime = secondprime; prime < numPrimes; ++prime) {
    /* (hv) is this possible? needs aind[0..s-1] + min >= secondprime */
    if (soln2[prime] == (unsigned long)-1)
      continue;

    p = factorBase[prime];
    size = primeSizes[prime];
    pos1 = sieve + soln1[prime];
    pos2 = sieve + soln2[prime];

    while (end - pos1 > 0) {
      flags[prime] |= (unsigned char)1 << ((pos1 - sieve) & 7);
      *pos1 += size;
      pos1 += p;
    }

    while (end - pos2 > 0) {
      flags[prime] |= (unsigned char)1 << ((pos2 - sieve) & 7);
      *pos2 += size;
      pos2 += p;
    }
  }
}

/*============================================================================

   random:

   Function: Generates a pseudo-random integer between 0 and n-1 inclusive

============================================================================*/
static unsigned long randval = 2994439072U;
static unsigned long silly_random(unsigned long upto) {
  randval = ((unsigned long)randval * 1025416097U + 286824428U)
          % (unsigned long)4294967291U;
  return randval % upto;
}

/*============================================================================

  danaj: added these routines to reduce the set of factors to co-primes.
  It's not the most efficient solution, but it's trivial in time compared
  to the loop it's in, much less the rest of the QS.  It gives us a nice
  set of factors back, which is much more useful than the essentially
  random combinations we discover.

============================================================================*/

/* Verify that the factor reduction hasn't broken anything */
static void verify_factor_array(mpz_t n, mpz_t *farray, int nfacs) {
  int i, j;
  mpz_t t;
  mpz_init_set_ui(t, 1);
  /* Assert we don't have duplicates */
  for (i = 0; i < nfacs; ++i)
    for (j = i + 1; j < nfacs; ++j)
      if (mpz_cmp(farray[i], farray[j]) == 0) {
        gmp_printf("duplicate: F[%d] = F[%d] = %Zd\n", i, j, farray[i]);
        croak("assert");
      }
  /* Assert that all factors multiply to n */
  for (i = 0; i < nfacs; ++i)
    mpz_mul(t, t, farray[i]);
  if (mpz_cmp(t, n) != 0) {
    gmp_printf("farray doesn't multiply: n=%Zd t=%Zd\n", n, t);
    croak("assert");
  }
  /* Assert that gcd of each non-identical factor is 1 */
  for (i = 0; i < nfacs; ++i)
    for (j = i + 1; j < nfacs; ++j)
      if (mpz_cmp(farray[i], farray[j]) != 0) {
        mpz_gcd(t, farray[i], farray[j]);
        if (mpz_cmp_ui(t, 1) != 0) {
          gmp_printf("gcd: farray[%d] = %Zd  farray[%d] = %Zd\n",
              i, farray[i], j, farray[j]);
          croak("assert");
        }
      }
  mpz_clear(t);
}

static int allprime_factor_array(mpz_t *farray, int nfacs) {
  int i;
  for (i = 0; i < nfacs; ++i)
    if (!mpz_probab_prime_p(farray[i], 5))   /* Be lazy */
      return 0;
  return 1;
}

static int insert_factor(mpz_t n, mpz_t *farray, int nfacs, mpz_t f) {
  int i, j;
  mpz_t t, t2;

  if (mpz_cmp_ui(f, 1) <= 0)
    return nfacs;

  /* skip duplicates */
  /* (hv) hmm, what if p^2 divides n? */
  for (i = 0; i < nfacs; ++i)
    if (mpz_cmp(farray[i], f) == 0)
      break;
  if (i != nfacs)
    return nfacs;

  /* Look for common factors in all the existing set */
  /* for (i = 0; i < nfacs; i++) gmp_printf("  F[%d] = %Zd\n", i, farray[i]); */
  mpz_init(t);
  mpz_init(t2);
  for (i = 0; i < nfacs; ++i) {
    mpz_gcd(t, farray[i], f);
    if (mpz_cmp_ui(t, 1) == 0) /* t=1:   F and f unchanged */
      continue;
    mpz_divexact(t2, farray[i], t);    /* t2 = F/t */
    mpz_divexact(f, f, t);             /* f  = f/t */
    /* Remove the old farray[i] */
    for (j = i + 1; j < nfacs; ++j)
      mpz_set(farray[j - 1], farray[j]);
    mpz_set_ui(farray[nfacs--], 0);
    /* Insert F/t, t, f/t */
    nfacs = insert_factor(n, farray, nfacs, t2);
    nfacs = insert_factor(n, farray, nfacs, t);
    nfacs = insert_factor(n, farray, nfacs, f);
    i = 0;
    break;
  }
  /* If nothing common, insert it. */
  if (i == nfacs)
    mpz_set(farray[nfacs++], f);
  mpz_clear(t);
  mpz_clear(t2);
  return nfacs;
}

/*============================================================================
   mainRoutine:

   Function: Generates the polynomials, initialises and calls the sieve,
             implementing cache blocking (breaking the sieve interval into
             small blocks for the small primes).

============================================================================*/
static int mainRoutine(
  unsigned long numPrimes,
  unsigned long Mdiv2,
  unsigned long relSought,
  mpz_t n,
  mpz_t *farray,
  unsigned long multiplier
) {
  mpz_t A, B, C, D, Bdivp2, nsqrtdiv, temp, temp2, temp3, temp4;
  int i, j, s, fact, span, min, nfactors, verbose;
  unsigned long u1, p, reps, M, Mq, Mr;
  unsigned long curves = 0;
  unsigned int   *primecount;
  unsigned char  *sieve;
  unsigned long  *aind;
  unsigned long  *amodp;
  unsigned long  *Ainv;
  unsigned long  *soln1;
  unsigned long  *soln2;
  unsigned char  *flags;
  unsigned long **Ainv2B;
  unsigned char **offsets;
  unsigned char **offsets2;
  la_col_t       *colarray;
  mpz_t          *Bterms;
  mpz_t          *sqrts;

  unsigned long next_cutoff = (relSought - 1) / 40 + 1;
  unsigned long next_inc = next_cutoff;

  arel_t *frels = new_arel();   /* all natural full relations */
  arel_t *rels = new_arel();    /* new natural full relations, to be merged */
  arel_t *lprels = new_arel();  /* all partial relations */
  arel_t *lpnew = new_arel();   /* new partial relations, to be merged */
  arel_t *comb = new_arel();    /* new full relations found from partials */
  arel_t *flprels = new_arel(); /* all combined full relations */

  verbose = get_verbose_level();
  s = mpz_sizeinbase(n, 2) / 28 + 1;

  Newz(0, aind,          s, unsigned long);
  Newz(0, amodp,         s, unsigned long);
  Newz(0, Ainv,  numPrimes, unsigned long);
  Newz(0, soln1, numPrimes, unsigned long);
  Newz(0, soln2, numPrimes, unsigned long);
  Newz(0, Ainv2B,        s, unsigned long *);
  New( 0, Bterms,        s, mpz_t);
  Newz(0, colarray, relSought, la_col_t);
  if (aind == 0 || amodp == 0 || Ainv == 0
    || soln1 == 0 || soln2 == 0 || Ainv2B == 0 || Bterms == 0
    || colarray == 0
  )
    croak("SIMPQS: Unable to allocate memory!\n");

  flags = 0;
  if (midprime < numPrimes) {
    Newz(0, flags, numPrimes, unsigned char);
    if (flags == 0)
      croak("SIMPQS: Unable to allocate memory!\n");
  }

  for (i = 0; i < s; ++i) {
    New(0, Ainv2B[i], numPrimes, unsigned long);
    if (Ainv2B[i] == NULL)
      croak("SIMPQS: Unable to allocate memory!\n");
    mpz_init(Bterms[i]);
  }

  /* one extra word for sentinel */
  Newz(0, sieve,     Mdiv2 * 2 + sizeof(unsigned long), unsigned char);
  New( 0, offsets,   midprime, unsigned char *);
  New( 0, offsets2,  midprime, unsigned char *);

  if (sieve == 0 || offsets == 0 || offsets2 == 0)
    croak("SIMPQS: Unable to allocate memory!\n");

  mpz_init(A); mpz_init(B); mpz_init(C); mpz_init(D);
  mpz_init(Bdivp2); mpz_init(nsqrtdiv);
  mpz_init(temp); mpz_init(temp2); mpz_init(temp3); mpz_init(temp4);

  /* Compute sqrt(n) mod factorbase[i] */
  New(0, sqrts, numPrimes, mpz_t);
  if (sqrts == 0)
    croak("SIMPQS: Unable to allocate memory!\n");
  for (i = 0; i < numPrimes; ++i)
    mpz_init(sqrts[i]);
  tonelliShanks(numPrimes, n, sqrts);

  /* Compute min A_prime and A_span */
  mpz_mul_ui(temp, n, 2);
  mpz_sqrt(temp, temp);
  mpz_tdiv_q_ui(nsqrtdiv, temp, Mdiv2);
  mpz_root(temp, nsqrtdiv, s);
  for (fact = 0; mpz_cmp_ui(temp, factorBase[fact]) >= 0; ++fact)
    ;
  span = numPrimes / s / s / 2;
  min = fact - span / 2;
  while (min > 0 && (fact * fact) / min - min < span)
    --min;

#ifdef ADETAILS
  printf("s = %d, fact = %d, min = %d, span = %d\n", s, fact, min, span);
#endif

  /* Compute first polynomial and adjustments */

  while (frels->count + flprels->count < relSought) {
    mpz_set_ui(A, 1);
    for (i = 0; i < s - 1; ) {
      unsigned long ran = span / 2 + silly_random(span / 2);
      j = -1;
      while (j != i) {
        ++ran;
        for (j = 0; j < i && aind[j] != ran; ++j)
          ;
      }
      aind[i] = ran;
      mpz_mul_ui(A, A, factorBase[ran + min]);
      ++i;
      if (i < s - 1) {
        j = -1;
        ran = ((min + span / 2) * (min + span / 2)) / (ran + min)
            - silly_random(10) - min;
        while (j != i) {
          ++ran;
          for (j = 0; j < i && aind[j] != ran; ++j)
            ;
        }
        aind[i] = ran;
        mpz_mul_ui(A, A, factorBase[ran + min]);
        ++i;
      }
    }

    mpz_fdiv_q(temp, nsqrtdiv, A);
    for (fact = 1; mpz_cmp_ui(temp, factorBase[fact]) >= 0; ++fact)
      ;
    fact -= min;
    do {
      for (j = 0; j < i && aind[j] != fact; ++j)
        ;
      ++fact;
    } while (j != i);
    --fact;
    aind[i] = fact;
    mpz_mul_ui(A, A, factorBase[fact + min]);

    for (i = 0; i < s; ++i) {
      p = factorBase[aind[i] + min];
      mpz_tdiv_q_ui(temp, A, p);
      amodp[i] = mpz_fdiv_r_ui(temp, temp, p);
      /* (hv) how do we know A/p^2 fits ui? */
      mpz_set_ui(temp, modinverse(mpz_get_ui(temp), p));
      mpz_mul(temp, temp, sqrts[aind[i] + min]);
      mpz_fdiv_r_ui(temp, temp, p);
      if (mpz_cmp_ui(temp, p / 2) > 0) {
        mpz_sub_ui(temp, temp, p);
        mpz_neg(temp, temp);
      }
      mpz_mul(temp, temp, A);
      mpz_tdiv_q_ui(Bterms[i], temp, p);
    }

    mpz_set(B, Bterms[0]);
    for (i = 1; i < s; ++i)
      mpz_add(B, B, Bterms[i]);

    for (i = 0; i < numPrimes; ++i) {
      p = factorBase[i];
      Ainv[i] = modinverse(mpz_fdiv_r_ui(temp, A, p), p);

      for (j = 0; j < s; ++j) {
        mpz_fdiv_r_ui(temp, Bterms[j], p);
        mpz_mul_ui(temp, temp, 2 * Ainv[i]);
        Ainv2B[j][i] = mpz_fdiv_r_ui(temp, temp, p);
      }

      mpz_fdiv_r_ui(temp, B, p);
      mpz_sub(temp, sqrts[i], temp);
      mpz_mul_ui(temp, temp, Ainv[i]);
      mpz_add_ui(temp, temp, Mdiv2);
      soln1[i] = mpz_fdiv_r_ui(temp, temp, p);
      mpz_sub_ui(temp, sqrts[i], p);
      mpz_neg(temp, temp);
      mpz_mul_ui(temp, temp, 2 * Ainv[i]);
      soln2[i] = mpz_fdiv_ui(temp, p) + soln1[i];
    }

    for (int polyindex = 1; polyindex < (1 << (s - 1)) - 1; ++polyindex) {
      int polyadd;
      unsigned long *polycorr;
      for (j = 0; j < s; ++j)
        if (((polyindex >> j) & 1) != 0)
          break;
      polyadd = (((polyindex >> j) & 2) != 0);
      if (polyadd) {
        mpz_add(B, B, Bterms[j]);
        mpz_add(B, B, Bterms[j]);
      } else {
        mpz_sub(B, B, Bterms[j]);
        mpz_sub(B, B, Bterms[j]);
      }
      polycorr = Ainv2B[j];

      for (j = 0; j < s; ++j) {
        int findex = aind[j] + min;
        p = factorBase[findex];
        mpz_fdiv_r_ui(D, n, p * p);
        mpz_fdiv_r_ui(Bdivp2, B, p * p);
        mpz_mul_ui(temp, Bdivp2, amodp[j]);
        u1 = modinverse(mpz_fdiv_ui(temp, p), p);
        mpz_mul(temp, Bdivp2, Bdivp2);
        mpz_sub(temp, temp, D);
        mpz_neg(temp, temp);
        mpz_tdiv_q_ui(temp, temp, p);
        mpz_mul_ui(temp, temp, u1);
        mpz_add_ui(temp, temp, Mdiv2);
        soln1[findex] = mpz_fdiv_ui(temp, p);
        soln2[findex] = (unsigned long)-1;
      }

      /* Count the number of polynomial curves used so far and compute
       * the C coefficient of our polynomial */
      ++curves;
      mpz_mul(C, B, B);
      mpz_sub(C, C, n);
      mpz_divexact(C, C, A);

      /* Do the sieving and relation collection */
      M = Mdiv2 * 2;
      Mq = M / CACHEBLOCKSIZE;
      Mr = M % CACHEBLOCKSIZE;

      /* set the solns1 and solns2 arrays */
      update_solns(1, numPrimes, soln1, soln2, polyadd, polycorr);
      /* Clear sieve and insert sentinel at end (used in evaluateSieve) */
      memset(sieve, 0, M * sizeof(unsigned char));
      sieve[M] = 255;
      /* Sieve [midprime, numPrimes) */
      if (midprime < numPrimes)
        sieve2(M, numPrimes, sieve, soln1, soln2, flags);
      /* Set the offsets and offsets2 arrays used for small sieve */
      set_offsets(sieve, soln1, soln2, offsets, offsets2);
      /* Sieve [firstprime, midprime) */
      sieveInterval(CACHEBLOCKSIZE, sieve, 1, offsets, offsets2);
      if (Mq > 0) {
        unsigned long maxreps = Mq - 1;
        for (reps = 1; reps < maxreps; ++reps)
          sieveInterval(
            CACHEBLOCKSIZE, sieve + CACHEBLOCKSIZE * reps, 1, offsets, offsets2
          );
        if (Mr == 0)
          sieveInterval(
            CACHEBLOCKSIZE, sieve + CACHEBLOCKSIZE * reps, 0, offsets, offsets2
          );
        else {
          sieveInterval(
            CACHEBLOCKSIZE, sieve + CACHEBLOCKSIZE * reps, 1, offsets, offsets2
          );
          ++reps;
          sieveInterval(
            Mr            , sieve + CACHEBLOCKSIZE * reps, 0, offsets, offsets2
          );
        }
      }

      evaluateSieve(
        numPrimes, Mdiv2,
        0, M, sieve, A, B, C,
        soln1, soln2, flags, colarray, aind,
        min, s,
        rels, lpnew, temp, temp2, temp3, temp4
      );

      if (2 * (rels->count + frels->count) >= next_cutoff) {
        sort_rel(lpnew);
        /* full relations found while merging are extracted into comb */
        merge_partial(lprels, lpnew, comb, numPrimes, n);
        reset_arel(lpnew);

        sort_rel(rels);
        merge_full(frels, rels);
        reset_arel(rels);

        /* (hv) why keep combined full relations separate from naturals? */
        sort_rel(comb);
        merge_full(flprels, comb);
        reset_arel(comb);
#ifdef COUNT
        printf("%u full, %u combined, %u partial\n",
            frels->count, flprels->count, lprels->count);
#endif
        if (
          next_cutoff < relSought
          && next_cutoff + next_inc / 2 >= relSought
        )
          next_inc = next_inc / 2;
        next_cutoff += next_inc;
      }
    }

#ifdef COUNT
    if ((curves % 20) == 0)
      printf("%lu curves.\n", curves);
#endif
  }

#ifdef CURPARTS
  printf("%lu curves, %u partials.\n", curves, lprels->count);
#endif

  if (verbose > 4)
    printf("# qs done sieving\n");

  /* Free everything we don't need for the linear algebra */

  for (i = 0; i < numPrimes; ++i)
    mpz_clear(sqrts[i]);
  Safefree(sqrts);
  for (i = 0; i < s; ++i) {
    Safefree(Ainv2B[i]);
    mpz_clear(Bterms[i]);
  }
  Safefree(aind);
  Safefree(amodp);
  Safefree(Ainv);
  Safefree(soln1);
  Safefree(soln2);
  Safefree(Ainv2B);
  Safefree(Bterms);
  if (flags)
    Safefree(flags);

  Safefree(sieve);    sieve = 0;
  Safefree(offsets);  offsets = 0;
  Safefree(offsets2); offsets2 = 0;

  mpz_clear(A);  mpz_clear(B);  mpz_clear(C);  mpz_clear(D);
  mpz_clear(Bdivp2); mpz_clear(nsqrtdiv);

  /* Do the matrix algebra step */

  unsigned long ncols = relSought;
  unsigned long nrows = numPrimes;

#ifdef ERRORS
  for (j = frels->count; j < relSought; ++j)
    if (colarray[j].weight != 0)
      printf("Dirty at col %d\n", j);
#endif

#ifdef COUNT
  printf("%u relations found in total!\n", frels->count + flprels->count);
#endif

  unsigned long relsFound = 0;
  relsFound = read_matrix(frels, colarray, 0, relSought, n);
  relsFound = read_matrix(flprels, colarray, relsFound, relSought, n);

/* temporary, while we store frels and flprels separately: return a (rel_t *)
 * pointer to the relation with the given index */
#define RELINDEX(ri) ((ri >= frels->count) ? flprels->r[ri - frels->count] : frels->r[ri])

#ifdef ERRORS
  for (j = 0; j < relSought; ++j)
    if (colarray[j].orig != j) {
      printf("Column numbering error, %d\n", j);
      colarray[j].orig = j;
    }

  for (j = 0; j < relSought; ++j)
    for (i = 0; i < colarray[j].weight; ++i)
      if (colarray[j].data[i] > numPrimes)
        printf("Error prime too large: %lu\n", colarray[j].data[i]);

  mpz_t test1, test2, test3;
  mpz_init(test1);
  mpz_init(test2);
  mpz_init(test3);
  unsigned int *exps = malloc(numPrimes * sizeof(unsigned int));
  for (j = 0; j < relSought; ++j) {
    for (i = 0; i < numPrimes; ++i)
      exps[i] = 0;
    rel_t *r = RELINDEX(j);
    mpz_set_ui(test1, 1);
    for (i = 0; i < r->f.count; ++i) {
      mpz_ui_pow_ui(test2, factorBase[ r->f.fact[i].p ], r->f.fact[i].e);
      mpz_mul(test1, test1, test2);
      exps[ r->f.fact[i].p ] += r->f.fact[i].e;
    }
    mpz_mod(test1, test1, n);
    mpz_mul(test2, r->X, r->X);
    mpz_mod(test2, test2, n);
    if (mpz_cmp(test1, test2) != 0) {
      mpz_add(test3, test1, test2);
      if (mpz_cmp(test3, n) != 0) {
        gmp_printf("%Zd !=\n%Zd\nin column %d and\n", test3, n, j);
        gmp_printf("%Zd !=\n%Zd\n\n", test1, test2);
      }
    }
    for (i = 0; i < colarray[j].weight; ++i) {
      if (exps[colarray[j].data[i]] % 2 != 1)
        printf("Col %d, row %d incorrect\n", j, i);
      exps[colarray[j].data[i]] = 0;
    }
    for (i = 0; i < numPrimes; ++i)
      if ((exps[i] % 2) == 1)
        printf("exps[%d] is not even in row %d\n", i, j);
  }
  free(exps);
  mpz_clear(test1);
  mpz_clear(test2);
  mpz_clear(test3);
#endif

  reduce_matrix(&nrows, &ncols, colarray);

#ifdef ERRORS
  exps = (unsigned int *)malloc(numPrimes * sizeof(unsigned int));
  for (j = 0; j < ncols; ++j) {
    for (i = 0; i < numPrimes; ++i)
      exps[i] = 0;
    unsigned int index = colarray[j].orig;
    rel_t *r = RELINDEX(index);
    for (i = 0; i < r->f.count; ++i)
      exps[r->f.fact[i].p] += r->f.fact[i].e;
    for (i = 0; i < colarray[j].weight; ++i) {
      for (int k = 0; k < i; ++k)
        if (colarray[j].data[i] == colarray[j].data[k])
          printf("Duplicate in column %d: %d, %d\n", j, i, k);
      if ((exps[colarray[j].data[i]] % 2) != 1)
        printf("Col %d, row %d incorrect\n", j, i);
      exps[colarray[j].data[i]] = 0;
    }
    for (i = 0; i < numPrimes; ++i)
      if ((exps[i] % 2) == 1)
        printf("exps[%d] is not even in row %d\n", i, j);
  }
  free(exps);
#endif

  u_int64_t *nullrows;
  long mask = 0;
  int fail_count = 0;
  while (mask == 0 && fail_count < BL_MAX_FAIL) {
    nullrows = block_lanczos(nrows, 0, ncols, colarray);
    if (nullrows == NULL) {
      ++fail_count;
      continue;
    }
    for (i = 0; i < ncols; ++i)
      mask |= nullrows[i];
    if (mask == 0)
      ++fail_count;
  }
  if (fail_count >= BL_MAX_FAIL) {
	gmp_printf(
	  "block_lanczos failed %d times on target %Zd (multiplier %d), giving up",
	  fail_count, n, multiplier
    );
	croak("assert");
  }

  if (verbose > 3) {
    for (i = j = 0; i < 64; ++i)
      if (mask & ((u_int64_t)1 << i))
        ++j;
    printf("%d nullspace vectors found.\n", j);
  }

#ifdef ERRORS
  exps = (unsigned int *)malloc(numPrimes * sizeof(unsigned int));
  for (j = 0; j < ncols; ++j) {
    for (i = 0; i < numPrimes; ++i)
      exps[i] = 0;
    unsigned int index = colarray[j].orig;
    rel_t *r = RELINDEX(index);
    for (i = 0; i < r->f.count; ++i)
      exps[ r->f.fact[i].p ] += r->f.fact[i].e;
    for (i = 0; i < colarray[j].weight; ++i) {
      if ((exps[ colarray[j].data[i] ] % 2) != 1)
        printf("Col %d, row %d incorrect\n", j, i);
      exps[colarray[j].data[i]] = 0;
    }
    for (i = 0; i < numPrimes; ++i)
      if ((exps[i] % 2) == 1)
        printf("exps[%d] is not even in row %d\n", i, j);
  }
  free(exps);
#endif

  /* We want factors of n, not kn, so divide out by the multiplier */
  mpz_divexact_ui(n, n, multiplier);

  /* Now find the factors via square root and gcd */
  mpz_set(farray[0], n);
  nfactors = 1;  /* We have one result -- n */
  New(0, primecount, numPrimes, unsigned int);
  if (primecount == 0)
    croak("SIMPQS: Unable to allocate memory!\n");
  for (int l = 0; l < 64; ++l) {
    while (!(mask & ((u_int64_t)1 << l)))
      ++l;
    mpz_set_ui(temp, 1);
    mpz_set_ui(temp2, 1);
    memset(primecount, 0, numPrimes * sizeof(unsigned int));
    for (i = 0; i < ncols; ++i) {
      if (getNullEntry(nullrows, i, l)) {
        unsigned int index = colarray[i].orig;
        rel_t *r = RELINDEX(index);
        mpz_mul(temp2, temp2, r->X);
        for (j = 0; j < r->f.count; ++j)
          primecount[r->f.fact[j].p] += r->f.fact[j].e;
      }
      if (((i + 1) % 16) == 0)
        mpz_mod(temp2, temp2, n);
    }
    for (j = 0; j < numPrimes; ++j) {
      mpz_set_ui(temp3, factorBase[j]);
      mpz_pow_ui(temp3, temp3, primecount[j] / 2);
      mpz_mul(temp, temp, temp3);
      if (((j + 1) % 16) == 0)
        mpz_mod(temp, temp, n);
    }
    mpz_sub(temp, temp2, temp);
    mpz_gcd(temp, temp, n);
#if 0
    /* (hv) shouldn't every hit give gcd > 1?
     * If so, failures represent a bug somewhere to investigate. */
    if (mpz_cmp_ui(temp, 1) == 0)
      gmp_printf("failed result for l = %d\n", l);
#endif
    /* only non-trivial factors */
    if (mpz_cmp_ui(temp, 1) && mpz_cmp(temp, n)) {
      if (verbose > 4)
        gmp_printf("# qs factor %Zd\n", temp);
      nfactors = insert_factor(n, farray, nfactors, temp);
      verify_factor_array(n, farray, nfactors);
      if (allprime_factor_array(farray, nfactors))
        break;
    }
  }

  /* Free everything remaining */
  free(nullrows);
  free_arel(frels);
  free_arel(rels);
  free_arel(lprels);
  free_arel(lpnew);
  free_arel(comb);
  free_arel(flprels);
  for (i = 0; i < relSought; ++i)
    free(colarray[i].data);
  Safefree(colarray);
  Safefree(primecount);
  mpz_clear(temp);  mpz_clear(temp2);  mpz_clear(temp3);  mpz_clear(temp4);

  return nfactors;
}

int _GMP_simpqs(mpz_t n, mpz_t *farray) {
  unsigned long numPrimes, Mdiv2, multiplier, decdigits, relSought;
  int result = 0;
  int verbose = get_verbose_level();

  mpz_set(farray[0], n);
  decdigits = mpz_sizeinbase(n, 10); /* often 1 too big */
  if (decdigits < MINDIG)
    return 0;

  if (verbose > 2)
    gmp_printf("# qs trying %Zd (%lu digits)\n", n, decdigits);

  /* It's important to remove small factors. */
  /* (hv) this isn't needed for the entry-point used by factor(), which
   * will already have done trial division - can it move to STANDALONE? */
  {
    UV p;
    PRIME_ITERATOR(iter);
    for (p = 2; p < 1000; p = prime_iterator_next(&iter)) {
      if (mpz_cmp_ui(n, p * p) < 0) break;
      while (mpz_divisible_ui_p(n, p)) {
        mpz_set_ui(farray[result++], p);
        mpz_divexact_ui(n, n, p);
      }
    }
    /* (hv) this duplicates work done above when no new factors found */
    decdigits = mpz_sizeinbase(n, 10);
    /* (hv) return here misses mpz_set(farray[result], n) for last factor */
    if (decdigits < MINDIG)
      return result;
    mpz_set(farray[result], n);
  }

  /* Get a preliminary number of primes, pick a multiplier, apply it */
  numPrimes = (decdigits <= 91) ? primesNo[decdigits - MINDIG] : 64000;
  multiplier = knuthSchroeppel(n, numPrimes);
  mpz_mul_ui(n, n, multiplier);
  decdigits = mpz_sizeinbase(n, 10);

  if (decdigits <= 91) {
    numPrimes = primesNo[decdigits - MINDIG];
    Mdiv2 = sieveSize[decdigits - MINDIG] / SIEVEDIV;
    if (Mdiv2 * 2 < CACHEBLOCKSIZE)
      Mdiv2 = CACHEBLOCKSIZE / 2;
    largeprime = 1000 * largeprimes[decdigits - MINDIG];
    secondprime = (numPrimes < SECONDPRIME) ? numPrimes : SECONDPRIME;
    midprime = (numPrimes < MIDPRIME) ? numPrimes : MIDPRIME;
    firstprime = firstPrimes[decdigits - MINDIG];
    errorbits = errorAmounts[decdigits - MINDIG];
    threshold = thresholds[decdigits - MINDIG];
  } else {
    /* (hv) this config makes no sense: it should surely be based on
     * that for 91 digits, which would be:
     *   numPrimes >= 80000
     *   largeprime >= 1000 * 500000 (not 1/10 of that)
     *   errorbits >= 33
     */
    numPrimes = 64000;
    Mdiv2 = 192000 / SIEVEDIV;
    largeprime = numPrimes * 10 * decdigits;
    secondprime = SECONDPRIME;
    midprime = MIDPRIME;
    firstprime = 30;
    errorbits = decdigits / 4 + 2;
    threshold = 43 + (7 * decdigits) / 10;
  }

  if (verbose > 2)
    gmp_printf("# qs    mult %lu, digits %lu, sieving %lu, primes %lu\n",
        multiplier, decdigits, Mdiv2 * 2, numPrimes);

  /* as numPrimes+64, was commented: we probably need fewer than this */
  relSought = numPrimes;
  initFactorBase();
  computeFactorBase(n, numPrimes, multiplier);

  result += mainRoutine(
    numPrimes, Mdiv2, relSought, n, farray + result, multiplier
  );

  clearFactorBase();
  if (verbose > 2) {
    int i;
    gmp_printf("# qs:");
    for (i = 0; i < result; ++i)
      gmp_printf(" %Zd", farray[i]);
    gmp_printf("%s\n", (result) ? "" : " no factors");
  }
  if (verbose > 2 && !result)
    gmp_printf("QS Fail: %Zd (%lu digits)\n", n, decdigits);
  return result;
}

#ifdef STANDALONE_SIMPQS
/* Main Program: factors a user specified number using a quadratic sieve */
int main(int argc, char **argv) {
  int i = 1, nfactors;
  mpz_t n;
  mpz_t *farray;

  while (i < argc && argv[i][0] == '-') {
    char *arg = argv[i++];
    if (arg[1] == '-')
      break;
    if (arg[1] == 'v') {
      set_verbose_level(atoi(&arg[2]));
      continue;
    }
    croak("unknown option '%s'\n", arg);
  }
  if (i + 1 == argc)
    mpz_init_set_str(n, argv[i], 10);
  else
    croak("usage: %s [options] n\n", argv[0]);

  if (mpz_sizeinbase(n, 10) < MINDIG)
    croak("Error in input or number has too few digits.\n");

  _GMP_init();
  New(0, farray, 64, mpz_t);
  for (i = 0; i < 64; ++i)
    mpz_init_set_ui(farray[i], 0);

  nfactors = _GMP_simpqs(n, farray);

  for (i = 0; i < nfactors; ++i)
    gmp_printf("  %Zd\n", farray[i]);
  for (i = 0; i < 64; ++i)
    mpz_clear(farray[i]);
  Safefree(farray);
  mpz_clear(n);
  _GMP_destroy();
}
#endif
