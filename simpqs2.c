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
   gcc -O2 -DSTANDALONE_SIMPQS -DSTANDALONE simpqs2.c lanczos.c utility.c \
     bls75.c ecm.c ecpp.c factor.c gmp_main.c isaac.c lucas_seq.c pbrent63.c \
     primality.c prime_iterator.c random_prime.c real.c rootmod.c squfof126.c \
     tinyqs.c -lgmp -lm

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
#include <limits.h>
#include <math.h>
#include <gmp.h>

#include "ptypes.h"
#ifdef STANDALONE_SIMPQS
# include "gmp_main.h"
# define UV_MAX ULONG_MAX
# define UVCONST(x) ((unsigned long)x##UL)
# define New(id, mem, size, type)  mem = (type*) malloc((size)*sizeof(type))
# define Newz(id, mem, size, type) mem = (type*) calloc(size, sizeof(type))
# define Safefree(mem)             free((void*)mem)
#else
# include "simpqs2.h"
#endif

#include "prime_iterator.h"
#include "utility.h"
#include "misc_ui.h"
#include "rootmod.h"
#include "lanczos.h"

typedef struct qs_factor_array_s qs_factor_array_t;
static void insert_factor(qs_factor_array_t *fa, const mpz_t f);

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

/* Trial division is useful for direct callers, but factor() has already
 * removed substantially more than this. */
#define SIMPQS_TRIAL_LIMIT 1000U

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

/* Poor man's random number generator. It satisfies no particularly good
 * randomness properties, but is good enough for this application
 */
static unsigned long randval = 4035456057U;
static unsigned long init_randval;   /* for diagnostics */
static unsigned long simple_random(unsigned long upto) {
  randval = (unsigned long)(((uint64_t)randval * 1025416097U + 286824428U)
                           % 4294967291U);
  return upto ? (randval % upto) : randval;
}

/* Insert an occupied row into a matrix column. */
static INLINE void insertColEntry(
  la_col_t *colarray, unsigned long colNum, unsigned long entry
) {
  unsigned long i;
  unsigned long *temp;

  if ((colarray[colNum].weight & 0x0f) == 0) {
    temp = colarray[colNum].data;
    colarray[colNum].data = (unsigned long *)malloc(
        (colarray[colNum].weight + 16) * sizeof(unsigned long));
    for (i = 0; i < colarray[colNum].weight; ++i)
      colarray[colNum].data[i] = temp[i];
    free(temp);
  }

  colarray[colNum].data[colarray[colNum].weight] = entry;
  ++colarray[colNum].weight;
  colarray[colNum].orig = colNum;
}

/* Toggle an occupied row in a matrix column. */
static INLINE void xorColEntry(
  la_col_t *colarray, unsigned long colNum, unsigned long entry
) {
  unsigned long i, j;

  for (i = 0; i < colarray[colNum].weight; ++i) {
    if (colarray[colNum].data[i] == entry) {
      for (j = i; j < colarray[colNum].weight - 1; ++j)
        colarray[colNum].data[j] = colarray[colNum].data[j + 1];
      --colarray[colNum].weight;
      return;
    }
  }
  insertColEntry(colarray, colNum, entry);
}

static INLINE void clearCol(la_col_t *colarray, unsigned long colNum) {
  colarray[colNum].weight = 0;
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
  unsigned int i;
  for (i = 0; i < ar->count; ++i)
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
  arel_t *arc;
  arel_t tmp;
  unsigned int ia = 0, ib = 0, ic = 0;
  rel_t *next_ra, *next_rb;
  mpz_t *lastp = NULL;

  arc = new_arel();
  resize_arel(arc, ara->count + arb->count);
  next_ra = ara->count ? ara->r[ia++] : NULL;
  next_rb = arb->count ? arb->r[ib++] : NULL;

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
  tmp = *ara;
  *ara = *arc;
  *arc = tmp;
  steal_arel(arc);
  return ara->count;
}

unsigned int combine_partial(arel_t *comb, rel_t *ra, rel_t *rb,
  unsigned long numPrimes, mpz_t n, qs_factor_array_t *factors
) {
  mpz_t X, X2;
  rel_t *r;
  fact_t *fa, *fb;
  unsigned int ia = 0, ib = 0;

  /* same partial, nothing to do */
  if (mpz_cmpabs(ra->X, rb->X) == 0)
    return 0;

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
    if (get_verbose_level() > 4)
      gmp_printf("# qs early factor %Zd\n", X);
    insert_factor(factors, X);
    mpz_clear(X);
    return 0;
  }
  mpz_mul(X, X, ra->X);
  mpz_mul(X, X, rb->X);
  mpz_mod(X, X, n);

  /* prefer the smaller of (X, n - X) */
  mpz_init(X2);
  mpz_sub(X2, n, X);
  if (mpz_cmp(X2, X) < 0)
    mpz_set(X, X2);
  mpz_clear(X2);

  r = new_rel();
  mpz_abs(r->X, X);
  fa = &ra->f;
  fb = &rb->f;
  resize_fact(&r->f, fa->count + fb->count);

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
  unsigned long numPrimes, mpz_t n, qs_factor_array_t *factors
) {
  arel_t *arc;
  arel_t tmp;
  unsigned int i, ia = 0, ib = 0, ic = 0;
  rel_t *ra, *rb, *rb2;

  arc = new_arel();
  /* make room for worst case, if no full relations are extracted */
  resize_arel(arc, ara->count + arb->count);
  ra = ara->count ? ara->r[ia++] : NULL;
  rb = arb->count ? arb->r[ib++] : NULL;

  while (rb) {
    if (ra) {
      if (ra->Q < rb->Q) {
        arc->r[ic++] = ra;
        ra = (ia < ara->count) ? ara->r[ia++] : NULL;
        continue;
      }
      if (ra->Q == rb->Q) {
        while (rb && ra->Q == rb->Q) {
          combine_partial(comb, ra, rb, numPrimes, n, factors);
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
      rb2 = arb->r[ib++];
      while (rb2 && rb->Q == rb2->Q) {
        combine_partial(comb, rb, rb2, numPrimes, n, factors);
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
  tmp = *ara;
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
  unsigned int i, j;
  rel_t *r;
  spp_t *f;
#ifdef ERRORS
  mpz_t test1, test2;
  mpz_init(test1);
  mpz_init(test2);
#endif

  for (i = 0; i < full->count && relsFound < relSought; ++i) {
    r = full->r[i];
    for (j = 0; j < r->f.count; ++j) {
      f = &r->f.fact[j];
      if (f->e & 1)
        xorColEntry(colarray, relsFound, (unsigned long)f->p);
    }

#ifdef ERRORS
    mpz_set_ui(test1,1);
    for (j = 0; j < r->f.count; ++j) {
      mpz_set_ui(test2, factorBase[r->f.fact[j].p]);
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
        for (j = 0; j < r->f.count; ++j) {
          f = &r->f.fact[j];
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
#ifdef ERRORS
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
  const unsigned char *sieve,
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
  int ai;
  unsigned int k;
  unsigned int exponent;
  unsigned int extra;
  unsigned int score;
  unsigned char vv;
  unsigned int modp;
  const unsigned long *sieve2 = (const unsigned long *)sieve;
  size_t bits;
  rel_t *rel = NULL;  /* new relation to add */

  i = 0;
  j = 0;
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
        if (factorBase[0] == 2) {
          exponent = mpz_scan1(res, 0);
          mpz_tdiv_q_2exp(res, res, exponent);
        } else {
          mpz_set_ui(temp, factorBase[0]);
          exponent = mpz_remove(res, res, temp);
        }
        add_factor(rel, 0, exponent);
        extra += exponent * primeSizes[0];
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
        extra += exponent * primeSizes[1];
      }

      for (k = 2; k < firstprime; ++k) {
        modp = (i + ctimesreps) % factorBase[k];
        if (soln2[k] != (unsigned long)-1) {
          if (modp == soln1[k] || modp == soln2[k]) {
            mpz_set_ui(temp, factorBase[k]);
            exponent = mpz_remove(res, res, temp);
            CHECK_EXPONENT(exponent, k);
            PRINT_FB(exponent, k);
            add_factor(rel, k, exponent);
            extra += exponent * primeSizes[k];
          }
        } else {
          mpz_set_ui(temp, factorBase[k]);
          exponent = mpz_remove(res, res, temp);
          if (exponent) {
            PRINT_FB(exponent, k);
            add_factor(rel, k, exponent);
            extra += exponent * primeSizes[k];
          }
        }
      }

      score = (unsigned int)sieve[i] + extra;

      if (score >= bits) {
        vv = (unsigned char)1 << (i & 7);
        for (k = firstprime; k < secondprime && extra < score; ++k) {
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

        for (k = secondprime; k < numPrimes && extra < score; ++k) {
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
            for (ai = 0; ai < s; ++ai)
              add_factor(rel, aind[ai] + min, 1);
            sort_fact(&rel->f);
            mpz_abs(rel->X, temp3);
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
              for (ai = 0; ai < s; ++ai)
                add_factor(rel, aind[ai] + min, 1);
              sort_fact(&rel->f);
              mpz_abs(rel->X, temp3);
              rel->Q = mpz_get_ui(res);
              save_rel(lpnew, rel);
              rel = NULL;
            }
          } else {
#ifdef RELPRINT
            printf("....R\n");
#endif
            for (ai = 0; ai < s; ++ai)
              add_factor(rel, aind[ai] + min, 1);
            sort_fact(&rel->f);
            mpz_abs(rel->X, temp3);
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
    /* Factors of poly A for our input sizes are normally far below midprime.
     * Keep this guard as it's otherwise harmless. */
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

  Keep the factors as a dynamically grown multiplicative partition of the
  original input.  A factor discovered by trial division or linear algebra
  refines existing entries; it is not itself additional multiplicative mass.
  This also handles non-square-free inputs by allowing repeated factors.

============================================================================*/

#define FACTOR_ARRAY_GROW 16

struct qs_factor_array_s {
  mpz_t *factors;
  uint32_t count;
  uint32_t alloc;
};

static void factor_array_init(qs_factor_array_t *fa, const mpz_t n) {
  fa->count = 1;
  fa->alloc = FACTOR_ARRAY_GROW;
  New(0, fa->factors, fa->alloc, mpz_t);
  if (fa->factors == NULL)
    croak("SIMPQS2: Unable to allocate factor array\n");
  mpz_init_set(fa->factors[0], n);
}

static void factor_array_append(qs_factor_array_t *fa, const mpz_t f) {
  if (fa->count >= fa->alloc) {
    fa->alloc += FACTOR_ARRAY_GROW;
    Renew(fa->factors, fa->alloc, mpz_t);
    if (fa->factors == NULL)
      croak("SIMPQS2: Unable to grow factor array\n");
  }
  mpz_init_set(fa->factors[fa->count++], f);
}

/* Refine every applicable entry using the known divisor f.  The loop includes
 * newly appended quotients, so one prime divisor separates all of its copies.
 */
static void insert_factor(qs_factor_array_t *fa, const mpz_t f) {
  uint32_t i;
  mpz_t g, q;

  if (mpz_cmp_ui(f, 1) <= 0)
    return;

  mpz_init(g);
  mpz_init(q);
  for (i = 0; i < fa->count; ++i) {
    mpz_gcd(g, fa->factors[i], f);
    if (mpz_cmp_ui(g, 1) <= 0 || mpz_cmp(g, fa->factors[i]) == 0)
      continue;
    mpz_divexact(q, fa->factors[i], g);
    mpz_set(fa->factors[i], g);
    factor_array_append(fa, q);
  }
  mpz_clear(g);
  mpz_clear(q);
}

/* Verify that factor refinement hasn't changed the represented product. */
static void verify_factor_array(const mpz_t n, const qs_factor_array_t *fa) {
  uint32_t i;
  mpz_t t;
  mpz_init_set_ui(t, 1);
  for (i = 0; i < fa->count; ++i) {
    if (mpz_cmp_ui(fa->factors[i], 1) <= 0)
      croak("SIMPQS2: invalid factor array entry");
    mpz_mul(t, t, fa->factors[i]);
  }
  if (mpz_cmp(t, n) != 0) {
    gmp_printf("factor array doesn't multiply: n=%Zd t=%Zd\n", n, t);
    croak("assert");
  }
  mpz_clear(t);
}

static int allprime_factor_array(const qs_factor_array_t *fa) {
  uint32_t i;
  for (i = 0; i < fa->count; ++i)
    if (!mpz_probab_prime_p(fa->factors[i], 5))   /* Be lazy */
      return 0;
  return 1;
}

static mpz_t *factor_array_release(qs_factor_array_t *fa,
                                   uint32_t *nfactors) {
  *nfactors = fa->count;
  return fa->factors;
}

void _GMP_simpqs2_free(mpz_t *factors, uint32_t nfactors) {
  uint32_t i;
  for (i = 0; i < nfactors; ++i)
    mpz_clear(factors[i]);
  Safefree(factors);
}

/*============================================================================
   mainRoutine:

   Function: Generates the polynomials, initialises and calls the sieve,
             implementing cache blocking (breaking the sieve interval into
             small blocks for the small primes).

============================================================================*/
static void mainRoutine(
  unsigned long numPrimes,
  unsigned long Mdiv2,
  unsigned long relSought,
  mpz_t n,
  qs_factor_array_t *factors,
  const mpz_t original_n,
  unsigned long multiplier
) {
  mpz_t A, B, C, D, Bdivp2, nsqrtdiv, temp, temp2, temp3, temp4;
  int i, j, k, l, s, fact, span, min, verbose, polyindex;
  uint64_t mask;
  uint32_t lanczos_seed1, lanczos_seed2;
  unsigned long u1, p, reps, M, Mq, Mr, ncols, nrows, relsFound;
  unsigned long curves = 0;
  uint64_t *nullrows;
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
#ifdef ERRORS
  mpz_t test1, test2, test3;
  unsigned int *exps;
  unsigned int index;
  rel_t *r;
#endif

  unsigned long next_cutoff = (relSought - 1) / 40 + 1;
  unsigned long next_inc = next_cutoff;

  arel_t *frels = new_arel();   /* all full relations */
  arel_t *rels = new_arel();    /* new full relations, to be merged */
  arel_t *lprels = new_arel();  /* all partial relations */
  arel_t *lpnew = new_arel();   /* new partial relations, to be merged */

  verbose = get_verbose_level();
  s = mpz_sizeinbase(n, 2) / 28 + 1;
  lanczos_seed1 = 11111111U ^ (uint32_t)randval
                ^ (uint32_t)mpz_fdiv_ui(n, 4294967291UL);
  lanczos_seed2 = 22222222U ^ ~(uint32_t)randval
                ^ (uint32_t)mpz_fdiv_ui(n, 4294967279UL);

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

  while (frels->count < relSought) {
    mpz_set_ui(A, 1);
    for (i = 0; i < s - 1; ) {
      unsigned long ran = span / 2 + simple_random(span / 2);
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
            - simple_random(10) - min;
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
      amodp[i] = mpz_fdiv_ui(temp, p);
      mpz_set_ui(temp, modinverse(amodp[i], p));
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

    for (polyindex = 1; polyindex < (1 << (s - 1)) - 1; ++polyindex) {
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
        /* full relations found while merging are extracted into rels */
        merge_partial(lprels, lpnew, rels, numPrimes, n, factors);
        reset_arel(lpnew);

        /* combine everything into the main list */
        sort_rel(rels);
        merge_full(frels, rels);
        reset_arel(rels);

#ifdef COUNT
        printf("%u full, %u partial\n", frels->count, lprels->count);
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

  ncols = relSought;
  nrows = numPrimes;

#ifdef ERRORS
  for (j = frels->count; j < relSought; ++j)
    if (colarray[j].weight != 0)
      printf("Dirty at col %d\n", j);
#endif

#ifdef COUNT
  printf("%u relations found in total!\n", frels->count);
#endif

  relsFound = read_matrix(frels, colarray, 0, relSought, n);

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

  mpz_init(test1);
  mpz_init(test2);
  mpz_init(test3);
  exps = (unsigned int *)malloc(numPrimes * sizeof(unsigned int));
  for (j = 0; j < relSought; ++j) {
    for (i = 0; i < numPrimes; ++i)
      exps[i] = 0;
    r = frels->r[j];
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
    index = colarray[j].orig;
    r = frels->r[index];
    for (i = 0; i < r->f.count; ++i)
      exps[r->f.fact[i].p] += r->f.fact[i].e;
    for (i = 0; i < colarray[j].weight; ++i) {
      for (k = 0; k < i; ++k)
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

  nullrows = block_lanczos(
    nrows, 0, ncols, colarray, lanczos_seed1, lanczos_seed2, &mask
  );
  if (nullrows == NULL) {
    gmp_printf(
      "block_lanczos failed repeatedly on target %Zd (multiplier %d) from randval %lu, giving up",
      n, multiplier, init_randval
    );
    croak("assert");
  }

  if (verbose > 3) {
    for (i = j = 0; i < 64; ++i)
      if (mask & ((uint64_t)1 << i))
        ++j;
    printf("%d nullspace vectors found.\n", j);
  }

#ifdef ERRORS
  exps = (unsigned int *)malloc(numPrimes * sizeof(unsigned int));
  for (j = 0; j < ncols; ++j) {
    for (i = 0; i < numPrimes; ++i)
      exps[i] = 0;
    index = colarray[j].orig;
    r = frels->r[index];
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

  /* Now refine the factor array via square root and gcd */
  New(0, primecount, numPrimes, unsigned int);
  if (primecount == 0)
    croak("SIMPQS: Unable to allocate memory!\n");
  for (l = 0; l < 64; ++l) {
    while (l < 64 && !(mask & ((uint64_t)1 << l)))
      ++l;
    if (l == 64)
      break;
    mpz_set_ui(temp, 1);
    mpz_set_ui(temp2, 1);
    memset(primecount, 0, numPrimes * sizeof(unsigned int));
    for (i = 0; i < ncols; ++i) {
      if (getNullEntry(nullrows, i, l)) {
        unsigned int index = colarray[i].orig;
        rel_t *r = frels->r[index];
        mpz_mul(temp2, temp2, r->X);
        for (j = 0; j < r->f.count; ++j)
          primecount[r->f.fact[j].p] += r->f.fact[j].e;
      }
      if (((i + 1) % 16) == 0)
        mpz_mod(temp2, temp2, n);
    }
    for (j = 0; j < numPrimes; ++j) {
      if (primecount[j]) {
        mpz_set_ui(temp3, factorBase[j]);
        mpz_pow_ui(temp3, temp3, primecount[j] / 2);
        mpz_mul(temp, temp, temp3);
      }
      if (((j + 1) % 16) == 0)
        mpz_mod(temp, temp, n);
    }
    mpz_sub(temp, temp2, temp);
    mpz_gcd(temp, temp, n);
    /* only non-trivial factors */
    if (mpz_cmp_ui(temp, 1) && mpz_cmp(temp, n)) {
      if (verbose > 4)
        gmp_printf("# qs factor %Zd\n", temp);
      insert_factor(factors, temp);
      verify_factor_array(original_n, factors);
      if (allprime_factor_array(factors))
        break;
    }
  }

  /* Free everything remaining */
  free(nullrows);
  free_arel(frels);
  free_arel(rels);
  free_arel(lprels);
  free_arel(lpnew);
  for (i = 0; i < relSought; ++i)
    free(colarray[i].data);
  Safefree(colarray);
  Safefree(primecount);
  mpz_clear(temp);  mpz_clear(temp2);  mpz_clear(temp3);  mpz_clear(temp4);

}

mpz_t *_GMP_simpqs2(const mpz_t n, uint32_t *nfactors,
                     uint32_t trial_start) {
  qs_factor_array_t factors;
  mpz_t nred;
  unsigned long numPrimes, Mdiv2, multiplier, decdigits, relSought;
  int found_trial_factor = 0;
  int verbose = get_verbose_level();

  if (nfactors == NULL)
    croak("SIMPQS2: Missing factor count output");
  factor_array_init(&factors, n);
  decdigits = mpz_sizeinbase(n, 10); /* often 1 too big */
  if (decdigits < MINDIG)
    return factor_array_release(&factors, nfactors);

  if (verbose > 2)
    gmp_printf("# qs trying %Zd (%lu digits)\n", n, decdigits);

  mpz_init_set(nred, n);

  /* Remove small factors not already checked by the caller. */
  if (trial_start < SIMPQS_TRIAL_LIMIT) {
    UV p, first = trial_start < 2 ? 2 : trial_start;
    mpz_t fp;
    PRIME_ITERATOR(iter);
    mpz_init(fp);
    prime_iterator_setprime(&iter, first - 1);
    for (p = prime_iterator_next(&iter);
         p < SIMPQS_TRIAL_LIMIT;
         p = prime_iterator_next(&iter)) {
      if (mpz_cmp_ui(nred, p * p) < 0) break;
      if (mpz_divisible_ui_p(nred, p)) {
        mpz_set_ui(fp, p);
        insert_factor(&factors, fp);
        mpz_remove(nred, nred, fp);
        found_trial_factor = 1;
      }
    }
    mpz_clear(fp);
    prime_iterator_destroy(&iter);
    verify_factor_array(n, &factors);
  }

  decdigits = mpz_sizeinbase(nred, 10);
  if (decdigits < MINDIG ||
      (found_trial_factor && allprime_factor_array(&factors))) {
    mpz_clear(nred);
    return factor_array_release(&factors, nfactors);
  }

  /* Get a preliminary number of primes, pick a multiplier, apply it */
  numPrimes = (decdigits <= 91) ? primesNo[decdigits - MINDIG] : 80000;
  multiplier = knuthSchroeppel(nred, numPrimes);
  mpz_mul_ui(nred, nred, multiplier);
  decdigits = mpz_sizeinbase(nred, 10);

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
    /* No configurations have been tuned beyond 91 digits.  The inherited
     * SIMPQS/FlintQS fallback reduced several parameters discontinuously;
     * retain the 91-digit endpoint instead. */
    numPrimes   = 80000;
    Mdiv2       = 192000 / SIEVEDIV;
    largeprime  = 500000000U;
    secondprime = SECONDPRIME;
    midprime    = MIDPRIME;
    firstprime  = 29;
    errorbits   = 33;
    threshold   = 102;
  }

  if (verbose > 2)
    gmp_printf("# qs    mult %lu, digits %lu, sieving %lu, primes %lu\n",
        multiplier, decdigits, Mdiv2 * 2, numPrimes);

  /* as numPrimes+64, was commented: we probably need fewer than this */
  relSought = numPrimes;
  initFactorBase();
  computeFactorBase(nred, numPrimes, multiplier);

  init_randval = randval;
  mainRoutine(
    numPrimes, Mdiv2, relSought, nred, &factors, n, multiplier
  );

  clearFactorBase();
  mpz_clear(nred);
  if (verbose > 2) {
    int i;
    gmp_printf("# qs:");
    for (i = 0; i < (int)factors.count; ++i)
      gmp_printf(" %Zd", factors.factors[i]);
    gmp_printf("%s\n", (factors.count > 1) ? "" : " no factors");
  }
  if (verbose > 2 && factors.count <= 1)
    gmp_printf("QS Fail: %Zd (%lu digits)\n", n, decdigits);
  return factor_array_release(&factors, nfactors);
}

#ifdef STANDALONE_SIMPQS
/* Main Program: factors a user specified number using a quadratic sieve */
int main(int argc, char **argv) {
  int i = 1;
  uint32_t nfactors;
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
    if (arg[1] == 'r') {
      randval = strtoul(&arg[2], NULL, 10);
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
  farray = _GMP_simpqs2(n, &nfactors, 2);

  for (i = 0; i < (int)nfactors; ++i)
    gmp_printf("  %Zd\n", farray[i]);
  _GMP_simpqs2_free(farray, nfactors);
  mpz_clear(n);
  _GMP_destroy();
}
#endif
