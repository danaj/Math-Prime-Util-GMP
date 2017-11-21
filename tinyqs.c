#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "tinyqs.h"

/******************************************************************************/
typedef signed char s8;
typedef unsigned char u8;
typedef signed short s16;
typedef unsigned short u16;
typedef signed int s32;
typedef unsigned int u32;

#ifdef _MSC_VER
typedef signed __int64 s64;
typedef unsigned __int64 u64;
/* #define INLINE _inline */
#else
typedef long long s64;
typedef unsigned long long u64;
/* #define INLINE inline */
#endif
static INLINE u64 gmp2u64(mpz_t src)
{
  /* mpz_export is terribly slow */
  u64 ans = mpz_getlimbn(src, 0);
#if GMP_LIMB_BITS == 32
  if (mpz_size(src) >= 2)
    ans |= (u64)mpz_getlimbn(src, 1) << 32;
#endif
  return ans;
}
static INLINE void u64_2gmp(u64 src, mpz_t dest)
{
#if GMP_LIMB_BITS == 64
  dest->_mp_d[0] = src;
  dest->_mp_size = (src ? 1 : 0);
#else
  /* mpz_import is terribly slow */
  mpz_set_ui(dest, (u32)(src >> 32));
  mpz_mul_2exp(dest, dest, 32);
  mpz_add_ui(dest, dest, (u32)src);
#endif
}
/******************************************************************************/

#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif

#ifndef M_SQRT2
#define M_SQRT2	1.41421356237309504880
#endif

/* High-throughput, low-overhead implementation of the self-
   initializing multiple polynomial quadratic sieve, optimized
   for small inputs (50-120 bits). Many of the ideas here are
   extensions of the remarkable MPQS code of F. Bahr, used 
   in the lattice sievers by Jens Franke */

/* TODO: lots of static mpz_t is an issue for threading */

/* seeds for random numbers */

static u32 rand_seed1 = 11111111;
static u32 rand_seed2 = 22222222;

#define RAND_MULT 2131995753

static u32 get_rand(u32 *seed1, u32 *seed2)
{
  /* A multiply-with-carry generator by George Marsaglia.
     The period is about 2^63. */

  u64 temp = (u64)(*seed1) * (u64)RAND_MULT + (u64)(*seed2);
  *seed1 = (u32)temp;
  *seed2 = (u32)(temp >> 32);
  return (u32)temp;
}

/* masks for picking out individual bits of 64-bit
   words, used for the linear algebra */

#define B(x) ((u64)(1) << (x))

static const u64 bitmask[] = {
  B( 0), B( 1), B( 2), B( 3), B( 4), B( 5), B( 6), B( 7),
  B( 8), B( 9), B(10), B(11), B(12), B(13), B(14), B(15),
  B(16), B(17), B(18), B(19), B(20), B(21), B(22), B(23),
  B(24), B(25), B(26), B(27), B(28), B(29), B(30), B(31),
  B(32), B(33), B(34), B(35), B(36), B(37), B(38), B(39),
  B(40), B(41), B(42), B(43), B(44), B(45), B(46), B(47),
  B(48), B(49), B(50), B(51), B(52), B(53), B(54), B(55),
  B(56), B(57), B(58), B(59), B(60), B(61), B(62), B(63),
};

/* maximum size pool of primes from which 
   factor base is constructed */
#define NUM_PRIMES_TINY 1024

/* the number of dependencies the linear algebra
   will find */
#define NUM_EXTRA_RELATIONS_TINY 16

/* largest number of relations that can go into the
   linear algebra (includes relations combined from
   pairs of partial relations */
#define MAX_RELATIONS_TINY 512

/* the largest possible factor base */
#define MAX_FB_SIZE_TINY (MAX_RELATIONS_TINY - \
                          NUM_EXTRA_RELATIONS_TINY)

/* offset of the first valid factor base prime */
#define MIN_FB_OFFSET_TINY 1

/* offset of the first factor base prime 
   actually contributing to the sieving */
#define MIN_FB_OFFSET_TO_SIEVE_TINY 7

/* number of primes used when testing multipliers */
#define NUM_TEST_PRIMES_TINY 30

/* fudge factor to the target sieve value to account
   for not sieving with the smallest factor base primes */
#define SMALL_PRIME_FUDGE_TINY 10

/* maximum number of MPQS polynomials to be computed */
#define MAX_POLY_TINY 256

/* maximum number of FB primes that contribute to 
   a single polynomial 'A' value */
#define MAX_POLY_FACTORS_TINY 5

/* the size of the sieve interval. Each polynomial will
   sieve over this many positive and negative values */
#define SIEVE_SIZE_TINY 16384

/* value of the sieve root used when sieving is not
   to be performed for a given FB prime. Since this is
   larger than SIEVE_SIZE_TINY no special-case code
   is needed in the core sieve code */
#define DO_NOT_SIEVE_TINY 65535

/* maximum number of factors a relation can have (the
   large prime is stored separately) */
#define MAX_FACTORS_TINY 20

/* partial relations are listed in the order 
   in which they occur, and a hashtable matches 
   up partial relations with the same large prime. */
#define LOG2_PARTIAL_TABLE_SIZE 10
#define LARGE_PRIME_HASH(x) (((u32)(x) * ((u32)40499 * 65543)) >> \
                                (32 - LOG2_PARTIAL_TABLE_SIZE))

/* number of collisions allowed in one hashtable entry */
#define LP_HASH_DEPTH_TINY 3

/* scale factor for all log values */
#define LOGPRIME_SCALE_TINY 2

/* maximum number of relations to be saved for 
   resieving, used in place of trial factoring */
#define SIEVE_BATCH_SIZE_TINY 128

/* maximum size of the pool of FB primes that 
   can appear in a polynomial 'A' value */
#define POLY_SELECT_BITS_TINY 12

#define POSITIVE 0
#define NEGATIVE 1

/* structure describing a single relation */

typedef struct {
  u32 large_prime;      /* the large prime (may be 1) */
  s16 sieve_offset;     /* the sieve offset of the relation */
  u8 poly_num;          /* ID of the poly that produce the relation */
  u8 num_factors;       /* number of factors from the factor base
                           (duplicates count) */
  u16 fb_offsets[MAX_FACTORS_TINY]; /* offsets into FB of primes that
                                       divide this relation */
} tiny_relation;

/* structure describing a factor base entry */

typedef struct {
  u16 prime;    /* the factor base prime */
  u16 modsqrt;  /* x that solves x^2 = N mod p */
  u32 recip;    /* integer reciprocal of 'prime' */
  u8 logprime;  /* log value used in sieve */
  u16 roots[2]; /* the two sieve roots for 'prime' */
} tiny_fb;

/* structure describing one SIQS polynomial */

typedef struct {
  u16 a_fb_offsets[MAX_POLY_FACTORS_TINY];  /* factors of 'A' value */
  mpz_t b;                                  /* B value */
} tiny_poly;

/* main structure controlling the factorization */

typedef struct {

  /* basic stuff */

  mpz_t n;                          /* number to be factored */
  u32 multiplier;                   /* small multiplier of n */
  u16 multiplier_fb[2];             /* fb offsets of factors of multiplier */

  /* polynomial selection stuff */

  double target_a;                  /* the optimal size of poly A values */
  s32 poly_num;                     /* ID of current polynomial */
  s32 num_a_factors;                /* # of factors in poly 'A' values */
  s32 poly_select_idx;              /* ID of the combination of primes
                                       that will make current A value */
  u16 poly_select_offsets[POLY_SELECT_BITS_TINY]; /* pool of primes for A */
  mpz_t poly_b_aux[MAX_POLY_FACTORS_TINY];      /* scratch values for com-
                                                   puting poly B values */
  tiny_poly poly_list[MAX_POLY_TINY];      /* list of SIQS polynomials */

  /* sieve stuff */

  double align_me;
  u8 sieve_block[SIEVE_SIZE_TINY];  /* the sieve interval (8-byte aligned) */

  /* factor base stuff */

  s32 fb_size;                      /* number of FB primes */
  u16 prime_list[NUM_PRIMES_TINY];  /* complete list of primes from which
                                       factor base is generated */
  float test_prime_contrib[NUM_TEST_PRIMES_TINY]; /* scratch space used in 
                                                     multiplier selection */
  tiny_fb factor_base[MAX_FB_SIZE_TINY];          /* the factor base */
  u16 root_aux[MAX_POLY_FACTORS_TINY * 
               MAX_FB_SIZE_TINY];      /* scratch value for initializing
                                          sieve roots */
  /* relation stuff */

  s32 num_full_relations;   /* where next full relation will go */
  s32 partial_idx;          /* where next partial relation will go */
  s32 large_prime_max;      /* max value of a large prime */
  s32 error_bits;           /* value used for trial factoring cutoff */
  tiny_relation sieve_batch[SIEVE_BATCH_SIZE_TINY]; /* resieved relations */

  /* all relations that survive sieving are put in relation_list.
     Full relations (and partial relations whose large prime has
     occurred more than once) are stored in a list that grows up
     from the beginning of the list, while partial relations that
     have not been matched up yet are stored in a list growing down
     from the end of relation_list. num_full_relations is the index
     of the first free space for full relations, and partial_idx 
     does the same for unmatched partial relations. */

  tiny_relation relation_list[4 * MAX_RELATIONS_TINY];

  /* a hashtable is used to match up partial relations, using the
     large prime as a hash key. The hashtable stores the index in
     relation_list of the partial relation that connects up all the
     other partial relations with the same large prime (those other
     relations are treated as full relations) */

  u16 partial_hash[1 << LOG2_PARTIAL_TABLE_SIZE][LP_HASH_DEPTH_TINY];

  /* linear algebra stuff */

  u16 null_vectors[MAX_RELATIONS_TINY];
  u64 matrix[MAX_FB_SIZE_TINY][(MAX_RELATIONS_TINY+63) / 64];
} tiny_qs_params;


/* the following is reused across factorizations */

static tiny_qs_params *g_params = NULL;

/* The following utility routines are not really
   a performance bottleneck, but since they always
   deal with 16-bit data at most their input 
   datatypes should really by u16's. This will 
   make all the division and remainder operations 
   a lot faster */

/***********************************/
static s32 legendre_16(s32 a, s32 p)
/***********************************
Compute the Legendre symbol (a/p)
************************************/
{
  s32 tmp;
  s32 x = a;
  s32 y = p;
  s32 out = 1;

  while (x) {
    while ((x & 1) == 0) {
      x = x / 2;
      if ( (y & 7) == 3 || (y & 7) == 5 )
        out = -out;
    }

    tmp = x;
    x = y;
    y = tmp;

    if ( (x & 3) == 3 && (y & 3) == 3 )
      out = -out;

    x = x % y;
  }
  if (y == 1)
    return out;
  return 0;
}

/***********************************/
static s32 powm_16(s32 a, s32 b, s32 n) 
/***********************************
Compute a^b mod n
************************************/
{
  s32 res = 1;
  while (b) {
    if (b & 1)
      res = res * a % n;
    a = a * a % n;
    b = b >> 1;
  }
  return res;
}

/***********************************/
static s32 modinv_16(s32 a, s32 p)
/***********************************
High-speed modular inverse of 'a' mod 'p'
Thanks to the folks at www.mersenneforum.com
for coming up with this
************************************/
{
  s32 ps1, ps2, parity, dividend, divisor, rem, q, t;

  q = 1;
  rem = a;
  dividend = p;
  divisor = a;
  ps1 = 1;
  ps2 = 0;
  parity = 0;

  while (divisor > 1) {
    rem = dividend - divisor;
    t = rem - divisor;
    if (t >= 0) {
      q += ps1; rem = t; t -= divisor;
      if (t >= 0) {
        q += ps1; rem = t; t -= divisor;
        if (t >= 0) {
          q += ps1; rem = t; t -= divisor;
          if (t >= 0) {
            q += ps1; rem = t; t -= divisor;
            if (t >= 0) {
              q += ps1; rem = t; t -= divisor;
              if (t >= 0) {
                q += ps1; rem = t; t -= divisor;
                if (t >= 0) {
                  q += ps1; rem = t; t -= divisor;
                  if (t >= 0) {
                    q += ps1; rem = t;
                    if (rem >= divisor) {
                      q = dividend / divisor;
                      rem = dividend % divisor;
                      q *= ps1;
                    } } } } } } } }
    }
    q += ps2;
    parity = ~parity;
    dividend = divisor;
    divisor = rem;
    ps2 = ps1;
    ps1 = q;
  }
  
  if (parity == 0)
    return ps1;
  else
    return p - ps1;
}

/***********************************/
static s32 sqrtModP_16(s32 a, s32 p)
/***********************************
Compute the square root of 'a' mod 'p'
This is Algorithm 2.3.8 from Crandall & 
Pomerance, "Prime Numbers: A Computational
Perspective"
************************************/
{
  if ( (p & 7) == 3 || (p & 7) == 7 ) {
    return powm_16(a, (p+1)/4, p);
  }
  else if ( (p & 7) == 5 ) {
#if 0
    s32 x, y;
    
    x = powm_16(a, (p+3)/8, p);
    if ((x * x % p) == a)
      return x;

    y = powm_16(2, (p-1)/4, p);
    return (s32)x * y % p;
#else
#define mulm_16(a, b, n)  (((a) * (b)) % (n))
    u32 a2, alpha, beta, b;
    a2 = (a+a) % p;
    alpha = powm_16(a2, (p-5)>>3, p);
    beta = mulm_16(a2, mulm_16(alpha,alpha,p), p);
    b = mulm_16(alpha, mulm_16(a, (beta ? beta-1 : p-1), p), p);
    return b;
  }
  else if ( (p & 16) == 9 ) {
    s32 a2, alpha, beta, b, d = 1;
    a2 = (a+a) % p;
    alpha = powm_16(a2, (p-9)>>4, p);
    beta = mulm_16(a2, mulm_16(alpha,alpha,p), p);
    if (((beta*beta) % p) != p-1) {
      do { d += 2; } while (legendre_16(d,p) != -1 && d < p);
      alpha = mulm_16(alpha, powm_16(d, (p-9)>>3, p), p);
      beta = mulm_16(a2, mulm_16(mulm_16(d,d,p),mulm_16(alpha,p,p),p), p);
    }
    b = mulm_16(alpha, mulm_16(a, mulm_16(d, (beta ? beta-1 : p-1), p), p), p);
    return b;
#endif
  }
  else {
    s32 i, d0, d1, a1, s, t, m;

    d0 = get_rand(&rand_seed1, &rand_seed2) % p;
    while (legendre_16(d0, p) != -1)
      d0 = get_rand(&rand_seed1, &rand_seed2) % p;

    t = p - 1;
    s = 0;
    while (!(t & 1)) {
      s++;
      t = t / 2;
    }

    a1 = powm_16(a, t, p);
    d1 = powm_16(d0, t, p);

    for (i = 0, m = 0; i < s; i++) {
      s32 ad = powm_16(d1, m, p);
      ad = ad * a1 % p;
      ad = powm_16(ad, (u16)(1) << (s-1-i), p);
      if (ad == (p - 1))
        m += (1 << i);
    }

    a1 = powm_16(a, (t+1)/2, p);
    d1 = powm_16(d1, m/2, p);
    return a1 * d1 % p;
  }
}


/***********************************/
static void init_tinyqs(void)
/***********************************/
{
  s32 i, j, k, rem;
  tiny_qs_params *p;

  if (g_params)
    return;

  /* allocate the main structure */

  p = g_params = (tiny_qs_params *)malloc(sizeof(tiny_qs_params));
  mpz_init(p->n);

  /* fill in the pool of primes */

  p->prime_list[0] = 2;
  p->prime_list[1] = 3;
  for (i = 2, j = 5; i < NUM_PRIMES_TINY; j += 2) {
    for (k = 1, rem = 0; k < i; k++) {
      s32 prime = p->prime_list[k];
      rem = j % prime;
      if (prime * prime > j || rem == 0)
        break;
    }
    if (rem != 0)
      p->prime_list[i++] = j;
  }

  /* init the scratch values for polynomial 'B'
     value computations */

  for (i = 0; i < MAX_POLY_FACTORS_TINY; i++) {
    mpz_init(p->poly_b_aux[i]);
  }

  /* set up the list of sieve polynomials */

  for (i = 0; i < MAX_POLY_TINY; i++) {
    mpz_init(p->poly_list[i].b);
  }

  /* see the next routine for an explanation of what
     these quantities are */

  for (i = 1; i < NUM_TEST_PRIMES_TINY; i++) {
    p->test_prime_contrib[i] = 2 * log((double)p->prime_list[i]) / 
                               (p->prime_list[i] - 1) / M_LN2;
  }
}

/* Implementation of the modified Knuth-Schroeppel multiplier
   algorithm. This borrows ideas from at least four different
   sources, and seems to choose multipliers that are better on
   average than many of the other methods available.
   
   There are many misconceptions about what this algorithm is
   supposed to do. We want to multiply the input number n by a
   small odd squarefree constant k, chosen so that the factor base 
   for k * n contains as many small primes as possible. Since small primes
   occur more often than big ones, this makes sieve values smaller
   on average and so more likely to be smooth. We quantify this
   by measuring the average contribution of the first NUM_TEST_PRIMES_TINY
   primes to sieve values. There are two constraints: first, larger 
   multipliers mean a larger number to factor. Second, we can't spend 
   all day testing multipliers, so the set of multipliers to test should 
   be small. 

   The list of available multipliers depends on the value of n mod
   8, 3, and 5; each row of the table below gives the multipliers
   to try, pre-sorted by how well they approximately optimize sieving
   (the routine below computes a better approximation). Note that a
   multiplier of 1 (i.e. no multiplier) is always possible. Experiments
   show that 90% of the time the optimal multiplier is in one of the 
   first four columns of the table */

#define MAX_MULTIPLIERS 13                           /* for residue classes: */
static u8 mult_list[32][MAX_MULTIPLIERS] = {         /* mod 8  mod 3  mod 5 */
{ 1, 19, 61, 31, 21, 13,  7,  3, 73, 41,  5, 33, 37 }, /*  1      1      1 */
{ 1, 13,  7,  3, 73, 33, 37, 17, 57, 43,  5, 19, 15 }, /*  1      1      2 */
{ 1, 13,  7,  3, 73, 33, 37, 17, 57, 43,  5, 19, 15 }, /*  1      1      3 */
{ 1, 19, 61, 31, 21, 13,  7,  3, 73, 41,  5, 33, 37 }, /*  1      1      4 */
{ 1, 41,  5, 17, 11, 89, 29, 65, 21,  3, 59, 33, 35 }, /*  1      2      1 */
{ 1, 17,  5,  3, 33, 65, 57, 23, 41, 53, 47, 11, 89 }, /*  1      2      2 */
{ 1, 17,  5,  3, 33, 65, 57, 23, 41, 53, 47, 11, 89 }, /*  1      2      3 */
{ 1, 41,  5, 17, 11, 89, 29, 65, 21,  3, 59, 33, 35 }, /*  1      2      4 */
{ 1, 19,  3, 11, 31,  7, 51, 43, 15, 39, 61, 55, 21 }, /*  3      1      1 */
{ 1,  3,  7, 43, 19, 13, 37, 15, 55, 11, 73, 31, 35 }, /*  3      1      2 */
{ 1,  3,  7, 43, 19, 13, 37, 15, 55, 11, 73, 31, 35 }, /*  3      1      3 */
{ 1, 19,  3, 11, 31,  7, 51, 43, 15, 39, 61, 55, 21 }, /*  3      1      4 */
{ 1, 11,  3, 59, 35,  5, 51, 19, 29, 41, 15, 23, 39 }, /*  3      2      1 */
{ 1,  3, 11, 35,  5, 23, 17, 47,  7, 59, 43, 15, 53 }, /*  3      2      2 */
{ 1,  3, 11, 35,  5, 23, 17, 47,  7, 59, 43, 15, 53 }, /*  3      2      3 */
{ 1, 11,  3, 59, 35,  5, 51, 19, 29, 41, 15, 23, 39 }, /*  3      2      4 */
{ 1, 61, 21, 13,  5, 19, 37, 31, 29,  7,  3, 11, 15 }, /*  5      1      1 */
{ 1, 13, 37,  7,  3,  5, 73, 61, 21, 43, 33, 53, 17 }, /*  5      1      2 */
{ 1, 13, 37,  7,  3,  5, 73, 61, 21, 43, 33, 53, 17 }, /*  5      1      3 */
{ 1, 61, 21, 13,  5, 19, 37, 31, 29,  7,  3, 11, 15 }, /*  5      1      4 */
{ 1,  5, 29, 21, 11, 41, 53, 17, 89,  3, 59, 61, 65 }, /*  5      2      1 */
{ 1,  5, 53, 17,  3, 13, 29, 23, 21, 37, 47, 33, 11 }, /*  5      2      2 */
{ 1,  5, 53, 17,  3, 13, 29, 23, 21, 37, 47, 33, 11 }, /*  5      2      3 */
{ 1,  5, 29, 21, 11, 41, 53, 17, 89,  3, 59, 61, 65 }, /*  5      2      4 */
{ 1, 31,  7, 19, 15, 39, 55,  3, 11, 61, 21, 13, 51 }, /*  7      1      1 */
{ 1,  7,  3, 15, 13, 55, 31, 43, 23, 37, 19, 47, 73 }, /*  7      1      2 */
{ 1,  7,  3, 15, 13, 55, 31, 43, 23, 37, 19, 47, 73 }, /*  7      1      3 */
{ 1, 31,  7, 19, 15, 39, 55,  3, 11, 61, 21, 13, 51 }, /*  7      1      4 */
{ 1, 11,  5, 15, 23, 39,  3, 29, 47, 59, 31, 35,  7 }, /*  7      2      1 */
{ 1, 23,  3, 47,  7,  5, 15, 17, 11, 35, 53, 39, 33 }, /*  7      2      2 */
{ 1, 23,  3, 47,  7,  5, 15, 17, 11, 35, 53, 39, 33 }, /*  7      2      3 */
{ 1, 11,  5, 15, 23, 39,  3, 29, 47, 59, 31, 35,  7 }, /*  7      2      4 */
};

/***********************************/
static void find_multiplier_tiny(void)
/***********************************/
{
  tiny_qs_params *params = g_params;
  s32 i, j;
  u16 *prime_list = params->prime_list;
  u16 test_nmodp[NUM_TEST_PRIMES_TINY];
  s32 best_mult = 1;
  s32 nmod8 = mpz_get_ui(params->n) % 8;
  float best_score;
  u8 *mult_row;
  s32 num_tests;

  /* precompute information that will be needed 
     for all multipliers */

  for (i = 1; i < NUM_TEST_PRIMES_TINY; i++)
    test_nmodp[i] = mpz_tdiv_ui(params->n, prime_list[i]);

  /* find the row of the table that is approriate for this
     value of n */

  mult_row = mult_list[ test_nmodp[2] - 1 +
                        4*(test_nmodp[1] - 1) + 
		        8*(nmod8 / 2) ];

  /* test less than the whole row if n is small */

  num_tests = mpz_sizeinbase(params->n, 2) / 10;
  if (num_tests > MAX_MULTIPLIERS)
    num_tests = MAX_MULTIPLIERS;

  best_score = 1000.0;
  for (i = 0; i < num_tests; i++) {
    s32 curr_mult = mult_row[i];
    s32 knmod8 = (nmod8 * curr_mult) % 8;
    float score;

    /* measure the contribution of 2 as a factor of sieve
       values. The multiplier itself must also be taken into
       account in the score. 'score' is the correction that
       is implicitly applied to the size of sieve values; a
       negative score makes sieve values smaller, and so is 
       better. */

    if (knmod8 == 1)
      score = 0.5 * log((double)curr_mult) / M_LN2 - 2;
    else if (knmod8 == 5)
      score = 0.5 * log((double)curr_mult) / M_LN2 - 1;
    else
      score = 0.5 * log((double)curr_mult) / M_LN2 - 0.5;

    for (j = 1; j < NUM_TEST_PRIMES_TINY; j++) {
      s32 prime = prime_list[j];
      s32 knmodp = (s32)test_nmodp[j] * curr_mult % prime;

      /* if prime j is actually in the factor base 
         for k * n ... */

      if (legendre_16(knmodp, prime) != -1) {

        /* ...add its contribution. A prime p con-
           tributes log(p) to 1 in p sieve values, plus
           log(p) to 1 in p^2 sieve values, etc. The
           average contribution of all multiples of p 
           to a random sieve value is thus

           log(p) * (1/p + 1/p^2 + 1/p^3 + ...)
           = (log(p) / p) * 1 / (1 - (1/p)) 
           = log(p) / (p-1)

           This contribution occurs once for each
           square root used for sieving. There are two
           roots for each factor base prime, unless
           the prime divides the multiplier. In that
           case there is only one root. The scores are
           premultiplied by 2.0, and logarithms are 
           in base 2 (though any base will do) */

        if (knmodp == 0)
          score -= 0.5 * params->test_prime_contrib[j];
        else
          score -= params->test_prime_contrib[j];
      }
    }
    if (score < best_score) {
      best_score = score;
      best_mult = curr_mult;
    }
  }

  /* from now on we will factor best_mult * n */

  params->multiplier = best_mult;
  mpz_mul_ui(params->n, params->n, best_mult);
}

/***********************************/
static s32 init_fb_tiny(s32 fb_size)
/***********************************/
{
  tiny_qs_params *params = g_params;
  u16 *prime_list = params->prime_list;
  s32 i, j, mult_idx;
  tiny_fb *factor_base = params->factor_base;

  i = MIN_FB_OFFSET_TINY;
  mult_idx = 0;
  factor_base[i].prime = 2;
  params->multiplier_fb[0] = 0;
  params->multiplier_fb[1] = 0;

  /* Keep setting up factor base primes until enough
     are found or the pool of primes runs out */

  for (i++, j = 1; i < fb_size && j < NUM_PRIMES_TINY; j++) {
    tiny_fb *fbptr = factor_base + i;
    s32 prime = prime_list[j];
    s32 nmodp = mpz_tdiv_ui(params->n, prime);

    if (legendre_16(nmodp, prime) != -1) {
      fbptr->prime = prime;
      fbptr->logprime = (u8)(LOGPRIME_SCALE_TINY *
                             log((double)prime) / M_LN2 + 0.5);
      fbptr->recip = (u32)(B(32) / (u64)prime);

      /* if the prime divides n, it is part of n's 
         multiplier and is treated separately */

      if (nmodp != 0) {
        fbptr->modsqrt = (u16)sqrtModP_16(nmodp, prime);
      }
      else {
        fbptr->modsqrt = DO_NOT_SIEVE_TINY;
        params->multiplier_fb[mult_idx++] = i;
      }
      i++;
    }
  }
  params->fb_size = i;
  return i;
}

/***********************************/
static void fill_sieve_block_tiny(void)
/***********************************
Core sieving routine
************************************/
{
  tiny_qs_params *params = g_params;
  s32 i;
  s32 fb_size = params->fb_size;
  u8 *sieve_block = params->sieve_block;
  tiny_fb *factor_base = params->factor_base;

  /* Note that since this code will only ever 
     factor small inputs, the sieve interval will
     always be ridiculously small and does not 
     need to be broken up into chunks. Further,
     the bottleneck with small inputs is the trial
     factoring of relations and not the sieving,
     so no crazy unrolling tricks are needed
     here either */
     
  for (i = MIN_FB_OFFSET_TO_SIEVE_TINY; i < fb_size; i++) {
    tiny_fb *fbptr = factor_base + i;
    s32 prime = fbptr->prime;
    u8 logprime = fbptr->logprime;
    s32 root = fbptr->roots[0];

    while (root < SIEVE_SIZE_TINY) {
      sieve_block[root] -= logprime;
      root += prime;
    }

    root = fbptr->roots[1];
    while (root < SIEVE_SIZE_TINY) {
      sieve_block[root] -= logprime;
      root += prime;
    }
  }
}

#define PACKED_SIEVE_MASK ((u64)0x80808080 << 32 | 0x80808080)

/***********************************/
static s32 mark_sieve_block_tiny(void)
/***********************************
Walk through a filled-in sieve block and find 
the offsets correspodning to relations that
are probably useful
************************************/
{
  tiny_qs_params *params = g_params;
  s32 i, j, k;
  u8 *sieve_block = params->sieve_block;
  u64 *packed_sieve_block = (u64 *)params->sieve_block;

  /* standard technique for testing sieve locations
     in parallel: initialize each byte to the target
     sieve value, and subtract logs of the factor base
     primes instead of adding them. Sieve offsets that
     accumulate enough log values become negative, 
     and it's easy to simultaneously test for the top 
     bit in several bytes being set */

  for (i = j = 0; i < SIEVE_SIZE_TINY / 8; i += 4) {

    /* handle 32 bytes at a time */

    u64 accum = packed_sieve_block[i] |
                packed_sieve_block[i+1] |
                packed_sieve_block[i+2] |
                packed_sieve_block[i+3];

    if ((accum & PACKED_SIEVE_MASK) == (u64)(0))
      continue;

    /* at least one byte is a hit; go back and search
       the list one at a time. We treat the sieve interval
       as a hashtable, and associate entry j in the list
       of relations to be resieved (params->sieve_batch[])
       with a byte that is negative. The high-order bit of
       the byte is set to indicate that the low-order bits 
       mean something */

    for (k = 0; k < 32; k++) {
      u32 val = sieve_block[8 * i + k];
      if (val & 0x80) {
        if (j < SIEVE_BATCH_SIZE_TINY) {
          tiny_relation *r = params->sieve_batch + j;
          r->sieve_offset = 8 * i + k;
          r->num_factors = 0;
          sieve_block[8 * i + k] = j | 0x80;
          j++;
        }
        else {
          sieve_block[8 * i + k] = 0;
        }
      }
    }
  }

  return j;
}

/***********************************/
static void resieve_tiny(void)
/***********************************
Just like fill_sieve_block_tiny(), except
sieving is used to avoid trial division
on all the relations previously found
************************************/
{
  tiny_qs_params *params = g_params;
  s32 i;
  s32 fb_size = params->fb_size;
  u8 *sieve_block = params->sieve_block;
  tiny_fb *factor_base = params->factor_base;

  /* Note that even though this routine does only
     a little more work than fill_sieve_block_tiny(),
     it runs almost 3x slower */

  for (i = MIN_FB_OFFSET_TO_SIEVE_TINY; i < fb_size; i++) {
    tiny_fb *fbptr = factor_base + i;
    s32 prime = fbptr->prime;
    s32 root = fbptr->roots[0];

    while (root < SIEVE_SIZE_TINY) {
      s32 val = sieve_block[root];
      if (val & 0x80) {
        tiny_relation *r = params->sieve_batch + (val & 0x7f);
        r->fb_offsets[r->num_factors++] = i;
      }
      root += prime;
    }

    root = fbptr->roots[1];
    while (root < SIEVE_SIZE_TINY) {
      s32 val = sieve_block[root];
      if (val & 0x80) {
        tiny_relation *r = params->sieve_batch + (val & 0x7f);
        r->fb_offsets[r->num_factors++] = i;
      }
      root += prime;
    }
  }
}

/***********************************/
static s32 check_sieve_val_tiny(mpz_t a, mpz_t b, mpz_t c, 
                                 tiny_relation *r,
				 s32 sign_of_index)
/***********************************
Trial factor a relation that survived sieving
************************************/
{
  tiny_qs_params *params = g_params;
  s32 i, j;
  tiny_fb *factor_base = params->factor_base;
  s32 num_factors = 0;
  s32 sieve_offset = r->sieve_offset;
  tiny_relation *relation = params->relation_list +
                            params->num_full_relations;
  u16 *fb_offsets = relation->fb_offsets;
  static u8 initialized = 0;
  static mpz_t res, res2;

  if (initialized == 0) {
    mpz_init(res);
    mpz_init(res2);
    initialized = 1;
  }

  /* form the polynomial value */

  mpz_mul_ui(res, a, sieve_offset);
  if (sign_of_index == POSITIVE)
    mpz_add(res, res, b);
  else
    mpz_sub(res, res, b);
  mpz_mul_ui(res, res, sieve_offset);
  mpz_add(res, res, c);
  if (mpz_sgn(res) < 0) {
    mpz_abs(res, res);
    fb_offsets[num_factors++] = 0;
  }

  /* extract powers of two */

  i = mpz_scan1(res, 0);
  if (i) {
    mpz_tdiv_q_2exp(res, res, i);
    do {
      if (num_factors >= MAX_FACTORS_TINY)
        return 0;
      fb_offsets[num_factors++] = MIN_FB_OFFSET_TINY;
    } while (--i);
  }

  /* divide out the unsieved factor base primes */

  for (i = MIN_FB_OFFSET_TINY + 1; 
             i < MIN_FB_OFFSET_TO_SIEVE_TINY; i++) {
    tiny_fb *fbptr = factor_base + i;
    s32 prime = fbptr->prime;
    s32 root1 = fbptr->roots[0];
    s32 root2 = fbptr->roots[1];
    u32 recip = fbptr->recip;

    if (root1 == DO_NOT_SIEVE_TINY)
      continue;

    j = (s32)(((u64)sieve_offset * (u64)recip) >> 32);
    j = sieve_offset - j * prime;
    if (j >= prime)
      j -= prime;

    if (j == root1 || j == root2) {
      while (mpz_tdiv_q_ui(res2, res, prime) == 0) {
        if (num_factors >= MAX_FACTORS_TINY)
          return 0;

        fb_offsets[num_factors++] = i;
        mpz_swap(res, res2);
      }
    }
  }

  /* divide out the factors of the multiplier, 
     if any */
     
  for (i = 0; i < 2; i++) {
    if (params->multiplier_fb[i]) {
      s32 prime;
      j = params->multiplier_fb[i];
      prime = factor_base[j].prime;
      while (mpz_tdiv_q_ui(res2, res, prime) == 0) {
        if (num_factors >= MAX_FACTORS_TINY)
          return 0;
  
        fb_offsets[num_factors++] = j;
        mpz_swap(res, res2);
      }
    }
  }

  /* We should probably have been adding log values
     to the log of this relation in the previous loops,
     and testing that the complete log value now
     exceeds the trial factoring cutoff. However, 
     resieving has already found the remaining factors, 
     so we wouldn't save much time bailing out at 
     this point */

  for (i = 0; i < r->num_factors; i++) {
    s32 prime;
    j = r->fb_offsets[i];
    prime = factor_base[j].prime;

    while (mpz_tdiv_q_ui(res2, res, prime) == 0) {
      if (num_factors >= MAX_FACTORS_TINY)
        return 0;

      fb_offsets[num_factors++] = j;
      mpz_swap(res, res2);
    }
  }

  /* start filling in the final relation */

  if (sign_of_index == NEGATIVE)
    sieve_offset = -sieve_offset;
  relation->sieve_offset = sieve_offset;
  relation->num_factors = num_factors;
  relation->poly_num = params->poly_num;

  if (mpz_cmp_ui(res, 1) == 0) {

    /* full relation; we're done */

    relation->large_prime = 1;
    params->num_full_relations++;
  }
  else if (mpz_cmp_ui(res, params->large_prime_max) < 0) {
    u32 lp = mpz_get_ui(res);
    u32 table_idx = LARGE_PRIME_HASH(lp);
    s32 partial_idx;

    /* partial relation; see if it has occurred already */

    relation->large_prime = lp;
    for (i = 0; i < LP_HASH_DEPTH_TINY; i++) {
      partial_idx = params->partial_hash[table_idx][i];
      if (partial_idx == 0xffff ||
          lp == params->relation_list[partial_idx].large_prime)
        break;
    }

    if (i == LP_HASH_DEPTH_TINY) {

      /* not found, and no room to store it */

      return 0;
    }
    else if (partial_idx == 0xffff) {

      /* not found, but the hashtable entry has
         room to keep it; transfer the relation to
         the partial list */

      params->relation_list[params->partial_idx] = *relation;
      params->partial_hash[table_idx][i] = params->partial_idx--;
    }
    else {

      /* large prime has matched, new relation can stay */

      params->num_full_relations++;
    }
  }

  /* make sure the 'heap' of full relations has not
     overflowed into the 'stack' of partial relations */

  if (params->num_full_relations >= params->partial_idx)
    return -1;
  return 0;
}

/***********************************/
static void init_siqs_tiny(void)
/***********************************
Initialize the subystem for forming SIQS
sieve polynomials
************************************/
{
  tiny_qs_params *params = g_params;
  u32 i, j;
  u32 plus_idx, minus_idx;
  u32 fb_size = params->fb_size;
  u32 num_factors = params->num_a_factors;
  tiny_fb *factor_base = params->factor_base;

  /* compute the optimal size of the factors of
     the polynomial 'A' value. We know how many
     primes it should have, and know the optimal
     A value that will minimize sieving time. 
     Assume further that all factors are the
     same size. First compute the factor size,
     then locate the factor base offset where
     it approximately occurs */

  j = (u32)(exp(log(params->target_a) / num_factors) + 0.5);
  for (i = MIN_FB_OFFSET_TINY + 1; i < fb_size - 1; i++) {
    if (factor_base[i].prime > j)
      break;
  }
  if (i == MIN_FB_OFFSET_TINY + 1)
    i++;

  /* polynomial A values are built by selecting from
     a pool of primes. There are POLY_SELECT_BITS_TINY
     primes in the pool, evenly distributed above and
     below the optimal factor base offset */

  memset(params->poly_select_offsets, 0,
         sizeof(params->poly_select_offsets));
  plus_idx = i;
  minus_idx = i-1;
  i = 0;
  while (1) {
    if (plus_idx < fb_size && 
        factor_base[plus_idx].modsqrt != DO_NOT_SIEVE_TINY) {
      params->poly_select_offsets[i] = plus_idx;
      if (++i == POLY_SELECT_BITS_TINY)
        break;
    }
    if (minus_idx > MIN_FB_OFFSET_TINY + 1 &&
        factor_base[minus_idx].modsqrt != DO_NOT_SIEVE_TINY) {
      params->poly_select_offsets[i] = minus_idx;
      if (++i == POLY_SELECT_BITS_TINY)
        break;
    }
    plus_idx++;
    minus_idx--;
  }

  /* polynomial selection will begin at offset
     zero of the tables below */

  params->poly_select_idx = 0;
}

/* A perpetual problem with SIQS is deciding which
   primes should make up the next polynomial A value.
   The selected set must multiply out to a value as
   close as possible to the optimal A value, but must
   be sufficiently different from previously selected
   sets that the odds of producing duplicate relations
   are low. And the set has to be computed quickly.

   Fortunately, we will only need a few polynomials
   so the sets to use can be precomputed. Each prime in
   the pool is assigned a bit in a bitfield. Consecutive
   bits in the bitfield refer to primes alternately 
   above and below the optimal factor size. The low-order
   bits correspond to factors near the optimal value, and
   the more significant bits march away from the optimal
   value. Hence, setting the bitfield to an integer 
   will select a unique set of primes, the number of which
   is the number of set bits in the integer. Small integer
   values of the bitfield will pick primes close to the
   optimal factor size, with later bitfield values selecting
   prime factors that march away from the optimal size.

   popcount[] gives the number of set bits for each value
   of the bitfield, and a_choice[] lists the bitfields
   themselves. A given factorization only uses one of the
   population count sizes from the table; bitfields are
   arranged so that low-order bits are set first, then
   higher-order bits are set. Every bitfield is different
   by at least two bits from all other bitfields with the
   same weight, and there are enough bitfields to generate
   256 polynomials, whether A values contain 3, 4, or 5 primes */
   
static u8 popcount[] = {
       3,     4,     3,     5,     3,     4,
       3,     4,     3,     3,     4,     4,     
       3,     4,     5,     4,     5,     4,     
       4,     4,     4,     5,     5,     4,     
       4,     5,     5,     4,     5,     5,     
       5,     5,     3,     5,     5,     5,     
       5,     5,     3,     4,     3,     4,     
       4,     4,     3,     3,     4,     4,     
       4,     4,     3,     4,     4,     4,     
       4,     3,     4,     4,     3,     4,     
       4,     4,     4,     3,     4,     3,     
};

static u16 a_choice[] = {
       0x007, 0x00f, 0x019, 0x01f, 0x02a, 0x033,
       0x034, 0x03c, 0x04c, 0x052, 0x055, 0x05a,
       0x061, 0x066, 0x067, 0x069, 0x079, 0x096,
       0x099, 0x0a5, 0x0aa, 0x0ab, 0x0b5, 0x0c3,
       0x0cc, 0x0cd, 0x0d3, 0x0f0, 0x12d, 0x133,
       0x14b, 0x155, 0x181, 0x187, 0x199, 0x1e1,
       0x22e, 0x256, 0x282, 0x303, 0x304, 0x30c,
       0x330, 0x3c0, 0x484, 0x502, 0x505, 0x50a,
       0x550, 0x5a0, 0x601, 0x606, 0x609, 0x660,
       0x690, 0x888, 0x906, 0x909, 0x910, 0x960,
       0x990, 0xa05, 0xa0a, 0xa20, 0xa50, 0xc40,
};


/***********************************/
static s32 find_poly_a(mpz_t a)
/***********************************
Compute the next polynomial A value
************************************/
{
  tiny_qs_params *params = g_params;
  u32 i, j, mask;
  u32 num_a_factors = params->num_a_factors;
  tiny_fb *factor_base = params->factor_base;
  tiny_poly *poly = params->poly_list + params->poly_num;

  /* choose the next bitfield representing
     primes to use */

  for (i = params->poly_select_idx; i < sizeof(popcount); i++) {
    if (popcount[i] == num_a_factors)
      break;
  }
  if (i >= sizeof(popcount))
    return -1;
  mask = a_choice[i];
  params->poly_select_idx = i + 1;

  /* gather the chosen primes */

  for (i = j = 0; i < POLY_SELECT_BITS_TINY; i++) {
    if (!(mask & (1 << i)))
      continue;

    if (params->poly_select_offsets[i] == 0)
      return -2;

    poly->a_fb_offsets[j] = params->poly_select_offsets[i];
    if (++j == num_a_factors)
      break;
  }

  /* multiply them together */

  mpz_set_ui(a, 1);
  for (i = 0; i < num_a_factors; i++) {
    j = poly->a_fb_offsets[i];
    mpz_mul_ui(a, a, factor_base[j].prime);
  }

  return 0;
}

/***********************************/
static void find_first_poly_b(mpz_t a, mpz_t b, mpz_t c)
/***********************************
Compute the first of a list of polynomial
B values
************************************/
{
  tiny_qs_params *params = g_params;
  u32 i, j;
  u32 num_a_factors = params->num_a_factors;
  u32 fb_size = params->fb_size;
  tiny_fb *factor_base = params->factor_base;
  tiny_poly *poly = params->poly_list + params->poly_num;

  mpz_set_ui(b, 0);

  /* fill in the auxiliary quantities needed to
     compute future B values */

  for (i = 0; i < num_a_factors; i++) {
    tiny_fb *fbptr = factor_base + poly->a_fb_offsets[i];
    s32 g, prime = fbptr->prime;

    mpz_divexact_ui(params->poly_b_aux[i], a, prime);
    g = mpz_tdiv_ui(params->poly_b_aux[i], prime);
    g = modinv_16(g, prime);
    g = (s32)g * fbptr->modsqrt % prime;
    if (g > prime/2)
      g = prime - g;
    mpz_mul_ui(params->poly_b_aux[i], 
               params->poly_b_aux[i], g);
    mpz_add(b, b, params->poly_b_aux[i]);
    mpz_add(params->poly_b_aux[i],
            params->poly_b_aux[i],
            params->poly_b_aux[i]);
  }
  /* This first B is the sum of the auxiliary 
     quantities computed previously */

  mpz_set(poly->b, b);
  
  /* Form C, a helper for computing the value
     of a polynomial before trial factoring */

  mpz_mul(c, b, b);
  mpz_sub(c, c, params->n);
  mpz_divexact(c, c, a);

  /* initialize the factor base for sieving */

  for (i = MIN_FB_OFFSET_TINY + 1; i < fb_size; i++) {
    tiny_fb *fbptr = factor_base + i;
    s32 prime = fbptr->prime;
    s32 modsqrt = fbptr->modsqrt;
    s32 amodp = mpz_tdiv_ui(a, prime);
    s32 bmodp = prime - mpz_tdiv_ui(b, prime);

    if (fbptr->modsqrt == DO_NOT_SIEVE_TINY) {

      /* factors of the multiplier never 
         contribute to sieving */

      fbptr->roots[0] = DO_NOT_SIEVE_TINY;
      fbptr->roots[1] = DO_NOT_SIEVE_TINY;
      continue;
    }
    else if (amodp == 0) {

      /* factor base primes that divide the A value
         get one sieve root and not two */

      amodp = prime - mpz_tdiv_ui(c, prime);
      fbptr->roots[0] = amodp * modinv_16(2 * bmodp % prime, prime) % prime;
      fbptr->roots[1] = DO_NOT_SIEVE_TINY;
    }
    else {

      /* handle all the other FB primes, including the
         initialization that allows the next 2^(num_a_factors-1)-1
         factor bases to initialize quickly */

      amodp = modinv_16(amodp, prime);
      fbptr->roots[0] = amodp * (bmodp + modsqrt) % prime;
      fbptr->roots[1] = amodp * (bmodp + prime - modsqrt) % prime;

      for (j = 0; j < num_a_factors; j++) {
        bmodp = mpz_tdiv_ui(params->poly_b_aux[j], prime);
        params->root_aux[j * fb_size + i] = bmodp * amodp % prime;
      }
    }
  }
}

/***********************************/
static void find_next_poly_b(mpz_t a, mpz_t b, mpz_t c)
/***********************************
Initialize B values beyond the first
************************************/
{
  tiny_qs_params *params = g_params;
  s32 i, j;
  s32 num_a_factors = params->num_a_factors;
  s32 fb_size = params->fb_size;
  tiny_fb *factor_base = params->factor_base;
  tiny_poly *poly = params->poly_list + params->poly_num;
  u32 mask = params->poly_num & ((1 << (num_a_factors-1)) - 1);
  u8 do_sub;
  u16 *row;

  /* current poly starts of with the previous poly */

  mpz_set(b, poly[-1].b);
  for (i = 0; i < num_a_factors; i++)
    poly[0].a_fb_offsets[i] = poly[-1].a_fb_offsets[i];

  /* determine the auxiliary B constant that comes next
     in Gray code order, and add to or subtract from
     the current B. This also determines which of the
     rows from the table of corrections are applied to
     the factor base */

  i = 0;
  while ((mask & (1 << i)) == 0)
    i++;

  row = params->root_aux + fb_size * i;

  do_sub = 0;
  if (mask & (1 << (i+1))) {
    mpz_add(b, b, params->poly_b_aux[i]);
    do_sub = 1;
  }
  else {
    mpz_sub(b, b, params->poly_b_aux[i]);
  }

  /* form the C helper value */

  mpz_mul(c, b, b);
  mpz_sub(c, c, params->n);
  mpz_divexact(c, c, a);

  /* set up the factor base for the next B */

  for (j = MIN_FB_OFFSET_TINY + 1; j < fb_size; j++) {
    tiny_fb *fbptr = factor_base + j;
    s32 prime = fbptr->prime;
    s32 root1 = fbptr->roots[0];
    s32 root2 = fbptr->roots[1];

    /* apply the correction to each sieve root */

    if (root2 != DO_NOT_SIEVE_TINY) {

      /* ordinary FB prime. Note that the pevious sieving
         operation negated the roots to use, so they have
         to be negated again before the correction is applied */

      if (root1)
        root1 = prime - root1;
      if (root2)
        root2 = prime - root2;

      if (do_sub) {
        root1 -= row[j];
        root2 -= row[j];
      }
      else {
        root1 += row[j] - prime;
        root2 += row[j] - prime;
      }

      if (root1 < 0)
        root1 += prime;
      if (root2 < 0)
        root2 += prime;
      fbptr->roots[0] = root1;
      fbptr->roots[1] = root2;
    }
    else if (root1 != DO_NOT_SIEVE_TINY) {

      /* sieving with root1 but not root 2 only
         happens if the prime divides 'A'. Compute
         the new sieve root manually */

      s32 cmodp = prime - mpz_tdiv_ui(c, prime);
      s32 bmodp = mpz_tdiv_ui(b, prime);
      if (mpz_sgn(b) > 0)
        bmodp = prime - bmodp;
      fbptr->roots[0] = cmodp * modinv_16(2 * bmodp % prime, prime) % prime;
    }
  }
  mpz_set(poly[0].b, b);
}


/***********************************/
static s32 sieve_next_poly_tiny(void)
/***********************************
Do all the sieving for one polynomial
************************************/
{
  tiny_qs_params *params = g_params;
  s32 i;
  s32 fb_size = params->fb_size;
  u8 *sieve_block = params->sieve_block;
  tiny_fb *factor_base = params->factor_base;
  s32 cutoff1, num_surviving;
  s32 poly_num = params->poly_num;
  s32 target_relations = params->fb_size + NUM_EXTRA_RELATIONS_TINY;
  static u8 initialized = 0;
  static mpz_t a, b, c;

  if (initialized == 0) {
    mpz_init(a); mpz_init(b); mpz_init(c);
    initialized = 1;
  }

  /* generate the polynomial */

  if (!(poly_num & ((1 << (params->num_a_factors-1))-1))) {
    i = find_poly_a(a);
    if (i)
      return i;
    find_first_poly_b(a, b, c);
  }
  else {
    find_next_poly_b(a, b, c);
  }

  /* compute the cutoff beyond which trial factoring
     will be used on sieve values. */

  cutoff1 = LOGPRIME_SCALE_TINY * (mpz_sizeinbase(c, 2) - 
                  params->error_bits - SMALL_PRIME_FUDGE_TINY - 1);

  /* the trial factoring code wants 2*B and not B */

  mpz_add(b, b, b);
  
  /* sieve over positive offsets, mark the most
     promising offsets, resieve to trial factor
     them all at once and then finish each in turn */

  memset(sieve_block, cutoff1 - 1, SIEVE_SIZE_TINY);
  fill_sieve_block_tiny();
  num_surviving = mark_sieve_block_tiny();
  if (num_surviving) {
    resieve_tiny();
    for (i = 0; i < num_surviving; i++) {
      if (check_sieve_val_tiny(a, b, c,
                               params->sieve_batch + i,
                               POSITIVE) != 0) {
        return -3;
      }
      if (params->num_full_relations >= target_relations)
        return 0;
    }
  }
    
  /* flip the sieve roots from positive to
     negative values */
  for (i = MIN_FB_OFFSET_TINY + 1; i < fb_size; i++) {
    tiny_fb *fbptr = factor_base + i;
    s32 prime = fbptr->prime;
    s32 root1 = fbptr->roots[0];
    s32 root2 = fbptr->roots[1];
    if (root1 != DO_NOT_SIEVE_TINY && root1)
      fbptr->roots[0] = prime - root1;
    if (root2 != DO_NOT_SIEVE_TINY && root2)
      fbptr->roots[1] = prime - root2;
  }
    
  /* repeat the sieve procedure for negative
     sieve offsets */

  memset(sieve_block, cutoff1 - 1, SIEVE_SIZE_TINY);
  fill_sieve_block_tiny();
  num_surviving = mark_sieve_block_tiny();
  if (num_surviving) {
    resieve_tiny();
    for (i = 0; i < num_surviving; i++) {
      if (check_sieve_val_tiny(a, b, c,
                               params->sieve_batch + i,
                               NEGATIVE) != 0) {
        return -3;
      }
      if (params->num_full_relations >= target_relations)
        return 0;
    }
  }
    
  return 0;
}

/***********************************/
static void solve_linear_system_tiny(void)
/***********************************
Find linear dependencies among a set of relations
************************************/
{
  tiny_qs_params *params = g_params;
  s32 i, j, k, start_row;
  s32 nrows = params->fb_size;
  s32 ncols = params->num_full_relations;
  s32 num_a_factors = params->num_a_factors;
  u16 rowperm[MAX_FB_SIZE_TINY];
  u16 pivot[MAX_FB_SIZE_TINY];
  s32 row = 0;

  memset(params->matrix, 0, sizeof(params->matrix));

  /* build the matrix; relations become columns, and
     pairs of matched partial relations fuse into 
     columns as well */

  for (i = 0; i < ncols; i++) {
    tiny_relation *r;
    tiny_poly *poly;
    for (j = 0; j < 2; j++) {

      r = params->relation_list + i;
      if (j == 1) {
        s32 hash_idx = LARGE_PRIME_HASH(r->large_prime);
        s32 partial_idx;
        for (k = 0; k < LP_HASH_DEPTH_TINY; k++) {
          partial_idx = params->partial_hash[hash_idx][k];
          if (params->relation_list[partial_idx].large_prime == r->large_prime)
            break;
        }
        r = params->relation_list + partial_idx;
      }
      poly = params->poly_list + r->poly_num;

      for (k = 0; k < r->num_factors; k++) {
        row = r->fb_offsets[k];
        params->matrix[row][i / 64] ^= bitmask[i % 64];
      }

      /* the factors in the polynomial A value
         figure into the matrix as well */

      for (k = 0; k < num_a_factors; k++) {
        row = poly->a_fb_offsets[k];
        params->matrix[row][i / 64] ^= bitmask[i % 64];
      }
      if (r->large_prime == 1)
        break;
    }
  }
  for (i = 0; i < nrows; i++)
    rowperm[i] = i;

  /* begin with a random vector of dependencies */

  for (i = 0; i < ncols; i++)
    params->null_vectors[i] = (u16)get_rand(
                      &rand_seed1, &rand_seed2);

  /* perform the elimination */

  for (i = start_row = 0; start_row < nrows && i < ncols; i++) {
    
    /* find the next pivot */

    for (j = start_row; j < nrows; j++) {
      row = rowperm[j];
      if (params->matrix[row][i / 64] & bitmask[i % 64])
        break;
    }
    if (j == nrows)
      continue;

    rowperm[j] = rowperm[start_row];
    rowperm[start_row] = row;
    pivot[start_row++] = i;

    /* eliminate it from the other rows */

    for (j++; j < nrows; j++) {
      s32 row2 = rowperm[j];
      if (params->matrix[row2][i / 64] & bitmask[i % 64]) {
        for (k = i / 64; k < (ncols + 63) / 64; k++) {
          params->matrix[row2][k] ^= params->matrix[row][k];
        }
      }
    }
  }

  /* perform back substitution */

  for (i = start_row - 1; i >= 0; i--) {
    u16 accum;
    row = rowperm[i];

    for (j = pivot[i] + 1, accum = 0; j < ncols; j++) {
      if (params->matrix[row][j / 64] & bitmask[j & 63])
        accum ^= params->null_vectors[j];
    }
    params->null_vectors[pivot[i]] = accum;
  }
}

/***********************************/
static u32 find_factors_tiny(mpz_t factor1, 
                             mpz_t factor2)
/***********************************
perform MPQS square root phase
************************************/
{
  tiny_qs_params *params = g_params;
  s32 i, j, k;
  u16 mask;
  u16 fb_counts[MAX_FB_SIZE_TINY];
  tiny_fb *factor_base = params->factor_base;
  static mpz_t x, y, t0, t1;
  static u8 initialized = 0;

  if (initialized == 0) {
    mpz_init(x); mpz_init(y);
    mpz_init(t0); mpz_init(t1);
    initialized = 1;
  }

  /* for each dependency */

  for (mask = 1; mask; mask <<= 1) {

    memset(fb_counts, 0, sizeof(fb_counts));
    mpz_set_ui(x, 1);
    mpz_set_ui(y, 1);

    /* for each relation allowed in the dependency */

    for (i = 0; i < params->num_full_relations; i++) {

      if (!(params->null_vectors[i] & mask))
        continue;

      for (j = 0; j < 2; j++) {
        tiny_relation *r = params->relation_list + i;
        tiny_poly *poly;

        /* match up partials with the same large prime */

        if (j == 1) {
          s32 hash_idx = LARGE_PRIME_HASH(r->large_prime);
          s32 partial_idx;
          for (k = 0; k < LP_HASH_DEPTH_TINY; k++) {
            partial_idx = params->partial_hash[hash_idx][k];
            if (params->relation_list[partial_idx].large_prime == 
                r->large_prime)
              break;
          }
          r = params->relation_list + partial_idx;
          mpz_mul_ui(t0, y, r->large_prime);
          mpz_mod(y, t0, params->n);
        }
        poly = params->poly_list + r->poly_num;
  
        /* add the factors of this relation to the table
           of factors. Include the factors of A as well */

        for (k = 0; k < r->num_factors; k++)
          fb_counts[r->fb_offsets[k]]++;
  
        mpz_set_ui(t1, 1);
        for (k = 0; k < params->num_a_factors; k++) {
          s32 idx = poly->a_fb_offsets[k];
          fb_counts[idx]++;
          mpz_mul_ui(t1, t1, factor_base[idx].prime);
        }

        /* multiply A * sieve_offset + B into the left 
           side of the congruence */

        if (r->sieve_offset < 0) {
          mpz_mul_ui(t1, t1, -(r->sieve_offset));
          mpz_sub(t1, t1, poly->b);
        }
        else {
          mpz_mul_ui(t1, t1, r->sieve_offset);
          mpz_add(t1, t1, poly->b);
        }
        mpz_mul(t0, x, t1);
        mpz_mod(x, t0, params->n);

        if (r->large_prime == 1)
          break;
      }
    }

    /* Form the right side of the congruence; given its
       prime factorization, cut the exponent of each prime
       in half and perform a modular exponentiation */

    for (i = MIN_FB_OFFSET_TINY; i < params->fb_size; i++) {
      u16 mask2 = 0x8000;
      u16 exponent = fb_counts[i] / 2;
      u32 prime = params->factor_base[i].prime;

      if (exponent == 0)
        continue;

      mpz_set_ui(t0, prime);
      while (!(exponent & mask2))
        mask2 >>= 1;

      for (mask2 >>= 1; mask2; mask2 >>= 1) {
        mpz_mul(t1, t0, t0);
        mpz_mod(t0, t1, params->n);
        if (exponent & mask2) {
          mpz_mul_ui(t1, t0, prime);
          mpz_mod(t0, t1, params->n);
        }
      }
      mpz_mul(t1, t0, y);
      mpz_mod(y, t1, params->n);
    }

    /* For x and y the halves of the congruence, 
       compute gcd(x+-y, n) */

    for (i = 0; i < 2; i++) {
      if (i == 0)
        mpz_add(t0, x, y);
      else
        mpz_sub(t0, x, y);

      mpz_gcd(t1, t0, params->n);
      if (mpz_cmp_ui(t1, 1) && mpz_cmp(t1, params->n)) {

        /* we've possibly found a nontrivial factor of n.
           Divide any factors of the multiplier out from
           both factors */

        u32 mult1 = 0;
        u32 mult2 = 0;

        if (params->multiplier_fb[0])
          mult1 = params->factor_base[params->multiplier_fb[0]].prime;
        if (params->multiplier_fb[1])
          mult2 = params->factor_base[params->multiplier_fb[1]].prime;

        mpz_divexact(t0, params->n, t1);
        if (mult1) {
          if (mpz_tdiv_ui(t0, mult1) == 0)
            mpz_divexact_ui(t0, t0, mult1);
          if (mpz_tdiv_ui(t1, mult1) == 0)
            mpz_divexact_ui(t1, t1, mult1);
        }
        if (mult2) {
          if (mpz_tdiv_ui(t0, mult2) == 0)
            mpz_divexact_ui(t0, t0, mult2);
          if (mpz_tdiv_ui(t1, mult2) == 0)
            mpz_divexact_ui(t1, t1, mult2);
        }

        /* If both remaining factors exceed unity, 
           we've factored n and can stop */
        if (mpz_cmp_ui(t0, 1) && mpz_cmp_ui(t1, 1)) {
          mpz_set(factor1, t0);
          mpz_set(factor2, t1);
          return 1;
        }
      }
    }

    /* otherwise try the next dependency */
  }

  return 0;
}

typedef struct {
  s32 fb_size;
  s32 num_poly_factors;
} tiny_qs_config;

/* factor base sizes for 50 to 120-bit 
   factorizations */

static tiny_qs_config static_config[] = {
 { 40, 3 },
 { 50, 3 },
 { 60, 3 },
 { 70, 3 },
 { 80, 3 },
 { 90, 3 },
 { 110, 3 },
 { 120, 3 },
 { 140, 3 },
 { 140, 3 },
 { 160, 3 },
 { 180, 4 },
 { 230, 4 },
 { 280, 4 },
 { 350, 4 },
 { 420, 4 },
 { 490, 5 },
};

/***********************************/
unsigned int tinyqs(mpz_t n, mpz_t factor)
/***********************************
Main driver for MPQS factorization
Returns 1 and sets factor if
successful, returns 0 otherwise
************************************/
{
  tiny_qs_params *params;
  s32 bits, status = 0;
  s32 fb_size;
  s32 bound;
  s32 large_prime_mult;
  tiny_qs_config *config;
  mpz_t tmp;

  mpz_init(tmp);

  /* make sure the input isn't a perfect square.
     We may also want to add a test for a perfect
     cube, but that's so unlikely it's probably
     not worth worrying about */

  if (mpz_root(tmp, n, 2) != 0) {
    mpz_set(factor, tmp);
    mpz_clear(tmp);
    return 1;
  }

  /* start the initialization */

  init_tinyqs();
  params = g_params;
  mpz_set(params->n, n);
  params->num_full_relations = 0;
  params->partial_idx = 4 * MAX_RELATIONS_TINY - 1;
  params->poly_num = 0;

  bits = mpz_sizeinbase(params->n, 2);
  find_multiplier_tiny();

  /* determine the factor base size and the
     number of primes in a polynomial A value */

  if (bits < 50)
    bits = 50;
  if (bits > 116)
    bits = 116;
  config = static_config + ((bits - 50) / 4);
  fb_size = config->fb_size;
  params->num_a_factors = config->num_poly_factors;

  /* build the factor base */

  fb_size = init_fb_tiny(fb_size);

  /* compute the optimal A value */

  mpz_sqrt(tmp, params->n);
  params->target_a = mpz_get_d(tmp) * M_SQRT2 / SIEVE_SIZE_TINY;
  init_siqs_tiny();

  /* compute the large prime cutoff and the
     size of the fudge factor needed to account
     for it in the sieving cutoff */

  large_prime_mult = 15;
  bound = params->factor_base[fb_size - 1].prime;
  bound *= large_prime_mult;
  params->large_prime_max = bound;
  params->error_bits = (u32)(log(bound) / M_LN2 + 1);

  /* empty out the hashtable for partial relations */

  memset(params->partial_hash, 0xff, sizeof(params->partial_hash));

  /* do the sieving! */

  while (params->poly_num < MAX_POLY_TINY &&
         params->num_full_relations < fb_size + NUM_EXTRA_RELATIONS_TINY) {
    if (sieve_next_poly_tiny() != 0) {
      mpz_clear(tmp);
      return 0;
    }
    params->poly_num++;
  }

  /* if enough relations were found, finish
     off the factorization */

  if (params->num_full_relations >= fb_size + NUM_EXTRA_RELATIONS_TINY) {
    solve_linear_system_tiny();
    status = find_factors_tiny(factor, tmp);
    /* Return smaller of two factors */
    if (status && mpz_cmp(factor,tmp) > 0)  mpz_swap(factor,tmp);
  }

  mpz_clear(tmp);
  return status;
}
