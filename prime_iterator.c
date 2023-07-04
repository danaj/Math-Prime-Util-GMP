#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ptypes.h"
#define FUNC_isqrt 1
#include "utility.h"
#include "prime_iterator.h"


/* Add this to a number and you'll ensure you're on a wheel location */
static const unsigned char distancewheel30[30] =
    {1,0,5,4,3,2,1,0,3,2,1,0,1,0,3,2,1,0,1,0,3,2,1,0,5,4,3,2,1,0};
/* The bit mask within a byte */
static const unsigned char masktab30[30] = {
    0,  1,  0,  0,  0,  0,  0,  2,  0,  0,  0,  4,  0,  8,  0,
    0,  0, 16,  0, 32,  0,  0,  0, 64,  0,  0,  0,  0,  0,128  };
static const unsigned char nextwheel30[30] = {
    1,  7,  7,  7,  7,  7,  7, 11, 11, 11, 11, 13, 13, 17, 17,
   17, 17, 19, 19, 23, 23, 23, 23, 29, 29, 29, 29, 29, 29,  1 };
#if 0
static const unsigned char prevwheel30[30] = {
   29, 29,  1,  1,  1,  1,  1,  1,  7,  7,  7,  7, 11, 11, 13,
   13, 13, 13, 17, 17, 19, 19, 19, 19, 23, 23, 23, 23, 23, 23 };
#endif

static INLINE UV next_prime_in_segment( const unsigned char* sieve, UV segment_start, UV segment_bytes, UV p)
{
  UV d, m;
  if (p < segment_start) return 0;
  d = (p-segment_start)/30;
  if (d >= segment_bytes) return 0;
  m = (p-segment_start) - d*30;
  do {
    if (m==29) {
      d++; m = 1;
      if (d >= segment_bytes) return 0;
    } else {
      m = nextwheel30[m];
    }
  } while (sieve[d] & masktab30[m]);
  return (segment_start + d*30 + m);
}
static INLINE int is_prime_in_segment( const unsigned char* sieve, UV segment_start, UV segment_bytes, UV p)
{
  UV d, m, mtab;
  if (p < segment_start) return -1;
  d = (p-segment_start)/30;
  if (d >= segment_bytes) return -1;
  m = (p-segment_start) - d*30;
  mtab = masktab30[m];
  if (mtab == 0)  return 0;
  return ((sieve[d] & mtab) == 0);
}

#define next_prime_in_sieve(sieve, p) next_prime_in_segment(sieve, 0, UV_MAX, p)


/* 1001 bytes of presieved mod-30 bytes.  If the area to be sieved is
 * appropriately filled with this data, then 7, 11, and 13 do not have
 * to be sieved.  It wraps, so multiple memcpy's can be used.  Do be
 * aware that if you start at 0, you'll have to correct the first byte.
 */
#define PRESIEVE_SIZE (7*11*13)
static const unsigned char presieve13[PRESIEVE_SIZE] =
{ 0x0e,0x20,0x10,0x81,0x49,0x24,0xc2,0x06,0x2a,0x90,0xa1,0x0c,0x14,
  0x58,0x02,0x61,0x11,0xc3,0x28,0x0c,0x44,0x22,0xa4,0x10,0x91,0x18,
  0x4d,0x40,0x82,0x21,0x58,0xa1,0x28,0x04,0x42,0x92,0x20,0x51,0x91,
  0x8a,0x04,0x48,0x03,0x60,0x34,0x81,0x1c,0x06,0xc1,0x02,0xa2,0x10,
  0x89,0x08,0x24,0x45,0x42,0x30,0x10,0xc5,0x0a,0x86,0x40,0x0a,0x30,
  0x38,0x85,0x08,0x15,0x40,0x63,0x20,0x96,0x83,0x88,0x04,0x60,0x16,
  0x28,0x10,0x81,0x49,0x44,0xe2,0x02,0x2c,0x12,0xa1,0x0c,0x04,0x50,
  0x0a,0x61,0x10,0x83,0x48,0x2c,0x40,0x26,0x26,0x90,0x91,0x08,0x55,
  0x48,0x82,0x20,0x19,0xc1,0x28,0x04,0x44,0x12,0xa0,0x51,0x81,0x9a,
  0x0c,0x48,0x02,0x21,0x54,0xa1,0x18,0x04,0x43,0x82,0xa2,0x10,0x99,
  0x08,0x24,0x44,0x03,0x70,0x30,0xc1,0x0c,0x86,0xc0,0x0a,0x20,0x30,
  0x8d,0x08,0x14,0x41,0x43,0x20,0x92,0x85,0x0a,0x84,0x60,0x06,0x30,
  0x18,0x81,0x49,0x05,0xc2,0x22,0x28,0x14,0xa3,0x8c,0x04,0x50,0x12,
  0x69,0x10,0x83,0x09,0x4c,0x60,0x22,0x24,0x12,0x91,0x08,0x45,0x50,
  0x8a,0x20,0x18,0x81,0x68,0x24,0x40,0x16,0x22,0xd1,0x81,0x8a,0x14,
  0x48,0x02,0x20,0x15,0xc1,0x38,0x04,0x45,0x02,0xa2,0x10,0x89,0x18,
  0x2c,0x44,0x02,0x31,0x50,0xe1,0x08,0x86,0x42,0x8a,0x20,0x30,0x95,
  0x08,0x14,0x40,0x43,0x60,0xb2,0x81,0x0c,0x06,0xe0,0x06,0x20,0x10,
  0x89,0x49,0x04,0xc3,0x42,0x28,0x10,0xa5,0x0e,0x84,0x50,0x02,0x71,
  0x18,0x83,0x08,0x0d,0x40,0x22,0x24,0x14,0x93,0x88,0x45,0x40,0x92,
  0x28,0x18,0x81,0x29,0x44,0x60,0x12,0x24,0x53,0x81,0x8a,0x04,0x58,
  0x0a,0x20,0x14,0x81,0x58,0x24,0x41,0x06,0xa2,0x90,0x89,0x08,0x34,
  0x4c,0x02,0x30,0x11,0xc1,0x28,0x86,0x44,0x0a,0xa0,0x30,0x85,0x18,
  0x1c,0x40,0x43,0x21,0xd2,0xa1,0x08,0x04,0x62,0x86,0x20,0x10,0x91,
  0x49,0x04,0xc2,0x03,0x68,0x30,0xa1,0x0c,0x06,0xd0,0x02,0x61,0x10,
  0x8b,0x08,0x0c,0x41,0x62,0x24,0x10,0x95,0x0a,0xc5,0x40,0x82,0x30,
  0x18,0x81,0x28,0x05,0x40,0x32,0x20,0x55,0x83,0x8a,0x04,0x48,0x12,
  0x28,0x14,0x81,0x19,0x44,0x61,0x02,0xa6,0x12,0x89,0x08,0x24,0x54,
  0x0a,0x30,0x10,0xc1,0x48,0xa6,0x40,0x0e,0x22,0xb0,0x85,0x08,0x14,
  0x48,0x43,0x20,0x93,0xc1,0x28,0x04,0x64,0x06,0xa0,0x10,0x81,0x59,
  0x0c,0xc2,0x02,0x29,0x50,0xa1,0x0c,0x04,0x52,0x82,0x61,0x10,0x93,
  0x08,0x0c,0x40,0x23,0x64,0x30,0x91,0x0c,0x47,0xc0,0x82,0x20,0x18,
  0x89,0x28,0x04,0x41,0x52,0x20,0x51,0x85,0x8a,0x84,0x48,0x02,0x30,
  0x1c,0x81,0x18,0x05,0x41,0x22,0xa2,0x14,0x8b,0x88,0x24,0x44,0x12,
  0x38,0x10,0xc1,0x09,0xc6,0x60,0x0a,0x24,0x32,0x85,0x08,0x14,0x50,
  0x4b,0x20,0x92,0x81,0x48,0x24,0x60,0x06,0x22,0x90,0x81,0x49,0x14,
  0xca,0x02,0x28,0x11,0xe1,0x2c,0x04,0x54,0x02,0xe1,0x10,0x83,0x18,
  0x0c,0x40,0x22,0x25,0x50,0xb1,0x08,0x45,0x42,0x82,0x20,0x18,0x91,
  0x28,0x04,0x40,0x13,0x60,0x71,0x81,0x8e,0x06,0xc8,0x02,0x20,0x14,
  0x89,0x18,0x04,0x41,0x42,0xa2,0x10,0x8d,0x0a,0xa4,0x44,0x02,0x30,
  0x18,0xc1,0x08,0x87,0x40,0x2a,0x20,0x34,0x87,0x88,0x14,0x40,0x53,
  0x28,0x92,0x81,0x09,0x44,0x60,0x06,0x24,0x12,0x81,0x49,0x04,0xd2,
  0x0a,0x28,0x10,0xa1,0x4c,0x24,0x50,0x06,0x63,0x90,0x83,0x08,0x1c,
  0x48,0x22,0x24,0x11,0xd1,0x28,0x45,0x44,0x82,0xa0,0x18,0x81,0x38,
  0x0c,0x40,0x12,0x21,0x51,0xa1,0x8a,0x04,0x4a,0x82,0x20,0x14,0x91,
  0x18,0x04,0x41,0x03,0xe2,0x30,0x89,0x0c,0x26,0xc4,0x02,0x30,0x10,
  0xc9,0x08,0x86,0x41,0x4a,0x20,0x30,0x85,0x0a,0x94,0x40,0x43,0x30,
  0x9a,0x81,0x08,0x05,0x60,0x26,0x20,0x14,0x83,0xc9,0x04,0xc2,0x12,
  0x28,0x10,0xa1,0x0d,0x44,0x70,0x02,0x65,0x12,0x83,0x08,0x0c,0x50,
  0x2a,0x24,0x10,0x91,0x48,0x65,0x40,0x86,0x22,0x98,0x81,0x28,0x14,
  0x48,0x12,0x20,0x51,0xc1,0xaa,0x04,0x4c,0x02,0xa0,0x14,0x81,0x18,
  0x0c,0x41,0x02,0xa3,0x50,0xa9,0x08,0x24,0x46,0x82,0x30,0x10,0xd1,
  0x08,0x86,0x40,0x0b,0x60,0x30,0x85,0x0c,0x16,0xc0,0x43,0x20,0x92,
  0x89,0x08,0x04,0x61,0x46,0x20,0x10,0x85,0x4b,0x84,0xc2,0x02,0x38,
  0x18,0xa1,0x0c,0x05,0x50,0x22,0x61,0x14,0x83,0x88,0x0c,0x40,0x32,
  0x2c,0x10,0x91,0x09,0x45,0x60,0x82,0x24,0x1a,0x81,0x28,0x04,0x50,
  0x1a,0x20,0x51,0x81,0xca,0x24,0x48,0x06,0x22,0x94,0x81,0x18,0x14,
  0x49,0x02,0xa2,0x11,0xc9,0x28,0x24,0x44,0x02,0xb0,0x10,0xc1,0x18,
  0x8e,0x40,0x0a,0x21,0x70,0xa5,0x08,0x14,0x42,0xc3,0x20,0x92,0x91,
  0x08,0x04,0x60,0x07,0x60,0x30,0x81,0x4d,0x06,0xc2,0x02,0x28,0x10,
  0xa9,0x0c,0x04,0x51,0x42,0x61,0x10,0x87,0x0a,0x8c,0x40,0x22,0x34,
  0x18,0x91,0x08,0x45,0x40,0xa2,0x20,0x1c,0x83,0xa8,0x04,0x40,0x12,
  0x28,0x51,0x81,0x8b,0x44,0x68,0x02,0x24,0x16,0x81,0x18,0x04,0x51,
  0x0a,0xa2,0x10,0x89,0x48,0x24,0x44,0x06,0x32,0x90,0xc1,0x08,0x96,
  0x48,0x0a,0x20,0x31,0xc5,0x28,0x14,0x44,0x43,0xa0,0x92,0x81,0x18,
  0x0c,0x60,0x06,0x21,0x50,0xa1,0x49,0x04,0xc2,0x82,0x28,0x10,0xb1,
  0x0c,0x04,0x50,0x03,0x61,0x30,0x83,0x0c,0x0e,0xc0,0x22,0x24,0x10,
  0x99,0x08,0x45,0x41,0xc2,0x20,0x18,0x85,0x2a,0x84,0x40,0x12,0x30,
  0x59,0x81,0x8a,0x05,0x48,0x22,0x20,0x14,0x83,0x98,0x04,0x41,0x12,
  0xaa,0x10,0x89,0x09,0x64,0x64,0x02,0x34,0x12,0xc1,0x08,0x86,0x50,
  0x0a,0x20,0x30,0x85,0x48,0x34,0x40,0x47,0x22,0x92,0x81,0x08,0x14,
  0x68,0x06,0x20,0x11,0xc1,0x69,0x04,0xc6,0x02,0xa8,0x10,0xa1,0x1c,
  0x0c,0x50,0x02,0x61,0x50,0xa3,0x08,0x0c,0x42,0xa2,0x24,0x10,0x91,
  0x08,0x45,0x40,0x83,0x60,0x38,0x81,0x2c,0x06,0xc0,0x12,0x20,0x51,
  0x89,0x8a,0x04,0x49,0x42,0x20,0x14,0x85,0x1a,0x84,0x41,0x02,0xb2,
  0x18,0x89,0x08,0x25,0x44,0x22,0x30,0x14,0xc3,0x88,0x86,0x40,0x1a,
  0x28,0x30,0x85,0x09,0x54,0x60,0x43,0x24,0x92,0x81,0x08,0x04,0x70};

static void sieve_prefill(unsigned char* mem, UV startd, UV endd)
{
  UV nbytes = endd - startd + 1;
  MPUassert( (mem != 0) && (endd >= startd), "sieve_prefill bad arguments");

  /* Walk the memory, tiling in the presieve area using memcpy.
   * This is pretty fast, but it might still benefit from using copy
   * doubling (where we copy to the memory, then copy memory to memory
   * doubling in size each time), as memcpy usually loves big chunks.
   */
  while (startd <= endd) {
    UV pstartd = startd % PRESIEVE_SIZE;
    UV sieve_bytes = PRESIEVE_SIZE - pstartd;
    UV bytes = (nbytes > sieve_bytes) ? sieve_bytes : nbytes;
    memcpy(mem, presieve13 + pstartd, bytes);
    if (startd == 0)  mem[0] = 0x01; /* Correct first byte */
    startd += bytes;
    mem += bytes;
    nbytes -= bytes;
  }
}

/* Marking primes is done the same way we used to do with tables, but
 * now uses heavily unrolled code based on Kim Walisch's mod-30 sieve.
 */
#define set_bit(s,n)  *(s) |= (1 << n);
static const unsigned char masknum30[30] =
    {0,0,0,0,0,0,0,1,0,0,0,2,0,3,0,0,0,4,0,5,0,0,0,6,0,0,0,0,0,7};
static const unsigned char qinit30[30] =
    {0,0,1,1,1,1,1,1,2,2,2,2,3,3,4,4,4,4,5,5,6,6,6,6,7,7,7,7,7,7};
static const UV max_sieve_prime = (BITS_PER_WORD==64) ? 4294967291U : 65521U;

#define CROSS_INDEX(v, b0,b1,b2,b3,b4,b5,b6,b7,  i0,i1,i2,i3,i4,i5,i6,i7, it) \
  while (1) { \
    case (v+0): if (s >= send) break;  set_bit(s,b0);  s += r*6+i0; \
    case (v+1): if (s >= send) break;  set_bit(s,b1);  s += r*4+i1; \
    case (v+2): if (s >= send) break;  set_bit(s,b2);  s += r*2+i2; \
    case (v+3): if (s >= send) break;  set_bit(s,b3);  s += r*4+i3; \
    case (v+4): if (s >= send) break;  set_bit(s,b4);  s += r*2+i4; \
    case (v+5): if (s >= send) break;  set_bit(s,b5);  s += r*4+i5; \
    case (v+6): if (s >= send) break;  set_bit(s,b6);  s += r*6+i6; \
    case (v+7): if (s >= send) break;  set_bit(s,b7);  s += r*2+i7; \
    while (s + r*28 + it-1 < send) { \
      set_bit(s + r *  0 +  0, b0); \
      set_bit(s + r *  6 + i0, b1); \
      set_bit(s + r * 10 + i0+i1, b2); \
      set_bit(s + r * 12 + i0+i1+i2, b3); \
      set_bit(s + r * 16 + i0+i1+i2+i3, b4); \
      set_bit(s + r * 18 + i0+i1+i2+i3+i4, b5); \
      set_bit(s + r * 22 + i0+i1+i2+i3+i4+i5, b6); \
      set_bit(s + r * 28 + i0+i1+i2+i3+i4+i5+i6, b7); \
      s += r*30 + it; \
    } \
  }
static void mark_primes(unsigned char* s, const unsigned char* send, UV startp, UV endp, UV prime)
{
  UV p2, q, r;
  int index;

  q = prime;
  p2 = prime * prime;
  if (p2 < startp) {
    q = 1+(startp-1)/prime;
    q += distancewheel30[q % 30];
    p2 = prime * q;
  }
  if (p2 > endp || p2 < startp) return;

  s += (p2-startp) / 30;
  r = prime / 30;
  index = qinit30[q % 30]  +  8*masknum30[prime % 30];

  switch (index) {
    CROSS_INDEX( 0, 0,1,2,3,4,5,6,7, 0,0,0,0,0,0,0,1,  1); break;
    CROSS_INDEX( 8, 1,5,4,0,7,3,2,6, 1,1,1,0,1,1,1,1,  7); break;
    CROSS_INDEX(16, 2,4,0,6,1,7,3,5, 2,2,0,2,0,2,2,1, 11); break;
    CROSS_INDEX(24, 3,0,6,5,2,1,7,4, 3,1,1,2,1,1,3,1, 13); break;
    CROSS_INDEX(32, 4,7,1,2,5,6,0,3, 3,3,1,2,1,3,3,1, 17); break;
    CROSS_INDEX(40, 5,3,7,1,6,0,4,2, 4,2,2,2,2,2,4,1, 19); break;
    CROSS_INDEX(48, 6,2,3,7,0,4,5,1, 5,3,1,4,1,3,5,1, 23); break;
    CROSS_INDEX(56, 7,6,5,4,3,2,1,0, 6,4,2,4,2,4,6,1, 29); break;
  }
}

/* Monolithic mod-30 wheel sieve */
static unsigned char* sieve_erat30(UV end)
{
  unsigned char *mem;
  UV max_buf, limit, prime;

  max_buf = (end/30) + ((end%30) != 0);
  /* Round up to a word */
  max_buf = ((max_buf + sizeof(UV) - 1) / sizeof(UV)) * sizeof(UV);
  New(0, mem, max_buf, unsigned char);

  /* Fill buffer with marked 7, 11, and 13 */
  sieve_prefill(mem, 0, max_buf-1);

  limit = isqrt(end);  /* prime*prime can overflow */
  for (prime = 17; prime <= limit; prime = next_prime_in_sieve(mem,prime)) {
    mark_primes(mem, mem+max_buf, 0, end, prime);
  }
  return mem;
}

/* Segmented mod-30 wheel sieve */
static int sieve_segment(unsigned char* mem, UV startd, UV endd,
                         const unsigned char* prim_sieve, UV prim_limit)
{
  const unsigned char* sieve;
  UV limit, p;
  UV startp = 30*startd;
  UV endp = (endd >= (UV_MAX/30))  ?  UV_MAX-2  :  30*endd+29;

  MPUassert( (mem != 0) && (endd >= startd) && (endp >= startp),
             "sieve_segment bad arguments");

  /* Fill buffer with marked 7, 11, and 13 */
  sieve_prefill(mem, startd, endd);

  limit = isqrt(endp);
  if (limit > max_sieve_prime)  limit = max_sieve_prime;

  sieve = (prim_sieve != 0 && limit <= prim_limit)
        ? prim_sieve
        : sieve_erat30(limit);
  MPUassert( sieve != 0, "Could not generate base sieve" );

  for (p = 17; p <= limit; p = next_prime_in_sieve(sieve,p)) {
    mark_primes(mem, mem+endd-startd+1, startp, endp, p);
  }

  if (sieve != prim_sieve)  Safefree(sieve);
  return 1;
}




/*****************************************************************************/
/*                            Prime iterator                                 */
/*****************************************************************************/

/* These sizes are a tradeoff.  For better memory use I think 16k,4k is good.
 * For performance, 32k,16k or 64k,16k is better.  To avoid threading hell,
 * this is just decided statically.  At 24k,16k we handle 736800 numbers in
 * the primary sieve and won't redo for segments until after 5*10^11.  Each
 * segment will store a range of 30*(16384-16) = 491040 numbers.
 */
#define PRIMARY_SIZE  (32768-16)
#define SEGMENT_SIZE  (24576-16)
#define NSMALL_PRIMES (83970-180)

static const unsigned char* primary_sieve = 0;
static const UV primary_limit = (30 * PRIMARY_SIZE)-1;
static const uint32_t* small_primes = 0;
static UV num_small_primes = 0;

void prime_iterator_global_startup(void)
{
  primary_sieve = sieve_erat30(primary_limit);
#ifdef NSMALL_PRIMES
  {
    UV p;
    uint32_t *primes32;
    UV *primes64 = sieve_to_n(NSMALL_PRIMES + 180, &num_small_primes);
    New(0, primes32, num_small_primes, uint32_t);
    for (p = 0; p < num_small_primes; p++)  primes32[p] = primes64[p];
    Safefree(primes64);
    small_primes = primes32;
  }
#endif
}

void prime_iterator_global_shutdown(void)
{
  if (primary_sieve != 0)  Safefree(primary_sieve);
  if (small_primes != 0)   Safefree(small_primes);
  primary_sieve = 0;
  small_primes = 0;
}

void prime_iterator_init(prime_iterator *iter)
{
  iter->p = 2;
  iter->segment_start = 0;
  iter->segment_bytes = 0;
  iter->segment_mem = 0;
}

#if 0
prime_iterator prime_iterator_default(void)
{
  prime_iterator iter = {2, 0, 0, 0};
  return iter;
}
#endif

void prime_iterator_destroy(prime_iterator *iter)
{
  if (iter->segment_mem != 0)  Safefree(iter->segment_mem);
  iter->segment_mem = 0;
  iter->segment_start = 0;
  iter->segment_bytes = 0;
  iter->p = 0;
}

#ifdef NSMALL_PRIMES
static UV pcount(UV n)
{
  UV lo = 0 + (n >> 4);
  UV hi = (n >> 3) - (n >> 6) + ( (n<503) ? 40 : (n<1669) ? 80 : 139 );
  if (hi > num_small_primes) hi = num_small_primes;
  while (lo < hi) {
    UV mid = lo + (hi-lo)/2;
    if (small_primes[mid] <= n) lo = mid+1;
    else                        hi = mid;
  }
  return lo;   /* Because 2 is stored at location 0 */
}
#endif

void prime_iterator_setprime(prime_iterator *iter, UV n) {
  /* Is it inside the current segment? */
  if (    (iter->segment_mem != 0)
       && (n >= iter->segment_start)
       && (n <= iter->segment_start + 30*iter->segment_bytes - 1) ) {
    iter->p = n;
    return;
  }
  prime_iterator_destroy(iter);
#ifdef NSMALL_PRIMES
  /* In small area? */
  if (n < NSMALL_PRIMES) {
    UV pc = pcount(n);
    iter->segment_start = pc-1;
    iter->p = (pc == 0)  ?  2  :  small_primes[pc-1];
  } else
#endif
  if (n <= primary_limit) { /* Is it inside the primary cache range? */
    iter->p = n;
  } else { /* Sieve this range */
    UV lod, hid;
    lod = n/30;
    hid = lod + SEGMENT_SIZE;
    New(0, iter->segment_mem, SEGMENT_SIZE, unsigned char );
    iter->segment_start = lod * 30;
    iter->segment_bytes = SEGMENT_SIZE;
    if (!sieve_segment((unsigned char*)iter->segment_mem, lod, hid, primary_sieve, primary_limit))
      croak("Could not segment sieve");
    iter->p = n;
  }
}

UV prime_iterator_next(prime_iterator *iter)
{
  UV lod, hid, seg_beg, seg_end;
  const unsigned char* sieve;
  UV n = iter->p;

#ifdef NSMALL_PRIMES
  if (n < NSMALL_PRIMES) {
    iter->p = small_primes[++iter->segment_start];
    return iter->p;
  }
#else
  if (n < 11) {
    switch (n) {
      case 0: case 1:  iter->p =  2; break;
      case 2:          iter->p =  3; break;
      case 3: case 4:  iter->p =  5; break;
      case 5: case 6:  iter->p =  7; break;
      default:         iter->p = 11; break;
    }
    return iter->p;
  }
#endif

  /* Primary sieve */
  if (primary_sieve != 0 && n < 30*PRIMARY_SIZE) {
    n = next_prime_in_segment(primary_sieve, 0, PRIMARY_SIZE, iter->p);
    if (n > 0) {
      iter->p = n;
      return n;
    }
  }

  sieve = iter->segment_mem;
  /* Current segment */
  if (sieve != 0) {
    seg_beg = iter->segment_start;
    seg_end = iter->segment_start + 30*iter->segment_bytes - 1;
    n = next_prime_in_segment(sieve, seg_beg, iter->segment_bytes, iter->p);
    if (n > 0) {
      iter->p = n;
      return n;
    }
    /* Not found in this segment */
    lod = (seg_end+1)/30;
  } else {
    lod = PRIMARY_SIZE;
    New(0, sieve, SEGMENT_SIZE, unsigned char );
  }

  hid = lod + SEGMENT_SIZE - 1;
  iter->segment_start = lod * 30;
  iter->segment_bytes = SEGMENT_SIZE;
  seg_beg = iter->segment_start;
  seg_end = iter->segment_start + 30*iter->segment_bytes - 1;

  if (!sieve_segment((unsigned char*)sieve, lod, hid, primary_sieve, primary_limit))
    croak("Could not segment sieve from %"UVuf" to %"UVuf, seg_beg, seg_end);
  iter->segment_mem = sieve;

  n = next_prime_in_segment(sieve, seg_beg, iter->segment_bytes, seg_beg);
  if (n > 0) {
    iter->p = n;
    return n;
  }
  croak("MPU: segment size too small, could not find prime\n");
}

static int _is_trial_prime(UV n)
{
  UV i = 7;
  UV limit = (UV)sqrt(n);
  while (1) {   /* trial division, skipping multiples of 2/3/5 */
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 4;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 2;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 4;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 2;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 4;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 6;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 2;
    if (i > limit) break;  if ((n % i) == 0) return 0;  i += 6;
  }
  return 1;
}

int prime_iterator_isprime(prime_iterator *iter, UV n)
{
  if (n < 11) {
    switch (n) {
      case 2: case 3: case 5: case 7:  return 1;  break;
      default: break;
    }
    return 0;
  }

  /* Primary sieve */
  if (primary_sieve != 0 && n <= primary_limit) {
    UV d = n/30;
    UV m = n - d*30;
    unsigned char mtab = masktab30[m];
    return mtab && !(primary_sieve[d] & mtab);
  }

  /* Current segment */
  if (iter->segment_mem != 0) {
    int isp = is_prime_in_segment(iter->segment_mem, iter->segment_start, iter->segment_bytes, n);
    if (isp >= 0) return isp;
  }

  /* Out of segment range, can't answer.  Try simple divisibility */
  {
    UV d = n/30;
    UV m = n - d*30;
    unsigned char mtab = masktab30[m];
    if (mtab == 0)  return 0;
    return _is_trial_prime(n);
  }
}

UV* sieve_to_n(UV n, UV* count)
{
  UV pi_max, max_buf, i, p, pi;
  const unsigned char* sieve;
  UV* primes;

#ifdef NSMALL_PRIMES
  if (small_primes != 0 && n < NSMALL_PRIMES) {
    pi = pcount(n);
    New(0, primes, pi, UV);
    for (i = 0; i < pi; i++)  primes[i] = small_primes[i];
    if (count != 0) *count = pi;
    return primes;
  }
#endif

  pi_max = (n < 67)     ? 18
           : (n < 355991) ? 15+(n/(log(n)-1.09))
           : (n/log(n)) * (1.0+1.0/log(n)+2.51/(log(n)*log(n)));
  New(0, primes, pi_max + 10, UV);

  pi = 0;
  primes[pi++] =  2; primes[pi++] =  3; primes[pi++] =  5; primes[pi++] =  7;
  primes[pi++] = 11; primes[pi++] = 13; primes[pi++] = 17; primes[pi++] = 19;
  primes[pi++] = 23; primes[pi++] = 29;

  if (primary_sieve != 0 && n < 30*PRIMARY_SIZE)
    sieve = primary_sieve;
  else
    sieve = sieve_erat30(n);
  max_buf = (n/30) + ((n%30) != 0);
  for (i = 1, p = 30;   i < max_buf;   i++, p += 30) {
    UV c = sieve[i];
    if (!(c &   1)) primes[pi++] = p+ 1;
    if (!(c &   2)) primes[pi++] = p+ 7;
    if (!(c &   4)) primes[pi++] = p+11;
    if (!(c &   8)) primes[pi++] = p+13;
    if (!(c &  16)) primes[pi++] = p+17;
    if (!(c &  32)) primes[pi++] = p+19;
    if (!(c &  64)) primes[pi++] = p+23;
    if (!(c & 128)) primes[pi++] = p+29;
  }
  while (pi > 0 && primes[pi-1] > n) pi--;
  if (sieve != primary_sieve) Safefree(sieve);
  if (count != 0) *count = pi;
  return primes;
}
