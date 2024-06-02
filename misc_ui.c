/* Native functions, either ui or UV.
 *
 * Utility as well as some more specialized.
 */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <gmp.h>
#include <math.h>

#include "ptypes.h"

#define FUNC_ipow 1
#define FUNC_isqrt 1
#define FUNC_icbrt 1
#include "misc_ui.h"


#if BITS_PER_WORD == 64
  #define MPU_MAX_POW3 40
  static const uint32_t root_max[1+MPU_MAX_POW3] = {0,0,4294967295U,2642245,65535,7131,1625,565,255,138,84,56,40,30,23,19,15,13,11,10,9,8,7,6,6,5,5,5,4,4,4,4,3,3,3,3,3,3,3,3,3};
  static const unsigned char _debruijn64[64] = {
    63, 0,58, 1,59,47,53, 2, 60,39,48,27,54,33,42, 3, 61,51,37,40,49,18,28,20,
    55,30,34,11,43,14,22, 4, 62,57,46,52,38,26,32,41, 50,36,17,19,29,10,13,21,
    56,45,25,31,35,16, 9,12, 44,24,15, 8,23, 7, 6, 5 };
  unsigned int log2_ui(UV n) {
    if (n == 0) return 0;
    n |= n >> 1;   n |= n >> 2;   n |= n >> 4;
    n |= n >> 8;   n |= n >> 16;  n |= n >> 32;
    return _debruijn64[((n-(n>>1))*UVCONST(0x07EDD5E59A4E28C2)) >> 58];
  }
#else
  #define MPU_MAX_POW3 20
  static const uint32_t root_max[1+MPU_MAX_POW3] = {0,0,65535,1625,255,84,40,23,15,11,9,7,6,5,4,4,3,3,3,3,3};
  static const unsigned char _lead_debruijn32[32] = {
    0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30,
    8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31 };
  unsigned int log2_ui(UV n) {
    if (n == 0) return 0;
    n |= n >> 1;   n |= n >> 2;   n |= n >> 4;   n |= n >> 8;   n |= n >> 16;
    return _lead_debruijn32[(n * UVCONST(0x07C4ACDD)) >> 27];
  }
#endif

UV rootint_ui(UV n, UV k)
{
  UV lo, hi, max;
  if (k == 0) return 0;
  if (k == 1) return n;
  if (k == 2) return isqrt(n);
  if (k == 3) return icbrt(n);

  /* Bracket between powers of 2, but never exceed max power so ipow works */
  max = 1 + ((k > MPU_MAX_POW3) ? 2 : root_max[k]);
  lo = UVCONST(1) << (log2_ui(n)/k);
  hi = ((lo*2) < max) ? lo*2 : max;

  /* Binary search */
  while (lo < hi) {
    UV mid = lo + (hi-lo)/2;
    if (ipow(mid,k) <= n) lo = mid+1;
    else                  hi = mid;
  }
  return lo-1;
}
UV logint_ui(UV n, UV b)
{
  /* UV e;  for (e=0; n; n /= b) e++;  return e-1; */
  UV v, e = 0;
  if (b == 2)
    return log2_ui(n);
  if (b > n)
    return 0;
  if (n > UV_MAX/b) {
    n /= b;
    e = 1;
  }
  for (v = b; v <= n; v *= b)
    e++;
  return e;
}

/* Modular inversion: invert a mod p.
 * This implementation from William Hart, using extended gcd.
 */
unsigned long modinverse(unsigned long a, unsigned long p)
{
  long u1 = 1;
  long u3 = a;
  long v1 = 0;
  long v3 = p;
  long t1 = 0;
  long t3 = 0;
  long quot;
  while (v3) {
    quot = u3 - v3;
    if (u3 < (v3<<2)) {
      if (quot < v3) {
        if (quot < 0) {
          t1 = u1; u1 = v1; v1 = t1;
          t3 = u3; u3 = v3; v3 = t3;
        } else {
          t1 = u1 - v1; u1 = v1; v1 = t1;
          t3 = u3 - v3; u3 = v3; v3 = t3;
        }
      } else if (quot < (v3<<1)) {
        t1 = u1 - (v1<<1); u1 = v1; v1 = t1;
        t3 = u3 - (v3<<1); u3 = v3; v3 = t3;
      } else {
        t1 = u1 - v1*3; u1 = v1; v1 = t1;
        t3 = u3 - v3*3; u3 = v3; v3 = t3;
      }
    } else {
      quot = u3 / v3;
      t1 = u1 - v1*quot; u1 = v1; v1 = t1;
      t3 = u3 - v3*quot; u3 = v3; v3 = t3;
    }
 }
 if (u1 < 0) u1 += p;
 return u1;
}
