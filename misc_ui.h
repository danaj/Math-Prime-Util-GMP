#ifndef MPU_MISCUI_H
#define MPU_MISCUI_H

#include <math.h>
#include "ptypes.h"

extern unsigned long modinverse(unsigned long a, unsigned long p);
extern UV rootint_ui(UV n, UV k);
extern UV logint_ui(UV n, UV b);
extern unsigned int log2_ui(UV n);


#if defined(FUNC_isqrt) || defined(FUNC_is_perfect_square)
static UV isqrt(UV n) {
  UV root;
#if BITS_PER_WORD == 32
  if (n >= UVCONST(4294836225)) return UVCONST(65535);
#else
  if (n >= UVCONST(18446744065119617025)) return UVCONST(4294967295);
#endif
  root = (UV) sqrt((double)n);
  while (root*root > n)  root--;
  while ((root+1)*(root+1) <= n)  root++;
  return root;
}
#endif

#if defined(FUNC_icbrt) || defined(FUNC_is_perfect_cube)
static UV icbrt(UV n) {
  UV b, root = 0;
#if BITS_PER_WORD == 32
  int s = 30;
  if (n >= UVCONST(4291015625)) return UVCONST(1625);
#else
  int s = 63;
  if (n >= UVCONST(18446724184312856125)) return UVCONST(2642245);
#endif
  for ( ; s >= 0; s -= 3) {
    root += root;
    b = 3*root*(root+1)+1;
    if ((n >> s) >= b) {
      n -= b << s;
      root++;
    }
  }
  return root;
}
#endif

#if defined(FUNC_ipow)
static UV ipow(UV n, UV k) {
  UV p = 1;
  while (k) {
    if (k & 1) p *= n;
    k >>= 1;
    if (k)     n *= n;
  }
  return p;
}
#endif

#if defined(FUNC_gcd_ui) || defined(FUNC_lcm_ui)
static UV gcd_ui(UV x, UV y) {
  UV t;
  if (y < x) { t = x; x = y; y = t; }
  while (y > 0) {
    t = y;  y = x % y;  x = t;  /* y1 <- x0 % y0 ; x1 <- y0 */
  }
  return x;
}
#endif

#ifdef FUNC_lcm_ui
static UV lcm_ui(UV x, UV y) {
  /* Can overflow if lcm(x,y) > 2^64 (e.g. two primes each > 2^32) */
  return x * (y / gcd_ui(x,y));
}
#endif

#ifdef FUNC_is_perfect_square
static int is_perfect_square(UV n)
{
  /* Step 1, reduce to 18% of inputs */
  uint32_t m = n & 127;
  if ((m*0x8bc40d7d) & (m*0xa1e2f5d1) & 0x14020a)  return 0;
  /* Step 2, reduce to 7% of inputs (mod 99 reduces to 4% but slower) */
  m = n %240; if ((m*0xfa445556) & (m*0x8021feb1) & 0x614aaa0f) return 0;
  /* m = n % 99; if ((m*0x5411171d) & (m*0xe41dd1c7) & 0x80028a80) return 0; */
  /* Step 3, do the square root instead of any more rejections */
  m = isqrt(n);
  return (UV)m*(UV)m == n;
}
#endif

#ifdef FUNC_is_perfect_cube
static int is_perfect_cube(UV n)
{
  uint32_t m;
  m = n % 117; if ((m*833230740) & (m*120676722) & 813764715) return 0;
  m = n % 133; if ((m*76846229) & (m*305817297) & 306336544) return 0;
  m = n % 43; if ((m*193635074) & (m*3653322805U) & 74401) return 0;
  m = n % 37; if ((m*919307198) & (m*3908849845U) & 6665) return 0;
  m = icbrt(n);
  return (UV)m*m*m == n;
}
#endif

#endif
