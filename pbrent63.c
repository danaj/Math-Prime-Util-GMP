/* Native 63-bit Pollard-Rho-Brent for x86-64. */

#include <gmp.h>
#include "ptypes.h"
#include "pbrent63.h"

#if BITS_PER_WORD == 64 && HAVE_STD_U64 && defined(__GNUC__) && defined(__x86_64__)

#define FUNC_gcd_ui 1
#include "utility.h"
#include "misc_ui.h"

static INLINE UV mpz_getuv(const mpz_t n) {
  UV v = mpz_getlimbn(n,0);
  if (GMP_LIMB_BITS < 64 || sizeof(mp_limb_t) < sizeof(UV))
    v |= ((UV)mpz_getlimbn(n,1)) << 32;
  return v;
}

int pbrent63(const mpz_t n, mpz_t f, UV rounds) {
  UV facs[2];
  int nfactors;

  if (mpz_sizeinbase(n,2) > 63) return 0;
  nfactors = uvpbrent63(mpz_getuv(n), facs, rounds, 1);
  if (nfactors < 2) return 0;
  /* Smallest factor of 64-bit n always fits in 32-bit */
  mpz_set_ui(f, (facs[0] < facs[1]) ? facs[0] : facs[1]);
  return 1;
}

/* Trimmed out all the extra stuff and the 64-bit variation */

#define mont_get1(n)              _u64div(1,n)
/* Must have npi = mont_inverse(n), mont1 = mont_get1(n) */
#define mont_geta(a,n)            mulmod(a,mont1,n)
#define mont_mulmod(a,b,n)        _mulredc63(a,b,n,npi)

static INLINE uint64_t mont_inverse(const uint64_t n) {
  uint64_t ret = (3*n) ^ 2;
  ret *= (uint64_t)2 - n * ret;
  ret *= (uint64_t)2 - n * ret;
  ret *= (uint64_t)2 - n * ret;
  ret *= (uint64_t)2 - n * ret;
  return (uint64_t)0 - ret;
}
/* MULREDC asm from Ben Buhrow */
static INLINE uint64_t _mulredc63(uint64_t a, uint64_t b, uint64_t n, uint64_t npi) {
    asm("mulq %2 \n\t"
        "movq %%rax, %%r10 \n\t"
        "movq %%rdx, %%r11 \n\t"
        "mulq %3 \n\t"
        "mulq %4 \n\t"
        "addq %%r10, %%rax \n\t"
        "adcq %%r11, %%rdx \n\t"
        "xorq %%rax, %%rax \n\t"
        "subq %4, %%rdx \n\t"
        "cmovc %4, %%rax \n\t"
        "addq %%rdx, %%rax \n\t"
        : "=a"(a)
        : "0"(a), "r"(b), "r"(npi), "r"(n)
        : "rdx", "r10", "r11", "cc");
  return a;
}
static INLINE uint64_t _u64div(uint64_t c, uint64_t n) {
  asm("divq %4"
      : "=a"(c), "=d"(n)
      : "1"(c), "0"(0), "r"(n));
  return n;
}
static INLINE UV mulmod(UV a, UV b, UV n) {
  UV d, dummy;                    /* d will get a*b mod c */
  asm ("mulq %3\n\t"              /* mul a*b -> rdx:rax */
       "divq %4\n\t"              /* (a*b)/c -> quot in rax remainder in rdx */
       :"=a"(dummy), "=&d"(d)     /* output */
       :"a"(a), "r"(b), "r"(n)    /* input */
       :"cc"                      /* mulq and divq can set conditions */
      );
  return d;
}
static INLINE UV addmod(UV a, UV b, UV n) {
  UV t = a-n;
  a += b;
  asm ("add %2, %1\n\t"    /* t := t + b */
       "cmovc %1, %0\n\t"  /* if (carry) a := t */
       :"+r" (a), "+&r" (t)
       :"r" (b)
       :"cc"
      );
  return a;
}

#define ABSDIFF(x,y)  (x>y) ? x-y : y-x
/* Brent's modifications to Pollard's Rho. */
int uvpbrent63(UV n, UV *factors, UV rounds, UV a)
{
  UV const nbits = BITS_PER_WORD - __builtin_ctzll(n);
  const UV inner = (nbits <= 31) ? 32 : (nbits <= 35) ? 64 : (nbits <= 40) ? 160 : (nbits <= 52) ? 256 : 320;
  UV f, m, r, rleft, Xi, Xm, Xs;
  int irounds, fails = 6;
  const uint64_t npi = mont_inverse(n),  mont1 = mont_get1(n);

  if (n <= 3) { factors[0] = n; return 1; }
  if (!(n&1)) { factors[0] = 2; factors[1] = n/2; return 2; }

  r = f = 1;
  Xi = Xm = Xs = mont1;
  a = mont_geta(a,n);

  while (rounds > 0) {
    rleft = (r > rounds) ? rounds : r;
    Xm = Xi;
    /* Do rleft rounds, inner at a time */
    while (rleft > 0) {
      irounds = (rleft > (UV)inner) ? inner : rleft;
      rleft -= irounds;
      rounds -= irounds;
      Xs = Xi;
      Xi = mont_mulmod(Xi,Xi+a,n);
      m = ABSDIFF(Xi,Xm);
      while (--irounds > 0) {
        Xi = mont_mulmod(Xi,Xi+a,n);
        f = ABSDIFF(Xi,Xm);
        m = mont_mulmod(m, f, n);
      }
      f = gcd_ui(m, n);
      if (f != 1)
        break;
    }
    /* If f == 1, then we didn't find a factor.  Move on. */
    if (f == 1) {
      r *= 2;
      continue;
    }
    if (f == n) {  /* back up, with safety */
      Xi = Xs;
      do {
        Xi = mont_mulmod(Xi,Xi+a,n);
        m = ABSDIFF(Xi,Xm);
        f = gcd_ui(m, n);
      } while (f == 1 && r-- != 0);
    }
    if (f == 0 || f == n) {
      if (fails-- <= 0) break;
      Xi = Xm = mont1;
      a = addmod(a, mont_geta(11,n), n);
      continue;
    }
    factors[0] = f;
    factors[1] = n / f;
    MPUassert( factors[0] * factors[1] == n , "incorrect factoring");
    return 2;
  }
  factors[0] = n;
  return 1;
}

#else /* no 64-bit gcc x86-64 */
int pbrent63(const mpz_t n, mpz_t f, UV rounds) { return 0; }
int uvpbrent63(UV n, UV *factors, UV rounds, UV a) { factors[0] = n; return 1; }
#endif
