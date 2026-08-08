#include <math.h>
#include <gmp.h>

#include "ptypes.h"
#include "legendre_phi.h"
#include "gmp_main.h"
#include "prime_iterator.h"
#include "utility.h"

#define FUNC_isqrt 1
#include "misc_ui.h"

#define FAST_DIV(x,y) \
  ( ((x) <= 4294967295U) ? (uint32_t)(x)/(uint32_t)(y) : (x)/(y) )

#define PHIC  6U
#define PHIS 15U
#define PHIS_XMIN (_snth[PHIS+1]-1U)
#define PHIR 35U
#define PHIZ 25U
#define PHIZ_ROUTE_MAXA 75U
#define PHIZ_ROUTE_MINN (UVCONST(1) << 61)

static const int8_t _coprime_idx210[210]={-1,0,0,0,0,0,0,0,0,0,0,1,1,2,2,2,2,3,3,4,4,4,4,5,5,5,5,5,5,6,6,7,7,7,7,7,7,8,8,8,8,9,9,10,10,10,10,11,11,11,11,11,11,12,12,12,12,12,12,13,13,14,14,14,14,14,14,15,15,15,15,16,16,17,17,17,17,17,17,18,18,18,18,19,19,19,19,19,19,20,20,20,20,20,20,20,20,21,21,21,21,22,22,23,23,23,23,24,24,25,25,25,25,26,26,26,26,26,26,26,26,27,27,27,27,27,27,28,28,28,28,29,29,29,29,29,29,30,30,31,31,31,31,32,32,32,32,32,32,33,33,34,34,34,34,34,34,35,35,35,35,35,35,36,36,36,36,37,37,38,38,38,38,39,39,39,39,39,39,40,40,41,41,41,41,41,41,42,42,42,42,43,43,44,44,44,44,45,45,46,46,46,46,46,46,46,46,46,46,47};
static UV _toindex210(UV x) {
  UV q = x / 210, r = x % 210;
  return 48 * q + _coprime_idx210[r];
}

static const unsigned char _snth[50+1] = {
  0,2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,
  101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,
  193,197,199,211,223,227,229
};

static const uint8_t _s3[30] = {0,1,1,1,1,1,1,2,2,2,2,3,3,4,4,4,4,5,5,6,6,6,6,7,7,7,7,7,7,8};
static const uint8_t _s4[210]= {0,1,1,1,1,1,1,1,1,1,1,2,2,3,3,3,3,4,4,5,5,5,5,6,6,6,6,6,6,7,7,8,8,8,8,8,8,9,9,9,9,10,10,11,11,11,11,12,12,12,12,12,12,13,13,13,13,13,13,14,14,15,15,15,15,15,15,16,16,16,16,17,17,18,18,18,18,18,18,19,19,19,19,20,20,20,20,20,20,21,21,21,21,21,21,21,21,22,22,22,22,23,23,24,24,24,24,25,25,26,26,26,26,27,27,27,27,27,27,27,27,28,28,28,28,28,28,29,29,29,29,30,30,30,30,30,30,31,31,32,32,32,32,33,33,33,33,33,33,34,34,35,35,35,35,35,35,36,36,36,36,36,36,37,37,37,37,38,38,39,39,39,39,40,40,40,40,40,40,41,41,42,42,42,42,42,42,43,43,43,43,44,44,45,45,45,45,46,46,47,47,47,47,47,47,47,47,47,47,48};

static UV tablephi(UV x, uint32_t a)
{
  switch (a) {
    case 0: return x;
    case 1: return x-x/2;
    case 2: return x-x/2-x/3+x/6;
    case 3: return (x/    30U) *     8U + _s3[x %     30U];
    case 4: return (x/   210U) *    48U + _s4[x %    210U];
    case 5: {
              UV xp  = x / 11U;
              return ((x /210) * 48 + _s4[x  % 210]) -
                     ((xp/210) * 48 + _s4[xp % 210]);
             }
    case 6:
    default:
            {
              UV xp  = x / 11U;
              UV x2  = x / 13U;
              UV x2p = x2 / 11U;
              return ((x  /210) * 48 + _s4[x  % 210]) -
                     ((xp /210) * 48 + _s4[xp % 210]) -
                     ((x2 /210) * 48 + _s4[x2 % 210]) +
                     ((x2p/210) * 48 + _s4[x2p% 210]);
            }
  }
}

static void tablephi4_z(mpz_t r, const mpz_t x, mpz_t q)
{
  unsigned long rem = mpz_fdiv_q_ui(q, x, 210);
  mpz_mul_ui(r, q, 48);
  mpz_add_ui(r, r, _s4[rem]);
}

static void tablephi_z(mpz_t r, const mpz_t x, uint32_t a,
                       mpz_t t, mpz_t u)
{
  unsigned long rem;
  switch (a) {
    case 0:
      mpz_set(r, x);
      break;
    case 1:
      mpz_fdiv_q_2exp(t, x, 1);
      mpz_sub(r, x, t);
      break;
    case 2:
      mpz_set(r, x);
      mpz_fdiv_q_2exp(t, x, 1);
      mpz_sub(r, r, t);
      mpz_fdiv_q_ui(t, x, 3);
      mpz_sub(r, r, t);
      mpz_fdiv_q_ui(t, x, 6);
      mpz_add(r, r, t);
      break;
    case 3:
      rem = mpz_fdiv_q_ui(t, x, 30);
      mpz_mul_ui(r, t, 8);
      mpz_add_ui(r, r, _s3[rem]);
      break;
    case 4:
      tablephi4_z(r, x, t);
      break;
    case 5:
      tablephi4_z(r, x, t);
      mpz_fdiv_q_ui(t, x, 11);
      tablephi4_z(u, t, t);
      mpz_sub(r, r, u);
      break;
    case 6:
    default:
      tablephi4_z(r, x, t);
      mpz_fdiv_q_ui(t, x, 11);
      tablephi4_z(u, t, t);
      mpz_sub(r, r, u);
      mpz_fdiv_q_ui(t, x, 13);
      tablephi4_z(u, t, t);
      mpz_sub(r, r, u);
      mpz_fdiv_q_ui(t, x, 143);
      tablephi4_z(u, t, t);
      mpz_add(r, r, u);
      break;
  }
}

static void phi_recurse_small_z(mpz_t r, const mpz_t x, uint32_t a,
                                mpz_t t, mpz_t u)
{
  mpz_t xp, sub;

  if (mpz_sgn(x) < 1) {
    mpz_set_ui(r, 0);
    return;
  }
  if (a <= PHIC) {
    tablephi_z(r, x, a, t, u);
    return;
  }
  if (mpz_cmp_ui(x, _snth[a]) < 0) {
    mpz_set_ui(r, 1);
    return;
  }

  mpz_init(xp);
  mpz_init(sub);
  phi_recurse_small_z(r, x, a-1, t, u);
  mpz_fdiv_q_ui(xp, x, _snth[a]);
  phi_recurse_small_z(sub, xp, a-1, t, u);
  mpz_sub(r, r, sub);
  mpz_clear(sub);
  mpz_clear(xp);
}

static UV phi_small(UV x, uint32_t a) {
  UV sum = 0, xpos[1025], xneg[1025];
  uint32_t i, npos, nneg;

  if (a < 4) {
    return (a==0) ? x :
           (a==1) ? x-x/2 :
           (a==2) ? x-x/2-x/3+x/6
                  : (x/30U) * 8U + _s3[x % 30U];
  }
  MPUassert(a <= PHIS, "phi_small: a too large");
  if (x < _snth[a+1]) return (x>0);

  for (npos = nneg = 0, xpos[npos++] = x;  a > 4U;  a--) {
    uint32_t oneg = nneg,  opos = npos;
    for (i = 0; i < opos; i++)
      if (xpos[i] >= _snth[a])
        xneg[nneg++] = xpos[i]/_snth[a];
    for (i = 0; i < oneg; i++)
      if (xneg[i] >= _snth[a])
        xpos[npos++] = xneg[i]/_snth[a];
  }
  for (i = 0; i < npos; i++)
    sum += (xpos[i]/210U)*48U + _s4[xpos[i] % 210U];
  for (i = 0; i < nneg; i++)
    sum -= (xneg[i]/210U)*48U + _s4[xneg[i] % 210U];
  return sum;
}

static UV phi_recurse_small(UV x, UV a) {
  UV sum, i, xp, p;

  if (x < 1 || a >= x) return (x > 0);
  if (a <= PHIS) return phi_small(x, a);
  if (x <= PHIS_XMIN) return 1;

  sum = phi_small(x, PHIS);
  p = _snth[PHIS];
  for (i = PHIS+1; i <= a; i++) {
    p = next_prime_ui(p);
    xp = FAST_DIV(x,p);
    if (xp > 0)
      sum -= phi_recurse_small(xp, i-1);
  }
  return sum;
}

#define PHICACHEA 512
typedef struct
{
  uint32_t  siz[PHICACHEA];
  uint16_t *val[PHICACHEA];
  uint32_t  xlim;
} phi_cache_t;

static phi_cache_t* phi_cache_create(uint32_t xlim) {
  phi_cache_t *cache;
  int a;
  New(0, cache, 1, phi_cache_t);
  for (a = 0; a < PHICACHEA; a++) {
    cache->val[a] = 0;
    cache->siz[a] = 0;
  }
  cache->xlim = (xlim < 0xFFFFFFFFU) ? xlim : xlim-1;
  return cache;
}

static void phi_cache_destroy(phi_cache_t* cache) {
  int a;
  for (a = 0; a < PHICACHEA; a++) {
    if (cache->val[a] != 0)
      Safefree(cache->val[a]);
  }
  Safefree(cache);
}

static void phi_cache_insert(uint32_t x, uint32_t a, IV sum, phi_cache_t* cache) {
  uint32_t i, newsize;
  if (sum < 0) sum = -sum;
  if (sum > 65535) return;
  if (x >= cache->siz[a]) {
    newsize = (x >= 0xFFFFFFFFUL-32)  ?  0xFFFFFFFFUL-1  :  x+32;
    if (cache->val[a] == 0) {
      Newz(0, cache->val[a], newsize, uint16_t);
    } else {
      Renew(cache->val[a], newsize, uint16_t);
      for (i = cache->siz[a]; i < newsize; i++)
        cache->val[a][i] = 0;
    }
    cache->siz[a] = newsize;
  }
  cache->val[a][x] = (uint16_t) sum;
}

typedef struct {
  const uint32_t* primes;
  uint32_t lastidx;
  phi_cache_t* cachephi;
} phidata_t;

static UV nth_prime_upper_ui(UV n)
{
  double fn, upper;
  if (n <= 50) return _snth[n];
  fn = (double)n;
  upper = fn * (log(fn) + log(log(fn)));
  if (upper >= (double)UV_MAX-32)
    croak("legendre_phi: nth prime upper overflow");
  return (UV)ceil(upper) + 32;
}

static UV prime_count_ui(UV n)
{
  mpz_t zn, pc;
  UV ret;
  mpz_init(zn);
  mpz_set_uv(zn, n);
  mpz_init(pc);
  prime_count(pc, zn);
  if (!mpz_fits_uv_p(pc))
    croak("legendre_phi: prime count overflow");
  ret = mpz_get_uv(pc);
  mpz_clear(pc);
  mpz_clear(zn);
  return ret;
}

static UV prime_count_upper_ui(UV n)
{
  mpz_t zn, pc;
  UV ret;
  mpz_init(zn);
  mpz_set_uv(zn, n);
  mpz_init(pc);
  prime_count_upper(pc, zn);
  ret = mpz_fits_uv_p(pc) ? mpz_get_uv(pc) : UV_MAX;
  mpz_clear(pc);
  mpz_clear(zn);
  return ret;
}

static uint32_t prime_count_with_primes(UV n, const uint32_t *primes, uint32_t lastidx)
{
  uint32_t lo = 0, hi = lastidx;

  if (lastidx > 0 && n <= primes[lastidx]) {
    while (lo < hi) {
      uint32_t mid = lo + (hi-lo+1)/2;
      if (primes[mid] <= n) lo = mid;
      else                  hi = mid-1;
    }
    return lo;
  }
  return (uint32_t)prime_count_ui(n);
}

static uint32_t make_prime_table(uint32_t **primes, UV n)
{
  UV npr, i, *uvprimes;
  uint32_t *p32;

  if (n > UINT32_MAX)
    croak("legendre_phi: prime table too large");
  uvprimes = sieve_to_n(n, &npr);
  if (npr > UINT32_MAX-1)
    croak("legendre_phi: too many primes");
  New(0, p32, npr+1, uint32_t);
  p32[0] = 0;
  for (i = 0; i < npr; i++)
    p32[i+1] = (uint32_t)uvprimes[i];
  Safefree(uvprimes);
  *primes = p32;
  return (uint32_t)npr;
}

static phidata_t* phidata_create(const uint32_t* primes, uint32_t lastidx, UV x)
{
  phidata_t *d;
  uint32_t xlim = (UV) pow((double)x, 1.0/2.70);
  if (xlim < 256) xlim = 256;

  New(0, d, 1, phidata_t);
  d->primes = primes;
  d->lastidx = lastidx;
  d->cachephi = phi_cache_create(xlim);
  return d;
}

static void phidata_destroy(phidata_t *d)
{
  phi_cache_destroy(d->cachephi);
  Safefree(d);
}

#define PHI_IS_X_SMALL(x, a) \
  ( ((x) <= primes[d->lastidx]) && ((x) < (UV)primes[a+1] * primes[a+1]) )
#define PHI_PRIMECOUNT(x) \
  prime_count_with_primes((x), primes, d->lastidx)

static void mpz_add_sign_ui(mpz_t r, UV v, int sign, mpz_t t)
{
  mpz_set_uv(t, v);
  if (sign > 0) mpz_add(r, r, t);
  else          mpz_sub(r, r, t);
}

static void prime_count_z(mpz_t r, const mpz_t n)
{
  if (mpz_fits_uv_p(n)) {
    mpz_set_uv(r, prime_count_ui(mpz_get_uv(n)));
  } else {
    prime_count(r, n);
  }
}

static int phi_z_is_small(const mpz_t x, UV a, phidata_t *d)
{
  const uint32_t* const primes = d->primes;
  UV p, p2, xu;

  if (!mpz_fits_uv_p(x) || a >= d->lastidx)
    return 0;
  p = primes[a+1];
  if (p > UV_MAX / p)
    return 0;
  p2 = p * p;
  xu = mpz_get_uv(x);
  return xu < p2;
}

static IV _phi3(UV x, UV a, int sign, phidata_t *d);

static void _phi3_z(mpz_t r, const mpz_t x, UV a, int sign, phidata_t *d)
{
  const uint32_t* const primes = d->primes;
  mpz_t xp, tmp, term, u;
  UV c, i, iters;

  if (a >= d->lastidx)
    croak("legendre_phi: failed to generate enough primes");

  if (mpz_cmp_ui(x, primes[a+1]) < 0) {
    mpz_set_si(r, sign);
    return;
  } else if (a <= PHIC) {
    mpz_init(tmp);
    mpz_init(term);
    tablephi_z(r, x, (uint32_t)a, tmp, term);
    if (sign < 0)
      mpz_neg(r, r);
    mpz_clear(term);
    mpz_clear(tmp);
    return;
  } else if (mpz_fits_uv_p(x)) {
    mpz_set_iv(r, _phi3(mpz_get_uv(x), a, sign, d));
    return;
  }

  iters = a;
  if (a <= UV_MAX / a) {
    UV a2 = a * a;
    if (mpz_cmp_ui(x, a2) < 0) {
      mpz_t sx, pc;
      mpz_init(sx);
      mpz_init(pc);
      mpz_sqrt(sx, x);
      prime_count_z(pc, sx);
      if (!mpz_fits_uv_p(pc))
        croak("legendre_phi: prime count overflow");
      iters = mpz_get_uv(pc);
      mpz_clear(pc);
      mpz_clear(sx);
    }
  }

  c = (iters > PHIC) ? PHIC : iters;
  mpz_init(xp);
  mpz_init(tmp);
  mpz_init(term);
  mpz_init(u);
  tablephi_z(r, x, (uint32_t)c, tmp, u);
  if (iters >= a) mpz_add_ui(r, r, iters - a);
  else            mpz_sub_ui(r, r, a - iters);
  if (sign < 0)
    mpz_neg(r, r);

  if (c < iters) {
    mpz_fdiv_q_ui(xp, x, primes[c+1]);
    tablephi_z(term, xp, (uint32_t)c, tmp, u);
    if (sign > 0) mpz_sub(r, r, term);
    else          mpz_add(r, r, term);
  }

  for (i = c+1; i < iters; i++) {
    mpz_fdiv_q_ui(xp, x, primes[i+1]);
    if (phi_z_is_small(xp, i, d))
      break;
    _phi3_z(term, xp, i, -sign, d);
    mpz_add(r, r, term);
  }
  for (; i < iters; i++) {
    mpz_fdiv_q_ui(xp, x, primes[i+1]);
    if (mpz_cmp_ui(xp, primes[i+1]) < 0)
      break;
    prime_count_z(term, xp);
    mpz_sub_ui(term, term, i-1);
    if (sign > 0) mpz_sub(r, r, term);
    else          mpz_add(r, r, term);
  }
  if (i < iters)
    mpz_add_sign_ui(r, iters - i, -sign, tmp);

  mpz_clear(u);
  mpz_clear(term);
  mpz_clear(tmp);
  mpz_clear(xp);
}

static IV _phi3(UV x, UV a, int sign, phidata_t *d)
{
  const uint32_t* const primes = d->primes;
  phi_cache_t* pcache = d->cachephi;
  UV mapx;

  if (x < primes[a+1])
    return sign;
  else if (a <= PHIC)
    return sign * tablephi(x,a);
  else if (PHI_IS_X_SMALL(x,a))
    return sign * (PHI_PRIMECOUNT(x) - a + 1);

  mapx = (a < PHICACHEA)  ?  _toindex210(x)  :  0;

  if (a < PHICACHEA && mapx < pcache->siz[a]) {
    IV v = pcache->val[a][mapx];
    if (v != 0)
      return sign * v;
  }
  {
    UV xp, i, iters = ((UV)a*a > x)  ?  PHI_PRIMECOUNT(isqrt(x))  :  a;
    UV c = (iters > PHIC) ? PHIC : iters;
    IV sum = sign * (iters - a + tablephi(x,c));

    if (c < iters)
      sum += -sign * tablephi(FAST_DIV(x,primes[c+1]), c);
    for (i = c+1; i < iters; i++) {
      xp = FAST_DIV(x,primes[i+1]);
      if (PHI_IS_X_SMALL(xp,i))
        break;
      sum += _phi3(xp, i, -sign, d);
    }
    for (; i < iters; i++) {
      xp = FAST_DIV(x,primes[i+1]);
      if (xp < primes[i+1])
        break;
      sum += -sign * (PHI_PRIMECOUNT(xp) - i + 1);
    }
    if (i < iters)
      sum += -sign * (iters - i);

    if (a < PHICACHEA && mapx <= pcache->xlim)
      phi_cache_insert(mapx, a, sum, pcache);
    return sum;
  }
}

static UV phi_recurse(UV x, UV a)
{
  uint32_t* primes;
  uint32_t lastidx;
  UV primes_to_n, sqrtx, sum = 1;

  if (x < 1 || a >= x) return (x > 0);
  if (a <= PHIS) return phi_small(x, a);
  if (x <= PHIS_XMIN) return 1;
  if (a > 203280221) croak("legendre_phi: 64-bit phi out of range");

  sqrtx = isqrt(x);
  primes_to_n = nth_prime_upper_ui(a);
  if (sqrtx > primes_to_n) primes_to_n = sqrtx;
  lastidx = make_prime_table(&primes, primes_to_n);
  if (lastidx < a)
    croak("legendre_phi: failed to generate enough primes");

  if (primes[a] < x) {
    phidata_t *d = phidata_create(primes, lastidx, x);
    sum = (UV) _phi3(x, a-1, 1, d) - (UV) _phi3(x/primes[a], a-1, 1, d);
    phidata_destroy(d);
  }

  Safefree(primes);
  return sum;
}

static void phi_recurse_z(mpz_t r, const mpz_t x, UV a)
{
  uint32_t* primes;
  uint32_t lastidx;
  UV primes_to_n;
  mpz_t xp, t;

  if (mpz_sgn(x) < 1) {
    mpz_set_ui(r, 0);
    return;
  }
  if (mpz_cmp_ui(x, a) <= 0) {
    mpz_set_ui(r, 1);
    return;
  }
  if (a <= PHIZ) {
    mpz_init(t);
    mpz_init(xp);
    phi_recurse_small_z(r, x, (uint32_t)a, t, xp);
    mpz_clear(xp);
    mpz_clear(t);
    return;
  }

  primes_to_n = nth_prime_upper_ui(a);
  lastidx = make_prime_table(&primes, primes_to_n);
  if (lastidx < a)
    croak("legendre_phi: failed to generate enough primes");

  if (mpz_cmp_ui(x, primes[a]) <= 0) {
    mpz_set_ui(r, 1);
  } else {
    phidata_t *d = phidata_create(primes, lastidx, mpz_fits_uv_p(x) ? mpz_get_uv(x) : UVCONST(100000000));
    mpz_init(xp);
    mpz_init(t);
    _phi3_z(r, x, a-1, 1, d);
    mpz_fdiv_q_ui(xp, x, primes[a]);
    _phi3_z(t, xp, a-1, 1, d);
    mpz_sub(r, r, t);
    mpz_clear(t);
    mpz_clear(xp);
    phidata_destroy(d);
  }

  Safefree(primes);
}

static UV legendre_phi_ui(UV x, UV a)
{
  UV sqrtx = isqrt(x);

  if (x < 1 || a >= x) return (x > 0);
  if (x <= PHIC || a <= PHIC)  return tablephi(x, (a > PHIC) ? PHIC : a);

  if (a > (x >> 1))
    return 1;
  if (a >= sqrtx || a > 203280221) {
    UV pc = prime_count_ui(x);
    return (a >= pc)  ?  1  :  pc - a + 1;
  }

  if (a <= PHIS)  return phi_small(x, a);
  if (a <= PHIR)  return phi_recurse_small(x, a);

  if (prime_count_upper_ui(x) <= a)
    return 1;
  if (prime_count_upper_ui(sqrtx) < a) {
    UV pc = prime_count_ui(x);
    return (a >= pc)  ?  1  :  pc - a + 1;
  }

  return phi_recurse(x, a);
}

void legendre_phi(mpz_t r, const mpz_t n, const mpz_t a)
{
  mpz_t sqrtn, t;
  UV au;

  if (r == n || r == a) {
    mpz_t N, A;
    mpz_init_set(N, n);
    mpz_init_set(A, a);
    legendre_phi(r, N, A);
    mpz_clear(A);
    mpz_clear(N);
    return;
  }

  if (mpz_sgn(n) == 0) {
    mpz_set_ui(r, 0);
    return;
  }
  if (mpz_sgn(a) == 0) {
    mpz_set(r, n);
    return;
  }
  if (mpz_cmp(a, n) >= 0) {
    mpz_set_ui(r, 1);
    return;
  }

  mpz_init(t);
  mpz_fdiv_q_2exp(t, n, 1);
  if (mpz_cmp(a, t) > 0) {
    mpz_set_ui(r, 1);
    mpz_clear(t);
    return;
  }

  if (mpz_fits_uv_p(n) && mpz_fits_uv_p(a)) {
    UV nu = mpz_get_uv(n);
    UV auv = mpz_get_uv(a);
#if BITS_PER_WORD >= 64
    /* The UV cached path can spend heavily near UV_MAX without phi_walk. */
    if (nu >= PHIZ_ROUTE_MINN && auv > PHIZ && auv <= PHIZ_ROUTE_MAXA) {
      phi_recurse_z(r, n, auv);
      mpz_clear(t);
      return;
    }
#endif
    mpz_set_uv(r, legendre_phi_ui(nu, auv));
    mpz_clear(t);
    return;
  }

  if (mpz_fits_uv_p(a))
    au = mpz_get_uv(a);
  else
    au = 0;

  mpz_init(sqrtn);
  mpz_sqrt(sqrtn, n);
  if (mpz_cmp(a, sqrtn) >= 0) {
    mpz_clear(sqrtn);
    mpz_clear(t);
    croak("legendre_phi: arguments too large for processing");
  }
  mpz_clear(sqrtn);
  mpz_clear(t);

  if (au > 0) {
    phi_recurse_z(r, n, au);
    return;
  }

  croak("legendre_phi: arguments too large for processing");
}
