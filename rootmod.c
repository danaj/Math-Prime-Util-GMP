#include <gmp.h>
#include "ptypes.h"

#include "rootmod.h"
#include "utility.h"
#include "factor.h"
#include "primality.h"
#include "znlog.h"



/******************************************************************************/
/*                               SQRT(N) MOD M                                */
/******************************************************************************/

static int _sqrtmod_return(mpz_t r, const mpz_t a, const mpz_t n, mpz_t t) {
  mpz_sub(t, n, r);
  if (mpz_cmp(t, r) < 0)
    mpz_set(r, t);
  mpz_mul(t, r, r);
  mpz_sub(t, t, a);
  mpz_mod(t, t, n);
  if (mpz_sgn(t) == 0)
    return 1;
  /* (r*r) mod n != a mod n : r is not a square root of a mod n */
  mpz_set_ui(r, 0);
  return 0;
}

/* No aliasing and 4 temp variables passed in. */
static int _sqrtmod_prime(mpz_t x, const mpz_t a, const mpz_t p,
                          mpz_t t, mpz_t q, mpz_t b, mpz_t z)
{
  int r, e, m;

#if 0
  if (mpz_perfect_square_p(a)) {
    mpz_sqrt(x, a);
    mpz_mod(x, x, p);
    return _sqrtmod_return(x, a, p, t);
  }
#endif

  /* Easy cases from page 31 (or Menezes 3.36, 3.37) */
  if (mpz_congruent_ui_p(p, 3, 4)) {
    mpz_add_ui(t, p, 1);
    mpz_tdiv_q_2exp(t, t, 2);
    mpz_powm(x, a, t, p);
    return _sqrtmod_return(x, a, p, t);
  }

  if (mpz_congruent_ui_p(p, 5, 8)) {
    mpz_sub_ui(t, p, 1);
    mpz_tdiv_q_2exp(t, t, 2);
    mpz_powm(q, a, t, p);
    if (mpz_cmp_si(q, 1) == 0) {  /* s = a^((p+3)/8) mod p */
      mpz_add_ui(t, p, 3);
      mpz_tdiv_q_2exp(t, t, 3);
      mpz_powm(x, a, t, p);
    } else {                      /* s = 2a * (4a)^((p-5)/8) mod p */
      mpz_sub_ui(t, p, 5);
      mpz_tdiv_q_2exp(t, t, 3);
      mpz_mul_ui(q, a, 4);
      mpz_powm(x, q, t, p);
      mpz_mul_ui(x, x, 2);
      mpz_mulmod(x, x, a, p, x);
    }
    return _sqrtmod_return(x, a, p, t);
  }

  if (mpz_kronecker(a, p) != 1) {
    /* Possible no solution exists.  Check Euler criterion. */
    mpz_sub_ui(t, p, 1);
    mpz_tdiv_q_2exp(t, t, 1);
    mpz_powm(x, a, t, p);
    if (mpz_cmp_si(x, 1) != 0) {
      mpz_set_ui(x, 0);
      return 0;
    }
  }

  mpz_sub_ui(q, p, 1);
  e = mpz_scan1(q, 0);                 /* Remove 2^e from q */
  mpz_tdiv_q_2exp(q, q, e);
  mpz_set_ui(t, 2);
  while (mpz_kronecker(t, p) != -1) {  /* choose t "at random" */
    mpz_add_ui(t, t, 1);
    if (!mpz_cmp_ui(t,133)) {
      /* If a root of p exists, then our chances are nearly 1/2 that
       * (t|p) = -1.  After 133 tries it seems dubious that a root
       * exists.  It's likely that p is not prime. */
      if (mpz_even_p(p)) { mpz_set_ui(x,0); return 0; }
      /* Euler probable prime test with base t.  (t|p) = 1 or t divides p */
      if (mpz_divisible_p(p, t)) { mpz_set_ui(x,0); return 0; }
      mpz_sub_ui(z, p, 1);  mpz_fdiv_q_2exp(b,z,1);  mpz_powm(z, t, b, p);
      if (mpz_cmp_ui(z,1)) { mpz_set_ui(x,0); return 0; }
      /* Fermat base 2 */
      mpz_set_ui(b,2);  mpz_sub_ui(z, p, 1);  mpz_powm(z, b, z, p);
      if (mpz_cmp_ui(z,1)) { mpz_set_ui(x,0); return 0; }
    }
    if (!mpz_cmp_ui(t,286)) {
      /* Another Euler probable prime test, p not even so t can't divide. */
      mpz_sub_ui(z, p, 1);  mpz_fdiv_q_2exp(b,z,1);  mpz_powm(z, t, b, p);
      if (mpz_cmp_ui(z,1)) { mpz_set_ui(x,0); return 0; }
    }
    if (!mpz_cmp_ui(t,20000)) { mpz_set_ui(x,0); return 0; }
  }
  mpz_powm(z, t, q, p);                     /* Step 1 complete */
  r = e;

  mpz_powm(b, a, q, p);
  mpz_add_ui(q, q, 1);
  mpz_divexact_ui(q, q, 2);
  mpz_powm(x, a, q, p);   /* Done with q, will use it for y now */

  while (mpz_cmp_ui(b, 1)) {
    /* calculate how many times b^2 mod p == 1 */
    mpz_set(t, b);
    m = 0;
    do {
      mpz_powm_ui(t, t, 2, p);
      m++;
    } while (m < r && mpz_cmp_ui(t, 1));
    if (m >= r) break;
    mpz_ui_pow_ui(t, 2, r-m-1);
    mpz_powm(t, z, t, p);
    mpz_mulmod(x, x, t, p, x);
    mpz_powm_ui(z, t, 2, p);
    mpz_mulmod(b, b, z, p, b);
    r = m;
  }
  return _sqrtmod_return(x, a, p, t);
}

/******************************************************************************/

static int _sqrtmod_prime_power(mpz_t r, const mpz_t a, const mpz_t p, int e, mpz_t t, mpz_t u, mpz_t v, mpz_t w) {
  mpz_t n, pk, s;
  int ret, ered;

  if (e == 1) {
    if (mpz_mod(r,a,p), (mpz_cmp_ui(p,2) == 0 || mpz_cmp_ui(r,0) == 0))
      return _sqrtmod_return(r, a, p, t);
    return _sqrtmod_prime(r, a, p,  t,u,v,w);
  }

  mpz_init(n); mpz_init(pk), mpz_init(s);

  mpz_pow_ui(n, p, e);
  mpz_mul(pk, p, p);

  if (mpz_mod(t,a,n), !mpz_cmp_ui(t,0)) {
    mpz_clear(s); mpz_clear(pk); mpz_clear(n);
    mpz_set_ui(r,0);
    return 1;
  }
  if (mpz_mod(t,a,pk), !mpz_cmp_ui(t,0)) {
    mpz_divexact(pk, a, pk);
    ret = _sqrtmod_prime_power(s, pk, p, e-2, t,u,v,w);
    if (ret) mpz_mul(r, s, p);
    mpz_clear(s); mpz_clear(pk); mpz_clear(n);
    return ret;   /* TODO: No verify? */
  }
  if (mpz_mod(t,a,p), !mpz_cmp_ui(t,0)) {
    mpz_clear(s); mpz_clear(pk); mpz_clear(n);
    mpz_set_ui(r,0);
    return 0;
  }

  ered = (mpz_cmp_ui(p,2) > 0 || e < 5)  ?  (e+1)>>1  :  (e+3)>>1;
  if (!_sqrtmod_prime_power(s, a, p, ered, t,u,v,w)) {
    mpz_clear(s); mpz_clear(pk); mpz_clear(n);
    mpz_set_ui(r,0);
    return 0;
  }

  /* This differs from the XS code, which always tries to use n*p */
  /* my $np  = ($p == 2)  ?  Mmulint($n,$p)  :  $n; */
  if (mpz_cmp_ui(p,2)==0) mpz_mul_ui(u,n,2);
  else                    mpz_set(u, n);
  /* my $t1  = Msubmod($a, Mmulmod($s,$s,$np), $np); */
  mpz_mul(v, s, s);
  mpz_sub(v, a, v);
  mpz_mod(v, v, u);
  /* my $t2  = Maddmod($s, $s, $np); */
  mpz_add(w, s, s);
  mpz_mod(w, w, u);
  /* my $gcd = Mgcd($t1, $t2); */
  mpz_gcd(t, v, w);
  /* $r = Maddmod($s, Mdivmod(Mdivint($t1,$gcd),Mdivint($t2,$gcd),$n), $n); */
  mpz_divexact(v, v, t);
  mpz_divexact(w, w, t);
  mpz_divmod(t, v, w, n, u);
  mpz_add(r, s, t);
  mpz_mod(r, r, n);
  /* return ((Mmulmod($r,$r,$n) == ($a % $n)) ? $r : undef); */
  ret = _sqrtmod_return(r, a, n, t);
  mpz_clear(s); mpz_clear(pk); mpz_clear(n);
  return ret;
}


/******************************************************************************/

static int _sqrtmod_composite(mpz_t r, const mpz_t a, const mpz_t n, mpz_t t, mpz_t u, mpz_t v, mpz_t w) {
  mpz_t N, s, fe, *fac;
  int i, nfactors, *exp;

  if (mpz_mod(t,a,n), mpz_perfect_square_p(t)) {
    mpz_sqrt(r, t);
    return _sqrtmod_return(r, a, n, t);
  }

  nfactors = factor(n, &fac, &exp);
  mpz_init_set_ui(N, 1);
  mpz_set_ui(r, 0);
  mpz_init(fe);
  mpz_init(s);

  for (i = 0; i < nfactors; i++) {
    if (!_sqrtmod_prime_power(s, a, fac[i], exp[i], t,u,v,w))
      break;
    mpz_pow_ui(fe, fac[i], exp[i]);
    mpz_sub(t, s, r);
    mpz_mod(t, t, fe);
    if (!mpz_divmod(t, t, N, fe, u))
      break;
    mpz_mul(t, t, N);
    mpz_add(r, r, t);
    mpz_mod(r, r, n);
    mpz_mul(N, N, fe);
  }
  clear_factors(nfactors, &fac, &exp);
  mpz_clear(s); mpz_clear(fe); mpz_clear(N);
  if (i < nfactors) {
    mpz_set_ui(r, 0);
    return 0;
  }
  return _sqrtmod_return(r, a, n, t);
}

/******************************************************************************/

static int sqrtmod_t(mpz_t r, const mpz_t a, const mpz_t n, int isprime,
                     mpz_t t, mpz_t u, mpz_t v, mpz_t w)
{
  if (mpz_cmp_ui(n,2) <= 0) {
    if (mpz_cmp_ui(n,0) <= 0) {
      mpz_set_ui(r,0);
      return 0;
    }
    mpz_mod(r, a, n);
    return _sqrtmod_return(r, a, n, t);
  }
  if (   (mpz_set_ui(r,0), mpz_congruent_p(a, r, n))
      || (mpz_set_ui(r,1), mpz_congruent_p(a, r, n)) )
    return _sqrtmod_return(r, a, n, t);
#if 0
  if (mpz_perfect_square_p(a)) {
    mpz_sqrt(r, a);
    mpz_mod(r, r, n);
    return _sqrtmod_return(r, a, n, t);
  }
#endif
  return (isprime)  ?  _sqrtmod_prime(    r,a,n, t,u,v,w)
                    :  _sqrtmod_composite(r,a,n, t,u,v,w);
}
#define NSMALL 16
static char _small[NSMALL-3+1][NSMALL-2+1] = {
  {0},
  {0,0},
  {0,0,2},
  {0,3,2,0},
  {3,0,2,0,0},
  {0,0,2,0,0,0},
  {0,0,2,0,0,4,0},
  {0,0,2,5,4,0,0,3},
  {0,5,2,4,0,0,0,3,0},
  {0,0,2,0,0,0,0,3,0,0},
  {0,4,2,0,0,0,0,3,6,0,5},
  {4,0,2,0,0,7,6,3,0,5,0,0},
  {0,0,2,0,6,0,0,3,5,0,0,0,0},
  {0,0,2,0,0,0,0,3,0,0,0,0,0,0},
};

/* No temps and r is allowed to alias a */
static int _sqrtmodi(mpz_t r, const mpz_t a, const mpz_t n, int isprime) {
  int res;
  mpz_t x, t1, t2, t3, t4;

  /* Accelerate tiny n as well as a = {0,1} */
  if (mpz_cmp_ui(n,NSMALL) <= 0 || (mpz_sgn(a) >= 0 && mpz_cmp_ui(a,1) <= 0)) {
    unsigned long ua, un = mpz_get_ui(n);
    if (un == 0) { mpz_set_ui(r,0); return 0; }
    ua = mpz_fdiv_ui(a, un);
    if (un > 2 && ua > 1) {
      ua = _small[un-3][ua-2];
      if (ua == 0) { mpz_set_ui(r,0); return 0; }
    }
    mpz_set_ui(r, ua);
    return 1;
  }
  mpz_init(x); mpz_init(t1), mpz_init(t2); mpz_init(t3); mpz_init(t4);
  res = sqrtmod_t(x, a, n, isprime, t1, t2, t3, t4);
  mpz_set(r, x);
  mpz_clear(t4); mpz_clear(t3); mpz_clear(t2); mpz_clear(t1); mpz_clear(x);
  return res;
}

int sqrtmod( mpz_t r, const mpz_t a, const mpz_t n) {return _sqrtmodi(r,a,n,0);}
int sqrtmodp(mpz_t r, const mpz_t a, const mpz_t n) {return _sqrtmodi(r,a,n,1);}

int sqrtmodp_t(mpz_t r, const mpz_t a, const mpz_t p,  mpz_t t1,mpz_t t2,mpz_t t3,mpz_t t4)
{ return sqrtmod_t(r, a, p, 1, t1, t2, t3, t4); }

/******************************************************************************/

#define MAX_ROOTS_RETURNED 600000000

static int _mpz_cmp_roots(const void *av, const void *bv)
  { return mpz_cmp(*(const mpz_t*)av, *(const mpz_t*)bv); }

static int _roots_count_from_mpz(const mpz_t n, UV *count)
{
  if (mpz_sgn(n) < 0 || mpz_cmp_ui(n, MAX_ROOTS_RETURNED) > 0)
    return 0;
  *count = mpz_get_uv(n);
  return 1;
}

static int _roots_count_add_ok(UV a, UV b, UV *sum)
{
  if (a > MAX_ROOTS_RETURNED || b > MAX_ROOTS_RETURNED) return 0;
  if (a > MAX_ROOTS_RETURNED - b) return 0;
  if (sum) *sum = a + b;
  return 1;
}

static int _roots_count_mul_ok(UV a, UV b, UV *prod)
{
  if (a == 0 || b == 0) {
    if (prod) *prod = 0;
    return 1;
  }
  if (a > MAX_ROOTS_RETURNED || b > MAX_ROOTS_RETURNED) return 0;
  if (a > MAX_ROOTS_RETURNED / b) return 0;
  if (prod) *prod = a * b;
  return 1;
}

static void _roots_count_or_croak(const mpz_t n, UV *count)
{
  if (!_roots_count_from_mpz(n, count))
    croak("Maximum returned roots exceeded");
}

static mpz_t* _roots_alloc(UV nroots)
{
  mpz_t *roots;
  UV i;

  if (nroots > MAX_ROOTS_RETURNED)
    croak("Maximum returned roots exceeded");
  if (nroots == 0)
    return 0;
  New(0, roots, nroots, mpz_t);
  for (i = 0; i < nroots; i++)
    mpz_init(roots[i]);
  return roots;
}

void clear_rootmod_list(mpz_t* roots, UV nroots)
{
  UV i;

  if (roots == 0)
    return;
  for (i = 0; i < nroots; i++)
    mpz_clear(roots[i]);
  Safefree(roots);
}

static mpz_t* _one_root(UV *nroots, const mpz_t r)
{
  mpz_t *roots = _roots_alloc(1);

  mpz_set(roots[0], r);
  *nroots = 1;
  return roots;
}

static mpz_t* _two_roots(UV *nroots, const mpz_t r, const mpz_t s)
{
  mpz_t *roots = _roots_alloc(2);

  mpz_set(roots[0], r);
  mpz_set(roots[1], s);
  *nroots = 2;
  return roots;
}

static UV _roots_sort_unique(mpz_t *roots, UV nroots)
{
  UV i, j;

  if (nroots <= 1)
    return nroots;
  qsort(roots, nroots, sizeof(mpz_t), _mpz_cmp_roots);
  for (j = 0, i = 1; i < nroots; i++) {
    if (mpz_cmp(roots[j], roots[i]) != 0) {
      j++;
      if (j != i)
        mpz_swap(roots[j], roots[i]);
    }
  }
  for (i = j+1; i < nroots; i++)
    mpz_clear(roots[i]);
  return j+1;
}

static mpz_t* _rootmod_cprod(UV *nroots,
                             UV nr1, mpz_t *roots1, const mpz_t p1,
                             UV nr2, mpz_t *roots2, const mpz_t p2)
{
  mpz_t *roots;
  mpz_t n, inv, t, u;
  UV nr, i, j;

  if (nr1 == 0 || nr2 == 0) {
    clear_rootmod_list(roots1, nr1);
    clear_rootmod_list(roots2, nr2);
    *nroots = 0;
    return 0;
  }
  if (!_roots_count_mul_ok(nr1, nr2, &nr)) {
    clear_rootmod_list(roots1, nr1);
    clear_rootmod_list(roots2, nr2);
    croak("Maximum returned roots exceeded");
  }

  roots = _roots_alloc(nr);
  mpz_init(n);
  mpz_init(inv);
  mpz_init(t);
  mpz_init(u);

  mpz_mul(n, p1, p2);
  if (!mpz_invert(inv, p1, p2)) {
    clear_rootmod_list(roots1, nr1);
    clear_rootmod_list(roots2, nr2);
    mpz_clear(u); mpz_clear(t); mpz_clear(inv); mpz_clear(n);
    clear_rootmod_list(roots, nr);
    croak("CRT has undefined inverse");
  }

  for (i = 0; i < nr1; i++) {
    for (j = 0; j < nr2; j++) {
      mpz_sub(t, roots2[j], roots1[i]);
      mpz_mod(t, t, p2);
      mpz_mul(t, t, inv);
      mpz_mod(t, t, p2);
      mpz_mul(u, p1, t);
      mpz_add(roots[i*nr2+j], roots1[i], u);
      mpz_mod(roots[i*nr2+j], roots[i*nr2+j], n);
    }
  }

  clear_rootmod_list(roots1, nr1);
  clear_rootmod_list(roots2, nr2);
  mpz_clear(u); mpz_clear(t); mpz_clear(inv); mpz_clear(n);
  *nroots = nr;
  return roots;
}


/******************************************************************************/
/*               SQRTMOD AND ROOTMOD RETURNING ALL RESULTS                    */
/******************************************************************************/

static mpz_t* _allsqrtmodpk(UV *nroots, const mpz_t a, const mpz_t p, int e,
                            mpz_t t, mpz_t u, mpz_t v, mpz_t w)
{
  mpz_t n, pk, low, high, q, q2, pj, idx, a2;
  mpz_t *roots, *roots2;
  UV numr, nr2, pu, i, j;

  *nroots = 0;
  mpz_init(n); mpz_init(pk); mpz_init(low); mpz_init(high);
  mpz_init(q); mpz_init(q2); mpz_init(pj); mpz_init(idx); mpz_init(a2);
  mpz_pow_ui(n, p, (unsigned long)e);

  if (mpz_divisible_p(a, p)) {
    if (mpz_divisible_p(a, n)) {
      mpz_pow_ui(low, p, (unsigned long)(e >> 1));
      if (e & 1) mpz_mul(high, low, p);
      else       mpz_set(high, low);
      _roots_count_or_croak(low, &numr);
      roots = _roots_alloc(numr);
      for (i = 0; i < numr; i++) {
        mpz_set_uv(idx, i);
        mpz_mul(roots[i], high, idx);
      }
      *nroots = numr;
      goto DONE_RETURN;
    }

    mpz_mul(pk, p, p);
    mpz_divexact(a2, a, p);
    if (!mpz_divisible_p(a2, p))
      goto DONE_NONE;

    mpz_divexact(pj, n, p);
    mpz_divexact(a2, a2, p);
    roots2 = _allsqrtmodpk(&nr2, a2, p, e-2, t,u,v,w);
    if (nr2 == 0)
      goto DONE_NONE;
    _roots_count_or_croak(p, &pu);
    if (!_roots_count_mul_ok(nr2, pu, &numr)) {
      clear_rootmod_list(roots2, nr2);
      croak("Maximum returned roots exceeded");
    }
    roots = _roots_alloc(numr);
    for (i = 0; i < nr2; i++) {
      mpz_mul(q, roots2[i], p);
      for (j = 0; j < pu; j++) {
        mpz_set_uv(idx, j);
        mpz_mul(roots[i*pu+j], idx, pj);
        mpz_add(roots[i*pu+j], roots[i*pu+j], q);
      }
    }
    clear_rootmod_list(roots2, nr2);
    *nroots = numr;
    goto DONE_RETURN;
  }

  if (!_sqrtmod_prime_power(q, a, p, e, t,u,v,w))
    goto DONE_NONE;

  roots = _roots_alloc(4);
  mpz_set(roots[0], q);
  mpz_sub(roots[1], n, q);
  if (mpz_cmp_ui(p,2) != 0) {
    *nroots = 2;
  } else if (e == 1) {
    mpz_clear(roots[3]); mpz_clear(roots[2]); mpz_clear(roots[1]);
    *nroots = 1;
  } else if (e == 2) {
    mpz_clear(roots[3]); mpz_clear(roots[2]);
    *nroots = 2;
  } else {
    mpz_divexact_ui(pj, n, 2);
    mpz_sub_ui(t, pj, 1);
    mpz_mul(q2, q, t);
    mpz_mod(q2, q2, n);
    mpz_set(roots[2], q2);
    mpz_sub(roots[3], n, q2);
    *nroots = 4;
  }
  goto DONE_RETURN;

DONE_NONE:
  roots = 0;
DONE_RETURN:
  mpz_clear(a2); mpz_clear(idx); mpz_clear(pj); mpz_clear(q2); mpz_clear(q);
  mpz_clear(high); mpz_clear(low); mpz_clear(pk); mpz_clear(n);
  return roots;
}

static mpz_t* _allsqrtmodfact(UV *nroots, const mpz_t a,
                              mpz_t *fac, int *exp, int nfactors,
                              mpz_t t, mpz_t u, mpz_t v, mpz_t w)
{
  mpz_t N, fe;
  mpz_t *roots, *roots2;
  UV nr, nr2;
  int i;

  *nroots = 0;
  if (nfactors <= 0)
    return 0;

  mpz_init(N);
  mpz_init(fe);
  mpz_pow_ui(N, fac[0], (unsigned long)exp[0]);
  roots = _allsqrtmodpk(&nr, a, fac[0], exp[0], t,u,v,w);
  if (nr == 0) {
    mpz_clear(fe); mpz_clear(N);
    return 0;
  }

  for (i = 1; i < nfactors; i++) {
    mpz_pow_ui(fe, fac[i], (unsigned long)exp[i]);
    roots2 = _allsqrtmodpk(&nr2, a, fac[i], exp[i], t,u,v,w);
    if (nr2 == 0) {
      clear_rootmod_list(roots, nr);
      mpz_clear(fe); mpz_clear(N);
      return 0;
    }
    roots = _rootmod_cprod(&nr, nr, roots, N, nr2, roots2, fe);
    mpz_mul(N, N, fe);
  }

  mpz_clear(fe);
  mpz_clear(N);
  *nroots = nr;
  return roots;
}

mpz_t* allsqrtmod(UV* nroots, const mpz_t a, const mpz_t n)
{
  mpz_t A, t, u, v, w, *fac;
  mpz_t *roots;
  int *exp, nfactors;

  *nroots = 0;
  if (mpz_sgn(n) == 0)
    return 0;

  mpz_init(A);
  mpz_mod(A, a, n);
  if (mpz_cmp_ui(n,2) <= 0) {
    roots = _one_root(nroots, A);
    mpz_clear(A);
    return roots;
  }

  mpz_init(t); mpz_init(u); mpz_init(v); mpz_init(w);
  nfactors = factor(n, &fac, &exp);
  roots = _allsqrtmodfact(nroots, A, fac, exp, nfactors, t,u,v,w);
  clear_factors(nfactors, &fac, &exp);
  if (*nroots > 0)
    *nroots = _roots_sort_unique(roots, *nroots);
  mpz_clear(w); mpz_clear(v); mpz_clear(u); mpz_clear(t);
  mpz_clear(A);
  return roots;
}

UV allsqrtmod_count(const mpz_t a, const mpz_t n)
{
  mpz_t *list;
  UV nroots;

  list = allsqrtmod(&nroots, a, n);
  clear_rootmod_list(list, nroots);
  return nroots;
}


/******************************************************************************/
/*                          K-TH ROOT OF N MOD M                              */
/******************************************************************************/

static int _rootmod_verify(mpz_t r, const mpz_t a, const mpz_t k,
                           const mpz_t n, mpz_t t, mpz_t u)
{
  if (mpz_cmp_ui(k,2) == 0) {
    mpz_sub(t, n, r);
    if (mpz_cmp(t, r) < 0)
      mpz_set(r, t);
  }
  mpz_powm(t, r, k, n);
  mpz_mod(u, a, n);
  if (mpz_cmp(t, u) == 0)
    return 1;
  mpz_set_ui(r, 0);
  return 0;
}

static void _find_ts_generator(mpz_t y, mpz_t m, const mpz_t k, const mpz_t p)
{
  mpz_t p1, r, ke1, x;
  unsigned long e = 0;

  mpz_init(p1); mpz_init(r); mpz_init(ke1); mpz_init_set_ui(x, 2);
  mpz_sub_ui(p1, p, 1);
  mpz_set(r, p1);
  while (mpz_divisible_p(r, k)) {
    mpz_divexact(r, r, k);
    e++;
  }
  if (e == 0)
    croak("bad Tonelli-Shanks input");
  mpz_pow_ui(ke1, k, e-1);

  mpz_set_ui(m, 1);
  while (mpz_cmp_ui(m, 1) == 0) {
    mpz_powm(y, x, r, p);
    if (mpz_cmp_ui(y, 1) != 0)
      mpz_powm(m, y, ke1, p);
    if (mpz_cmp(x, p) >= 0)
      croak("bad Tonelli-Shanks input");
    mpz_add_ui(x, x, 1);
  }

  mpz_clear(x); mpz_clear(ke1); mpz_clear(r); mpz_clear(p1);
}

static int _ts_rootmod(mpz_t x, const mpz_t a, const mpz_t k, const mpz_t p,
                       const mpz_t yin, const mpz_t min)
{
  mpz_t ain, p1, r, inv, A, y, m, T, z, kz, exp, tmp;
  unsigned long e = 0, l;
  int ret = 0;

  mpz_init_set(ain, a);
  mpz_init(p1); mpz_init(r); mpz_init(inv); mpz_init(A);
  mpz_init(y); mpz_init(m); mpz_init(T); mpz_init(z);
  mpz_init(kz); mpz_init(exp); mpz_init(tmp);

  mpz_sub_ui(p1, p, 1);
  mpz_set(r, p1);
  while (mpz_divisible_p(r, k)) {
    mpz_divexact(r, r, k);
    e++;
  }
  if (e == 0)
    goto DONE;

  mpz_mod(tmp, k, r);
  if (!mpz_invert(inv, tmp, r))
    goto DONE;
  mpz_powm(x, ain, inv, p);
  mpz_powm(A, x, k, p);
  if (!mpz_invert(inv, ain, p))
    goto DONE;
  mpz_mul(A, A, inv);
  mpz_mod(A, A, p);

  mpz_set(y, yin);
  mpz_set(m, min);
  if (mpz_cmp_ui(y, 0) == 0 && mpz_cmp_ui(A, 1) != 0)
    _find_ts_generator(y, m, k, p);

  while (mpz_cmp_ui(A, 1) != 0) {
    mpz_set(T, A);
    for (l = 1; mpz_cmp_ui(T, 1) != 0; l++) {
      if (l >= e)
        goto DONE;
      mpz_set(z, T);
      mpz_powm(T, T, k, p);
    }
    if (!_GMP_znlog(kz, z, m, p))
      goto DONE;
    mpz_mod(kz, kz, k);
    if (mpz_sgn(kz) != 0)
      mpz_sub(kz, k, kz);

    mpz_powm(m, m, kz, p);
    mpz_pow_ui(exp, k, e-l);
    mpz_mul(exp, exp, kz);
    mpz_powm(T, y, exp, p);
    e = l-1;
    mpz_mul(x, x, T);
    mpz_mod(x, x, p);
    mpz_powm(y, T, k, p);
    if (mpz_cmp_ui(y, 1) <= 0)
      goto DONE;
    mpz_mul(A, A, y);
    mpz_mod(A, A, p);
  }
  ret = 1;

DONE:
  mpz_clear(tmp); mpz_clear(exp); mpz_clear(kz); mpz_clear(z);
  mpz_clear(T); mpz_clear(m); mpz_clear(y); mpz_clear(A);
  mpz_clear(inv); mpz_clear(r); mpz_clear(p1); mpz_clear(ain);
  return ret;
}

static void _compute_generator(mpz_t y, const mpz_t l, unsigned long e,
                               const mpz_t r, const mpz_t p)
{
  mpz_t x, m, lem1;

  mpz_init_set_ui(x, 2);
  mpz_init_set_ui(m, 1);
  mpz_init(lem1);
  mpz_pow_ui(lem1, l, e-1);
  while (mpz_cmp_ui(m, 1) == 0) {
    mpz_powm(y, x, r, p);
    if (mpz_cmp_ui(y, 1) != 0)
      mpz_powm(m, y, lem1, p);
    mpz_add_ui(x, x, 1);
  }
  mpz_clear(lem1); mpz_clear(m); mpz_clear(x);
}

static int _rootmod_prime_splitk(mpz_t r, mpz_t zeta, const mpz_t a,
                                 const mpz_t k, const mpz_t p)
{
  mpz_t A, p1, g, kg, pg, inv, t, u, y, m, rem, gen, *fac;
  int *exp, nfactors, i, j, ret = 0;

  mpz_init(A); mpz_init(p1); mpz_init(g); mpz_init(kg); mpz_init(pg);
  mpz_init(inv); mpz_init(t); mpz_init(u); mpz_init(y); mpz_init(m);
  mpz_init(rem); mpz_init(gen);

  if (zeta != 0)
    mpz_set_ui(zeta, 1);
  mpz_mod(A, a, p);
  if (mpz_sgn(A) == 0 || (mpz_cmp_ui(A, 1) == 0 && zeta == 0)) {
    mpz_set(r, A);
    ret = 1;
    goto DONE;
  }

  mpz_sub_ui(p1, p, 1);
  if (mpz_cmp_ui(k, 2) == 0) {
    ret = sqrtmodp(r, A, p);
    if (zeta != 0) {
      if (ret) mpz_set(zeta, p1);
      else     mpz_set_ui(zeta, 0);
    }
    goto DONE_VERIFY;
  }

  mpz_gcd(g, k, p1);
  mpz_set(r, A);

  if (mpz_cmp_ui(g, 1) != 0) {
    nfactors = factor(g, &fac, &exp);
    for (i = 0; i < nfactors; i++) {
      if (mpz_sgn(r) == 0)
        break;
      if (zeta != 0) {
        unsigned long V = mpz_remove(rem, p1, fac[i]);
        _compute_generator(gen, fac[i], V, rem, p);
        mpz_pow_ui(t, fac[i], V - (unsigned long)exp[i]);
        mpz_powm(t, gen, t, p);
        mpz_mul(zeta, zeta, t);
        mpz_mod(zeta, zeta, p);
      }
      _find_ts_generator(y, m, fac[i], p);
      for (j = 0; j < exp[i]; j++) {
        if (!_ts_rootmod(r, r, fac[i], p, y, m)) {
          clear_factors(nfactors, &fac, &exp);
          goto DONE;
        }
      }
    }
    clear_factors(nfactors, &fac, &exp);
  }

  if (mpz_cmp(g, k) != 0) {
    mpz_divexact(kg, k, g);
    mpz_divexact(pg, p1, g);
    mpz_mod(t, kg, pg);
    if (!mpz_invert(inv, t, pg))
      goto DONE;
    mpz_powm(r, r, inv, p);
  }
  ret = 1;

DONE_VERIFY:
  if (ret && !_rootmod_verify(r, A, k, p, t, u)) {
    if (zeta != 0)
      mpz_set_ui(zeta, 0);
    ret = 0;
  }

DONE:
  mpz_clear(gen); mpz_clear(rem); mpz_clear(m); mpz_clear(y); mpz_clear(u);
  mpz_clear(t); mpz_clear(inv); mpz_clear(pg); mpz_clear(kg); mpz_clear(g);
  mpz_clear(p1); mpz_clear(A);
  return ret;
}

static mpz_t* _allrootmod_prime(UV *nroots, const mpz_t a,
                                const mpz_t k, const mpz_t p)
{
  mpz_t A, p1, g, exp, inv, r, zeta, r2;
  mpz_t *roots;
  UV numr, i, j;

  *nroots = 0;
  mpz_init(A); mpz_init(p1); mpz_init(g); mpz_init(exp); mpz_init(inv);
  mpz_init(r); mpz_init(zeta); mpz_init(r2);
  mpz_mod(A, a, p);

  if (mpz_cmp_ui(p, 2) == 0 || mpz_sgn(A) == 0) {
    roots = _one_root(nroots, A);
    goto DONE_RETURN;
  }

  mpz_sub_ui(p1, p, 1);
  mpz_gcd(g, k, p1);
  if (mpz_cmp_ui(g, 1) == 0) {
    mpz_mod(exp, k, p1);
    if (!mpz_invert(inv, exp, p1))
      goto DONE_NONE;
    mpz_powm(r, A, inv, p);
    roots = _one_root(nroots, r);
    goto DONE_RETURN;
  }

  mpz_divexact(exp, p1, g);
  mpz_powm(r, A, exp, p);
  if (mpz_cmp_ui(r, 1) != 0)
    goto DONE_NONE;

  if (mpz_cmp_ui(p, 3) == 0) {
    mpz_set_ui(r, 1);
    mpz_set_ui(r2, 2);
    roots = _two_roots(nroots, r, r2);
    goto DONE_RETURN;
  }

  if (!_rootmod_prime_splitk(r, zeta, A, k, p) || mpz_sgn(zeta) == 0)
    croak("allrootmod: failed to find root");
  _roots_count_or_croak(k, &numr);
  roots = _roots_alloc(numr);
  mpz_set(roots[0], r);
  mpz_mul(r2, r, zeta);
  mpz_mod(r2, r2, p);
  for (i = 1; i < numr && mpz_cmp(r2, r) != 0; i++) {
    mpz_set(roots[i], r2);
    mpz_mul(r2, r2, zeta);
    mpz_mod(r2, r2, p);
  }
  if (mpz_cmp(r2, r) != 0) {
    clear_rootmod_list(roots, numr);
    croak("allrootmod: excess roots found");
  }
  for (j = i; j < numr; j++)
    mpz_clear(roots[j]);
  *nroots = i;
  goto DONE_RETURN;

DONE_NONE:
  roots = 0;
DONE_RETURN:
  mpz_clear(r2); mpz_clear(zeta); mpz_clear(r); mpz_clear(inv);
  mpz_clear(exp); mpz_clear(g); mpz_clear(p1); mpz_clear(A);
  return roots;
}

static mpz_t* _allrootmod_prime_power(UV *nroots, const mpz_t a,
                                      const mpz_t k, const mpz_t p, int e)
{
  mpz_t n, pk, step, countz, idx, apk, pe1, pek, s, tt, t1, t2, gcd, r, np, ndivp;
  mpz_t *roots, *roots2;
  UV numr, nr2, pu, pe1u, i, j;
  unsigned long ku = 0, texp;
  int k_le_e;

  *nroots = 0;
  if (e == 1)
    return _allrootmod_prime(nroots, a, k, p);

  mpz_init(n); mpz_init(pk); mpz_init(step); mpz_init(countz); mpz_init(idx);
  mpz_init(apk); mpz_init(pe1); mpz_init(pek); mpz_init(s); mpz_init(tt);
  mpz_init(t1); mpz_init(t2); mpz_init(gcd); mpz_init(r); mpz_init(np);
  mpz_init(ndivp);
  mpz_pow_ui(n, p, (unsigned long)e);
  k_le_e = mpz_fits_ulong_p(k) && (ku = mpz_get_ui(k)) <= (unsigned long)e;

  if (mpz_divisible_p(a, n)) {
    texp = k_le_e ? (((unsigned long)e + ku - 1) / ku) : 1;
    mpz_pow_ui(step, p, texp);
    mpz_pow_ui(countz, p, (unsigned long)e - texp);
    _roots_count_or_croak(countz, &numr);
    roots = _roots_alloc(numr);
    for (i = 0; i < numr; i++) {
      mpz_set_uv(idx, i);
      mpz_mul(roots[i], idx, step);
      mpz_mod(roots[i], roots[i], n);
    }
    *nroots = numr;
    goto DONE_RETURN;
  }

  if (k_le_e) {
    mpz_pow_ui(pk, p, ku);
    if (mpz_divisible_p(a, pk)) {
      mpz_divexact(apk, a, pk);
      mpz_pow_ui(pe1, p, ku-1);
      mpz_pow_ui(pek, p, (unsigned long)e - ku + 1);
      roots2 = _allrootmod_prime_power(&nr2, apk, k, p, e-(int)ku);
      _roots_count_or_croak(pe1, &pe1u);
      if (!_roots_count_mul_ok(nr2, pe1u, &numr)) {
        clear_rootmod_list(roots2, nr2);
        croak("Maximum returned roots exceeded");
      }
      roots = _roots_alloc(numr);
      for (i = 0; i < nr2; i++) {
        mpz_mul(s, roots2[i], p);
        mpz_mod(s, s, n);
        for (j = 0; j < pe1u; j++) {
          mpz_set_uv(idx, j);
          mpz_mul(roots[i*pe1u+j], idx, pek);
          mpz_add(roots[i*pe1u+j], roots[i*pe1u+j], s);
          mpz_mod(roots[i*pe1u+j], roots[i*pe1u+j], n);
        }
      }
      clear_rootmod_list(roots2, nr2);
      *nroots = numr;
      goto DONE_RETURN;
    }
  }

  if (mpz_divisible_p(a, p))
    goto DONE_NONE;

  {
    int ered = (mpz_cmp_ui(p,2) > 0 || e < 5) ? (e+1)>>1 : (e+3)>>1;
    roots2 = _allrootmod_prime_power(&nr2, a, k, p, ered);
  }
  if (nr2 == 0)
    goto DONE_NONE;

  if (mpz_cmp(k, p) != 0) {
    for (j = 0; j < nr2; j++) {
      mpz_set(s, roots2[j]);
      mpz_sub_ui(tt, k, 1);
      mpz_powm(tt, s, tt, n);
      mpz_mul(t1, tt, s);
      mpz_sub(t1, a, t1);
      mpz_mod(t1, t1, n);
      mpz_mul(t2, k, tt);
      mpz_mod(t2, t2, n);
      mpz_gcd(gcd, t1, t2);
      mpz_divexact(t1, t1, gcd);
      mpz_divexact(t2, t2, gcd);
      if (!mpz_divmod(r, t1, t2, n, tt)) {
        clear_rootmod_list(roots2, nr2);
        goto DONE_NONE;
      }
      mpz_add(roots2[j], s, r);
      mpz_mod(roots2[j], roots2[j], n);
    }
    roots = roots2;
    *nroots = nr2;
    goto DONE_RETURN;
  }

  mpz_mul(np, n, p);
  for (j = 0, numr = 0; j < nr2; j++) {
    mpz_set(s, roots2[j]);
    mpz_sub_ui(tt, k, 1);
    mpz_powm(tt, s, tt, np);
    mpz_mul(t1, tt, s);
    mpz_sub(t1, a, t1);
    mpz_mod(t1, t1, np);
    mpz_mul(t2, k, tt);
    mpz_mod(t2, t2, np);
    mpz_gcd(gcd, t1, t2);
    mpz_divexact(t1, t1, gcd);
    mpz_divexact(t2, t2, gcd);
    if (!mpz_divmod(r, t1, t2, n, tt))
      continue;
    mpz_add(r, s, r);
    mpz_mod(r, r, n);
    mpz_powm(tt, r, k, n);
    mpz_mod(t1, a, n);
    if (mpz_cmp(tt, t1) == 0)
      mpz_set(roots2[numr++], r);
  }
  for (j = numr; j < nr2; j++)
    mpz_clear(roots2[j]);
  nr2 = numr;
  if (nr2 == 0) {
    Safefree(roots2);
    goto DONE_NONE;
  }

  _roots_count_or_croak(p, &pu);
  if (!_roots_count_mul_ok(nr2, pu, &numr)) {
    clear_rootmod_list(roots2, nr2);
    croak("Maximum returned roots exceeded");
  }
  roots = _roots_alloc(numr);
  mpz_divexact(ndivp, n, p);
  for (j = 0; j < nr2; j++) {
    for (i = 0; i < pu; i++) {
      mpz_set_uv(idx, i);
      mpz_mul(tt, idx, ndivp);
      mpz_add_ui(tt, tt, 1);
      mpz_mul(roots[j*pu+i], roots2[j], tt);
      mpz_mod(roots[j*pu+i], roots[j*pu+i], n);
    }
  }
  clear_rootmod_list(roots2, nr2);
  *nroots = _roots_sort_unique(roots, numr);
  goto DONE_RETURN;

DONE_NONE:
  roots = 0;
DONE_RETURN:
  mpz_clear(ndivp); mpz_clear(np); mpz_clear(r); mpz_clear(gcd); mpz_clear(t2);
  mpz_clear(t1); mpz_clear(tt); mpz_clear(s); mpz_clear(pek); mpz_clear(pe1);
  mpz_clear(apk); mpz_clear(idx); mpz_clear(countz); mpz_clear(step);
  mpz_clear(pk); mpz_clear(n);
  return roots;
}

static mpz_t* _allrootmod_kprime(UV *nroots, const mpz_t a,
                                 const mpz_t k, const mpz_t n,
                                 mpz_t *fac, int *exp, int nfactors)
{
  mpz_t N, fe, t, u, v, w;
  mpz_t *roots, *roots2;
  UV nr, nr2;
  int i;

  if (mpz_cmp_ui(k, 2) == 0) {
    mpz_init(t); mpz_init(u); mpz_init(v); mpz_init(w);
    roots = _allsqrtmodfact(nroots, a, fac, exp, nfactors, t,u,v,w);
    mpz_clear(w); mpz_clear(v); mpz_clear(u); mpz_clear(t);
    return roots;
  }

  *nroots = 0;
  if (nfactors <= 0)
    return 0;
  mpz_init(N);
  mpz_init(fe);
  mpz_pow_ui(N, fac[0], (unsigned long)exp[0]);
  roots = _allrootmod_prime_power(&nr, a, k, fac[0], exp[0]);
  if (nr == 0) {
    mpz_clear(fe); mpz_clear(N);
    return 0;
  }

  for (i = 1; i < nfactors; i++) {
    mpz_pow_ui(fe, fac[i], (unsigned long)exp[i]);
    roots2 = _allrootmod_prime_power(&nr2, a, k, fac[i], exp[i]);
    if (nr2 == 0) {
      clear_rootmod_list(roots, nr);
      mpz_clear(fe); mpz_clear(N);
      return 0;
    }
    roots = _rootmod_cprod(&nr, nr, roots, N, nr2, roots2, fe);
    mpz_mul(N, N, fe);
  }

  mpz_clear(fe);
  mpz_clear(N);
  *nroots = nr;
  return roots;
}

mpz_t* allrootmod(UV* nroots, const mpz_t a, const mpz_t k, const mpz_t n)
{
  mpz_t A, K, i, *fac, *kfac;
  mpz_t *roots, *roots2, *roots3;
  int *exp, *kexp, nfactors, nkfactors, fi, ei;
  UV numr, nr2, nr3, need, j, t;

  *nroots = 0;
  if (mpz_sgn(n) == 0)
    return 0;

  mpz_init(A);
  mpz_init(K);
  mpz_mod(A, a, n);
  mpz_set(K, k);

  if (mpz_cmp_ui(n, 1) == 0) {
    mpz_set_ui(A, 0);
    roots = _one_root(nroots, A);
    mpz_clear(K); mpz_clear(A);
    return roots;
  }

  if (mpz_sgn(K) < 0) {
    if (mpz_sgn(A) == 0 || !mpz_invert(A, A, n)) {
      mpz_clear(K); mpz_clear(A);
      return 0;
    }
    mpz_abs(K, K);
  }

  if (mpz_sgn(K) == 0) {
    if (mpz_cmp_ui(A, 1) != 0) {
      mpz_clear(K); mpz_clear(A);
      return 0;
    }
    _roots_count_or_croak(n, &numr);
    roots = _roots_alloc(numr);
    mpz_init(i);
    for (j = 0; j < numr; j++) {
      mpz_set_uv(i, j);
      mpz_set(roots[j], i);
    }
    mpz_clear(i);
    *nroots = numr;
    mpz_clear(K); mpz_clear(A);
    return roots;
  }

  if (mpz_cmp_ui(n, 2) == 0 || mpz_cmp_ui(K, 1) == 0) {
    roots = _one_root(nroots, A);
    mpz_clear(K); mpz_clear(A);
    return roots;
  }

  nfactors = factor(n, &fac, &exp);
  if (_GMP_is_prime(K)) {
    roots = _allrootmod_kprime(&numr, A, K, n, fac, exp, nfactors);
  } else {
    roots = _one_root(&numr, A);
    nkfactors = factor(K, &kfac, &kexp);
    for (fi = 0; numr > 0 && fi < nkfactors; fi++) {
      for (ei = 0; numr > 0 && ei < kexp[fi]; ei++) {
        roots3 = 0;
        nr3 = 0;
        for (j = 0; j < numr; j++) {
          roots2 = _allrootmod_kprime(&nr2, roots[j], kfac[fi], n, fac, exp, nfactors);
          if (nr2 == 0)
            continue;
          if (!_roots_count_add_ok(nr3, nr2, &need)) {
            clear_rootmod_list(roots2, nr2);
            clear_rootmod_list(roots3, nr3);
            clear_rootmod_list(roots, numr);
            clear_factors(nkfactors, &kfac, &kexp);
            clear_factors(nfactors, &fac, &exp);
            mpz_clear(K); mpz_clear(A);
            croak("Maximum returned roots exceeded");
          }
          if (roots3 == 0) New(0, roots3, need, mpz_t);
          else             Renew(roots3, need, mpz_t);
          for (t = nr3; t < need; t++)
            mpz_init(roots3[t]);
          for (t = 0; t < nr2; t++)
            mpz_set(roots3[nr3+t], roots2[t]);
          nr3 = need;
          clear_rootmod_list(roots2, nr2);
        }
        clear_rootmod_list(roots, numr);
        roots = roots3;
        numr = nr3;
      }
    }
    clear_factors(nkfactors, &kfac, &kexp);
  }
  clear_factors(nfactors, &fac, &exp);

  if (numr > 0)
    numr = _roots_sort_unique(roots, numr);
  *nroots = numr;
  mpz_clear(K); mpz_clear(A);
  return roots;
}

UV allrootmod_count(const mpz_t a, const mpz_t k, const mpz_t n)
{
  mpz_t *list;
  UV nroots;

  if (mpz_sgn(n) != 0 && mpz_sgn(k) == 0) {
    mpz_t A;
    if (mpz_cmp_ui(n, 1) == 0)
      return 1;
    mpz_init(A);
    mpz_mod(A, a, n);
    if (mpz_cmp_ui(A, 1) == 0)
      _roots_count_or_croak(n, &nroots);
    else
      nroots = 0;
    mpz_clear(A);
    return nroots;
  }

  list = allrootmod(&nroots, a, k, n);
  clear_rootmod_list(list, nroots);
  return nroots;
}

int rootmod(mpz_t r, const mpz_t a, const mpz_t k, const mpz_t n)
{
  mpz_t A, K, t, u;
  mpz_t *roots;
  UV nroots;
  int ret = 0;

  if (mpz_sgn(n) == 0) {
    mpz_set_ui(r, 0);
    return 0;
  }
  if (mpz_cmp_ui(n, 1) == 0) {
    mpz_set_ui(r, 0);
    return 1;
  }

  mpz_init(A); mpz_init(K); mpz_init(t); mpz_init(u);
  mpz_mod(A, a, n);
  mpz_set(K, k);

  if (mpz_sgn(A) == 0) {
    if (mpz_sgn(K) > 0) {
      mpz_set_ui(r, 0);
      ret = 1;
    }
    goto DONE;
  }
  if (mpz_sgn(K) < 0) {
    if (!mpz_invert(A, A, n))
      goto DONE;
    mpz_abs(K, K);
  }
  if (mpz_sgn(K) == 0) {
    if (mpz_cmp_ui(A, 1) == 0) {
      mpz_set_ui(r, 1);
      ret = 1;
    }
    goto DONE;
  }
  if (mpz_cmp_ui(A, 1) == 0) {
    mpz_set_ui(r, 1);
    ret = 1;
    goto DONE;
  }
  if (mpz_cmp_ui(K, 1) == 0) {
    mpz_set(r, A);
    ret = 1;
    goto DONE;
  }
  if (mpz_cmp_ui(K, 2) == 0) {
    ret = sqrtmod(r, A, n);
    goto DONE;
  }
  if (_GMP_is_prime(n)) {
    ret = _rootmod_prime_splitk(r, 0, A, K, n);
    goto DONE;
  }

  roots = allrootmod(&nroots, A, K, n);
  if (nroots > 0) {
    mpz_set(r, roots[0]);
    ret = _rootmod_verify(r, A, K, n, t, u);
  }
  clear_rootmod_list(roots, nroots);

DONE:
  mpz_clear(u); mpz_clear(t); mpz_clear(K); mpz_clear(A);
  if (!ret)
    mpz_set_ui(r, 0);
  return ret;
}
