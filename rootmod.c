#include <gmp.h>
#include "ptypes.h"

#include "rootmod.h"
#include "utility.h"
#include "factor.h"


static int _sqrtmod_return(mpz_t r, mpz_t a, mpz_t n, mpz_t t) {
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
static int _sqrtmod_prime(mpz_t x, mpz_t a, mpz_t p,
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

static int _sqrtmod_prime_power(mpz_t r, mpz_t a, mpz_t p, int e, mpz_t t, mpz_t u, mpz_t v, mpz_t w) {
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

static int _sqrtmod_composite(mpz_t r, mpz_t a, mpz_t n, mpz_t t, mpz_t u, mpz_t v, mpz_t w) {
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

static int sqrtmod_t(mpz_t r, mpz_t a, mpz_t n, int isprime,
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

/* No temps and s is allowed to alias a */
static int _sqrtmodi(mpz_t r, mpz_t a, mpz_t n, int isprime) {
  int res;
  mpz_t x, t1, t2, t3, t4;
  mpz_init(x); mpz_init(t1), mpz_init(t2); mpz_init(t3); mpz_init(t4);
  res = sqrtmod_t(x, a, n, isprime, t1, t2, t3, t4);
  mpz_set(r, x);
  mpz_clear(t4); mpz_clear(t3); mpz_clear(t2); mpz_clear(t1); mpz_clear(x);
  return res;
}

int sqrtmod( mpz_t r, mpz_t a, mpz_t n)  { return _sqrtmodi(r,a,n,0); }
int sqrtmodp(mpz_t r, mpz_t a, mpz_t n)  { return _sqrtmodi(r,a,n,1); }

int sqrtmodp_t(mpz_t r, mpz_t a, mpz_t p,  mpz_t t1,mpz_t t2,mpz_t t3,mpz_t t4)
{ return sqrtmod_t(r, a, p, 1, t1, t2, t3, t4); }

/******************************************************************************/

/* TODO: rootmod, rootmodp */
/* TODO: allsqrtmod */
/* TODO: allrootmod */
