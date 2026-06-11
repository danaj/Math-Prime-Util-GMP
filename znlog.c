#include <gmp.h>

#include "ptypes.h"
#include "znlog.h"
#include "factor.h"
#include "utility.h"

#if BITS_PER_WORD == 32
  #define ZNLOG_BSGS_MAX       65535UL
  #define ZNLOG_BSGS_HARD_MAX  65535UL
#else
  #define ZNLOG_BSGS_MAX       1000000UL
  #define ZNLOG_BSGS_HARD_MAX  12000000UL
#endif
#define ZNLOG_RHO_ATTEMPTS     24
#define ZNLOG_RHO_MIN_ITERS    10000UL
#define ZNLOG_RHO_MAX_ITERS    5000000UL

/* Discrete log via coprime reduction, Pohlig-Hellman, BSGS, and rho.
 * The subproblem solvers are bounded: very hard prime-order logs can fail.
 */

typedef struct znlog_bsgs_entry_t {
  mpz_t value;
  unsigned long exp;
  struct znlog_bsgs_entry_t *next;
} znlog_bsgs_entry_t;

static void _mulmod(mpz_t r, const mpz_t a, const mpz_t b, const mpz_t n, mpz_t t)
{
  mpz_mul(t, a, b);
  mpz_mod(r, t, n);
}

static int _pow_equal(const mpz_t g, const mpz_t k, const mpz_t n, const mpz_t a)
{
  mpz_t t;
  int eq;

  mpz_init(t);
  mpz_powm(t, g, k, n);
  eq = (mpz_cmp(t, a) == 0);
  mpz_clear(t);
  return eq;
}

static unsigned long _rho_max_iters(const mpz_t order)
{
  mpz_t s;
  unsigned long max_iters;

  mpz_init(s);
  mpz_sqrt(s, order);
  if (mpz_fits_ulong_p(s)) {
    unsigned long su = mpz_get_ui(s);
    if (su > (ZNLOG_RHO_MAX_ITERS - 1000UL) / 8UL)
      max_iters = ZNLOG_RHO_MAX_ITERS;
    else
      max_iters = 8UL * su + 1000UL;
  } else {
    max_iters = ZNLOG_RHO_MAX_ITERS;
  }
  if (max_iters < ZNLOG_RHO_MIN_ITERS)
    max_iters = ZNLOG_RHO_MIN_ITERS;
  if (max_iters > ZNLOG_RHO_MAX_ITERS)
    max_iters = ZNLOG_RHO_MAX_ITERS;
  mpz_clear(s);
  return max_iters;
}

static unsigned long _bsgs_steps(const mpz_t order)
{
  mpz_t root;
  unsigned long m;

  mpz_init(root);
  mpz_sqrt(root, order);
  mpz_add_ui(root, root, 1);
  m = mpz_fits_ulong_p(root) ? mpz_get_ui(root) : 0;
  mpz_clear(root);
  return m;
}

static void _bsgs_clear(znlog_bsgs_entry_t **table, unsigned long hsize)
{
  unsigned long i;

  for (i = 0; i < hsize; i++) {
    znlog_bsgs_entry_t *entry = table[i];
    while (entry != 0) {
      znlog_bsgs_entry_t *next = entry->next;
      mpz_clear(entry->value);
      Safefree(entry);
      entry = next;
    }
  }
  Safefree(table);
}

static int _bsgs_lookup(znlog_bsgs_entry_t **table, unsigned long hsize,
                        const mpz_t value, unsigned long *exp)
{
  unsigned long idx;
  znlog_bsgs_entry_t *entry;

  idx = mpz_fdiv_ui(value, hsize);
  entry = table[idx];
  while (entry != 0) {
    if (mpz_cmp(entry->value, value) == 0) {
      *exp = entry->exp;
      return 1;
    }
    entry = entry->next;
  }
  return 0;
}

static void _bsgs_insert(znlog_bsgs_entry_t **table, unsigned long hsize,
                         const mpz_t value, unsigned long exp)
{
  unsigned long idx;
  znlog_bsgs_entry_t *entry;

  idx = mpz_fdiv_ui(value, hsize);
  entry = table[idx];
  while (entry != 0) {
    if (mpz_cmp(entry->value, value) == 0)
      return;
    entry = entry->next;
  }

  New(0, entry, 1, znlog_bsgs_entry_t);
  mpz_init_set(entry->value, value);
  entry->exp = exp;
  entry->next = table[idx];
  table[idx] = entry;
}

static int _znlog_trial(mpz_t r, const mpz_t h, const mpz_t g,
                        const mpz_t n, const mpz_t order,
                        unsigned long limit)
{
  mpz_t x, k;
  unsigned long i;
  int found;

  found = 0;
  mpz_init_set_ui(x, 1);
  mpz_init(k);
  for (i = 0; i <= limit; i++) {
    mpz_set_ui(k, i);
    if (mpz_cmp(k, order) >= 0)
      break;
    if (mpz_cmp(x, h) == 0) {
      mpz_set(r, k);
      found = 1;
      break;
    }
    _mulmod(x, x, g, n, k);
  }
  mpz_clear(k);
  mpz_clear(x);
  return found;
}

static int _znlog_bsgs(mpz_t r, const mpz_t h, const mpz_t g,
                       const mpz_t n, const mpz_t order, unsigned long m)
{
  znlog_bsgs_entry_t **table;
  mpz_t x, y, invg, gm, cand, tmp;
  unsigned long i, j, hsize;
  int found;

  found = 0;
  if (m == 0)
    return 0;

  hsize = 2UL * m + 1UL;
  if (hsize < 1021UL)
    hsize = 1021UL;
  Newz(0, table, hsize, znlog_bsgs_entry_t*);

  mpz_init_set_ui(x, 1);
  mpz_init_set(y, h);
  mpz_init(invg);
  mpz_init(gm);
  mpz_init(cand);
  mpz_init(tmp);

  for (i = 0; i < m; i++) {
    _bsgs_insert(table, hsize, x, i);
    _mulmod(x, x, g, n, tmp);
  }

  if (mpz_invert(invg, g, n)) {
    mpz_powm_ui(gm, invg, m, n);
    for (i = 0; i <= m; i++) {
      if (_bsgs_lookup(table, hsize, y, &j)) {
        mpz_set_ui(cand, i);
        mpz_mul_ui(cand, cand, m);
        mpz_add_ui(cand, cand, j);
        mpz_mod(cand, cand, order);
        if (_pow_equal(g, cand, n, h)) {
          mpz_set(r, cand);
          found = 1;
          break;
        }
      }
      _mulmod(y, y, gm, n, tmp);
    }
  }

  mpz_clear(tmp);
  mpz_clear(cand);
  mpz_clear(gm);
  mpz_clear(invg);
  mpz_clear(y);
  mpz_clear(x);
  _bsgs_clear(table, hsize);

  return found;
}

static void _rho_step(mpz_t x, mpz_t a, mpz_t b,
                      const mpz_t g, const mpz_t h,
                      const mpz_t n, const mpz_t order, mpz_t tmp)
{
  unsigned long bucket;

  bucket = mpz_fdiv_ui(x, 3);
  if (bucket == 0) {
    _mulmod(x, x, g, n, tmp);
    mpz_add_ui(a, a, 1);
    if (mpz_cmp(a, order) >= 0)
      mpz_sub(a, a, order);
  } else if (bucket == 1) {
    _mulmod(x, x, x, n, tmp);
    mpz_mul_ui(a, a, 2);
    mpz_mod(a, a, order);
    mpz_mul_ui(b, b, 2);
    mpz_mod(b, b, order);
  } else {
    _mulmod(x, x, h, n, tmp);
    mpz_add_ui(b, b, 1);
    if (mpz_cmp(b, order) >= 0)
      mpz_sub(b, b, order);
  }
}

static int _znlog_rho(mpz_t r, const mpz_t h, const mpz_t g,
                      const mpz_t n, const mpz_t order)
{
  mpz_t x1, x2, a1, a2, b1, b2, da, db, inv, tmp;
  unsigned long attempt, iter, max_iters;
  int found;

  found = 0;
  max_iters = _rho_max_iters(order);

  mpz_init(x1); mpz_init(x2);
  mpz_init(a1); mpz_init(a2);
  mpz_init(b1); mpz_init(b2);
  mpz_init(da); mpz_init(db);
  mpz_init(inv); mpz_init(tmp);

  for (attempt = 0; attempt < ZNLOG_RHO_ATTEMPTS && !found; attempt++) {
    mpz_isaac_urandomm(a1, order);
    mpz_isaac_urandomm(b1, order);
    if (mpz_sgn(a1) == 0 && mpz_sgn(b1) == 0)
      mpz_add_ui(a1, a1, 1);
    mpz_powm(x1, g, a1, n);
    mpz_powm(tmp, h, b1, n);
    mpz_mul(x1, x1, tmp);
    mpz_mod(x1, x1, n);

    mpz_set(x2, x1);
    mpz_set(a2, a1);
    mpz_set(b2, b1);

    for (iter = 0; iter < max_iters; iter++) {
      _rho_step(x1, a1, b1, g, h, n, order, tmp);
      _rho_step(x2, a2, b2, g, h, n, order, tmp);
      _rho_step(x2, a2, b2, g, h, n, order, tmp);

      if (mpz_cmp(x1, x2) == 0) {
        mpz_sub(da, a1, a2);
        mpz_mod(da, da, order);
        mpz_sub(db, b2, b1);
        mpz_mod(db, db, order);
        if (mpz_sgn(db) != 0 && mpz_invert(inv, db, order)) {
          mpz_mul(r, da, inv);
          mpz_mod(r, r, order);
          if (_pow_equal(g, r, n, h)) {
            found = 1;
            break;
          }
        }
        break;
      }
    }
  }

  mpz_clear(tmp); mpz_clear(inv);
  mpz_clear(db); mpz_clear(da);
  mpz_clear(b2); mpz_clear(b1);
  mpz_clear(a2); mpz_clear(a1);
  mpz_clear(x2); mpz_clear(x1);

  return found;
}

static int _znlog_prime_order(mpz_t r, const mpz_t h, const mpz_t g,
                              const mpz_t n, const mpz_t order)
{
  unsigned long m, VL;

  if (mpz_cmp_ui(h, 1) == 0) {
    mpz_set_ui(r, 0);
    return 1;
  }
  if (mpz_cmp(h, g) == 0) {
    mpz_set_ui(r, 1);
    return 1;
  }

  VL = get_verbose_level();

  if (_znlog_trial(r, h, g, n, order, 200UL))
    return 1;

  m = _bsgs_steps(order);
  if (m <= ZNLOG_BSGS_MAX && _znlog_bsgs(r, h, g, n, order, m)) {
    if (VL) gmp_printf("znlog: BSGS found order %Zd residue %Zd\n",order,r);
    return 1;
  }

  if (_znlog_rho(r, h, g, n, order)) {
    if (VL) gmp_printf("znlog: rho found order %Zd residue %Zd\n",order,r);
    return 1;
  }

  /* BSGS is very memory-heavy here, but it gives a deterministic fallback
   * after rho for large subproblems that still fit our hard cap.
   */
  if (m > ZNLOG_BSGS_MAX && m <= ZNLOG_BSGS_HARD_MAX &&
      _znlog_bsgs(r, h, g, n, order, m)) {
    if (VL) gmp_printf("znlog: large BSGS found order %Zd residue %Zd\n",order,r);
    return 1;
  }

  return 0;
}

static int _znlog_prime_power(mpz_t r, const mpz_t delta, const mpz_t gamma,
                              const mpz_t n, const mpz_t p, int e)
{
  mpz_t pe, pj, exp, gamma0, k, t, h, digit, tmp;
  int j, ok;

  ok = 1;
  mpz_init(pe); mpz_init(pj); mpz_init(exp);
  mpz_init(gamma0); mpz_init(k);
  mpz_init(t); mpz_init(h); mpz_init(digit); mpz_init(tmp);

  mpz_pow_ui(pe, p, (unsigned long)e);
  mpz_pow_ui(exp, p, (unsigned long)(e-1));
  mpz_powm(gamma0, gamma, exp, n);
  mpz_set_ui(k, 0);
  mpz_set_ui(pj, 1);

  for (j = 0; j < e; j++) {
    if (mpz_sgn(k) == 0) {
      mpz_set(t, delta);
    } else {
      mpz_sub(exp, pe, k);
      mpz_powm(t, gamma, exp, n);
      _mulmod(t, t, delta, n, tmp);
    }

    if (j == e-1) {
      mpz_set(h, t);
    } else {
      mpz_pow_ui(exp, p, (unsigned long)(e-1-j));
      mpz_powm(h, t, exp, n);
    }

    if (!_znlog_prime_order(digit, h, gamma0, n, p)) {
      ok = 0;
      break;
    }
    mpz_mul(tmp, digit, pj);
    mpz_add(k, k, tmp);
    if (j != e-1)
      mpz_mul(pj, pj, p);
  }

  if (ok)
    mpz_set(r, k);

  mpz_clear(tmp); mpz_clear(digit); mpz_clear(h); mpz_clear(t);
  mpz_clear(k); mpz_clear(gamma0);
  mpz_clear(exp); mpz_clear(pj); mpz_clear(pe);

  return ok;
}

static int _znlog_coprime(mpz_t r, const mpz_t a, const mpz_t g, const mpz_t n)
{
  mpz_t order, lcm;
  mpz_t *factors, *residues, *moduli;
  int *exponents;
  int nfactors, i, ninit, ok;

  ok = 0;
  nfactors = 0;
  factors = residues = moduli = 0;
  exponents = 0;
  ninit = 0;
  mpz_init(order);

  if (!znorder(order, g, n) || mpz_sgn(order) == 0)
    goto DONE_ORDER;

  {
    mpz_t t;
    mpz_init(t);
    mpz_powm(t, a, order, n);
    if (mpz_cmp_ui(t, 1) != 0) {
      mpz_clear(t);
      goto DONE_ORDER;
    }
    mpz_clear(t);
  }

  if (mpz_cmp_ui(order, 1) == 0) {
    if (mpz_cmp_ui(a, 1) == 0) {
      mpz_set_ui(r, 0);
      ok = 1;
    }
    goto DONE_ORDER;
  }

  nfactors = factor(order, &factors, &exponents);
  if (nfactors <= 0)
    goto DONE_FACTORS;

  New(0, residues, nfactors, mpz_t);
  New(0, moduli, nfactors, mpz_t);
  for (i = 0; i < nfactors; i++) {
    mpz_t P, G, delta, gamma;
    int subok;

    mpz_init(residues[i]);
    mpz_init(moduli[i]);
    ninit++;

    mpz_init(P); mpz_init(G); mpz_init(delta); mpz_init(gamma);
    mpz_pow_ui(P, factors[i], (unsigned long)exponents[i]);
    mpz_set(moduli[i], P);
    mpz_divexact(G, order, P);
    mpz_powm(delta, a, G, n);

    if (mpz_cmp_ui(delta, 1) == 0) {
      mpz_set_ui(residues[i], 0);
      subok = 1;
    } else {
      mpz_powm(gamma, delta, P, n);
      if (mpz_cmp_ui(gamma, 1) != 0) {
        subok = 0;
      } else {
        mpz_powm(gamma, g, G, n);
        subok = _znlog_prime_power(residues[i], delta, gamma, n,
                                   factors[i], exponents[i]);
      }
    }

    mpz_clear(gamma); mpz_clear(delta); mpz_clear(G); mpz_clear(P);
    if (!subok)
      goto DONE_FACTORS;
  }

  mpz_init(lcm);
  if (chinese(r, lcm, residues, moduli, nfactors)) {
    mpz_mod(r, r, order);
    ok = _pow_equal(g, r, n, a);
  }
  mpz_clear(lcm);

DONE_FACTORS:
  for (i = 0; i < ninit; i++) {
    mpz_clear(residues[i]);
    mpz_clear(moduli[i]);
  }
  if (residues != 0) Safefree(residues);
  if (moduli != 0) Safefree(moduli);
  if (factors != 0) clear_factors(nfactors, &factors, &exponents);

DONE_ORDER:
  mpz_clear(order);
  return ok;
}

static int _znlog_general(mpz_t r, const mpz_t a, const mpz_t g, const mpz_t n)
{
  mpz_t d, nr, ar, gp, dacc, tmp, inv, y;
  unsigned long k;
  int ok;

  mpz_init(d);
  mpz_gcd(d, g, n);
  if (mpz_cmp_ui(d, 1) == 0) {
    mpz_clear(d);
    return _znlog_coprime(r, a, g, n);
  }

  ok = 0;
  k = 0;
  mpz_init_set(nr, n);
  mpz_init_set(ar, a);
  mpz_init_set_ui(gp, 1);
  mpz_init_set_ui(dacc, 1);
  mpz_init(tmp); mpz_init(inv); mpz_init(y);

  /* Strip common factors until the remaining base is coprime to the
   * remaining modulus, checking the removed prefix powers as we go.
   */
  while (mpz_cmp_ui(d, 1) > 0) {
    if (mpz_cmp(gp, a) == 0) {
      mpz_set_ui(r, k);
      ok = 1;
      goto DONE;
    }
    if (!mpz_divisible_p(ar, d))
      goto DONE;

    mpz_divexact(nr, nr, d);
    mpz_divexact(ar, ar, d);

    mpz_divexact(tmp, g, d);
    mpz_mul(dacc, dacc, tmp);
    if (mpz_cmp_ui(nr, 1) > 0)
      mpz_mod(dacc, dacc, nr);
    else
      mpz_set_ui(dacc, 0);

    _mulmod(gp, gp, g, n, tmp);
    k++;

    mpz_gcd(d, g, nr);
  }

  if (mpz_cmp(gp, a) == 0) {
    mpz_set_ui(r, k);
    ok = 1;
    goto DONE;
  }
  if (mpz_cmp_ui(nr, 1) == 0)
    goto DONE;
  if (!mpz_invert(inv, dacc, nr))
    goto DONE;

  mpz_mul(ar, ar, inv);
  mpz_mod(ar, ar, nr);
  mpz_mod(tmp, g, nr);

  if (_znlog_general(y, ar, tmp, nr)) {
    mpz_add_ui(r, y, k);
    ok = 1;
  }

DONE:
  mpz_clear(y); mpz_clear(inv); mpz_clear(tmp);
  mpz_clear(dacc); mpz_clear(gp); mpz_clear(ar); mpz_clear(nr);
  mpz_clear(d);
  return ok;
}

int _GMP_znlog(mpz_t r, const mpz_t a, const mpz_t g, const mpz_t n)
{
  mpz_t A, G;
  int ok;

  if (mpz_sgn(n) == 0)
    return 0;
  if (mpz_cmp_ui(n, 1) == 0) {
    mpz_set_ui(r, 0);
    return 1;
  }

  mpz_init(A); mpz_init(G);
  mpz_mod(A, a, n);
  mpz_mod(G, g, n);

  if (mpz_cmp_ui(G, 0) == 0) {
    if (mpz_cmp_ui(A, 1) == 0) {
      mpz_set_ui(r, 0);
      ok = 1;
    } else if (mpz_cmp_ui(A, 0) == 0) {
      mpz_set_ui(r, 1);
      ok = 1;
    } else {
      ok = 0;
    }
  } else if (mpz_cmp_ui(G, 1) == 0) {
    if (mpz_cmp_ui(A, 1) == 0) {
      mpz_set_ui(r, 0);
      ok = 1;
    } else {
      ok = 0;
    }
  } else if (mpz_cmp_ui(A, 1) == 0) {
    mpz_set_ui(r, 0);
    ok = 1;
  } else {
    ok = _znlog_general(r, A, G, n);
    if (ok && !_pow_equal(G, r, n, A))
      ok = 0;
  }

  mpz_clear(G); mpz_clear(A);
  return ok;
}
