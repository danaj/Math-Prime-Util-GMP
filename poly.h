#ifndef MPU_POLY_H
#define MPU_POLY_H

#include <math.h>
#include <gmp.h>
#include "ptypes.h"

extern void poly_mod_mul(mpz_t* px, mpz_t* py, UV r, const mpz_t mod, mpz_t t1, mpz_t t2, mpz_t t3);
extern void poly_mod_pow(mpz_t *pres, mpz_t *pn, const mpz_t power, UV r, const mpz_t mod);

extern void poly_mod(mpz_t *pres, mpz_t *pn, UV *dn, const mpz_t mod);


extern void polyz_mod(mpz_t *pres, mpz_t *pn, long *dn, const mpz_t mod);

extern void polyz_set(mpz_t* pr, long* dr, mpz_t* ps, long ds);
extern void polyz_print(const char* header, mpz_t* pn, long dn);
extern void polyz_mulmod(mpz_t* pr, mpz_t* px, mpz_t *py, long *dr, long dx, long dy, const mpz_t mod);
extern void polyz_div(mpz_t *pq, mpz_t *pr, mpz_t *pn, mpz_t *pd,
                      long *dq, long *dr, long dn, long dd, const mpz_t NMOD);
extern void polyz_pow_polymod(mpz_t* pres,  mpz_t* pn,  mpz_t* pmod,
                              long *dres,   long   dn,  long   dmod,
                              const mpz_t power, const mpz_t NMOD);
extern void polyz_gcd(mpz_t* pres, mpz_t* pa, mpz_t* pb, long* dres, long da, long db, const mpz_t MODN);

extern void polyz_root_deg1(mpz_t root, mpz_t* pn, const mpz_t NMOD);
extern void polyz_root_deg2(mpz_t root1, mpz_t root2, mpz_t* pn, const mpz_t NMOD);

/* Find roots of a polynomial mod a prime, slightly modified. */
/* We will stop if we've found at least maxroots unique roots. */
extern void polyz_roots_modp(mpz_t** roots, long *nroots, long maxroots,
                             mpz_t *pP, long dP, const mpz_t NMOD);

extern const char* poly_class_type_name(int type);

/* List of class polynomial indices in order */
extern int* poly_class_nums(void);

/* Given a class poly index, return the degree and fill in (if not null):
 *   D     the discriminant number
 *   T     the polynomial coefficients
 *   type  the poly type:  1 Hilber, 2 Weber
 */
extern UV poly_class_poly_num(int i, int *D, mpz_t**T, int* type);

#endif
