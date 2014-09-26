#ifndef MPU_UTILITY_H
#define MPU_UTILITY_H

#include <gmp.h>
#include "ptypes.h"

extern int get_verbose_level(void);
extern void set_verbose_level(int level);

extern gmp_randstate_t* get_randstate(void);
extern void init_randstate(unsigned long seed);
extern void clear_randstate(void);

/* tdiv_r is faster, but we'd need to guarantee the input is positive */
#define mpz_mulmod(r, a, b, n, t)  \
  do { mpz_mul(t, a, b); mpz_mod(r, t, n); } while (0)

#undef mpz_divmod
extern int mpz_divmod(mpz_t r, mpz_t a, mpz_t b, mpz_t n, mpz_t t);

/* s = sqrt(a) mod p */
extern int sqrtmod(mpz_t s, mpz_t a, mpz_t p,
                   mpz_t t, mpz_t t2, mpz_t b, mpz_t g); /* 4 temp variables */

extern unsigned long modinverse(unsigned long a, unsigned long p);

extern UV mpz_order_ui(UV r, mpz_t n, UV limit);

extern void mpz_arctan(mpz_t r, unsigned long base, mpz_t pow, mpz_t t1, mpz_t t2);
extern void mpz_product(mpz_t* A, UV a, UV b);

extern void poly_mod_mul(mpz_t* px, mpz_t* py, UV r, mpz_t mod, mpz_t t1, mpz_t t2, mpz_t t3);
extern void poly_mod_pow(mpz_t *pres, mpz_t *pn, mpz_t power, UV r, mpz_t mod);

extern void poly_mod(mpz_t *pres, mpz_t *pn, UV *dn, mpz_t mod);
extern void polyz_mod(mpz_t *pres, mpz_t *pn, long *dn, mpz_t mod);

extern void polyz_set(mpz_t* pr, long* dr, mpz_t* ps, long ds);
extern void polyz_print(const char* header, mpz_t* pn, long dn);
extern void polyz_mulmod(mpz_t* pr, mpz_t* px, mpz_t *py, long *dr, long dx, long dy, mpz_t mod);
extern void polyz_div(mpz_t *pq, mpz_t *pr, mpz_t *pn, mpz_t *pd,
                      long *dq, long *dr, long dn, long dd, mpz_t NMOD);
extern void polyz_pow_polymod(mpz_t* pres,  mpz_t* pn,  mpz_t* pmod,
                              long *dres,   long   dn,  long   dmod,
                              mpz_t power, mpz_t NMOD);
extern void polyz_gcd(mpz_t* pres, mpz_t* pa, mpz_t* pb, long* dres, long da, long db, mpz_t MODN);

extern void polyz_root_deg1(mpz_t root, mpz_t* pn, mpz_t NMOD);
extern void polyz_root_deg2(mpz_t root1, mpz_t root2, mpz_t* pn, mpz_t NMOD);

/* Find roots of a polynomial mod a prime, slightly modified. */
/* We will stop if we've found at least maxroots unique roots. */
extern void polyz_roots_modp(mpz_t** roots, long *nroots, long maxroots,
                             mpz_t *pP, long dP, mpz_t NMOD,
                             gmp_randstate_t* p_randstate);

/* Solve x^2 + |D|y^2 = p */
extern int cornacchia(mpz_t x, mpz_t y, mpz_t D, mpz_t p);
/* Solve x^2 + |D|y^2 = 4p */
extern int modified_cornacchia(mpz_t x, mpz_t y, mpz_t D, mpz_t p);

/* return a class poly (Hilbert [type 1] or Weber [type 2]) */
extern UV poly_class_poly(IV D, mpz_t**T, int* type);

/* return a 0 terminated list of all D's sorted by degree */
extern IV* poly_class_degrees(int insert_1s);

/* List of class polynomial indices in order */
extern int* poly_class_nums(void);
/* Given a class poly index, return the degree and fill in (if not null):
 *   D     the discriminant number
 *   T     the polynomial coefficients
 *   type  the poly type:  1 Hilber, 2 Weber
 */
extern UV poly_class_poly_num(int i, int *D, mpz_t**T, int* type);

#endif
