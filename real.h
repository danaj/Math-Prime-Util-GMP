#ifndef MPU_REAL_H
#define MPU_REAL_H

#include <gmp.h>
#include "ptypes.h"

extern void bernfrac(mpz_t num, mpz_t den, mpz_t n);
extern void harmfrac(mpz_t num, mpz_t den, mpz_t n);
extern void li(mpf_t li, mpf_t x, unsigned long prec);
extern void ei(mpf_t li, mpf_t x, unsigned long prec);

extern void const_euler(mpf_t gamma, unsigned long prec);
extern void const_pi(mpf_t pi, unsigned long prec);
extern void const_log2(mpf_t logn, unsigned long prec);
extern void free_float_constants(void);

extern char* bernreal(mpz_t zn, unsigned long prec);
extern char* harmreal(mpz_t zn, unsigned long prec);
extern char* zetareal(mpf_t r, unsigned long prec);
extern char* eireal(mpf_t r, unsigned long prec);
extern char* lireal(mpf_t r, unsigned long prec);
extern char* riemannrreal(mpf_t r, unsigned long prec);
extern char* lambertwreal(mpf_t r, unsigned long prec);
extern char* logreal(mpf_t r, unsigned long prec);
extern char* expreal(mpf_t r, unsigned long prec);
extern char* powreal(mpf_t r, mpf_t x, unsigned long prec);
extern char* agmreal(mpf_t a, mpf_t b, unsigned long prec);
extern char* eulerconst(unsigned long n);
extern char* piconst(unsigned long n);

#endif
