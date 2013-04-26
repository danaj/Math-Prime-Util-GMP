/*
 * Utility functions, such as sqrt mod p, polynomial manipulation, etc.
 */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <gmp.h>

#ifdef STANDALONE
  typedef unsigned long UV;
  typedef   signed long IV;
  #define INLINE
  #define UV_MAX ULONG_MAX
  #define UVCONST(x) ((unsigned long)x##UL)
  #define croak(fmt,...)            { printf(fmt,##__VA_ARGS__); exit(1); }
  #define New(id, mem, size, type)  mem = (type*) malloc((size)*sizeof(type))
  #define Newz(id, mem, size, type) mem = (type*) calloc(size, sizeof(type))
  #define Safefree(mem)             free((void*)mem)
  #define PRIME_ITERATOR(i) mpz_t i; mpz_init_set_ui(i, 2)
  /*
  static UV prime_iterator_next(mpz_t *iter) { mpz_nextprime(*iter, *iter); return mpz_get_ui(*iter); }
  static void prime_iterator_destroy(mpz_t *iter) { mpz_clear(*iter); }
  static void prime_iterator_setprime(mpz_t *iter, UV n) {mpz_set_ui(*iter, n);}
  static int prime_iterator_isprime(mpz_t *iter, UV n) {int isp; mpz_t t; mpz_init_set_ui(t, n); isp = mpz_probab_prime_p(t, 10); mpz_clear(t); return isp;}
  */
#else
  #include "EXTERN.h"
  #include "perl.h"
  #include "XSUB.h"
#endif

#include "utility.h"


/* set x to sqrt(a) mod p.  Returns 0 if a is not a square root mod p */
/* See Cohen section 1.5.
 * See http://www.math.vt.edu/people/brown/doc/sqrts.pdf
 */
int sqrtmod(mpz_t x, mpz_t a, mpz_t p,
            mpz_t t, mpz_t q, mpz_t b, mpz_t z) /* 4 temp variables */
{
  int r, e, m;

  if (mpz_kronecker(a, p) != 1) { /* No solution exists */
    mpz_set_ui(x, 0);
    return 0;
  }

  /* Easy cases from page 31 */
  if (mpz_congruent_ui_p(p, 3, 4)) {
    mpz_add_ui(t, p, 1);
    mpz_divexact_ui(t, t, 4);
    mpz_powm(x, a, t, p);
    return 1;
  }

  if (mpz_congruent_ui_p(p, 5, 8)) {
    mpz_sub_ui(t, p, 1);
    mpz_divexact_ui(t, t, 4);
    mpz_powm(q, a, t, p);
    if (mpz_cmp_si(q, 1) == 0) { /* s = a^((p+3)/8) mod p */
      mpz_add_ui(t, p, 3);
      mpz_divexact_ui(t, t, 8);
      mpz_powm(x, a, t, p);
    } else {                      /* s = 2a * (4a)^((p-5)/8) mod p */
      mpz_sub_ui(t, p, 5);
      mpz_divexact_ui(t, t, 8);
      mpz_mul_ui(q, a, 4);
      mpz_powm(x, q, t, p);
      mpz_mul_ui(x, x, 2);
      mpz_mulmod(x, x, a, p, x);
    }
    return 1;
  }

  mpz_sub_ui(q, p, 1);
  e = mpz_scan1(q, 0);              /* Remove 2^e from q */
  mpz_tdiv_q_2exp(q, q, e);
  mpz_set_ui(t, 2);
  while (mpz_legendre(t, p) != -1)  /* choose t "at random" */
    mpz_add_ui(t, t, 1);
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
    if (m == r) return 0;
    mpz_ui_pow_ui(t, 2, r-m-1);
    mpz_powm(t, z, t, p);
    mpz_mulmod(x, x, t, p, x);
    mpz_powm_ui(z, t, 2, p);
    mpz_mulmod(b, b, z, p, b);
    r = m;
  }
  return 1;
}

/* Smith-Cornacchia: Solve x,y for x^2 + |D|y^2 = p given prime p */
/* See Cohen 1.5.2 */
int cornacchia(mpz_t x, mpz_t y, mpz_t D, mpz_t p)
{
  int result = 0;
  mpz_t a, b, c, d;

  if (mpz_jacobi(D, p) < 0)     /* No solution */
    return 0;

  mpz_init(a); mpz_init(b); mpz_init(c); mpz_init(d);

  sqrtmod(x, D, p, a, b, c, d);
  mpz_set(a, p);
  mpz_set(b, x);
  mpz_sqrt(c, p);

  while (mpz_cmp(b,c) > 0) {
    mpz_set(d, a);
    mpz_set(a, b);
    mpz_mod(b, d, b);
  }

  mpz_mul(a, b, b);
  mpz_sub(a, p, a);   /* a = p - b^2 */
  mpz_abs(d, D);      /* d = |D| */

  if (mpz_divisible_p(a, d)) {
    mpz_divexact(c, a, d);
    if (mpz_perfect_square_p(c)) {
      mpz_set(x, b);
      mpz_sqrt(y, c);
      result = 1;
    }
  }

  mpz_clear(a); mpz_clear(b); mpz_clear(c); mpz_clear(d);

  return result;
}

/* Modified Cornacchia, Solve x,y for x^2 + |D|y^2 = 4p given prime p */
/* See Cohen 1.5.3 */
int modified_cornacchia(mpz_t x, mpz_t y, mpz_t D, mpz_t p)
{
  int result = 0;
  mpz_t a, b, c, d;

  if (mpz_cmp_ui(p, 2) == 0) {
    mpz_add_ui(x, D, 8);
    if (mpz_perfect_square_p(x)) {
      mpz_sqrt(x, x);
      mpz_set_ui(y, 1);
      result = 1;
    }
    return result;
  }
  if (mpz_jacobi(D, p) == -1)     /* No solution */
    return 0;

  mpz_init(a); mpz_init(b); mpz_init(c); mpz_init(d);

  sqrtmod(x, D, p, a, b, c, d);
  if ( (mpz_even_p(D) && mpz_odd_p(x)) || (mpz_odd_p(D) && mpz_even_p(x)) )
    mpz_sub(x, p, x);

  mpz_mul_ui(a, p, 2);
  mpz_set(b, x);
  mpz_sqrt(c, p);
  mpz_mul_ui(c, c, 2);

  /* Euclidean algorithm */
  while (mpz_cmp(b, c) > 0) {
    mpz_set(d, a);
    mpz_set(a, b);
    mpz_mod(b, d, b);
  }

  mpz_mul_ui(c, p, 4);
  mpz_mul(a, b, b);
  mpz_sub(a, c, a);   /* a = 4p - b^2 */
  mpz_abs(d, D);      /* d = |D| */

  if (mpz_divisible_p(a, d)) {
    mpz_divexact(c, a, d);
    if (mpz_perfect_square_p(c)) {
      mpz_set(x, b);
      mpz_sqrt(y, c);
      result = 1;
    }
  }

  mpz_clear(a); mpz_clear(b); mpz_clear(c); mpz_clear(d);

  return result;
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



UV mpz_order_ui(UV r, mpz_t n, UV limit) {
  UV j;
  mpz_t t;

  /* If n < limit, set limit to n */
  if (mpz_cmp_ui(n, limit) < 0)
    limit = mpz_get_ui(n);
  mpz_init_set_ui(t, 1);
  for (j = 1; j <= limit; j++) {
    mpz_mul(t, t, n);
    mpz_mod_ui(t, t, r);
    if (!mpz_cmp_ui(t, 1))
      break;
  }
  mpz_clear(t);
  return j;
}



#if 0
/* Simple polynomial multiplication */
void poly_mod_mul(mpz_t* px, mpz_t* py, mpz_t* ptmp, UV r, mpz_t mod)
{
  UV i, j, prindex;

  for (i = 0; i < r; i++)
    mpz_set_ui(ptmp[i], 0);
  for (i = 0; i < r; i++) {
    if (!mpz_sgn(px[i])) continue;
    for (j = 0; j < r; j++) {
      if (!mpz_sgn(py[j])) continue;
      prindex = (i+j) % r;
      mpz_addmul( ptmp[prindex], px[i], py[j] );
    }
  }
  /* Put ptmp into px and mod n */
  for (i = 0; i < r; i++)
    mpz_mod(px[i], ptmp[i], mod);
}
void poly_mod_sqr(mpz_t* px, mpz_t* ptmp, UV r, mpz_t mod)
{
  UV i, d, s;
  UV degree = r-1;

  for (i = 0; i < r; i++)
    mpz_set_ui(ptmp[i], 0);
  for (d = 0; d <= 2*degree; d++) {
    UV prindex = d % r;
    for (s = (d <= degree) ? 0 : d-degree; s <= (d/2); s++) {
      if (s*2 == d) {
        mpz_addmul( ptmp[prindex], px[s], px[s] );
      } else {
        mpz_addmul( ptmp[prindex], px[s], px[d-s] );
        mpz_addmul( ptmp[prindex], px[s], px[d-s] );
      }
    }
  }
  /* Put ptmp into px and mod n */
  for (i = 0; i < r; i++)
    mpz_mod(px[i], ptmp[i], mod);
}
#endif

#if 0
/* Binary segmentation, using simple shift+add method for processing p.
 * Faster than twiddling bits, but not nearly as fast as import/export.
 */
void poly_mod_mul(mpz_t* px, mpz_t* py, UV r, mpz_t mod, mpz_t p, mpz_t p2, mpz_t t)
{
  UV i, d, bits;
  UV degree = r-1;

  mpz_mul(t, mod, mod);
  mpz_mul_ui(t, t, r);
  bits = mpz_sizeinbase(t, 2);

  mpz_set_ui(p, 0);
  for (i = 0; i < r; i++) {
    mpz_mul_2exp(p, p, bits);
    mpz_add(p, p, px[r-i-1]);
  }

  if (px == py) {
    mpz_mul(p, p, p);
  } else {
    mpz_set_ui(p2, 0);
    for (i = 0; i < r; i++) {
      mpz_mul_2exp(p2, p2, bits);
      mpz_add(p2, p2, py[r-i-1]);
    }
    mpz_mul(p, p, p2);
  }

  for (d = 0; d <= 2*degree; d++) {
    mpz_tdiv_r_2exp(t, p, bits);
    mpz_tdiv_q_2exp(p, p, bits);
    if (d < r)
      mpz_set(px[d], t);
    else
      mpz_add(px[d-r], px[d-r], t);
  }
  for (i = 0; i < r; i++)
    mpz_mod(px[i], px[i], mod);
}
#endif

/* Binary segmentation, using import/export method for processing p.
 * Thanks to Dan Bernstein's 2007 Quartic paper.
 */
void poly_mod_mul(mpz_t* px, mpz_t* py, UV r, mpz_t mod, mpz_t p, mpz_t p2, mpz_t t)
{
  UV i, bytes;
  char* s;

  mpz_mul(t, mod, mod);
  mpz_mul_ui(t, t, r);
  bytes = mpz_sizeinbase(t, 256);
  mpz_set_ui(p, 0);
  mpz_set_ui(p2, 0);

  /* 1. Create big integers from px and py with padding. */
  {
    Newz(0, s, r*bytes, char);
    for (i = 0; i < r; i++)
      mpz_export(s + i*bytes, NULL, -1, 1, 0, 0, px[i]);
    mpz_import(p, r*bytes, -1, 1, 0, 0, s);
    Safefree(s);
  }
  if (px != py) {
    Newz(0, s, r*bytes, char);
    for (i = 0; i < r; i++)
      mpz_export(s + i*bytes, NULL, -1, 1, 0, 0, py[i]);
    mpz_import(p2, r*bytes, -1, 1, 0, 0, s);
    Safefree(s);
  }

  /* 2. Multiply using the awesomeness that is GMP. */
  mpz_mul( p, p, (px == py) ? p : p2 );

  /* 3. Pull out the parts of the result, add+mod, and put in px. */
  {
    Newz(0, s, 2*r*bytes, char);
    /* fill s with data from p */
    mpz_export(s, NULL, -1, 1, 0, 0, p);
    for (i = 0; i < r; i++) {
      /* Set px[i] to the top part, t to the bottom. */
      mpz_import(px[i], bytes, -1, 1, 0, 0, s + (i+r)*bytes);
      mpz_import(t,     bytes, -1, 1, 0, 0, s +     i*bytes);
      /* Add and mod */
      mpz_add(px[i], px[i], t);
      mpz_mod(px[i], px[i], mod);
    }
    Safefree(s);
  }
}

void poly_mod_pow(mpz_t *pres, mpz_t *pn, mpz_t power, UV r, mpz_t mod)
{
  UV i;
  mpz_t mpow, t1, t2, t3;

  for (i = 0; i < r; i++)
    mpz_set_ui(pres[i], 0);
  mpz_set_ui(pres[0], 1);

  mpz_init_set(mpow, power);
  mpz_init(t1);  mpz_init(t2);  mpz_init(t3);

  while (mpz_cmp_ui(mpow, 0) > 0) {
    if (mpz_odd_p(mpow))            poly_mod_mul(pres, pn, r, mod, t1, t2, t3);
    mpz_tdiv_q_2exp(mpow, mpow, 1);
    if (mpz_cmp_ui(mpow, 0) > 0)    poly_mod_mul(pn, pn, r, mod, t1, t2, t3);
  }
  mpz_clear(t1);  mpz_clear(t2);  mpz_clear(t3);
  mpz_clear(mpow);
}

/* TODO: (1) store these more efficiently, (2) store without variable len */
struct _hilbert_poly {
  signed long  D;
  unsigned long  degree;
  const char *coef[10];
};

static const struct _hilbert_poly _hilbert_data[] =
{
  {  -3, 1, {"0"}},
  {  -4, 1, {"-1728"}},
  {  -7, 1, {"3375"}},
  {  -8, 1, {"-8000"}},
  { -11, 1, {"32768"}},
  { -12, 2, {"-54000","0"}},
  { -15, 2, {"191025","-121287375"}},
  { -16, 2, {"-289224","496793088"}},
  { -19, 1, {"884736"}},
  { -20, 2, {"-1264000","-681472000"}},
  { -23, 3, {"3491750","-5151296875","12771880859375"}},
  { -24, 2, {"-4834944","14670139392"}},
  { -27, 2, {"12288000","0"}},
  { -28, 2, {"-16578000","-55962140625"}},
  { -31, 3, {"39491307","-58682638134","1566028350940383"}},
  { -32, 3, {"-52258000","430167000000","-97336000000000"}},
  { -35, 2, {"117964800","-134217728000"}},
  { -36, 3, {"-153543744","-1525636878336","3094774528868352"}},
  { -39, 4, {"331531596","-429878960946","109873509788637459","20919104368024767633"}},
  { -40, 2, {"-425692800","9103145472000"}},
  { -43, 1, {"884736000"}},
  { -44, 4, {"-1122629840","-36516994456832","8207673077469184","-21405663611349630976"}},
  { -47, 5, {"2257834125","-9987963828125","5115161850595703125","-14982472850828613281250","16042929600623870849609375"}},
  { -48, 4, {"-2835864000","159683258250000","-353673985500000000","0"}},
  { -51, 2, {"5541101568","6262062317568"}},
  { -52, 2, {"-6896880000","-567663552000000"}},
  { -55, 4, {"13136684625","-20948398320757","172576736358790948487","-18577987020116397699793"}},
  { -56, 4, {"-16220384512","2059647197077504","2257767342088912896","10064086044321563803648"}},
  { -59, 3, {"30197678080","-140811576541184","374643194001883136"}},
  { -60, 4, {"-37017885600","-6918204895815375","33749757411875550000","-18577989025032784359375"}},
  { -63, 5, {"67515203250","34794956993579","4557799635952495095663","9127868937627668308696475","-21117051028677601907836561882"}},
  { -64, 4, {"-82226605464","23774457518370936","-38718733083656720832","-3659907775607804768256"}},
  { -67, 1, {"147197952000"}},
  { -68, 4, {"-178211040000","-75843692160000000","-318507038720000000000","-2089297506304000000000000"}},
  { -71, 7, {"313645809715","-3091990136606714","98394038810025431089374","-823534262812956038788256780","5138800360365728878989759011724","-425319434303408151163797366739686","737707041351808575571637718805297513"}},
  { -72, 3, {"-377674776000","235402911936000000","-1859052110336000000000"}},
  { -75, 3, {"654403829760","5209253090426880","0"}},
  { -76, 4, {"-784073554128","-692570204076416256","997755410671625244672","-731887143867097773244416"}},
  { -79, 5, {"1339190283240","-6366718457553016","1793441424178094235755639","-5859423012842750197608704199","5458041018770346219836712950203"}},
  { -80, 6, {"-1597178436000","2005803389486704000","17385261375717088000000","225778539520890715712000000","-417059021735326620672000000000","-287776687393162828546048000000000"}},
  { -83, 3, {"2691907584000","-41490055168000000","549755813888000000000"}},
  { -84, 4, {"-3196800946944","-5663679223085309952","88821246589810089394176","-5133201653210986057826304"}},
  { -87, 6, {"5321761711875","85585228383842544","28321090578679369619067332","497577733930264324271280712821","432181202369154405292219074520147","549806419230953674241390560606228832"}},
  { -88, 2, {"-6294842640000","15798135578688000000"}},
  { -91, 2, {"10359073013760","-3845689020776448"}},
  { -92, 6, {"-12207820358000","-42889702199368156250","-861827826871539375000000","-20685644905048757436279296875","28926541005443767388183593750000","-80048302266801489022434234619140625"}},
  { -95, 8, {"19874477919500","-688170786043879759","395013575867144676284045625","-13089776537013942054366660169906","352163322861023761897441060595178987","-1437415939907503861365582541806826074117","2110631629876441781477723836302339952041537","107789694120266333783482104311412861034259947"}},
  { -96, 6, {"-23340149131680","113518711696500238464","-3583403623173383855142912","7669072309395137266754260992","11327741189564786709352725086208","-14437811688623350992360321646067712"}},
  { -99, 3, {"37616060989440","1176431759374417920","-1840622012131251847168"}},
  { -100, 3, {"-44031499228224","-292067672456279052288","504824415356636531785728"}},
  { -103, 5, {"70292286280125","85475283688254438","4941005649165510742913088361","13355527722149663832413204182160","28826612701416548955034029148929165"}},
  { -104, 6, {"-82028232174464","739545196164376195072","31013571054009020830449664","1378339984770204584193868955648","-25735039642229334200564710375424","65437179730333545242323676123103232"}},
  { -107, 3, {"129783279616000","-6764523159552000000","337618789203968000000000"}},
  { -108, 6, {"-151013216472000","-1847271661190568000000","102945861480026016000000000","-171754585723158768000000000000","1247474246948683776000000000000000","0"}},
  { -111, 8, {"236917342626795","12257744369666687831","56129700127461606440846272378","2987537813843061046262473096622198","-25675269519969449428856645579140830029","88953282142820514438260212803048864368054","-64773995228485734662840257653334874959664148","27524793728233755290201499581457721229607111305"}},
  { -112, 4, {"-274917340548000","4558917032465838750000","-6790363471843861500000000","-74856959786525462843994140625"}},
  { -115, 2, {"427864611225600","130231327260672000"}},
};
#define NUM_HILBERT_POLYS (sizeof(_hilbert_data)/sizeof(_hilbert_data[0]))

UV poly_hilbert(IV D, mpz_t**T)
{
  UV i, j;
  if (T == 0) return 0;
  *T = 0;
  for (i = 0; i < NUM_HILBERT_POLYS; i++) {
    if (_hilbert_data[i].D == D) {
      UV degree = _hilbert_data[i].degree;
      New(0, *T, degree+1, mpz_t);
      mpz_init_set_ui( (*T)[degree], 1 );
      for (j = 1; j <= degree; j++) {
        mpz_init_set_str( (*T)[degree-j], _hilbert_data[i].coef[j-1], 10 );
      }
      return degree;
    }
  }
  return 0;
}
