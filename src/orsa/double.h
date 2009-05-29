#ifndef _ORSA_DOUBLE_
#define _ORSA_DOUBLE_

#include <gmpxx.h>

#include <complex>

// while debugging...
#include <orsa/debug.h>
#include <iostream>
#include <string>

namespace orsa {
  
  // standard double 
  
  // typedef double double;
  
  // GMP mpf_class
  
  // typedef mpf_class double;
  // the double constructor should call mpf_set_default_prec()
  // to set the default number of bits (GMP default: 64)
  
  // complex type
  // typedef std::complex<double> Complex;
  
  /* 
     inline const double & 0 {
     static double _zero;
     static unsigned int old_prec = 0;
     if (old_prec != mpf_get_default_prec()) {
     _zero = double("0");
     old_prec = mpf_get_default_prec();
     }	
     return _zero;
     }
  */
  //
  // const double & 0;
  
  /* 
     inline const double & 1 {
     static double _one;
     static unsigned int old_prec = 0;
     if (old_prec != mpf_get_default_prec()) {
     _one = double("1");
     old_prec = mpf_get_default_prec();
     }	
     return _one;
     }
  */
  //
  // const double & 1;
  
  /* 
     inline const double & two() {
     static double _two;
     static unsigned int old_prec = 0;
     if (old_prec != mpf_get_default_prec()) {
     _two = double("2");
     old_prec = mpf_get_default_prec();
     }	
     return _two;
     }
  */
  //
  // const double & two();
  
  // const double & three();
  
  // const double & six();
  
  /* inline double fabs(const double & x) {
     return abs(x);
     }
  */
  
  // computes i! (factorial)
  /* 
     inline mpz_class factorial(const mpz_class & i) {
     if (i < 1) {
     return 1;
     }
     return (i*factorial(i-1));
     }
  */
  //
  mpz_class factorial(const mpz_class & i);
  
  // computes i!! (bi-factorial)
  /* 
     inline mpz_class bi_factorial(const mpz_class & i) {
     if (i < 2) {
     return 1;
     }
     return (i*bi_factorial(i-2));
     }
  */
  // 
  mpz_class bi_factorial(const mpz_class & i);
  
  mpz_class binomial(const mpz_class & n, const mpz_class & k);
  
  // (-1)^l
  /* 
     inline mpz_class power_sign(const mpz_class & l) {
     if ((l%2)==1) {
     return -1.0;
     } else {
     return 1.0;
     }
     }
  */
  //
  int power_sign(const mpz_class & l);
  
  /* 
     inline double int_pow(const double & x, const mpz_class & p) {
     if (p == 0) return 1;
     double _pow = x;
     const mpz_class max_k = abs(p);
     for (mpz_class k=1; k < max_k; ++k) {
     _pow *= x;
     }
     if (p < 0) _pow = 1/_pow;
     return _pow;
     }
  */
  //
  // double int_pow(const double & x, const mpz_class & p);
  // double int_pow(const double & x, const mpz_class & p, const bool can_call_epsilon=true);
  // double int_pow(const double & x, const mpz_class & p);
  double int_pow(const double & x, const int & p);
  
  // square and cube
  inline double square(const double & x) { return (x*x); }
  inline double cube(const double & x) { return (x*x*x); }
  
  //! x^y
  // double pow(const double & x, const double & y);
  
  // double copysign(const double & x, const double & y);
  
  /* 
     inline const double & epsilon() {
     // approx 15 digits (base 10) every 64 bits... is this correct?
     // return int_pow(double("10.0"),mpz_class(-15.0*(mpf_get_default_prec()/64.0))); 
     // ...or a bit better: 16 digits for 64 bits:
     // return int_pow(double("10.0"),mpz_class(-(int)mpf_get_default_prec()/4)); 
     // ...or more... (0.3125 = 1/3.2 -> 20 digits per 64 bits)
     // return int_pow(double("10.0"),mpz_class(-0.3125*mpf_get_default_prec())); 
     // ...or simply 0.3, which seems to be the one working better...
     // return 
     // int_pow(double("10.0"),mpz_class(-0.3*mpf_get_default_prec())); 
     // ... OK, 0.28125 = 18/64 works great.
     static double _eps;
     static unsigned int old_prec = 0;
     if (old_prec != mpf_get_default_prec()) {
     _eps = int_pow(double("10.0"),mpz_class(-0.28125*mpf_get_default_prec())); 
     old_prec = mpf_get_default_prec();
     } 
     return _eps;
     }
  */
  //
  const double & epsilon();
  
  /* 
     inline const double & pi() {
     static double _pi;
     static unsigned int old_prec = 0;
     if (old_prec != mpf_get_default_prec()) {
     const double _eps = epsilon();
     _pi = double("0");
     double _add;
     mpz_class k("0");
     mpz_class eight_k;
     do {
     eight_k = 8*k;
     _add = double("1")/int_pow(mpz_class("16"),k)*
     ( (double("4")/double(eight_k+1)) -
     (double("2")/double(eight_k+4)) -
     (double("1")/double(eight_k+5)) -
     (double("1")/double(eight_k+6)) );
     _pi += _add;
     ++k;
     } while (fabs(_add/_pi) > _eps);
      old_prec = mpf_get_default_prec();
      }
      return _pi;
      }
  */
  //
  const double & pi();
  
  /* 
     inline const double & halfpi() {
     static double _halfpi;
     static unsigned int old_prec = 0;
     if (old_prec != mpf_get_default_prec()) {
     _halfpi = double("0.5")*pi();
     old_prec = mpf_get_default_prec();
     }	
     return _halfpi;
     }
  */
  const double & halfpi();
  
  /* 
     inline const double & twopi() {
     static double _twopi;
     static unsigned int old_prec = 0;
     if (old_prec != mpf_get_default_prec()) {
     _twopi = double("2")*pi();
     old_prec = mpf_get_default_prec();
     }	
     return _twopi;
     }
  */
  //
  const double & twopi();
  
  const double & pisquared();
  
  const double & radToDeg();
  
  const double & degToRad();
  
  const double & radToArcmin();
  
  const double & arcminToRad();
  
  const double & radToArcsec();
  
  const double & arcsecToRad();
  
  /* 
     inline mpz_class kronecker(const mpz_class & i,
     const mpz_class & j) {
     if (i==j) {
     return mpz_class("1");
     } else {
     return mpz_class("0");
     }
     }
  */
  //
  int kronecker(const mpz_class & i,
		const mpz_class & j);
  
  /* 
     inline double fmod(const double & x,
     const double & y) {
     double _ratio = x/y;
     mpz_class _n(_ratio);
     if (_ratio < 0) --_n;
     return (x-_n*y);
     }
  */
  //
  /* double fmod(const double & x,
     const double & y);
  */
  
  // cubic root
  /* 
     inline double cbrt(const double & x) {
     double _c = x;
     double _old_c;
     do {
     _old_c = _c;
     _c = (double("2")*_c+x/(_c*_c))/double("3");
     } while (fabs(_c-_old_c)/fabs(_c) > epsilon());
     return _c;
     }
  */
  //
  // double cbrt(const double & x);
  
  /* 
     inline double pow(const double & x,
     const double & y) {
     ORSA_DEBUG("code needed here!");
     return 0;
     }
  */
  
  // (a)_n = Product[a+k,{k,0,n-1}] with n in N+
  /* 
     inline double pochhammer(const double & a, const mpz_class & n) {
     if (n == 0) {
     return 1;
     }
     double _result = 1;
     mpz_class _k = 0;
     if (n > 0) {
     do { 
     // ORSA_DEBUG("PH+: k: %Zi n: %Zi",_k.get_mpz_t(),n.get_mpz_t());
     _result *= (a+_k);
     ++_k;
     } while (_k != n);
     } else {
     do { 
     // ORSA_DEBUG("PH-: k: %Zi n: %Zi",_k.get_mpz_t(),n.get_mpz_t());
     _result *= (a+_k);
     --_k;
     } while (_k != n);
     }
     return _result;
     }
  */
  //
  double pochhammer(const double    & a, 
		    const mpz_class & n);
  
  /* 
     inline double sin(const double & x) {
     const double _local_x = fmod(x,twopi());
     if (fabs(_local_x) < epsilon()) {
     return 0;
     }
     const double _eps = epsilon();
     double _s("0");
     double _add;
     mpz_class k("0");
     mpz_class twice_k_plus_one;
     do {
     twice_k_plus_one = 2*k+1;
     _add = power_sign(k)*int_pow(_local_x,twice_k_plus_one)/factorial(twice_k_plus_one);
     _s += _add;
     ++k;
     } while (fabs(_add/_s) > _eps);
     return _s;
     }
  */
  //
  // double sin(const double & x);
  
  /* 
     inline double cos(const double & x) {
     const double _local_x = fmod(x,twopi());
     if (fabs(_local_x) < epsilon()) {
     return 1;
     }
     const double _eps = epsilon();
     double _c("0");
     double _add;
     mpz_class k("0");
     mpz_class twice_k;
     do {
     twice_k = 2*k;
     _add = power_sign(k)*int_pow(_local_x,twice_k)/factorial(twice_k);
     _c += _add;
     ++k;
     } while (fabs(_add/_c) > _eps);
     return _c;
     }
  */
  //
  // double cos(const double & x);
  
  // double tan(const double & x);
  
  /* 
     inline void sincos(const double & x, double & s, double & c) {
     s = sin(x);
     // c = sqrt(1.0-s*s); // Warning: sign problem, need value of pi...
     c = cos(x);
     }
  */
  //
  /* void sincos(const double & x, 
     double & s, 
     double & c);
  */
  
  /* 
     inline double asin(const double & x) {
     ORSA_DEBUG("write 3 different cases!!");
     // ORSA_DEBUG("asin(%20.12Fg) called...",x());
     if (fabs(x) < epsilon()) {
     return 0;
     }
     if (fabs(x) > (1+epsilon())) {
     ORSA_ERROR("out-of-domain error: x = %Fg",x());
     return 0;
     }
     const double _eps = epsilon();
     double _as("0");
     double _add;
     const double _one_over_two("0.5");
     mpz_class k("0");
     mpz_class twice_k_plus_one;
     do {
     // ORSA_DEBUG("asin: k: %Zi",k.get_mpz_t());
     twice_k_plus_one = 2*k+1;
     _add = pochhammer(_one_over_two,k)*int_pow(x,twice_k_plus_one)/(twice_k_plus_one*factorial(k));
     _as += _add;
     ++k;
     if (k > 10*mpf_get_default_prec()) {
     ORSA_WARNING("max number of iterations hit while computing asin(%20.12Fg)",x());
     break;
     }
     } while (fabs(_add/_as) > _eps);
     return _as;
     }
  */
  //
  // double asin(const double & x);
  
  /* 
     inline double acos(const double & x) {
     return (halfpi()-asin(x));
     }
  */
  //
  // double acos(const double & x);
  
  /* 
     inline double atan2(const double & y, const double & x) {
     if (fabs(y) < epsilon()) {
     return 0;
     }
     if (fabs(x) < epsilon()) {
     return halfpi();
     }
     
     const double _l = sqrt(x*x+y*y);
     //
     if (fabs(y) < fabs(x)) {
     const double _phi = asin(y/_l);
     //
     if (x > 0) {
     return (_phi);
     } else {
     return (pi()-_phi);
     }
     } else {
     const double _phi = acos(x/_l);
     //
     if (y > 0) {
     return (_phi);
     } else { 
     return (-_phi);
     }
     }
     }
  */
  //
  /* 
     double atan2(const double & y, 
     const double & x);
     
     double atan(const double & x);
     
     double log(const double & x);
  
     double log(const double & base,
     const double & x);
     
     double log10(const double & x);
     
     double exp(const double & x);
     
     double cosh(const double & x);
     double sinh(const double & x);
  */
  
} // namespace orsa

#endif // _ORSA_DOUBLE_
