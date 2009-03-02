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
  
  // typedef double Double;
  
  // GMP mpf_class
  
  typedef mpf_class Double;
  // the Double constructor should call mpf_set_default_prec()
  // to set the default number of bits (GMP default: 64)
  
  // complex type
  typedef std::complex<orsa::Double> Complex;
  
  /* 
     inline const Double & zero() {
     static Double _zero;
     static unsigned int old_prec = 0;
     if (old_prec != mpf_get_default_prec()) {
     _zero = Double("0");
     old_prec = mpf_get_default_prec();
     }	
     return _zero;
     }
  */
  //
  const Double & zero();
  
  /* 
     inline const Double & one() {
     static Double _one;
     static unsigned int old_prec = 0;
     if (old_prec != mpf_get_default_prec()) {
     _one = Double("1");
     old_prec = mpf_get_default_prec();
     }	
     return _one;
     }
  */
  //
  const Double & one();
  
  /* 
     inline const Double & two() {
     static Double _two;
     static unsigned int old_prec = 0;
     if (old_prec != mpf_get_default_prec()) {
     _two = Double("2");
     old_prec = mpf_get_default_prec();
     }	
     return _two;
     }
  */
  //
  const Double & two();
  
  const Double & three();
  
  const Double & six();
  
  inline Double fabs(const Double & x) {
    return abs(x);
  }
  
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
  mpz_class power_sign(const mpz_class & l);
  
  /* 
     inline Double int_pow(const Double & x, const mpz_class & p) {
     if (p == 0) return one();
     Double _pow = x;
     const mpz_class max_k = abs(p);
     for (mpz_class k=1; k < max_k; ++k) {
     _pow *= x;
     }
     if (p < 0) _pow = one()/_pow;
     return _pow;
     }
  */
  //
  Double int_pow(const Double & x, const mpz_class & p);
  // Double int_pow(const Double & x, const mpz_class & p, const bool can_call_epsilon=true);
  
  //! x^y
  Double pow(const Double & x, const Double & y);
  
  Double copysign(const Double & x, const Double & y);
  
  /* 
     inline const Double & epsilon() {
     // approx 15 digits (base 10) every 64 bits... is this correct?
     // return int_pow(Double("10.0"),mpz_class(-15.0*(mpf_get_default_prec()/64.0))); 
     // ...or a bit better: 16 digits for 64 bits:
     // return int_pow(Double("10.0"),mpz_class(-(int)mpf_get_default_prec()/4)); 
     // ...or more... (0.3125 = 1/3.2 -> 20 digits per 64 bits)
     // return int_pow(Double("10.0"),mpz_class(-0.3125*mpf_get_default_prec())); 
     // ...or simply 0.3, which seems to be the one working better...
     // return 
     // int_pow(Double("10.0"),mpz_class(-0.3*mpf_get_default_prec())); 
     // ... OK, 0.28125 = 18/64 works great.
     static Double _eps;
     static unsigned int old_prec = 0;
     if (old_prec != mpf_get_default_prec()) {
     _eps = int_pow(Double("10.0"),mpz_class(-0.28125*mpf_get_default_prec())); 
     old_prec = mpf_get_default_prec();
     } 
     return _eps;
     }
  */
  //
  const Double & epsilon();
  
  /* 
     inline const Double & pi() {
     static Double _pi;
     static unsigned int old_prec = 0;
     if (old_prec != mpf_get_default_prec()) {
     const Double _eps = epsilon();
     _pi = Double("0");
     Double _add;
     mpz_class k("0");
     mpz_class eight_k;
     do {
     eight_k = 8*k;
     _add = Double("1")/int_pow(mpz_class("16"),k)*
     ( (Double("4")/Double(eight_k+1)) -
     (Double("2")/Double(eight_k+4)) -
     (Double("1")/Double(eight_k+5)) -
     (Double("1")/Double(eight_k+6)) );
     _pi += _add;
     ++k;
     } while (fabs(_add/_pi) > _eps);
      old_prec = mpf_get_default_prec();
      }
      return _pi;
      }
  */
  //
  const Double & pi();
  
  /* 
     inline const Double & halfpi() {
     static Double _halfpi;
     static unsigned int old_prec = 0;
     if (old_prec != mpf_get_default_prec()) {
     _halfpi = Double("0.5")*pi();
     old_prec = mpf_get_default_prec();
     }	
     return _halfpi;
     }
  */
  const Double & halfpi();
  
  /* 
     inline const Double & twopi() {
     static Double _twopi;
     static unsigned int old_prec = 0;
     if (old_prec != mpf_get_default_prec()) {
     _twopi = Double("2")*pi();
     old_prec = mpf_get_default_prec();
     }	
     return _twopi;
     }
  */
  //
  const Double & twopi();
  
  const Double & pisquared();
  
  const Double & radToDeg();
  
  const Double & degToRad();
  
  const Double & radToArcmin();
  
  const Double & arcminToRad();
  
  const Double & radToArcsec();
  
  const Double & arcsecToRad();
  
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
  mpz_class kronecker(const mpz_class & i,
		      const mpz_class & j);
  
  /* 
     inline Double fmod(const Double & x,
     const Double & y) {
     Double _ratio = x/y;
     mpz_class _n(_ratio);
     if (_ratio < 0) --_n;
     return (x-_n*y);
     }
  */
  //
  Double fmod(const Double & x,
	      const Double & y);
  
  // cubic root
  /* 
     inline Double cbrt(const Double & x) {
     Double _c = x;
     Double _old_c;
     do {
     _old_c = _c;
     _c = (Double("2")*_c+x/(_c*_c))/Double("3");
     } while (fabs(_c-_old_c)/fabs(_c) > epsilon());
     return _c;
     }
  */
  //
  Double cbrt(const Double & x);
  
  /* 
     inline Double pow(const Double & x,
     const Double & y) {
     ORSA_DEBUG("code needed here!");
     return zero();
     }
  */
  
  // (a)_n = Product[a+k,{k,0,n-1}] with n in N+
  /* 
     inline Double pochhammer(const Double & a, const mpz_class & n) {
     if (n == 0) {
     return one();
     }
     Double _result = one();
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
  Double pochhammer(const Double & a, 
		    const mpz_class & n);
  
  /* 
     inline Double sin(const Double & x) {
     const Double _local_x = fmod(x,twopi());
     if (fabs(_local_x) < epsilon()) {
     return zero();
     }
     const Double _eps = epsilon();
     Double _s("0");
     Double _add;
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
  Double sin(const Double & x);
  
  /* 
     inline Double cos(const Double & x) {
     const Double _local_x = fmod(x,twopi());
     if (fabs(_local_x) < epsilon()) {
     return one();
     }
     const Double _eps = epsilon();
     Double _c("0");
     Double _add;
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
  Double cos(const Double & x);
  
  Double tan(const Double & x);
  
  /* 
     inline void sincos(const Double & x, Double & s, Double & c) {
     s = sin(x);
     // c = sqrt(1.0-s*s); // Warning: sign problem, need value of pi...
     c = cos(x);
     }
  */
  //
  void sincos(const Double & x, 
	      Double & s, 
	      Double & c);
  
  /* 
     inline Double asin(const Double & x) {
     ORSA_DEBUG("write 3 different cases!!");
     // ORSA_DEBUG("asin(%20.12Fg) called...",x.get_mpf_t());
     if (fabs(x) < epsilon()) {
     return zero();
     }
     if (fabs(x) > (one()+epsilon())) {
     ORSA_ERROR("out-of-domain error: x = %Fg",x.get_mpf_t());
     return zero();
     }
     const Double _eps = epsilon();
     Double _as("0");
     Double _add;
     const Double _one_over_two("0.5");
     mpz_class k("0");
     mpz_class twice_k_plus_one;
     do {
     // ORSA_DEBUG("asin: k: %Zi",k.get_mpz_t());
     twice_k_plus_one = 2*k+1;
     _add = pochhammer(_one_over_two,k)*int_pow(x,twice_k_plus_one)/(twice_k_plus_one*factorial(k));
     _as += _add;
     ++k;
     if (k > 10*mpf_get_default_prec()) {
     ORSA_WARNING("max number of iterations hit while computing asin(%20.12Fg)",x.get_mpf_t());
     break;
     }
     } while (fabs(_add/_as) > _eps);
     return _as;
     }
  */
  //
  Double asin(const Double & x);
  
  /* 
     inline Double acos(const Double & x) {
     return (halfpi()-asin(x));
     }
  */
  //
  Double acos(const Double & x);
  
  /* 
     inline Double atan2(const Double & y, const Double & x) {
     if (fabs(y) < epsilon()) {
     return zero();
     }
     if (fabs(x) < epsilon()) {
     return halfpi();
     }
     
     const Double _l = sqrt(x*x+y*y);
     //
     if (fabs(y) < fabs(x)) {
     const Double _phi = asin(y/_l);
     //
     if (x > zero()) {
     return (_phi);
     } else {
     return (pi()-_phi);
     }
     } else {
     const Double _phi = acos(x/_l);
     //
     if (y > zero()) {
     return (_phi);
     } else { 
     return (-_phi);
     }
     }
     }
  */
  //
  Double atan2(const Double & y, 
	       const Double & x);
  
  Double atan(const Double & x);
  
  Double log(const Double & x);
  
  Double log(const Double & base,
	     const Double & x);
  
  Double log10(const Double & x);
  
  Double exp(const Double & x);
  
  Double cosh(const Double & x);
  Double sinh(const Double & x);
  
} // namespace orsa

#endif // _ORSA_DOUBLE_
