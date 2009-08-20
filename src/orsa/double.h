#ifndef _ORSA_DOUBLE_
#define _ORSA_DOUBLE_

#include <gmpxx.h>

#include <complex>

// while debugging...
#include <orsa/debug.h>
#include <iostream>
#include <string>

namespace orsa {
  
  // computes i! (factorial)
  mpz_class factorial(const mpz_class & i);
  
  // computes i!! (bi-factorial)
  mpz_class bi_factorial(const mpz_class & i);
  
  mpz_class binomial(const mpz_class & n, const mpz_class & k);
  
  // (-1)^l
  int power_sign(const mpz_class & l);
  
  double int_pow(const double & x, const int & p);
  
  inline double square(const double & x) { return (x*x); }
  
  inline double cube(const double & x) { return (x*x*x); }
  
  const double & epsilon();
  
  const double & pi();
  
  const double & halfpi();
  
  const double & twopi();
  
  const double & pisquared();
  
  const double & radToDeg();
  
  const double & degToRad();
  
  const double & radToArcmin();
  
  const double & arcminToRad();
  
  const double & radToArcsec();
  
  const double & arcsecToRad();
  
  int kronecker(const mpz_class & i,
		const mpz_class & j);
  
  double pochhammer(const double    & a, 
		    const mpz_class & n);
  
} // namespace orsa

#endif // _ORSA_DOUBLE_
