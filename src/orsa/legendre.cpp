#include <orsa/legendre.h>

using namespace orsa;

double Legendre::norm(const mpz_class & l,
		      const mpz_class & m) {
  // return sqrt((__fact__(l-m)*(2.0*l+1.0)*(2.0-__kronecker__(0,m)))/(__fact__(l+m)));
  // return sqrt((factorial(l-m)*(2.0*l+1.0)*(2.0-__kronecker__(0,m)))/(__fact__(l+m)));
  //
  // return sqrt(double(factorial(l-m)*(2-kronecker(0,m)))/double(factorial(l+m)));
  //
  // return sqrt(double(factorial(l-m)*(2-kronecker(0,m))*(2*l-1))/double(factorial(l+m)));
  //
  // ORSA_DEBUG("pre-sqrt: %Fg",double(double(factorial(l+m))/double(factorial(l-m)*(2-kronecker(0,m))*(2*l+1)))());
  //
  // almost same as paper
  // return sqrt(double(factorial(l+m))/double(factorial(l-m)*(2-kronecker(0,m))*(2*l+1)));
  //
  // fully norm?
  // return sqrt(1/double((2*l+1)*(2-kronecker(0,m))));
  // return sqrt(factorial(l-m)/double(factorial(l+m)*(2*l+1)*(2-kronecker(0,m)))); // this almost, except for a factor 2 for m != 0
  //
  //
  // QUESTO SEMBRA IN OTTIMO ACCORDO CON L'ARTICOLO! (forse manca ancora una sqrt(2) per m != 0)
  // return sqrt(factorial(l-m)/double(factorial(l+m)*(2*l+1)));
  // prova con sqrt(2) per m != 0
  // return sqrt(double((2-kronecker(0,m))*factorial(l-m))/double(factorial(l+m)*(2*l+1)));
  //





  
  // PERFECT!! SAME AS ARTICLE!!
  return sqrt(mpz_class((2-kronecker(0,m))*factorial(l-m)).get_d()/mpz_class(factorial(l+m)*(2*l+1)).get_d());
  //
  
  
  
  

  //
  // return sqrt(double(factorial(l+m))/double(factorial(l-m)*(2*l+1)));
  //
  // another test
  // return sqrt(double(factorial(l-m)*(2-kronecker(0,m))*(2*l+1))/double(factorial(l+m)));
  //
  // another test
  // return sqrt(double(factorial(l+m)*(2-kronecker(0,m))*(2*l+1))/double(factorial(l-m)));
  // return 1;
}
