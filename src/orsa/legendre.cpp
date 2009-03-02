#include <orsa/legendre.h>

using namespace orsa;

Double Legendre::norm(const mpz_class & l,
		      const mpz_class & m) {
  // return sqrt((__fact__(l-m)*(2.0*l+1.0)*(2.0-__kronecker__(0,m)))/(__fact__(l+m)));
  // return sqrt((factorial(l-m)*(2.0*l+1.0)*(2.0-__kronecker__(0,m)))/(__fact__(l+m)));
  //
  // return sqrt(Double(factorial(l-m)*(2-kronecker(0,m)))/Double(factorial(l+m)));
  //
  // return sqrt(Double(factorial(l-m)*(2-kronecker(0,m))*(2*l-1))/Double(factorial(l+m)));
  //
  // ORSA_DEBUG("pre-sqrt: %Fg",Double(Double(factorial(l+m))/Double(factorial(l-m)*(2-kronecker(0,m))*(2*l+1))).get_mpf_t());
  //
  // almost same as paper
  // return sqrt(Double(factorial(l+m))/Double(factorial(l-m)*(2-kronecker(0,m))*(2*l+1)));
  //
  // fully norm?
  // return sqrt(one()/Double((2*l+1)*(2-kronecker(0,m))));
  // return sqrt(factorial(l-m)/Double(factorial(l+m)*(2*l+1)*(2-kronecker(0,m)))); // this almost, except for a factor 2 for m != 0
  //
  //
  // QUESTO SEMBRA IN OTTIMO ACCORDO CON L'ARTICOLO! (forse manca ancora una sqrt(2) per m != 0)
  // return sqrt(factorial(l-m)/Double(factorial(l+m)*(2*l+1)));
  // prova con sqrt(2) per m != 0
  // return sqrt(Double((2-kronecker(0,m))*factorial(l-m))/Double(factorial(l+m)*(2*l+1)));
  //





  
  // PERFECT!! SAME AS ARTICLE!!
  return sqrt(Double((2-kronecker(0,m))*factorial(l-m))/Double(factorial(l+m)*(2*l+1)));
  //





  //
  // return sqrt(Double(factorial(l+m))/Double(factorial(l-m)*(2*l+1)));
  //
  // another test
  // return sqrt(Double(factorial(l-m)*(2-kronecker(0,m))*(2*l+1))/Double(factorial(l+m)));
  //
  // another test
  // return sqrt(Double(factorial(l+m)*(2-kronecker(0,m))*(2*l+1))/Double(factorial(l-m)));
  // return one();
}
