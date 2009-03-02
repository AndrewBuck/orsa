#include <orsa/double.h>
#include <orsa/cache.h>

#include <vector>

using namespace orsa;

const Double & orsa::zero() {
  static Double _zero;
  static unsigned int old_prec = 0;
  if (old_prec != mpf_get_default_prec()) {
    _zero = Double("0");
    old_prec = mpf_get_default_prec();
  }	
  return _zero;
}

const Double & orsa::one() {
  static Double _one;
  static unsigned int old_prec = 0;
  if (old_prec != mpf_get_default_prec()) {
    _one = Double("1");
    old_prec = mpf_get_default_prec();
  }	
  return _one;
}

const Double & orsa::two() {
  static Double _two;
  static unsigned int old_prec = 0;
  if (old_prec != mpf_get_default_prec()) {
    _two = Double("2");
    old_prec = mpf_get_default_prec();
  }	
  return _two;
}

const Double & orsa::three() {
  static Double _three;
  static unsigned int old_prec = 0;
  if (old_prec != mpf_get_default_prec()) {
    _three = Double("3");
    old_prec = mpf_get_default_prec();
  }	
  return _three;
}

const Double & orsa::six() {
  static Double _six;
  static unsigned int old_prec = 0;
  if (old_prec != mpf_get_default_prec()) {
    _six = Double("6");
    old_prec = mpf_get_default_prec();
  }	
  return _six;
}

/* 
   mpz_class orsa::factorial(const mpz_class & i) {
   if (i <= 1) {
   return 1;
   }
   return (i*factorial(i-1));
   }
*/

mpz_class orsa::factorial(const mpz_class & i) {
  // ORSA_DEBUG("f: %Zi",i.get_mpz_t());
  static std::vector< orsa::Cache<mpz_class> > _factorial_table;
  const unsigned long int index = i.get_ui();
  static mpz_class _mpz_one("1");
  if (i <= _mpz_one) {
    return _mpz_one;
  } else if (_factorial_table.size() > i) {
    if (!_factorial_table[index].isSet()) {
      _factorial_table[index].set(i*factorial(i-_mpz_one));
    }
    return _factorial_table[index].get();
  } else {
    _factorial_table.resize(index+1);
    _factorial_table[index].set(i*factorial(i-_mpz_one));
    return _factorial_table[index].get();
  }
}

/* 
   mpz_class orsa::bi_factorial(const mpz_class & i) {
   ORSA_DEBUG("bf: %Zi",i.get_mpz_t());
   if (i <= mpz_class("1")) {
   return mpz_class("1");
   }
   return (i*bi_factorial(i-mpz_class("2")));
   }
*/

mpz_class orsa::bi_factorial(const mpz_class & i) {
  // ORSA_DEBUG("bf: %Zi",i.get_mpz_t());
  static std::vector< orsa::Cache<mpz_class> > _bi_factorial_table;
  const unsigned long int index = i.get_ui();
  static mpz_class _mpz_one("1");
  static mpz_class _mpz_two("2");
  if (i <= _mpz_one) {
    return _mpz_one;
  } else if (_bi_factorial_table.size() > i) {
    if (!_bi_factorial_table[index].isSet()) {
      _bi_factorial_table[index].set(i*bi_factorial(i-_mpz_two));
    }
    return _bi_factorial_table[index].get();
  } else {
    _bi_factorial_table.resize(index+1);
    _bi_factorial_table[index].set(i*bi_factorial(i-_mpz_two));
    return _bi_factorial_table[index].get();
  }
}

mpz_class orsa::binomial(const mpz_class & n, const mpz_class & k) {
  const mpz_class retVal = ( (factorial(n)) / 
			     (factorial(k)*factorial(n-k)) );

  /* 
     ORSA_DEBUG("binomial(%Zi,%Zi) = %Zi",
     n.get_mpz_t(),
     k.get_mpz_t(),
     retVal.get_mpz_t());
  */
  
  return retVal;
}

mpz_class orsa::power_sign(const mpz_class & l) {
  if ((l%2)==1) {
    return -1.0;
  } else {
    return 1.0;
  }
}

// call this from int_pow, to avoid infinite loops...
/* 
   static Double __rough_epsilon__() {
   // keep this number in sync with epsilon()
   // positive!
   // const mpz_class (0.28125*mpf_get_default_prec())
   static Double _r_eps;
   static unsigned int old_prec = 0;
   if (old_prec != mpf_get_default_prec()) {
   _r_eps = one();
   const Double _one_over_ten("0.1");
   const mpz_class _pow = abs(mpz_class(0.28125*mpf_get_default_prec()));
   for (unsigned int p=0; p<_pow; ++p) {
   _r_eps *= _one_over_ten;
   }
   old_prec = mpf_get_default_prec();
   } 
   return _r_eps;
   }
*/

Double orsa::int_pow(const Double & x, 
		     const mpz_class & p) {
  // ORSA_DEBUG("int_pow(%Ff,%Zi)",x.get_mpf_t(),p.get_mpz_t());
  if (p ==  2) return x*x;
  if (p ==  1) return x;
  if (p ==  0) return one();
  if (p == -1) return one()/x;
  if (fabs(x) < epsilon()) {
    return zero();
  }
  Double _pow = x;
  const mpz_class max_k = abs(p);
  for (mpz_class k=1; k < max_k; ++k) {
    _pow *= x;
  }
  if (p < 0) _pow = one()/_pow;
  return _pow;
}

Double orsa::pow(const Double & x, const Double & y) {
  if ((fabs(x) < epsilon()) && (y > zero())) {
    return zero();
  }
  if (x < zero()) {
    ORSA_ERROR("called pow(%Fg,%Fg) with non-positive argument; if the power is integer, then call int_pow(...); returning zero.",
	       x.get_mpf_t(),
	       y.get_mpf_t());
    return zero();
  }
  //
  if (fabs(x-one()) < epsilon()) {
    return one();
  }
  //
  if (fabs(y) < epsilon()) {
    return one();
  }
  //
  mpz_class _q("1");
  Double _val(fabs(x));
  // domain limits are 0->2, but convergence is faster closer to 1.0
  while (_val > Double("1.1")) {
    _val = sqrt(_val);
    _q *= mpz_class("2");
  }
  //
  while (_val < Double("0.9")) {
    _val = sqrt(_val);
    _q *= mpz_class("2");
  }
  //
  const Double _val_minus_one = _val - one();
  // const Double _val_minus_one = copysign(_val,x) - one();
  //
  const Double _eps = epsilon();
  Double _pow("0");
  Double _add;
  mpz_class k("0");
  do {
    _add = power_sign(k)*pochhammer(-(_q*y),k)*int_pow(_val_minus_one,k)/factorial(k);
    _pow += _add;
    //
    /* 
       ORSA_DEBUG("k: %Zi   _pow: %Fg   _add: %Fg",
       k.get_mpz_t(),
       _pow.get_mpf_t(),
       _add.get_mpf_t());
    */
    //
    ++k;
  } while (fabs(_add/_pow) > _eps);
  /* 
     ORSA_DEBUG("pow(%Ff,%Ff) =  %Ff",
     x.get_mpf_t(),
     y.get_mpf_t(),
     _pow.get_mpf_t());
  */
  return (_pow);
}

Double orsa::copysign(const Double & x, 
		      const Double & y) {
  Double _result = fabs(x);
  if (y < zero()) _result *= -one();
  return _result;
}

/* 
   const Double & orsa::epsilon() {
   // ORSA_DEBUG("epsilon power: %Zi",mpz_class(-0.28125*mpf_get_default_prec()).get_mpz_t());
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
   _eps = int_pow(Double("10.0"),mpz_class(-0.28125*mpf_get_default_prec()),false); 
   old_prec = mpf_get_default_prec();
   } 
   return _eps;
   }
*/

const Double & orsa::epsilon() {
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
    _eps = one();
    // const Double _one_over_ten("0.1");
    const Double _one_over_ten_to_the_9th("1.0e-9");
    // _eps = int_pow(Double("10.0"),mpz_class(-0.28125*mpf_get_default_prec()),false); 
    // const mpz_class _pow = abs(mpz_class(0.28125*mpf_get_default_prec()));
    // 18/64 = 9/32... don't change parenthesis!!
    // const mpz_class _pow = mpz_class("9")*(mpz_class(mpf_get_default_prec()+mpz_class("31"))/mpz_class("32"));
    //
    const mpz_class _pow9 = mpz_class(mpf_get_default_prec()+mpz_class("31"))/mpz_class("32");
    /* 
       for (unsigned int p=0; p<_pow; ++p) {
       _eps *= _one_over_ten;
       }
    */
    for (unsigned int p9=0; p9<_pow9; ++p9) {
      _eps *= _one_over_ten_to_the_9th;
    }
    //
    old_prec = mpf_get_default_prec();
  } 
  //
  /* 
     ORSA_DEBUG("default_prec.: %i",mpf_get_default_prec());
     // ORSA_DEBUG("epsilon power: %Zi",mpz_class(-0.28125*mpf_get_default_prec()).get_mpz_t());
     ORSA_DEBUG("epsilon power: %Zi",mpz_class(mpz_class("9")*(mpz_class(mpf_get_default_prec()+mpz_class("31"))/mpz_class("32"))).get_mpz_t());
     ORSA_DEBUG("epsilon......: %Fg",_eps.get_mpf_t());
  */
  //
  return _eps;
}

const Double & orsa::pi() {
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
      /*  
	  {
	  mp_exp_t t;
	  std::string str = _pi.get_str(&t,10,0);
	  std::cout << "k: " << k << "  t: " << t << "  _pi: " << str << std::endl;
	  std::cerr << "_add: " << _add << "   _eps: " << _eps << std::endl;
	  }
      */
      ++k;
    } while (fabs(_add/_pi) > _eps);
    old_prec = mpf_get_default_prec();
  }
  return _pi;
}

const Double & orsa::halfpi() {
  static Double _halfpi;
  static unsigned int old_prec = 0;
  if (old_prec != mpf_get_default_prec()) {
    _halfpi = Double("0.5")*pi();
    old_prec = mpf_get_default_prec();
  }	
  return _halfpi;
}

const Double & orsa::twopi() {
  static Double _twopi;
  static unsigned int old_prec = 0;
  if (old_prec != mpf_get_default_prec()) {
    _twopi = Double("2")*pi();
    old_prec = mpf_get_default_prec();
  }	
  return _twopi;
}

const Double & orsa::pisquared() {
  static Double _pisquared;
  static unsigned int old_prec = 0;
  if (old_prec != mpf_get_default_prec()) {
    _pisquared = pi()*pi();
    old_prec = mpf_get_default_prec();
  }	
  return _pisquared;
}

const orsa::Double & orsa::radToDeg() {
  static orsa::Double _radToDeg;
  static unsigned int old_prec = 0;
  if (old_prec != mpf_get_default_prec()) {
    _radToDeg = orsa::Double("180")/pi();
    old_prec = mpf_get_default_prec();
  }	
  return _radToDeg;
}

const orsa::Double & orsa::degToRad() {
  static orsa::Double _degToRad;
  static unsigned int old_prec = 0;
  if (old_prec != mpf_get_default_prec()) {
    _degToRad = pi()/orsa::Double("180");
    old_prec = mpf_get_default_prec();
  }	
  return _degToRad;
}

const orsa::Double & orsa::radToArcmin() {
  static orsa::Double _radToArcmin;
  static unsigned int old_prec = 0;
  if (old_prec != mpf_get_default_prec()) {
    _radToArcmin = orsa::Double("60")*orsa::Double("180")/pi();
    old_prec = mpf_get_default_prec();
  }	
  return _radToArcmin;
}

const orsa::Double & orsa::arcminToRad() {
  static orsa::Double _arcminToRad;
  static unsigned int old_prec = 0;
  if (old_prec != mpf_get_default_prec()) {
    _arcminToRad = pi()/(orsa::Double("60")*orsa::Double("180"));
    old_prec = mpf_get_default_prec();
  }	
  return _arcminToRad;
}

const orsa::Double & orsa::radToArcsec() {
  static orsa::Double _radToArcsec;
  static unsigned int old_prec = 0;
  if (old_prec != mpf_get_default_prec()) {
    _radToArcsec = orsa::Double("3600")*orsa::Double("180")/pi();
    old_prec = mpf_get_default_prec();
  }	
  return _radToArcsec;
}

const orsa::Double & orsa::arcsecToRad() {
  static orsa::Double _arcsecToRad;
  static unsigned int old_prec = 0;
  if (old_prec != mpf_get_default_prec()) {
    _arcsecToRad = pi()/(orsa::Double("3600")*orsa::Double("180"));
    old_prec = mpf_get_default_prec();
  }	
  return _arcsecToRad;
}

mpz_class orsa::kronecker(const mpz_class & i,
			  const mpz_class & j) {
  if (i==j) {
    return mpz_class("1");
  } else {
    return mpz_class("0");
  }
}

Double orsa::fmod(const Double & x,
		  const Double & y) {
  Double _ratio = x/y;
  if (fabs(_ratio*orsa::epsilon()) > orsa::one()) {
    ORSA_DEBUG("problems... x: %Fg   y: %Fg   ratio: %Fg",
	       x.get_mpf_t(),
	       y.get_mpf_t(),
	       _ratio.get_mpf_t());
    return orsa::zero();
  }
  mpz_class _n(_ratio);
  if (_ratio < 0) --_n;
  return (x-_n*y);
  // return ((_ratio-_n)*y);
}

Double orsa::cbrt(const Double & x) {
  const Double _three("3");
  Double _c = x;
  Double _old_c;
  do {
    _old_c = _c;
    // _c = (Double("2")*_c+x/(_c*_c))/Double("3");
    // _c = (two()*_c+x/(_c*_c))/_three;
    _c = (two()*_c+x/(_c*_c+epsilon()))/_three;
   // } while (fabs(_c-_old_c)/fabs(_c) > epsilon());
  } while (fabs(_c-_old_c)/fabs(_c+epsilon()) > epsilon());
  return _c;
}

Double orsa::pochhammer(const Double & a, const mpz_class & n) {
  if (n == 0) {
    return one();
  }
  if (n == 1) {
    return a;
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

Double orsa::sin(const Double & x) {
  // ORSA_DEBUG("sin(%Fg)",x.get_mpf_t());
  const Double _local_x = fmod(x,twopi());
  if (fabs(_local_x) < epsilon()) {
    return zero();
  }
  // ORSA_DEBUG("local_x: %Fg",_local_x.get_mpf_t());
  const Double _eps = epsilon();
  Double _s("0");
  Double _add;
  mpz_class k("0");
  mpz_class twice_k_plus_one;
  do {
    twice_k_plus_one = 2*k+1;
    _add = power_sign(k)*int_pow(_local_x,twice_k_plus_one)/factorial(twice_k_plus_one);
    _s += _add;
    /* 
       ORSA_DEBUG("k: %Zi   s: %Fg",
       k.get_mpz_t(),
       _s.get_mpf_t());
    */
    ++k;
  } while (fabs(_add/_s) > _eps);
  return _s;
}

Double orsa::cos(const Double & x) {
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
    /* {
       mp_exp_t t;
       std::string str = _c.get_str(&t,10,0);
       std::cout << "k: " << k << "  t: " << t << "  _c: " << str << std::endl;
       }
    */
    ++k;
  } while (fabs(_add/_c) > _eps);
  return _c;
}

void orsa::sincos(const Double & x, Double & s, Double & c) {
  s = sin(x);
  // c = sqrt(1.0-s*s); // Warning: sign problem, need value of pi...
  c = cos(x);
}

Double orsa::tan(const Double & x) {
  return (sin(x)/cos(x));
}

Double orsa::asin(const Double & x) {
  // ORSA_DEBUG("asin(%20.12Fg) called...",x.get_mpf_t());
  //
  if (fabs(x) < epsilon()) {
    return zero();
  }
  //
  if (x > one()) {
    if (x > 1.01) {
      ORSA_ERROR("out-of-domain error: x = %Fg > 1.0",x.get_mpf_t());
    }
    return halfpi();
  } else if (x < (-one())) {
    if (x < -1.01) {
      ORSA_ERROR("out-of-domain error: x = %Fg < -1.0",x.get_mpf_t());
    }
    return -halfpi();
  }
  //
  const Double _eps = epsilon();
  const Double _one_over_two("0.5");
  Double _as;
  Double _add;
  mpz_class k;
  mpz_class twice_k_plus_one;
  //
  if (x > Double("0.5")) {
    const Double _one_minus_x(one()-x);
    _as = zero();
    k = 0;
    do {
      // ORSA_DEBUG("asin: k: %Zi",k.get_mpz_t());
      twice_k_plus_one = 2*k+1;
      _add = pochhammer(_one_over_two,k)*int_pow(_one_minus_x,k)/(int_pow(two(),k)*twice_k_plus_one*factorial(k));
      _as += _add;
      ++k;
      if (k > mpf_get_default_prec()) {
	ORSA_WARNING("max number of iterations hit while computing asin(%20.12Fg)",x.get_mpf_t());
	break;
      }
    } while (fabs(_add/_as) > _eps);
    _as = halfpi() - sqrt(two())*sqrt(_one_minus_x)*_as;
  } else if (x < Double("-0.5")) {
    const Double _one_plus_x(one()+x);
    _as = zero();
    k = 0;
    do {
      // ORSA_DEBUG("asin: k: %Zi",k.get_mpz_t());
      twice_k_plus_one = 2*k+1;
      _add = pochhammer(_one_over_two,k)*int_pow(_one_plus_x,k)/(int_pow(two(),k)*twice_k_plus_one*factorial(k));
      _as += _add;
      ++k;
      if (k > mpf_get_default_prec()) {
	ORSA_WARNING("max number of iterations hit while computing asin(%20.12Fg)",x.get_mpf_t());
	break;
      }
    } while (fabs(_add/_as) > _eps);
    _as = sqrt(two())*sqrt(_one_plus_x)*_as - halfpi();
  } else {
    _as = zero();
    k = 0;
    do {
      // ORSA_DEBUG("asin: k: %Zi",k.get_mpz_t());
      twice_k_plus_one = 2*k+1;
      _add = pochhammer(_one_over_two,k)*int_pow(x,twice_k_plus_one)/(twice_k_plus_one*factorial(k));
      _as += _add;
      ++k;
      if (k > mpf_get_default_prec()) {
	ORSA_WARNING("max number of iterations hit while computing asin(%20.12Fg)",x.get_mpf_t());
	break;
      }
    } while (fabs(_add/_as) > _eps);
  }
  
  return _as;
}
   
Double orsa::acos(const Double & x) {
  return (halfpi()-asin(x));
}

Double orsa::atan2(const Double & y,
		   const Double & x) {
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

Double orsa::atan(const Double & x) {
  orsa::Double tmpVal = fmod(atan2(x,orsa::one()),pi());
  if (tmpVal > halfpi()) {
    tmpVal -= pi();
  }	
  return tmpVal;
}

Double orsa::log(const Double & x) {
  // ORSA_DEBUG("called log(%Ff)",x.get_mpf_t());
  if (x <= zero()) {
    ORSA_ERROR("called log() with non-positive argument");
    return zero();
  }
  //
  if (fabs(x-one()) < epsilon()) {
    return zero();
  }
  //
  mpz_class _q("1");
  Double _val(x);
  // domain limits are 0->2, but convergence is faster closer to 1.0
  while (_val > Double("1.1")) {
    _val = sqrt(_val);
    _q *= mpz_class("2");
  }
  //
  while (_val < Double("0.9")) {
    _val = sqrt(_val);
    _q *= mpz_class("2");
  }
  //
  const Double _val_minus_one = _val - one();
  // std::cerr << "_q: " << _q << "  _val: " << _val << std::endl;
  //
  const Double _eps = epsilon();
  Double _log("0");
  Double _add;
  mpz_class k("1");
  do {
    _add = power_sign(k+1)*int_pow(_val_minus_one,k)/Double(k);
    _log += _add;
    //
    /* 
       {
       mp_exp_t t;
       std::string str = _log.get_str(&t,10,0);
       std::cout << "k: " << k << "  t: " << t << "  _log: " << str << std::endl;
       }
    */
    //
    ++k;
  } while (fabs(_add/_log) > _eps);
  // ORSA_DEBUG("log(%Ff) = %Ff",x.get_mpf_t(),orsa::Double(_q*_log).get_mpf_t());
  return (_q*_log);
}

Double orsa::log(const Double & base,
		 const Double & x) {
  return (log(x)/log(base));
}

Double orsa::log10(const Double & x) {
  return log(10,x);
}

Double orsa::exp(const Double & x) {
  if (fabs(x) < epsilon()) {
    return one();
  }
  // 
  // ix = an integer close to x, so that x/ix is close to one
  mpz_class ix(1+fabs(x)); 
  const Double val(x/ix);
  //
  /* 
     ORSA_DEBUG("ix: %Zi   val: %Ff",
     ix.get_mpz_t(),
     val.get_mpf_t());
  */
  //
  const Double _eps = epsilon();
  Double _exp("0");
  Double _add;
  mpz_class k("0");
  do {
    _add = int_pow(val,k)/factorial(k);
    _exp += _add;
    //
    /* 
       ORSA_DEBUG("k: %Zi   _exp: %Fg   _add: %Fg",
       k.get_mpz_t(),
       _exp.get_mpf_t(),
       _add.get_mpf_t());
    */
    //
    ++k;
  } while (fabs(_add/_exp) > _eps);
  /* 
     ORSA_DEBUG("exp(%Ff) = %Ff",
     x.get_mpf_t(),
     int_pow(_exp,ix).get_mpf_t());
  */
  return (int_pow(_exp,ix));
}

Double orsa::cosh(const Double & x) {
  return (exp(x)+exp(-x))/two();
}

Double orsa::sinh(const Double & x) {
  return (exp(x)-exp(-x))/two();
}
