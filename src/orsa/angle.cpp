#include <orsa/angle.h>

using namespace orsa;

void Angle::setRad(const orsa::Double & a) {
  // ORSA_DEBUG("args: a=%Ff",a.get_mpf_t());
  _rad = a;
}

const orsa::Double & Angle::getRad(orsa::Double & a) const {
  a = _rad;
  return _rad;
}

const orsa::Double & Angle::getRad() const {
  return _rad;
}
  
/* 
   void Angle::setDPS(const orsa::Double d, const orsa::Double p, const orsa::Double s) {
   if ( (d > 0) ||
   ( (d == 0) && (p > 0) ) ||
   ( (d == 0) && (p == 0) && (s >= 0) ) ) {
   _rad = (pi()/180)*(d+p/60.0+s/3600.0);
   } else {
   _rad = (pi()/180)*(d-p/60.0-s/3600.0);
   }
   }
*/
  
int Angle::check_sign(int sign) {
  if (sign == 0) {
    ORSA_ERROR("Hmmm, sign equal to zero...");
    return 1;
  } else {
    return (sign/abs(sign));
  }
}
  
void Angle::setDPS(const orsa::Double & d, 
		   const orsa::Double & p, 
		   const orsa::Double & s,
		   const int sign) {   
  // ORSA_DEBUG("args: d=%Ff   p=%Ff   s=%Ff   sign: %i",d.get_mpf_t(),p.get_mpf_t(),s.get_mpf_t(),sign);
  _rad = check_sign(sign)*degToRad()*(d+p/60.0+s/3600.0);
}

/* 
   void Angle::getDPS(orsa::Double & d, orsa::Double & p, orsa::Double & s) const {
   orsa::Double frac;
   orsa::Double fdeg = (180/pi())*_rad;
   if (fdeg < 0.0) {
   d    = -floor(-fdeg);
   frac = d - fdeg;
   } else {
   d    = floor(fdeg);
   frac = fdeg - d;
   }
   //
   p = floor(frac*60.0);  
   s = frac*3600.0 - p*60.0;
   } 
*/
  
void Angle::getDPS(orsa::Double & d, 
		   orsa::Double & p, 
		   orsa::Double & s,
		   int & sign) const {
  // const orsa::Double fdeg = (180.0/pi())*_rad;
  /* 
     if (fdeg < 0.0) {
     sign = -1;
     d    = -floor(-fdeg);
     frac = d - fdeg;
     } else {
     sign = 1;
     d    = floor(fdeg);
     frac = fdeg - d;
     }
  */
  //
  if (_rad < 0.0) {
    sign = -1;
  } else {
    sign = 1;
  }
  //
  const orsa::Double abs_fdeg = fabs(radToDeg()*_rad);
  d = floor(abs_fdeg);
  const orsa::Double frac = abs_fdeg - d;
  p = floor(frac*60.0);  
  s = frac*3600.0 - p*60.0;
} 

/* 
   void Angle::setHMS(const orsa::Double h, const orsa::Double m, const orsa::Double s) {
   if (h >= 0) 
   _rad = 15*(pi()/180)*(h+m/60.0+s/3600.0);
   else
   _rad = 15*(pi()/180)*(h-m/60.0-s/3600.0);
   }
*/
  
void Angle::setHMS(const orsa::Double & h, 
		   const orsa::Double & m,
		   const orsa::Double & s, 
		   const int sign) { 
  // ORSA_DEBUG("args: h=%Ff   m=%Ff   s=%Ff   sign: %i",h.get_mpf_t(),m.get_mpf_t(),s.get_mpf_t(),sign);
  _rad = check_sign(sign)*15.0*degToRad()*(h+m/60.0+s/3600.0);
}

void Angle::getHMS(orsa::Double & h, 
		   orsa::Double & m, 
		   orsa::Double & s, 
		   int & sign) const {
  /* 
     orsa::Double frac;
     orsa::Double fh = (180/pi())*_rad/15.0;
     if (fh < 0.0) {
     sign = -1;
     h    = -floor(-fh);
     frac = h - fh;
     } else {
     sign = 1;
     h    = floor(fh);
     frac = fh - h;
     }
  */
  //
  if (_rad < 0.0) {
    sign = -1;
  } else {
    sign = 1;
  }
  //
  const orsa::Double abs_fh = fabs(radToDeg()*_rad/15.0);
  h = floor(abs_fh);
  const orsa::Double frac = abs_fh - h;
  m = floor(frac*60.0);  
  s = frac*3600.0 - m*60.0; 
}    

