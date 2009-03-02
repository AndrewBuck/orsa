#include <orsa/angle.h>

using namespace orsa;

void Angle::setRad(const double & a) {
  // ORSA_DEBUG("args: a=%f",a());
  _rad = a;
}

const double & Angle::getRad(double & a) const {
  a = _rad;
  return _rad;
}

const double & Angle::getRad() const {
  return _rad;
}
  
/* 
   void Angle::setDPS(const double d, const double p, const double s) {
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
  
void Angle::setDPS(const double & d, 
		   const double & p, 
		   const double & s,
		   const int sign) {   
  // ORSA_DEBUG("args: d=%f   p=%f   s=%f   sign: %i",d(),p(),s(),sign);
  _rad = check_sign(sign)*degToRad()*(d+p/60.0+s/3600.0);
}

/* 
   void Angle::getDPS(double & d, double & p, double & s) const {
   double frac;
   double fdeg = (180/pi())*_rad;
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
  
void Angle::getDPS(double & d, 
		   double & p, 
		   double & s,
		   int & sign) const {
  // const double fdeg = (180.0/pi())*_rad;
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
  const double abs_fdeg = fabs(radToDeg()*_rad);
  d = floor(abs_fdeg);
  const double frac = abs_fdeg - d;
  p = floor(frac*60.0);  
  s = frac*3600.0 - p*60.0;
} 

/* 
   void Angle::setHMS(const double h, const double m, const double s) {
   if (h >= 0) 
   _rad = 15*(pi()/180)*(h+m/60.0+s/3600.0);
   else
   _rad = 15*(pi()/180)*(h-m/60.0-s/3600.0);
   }
*/
  
void Angle::setHMS(const double & h, 
		   const double & m,
		   const double & s, 
		   const int sign) { 
  // ORSA_DEBUG("args: h=%f   m=%f   s=%f   sign: %i",h(),m(),s(),sign);
  _rad = check_sign(sign)*15.0*degToRad()*(h+m/60.0+s/3600.0);
}

void Angle::getHMS(double & h, 
		   double & m, 
		   double & s, 
		   int & sign) const {
  /* 
     double frac;
     double fh = (180/pi())*_rad/15.0;
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
  const double abs_fh = fabs(radToDeg()*_rad/15.0);
  h = floor(abs_fh);
  const double frac = abs_fh - h;
  m = floor(frac*60.0);  
  s = frac*3600.0 - m*60.0; 
}    

