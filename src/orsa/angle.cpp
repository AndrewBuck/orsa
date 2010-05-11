#include <orsa/angle.h>

using namespace orsa;

int Angle::check_sign(int sign) {
    if (sign == 0) {
        ORSA_ERROR("sign equal to zero...");
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
    _rad = check_sign(sign)*degToRad()*(d+p/60+s/3600);
}

void Angle::getDPS(double & d, 
                   double & p, 
                   double & s,
                   int & sign) const {
    if (_rad < 0) {
        sign = -1;
    } else {
        sign = 1;
    }
    //
    const double abs_fdeg = fabs(radToDeg()*_rad);
    d = floor(abs_fdeg);
    const double frac = abs_fdeg - d;
    p = floor(frac*60);  
    s = frac*3600 - p*60;
} 

void Angle::setHMS(const double & h, 
                   const double & m,
                   const double & s, 
                   const int sign) { 
    // ORSA_DEBUG("args: h=%f   m=%f   s=%f   sign: %i",h(),m(),s(),sign);
    _rad = check_sign(sign)*15*degToRad()*(h+m/60+s/3600);
}

void Angle::getHMS(double & h, 
                   double & m, 
                   double & s, 
                   int & sign) const {
    if (_rad < 0) { 
        sign = -1; 
    } else { 
        sign = 1; 
    } 
    //
    const double abs_fh = fabs(radToDeg()*_rad/15);
    h = floor(abs_fh);
    const double frac = abs_fh - h;
    m = floor(frac*60);  
    s = frac*3600 - m*60; 
}    
