#ifndef _ORSA_ANGLE_
#define _ORSA_ANGLE_

#include <orsa/debug.h>
#include <orsa/double.h>

namespace orsa {
  
  class Angle {
  public:
    Angle() : _rad(zero()) { }
  public:
    Angle(const orsa::Double & x) : _rad(x) { }
      
  public:
    void  setRad(const orsa::Double &);
    const orsa::Double & getRad(orsa::Double &) const;
    const orsa::Double & getRad() const;
    //
    void setDPS(const orsa::Double & d, 
		const orsa::Double & p, 
		const orsa::Double & s, 
		const int sign = 1);
    
    void getDPS(orsa::Double & d, 
		orsa::Double & p, 
		orsa::Double & s, 
		int & sign) const;
    //
    void setHMS(const orsa::Double & h, 
		const orsa::Double & m, 
		const orsa::Double & s, 
		const int sign = 1);
   
    void getHMS(orsa::Double & h, 
		orsa::Double & m, 
		orsa::Double & s,
		int & sign) const;
    
  private:
    int check_sign(int sign);
    
  private:
    orsa::Double _rad;
  };
  
  inline orsa::Double sin(const Angle & alpha) {
    return orsa::sin(alpha.getRad());
  }
  
  inline orsa::Double cos(const Angle & alpha) {
    return orsa::cos(alpha.getRad());
  }
  
  inline orsa::Double tan(const Angle & alpha) {
    return orsa::tan(alpha.getRad());
  }
  
  inline void sincos(const Angle & alpha, orsa::Double & s, orsa::Double & c) {
    orsa::sincos(alpha.getRad(),s,c); 
  }
  
} // namespace orsa

#endif // _ORSA_ANGLE_
