#ifndef _ORSA_ANGLE_
#define _ORSA_ANGLE_

#include <orsa/debug.h>
#include <orsa/double.h>

namespace orsa {
  
  class Angle {
  public:
    Angle() : _rad(0) { }
  public:
    Angle(const double & x) : _rad(x) { }
      
  public:
    void  setRad(const double &);
    const double & getRad(double &) const;
    const double & getRad() const;
    //
    void setDPS(const double & d, 
		const double & p, 
		const double & s, 
		const int sign = 1);
    
    void getDPS(double & d, 
		double & p, 
		double & s, 
		int & sign) const;
    //
    void setHMS(const double & h, 
		const double & m, 
		const double & s, 
		const int sign = 1);
   
    void getHMS(double & h, 
		double & m, 
		double & s,
		int & sign) const;
    
  private:
    int check_sign(int sign);
    
  private:
    double _rad;
  };
  
  inline double sin(const Angle & alpha) {
    return orsa::sin(alpha.getRad());
  }
  
  inline double cos(const Angle & alpha) {
    return orsa::cos(alpha.getRad());
  }
  
  inline double tan(const Angle & alpha) {
    return orsa::tan(alpha.getRad());
  }
  
  inline void sincos(const Angle & alpha, double & s, double & c) {
    orsa::sincos(alpha.getRad(),s,c); 
  }
  
} // namespace orsa

#endif // _ORSA_ANGLE_
