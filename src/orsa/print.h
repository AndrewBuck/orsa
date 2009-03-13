#ifndef _ORSA_PRINT_
#define _ORSA_PRINT_

#include <orsa/angle.h>
#include <orsa/datetime.h>
#include <orsa/matrix.h>
#include <orsa/quaternion.h>

#include <orsaSolarSystem/datetime.h>

namespace orsa {
  
  inline void print(const double & d) {
    ORSA_DEBUG("d: %+20.12g",d);
  }
  
  inline void print(const orsa::Angle & a) {
    double H,M,S;
    int sign_HMS;
    double d,p,s;
    int sign_dps;
    a.getHMS(H,M,S,sign_HMS);
    a.getDPS(d,p,s,sign_dps);
    ORSA_DEBUG("a: %+17.12f [rad] = %+17.12f [deg] = %+3.0fH %02.0fM %06.3fS = %+4.0fd %02.0fp %06.3fs",
	       a.getRad(),
	       a.getRad()*orsa::radToDeg(),
	       sign_HMS*H,M,S,
	       sign_dps*d,p,s);
  }
  
  inline void print(const orsa::Vector & v) {
    ORSA_DEBUG("v: length: %.12e [%+.12e,%+.12e,%+.12e]",
	       v.length(),
	       v.getX(),
	       v.getY(),
	       v.getZ());
  }
  
  inline void print(const orsa::Matrix & m) {
    ORSA_DEBUG("m: det(m) = %.12e\n"
	       "[%+.12e,%+.12e,%+.12e]\n"
	       "[%+.12e,%+.12e,%+.12e]\n"
	       "[%+.12e,%+.12e,%+.12e]",
	       m.determinant(),
	       m.getM11(),  m.getM12(),  m.getM13(), 
	       m.getM21(),  m.getM22(),  m.getM23(), 
	       m.getM31(),  m.getM32(),  m.getM33());  
  }
  
  inline void print(const orsa::Quaternion & q) {
    ORSA_DEBUG("q: s: %+9.6e   v: %9.6e x [%6.3f,%6.3f,%6.3f]  l: %9.6e",
	       q.getScalar(),
	       q.getVector().length(),
	       q.getVector().normalized().getX(),
	       q.getVector().normalized().getY(),
	       q.getVector().normalized().getZ(),
	       q.length());
  }
  
  inline void print(const orsa::Time & t) {
    
    // ORSA_DEBUG("t: %20.14f [day]",orsa::FromUnits(t.get_d(),orsa::Unit::DAY,-1));
    
    int y,m,d,H,M,S,ms;
    // double fd;
    int y_UTC,m_UTC,d_UTC,H_UTC,M_UTC,S_UTC,ms_UTC;
    orsaSolarSystem::gregorDay(t,y,m,d,H,M,S,ms);
    // orsaSolarSystem::gregorDay(t,y,m,d,fd);
    orsaSolarSystem::gregorDay(orsaSolarSystem::ToTimeScale(t,orsaSolarSystem::TS_UTC),
			       y_UTC,m_UTC,d_UTC,H_UTC,M_UTC,S_UTC,ms_UTC);
    ORSA_DEBUG("t: %.12f [day] = JD %.5f (TDT) = JD %.5f (UTC) = %i %i %i %02i:%02i:%02i.%03i (TDT) = %i %i %i %02i:%02i:%02i.%03i (UTC)",
	       orsa::FromUnits(t.get_d(),orsa::Unit::DAY,-1),
	       orsaSolarSystem::timeToJulian(t),
	       orsaSolarSystem::timeToJulian(orsaSolarSystem::ToTimeScale(t,orsaSolarSystem::TS_UTC)),
	       y,m,d,H,M,S,ms,
	       y_UTC,m_UTC,d_UTC,H_UTC,M_UTC,S_UTC,ms_UTC);
  }
  
} // namespace orsa

#endif // _ORSA_PRINT_
