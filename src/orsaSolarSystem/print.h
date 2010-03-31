#ifndef _ORSA_SOLAR_SYSTEM_PRINT_
#define _ORSA_SOLAR_SYSTEM_PRINT_

#include <orsaSolarSystem/data.h>
#include <orsaSolarSystem/datetime.h>

namespace orsaSolarSystem {
  
  inline void print(const orsa::Time & t) {
    
    // ORSA_DEBUG("t: %20.14f [day]",orsa::FromUnits(t.get_d(),orsa::Unit::DAY,-1));
    
    int y,m,d,H,M,S,ms;
    // double fd;
    int y_UTC,m_UTC,d_UTC,H_UTC,M_UTC,S_UTC,ms_UTC;
    orsaSolarSystem::gregorDay(t,y,m,d,H,M,S,ms);
    orsaSolarSystem::gregorDay(orsaSolarSystem::ToTimeScale(t,orsaSolarSystem::TS_UTC),
			       y_UTC,m_UTC,d_UTC,H_UTC,M_UTC,S_UTC,ms_UTC);
    ORSA_DEBUG("t: %Zi [musec] = %.12f [day] = JD %.5f (TDT) = JD %.5f (UTC) = %i %2i %2i %02i:%02i:%02i.%03i (TDT) = %i %2i %2i %02i:%02i:%02i.%03i (UTC)",
	       t.getMuSec().get_mpz_t(),
	       orsa::FromUnits(t.get_d(),orsa::Unit::DAY,-1),
	       orsaSolarSystem::timeToJulian(t),
	       orsaSolarSystem::timeToJulian(orsaSolarSystem::ToTimeScale(t,orsaSolarSystem::TS_UTC)),
	       y,m,d,H,M,S,ms,
	       y_UTC,m_UTC,d_UTC,H_UTC,M_UTC,S_UTC,ms_UTC);
  }
  
  inline void print(const orsa::Orbit & o) {
    const double G    = orsa::Unit::G();
    const double MSun = orsaSolarSystem::Data::MSun();
    ORSA_DEBUG("orbit:\n"
	       "   a: %g [AU] = %g [km]\n"
	       "   e: %g\n"
	       "   i: %g [deg]\n"
	       "node: %g [deg]\n"
	       "peri: %g [deg]\n"
	       "   M: %g [deg]\n"
	       "mu/G: %g [MSun] = %g [kg]",
	       orsa::FromUnits(o.a,orsa::Unit::AU,-1),
	       orsa::FromUnits(o.a,orsa::Unit::KM,-1),
	       o.e,
	       orsa::radToDeg()*o.i,
	       orsa::radToDeg()*o.omega_node,
	       orsa::radToDeg()*o.omega_pericenter,
	       orsa::radToDeg()*o.M,
	       o.mu/G/MSun,
	       orsa::FromUnits(o.mu/G,orsa::Unit::KG,-1));
  }
  
} // namespace orsaSolarSystem

#endif // _ORSA_SOLAR_SYSTEM_PRINT_
