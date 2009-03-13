#ifndef _ORSA_SOLAR_SYSTEM_OBLEQ_
#define _ORSA_SOLAR_SYSTEM_OBLEQ_

#include <orsa/angle.h>
#include <orsa/datetime.h>
#include <orsa/matrix.h>

namespace orsaSolarSystem {
  
  orsa::Angle obleq(const orsa::Time &);
  
  orsa::Angle obleqJ2000();
  
  orsa::Matrix eclipticToEquatorial();
  orsa::Matrix equatorialToEcliptic();
  
}; // namespace orsaSolarSystem

#endif // _ORSA_SOLAR_SYSTEM_OBLEQ_
