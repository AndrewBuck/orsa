#ifndef _ORSA_SOLAR_SYSTEM_ORBIT_
#define _ORSA_SOLAR_SYSTEM_ORBIT_

#include <orsa/bodygroup.h>
#include <orsa/double.h>
#include <orsa/integrator_radau.h>
#include <orsa/multifit.h>
#include <orsa/orbit.h>
#include <orsa/print.h>

#include <orsaSolarSystem/data.h>
#include <orsaSolarSystem/observation.h>
#include <orsaSolarSystem/observatory.h>
#include <orsaSolarSystem/obleq.h>

#include <orsaSPICE/spice.h>
#include <orsaSPICE/spiceBodyTranslationalCallback.h>

namespace orsaSolarSystem {
  
  class OrbitWithEpoch : public orsa::Orbit {
  public:
    orsa::Cache<orsa::Time> epoch;
  };
  
}; // namespace orsaSolarSystem

#endif // _ORSA_SOLAR_SYSTEM_ORBIT_
