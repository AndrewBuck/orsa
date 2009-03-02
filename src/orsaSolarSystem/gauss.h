#ifndef ORSA_SOLAR_SYSTEM_GAUSS_H
#define ORSA_SOLAR_SYSTEM_GAUSS_H

#include <vector>

#include <orsaSolarSystem/observation.h>
#include <orsaSolarSystem/observatory.h>
#include <orsaSolarSystem/orbit.h>

namespace orsaSolarSystem {
  
  void GaussMethod(std::vector<orsaSolarSystem::OrbitWithEpoch> & preliminaryOrbitVector,
		   const orsaSolarSystem::ObservationVector & observationVector,
		   const orsaSolarSystem::ObservatoryPositionCallback * obsPosCB,
		   const orsa::Body * refBody,
		   orsa::BodyGroup * bg);
  
} // namespace orsaSolarSystem

#endif // ORSA_SOLAR_SYSTEM_GAUSS_H

