#ifndef _ORSA_SOLAR_SYSTEM_OBSERVATION_
#define _ORSA_SOLAR_SYSTEM_OBSERVATION_

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/angle.h>
#include <orsa/cache.h>
#include <orsa/datetime.h>
#include <orsa/double.h>

#include <string>

namespace orsaSolarSystem {
  
  class Observation : public osg::Referenced {
  public:
    Observation() : Referenced(true) { }
  protected:
    virtual ~Observation() { }
    
  public:
    orsa::Cache<unsigned int> number; // code needed, and MPC->real number code...
    orsa::Cache<std::string>  designation, obsCode, magCode;
    orsa::Cache<bool>         discovery;
    orsa::Cache<orsa::Time>   epoch;
    orsa::Cache<orsa::Angle>  ra, dec;
    orsa::Cache<double> mag;
  };
  
  class RadarObservation : public osg::Referenced {
  public:
    RadarObservation() : Referenced(true) { }
  protected:
    virtual ~RadarObservation() { }
    // TO BE CONTINUED...
  };
  
  typedef std::vector< osg::ref_ptr<orsaSolarSystem::Observation> > ObservationVector;
  
} // namespace orsaSolarSystem

#endif // _ORSA_SOLAR_SYSTEM_OBSERVATION_
