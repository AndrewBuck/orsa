#ifndef _ORSA_SOLAR_SYSTEM_OBSERVATORY_
#define _ORSA_SOLAR_SYSTEM_OBSERVATORY_

#include <orsa/cache.h>
#include <orsa/datetime.h>
#include <orsa/double.h>
#include <orsa/vector.h>

#include <string>

#include <osg/Referenced>
#include <osg/ref_ptr>

namespace orsaSolarSystem {
  
  class Observatory {
  public:
    orsa::Cache<double> lon, pxy, pz;
    orsa::Cache<std::string> obsCode;
    orsa::Cache<std::string> name;
  };
  
  class ObservatoryPositionCallback : public osg::Referenced {
  public:
    ObservatoryPositionCallback() : osg::Referenced(true) { }
  protected:
    virtual ~ObservatoryPositionCallback() { }
  public:
    // absolute position
    virtual bool getPosition(orsa::Vector & position,
			     const std::string & obsCode,
			     const orsa::Time &) const = 0;
  };
  
}; // namespace orsaSolarSystem

#endif // _ORSA_SOLAR_SYSTEM_OBSERVATORY_
