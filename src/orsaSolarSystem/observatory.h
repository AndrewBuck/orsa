#ifndef _ORSA_SOLAR_SYSTEM_OBSERVATORY_
#define _ORSA_SOLAR_SYSTEM_OBSERVATORY_

#include <orsa/body.h>
#include <orsa/bodygroup.h>
#include <orsa/cache.h>
#include <orsa/datetime.h>
#include <orsa/double.h>
#include <orsa/vector.h>

#include <orsaSolarSystem/observation.h>

#include <string>

#include <osg/Referenced>
#include <osg/ref_ptr>

namespace orsaInputOutput {
  class MPCObsCodeFile;
}

namespace orsaSolarSystem {
  
  class Observatory {
  public:
    orsa::Cache<double> lon, pxy, pz;
    orsa::Cache<std::string> obsCode;
    orsa::Cache<std::string> name;
  public:
    // satellite or roving observer
    bool moving() const {
      return ((!lon.isSet()) && (!pxy.isSet()) && (!pz.isSet()));
    }
  };
  
  class ObservatoryPositionCallback : public osg::Referenced {
  public:
    ObservatoryPositionCallback() : osg::Referenced(true) { }
  protected:
    virtual ~ObservatoryPositionCallback() { }
  public:
    // absolute position
    virtual bool getPosition(orsa::Vector & position,
			     const orsaSolarSystem::Observation * obs) const = 0;
  public:
    // useful wrapper, equivalent if the observatory is not "moving"
    virtual bool getPosition(orsa::Vector & position,
			     const std::string & obsCode,
			     const orsa::Time  & t) const = 0;
  public:
    // absolute position and velocity
    virtual bool getPosVel(orsa::Vector & position,
			   orsa::Vector & velocity,
			   const orsaSolarSystem::Observation * obs) const = 0;
  public:
    // useful wrapper, equivalent if the observatory is not "moving"
    virtual bool getPosVel(orsa::Vector & position,
			   orsa::Vector & velocity,
			   const std::string & obsCode,
			   const orsa::Time  & t) const = 0;
  };
  
  class StandardObservatoryPositionCallback : public ObservatoryPositionCallback {
  public:
    StandardObservatoryPositionCallback(orsaInputOutput::MPCObsCodeFile * ocf);
  public:
    bool getPosition(orsa::Vector & position,
		     const orsaSolarSystem::Observation * obs) const; 
  public:
    bool getPosition(orsa::Vector & position,
		     const std::string & obsCode,
		     const orsa::Time  & t) const;
  public:
    bool getPosVel(orsa::Vector      & position,
		   orsa::Vector      & velocity,
		   const orsaSolarSystem::Observation * obs) const;
  public:
    bool getPosVel(orsa::Vector & position,
		   orsa::Vector & velocity,
		   const std::string & obsCode,
		   const orsa::Time  & t) const;
  protected:
    osg::ref_ptr<orsa::BodyGroup> bg;   
    osg::ref_ptr<orsa::Body> earth; 
    osg::ref_ptr<orsaInputOutput::MPCObsCodeFile> obsCodeFile;
  };
  
}; // namespace orsaSolarSystem

#endif // _ORSA_SOLAR_SYSTEM_OBSERVATORY_
