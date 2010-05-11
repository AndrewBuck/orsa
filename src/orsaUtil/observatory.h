#ifndef _ORSA_UTIL_OBSERVATORY_
#define _ORSA_UTIL_OBSERVATORY_

#include <orsa/body.h>
#include <orsa/bodygroup.h>
#include <orsa/cache.h>
#include <orsa/datetime.h>
#include <orsa/double.h>
#include <orsa/vector.h>

#include <orsaSolarSystem/observation.h>
#include <orsaSolarSystem/observatory.h>

#include <string>

#include <osg/Referenced>
#include <osg/ref_ptr>

namespace orsaInputOutput {
    class MPCObsCodeFile;
}

namespace orsaUtil {
  
    class StandardObservatoryPositionCallback : public orsaSolarSystem::ObservatoryPositionCallback {
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
    public:   
        const orsaSolarSystem::Observatory & getObservatory(const std::string & obsCode) const;
        const orsaSolarSystem::Observatory & getObservatory(const orsaSolarSystem::Observation * obs) const;
    protected:
        osg::ref_ptr<orsa::BodyGroup> bg;   
        osg::ref_ptr<orsa::Body> earth; 
        osg::ref_ptr<orsaInputOutput::MPCObsCodeFile> obsCodeFile;
    };
  
}; // namespace orsaUtil

#endif // _ORSA_UTIL_OBSERVATORY_
