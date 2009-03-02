#ifndef _ORSA_SPICE_SPICE_BODY_INTERPOLATED_POS_VEL_CALLBACK_
#define _ORSA_SPICE_SPICE_BODY_INTERPOLATED_POS_VEL_CALLBACK_

#include <orsa/body.h>
#include <orsa/bodygroup.h>

namespace orsaSPICE {
  
  /* 
     class SpiceBodyInterpolatedPosVelCallback : public orsa::BodyPosVelCallback {
     public:    
     SpiceBodyInterpolatedPosVelCallback(const std::string & name,
     const orsa::Time  & start,
     const orsa::Time  & stop,
     const orsa::Time  & samplingPeriod);
     public:
     void getPosVel(const orsa::Time  & t,
     orsa::Vector      & relativePosition,
     orsa::Vector      & relativeVelocity) const;
     private:
     osg::ref_ptr<orsa::Interval<orsa::BodyGroup::TRV> > _interval;
     };
  */
  
}; // namespace orsaSPICE

#endif
