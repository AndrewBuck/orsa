#ifndef _ORSA_PROPULSION_
#define _ORSA_PROPULSION_

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/cache.h>
#include <orsa/datetime.h>
#include <orsa/double.h>
#include <orsa/matrix.h>
#include <orsa/vector.h>

#include <vector>

namespace orsa {
  
  // class Body;
  // class BodyGroup;
  
  // NOTE: be sure to include the reaction mass into the 
  //       total body mass computation at a given time t.
  
  // #warning "when using a Propulsion, in general you need to set orsa::Interaction::dependsOnVelocity() as TRUE"
  
  class Propulsion : public osg::Referenced {
  public:	
    Propulsion() : osg::Referenced(true) { }
  protected:
    virtual ~Propulsion() { }
  public:	
    virtual orsa::Vector getThrust(const orsa::Time & t) const = 0;
  public:
    // Thrust is ON at the "start" event time, and OFF at the "stop" event time (semi-open interval)
    // This function must always set t to a time different from the input t, or return false
    // The value of sign (+/- 1) determines the direction of the new t.
    virtual bool nextEventTime(orsa::Time      & t,
			       const mpz_class & sign) const = 0;
  };
  
} // namespace orsa

#endif // _ORSA_PROPULSION_
