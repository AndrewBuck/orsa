#ifndef __MULTIMIN_PHASE__
#define __MULTIMIN_PHASE__

#include <orsa/double.h>
#include <orsa/multimin.h> 
#include <orsa/vector.h>

class MultiminPhase : public orsa::Multimin {
 public:
  orsa::Double fun(const orsa::MultiminParameters * par) const;
 public:
  orsa::Double getAlpha(const orsa::Double & phaseAngle,
			const orsa::Vector & uSun,
			const orsa::Vector & uInclination);
 protected:
  orsa::Cache<orsa::Double> phi;
  orsa::Cache<orsa::Vector> uS, uI;  
};

#endif // __MULTIMIN_PHASE__
