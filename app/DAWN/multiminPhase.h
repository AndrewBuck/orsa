#ifndef __MULTIMIN_PHASE__
#define __MULTIMIN_PHASE__

#include <orsa/double.h>
#include <orsa/multimin.h> 
#include <orsa/vector.h>

class MultiminPhase : public orsa::Multimin {
 public:
  double fun(const orsa::MultiminParameters * par) const;
 public:
  double getAlpha(const double & phaseAngle,
		  const orsa::Vector & uSun,
		  const orsa::Vector & uInclination);
 protected:
  orsa::Cache<double> phi;
  orsa::Cache<orsa::Vector> uS, uI;  
};

#endif // __MULTIMIN_PHASE__
