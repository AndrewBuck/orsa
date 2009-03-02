#ifndef _ORSA_EULER_H_
#define _ORSA_EULER_H_

#include <orsa/double.h>
#include <orsa/matrix.h>
#include <orsa/vector.h>

namespace orsa {
  
  bool Euler(orsa::Vector       & omegaDot,
	     const orsa::Vector & omega,
	     const orsa::Matrix & I,
	     const orsa::Vector & T,
	     const bool usePrincipalAxisRotationTransformation = true);
  
} // namespace orsa

#endif // _ORSA_EULER_H_
