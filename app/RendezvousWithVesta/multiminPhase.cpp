#include "multiminPhase.h"

#include <orsa/double.h>
#include <orsa/matrix.h>

#include <algorithm>

using namespace orsa;

orsa::Double MultiminPhase::fun(const orsa::MultiminParameters * par) const {
  
  if (0) {
    // debug output 
    for (unsigned int k=0; k<par->size(); ++k) {
      ORSA_DEBUG("par[%02i] = [%20s] = %18.8Fg   step: %18.8Fg",
		 k,
		 par->name(k).c_str(),
		 par->get(k).get_mpf_t(),
		 par->getStep(k).get_mpf_t());
    }
  }
  
  Matrix m = Matrix::identity();
  m.rotZ(par->get("alpha"));
  
  const orsa::Double retVal = fabs(phi.getRef() - 
				   fabs(acos((m*uI.getRef()).normalized() * uS.getRef())-halfpi()));
  
  // ORSA_DEBUG("fun: %Fg",retVal.get_mpf_t());
  
  return retVal;
}

orsa::Double MultiminPhase::getAlpha(const orsa::Double & phaseAngle,
				     const orsa::Vector & uSun,
				     const orsa::Vector & uInclination) {
  
  phi = std::min(phaseAngle,fabs(twopi()-phaseAngle));
  uS  = uSun.normalized();
  uI  = uInclination.normalized();
  
  osg::ref_ptr<orsa::MultiminParameters> par = new orsa::MultiminParameters;
  //
  // par->insert("alpha",halfpi(),halfpi());
  // par->insert("alpha",phaseAngle,halfpi());
  par->insert("alpha",phaseAngle,0.01);
  //
  // test range
  // par->setRange("alpha",-pi(),pi());
  //
  setMultiminParameters(par.get());
  
  if (!run_nmsimplex(256,1.0e-14)) {
    ORSA_WARNING("the phase search did not converge.");
  }
  
  osg::ref_ptr<const orsa::MultiminParameters> parFinal = getMultiminParameters();
  //
  return (parFinal->get("alpha"));
} 
