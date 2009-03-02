#include <orsaSolarSystem/attitude.h>

#include <orsaSolarSystem/obleq.h>

#include <orsa/print.h>
#include <orsa/util.h>

using namespace orsa;
using namespace orsaSolarSystem;

// ConstantZRotationEcliptic_RotationalBodyProperty

ConstantZRotationEcliptic_RotationalBodyProperty::ConstantZRotationEcliptic_RotationalBodyProperty (const Time   & t0,
												    const Double & phi0,
												    const Double & omega,
												    const Double & lambda,
												    const Double & beta) : 
  orsa::PrecomputedRotationalBodyProperty(), 
  _t0(t0), 
  _phi0(phi0),
  _omega(omega),
  _lambda(lambda),
  _beta(beta) { 
  // ORSA_DEBUG("NOTE: check lambda and beta definitions!");
  // ORSA_DEBUG("constructing %x",this);
}

ConstantZRotationEcliptic_RotationalBodyProperty::ConstantZRotationEcliptic_RotationalBodyProperty(const ConstantZRotationEcliptic_RotationalBodyProperty & bp) : 
  orsa::PrecomputedRotationalBodyProperty(),
  _t0(bp._t0), 
  _phi0(bp._phi0),
  _omega(bp._omega),
  _lambda(bp._lambda),
  _beta(bp._beta),
  _q(bp._q),
  _omegaVector(bp._omegaVector),
  _previousTime(bp._previousTime) {
  // ORSA_DEBUG("copy-constructor: %x copy of %x",this,&bp);
  // ORSA_DEBUG("orig: %i   copy: %i",bp._previousTime.isSet(),_previousTime.isSet());
}

bool ConstantZRotationEcliptic_RotationalBodyProperty::update(const orsa::Time & t) {
  
  /* 
     if (_previousTime.isSet()) {
     ORSA_DEBUG("this: %x   t: %Zi   prev.t: %Zi !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",
     this,
     t.getMuSec().get_mpz_t(),
     _previousTime.getRef().getMuSec().get_mpz_t());
     } else {
     ORSA_DEBUG("this: %x   t: %Zi   prev.t: not set",
     this,
     t.getMuSec().get_mpz_t());
     }
  */
  
  if (_previousTime.isSet()) {
    if (_previousTime.getRef() == t) {
      ORSA_DEBUG("cached...");
      return true;
    }    
  }
  
  _previousTime = t;
  
  // const Double _phi = _phi0 + _omega*(t-_t0).asDouble();
  // 
  const Double _phi = fmod(_phi0 + _omega*(t-_t0).asDouble(), orsa::twopi());
  
  Matrix _m = Matrix::identity();
  
  _m.rotZ(_phi);
  
  _m.rotY(halfpi()-_beta);
  
  _m.rotZ(_lambda);
  
  _q = MatrixToQuaternion(_m);
  
  // debug
  // orsa::print(_q.getRef());
  
#warning "test this definition for omegaVector"
  _omegaVector = _omega * (_m*orsa::Vector(0,0,1)).normalized();
  
  // ORSA_DEBUG("following: omegaVector...");
  // print(_omegaVector.getRef());
  
  /* 
     ORSA_DEBUG("phi: %Fg   t: %Zi   prev.t: %Zi",
     _phi.get_mpf_t(),
     t.getMuSec().get_mpz_t(),
     _previousTime.getRef().getMuSec().get_mpz_t());
  */
  
  return true;
}

// ConstantZRotationEquatorial_RotationalBodyProperty

ConstantZRotationEquatorial_RotationalBodyProperty::ConstantZRotationEquatorial_RotationalBodyProperty (const Time   & t0,
													const Double & phi0,
													const Double & omega,
													const Double & alpha,
													const Double & delta) : 
  orsa::PrecomputedRotationalBodyProperty(), 
  _t0(t0), 
  _phi0(phi0),
  _omega(omega),
  _alpha(alpha),
  _delta(delta) { 
  // ORSA_DEBUG("NOTE: check alpha and delta definitions!");
}

bool ConstantZRotationEquatorial_RotationalBodyProperty::update(const orsa::Time & t) {
  
  if (_previousTime.isSet()) {
    if (_previousTime.getRef() == t) {
      // ORSA_DEBUG("cached...");
      return true;
    }    
  }
  
  _previousTime = t;
  
  const Double _phi = _phi0 + _omega*(t-_t0).asDouble();
  
  Matrix _m = Matrix::identity();
  
  _m.rotZ(_phi);
  
  _m.rotY(halfpi()-_delta);
  
  _m.rotZ(_alpha);
  
#warning "check this!!"
  // _m.rotX(-orsaSolarSystem::obleqJ2000());
  _m = orsaSolarSystem::equatorialToEcliptic()*_m;
  
  _q = MatrixToQuaternion(_m);
  
#warning "test this definition for omegaVector"
  _omegaVector = _omega * (_m*orsa::Vector(0,0,1)).normalized();
  
  return true;
}
