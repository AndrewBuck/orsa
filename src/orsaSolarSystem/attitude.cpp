#include <orsaSolarSystem/attitude.h>

#include <orsaSolarSystem/obleq.h>

#include <orsa/print.h>
#include <orsa/util.h>

using namespace orsa;
using namespace orsaSolarSystem;

// ConstantZRotationEcliptic_RotationalBodyProperty

ConstantZRotationEcliptic_RotationalBodyProperty::ConstantZRotationEcliptic_RotationalBodyProperty (const Time   & t0,
                                                                                                    const double & phi0,
                                                                                                    const double & omega,
                                                                                                    const double & lambda,
                                                                                                    const double & beta) : 
    orsa::PrecomputedRotationalBodyProperty(), 
    _t0(t0), 
    _phi0(phi0),
    _omega(omega),
    _lambda(lambda),
    _beta(beta) { }

ConstantZRotationEcliptic_RotationalBodyProperty::ConstantZRotationEcliptic_RotationalBodyProperty(const ConstantZRotationEcliptic_RotationalBodyProperty & bp) : 
    orsa::PrecomputedRotationalBodyProperty(),
    _t0(bp._t0), 
    _phi0(bp._phi0),
    _omega(bp._omega),
    _lambda(bp._lambda),
    _beta(bp._beta),
    _q(bp._q),
    _omegaVector(bp._omegaVector),
    _previousTime(bp._previousTime) { }

bool ConstantZRotationEcliptic_RotationalBodyProperty::update(const orsa::Time & t) {
  
    if (_previousTime.isSet()) {
        if (_previousTime.getRef() == t) {
            // ORSA_DEBUG("cached...");
            return true;
        }    
    }
  
    _previousTime = t;
  
    const double _phi = fmod(_phi0 + _omega*(t-_t0).get_d(), orsa::twopi());
  
    Matrix _m = Matrix::identity();
  
    _m.rotZ(_phi);
  
    _m.rotY(halfpi()-_beta);
  
    _m.rotZ(_lambda);
  
    _q = MatrixToQuaternion(_m);
  
    _omegaVector = _omega * (_m*orsa::Vector(0,0,1)).normalized();
  
    return true;
}

// ConstantZRotationEquatorial_RotationalBodyProperty

ConstantZRotationEquatorial_RotationalBodyProperty::ConstantZRotationEquatorial_RotationalBodyProperty (const Time   & t0,
                                                                                                        const double & phi0,
                                                                                                        const double & omega,
                                                                                                        const double & alpha,
                                                                                                        const double & delta) : 
    orsa::PrecomputedRotationalBodyProperty(), 
    _t0(t0), 
    _phi0(phi0),
    _omega(omega),
    _alpha(alpha),
    _delta(delta) { }

bool ConstantZRotationEquatorial_RotationalBodyProperty::update(const orsa::Time & t) {
  
    if (_previousTime.isSet()) {
        if (_previousTime.getRef() == t) {
            // ORSA_DEBUG("cached...");
            return true;
        }    
    }
  
    _previousTime = t;
  
    const double _phi = fmod(_phi0 + _omega*(t-_t0).get_d(), orsa::twopi());
  
    Matrix _m = Matrix::identity();
  
    _m.rotZ(_phi);
  
    _m.rotY(halfpi()-_delta);
  
    _m.rotZ(_alpha);
  
    _m = orsaSolarSystem::equatorialToEcliptic()*_m;
  
    _q = MatrixToQuaternion(_m);
  
    _omegaVector = _omega * (_m*orsa::Vector(0,0,1)).normalized();
  
    return true;
}
