#include <orsaSPICE/spiceBodyRotationalCallback.h>

#include <orsa/util.h>

#include <orsaSPICE/spice.h>

using namespace orsa;
using namespace orsaSPICE;

SpiceBodyRotationalCallback::SpiceBodyRotationalCallback(const std::string & name) :
    orsa::PrecomputedRotationalBodyProperty(),
    _name(name) { }

SpiceBodyRotationalCallback::SpiceBodyRotationalCallback(const SpiceBodyRotationalCallback & sbrc) :
    orsa::PrecomputedRotationalBodyProperty(),
    _name(sbrc._name) {
    if (sbrc._previousTime.isSet()) {
        _q            = sbrc._q.getRef();
        _omega        = sbrc._omega.getRef();
        _previousTime = sbrc._previousTime.getRef();
    }
}

bool SpiceBodyRotationalCallback::get(orsa::Quaternion & q,
                                      orsa::Vector     & omega) const {
    q = _q.getRef();
    omega = _omega.getRef();
    return true;
}

orsa::Quaternion SpiceBodyRotationalCallback::getQ() const { return _q.getRef(); }

orsa::Vector SpiceBodyRotationalCallback::getOmega() const { return _omega.getRef(); }

bool SpiceBodyRotationalCallback::update(const orsa::Time & t) {
  
    if (_previousTime.isSet()) {
        if (_previousTime.getRef() == t) {
            // ORSA_DEBUG("cached...");
            return true;
        }    
    }
  
    _previousTime = t;
  
    const orsa::Matrix l2g = orsaSPICE::SPICE::instance()->localToGlobal(_name,t);
    _q = orsa::MatrixToQuaternion(l2g);
    // note: _omega is set to zero for now, since the rotation of the body is precomputed anyway
    // this is wrong in general, an a correct value for _omega should be computed
    _omega = orsa::Vector(0,0,0);
  
    return true;
}
