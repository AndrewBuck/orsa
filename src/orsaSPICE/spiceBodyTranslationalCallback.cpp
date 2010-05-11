#include <orsaSPICE/spiceBodyTranslationalCallback.h>

#include <orsaSPICE/spice.h>

using namespace orsa;
using namespace orsaSPICE;

SpiceBodyTranslationalCallback::SpiceBodyTranslationalCallback(const std::string & name) : 
    orsa::PrecomputedTranslationalBodyProperty(),
    _name(name) { }


SpiceBodyTranslationalCallback::SpiceBodyTranslationalCallback(const SpiceBodyTranslationalCallback & sbtc) :
    orsa::PrecomputedTranslationalBodyProperty(),
    _name(sbtc._name) { 
    if (sbtc._previousTime.isSet()) {
        _position     = sbtc._position.getRef();
        _velocity     = sbtc._velocity.getRef();
        _previousTime = sbtc._previousTime.getRef();
    }
}

orsa::Vector SpiceBodyTranslationalCallback::position() const { return _position.getRef(); }

orsa::Vector SpiceBodyTranslationalCallback::velocity() const { return _velocity.getRef(); }

bool SpiceBodyTranslationalCallback::update(const orsa::Time & t) {
  
    if (_previousTime.isSet()) {
        if (_previousTime.getRef() == t) {
            // ORSA_DEBUG("cached...");
            return true;
        }    
    }
  
    _previousTime = t;
  
    orsa::Vector r, v;
  
    SPICE::instance()->getPosVel(_name,
                                 t,
                                 r,
                                 v);
    _position = r;
    _velocity = v;
  
    return true;
}
