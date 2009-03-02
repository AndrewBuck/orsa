#include <orsaSPICE/spiceBodyPosVelCallback.h>

#include <orsaSPICE/spice.h>

using namespace orsa;
using namespace orsaSPICE;

SpiceBodyPosVelCallback::SpiceBodyPosVelCallback(const std::string & name) : 
  orsa::PrecomputedTranslationalBodyProperty(),
  _name(name) { }


SpiceBodyPosVelCallback::SpiceBodyPosVelCallback(const SpiceBodyPosVelCallback & sbpvc) :
  orsa::PrecomputedTranslationalBodyProperty(),
  _name(sbpvc._name) { 
  // copying these members is necessary, otherwise the 'cache' never gets really used
  if (sbpvc._previousTime.isSet()) {
    _position     = sbpvc._position.getRef();
    _velocity     = sbpvc._velocity.getRef();
    _previousTime = sbpvc._previousTime.getRef();
  }
}

orsa::Vector SpiceBodyPosVelCallback::position() const { return _position.getRef(); }

orsa::Vector SpiceBodyPosVelCallback::velocity() const { return _velocity.getRef(); }

bool SpiceBodyPosVelCallback::update(const orsa::Time & t) {
  
  /* 
     ORSA_DEBUG("--MARK--  pt.set: %i  this: 0x%x  t: %Zi  [%s]",
     _previousTime.isSet(),
     this,
     t.getMuSec().get_mpz_t(),
     _name.c_str());
  */
  
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
  
  /* 
     ORSA_DEBUG("_position.length(): %.20e   [%s]",
     _position.getRef().length(),
     _name.c_str());
  */
  
  return true;
}
