#ifndef _ORSA_SPICE_SPICE_BODY_POS_VEL_CALLBACK_
#define _ORSA_SPICE_SPICE_BODY_POS_VEL_CALLBACK_

#include <orsa/body.h>

namespace orsaSPICE {
  
  class SpiceBodyPosVelCallback : public orsa::PrecomputedTranslationalBodyProperty {
  public:    
    SpiceBodyPosVelCallback(const std::string & name);
  public:
    SpiceBodyPosVelCallback(const SpiceBodyPosVelCallback & sbpvc);
  public:
    orsa::Vector position() const;
    orsa::Vector velocity() const;
  public:
    orsa::TranslationalBodyProperty * clone() const {
      return new SpiceBodyPosVelCallback(*this);
    }
  public:
    bool update(const orsa::Time &);
  protected:
    const std::string _name;
  protected:
    orsa::Cache<orsa::Vector> _position, _velocity;
  protected:
    orsa::Cache<orsa::Time> _previousTime;
  };
  
}; // namespace orsaSPICE

#endif // _ORSA_SPICE_SPICE_BODY_POS_VEL_CALLBACK_
