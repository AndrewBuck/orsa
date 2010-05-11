#ifndef _ORSA_SPICE_SPICE_BODY_TRANSLATIONAL_CALLBACK_
#define _ORSA_SPICE_SPICE_BODY_TRANSLATIONAL_CALLBACK_

#include <orsa/body.h>

namespace orsaSPICE {
  
    class SpiceBodyTranslationalCallback : public orsa::PrecomputedTranslationalBodyProperty {
    public:    
        SpiceBodyTranslationalCallback(const std::string & name);
    public:
        SpiceBodyTranslationalCallback(const SpiceBodyTranslationalCallback & sbtc);
    public:
        orsa::Vector position() const;
        orsa::Vector velocity() const;
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

#endif // _ORSA_SPICE_SPICE_BODY_TRANSLATIONAL_CALLBACK_
