#ifndef _ORSA_SPICE_BODY_ROTATIONAL_CALLBACK_
#define _ORSA_SPICE_BODY_ROTATIONAL_CALLBACK_

#include <orsa/body.h>

namespace orsaSPICE {
  
    class SpiceBodyRotationalCallback : public orsa::PrecomputedRotationalBodyProperty {
    public:
        SpiceBodyRotationalCallback(const std::string & name);
    public:
        SpiceBodyRotationalCallback(const SpiceBodyRotationalCallback & sbrc);
    public:
        bool get(orsa::Quaternion & q,
                 orsa::Vector     & omega) const;
    public:
        orsa::Quaternion getQ() const;
        orsa::Vector     getOmega() const;
    public:
        bool update(const orsa::Time &);
    protected:
        const std::string _name;
    protected:
        orsa::Cache<orsa::Quaternion> _q;
        orsa::Cache<orsa::Vector> _omega;
    protected:
        orsa::Cache<orsa::Time> _previousTime;
    };
  
}; // namespace orsaSPICE

#endif // _ORSA_SPICE_BODY_ROTATIONAL_CALLBACK_
