#ifndef _ORSA_SOLAR_SYSTEM_ATTITUDE_
#define _ORSA_SOLAR_SYSTEM_ATTITUDE_

// #include <orsa/attitude.h>
#include <orsa/body.h>
#include <orsa/print.h>

namespace orsaSolarSystem {
  
    class ConstantZRotationEcliptic_RotationalBodyProperty : 
        public orsa::PrecomputedRotationalBodyProperty {
  
    public:
        ConstantZRotationEcliptic_RotationalBodyProperty(const orsa::Time & t0,
                                                         const double     & phi0,
                                                         const double     & omega,
                                                         const double     & lambda,
                                                         const double     & beta);
    
    public:
        ConstantZRotationEcliptic_RotationalBodyProperty(const ConstantZRotationEcliptic_RotationalBodyProperty &);
    
    protected:
        ~ConstantZRotationEcliptic_RotationalBodyProperty() { }
    
    public:
        bool get(orsa::Quaternion & q,
                 orsa::Vector     & omega) const { 
            q     = _q.getRef();
            omega = _omegaVector.getRef();
            return true;
        }
    public:
        orsa::Quaternion getQ()    const  { return _q.getRef(); } 
        orsa::Vector     getOmega() const { return _omegaVector.getRef(); }
    public:
        bool update(const orsa::Time &);
    private:
        const orsa::Time   _t0;
        const double _phi0;
        const double _omega;
        const double _lambda, _beta;    
    protected:
        orsa::Cache<orsa::Quaternion> _q;
        orsa::Cache<orsa::Vector>     _omegaVector;
    protected:
        orsa::Cache<orsa::Time> _previousTime;
    };
  
    class ConstantZRotationEquatorial_RotationalBodyProperty : 
        public orsa::PrecomputedRotationalBodyProperty {
    public:
        ConstantZRotationEquatorial_RotationalBodyProperty(const orsa::Time   & t0,
                                                           const double & phi0,
                                                           const double & omega,
                                                           const double & alpha,
                                                           const double & delta);
    public:
        bool get(orsa::Quaternion & q,
                 orsa::Vector     & omega) const { 
            q     = _q.getRef();
            omega = _omegaVector.getRef();
            return true;
        }
    public:
        orsa::Quaternion getQ()    const  { return _q.getRef(); } 
        orsa::Vector     getOmega() const { return _omegaVector.getRef(); }
    public:
        bool update(const orsa::Time &);
    private:
        const orsa::Time   _t0;
        const double _phi0;
        const double _omega;
        const double _alpha, _delta;    
    protected:
        orsa::Cache<orsa::Quaternion> _q;
        orsa::Cache<orsa::Vector>     _omegaVector;
    protected:
        orsa::Cache<orsa::Time> _previousTime;
    };
  
}; // namespace orsaSolarSystem

#endif // _ORSA_SOLAR_SYSTEM_ATTITUDE_
