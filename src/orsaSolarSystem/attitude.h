#ifndef _ORSA_SOLAR_SYSTEM_ATTITUDE_
#define _ORSA_SOLAR_SYSTEM_ATTITUDE_

// #include <orsa/attitude.h>
#include <orsa/body.h>

namespace orsaSolarSystem {
  
  class ConstantZRotationEcliptic_RotationalBodyProperty : 
  public orsa::PrecomputedRotationalBodyProperty {
  
  public:
    ConstantZRotationEcliptic_RotationalBodyProperty(const orsa::Time   & t0,
						     const orsa::Double & phi0,
						     const orsa::Double & omega,
						     const orsa::Double & lambda,
						     const orsa::Double & beta);
    
  public:
    ConstantZRotationEcliptic_RotationalBodyProperty(const ConstantZRotationEcliptic_RotationalBodyProperty &);
    
  protected:
    ~ConstantZRotationEcliptic_RotationalBodyProperty() { 
      // ORSA_DEBUG("destroying %x",this);
    }
    
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
    orsa::RotationalBodyProperty * clone() const {
#warning "IMPORTANT: cloning prevents caching and optimizations... is cloning needed? where does it gets called?"	
      // ORSA_DEBUG("cloning %x",this);
      return new ConstantZRotationEcliptic_RotationalBodyProperty(*this);
    }
  public:
    bool update(const orsa::Time &);
  private:
    const orsa::Time   _t0;
    const orsa::Double _phi0;
    const orsa::Double _omega;
    const orsa::Double _lambda, _beta;    
  protected:
    orsa::Cache<orsa::Quaternion> _q;
    orsa::Cache<orsa::Vector>     _omegaVector;
  protected:
    orsa::Cache<orsa::Time> _previousTime;
  };
  
  //! this will be removed in future;
  /* class ConstantZRotationEclipticAttitude : public orsa::Attitude {
     public:	
     ConstantZRotationEclipticAttitude(const orsa::Time   & t0,
     const orsa::Double & phi0,
     const orsa::Double & omega,
     const orsa::Double & lambda,
     const orsa::Double & beta);
     
     protected:
     virtual ~ConstantZRotationEclipticAttitude();
     
     public:
     orsa::Matrix localToGlobal(const orsa::Time & t) const;
     orsa::Matrix globalToLocal(const orsa::Time & t) const;
     public:
     bool dynamic() const { return false; }
     
     private:
     const orsa::Time   _t0;
     const orsa::Double _phi0;
     const orsa::Double _omega;
     const orsa::Double _lambda, _beta;    
     };
  */
  
  class ConstantZRotationEquatorial_RotationalBodyProperty : 
  public orsa::PrecomputedRotationalBodyProperty {
  public:
    ConstantZRotationEquatorial_RotationalBodyProperty(const orsa::Time   & t0,
						       const orsa::Double & phi0,
						       const orsa::Double & omega,
						       const orsa::Double & alpha,
						       const orsa::Double & delta);
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
    orsa::RotationalBodyProperty * clone() const {
      return new ConstantZRotationEquatorial_RotationalBodyProperty(*this);
    }
  public:
    bool update(const orsa::Time &);
  private:
    const orsa::Time   _t0;
    const orsa::Double _phi0;
    const orsa::Double _omega;
    const orsa::Double _alpha, _delta;    
  protected:
    orsa::Cache<orsa::Quaternion> _q;
    orsa::Cache<orsa::Vector>     _omegaVector;
  protected:
    orsa::Cache<orsa::Time> _previousTime;
  };
  
  //! this will be removed in future;
  /* class ConstantZRotationEquatorialAttitude : public orsa::Attitude {
     public:	
     ConstantZRotationEquatorialAttitude(const orsa::Time   & t0,
     const orsa::Double & phi0,
     const orsa::Double & omega,
     const orsa::Double & alpha,
     const orsa::Double & delta);
     
     protected:
     virtual ~ConstantZRotationEquatorialAttitude();
     
     public:
     orsa::Matrix localToGlobal(const orsa::Time & t) const;
     orsa::Matrix globalToLocal(const orsa::Time & t) const;
     public:
     bool dynamic() const { return false; }
     
     private:
     const orsa::Time   _t0;
     const orsa::Double _phi0;
     const orsa::Double _omega;
     const orsa::Double _alpha, _delta;    
     };
  */
  
}; // namespace orsaSolarSystem

#endif // _ORSA_SOLAR_SYSTEM_ATTITUDE_
