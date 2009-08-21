#ifndef _ORSA_SLERP_
#define _ORSA_SLERP_

#include <osg/Referenced>

#include <orsa/datetime.h>
#include <orsa/double.h>
#include <orsa/quaternion.h>

namespace orsa {
  
  //! SLERP = Spherical Linear intERPolation
  
  class Slerp : public osg::Referenced {
  public:
    Slerp() : Referenced(true), _set(false) { }
  protected:
    ~Slerp() { }
    
  public:
    inline bool set(const orsa::Quaternion & q0, const orsa::Time & t0, 
		    const orsa::Quaternion & q1, const orsa::Time & t1) {
      _q0 = q0;
      _t0 = t0;
      //
      _q1 = q1;
      _t1 = t1;
      //
      _set = true;
      
      /* 
	 ORSA_DEBUG("t0: %f",_t0.get_d());
	 ORSA_DEBUG("t1: %f",_t1.get_d());
      */
      
      if (!valid()) {
	return false;
      }
      
      // update derived vars
      _period = (_t1 - _t0).get_d();
      //
      if (fabs(_q0.getVector()*_q1.getVector()) > 1) { ORSA_ERROR("acos() argument larger than 1 or smaller than -1"); }
      _Omega = acos(_q0.getVector()*_q1.getVector());
      //
      // _oneOverSinOmega = 1/sin(_Omega);
      
      return valid();
    }
    
  public:
    inline bool valid() const {
      return (_set && (_t0 != _t1));
    }
  private:
    bool _set;
    
  public:
    inline virtual bool get(orsa::Quaternion & q, const orsa::Time & t) const {
      
      if (!valid()) {
	return false;
      }
      
      // Add a check: t MUST be between _t0 and _t1 included
      if (t < _t0) {
	ORSA_ERROR("t < t0");
	ORSA_DEBUG("t.: %f",t.get_d());
	ORSA_DEBUG("t0: %f",_t0.get_d());
	ORSA_DEBUG("t1: %f",_t1.get_d());
	q = _q0;
      	return false;
      }
      //
      if (t > _t1) {
	ORSA_ERROR("t > t1");
	ORSA_DEBUG("t.: %f",t.get_d());
	ORSA_DEBUG("t0: %f",_t0.get_d());
	ORSA_DEBUG("t1: %f",_t1.get_d());
	q = _q1;
	return false;
      }
      
      const double u = (t - _t0).get_d() / _period;
      
      // q = (sin((1-u)*_Omega)*_q0 + sin(u*_Omega)*_q1)*_oneOverSinOmega;
      //
      q = unitQuaternion(sin((1-u)*_Omega)*_q0 + sin(u*_Omega)*_q1);
      
      return true;
    }
    
  protected:    
    orsa::Quaternion _q0, _q1;
    orsa::Time _t0, _t1;
    
    // derived vars
  protected:
    double _period;
    double _Omega;
    // double _oneOverSinOmega;
  };
  
} // namespace orsa

#endif // _ORSA_SLERP_
