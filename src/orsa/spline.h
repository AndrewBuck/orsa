#ifndef _ORSA_SPLINE_
#define _ORSA_SPLINE_

#include <osg/Referenced>

#include <orsa/datetime.h>
#include <orsa/double.h>
// #include <orsa/vector.h>

namespace orsa {
  
    //! Spline with exact position and velocity at initial and final time,
    //! and with acceleration = a_c + (1-2*tau) * a_p, tau = (t-t0)/(t1-t0),
    //! p position, v velocity at times t0 and t1
    template <class T> class PhysicalSpline : public osg::Referenced {
    public:
        PhysicalSpline() : Referenced(true), _set(false) { }
    protected:
        ~PhysicalSpline() { }
    
    public:
        bool set(const T & p0, const T & v0, const orsa::Time & t0, 
                 const T & p1, const T & v1, const orsa::Time & t1) {
            _p0 = p0;
            _v0 = v0;
            _t0 = t0;
            //
            _p1 = p1;
            _v1 = v1;
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
            _period     = (_t1 - _t0).get_d();
            _periodsq   = _period*_period;
            _halfperiod = _period*0.5;
            //
            _a_c =  (_v1 - _v0) / _period;
            _a_p = ((_p1 - _p0) - _halfperiod*(_v0 + _v1))*6/_periodsq;
      
            return valid();
        }
    
    public:
        bool valid() const {
            return (_set && (_t0 != _t1));
        }
    private:
        bool _set;
    
    public:
        virtual bool get(T & p, T & v, const orsa::Time & t) const {
      
            if (!valid()) {
                return false;
            }
      
            /* 
               ORSA_DEBUG("t.: %f",t.get_d());
               ORSA_DEBUG("t0: %f",_t0.get_d());
               ORSA_DEBUG("t1: %f",_t1.get_d());
            */
      
            // Add a check: t MUST be between _t0 and _t1 included
            if (t < _t0) {
                ORSA_ERROR("t < t0");
                ORSA_DEBUG("t.: %f",t.get_d());
                ORSA_DEBUG("t0: %f",_t0.get_d());
                ORSA_DEBUG("t1: %f",_t1.get_d());
                p = _p0;
                v = _v0;
                return false;
            }
            //
            if (t > _t1) {
                ORSA_ERROR("t > t1");
                ORSA_DEBUG("t.: %f",t.get_d());
                ORSA_DEBUG("t0: %f",_t0.get_d());
                ORSA_DEBUG("t1: %f",_t1.get_d());
                p = _p1;
                v = _v1;
                return false;
            }
      
            const double tau   = (t - _t0).get_d() / _period;
            const double tausq = tau*tau;
            const double taucb = tau*tausq;
      
            // ORSA_DEBUG("spline tau = %Fg",tau());
      
            p = _periodsq*(_a_c*tausq/2 + _a_p*(tausq/2 - taucb/3)) + _period*_v0*tau + _p0;
      
            v = _period*(_a_c*tau + _a_p*(tau-tausq)) + _v0;
      
            return true;
        }
    
    protected:    
        T _p0, _v0, _p1, _v1;
        orsa::Time _t0, _t1;

        // derived vars
    protected:
        double _period, _periodsq, _halfperiod;
        T _a_c, _a_p;
    };
  
} // namespace orsa

#endif // _ORSA_SPLINE_
