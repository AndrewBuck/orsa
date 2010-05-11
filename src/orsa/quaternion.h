#ifndef _ORSA_QUATERNION_
#define _ORSA_QUATERNION_

#include <orsa/matrix.h>
#include <orsa/vector.h>

namespace orsa {
  
    class Quaternion {
    
    public:
        Quaternion() : _s(0), _v(0,0,0) { }
    public:
        Quaternion(const orsa::Vector & v) : _s(0), _v(v) { check(); }
    public:
        Quaternion(const double & s,
                   const orsa::Vector & v) : _s(s), _v(v) { check(); }
    public:
        virtual ~Quaternion() { }
    
    public:
        inline const double       & getScalar() const { return _s; }
        inline const orsa::Vector & getVector() const { return _v; }
    
    public:
        inline void setScalar(const double & s) { _s += s; check(); }
        inline void setVector(const orsa::Vector & v) { _v += v; check(); }
        inline void set(const double & s,
                        const orsa::Vector & v) { 
            _s = s;
            _v = v;
            //
            check();
        }
    
        // unary operators
    public:
        inline Quaternion & operator += (const Quaternion & q) {
            _s += q._s;
            _v += q._v;
            return *this;
        }
    public:
        inline Quaternion & operator -= (const Quaternion & q) {
            _s -= q._s;
            _v -= q._v;
            return *this;
        }
    public:
        inline Quaternion & operator *= (const double & f) {
            orsa::check(f);
            _s *= f;
            _v *= f;
            return *this;
        }
    public:
        inline Quaternion & operator /= (const double & f) {
            orsa::check(f);
            _s /= f;
            _v /= f;
            return *this;
        }
    
        // sign
    public:
        inline Quaternion operator + () const { return Quaternion( _s, _v); }
        inline Quaternion operator - () const { return Quaternion(-_s,-_v); }
    
        // binary operators
    public:
        Quaternion operator + (const Quaternion & rhs) const {
            return Quaternion(_s+rhs._s,
                              _v+rhs._v);
        }
    public:
        Quaternion operator - (const Quaternion & rhs) const {
            return Quaternion(_s-rhs._s,
                              _v-rhs._v);
        }
    
        // the only product defined for Quaternions
    public:
        Quaternion operator * (const Quaternion & rhs) const;
    
        // length
    public:
        double length() const {
            return sqrt(lengthSquared());
        }
    public:
        double lengthSquared() const {
            return (_s*_s + _v.lengthSquared());
        }
    
    protected:
        void check() const;
    
    protected:
        double       _s; // scalar component
        orsa::Vector _v; // vector component
    };
  
    inline Quaternion operator * (const double & f, const Quaternion & q) {
        Quaternion _q(q);
        _q *= f;
        return _q;
    }
  
    inline Quaternion operator * (const Quaternion & q, const double & f) {
        Quaternion _q(q);
        _q *= f;
        return _q;
    }
  
    inline Quaternion operator / (const Quaternion & q, const double & f) {
        Quaternion _q(q);
        _q /= f;
        return _q;
    }
  
    inline bool operator == (const Quaternion & q1, const Quaternion & q2) {
        if (q1.getScalar() != q2.getScalar()) return false;
        if (q1.getVector() != q2.getVector()) return false;
        return true;
    }
  
    inline bool operator != (const Quaternion & q1, const Quaternion & q2) {
        return !(q1 == q2);
    }
  
    inline Quaternion conjugate (const Quaternion & q) {
        return Quaternion(q.getScalar(), -q.getVector());
    }
  
    inline Quaternion inverse (const Quaternion & q) {
        return Quaternion(conjugate(q)/q.lengthSquared());
    }
  
    inline Quaternion unitQuaternion (const Quaternion & q) {
        return (q/q.length());
    }
  
} // namespace orsa

#endif // _ORSA_QUATERNION_
