#ifndef _ORSA_QUATERNION_
#define _ORSA_QUATERNION_

#include <orsa/matrix.h>
#include <orsa/vector.h>

namespace orsa {
  
  class Quaternion {
    
  public:
    Quaternion() : _s(zero()), _v(zero(),zero(),zero()) { }
  public:
    Quaternion(const orsa::Vector & v) : _s(zero()), _v(v) { }
  public:
    Quaternion(const orsa::Double & s,
	       const orsa::Vector & v) : _s(s), _v(v) { }
  public:
    virtual ~Quaternion() { }
    
  public:
    inline const orsa::Double & getScalar() const { return _s; }
    inline const orsa::Vector & getVector() const { return _v; }
    
  public:
    inline void setScalar(const orsa::Double & s) { _s += s; }
    inline void setVector(const orsa::Vector & v) { _v += v; }
    inline void set(const orsa::Double & s,
		    const orsa::Vector & v) { 
      _s = s;
      _v = v;
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
    inline Quaternion & operator *= (const orsa::Double & f) {
      _s *= f;
      _v *= f;
      return *this;
    }
  public:
    inline Quaternion & operator /= (const orsa::Double & f) {
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
    const orsa::Double length() const {
      return sqrt(lengthSquared());
    }
  public:
    const orsa::Double lengthSquared() const {
      return (_s*_s + _v.lengthSquared());
    }
    
  protected:
    orsa::Double _s; // scalar component
    orsa::Vector _v; // vector component
  };
  
  inline Quaternion operator * (const orsa::Double & f, const Quaternion & q) {
    Quaternion _q(q);
    _q *= f;
    return _q;
  }
  
  inline Quaternion operator * (const Quaternion & q, const orsa::Double & f) {
    Quaternion _q(q);
    _q *= f;
    return _q;
  }
  
  inline Quaternion operator / (const Quaternion & q, const orsa::Double & f) {
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
