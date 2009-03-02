#include <cmath>
#include <limits>

#include <orsa/debug.h>
#include <orsa/vector.h>

using namespace orsa;

osg::Vec3d Vector::getVec3d() const {
  return osg::Vec3d(getX(),
		    getY(),
		    getZ());
}

osg::Vec3f Vector::getVec3f() const {
  return osg::Vec3f(getX(),
		    getY(),
		    getZ());
}

void Vector::_update_l() const {
  if (!_l2.isSet()) {
    _update_l2();
  }
  _l.set(sqrt(_l2.getRef()));
}

void Vector::_update_l2() const {
  _l2.set((_x*_x) +
	  (_y*_y) +
	  (_z*_z) );
}

const double & Vector::length() const {
  if (_l.isSet()) {
    return _l.getRef();
  }
  _update_l();
  return _l.getRef();
}

const double & Vector::lengthSquared() const {
  if (_l2.isSet()) {
    return _l2.getRef();
  }
  _update_l2();
  return _l2.getRef();
}

bool Vector::isZero() const {
  return (lengthSquared() < (epsilon()*epsilon()));
}

// normalization
Vector Vector::normalized() const {
  if (isZero()) {
    return Vector(0,0,0);
  } else {
    const double _one_over_l = 1/length();
    return Vector(_x*_one_over_l,
		  _y*_one_over_l,
		  _z*_one_over_l);
  }  
}

/* 
   Vector & Vector::normalize() {
   ORSA_DEBUG("rewrite Vector::normalize() for double!!");
   double l = length();
   if (l > (std::numeric_limits<double>::min() * 1.0e3)) {
   _x /= l;
   _y /= l;
   _z /= l;
   } else {
   _z = 0.0;
   _y = 0.0;
   _z = 0.0;
   }
   _reset_cache(); 
   return *this;
   }
*/

Vector & Vector::normalize() {
  if (isZero()) {
    _x = _y = _z = 0;
  } else {
    const double _one_over_l = 1/length();
    _x *= _one_over_l;
    _y *= _one_over_l;
    _z *= _one_over_l;
  }
  _reset_cache(); 
  return *this;
}

/* 
   Vector operator + (const Vector & u, const Vector & v) {
   return Vector(u.getX()+v.getX(),
   u.getY()+v.getY(),
   u.getZ()+v.getZ());
   }
*/

/* 
   Vector operator - (const Vector & u, const Vector & v) {
   return Vector(u.getX()-v.getX(),
   u.getY()-v.getY(),
   u.getZ()-v.getZ());
   }
*/

/* 
   Vector externalProduct (const Vector & lhs, const Vector & rhs) {
   return Vector (lhs.getY()*rhs.getZ()-lhs.getZ()*rhs.getY(),
   lhs.getZ()*rhs.getX()-lhs.getX()*rhs.getZ(),
   lhs.getX()*rhs.getY()-lhs.getY()*rhs.getX()); 
   }
*/

/* 
   double operator * (const Vector & u, const Vector & v) {
   return (u.getX()*v.getX()+
   u.getY()*v.getY()+
   u.getZ()*v.getZ());
   }  
*/

bool orsa::operator == (const Vector & v1, const Vector & v2) {
  if (v1.getX() != v2.getX()) return false;
  if (v1.getY() != v2.getY()) return false;
  if (v1.getZ() != v2.getZ()) return false;
  return true;
}

