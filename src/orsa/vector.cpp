#include <cmath>
#include <limits>

#include <orsa/debug.h>
#include <orsa/vector.h>

using namespace orsa;

osg::Vec3d Vector::getVec3d() const {
  return osg::Vec3d(_x,
		    _y,
		    _z);
}

osg::Vec3f Vector::getVec3f() const {
  return osg::Vec3f(_x,
		    _y,
		    _z);
}

void Vector::_update_l() const {
  if (!_l2.isSet()) {
    _update_l2();
  }
  _l = sqrt(_l2.getRef());
}

void Vector::_update_l2() const {
  _l2 = (_x*_x) + (_y*_y) + (_z*_z);
}

const double & Vector::length() const {
  if (_l.isSet()) {
    return _l.getRef();
  } else {
    _update_l();
    return _l.getRef();
  }
}

const double & Vector::lengthSquared() const {
  if (_l2.isSet()) {
    return _l2.getRef();
  } else {
    _update_l2();
    return _l2.getRef();
  }
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
   return Vector(u._x+v._x,
   u._y+v._y,
   u._z+v._z);
   }
*/

/* 
   Vector operator - (const Vector & u, const Vector & v) {
   return Vector(u._x-v._x,
   u._y-v._y,
   u._z-v._z);
   }
*/

/* 
   Vector externalProduct (const Vector & lhs, const Vector & rhs) {
   return Vector (lhs._y*rhs._z-lhs._z*rhs._y,
   lhs._z*rhs._x-lhs._x*rhs._z,
   lhs._x*rhs._y-lhs._y*rhs._x); 
   }
*/

/* 
   double operator * (const Vector & u, const Vector & v) {
   return (u._x*v._x+
   u._y*v._y+
   u._z*v._z);
   }  
*/

bool orsa::operator == (const Vector & v1, const Vector & v2) {
  if (v1.getX() != v2.getX()) return false;
  if (v1.getY() != v2.getY()) return false;
  if (v1.getZ() != v2.getZ()) return false;
  return true;
}

