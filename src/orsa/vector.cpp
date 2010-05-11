#include <cmath>
#include <limits>

#include <orsa/debug.h>
#include <orsa/vector.h>

#include <osg/Vec3f>
#include <osg/Vec3d>

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

// normalization
Vector Vector::normalized() const {
    if (length() > 0) {
        const double _one_over_l = 1/length();
        return Vector(_x*_one_over_l,
                      _y*_one_over_l,
                      _z*_one_over_l);
    } else {
        ORSA_DEBUG("cannot normalize zero vector, time to die!");
        orsa::crash();
        // placeholder, to keep compiler happy
        return (*this);
    }
}

Vector & Vector::normalize() {
    if (length() > 0) {
        const double _one_over_l = 1/length();
        _x *= _one_over_l;
        _y *= _one_over_l;
        _z *= _one_over_l;
        check();
    } else {
        orsa::crash();
    } 
    _reset_cache(); 
    return *this;
}

void Vector::check() const {
    orsa::check(_x);
    orsa::check(_y);
    orsa::check(_z);
}

bool orsa::operator == (const Vector & v1, const Vector & v2) {
    if (v1.getX() != v2.getX()) return false;
    if (v1.getY() != v2.getY()) return false;
    if (v1.getZ() != v2.getZ()) return false;
    return true;
}

orsa::Vector orsa::externalProduct (const orsa::Vector & lhs, const orsa::Vector & rhs) {
    return orsa::Vector (lhs.getY()*rhs.getZ()-lhs.getZ()*rhs.getY(),
                         lhs.getZ()*rhs.getX()-lhs.getX()*rhs.getZ(),
                         lhs.getX()*rhs.getY()-lhs.getY()*rhs.getX()); 
}
