#include <orsa/box.h>

#include <osg/BoundingBox>

using namespace orsa;

Box::Box() { }

Box::Box(const double & xMin,
         const double & xMax,
         const double & yMin,
         const double & yMax,
         const double & zMin,
         const double & zMax) {
    set(xMin, xMax, yMin, yMax, zMin, zMax);
}

Box::Box(const orsa::Vector & v1,
         const orsa::Vector & v2) {
    set(v1, v2);
}

void Box::set(const double & xMin,
              const double & xMax,
              const double & yMin,
              const double & yMax,
              const double & zMin,
              const double & zMax) {
    _xMin.set(xMin);
    _xMax.set(xMax);
    _yMin.set(yMin);
    _yMax.set(yMax);
    _zMin.set(zMin);
    _zMax.set(zMax);
}

void Box::set(const orsa::Vector & v1,
              const orsa::Vector & v2) {
    if (v1.getX() < v2.getX()) {
        _xMin.set(v1.getX());
        _xMax.set(v2.getX());
    } else {
        _xMin.set(v2.getX());
        _xMax.set(v1.getX());
    }
  
    if (v1.getY() < v2.getY()) {
        _yMin.set(v1.getY());
        _yMax.set(v2.getY());
    } else {
        _yMin.set(v2.getY());
        _yMax.set(v1.getY());
    }
  
    if (v1.getZ() < v2.getZ()) {
        _zMin.set(v1.getZ());
        _zMax.set(v2.getZ());
    } else {
        _zMin.set(v2.getZ());
        _zMax.set(v1.getZ());
    }
}

bool Box::isSet() const {
    return (_xMin.isSet() &&
            _xMax.isSet() &&
            _yMin.isSet() &&
            _yMax.isSet() &&
            _zMin.isSet() &&
            _zMax.isSet());
}

void Box::reset() {
    _xMin.reset();
    _xMax.reset();
    _yMin.reset();
    _yMax.reset();
    _zMin.reset();
    _zMax.reset();
}

double Box::volume() const {
    return fabs((_xMax.getRef()-_xMin.getRef())*
                (_yMax.getRef()-_yMin.getRef())*
                (_zMax.getRef()-_zMin.getRef()));
}

bool Box::isInside(const orsa::Vector & v) const {
    if (v.getX() < _xMin.getRef()) return false;
    if (v.getX() > _xMax.getRef()) return false;
    if (v.getY() < _yMin.getRef()) return false;
    if (v.getY() > _yMax.getRef()) return false;
    if (v.getZ() < _zMin.getRef()) return false;
    if (v.getZ() > _zMax.getRef()) return false;
    return true;
}

osg::BoundingBox Box::getOSGBoundingBox() const {
    return osg::BoundingBox(_xMin.getRef(),
                            _yMin.getRef(),
                            _zMin.getRef(),
                            _xMax.getRef(),
                            _yMax.getRef(),
                            _zMax.getRef());
}

