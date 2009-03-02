#ifndef _ORSA_BOX_
#define _ORSA_BOX_

#include <orsa/cache.h>
#include <orsa/double.h>
#include <orsa/vector.h>

#include <osg/BoundingBox>

namespace orsa {
  
  class Box {
  public:
    Box();
  public:
    Box(const orsa::Double & xMin,
	const orsa::Double & xMax,
	const orsa::Double & yMin,
	const orsa::Double & yMax,
	const orsa::Double & zMin,
	const orsa::Double & zMax);
  public:	    
    Box(const orsa::Vector & v1,
	const orsa::Vector & v2);
    
  public:
    void set(const orsa::Double & xMin,
	     const orsa::Double & xMax,
	     const orsa::Double & yMin,
	     const orsa::Double & yMax,
	     const orsa::Double & zMin,
	     const orsa::Double & zMax);
  public:
    void set(const orsa::Vector & v1,
	     const orsa::Vector & v2);
  public:
    bool isSet() const;
  public:
    void reset();
    
  public:
    const orsa::Double & getXMin() const { return _xMin.getRef(); }
    const orsa::Double & getXMax() const { return _xMax.getRef(); }
    const orsa::Double & getYMin() const { return _yMin.getRef(); }
    const orsa::Double & getYMax() const { return _yMax.getRef(); }
    const orsa::Double & getZMin() const { return _zMin.getRef(); }
    const orsa::Double & getZMax() const { return _zMax.getRef(); }

  public:
    osg::BoundingBox getOSGBoundingBox() const {
      return osg::BoundingBox(_xMin.getRef().get_d(),
			      _yMin.getRef().get_d(),
			      _zMin.getRef().get_d(),
			      _xMax.getRef().get_d(),
			      _yMax.getRef().get_d(),
			      _zMax.getRef().get_d());
    }
    
  public:
    Double volume() const;
    
  public:
    bool isInside(const orsa::Vector &) const;
    
  protected:
    orsa::Cache<orsa::Double> _xMin, _xMax, _yMin, _yMax, _zMin, _zMax;
  };
  
} // namespace orsa

#endif // _ORSA_BOX_
