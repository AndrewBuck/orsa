#ifndef _ORSA_BOX_
#define _ORSA_BOX_

#include <orsa/cache.h>
#include <orsa/double.h>
#include <orsa/vector.h>

// keep this in sync with include/osg/BoundingBox in OSG source
// this avoids including directly the OSG header, makes code lighter
#include <osg/Config>
namespace osg {
    template<typename VT> class BoundingBoxImpl;
    typedef BoundingBoxImpl<Vec3f> BoundingBoxf;
    typedef BoundingBoxImpl<Vec3d> BoundingBoxd;
#ifdef OSG_USE_FLOAT_BOUNDINGBOX
    typedef BoundingBoxf BoundingBox;
#else
    typedef BoundingBoxd BoundingBox;
#endif
}

namespace orsa {
  
    class Box {
    public:
        Box();
    public:
        Box(const double & xMin,
            const double & xMax,
            const double & yMin,
            const double & yMax,
            const double & zMin,
            const double & zMax);
    public:	    
        Box(const orsa::Vector & v1,
            const orsa::Vector & v2);
    
    public:
        void set(const double & xMin,
                 const double & xMax,
                 const double & yMin,
                 const double & yMax,
                 const double & zMin,
                 const double & zMax);
    public:
        void set(const orsa::Vector & v1,
                 const orsa::Vector & v2);
    public:
        bool isSet() const;
    public:
        void reset();
    
    public:
        const double & getXMin() const { return _xMin.getRef(); }
        const double & getXMax() const { return _xMax.getRef(); }
        const double & getYMin() const { return _yMin.getRef(); }
        const double & getYMax() const { return _yMax.getRef(); }
        const double & getZMin() const { return _zMin.getRef(); }
        const double & getZMax() const { return _zMax.getRef(); }
    
    public:
        osg::BoundingBox getOSGBoundingBox() const;
    
    public:
        double volume() const;
    
    public:
        bool isInside(const orsa::Vector &) const;
    
    protected:
        orsa::Cache<double> _xMin, _xMax, _yMin, _yMax, _zMin, _zMax;
    };
  
} // namespace orsa

#endif // _ORSA_BOX_
