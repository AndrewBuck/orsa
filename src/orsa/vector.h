#ifndef _ORSA_VECTOR_
#define _ORSA_VECTOR_

#include <orsa/cache.h>
#include <orsa/double.h>

namespace osg {
    class Vec3f;
    class Vec3d;
}

namespace orsa {
  
    class Vector {
    
        // constructors
    public:
        Vector() { }
    public:
        Vector(const Vector & v) : _x(v._x), _y(v._y), _z(v._z), _l(v._l), _l2(v._l2) { 
            check();
        }
    public:	
        Vector(const double & x, 
               const double & y, 
               const double & z) : _x(x), _y(y), _z(z) {
            check();
        }
    public:     
        ~Vector() { }
    
    public:
        // copy operator
        Vector & operator = (const Vector & v) {
            _x = v._x;
            _y = v._y;
            _z = v._z;
            //
            _l  = v._l;
            _l2 = v._l2; 
            return (*this);
        }
    
    public:
        inline const double & getX() const { return _x; }
        inline const double & getY() const { return _y; }
        inline const double & getZ() const { return _z; }
    public:    
        inline void setX(const double & x) { orsa::check(x); _x = x; _reset_cache(); }
        inline void setY(const double & y) { orsa::check(y); _y = y; _reset_cache(); }
        inline void setZ(const double & z) { orsa::check(z); _z = z; _reset_cache(); }
    public:
        osg::Vec3f getVec3f() const;
        osg::Vec3d getVec3d() const;
    
        // unary operators
        Vector & operator += (const Vector & v) {
            _x += v._x;
            _y += v._y;
            _z += v._z;
            _reset_cache(); 
            return *this;
        }
    
        Vector & operator -= (const Vector & v) {
            _x -= v._x;
            _y -= v._y;
            _z -= v._z;
            _reset_cache(); 
            return *this;
        }
    
        Vector & operator *= (const double & f) {
            orsa::check(f);
            _x *= f;
            _y *= f;
            _z *= f;
            _reset_cache(); 
            return *this;
        }
    
        Vector & operator /= (const double & f) {
            orsa::check(f);
            if (f==0) { 
                ORSA_DEBUG("catched attempt of divide by zero");
                orsa::crash();
            }
            _x /= f;
            _y /= f;
            _z /= f;
            _reset_cache(); 
            return *this;
        }
    
        // sign
        inline Vector operator + () const { return Vector( _x, _y, _z); }   
        inline Vector operator - () const { return Vector(-_x,-_y,-_z); }    
    
        // binary operators
        inline Vector operator + (const Vector & rhs) const {
            return Vector(_x+rhs._x,
                          _y+rhs._y,
                          _z+rhs._z);
        }
    
        inline Vector operator - (const Vector & rhs) const {
            // ORSA_DEBUG("Vector -> binary minus...");
            return Vector(_x-rhs._x,
                          _y-rhs._y,
                          _z-rhs._z);
        }
    
        // scalar product
        inline double operator * (const Vector & rhs) const {
            return double(_x*rhs._x+
                          _y*rhs._y+
                          _z*rhs._z);
        }
    
        inline void set(const double & x, const double & y, const double & z) {
            _x = x;
            _y = y;
            _z = z;
            //
            check();
            //
            _reset_cache(); 
        }
    
        // metrics
        const double & length() const;    
        const double & lengthSquared() const;
    
        double manhattanLength() const {
            return (fabs(_x)+fabs(_y)+fabs(_z));
        }
        //
        // double manhattanLength() const;    
    
        // normalization
        /* inline Vector normalized() const {
           double l = Length();
           if (l > (std::numeric_limits<double>::min() * 1.0e3))
           return Vector(_x/l, _y/l, _z/l);
           else
           return Vector(0.0, 0.0, 0.0);
           }
        */
        //
        Vector normalized() const;
    
        /* 
           inline Vector & normalize() {
           double l = Length();
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
        //
        Vector & normalize();
    
    protected:
        void check() const;
    
    protected:
        double _x, _y, _z;
    
        // cache
    private:
        // mutable Cache<double> _l, _l2, _ml;
        mutable Cache<double> _l, _l2;
    private:
        void _reset_cache() const {
            _l.reset();
            _l2.reset();
            // _ml.reset();
            // changed();
        }
    public:
        // inline virtual void changed() const { }
    private:
        void _update_l() const;
        void _update_l2() const;
    };
  
    inline Vector operator * (const double & f, const Vector & v) {
        Vector p(v);
        p *= f;
        return p;
    }
  
    inline Vector operator * (const Vector & v, const double & f) {
        orsa::check(f);
        Vector p(v);
        p *= f;
        return p;
    }
  
    inline Vector operator / (const Vector & v, const double & f) {
        orsa::check(f);
        Vector p(v);
        p /= f;
        return p;
    }
  
    Vector externalProduct (const Vector & lhs, const Vector & rhs);
  
    bool operator == (const Vector &, const Vector &);
  
    inline bool operator != (const Vector & v1, const Vector & v2) {
        return !(v1 == v2);
    }
  
} // namespace orsa

#endif // _ORSA_VECTOR_
