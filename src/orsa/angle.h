#ifndef _ORSA_ANGLE_
#define _ORSA_ANGLE_

#include <orsa/debug.h>
#include <orsa/double.h>

namespace orsa {
  
    class Angle {
    public:
        Angle() { }
    public:
        Angle(const double & angle) : _rad(angle) { }
    public:
        Angle operator + () const { return Angle( _rad); }
        Angle operator - () const { return Angle(-_rad); }
    public:
        void   setRad(const double & angle) { _rad = angle; }
        double getRad() const { return _rad; }
        //
        void setDPS(const double & d, 
                    const double & p, 
                    const double & s, 
                    const int sign = 1);
    
        void getDPS(double & d, 
                    double & p, 
                    double & s, 
                    int & sign) const;
        //
        void setHMS(const double & h, 
                    const double & m, 
                    const double & s, 
                    const int sign = 1);
   
        void getHMS(double & h, 
                    double & m, 
                    double & s,
                    int & sign) const;
    
    private:
        int check_sign(int sign);
    
    private:
        double _rad;
    };
  
    inline double sin(const Angle & angle) {
        return ::sin(angle.getRad());
    }
  
    inline double cos(const Angle & angle) {
        return ::cos(angle.getRad());
    }
  
    inline double tan(const Angle & angle) {
        return ::tan(angle.getRad());
    }
  
    inline void sincos(const Angle & angle, double * s, double * c) {
#if defined(__APPLE__) || defined(__MINGW32__)
        (*s) = sin(angle.getRad());
        (*c) = cos(angle.getRad());
#else
        ::sincos(angle.getRad(),s,c); 
#endif
    }
  
} // namespace orsa

#endif // _ORSA_ANGLE_
