#ifndef _ORSA_SOLAR_SYSTEM_DATETIME_
#define _ORSA_SOLAR_SYSTEM_DATETIME_

#include <orsa/datetime.h>
#include <orsa/double.h>

namespace orsaSolarSystem {
  
    /*! A note on Gregorian and Julian calendars:
     * since its introduction in 1582, the Gregorian calendar
     * has been adopted in different years from different countries,
     * so the date obtained from the Date class can be different from 
     * the real 'old' one which used the Julian date, and usually
     * the difference is of a few days.
     *
     * In particular, ORSA applies the Gregorian calendar in any epoch,
     * i.e. even before 1582. This appears to be the simplest solution,
     * at least for the moment.
     * 
     * For more info, check out i.e. http://www.dome-igm.com/convers.htm
     */

    /* From JPL's Horizons system:
     *
     *  The three time systems are described as follows: 
     *
     *   CT ("Coordinate Time"); typically for cartesian and osculating element 
     *       tables. The uniform time scale and independent variable of the 
     *       ephemerides.
     *
     *   TT  ("Terrestrial (Dynamic) Time"), called TDT prior to 1991, used for
     *       observer quantity tables. This is proper time as measured by an 
     *       Earth-bound observer and is directly related to atomic time, TAI.
     *       TT periodically differs from CT by, at most, 0.002 seconds.
     *
     *   UT  is Universal Time. This can mean one of two non-uniform time-scales 
     *       based on the rotation of the Earth. For this program, prior to 1962, 
     *       UT means UT1.  After 1962, UT means UTC or "Coordinated Universal 
     *       Time". Future UTC leap-seconds are not known yet, so the closest 
     *       known leap-second correction is used over future time-spans.
     *
     */
  
    /*
     * Also read: http://www.cv.nrao.edu/~rfisher/Ephemerides/times.html
     */
  
    enum TimeScale {
        TS_UTC = 1,
        TS_UT  = 2,
        TS_TAI = 3,
        TS_TDT = 4,
        TS_GPS = 5,
        // aliases
        TS_UT1 = TS_UT,
        TS_ET  = TS_TDT,
        TS_TT  = TS_TDT
    };
  
    orsa::Time FromTimeScale(const orsa::Time & t,
                             const orsaSolarSystem::TimeScale ts);
  
    orsa::Time ToTimeScale(const orsa::Time & t,
                           const orsaSolarSystem::TimeScale ts);
  
    // conversion
    inline double MJD2JD(const double & MJD) { return (MJD+2400000.5); }
    inline double JD2MJD(const double &  JD) { return ( JD-2400000.5); }
  
    orsa::Time gregorTime(int y, 
                          int m, 
                          int d, 
                          int H, 
                          int M, 
                          int S, 
                          int ms);
  
    orsa::Time gregorTime(int y, 
                          int m, 
                          double d);
  
    void gregorDay(const orsa::Time & t,
                   int & y,
                   int & m,
                   int & d,
                   int & H,
                   int & M,
                   int & S,
                   int & ms);
  
    void gregorDay(const orsa::Time & t,
                   int & y,
                   int & m,
                   int & d,
                   double & fd);
  
    double     timeToJulian(const orsa::Time &);
    orsa::Time julianToTime(const double &);
  
    orsa::Time J2000();
  
    orsa::Time now();
  
    //! approximated period between two full moons
    //! real full moon times can be off by up to about 15 hours
    //! only really useful to divide time in lunation periods of (almost) constant length (very small quadratic term)
    void ApproximatedLunation(const orsa::Time & t,
                              orsa::Time & begin,
                              orsa::Time & end,
                              int & lunationID);
  
    // utility: returns the "fractional year" i.e. 2009.121233, useful for plots...
    double fractionalYear(const orsa::Time & t);
  
}; // namespace orsaSolarSystem

#endif // _ORSA_SOLAR_SYSTEM_DATETIME_
