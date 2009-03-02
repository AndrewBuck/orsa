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
  
  //! More information can be obtained here: http://www.hartrao.ac.za/nccsdoc/slalib/sun67.htx/node217.html
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
  
  double timeToJulian(const orsa::Time &);
  orsa::Time   julianToTime(const double &);
  
  orsa::Time J2000();
  
  orsa::Time now();
  
  double dayFraction(const orsa::Time &);
  
}; // namespace orsaSolarSystem

#endif // _ORSA_SOLAR_SYSTEM_DATETIME_
