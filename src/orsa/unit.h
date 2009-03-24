#ifndef _ORSA_UNIT_
#define _ORSA_UNIT_

#include <orsaTBB/malloc.h>

#include <orsa/cache.h>
#include <orsa/double.h>
#include <orsa/debug.h>

namespace orsa {
  
  //! All units in MKS
  
  class Unit {
    
  public:
    enum TimeUnit {
      YEAR,
      DAY,
      HOUR,
      MINUTE,
      SECOND,
      MILLISECOND,
      MICROSECOND
    };
  public:
    enum LengthUnit {
      MPARSEC,
      KPARSEC,
      PARSEC,
      LY,
      AU,
      KM,
      M,
      CM,
      // aliases
      METER=M
    };
  public:
    enum MassUnit {
      KG,
      GRAM
    };
    
  protected:
    static double getTimeScale(const TimeUnit &);
    static double getLengthScale(const LengthUnit &); 
    static double getMassScale(const MassUnit &); 
    
  public:
    inline static double FromUnits(const double   & x,
				   const TimeUnit & tu, 
				   const int      & power = 1) { 
      return (x*orsa::int_pow(getTimeScale(tu),power));
    }
  public:
    inline static double FromUnits(const double     & x,
				   const LengthUnit & lu,
				   const int        & power = 1) { 
      return (x*orsa::int_pow(getLengthScale(lu),power)); 
    }
  public:
    inline static double FromUnits(const double   & x, 
				   const MassUnit & mu, 
				   const int      & power = 1) {
      return (x*orsa::int_pow(getMassScale(mu),power));  
    }
    
  public:
    inline static const double & G() {
      static const double _G = 6.67259e-11;
      return _G;
    }
  public:
    inline static const double & c() {
      static const double _c = 299792458;
      return _c;
    }
    
  };
  
  inline double FromUnits(const double         & x, 
			  const Unit::TimeUnit & tu, 
			  const int            & power = 1) {
    return Unit::FromUnits(x, tu, power);
  }
  
  inline double FromUnits(const double           & x,
			  const Unit::LengthUnit & lu, 
			  const int              & power = 1) {
    return Unit::FromUnits(x, lu, power);
  }
  
  inline double FromUnits(const double         & x, 
			  const Unit::MassUnit & mu, 
			  const int            & power = 1) {
    return Unit::FromUnits(x, mu, power);
  }
  
} // namespace orsa

#endif // _ORSA_UNIT_
