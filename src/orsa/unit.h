#ifndef _ORSA_UNIT_
#define _ORSA_UNIT_

#include <orsaTBB/malloc.h>

#include <orsa/double.h>

#include <orsa/debug.h>

namespace orsa {
  
  template <class UNIT> class UnitBaseScale {
  public:
    UnitBaseScale() { }
  public:
    UnitBaseScale(UNIT unit) { 
      _base_unit = unit;
    }
  public:
    void set(UNIT unit) { 
      _base_unit = unit;
    }
  public:
    const UNIT & get() const {
      return _base_unit;
    }
    
  private:
    UNIT _base_unit;
  };
  
  class Unit {
  public:
    static Unit * instance() {
      if (_instance == 0) {
	_instance = new Unit;
      }
      return _instance;
    }
  public:
    static bool instanciated() {
      return (_instance != 0);
    }
  protected:
    Unit();
  public:
    virtual ~Unit() {
      _instance = 0;
    }
  protected:
    static Unit * _instance;
  private:
    void _init();
    
  public:
    enum TimeUnit {
      YEAR=1,
      DAY=2,
      HOUR=3,
      MINUTE=4,
      SECOND=5,
      MILLISECOND=6,
      MICROSECOND=7
    };
  public:
    enum LengthUnit {
      MPARSEC=1,
      KPARSEC=2,
      PARSEC=3,
      LY=4,
      AU=5,
      EARTHMOON=6,
      REARTH=7,
      RMOON=8,
      KM=9,
      M=10,
      CM=11,
      // aliases
      LD=EARTHMOON,
      ER=REARTH,
      MR=RMOON,
      METER=M
    };
  public:
    enum MassUnit {
      MSUN=1,
      MJUPITER=2,
      MEARTH=3,
      MMOON=4,
      KG=5,
      GRAM=6
    };
    
  private:
    UnitBaseScale<TimeUnit>   _time;
    UnitBaseScale<LengthUnit> _length;
    UnitBaseScale<MassUnit>   _mass;
    
  protected:
    orsa::Double getTimeScale(const TimeUnit & tu) const;
    orsa::Double getLengthScale(const LengthUnit &) const; 
    orsa::Double getMassScale(const MassUnit &) const; 
    
  protected:
    orsa::Double getTimeScale() const;    
    orsa::Double getLengthScale() const; 
    orsa::Double getMassScale() const; 
    
  protected:
    void recompute();
    
  public:
    inline orsa::Double FromUnits(const orsa::Double & x,
				  const TimeUnit & tu, 
				  const mpz_class & power = mpz_class("1")) const { 
      return (x*int_pow(getTimeScale(tu)/getTimeScale(),power));
    }
  public:
    inline orsa::Double FromUnits(const orsa::Double & x,
				  const LengthUnit & lu,
				  const mpz_class & power = mpz_class("1")) const { 
      return (x*int_pow(getLengthScale(lu)/getLengthScale(),power)); 
    }
  public:
    inline orsa::Double FromUnits(const orsa::Double & x, 
				  const MassUnit & mu, 
				  const mpz_class & power = mpz_class("1")) const {
      
      /* 
	 ORSA_DEBUG("called FromUnits(%Fg,%Zi) [mass] [inner]",x.get_mpf_t(),power.get_mpz_t());
	 ORSA_DEBUG("getMassScale(mu): %Fg",getMassScale(mu).get_mpf_t());
	 ORSA_DEBUG("getMassScale(): %Fg",getMassScale().get_mpf_t());
	 const orsa::Double retVal = x*int_pow(getMassScale(mu)/getMassScale(),power);
	 ORSA_DEBUG("retVal: %Fg [inner]",retVal.get_mpf_t());
	 return retVal;
      */
      
      return (x*int_pow(getMassScale(mu)/getMassScale(),power));  
    }
    
  public:
    inline const TimeUnit & getTimeBaseUnit() const {
      return _time.get();
    };
  public:    
    inline const LengthUnit & getLengthBaseUnit() const { 
      return _length.get();
    };
  public:
    inline const MassUnit & getMassBaseUnit() const {
      return _mass.get();
    };
    
  public:
    inline orsa::Double getG()    const { return G; };
    inline orsa::Double getMSun() const { return MSun; };
    inline orsa::Double getC()    const { return c; }; // c = speed of light
    inline orsa::Double getC2()   const { return c2; }; // c2 = speed of light squared
  public:
    orsa::Double getG_MKS() const;
    
  private:
    orsa::Double G,G_base;
    orsa::Double MSun,MSun_base;
    orsa::Double MJupiter_base, MEarth_base, MMoon_base;
    orsa::Double AU_base;
    orsa::Double c,c2,c_base;
    orsa::Double r_earth_base;
    orsa::Double r_moon_base;
    orsa::Double parsec_base;
 
  private:
    const orsa::Double        G_MKS;
    const orsa::Double     MSUN_MKS;
    const orsa::Double MJUPITER_MKS;
    const orsa::Double   MEARTH_MKS;
    const orsa::Double    MMOON_MKS;
    const orsa::Double       AU_MKS;
    const orsa::Double        c_MKS;
    const orsa::Double  R_EARTH_MKS;
    const orsa::Double   R_MOON_MKS;
  };
  
  
  // _99%_ the user interface follows 
  inline orsa::Double FromUnits(const orsa::Double & x, const Unit::TimeUnit & tu, const mpz_class & power = mpz_class("1")) {
    return Unit::instance()->FromUnits(x, tu, power);
  }
  
  inline orsa::Double FromUnits(const orsa::Double & x, const Unit::LengthUnit & lu, const mpz_class & power = mpz_class("1")) {
    return Unit::instance()->FromUnits(x, lu, power);
  }
  
  inline orsa::Double FromUnits(const orsa::Double & x, const Unit::MassUnit & mu, const mpz_class & power = mpz_class("1")) {
    
    /* 
       ORSA_DEBUG("called FromUnits(%Fg,%Zi) [mass]",x.get_mpf_t(),power.get_mpz_t());
       const orsa::Double retVal = Unit::instance()->FromUnits(x, mu, power);
       ORSA_DEBUG("retVal: %Fg",retVal.get_mpf_t());
       return retVal;
    */
    
    return Unit::instance()->FromUnits(x, mu, power);
  }
  
} // namespace orsa

#endif // _ORSA_UNIT_
