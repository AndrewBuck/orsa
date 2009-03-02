#include <orsa/unit.h>

using namespace orsa;

Unit * Unit::_instance = 0;

Unit::Unit() :
  // default constants
  G_MKS("6.67428e-11"),
  MSUN_MKS("1.9884e30"),
  MJUPITER_MKS("1.8985e27"),
  MEARTH_MKS("5.97214e24"),
  MMOON_MKS("7.34575e22"),
  AU_MKS("1.49597870691e11"),
  c_MKS("299792458"),
  R_EARTH_MKS("6378136.6"),
  R_MOON_MKS("1737400.0") {
  _init();
  _time.set(SECOND);
  _length.set(M);
  _mass.set(KG);
  recompute();
}

/* 
   Unit::Unit() :
   // default constants
   G_MKS("6.67259e-11"),
   MSUN_MKS("1.9889e30"),
   MJUPITER_MKS("1.8989e27"),
   MEARTH_MKS("5.9742e24"),
   MMOON_MKS("7.3483e22"),
   AU_MKS("1.49597870660e11"),
   c_MKS("299792458.0"),
   R_EARTH_MKS("6378137.0"),
   R_MOON_MKS("1737400.0") {
   _init();
   _time.set(DAY);
   _length.set(AU);
   _mass.set(MSUN);
   recompute();
   }
*/

void Unit::_init() {
  G_base        =        G_MKS;
  MSun_base     =     MSUN_MKS;
  MJupiter_base = MJUPITER_MKS;
  MEarth_base   =   MEARTH_MKS;
  MMoon_base    =    MMOON_MKS;
  AU_base       =       AU_MKS;
  c_base        =        c_MKS;
  r_earth_base  =  R_EARTH_MKS;
  r_moon_base   =   R_MOON_MKS;
}

Double Unit::getTimeScale(const TimeUnit & tu) const {
  // all in seconds
  switch(tu) {
  case YEAR:        return Double("31557600.0"); break;
  case DAY:         return Double("86400.0");    break;
  case HOUR:        return Double("3600.0");     break;
  case MINUTE:      return Double("60.0");       break;
  case SECOND:      return Double("1.0");        break;
  case MILLISECOND: return Double("1.0e-3");     break;
  case MICROSECOND: return Double("1.0e-6");     break;
  };
  ORSA_ERROR("returning dummy value");
  return Double("1.0");
}

Double Unit::getLengthScale(const LengthUnit & lu) const {
  // all in meters
  Double ls("-1.0");
  switch(lu) {
  case   MPARSEC: ls = parsec_base*Double("1.0e6"); break;
  case   KPARSEC: ls = parsec_base*Double("1.0e3"); break;
  case    PARSEC: ls = parsec_base;                 break;
  case        LY: ls = c_base*Double("31557600.0"); break;
  case        AU: ls = AU_base;                     break;
  case EARTHMOON: ls = Double("3.844e8");           break;
  case    REARTH: ls = r_earth_base;                break;
  case     RMOON: ls = r_moon_base;                 break;
  case        KM: ls = Double("1000.0");            break;
  case         M: ls = Double("1.0");               break;
  case        CM: ls = Double("0.01");              break;
  }
  if (ls < Double("0.0")) {
    ORSA_ERROR("returning dummy value");
  }
  return ls;
}

Double Unit::getMassScale(const MassUnit & mu) const {
  // all in kilo-grams
  Double ms = Double("-1.0");
  switch(mu) {
  case     MSUN: ms = MSun_base;       break;
  case MJUPITER: ms = MJupiter_base;   break;
  case   MEARTH: ms = MEarth_base;     break;
  case    MMOON: ms = MMoon_base;      break;
  case       KG: ms = Double("1.0");   break;
  case     GRAM: ms = Double("0.001"); break;
  }
  if (ms < Double("0.0")) {
    ORSA_ERROR("returning dummy value");
  }
  return ms;
}

Double Unit::getTimeScale() const { 
  return Unit::getTimeScale(_time.get()); 
} 

Double Unit::getLengthScale() const { 
  return Unit::getLengthScale(_length.get()); 
} 

Double Unit::getMassScale() const { 
  return Unit::getMassScale(_mass.get());
} 

void Unit::recompute() { 
  G = G_base*int_pow(getTimeScale(),2)*int_pow(getLengthScale(),-3)*int_pow(getMassScale(),1); 
  MSun = MSun_base/getMassScale();
  c = c_base*getTimeScale()/getLengthScale();
  c2 = c*c;
  parsec_base = AU_base/(Double("2.0")*sin((pi()/Double("180.0"))/Double("3600.0")/Double("2.0")));
}

Double Unit::getG_MKS() const {
  return G_MKS; 
}

