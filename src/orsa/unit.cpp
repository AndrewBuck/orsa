#include <orsa/unit.h>

using namespace orsa;

Unit * Unit::_instance = 0;

Unit::Unit() :
  // default constants
  G_MKS(6.67428e-11),
  MSUN_MKS(1.9884e30),
  MJUPITER_MKS(1.8985e27),
  MEARTH_MKS(5.97214e24),
  MMOON_MKS(7.34575e22),
  AU_MKS(1.49597870691e11),
  c_MKS(299792458),
  R_EARTH_MKS(6378136.6),
  R_MOON_MKS(1737400.0) {
  _init();
  _time.set(SECOND);
  _length.set(M);
  _mass.set(KG);
  recompute();
}

/* 
   Unit::Unit() :
   // default constants
   G_MKS(6.67259e-11),
   MSUN_MKS(1.9889e30),
   MJUPITER_MKS(1.8989e27),
   MEARTH_MKS(5.9742e24),
   MMOON_MKS(7.3483e22),
   AU_MKS(1.49597870660e11),
   c_MKS(299792458.0),
   R_EARTH_MKS(6378137.0),
   R_MOON_MKS(1737400.0) {
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

double Unit::getTimeScale(const TimeUnit & tu) const {
  // all in seconds
  switch(tu) {
  case YEAR:        return 31557600; break;
  case DAY:         return 86400;    break;
  case HOUR:        return 3600;     break;
  case MINUTE:      return 60;       break;
  case SECOND:      return 1;        break;
  case MILLISECOND: return 1.0e-3;   break;
  case MICROSECOND: return 1.0e-6;   break;
  };
  ORSA_ERROR("returning dummy value");
  return 1;
}

double Unit::getLengthScale(const LengthUnit & lu) const {
  // all in meters
  double ls = -1;
  switch(lu) {
  case   MPARSEC: ls = parsec_base*1.0e6; break;
  case   KPARSEC: ls = parsec_base*1.0e3; break;
  case    PARSEC: ls = parsec_base;       break;
  case        LY: ls = c_base*31557600;   break;
  case        AU: ls = AU_base;           break;
  case EARTHMOON: ls = 3.844e8;           break;
  case    REARTH: ls = r_earth_base;      break;
  case     RMOON: ls = r_moon_base;       break;
  case        KM: ls = 1000;              break;
  case         M: ls = 1;                 break;
  case        CM: ls = 0.01;              break;
  }
  if (ls < 0) {
    ORSA_ERROR("returning dummy value");
  }
  return ls;
}

double Unit::getMassScale(const MassUnit & mu) const {
  // all in kilo-grams
  double ms = -1;
  switch(mu) {
  case     MSUN: ms = MSun_base;     break;
  case MJUPITER: ms = MJupiter_base; break;
  case   MEARTH: ms = MEarth_base;   break;
  case    MMOON: ms = MMoon_base;    break;
  case       KG: ms = 1;             break;
  case     GRAM: ms = 0.001;         break;
  }
  if (ms < 0) {
    ORSA_ERROR("returning dummy value");
  }
  return ms;
}

double Unit::getTimeScale() const { 
  return Unit::getTimeScale(_time.get()); 
} 

double Unit::getLengthScale() const { 
  return Unit::getLengthScale(_length.get()); 
} 

double Unit::getMassScale() const { 
  return Unit::getMassScale(_mass.get());
} 

void Unit::recompute() { 
  G = G_base*int_pow(getTimeScale(),2)*int_pow(getLengthScale(),-3)*int_pow(getMassScale(),1); 
  MSun = MSun_base/getMassScale();
  c = c_base*getTimeScale()/getLengthScale();
  c2 = c*c;
  parsec_base = AU_base/(2*sin(pi()/(180*3600*2)));
}

double Unit::getG_MKS() const {
  return G_MKS; 
}

