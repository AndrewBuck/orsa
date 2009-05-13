#include <orsa/unit.h>

using namespace orsa;

// Unit * Unit::_instance = 0;

double Unit::getTimeScale(const TimeUnit & tu) {
  switch(tu) {
  case YEAR:        return 31557600; break;
  case DAY:         return 86400;    break;
  case HOUR:        return 3600;     break;
  case MINUTE:      return 60;       break;
  case SECOND:      return 1;        break;
  case MILLISECOND: return 1.0e-3;   break;
  case MICROSECOND: return 1.0e-6;   break;
  }
  ORSA_DEBUG("problems...");
  return 1;
}

double Unit::getLengthScale(const LengthUnit & lu) {
  static const double _AU     = 1.49597870691e11;
  static const double _parsec = _AU/(2*sin(pi()/(180*3600*2)));
  switch(lu) {
  case   MPARSEC: return _parsec*1.0e6; break;
  case   KPARSEC: return _parsec*1.0e3; break;
  case    PARSEC: return _parsec;       break;
  case        LY: return c()*31557600;  break;
  case        AU: return _AU;           break;
  case        KM: return 1000;          break;
  case         M: return 1;             break;
  case        CM: return 0.01;          break;
  }
  ORSA_DEBUG("problems...");
  return 1;
}

double Unit::getMassScale(const MassUnit & mu) {
  switch(mu) {
  case   KG: return 1;     break;
  case GRAM: return 0.001; break;
  }
  ORSA_DEBUG("problems...");
  return 1;
}
