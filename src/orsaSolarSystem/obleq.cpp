#include <orsaSolarSystem/obleq.h>

#include <orsa/datetime.h>
#include <orsaSolarSystem/datetime.h>

using namespace orsa;
using namespace orsaSolarSystem;

double orsaSolarSystem::obleq(const orsa::Time & t) {
  // T in centuries from J2000 
  // const double T = (d.GetJulian() - 2451545.0)/36525.0;
  // DOUBLE-CHECK this "UT"!!!
  // Time J2000; J2000.setJ2000();
  const double T = FromUnits((t - J2000()).get_d(),Unit::YEAR,-1)/36525.0;
  // double a;
  // updated Feb 2004
  // a.SetDPS(23,26,21.448+((0.001813*T-0.00059)*T-46.8150)*T);
  //
  const double oneOverSixty = 1.0/60.0;
  const double obleqDEG = (23+oneOverSixty*(26+oneOverSixty*(21.448+((0.001813*T-0.00059)*T-46.8150)*T)));
  //
  /* 
     ORSA_DEBUG("T: %f   obleq: %f [deg]",
     T(),
     obleqDEG());
  */
  //
  // return (degToRad()*(23+oneOverSixty*(26+oneOverSixty*(21.448+((0.001813*T-0.00059)*T-46.8150)*T))));
  return (degToRad()*obleqDEG);
}

double orsaSolarSystem::obleqJ2000() {
  static double _obleqJ2000 = orsaSolarSystem::obleq(J2000());
  return _obleqJ2000;
}

orsa::Matrix orsaSolarSystem::eclipticToEquatorial() {
  static orsa::Matrix _m;
  static bool firstCall = true;
  if (firstCall) {
    _m = Matrix::identity();
    _m.rotX(orsaSolarSystem::obleqJ2000());
    firstCall = false;
  }
  return _m;
}

orsa::Matrix orsaSolarSystem::equatorialToEcliptic() {
  static orsa::Matrix _m;
  static bool firstCall = true;
  if (firstCall) {
    _m = Matrix::identity();
    _m.rotX(-orsaSolarSystem::obleqJ2000());
    firstCall = false;
  }	
  return _m;
}
