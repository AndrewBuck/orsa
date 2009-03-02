#include <orsaSolarSystem/obleq.h>

#include <orsa/datetime.h>
#include <orsaSolarSystem/datetime.h>

using namespace orsa;
using namespace orsaSolarSystem;

orsa::Double orsaSolarSystem::obleq(const orsa::Time & t) {
  // T in centuries from J2000 
  // const double T = (d.GetJulian() - 2451545.0)/36525.0;
  // DOUBLE-CHECK this "UT"!!!
  // Time J2000; J2000.setJ2000();
  const Double T = FromUnits((t - J2000()).asDouble(),Unit::YEAR,-1)/36525.0;
  // Double a;
  // updated Feb 2004
  // a.SetDPS(23,26,21.448+((0.001813*T-0.00059)*T-46.8150)*T);
  //
  const Double oneOverSixty = one()/Double("60.0");
  const Double obleqDEG = (23+oneOverSixty*(26+oneOverSixty*(21.448+((0.001813*T-0.00059)*T-46.8150)*T)));
  //
  /* 
     ORSA_DEBUG("T: %Ff   obleq: %Ff [deg]",
     T.get_mpf_t(),
     obleqDEG.get_mpf_t());
  */
  //
  // return (degToRad()*(23+oneOverSixty*(26+oneOverSixty*(21.448+((0.001813*T-0.00059)*T-46.8150)*T))));
  return (degToRad()*obleqDEG);
}

orsa::Double orsaSolarSystem::obleqJ2000() {
  static Double _obleqJ2000;
  static unsigned int old_prec = 0;
  if (old_prec != mpf_get_default_prec()) {
    // Time J2000; J2000.setJ2000();
    _obleqJ2000 = orsaSolarSystem::obleq(J2000());
    old_prec = mpf_get_default_prec();
  }	
  return _obleqJ2000;
}

orsa::Matrix orsaSolarSystem::eclipticToEquatorial() {
  static orsa::Matrix _m;
  static unsigned int old_prec = 0;
  if (old_prec != mpf_get_default_prec()) {
    _m = Matrix::identity();
    _m.rotX(orsaSolarSystem::obleqJ2000());
    old_prec = mpf_get_default_prec();
  }	
  return _m;
}

orsa::Matrix orsaSolarSystem::equatorialToEcliptic() {
  static orsa::Matrix _m;
  static unsigned int old_prec = 0;
  if (old_prec != mpf_get_default_prec()) {
    _m = Matrix::identity();
    _m.rotX(-orsaSolarSystem::obleqJ2000());
    old_prec = mpf_get_default_prec();
  }	
  return _m;
}
