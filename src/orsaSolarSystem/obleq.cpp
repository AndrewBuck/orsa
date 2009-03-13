#include <orsaSolarSystem/obleq.h>

#include <orsa/datetime.h>
#include <orsaSolarSystem/datetime.h>

using namespace orsa;
using namespace orsaSolarSystem;

orsa::Angle orsaSolarSystem::obleq(const orsa::Time & t) {
  const orsa::Time t_UT = orsaSolarSystem::ToTimeScale(t,orsaSolarSystem::TS_UT);
  const double T = (timeToJulian(t_UT) - 2451545)/36525;
  orsa::Angle obleq_angle;
  obleq_angle.setDPS(23,26,21.448+((0.001813*T-0.00059)*T-46.8150)*T);
  return obleq_angle;
}

orsa::Angle orsaSolarSystem::obleqJ2000() {
  static orsa::Angle _obleqJ2000 = orsaSolarSystem::obleq(J2000());
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
