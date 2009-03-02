#include <orsaSolarSystem/gmst.h>

#include <orsaSolarSystem/datetime.h>

using namespace orsa;
using namespace orsaSolarSystem;

double orsaSolarSystem::gmst(const orsa::Time & t) {
  const double T = (timeToJulian(orsaSolarSystem::ToTimeScale(t,orsaSolarSystem::TS_UT)) - 2451545.0)/36525.0;
  const double gmstARCSEC = 15.0*(((orsaSolarSystem::dayFraction(t)*24.0+6.0)*60.0+41.0)*60.0+50.54841+((-0.0000062*T+0.093104)*T+8640184.812866)*T);
  return ((gmstARCSEC/3600.0)*degToRad());
}
