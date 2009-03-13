#include <orsaSolarSystem/gmst.h>

// #include <orsa/print.h>
#include <orsaSolarSystem/datetime.h>

using namespace orsa;
using namespace orsaSolarSystem;


/* DEBUG: Using our example time of 1994 June 16 at 18h UT, we get:
 *
 *            GMST = 174.7711135
 *                 = 11h 39m 5.0672s
 */

orsa::Angle orsaSolarSystem::gmst(const orsa::Time & t) {
  const orsa::Time t_UT = orsaSolarSystem::ToTimeScale(t,orsaSolarSystem::TS_UT);
  const double T = (timeToJulian(t_UT) - 2451545)/36525;
  const double dayFraction = mpz_class(t_UT.getMuSec() % mpz_class("86400000000")).get_d()/86400000000.0;
  const orsa::Angle gmst_angle = fmod(fmod(orsa::arcsecToRad()*15*(((dayFraction*24+6)*60+41)*60+50.54841+((-0.0000062*T+0.093104)*T+8640184.812866)*T),orsa::twopi())+orsa::twopi(),orsa::twopi());
  /* 
     orsa::print(t);
     orsa::print(gmst_angle);
     ORSA_DEBUG("gmst: %.9f [deg]   T: %.12f   dayFraction: %f",
     gmst_angle.getRad()*orsa::radToDeg(),
     T,
     dayFraction);
  */
  return gmst_angle;
}
