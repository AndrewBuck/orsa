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

// equatorial coordinates of the north galactic pole
// data relative to B1950 epoch
static const double  ra_g = 192.25*orsa::degToRad();
static const double dec_g =  27.40*orsa::degToRad();
static const double   l_g =  33.00*orsa::degToRad();

void orsaSolarSystem::equatorialToGalactic(double & l,
                                           double & b,
                                           const double & ra,
                                           const double & dec) {
    double c_ra, s_ra;
    orsa::sincos(ra,&s_ra,&c_ra);
    double c_dec, s_dec;
    orsa::sincos(dec,&s_dec,&c_dec);
    b = asin(c_dec*cos(dec_g)*cos(ra-ra_g)+s_dec*sin(dec_g));
    l = l_g+atan2((s_dec/c_dec)*cos(dec_g)-cos(ra-ra_g)*sin(dec_g),sin(ra-ra_g));
}

void orsaSolarSystem::galacticToEquatorial(double & ra,
                                           double & dec,
                                           const double & l,
                                           const double & b) {
    double c_b, s_b;
    orsa::sincos(b,&s_b,&c_b);
    double c_l, s_l;
    orsa::sincos(l,&s_l,&c_l);
    dec = asin(c_b*cos(dec_g)*sin(l-l_g)+s_b*sin(dec_g));
    ra = ra_g+atan2(cos(l-l_g),(s_b/c_b)*cos(dec_g)-sin(dec_g)*sin(l-l_g));
}
