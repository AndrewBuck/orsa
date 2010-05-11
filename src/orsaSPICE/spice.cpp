#include <orsaSPICE/spice.h>

#include <orsaSolarSystem/datetime.h>

// CSPICE prototypes and definitions.      
#include "SpiceUsr.h"

using namespace orsa;
using namespace orsaSolarSystem;
using namespace orsaSPICE;

SPICE * SPICE::_instance = 0;

void SPICE::loadKernel(const std::string & filename) {
    mutex.lock();
    furnsh_c(filename.c_str());
    mutex.unlock();
}

void SPICE::unloadKernel(const std::string & filename) {
    mutex.lock();
    unload_c(filename.c_str());
    mutex.unlock();
}

void SPICE::getPosVel(const std::string & target,
                      const orsa::Time  & ephemerisTime,
                      const std::string & observer,
                      orsa::Vector      & relativePosition,
                      orsa::Vector      & relativeVelocity) {
    SpiceDouble lt;
    SpiceDouble state[6];
    mutex.lock();
    spkezr_c(target.c_str(),
             SPICETime(ephemerisTime), 
             _global.getRef().c_str(), // "J2000", //! note: J2000 is the equatorial! eclipj2000 is the ecliptic one!
             "NONE",
             observer.c_str(),
             state,
             &lt);
    mutex.unlock();
    relativePosition.set(FromUnits(state[0],Unit::KM),
                         FromUnits(state[1],Unit::KM),
                         FromUnits(state[2],Unit::KM));
    relativeVelocity.set(FromUnits(FromUnits(state[3],Unit::KM),Unit::SECOND,-1),
                         FromUnits(FromUnits(state[4],Unit::KM),Unit::SECOND,-1),
                         FromUnits(FromUnits(state[5],Unit::KM),Unit::SECOND,-1));
}

void SPICE::getPosVel(const std::string & target,
                      const orsa::Time  & ephemerisTime,
                      orsa::Vector      & relativePosition,
                      orsa::Vector      & relativeVelocity) {
    getPosVel(target,
              ephemerisTime,
              _observer.getRef(),
              relativePosition,
              relativeVelocity);
}

orsa::Matrix SPICE::localToGlobal(const std::string & local,
                                  const orsa::Time  & ephemerisTime) {
    SpiceDouble rotate[3][3];
    mutex.lock();
    pxform_c(local.c_str(),
             _global.getRef().c_str(),
             SPICETime(ephemerisTime),
             rotate);
    mutex.unlock();
    orsa::Matrix m;
    m.set(rotate[0][0],rotate[0][1],rotate[0][2],
          rotate[1][0],rotate[1][1],rotate[1][2],
          rotate[2][0],rotate[2][1],rotate[2][2]);
    return m;
}

double SPICE::SPICETime(const orsa::Time & t) {
    // Corrected and verified: no TimeScale needed
    return FromUnits((t-J2000()).get_d(),Unit::SECOND,-1);
}

orsa::Matrix SPICE::globalToLocal(const std::string & local,
                                  const orsa::Time  & ephemerisTime) {
    // SINTAX: pxform_c ( from, to,  et,  rotate );
    SpiceDouble rotate[3][3];
    mutex.lock();
    pxform_c(_global.getRef().c_str(),
             local.c_str(),
             SPICETime(ephemerisTime),
             rotate);
    mutex.unlock();
    orsa::Matrix m;
    m.set(rotate[0][0],rotate[0][1],rotate[0][2],
          rotate[1][0],rotate[1][1],rotate[1][2],
          rotate[2][0],rotate[2][1],rotate[2][2]);
    return m;
}
