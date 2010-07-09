#include <orsaSolarSystem/obleq.h>
#include <orsaSolarSystem/gmst.h>

#include <orsaUtil/observatory.h>

#include <orsaInputOutput/MPC_obscode.h>

#include <orsaSPICE/spiceBodyTranslationalCallback.h>

using namespace orsa;
using namespace orsaUtil;

StandardObservatoryPositionCallback::StandardObservatoryPositionCallback(orsaInputOutput::MPCObsCodeFile * ocf) :
    orsaSolarSystem::ObservatoryPositionCallback(),
    obsCodeFile(ocf) {
  
    bg = new orsa::BodyGroup;
  
    earth = new orsa::Body;
  
    earth->setName("EARTH");
    orsaSPICE::SpiceBodyTranslationalCallback * sbtc = new orsaSPICE::SpiceBodyTranslationalCallback(earth->getName());
    orsa::IBPS ibps;
    ibps.inertial = new orsa::PointLikeConstantInertialBodyProperty(orsaSolarSystem::Data::MEarth());
    ibps.translational = sbtc;
    earth->setInitialConditions(ibps);
  
    bg->addBody(earth.get());
}

bool StandardObservatoryPositionCallback::getPosition(orsa::Vector & position,
                                                      const orsaSolarSystem::Observation * obs) const {
    orsa::Vector velocity;
    return getPosVel(position,velocity,obs);
}

bool  StandardObservatoryPositionCallback::getPosition(orsa::Vector & position,
                                                       const std::string & obsCode,
                                                       const orsa::Time  & t) const {
    const orsaSolarSystem::Observatory & observatory = obsCodeFile->_data.observatory[obsCode];
    if (observatory.moving()) {
        ORSA_DEBUG("problems...");
        return false;
    }
    osg::ref_ptr<orsaSolarSystem::Observation> obs = new orsaSolarSystem::Observation;
    obs->obsCode=obsCode;
    obs->epoch=t;
    return getPosition(position,obs.get());
}

bool StandardObservatoryPositionCallback::getPosVel(orsa::Vector      & position,
                                                    orsa::Vector      & velocity,
                                                    const orsaSolarSystem::Observation * obs) const {
  
    const orsa::Time t = obs->epoch.getRef();
  
    orsa::Vector rEarth, vEarth;
    if (!bg->getInterpolatedPosVel(rEarth,vEarth,earth.get(),t)) { 
        ORSA_DEBUG("problems...");
        return false;
    }
  
    const orsaSolarSystem::Observatory & observatory = obsCodeFile->_data.observatory[obs->obsCode.getRef()];
  
    orsa::Vector obsPos;
  
    if (observatory.moving()) {
    
        const orsaSolarSystem::SatelliteObservation * satelliteObservation =  
            dynamic_cast<const orsaSolarSystem::SatelliteObservation *>(obs);
        if (satelliteObservation) {
            obsPos = satelliteObservation->obsPos.getRef();
        }
    
    } else {
    
        // LMST = GMST + longitude
        double s, c;
        orsa::sincos(orsaSolarSystem::gmst(t).getRad()+observatory.lon.getRef(),&s,&c);
        obsPos = orsa::Vector(observatory.pxy.getRef()*c,
                              observatory.pxy.getRef()*s,
                              observatory.pz.getRef());
    }
  
    obsPos = orsaSolarSystem::equatorialToEcliptic()*obsPos;
  
    position = rEarth + obsPos;
    velocity = vEarth; // should correct for obsVel...
  
    return true;
}

bool StandardObservatoryPositionCallback::getPosVel(orsa::Vector & position,
                                                    orsa::Vector & velocity,
                                                    const std::string & obsCode,
                                                    const orsa::Time  & t) const {
    const orsaSolarSystem::Observatory & observatory = obsCodeFile->_data.observatory[obsCode];
    if (observatory.moving()) {
        ORSA_DEBUG("problems...");
        return false;
    }
    osg::ref_ptr<orsaSolarSystem::Observation> obs = new orsaSolarSystem::Observation;
    obs->obsCode=obsCode;
    obs->epoch=t;
    return getPosVel(position,velocity,obs.get());
}

const orsaSolarSystem::Observatory & StandardObservatoryPositionCallback::getObservatory(const std::string & obsCode) const {
    return obsCodeFile->_data.observatory[obsCode];
}

const orsaSolarSystem::Observatory & StandardObservatoryPositionCallback::getObservatory(const orsaSolarSystem::Observation * obs) const {
    return obsCodeFile->_data.observatory[obs->obsCode.getRef()];
}
