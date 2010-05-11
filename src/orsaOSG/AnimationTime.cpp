#include <orsaOSG/AnimationTime.h>

#include <orsaSolarSystem/datetime.h>

using namespace orsaOSG;

orsa::Time AnimationTime::getSimulationTime(const int frameID) const {
  
    if (!(_initialTime.isSet())) {
        _initialTime = QTime::currentTime();
    }
  
    if ( (!(_lastFrameID.isSet())) ||
         ((_lastFrameID.isSet()) && (frameID != _lastFrameID.getRef())) ) {
    
        orsa::Time t_start, t_stop;
        _bg->getGlobalInterval(t_start,t_stop,false);
        const orsa::Time period = t_stop - t_start;
    
        /* 
           ORSA_DEBUG("period.getMuSec(): %Zi",
           period.getMuSec().get_mpz_t());
        */
    
        if (period.getMuSec() == 0) {
            return t_start;
        }
    
        if (_realTime) {
            _lastSimulationTime = orsaSolarSystem::now(); // accuracy?
            if (_lastSimulationTime < t_start) _lastSimulationTime = t_start; 
            if (_lastSimulationTime > t_stop)  _lastSimulationTime = t_stop; 
            _lastFrameID = frameID;
        } else {
            _elapsedMSec = _initialTime.getRef().elapsed();
            _lastFrameID = frameID;
      
            if (_timeMultiplier < 0) {
                _lastSimulationTime = 
                    t_stop + 
                    orsa::Time((_elapsedMSec*mpz_class(1000*_timeMultiplier)) % period.getMuSec());
            } else {
                _lastSimulationTime = 
                    t_start + 
                    orsa::Time((_elapsedMSec*mpz_class(1000*_timeMultiplier)) % period.getMuSec());
            }
        }
    
        emit simulationTimeChanged(_lastSimulationTime);
    }
  
    /* 
       ORSA_DEBUG("frameID: %03i   simTime: %f   t: %f this: %x",
       frameID,
       orsa::FromUnits(FromUnits(_lastSimulationTime.getMuSec(), 
       orsa::Unit::MICROSECOND), 
       orsa::Unit::DAY,-1)(),
       _lastSimulationTime.get_d(),
       this);	     
    */
  
    return _lastSimulationTime;
}

const orsa::Vector AnimationTime::centralBodyPosition(const int frameID) const {
    return centralBodyPosition(getSimulationTime(frameID));
}

// to be optimized / cached for effective use of _lastCentralBodyPosition
const orsa::Vector AnimationTime::centralBodyPosition(const orsa::Time & t) const {
  
    if (_centralBody.get()) {
    
        if ( (!(_lastCentralBodyTime.isSet())) ||
             ((_lastCentralBodyTime.isSet()) && (t != _lastCentralBodyTime.getRef())) ) {
      
            _bg->getInterpolatedPosition(_lastCentralBodyPosition,
                                         _centralBody.get(),
                                         t);
            _lastCentralBodyTime = t;
        }
    
        return _lastCentralBodyPosition;
    
    } else if (_followCenterOfMass) {
    
        if ( (!(_lastCentralBodyTime.isSet())) ||
             ((_lastCentralBodyTime.isSet()) && (t != _lastCentralBodyTime.getRef())) ) {
      
            // _lastCentralBodyPosition = _bg->centerOfMassPosition(t);
            //
            orsa::Vector rcm,vcm;
            _bg->centerOfMassPosVel(rcm,vcm,t);
            _lastCentralBodyPosition = rcm;
      
            _lastCentralBodyTime = t;
        }
    
        return _lastCentralBodyPosition;

    } else {
    
        return orsa::Vector(0,0,0);

    }
}
