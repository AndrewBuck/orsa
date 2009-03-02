#ifndef _ORSA_OSG_ANIMATION_TIME_
#define _ORSA_OSG_ANIMATION_TIME_

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/bodygroup.h>
#include <orsa/datetime.h>
#include <orsa/double.h>
#include <orsa/cache.h>

#include <QObject>
#include <QTime>

namespace orsaOSG {
  
  // should rename to something as AnimationData... because I'm adding here the data relative to the "central" body (only used for position...).
  
  class AnimationTime : public QObject, public osg::Referenced {
    
    Q_OBJECT;
    
  public:    
    AnimationTime(orsa::BodyGroup * bg) : 
      QObject(),
      osg::Referenced(),
      _bg(bg) { 
      
      _timeMultiplier = 1;
      _realTime       = false;
    }
  protected:
    ~AnimationTime() { }
    
  public:
    void setTimeMultiplier(const double & m) { _timeMultiplier = m; }
  public:
    const double & getTimeMultiplier() const { return _timeMultiplier; }
  protected:
    double _timeMultiplier;
    
  public:
    void setRealTime(const bool b) { _realTime = b; }
  public:
    bool getRealTime() const { return _realTime; }
  protected:
    bool _realTime;
    
  public:
    orsa::Time getSimulationTime(const int frameID) const;
    
  signals:
    void simulationTimeChanged(const orsa::Time & t) const;
    
  public:
    void setCentralBody(const orsa::Body * b) {
      _centralBody = b;
      if (_centralBody.get()) {
	_followCenterOfMass = false;
      }
      _lastCentralBodyTime.reset();
    }
  public:
    const orsa::Body * getCentralBody() const {
      return _centralBody.get();
    }
  protected:
    osg::ref_ptr<const orsa::Body> _centralBody;
    
  public:
    void followCenterOfMass(const bool b = true) {
      _followCenterOfMass = b;
      if (_followCenterOfMass) {
	_centralBody = 0;
      }
      _lastCentralBodyTime.reset();
    }
  public:
    bool followingCenterOfMass() const {
      return _followCenterOfMass;
    }
  protected:
    bool _followCenterOfMass;
    
  public:
    const orsa::Vector centralBodyPosition(const int frameID) const;
  public:
    const orsa::Vector centralBodyPosition(const orsa::Time & t) const;
  protected:
    mutable orsa::Vector            _lastCentralBodyPosition;
    mutable orsa::Cache<orsa::Time> _lastCentralBodyTime;
    
  protected:
    osg::ref_ptr<orsa::BodyGroup> _bg;
  protected:
    mutable orsa::Cache<QTime> _initialTime;
  protected:
    mutable int _elapsedMSec;
  protected:    
    mutable orsa::Cache<int> _lastFrameID;
  protected:
    mutable orsa::Time       _lastSimulationTime;
  };
  
}; // namespace orsaOSG

#endif // _ORSA_OSG_ANIMATION_TIME_
