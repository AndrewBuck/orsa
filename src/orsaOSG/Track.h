#ifndef __ORSA_OSG_TRACK_H__
#define __ORSA_OSG_TRACK_H__

#include <osg/Geode>
#include <osg/Geometry>

#include <orsaOSG/AnimationTime.h>

#include <orsa/bodygroup.h>
#include <orsa/unit.h>

#include <QMutex>
#include <QThread>
#include <QTimer>

namespace orsaOSG {
  
  class Track;
  
  class TrackCallback : public osg::NodeCallback {
  public: 
    TrackCallback(const Track * track,
		  const orsaOSG::AnimationTime * at) :
      osg::NodeCallback(),
      _track(track),
      _at (at) {
      oldTrackVectorSize=0; 
    }  
  public:
    void operator () (osg::Node *, 
		      osg::NodeVisitor *);
  protected:
    osg::ref_ptr<const Track> _track; 
  protected:
    unsigned int oldTrackVectorSize;
  protected:
    osg::ref_ptr<const orsaOSG::AnimationTime> _at;
  protected:
    osg::ref_ptr<osg::Geometry> _trackGeometry;
  };
  
  class TrackElement {
  public:
    orsa::Time   t;
    orsa::Vector r;
  };
  
  class TrackFillThread : public QThread {  

    Q_OBJECT;
  
  public:
    TrackFillThread(orsa::BodyGroup  * bg,
		    const orsa::Body * b,
		    const orsa::Body * ref_b,
		    const orsa::Time & dt,
		    const bool         groundTrack);
  public:
    virtual ~TrackFillThread() { }
    
  protected:
    void run();
    
  public slots:
    void abort() {
      doAbort = true;
    }
  private:
    bool doAbort;
    
  protected:    
    osg::ref_ptr<orsa::BodyGroup>  _bg;
    osg::ref_ptr<const orsa::Body> _b;
    osg::ref_ptr<const orsa::Body> _ref_b;
  protected:    
    const orsa::Time               _dt;
  protected:
    const bool _groundTrack;
    
  public:
    std::vector<TrackElement> data;
    QMutex                    dataMutex;
  };
  
  class Track : public QObject, public osg::Geode {
    
    Q_OBJECT;

  protected:
    ~Track() {
      if (trackFillThread) {
	trackFillThread->abort();
	trackFillThread->wait();
      }
    }
    
  public:
    Track(orsa::BodyGroup  * bg,
	  const orsa::Body * b,
	  const orsa::Body * ref_b,
	  const orsa::Time & dt,
	  const orsaOSG::AnimationTime * at,
	  const bool groundTrack) :
      QObject(),
      osg::Geode(),
      _bg(bg),
      _b(b),
      _ref_b(ref_b),
      _dt(dt),
      _at(at),
      _groundTrack(groundTrack) { 
      
      setUpdateCallback(new TrackCallback(this,
					  _at.get()));
      
      trackVector.clear();
      
      {
	const orsa::Time dt(0,0,15,0,0);
	trackFillThread = new TrackFillThread(_bg.get(),
					      _b.get(),
					      _ref_b.get(),
					      _dt,
					      _groundTrack);
	//
	// this returns instantly!
	trackFillThread->start(QThread::LowestPriority);
      }
      
      trackFillTimer.start(200);
      
      connect(&trackFillTimer,
	      SIGNAL(timeout()),
	      this,
	      SLOT(trackFill()));
      
      connect(trackFillThread,
	      SIGNAL(finished()),
	      this,
	      SLOT(trackFill()));
      
      connect(trackFillThread,
	      SIGNAL(finished()),
	      &trackFillTimer,
	      SLOT(stop()));
      
      // one call, to fill a bit...
      trackFill();
      
    }
      
  protected slots:
    void trackFill();
    
  protected:
    TrackFillThread * trackFillThread;
    mutable QTimer    trackFillTimer;
    
  public:
    const std::vector<TrackElement> & getTrackVector() const {
      return trackVector;
    }	
  protected:
    std::vector<TrackElement> trackVector;
  public:
    mutable QMutex            trackVectorMutex;
    
  protected:
    osg::ref_ptr<orsa::BodyGroup>  _bg;
    osg::ref_ptr<const orsa::Body> _b;
    osg::ref_ptr<const orsa::Body> _ref_b;
  protected:
    const orsa::Time               _dt;
  protected:
    osg::ref_ptr<const orsaOSG::AnimationTime> _at;
  protected:
    const bool _groundTrack;
    
  };
  
}; // namespace orsaOSG

#endif // __ORSA_OSG_TRACK_H__
