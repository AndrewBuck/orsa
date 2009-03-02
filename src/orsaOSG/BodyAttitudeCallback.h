#ifndef _ORSA_OSG_BODY_ATTITUDE_CALLBACK_
#define _ORSA_OSG_BODY_ATTITUDE_CALLBACK_

#include <osg/NodeCallback>
#include <osg/NodeVisitor>
#include <osg/PositionAttitudeTransform>

#include <orsaOSG/AnimationTime.h>

// #include <orsa/attitude.h>
#include <orsa/bodygroup.h>
#include <orsa/print.h>
#include <orsa/unit.h>
#include <orsa/util.h>

namespace orsaOSG {
  
  class BodyAttitudeCallback : public osg::NodeCallback {
  public:
    BodyAttitudeCallback(orsa::BodyGroup  * bg,
			 const orsa::Body * b,
			 orsaOSG::AnimationTime * at) : 
      osg::NodeCallback(),
      _bg(bg),
      _b(b),
      _at(at) { }
  protected:
    ~BodyAttitudeCallback() { }
    
  public:
    void operator () (osg::Node * node, osg::NodeVisitor * nv) {
      
      /* 
	 ORSA_DEBUG("<<<body [%s]>>>",
	 _b->getName().c_str());
      */
      
      if (nv->getVisitorType()==osg::NodeVisitor::UPDATE_VISITOR && nv->getFrameStamp()) {
        
	osg::PositionAttitudeTransform * pat = dynamic_cast<osg::PositionAttitudeTransform * > (node);
	//
	if (pat) {
	  
       	  const orsa::Time simulationTime = _at->getSimulationTime(nv->getFrameStamp()->getFrameNumber());
	  
	  /* 
	     ORSA_DEBUG("simTime: %Ff  body [%s]",
	     simulationTime.asDouble().get_mpf_t(),
	     _b->getName().c_str());
	  */
	  
	  if (_b->alive(simulationTime)) {
	    
	    if ( (!(_lastSimulationTime.isSet())) ||
		 ((_lastSimulationTime.isSet()) && (simulationTime != _lastSimulationTime.getRef())) ) {
	      
	      _lastSimulationTime = simulationTime;
	      
	      orsa::Vector position;
	      
	      /* 
		 if (_bg->getInterpolatedPosition(position,
		 _b.get(),
		 simulationTime)) {
	      */
	      
	      // pat->setPosition(position.getVec3d());
	      // pat->setPosition((1.0e9*position).getVec3d());
	      // pat->setPosition(osg::Vec3d(0,0,0));
	      
	      // pat->setPosition((position-_at->centralBodyPosition(simulationTime)).getVec3d());
	      //
	      pat->setPosition(osg::Vec3d(0,0,0));
	      
	      // if (_b->getAttitude()) {
	      /* 
		 if (_bg->getInterpolatedAttitude(_b.get(),simulationTime)) {
		 osg::Quat q;
		 // q.set(_b->getAttitude()->globalToLocal(simulationTime).getMatrixd());
		 q.set(_bg->getInterpolatedAttitude(_b.get(),simulationTime)->globalToLocal(simulationTime).getMatrixd());
		 pat->setAttitude(q);
		 } else {
		 pat->setAttitude(osg::Quat(0,0,0,1));
		 }
	      */
	      //
	      {
		// osg::ref_ptr<orsa::Attitude> attitude = new orsa::BodyAttitude(_b.get(),_bg.get());
		osg::Quat q;
	    	// q.set(attitude->globalToLocal(simulationTime).getMatrixd());
		q.set(orsa::globalToLocal(_b.get(),
					  _bg.get(),
					  simulationTime).getMatrixd());
		pat->setAttitude(q);
	      }
	      
	      // ORSA_DEBUG("used...");
	      
	      
	      
	      if (1) {
		// debug
		
		// osg::ref_ptr<orsa::Attitude> attitude = new orsa::BodyAttitude(_b.get(),_bg.get());
		
		// const orsa::Matrix g2l = attitude->globalToLocal(simulationTime);
		// const orsa::Matrix l2g = attitude->localToGlobal(simulationTime);
		
		const orsa::Matrix g2l = orsa::globalToLocal(_b.get(),_bg.get(),simulationTime);
		const orsa::Matrix l2g = orsa::localToGlobal(_b.get(),_bg.get(),simulationTime);
		
		/* 
		   print(l2g);
		   print(g2l);
		   print(g2l*l2g);
		   print(l2g*g2l);
		*/
		
		orsa::IBPS ibps;	
		_bg->getInterpolatedIBPS(ibps,
					 _b.get(),
					 simulationTime);
		
		/* 
		   {
		   // test
		   const orsa::Vector r(-1,1,-3);
		   ORSA_DEBUG("test...");
		   const orsa::Quaternion q = ibps.rotational->getQ();
		   print(q*r*conjugate(q));
		   print(QuaternionToMatrix(q)*r);
		   }
		*/
		
		/* 
		   orsa::Matrix inertiaMoment = 
		   _b->getMass() * 
		   _b->getPaulMoment()->getInertiaMoment();
		   
		   orsa::Vector L = l2g * inertiaMoment * g2l * ibps.rotational->getOmega();
		   
		   // print(inertiaMoment);
		   // print(l2g * inertiaMoment * g2l);
		   
		   const orsa::Vector uL = L.normalized();
		   const orsa::Vector uO = ibps.rotational->getOmega().normalized();
		   
		   ORSA_DEBUG("[%x|t:%12.3Ff][%s] L: %12.9Fe x [%+6.3Ff,%+6.3Ff,%+6.3Ff]   omega: %12.9Fe x [%+6.3Ff,%+6.3Ff,%+6.3Ff]   sp: %12.9Ff",
		   &ibps,
		   simulationTime.asDouble().get_mpf_t(),
		   _b->getName().c_str(),
		   L.length().get_mpf_t(),
		   uL.getX().get_mpf_t(), uL.getY().get_mpf_t(), uL.getZ().get_mpf_t(),
		   ibps.rotational->getOmega().length().get_mpf_t(),
		   uO.getX().get_mpf_t(), uO.getY().get_mpf_t(), uO.getZ().get_mpf_t(),
		   (uL*uO).get_mpf_t());
		*/
		
		/* 
		   ORSA_DEBUG("g2l: m11 = %Fg",
		   g2l.getM11().get_mpf_t());
		*/
	      }
	      
	      
	      
	      
	      // pat->setScale(osg::Vec3d(1,1,1));
	      // pat->setPivotPoint(osg::Vec3d(0,0,0));
	      
	      if (0) {
		ORSA_DEBUG("simTime: %Ff  body [%s] position: %Ff %Ff %Ff",
			   orsa::FromUnits(FromUnits(simulationTime.getMuSec(), orsa::Unit::MICROSECOND), orsa::Unit::DAY,-1).get_mpf_t(),
			   _b->getName().c_str(),
			   orsa::FromUnits(position.getX(), orsa::Unit::AU,-1).get_mpf_t(),
			   orsa::FromUnits(position.getY(), orsa::Unit::AU,-1).get_mpf_t(),
			   orsa::FromUnits(position.getZ(), orsa::Unit::AU,-1).get_mpf_t());
	      }
	      
	      /* 
		 } else {
		 
		 ORSA_DEBUG("problems... simTime: %Ff   body [%s]",
		 orsa::FromUnits(FromUnits(simulationTime.getMuSec(), orsa::Unit::MICROSECOND), orsa::Unit::DAY,-1).get_mpf_t(),
		 _b->getName().c_str());
		 
		 }
	      */
	      
	    } else {
	      
	      // nothing to do...
	      // ORSA_DEBUG("skipping...");
	      
	    }
	    
	  } else {
	    
	    /* 
	       ORSA_DEBUG("body [%s] not alive at time %Ff",
	       _b->getName().c_str(),
	       simulationTime.asDouble().get_mpf_t());
	    */
	  }
	  
	} else {
	  
	  ORSA_DEBUG("pat == 0...");
	  
	}
	
      }
      
      // must call any nested node callbacks and continue subgraph traversal.
      NodeCallback::traverse(node,nv);
    }
    
  protected:    
    osg::ref_ptr<orsa::BodyGroup> _bg;
    osg::ref_ptr<const orsa::Body> _b;
  protected:
    osg::ref_ptr<AnimationTime> _at;
  protected:
    orsa::Cache<orsa::Time> _lastSimulationTime;
  };
  
}; // namespace orsaOSG

#endif // _ORSA_OSG_BODY_ATTITUDE_CALLBACK_
