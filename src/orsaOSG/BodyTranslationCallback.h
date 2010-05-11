#ifndef _ORSA_OSG_BODY_TRANSLATION_CALLBACK_
#define _ORSA_OSG_BODY_TRANSLATION_CALLBACK_

#include <osg/NodeCallback>
#include <osg/NodeVisitor>
#include <osg/PositionAttitudeTransform>

#include <orsaOSG/AnimationTime.h>

#include <orsa/bodygroup.h>
#include <orsa/unit.h>

namespace orsaOSG {
  
    class BodyTranslationCallback : public osg::NodeCallback {
    public:
        BodyTranslationCallback(orsa::BodyGroup  * bg,
                                const orsa::Body * b,
                                orsaOSG::AnimationTime * at) : 
            osg::NodeCallback(),
            _bg(bg),
            _b(b),
            _at(at) { }
    protected:
        ~BodyTranslationCallback() { }
    
    public:
        void operator () (osg::Node * node, osg::NodeVisitor * nv) {
      
            if (nv->getVisitorType()==osg::NodeVisitor::UPDATE_VISITOR && nv->getFrameStamp()) {
        
                osg::PositionAttitudeTransform * pat = dynamic_cast<osg::PositionAttitudeTransform * > (node);
                //
                if (pat) {
	  
                    const orsa::Time simulationTime = _at->getSimulationTime(nv->getFrameStamp()->getFrameNumber());
	  
                    /* 
                       ORSA_DEBUG("simTime: %f  body [%s]",
                       simulationTime.get_d(),
                       _b->getName().c_str());
                    */
	  
                    if (_b->alive(simulationTime)) {
	    
                        if ( (!(_lastSimulationTime.isSet())) ||
                             ((_lastSimulationTime.isSet()) && (simulationTime != _lastSimulationTime.getRef())) ) {
	      
                            _lastSimulationTime = simulationTime;
	      
                            orsa::Vector position;
	      
                            if (_bg->getInterpolatedPosition(position,
                                                             _b.get(),
                                                             simulationTime)) {
		
                                // pat->setPosition(position.getVec3d());
                                // pat->setPosition((1.0e9*position).getVec3d());
                                // pat->setPosition(osg::Vec3d(0,0,0));
		
                                pat->setPosition((position-_at->centralBodyPosition(simulationTime)).getVec3d());
		
                                /* 
                                   if (_b->getAttitude()) {
                                   osg::Quat q;
                                   q.set(_b->getAttitude()->globalToLocal(simulationTime).getMatrixd());
                                   pat->setAttitude(q);
                                   } else {
                                   pat->setAttitude(osg::Quat(0,0,0,1));
                                   }
                                */
		
                                // pat->setScale(osg::Vec3d(1,1,1));
                                // pat->setPivotPoint(osg::Vec3d(0,0,0));
		
                                if (0) {
                                    ORSA_DEBUG("simTime: %f  body [%s] position: %f %f %f",
                                               orsa::FromUnits(FromUnits(simulationTime.getMuSec().get_d(), orsa::Unit::MICROSECOND), orsa::Unit::DAY,-1),
                                               _b->getName().c_str(),
                                               orsa::FromUnits(position.getX(), orsa::Unit::AU,-1),
                                               orsa::FromUnits(position.getY(), orsa::Unit::AU,-1),
                                               orsa::FromUnits(position.getZ(), orsa::Unit::AU,-1));
                                }
		
                            } else {
		
                                ORSA_DEBUG("problems... simTime: %f   body [%s]",
                                           orsa::FromUnits(FromUnits(simulationTime.getMuSec().get_d(), orsa::Unit::MICROSECOND), orsa::Unit::DAY,-1),
                                           _b->getName().c_str());
		
                            }
	      
                        } else {
	      
                            // nothing to do...
                            // ORSA_DEBUG("skipping...");
	      
                        }
	    
                    } else {
	    
                        /* 
                           ORSA_DEBUG("body [%s] not alive at time %f",
                           _b->getName().c_str(),
                           simulationTime.get_d());
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

#endif // _ORSA_OSG_BODY_TRANSLATION_CALLBACK_
