#ifndef _ORSA_OSG_CENTER_OF_MASS_POSITION_CALLBACK_
#define _ORSA_OSG_CENTER_OF_MASS_POSITION_CALLBACK_

#include <osg/NodeCallback>
#include <osg/NodeVisitor>
#include <osg/PositionAttitudeTransform>

#include <orsaOSG/AnimationTime.h>

#include <orsa/bodygroup.h>
#include <orsa/unit.h>

namespace orsaOSG {
  
    class CenterOfMassPositionCallback : public osg::NodeCallback {
    public:
        CenterOfMassPositionCallback(orsa::BodyGroup  * bg,
                                     orsaOSG::AnimationTime * at) : 
            osg::NodeCallback(),
            _bg(bg),
            _at(at) { }
    protected:
        ~CenterOfMassPositionCallback() { }
    
    public:
        void operator () (osg::Node * node, osg::NodeVisitor * nv) {
      
            if (nv->getVisitorType()==osg::NodeVisitor::UPDATE_VISITOR && nv->getFrameStamp()) {
        
                osg::PositionAttitudeTransform * pat = dynamic_cast<osg::PositionAttitudeTransform * > (node);
                //
                if (pat) {
	  
                    const orsa::Time simulationTime = _at->getSimulationTime(nv->getFrameStamp()->getFrameNumber());
	  
                    if ( (!(_lastSimulationTime.isSet())) ||
                         ((_lastSimulationTime.isSet()) && (simulationTime != _lastSimulationTime.getRef())) ) {
	    
                        _lastSimulationTime = simulationTime;
	    
                        orsa::Vector rcm, vcm;
                        _bg->centerOfMassPosVel(rcm,vcm,simulationTime);
	    
                        pat->setPosition((rcm - _at->centralBodyPosition(simulationTime)).getVec3d());
	    
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
    protected:
        osg::ref_ptr<AnimationTime> _at;
    protected:
        orsa::Cache<orsa::Time> _lastSimulationTime;
    };
  
}; // namespace orsaOSG

#endif // _ORSA_OSG_CENTER_OF_MASS_POSITION_CALLBACK_
