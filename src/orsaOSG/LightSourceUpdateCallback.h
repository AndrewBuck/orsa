#ifndef _ORSA_OSG_LIGHTSOURCE_UPDATE_CALLBACK_
#define _ORSA_OSG_LIGHTSOURCE_UPDATE_CALLBACK_

#include <osg/NodeCallback>
#include <osg/NodeVisitor>

#include <orsaOSG/AnimationTime.h>

#include <orsa/bodygroup.h>
#include <orsa/unit.h>

namespace orsaOSG {
  
    class LightSourceUpdateCallback : public osg::NodeCallback {
    public:
        LightSourceUpdateCallback(orsa::BodyGroup  * bg,
                                  const orsa::Body * b,
                                  orsaOSG::AnimationTime * at) : 
            osg::NodeCallback(),
            _bg(bg),
            _b(b),
            _at(at) { }
    protected:
        ~LightSourceUpdateCallback() { }
    
    public:
        void operator () (osg::Node * node, osg::NodeVisitor * nv) {
     
            if (nv->getVisitorType()==osg::NodeVisitor::UPDATE_VISITOR && nv->getFrameStamp()) {
	
                osg::LightSource * lightSource = dynamic_cast<osg::LightSource * > (node);
                //
                if (lightSource) {
	  
                    const orsa::Time simulationTime = _at->getSimulationTime(nv->getFrameStamp()->getFrameNumber());
	  
                    if (_b->alive(simulationTime)) {
	    
                        if ( (!(_lastSimulationTime.isSet())) ||
                             ((_lastSimulationTime.isSet()) && (simulationTime != _lastSimulationTime.getRef())) ) {
	      
                            _lastSimulationTime = simulationTime;
	      
                            orsa::Vector position;
	      
                            if (_bg->getInterpolatedPosition(position,
                                                             _b.get(),
                                                             simulationTime)) {
		
                                osg::Light * light = lightSource->getLight();
                                //
                                light->setPosition(osg::Vec4((position-_at->centralBodyPosition(simulationTime)).getVec3d(), 
                                                             1.0));
                                //
                                lightSource->setLight(light);
		
                                lightSource->setLocalStateSetModes( osg::StateAttribute::ON );
                                lightSource->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::ON);
		
                            }
	      
                        }
	    
                    } else {
	    
                        lightSource->setLocalStateSetModes( osg::StateAttribute::OFF );
                        lightSource->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::OFF);
	    
                    }
	  
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

#endif // _ORSA_OSG_LIGHTSOURCE_UPDATE_CALLBACK_
