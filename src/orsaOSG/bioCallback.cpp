#include <orsaOSG/bioCallback.h>

#include <osg/Geometry>
#include <osg/Switch>

using namespace orsa;
using namespace orsaOSG;
using namespace osg;

bioCallback::bioCallback(orsa::BodyGroup  * bg,
			 const orsa::Body * b,
			 orsaOSG::AnimationTime * at) :
  NodeCallback(),
  _bg(bg),
  _b(b),
  _at(at) {
  /* 
     orsa::Time t_start, t_stop;
     // bg->getCommonInterval(t_start,t_stop,false);
     bg->getGlobalInterval(t_start,t_stop,false);
     _timing = new orsaOSG::AnimationTiming(-t_start.asDouble().get_d()/time_multiplier,
     time_multiplier);
  */
}

void bioCallback::operator()(osg::Node * node, osg::NodeVisitor * nv) {
  
  if ((_bg.get() == 0) || (_b.get() == 0)) {
    
    ORSA_ERROR("zero pointers...");
    
    traverse(node,nv);
    
    return;
  }
  
  const orsa::Time & simulationTime = _at->getSimulationTime(nv->getFrameStamp()->getFrameNumber());
  
  osg::Switch * s = dynamic_cast<osg::Switch *>(node);
  if (s) {
    if (_b->alive(simulationTime)) {
      s->setAllChildrenOn();
    } else {
      s->setAllChildrenOff();
    }
  }
  
  // note, callback is responsible for scenegraph traversal so
  // should always include call traverse(node,nv) to ensure 
  // that the rest of cullbacks and the scene graph are traversed.
  traverse(node,nv);
}
