#include <orsaOSG/HUDCallback.h>

using namespace orsaOSG;

#include <orsaOSG/HUD.h>

#include <orsa/unit.h>

#include <orsaSolarSystem/datetime.h>

void HUDCallback::operator () (osg::Node * node, osg::NodeVisitor * nv) {
  
  if (nv->getVisitorType()==osg::NodeVisitor::UPDATE_VISITOR && nv->getFrameStamp()) {
    
    orsaOSG::HUD * h = dynamic_cast <orsaOSG::HUD * > (_HUD.get());
    if (h != 0) {
      
      const orsa::Time & simulationTime = _at->getSimulationTime(nv->getFrameStamp()->getFrameNumber());
      
      {
	char label[1024];
	// int y, m, d, H, M, S, ms;
	// gmp_sprintf(label,"%Ff",orsa::FromUnits(FromUnits(simulationTime.getMuSec(), orsa::Unit::MICROSECOND), orsa::Unit::DAY,-1).get_mpf_t());
	
	gmp_sprintf(label,"JD %Ff",orsaSolarSystem::timeToJulian(simulationTime).get_mpf_t());
	//
	/* 
	   const double tf = (double)((float)(orsaSolarSystem::julianTime(simulationTime).get_d()));
	   gmp_sprintf(label,"JD %Ff   float: %f",
	   orsaSolarSystem::julianTime(simulationTime).get_mpf_t(),
	   tf);
	*/
	
	h->_timeLabel->setText(label);
      } 
      
    }
    
  }
  
  // must call any nested node callbacks and continue subgraph traversal.
  NodeCallback::traverse(node,nv);
}
