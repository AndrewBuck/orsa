#include "NEAR.h"

#include <QApplication>

#include <orsaQt/debug.h>

#include <orsaOSG/viz.h>
#include <orsaOSG/FindNamedNodeVisitor.h>
#include <orsaOSG/Track.h>
#include <orsaOSG/DepthPartitionNode.h>
#include <orsaOSG/DistanceAccumulator.h>

#include <iostream>

#include <osgGA/KeySwitchMatrixManipulator>
#include <osgGA/NodeTrackerManipulator>

#include "viz.h"

// main purpose: disable the rotation after the release of the mouse button
class VestaNodeTrackerManipulator : public osgGA::NodeTrackerManipulator {
 public:
  bool handle(const osgGA::GUIEventAdapter & ea,
	      osgGA::GUIActionAdapter      & us) {
    if (ea.getEventType() == osgGA::GUIEventAdapter::RELEASE) {
      osg::ref_ptr<osgGA::GUIEventAdapter> mod_ea = new osgGA::GUIEventAdapter(ea);
      mod_ea->setButtonMask(99);
      return osgGA::NodeTrackerManipulator::handle((*mod_ea.get()),us);
    } else {
      return osgGA::NodeTrackerManipulator::handle(ea,us);
    }
  }
};


int main(int argc, char ** argv) {
  
  QApplication app(argc, argv);
  
  // orsaQt::Debug::instance()->initTimer();
  //
  orsa::Debug::instance()->initTimer();
  
  osg::ref_ptr<orsa::BodyGroup> bg = run();
  
  if (!bg.get()) {
    exit(0);
  }	
  
  if (0) {
    
    // graphics...
    
    osg::ref_ptr<orsaOSG::Viz> viz = new orsaOSG::Viz(bg.get(),
						      600.0);
    
    ViewerQT * viewerWindow = new ViewerQT;
  
    osg::Group * rootNode = viz->createRoot();
    
    viz->_at->setCentralBody(bg->getBody("EROS"));
    
    // viewerWindow->setSceneData(rootNode);
    //
    osg::ref_ptr<DepthPartitionNode> dpn = new DepthPartitionNode;
    dpn->addChild(rootNode);
    dpn->setActive(true);
    //
    viewerWindow->setSceneData(dpn.get());
    //
    // depth partion node only supports single window/single threaded at present.
    viewerWindow->setThreadingModel(osgViewer::Viewer::SingleThreaded);
    
    orsaOSG::FindNamedNodeVisitor fnnv("EROS");
    rootNode->accept(fnnv);
    if (!fnnv._foundNodes.empty()) {
      // set up the node tracker.
      // osgGA::NodeTrackerManipulator * tm = new osgGA::NodeTrackerManipulator;
      VestaNodeTrackerManipulator * tm = new VestaNodeTrackerManipulator;
      
      osgGA::NodeTrackerManipulator::TrackerMode   trackerMode =
	osgGA::NodeTrackerManipulator::NODE_CENTER; // NODE_CENTER_AND_ROTATION
      
      osgGA::NodeTrackerManipulator::RotationMode rotationMode =
	osgGA::NodeTrackerManipulator::TRACKBALL;   // TRACKBALL ELEVATION_AZIM
      
      tm->setTrackerMode(  trackerMode);
      tm->setRotationMode(rotationMode);
      //
      tm->setTrackNode(fnnv._foundNodes.front().get());
      viewerWindow->setCameraManipulator(tm);
    }
    
    viewerWindow->show();
    
    app.connect(&app,
		SIGNAL(lastWindowClosed()), 
		&app, 
		SLOT(quit()));
    
    return app.exec();
  
  } else {
    
    // no graphics
    
    return 0;
    
  } 
  
}
