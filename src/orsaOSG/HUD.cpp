#include <orsaOSG/HUD.h>

using namespace orsaOSG;

#include <osg/Geode>

#include <orsaOSG/HUDCallback.h>

HUD::HUD(orsaOSG::AnimationTime * at) : 
  osg::CameraNode(),
  _at(at) { 
  
  osg::Geode * geode = new osg::Geode();
  
  std::string font("fonts/arial.ttf");
  
  // turn lighting off for the text and disable depth test to ensure its always ontop.
  osg::StateSet * stateset = geode->getOrCreateStateSet();
  stateset->setMode(GL_LIGHTING,osg::StateAttribute::OFF);
  
  // osg::Vec3 position(150.0f,800.0f,0.0f);
  osg::Vec3 position(20.0f,20.0f,0.0f);
  
  osg::Vec3 delta(0.0f,-120.0f,0.0f);
  
  {
    // osgText::Text * text = new  osgText::Text;
    _timeLabel = new osgText::Text;
    
    geode->addDrawable(_timeLabel.get());
    
    _timeLabel->setFont(font);
    _timeLabel->setPosition(position);
    _timeLabel->setText("time here...");
    
    // text->setName("timeLabel");
    
    position += delta;
  }    
  
  //  osg::CameraNode * camera = new osg::CameraNode;
  
  // set the projection matrix
  // camera->setProjectionMatrix(osg::Matrix::ortho2D(0,1280,0,1024));
  setProjectionMatrix(osg::Matrix::ortho2D(0,1280,0,1024));
  
  // set the view matrix    
  // camera->setReferenceFrame(osg::Transform::ABSOLUTE_RF);
  setReferenceFrame(osg::Transform::ABSOLUTE_RF);
  // camera->setViewMatrix(osg::Matrix::identity());
  setViewMatrix(osg::Matrix::identity());
  
  // only clear the depth buffer
  // camera->setClearMask(GL_DEPTH_BUFFER_BIT);
  setClearMask(GL_DEPTH_BUFFER_BIT);
  
  // draw subgraph after main camera view.
  // camera->setRenderOrder(osg::CameraNode::POST_RENDER);
  setRenderOrder(osg::CameraNode::POST_RENDER);
  
  // camera->addChild(geode);
  addChild(geode);

  setUpdateCallback(new HUDCallback(this, _at.get()));
}
