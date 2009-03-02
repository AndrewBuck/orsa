#ifndef _ORSA_OSG_HUD_CALLBACK_
#define _ORSA_OSG_HUD_CALLBACK_

#include <osg/CameraNode>
#include <osg/NodeCallback>

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsaOSG/AnimationTime.h>

namespace orsaOSG {
  
  class HUD;
  
  class HUDCallback : public osg::NodeCallback {
  public:
    HUDCallback(osg::CameraNode * HUD,
		orsaOSG::AnimationTime * at) : 
      osg::NodeCallback(),
      _HUD(HUD),
      _at(at) { }
  protected:    
    ~HUDCallback() { }
    
  public:
    void operator () (osg::Node *, osg::NodeVisitor *);
    
  protected:
    osg::ref_ptr<osg::CameraNode> _HUD;
  protected:
    osg::ref_ptr<AnimationTime> _at;
 };
  
}; // namespace orsaOSG

#endif // _ORSA_OSG_HUD_CALLBACK_
