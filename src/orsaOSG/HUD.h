#ifndef _ORSA_OSG_HUD_H_
#define _ORSA_OSG_HUD_H_

#include <osg/ref_ptr>

#include <osg/CameraNode>

#include <osgText/Text>

#include <orsaOSG/AnimationTime.h>

namespace orsaOSG {
  
  class HUD : public osg::CameraNode {
  public:
    HUD(orsaOSG::AnimationTime *);
  protected:	
    ~HUD() { }
    
  public:	
    osg::ref_ptr<osgText::Text> _timeLabel;
  protected:
    osg::ref_ptr<AnimationTime> _at;
  };
  
} // namespace orsaOSG

#endif // _ORSA_OSG_HUD_H_
