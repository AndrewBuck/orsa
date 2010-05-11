#ifndef _ORSA_OSG_BIO_CALLBACK_H_
#define _ORSA_OSG_BIO_CALLBACK_H_

#include <osg/Node>
#include <osg/NodeCallback>

#include <orsa/bodygroup.h>

// #include <orsaOSG/animation_timing.h>
#include <orsaOSG/AnimationTime.h>

namespace orsaOSG {
  
    class bioCallback : public osg::NodeCallback {
    public:
        bioCallback(orsa::BodyGroup  * bg,
                    const orsa::Body * b,
                    orsaOSG::AnimationTime * at);
    
    public:
        void operator () (osg::Node *, osg::NodeVisitor *);
    
    protected:
        osg::ref_ptr<orsa::BodyGroup> _bg;
        osg::ref_ptr<const orsa::Body> _b;
    
    protected:
        osg::ref_ptr<AnimationTime> _at;
    };
  
} // namespace orsaOSG

#endif // _ORSA_OSG_BIO_CALLBACK_H_
