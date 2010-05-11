#ifndef _ORSA_OSG_VIZ_
#define _ORSA_OSG_VIZ_

#include <QHash>

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <osg/Group>
#include <osg/Node>
#include <osg/PositionAttitudeTransform>
#include <osg/Switch>

#include <orsa/bodygroup.h>

#include <orsaOSG/AnimationTime.h>

namespace orsaOSG {
  
    class Viz : public osg::Referenced {
    public:
        Viz(orsa::BodyGroup *, 
            const double & timeMultiplier = 1.0);
    protected:
        virtual ~Viz() { }
    public:
        osg::Group * createRoot();
    public:    
        osg::Node * createHUD();
    private:
        osg::Group * oldCode();
    protected:
        osg::ref_ptr<orsa::BodyGroup> _bg;
        const double _timeMultiplier;
    
    public:
        osg::PositionAttitudeTransform * getBodyPositionTransform(const orsa::Body * b) const {
            return _bodyPositionTransform[b].get();
        }
    public:
        osg::PositionAttitudeTransform * getBodyAttitudeTransform(const orsa::Body * b) const {
            return _bodyAttitudeTransform[b].get();
        }
    private: 
        QHash < const orsa::Body * , osg::ref_ptr<osg::Switch> >                    _bodySwitch;
        QHash < const orsa::Body * , osg::ref_ptr<osg::PositionAttitudeTransform> > _bodyPositionTransform;
        QHash < const orsa::Body * , osg::ref_ptr<osg::PositionAttitudeTransform> > _bodyAttitudeTransform;
    
    private:
        // QHash < const orsa::Body * , osg::ref_ptr<osg::Switch> > _orbitSwitch;
    
    public:
        osg::ref_ptr<orsaOSG::AnimationTime> _at;
    };
  
} // namespace orsaOSG

#endif // _ORSA_OSG_VIZ_
