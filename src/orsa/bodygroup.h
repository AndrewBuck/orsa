#ifndef _ORSA_BODYGROUP_
#define _ORSA_BODYGROUP_

// #include <list>
#include <string>
#include <vector>

#include <orsa/body.h>
#include <orsa/interaction.h>
#include <orsa/interval.h>

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <QHash>

namespace orsa {
  
    // class Attitude;
  
    class BodyGroup : public osg::Referenced {
    public:    
        BodyGroup();
    protected:
        virtual ~BodyGroup();
    
    public:
        virtual bool setName(const std::string &);
        virtual const std::string & getName() const;
    protected:
        std::string _name;
    
    public:
        void clear();
        void clearIntegration(const bool restoreInitialConditions=true);
    
    public:
        bool addBody(const Body *);
        bool removeBody(const Body *);
        inline size_t size() const { return _b_list.size(); } 
        
    public:
        typedef std::vector < osg::ref_ptr<const Body> > BodyList;
    protected:
        BodyList _b_list;
        // NOTE: some (all?) non-const access member should be protected
        // and the access allowet to friend classes only
    public:
        const BodyList & getBodyList() const { 
            /* for (unsigned int k=0; k<_b_list.size(); ++k) {
               ORSA_DEBUG("member[%04i] = %x [%s]",k,_b_list[k].get(),_b_list[k]->getName().c_str());
               }
            */
            //
            return _b_list; 
        }
        // BodyList & getBodyList() { return _b_list; }
    
    public:
        const Body * getBody(const orsa::Body::BodyID & bodyID) const;
        const Body * getBody(const std::string & bodyName) const;
        
        /* 
           public:
           bool setInteractionType(const orsa::Interaction::InteractionType);
           public:
           orsa::Interaction::InteractionType getInteractionType() const;
        */
    public:
        const orsa::Interaction * getInteraction() const {
            return _itr.get();
        }	
    protected:
        osg::ref_ptr<orsa::Interaction> _itr;
    
        /* 
           public:
           class TRV {
           public:
           TRV() : tmp(false) { }
           public:
           Time   t;
           Vector r;
           Vector v;
           public:	
           // osg::ref_ptr<orsa::RotationAngles> rot;
           public:
           bool tmp;
           public:
           inline bool operator == (const TRV & rhs) const {
           return (t == rhs.t);
           }
           public:
           inline bool operator != (const TRV & rhs) const {
           return (t != rhs.t);
           }
           public:
           inline bool operator < (const TRV & rhs) const {
           return (t < rhs.t);
           }
           public:
           inline bool operator > (const TRV & rhs) const {
           return (t > rhs.t);
           }
           public:
           inline bool operator <= (const TRV & rhs) const {
           return (t <= rhs.t);
           }
           public:
           inline bool operator >= (const TRV & rhs) const {
           return (t >= rhs.t);
           }
           };
        */
        //
        /* 
           public:
           class TRVB : public TRV {
           public:
           osg::ref_ptr<Body> b;
           };
        */
        //
    public:	
        typedef orsa::Interval<orsa::IBPS> BodyInterval;
    public:
        // typedef std::map<const orsa::Body *, osg::ref_ptr<orsa::Interval<TRV> > > BodyIntervalMap;
        // typedef QHash<const orsa::Body *, osg::ref_ptr<orsa::Interval<TRV> > > BodyIntervalMap;
        typedef QHash<const orsa::Body *, osg::ref_ptr<BodyInterval> > BodyIntervalMap;
    protected:
        // std::map<osg::ref_ptr<orsa::BodyGroup>, BodyInterval> _interval;
    protected:
        BodyIntervalMap _b_interval;
    public:
        // const osg::ref_ptr<orsa::Interval<TRV> > & getBodyInterval(const orsa::Body *) const;
        // orsa::Interval<TRV> * getBodyInterval(const orsa::Body *);
    public:
        const BodyInterval * getBodyInterval(const orsa::Body *) const;
        BodyInterval * getBodyInterval(const orsa::Body *);
        // BodyInterval * getBodyInterval(const orsa::Body *);
        /* 
           inline orsa::Interval<TRV> * getBodyInterval(osg::ref_ptr<const orsa::Body> b) {
           return getBodyInterval(b.get());
           }
        */
    
    public:
        orsa::Time getTimeInterval(const orsa::Body * b) const {
            osg::ref_ptr<const BodyInterval> i = getBodyInterval(b);
            return (i->max().time.getRef() - i->min().time.getRef());
        }
    
    public:
        orsa::Time getLongestTimeInterval() const {
            orsa::Time _t(0);
            BodyGroup::BodyList::const_iterator _b_it = getBodyList().begin();
            while (_b_it != getBodyList().end()) {
                osg::ref_ptr<const BodyInterval> i = getBodyInterval((*_b_it).get());
                const Time _tmp_t = (i->max().time.getRef() - i->min().time.getRef());
                if (_tmp_t > _t) {
                    _t = _tmp_t;
                }
                ++_b_it;
            }
            return _t;
        }
        
    public:
        bool insertIBPS(const orsa::IBPS & ibps,
                        const orsa::Body * b,
                        const bool         onlyIfExtending,
                        const bool         replaceIfDouble);
        
    public:
        bool getIBPS(orsa::IBPS       & ibps,
                     const orsa::Body * b,
                     const Time       & t) const;
        
    public:
        bool getInterpolatedIBPS(orsa::IBPS       & ibps,
                                 const orsa::Body * b,
                                 const Time       & t) const;
    
    public:
        bool getClosestIBPS(orsa::IBPS       & ibps,
                            const orsa::Body * b,
                            const Time       & t) const;
    
    public:
        bool getInterpolatedPosition(Vector     & position,
                                     const Body * b,
                                     const Time & t) const;
    
    public:
        bool getInterpolatedVelocity(Vector     & velocity,
                                     const Body * b,
                                     const Time & t) const;
    
    public:
        bool getInterpolatedPosVel(Vector     & position,
                                   Vector     & velocity,
                                   const Body * b,
                                   const Time & t) const;

    public:
        bool getInterpolatedMass(double     & mass,
                                 const Body * b,
                                 const Time & t) const;
        
    public:
        bool haveDynamicBodies(const bool massive_only) const;
        
    public:
        bool getClosestCommonTime(orsa::Time       & t,
                                  const orsa::Time & refTime,
                                  const bool         massive_only) const;
    
    public:
        bool getCommonInterval(orsa::Time & start,
                               orsa::Time & stop,
                               const bool   massive_only) const;
    
    public:
        bool getGlobalInterval(orsa::Time & start,
                               orsa::Time & stop,
                               const bool   massive_only) const;
    
    public:
        void centerOfMassPosVel(orsa::Vector & r,
                                orsa::Vector & v,
                                const orsa::Time & t) const;
  
    public:
        //! total energy with respect to the barycenter
        double totalEnergy(const orsa::Time & t) const;
  
    public:
        //! total angular momentum with respect to the barycenter
        orsa::Vector totalAngularMomentum(const orsa::Time & t) const;
    };
  
} // namespace orsa

#endif // _ORSA_BODYGROUP_
