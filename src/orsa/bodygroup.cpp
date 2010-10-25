#include <orsa/bodygroup.h>

#include <orsa/print.h>
#include <orsa/paul.h>
#include <orsa/slerp.h>
#include <orsa/spline.h>
#include <orsa/util.h>
#include <orsa/vector.h>

#include <algorithm>
#include <map>

using namespace std;
using namespace orsa;

BodyGroup::BodyGroup() : Referenced(true) {
    _itr = new Interaction;
}

BodyGroup::~BodyGroup() {
  
}

void BodyGroup::clear() {
    _name.clear();
    _b_list.clear();
    _b_interval.clear();
}

void BodyGroup::clearIntegration(const bool restoreInitialConditions) {
    BodyList::const_iterator it = _b_list.begin();
    while (it != _b_list.end()) {
        getBodyInterval(*it)->reset();
        if (restoreInitialConditions) {
            getBodyInterval(*it)->insert((*it)->getInitialConditions(),false,false);
        }
        ++it;
    }
}

bool BodyGroup::setName(const std::string & s) {
    _name = s;
    return true;
}

const std::string & BodyGroup::getName() const {
    return _name;
}

bool BodyGroup::addBody(const Body * b) {
    
    const BodyList::const_iterator _b = find(_b_list.begin(),_b_list.end(),b);
    
    if (_b != _b_list.end()) {
        ORSA_ERROR("Body already included in BodyGroup.");
        return false;
    } 
    
    _b_list.push_back(b);
    
    if (b->getInitialConditions().dynamic()) {
        const IBPS ibps = b->getInitialConditions();
        insertIBPS(ibps, b, false, false);
    }
    
    return true;
}

bool BodyGroup::removeBody(const Body * b) {
    
    BodyList::iterator _bl = find(_b_list.begin(),_b_list.end(),b);
    if (_bl != _b_list.end()) {
        _bl = _b_list.erase(_bl);
        BodyIntervalMap::iterator _bi = _b_interval.find(b);
        if (_bi != _b_interval.end()) {
            _bi = _b_interval.erase(_bi);
            return true;
        } else {
            ORSA_DEBUG("Body found in list but not in map");
            return false;
        }
    } else {
        ORSA_ERROR("Body not included in BodyGroup");
        return false;
    }
    return false;
}

bool BodyGroup::insertIBPS(const orsa::IBPS & ibps,
                           const orsa::Body * b,
                           const bool         onlyIfExtending,
                           const bool         replaceIfDouble) {
    return getBodyInterval(b)->insert(ibps,onlyIfExtending,replaceIfDouble);
}

bool BodyGroup::getIBPS(orsa::IBPS       & ibps,
                        const orsa::Body * b,
                        const orsa::Time & t) const {
    // if (b->getBodyPosVelCallback()) {
    if (!(b->getInitialConditions().dynamic())) { 
        // trv.t = t;
        // b->getBodyPosVelCallback()->getPosVel(t,trv.r,trv.v);
        ibps = b->getInitialConditions();
        ibps.update(t);
        return true;
    } else {
        // BodyGroup::TRV _trv1, _trv2;
        IBPS ibps1, ibps2;
        // orsa::Interval<BodyGroup::TRV> * _b_interval = getBodyInterval(b);
        osg::ref_ptr<const orsa::BodyGroup::BodyInterval> bi = getBodyInterval(b);
        if (bi.get() == 0) {
            return false;
        }
        // trv.t = t;
        ibps.time = t;
        // if (_b_interval->getSubInterval(trv, _trv1, _trv2)) {
        if (bi->getSubInterval(ibps,ibps1,ibps2)) {
            if ((t == ibps1.time.getRef()) && 
                (t == ibps2.time.getRef())) {
                ibps = ibps1;
                ibps.update(t);
                return true;
            } else {
                // ORSA_ERROR("point not present in interval...");
                return false;
            }
        } else {
            // ORSA_ERROR("problems encountered in getSubInterval(...) call");
            return false;
        }	
    }
}

bool BodyGroup::getClosestIBPS(orsa::IBPS       & ibps,
                               const orsa::Body * b,
                               const orsa::Time & t) const {
    if (!(b->getInitialConditions().dynamic())) { 
        ibps = b->getInitialConditions();
        ibps.update(t);
        return true;
    } else {
        // BodyGroup::TRV _trv1, _trv2;
        IBPS ibps1, ibps2;
        // orsa::Interval<BodyGroup::TRV> * _b_interval = getBodyInterval(b);
        osg::ref_ptr<const orsa::BodyGroup::BodyInterval> bi = getBodyInterval(b);
        // trv.t = t;
        ibps.time = t; 
        // if (_b_interval->getSubInterval(trv, _trv1, _trv2)) {
        if (bi->getSubInterval(ibps,ibps1,ibps2)) {
            // if (fabs((_trv1.t-t).get_d()) < fabs((_trv2.t-t).get_d())) {
            if (fabs((ibps1.time.getRef()-t).get_d()) < 
                fabs((ibps2.time.getRef()-t).get_d())) {
                ibps = ibps1;
            } else {
                ibps = ibps2;
            }  
            ibps.update(t);
            return true;
        } else {
            return false;
        }	
    }
}

bool BodyGroup::getInterpolatedIBPS(orsa::IBPS       & ibps,
                                    const orsa::Body * b,
                                    const orsa::Time & t) const {
  
    if (!b->alive(t)) {
        // ORSA_DEBUG("out, body [%s]",b->getName().c_str());
        return false;
    }
  
    /* 
       ORSA_DEBUG("body: [%s]  bodies: %i",
       b->getName().c_str(),
       size());
    */
  
    if (b->getInitialConditions().dynamic()) {
        IBPS ibps1, ibps2;
        osg::ref_ptr<const orsa::BodyGroup::BodyInterval> bi = getBodyInterval(b);
        ibps.time = t; 
        if (!bi) {
            ORSA_DEBUG("problems... body: [%s]",b->getName().c_str());
            return false;
        }
        if (bi->getSubInterval(ibps,ibps1,ibps2)) {
            if ((t == ibps1.time.getRef()) && (t == ibps2.time.getRef())) {
                ibps = ibps1;	
                //
                // VERY IMPORTANT call to update(t);
                ibps.update(t);
                // ORSA_DEBUG("out, body [%s]",b->getName().c_str());
                return true;
            } else {
                if ( (t < ibps1.time.getRef()) || 
                     (t > ibps2.time.getRef())) {
                    ORSA_WARNING("outside boundaries!!");
                    getClosestIBPS(ibps,b,t);
                    ORSA_DEBUG("out, body [%s]",b->getName().c_str());
                    return false;
                }
                if (ibps1.time.getRef() == ibps2.time.getRef()) {
                    ibps = ibps1;
                    // ORSA_DEBUG("out, body [%s]",b->getName().c_str());
                    ibps.update(t);
                    return true;	  
                } else {
                    // to copy pointers...
                    ibps = ibps1;	
                    // set time again...
                    ibps.time = t; 
	  
	  
                    if (ibps.inertial.get()) {
                        if (ibps.inertial->dynamic()) {
                            // linear interpolation (enough?)
                            const double m1 = ibps1.inertial->mass();
                            const double m2 = ibps2.inertial->mass();
                            //
                            const orsa::Time & t1 = ibps1.time.getRef();
                            const orsa::Time & t2 = ibps2.time.getRef();
                            //
                            const double mt = ((m2-m1)/(t2-t1).get_d())*(t-t1).get_d();
                            ibps.inertial->setMass(mt);
                        } else {
                            ibps.inertial->update(t);
                        }
                    }
	  
	  
                    if (ibps.translational.get()) {
                        if (ibps.translational->dynamic()) {
                            osg::ref_ptr<orsa::PhysicalSpline<orsa::Vector> > s = new orsa::PhysicalSpline<orsa::Vector>;
                            //
                            if (s->set(ibps1.translational->position(),
                                       ibps1.translational->velocity(),
                                       ibps1.time.getRef(),
                                       ibps2.translational->position(),
                                       ibps2.translational->velocity(),
                                       ibps2.time.getRef())) {
                                orsa::Vector r,v;
                                if (s->get(r,v,t)) {
                                    ibps.translational->setPosition(r);
                                    ibps.translational->setVelocity(v);
		  
                                    /* 
                                       ORSA_DEBUG("spline, r: %f   r1: %f   r2: %f",
                                       r.length(),
                                       ibps1.translational->position().length(),
                                       ibps2.translational->position().length());
                                    */
		  
                                }
                            } else {
                                ORSA_DEBUG("problems...");
                                return false;
                            }
                        } else {
                            ibps.translational->update(t);
                        }
                    } 
	  
	  
                    if (ibps.rotational.get()) {
	    
                        if (ibps.rotational->dynamic()) {
	      
                            orsa::Quaternion sQ, sQDot;
	      
                            osg::ref_ptr<orsa::Slerp> slerp = new Slerp;
	      
                            if (slerp->set(ibps1.rotational->getQ(),
                                           ibps1.time.getRef(),
                                           ibps2.rotational->getQ(),
                                           ibps2.time.getRef())) {
                                if (!slerp->get(sQ,t)) {
                                    ORSA_DEBUG("problems...");
                                    return false;
                                }
                            } else {
                                ORSA_DEBUG("problems...");
                                return false;   
                            }	
	      
                            sQ = unitQuaternion(sQ);
	      
	      
                            ibps.rotational->set(sQ,
                                                 ibps1.rotational->getOmega());
	      
                        } else {
                            ibps.rotational->update(t); 
                        }
                    }
	  
	  
                    // ORSA_DEBUG("out, body [%s]",b->getName().c_str());
                    ibps.update(t);
                    return true;
                }
            }	
        } else {
            ORSA_DEBUG("problems... body: [%s] dynamic: %i time: [below]",
                       b->getName().c_str(),
                       b->getInitialConditions().dynamic());
            orsa::print(t);
      
            ORSA_DEBUG("body initial conditions time: [below]");
            orsa::print(b->getInitialConditions().time.getRef());
      
            // ORSA_DEBUG("..out..");
            return false;
        }	
    } else {
        ibps = b->getInitialConditions();
        ibps.update(t);
        // ORSA_DEBUG("out, body [%s]",b->getName().c_str());
        return true;
    } 
}

bool BodyGroup::getInterpolatedPosition(Vector     & position,
                                        const Body * b,
                                        const Time & t) const {
    orsa::IBPS ibps;
    if (getInterpolatedIBPS(ibps,b,t)) {
        if (ibps.translational.get()) {
            position = ibps.translational->position();
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

bool BodyGroup::getInterpolatedVelocity(Vector     & velocity,
                                        const Body * b,
                                        const Time & t) const {
    orsa::IBPS ibps;
    if (getInterpolatedIBPS(ibps,b,t)) {
        if (ibps.translational.get()) {
            velocity = ibps.translational->velocity();
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

bool BodyGroup::getInterpolatedPosVel(Vector     & position,
                                      Vector     & velocity,
                                      const Body * b,
                                      const Time & t) const {
    orsa::IBPS ibps;
    if (getInterpolatedIBPS(ibps,b,t)) {
        if (ibps.translational.get()) {
            position = ibps.translational->position();
            velocity = ibps.translational->velocity();
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

bool BodyGroup::getInterpolatedMass(double     & mass,
                                    const Body * b,
                                    const Time & t) const {
    orsa::IBPS ibps;
    if (getInterpolatedIBPS(ibps,b,t)) {
        if (ibps.inertial.get()) {
            mass = ibps.inertial->mass();
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

class TimeBool {
public:
    TimeBool(const orsa::Time & t,
             const orsa::Time & t0) : 
        time(t), 
        refTime(t0) { }
    // check(true) { }
public:
    bool operator < (const TimeBool & rhs) const {
        if (refTime != rhs.refTime) {
            ORSA_ERROR("inconsistent reference time");
        }
        return (abs((time-refTime).getMuSec()) < 
                abs((rhs.time-rhs.refTime).getMuSec()));
    }
public:
    const orsa::Time time, refTime;
public:
    // bool check;
};

bool BodyGroup::getClosestCommonTime(orsa::Time       & t,
                                     const orsa::Time & refTime,
                                     const bool         massive_only) const {
  
    typedef std::list<TimeBool> TimeList;
    TimeList timeList;
  
    BodyGroup::BodyList::const_iterator b_it = getBodyList().begin();
    while (b_it != getBodyList().end()) { 
    
        if (massive_only && (*b_it)->getInitialConditions().inertial->mass() == 0.0) { 
            ++b_it;
            continue;
        }
    
        if (!(*b_it)->getInitialConditions().dynamic()) { 
            ++b_it;
            continue;
        }
    
        if (timeList.size() == 0) {
      
            const osg::ref_ptr<const orsa::BodyGroup::BodyInterval> i_b_it = getBodyInterval((*b_it).get());
            const orsa::BodyGroup::BodyInterval::DataType & interval_data = i_b_it->getData();
            orsa::BodyGroup::BodyInterval::DataType::const_iterator interval_data_it = interval_data.begin();
            while (interval_data_it != interval_data.end()) {
                timeList.push_back(TimeBool((*interval_data_it).time.getRef(),
                                            refTime));
                ++interval_data_it;
            }
      
            // sort for candidate closest time...
            timeList.sort();
      
            TimeList::iterator tl_it = timeList.begin();
            while (tl_it != timeList.end()) {
	
                IBPS ibps;
	
                bool found=true;
	
                unsigned int trials=0;
	
                BodyGroup::BodyList::const_iterator b2_it = b_it; ++b2_it;
                while (b2_it != getBodyList().end()) { 
	  
                    ++trials;
	  
                    if (b2_it == b_it) {
                        ++b2_it;
                        continue;
                    }
	  
                    if (massive_only && (*b2_it)->getInitialConditions().inertial->mass() == 0.0) { 
                        ++b2_it;
                        continue;
                    }
	  
                    if (!(*b2_it)->getInitialConditions().dynamic()) { 
                        ++b2_it;
                        continue;
                    }
	  
                    if (!getIBPS(ibps,(*b2_it).get(),(*tl_it).time)) {
                        found=false;
                        break;
                    }
	  
                    ++b2_it;
                }
	
                if (found) {
                    // ORSA_DEBUG("FOUND in %i trials",trials);
                    t = (*tl_it).time;
                    return true;
                }
	
                ++tl_it;
            }	
      
            return false;
        }
        
        ++b_it;
    }
    
    return false;
}

bool BodyGroup::haveDynamicBodies(const bool massive_only) const {
    
    BodyIntervalMap::const_iterator _bi = _b_interval.begin();
    
    while (_bi != _b_interval.end()) {
        
        if (massive_only && _bi.key()->getInitialConditions().inertial->mass() == 0.0) { 
            ++_bi;
            continue;
        }
        
        if (!(_bi.key()->getInitialConditions().dynamic())) { 
            ++_bi;
            continue;
        }
        
        return true;
        
        ++_bi;
    }
    
    return false;
}

bool BodyGroup::getCommonInterval(orsa::Time & start, orsa::Time & stop, const bool massive_only) const {
    
    Cache<Time> _min, _max;
    
    BodyIntervalMap::const_iterator _bi = _b_interval.begin();
    
    while (_bi != _b_interval.end()) {
        
        if (massive_only && _bi.key()->getInitialConditions().inertial->mass() == 0.0) { 
            ++_bi;
            continue;
        }
        
        if (!(_bi.key()->getInitialConditions().dynamic())) { 
            ++_bi;
            continue;
        }
        
        const osg::ref_ptr<orsa::BodyGroup::BodyInterval> _interval = _bi.value();
        
        /* ORSA_DEBUG("considering body [%s] id: %i  interval: ",
           _bi.key()->getName().c_str(),
           _bi.key()->id());
           orsa::print(_interval->min().time.getRef());
           orsa::print(_interval->max().time.getRef());
        */
        
        if ( (!_min.isSet()) && 
             (!_max.isSet()) ) {
            _min.set(_interval->min().time.getRef());
            _max.set(_interval->max().time.getRef());
        } else if ( (_min.isSet()) && 
                    (_max.isSet()) ) {
            if ( (_interval->min().time.getRef() > _max.getRef()) ||
                 (_interval->max().time.getRef() < _min.getRef()) ) {
                // ORSA_DEBUG("no common interval found");
                return false;
            } else {
                if (_interval->min().time.getRef() > _min.getRef()) {
                    _min.set(_interval->min().time.getRef());
                }
                if (_interval->max().time.getRef() < _max.getRef()) {
                    _max.set(_interval->max().time.getRef());
                }
            }
        } else {
            // one only is set, quite odd...
            ORSA_ERROR("logic error: we shouldn't be here");
        }
        
        ++_bi;
    }
    
    if ( (_min.isSet()) && 
         (_max.isSet()) && 
         (_min.getRef() <= _max.getRef()) ) {
        start = _min.getRef();
        stop  = _max.getRef();
        return true;
    }
    
    return false;  
}

bool BodyGroup::getGlobalInterval(orsa::Time & start, orsa::Time & stop, const bool massive_only) const {
    Cache<Time> _min, _max;
    BodyIntervalMap::const_iterator _bi = _b_interval.begin();
    while (_bi != _b_interval.end()) {
        if (massive_only && (_bi.key()->getInitialConditions().inertial->mass() == 0.0)) {
            ++_bi;
            continue;
        }	
        /* 
           if (_bi.key()->getBodyPosVelCallback() != 0) {
           ++_bi;
           continue;
           }
        */
        //
        if (!(_bi.key()->getInitialConditions().dynamic())) { 
            ++_bi;
            continue;
        }
    
        const osg::ref_ptr<orsa::BodyGroup::BodyInterval> _interval = _bi.value();
        if ( (!_min.isSet()) && 
             (!_max.isSet()) ) {
            // if ( (_bi.key()->alive(_interval->min().t)) &&
            // (_bi.key()->alive(_interval->max().t)) ) {
            _min.set(_interval->min().time.getRef());
            _max.set(_interval->max().time.getRef());
            // ORSA_DEBUG("min: %f   max: %f",_min.getRef().get_d(),_max.getRef().get_d());
            // ORSA_DEBUG("//1//");
            // } 
        } else if ( (_min.isSet()) && 
                    (_max.isSet()) ) {
            if (_interval->min().time.getRef() < _min.getRef()) {
                _min.set(_interval->min().time.getRef());
                // ORSA_DEBUG("//2//");
            }
            if (_interval->max().time.getRef() > _max.getRef()) {
                _max.set(_interval->max().time.getRef());
                // ORSA_DEBUG("//3//");
            }
        } else {
            // one only is set, quite odd...
            ORSA_ERROR("logic error: we shouldn't be here");
        }
        ++_bi;
    }
    if ( (_min.isSet()) && 
         (_max.isSet()) && 
         (_min.getRef() <= _max.getRef()) ) {
        start = _min.getRef();
        stop  = _max.getRef();
        /* 
           ORSA_DEBUG("start: %f   stop: %f",
           start.get_d(),
           stop.get_d());
        */
        return true;
    } else {
        ORSA_WARNING("start and stop are unset...");
        return false;  
    }
}

const orsa::BodyGroup::BodyInterval * BodyGroup::getBodyInterval(const orsa::Body * b) const {
    // ORSA_DEBUG("body: [%s]",b->getName().c_str());
    if (_b_interval[b].get()) {
        return (_b_interval[b].get());
    } else {
        ORSA_ERROR("cannot use const version of getBodyInterval");
        return (0);
    }
}

orsa::BodyGroup::BodyInterval * BodyGroup::getBodyInterval(const orsa::Body * b) {
    if (_b_interval[b].get()) {
        return (_b_interval[b].get());
    } else {
        /* 
           osg::ref_ptr<orsa::BodyGroup::BodyInterval> bi = new orsa::BodyGroup::BodyInterval;
           bi->enableDataStoring();
           _b_interval[b] = bi;
           // ORSA_DEBUG("creating new interval for body [%s] ... done.",b->getName().c_str());
           return (bi.get());
        */
        //
        // ORSA_DEBUG("creating new interval for body [%s]",b->getName().c_str());
        orsa::BodyGroup::BodyInterval * bi = new orsa::BodyGroup::BodyInterval;
        bi->enableDataStoring();
        _b_interval[b] = bi;
        // ORSA_DEBUG("creating new interval for body [%s] ... done.",b->getName().c_str());
        return (bi);
    }
}

const Body * BodyGroup::getBody(const orsa::Body::BodyID & bodyID) const {
    BodyList::const_iterator it = _b_list.begin();
    while (it != _b_list.end()) {
        if ((*it)->id() == bodyID) {
            return (*it).get();
        }
        ++it;
    }
    return 0;
}

const Body * BodyGroup::getBody(const std::string & bodyName) const {
    // returns the "first" body found with this name, or zero
    BodyList::const_iterator it = _b_list.begin();
    while (it != _b_list.end()) {
        if ((*it)->getName() == bodyName) {
            return (*it).get();
        }
        ++it;
    }
    return 0;
}

/* 
   orsa::Vector BodyGroup::centerOfMassPosition(const orsa::Time & t) {
   double sum_mb(0);
   orsa::Vector sum_mb_rb(0,0,0);
   //
   double mb;
   orsa::Vector rb;
   BodyList::const_iterator it = _b_list.begin();
   while (it != _b_list.end()) {
   mb = (*it)->getMass();
   if (mb > 0) {
   ORSA_DEBUG("more checks needed for extended bodies");
   if (getInterpolatedPosition(rb,(*it).get(),t)) {
   sum_mb    += mb;
   sum_mb_rb += mb*rb;
   }
   }
   ++it;
   }
   return (sum_mb_rb/sum_mb);
   }
*/

void BodyGroup::centerOfMassPosVel(orsa::Vector     & r,
                                   orsa::Vector     & v,
                                   const orsa::Time & t) const {
  
    double sum_mb(0);
    orsa::Vector sum_mb_rb(0,0,0);
    orsa::Vector sum_mb_vb(0,0,0);
  
    double mb;
    orsa::Vector rb;
    orsa::Vector vb;
  
    BodyList::const_iterator it = _b_list.begin();
    while (it != _b_list.end()) {
        // mb = (*it)->getMass();
        if (getInterpolatedMass(mb,(*it).get(),t)) {
            if (mb > 0) {
                // ORSA_DEBUG("more checks needed for extended bodies");
#warning "more checks needed for extended bodies"
	
                if (getInterpolatedPosVel(rb,vb,(*it).get(),t)) {
                    sum_mb    += mb;
                    sum_mb_rb += mb*rb;
                    sum_mb_vb += mb*vb;
                }
            }
        }
        ++it;
    }
  
    r = sum_mb_rb/sum_mb;
    v = sum_mb_vb/sum_mb;
}

double BodyGroup::totalEnergy(const orsa::Time & t) const {
  
    orsa::Vector r,v;
    centerOfMassPosVel(r,v,t);
    const orsa::Vector rcm = r;
    const orsa::Vector vcm = v;
  
    osg::ref_ptr<PaulMoment> dummy_pm = new PaulMoment(0);
    dummy_pm->setM(1,0,0,0);
    // dummy_pm->setCenterOfMass(orsa::Vector(0,0,0));
    // dummy_pm->setInertiaMoment(orsa::Matrix::identity());
  
    double E(0);
  
    // kinetic energy contribution
    {
        double m;
        BodyList::const_iterator it = _b_list.begin();
        while (it != _b_list.end()) {
      
            if (!(*it)->alive(t)) {
                ++it;
                continue;
            }
      
            if (!getInterpolatedMass(m,(*it).get(),t)) {
                ORSA_DEBUG("problems...");
            }	
      
            if (getInterpolatedPosVel(r,v,(*it).get(),t)) {
                v -= vcm;
                E += m*v*v/2;
            }
      
            ++it;
        }
    }
  
    // rotational energy contribution 
    {
        double m;
        BodyList::const_iterator it = _b_list.begin();
        while (it != _b_list.end()) {
      
            if (!(*it)->alive(t)) {
                ++it;
                continue;
            }
      
            if (!getInterpolatedMass(m,(*it).get(),t)) {
                ORSA_DEBUG("problems...");
            }	
      
            /* osg::ref_ptr<orsa::Attitude> attitude =
               new BodyAttitude((*it).get(),this);
            */
            //
            const orsa::Matrix g2l = orsa::globalToLocal((*it).get(),this,t);
      
            /* 
               if ((*it)->getPaulMoment()) {
               orsa::IBPS ibps;
               if (getInterpolatedIBPS(ibps,(*it).get(),t)) {
               if (ibps.rotational.get()) {
               const orsa::Vector omega = 
               g2l * ibps.rotational->getOmega();
               E += 
               omega * 
               (*it)->getPaulMoment()->getInertiaMoment() * m *
               omega / 2;
               }
               }
               }
            */
            //
            // if ((*it)->getPaulMoment()) {
            orsa::IBPS ibps;
            if (getInterpolatedIBPS(ibps,(*it).get(),t)) {
                if (ibps.rotational.get()) {
                    if (ibps.inertial.get()) {
                        const orsa::Vector omega = 
                            g2l * ibps.rotational->getOmega();
                        E += 
                            omega * 
                            ibps.inertial->inertiaMatrix() * m *
                            omega / 2;
                    }
                }
            }
      
            ++it;
        }
    }
  
    // potential energy contribution
    {
        orsa::Vector r1, v1, r2, v2;
    
        BodyList::const_iterator it1 = _b_list.begin();
        while (it1 != _b_list.end()) {
      
            if (!(*it1)->alive(t)) {
                ++it1;
                continue;
            }
      
            if (!getInterpolatedPosVel(r1,v1,(*it1).get(),t)) {
                ORSA_DEBUG("problems...");
            }
      
            /* osg::ref_ptr<orsa::Attitude> b1_attitude = 
               new orsa::BodyAttitude((*it1).get(),this);
            */
            //
            const orsa::Matrix b1_l2g = orsa::localToGlobal((*it1).get(),this,t);
            const orsa::Matrix b1_g2l = orsa::globalToLocal((*it1).get(),this,t);
      
            orsa::IBPS ibps1;
            if (!getInterpolatedIBPS(ibps1,(*it1).get(),t)) {
                ORSA_DEBUG("problems...");
            }
      
            osg::ref_ptr<const PaulMoment> b1_pm = 
                (ibps1.inertial->paulMoment()) ? 
                (ibps1.inertial->paulMoment()) :
                (dummy_pm.get());
      
            BodyList::const_iterator it2 = _b_list.begin();
            while (it2 != _b_list.end()) {
	
                if (!(*it2)->alive(t)) {
                    ++it2;
                    continue;
                }
	
                if ((*it1).get() == (*it2).get()) {
                    break;
                }
	
                if (!getInterpolatedPosVel(r2,v2,(*it2).get(),t)) {
                    ORSA_DEBUG("problems...");
                }
	
                /* osg::ref_ptr<orsa::Attitude> b2_attitude = 
                   new orsa::BodyAttitude((*it2).get(),this);
                */
                //
                const orsa::Matrix b2_l2g = orsa::localToGlobal((*it2).get(),this,t);
                const orsa::Matrix b2_g2l = orsa::globalToLocal((*it2).get(),this,t);
	
                orsa::IBPS ibps2;
                if (!getInterpolatedIBPS(ibps2,(*it2).get(),t)) {
                    ORSA_DEBUG("problems...");
                }
	
                osg::ref_ptr<const PaulMoment> b2_pm = 
                    (ibps2.inertial->paulMoment()) ? 
                    (ibps2.inertial->paulMoment()) :
                    (dummy_pm.get());
	
                // dr = r2-r1
                /* const orsa::Vector dr =
                   (r2 + b2_l2g*b2_pm->getCenterOfMass()) - 
                   (r1 + b1_l2g*b1_pm->getCenterOfMass());
                */
                //
                const orsa::Vector dr = r2 - r1;
	
                double m1, m2;
	
                if ( (!getInterpolatedMass(m1,(*it1).get(),t)) ||
                     (!getInterpolatedMass(m2,(*it2).get(),t)) ) {
                    ORSA_DEBUG("problems...");
                }	
	
                E -= 
                    orsa::Unit::G() *
                    m1 * 
                    m2 *
                    Paul::gravitationalPotential(b1_pm.get(),
                                                 b1_g2l,
                                                 b2_pm.get(),
                                                 b2_g2l,
                                                 dr);
	
                ++it2;
            }
      
            ++it1;
        }
    }
  
    return E;
}

orsa::Vector BodyGroup::totalAngularMomentum(const orsa::Time & t) const {
  
    orsa::Vector r, v;
    centerOfMassPosVel(r, v, t);
    const orsa::Vector rcm = r;
    const orsa::Vector vcm = v;
  
    osg::ref_ptr<PaulMoment> dummy_pm = new PaulMoment(0);
    dummy_pm->setM(1,0,0,0);
    // dummy_pm->setCenterOfMass(orsa::Vector(0,0,0));
    // dummy_pm->setInertiaMoment(orsa::Matrix::identity());
  
    orsa::Vector L(0,0,0);
  
    // contribution from body's rotation
    {
        double m;
        BodyList::const_iterator it = _b_list.begin();
        while (it != _b_list.end()) {
      
            if (!(*it)->alive(t)) {
                ++it;
                continue;
            }
      
            if (!getInterpolatedMass(m,(*it).get(),t)) {
                ORSA_DEBUG("problems...");
            }	
      
            /* osg::ref_ptr<orsa::Attitude> attitude =
               new BodyAttitude((*it).get(),this);
            */
            //
            const orsa::Matrix g2l = orsa::globalToLocal((*it).get(),this,t);
      
            // if ((*it)->getPaulMoment()) {
            orsa::IBPS ibps;
            if (getInterpolatedIBPS(ibps,(*it).get(),t)) {
                if (ibps.rotational.get()) {
                    if (ibps.inertial.get()) {
                        const orsa::Vector omega = 
                            g2l * ibps.rotational->getOmega();
                        L += 
                            ibps.inertial->inertiaMatrix() * m *
                            omega;
                    }
                }
            }
      
            ++it;
        }
    }
  
    // contribution from rotation about the barycenter
    {
        double m;
        BodyList::const_iterator it = _b_list.begin();
        while (it != _b_list.end()) {
      
            if (!(*it)->alive(t)) {
                ++it;
                continue;
            }
      
            if (!getInterpolatedMass(m,(*it).get(),t)) {
                ORSA_DEBUG("problems...");
            }	
      
            if (!getInterpolatedPosVel(r,v,(*it).get(),t)) {
                ORSA_DEBUG("problems...");
            }
      
            /* osg::ref_ptr<orsa::Attitude> b_attitude = 
               new orsa::BodyAttitude((*it).get(),this);
            */
            //
            const orsa::Matrix l2g = orsa::localToGlobal((*it).get(),this,t);
      
            orsa::IBPS ibps;
            if (!getInterpolatedIBPS(ibps,(*it).get(),t)) {
                ORSA_DEBUG("problems...");
            }
      
            osg::ref_ptr<const PaulMoment> b_pm = 
                (ibps.inertial->paulMoment()) ? 
                (ibps.inertial->paulMoment()) :
                (dummy_pm.get());
      
            // const orsa::Vector rb = r + l2g*b_pm->getCenterOfMass() - rcm;
            const orsa::Vector rb = r - rcm;
            const orsa::Vector vb = v + vcm;
      
            L += m * orsa::externalProduct(rb,vb);
      
            ++it;
        }
    }
  
    return L;
}
