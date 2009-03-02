#include <orsa/interaction.h>

// #include <orsa/attitude.h>
#include <orsa/bodygroup.h>
// #include <orsa/legendre.h>
#include <orsa/paul.h>
#include <orsa/print.h>
#include <orsa/util.h>

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

#include <tbb/tick_count.h>

using namespace orsa;

class IBPSCacheTBB;

namespace orsa {
  class BodyPair {
  public:
    // orsa::BodyGroup::BodyList::const_iterator ref_b_it;
    // orsa::BodyGroup::BodyList::const_iterator     b_it;
    //
    // const orsa::Body * ref_b;
    // const orsa::Body * b;
    unsigned int ref_b_index;
    unsigned int     b_index;
    //
    const IBPS * ref_b_ibps;
    const IBPS *     b_ibps;
  };	
} // namespace orsa;

#ifdef ORSA_USE_TBB
typedef std::vector< BodyPair, tbb::cache_aligned_allocator<BodyPair> > BodyPairVector;
#else
typedef std::vector< BodyPair> BodyPairVector;
#endif

class InteractionTBB {
public:
  InteractionTBB(orsa::BodyGroup  * bodyGroup,
		 const orsa::Time & time,
		 const BodyPairVector & bpv) : 
    interaction(new orsa::Interaction),
    bg(bodyGroup), 
    t(time), 
    bodyPairVector(bpv) {
    // each component is Vector(0,0,0) by default
    // a.resize(bg->largestBodyID()+1);
    a.resize(bg->getBodyList().size());
  }
public:
  InteractionTBB(InteractionTBB & i, tbb::split) : 
    interaction(i.interaction),
    bg(i.bg),
    t(i.t),
    bodyPairVector(i.bodyPairVector) {
    // each component is Vector(0,0,0) by default
    // a.resize(bg->largestBodyID()+1);
    a.resize(bg->getBodyList().size());
  } 
protected:
  osg::ref_ptr<orsa::Interaction> interaction;
protected:
  orsa::BodyGroup * bg;
  orsa::Time        t;
protected:
  const BodyPairVector & bodyPairVector;
public:
  void operator()(const tbb::blocked_range<unsigned int> & range) {
    // ORSA_DEBUG("range: %04i - %04i",range.begin(),range.end());
    orsa::Vector a_ref_b, a_b;
    for (unsigned int k=range.begin(); k!=range.end(); ++k) {
      const orsa::BodyPair & bp = bodyPairVector[k];
      if (interaction->bodyPairAccelerationTerm(a_ref_b,
						a_b,
						bg,
						&bp,
						t)) {
	// a[bp.ref_b->id()] += a_ref_b;
	// a[bp.b->id()]     += a_b;
	a[bp.ref_b_index] += a_ref_b;
	a[bp.b_index]     += a_b;
      } else {
	ORSA_DEBUG("problems...");
      }	
    }
  }
public:
  /* 
     void operator()(const tbb::blocked_range<unsigned int> & range) {
     // ORSA_DEBUG("range: %04i - %04i",range.begin(),range.end());
     orsa::Vector a_ref_b, a_b;
     const orsa::Double G = orsa::Unit::instance()->getG();
     for (unsigned int k=range.begin(); k!=range.end(); ++k) {
     
     const orsa::BodyPair & bp = bodyPairVector[k];
     
     const orsa::Double m_ref_b = bp.ref_b_ibps->inertial->mass();
     const orsa::Double m_b     = bp.b_ibps->inertial->mass();
     
     const orsa::Body * ref_b = bp.ref_b;
     const orsa::Body * b     = bp.b;
     
     orsa::Vector _d =
     bp.b_ibps->translational->position() - 
     bp.ref_b_ibps->translational->position();
     
     const Double _l = _d.length();
     
     _d /= (_l*_l*_l);
     
     const orsa::Vector accTerm = _d;
     
     a_ref_b =   G * m_b     * accTerm;
     a_b     = - G * m_ref_b * accTerm;
     
     a[bp.ref_b->id()] += a_ref_b;
     a[bp.b->id()]     += a_b;
     
     }
     }
  */
public:
  void join(InteractionTBB & rhs) {
    // const tbb::tick_count tick_in = tbb::tick_count::now();
    for (unsigned int k=0; k<a.size(); ++k) {
      a[k] += rhs.a[k];
    }
    // const tbb::tick_count tick_out = tbb::tick_count::now();
    // ORSA_DEBUG("join call took %g [s]",(tick_out-tick_in).seconds());
  }
public:
  orsa::Interaction::InteractionVector a;
  // orsa::Interaction::VectorHash N;
};

/* 
   class IBPSCache {
   public:
   IBPSCache(orsa::BodyGroup  * bodyGroup,
   const orsa::Time & time) :
   bg(bodyGroup),
   t(time) {
   cachedIBPS.resize(bg->largestBodyID()+1);
   } 
   public:
   bool getInterpolatedIBPS(const orsa::IBPS * * ibps_ptr,
   const orsa::Body * b) {
   if (!cachedIBPS[b->id()].isSet()) {
   // ORSA_DEBUG("computing...");
   orsa::IBPS ibps_tmp;
   if (bg->getInterpolatedIBPS(ibps_tmp,
   b,
   t)) {
   cachedIBPS[b->id()] = ibps_tmp;
   (*ibps_ptr) = cachedIBPS[b->id()].getPtr();
   return true;
   } else {
   return false;
   }
   } else {
   // ORSA_DEBUG("cached...");
   (*ibps_ptr) = cachedIBPS[b->id()].getPtr();
   return true;
   }
   }
   protected:
   orsa::BodyGroup * bg;
   const orsa::Time t;
   protected:
   typedef std::vector< orsa::Cache<orsa::IBPS> > BodyIDIBPS;
   BodyIDIBPS cachedIBPS;
   };
*/

class IBPSCacheTBB { 
public:
#ifdef ORSA_USE_TBB
  typedef std::vector< orsa::IBPS, tbb::cache_aligned_allocator<orsa::IBPS> > BodyIDIBPS;
#else 
  typedef std::vector< orsa::IBPS> BodyIDIBPS;
#endif
  
public:
  IBPSCacheTBB(orsa::BodyGroup  * bodyGroup,
	       const orsa::Time & time,
	       BodyIDIBPS * cIBPS) :
    bg(bodyGroup),
    t(time),
    cachedIBPS(cIBPS) {
    // cachedIBPS->resize(bg->largestBodyID()+1);
    // bodyPtrVector.resize(bg->largestBodyID()+1);
    cachedIBPS->resize(bg->getBodyList().size());
    // bodyPtrVector.resize(bg->getBodyList().size());
    //
    /* for (unsigned int k=0; k<bodyPtrVector.size(); ++k) {
       bodyPtrVector[k] = 0;
       }
    */
    //
    /* BodyGroup::BodyList bl = bg->getBodyList();
       BodyGroup::BodyList::const_iterator b_it = bl.begin();
       while (b_it != bl.end()) { 
       bodyPtrVector[(*b_it)->id()] = (*b_it).get();
       ++b_it;
       } 
    */
    // debug
    /* for (unsigned int k=0; k<bodyPtrVector.size(); ++k) {
       ORSA_DEBUG("bodyPtrVector[%i] = %x",
       k, bodyPtrVector[k]);
       }
    */
  } 
public:
  // copy constructor important when using tbb::parallel_for
  IBPSCacheTBB(const IBPSCacheTBB & c) : 
    bg(c.bg),
    t(c.t),
    cachedIBPS(c.cachedIBPS) // , bodyPtrVector(c.bodyPtrVector)
  { }
public:
  /* 
     ~IBPSCacheTBB() {
     // for debugging only...
     BodyGroup::BodyList bl = bg->getBodyList();
     BodyGroup::BodyList::const_iterator b_it = bl.begin();
     while (b_it != bl.end()) { 
     const unsigned int id = (*b_it)->id();
     ORSA_DEBUG("id: %i   ptr: %x   translational.get(): %x",
     id, 
     &cachedIBPS[id],
     cachedIBPS[id].translational.get());
     ++b_it;
     } 
     }
  */
public:
  void operator() (const tbb::blocked_range<unsigned int> & range) const {
    // ORSA_DEBUG("range: %04i - %04i",range.begin(),range.end());
    for (unsigned int k=range.begin(); k!=range.end(); ++k) {
      // if (bodyPtrVector[k] != 0) {
      if (bg->getInterpolatedIBPS((*cachedIBPS)[k],
				  // bodyPtrVector[k],
				  bg->getBodyList()[k].get(),
				  t)) {
	/* if (!(*cachedIBPS)[k].translational.get()) {
	   ORSA_DEBUG("problems...");
	   }
	*/
      } else {
	ORSA_DEBUG("problems...");
      }
      // }
    }
  }
protected:
  orsa::BodyGroup * bg;
  const orsa::Time t;
public:
  BodyIDIBPS * cachedIBPS;
  /* 
     protected:
     #ifdef ORSA_USE_TBB
     std::vector<const orsa::Body *, tbb::cache_aligned_allocator<const orsa::Body *> > bodyPtrVector;
     #else 
     std::vector<const orsa::Body *> bodyPtrVector;
     #endif
  */
};


Interaction::Interaction() : osg::Referenced(true) { 
  dummyPaulMoment = new PaulMoment(0);
  dummyPaulMoment->setM(one(),0,0,0);
  dummyPaulMoment->setCenterOfMass(orsa::Vector(0,0,0));
  dummyPaulMoment->setInertiaMoment(orsa::Matrix::identity());
}	

bool Interaction::acceleration(InteractionVector & a,  
			       orsa::BodyGroup   * bg,
			       const orsa::Time  & t) const {
  
  // ORSA_DEBUG("called...");
  
  // const tbb::tick_count tick_in = tbb::tick_count::now();
  
  BodyGroup::BodyList bl = bg->getBodyList();
  
  /* 
     {
     // reset accel. vector
     BodyGroup::BodyList::const_iterator ref_b_it = bl.begin();
     while (ref_b_it != bl.end()) {
     a[(*ref_b_it).get()] = orsa::Vector(0,0,0);
     ++ref_b_it;
     }
     }
  */
  //
  // a.resize(bg->largestBodyID()+1);
  
  IBPSCacheTBB::BodyIDIBPS cachedIBPS;
  
  IBPSCacheTBB cachedIBPSTBB(bg,t,&cachedIBPS);
  //
  tbb::parallel_for(tbb::blocked_range<unsigned int>(0,bl.size()),
		    cachedIBPSTBB,
		    tbb::auto_partitioner());
  
  // std::vector<BodyPair> bodyPairVector;
  BodyPairVector bodyPairVector;
  //
  /* 
     {
     // create body pairs list
     // ORSA_DEBUG("pairs creation started...");
     // orsa::Double m_ref_b, m_b;
     BodyGroup::BodyList::const_iterator ref_b_it = bl.begin();
     while (ref_b_it != bl.end()) {
     if (!(*ref_b_it)->alive(t)) {
     ++ref_b_it;
     continue;
     }
     if ((*ref_b_it)->getInitialConditions().translational.get()) {
     if (!((*ref_b_it)->getInitialConditions().translational->dynamic())) {
     ++ref_b_it;
     continue;
     }
     } else {
     ++ref_b_it;
     continue;
     }
     BodyGroup::BodyList::const_iterator b_it = bl.begin();
     while (b_it != ref_b_it) {
     if ((*ref_b_it)->nonInteractingGroup && 
     (*b_it)->nonInteractingGroup) {
     ++b_it;
     continue;
     }
     if (!(*b_it)->alive(t)) {
     ++b_it;
     continue;
     }
     //
     BodyPair bp;
     bp.ref_b      = (*ref_b_it).get();
     bp.b          = (*b_it).get();
     // bp.ref_b_ibps = cachedIBPSTBB.cachedIBPS[bp.ref_b->id()].getPtr();
     // bp.b_ibps     = cachedIBPSTBB.cachedIBPS[bp.b->id()].getPtr();
     bp.ref_b_ibps = &(cachedIBPS[bp.ref_b->id()]);
     bp.b_ibps     = &(cachedIBPS[bp.b->id()]);
     bodyPairVector.push_back(bp);
     ++b_it;
     }
     ++ref_b_it;
     }
     // ORSA_DEBUG("pairs creation done.");
     }
  */
  //
  {
    for (unsigned int j=0; j<bl.size(); ++j) {
      const orsa::Body * ref_b = bl[j].get();
      if (!ref_b->alive(t)) {
     	continue;
      }
      if (ref_b->getInitialConditions().translational.get()) {
	if (!(ref_b->getInitialConditions().translational->dynamic())) {
	  continue;
	}
      } else {
     	continue;
      }
      for (unsigned int k=0; k<j; ++k) {
	const orsa::Body * b = bl[k].get();
	if (ref_b->nonInteractingGroup && 
	    b->nonInteractingGroup) {
	  /* ORSA_DEBUG("skipping two bodies belonging to nonInteractingGroup: [%s] & [%s]",
	     (*ref_b_it)->getName().c_str(),
	     (*b_it)->getName().c_str());
	  */
	  continue;
	}
	if (!b->alive(t)) {
	  continue;
	}
	BodyPair bp;
	// bp.ref_b      = (*ref_b_it).get();
	// bp.b          = (*b_it).get();
	bp.ref_b_index = j;
	bp.b_index     = k;
	// bp.ref_b_ibps = cachedIBPSTBB.cachedIBPS[bp.ref_b->id()].getPtr();
	// bp.b_ibps     = cachedIBPSTBB.cachedIBPS[bp.b->id()].getPtr();
	bp.ref_b_ibps = &(cachedIBPS[j]);
	bp.b_ibps     = &(cachedIBPS[k]);
	bodyPairVector.push_back(bp);
      }
    }
  }
  
  // ORSA_DEBUG("bodyPairVector.size() %i   bg->getBodyList().size() %i",bodyPairVector.size(),bg->getBodyList().size());
  
  // const tbb::tick_count tick_pairs = tbb::tick_count::now();
  // ORSA_DEBUG("pairs done, so far call took %g [s]",(tick_pairs-tick_in).seconds());
  
  /* 
     IBPSCache cachedIBPS(bg,t);
     {
     // create body pairs list
     // ORSA_DEBUG("pairs creation started...");
     BodyGroup::BodyList::const_iterator ref_b_it = bl.begin();
     while (ref_b_it != bl.end()) {
     BodyGroup::BodyList::const_iterator b_it = bl.begin();
     while (b_it != ref_b_it) {
     if ((*ref_b_it)->nonInteractingGroup && 
     (*b_it)->nonInteractingGroup) {
     ++b_it;
     continue;
     }
     BodyPair bp;
     bp.ref_b = (*ref_b_it).get();
     bp.b     = (*b_it).get();
     if (!(cachedIBPS.getInterpolatedIBPS(&bp.ref_b_ibps,
     bp.ref_b))) {
     ORSA_DEBUG("problems...");
     return false;
     }
     if (!(cachedIBPS.getInterpolatedIBPS(&bp.b_ibps,
     bp.b))) {
     ORSA_DEBUG("problems...");
     return false;
     }
     bodyPairVector.push_back(bp);
     ++b_it;
     }
     ++ref_b_it;
     }
     // ORSA_DEBUG("pairs creation done.");
     }
  */
  
  // ORSA_DEBUG("bodyPairVector.size(): %i",bodyPairVector.size());
  
  /* 
     a.resize(bg->largestBodyID()+1);
     {
     // reset accel. vector
     BodyGroup::BodyList::const_iterator ref_b_it = bl.begin();
     while (ref_b_it != bl.end()) {
     a[(*ref_b_it)->id()] = orsa::Vector(0,0,0);
     ++ref_b_it;
     }
     }
     orsa::Vector a_ref_b, a_b;
     for (unsigned int k=0; k<bodyPairVector.size(); ++k) {
     const BodyPair bp = bodyPairVector[k];
     if (bodyPairAccelerationTerm(a_ref_b,
     a_b,
     bg,
     &bp,
     t)) {
     a[bp.ref_b->id()] += a_ref_b;
     a[bp.b->id()]     += a_b;
     } else {
     ORSA_DEBUG("problems...");
     }	
     }
  */
  //
  InteractionTBB interactionTBB(bg,t,bodyPairVector);
  // parallel_reduce(tbb::blocked_range<unsigned int>(0,bodyPairVector.size(),(100+bodyPairVector.size()/2)),interactionTBB);
  parallel_reduce(tbb::blocked_range<unsigned int>(0,bodyPairVector.size()),
		  interactionTBB,
		  tbb::auto_partitioner());
  // save computed value
  a = interactionTBB.a;
  
  // const tbb::tick_count tick_out = tbb::tick_count::now();
  // ORSA_DEBUG("done.... call took %g [s]",(tick_out-tick_in).seconds());
  
  return true;
}

bool Interaction::bodyPairAccelerationTerm(orsa::Vector     & a_ref_b,
					   orsa::Vector     & a_b,
					   orsa::BodyGroup  * bg,
					   const BodyPair   * bp,
					   const orsa::Time & t) const {
  
  // const tbb::tick_count tick_in = tbb::tick_count::now();
  
#warning "keep checking if reset needed..."
  // a_ref_b = orsa::Vector(0,0,0);
  // a_b     = orsa::Vector(0,0,0);
  
  /* ORSA_DEBUG("bp->ref_b_ibps: %x   bp->ref_b_ibps->inertial.get(): %x   [%s]",
     bp->ref_b_ibps,bp->ref_b_ibps->inertial.get(),bp->ref_b->getName().c_str());
  */
  
  const orsa::Double m_ref_b = bp->ref_b_ibps->inertial->mass();
  
  const orsa::Double m_b     = bp->b_ibps->inertial->mass();
  
  // orsa::Double m_ref_b, m_b;
  
  BodyGroup::BodyList bl = bg->getBodyList();
  
  // BodyGroup::BodyList::const_iterator ref_b_it = bp->ref_b_it;
  
  // const orsa::Body * ref_b = bp->ref_b;
  const orsa::Body * ref_b = bl[bp->ref_b_index].get();
  
  /* if (!bg->getInterpolatedMass(m_ref_b,ref_b,t)) {
     ORSA_DEBUG("problems...");
     }
  */
  
  // #warning "make this check more generic/accurate..."
  
  /* 
     if (ref_b->getInitialConditions().translational.get()) {
     if (!(ref_b->getInitialConditions().translational->dynamic())) {
     return true;
     }
     } else {
     return true;
     }
  */
  
  /* 
     if (!ref_b->alive(t)) {
     return true;
     }
  */
  
  orsa::Vector thrust;
  //
  if (ref_b->propulsion.get()) {
    if (!dependsOnVelocity()) {
      ORSA_DEBUG("PROBLEMS: when using a Propulsion, in general you need to set orsa::Interaction::dependsOnVelocity() as TRUE");
    }
    thrust = ref_b->propulsion->getThrust(t);
  }
  //
  /* The problem with using FrenetSerret is that positions and velocities
     are not well defined at the time t, because the relative IBPS is temporary.
  */
  /* 
     if (ref_b->propulsion.get()) {
     const orsa::Vector propulsionThrust = ref_b->propulsion->getThrust(t);
     if (propulsionThrust.lengthSquared() > orsa::epsilon()*orsa::epsilon()) {
     const orsa::BodyGroup::BodyInterval * bi =
     bg->getBodyInterval(ref_b);
     const orsa::Time dt = 
     std::min(orsa::Time(0,0,0,5,0),
     std::min(bi->max().time.getRef()-t,
     t-bi->min().time.getRef()));
     // REMEMBER: now thrust is in FrenetSerret components
     orsa::Vector T, N, B;
     ORSA_DEBUG("bi->min().time.getRef().tmp: %i",bi->min().tmp);
     ORSA_DEBUG("bi->max().time.getRef().tmp: %i",bi->max().tmp);
     if (t > bi->min().time.getRef()) {
     FrenetSerret(ref_b, bg,
     t,
     -dt,
     T, N, B);
     } else if (t < bi->max().time.getRef()) {
     FrenetSerret(ref_b, bg,
     t,
     dt,
     T, N, B);
     } else {
     ORSA_DEBUG("interval smaller than dt");
     ORSA_DEBUG("--BODY-INTERVAL-TIME-RANGE--");
     print(bi->min().time.getRef());
     print(bi->max().time.getRef());
     ORSA_DEBUG("call time:");
     print(t);
     //
     T = orsa::Vector(1,0,0);
     N = orsa::Vector(0,1,0);
     B = orsa::Vector(0,0,1);
     }
     thrust = 
     propulsionThrust.getX()*T +	
     propulsionThrust.getY()*N +
     propulsionThrust.getZ()*B;
     }
     }
  */
  
  // ORSA_DEBUG("acc thrust");
  // orsa::print(thrust);
  
  // propulsion
  // #warning "propulsion code disabled for now"
  /* 
     if (ref_b->getPropulsion() != 0) {
     osg::ref_ptr<PropulsionEvent> pe = ref_b->getPropulsion()->getPropulsionEvent(t);
     // if (pe.thrustMagnitude.getRef() > orsa::zero()) {
     // ORSA_DEBUG("using propulsion!! time: %Ff",t.asDouble().get_mpf_t());
     
     // IMPORTANT!!!!
     // ORSA_DEBUG("a=F/m and F=thrust, so I have to divide by the instantaneous mass...");
     
     // a[(*ref_b_it).get()] += pe.thrustMagnitude.getRef() * pe.thrustDirection();
     
     const orsa::Double bodyMass = ref_b->getMass() - ref_b->getPropulsion()->massLost(t);
     
     if (ref_b->getPropulsion()->thrustToMass.isSet()) {
     ORSA_DEBUG("t: %Ff   bodyMass: %Fg",
     t.asDouble().get_mpf_t(),
     bodyMass.get_mpf_t());
     }
     
     if (bodyMass > orsa::zero()) {
     a_ref_b += pe->thrust.getRef() / bodyMass;
     } else {
     ORSA_DEBUG("problem: non-positive mass!! (m=%Fe)",bodyMass.get_mpf_t());
     }
     
     // } else {
     // ORSA_DEBUG("not using propulsion, time: %Ff",t.asDouble().get_mpf_t());
     // }
     }
  */
  
  /* 
     if (!(bg->getInterpolatedIBPS(ref_b_ibps,
     (*ref_b_it).get(),
     t))) {
     ORSA_DEBUG("problems...");
     return false;
     }
  */
  
  // BodyGroup::BodyList::const_iterator b_it = bp->b_it;
  // const orsa::Body * b = bp->b;
  const orsa::Body * b = bl[bp->b_index].get();
  
  /*
    if (!bg->getInterpolatedMass(m_b,b,t)) {
    ORSA_DEBUG("problems...");
    }
  */
  
  /* 
     ORSA_DEBUG("ref_b: [%s]   b: [%s]",
     (*ref_b_it)->getName().c_str(),
     (*b_it)->getName().c_str());
  */
  
  // already checked...
  /* if (b == ref_b) {
     return true;
     }
  */
  
  // this check already included in pairs creation code
  /* if (ref_b->nonInteractingGroup && 
      b->nonInteractingGroup) {
      return true;
      }
  */
  
  /* 
     ORSA_DEBUG("ref_b: [%s]   b: [%s]",
     (*ref_b_it)->getName().c_str(),
     (*b_it)->getName().c_str());
  */
  
  // cannot do this!
  /* 
     if ((*b_it)->getInitialConditions().translational.get()) {
     if (!((*b_it)->getInitialConditions().translational->dynamic())) {
     ++b_it;
     continue;
     }
     }
  */
  
  /* 
     ORSA_DEBUG("ref_b: [%s]   b: [%s]",
     (*ref_b_it)->getName().c_str(),
     (*b_it)->getName().c_str());
  */
  
  /* if (!b->alive(t)) {
     return true;
     }
  */
  /* 
     ORSA_DEBUG("ref_b: [%s]   b: [%s]",
     (*ref_b_it)->getName().c_str(),
     (*b_it)->getName().c_str());
  */
  
  /* if ((b->getMass() == zero()) && 
     (ref_b->getMass() == zero())) {
     return true;
     }
  */
  
  if ((m_b == zero()) && 
      (m_ref_b == zero())) {
    return true;
  }
  
  /* 
     ORSA_DEBUG("ref_b: [%s]   b: [%s]",
     (*ref_b_it)->getName().c_str(),
     (*b_it)->getName().c_str());
  */
  
  /* 
     if (!(bg->getInterpolatedIBPS(b_ibps,
     (*b_it).get(),
     t))) {
     ORSA_DEBUG("problems...");
     return false;
     }
  */
  
  /* 
     ORSA_DEBUG("ref_b: [%s]   b: [%s]",
     (*ref_b_it)->getName().c_str(),
     (*b_it)->getName().c_str());
  */
  
  if (ref_b->getPaulMoment() || 
      b->getPaulMoment()) {
    
    // ORSA_DEBUG("--MARK--");
    
    osg::ref_ptr<const PaulMoment> ref_b_pm = 
      (ref_b->getPaulMoment()) ? 
      (ref_b->getPaulMoment()) :
      (dummyPaulMoment.get());
    
    osg::ref_ptr<const PaulMoment> b_pm = 
      (b->getPaulMoment()) ? 
      (b->getPaulMoment()) :
      (dummyPaulMoment.get());
    
    /* 
       if (b_pm.get() == dummyPaulMoment.get()) {
       ORSA_DEBUG("using dummy PaulMoment, body [%s]",(*b_it)->getName().c_str());
       }
    */
    
    /*
      if (bp->b_ibps->translational.get() &&
      bp->ref_b_ibps->translational.get()) {
    */
    
    /* osg::ref_ptr<orsa::Attitude> ref_b_attitude = 
       new orsa::BodyAttitude(ref_b,bg);
    */
    //
    const orsa::Matrix ref_b_l2g = orsa::localToGlobal(ref_b,bg,t);
    const orsa::Matrix ref_b_g2l = orsa::globalToLocal(ref_b,bg,t);
    
    /* osg::ref_ptr<orsa::Attitude> b_attitude = 
       new orsa::BodyAttitude(b,bg);
    */
    //
    const orsa::Matrix b_l2g = orsa::localToGlobal(b,bg,t);
    const orsa::Matrix b_g2l = orsa::globalToLocal(b,bg,t);
    
    const orsa::Vector R =
      (bp->b_ibps->translational->position() + b_l2g*b_pm->getCenterOfMass()) - 
      (bp->ref_b_ibps->translational->position() + ref_b_l2g*ref_b_pm->getCenterOfMass());
    
    /* 
       const orsa::Vector accTerm =
       Paul::gravitationalForce(ref_b_pm.get(),
       new orsa::BodyAttitude((*ref_b_it).get(),bg),
       b_pm.get(),
       new orsa::BodyAttitude((*b_it).get(),bg),
       R,
       t);
    */
    //
    const orsa::Vector accTerm =
      Paul::gravitationalForce(ref_b_pm.get(),
			       ref_b_g2l,
			       b_pm.get(),
			       b_g2l,
			       R);
    
    // a[(*ref_b_it).get()] += (*b_it)->getMu()     * accTerm;
    // a[    (*b_it).get()] -= (*ref_b_it)->getMu() * accTerm;
    
    // a_ref_b += b->getMu()     * accTerm;
    // a_b     -= ref_b->getMu() * accTerm;
    
    a_ref_b =   orsa::Unit::instance()->getG() * m_b     * accTerm;
    a_b     = - orsa::Unit::instance()->getG() * m_ref_b * accTerm;
    
    /*   
	 } else {
	 if (!bp->ref_b_ibps->translational.get()) ORSA_DEBUG("problems: [%s].translational.get() = %x",
	 ref_b->getName().c_str(),
	 bp->ref_b_ibps->translational.get());
	 if (!bp->b_ibps->translational.get()) ORSA_DEBUG("problems: [%s].translational.get() = %x",
	 b->getName().c_str(),
	 bp->b_ibps->translational.get());
	 }
    */
    
  } else {
    
    // ORSA_DEBUG("--MARK--");
    
    /* 
       if (bp->b_ibps->translational.get() &&
       bp->ref_b_ibps->translational.get()) {
    */
    
    orsa::Vector _d =
      bp->b_ibps->translational->position() - 
      bp->ref_b_ibps->translational->position();
    
    const Double _l = _d.length();
    
    /* 
       if (_l > epsilon()) {
    */
    
    _d /= (_l*_l*_l);
    
    if (ref_b->betaSun == b) {
      const orsa::Vector accTerm = (one() - ref_b->beta.getRef()) * _d;
      // a[(*ref_b_it).get()] += (*b_it)->getMu()     * accTerm;
      // a[    (*b_it).get()] -= (*ref_b_it)->getMu() * accTerm;
      // a_ref_b += b->getMu()     * accTerm;
      // a_b     -= ref_b->getMu() * accTerm;
      a_ref_b =   orsa::Unit::instance()->getG() * m_b     * accTerm;
      a_b     = - orsa::Unit::instance()->getG() * m_ref_b * accTerm;
    } else {
      const orsa::Vector accTerm = _d;
      // a[(*ref_b_it).get()] += (*b_it)->getMu()     * accTerm;
      // a[    (*b_it).get()] -= (*ref_b_it)->getMu() * accTerm;
      // a_ref_b += b->getMu()     * accTerm;
      // a_b     -= ref_b->getMu() * accTerm;
      a_ref_b =   orsa::Unit::instance()->getG() * m_b     * accTerm;
      a_b     = - orsa::Unit::instance()->getG() * m_ref_b * accTerm;
    }
    
    /* 
       } else {
       
       ORSA_WARNING("skipping: zero distance between bodies [%s] and [%s]",
       ref_b->getName().c_str(),
       b->getName().c_str());
       
       }
    */
    
    /* 
       } else {
       if (!bp->ref_b_ibps->translational.get()) ORSA_DEBUG("problems: [%s].translational.get() = %x",
       ref_b->getName().c_str(),
       bp->ref_b_ibps->translational.get());
       if (!bp->b_ibps->translational.get()) ORSA_DEBUG("problems: [%s].translational.get() = %x",
       b->getName().c_str(),
       bp->b_ibps->translational.get());  
       }
    */
    
  }
  
  // const tbb::tick_count tick_out = tbb::tick_count::now();
  /* ORSA_DEBUG("done.... call took %g [s]   bodies: [%s] & [%s]",
     (tick_out-tick_in).seconds(),
     ref_b->getName().c_str(),
     b->getName().c_str());
  */
  
  // thrust
  if (ref_b->propulsion.get()) {
    if (m_ref_b < orsa::epsilon()) {
      ORSA_DEBUG("Propulsion problems: non-positive mass...");
    } else {
      a_ref_b += thrust / m_ref_b;
    }
  }
  
  return true;
}

bool Interaction::torque(InteractionVector & N,  
			 orsa::BodyGroup   * bg,
			 const orsa::Time  & t) const {
  
  // ORSA_DEBUG("called...");
  
  BodyGroup::BodyList bl = bg->getBodyList();
  
  {
    // reset torque vector
    // N.resize(bg->largestBodyID()+1);
    N.resize(bl.size());
    // BodyGroup::BodyList::const_iterator ref_b_it = bl.begin();
    // while (ref_b_it != bl.end()) {
    // N[(*ref_b_it)->id()] = orsa::Vector(0,0,0);
    // ++ref_b_it;
    // }
    for (unsigned int k=0; k<bl.size(); ++k) {
      N[k] = orsa::Vector(0,0,0);
    }	
  }
  
  {
    // main loop
    orsa::Double m_ref_b, m_b;
    
    // BodyGroup::BodyList::const_iterator ref_b_it = bl.begin();
    // while (ref_b_it != bl.end()) {
    
    for (unsigned int j=0; j<bl.size(); ++j) {
      
      const orsa::Body * ref_b = bl[j].get();
      
      // should this be here??
      /* if ((*ref_b_it)->getInitialConditions().translational.get()) {
	 if (!((*ref_b_it)->getInitialConditions().translational->dynamic())) {
	 ORSA_DEBUG("--MARK--");
	 ++ref_b_it;
	 continue;
	 }
	 }
      */
      
      if (ref_b->getInitialConditions().rotational.get()) {
	if (!(ref_b->getInitialConditions().rotational->dynamic())) {
	  // ++ref_b_it;
	  continue;
	}
      } else {
	// ++ref_b_it;
	continue;
      }
      
      if (!ref_b->alive(t)) {
       	// ++ref_b_it;
	continue;
      }
      
      // propulsion torque?
      
      IBPS ref_b_ibps;
      
      if (!(bg->getInterpolatedIBPS(ref_b_ibps,
				    ref_b,
				    t))) {
	ORSA_DEBUG("problems...");
	return false;
      }
      
      osg::ref_ptr<const PaulMoment> ref_b_pm = 
	(ref_b->getPaulMoment()) ? 
	(ref_b->getPaulMoment()) :
	(dummyPaulMoment.get());
      
      if (!bg->getInterpolatedMass(m_ref_b,ref_b,t)) {
	ORSA_DEBUG("problems...");
      }
      
      // BodyGroup::BodyList::const_iterator b_it = bl.begin();
      // while (b_it != bl.end()) {
      
      for (unsigned int k=0; k<j; ++k) {
	
	const orsa::Body * b = bl[k].get();
	
	/* 
	   if ((*b_it).get() == 
	   (*ref_b_it).get()) {
	   break;
	   }
	*/
	
	// cannot do this!
	/* 
	   if ((*b_it)->getInitialConditions().translational.get()) {
	   if (!((*b_it)->getInitialConditions().translational->dynamic())) {
	   ++b_it;
	   continue;
	   }
	   }
	*/
	//
	/* 
	   if ((*b_it)->getInitialConditions().rotational.get()) {
	   if (!((*b_it)->getInitialConditions().rotational->dynamic())) {
	   ++b_it;
	   continue;
	   }
	   }
	*/
	
	/* 
	   if ((*ref_b_it)->nonInteractingGroupID.isSet() && 
	   (*b_it)->nonInteractingGroupID.isSet()) {
	   if ((*ref_b_it)->nonInteractingGroupID.get() == 
	   (*b_it)->nonInteractingGroupID.get()) {
	   ++b_it;
	   continue;
	   }
	   }
	*/
	//
	if (ref_b->nonInteractingGroup && 
	    b->nonInteractingGroup) {
	  // ++b_it;
	  continue;
	}	
	
	if (!b->alive(t)) {
	  // ++b_it;
	  continue;
	}
	
	if (!bg->getInterpolatedMass(m_b,b,t)) {
	  ORSA_DEBUG("problems...");
	}
	
	if ((m_b == zero()) && 
	    (m_ref_b == zero())) {
	  // ++b_it;
	  continue;
	}
	
	IBPS b_ibps;
	
	if (!(bg->getInterpolatedIBPS(b_ibps,
				      b,
				      t))) {
	  ORSA_DEBUG("problems...");
	  return false;
	}
	
	if (ref_b->getPaulMoment() || 
	    b->getPaulMoment()) {
	  
	  osg::ref_ptr<const PaulMoment> b_pm = 
	    (b->getPaulMoment()) ? 
	    (b->getPaulMoment()) :
	    (dummyPaulMoment.get());
	  
	  if ( (b_ibps.translational.get() &&
		ref_b_ibps.translational.get()) && 
	       (b_ibps.rotational.get() &&
		ref_b_ibps.rotational.get()) ) {
	    
	    /* 
	       const orsa::Vector R =
	       b_ibps.translational->position() - 
	       ref_b_ibps.translational->position();
	    */
	    
	    /* osg::ref_ptr<orsa::Attitude> ref_b_attitude = 
	       new orsa::BodyAttitude((*ref_b_it).get(),bg);
	    */
	    //
	    const orsa::Matrix ref_b_l2g = orsa::localToGlobal(ref_b,bg,t);
	    const orsa::Matrix ref_b_g2l = orsa::globalToLocal(ref_b,bg,t);
	    
	    /* osg::ref_ptr<orsa::Attitude> b_attitude = 
	       new orsa::BodyAttitude((*b_it).get(),bg);
	    */
	    //
	    const orsa::Matrix b_l2g = orsa::localToGlobal(b,bg,t);
	    const orsa::Matrix b_g2l = orsa::globalToLocal(b,bg,t);
	    
	    const orsa::Vector R =
	      (b_ibps.translational->position() + b_l2g*b_pm->getCenterOfMass()) - 
	      (ref_b_ibps.translational->position() + ref_b_l2g*ref_b_pm->getCenterOfMass());
	    
	    const orsa::Vector torqueTerm =
	      Paul::gravitationalTorque(ref_b_pm.get(),
					ref_b_g2l,
					b_pm.get(),
					b_g2l,
					R);
	    
	    // N[(*ref_b_it)->id()] += orsa::Unit::instance()->getG() * m_b * torqueTerm;
	    // N[    (*b_it)->id()] -= orsa::Unit::instance()->getG() * m_ref_b * torqueTerm;
	    //
	    N[j] += orsa::Unit::instance()->getG() * m_b     * torqueTerm;
	    N[k] -= orsa::Unit::instance()->getG() * m_ref_b * torqueTerm;
	    
	  } else {
	    if (!ref_b_ibps.translational.get()) ORSA_DEBUG("problems: [%s].translational.get() = %x",
							    ref_b->getName().c_str(),
							    ref_b_ibps.translational.get());
	    if (!b_ibps.translational.get()) ORSA_DEBUG("problems: [%s].translational.get() = %x",
							b->getName().c_str(),
							b_ibps.translational.get());
	  }
	  
	}
	
	//++b_it;
      }
      
      // ++ref_b_it;
    }
  }
  
  // ORSA_DEBUG("done.");
  
  return true;
}
