#include <orsaSPICE/spiceBodyInterpolatedPosVelCallback.h>

#include <orsaSPICE/spice.h>

#include <orsa/spline.h>

using namespace orsa;
using namespace orsaSPICE;

/* 
   SpiceBodyInterpolatedPosVelCallback::SpiceBodyInterpolatedPosVelCallback(const std::string & name,
   const orsa::Time  & start,
   const orsa::Time  & stop,
   const orsa::Time  & samplingPeriod) :
   BodyPosVelCallback() {
   
   _interval = new Interval<BodyGroup::TRV>;
   _interval->enableDataStoring();
   
   BodyGroup::TRV trv;
   trv.tmp = false;
   orsa::Time t = start;
   while (t <= stop) {
   // ORSA_DEBUG("t: %Ff",t.asDouble().get_mpf_t());
   trv.t = t;
   SPICE::instance()->getPosVel(name,
   t,
   trv.r,
   trv.v);
   _interval->insert(trv);
   t += samplingPeriod;
   }
   }
   
   void SpiceBodyInterpolatedPosVelCallback::getPosVel(const orsa::Time  & t,
   orsa::Vector      & relativePosition,
   orsa::Vector      & relativeVelocity) const {
   
   BodyGroup::TRV trv, _trv1, _trv2;
   // orsa::Interval<BodyGroup::TRV> * _b_interval = getBodyInterval(b);
   trv.t = t;
   if (_interval->getSubInterval(trv, _trv1, _trv2)) {
   if ((_trv1.t == _trv2.t) && (t == _trv1.t)) {
   // trv = _trv1;
   relativePosition = _trv1.r;
   relativeVelocity = _trv1.v;
   
   return;
   } else {
   
   // PhysicalSpline
   osg::ref_ptr<orsa::PhysicalSpline<orsa::Vector> > s = new orsa::PhysicalSpline<orsa::Vector>;
   // s->set(p0,v0,t0,p1,v1,t1);
   s->set(_trv1.r,
   _trv1.v,
   _trv1.t,
   _trv2.r,
   _trv2.v,
   _trv2.t);
   //
   s->get(relativePosition,
   relativeVelocity,
   t);
   
   return;
   }
   } else {
   
   ORSA_DEBUG("problems...");
   
   return;
   }	
   }
*/
