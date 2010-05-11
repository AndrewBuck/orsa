#include <orsa/interaction.h>

// #include <orsa/attitude.h>
#include <orsa/bodygroup.h>
// #include <orsa/legendre.h>
#include <orsa/paul.h>
#include <orsa/print.h>
#include <orsa/util.h>

using namespace orsa;

Interaction::Interaction() : osg::Referenced(true) { 
    dummyPaulMoment = new PaulMoment(0);
    dummyPaulMoment->setM(1,0,0,0);
    // dummyPaulMoment->setCenterOfMass(orsa::Vector(0,0,0));
    // dummyPaulMoment->setInertiaMoment(orsa::Matrix::identity());
}	

bool Interaction::acceleration(InteractionVector & a,  
                               orsa::BodyGroup   * bg,
                               const orsa::Time  & t) const {
  
    // ORSA_DEBUG("called...");
  
    BodyGroup::BodyList bl = bg->getBodyList();
  
    {
        // reset accel. vector
        a.resize(bl.size());
        for (unsigned int k=0; k<bl.size(); ++k) {
            a[k] = orsa::Vector(0,0,0);
        }
        /* 
           BodyGroup::BodyList::const_iterator ref_b_it = bl.begin();
           while (ref_b_it != bl.end()) {
           a[(*ref_b_it).get()] = orsa::Vector(0,0,0);
           ++ref_b_it;
           }
        */
    }
    //
    // a.resize(bg->largestBodyID()+1);
  
    {
        // main loop
        double m_ref_b, m_b;
    
        // BodyGroup::BodyList::const_iterator ref_b_it = bl.begin();
        // while (ref_b_it != bl.end()) {
    
        for (unsigned int j=0; j<bl.size(); ++j) {
      
            const orsa::Body * ref_b = bl[j].get();
      
            if (!ref_b->alive(t)) {
                continue;
            }
      
            // ORSA_DEBUG("ref_b is [%s]",ref_b->getName().c_str());
      
            if (ref_b->getInitialConditions().translational.get()) {
                if (!(ref_b->getInitialConditions().translational->dynamic())) {
                    continue;
                }
            } else {
                continue;
            }
      
            orsa::Vector thrust(0,0,0);
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
               // if (pe.thrustMagnitude.getRef() > 0) {
               // ORSA_DEBUG("using propulsion!! time: %f",t.get_d());
     
               // IMPORTANT!!!!
               // ORSA_DEBUG("a=F/m and F=thrust, so I have to divide by the instantaneous mass...");
	 
               // a[(*ref_b_it).get()] += pe.thrustMagnitude.getRef() * pe.thrustDirection();
	 
               const double bodyMass = ref_b->getMass() - ref_b->getPropulsion()->massLost(t);
	 
               if (ref_b->getPropulsion()->thrustToMass.isSet()) {
               ORSA_DEBUG("t: %f   bodyMass: %Fg",
               t.get_d(),
               bodyMass());
               }
	 
               if (bodyMass > 0) {
               a_ref_b += pe->thrust.getRef() / bodyMass;
               } else {
               ORSA_DEBUG("problem: non-positive mass!! (m=%Fe)",bodyMass());
               }
	 
               // } else {
               // ORSA_DEBUG("not using propulsion, time: %f",t.get_d());
               // }
               }
            */
      
            IBPS ref_b_ibps;
      
            if (!(bg->getInterpolatedIBPS(ref_b_ibps,
                                          ref_b,
                                          t))) {
                ORSA_DEBUG("problems...");
                return false;
            }
      
            if (!bg->getInterpolatedMass(m_ref_b,ref_b,t)) {
                ORSA_DEBUG("problems...");
            }
      
            for (unsigned int k=0; k<j; ++k) {
	
                const orsa::Body * b = bl[k].get();
	
                if (!b->alive(t)) {
                    continue;
                }
	
                // ORSA_DEBUG("b is [%s]",b->getName().c_str());
	
                if (!bg->getInterpolatedMass(m_b,b,t)) {
                    ORSA_DEBUG("problems...");
                }
	
                if (ref_b->nonInteractingGroup && 
                    b->nonInteractingGroup) {
                    continue;
                }
	
                // cannot do this!
                /* 
                   if ((*b_it)->getInitialConditions().translational.get()) {
                   if (!((*b_it)->getInitialConditions().translational->dynamic())) {
                   ++b_it;
                   continue;
                   }
                   }
                */
	
                if ((m_b == 0) && 
                    (m_ref_b == 0)) {
                    continue;
                }
	
                IBPS b_ibps;
	
                if (!(bg->getInterpolatedIBPS(b_ibps,
                                              b,
                                              t))) {
                    ORSA_DEBUG("problems...");
                    return false;
                }
	
                /* if (ref_b->getPaulMoment() || 
                   b->getPaulMoment()) {
                */
                //
                if (ref_b_ibps.inertial->paulMoment() || b_ibps.inertial->paulMoment()) {
	  
                    // ORSA_DEBUG("--MARK--");
	  
                    osg::ref_ptr<const PaulMoment> ref_b_pm = 
                        (ref_b_ibps.inertial->paulMoment()) ? 
                        (ref_b_ibps.inertial->paulMoment()) :
                        (dummyPaulMoment.get());
	  
                    osg::ref_ptr<const PaulMoment> b_pm = 
                        (b_ibps.inertial->paulMoment()) ? 
                        (b_ibps.inertial->paulMoment()) :
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
	  
                    /* const orsa::Vector R =
                       (b_ibps.translational->position() + b_l2g*b_pm->getCenterOfMass()) - 
                       (ref_b_ibps.translational->position() + ref_b_l2g*ref_b_pm->getCenterOfMass());
                    */
                    //
                    const orsa::Vector R =
                        b_ibps.translational->position() - 
                        ref_b_ibps.translational->position();
	  
                    /* 
                       ORSA_DEBUG("---R---");
                       orsa::print(R);
                       orsa::print(b_ibps.translational->position());
                       orsa::print(b_l2g*b_pm->getCenterOfMass());
                       orsa::print(ref_b_ibps.translational->position());
                       orsa::print(ref_b_l2g*ref_b_pm->getCenterOfMass());
                       orsa::print(t);
                    */
	  
                    const double     b_radius =     b_ibps.inertial->localShape() ?     b_ibps.inertial->localShape()->boundingRadius() : 0;
                    const double ref_b_radius = ref_b_ibps.inertial->localShape() ? ref_b_ibps.inertial->localShape()->boundingRadius() : 0;
	  
                    // if (b->getRadius()+ref_b->getRadius() > R.length()) {
                    //
                    if ((b_radius+ref_b_radius) > R.length()) {
	    
                        ORSA_DEBUG("bodies too close: R<R1+R2, R=%g R1=%g R2=%g [km] [b1:%s] [b2:%s]",
                                   orsa::FromUnits(R.length(),orsa::Unit::KM,-1),
                                   orsa::FromUnits(b_radius,orsa::Unit::KM,-1),
                                   orsa::FromUnits(ref_b_radius,orsa::Unit::KM,-1),
                                   b->getName().c_str(),
                                   ref_b->getName().c_str());
                        ORSA_DEBUG("reverting to pointlike...");
#warning should handle this better....
	   
                        orsa::Vector _d =
                            b_ibps.translational->position() - 
                            ref_b_ibps.translational->position();
	    
                        const double _l = _d.length();
	    
                        /* 
                           if (_l > epsilon()) {
                        */
	    
                        _d /= (_l*_l*_l);
	    
                        if (ref_b->betaSun == b) {
                            const orsa::Vector accTerm = (1 - ref_b->beta.getRef()) * _d;
                            // a[(*ref_b_it).get()] += (*b_it)->getMu()     * accTerm;
                            // a[    (*b_it).get()] -= (*ref_b_it)->getMu() * accTerm;
                            // a_ref_b += b->getMu()     * accTerm;
                            // a_b     -= ref_b->getMu() * accTerm;
                            a[j] += orsa::Unit::G() * m_b     * accTerm;
                            a[k] -= orsa::Unit::G() * m_ref_b * accTerm;
                        } else {
	     
                            // ORSA_DEBUG("--MARK--");
	      
                            const orsa::Vector accTerm = _d;
                            // a[(*ref_b_it).get()] += (*b_it)->getMu()     * accTerm;
                            // a[    (*b_it).get()] -= (*ref_b_it)->getMu() * accTerm;
                            // a_ref_b += b->getMu()     * accTerm;
                            // a_b     -= ref_b->getMu() * accTerm;
                            a[j] += orsa::Unit::G() * m_b     * accTerm;
                            a[k] -= orsa::Unit::G() * m_ref_b * accTerm;
                        }	
	    
                        /* ORSA_DEBUG("--tmp--");
                           orsa::print(a[j]);
                        */
	    
                    } else {
	    
                        // ORSA_DEBUG("--MARK--");
	    
                        const orsa::Vector accTerm =
                            Paul::gravitationalForce(ref_b_pm.get(),
                                                     ref_b_g2l,
                                                     b_pm.get(),
                                                     b_g2l,
                                                     R);
	    
                        a[j] += orsa::Unit::G() * m_b     * accTerm;
                        a[k] -= orsa::Unit::G() * m_ref_b * accTerm;
	    
                        /* 
                           {
                           ORSA_DEBUG("--tmp--");
                           orsa::print(Paul::gravitationalPotential(ref_b_pm.get(),
                           ref_b_g2l,
                           b_pm.get(),
                           b_g2l,
                           R));
                           orsa::print(Paul::gravitationalPotential(b_pm.get(),
                           b_g2l,
                           ref_b_pm.get(),
                           ref_b_g2l,
                           -R));
                           orsa::print(Paul::gravitationalForce(ref_b_pm.get(),
                           ref_b_g2l,
                           b_pm.get(),
                           b_g2l,
                           R));
                           orsa::print(Paul::gravitationalForce(b_pm.get(),
                           b_g2l,
                           ref_b_pm.get(),
                           ref_b_g2l,
                           -R));
                           orsa::print(Paul::gravitationalTorque(ref_b_pm.get(),
                           ref_b_g2l,
                           b_pm.get(),
                           b_g2l,
                           R));
                           orsa::print(Paul::gravitationalTorque(b_pm.get(),
                           b_g2l,
                           ref_b_pm.get(),
                           ref_b_g2l,
                           -R));
                           }
                        */
	    
                    }
	  
                } else {
	  
                    // ORSA_DEBUG("--MARK--");
	  
                    /* 
                       if (bp->b_ibps->translational.get() &&
                       bp->ref_b_ibps->translational.get()) {
                    */
	  
                    orsa::Vector _d =
                        b_ibps.translational->position() - 
                        ref_b_ibps.translational->position();
	  
                    const double _l = _d.length();
	  
                    if (_l > epsilon()) {
	    
                        _d /= (_l*_l*_l);
	    
                        if (ref_b->betaSun == b) {
                            const orsa::Vector accTerm = (1 - ref_b->beta.getRef()) * _d;
                            // a[(*ref_b_it).get()] += (*b_it)->getMu()     * accTerm;
                            // a[    (*b_it).get()] -= (*ref_b_it)->getMu() * accTerm;
                            // a_ref_b += b->getMu()     * accTerm;
                            // a_b     -= ref_b->getMu() * accTerm;
                            a[j] += orsa::Unit::G() * m_b     * accTerm;
                            a[k] -= orsa::Unit::G() * m_ref_b * accTerm;
                        } else {
                            const orsa::Vector accTerm = _d;
                            // a[(*ref_b_it).get()] += (*b_it)->getMu()     * accTerm;
                            // a[    (*b_it).get()] -= (*ref_b_it)->getMu() * accTerm;
                            // a_ref_b += b->getMu()     * accTerm;
                            // a_b     -= ref_b->getMu() * accTerm;
                            a[j] += orsa::Unit::G() * m_b     * accTerm;
                            a[k] -= orsa::Unit::G() * m_ref_b * accTerm;
                        }
	    
                        /* ORSA_DEBUG("--tmp--");
                           orsa::print(a[j]);
                        */
	    
                    } else {
	    
                        ORSA_WARNING("skipping: zero distance between bodies [%s] and [%s]",
                                     ref_b->getName().c_str(),
                                     b->getName().c_str());
                        ORSA_DEBUG("[%s] position:",b->getName().c_str());
                        orsa::print(b_ibps.translational->position());
                        ORSA_DEBUG("[%s] position:",ref_b->getName().c_str());
                        orsa::print(ref_b_ibps.translational->position());
                    }
	  
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
	
            }
      
            // thrust
            if (ref_b->propulsion.get()) {
                if (m_ref_b < orsa::epsilon()) {
                    ORSA_DEBUG("Propulsion problems: non-positive mass...");
                } else {
                    /* ORSA_DEBUG("adding propulsion thrust for body [%s]",bl[j]->getName().c_str());
                       print(j);
                       print(a[j]);
                       print(thrust);
                       print(m_ref_b);
                    */
	  
                    a[j] += thrust / m_ref_b;
	  
                    // print(a[j]);
                }
	
            }
      
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
        double m_ref_b, m_b;
    
        // BodyGroup::BodyList::const_iterator ref_b_it = bl.begin();
        // while (ref_b_it != bl.end()) {
    
        for (unsigned int j=0; j<bl.size(); ++j) {
      
            const orsa::Body * ref_b = bl[j].get();
      
            if (!ref_b->alive(t)) {
                // ++ref_b_it;
                continue;
            }
      
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
      
            // propulsion torque?
      
            IBPS ref_b_ibps;
      
            if (!(bg->getInterpolatedIBPS(ref_b_ibps,
                                          ref_b,
                                          t))) {
                ORSA_DEBUG("problems...");
                return false;
            }
      
            osg::ref_ptr<const PaulMoment> ref_b_pm = 
                (ref_b_ibps.inertial->paulMoment()) ? 
                (ref_b_ibps.inertial->paulMoment()) :
                (dummyPaulMoment.get());
      
            if (!bg->getInterpolatedMass(m_ref_b,ref_b,t)) {
                ORSA_DEBUG("problems...");
            }
      
            // BodyGroup::BodyList::const_iterator b_it = bl.begin();
            // while (b_it != bl.end()) {
      
            for (unsigned int k=0; k<j; ++k) {
	
                const orsa::Body * b = bl[k].get();
	
                if (!b->alive(t)) {
                    // ++b_it;
                    continue;
                }
	
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
	
                if (!bg->getInterpolatedMass(m_b,b,t)) {
                    ORSA_DEBUG("problems...");
                }
	
                if ((m_b == 0) && 
                    (m_ref_b == 0)) {
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
	
                /* if (ref_b->getPaulMoment() || 
                   b->getPaulMoment()) {
                */
                //
                if (ref_b_ibps.inertial->paulMoment() || b_ibps.inertial->paulMoment()) {
	  
                    osg::ref_ptr<const PaulMoment> b_pm = 
                        (b_ibps.inertial->paulMoment()) ? 
                        (b_ibps.inertial->paulMoment()) :
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
	    
                        /* const orsa::Vector R =
                           (b_ibps.translational->position() + b_l2g*b_pm->getCenterOfMass()) - 
                           (ref_b_ibps.translational->position() + ref_b_l2g*ref_b_pm->getCenterOfMass());
                        */
                        //
                        const orsa::Vector R =
                            b_ibps.translational->position() - 
                            ref_b_ibps.translational->position();
	    
                        const orsa::Vector torqueTerm =
                            Paul::gravitationalTorque(ref_b_pm.get(),
                                                      ref_b_g2l,
                                                      b_pm.get(),
                                                      b_g2l,
                                                      R);
	    
                        // N[(*ref_b_it)->id()] += orsa::Unit::instance()->getG() * m_b * torqueTerm;
                        // N[    (*b_it)->id()] -= orsa::Unit::instance()->getG() * m_ref_b * torqueTerm;
                        //
                        N[j] += orsa::Unit::G() * m_b     * torqueTerm;
                        N[k] -= orsa::Unit::G() * m_ref_b * torqueTerm;
	    
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
