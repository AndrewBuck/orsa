#include <orsa/integrator_leapfrog.h>

// #include <orsa/attitude.h>
#include <orsa/body.h>
#include <orsa/bodygroup.h>
#include <orsa/euler.h>
#include <orsa/interaction.h>
#include <orsa/print.h>
#include <orsa/util.h>

using namespace orsa;

IntegratorLeapFrog::IntegratorLeapFrog() : Integrator() {
  
}

IntegratorLeapFrog::~IntegratorLeapFrog() {
  
}

bool IntegratorLeapFrog::step(orsa::BodyGroup  * bg,
                              const orsa::Time & start,
                              const orsa::Time & timestep,
                              orsa::Time       & next_timestep) {
  
    /* 
       ORSA_DEBUG("start:     %12.6f [day]",
       FromUnits(FromUnits(start.getMuSec(),Unit::MICROSECOND),Unit::DAY,-1)());
       ORSA_DEBUG("timestep:  %12.6f [day]",
       FromUnits(FromUnits(timestep.getMuSec(),Unit::MICROSECOND),Unit::DAY,-1)());
    */
  
    const Time _h  = timestep;
    const Time _h2 = _h/2;
    //
    orsa::IBPS ibps;
  
    // drift h/2
    {
        BodyGroup::BodyList::const_iterator _b_it = bg->getBodyList().begin();
        while (_b_it != bg->getBodyList().end()) {
      
            if (!(*_b_it)->alive(start)) {
                ++_b_it;
                continue;
            }
      
            if (bg->getInterpolatedIBPS(ibps,(*_b_it).get(),start)) {
	
	
                if (!ibps.translational->dynamic()) {
                    ++_b_it;
                    continue;
                }
	
                ibps.time += _h2;
	
                if ((*_b_it)->getInitialConditions().translational.get()) {
                    if ((*_b_it)->getInitialConditions().translational->dynamic()) {
                        ibps.translational->setPosition(ibps.translational->position() +
                                                        ibps.translational->velocity() * _h2.get_d());
                    }
                }
	
                if ((*_b_it)->getInitialConditions().rotational.get()) {
                    if ((*_b_it)->getInitialConditions().rotational->dynamic()) {
	    
                        // old code
                        /* 
                           ibps.rotational->set(RotationalBodyProperty::qFiniteRotation(ibps.rotational->getQ(),
                           ibps.rotational->getOmega(),
                           _h2),
                           ibps.rotational->getOmega());
                        */
	    
                        // new code: nothing done now...
                        ibps.rotational->set(ibps.rotational->getQ(),
                                             ibps.rotational->getOmega());
	    
                    }
                }
	
                ibps.tmp = true;
	
                if (!(bg->getBodyInterval((*_b_it).get())->insert(ibps,onlyIfExtending.getRef(),false))) {
                    ORSA_DEBUG("problems with insert, body [%s]",
                               (*_b_it)->getName().c_str());
                }
	
            } else {
                ORSA_ERROR("point not present in interval, body [%s]",
                           (*_b_it)->getName().c_str());
                return false;
            }
      
            ++_b_it;
        }
    }
  
    // kick h
    // orsa::Interaction::VectorHash a;
    orsa::Interaction::InteractionVector a;
    //
    if (!bg->getInteraction()->acceleration(a,  
                                            bg,
                                            start+_h2)) {
        ORSA_DEBUG("problems...");
        return false;
    }
  
    // orsa::Interaction::VectorHash N;
    orsa::Interaction::InteractionVector N;
    //
    if (!bg->getInteraction()->torque(N,  
                                      bg,
                                      start+_h2)) {
        ORSA_DEBUG("problems...");
        return false;
    }

    BodyGroup::BodyList bl = bg->getBodyList();
  
    // _b_it = bg->getBodyList().begin();
    //while (_b_it != bg->getBodyList().end()) {
    //
    for (unsigned int j=0; j<bg->getBodyList().size(); ++j) {
    
        const orsa::Body * b = bl[j].get();
    
        if (!b->alive(start+_h2)) {
            continue;
        }
    
        if (bg->getInterpolatedIBPS(ibps,b,start+_h2)) {
      
            if (!ibps.translational->dynamic()) {
                continue;
            }
      
            if (b->getInitialConditions().translational.get()) {
                if (b->getInitialConditions().translational->dynamic()) {
                    ibps.translational->setVelocity(ibps.translational->velocity() +
                                                    a[j] * _h.get_d());
                }
            }
      
            if (b->getInitialConditions().rotational.get()) {
                if (b->getInitialConditions().rotational->dynamic()) {
	  
#warning "inertia moments needed!! plus rotation to get to the principal axis..."
	  
                    // if (b->getPaulMoment()) {
                    if (ibps.inertial->paulMoment()) {
	    
                        double mass;
                        if (!bg->getInterpolatedMass(mass,b,start+_h2)) {
                            ORSA_DEBUG("problems...");
                        }	
	    
                        orsa::Matrix inertiaMoment = 
                            mass * 
                            ibps.inertial->inertiaMatrix();
                        // b->getPaulMoment()->getInertiaMoment();
	    
                        // new approach, following I.P. Omelyan (1999)
                        if (1) {
	      
                            // const orsa::Matrix g2l = BodyAttitude(b,bg).globalToLocal(start+_h2);
                            // const orsa::Matrix l2g = BodyAttitude(b,bg).localToGlobal(start+_h2);
	      
                            // osg::ref_ptr<orsa::Attitude> attitude = new orsa::BodyAttitude(b,bg);
	      
                            // const orsa::Matrix g2l = attitude->globalToLocal(start+_h2);
                            // const orsa::Matrix l2g = attitude->localToGlobal(start+_h2);
	      
                            const orsa::Matrix g2l = orsa::globalToLocal(b,bg,start+_h2);
                            const orsa::Matrix l2g = orsa::localToGlobal(b,bg,start+_h2);
	      
                            const orsa::Vector oldOmega = ibps.rotational->getOmega();
                            const orsa::Matrix I        = inertiaMoment;
                            const orsa::Vector T        = N[j];
	      
                            orsa::Vector omegaDot;
                            orsa::Vector omegaIter = oldOmega;
                            orsa::Vector oldOmegaIter;
	      
                            unsigned int niter = 0;
                            //
                            do {
		
                                ++niter;
                                // ORSA_DEBUG("iteration: %i",niter);
		
                                oldOmegaIter = omegaIter;
		
                                // we need to implement Euler manually here, to use Eq.(5)
                                {
                                    orsa::Matrix Ip, genericToPrincipal, principalToGeneric;
		  
                                    orsa::principalAxis(genericToPrincipal,
                                                        Ip,
                                                        I);
		  
                                    if (!orsa::Matrix::invert(genericToPrincipal, principalToGeneric)) {
                                        ORSA_DEBUG("problems...");
                                    }
		  
                                    /* 
                                       print(I);
                                       print(Ip);
                                       print(genericToPrincipal*Ip*principalToGeneric);
                                       print(principalToGeneric*I *genericToPrincipal);
		     
                                       print(genericToPrincipal);
                                       print(principalToGeneric);
                                       print(principalToGeneric*genericToPrincipal);
                                       print(genericToPrincipal*principalToGeneric);
                                    */
		  
                                    const orsa::Matrix globalToPrincipal = genericToPrincipal * g2l;
                                    const orsa::Matrix principalToGlobal = l2g * principalToGeneric;
		  
                                    const orsa::Matrix & g2p = globalToPrincipal;
                                    const orsa::Matrix & p2g = principalToGlobal;
		  
                                    const orsa::Vector T_p  = g2p * T;
		  
                                    const orsa::Vector oldOmega_p  = g2p * oldOmega;
                                    const orsa::Vector omegaIter_p = g2p * omegaIter; 
		  
                                    // Eq.(5)
                                    const double omegaYZn_p = (oldOmega_p.getY()*oldOmega_p.getZ() + omegaIter_p.getY()*omegaIter_p.getZ()) / 2;
                                    const double omegaZXn_p = (oldOmega_p.getZ()*oldOmega_p.getX() + omegaIter_p.getZ()*omegaIter_p.getX()) / 2;
                                    const double omegaXYn_p = (oldOmega_p.getX()*oldOmega_p.getY() + omegaIter_p.getX()*omegaIter_p.getY()) / 2;
		  
                                    const double Jx = Ip.getM11();
                                    const double Jy = Ip.getM22();
                                    const double Jz = Ip.getM33();
		  
                                    // Eq.(4)
                                    omegaIter.setX(oldOmega_p.getX() - (_h.get_d()/Jx) * (T_p.getX() + (Jy-Jz)*omegaYZn_p));
                                    omegaIter.setY(oldOmega_p.getY() - (_h.get_d()/Jy) * (T_p.getY() + (Jz-Jx)*omegaZXn_p));
                                    omegaIter.setZ(oldOmega_p.getZ() - (_h.get_d()/Jz) * (T_p.getZ() + (Jx-Jy)*omegaXYn_p));
                                    //
                                    omegaIter = p2g * omegaIter;
                                }
		
                                /* 
                                   ORSA_DEBUG("--prints--");
                                   print(oldOmega);
                                   print(omegaIter);
                                   print(omegaIter-oldOmegaIter);
                                */
		
                                // ORSA_DEBUG("use better/normalized check...");
                            } while ((omegaIter-oldOmegaIter).length() > orsa::epsilon());
	      
                            ibps.rotational->set(ibps.rotational->getQ(),
                                                 omegaIter);
	      
                        }
	    
                    } else {
                        ORSA_DEBUG("PaulMoment not available for body [%s]",
                                   b->getName().c_str());
                    }
                }
            }
      
            ibps.tmp = true;
      
            if (!(bg->getBodyInterval(b)->insert(ibps,onlyIfExtending.getRef(),true))) {
                ORSA_DEBUG("problems with insert, body [%s]",
                           b->getName().c_str());
            }
      
        } else {
            ORSA_ERROR("point not present in interval, body [%s]",
                       b->getName().c_str());
            return false;
        }
    
        // ++_b_it;
    }
  
    // drift h/2
    {
        BodyGroup::BodyList::const_iterator _b_it = bg->getBodyList().begin();
        while (_b_it != bg->getBodyList().end()) {
      
            if (!(*_b_it)->alive(start+_h2)) {
                ++_b_it;
                continue;
            }
      
            if (bg->getInterpolatedIBPS(ibps,(*_b_it).get(),start+_h2)) {
	
                if (!ibps.translational->dynamic()) {
                    ++_b_it;
                    continue;
                }
	
                ibps.time += _h2;
	
                if ((*_b_it)->getInitialConditions().translational.get()) {
                    if ((*_b_it)->getInitialConditions().translational->dynamic()) {
                        ibps.translational->setPosition(ibps.translational->position() + 
                                                        ibps.translational->velocity() * _h2.get_d());
                    }
                }
	
                if ((*_b_it)->getInitialConditions().rotational.get()) {
                    if ((*_b_it)->getInitialConditions().rotational->dynamic()) {
	    
                        // old code
                        /* 
                           ibps.rotational->set(RotationalBodyProperty::qFiniteRotation(ibps.rotational->getQ(),
                           ibps.rotational->getOmega(),
                           _h2),
                           ibps.rotational->getOmega());
                        */
	    
                        // new code
                        {
                            const double omegaSq = ibps.rotational->getOmega().lengthSquared();
                            const double     hSq = _h.get_d()*_h.get_d();
	      
                            const orsa::Quaternion qDot = RotationalBodyProperty::qDot(ibps.rotational->getQ(),
                                                                                       ibps.rotational->getOmega());
	      
                            const orsa::Quaternion q = ( (1 - hSq*omegaSq/16) * ibps.rotational->getQ() + 
                                                         _h.get_d() * qDot ) / 
                                (1 + hSq*omegaSq/16);
	      
                            ibps.rotational->set(unitQuaternion(q),
                                                 ibps.rotational->getOmega());    
                        }
	    
                    }
                }
	
                ibps.tmp = false;
	
                if (!(bg->getBodyInterval((*_b_it).get())->insert(ibps,onlyIfExtending.getRef(),false))) {
                    ORSA_DEBUG("problems with insert, body [%s]",
                               (*_b_it)->getName().c_str());
                }
	
            } else {
                ORSA_ERROR("point not present in interval, body [%s]",
                           (*_b_it)->getName().c_str());
                return false;
            }
      
            ++_b_it;
        }
    }
  
    // constant timestep integrator
    next_timestep = timestep;
  
    return true;
}

