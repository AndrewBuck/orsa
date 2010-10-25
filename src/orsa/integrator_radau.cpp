#include <orsa/integrator_radau.h>

// #include <orsa/attitude.h>
#include <orsa/bodygroup.h>
#include <orsa/crash.h>
#include <orsa/euler.h>
#include <orsa/interaction.h>
#include <orsa/print.h>
#include <orsa/util.h>

#include <algorithm>

using namespace orsa;

// #warning "remember: constraint on qDot too!!"

IntegratorRadau::IntegratorRadau() : Integrator() {
    _init();
}

IntegratorRadau::~IntegratorRadau() {
  
}

bool IntegratorRadau::step(orsa::BodyGroup  * bg,
                           const orsa::Time & start,
                           const orsa::Time & timestep,
                           orsa::Time       & next_timestep) {
  
    // ORSA_DEBUG("called...");
  
    unsigned int niter = 2;
  
    const BodyGroup::BodyList & bl = bg->getBodyList();
  
    if (bg->size() != size) {
        _body_mass_or_number_changed(bg,start);
        niter = 6;
    } else {
        double m;
        BodyGroup::BodyList::const_iterator bl_it = bl.begin();
        while (bl_it != bl.end()) {
            if (!(*bl_it)->alive(start)) {
                ++bl_it;
                continue;
            }
            if (!bg->getInterpolatedMass(m,(*bl_it).get(),start)) {
                ORSA_DEBUG("problems...");
            }
            if (m != mass[(*bl_it).get()]) {
                _body_mass_or_number_changed(bg,start);
                niter = 6;
                break;
            }
            ++bl_it;
        }
    }
  
    if (_lastCallRejected.isSet()) {
        if (_lastCallRejected.get()) {
            niter = 6;
        }
    }
  
    // cerr << "niter: " << niter << endl;
    
    // interaction->Acceleration(frame_out,acc);
    /* {
       BodyGroup::BodyList::const_iterator bl_it = bl.begin();
       while (bl_it != bl.end()) {
       // if ((*bl_it)->getBodyPosVelCallback() != 0) {
       if (!((*bl_it)->getInitialConditions().translational->dynamic())) {
       // ORSA_DEBUG("skipping body [%s]",(*bl_it)->getName().c_str());
       ++bl_it;
       continue;
       }
       //
       bg->getInteraction()->acceleration(acc[(*bl_it).get()],  
       (*bl_it).get(),
       bg,
       start);
       ++bl_it;
       }
       }
    */
    //
    if (!bg->getInteraction()->acceleration(acc,  
                                            bg,
                                            start)) {
        ORSA_DEBUG("problems...");
        return false;
    }
    //
    if (!bg->getInteraction()->torque(torque,  
                                      bg,
                                      start)) {
        ORSA_DEBUG("problems...");
        return false;
    }
  
    // ORSA_DEBUG("--MARK--");
  
    // unsigned int j,k;
    /* for(k=0;k<frame_in.size();++k) {
       if (interaction->IsSkippingJPLPlanets() && frame_in[k].JPLPlanet() != NONE) continue;
     
       x1[3*k]   = frame_in[k].position().x;
       x1[3*k+1] = frame_in[k].position().y;
       x1[3*k+2] = frame_in[k].position().z;
       //
       v1[3*k]   = frame_in[k].velocity().x;
       v1[3*k+1] = frame_in[k].velocity().y;
       v1[3*k+2] = frame_in[k].velocity().z;
       //
       a1[3*k]   = acc[k].x;
       a1[3*k+1] = acc[k].y;  
       a1[3*k+2] = acc[k].z;
       }
    */
    //
    {
        // BodyGroup::BodyList::const_iterator bl_it = bl.begin();
        // while (bl_it != bl.end()) {
        for (unsigned int bodyIndex=0; bodyIndex<bl.size(); ++bodyIndex) {
      
            const orsa::Body * b = bl[bodyIndex].get();
      
            if (!b->alive(start)) {
                continue;
            }
      
            //
            // if (bg->getTRV(trv,(*bl_it).get(),start)) {
            // IBPS ibps;
            // if (bg->getInterpolatedIBPS(ibps,(*bl_it).get(),start)) {
      
            IBPS ibps;
            if (!bg->getInterpolatedIBPS(ibps,b,start)) {
                ORSA_DEBUG("problem, body: [%s]",b->getName().c_str());
            }
      
            const orsa::Body * k = b;
      
            if (b->getInitialConditions().translational.get()) {
                if (b->getInitialConditions().translational->dynamic()) {
                    // if (bg->getInterpolatedIBPS(ibps,b,start)) {
                    x1[k] = ibps.translational->position();
                    v1[k] = ibps.translational->velocity();
                    a1[k] = acc[bodyIndex];
                    /* } else {
                       ORSA_DEBUG("problem, body: [%s]",b->getName().c_str());
                       } */
                }
            }
            //
            /* 
               x1[k] = ibps.translational->position();
               v1[k] = ibps.translational->velocity();
               a1[k] = acc[(*bl_it).get()].getRef();
            */
      
            // should this be here??
            // if (b->getInitialConditions().translational.get()) {
            // if (b->getInitialConditions().translational->dynamic()) {
      
            if (b->getInitialConditions().rotational.get()) {
                if (b->getInitialConditions().rotational->dynamic()) {
	  
                    // #warning "inertia moments needed!! plus rotation to get to the principal axis..."
	  
                    // if (b->getPaulMoment()) {
                    if (ibps.inertial->paulMoment()) {
	    
                        double m;
                        if (!bg->getInterpolatedMass(m,b,start)) {
                            ORSA_DEBUG("problems...");
                        }	
	    
                        orsa::Matrix inertiaMoment = 
                            m * 
                            ibps.inertial->inertiaMatrix();
                        // b->getPaulMoment()->getInertiaMoment();
	    
                        //
                        {
	      
                            // osg::ref_ptr<orsa::Attitude> attitude = new orsa::BodyAttitude((*bl_it).get(),bg);
	      
                            // const orsa::Matrix g2l = attitude->globalToLocal(start);
                            // const orsa::Matrix l2g = attitude->localToGlobal(start);
	      
                            const orsa::Matrix g2l = orsa::globalToLocal(b,bg,start);
                            const orsa::Matrix l2g = orsa::localToGlobal(b,bg,start);
	      
                            const orsa::Vector omega = ibps.rotational->getOmega();
                            const orsa::Matrix I     = inertiaMoment;
                            const orsa::Vector T     = torque[bodyIndex];
	      
                            orsa::Vector omegaDot;
	      
                            Euler(omegaDot,
                                  g2l*omega,
                                  I,
                                  g2l*T);
                            //
                            omegaDot = l2g * omegaDot;
	      
                            Q1[k]       = ibps.rotational->getQ();
                            Q1Dot[k]    = RotationalBodyProperty::qDot(Q1[k], omega);
                            Q1DotDot[k] = RotationalBodyProperty::qDotDot(Q1[k], omega, omegaDot);
                        }
	    
                    } else {
                        ORSA_DEBUG("PaulMoment not available for body [%s]",
                                   b->getName().c_str());
                    }
                }
            }
      
            // }
            // }
      
            // } else {
            // ORSA_DEBUG("problem, body: [%s]",(*bl_it).get()->getName().c_str());
            // }
      
            // ++bl_it;
        }
    }
  
    // ORSA_DEBUG("--MARK--");
  
    /* 
       for(k=0;k<nv;++k) {
       if (interaction->IsSkippingJPLPlanets() && frame_in[k/3].JPLPlanet() != NONE) continue;
     
       g[0][k] = b[6][k]*d[15] + b[5][k]*d[10] + b[4][k]*d[6] + b[3][k]*d[3]  + b[2][k]*d[1]  + b[1][k]*d[0]  + b[0][k];
       g[1][k] = b[6][k]*d[16] + b[5][k]*d[11] + b[4][k]*d[7] + b[3][k]*d[4]  + b[2][k]*d[2]  + b[1][k];
       g[2][k] = b[6][k]*d[17] + b[5][k]*d[12] + b[4][k]*d[8] + b[3][k]*d[5]  + b[2][k];
       g[3][k] = b[6][k]*d[18] + b[5][k]*d[13] + b[4][k]*d[9] + b[3][k];
       g[4][k] = b[6][k]*d[19] + b[5][k]*d[14] + b[4][k];
       g[5][k] = b[6][k]*d[20] + b[5][k];
       g[6][k] = b[6][k];
       }
    */
    //
    {
        BodyGroup::BodyList::const_iterator bl_it = bl.begin();
        while (bl_it != bl.end()) {
            // if ((*bl_it)->getBodyPosVelCallback() != 0) {
            /* 
               if (!((*bl_it)->getInitialConditions().translational->dynamic())) {
               // ORSA_DEBUG("skipping body [%s]",(*bl_it)->getName().c_str());
               ++bl_it;
               continue;
               }
            */

            const orsa::Body * k = (*bl_it).get();
      
            if ((*bl_it)->getInitialConditions().translational.get()) {
                if ((*bl_it)->getInitialConditions().translational->dynamic()) {
	  
                    g[0][k] = b[6][k]*d[15] + b[5][k]*d[10] + b[4][k]*d[6] + b[3][k]*d[3]  + b[2][k]*d[1]  + b[1][k]*d[0]  + b[0][k];
                    g[1][k] = b[6][k]*d[16] + b[5][k]*d[11] + b[4][k]*d[7] + b[3][k]*d[4]  + b[2][k]*d[2]  + b[1][k];
                    g[2][k] = b[6][k]*d[17] + b[5][k]*d[12] + b[4][k]*d[8] + b[3][k]*d[5]  + b[2][k];
                    g[3][k] = b[6][k]*d[18] + b[5][k]*d[13] + b[4][k]*d[9] + b[3][k];
                    g[4][k] = b[6][k]*d[19] + b[5][k]*d[14] + b[4][k];
                    g[5][k] = b[6][k]*d[20] + b[5][k];
                    g[6][k] = b[6][k];
	  
                }
            }
      
            /* 
               if ((*bl_it)->getInitialConditions().rotational.get()) {
               if ((*bl_it)->getInitialConditions().rotational->dynamic()) {
	 
               gPhi[0][k] = bPhi[6][k]*d[15] + bPhi[5][k]*d[10] + bPhi[4][k]*d[6] + bPhi[3][k]*d[3]  + bPhi[2][k]*d[1]  + bPhi[1][k]*d[0]  + bPhi[0][k];
               gPhi[1][k] = bPhi[6][k]*d[16] + bPhi[5][k]*d[11] + bPhi[4][k]*d[7] + bPhi[3][k]*d[4]  + bPhi[2][k]*d[2]  + bPhi[1][k];
               gPhi[2][k] = bPhi[6][k]*d[17] + bPhi[5][k]*d[12] + bPhi[4][k]*d[8] + bPhi[3][k]*d[5]  + bPhi[2][k];
               gPhi[3][k] = bPhi[6][k]*d[18] + bPhi[5][k]*d[13] + bPhi[4][k]*d[9] + bPhi[3][k];
               gPhi[4][k] = bPhi[6][k]*d[19] + bPhi[5][k]*d[14] + bPhi[4][k];
               gPhi[5][k] = bPhi[6][k]*d[20] + bPhi[5][k];
               gPhi[6][k] = bPhi[6][k];
	 
               gTheta[0][k] = bTheta[6][k]*d[15] + bTheta[5][k]*d[10] + bTheta[4][k]*d[6] + bTheta[3][k]*d[3]  + bTheta[2][k]*d[1]  + bTheta[1][k]*d[0]  + bTheta[0][k];
               gTheta[1][k] = bTheta[6][k]*d[16] + bTheta[5][k]*d[11] + bTheta[4][k]*d[7] + bTheta[3][k]*d[4]  + bTheta[2][k]*d[2]  + bTheta[1][k];
               gTheta[2][k] = bTheta[6][k]*d[17] + bTheta[5][k]*d[12] + bTheta[4][k]*d[8] + bTheta[3][k]*d[5]  + bTheta[2][k];
               gTheta[3][k] = bTheta[6][k]*d[18] + bTheta[5][k]*d[13] + bTheta[4][k]*d[9] + bTheta[3][k];
               gTheta[4][k] = bTheta[6][k]*d[19] + bTheta[5][k]*d[14] + bTheta[4][k];
               gTheta[5][k] = bTheta[6][k]*d[20] + bTheta[5][k];
               gTheta[6][k] = bTheta[6][k];
	 
               gPsi[0][k] = bPsi[6][k]*d[15] + bPsi[5][k]*d[10] + bPsi[4][k]*d[6] + bPsi[3][k]*d[3]  + bPsi[2][k]*d[1]  + bPsi[1][k]*d[0]  + bPsi[0][k];
               gPsi[1][k] = bPsi[6][k]*d[16] + bPsi[5][k]*d[11] + bPsi[4][k]*d[7] + bPsi[3][k]*d[4]  + bPsi[2][k]*d[2]  + bPsi[1][k];
               gPsi[2][k] = bPsi[6][k]*d[17] + bPsi[5][k]*d[12] + bPsi[4][k]*d[8] + bPsi[3][k]*d[5]  + bPsi[2][k];
               gPsi[3][k] = bPsi[6][k]*d[18] + bPsi[5][k]*d[13] + bPsi[4][k]*d[9] + bPsi[3][k];
               gPsi[4][k] = bPsi[6][k]*d[19] + bPsi[5][k]*d[14] + bPsi[4][k];
               gPsi[5][k] = bPsi[6][k]*d[20] + bPsi[5][k];
               gPsi[6][k] = bPsi[6][k];
	 
               }
               }
            */
            //
     
            // if ((*bl_it)->getInitialConditions().translational.get()) {
            // if ((*bl_it)->getInitialConditions().translational->dynamic()) {
      
            if ((*bl_it)->getInitialConditions().rotational.get()) {
                if ((*bl_it)->getInitialConditions().rotational->dynamic()) {
	  
                    gQ[0][k] = bQ[6][k]*d[15] + bQ[5][k]*d[10] + bQ[4][k]*d[6] + bQ[3][k]*d[3]  + bQ[2][k]*d[1]  + bQ[1][k]*d[0]  + bQ[0][k];
                    gQ[1][k] = bQ[6][k]*d[16] + bQ[5][k]*d[11] + bQ[4][k]*d[7] + bQ[3][k]*d[4]  + bQ[2][k]*d[2]  + bQ[1][k];
                    gQ[2][k] = bQ[6][k]*d[17] + bQ[5][k]*d[12] + bQ[4][k]*d[8] + bQ[3][k]*d[5]  + bQ[2][k];
                    gQ[3][k] = bQ[6][k]*d[18] + bQ[5][k]*d[13] + bQ[4][k]*d[9] + bQ[3][k];
                    gQ[4][k] = bQ[6][k]*d[19] + bQ[5][k]*d[14] + bQ[4][k];
                    gQ[5][k] = bQ[6][k]*d[20] + bQ[5][k];
                    gQ[6][k] = bQ[6][k];
                }
            }
      
            // }
            // }
      
            ++bl_it;
        }
    }
  
    // ORSA_DEBUG("--MARK--");
  
    // orsa::Vector tmp,gk;
    // double q1,q2,q3,q4,q5,q6,q7;
  
    // unsigned int main_loop_counter;
  
    for (unsigned int main_loop_counter=0; main_loop_counter<niter; ++main_loop_counter) {
    
        double s[9];
    
        for(unsigned int j=1; j<8; ++j) {
      
            // s[0] = timestep * h[j];
            s[0] = timestep.get_d() * h[j];
            s[1] = s[0] * s[0] * 0.5;
            s[2] = s[1] * h[j] * 0.3333333333333333;
            s[3] = s[2] * h[j] * 0.5;
            s[4] = s[3] * h[j] * 0.6;
            s[5] = s[4] * h[j] * 0.6666666666666667;
            s[6] = s[5] * h[j] * 0.7142857142857143;
            s[7] = s[6] * h[j] * 0.75;
            s[8] = s[7] * h[j] * 0.7777777777777778;
      
            /* 
               for(k=0;k<nv;++k) {
               if (interaction->IsSkippingJPLPlanets() && frame_in[k/3].JPLPlanet() != NONE) continue;
	
               x[k] = ( s[8]*b[6][k] +
               s[7]*b[5][k] + 
               s[6]*b[4][k] + 
               s[5]*b[3][k] + 
               s[4]*b[2][k] + 
               s[3]*b[1][k] + 
               s[2]*b[0][k] ) +
               s[1]*a1[k] + 
               s[0]*v1[k] + 
               x1[k];
               }
            */
            //
            {
                BodyGroup::BodyList::const_iterator bl_it = bl.begin();
                while (bl_it != bl.end()) {
	  
                    const orsa::Body * k = (*bl_it).get();
	
                    if ((*bl_it)->getInitialConditions().translational.get()) {
                        if ((*bl_it)->getInitialConditions().translational->dynamic()) {
	      
                            x[k] = ( s[8]*b[6][k] +
                                     s[7]*b[5][k] + 
                                     s[6]*b[4][k] + 
                                     s[5]*b[3][k] + 
                                     s[4]*b[2][k] + 
                                     s[3]*b[1][k] + 
                                     s[2]*b[0][k] ) +
                                s[1]*a1[k] + 
                                s[0]*v1[k] + 
                                x1[k];
	      
                            /* 
                               ORSA_DEBUG("--MARK-- x[k] , k=%s",(*bl_it)->getName().c_str());
                               orsa::print(x[k]);
		 
                               ORSA_DEBUG("--MARK-- x1[k] , k=%s",(*bl_it)->getName().c_str());
                               orsa::print(x1[k]);
		 
                               for (unsigned int dd=0; dd<=6; ++dd) {
                               ORSA_DEBUG("b[%i][k]:",dd);
                               orsa::print(b[dd][k]);
                               }	
		 
                               for (unsigned int sj=0; sj<=8; ++sj) {
                               ORSA_DEBUG("s[%i]:",sj);
                               orsa::print(s[sj]);	
                               }
		 
                               ORSA_DEBUG("a1:");
                               orsa::print(a1[k]);
		 
                               ORSA_DEBUG("v1:");
                               orsa::print(v1[k]);
                            */
	      
                        }
                    }
	  
                    /* 
                       x[k] = ( s[8]*b[6][k] +
                       s[7]*b[5][k] + 
                       s[6]*b[4][k] + 
                       s[5]*b[3][k] + 
                       s[4]*b[2][k] + 
                       s[3]*b[1][k] + 
                       s[2]*b[0][k] ) +
                       s[1]*a1[k] + 
                       s[0]*v1[k] + 
                       x1[k];
                    */

                    /* 
                       if ((*bl_it)->getInitialConditions().rotational.get()) {
                       if ((*bl_it)->getInitialConditions().rotational->dynamic()) {
	     
                       phi[k] = ( s[8]*bPhi[6][k] +
                       s[7]*bPhi[5][k] + 
                       s[6]*bPhi[4][k] + 
                       s[5]*bPhi[3][k] + 
                       s[4]*bPhi[2][k] + 
                       s[3]*bPhi[1][k] + 
                       s[2]*bPhi[0][k] ) +
                       s[1]*phi1DotDot[k] + 
                       s[0]*phi1Dot[k] + 
                       phi1[k];
	     
                       theta[k] = ( s[8]*bTheta[6][k] +
                       s[7]*bTheta[5][k] + 
                       s[6]*bTheta[4][k] + 
                       s[5]*bTheta[3][k] + 
                       s[4]*bTheta[2][k] + 
                       s[3]*bTheta[1][k] + 
                       s[2]*bTheta[0][k] ) +
                       s[1]*theta1DotDot[k] + 
                       s[0]*theta1Dot[k] + 
                       theta1[k];
	     
                       psi[k] = ( s[8]*bPsi[6][k] +
                       s[7]*bPsi[5][k] + 
                       s[6]*bPsi[4][k] + 
                       s[5]*bPsi[3][k] + 
                       s[4]*bPsi[2][k] + 
                       s[3]*bPsi[1][k] + 
                       s[2]*bPsi[0][k] ) +
                       s[1]*psi1DotDot[k] + 
                       s[0]*psi1Dot[k] + 
                       psi1[k];
                       }
                       }
                    */
                    //
                    // if ((*bl_it)->getInitialConditions().translational.get()) {
                    // if ((*bl_it)->getInitialConditions().translational->dynamic()) {
	  
                    if ((*bl_it)->getInitialConditions().rotational.get()) {
                        if ((*bl_it)->getInitialConditions().rotational->dynamic()) {
	      
                            // __THERE__
                            // moved part of this to __HERE__
                            // this can be done because s0 and s1 are the same here and there
                            // 
                            /* 
                               Q[k] = ( s[8]*bQ[6][k] +
                               s[7]*bQ[5][k] + 
                               s[6]*bQ[4][k] + 
                               s[5]*bQ[3][k] + 
                               s[4]*bQ[2][k] + 
                               s[3]*bQ[1][k] + 
                               s[2]*bQ[0][k] ) +
                               s[1]*Q1DotDot[k] + 
                               s[0]*Q1Dot[k] + 
                               Q1[k];
                            */		 
                            // need to do this here because s[j] changes later
                            Q[k] = ( s[8]*bQ[6][k] +
                                     s[7]*bQ[5][k] + 
                                     s[6]*bQ[4][k] + 
                                     s[5]*bQ[3][k] + 
                                     s[4]*bQ[2][k] + 
                                     s[3]*bQ[1][k] + 
                                     s[2]*bQ[0][k] ) +
                                // s[1]*Q1DotDot[k] + 
                                // s[0]*Q1Dot[k] + 
                                Q1[k];
	      
                        }
                    }
	  
                    // }
                    // }
	  
                    ++bl_it;
                }
            }
      
            // ORSA_DEBUG("--MARK--");
      
            // needed only if using a velocity-dependent interaction...
            /* 
               if (interaction->depends_on_velocity()) {	
               // s[0] = timestep * h[j];
               s[0] = timestep.Getdouble() * h[j];
               s[1] = s[0] * h[j] * 0.5;
               s[2] = s[1] * h[j] * 0.6666666666666667;
               s[3] = s[2] * h[j] * 0.75;
               s[4] = s[3] * h[j] * 0.8;
               s[5] = s[4] * h[j] * 0.8333333333333333;
               s[6] = s[5] * h[j] * 0.8571428571428571;
               s[7] = s[6] * h[j] * 0.875;
	 
               for(k=0;k<nv;++k) {
               if (interaction->IsSkippingJPLPlanets() && frame_in[k/3].JPLPlanet() != NONE) continue;
	 
               v[k] = ( s[7]*b[6][k] + 
               s[6]*b[5][k] + 
               s[5]*b[4][k] + 
               s[4]*b[3][k] + 
               s[3]*b[2][k] + 
               s[2]*b[1][k] +
               s[1]*b[0][k] ) +
               s[0]*a1[k] + 
               v1[k];
               }
               }
            */
            //
            {
                // if (bg->getInteraction()->dependsOnVelocity()) {	
	
                s[0] = timestep.get_d() * h[j];
                s[1] = s[0] * h[j] * 0.5;
                s[2] = s[1] * h[j] * 0.6666666666666667;
                s[3] = s[2] * h[j] * 0.75;
                s[4] = s[3] * h[j] * 0.8;
                s[5] = s[4] * h[j] * 0.8333333333333333;
                s[6] = s[5] * h[j] * 0.8571428571428571;
                s[7] = s[6] * h[j] * 0.875;
	
                BodyGroup::BodyList::const_iterator bl_it = bl.begin();
                while (bl_it != bl.end()) {
	  
                    /* 
                       if (!((*bl_it)->getInitialConditions().translational->dynamic())) {
                       // ORSA_DEBUG("skipping body [%s]",(*bl_it)->getName().c_str());
                       ++bl_it;
                       continue;
                       }
                    */
	  
                    const orsa::Body * k = (*bl_it).get();
	  
                    if (bg->getInteraction()->dependsOnVelocity()) {
                        if ((*bl_it)->getInitialConditions().translational.get()) {
                            if ((*bl_it)->getInitialConditions().translational->dynamic()) {
		
                                v[k] = ( s[7]*b[6][k] + 
                                         s[6]*b[5][k] + 
                                         s[5]*b[4][k] + 
                                         s[4]*b[3][k] + 
                                         s[3]*b[2][k] + 
                                         s[2]*b[1][k] +
                                         s[1]*b[0][k] ) +
                                    s[0]*a1[k] + 
                                    v1[k];
                            }
                        }
                    } 
	  
                    if ((*bl_it)->getInitialConditions().rotational.get()) {
                        if ((*bl_it)->getInitialConditions().rotational->dynamic()) {
	      
                            // moved __HERE__ from __THERE__
                            /* 
                               Q[k] = ( s[8]*bQ[6][k] +
                               s[7]*bQ[5][k] + 
                               s[6]*bQ[4][k] + 
                               s[5]*bQ[3][k] + 
                               s[4]*bQ[2][k] + 
                               s[3]*bQ[1][k] + 
                               s[2]*bQ[0][k] ) +
                               s[1]*Q1DotDot[k] + 
                               s[0]*Q1Dot[k] + 
                               Q1[k];
                            */
                            //
                            /* 
                               QDot[k] = ( s[7]*bQ[6][k] + 
                               s[6]*bQ[5][k] + 
                               s[5]*bQ[4][k] + 
                               s[4]*bQ[3][k] + 
                               s[3]*bQ[2][k] + 
                               s[2]*bQ[1][k] +
                               s[1]*bQ[0][k] ) +
                               s[0]*Q1DotDot[k] + 
                               Q1Dot[k];
                            */
                            //
                            {
		
                                const orsa::Quaternion local_Q = unitQuaternion(Q[k]);
		
                                const orsa::Quaternion tmp_QDot = 
                                    ( s[7]*bQ[6][k] + 
                                      s[6]*bQ[5][k] + 
                                      s[5]*bQ[4][k] + 
                                      s[4]*bQ[3][k] + 
                                      s[3]*bQ[2][k] + 
                                      s[2]*bQ[1][k] +
                                      s[1]*bQ[0][k] ) +
                                    // s[0]*Q1DotDot[k] + 
                                    Q1Dot[k];
		
                                // important constraint on qDot!
                                const double delta = 
                                    local_Q.getScalar()*tmp_QDot.getScalar() +
                                    local_Q.getVector()*tmp_QDot.getVector();
		
                                const orsa::Quaternion local_QDot    = tmp_QDot - local_Q*delta;
		
                                const orsa::Quaternion local_QDotDot = Q1DotDot[k];
		
                                const orsa::Vector omega    = RotationalBodyProperty::omega(local_Q,
                                                                                            local_QDot);
		
                                const orsa::Vector omegaDot = RotationalBodyProperty::omegaDot(local_Q,
                                                                                               local_QDot,
                                                                                               local_QDotDot);
                                //
                                // const orsa::Vector omegaDot(0,0,0);
		
                                const orsa::Vector newOmega = RotationalBodyProperty::newOmega(omega,
                                                                                               omegaDot,
                                                                                               timestep);
		
                                Q[k]    = unitQuaternion(RotationalBodyProperty::qFiniteRotation(local_Q,
                                                                                                 newOmega,
                                                                                                 timestep));
		
                                QDot[k] = RotationalBodyProperty::qDot(Q[k],
                                                                       newOmega);
		
                                // important constraint on qDot!
                                const double finalDelta = 
                                    Q[k].getScalar()*QDot[k].getScalar() +
                                    Q[k].getVector()*QDot[k].getVector();	
                                // 
                                QDot[k] -= Q[k]*finalDelta;
		
                            }
	      
                        }
                    }
	  
                    ++bl_it;
	  
                }
            }
      
            // ORSA_DEBUG("--MARK--");
      
            /* 
               {
               Vector rr,vv,drr,dvv;
               for(k=0;k<frame_out.size();++k) {
               if (interaction->IsSkippingJPLPlanets() && frame_in[k].JPLPlanet() != NONE) continue;
	 
               frame_out[k] = frame_in[k];
	 
               rr.x = x[3*k];
               rr.y = x[3*k+1];
               rr.z = x[3*k+2];
	 
               drr = rr - frame_in[k].position();
               frame_out[k].AddToPosition(drr);
	 
               vv.x = v[3*k];
               vv.y = v[3*k+1];
               vv.z = v[3*k+2];
	 
               dvv = vv - frame_in[k].velocity();
               frame_out[k].AddToVelocity(dvv);
               }
               }
            */
            //
            {
                BodyGroup::BodyList::const_iterator bl_it = bl.begin();
                while (bl_it != bl.end()) {
	  
                    const orsa::Body * k = (*bl_it).get();
	  
                    IBPS ibps = k->getInitialConditions();
	  
                    if (!ibps.dynamic()) {
                        ++bl_it;
                        continue;
                    }
	  
                    ibps.time = start + orsa::Time(FromUnits(h[j]*timestep.get_d(),Unit::MICROSECOND,-1));
	  
                    if (!k->alive(ibps.time.getRef())) {
                        ++bl_it;
                        continue;
                    }
	  
                    if ((*bl_it)->getInitialConditions().translational.get()) {
                        if ((*bl_it)->getInitialConditions().translational->dynamic()) {
	      
                            // ORSA_DEBUG("--MARK--");
                            // orsa::print(x[k]);
	      
                            ibps.translational->setPosition(x[k]);
                            ibps.translational->setVelocity(v[k]);
	      
                        }
                    }
	  
                    if ((*bl_it)->getInitialConditions().rotational.get()) {
                        if ((*bl_it)->getInitialConditions().rotational->dynamic()) {
	     
                            // ORSA_DEBUG("CODE NEEDED HERE!!");
                            //
                            /* 
                               ibps.rotational->set(phi[k],
                               theta[k],
                               psi[k],
                               phiDot[k],
                               thetaDot[k],
                               psiDot[k]);
                            */
                            //
                            // ORSA_DEBUG("omega...");
                            // print(RotationalBodyProperty::omega(Q[k],QDot[k]));
                            //
                            Q[k] = unitQuaternion(Q[k]);
	      
                            // important constraint on qDot!
                            const double delta = 
                                Q[k].getScalar()*QDot[k].getScalar() +
                                Q[k].getVector()*QDot[k].getVector();
                            //
                            QDot[k] -= Q[k]*delta;
	      
                            // print(RotationalBodyProperty::omega(Q[k],QDot[k]));
                            //
                            ibps.rotational->set(Q[k],
                                                 RotationalBodyProperty::omega(Q[k],QDot[k]));
	      
                        }
                    }
	  
                    ibps.tmp = true;
	  
                    bg->insertIBPS(ibps,k,onlyIfExtending.getRef(),true);
	  
                    ++bl_it;
                }
            }
      
            // ORSA_DEBUG("--MARK--");
      
            /* 
               if (interaction->IsSkippingJPLPlanets()) {
               frame_out.SetTime(frame_in+timestep*h[j]);
               frame_out.ForceJPLEphemerisData();
               }
            */
            //
            // interaction->Acceleration(frame_out,acc);
            //
            /* 
               {
               BodyGroup::BodyList::const_iterator bl_it = bl.begin();
               while (bl_it != bl.end()) {
               // if ((*bl_it)->getBodyPosVelCallback() != 0) {
               if (!((*bl_it)->getInitialConditions().translational->dynamic())) {
               // ORSA_DEBUG("skipping body [%s]",(*bl_it)->getName().c_str());
               ++bl_it;
               continue;
               }
               //
               bg->getInteraction()->acceleration(acc[(*bl_it).get()],  
               (*bl_it).get(),
               bg,
               start+orsa::Time(FromUnits(h[j]*timestep.get_d(),Unit::MICROSECOND,-1)));
               ++bl_it;
               }
               }
            */
            //
            if (!bg->getInteraction()->acceleration(acc,  
                                                    bg,
                                                    start+orsa::Time(FromUnits(h[j]*timestep.get_d(),Unit::MICROSECOND,-1)))) {
                ORSA_DEBUG("problems...");
                return false;
            }
            //
            if (!bg->getInteraction()->torque(torque,  
                                              bg,
                                              start+orsa::Time(FromUnits(h[j]*timestep.get_d(),Unit::MICROSECOND,-1)))) {
                ORSA_DEBUG("problems...");
                return false;
            }
      
            // ORSA_DEBUG("--MARK--");
      
            /* 
               for(k=0;k<frame_out.size();++k) {
               if (interaction->IsSkippingJPLPlanets() && frame_in[k].JPLPlanet() != NONE) continue;
	 
               a[3*k]   = acc[k].x;
               a[3*k+1] = acc[k].y;  
               a[3*k+2] = acc[k].z;
               }
            */
            //
            {
                // BodyGroup::BodyList::const_iterator bl_it = bl.begin();
                // while (bl_it != bl.end()) {
	
                for (unsigned int bodyIndex=0; bodyIndex<bl.size(); ++bodyIndex) {
	  
                    const orsa::Body * b = bl[bodyIndex].get();
	  
                    // should check if alive?
	  
                    /* 
                       if (!((*bl_it)->getInitialConditions().translational->dynamic())) {
                       // ORSA_DEBUG("skipping body [%s]",(*bl_it)->getName().c_str());
                       ++bl_it;
                       continue;
                       }
                    */
	  
                    if (!b->alive(start)) {
                        // ++bl_it;
                        continue;
                    }
	  
                    IBPS ibps;
                    if (!bg->getInterpolatedIBPS(ibps,b,start+orsa::Time(FromUnits(h[j]*timestep.get_d(),Unit::MICROSECOND,-1)))) {
                        ORSA_DEBUG("problem, body: [%s]",b->getName().c_str());
                    }
	  
                    const orsa::Body * k = b;
	  
                    if (b->getInitialConditions().translational.get()) {
                        if (b->getInitialConditions().translational->dynamic()) {
                            a[k] = acc[bodyIndex];
                        }
                    }
	  
                    if (b->getInitialConditions().rotational.get()) {
                        if (b->getInitialConditions().rotational->dynamic()) {
	      
                            // #warning "inertia moments needed!! plus rotation to get to the principal axis..."
	      
                            // if (b->getPaulMoment()) {
                            if (ibps.inertial->paulMoment()) {
		
                                double m;
                                if (!bg->getInterpolatedMass(m,b,start+orsa::Time(FromUnits(h[j]*timestep.get_d(),Unit::MICROSECOND,-1)))) {
                                    ORSA_DEBUG("problems...");
                                }		
		
                                orsa::Matrix inertiaMoment = 
                                    m * 
                                    ibps.inertial->inertiaMatrix();
                                // b->getPaulMoment()->getInertiaMoment();
		
                                // ORSA_DEBUG("CODE NEEDED HERE!!");
                                //
                                /* 
                                   double localPhiDotDot, localThetaDotDot, localPsiDotDot;
		   
                                   orsa::Euler(localPhiDotDot,
                                   localThetaDotDot,
                                   localPsiDotDot,
                                   phi[k],
                                   theta[k],
                                   psi[k],
                                   phiDot[k],
                                   thetaDot[k],
                                   psiDot[k],
                                   inertiaMoment.getM11(),
                                   inertiaMoment.getM22(),
                                   inertiaMoment.getM33(),
                                   torque[k].getRef().getX(),
                                   torque[k].getRef().getY(),
                                   torque[k].getRef().getZ());
		   
                                   phiDotDot[k]   = localPhiDotDot;
                                   thetaDotDot[k] = localThetaDotDot;
                                   psiDotDot[k]   = localPsiDotDot;
                                */
                                //
                                {
		  
                                    // const orsa::Matrix g2l = BodyAttitude((*bl_it).get(),bg).globalToLocal(start+orsa::Time(FromUnits(h[j]*timestep.get_d(),Unit::MICROSECOND,-1)));
                                    // const orsa::Matrix l2g = BodyAttitude((*bl_it).get(),bg).localToGlobal(start+orsa::Time(FromUnits(h[j]*timestep.get_d(),Unit::MICROSECOND,-1)));
		  
                                    // osg::ref_ptr<orsa::Attitude> attitude = new orsa::BodyAttitude((*bl_it).get(),bg);
		  
                                    // const orsa::Matrix g2l = attitude->globalToLocal(start+orsa::Time(FromUnits(h[j]*timestep.get_d(),Unit::MICROSECOND,-1)));
                                    // const orsa::Matrix l2g = attitude->localToGlobal(start+orsa::Time(FromUnits(h[j]*timestep.get_d(),Unit::MICROSECOND,-1)));
		  
                                    const orsa::Matrix g2l = orsa::globalToLocal(b,bg,start+orsa::Time(FromUnits(h[j]*timestep.get_d(),Unit::MICROSECOND,-1)));
                                    const orsa::Matrix l2g = orsa::localToGlobal(b,bg,start+orsa::Time(FromUnits(h[j]*timestep.get_d(),Unit::MICROSECOND,-1)));
		  
                                    const orsa::Vector omega = RotationalBodyProperty::omega(Q[k],QDot[k]);
                                    const orsa::Matrix I     = inertiaMoment;
                                    const orsa::Vector T     = torque[bodyIndex];
		  
                                    orsa::Vector omegaDot;
		  
                                    Euler(omegaDot,
                                          g2l*omega,
                                          I,
                                          g2l*T);
                                    //
                                    omegaDot = l2g * omegaDot;
		  
                                    QDotDot[k] = RotationalBodyProperty::qDotDot(Q[k], omega, omegaDot);
                                }
		
                            } else {
                                ORSA_DEBUG("PaulMoment not available for body [%s]",
                                           b->getName().c_str());
                            }
	      
                        }
                    }
	  
                    // ++bl_it;
                }
            }
      
            orsa::Quaternion tmpQ, gkQ;
            orsa::Vector     tmpV, gkV;
            // double     tmpD, gkD;
      
            switch (j) {
                case 1: 
                    /* 
                       for(k=0;k<nv;++k) {
                       if (interaction->IsSkippingJPLPlanets() && frame_in[k/3].JPLPlanet() != NONE) continue;
                       tmp = g[0][k];
                       g[0][k]  = (a[k] - a1[k]) * r[0];
                       b[0][k] += g[0][k] - tmp;
                       }
                    */
                    //
                {
                    BodyGroup::BodyList::const_iterator bl_it = bl.begin();
                    while (bl_it != bl.end()) {
                        // if ((*bl_it)->getBodyPosVelCallback() != 0) {
                        /* 
                           if (!((*bl_it)->getInitialConditions().translational->dynamic())) {
                           // ORSA_DEBUG("skipping body [%s]",(*bl_it)->getName().c_str());
                           ++bl_it;
                           continue;
                           }
                        */
	    
                        const orsa::Body * k = (*bl_it).get();
	    
                        if ((*bl_it)->getInitialConditions().translational.get()) {
                            if ((*bl_it)->getInitialConditions().translational->dynamic()) {
		
                                tmpV     = g[0][k];
                                g[0][k]  = (a[k] - a1[k]) * r[0];
                                b[0][k] += g[0][k] - tmpV;
                            }
                        }
	    
                        if ((*bl_it)->getInitialConditions().rotational.get()) {
                            if ((*bl_it)->getInitialConditions().rotational->dynamic()) {
		
                                /* 
                                   tmpD        = gPhi[0][k];
                                   gPhi[0][k]  = (phiDotDot[k] - phi1DotDot[k]) * r[0];
                                   bPhi[0][k] += gPhi[0][k] - tmpD;
                                */
                                // 
                                tmpQ      = gQ[0][k];
                                gQ[0][k]  = (QDotDot[k] - Q1DotDot[k]) * r[0];
                                bQ[0][k] += gQ[0][k] - tmpQ;
		
                            }
                        }
	    
                        ++bl_it;
                    }
                }
                //
                break;
                case 2: 
                    /* 
                       for(k=0;k<nv;++k) {
                       if (interaction->IsSkippingJPLPlanets() && frame_in[k/3].JPLPlanet() != NONE) continue;
                       tmp = g[1][k];
                       gk = a[k] - a1[k];
                       g[1][k] = (gk*r[1] - g[0][k])*r[2];
                       tmp = g[1][k] - tmp;
                       b[0][k] += tmp * c[0];
                       b[1][k] += tmp;
                       }
                    */
                    //
                {
                    BodyGroup::BodyList::const_iterator bl_it = bl.begin();
                    while (bl_it != bl.end()) {
	
                        /* 
                           if (!((*bl_it)->getInitialConditions().translational->dynamic())) {
                           // ORSA_DEBUG("skipping body [%s]",(*bl_it)->getName().c_str());
                           ++bl_it;
                           continue;
                           }
                        */
	    
                        const orsa::Body * k = (*bl_it).get();
	    
                        /* 
                           tmpV      = g[1][k];
                           gkV       = a[k] - a1[k];
                           g[1][k]  = (gkV*r[1] - g[0][k])*r[2];
                           tmpV      = g[1][k] - tmpV;
                           b[0][k] += tmpV * c[0];
                           b[1][k] += tmpV;
                        */
	    
                        if ((*bl_it)->getInitialConditions().translational.get()) {
                            if ((*bl_it)->getInitialConditions().translational->dynamic()) {
		
                                tmpV      = g[1][k];
                                gkV       = a[k] - a1[k];
                                g[1][k]  = (gkV*r[1] - g[0][k])*r[2];
                                tmpV      = g[1][k] - tmpV;
                                b[0][k] += tmpV * c[0];
                                b[1][k] += tmpV;
                            }
                        }
	    
                        if ((*bl_it)->getInitialConditions().rotational.get()) {
                            if ((*bl_it)->getInitialConditions().rotational->dynamic()) {
		
                                /* 
                                   tmpD        = gPhi[1][k];
                                   gkD         = phiDotDot[k] - phi1DotDot[k];
                                   gPhi[1][k]  = (gkD*r[1] - gPhi[0][k])*r[2];
                                   tmpD        = gPhi[1][k] - tmpD;
                                   bPhi[0][k] += tmpD * c[0];
                                   bPhi[1][k] += tmpD;
                                */
                                //
                                tmpQ      = gQ[1][k];
                                gkQ       = QDotDot[k] - Q1DotDot[k];
                                gQ[1][k]  = (gkQ*r[1] - gQ[0][k])*r[2];
                                tmpQ      = gQ[1][k] - tmpQ;
                                bQ[0][k] += tmpQ * c[0];
                                bQ[1][k] += tmpQ;
		
                            }
                        }
	    
                        ++bl_it;
                    }
                }
                //
                break;
                case 3: 
                    /* 
                       for(k=0;k<nv;++k) {
                       if (interaction->IsSkippingJPLPlanets() && frame_in[k/3].JPLPlanet() != NONE) continue;
                       tmp = g[2][k];
                       gk = a[k] - a1[k];
                       g[2][k] = ((gk*r[3] - g[0][k])*r[4] - g[1][k])*r[5];
                       tmp = g[2][k] - tmp;
                       b[0][k] += tmp * c[1];
                       b[1][k] += tmp * c[2];
                       b[2][k] += tmp;
                       } 
                    */
                    //
                {
                    BodyGroup::BodyList::const_iterator bl_it = bl.begin();
                    while (bl_it != bl.end()) {
	   
                        /* 
                           if (!((*bl_it)->getInitialConditions().translational->dynamic())) {
                           // ORSA_DEBUG("skipping body [%s]",(*bl_it)->getName().c_str());
                           ++bl_it;
                           continue;
                           }
                        */
	    
                        const orsa::Body * k = (*bl_it).get();
	    
                        if ((*bl_it)->getInitialConditions().translational.get()) {
                            if ((*bl_it)->getInitialConditions().translational->dynamic()) {
		
                                tmpV      = g[2][k];
                                gkV       = a[k] - a1[k];
                                g[2][k]   = ((gkV*r[3] - g[0][k])*r[4] - g[1][k])*r[5];
                                tmpV      = g[2][k] - tmpV;
                                b[0][k] += tmpV * c[1];
                                b[1][k] += tmpV * c[2];
                                b[2][k] += tmpV;
                            }
                        }
	    
                        if ((*bl_it)->getInitialConditions().rotational.get()) {
                            if ((*bl_it)->getInitialConditions().rotational->dynamic()) {
		
                                /* 
                                   tmpD        = gPhi[2][k];
                                   gkD         = phiDotDot[k] - phi1DotDot[k];
                                   gPhi[2][k]  = ((gkD*r[3] - gPhi[0][k])*r[4] - gPhi[1][k])*r[5];
                                   tmpD        = gPhi[2][k] - tmpD;
                                   bPhi[0][k] += tmpD * c[1];
                                   bPhi[1][k] += tmpD * c[2];
                                   bPhi[2][k] += tmpD;
                                */
                                //
                                tmpQ      = gQ[2][k];
                                gkQ       = QDotDot[k] - Q1DotDot[k];
                                gQ[2][k]  = ((gkQ*r[3] - gQ[0][k])*r[4] - gQ[1][k])*r[5];
                                tmpQ      = gQ[2][k] - tmpQ;
                                bQ[0][k] += tmpQ * c[1];
                                bQ[1][k] += tmpQ * c[2];
                                bQ[2][k] += tmpQ;
		
                            }
                        }
	    
                        ++bl_it;
                    }
                }
                //
                break;
                case 4:
                    /* 
                       for(k=0;k<nv;++k) {
                       if (interaction->IsSkippingJPLPlanets() && frame_in[k/3].JPLPlanet() != NONE) continue;
                       tmp = g[3][k];
                       gk = a[k] - a1[k];
                       g[3][k] = (((gk*r[6] - g[0][k])*r[7] - g[1][k])*r[8] - g[2][k])*r[9];
                       tmp = g[3][k] - tmp;
                       b[0][k] += tmp * c[3];
                       b[1][k] += tmp * c[4];
                       b[2][k] += tmp * c[5];
                       b[3][k] += tmp;
                       }
                    */
                    //
                {
                    BodyGroup::BodyList::const_iterator bl_it = bl.begin();
                    while (bl_it != bl.end()) {
	    
                        /* 
                           if (!((*bl_it)->getInitialConditions().translational->dynamic())) {
                           // ORSA_DEBUG("skipping body [%s]",(*bl_it)->getName().c_str());
                           ++bl_it;
                           continue;
                           }
                        */
	    
                        const orsa::Body * k = (*bl_it).get();
	    
                        if ((*bl_it)->getInitialConditions().translational.get()) {
                            if ((*bl_it)->getInitialConditions().translational->dynamic()) {
		
                                tmpV     = g[3][k];
                                gkV      = a[k] - a1[k];
                                g[3][k]  = (((gkV*r[6] - g[0][k])*r[7] - g[1][k])*r[8] - g[2][k])*r[9];
                                tmpV     = g[3][k] - tmpV;
                                b[0][k] += tmpV * c[3];
                                b[1][k] += tmpV * c[4];
                                b[2][k] += tmpV * c[5];
                                b[3][k] += tmpV;
                            }
                        }
	    
                        if ((*bl_it)->getInitialConditions().rotational.get()) {
                            if ((*bl_it)->getInitialConditions().rotational->dynamic()) {
		
                                /* 
                                   tmpD        = gPhi[3][k];
                                   gkD         = phiDotDot[k] - phi1DotDot[k];
                                   gPhi[3][k]  = (((gkD*r[6] - gPhi[0][k])*r[7] - gPhi[1][k])*r[8] - gPhi[2][k])*r[9];
                                   tmpD        = gPhi[3][k] - tmpD;
                                   bPhi[0][k] += tmpD * c[3];
                                   bPhi[1][k] += tmpD * c[4];
                                   bPhi[2][k] += tmpD * c[5];
                                   bPhi[3][k] += tmpD;
                                */
                                //
                                tmpQ      = gQ[3][k];
                                gkQ       = QDotDot[k] - Q1DotDot[k];
                                gQ[3][k]  = (((gkQ*r[6] - gQ[0][k])*r[7] - gQ[1][k])*r[8] - gQ[2][k])*r[9];
                                tmpQ      = gQ[3][k] - tmpQ;
                                bQ[0][k] += tmpQ * c[3];
                                bQ[1][k] += tmpQ * c[4];
                                bQ[2][k] += tmpQ * c[5];
                                bQ[3][k] += tmpQ;
		
                            }
                        }
	    
                        ++bl_it;
                    }
                }
                //
                break;
                case 5:
                    /* 
                       for(k=0;k<nv;++k) {
                       if (interaction->IsSkippingJPLPlanets() && frame_in[k/3].JPLPlanet() != NONE) continue;
                       tmp = g[4][k];
                       gk = a[k] - a1[k];
                       g[4][k] = ((((gk*r[10] - g[0][k])*r[11] - g[1][k])*r[12] - g[2][k])*r[13] - g[3][k])*r[14];
                       tmp = g[4][k] - tmp;
                       b[0][k] += tmp * c[6];
                       b[1][k] += tmp * c[7];
                       b[2][k] += tmp * c[8];
                       b[3][k] += tmp * c[9];
                       b[4][k] += tmp;
                       } 
                    */
                    //
                {
                    BodyGroup::BodyList::const_iterator bl_it = bl.begin();
                    while (bl_it != bl.end()) {
	    
                        /* 
                           if (!((*bl_it)->getInitialConditions().translational->dynamic())) {
                           // ORSA_DEBUG("skipping body [%s]",(*bl_it)->getName().c_str());
                           ++bl_it;
                           continue;
                           }
                        */
	    
                        const orsa::Body * k = (*bl_it).get();
	    
                        if ((*bl_it)->getInitialConditions().translational.get()) {
                            if ((*bl_it)->getInitialConditions().translational->dynamic()) {

                                tmpV     = g[4][k];
                                gkV      = a[k] - a1[k];
                                g[4][k]  = ((((gkV*r[10] - g[0][k])*r[11] - g[1][k])*r[12] - g[2][k])*r[13] - g[3][k])*r[14];
                                tmpV     = g[4][k] - tmpV;
                                b[0][k] += tmpV * c[6];
                                b[1][k] += tmpV * c[7];
                                b[2][k] += tmpV * c[8];
                                b[3][k] += tmpV * c[9];
                                b[4][k] += tmpV;
                            }
                        }
	    
                        if ((*bl_it)->getInitialConditions().rotational.get()) {
                            if ((*bl_it)->getInitialConditions().rotational->dynamic()) {
		
                                /* 
                                   tmpD        = gPhi[4][k];
                                   gkD         = phiDotDot[k] - phi1DotDot[k];
                                   gPhi[4][k]  = ((((gkD*r[10] - gPhi[0][k])*r[11] - gPhi[1][k])*r[12] - gPhi[2][k])*r[13] - gPhi[3][k])*r[14];
                                   tmpD        = gPhi[4][k] - tmpD;
                                   bPhi[0][k] += tmpD * c[6];
                                   bPhi[1][k] += tmpD * c[7];
                                   bPhi[2][k] += tmpD * c[8];
                                   bPhi[3][k] += tmpD * c[9];
                                   bPhi[4][k] += tmpD;
                                */
                                //
                                tmpQ      = gQ[4][k];
                                gkQ       = QDotDot[k] - Q1DotDot[k];
                                gQ[4][k]  = ((((gkQ*r[10] - gQ[0][k])*r[11] - gQ[1][k])*r[12] - gQ[2][k])*r[13] - gQ[3][k])*r[14];
                                tmpQ      = gQ[4][k] - tmpQ;
                                bQ[0][k] += tmpQ * c[6];
                                bQ[1][k] += tmpQ * c[7];
                                bQ[2][k] += tmpQ * c[8];
                                bQ[3][k] += tmpQ * c[9];
                                bQ[4][k] += tmpQ;
		
                            }
                        }
	    
                        ++bl_it;
                    }
                }
                //
                break;
                case 6:
                    /* 
                       for(k=0;k<nv;++k) {
                       if (interaction->IsSkippingJPLPlanets() && frame_in[k/3].JPLPlanet() != NONE) continue;
                       tmp = g[5][k];
                       gk = a[k] - a1[k];
                       g[5][k] = (((((gk*r[15] - g[0][k])*r[16] - g[1][k])*r[17] - g[2][k])*r[18] - g[3][k])*r[19] - g[4][k])*r[20];
                       tmp = g[5][k] - tmp;
                       b[0][k] += tmp * c[10];
                       b[1][k] += tmp * c[11];
                       b[2][k] += tmp * c[12];
                       b[3][k] += tmp * c[13];
                       b[4][k] += tmp * c[14];
                       b[5][k] += tmp;
                       }
                    */
                    //
                {
                    BodyGroup::BodyList::const_iterator bl_it = bl.begin();
                    while (bl_it != bl.end()) {
	   
                        /* 
                           if (!((*bl_it)->getInitialConditions().translational->dynamic())) {
                           // ORSA_DEBUG("skipping body [%s]",(*bl_it)->getName().c_str());
                           ++bl_it;
                           continue;
                           }
                        */
	    
                        const orsa::Body * k = (*bl_it).get();
	    
                        if ((*bl_it)->getInitialConditions().translational.get()) {
                            if ((*bl_it)->getInitialConditions().translational->dynamic()) {
		
                                tmpV     = g[5][k];
                                gkV      = a[k] - a1[k];
                                g[5][k]  = (((((gkV*r[15] - g[0][k])*r[16] - g[1][k])*r[17] - g[2][k])*r[18] - g[3][k])*r[19] - g[4][k])*r[20];
                                tmpV     = g[5][k] - tmpV;
                                b[0][k] += tmpV * c[10];
                                b[1][k] += tmpV * c[11];
                                b[2][k] += tmpV * c[12];
                                b[3][k] += tmpV * c[13];
                                b[4][k] += tmpV * c[14];
                                b[5][k] += tmpV;
                            }
                        }
	    
                        if ((*bl_it)->getInitialConditions().rotational.get()) {
                            if ((*bl_it)->getInitialConditions().rotational->dynamic()) {
		
                                /* 
                                   tmpD        = gPhi[5][k];
                                   gkD         = phiDotDot[k] - phi1DotDot[k];
                                   gPhi[5][k]  = (((((gkD*r[15] - gPhi[0][k])*r[16] - gPhi[1][k])*r[17] - gPhi[2][k])*r[18] - gPhi[3][k])*r[19] - gPhi[4][k])*r[20];
                                   tmpD        = gPhi[5][k] - tmpD;
                                   bPhi[0][k] += tmpD * c[10];
                                   bPhi[1][k] += tmpD * c[11];
                                   bPhi[2][k] += tmpD * c[12];
                                   bPhi[3][k] += tmpD * c[13];
                                   bPhi[4][k] += tmpD * c[14];
                                   bPhi[5][k] += tmpD;
                                */
                                //
                                tmpQ      = gQ[5][k];
                                gkQ       = QDotDot[k] - Q1DotDot[k];
                                gQ[5][k]  = (((((gkQ*r[15] - gQ[0][k])*r[16] - gQ[1][k])*r[17] - gQ[2][k])*r[18] - gQ[3][k])*r[19] - gQ[4][k])*r[20];
                                tmpQ      = gQ[5][k] - tmpQ;
                                bQ[0][k] += tmpQ * c[10];
                                bQ[1][k] += tmpQ * c[11];
                                bQ[2][k] += tmpQ * c[12];
                                bQ[3][k] += tmpQ * c[13];
                                bQ[4][k] += tmpQ * c[14];
                                bQ[5][k] += tmpQ;
		
                            }
                        }
	    
                        ++bl_it;
                    }
                }
                //
                break;
                case 7:
                    /* 
                       for(k=0;k<nv;++k) {
                       if (interaction->IsSkippingJPLPlanets() && frame_in[k/3].JPLPlanet() != NONE) continue;
                       tmp = g[6][k];
                       gk = a[k] - a1[k];
                       g[6][k] = ((((((gk*r[21] - g[0][k])*r[22] - g[1][k])*r[23] - g[2][k])*r[24] - g[3][k])*r[25] - g[4][k])*r[26] - g[5][k])*r[27];
                       tmp = g[6][k] - tmp;
                       b[0][k] += tmp * c[15];
                       b[1][k] += tmp * c[16];
                       b[2][k] += tmp * c[17];
                       b[3][k] += tmp * c[18];
                       b[4][k] += tmp * c[19];
                       b[5][k] += tmp * c[20];
                       b[6][k] += tmp;
                       } 
                    */
                    //
                {
                    BodyGroup::BodyList::const_iterator bl_it = bl.begin();
                    while (bl_it != bl.end()) {
	  
                        /* 
                           if (!((*bl_it)->getInitialConditions().translational->dynamic())) {
                           // ORSA_DEBUG("skipping body [%s]",(*bl_it)->getName().c_str());
                           ++bl_it;
                           continue;
                           }
                        */
	    
                        const orsa::Body * k = (*bl_it).get();
	    
                        if ((*bl_it)->getInitialConditions().translational.get()) {
                            if ((*bl_it)->getInitialConditions().translational->dynamic()) {
		
                                tmpV     = g[6][k];
                                gkV      = a[k] - a1[k];
                                g[6][k]  = ((((((gkV*r[21] - g[0][k])*r[22] - g[1][k])*r[23] - g[2][k])*r[24] - g[3][k])*r[25] - g[4][k])*r[26] - g[5][k])*r[27];
                                tmpV     = g[6][k] - tmpV;
                                b[0][k] += tmpV * c[15];
                                b[1][k] += tmpV * c[16];
                                b[2][k] += tmpV * c[17];
                                b[3][k] += tmpV * c[18];
                                b[4][k] += tmpV * c[19];
                                b[5][k] += tmpV * c[20];
                                b[6][k] += tmpV;
                            }
                        }
	    
                        if ((*bl_it)->getInitialConditions().rotational.get()) {
                            if ((*bl_it)->getInitialConditions().rotational->dynamic()) {
		
                                /* 
                                   tmpD     = gPhi[6][k];
                                   gkD      = phiDotDot[k] - phi1DotDot[k];
                                   gPhi[6][k]  = ((((((gkD*r[21] - gPhi[0][k])*r[22] - gPhi[1][k])*r[23] - gPhi[2][k])*r[24] - gPhi[3][k])*r[25] - gPhi[4][k])*r[26] - gPhi[5][k])*r[27];
                                   tmpD     = gPhi[6][k] - tmpD;
                                   bPhi[0][k] += tmpD * c[15];
                                   bPhi[1][k] += tmpD * c[16];
                                   bPhi[2][k] += tmpD * c[17];
                                   bPhi[3][k] += tmpD * c[18];
                                   bPhi[4][k] += tmpD * c[19];
                                   bPhi[5][k] += tmpD * c[20];
                                   bPhi[6][k] += tmpD;
                                */
                                //
                                tmpQ      = gQ[6][k];
                                gkQ       = QDotDot[k] - Q1DotDot[k];
                                gQ[6][k]  = ((((((gkQ*r[21] - gQ[0][k])*r[22] - gQ[1][k])*r[23] - gQ[2][k])*r[24] - gQ[3][k])*r[25] - gQ[4][k])*r[26] - gQ[5][k])*r[27];
                                tmpQ      = gQ[6][k] - tmpQ;
                                bQ[0][k] += tmpQ * c[15];
                                bQ[1][k] += tmpQ * c[16];
                                bQ[2][k] += tmpQ * c[17];
                                bQ[3][k] += tmpQ * c[18];
                                bQ[4][k] += tmpQ * c[19];
                                bQ[5][k] += tmpQ * c[20];
                                bQ[6][k] += tmpQ;
		
                            }
                        }
	    
                        ++bl_it;
                    }
                }
                //
                break;
                default:
                    ORSA_ERROR("aieeee!!!");
            }
        }
    }
  
    // ORSA_DEBUG("--MARK--");
  
    // const orsa::Time timestep_done = timestep;
  
    // Estimate suitable sequence size for the next call
    double tmp = 0;
    /* 
       for(k=0;k<nv;++k) {
       if (interaction->IsSkippingJPLPlanets() && frame_in[k/3].JPLPlanet() != NONE) continue;
       tmp = MAX(tmp,fabs(b[6][k]));
       }
    */
    {
        BodyGroup::BodyList::const_iterator bl_it = bl.begin();
        while (bl_it != bl.end()) {
      
            /* 
               if (!((*bl_it)->getInitialConditions().translational->dynamic())) {
               // ORSA_DEBUG("skipping body [%s]",(*bl_it)->getName().c_str());
               ++bl_it;
               continue;
               }
            */
      
            const orsa::Body * k = (*bl_it).get();
      
            if ((*bl_it)->getInitialConditions().translational.get()) {
                if ((*bl_it)->getInitialConditions().translational->dynamic()) {
	  
                    tmp = std::max(tmp,fabs(b[6][k].getX()));
                    tmp = std::max(tmp,fabs(b[6][k].getY()));
                    tmp = std::max(tmp,fabs(b[6][k].getZ()));
                }
            }
      
            if ((*bl_it)->getInitialConditions().rotational.get()) {
                if ((*bl_it)->getInitialConditions().rotational->dynamic()) {
	  
                    /* 
                       tmp = std::max(tmp,fabs(bPhi[6][k]));
                       tmp = std::max(tmp,fabs(bTheta[6][k]));
                       tmp = std::max(tmp,fabs(bPsi[6][k]));
                    */
                    //
                    // #warning "check this length()..."
                    
                    /* ORSA_DEBUG("tmp test...");
                       print(tmp);
                       print(fabs(bQ[6][k].length()));
                    */
                    //
                    tmp = std::max(tmp,fabs(bQ[6][k].length()));
                }
            }
      
            ++bl_it;
        }
    }
  
    // ORSA_DEBUG("1: tmp: %Fg",tmp());
  
    // if (tmp!=0.0) tmp /= (72.0 * secure_pow(fabs(timestep),7));
    // if (tmp!=0.0) tmp /= (72.0 * secure_pow(fabs(timestep.Getdouble()),7));
    //
    // if (tmp!=0.0) tmp /= (72.0 * pow(fabs(timestep.Getdouble()),7));
    //
    if (tmp != 0) tmp /= (72*int_pow(fabs(timestep.get_d()),7));
  
    // ORSA_DEBUG("2: tmp: %Fg",tmp());
  
    // if (tmp < 1.0e-50) { // is equal to zero?
    // if (tmp < epsilon()) { // is equal to zero?
    if (tmp == 0) {
        // ORSA_DEBUG("zero test...");
        // timestep = timestep_done * 1.4;
        next_timestep = orsa::Time(FromUnits(timestep.get_d()*1.4,Unit::MICROSECOND,-1));
    } else {
        // old rule...
        // timestep = copysign(secure_pow(accuracy/tmp,0.1111111111111111),timestep_done); // 1/9=0.111...
        // timestep = copysign(secure_pow(accuracy/tmp,0.1111111111111111),timestep_done.Getdouble()); // 1/9=0.111...
        //
        // timestep = copysign(pow(accuracy/tmp,0.1111111111111111),timestep_done.Getdouble()); // 1/9=0.111...
        //
    
        next_timestep = orsa::Time(FromUnits(copysign(pow(_accuracy.getRef()/tmp, 1.0/9.0), timestep.get_d()),Unit::MICROSECOND,-1)); // 1/9=0.111...
    
        // next_timestep = orsa::Time(copysign(pow(_accuracy.getRef()/tmp, 1/double("9.0")), FromUnits(timestep.get_d(),Unit::MICROSECOND,-1))); // 1/9=0.111...
    }
  
    /* 
       ORSA_DEBUG("proposed next_timestep: %f [day]     old timestep:  %f [day]",
       FromUnits(FromUnits(next_timestep.getMuSec(),Unit::MICROSECOND),Unit::DAY,-1)(),
       FromUnits(FromUnits(timestep.getMuSec(),Unit::MICROSECOND),Unit::DAY,-1)());
    */
  
    // check this only on "initialization" steps, i.e. when niter > 2
    if ((niter > 2) && (fabs(next_timestep.get_d()/timestep.get_d()) < 1.0)) {
        next_timestep = orsa::Time(FromUnits(timestep.get_d()*0.8,Unit::MICROSECOND,-1));
        // std::cerr << "Radau: step rejected! New proposed timestep: " << timestep.Getdouble() << std::endl;
        // frame_out = frame_in;
        _lastCallRejected = true;
        // niter = 6;
        // std::cerr << "[rej]" << std::endl;
    
        /* 
           ORSA_DEBUG("REJECTED, next_timestep: %20.12f [day]",
           FromUnits(FromUnits(next_timestep.getMuSec(),Unit::MICROSECOND),Unit::DAY,-1)());
        */
    
        // ORSA_DEBUG("done.");
    
        return true;
    
    } else {
        // ORSA_DEBUG("->>>>>>>> CALL NOT REJECTED!!!");
        _lastCallRejected = false;
    }
  
    if (fabs(timestep.get_d()/timestep.get_d()) > 1.4) {
        // timestep = timestep_done * 1.4;
        next_timestep = orsa::Time(FromUnits(timestep.get_d()*1.4,Unit::MICROSECOND,-1));
    }
  
    // std::cerr << "RA15: new timestep: " << timestep.Getdouble() << std::endl;
  
    // Find new position and velocity values at end of the sequence
    tmp = timestep.get_d() * timestep.get_d();
    //
    /* 
       for(k=0;k<nv;++k) {
       if (interaction->IsSkippingJPLPlanets() && frame_in[k/3].JPLPlanet() != NONE) continue;
     
       x1[k] = ( xc[7]*b[6][k] +
       xc[6]*b[5][k] + 
       xc[5]*b[4][k] + 
       xc[4]*b[3][k] + 
       xc[3]*b[2][k] + 
       xc[2]*b[1][k] + 
       xc[1]*b[0][k] + 
       xc[0]*a1[k]   ) * tmp + v1[k]*timestep_done.Getdouble() + x1[k];
     
       v1[k] = ( vc[6]*b[6][k] + 
       vc[5]*b[5][k] + 
       vc[4]*b[4][k] +
       vc[3]*b[3][k] + 
       vc[2]*b[2][k] + 
       vc[1]*b[1][k] +
       vc[0]*b[0][k] + 
       a1[k])        * timestep_done.Getdouble() + v1[k];
       }
    */
    //
    {
        BodyGroup::BodyList::const_iterator bl_it = bl.begin();
        while (bl_it != bl.end()) {
      
            /* 
               if (!((*bl_it)->getInitialConditions().translational->dynamic())) {
               // ORSA_DEBUG("skipping body [%s]",(*bl_it)->getName().c_str());
               ++bl_it;
               continue;
               }
            */
      
            const orsa::Body * k = (*bl_it).get();
      
            if ((*bl_it)->getInitialConditions().translational.get()) {
                if ((*bl_it)->getInitialConditions().translational->dynamic()) {
	  
                    x1[k] = ( xc[7]*b[6][k] +
                              xc[6]*b[5][k] + 
                              xc[5]*b[4][k] + 
                              xc[4]*b[3][k] + 
                              xc[3]*b[2][k] + 
                              xc[2]*b[1][k] + 
                              xc[1]*b[0][k] + 
                              xc[0]*a1[k]) * tmp + v1[k]*timestep.get_d() + x1[k];
	  
                    v1[k] = ( vc[6]*b[6][k] + 
                              vc[5]*b[5][k] + 
                              vc[4]*b[4][k] +
                              vc[3]*b[3][k] + 
                              vc[2]*b[2][k] + 
                              vc[1]*b[1][k] +
                              vc[0]*b[0][k] + 
                              a1[k]) * timestep.get_d() + v1[k];
                }
            }
      
            if ((*bl_it)->getInitialConditions().rotational.get()) {
                if ((*bl_it)->getInitialConditions().rotational->dynamic()) {
	  
                    /* 
                       phi1[k] = ( xc[7]*bPhi[6][k] +
                       xc[6]*bPhi[5][k] + 
                       xc[5]*bPhi[4][k] + 
                       xc[4]*bPhi[3][k] + 
                       xc[3]*bPhi[2][k] + 
                       xc[2]*bPhi[1][k] + 
                       xc[1]*bPhi[0][k] + 
                       xc[0]*phi1DotDot[k]) * tmp + phi1Dot[k]*timestep.get_d() + phi1[k];
                    */
                    //
                    /* 
                       Q1[k] = ( xc[7]*bQ[6][k] +
                       xc[6]*bQ[5][k] + 
                       xc[5]*bQ[4][k] + 
                       xc[4]*bQ[3][k] + 
                       xc[3]*bQ[2][k] + 
                       xc[2]*bQ[1][k] + 
                       xc[1]*bQ[0][k] + 
                       xc[0]*Q1DotDot[k]) * tmp + Q1Dot[k]*timestep.get_d() + Q1[k];
                    */
	  
                    // dot
	  
                    /* 
                       phi1Dot[k] = ( vc[6]*bPhi[6][k] + 
                       vc[5]*bPhi[5][k] + 
                       vc[4]*bPhi[4][k] +
                       vc[3]*bPhi[3][k] + 
                       vc[2]*bPhi[2][k] + 
                       vc[1]*bPhi[1][k] +
                       vc[0]*bPhi[0][k] + 
                       phi1DotDot[k]) * timestep.get_d() + phi1Dot[k];
                    */
                    //
                    /* 
                       Q1Dot[k] = ( vc[6]*bQ[6][k] + 
                       vc[5]*bQ[5][k] + 
                       vc[4]*bQ[4][k] +
                       vc[3]*bQ[3][k] + 
                       vc[2]*bQ[2][k] + 
                       vc[1]*bQ[1][k] +
                       vc[0]*bQ[0][k] + 
                       Q1DotDot[k]) * timestep.get_d() + Q1Dot[k];
                    */
                    //
                    {
	    
                        Q1[k] = ( xc[7]*bQ[6][k] +
                                  xc[6]*bQ[5][k] + 
                                  xc[5]*bQ[4][k] + 
                                  xc[4]*bQ[3][k] + 
                                  xc[3]*bQ[2][k] + 
                                  xc[2]*bQ[1][k] + 
                                  xc[1]*bQ[0][k] + 
                                  xc[0]*Q1DotDot[k]) * tmp + 
                            // Q1Dot[k]*timestep.get_d() +
                            Q1[k];
	    
                        const orsa::Quaternion local_Q1       = unitQuaternion(Q1[k]);
	    
                        const orsa::Quaternion tmp_Q1Dot      = Q1Dot[k];
	
                        // important constraint on qDot!
                        const double delta = 
                            local_Q1.getScalar()*tmp_Q1Dot.getScalar() +
                            local_Q1.getVector()*tmp_Q1Dot.getVector();
	    
                        const orsa::Quaternion local_Q1Dot    = tmp_Q1Dot - local_Q1*delta;
                        /* 
                           const orsa::Quaternion local_Q1DotDot = ( vc[6]*bQ[6][k] + 
                           vc[5]*bQ[5][k] + 
                           vc[4]*bQ[4][k] +
                           vc[3]*bQ[3][k] + 
                           vc[2]*bQ[2][k] + 
                           vc[1]*bQ[1][k] +
                           vc[0]*bQ[0][k] + 
                           Q1DotDot[k]);
                        */
                        const orsa::Quaternion local_Q1DotDot = Q1DotDot[k];
	    
                        const orsa::Vector omega    = RotationalBodyProperty::omega(local_Q1,
                                                                                    local_Q1Dot);
	    
                        const orsa::Vector omegaDot = RotationalBodyProperty::omegaDot(local_Q1,
                                                                                       local_Q1Dot,
                                                                                       local_Q1DotDot);
                        //
                        // const orsa::Vector omegaDot(0,0,0);
	    
                        const orsa::Vector newOmega = RotationalBodyProperty::newOmega(omega,
                                                                                       omegaDot,
                                                                                       timestep);
	    
                        Q1[k] = unitQuaternion(RotationalBodyProperty::qFiniteRotation(local_Q1,
                                                                                       omega,
                                                                                       timestep));
	    
                        Q1Dot[k] = RotationalBodyProperty::qDot(Q1[k],
                                                                newOmega);
	    
                        // important constraint on qDot!
                        const double finalDelta = 
                            Q1[k].getScalar()*Q1Dot[k].getScalar() +
                            Q1[k].getVector()*Q1Dot[k].getVector();	
                        // 
                        Q1Dot[k] -= Q1[k]*finalDelta;  
	  
                    }
	  
                }
            }
      
            ++bl_it;
        }
    }
  
    // ORSA_DEBUG("--MARK--");
  
    /* 
       {
       Vector rr,vv,drr,dvv;
       for(k=0;k<frame_out.size();++k) {
       if (interaction->IsSkippingJPLPlanets() && frame_in[k].JPLPlanet() != NONE) continue;
     
       frame_out[k] = frame_in[k];
     
       rr.x = x1[3*k];
       rr.y = x1[3*k+1];
       rr.z = x1[3*k+2];
     
       drr = rr - frame_in[k].position();  
       frame_out[k].AddToPosition(drr);
     
       vv.x = v1[3*k];
       vv.y = v1[3*k+1];
       vv.z = v1[3*k+2];
     
       dvv = vv - frame_in[k].velocity();
       frame_out[k].AddToVelocity(dvv);
       }
       }
    */
    //
    {
        BodyGroup::BodyList::const_iterator bl_it = bl.begin();
        while (bl_it != bl.end()) {
      
            const orsa::Body * k = (*bl_it).get();
      
            if (!k->alive(start+timestep)) {
                ++bl_it;
                continue;
            }
      
            IBPS ibps = k->getInitialConditions();
      
            if (!ibps.dynamic()) {
                ++bl_it;
                continue;
            }
      
            ibps.time = start + timestep;
      
            if ((*bl_it)->getInitialConditions().translational.get()) {
                if ((*bl_it)->getInitialConditions().translational->dynamic()) {
	  
                    // ORSA_DEBUG("--MARK--");
                    // orsa::print(x1[k]);
	  
                    ibps.translational->setPosition(x1[k]);
                    ibps.translational->setVelocity(v1[k]);
	  
                }
            }
      
            if ((*bl_it)->getInitialConditions().rotational.get()) {
                if ((*bl_it)->getInitialConditions().rotational->dynamic()) {
	  
                    // ORSA_DEBUG("CODE NEEDED HERE!!");
                    //
                    /* 
                       ibps.rotational->set(phi1[k],
                       theta1[k],
                       psi1[k],
                       phi1Dot[k],
                       theta1Dot[k],
                       psi1Dot[k]); 
                    */
                    //
                    // ORSA_DEBUG("omega...");
                    // print(RotationalBodyProperty::omega(Q1[k],Q1Dot[k]));
                    //	  
                    Q1[k] = unitQuaternion(Q1[k]);
	  
                    // important constraint on qDot!
                    const double delta = 
                        Q1[k].getScalar()*Q1Dot[k].getScalar() +
                        Q1[k].getVector()*Q1Dot[k].getVector();
                    //
                    Q1Dot[k] -= Q1[k]*delta;
	  
                    // print(RotationalBodyProperty::omega(Q1[k],Q1Dot[k]));
                    //
                    ibps.rotational->set(Q1[k],
                                         RotationalBodyProperty::omega(Q1[k],Q1Dot[k]));
	  
                }
            }
      
            ibps.tmp = false;
      
            if (!bg->insertIBPS(ibps,k,onlyIfExtending.getRef(),false)) {
                ORSA_DEBUG("problem, body: [%s]",(*bl_it).get()->getName().c_str());
            }
      
            ++bl_it;
        }
    }
  
    // ORSA_DEBUG("--MARK--");
  
    // frame_out += timestep_done;
    // frame_out.SetTime(frame_in + timestep_done);
    //
    /* 
       frame_out.SetTime(frame_in);
       frame_out += timestep_done;
    */
  
    // Predict new B values to use at the start of the next sequence. The predicted
    // values from the last call are saved as E. The correction, BD, between the
    // actual and predicted values of B is applied in advance as a correction.
    //
    double q1,q2,q3,q4,q5,q6,q7;
    //
    q1 = next_timestep.get_d() / timestep.get_d();
    q2 = q1 * q1;
    q3 = q1 * q2;
    q4 = q2 * q2;
    q5 = q2 * q3;
    q6 = q3 * q3;
    q7 = q3 * q4;
  
    /* 
       for(k=0;k<nv;++k) {
       if (interaction->IsSkippingJPLPlanets() && frame_in[k/3].JPLPlanet() != NONE) continue;
    */
    //
    {
        BodyGroup::BodyList::const_iterator bl_it = bl.begin();
        while (bl_it != bl.end()) {
      
            /* 
               if (!((*bl_it)->getInitialConditions().translational->dynamic())) {
               // ORSA_DEBUG("skipping body [%s]",(*bl_it)->getName().c_str());
               ++bl_it;
               continue;
               }
            */
            
            // ORSA_DEBUG("body [%s]",(*bl_it)->getName().c_str());
            
            const orsa::Body * k = (*bl_it).get();
      
            if ((*bl_it)->getInitialConditions().translational.get()) {
                if ((*bl_it)->getInitialConditions().translational->dynamic()) {
	  
                    orsa::Vector local_s[7];
	  
                    local_s[0] = b[0][k] - e[0][k];
                    local_s[1] = b[1][k] - e[1][k];
                    local_s[2] = b[2][k] - e[2][k];
                    local_s[3] = b[3][k] - e[3][k];
                    local_s[4] = b[4][k] - e[4][k];
                    local_s[5] = b[5][k] - e[5][k];
                    local_s[6] = b[6][k] - e[6][k];
	  
                    // Estimate B values for the next sequence
	  
                    e[0][k] = q1*(b[6][k]* 7.0 + b[5][k]* 6.0 + b[4][k]* 5.0 + b[3][k]* 4.0 + b[2][k]* 3.0 + b[1][k]*2.0 + b[0][k]);
                    e[1][k] = q2*(b[6][k]*21.0 + b[5][k]*15.0 + b[4][k]*10.0 + b[3][k]* 6.0 + b[2][k]* 3.0 + b[1][k]);
                    e[2][k] = q3*(b[6][k]*35.0 + b[5][k]*20.0 + b[4][k]*10.0 + b[3][k]* 4.0 + b[2][k]);
                    e[3][k] = q4*(b[6][k]*35.0 + b[5][k]*15.0 + b[4][k]* 5.0 + b[3][k]);
                    e[4][k] = q5*(b[6][k]*21.0 + b[5][k]* 6.0 + b[4][k]);
                    e[5][k] = q6*(b[6][k]* 7.0 + b[5][k]);
                    e[6][k] = q7* b[6][k];
	  
                    b[0][k] = e[0][k] + local_s[0];
                    b[1][k] = e[1][k] + local_s[1];
                    b[2][k] = e[2][k] + local_s[2];
                    b[3][k] = e[3][k] + local_s[3];
                    b[4][k] = e[4][k] + local_s[4];
                    b[5][k] = e[5][k] + local_s[5];
                    b[6][k] = e[6][k] + local_s[6];
                } 
            }
      
            if ((*bl_it)->getInitialConditions().rotational.get()) {
                if ((*bl_it)->getInitialConditions().rotational->dynamic()) {
	  
                    orsa::Quaternion local_s[7];
	  
                    // 
                    /* 
                       local_s[0] = bPhi[0][k] - ePhi[0][k];
                       local_s[1] = bPhi[1][k] - ePhi[1][k];
                       local_s[2] = bPhi[2][k] - ePhi[2][k];
                       local_s[3] = bPhi[3][k] - ePhi[3][k];
                       local_s[4] = bPhi[4][k] - ePhi[4][k];
                       local_s[5] = bPhi[5][k] - ePhi[5][k];
                       local_s[6] = bPhi[6][k] - ePhi[6][k];
	     
                       // Estimate B values for the next sequence
	     
                       ePhi[0][k] = q1*(bPhi[6][k]* 7.0 + bPhi[5][k]* 6.0 + bPhi[4][k]* 5.0 + bPhi[3][k]* 4.0 + bPhi[2][k]* 3.0 + bPhi[1][k]*2.0 + bPhi[0][k]);
                       ePhi[1][k] = q2*(bPhi[6][k]*21.0 + bPhi[5][k]*15.0 + bPhi[4][k]*10.0 + bPhi[3][k]* 6.0 + bPhi[2][k]* 3.0 + bPhi[1][k]);
                       ePhi[2][k] = q3*(bPhi[6][k]*35.0 + bPhi[5][k]*20.0 + bPhi[4][k]*10.0 + bPhi[3][k]* 4.0 + bPhi[2][k]);
                       ePhi[3][k] = q4*(bPhi[6][k]*35.0 + bPhi[5][k]*15.0 + bPhi[4][k]* 5.0 + bPhi[3][k]);
                       ePhi[4][k] = q5*(bPhi[6][k]*21.0 + bPhi[5][k]* 6.0 + bPhi[4][k]);
                       ePhi[5][k] = q6*(bPhi[6][k]* 7.0 + bPhi[5][k]);
                       ePhi[6][k] = q7* bPhi[6][k];
	     
                       bPhi[0][k] = ePhi[0][k] + local_s[0];
                       bPhi[1][k] = ePhi[1][k] + local_s[1];
                       bPhi[2][k] = ePhi[2][k] + local_s[2];
                       bPhi[3][k] = ePhi[3][k] + local_s[3];
                       bPhi[4][k] = ePhi[4][k] + local_s[4];
                       bPhi[5][k] = ePhi[5][k] + local_s[5];
                       bPhi[6][k] = ePhi[6][k] + local_s[6];
                    */
	  
                    local_s[0] = bQ[0][k] - eQ[0][k];
                    local_s[1] = bQ[1][k] - eQ[1][k];
                    local_s[2] = bQ[2][k] - eQ[2][k];
                    local_s[3] = bQ[3][k] - eQ[3][k];
                    local_s[4] = bQ[4][k] - eQ[4][k];
                    local_s[5] = bQ[5][k] - eQ[5][k];
                    local_s[6] = bQ[6][k] - eQ[6][k];
	  
                    // Estimate B values for the next sequence
	  
                    eQ[0][k] = q1*(bQ[6][k]* 7.0 + bQ[5][k]* 6.0 + bQ[4][k]* 5.0 + bQ[3][k]* 4.0 + bQ[2][k]* 3.0 + bQ[1][k]*2.0 + bQ[0][k]);
                    eQ[1][k] = q2*(bQ[6][k]*21.0 + bQ[5][k]*15.0 + bQ[4][k]*10.0 + bQ[3][k]* 6.0 + bQ[2][k]* 3.0 + bQ[1][k]);
                    eQ[2][k] = q3*(bQ[6][k]*35.0 + bQ[5][k]*20.0 + bQ[4][k]*10.0 + bQ[3][k]* 4.0 + bQ[2][k]);
                    eQ[3][k] = q4*(bQ[6][k]*35.0 + bQ[5][k]*15.0 + bQ[4][k]* 5.0 + bQ[3][k]);
                    eQ[4][k] = q5*(bQ[6][k]*21.0 + bQ[5][k]* 6.0 + bQ[4][k]);
                    eQ[5][k] = q6*(bQ[6][k]* 7.0 + bQ[5][k]);
                    eQ[6][k] = q7* bQ[6][k];
	  
                    bQ[0][k] = eQ[0][k] + local_s[0];
                    bQ[1][k] = eQ[1][k] + local_s[1];
                    bQ[2][k] = eQ[2][k] + local_s[2];
                    bQ[3][k] = eQ[3][k] + local_s[3];
                    bQ[4][k] = eQ[4][k] + local_s[4];
                    bQ[5][k] = eQ[5][k] + local_s[5];
                    bQ[6][k] = eQ[6][k] + local_s[6];
	  
                }
            }
      
            ++bl_it;
        }
    }
  
    // cerr << "-> out of Radau15::Step()..." << endl;
  
    // ORSA_DEBUG("done.");
  
    return true;
}

void IntegratorRadau::_init() {
  
    progressiveCleaningSteps = 32;
  
    _accuracy.set(1.0e-8);
  
    // for h and xc should use higher accuracy, since the double class has arbitrary precision!
  
    h[0] = 0.0;
    h[1] = 0.05626256053692215;
    h[2] = 0.18024069173689236;
    h[3] = 0.35262471711316964;
    h[4] = 0.54715362633055538;
    h[5] = 0.73421017721541053;
    h[6] = 0.88532094683909577;
    h[7] = 0.97752061356128750;
  
    xc[0] = 0.5;
    xc[1] = 0.16666666666666667;
    xc[2] = 0.08333333333333333;
    xc[3] = 0.05;
    xc[4] = 0.03333333333333333;
    xc[5] = 0.02380952380952381;
    xc[6] = 0.01785714285714286;
    xc[7] = 0.01388888888888889;
  
    vc[0] = 0.5;
    vc[1] = 0.3333333333333333;
    vc[2] = 0.25;
    vc[3] = 0.2;
    vc[4] = 0.1666666666666667;
    vc[5] = 0.1428571428571429;
    vc[6] = 0.125;
  
    // r.resize(28);
    //
    int j,k,l;
    l=0;
    for (j=1;j<8;++j) {
        for(k=0;k<j;++k) {
            r[l] = 1.0 / (h[j] - h[k]);
            ++l;
        }
    }
  
    /* 
       for(k=0;k<28;++k) {
       printf("r[%02i]= %e\n",k,r[k]);
       }
    */
  
    // c.resize(21);
    // d.resize(21);
    //
    c[0] = -h[1];
    d[0] =  h[1];
    l=0;
    for (j=2;j<7;++j) {
        ++l;
        c[l] = -h[j] * c[l-j+1];
        d[l] =  h[1] * d[l-j+1];
        for(k=2;k<j;++k) {
            ++l;
            c[l] = c[l-j] - h[j] * c[l-j+1];
            d[l] = d[l-j] + h[k] * d[l-j+1];
        }
        ++l;
        c[l] = c[l-j] - h[j];
        d[l] = d[l-j] + h[j]; 
    }
  
    /*
      for(k=0;k<21;++k) {
      printf("c[%02i]= %e    d[%02i]= %e\n",k,c[k],k,d[k]);
      }
    */
  
    // nv    = 0;
  
    // niter = 6;
  
    // s.resize(9);
  
    // _lastCallRejected = false;
  
    size = 0;
}

void IntegratorRadau::_body_mass_or_number_changed(orsa::BodyGroup  * bg,
                                                   const orsa::Time & t) {
  
    if (bg->size() > size) {
    
        // translational
    
        g.resize(7);
        b.resize(7);
        e.resize(7);
        //
        for (unsigned int l=0;l<7;++l) {
            g[l].clear();
            b[l].clear();
            e[l].clear();
        }
        x.clear();
        v.clear();
        a.clear();
        //
        x1.clear();
        v1.clear();
        a1.clear();
        //
        acc.clear();
    
        // rotational
    
        /* 
           gPhi.resize(7);
           bPhi.resize(7);
           ePhi.resize(7);
           //
           for (unsigned int l=0;l<7;++l) {
           gPhi[l].clear();
           bPhi[l].clear();
           ePhi[l].clear();
           }
           phi.clear();
           phiDot.clear();
           phiDotDot.clear();
           //
           phi1.clear();
           phi1Dot.clear();
           phi1DotDot.clear();
        */
        //
        gQ.resize(7);
        bQ.resize(7);
        eQ.resize(7);
        //
        for (unsigned int l=0;l<7;++l) {
            gQ[l].clear();
            bQ[l].clear();
            eQ[l].clear();
        }
        Q.clear();
        QDot.clear();
        QDotDot.clear();
        //
        Q1.clear();
        Q1Dot.clear();
        Q1DotDot.clear();
    
        torque.clear();
    }
  
    {   
        mass.clear();
        double m;
        const BodyGroup::BodyList & bl = bg->getBodyList();
        BodyGroup::BodyList::const_iterator bl_it = bl.begin();
        while (bl_it != bl.end()) {
            if (!(*bl_it)->alive(t)) {
                ++bl_it;
                continue;
            }
            if (!bg->getInterpolatedMass(m,(*bl_it).get(),t)) {
                ORSA_DEBUG("problems...");
            }	
            mass[(*bl_it).get()] = m;
            ++bl_it;
        }
    }
  
    size = bg->size();
}

void IntegratorRadau::reset() {
    // ORSA_DEBUG("RESET called");
    // this should be sufficient
    size = 0;
}

