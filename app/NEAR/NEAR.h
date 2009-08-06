#ifndef _NEAR_H_
#define _NEAR_H_

#include <orsa/bodygroup.h>
#include <orsa/print.h>
#include <orsa/integrator.h>
#include <orsa/integrator_radau.h>
#include <orsa/multimin.h>
#include <orsa/util.h>

orsa::BodyGroup * run();

bool processGravityFile(unsigned int & order,
			double       & mu, /* G*m */
			double       & R0,
			std::vector< std::vector<double> > & norm_C, 
			std::vector< std::vector<double> > & norm_S,
			const std::string & fileName);

class CustomIntegrator : public orsa::IntegratorRadau {
 public:
  CustomIntegrator(const orsa::Time & t0_in) :
    orsa::IntegratorRadau(),
    t0(t0_in) { }
 protected:
  const orsa::Time t0;
 public:
  void singleStepDone(orsa::BodyGroup  * bg,
		      const orsa::Time & call_t,
		      const orsa::Time & call_dt,
		      orsa::Time       & /* next_dt */ ) const {
    
    const orsa::Time t = call_t+call_dt;
    
    const orsa::Body * eros  = bg->getBody("EROS");
    const orsa::Body * near  = bg->getBody("NEAR");
    const orsa::Body * clone = bg->getBody("CLONE");
    
    orsa::IBPS eros_ibps;
    bg->getInterpolatedIBPS(eros_ibps,eros,t);
    //
    const double eros_mu = orsa::Unit::G()*eros_ibps.inertial->mass();
    
    orsa::Vector eros_position;
    orsa::Vector eros_velocity;
    bg->getInterpolatedPosVel(eros_position,eros_velocity,eros,t);
    
    orsa::Vector near_position;
    orsa::Vector near_velocity;
    bg->getInterpolatedPosVel(near_position,near_velocity,near,t);
    
    orsa::Vector clone_position;
    orsa::Vector clone_velocity;
    bg->getInterpolatedPosVel(clone_position,clone_velocity,clone,t);
    
    near_position  -= eros_position;
    near_velocity  -= eros_velocity;
    
    clone_position -= eros_position;
    clone_velocity -= eros_velocity;
    
    const orsa::Matrix eros_g2l = orsa::globalToLocal(eros,bg,t);
    //
    near_position  = eros_g2l*near_position;
    near_velocity  = eros_g2l*near_velocity;
    //
    clone_position = eros_g2l*clone_position;
    clone_velocity = eros_g2l*clone_velocity;
    
    orsa::Orbit near_orbit;
    near_orbit.compute(near_position,near_velocity,eros_mu);

    orsa::Orbit clone_orbit;
    clone_orbit.compute(clone_position,clone_velocity,eros_mu);
    
    ORSA_DEBUG("STEPDONE: %f %e %e %e %e %e %e %e %e %e",
	       (t-t0).get_d(),
	       near_position.length(),
	       clone_position.length(),
	       (near_position-clone_position).length(),
	       near_orbit.a,
	       near_orbit.e,
	       near_orbit.i*orsa::radToDeg(),
	       clone_orbit.a,
	       clone_orbit.e,
	       clone_orbit.i*orsa::radToDeg());
  }

};

class InitialConditionsMultimin : public orsa::Multimin {
 public:
  InitialConditionsMultimin() : 
    orsa::Multimin() { }
 public:
  orsa::Cache<orsa::Time>        t0;
  orsa::Cache<orsa::Time>        duration;
  orsa::Cache<double>            accuracy;
  osg::ref_ptr<const orsa::Body> near;
  osg::ref_ptr<orsa::Body>       clone;
  osg::ref_ptr<orsa::BodyGroup>  bg;
 public:
  double fun(const orsa::MultiminParameters * par) const {
    
    bg->clearIntegration();
    
    const orsa::Vector r(par->get("rx"),
			 par->get("ry"),
			 par->get("rz"));
    
    const orsa::Vector v(par->get("vx"),
			 par->get("vy"),
			 par->get("vz")); 
    
    orsa::IBPS ibps = clone->getInitialConditions();
    ibps.translational->setPosition(r);
    ibps.translational->setVelocity(v);
    clone->setInitialConditions(ibps);
    
    osg::ref_ptr<CustomIntegrator> radau = new CustomIntegrator(t0.getRef());
    radau->_accuracy = accuracy;
    radau->keepOnlyLastStep = false;
    const orsa::Time samplingPeriod = orsa::Time(0,0,0,10,0);  
    // first call to output function
    // radau->singleStepDone(bg,t0,orsa::Time(0),orsa::Time(0));
    const bool goodIntegration = radau->integrate(bg,
						  t0.getRef(),
						  t0.getRef()+duration.getRef(),
						  samplingPeriod);
    
    if (!goodIntegration) {
      ORSA_DEBUG("problems...");
      exit(0);
    }
    
    double retVal = 0;
    const unsigned int npoints = 32;
    const orsa::Time dt = duration.getRef()/(npoints-1);
    for (unsigned int k=0; k<npoints; ++k) {
      
      const orsa::Time t = t0.getRef()+k*dt;
      
      orsa::Vector near_position;
      orsa::Vector near_velocity;
      bg->getInterpolatedPosVel(near_position,near_velocity,near.get(),t);      
    
      orsa::Vector clone_position;
      orsa::Vector clone_velocity;
      bg->getInterpolatedPosVel(clone_position,clone_velocity,clone.get(),t);
      
      const orsa::Vector dr = clone_position - near_position;
      
      ORSA_DEBUG("dr[%03i] = %g",k,dr.length());
      
      retVal += dr.lengthSquared();
    }
    
    ORSA_DEBUG("retVal: %g",retVal);
    
    bg->clearIntegration();
    
    return retVal;
  }
};

inline bool betterInitialConditions(orsa::Vector     & r,
				    orsa::Vector     & v,
				    const orsa::Time & t0,
				    const orsa::Time & duration,
				    const double     & accuracy,   
				    const orsa::Body * near,
				    orsa::Body       * clone,
				    orsa::BodyGroup  * bg) {
  
  // standard uncertainty
  const double dr = orsa::FromUnits(10,orsa::Unit::METER);
  const double dv = orsa::FromUnits(orsa::FromUnits(0.01,orsa::Unit::METER),orsa::Unit::SECOND,-1);
  
  // initial values, same as NEAR (not clone...)
  bg->getInterpolatedPosVel(r,v,near,t0);
  
  osg::ref_ptr<orsa::MultiminParameters> par = new orsa::MultiminParameters;
  //
  par->insert("rx",r.getX(),dr);
  par->insert("ry",r.getY(),dr);
  par->insert("rz",r.getZ(),dr);
  //
  par->insert("vx",v.getX(),dv);
  par->insert("vy",v.getY(),dv);
  par->insert("vz",v.getZ(),dv);
  
  osg::ref_ptr<InitialConditionsMultimin> multimin = new InitialConditionsMultimin;
  //
  multimin->t0       = t0;
  multimin->duration = duration;
  multimin->accuracy = accuracy;
  multimin->near     = near;
  multimin->clone    = clone;
  multimin->bg       = bg;
  //
  multimin->setMultiminParameters(par.get());
  //
  multimin->run_nmsimplex(1024,1.0e-6);
  
  r.set(par->get("rx"),
	par->get("ry"),
	par->get("rz"));
  
  v.set(par->get("vx"),
	par->get("vy"),
	par->get("vz"));
  
  bg->clearIntegration();
  
  return true;
}


// for a constant mass body
/* class SolarRadiationPressure : public orsa::Propulsion {
   public:
   SolarRadiationPressure(const double     & bodyMass_in,
   orsa::BodyGroup  * bg_in,
   const orsa::Body * sun_in,
   const orsa::Body * body_in) :
   orsa::Propulsion(),
   bodyMass(bodyMass_in),
   bg(bg_in),
   sun(sun_in),
   body(body_in) {
   }
   protected:
   const double bodyMass;
   osg::ref_ptr<orsa::BodyGroup> bg;
   osg::ref_ptr<const orsa::Body> sun;
   osg::ref_ptr<const orsa::Body> body;
   public:
   orsa::Vector getThrust(const orsa::Time & t) const {
   // compute it once only, use forever...
   if (!thrust.isSet()) {
   orsa::Vector rSun;
   if (!bg->getInterpolatedPosition(rSun,sun.get(),t)) {
   ORSA_DEBUG("problems...");
   }	
   orsa::print(rSun);
   orsa::Vector rBody;
   if (!bg->getInterpolatedPosition(rBody,body.get(),t)) {
   ORSA_DEBUG("problems...");
   }	
   orsa::print(rBody);
   // Scheeres formalism
   const orsa::Vector   s2b = (rBody-rSun);
   const orsa::Vector u_s2b = s2b.normalized();
   const double G1 = orsa::FromUnits(orsa::FromUnits(orsa::FromUnits(orsa::FromUnits(1e8,orsa::Unit::KG),orsa::Unit::KM,3),orsa::Unit::SECOND,-2),orsa::Unit::METER,-2);
   const double B  = orsa::FromUnits(487,orsa::Unit::KG)/orsa::FromUnits(10,orsa::Unit::METER,2); // Spacecraft mass to projected area ratio
   const double R  = s2b.length();
   const double acc = G1/(B*R*R);
   ORSA_DEBUG("acc: %g",acc);
   thrust = (bodyMass*acc*u_s2b);
   }
   return thrust.getRef();
   }
   protected:
   mutable orsa::Cache<orsa::Vector> thrust;
   public:
   bool nextEventTime(orsa::Time      &,
   const mpz_class &) const {
   return false; // always ON
   }
   };
*/

#endif // _NEAR_H_
