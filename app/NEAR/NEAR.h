#ifndef _NEAR_H_
#define _NEAR_H_

#include <orsa/body.h>
#include <orsa/bodygroup.h>
#include <orsa/print.h>
#include <orsa/integrator.h>
#include <orsa/integrator_radau.h>
#include <orsa/multimin.h>
#include <orsa/util.h>

// for a constant mass body
class SolarRadiationPressure : public orsa::Propulsion {
 public:
  SolarRadiationPressure(const double     & bodyMass_in,
			 const double     & B_in,
			 orsa::BodyGroup  * bg_in,
			 const orsa::Body * sun_in,
			 const orsa::Body * body_in) :
    orsa::Propulsion(),
    bodyMass(bodyMass_in),
    B(B_in),
    bg(bg_in),
    sun(sun_in),
    body(body_in) {
  }
 protected:
  const double bodyMass;
  const double B;
  osg::ref_ptr<orsa::BodyGroup> bg;
  osg::ref_ptr<const orsa::Body> sun;
  osg::ref_ptr<const orsa::Body> body;
 public:
  orsa::Vector getThrust(const orsa::Time & t) const {
    // compute it once only, use forever...
    // if (!thrust.isSet()) {
    {
      orsa::Vector rSun;
      if (!bg->getInterpolatedPosition(rSun,sun.get(),t)) {
	ORSA_DEBUG("problems...");
      }	
      // orsa::print(rSun);
      orsa::Vector rBody;
      if (!bg->getInterpolatedPosition(rBody,body.get(),t)) {
	ORSA_DEBUG("problems...");
      }	
      // orsa::print(rBody);
      // Scheeres (1999) formalism
      const orsa::Vector   s2b = (rBody-rSun);
      // orsa::print(s2b);
      const orsa::Vector u_s2b = s2b.normalized();
      const double G1 = orsa::FromUnits(orsa::FromUnits(orsa::FromUnits(orsa::FromUnits(1.0e8,orsa::Unit::KG),orsa::Unit::KM,3),orsa::Unit::SECOND,-2),orsa::Unit::METER,-2);
      
      // now B is an external parameter
      // const double B  = orsa::FromUnits(487.0,orsa::Unit::KG)/orsa::FromUnits(10.0,orsa::Unit::METER,2); // Spacecraft mass to projected area ratio
      
      const double R  = s2b.length();
      const double acc = G1/(B*R*R);
      // ORSA_DEBUG("acc: %g",acc);
      const orsa::Vector thrust = (bodyMass*acc*u_s2b);
      return thrust;
    }
    // return thrust.getRef();
  }
 protected:
  // mutable orsa::Cache<orsa::Vector> thrust;
 public:
  bool nextEventTime(orsa::Time      &,
		     const mpz_class &) const {
    return false; // always ON
  }
};

orsa::BodyGroup * run();

bool processGravityFile(unsigned int & order,
			double       & mu, /* G*m */
			double       & R0,
			std::vector< std::vector<double> > & norm_C, 
			std::vector< std::vector<double> > & norm_S,
			const std::string & fileName);

class CustomIntegrator : public orsa::IntegratorRadau {
 public:
  CustomIntegrator(const orsa::Time & t0_in, const bool _verbose=true) :
    orsa::IntegratorRadau(),
    t0(t0_in),
    verbose(_verbose) { }
 protected:
  const orsa::Time t0;
  const bool verbose;
  mutable orsa::Cache<orsa::Time> lastOutput;
 public:
  void singleStepDone(orsa::BodyGroup  * bg,
		      const orsa::Time & call_t,
		      const orsa::Time & call_dt,
		      orsa::Time       & /* next_dt */ ) const {
    
    if (!verbose) return;
    
    const orsa::Time t = call_t+call_dt;
    
    const orsa::Time minOutputInterval = orsa::Time(0,0,1,0,0);
    bool print=false;
    if (lastOutput.isSet()) {
      if (t>lastOutput.getRef()) {
	if ((t-lastOutput.getRef())>minOutputInterval) {
	  print=true;
	}
      } else {
	if ((lastOutput.getRef()-t)>minOutputInterval) {
	  print=true;
	}
      }
    } else {
      print=true;
    }
    
    if (!print) return;
    
    lastOutput=t;
    
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
    
    // global dr
    const orsa::Vector dr = near_position - clone_position;
    //
    const orsa::Vector u_r = near_position.normalized();
    const orsa::Vector u_v = near_velocity.normalized();
    const orsa::Vector u_t = orsa::externalProduct(u_r,u_v).normalized();
    
    // rotate, for orbital elements
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
    
    ORSA_DEBUG("STEPDONE: %f %e %e %e %e %e %e %e %e %e %e %e",
	       (t-t0).get_d(),
	       near_position.length(),
	       clone_position.length(),
	       dr*u_r,
	       dr*u_v,
	       dr*u_t,
	       near_orbit.a,
	       near_orbit.e,
	       near_orbit.i*orsa::radToDeg(),
	       clone_orbit.a,
	       clone_orbit.e,
	       clone_orbit.i*orsa::radToDeg());
  }

};

// initial conditions multimin

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
    
    osg::ref_ptr<CustomIntegrator> radau = new CustomIntegrator(t0.getRef(),false);
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
  multimin->run_nmsimplex(256,1.0e-3);
  
  r.set(par->get("rx"),
	par->get("ry"),
	par->get("rz"));
  
  v.set(par->get("vx"),
	par->get("vy"),
	par->get("vz"));
  
  bg->clearIntegration();
  
  return true;
}

// eros mass multimin

class ErosMassMultimin : public orsa::Multimin {
 public:
  ErosMassMultimin() : 
    orsa::Multimin() { }
 public:
  orsa::Cache<orsa::Time>        t0;
  orsa::Cache<orsa::Time>        duration;
  orsa::Cache<double>            accuracy;
  osg::ref_ptr<const orsa::Body> sun;
  osg::ref_ptr<orsa::Body>       eros;
  osg::ref_ptr<const orsa::Body> near;
  osg::ref_ptr<orsa::Body>       clone;
  osg::ref_ptr<orsa::BodyGroup>  bg;
 public:
  double fun(const orsa::MultiminParameters * par) const {
    
    bg->clearIntegration();
    
    const double erosMass = par->get("erosMass");
    
    const double solarRadiationPressure_B = par->get("solarRadiationPressure_B");
    
    orsa::IBPS ibps = eros->getInitialConditions();
    osg::ref_ptr<orsa::ConstantInertialBodyProperty> inertial = new
      orsa::ConstantInertialBodyProperty(erosMass,
					 ibps.inertial->originalShape(),
					 ibps.inertial->centerOfMass(),
					 ibps.inertial->shapeToLocal(),
					 ibps.inertial->localToShape(),
					 ibps.inertial->inertiaMatrix(),
					 ibps.inertial->paulMoment());
    ibps.inertial = inertial.get();
    eros->setInitialConditions(ibps);

    double cloneMass;
    bg->getInterpolatedMass(cloneMass,clone.get(),t0.getRef());
    clone->propulsion = new SolarRadiationPressure(cloneMass,
						   solarRadiationPressure_B,
						   bg,
						   sun.get(),
						   clone.get());
    
    osg::ref_ptr<CustomIntegrator> radau = new CustomIntegrator(t0.getRef(),false);
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

inline bool betterErosMass(double           & mass,
			   double           & solarRadiationPressure_B,
			   const orsa::Time & t0,
			   const orsa::Time & duration,
			   const double     & accuracy, 
			   const orsa::Body * sun,  
			   orsa::Body       * eros,
			   const orsa::Body * near,
			   orsa::Body       * clone,
			   orsa::BodyGroup  * bg) {
  
  // initial values
  bg->getInterpolatedMass(mass,eros,t0);
  solarRadiationPressure_B = 
    orsa::FromUnits(487.0,orsa::Unit::KG)/orsa::FromUnits(10.0,orsa::Unit::METER,2); // Spacecraft mass to projected area ratio
  
  // mass uncertainty
  const double dm = 0.01*mass;
  const double dB = 0.25*solarRadiationPressure_B;
  
  osg::ref_ptr<orsa::MultiminParameters> par = new orsa::MultiminParameters;
  //
  par->insert("erosMass",mass,dm);
  par->insert("solarRadiationPressure_B",solarRadiationPressure_B,dB);
  
  osg::ref_ptr<ErosMassMultimin> multimin = new ErosMassMultimin;
  //
  multimin->t0       = t0;
  multimin->duration = duration;
  multimin->accuracy = accuracy;
  multimin->sun      = sun;
  multimin->eros     = eros;
  multimin->near     = near;
  multimin->clone    = clone;
  multimin->bg       = bg;
  //
  multimin->setMultiminParameters(par.get());
  //
  multimin->run_nmsimplex(64,1.0e-3);
  
  mass = par->get("erosMass");
  solarRadiationPressure_B = par->get("solarRadiationPressure_B");
  
  bg->clearIntegration();
  
  return true;
}

#endif // _NEAR_H_
