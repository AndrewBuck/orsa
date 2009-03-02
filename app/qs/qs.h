#ifndef _QS_H_
#define _QS_H_

#include <orsaTBB/malloc.h>

#include <orsa/body.h>
#include <orsa/bodygroup.h>
#include <orsa/integrator_leapfrog.h>
// #include <orsa/integrator_radau.h>
#include <orsa/orbit.h>
#include <orsa/unit.h>
#include <orsa/util.h>

#include <QHash>

class QSIntegrator : public orsa::IntegratorLeapFrog {
  
 public:
  QSIntegrator(const int rs) : orsa::IntegratorLeapFrog(), randomSeed(rs) {
    rnd = new orsa::RNG(randomSeed);
    outputInitialized = false;
  }
    
 protected:
  ~QSIntegrator() { }
  
 protected:
  const int randomSeed;
  
 public:
  mutable osg::ref_ptr<orsa::Body> sun;
  mutable osg::ref_ptr<orsa::Body> planet;
  
 private:
  mutable bool outputInitialized;
  
 private:
  mutable orsa::Cache<orsa::Time> lastPerturbTime;
  
 public:
  typedef QHash <
    const orsa::Body * , 
    std::FILE * > BodyOutputFile;
 public:
  mutable BodyOutputFile bodyOutputFile;
  
 private:
  std::string outputFileName(const std::string & bodyName, const int randomSeed) const {
    char filename[1024];
    snprintf(filename,1024,"%i.%s.dat",randomSeed,bodyName.c_str());
    return filename;
  }
  
 protected:
  void output(orsa::BodyGroup  * bg,
	      const orsa::Time & t) const {
    
    if (!outputInitialized) { 
      orsa::BodyGroup::BodyList::iterator b_it = bg->getBodyList().begin();
      while (b_it != bg->getBodyList().end()) {
	if ((*b_it) == sun.get()) {
	  ++b_it;
	  continue;
	}
	
	bodyOutputFile.insert((*b_it).get(),fopen(outputFileName((*b_it)->getName(),randomSeed).c_str(),"w"));
	
	// ORSA_DEBUG("%x",bodyOutputFile[(*b_it).get()]);
	
	++b_it;
      }
      outputInitialized = true;
    }
    
    orsa::Orbit orbitPlanet;
    orsa::Vector rp, vp;
    if (bg->getInterpolatedPosVel(rp,vp,planet.get(),t)) {
      if (!orbitPlanet.compute(planet.get(),sun.get(),bg,t)) {
	ORSA_DEBUG("problems...");
      }
    } else {
      ORSA_DEBUG("problems...");
    }
    
    orsa::Vector r, v;
    orsa::BodyGroup::BodyList::iterator b_it = bg->getBodyList().begin();
    while (b_it != bg->getBodyList().end()) {
      
      if ((*b_it) == sun.get()) {
	++b_it;
	continue;
      }
      
      if (!(*b_it)->alive(t)) {
	++b_it;
	continue;
      }
      
      if (bg->getInterpolatedPosVel(r,v,(*b_it).get(),t)) {
	
	orsa::Orbit orbit;
	if (!orbit.compute((*b_it).get(),sun.get(),bg,t)) {
	  ORSA_DEBUG("problems");
	}
	
	orsa::Double criticalArgument = 
	  (      orbit.omega_node +       orbit.omega_pericenter +       orbit.M) - 
	  (orbitPlanet.omega_node + orbitPlanet.omega_pericenter + orbitPlanet.M);
	criticalArgument = orsa::fmod(6*orsa::twopi() + criticalArgument, orsa::twopi());
	if (criticalArgument > orsa::pi()) criticalArgument -= orsa::twopi(); 
	
	gmp_fprintf(bodyOutputFile[(*b_it).get()],
		    "%14.3Ff %12.9Fe %12.9Fe %10.6Ff %10.6Ff %10.6Ff %10.6Ff %+12.6Ff %12.6Fe %12.6Fe\n",
		    FromUnits(t.asDouble(),orsa::Unit::YEAR,-1).get_mpf_t(),
		    orsa::FromUnits(orbit.a,orsa::Unit::AU,-1).get_mpf_t(),
		    orbit.e.get_mpf_t(),
		    orsa::Double(orsa::radToDeg()*orbit.i).get_mpf_t(),
		    orsa::Double(orsa::radToDeg()*orbit.omega_node).get_mpf_t(),  
		    orsa::Double(orsa::radToDeg()*orbit.omega_pericenter).get_mpf_t(),
		    orsa::Double(orsa::radToDeg()*orbit.M).get_mpf_t(),
		    orsa::Double(orsa::radToDeg()*criticalArgument).get_mpf_t(),
		    orsa::FromUnits(orbit.a-orbitPlanet.a,orsa::Unit::AU,-1).get_mpf_t(),
		    orsa::FromUnits((r-rp).length(),orsa::Unit::AU,-1).get_mpf_t());
	// fflush(bodyOutputFile[(*b_it).get()]);
	
      } else {
	ORSA_DEBUG("problems");
	return;
      }
      
      ++b_it;
    }
    
  }
  
 protected:
  
  void cleanup(orsa::BodyGroup  * bg,
	       const orsa::Time & t) const {
    
    // remove all data other not at time = call_t+call_dt
    {
      // ORSA_DEBUG("remember, you're cleaning the BodyGroup intervals...");
      
      // const orsa::Time t = call_t + call_dt;
      
      orsa::BodyGroup::BodyList::iterator _b_it = bg->getBodyList().begin();
      while (_b_it != bg->getBodyList().end()) { 
	if (!((*_b_it)->getInitialConditions().translational->dynamic())) { 
	  ++_b_it;
	  continue;
	}
	orsa::BodyGroup::BodyInterval * _b_interval = bg->getBodyInterval((*_b_it).get());
	orsa::BodyGroup::BodyInterval::DataType & _b_interval_data = _b_interval->getData();
	orsa::BodyGroup::BodyInterval::DataType::iterator _b_interval_data_it = _b_interval_data.begin();
	//
	while (_b_interval_data_it != _b_interval_data.end()) {
	  if ((*_b_interval_data_it).time.getRef() != t) {
	    _b_interval_data_it = _b_interval_data.erase(_b_interval_data_it);
	  } else {
	    ++_b_interval_data_it;
	  }
	}     
	// IMPORTANT!
	_b_interval->update();
	
	++_b_it;
      }
    }
    
    // "ejected" bodies?
    {
      
      orsa::Orbit orbitPlanet;
      orsa::Vector rp, vp;
      if (bg->getInterpolatedPosVel(rp,vp,planet.get(),t)) {
	if (!orbitPlanet.compute(planet.get(),sun.get(),bg,t)) {
	  ORSA_DEBUG("problems...");
	}
      } else {
	ORSA_DEBUG("problems...");
      }
      
      // orsa::Vector r, v;
      orsa::BodyGroup::BodyList::iterator b_it = bg->getBodyList().begin();
      while (b_it != bg->getBodyList().end()) {
	
	if ((*b_it) == sun.get()) {
	  ++b_it;
	  continue;
	}
	
	if ((*b_it) == planet.get()) {
	  ++b_it;
	  continue;
	}
	
	if (!(*b_it)->alive(t)) {
	  ++b_it;
	  continue;
	}
	
	// if (bg->getInterpolatedPosVel(r,v,(*b_it).get(),t)) {
	
	orsa::Orbit orbit;
	if (!orbit.compute((*b_it).get(),sun.get(),bg,t)) {
	  ORSA_DEBUG("problems");
	} else {
	
	  // 20% difference max in semi-major axis, or big eccentricity
	  if ( (orsa::fabs((orbit.a-orbitPlanet.a)/orbitPlanet.a) > 0.2) ||
	       (orsa::fabs(orbit.e) > 0.8) ) {
	    ORSA_DEBUG("ejecting body [%s]",
		       (*b_it)->getName().c_str());
	    (*b_it)->deathTime = t;
	    // flush relative output file
	    fflush(bodyOutputFile[(*b_it).get()]);
	  }
	}
	
	// } else {
	// ORSA_DEBUG("problems");
	// return;
	// }
	
	++b_it;
      }
      
    }
    
  }
  
 public:
  osg::ref_ptr<orsa::RNG> rnd;
  
 protected:
  mutable orsa::Cache<bool> resetNeeded;
  
 protected:    
  bool step(orsa::BodyGroup  * bg,
	    const orsa::Time & start,
	    const orsa::Time & timestep,
	    orsa::Time       & next_timestep) {
    if (resetNeeded.isSet()) {
      if (resetNeeded.getRef()) {
	// ORSA_DEBUG("calling reset()...");
	reset();
	resetNeeded = false;
      }
    }
    return orsa::IntegratorLeapFrog::step(bg,start,timestep,next_timestep);
  }
  
 protected:
  void singleStepDone(orsa::BodyGroup  * bg,
		      const orsa::Time & call_t,
		      const orsa::Time & call_dt,
		      orsa::Time       & ) const {
    
    const orsa::Time perturbPeriod(100,0,0,0,0);
    
    bool doPerturb=false;
    // comment all this block to disable perturbations
    if (lastPerturbTime.isSet()) {
      if ((call_t+call_dt-lastPerturbTime.getRef()) >= perturbPeriod) {
	doPerturb=true;
      }
    } else {
      doPerturb=true;
    }	
    
    if (doPerturb) {
      
      lastPerturbTime = call_t+call_dt;
      
      const orsa::Double minimumPerturberRelativeDistance = orsa::FromUnits(0.05,orsa::Unit::AU);
      const orsa::Double maximumPerturberRelativeDistance = orsa::FromUnits(2.00,orsa::Unit::AU);
      
      const orsa::Double maxPerturberEccentricity = 0.3;
      
      const orsa::Double    maxPerturberInclination = 30 * orsa::degToRad();
      const orsa::Double cosMaxPerturberInclination = orsa::cos(maxPerturberInclination);
      const orsa::Double sinMaxPerturberInclination = orsa::sin(maxPerturberInclination);
      
      orsa::Vector r_sun, v_sun;
      if (!bg->getInterpolatedPosVel(r_sun, v_sun, sun.get(), call_t+call_dt)) {
	ORSA_DEBUG("problems...");
      }
      
      double x,y,z;
      orsa::Vector rb, vb;
      
      orsa::BodyGroup::BodyList::iterator b_it = bg->getBodyList().begin();
      while (b_it != bg->getBodyList().end()) {
	
	if ((*b_it) == sun.get()) {
	  ++b_it;
	  continue;
	}
	
	if (!(*b_it)->alive(call_t+call_dt)) {
	  ++b_it;
	  continue;
	}
	
	if (!bg->getInterpolatedPosVel(rb, vb, (*b_it).get(), call_t+call_dt)) {
	  ORSA_DEBUG("problems...");
	}
	
	// coarse inclination check: if body too high, the check will always fail, and we detect this here
	{
	  const orsa::Vector r1 = rb + orsa::Vector(0,0,+1) * maximumPerturberRelativeDistance;
	  const orsa::Vector r2 = rb + orsa::Vector(0,0,-1) * maximumPerturberRelativeDistance;
	  if (std::min(orsa::fabs(r1.getZ()/r1.length()),
		       orsa::fabs(r2.getZ()/r2.length())) > sinMaxPerturberInclination) {
	    // ORSA_DEBUG("rejecting at coarse inclination check...");
	    ++b_it;
	    continue;
	  }
	}
	
	const orsa::Double limit = orsa::FromUnits(30,orsa::Unit::AU);
	
	const orsa::Double strength = 
	  std::min(orsa::Double(1),
		   std::max(orsa::Double(-1),
			    orsa::Double(4.0 * (limit - (rb-r_sun).length()) / limit)));
	
	// ORSA_DEBUG("strength: %Ff",strength.get_mpf_t());
	
	const orsa::Vector u_offset = strength * (vb-v_sun).normalized();
	
	// const orsa::Double perturberMass = orsa::FromUnits(1.0e22,orsa::Unit::KG);
	//
	/* const orsa::Double perturberMass = 
	   std::min(orsa::FromUnits(1.0e22,orsa::Unit::KG),
	   orsa::Double(strength*orsa::FromUnits(1.0e18,orsa::Unit::KG)*pow(10,fabs(rnd->gsl_ran_laplace(0.4)))));
	*/
	//
	orsa::Double perturberMass;
	do {
	  // perturberMass = orsa::fabs(strength)*orsa::FromUnits(1.0e18,orsa::Unit::KG)*pow(10,fabs(rnd->gsl_ran_laplace(0.4)));
	  //
	  perturberMass = orsa::FromUnits(1.0e18,orsa::Unit::KG)*pow(10,fabs(rnd->gsl_ran_laplace(0.4)));
	  if (perturberMass < orsa::FromUnits(1.0e22,orsa::Unit::KG)) break;
	} while (1);
	
	// const orsa::Double perturberRelativeDistance = orsa::FromUnits(0.3,orsa::Unit::AU);
	//
	/* const orsa::Double perturberRelativeDistance = 
	   std::min(orsa::FromUnits(2.0,orsa::Unit::AU),
	   std::max(orsa::FromUnits(0.1,orsa::Unit::AU),
	   strength * orsa::FromUnits(0.1,orsa::Unit::AU) * pow(10,fabs(rnd->gsl_ran_laplace(0.3)))));
	*/
	//
	orsa::Double perturberRelativeDistance;
	do {
	  perturberRelativeDistance = orsa::FromUnits(0.05,orsa::Unit::AU) * pow(10,fabs(rnd->gsl_ran_laplace(0.3)));
	  //
	  /* ORSA_DEBUG("perturberRelativeDistance: %Fg [AU]",
	     FromUnits(perturberRelativeDistance,orsa::Unit::AU,-1).get_mpf_t());
	  */
	  if ( (perturberRelativeDistance < maximumPerturberRelativeDistance) && 
	       (perturberRelativeDistance > minimumPerturberRelativeDistance) ) break;
	} while (1);
	
	orsa::Vector u;
	orsa::Double perturberRelativeVelocity;
	bool good = false;
	for (unsigned int trial=0; trial<32; ++trial) {
	  
	  rnd->gsl_ran_dir_3d(&x,&y,&z);
	  u = (u_offset + orsa::Vector(x,y,z).normalized()).normalized();
	  
	  rnd->gsl_ran_dir_3d(&x,&y,&z);
	  const orsa::Vector tmp_w = orsa::Vector(x,y,z).normalized();
	  
	  // orthogonal to u
	  const orsa::Vector v = externalProduct(u,tmp_w).normalized();
	  
	  const orsa::Vector perturberPosition = 
	    rb + u * perturberRelativeDistance;
	  
	  // inclination test
	  if (externalProduct((perturberPosition-r_sun),v).normalized().getZ() < cosMaxPerturberInclination) {
	    // ORSA_DEBUG("rejecting at inclination check...");
	    continue;
	  }
	  
	  const orsa::Double perturberRadialDistance = 
	    (perturberPosition-r_sun).length();
	  
	  const orsa::Double perturberSemiMajorAxis = 
	    orsa::FromUnits((25+5*rnd->gsl_rng_uniform()),orsa::Unit::AU);
	  
	  // ORSA_DEBUG("perturber a: %20.3Fg [AU]",
	  // FromUnits(perturberSemiMajorAxis,orsa::Unit::AU,-1).get_mpf_t());
	  
	  if (2*perturberSemiMajorAxis < perturberRadialDistance) {
	    /* 
	       ORSA_DEBUG("square root check failed:   a: %8.3Ff   R: %8.3Ff   body: [%s]",
	       FromUnits(perturberSemiMajorAxis,orsa::Unit::AU,-1).get_mpf_t(),
	       FromUnits(perturberRadialDistance,orsa::Unit::AU,-1).get_mpf_t(),
	       (*b_it)->getName().c_str());
	    */
	    //
	    // ORSA_DEBUG("rejecting at square root check...");
	    continue;
	  }
	  
	  const orsa::Double v_mu = 
	    sqrt(orsa::Unit::instance()->getG() * FromUnits(orsa::one(),orsa::Unit::MSUN) * 
		 (orsa::two()/perturberRadialDistance - orsa::one()/perturberSemiMajorAxis));
	  
	  // const orsa::Vector v_mu_vec = v_mu * v;
	  
	  // rnd->gsl_ran_dir_2d(&x,&y);
	  // const orsa::Double v_mu_x = v_mu * x;
	  // const orsa::Double v_mu_y = v_mu * y;
	  
	  // perturberRelativeVelocity = orsa::fabs(v_mu_x*v.getX()+v_mu_y*v.getY()-(vb-v_sun)*v);
	  
	  {
	    const orsa::Double RV = (perturberPosition-r_sun)*v*v_mu;
	    const orsa::Double maxRVsquared = 
	      orsa::Unit::instance()->getG() * FromUnits(orsa::one(),orsa::Unit::MSUN) * perturberSemiMajorAxis *
	      (maxPerturberEccentricity*maxPerturberEccentricity - orsa::int_pow(1-perturberRadialDistance/perturberSemiMajorAxis,2));
	    if (RV*RV > maxRVsquared) {
	      // ORSA_DEBUG("rejecting at eccentricity check...");
	      continue;
	    }
	  }
	  
	  perturberRelativeVelocity = orsa::fabs(v_mu-(vb-v_sun)*v);
	  
	  // maximum 1/1000 velocity kick
	  if (perturberRelativeVelocity < (1000*2*orsa::Unit::instance()->getG()*perturberMass/(perturberRelativeDistance*(vb-v_sun).length()))) {
	    // ORSA_DEBUG("rejecting at kick threshold check...");
	    continue;
	  }
	  
	  /* 
	     ORSA_DEBUG("perturberRelativeVelocity: %20.3Ff [km/s]   vb: %20.3Ff [km/s]   v_mu: %20.3Ff [km/s]",
	     orsa::FromUnits(orsa::FromUnits(perturberRelativeVelocity,orsa::Unit::KM,-1),orsa::Unit::SECOND).get_mpf_t(),
	     orsa::FromUnits(orsa::FromUnits(              vb.length(),orsa::Unit::KM,-1),orsa::Unit::SECOND).get_mpf_t(),
	     orsa::FromUnits(orsa::FromUnits(                     v_mu,orsa::Unit::KM,-1),orsa::Unit::SECOND).get_mpf_t());
	  */
	  
	  good = true;
	  break;
	}
	
	if (!good) {
	  
	  // ORSA_DEBUG("could not find appropriate perturbation for body [%s]",(*b_it)->getName().c_str());
	  
	} else {
	  
	  // ORSA_DEBUG("-- kick accepted --");
	  
	  const orsa::Double kickMagnitude = 
	    2 * orsa::Unit::instance()->getG() * perturberMass / 
	    (perturberRelativeDistance*perturberRelativeVelocity);
	
	  // ORSA_DEBUG("******************* kickMagnitude: %Fg",kickMagnitude.get_mpf_t());
	
	  const orsa::Vector r = rb;
	  const orsa::Vector v = vb + u * kickMagnitude;
	
	  orsa::IBPS ibps;
	  //
	  if (bg->getInterpolatedIBPS(ibps,(*b_it).get(),call_t+call_dt)) {
	  
	    if (!ibps.translational->dynamic()) {
	      ORSA_DEBUG("problems...");
	    }
	  
	    ibps.time = call_t+call_dt;
	    //
	    if ((*b_it)->getInitialConditions().translational.get()) {
	      if ((*b_it)->getInitialConditions().translational->dynamic()) {
		ibps.translational->setPosition(r);
		ibps.translational->setVelocity(v);
	      }
	    }
	    //
	    if (!(bg->getBodyInterval((*b_it).get())->insert(ibps,true))) {
	      ORSA_DEBUG("problems with insert, body [%s]",
			 (*b_it)->getName().c_str());
	    }
	  }
      
	}
      
	++b_it;
      }
    }    
    
    // output + cleanup
    {
      const unsigned int sample = 256;
      static unsigned int iter = sample;
      if (iter == sample) {
	output(bg, call_t+call_dt);
	cleanup(bg, call_t+call_dt);
	iter = 0;
      }
      ++iter;
    }
    
  }
  
};

#endif // _QS_H_
