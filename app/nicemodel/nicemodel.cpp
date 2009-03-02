#include "nicemodel.h"

#include <tbb/task_scheduler_init.h>

using namespace std;
using namespace orsa;

static void addBody(orsa::BodyGroup   * bg, 
		    const std::string & planet_name,
		    const orsa::Double  planet_a,
		    const orsa::Double  planet_e,
		    const orsa::Double  planet_i, 
		    const orsa::Double  planet_mass,
		    const orsa::Double  sun_mass,
		    const orsa::Vector  rSun,
		    const orsa::Vector  vSun,
		    const orsa::Time  & tStart,
		    orsa::RNG         * rnd) {
  osg::ref_ptr<orsa::Body> b = new orsa::Body;
  b->setName(planet_name);
  // b->setMass(planet_mass);
  //
  orsa::Orbit orbit;
  orbit.mu = orsa::Unit::instance()->getG() * (sun_mass + planet_mass);
  orbit.a  = planet_a;
  orbit.e  = planet_e;
  orbit.i  = planet_i;
  orbit.omega_node       = orsa::twopi() * rnd->gsl_rng_uniform();
  orbit.omega_pericenter = orsa::twopi() * rnd->gsl_rng_uniform();
  orbit.M                = orsa::twopi() * rnd->gsl_rng_uniform();
  // orbit.omega_node = orbit.omega_pericenter = orbit.M = orsa::zero();
  //
  orsa::Vector r, v;
  orbit.relativePosVel(r,v);
  r += rSun;
  v += vSun;
  //
  orsa::IBPS ibps;
  ibps.time = tStart;
  ibps.inertial = new ConstantMassBodyProperty(planet_mass);
  ibps.translational = new orsa::DynamicTranslationalBodyProperty;
  ibps.translational->setPosition(r);
  ibps.translational->setVelocity(v);
  b->setInitialConditions(ibps);
  bg->addBody(b.get());
  // ORSA_DEBUG("added body [%s]",planet_name.c_str());
}

int main(int argc, char ** argv) {
  
  // TBB
  tbb::task_scheduler_init init;
  
  if (argc != 2) {
    cerr << "Usage: " << argv[0] << " <randomSeed>" << endl;
    exit(0);
  }
  
  orsa::Debug::instance()->initTimer();
  
  ORSA_DEBUG("process ID: %i",getpid());
  
  const int randomSeed = atoi(argv[1]);
  
  ORSA_DEBUG("randomSeed: %i",randomSeed);
  
  const orsa::Double sun_mass    = FromUnits(one(),Unit::MSUN);
  const orsa::Vector rSun        = orsa::Vector(0,0,0);
  const orsa::Vector vSun        = orsa::Vector(0,0,0);
  
  const orsa::Double jupiter_a    = FromUnits(5.45,orsa::Unit::AU);
  const orsa::Double jupiter_e    = 0.001; 
  const orsa::Double jupiter_i    = 0.001; 
  const orsa::Double jupiter_mass = FromUnits(1.90e27,orsa::Unit::KG);
  
  const orsa::Double saturn_a    = FromUnits(8.45,orsa::Unit::AU);
  const orsa::Double saturn_e    = 0.001;
  const orsa::Double saturn_i    = 0.001;
  const orsa::Double saturn_mass = FromUnits(5.68e26,orsa::Unit::KG);
  
  const orsa::Double uranus_a    = FromUnits(16,orsa::Unit::AU);
  const orsa::Double uranus_e    = 0.001;
  const orsa::Double uranus_i    = 0.001;
  const orsa::Double uranus_mass = FromUnits(8.68e25,orsa::Unit::KG);
  
  const orsa::Double neptune_a    = FromUnits(12,orsa::Unit::AU);
  const orsa::Double neptune_e    = 0.001;
  const orsa::Double neptune_i    = 0.001;
  const orsa::Double neptune_mass = FromUnits(1.02e26,orsa::Unit::KG);
  
  const unsigned int numParticles  = 4096;
  const orsa::Double totalDiskMass = 35*orsa::FromUnits(6e24,orsa::Unit::KG);
  const orsa::Double aMinDisk      = std::max(uranus_a,neptune_a);
  const orsa::Double aMaxDisk      = orsa::FromUnits(30,orsa::Unit::AU);
  
  const orsa::Time tStart = Time(                       0,0,0,0,0);
  const orsa::Time tStop  = Time(mpz_class("20000000000"),0,0,0,0);
  // const orsa::Time tStop  = Time(mpz_class("10000"),0,0,0,0);
  //
  const orsa::Time samplingPeriod = Time(20,0,0,0,0);
  
  osg::ref_ptr<QSIntegrator> qsIntegrator = new QSIntegrator(randomSeed);
  
  BodyGroup * bg = new BodyGroup;
  
  osg::ref_ptr<orsa::Body> b;
  
  osg::ref_ptr<orsa::RNG> rnd = new RNG(randomSeed);
  
  {
    b = new orsa::Body;
    b->setName("Sun");
    // b->setMass(sun_mass);
    orsa::IBPS ibps;
    ibps.time = tStart;
    ibps.inertial = new ConstantMassBodyProperty(sun_mass);
    ibps.translational = new orsa::DynamicTranslationalBodyProperty;
    ibps.translational->setPosition(rSun);
    ibps.translational->setVelocity(vSun);
    b->setInitialConditions(ibps);
    qsIntegrator->sun = b;    
    bg->addBody(b.get());
  }
  
  /* 
     {
     b = new orsa::Body;
     b->setName("Jupiter");
     b->setMass(jupiter_mass);
     //
     orsa::Orbit orbit;
     orbit.mu = orsa::Unit::instance()->getG() * (sun_mass + planet_mass);
     orbit.a  = planet_a;
     orbit.e  = planet_e;
     orbit.i  = planet_i;
     orbit.omega_node = orbit.omega_pericenter = orbit.M = orsa::zero();
     orsa::Vector r, v;
     orbit.relativePosVel(r,v);
     r += sun_r;
     v += sun_v;
     //
     orsa::IBPS ibps;
     ibps.time = tStart;
     ibps.translational = new orsa::DynamicTranslationalBodyProperty;
     ibps.translational->setPosition(r);
     ibps.translational->setVelocity(v);
     b->setInitialConditions(ibps);
     planet_r = r;
     planet_v = v;
     bg->addBody(b.get());
     }
  */
  
  addBody(bg, 
	  "Jupiter",
	  jupiter_a,
	  jupiter_e,
	  jupiter_i, 
	  jupiter_mass,
	  sun_mass,
	  rSun,
	  vSun,
	  tStart,
	  rnd.get());
  
  addBody(bg, 
	  "Saturn",
	  saturn_a,
	  saturn_e,
	  saturn_i, 
	  saturn_mass,
	  sun_mass,
	  rSun,
	  vSun,
	  tStart,
	  rnd.get());
  
  addBody(bg, 
	  "Uranus",
	  uranus_a,
	  uranus_e,
	  uranus_i, 
	  uranus_mass,
	  sun_mass,
	  rSun,
	  vSun,
	  tStart,
	  rnd.get());
  
  addBody(bg, 
	  "Neptune",
	  neptune_a,
	  neptune_e,
	  neptune_i, 
	  neptune_mass,
	  sun_mass,
	  rSun,
	  vSun,
	  tStart,
	  rnd.get());
  
  ORSA_DEBUG("to-be-fixed: mass surface density should fall linearly with distance");
  
  if (numParticles > 0) {
    const orsa::Double singleMass = totalDiskMass / numParticles;
    for (unsigned int id=0; id<numParticles; ++id) {
      b = new orsa::Body;
      char name[1024];
      snprintf(name,1024,"QS%04i",id);
      b->setName(name);
      orsa::Orbit orbit;
      orbit.mu = orsa::Unit::instance()->getG() * (sun_mass);
      orbit.a  = aMinDisk + (aMaxDisk-aMinDisk)*rnd->gsl_rng_uniform();
      orbit.e  = 0.001;
      orbit.i  = 0.001;
      orbit.omega_node       = orsa::twopi() * rnd->gsl_rng_uniform();
      orbit.omega_pericenter = orsa::twopi() * rnd->gsl_rng_uniform();
      orbit.M                = orsa::twopi() * rnd->gsl_rng_uniform();
      //
      orsa::Vector r, v;
      orbit.relativePosVel(r,v);
      r += rSun;
      v += vSun;
      //
      orsa::IBPS ibps;
      ibps.time = tStart;
      ibps.inertial = new ConstantMassBodyProperty(singleMass);
      ibps.translational = new orsa::DynamicTranslationalBodyProperty;
      ibps.translational->setPosition(r);
      ibps.translational->setVelocity(v);
      b->setInitialConditions(ibps);
      b->nonInteractingGroup = true;
      bg->addBody(b.get());
    }
  }
  
  // print initial state
  qsIntegrator->output(bg,tStart);
  
  qsIntegrator->integrate(bg,
			  tStart,
			  tStop,
			  samplingPeriod);
  
  return 0;
}

