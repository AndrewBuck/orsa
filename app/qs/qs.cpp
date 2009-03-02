#include "qs.h"

using namespace std;
using namespace orsa;

int main(int argc, char ** argv) {
  
  if (argc != 2) {
    cerr << "Usage: " << argv[0] << " <randomSeed>" << endl;
    exit(0);
  }
  
  orsa::Debug::instance()->initTimer();
  
  ORSA_DEBUG("process ID: %i",getpid());
  
  const int randomSeed = atoi(argv[1]);
  
  ORSA_DEBUG("randomSeed: %i",randomSeed);
  
  const orsa::Double sun_mass    = FromUnits(one(),Unit::MSUN);
  const orsa::Vector sun_r       = orsa::Vector(0,0,0);
  const orsa::Vector sun_v       = orsa::Vector(0,0,0);
  
  const orsa::Double planet_a    = FromUnits(23,orsa::Unit::AU);
  const orsa::Double planet_e    = 0.00; // very small eccentricities should yeld more stable simulations
  const orsa::Double planet_mass = FromUnits(1.0e26,orsa::Unit::KG);
  
  // const orsa::Double discardFactor = 10.0;
  //
  /* const orsa::Double planetHillRadius = orsa::HillRadius(planet_a,
     planet_mass,
     sun_mass);
  */
  //
  // const orsa::Double discardDistance = discardFactor * planetHillRadius;
  
  /* ORSA_DEBUG("planet Hill radius = %Ff AU",
     FromUnits(planetHillRadius,orsa::Unit::AU,-1).get_mpf_t());
  */
  
  /* 
     ORSA_DEBUG("discard distance = %Ff AU",
     FromUnits(discardDistance,orsa::Unit::AU,-1).get_mpf_t());
  */
  
  const orsa::Time tStart = Time(                       0,0,0,0,0);
  const orsa::Time tStop  = Time(mpz_class("20000000000"),0,0,0,0);
  //
  const orsa::Time samplingPeriod = Time(25,0,0,0,0);
  
  osg::ref_ptr<QSIntegrator> qsIntegrator = new QSIntegrator(randomSeed);
  
  BodyGroup * bg = new BodyGroup;
  
  osg::ref_ptr<orsa::Body> b;
  
  orsa::Cache<orsa::Vector> planet_r, planet_v;
  
  {
    b = new orsa::Body;
    b->setName("Sun");
    b->setMass(sun_mass);
    orsa::IBPS ibps;
    ibps.time = tStart;
    ibps.translational = new orsa::DynamicTranslationalBodyProperty;
    ibps.translational->setPosition(sun_r);
    ibps.translational->setVelocity(sun_v);
    b->setInitialConditions(ibps);
    qsIntegrator->sun = b;    
    bg->addBody(b.get());
  }
  
  {
    b = new orsa::Body;
    b->setName("Planet");
    b->setMass(planet_mass);
    //
    orsa::Orbit orbit;
    orbit.mu = orsa::Unit::instance()->getG() * (sun_mass + planet_mass);
    orbit.a  = planet_a;
    orbit.e  = planet_e;
    orbit.i  = orsa::zero();
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
    qsIntegrator->planet = b;    
    planet_r = r;
    planet_v = v;
    bg->addBody(b.get());
  }
  
  /* 
     {
     b = new orsa::Body;
     b->setName("Moon001");
     b->setMass(orsa::zero());
     //
     orsa::Orbit orbit;
     orbit.mu = orsa::Unit::instance()->getG() * (planet_mass);
     orbit.a  = FromUnits(0.01,orsa::Unit::AU);
     orbit.e  = 0.04;
     orbit.i  = orsa::degToRad()*8.0;
     orbit.omega_node = orbit.omega_pericenter = orbit.M = orsa::zero();
     orsa::Vector r, v;
     orbit.relativePosVel(r,v);
     r += planet_r.getRef();
     v += planet_v.getRef();
     //
     orsa::IBPS ibps;
     ibps.time = tStart;
     ibps.translational = new orsa::DynamicTranslationalBodyProperty;
     ibps.translational->setPosition(r);
     ibps.translational->setVelocity(v);
     b->setInitialConditions(ibps);
     bg->addBody(b.get());
     }
  */
  
  osg::ref_ptr<orsa::RNG> rnd = new RNG(randomSeed);
  for (unsigned int id=0; id<16; ++id) {
    b = new orsa::Body;
    char name[1024];
    snprintf(name,1024,"QS%04i",id);
    b->setName(name);
    b->setMass(orsa::zero());
    //
    orsa::Orbit orbit;
    orbit.mu = orsa::Unit::instance()->getG() * (sun_mass);
    orbit.a  = planet_a + FromUnits(0.3*(2*rnd->gsl_rng_uniform()-1),orsa::Unit::AU);
    orbit.e  = 0.5*rnd->gsl_rng_uniform();
    orbit.i  = orsa::degToRad()*40.0*rnd->gsl_rng_uniform();
    orbit.omega_node       = orsa::twopi() * rnd->gsl_rng_uniform();
    orbit.omega_pericenter = orsa::twopi() * rnd->gsl_rng_uniform();
    orbit.M                = orsa::twopi() * rnd->gsl_rng_uniform();
    //
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
    bg->addBody(b.get());
  }
  
  qsIntegrator->integrate(bg,
			  tStart,
			  tStop,
			  samplingPeriod);
  
  return 0;
}

