#include "SurveyReview.h"

#include <orsaSPICE/spice.h>
#include <orsaSPICE/spiceBodyAttitudeCallback.h>
#include <orsaSPICE/spiceBodyPosVelCallback.h>

orsa::Body * SPICEBody (const std::string  & bodyName,
			const double & bodyMass) {
  
  orsa::Body * localBody = new orsa::Body;
  localBody->setName(bodyName);
  // localBody->setMass(bodyMass);  
  
  osg::ref_ptr<orsaSPICE::SpiceBodyPosVelCallback> sbpvc = new orsaSPICE::SpiceBodyPosVelCallback(bodyName);
  orsa::IBPS ibps;
  ibps.inertial = new orsa::ConstantMassBodyProperty(bodyMass);
  ibps.translational = sbpvc.get();
  localBody->setInitialConditions(ibps);
  
  return localBody;
}

// Taking definitons from: http://neo.jpl.nasa.gov/neo/groups.html

// it looks like only q<1.3 AU is required in the definition
bool OrbitID::isNEO() const {
  if (a*(1-e) < NEO_max_q) {
    return true;
  } else {
    return false;
  }
}

bool OrbitID::isIEO() const {
  if (a*(1+e) < EARTH_q) {
    return true;
  } else {
    return false;
  }
}

bool OrbitID::isAten() const {
  if ( (a < ONE_AU) && (a*(1+e) > EARTH_q) ) {
    return true;
  } else {
    return false;
  }
}

bool OrbitID::isApollo() const {
  if ( (a > ONE_AU) && (a*(1-e) < EARTH_Q) ) {
    return true;
  } else {
    return false;
  }
}

bool OrbitID::isAmor() const {
  if ( (a > ONE_AU) && (a*(1-e) > EARTH_Q) && (a*(1-e) < NEO_max_q) ) {
    return true;
  } else {
    return false;
  }
}

// Earth MOID < 0.05 AU
bool OrbitID::isPHO() const {
  
  // approximate elements
  orsa::Orbit earthOrbit;
  earthOrbit.mu = orsaSolarSystem::Data::GMSun();
  earthOrbit.a  = FromUnits(1,orsa::Unit::AU);
  earthOrbit.e  = 0.017;
  earthOrbit.i                =   0.002*orsa::degToRad();
  earthOrbit.omega_node       =   0.000*orsa::degToRad();
  earthOrbit.omega_pericenter = 100.000*orsa::degToRad();
  earthOrbit.M                =   0.000*orsa::degToRad(); // M does not matter when computing the MOID
  
  double moid, M1, M2;
  if (!orsa::MOID(moid, 
		  M1,
		  M2,
		  earthOrbit,
		  (*this),
		  19923,
		  16,
		  1e-6)) {
    ORSA_DEBUG("problems while computing MOID...");
  }
  
  // ORSA_DEBUG("moid: %f [AU]",orsa::FromUnits(moid,orsa::Unit::AU,-1));
  
  return (moid < FromUnits(0.05,orsa::Unit::AU));
}

OrbitID * OrbitFactory::sample() const {
  
  // OrbitID * orbit = new OrbitID(idCounter++,rnd->randomSeed);
  OrbitID * orbit = new OrbitID(idCounter++);
  
  orbit->a = FromUnits(a_AU_min+(a_AU_max-a_AU_min)*rnd->gsl_rng_uniform(),orsa::Unit::AU);  
  orbit->e = e_min+(e_max-e_min)*rnd->gsl_rng_uniform();
  orbit->i = orsa::degToRad()*(i_DEG_min+(i_DEG_max-i_DEG_min)*rnd->gsl_rng_uniform());
  //
  orbit->H = H_min+(H_max-H_min)*rnd->gsl_rng_uniform();
  
  orbit->omega_node       = orsa::twopi()*rnd->gsl_rng_uniform();
  orbit->omega_pericenter = orsa::twopi()*rnd->gsl_rng_uniform();
  orbit->M                = orsa::twopi()*rnd->gsl_rng_uniform();
  
  orbit->mu = GMSun;
  
  return orbit;
}
