#define RUN_IN_BOINC_CLIENT

#include "parameter.h"
#include "telescope.h"
#include "SurveySimulator.h"

#include <list>

#include <fcntl.h>

#include <orsa/body.h>
#include <orsa/bodygroup.h>
#include <orsa/double.h>
#include <orsa/matrix.h>
#include <orsa/unit.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_roots.h>

#include <orsaSPICE/spice.h>
#include <orsaSPICE/spiceBodyAttitudeCallback.h>
#include <orsaSPICE/spiceBodyPosVelCallback.h>

#include <orsaSolarSystem/datetime.h>
#include <orsaSolarSystem/gmst.h>
#include <orsaSolarSystem/obleq.h>
#include <orsaSolarSystem/observatory.h>

// BOINC API
#include <boinc_api.h>
#include <filesys.h>

// CSPICE prototypes and definitions.      
#include "SpiceUsr.h"

// globaly accessible proxy
osg::ref_ptr<PhaseComponentProxy> phaseComponentProxy = new PhaseComponentProxy(0.2);

// globaly accessible proxy
osg::ref_ptr<Log10Proxy> log10Proxy = new Log10Proxy(0.2);

osg::ref_ptr<NEO> NEOFactory::sampleNEO(const orsa::Time & orbitEpoch) const {
  
  osg::ref_ptr<NEO> newNEO = createNEO(idCounter++,
				       rnd->randomSeed);
  
  const orsa::Double sample_a_AU = ddu_a_AU->sample(rnd.get(), a_AU_min.get_d(), a_AU_max.get_d());
  if ( (sample_a_AU < a_AU_min) || 
       (sample_a_AU > a_AU_max) ) {
    ORSA_ERROR("problems");
  }
  newNEO->orbit.a = FromUnits(sample_a_AU,orsa::Unit::AU);
  
  const orsa::Double sample_e = ddu_e->sample(rnd.get(), e_min.get_d(), e_max.get_d());
  if ( (sample_e < e_min) || 
       (sample_e > e_max) ) {
    ORSA_ERROR("problems");
  }
  newNEO->orbit.e = sample_e;
  
  const orsa::Double sample_i_DEG = ddu_i_DEG->sample(rnd.get(), i_DEG_min.get_d(), i_DEG_max.get_d());
  if ( (sample_i_DEG < i_DEG_min) || 
       (sample_i_DEG > i_DEG_max) ) {
    ORSA_ERROR("problems");
  }
  newNEO->orbit.i = sample_i_DEG * orsa::degToRad();
  
  newNEO->orbit.omega_node       = orsa::twopi()*rnd->gsl_rng_uniform();
  
  newNEO->orbit.omega_pericenter = orsa::twopi()*rnd->gsl_rng_uniform();
  
  newNEO->orbit.M                = orsa::twopi()*rnd->gsl_rng_uniform();
  
  newNEO->orbit.mu = orsa::Unit::instance()->getG()*orsa::FromUnits(orsa::one(),orsa::Unit::MSUN); 
  
  newNEO->orbitEpoch = orbitEpoch;
  
  /* 
     ORSA_DEBUG("id: %6i   a: %Ff   e: %Ff   i: %Ff",
     newNEO->id,
     FromUnits(newNEO->orbit.a,orsa::Unit::AU,-1).get_mpf_t(),
     newNEO->orbit.e.get_mpf_t(),
     orsa::Double(orsa::radToDeg()*newNEO->orbit.i).get_mpf_t());
  */
  
  return newNEO;
}


/* 
   void writeDataFiles(int                  fd,
   const NEO          * neo,
   const orsa::Time   & t,
   const orsa::Double & V,
   Telescope          * telescope) {
   char line[1024];
   gmp_snprintf(line,1024,"%.5Ff %Zi %6i %.3Ff %.4Ff %6.2Ff %5.2Ff %5.2Ff [%s]",
   orsaSolarSystem::julianTime(t).get_mpf_t(),
   t.getMuSec().get_mpz_t(),
   neo->id,
   FromUnits(neo->orbit.a,orsa::Unit::AU,-1).get_mpf_t(),
   neo->orbit.e.get_mpf_t(),
   orsa::Double(neo->orbit.i*orsa::radToDeg()).get_mpf_t(),
   neo->getH(t).get_mpf_t(),
   V.get_mpf_t(),
   telescope->name.c_str());
   boinc_begin_critical_section();
   strcat(line,MODEOL);
   write(fd,line,strlen(line));
   boinc_end_critical_section();
   }
*/


void dumpSampledNEOs (const NEOList      & neoList,
		      const std::string  & fileName,
		      const orsa::Double & detectionProbabilityThreshold) {
  
  FILE * fp = fopen(fileName.c_str(),"w");
  if (!fp) {
    ORSA_DEBUG("cannot open file [%s] for writing...",fileName.c_str());
    return;
  }
  NEOList::const_iterator it = neoList.begin();
  while (it != neoList.end()) {
    
    const RealNEO      *      realNEO = dynamic_cast<const      RealNEO *> ((*it).get());
    const SyntheticNEO * syntheticNEO = dynamic_cast<const SyntheticNEO *> ((*it).get());
    
    if (realNEO) {      
      gmp_fprintf(fp,"%10i %12i %12.9Ff %12.10Ff %12.8Ff %5.2Ff"
		  MODEOL,
		  realNEO->randomSeed,
		  realNEO->id,
		  orsa::FromUnits(realNEO->orbit.a,orsa::Unit::AU,-1).get_mpf_t(),
		  realNEO->orbit.e.get_mpf_t(),
		  orsa::Double(realNEO->orbit.i*orsa::radToDeg()).get_mpf_t(),
		  realNEO->getH(orsaSolarSystem::J2000(),detectionProbabilityThreshold).get_mpf_t());      
    } else if (syntheticNEO) {
      gmp_fprintf(fp,"%10i %12i %12.9Ff %12.10Ff %12.8Ff"
		  MODEOL,
		  syntheticNEO->randomSeed,
		  syntheticNEO->id,
		  orsa::FromUnits(syntheticNEO->orbit.a,orsa::Unit::AU,-1).get_mpf_t(),
		  syntheticNEO->orbit.e.get_mpf_t(),
		  orsa::Double(syntheticNEO->orbit.i*orsa::radToDeg()).get_mpf_t());
    } else {
      ORSA_DEBUG("case not handled");
    }
    
    ++it;
  }
  fclose(fp);
}


osg::ref_ptr<orsa::Body> SPICEBody (const std::string  & bodyName,
				    const orsa::Double & bodyMass) {
  
  osg::ref_ptr<orsa::Body> localBody = new orsa::Body;
  localBody->setName(bodyName);
  // localBody->setMass(bodyMass);  
  
  osg::ref_ptr<orsaSPICE::SpiceBodyPosVelCallback> sbpvc = new orsaSPICE::SpiceBodyPosVelCallback(bodyName);
  orsa::IBPS ibps;
  ibps.inertial = new orsa::ConstantMassBodyProperty(bodyMass);
  ibps.translational = sbpvc.get();
  localBody->setInitialConditions(ibps);
  
  return localBody.get();
}

void simpleNEOVector(orsa::Vector & u_obs2neo,
		     orsa::Double & V,
		     orsa::Vector & neo2obs,
		     orsa::Vector & neo2sun,
		     orsa::Double & phaseAngle,
		     const NEO * neo,
		     const orsa::Time   & epoch,
	    	     const orsa::Vector & sunPosition,
		     const orsa::Vector & obsPosition,
		     const orsa::Vector & tp_u,
		     const orsa::Double & apertureAngle,
		     const orsa::Double & cos_apertureAngle,
		     const orsa::Double & detectionProbabilityThreshold,
		     const bool           cacheON) {
  
  if (cacheON && (neo->r_cache.isSet())) {
    const r_Cache      & rC   = neo->r_cache.getRef();
    const orsa::Time     dt   = epoch - rC.epoch.getRef();
    const orsa::Double   dMax = abs(rC.vMax.getRef() * dt.asDouble());
    const orsa::Double   dr   = (rC.r.getRef() - obsPosition).length();
    if (dMax < dr) {
      // good: NEO far enough to make cache still useful
      
      // quick test: arcsin(z)<=1.1*z for z<=0.68
      const orsa::Double ratio = dMax / dr;
      if (ratio <= 0.68) {
	const orsa::Double z = 1.1*ratio + apertureAngle;
	const orsa::Double mod_cos_scalarProduct = 1 - z*z/2;
	u_obs2neo = (rC.r.getRef() - obsPosition).normalized();
	if ((u_obs2neo*tp_u) < mod_cos_scalarProduct) {    
	  // ORSA_DEBUG("saved time by using the cache...");
	  V = 100.0;
	  return;
	}
      }
    } 
  }
  
  orsa::Orbit & localOrbit = neo->orbit;
  
  const orsa::Double orbitPeriod = localOrbit.period();
  
  const orsa::Double original_M  = localOrbit.M;
  
  localOrbit.M += orsa::fmod(orsa::twopi() * (epoch - neo->orbitEpoch).asDouble() / orbitPeriod, orsa::twopi());
  //
  orsa::Vector r;
  localOrbit.relativePosition(r);
  r += sunPosition;
  const orsa::Vector neoPosition = r;
  //
  localOrbit.M = original_M;
  
  /* 
     ORSA_DEBUG("OBS position: %Ff %Ff %Ff",
     orsa::FromUnits(obsPosition.getX(),orsa::Unit::AU,-1).get_mpf_t(),
     orsa::FromUnits(obsPosition.getY(),orsa::Unit::AU,-1).get_mpf_t(),
     orsa::FromUnits(obsPosition.getZ(),orsa::Unit::AU,-1).get_mpf_t());
  */
  //
  /* 
     ORSA_DEBUG("NEO position: %Ff %Ff %Ff",
     orsa::FromUnits(neoPosition.getX(),orsa::Unit::AU,-1).get_mpf_t(),
     orsa::FromUnits(neoPosition.getY(),orsa::Unit::AU,-1).get_mpf_t(),
     orsa::FromUnits(neoPosition.getZ(),orsa::Unit::AU,-1).get_mpf_t());
  */
  
  u_obs2neo = (neoPosition - obsPosition).normalized();
  
  if (cacheON) {
    r_Cache rC;
    rC.epoch = epoch;
    rC.r = r;
    if (neo->r_cache.isSet() && 
	neo->r_cache.getRef().vMax.isSet()) {
      rC.vMax = neo->r_cache.getRef().vMax.getRef();
    } else {
      rC.update_vMax(neo->orbit);
    }	
    neo->r_cache = rC;
  }
  
  if ((u_obs2neo*tp_u) > cos_apertureAngle) {
    
    // save these values in any case, because they're returned
    neo2obs    = obsPosition - neoPosition;
    neo2sun    = sunPosition - neoPosition;
    phaseAngle = orsa::acos((neo2obs.normalized())*(neo2sun.normalized()));
    
    /* 
       V = neo->getH(epoch,detectionProbabilityThreshold) + 
       P(phaseAngle) +
       5*orsa::log10(FromUnits(neo2obs.length(),orsa::Unit::AU,-1)*
       FromUnits(neo2sun.length(),orsa::Unit::AU,-1));
    */
    //
    /* 
       orsa::Double proxyP;
       if (phaseComponentProxy->get(proxyP,phaseAngle)) {
       V = neo->getH(epoch,detectionProbabilityThreshold) + 
       proxyP +
       5*orsa::log10(FromUnits(neo2obs.length(),orsa::Unit::AU,-1)*
       FromUnits(neo2sun.length(),orsa::Unit::AU,-1));
       } else {
       ORSA_DEBUG("problems: slow version");
       V = neo->getH(epoch,detectionProbabilityThreshold) + 
       P(phaseAngle) +
       5*orsa::log10(FromUnits(neo2obs.length(),orsa::Unit::AU,-1)*
       FromUnits(neo2sun.length(),orsa::Unit::AU,-1));
       }
    */
    //
    V = apparentMagnitude(neo->getH(epoch,detectionProbabilityThreshold),
			  phaseAngle,
			  neo2obs.length(),
			  neo2sun.length());
    
    // debug
    /* {
       const orsa::Double nominalP = P(phaseAngle);
       orsa::Double proxyP;
       if (phaseComponentProxy->get(proxyP,phaseAngle)) {
       ORSA_DEBUG("nominal: %Ff   proxy: %Ff   diff: %Ff",
       nominalP.get_mpf_t(),
       proxyP.get_mpf_t(),
       orsa::fabs(proxyP-nominalP).get_mpf_t());
       } else {
       ORSA_DEBUG("problems...");
       }
       }
    */
    
  } else {
    V = 100;
  }	
  
}


/*** SyntheticNEO::getH related code ***/ 

class synthetic_getH_parameters {
public:
  SyntheticNEO::detectionLog log;
  orsa::Time                 t;
  orsa::Double               detectionProbabilityThreshold;
};

double synthetic_getH_f (double H, void * p) {
  
  synthetic_getH_parameters * parameters = (synthetic_getH_parameters *) p;
  
  orsa::Double prob = orsa::one();  
  SyntheticNEO::detectionLog::const_iterator it = parameters->log.begin();
  while (it != parameters->log.end()) {
    if ((*it).first <= parameters->t) {
      /* 
	 const orsa::Double V =
	 H + 
	 P((*it).second->phaseAngle.getRef()) +
	 5*orsa::log10(FromUnits((*it).second->neo2obs.getRef(),orsa::Unit::AU,-1)*
	 FromUnits((*it).second->neo2sun.getRef(),orsa::Unit::AU,-1));
      */
      //
      /* 
	 orsa::Double proxyP;
	 if (phaseComponentProxy->get(proxyP,(*it).second->phaseAngle.getRef())) {
	 const orsa::Double V =
	 H + 
	 proxyP +
	 5*orsa::log10(FromUnits((*it).second->neo2obs.getRef(),orsa::Unit::AU,-1)*
	 FromUnits((*it).second->neo2sun.getRef(),orsa::Unit::AU,-1));
	 prob *= (orsa::one() - Telescope::detectionProbability(V,(*it).second->limitingMagnitude.getRef()));
	 } else {
	 ORSA_DEBUG("problems: slow version");
	 const orsa::Double V =
	 H + 
	 P((*it).second->phaseAngle.getRef()) +
	 5*orsa::log10(FromUnits((*it).second->neo2obs.getRef(),orsa::Unit::AU,-1)*
	 FromUnits((*it).second->neo2sun.getRef(),orsa::Unit::AU,-1));
	 prob *= (orsa::one() - Telescope::detectionProbability(V,(*it).second->limitingMagnitude.getRef()));
	 }
      */
      //
      const orsa::Double V = apparentMagnitude(H,
					       (*it).second->phaseAngle.getRef(),
					       (*it).second->neo2obs.getRef(),
					       (*it).second->neo2sun.getRef());
      prob *= (orsa::one() - Telescope::detectionProbability(V,(*it).second->limitingMagnitude.getRef()));
      
      /* 
	 if (boinc_is_standalone()) {
	 ORSA_DEBUG("V: %Ff",V.get_mpf_t());
	 }
      */
      // prob *= (orsa::one() - Telescope::detectionProbability(V,(*it).second->limitingMagnitude.getRef()));
    }
    ++it;
  }
  
  prob = (orsa::one() - prob);
  
  // ORSA_DEBUG("prob: %Ff",prob.get_mpf_t());
  
  return orsa::Double(prob - parameters->detectionProbabilityThreshold).get_d();
}


orsa::Double SyntheticNEO::getH(const orsa::Time   & t,
				const orsa::Double & detectionProbabilityThreshold) const {
  
  // debug
  /* 
     if (boinc_is_standalone()) {
     ORSA_DEBUG("id: %i   rS: %i",id,randomSeed);
     SyntheticNEO::detectionLog::const_iterator it = log.begin();
     while (it != log.end()) {
     ORSA_DEBUG("%Zi %12.9Ff %12.9Ff %13.9Ff %6.3Ff %s",
     (*it).first.getMuSec().get_mpz_t(),
     FromUnits((*it).second->neo2obs.getRef(),orsa::Unit::AU,-1).get_mpf_t(),
     FromUnits((*it).second->neo2sun.getRef(),orsa::Unit::AU,-1).get_mpf_t(),
     orsa::Double(orsa::radToDeg()*(*it).second->phaseAngle.getRef()).get_mpf_t(),
     (*it).second->limitingMagnitude.getRef().get_mpf_t(),
     (*it).second->telescopeName.getRef().c_str());
     ++it;
     }
     }
  */
  
  /* 
     if (getH_cache_set) {
     if ((getH_cache_t == t) &&
     (getH_cache_detectionProbabilityThreshold == detectionProbabilityThreshold)) {
     return getH_cache_H;
     }
     }
  */
  
  /* 
     if (boinc_is_standalone()) {
     ORSA_DEBUG("log.size(): %i   detectionInterval->size(): %i",log.size(),detectionInterval->size());
     }
  */
  
  bool activeLog=false;
  SyntheticNEO::detectionLog::const_iterator it = log.begin();
  while (it != log.end()) {
    if ((*it).first <= t) {
      activeLog=true;
    }
    ++it;
  }
  
  /* 
     if (boinc_is_standalone()) {
     ORSA_DEBUG("activeLog = %i",activeLog);
     }
  */
  
  if (!activeLog) {
    // returning "very bright" H, because a synthetic NEO that has never been observed 
    // does not have any constraint on H, so can be considered arbitrarily bright.
    return -100;
  }
  
  // try to use detectionInterval cache
  if (detectionInterval->size()) {
    
    if (t < detectionInterval->min().t.getRef()) {
      return -100;
    }
    
    if ( (t >= detectionInterval->min().t.getRef()) && 
	 (t <= detectionInterval->max().t.getRef()) ) {
      detectionIntervalEntry e, eMin, eMax;
      e.t = t;
      if (detectionInterval->getSubInterval(e, eMin, eMax)) {
	return eMin.H.getRef();
      } else {
	ORSA_DEBUG("problems");
      }
    }
    
    if (t > detectionInterval->max().t.getRef()) {
      if ((*log.rbegin()).first <= detectionInterval->max().t.getRef()) {
	return detectionInterval->max().H.getRef();
      }
    }
    
  }
  
  synthetic_getH_parameters parameters;
  //
  parameters.log                           = log;
  parameters.t                             = t;
  parameters.detectionProbabilityThreshold = detectionProbabilityThreshold;
  
  gsl_function F;
  //
  F.function = &synthetic_getH_f;
  F.params   = &parameters;
  
  gsl_root_fsolver * s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
  
  double H_min = detectionInterval->size() ? detectionInterval->max().H.getRef().get_d()   : 15.0;
  double H_max = detectionInterval->size() ? detectionInterval->max().H.getRef().get_d()+1 : 16.0;
  // make sure the H range above encloses the root
  {
    double f;
    unsigned int iter;
    const unsigned int maxIter = 16;
    
    /* 
       iter=0;
       do {
       ++iter;
       H_min -= 5.0;
       f = synthetic_getH_f(H_min, &parameters);
       // ORSA_DEBUG("iter: %i   H_min: %f   f: %f",iter,H_min,f);
       } while ((f < 0) && (iter < maxIter));
    */
    
    /* 
       iter=0;
       f = synthetic_getH_f(H_min, &parameters);
       while ((f < 0) && (iter < maxIter)) {
       ++iter;
       H_min -= 1.0;
       f = synthetic_getH_f(H_min, &parameters);
       // ORSA_DEBUG("iter: %i   H_min: %f   f: %f",iter,H_min,f);
       } 
    */
    
    /*  
	iter=0;
	do {
	++iter;
	H_max += 5.0;
	f = synthetic_getH_f(H_max, &parameters);
	// ORSA_DEBUG("iter: %i   H_max: %f   f: %f",iter,H_max,f);
	} while ((f > 0) && (iter < maxIter));
    */
    
    iter=0;
    f = synthetic_getH_f(H_max, &parameters);
    while ((f > 0) && (iter < maxIter)) {
      ++iter;
      H_max += 1.0;
      f = synthetic_getH_f(H_max, &parameters);
      // ORSA_DEBUG("iter: %i   H_max: %f   f: %f",iter,H_max,f);
    }
    
    H_min = H_max - 1;
    iter=0;
    f = synthetic_getH_f(H_min, &parameters);
    while ((f < 0) && (iter < maxIter)) {
      ++iter;
      H_min -= 1.0;
      f = synthetic_getH_f(H_min, &parameters);
      // ORSA_DEBUG("iter: %i   H_min: %f   f: %f",iter,H_min,f);
    } 
   
  }
  
  gsl_root_fsolver_set (s, &F, H_min, H_max);
  
  double H;
  const unsigned int maxIter = 128;
  int status;
  unsigned int iter = 0;
  do {
    iter++;
    status = gsl_root_fsolver_iterate(s);
    H      = gsl_root_fsolver_root(s);
    H_min  = gsl_root_fsolver_x_lower(s);
    H_max  = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(H_min, H_max, 0.05, 0);
    /* 
       if (boinc_is_standalone()) {
       ORSA_DEBUG("randomSeed: %10i   id: %10i   iter: %3i   H: %6.3f   [%f:%f]", 
       randomSeed, id, iter, H, H_min, H_max);
       }
    */
  } while ((status == GSL_CONTINUE) && (iter < maxIter));
  
  gsl_root_fsolver_free(s);
  
  // cache update
  // getH_cache_set = true;
  // getH_cache_t   = t;
  // getH_cache_detectionProbabilityThreshold = detectionProbabilityThreshold;
  // getH_cache_H   = orsa::Double(H);
  
  // cache update
  detectionIntervalEntry e;
  e.t = t;
  e.H = H;
  detectionInterval->insert(e);
  
  return orsa::Double(H);
}

/********/

void SyntheticNEO::logInsert(const orsa::Time   & t, 
			     NEO::LogEntry      * e,
			     const bool           writeFile) const {
  if (!e) return;
  SyntheticNEO::LogEntry * se = dynamic_cast<SyntheticNEO::LogEntry *>(e);
  if (!se) return;
  log[t] = se;
  //
  {
    DetectionIntervalType::DataType::iterator it = detectionInterval->getData().begin();
    while (it != detectionInterval->getData().end()) {
      if ((*it).t.getRef() >= t) {
	it = detectionInterval->getData().erase(it);
      } else {
	++it;
      }
    }
    detectionInterval->update();
  }
  //
  if (writeFile) {
    if (trustLogFile) {
      FILE * fp = fopen(logFileName().c_str(),"a");
      if (fp) {
	boinc_begin_critical_section();
	gmp_fprintf(fp,"%Zi %12.9Ff %12.9Ff %13.9Ff %6.3Ff %s"
		    MODEOL,
		    t.getMuSec().get_mpz_t(),
		    FromUnits(se->neo2obs.getRef(),orsa::Unit::AU,-1).get_mpf_t(),
		    FromUnits(se->neo2sun.getRef(),orsa::Unit::AU,-1).get_mpf_t(),
		    orsa::Double(orsa::radToDeg()*se->phaseAngle.getRef()).get_mpf_t(),
		    se->limitingMagnitude.getRef().get_mpf_t(),
		    se->telescopeName.getRef().c_str());
	fflush(fp);
	boinc_end_critical_section();
	fclose(fp);
      }
    } else {
      FILE * fp = fopen(logFileName().c_str(),"w");
      if (fp) {
	detectionLog::const_iterator it = log.begin();
	boinc_begin_critical_section();
	while (it != log.end()) {
	  gmp_fprintf(fp,"%Zi %12.9Ff %12.9Ff %13.9Ff %6.3Ff %s"
		      MODEOL,
		      (*it).first.getMuSec().get_mpz_t(),
		      FromUnits((*it).second->neo2obs.getRef(),orsa::Unit::AU,-1).get_mpf_t(),
		      FromUnits((*it).second->neo2sun.getRef(),orsa::Unit::AU,-1).get_mpf_t(),
		      orsa::Double(orsa::radToDeg()*(*it).second->phaseAngle.getRef()).get_mpf_t(),
		      (*it).second->limitingMagnitude.getRef().get_mpf_t(),
		      (*it).second->telescopeName.getRef().c_str());
	  ++it;
	}
	fflush(fp);
	boinc_end_critical_section();
	fclose(fp);
	trustLogFile=true;
      }
    }
  }
}

void SyntheticNEO::logPurify(const orsa::Double & detectionProbabilityThreshold) const {
  // ORSA_DEBUG("initial log size: %i",log.size());
  detectionInterval->reset();
  trustLogFile=false;
  const orsa::Double H = getH((*log.rbegin()).first, 
			      detectionProbabilityThreshold);
  detectionLog::iterator it = log.begin();
  while (it != log.end()) {
    const orsa::Double V = 
      apparentMagnitude(H,
			(*it).second->phaseAngle.getRef(),
			(*it).second->neo2obs.getRef(),
			(*it).second->neo2sun.getRef());
    const orsa::Double dp = 
      Telescope::detectionProbability(V,(*it).second->limitingMagnitude.getRef());
    if (dp <= 0) {
      log.erase(it++);
    } else {
      ++it;
    }
  }
  // ORSA_DEBUG("  final log size: %i",log.size());
}
