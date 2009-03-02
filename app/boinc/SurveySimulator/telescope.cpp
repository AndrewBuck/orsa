#include "telescope.h"
#include "SurveySimulator.h"
#include "boinc_util.h"

#include <orsaInputOutput/MPC_obscode.h>

#include <orsaSolarSystem/datetime.h>
#include <orsaSolarSystem/gmst.h>
#include <orsaSolarSystem/obleq.h>
#include <orsaSolarSystem/observatory.h>

#include <fcntl.h>

Telescope::Telescope(const orsa::Time   & start,
		     const orsa::Time   & stop,
		     const orsa::Double & limitingMagnitude_in,
		     const orsa::Double & FOV_DEG_in,
		     const orsa::Double & maxZenithDistanceAngle_DEG_in,
		     const orsa::Double & minMoonDistanceAngle_DEG_in,
		     const orsa::Double & minMoonPhase_DEG_in,
		     const int            recycleTime_DAY_in,
		     const int            dutyCycle_SEC_in,
		     const int            dutyCycleMultiplicity_in,
		     const std::string  & obscode_in,
		     const std::string  & name_in) : 
  osg::Referenced(), 
  tStart(start),
  tStop(stop), 
  limitingMagnitude(limitingMagnitude_in),
  FOV_DEG(FOV_DEG_in),
  maxZenithDistanceAngle_DEG(maxZenithDistanceAngle_DEG_in),
  minMoonDistanceAngle_DEG(minMoonDistanceAngle_DEG_in),
  minMoonPhase_DEG(minMoonPhase_DEG_in),
  recycleTime_DAY(recycleTime_DAY_in),
  dutyCycle_SEC(dutyCycle_SEC_in),
  dutyCycleMultiplicity(dutyCycleMultiplicity_in),
  obscode(obscode_in),
  name(name_in),
  effectiveDutyCycle(orsa::Time(0,0,0,dutyCycle_SEC,0)*dutyCycleMultiplicity),
  recycleTime(orsa::Time(recycleTime_DAY,0,0,0,0)),
  FOV(FOV_DEG*orsa::degToRad()),
  // cos_FOV(orsa::cos(FOV)),
  effectiveHalfFOV(1.15*(0.5*FOV)),
  cos_effectiveHalfFOV(orsa::cos(effectiveHalfFOV)),
  maxZenithDistanceAngle(maxZenithDistanceAngle_DEG*orsa::degToRad()),
  cos_maxZenithDistanceAngle(orsa::cos(maxZenithDistanceAngle_DEG*orsa::degToRad())),
  cos_minMoonDistanceAngle(orsa::cos(minMoonPhase_DEG*orsa::degToRad())),
  eclipticToEquatorial(orsa::Matrix::identity().rotX(+orsaSolarSystem::obleqJ2000())),
  equatorialToEcliptic(orsa::Matrix::identity().rotX(-orsaSolarSystem::obleqJ2000()))
{ 
  std::string resolvedFileName;
  boinc_resolve_filename_s(logFileName().c_str(),resolvedFileName);
  fd_tp = open(resolvedFileName.c_str(), open_flag, open_mode);
}


void Telescope::writeTP(const TelescopePointing & localTP,
			const orsa::Vector      & sunPosition,
			const orsa::Vector      & obsPosition) const {
  
  const orsa::Double lambda = 
    orsa::fmod(orsa::twopi() + 
	       orsa::atan2(localTP.u.getY(),
			   localTP.u.getX()),
	       orsa::twopi());
  const orsa::Double beta = orsa::halfpi() - orsa::acos(localTP.u.getZ());
  
  const orsa::Vector u_opposition = -(sunPosition-obsPosition).normalized();
  const orsa::Double lambdaOpposition = 
    orsa::fmod(orsa::twopi() + 
	       orsa::atan2(u_opposition.getY(),
			   u_opposition.getX()),
	       orsa::twopi());
  
  const orsa::Double dLambda = orsa::fmod(orsa::twopi() + 
					  lambda-lambdaOpposition,
					  orsa::twopi());
  
  char line[1024];
  gmp_snprintf(line,1024,"%.5Ff %Zi %+.12Ff %+.12Ff %+.12Ff %7.3Ff %+7.3Ff %7.3Ff %.3Ff [%s]",
	       orsaSolarSystem::timeToJulian(localTP.epoch).get_mpf_t(),
	       localTP.epoch.getMuSec().get_mpz_t(),
	       localTP.u.getX().get_mpf_t(),
	       localTP.u.getY().get_mpf_t(),
	       localTP.u.getZ().get_mpf_t(),
	       orsa::Double(dLambda*orsa::radToDeg()).get_mpf_t(),
	       orsa::Double(beta*orsa::radToDeg()).get_mpf_t(),
	       orsa::Double(orsa::acos(u_opposition*localTP.u)*orsa::radToDeg()).get_mpf_t(),
	       FromUnits((sunPosition-obsPosition).length(),orsa::Unit::AU,-1).get_mpf_t(),
	       name.c_str());
  boinc_begin_critical_section();
  strcat(line,MODEOL);
  write(fd_tp,line,strlen(line));
  boinc_end_critical_section();
}


const orsaSolarSystem::Observatory & Telescope::getObservatory() const {
  if (observatoryCache.isSet()) {
    return observatoryCache.getRef();
  }
  std::string resolvedFileName;
  boinc_resolve_filename_s("obscode.dat",resolvedFileName);
  osg::ref_ptr<orsaInputOutput::MPCObsCodeFile> obscodeFile = new orsaInputOutput::MPCObsCodeFile;
  obscodeFile->setFileName(resolvedFileName);
  obscodeFile->read();
  //
  observatoryCache = obscodeFile->_data.observatory[obscode];
  //
  return observatoryCache.getRef();
}


bool OppositionTelescope::sampleTP(TelescopePointing   & localTP,
				   const orsa::Time    & t,
				   const TelescopeList & telescopeList,
				   const NEOList       & ,
				   const orsa::Double  & ,
				   const orsa::Vector  &  sunPosition,
				   const orsa::Vector  & moonPosition,
				   const orsa::Vector  &  obsPosition,
				   const orsa::Vector  &  obsNormal,
				   const bool            ) {
  
  if (!isRunning(t)) {
    return false;
  }
  
  if (lastObservationTime.isSet()) {
    if ((t-lastObservationTime.getRef()) < effectiveDutyCycle) {
      // ORSA_DEBUG("--MARK-- name: %s [early exit]",name.c_str());
      return false;
    } 
  } else {
    // ORSA_DEBUG("--MARK-- name: %s [--observing--]",name.c_str());
  }	
  
  localTP.epoch = t;
  
  const orsa::Vector u_opposition = -(sunPosition-obsPosition).normalized();
  // const orsa::Vector u_zenith     = (obsNormal).normalized();
  const orsa::Vector u_northPole  = (equatorialToEcliptic*orsa::Vector(0,0,1)).normalized();
  const orsa::Vector u_ortho      = (orsa::externalProduct(u_northPole,u_opposition)).normalized();
  
  for (unsigned int m=0; m<1000000; ++m) {
    localTP.u = u_opposition;
    const orsa::Double sqrtM = sqrt(m);
    localTP.u = orsa::Matrix::axisRotation(u_ortho,
					   FOV*sqrtM*(2.0*rnd->gsl_rng_uniform()-1.0)) * localTP.u;
    localTP.u = orsa::Matrix::axisRotation(u_northPole,
					   2.0*FOV*sqrtM*(2.0*rnd->gsl_rng_uniform()-1.0)) * localTP.u;
    localTP.u.normalize();
    
    // const orsa::Double moonDistanceAngle = orsa::acos(((moonPosition-obsPosition).normalized())*localTP.u);
    //
    const orsa::Double cos_moonDistanceAngle = ((moonPosition-obsPosition).normalized()) * localTP.u;
    
    // const orsa::Double zenithDistanceAngle = orsa::acos(obsNormal*localTP.u);
    // const orsa::Double elevation = orsa::halfpi() - orsa::acos(obsNormal*localTP.u);
    //
    const orsa::Double cos_zenithDistanceAngle = obsNormal*localTP.u;
    
    /* 
       if ( ((moonDistanceAngle*orsa::radToDeg()) > 45.0) &&
       ((elevation*orsa::radToDeg()) > 45.0) ) {
    */
    //
    if ( (cos_zenithDistanceAngle > cos_maxZenithDistanceAngle) &&
	 (cos_moonDistanceAngle   < cos_minMoonDistanceAngle) ) {
      
      bool tooClose=false;
      //
      {
	TelescopeList::const_iterator tl_it = telescopeList.begin();
	while (tl_it != telescopeList.end()) {
	  if ((*tl_it)->limitingMagnitude+1.05 <= limitingMagnitude) {
	    /* 
	       ORSA_DEBUG("NOT using external TPs, this Vlim: %Ff   other Vlim: %Ff",
	       limitingMagnitude.get_mpf_t(),
	       (*tl_it)->limitingMagnitude.get_mpf_t());
	    */
	    ++tl_it;
	    continue;
	  } 
	  /* 
	     ORSA_DEBUG("using external TPs, this Vlim: %Ff   other Vlim: %Ff",
	     limitingMagnitude.get_mpf_t(),
	     (*tl_it)->limitingMagnitude.get_mpf_t());
	  */
	  tpList::const_iterator tp_it = (*tl_it)->getTPList().begin();
	  while (tp_it != (*tl_it)->getTPList().end()) {
	    if (((*tp_it).u*localTP.u) > cos_effectiveHalfFOV) {
	      tooClose=true;
	    }
	    if (tooClose) break;
	    ++tp_it;
	  }
	  if (tooClose) break;
	  ++tl_it;
	}
      }
      
      if (!tooClose) {
	
	alltp.push_back(localTP);
	
	// ORSA_DEBUG("setting observation time for %s",name.c_str());
	lastObservationTime = localTP.epoch;
	
	writeTP(localTP,
		sunPosition,
		obsPosition);
	
	return true;
      }
    }
  }
  
  ORSA_DEBUG("telescope [%s] did not find any good TP...",
	     name.c_str());
  
  {
    ORSA_DEBUG("no good TP found, aborting.");
    boinc_finish(0);
  }
  
  return false;
}


bool SieveTelescope::sampleTP(TelescopePointing   & localTP,
			      const orsa::Time    & t,
			      const TelescopeList & telescopeList,
			      const NEOList       & syntheticNEO,
			      const orsa::Double  & detectionProbabilityThreshold,
			      const orsa::Vector  &  sunPosition,
			      const orsa::Vector  & moonPosition,
			      const orsa::Vector  &  obsPosition,	
			      const orsa::Vector  &  obsNormal,
			      const bool            cacheON) {
  
  if (!isRunning(t)) {
    return false;
  }
  
  if (lastObservationTime.isSet()) {
    if ((t-lastObservationTime.getRef()) < effectiveDutyCycle) {
      // ORSA_DEBUG("--MARK-- name: %s [early exit]",name.c_str());
      return false;
    } 
  } else {
    // ORSA_DEBUG("--MARK-- name: %s [--observing--]",name.c_str());
  }	
  
  localTP.epoch = t;
  
  // const orsa::Vector u_opposition = -(sunPosition-obsPosition).normalized();
  const orsa::Vector u_zenith     = (obsNormal).normalized();
  const orsa::Vector u_northPole  = (equatorialToEcliptic*orsa::Vector(0,0,1)).normalized();
  const orsa::Vector u_ortho      = (orsa::externalProduct(u_northPole,u_zenith)).normalized();
  
  std::list<uVit> synthetic_uVit;
  //
  {
    orsa::Vector u_obs2neo;
    orsa::Double V; // apparent magnitude
    orsa::Vector neo2obs;
    orsa::Vector neo2sun;
    orsa::Double phaseAngle;
    
    uVit tmp_uVit;
    
    NEOList::const_iterator it = syntheticNEO.begin();
    while (it != syntheticNEO.end()) {
      
      // for now, skip all synthetic NEOs with a long log
      /* 
	 {
	 const SyntheticNEO * syntheticNEO = dynamic_cast<const SyntheticNEO *> ((*it).get());
	 if (!syntheticNEO) {
	 ORSA_DEBUG("problems...");
	 } else {
	 if (syntheticNEO->getLog().size() >= 3) {
	 ++it;
	 continue;
	 }
	 }
	 }
      */
      
      simpleNEOVector(u_obs2neo,
		      V,
		      neo2obs,
		      neo2sun,
		      phaseAngle,
		      (*it).get(),
		      t,
		      sunPosition,
		      obsPosition,
		      u_zenith,
		      maxZenithDistanceAngle,
		      cos_maxZenithDistanceAngle,
		      detectionProbabilityThreshold,
		      cacheON);
      
      if ((u_obs2neo*u_zenith) > cos_maxZenithDistanceAngle) {  
	if (detectionProbability(V) > 0) {
	  tmp_uVit.u  = u_obs2neo;
	  tmp_uVit.V  = V;
	  tmp_uVit.it = it;
	  //
	  synthetic_uVit.push_back(tmp_uVit);
	}
      }
      
      ++it;
    }
  }
  
  orsa::Cache<orsa::Vector> selected_u;
  std::list<uVit>           selected_synthetic_uVit;
  orsa::Double              selected_synthetic_total_prob = orsa::zero();
  //
  // const unsigned int num_iter = synthetic_uVit.size();
  unsigned int count=0;
  std::list<uVit>::const_iterator synthetic_uVit_it = synthetic_uVit.begin();
  while (synthetic_uVit_it != synthetic_uVit.end()) {
    
    if (count*count >= synthetic_uVit.size()) {
      break;
    }
    
    // ORSA_DEBUG("num_iter: %i",num_iter);
    // for (unsigned int m=0; m<num_iter; ++m) {
    
    // random sampling...
    {
      const long unsigned int forward = rnd->gsl_rng_uniform_int(synthetic_uVit.size());
      for (long unsigned int j=0; j<forward; ++j) {
	++synthetic_uVit_it;
	if (synthetic_uVit_it == synthetic_uVit.end()) {
	  synthetic_uVit_it = synthetic_uVit.begin();
	}
      }
    }
    
    localTP.u = (*synthetic_uVit_it).u;
    
    /* 
       if (synthetic_uVit_it != synthetic_uVit.end()) {
       localTP.u = (*synthetic_uVit_it).u;
       ++synthetic_uVit_it;
       } else {
       break;
       }
    */
    //
    /* 
       } else {
       localTP.u = u_zenith;
       localTP.u = orsa::Matrix::axisRotation(u_ortho,
       maxZenithDistanceAngle*rnd->gsl_rng_uniform()) * localTP.u;
       localTP.u = orsa::Matrix::axisRotation(u_northPole,
       orsa::twopi()*rnd->gsl_rng_uniform()) * localTP.u;
       localTP.u.normalize();
       }
    */
    
    const orsa::Double cos_moonDistanceAngle   = ((moonPosition-obsPosition).normalized())*localTP.u;
    const orsa::Double cos_zenithDistanceAngle = obsNormal*localTP.u;
    
    if ( (cos_zenithDistanceAngle > cos_maxZenithDistanceAngle) &&
	 (cos_moonDistanceAngle < cos_minMoonDistanceAngle) ) {
      
      bool tooClose=false;
      //
      {
	TelescopeList::const_iterator tl_it = telescopeList.begin();
	while (tl_it != telescopeList.end()) {
	  if ((*tl_it)->limitingMagnitude+1.05 <= limitingMagnitude) {
	    ++tl_it;
	    continue;
	  } 
	  tpList::const_iterator tp_it = (*tl_it)->getTPList().begin();
	  while (tp_it != (*tl_it)->getTPList().end()) {
	    if (((*tp_it).u*localTP.u) > cos_effectiveHalfFOV) {
	      tooClose=true;
	    }
	    if (tooClose) break;
	    ++tp_it;
	  }
	  if (tooClose) break;
	  ++tl_it;
	}
      }
      
      if (!tooClose) {
	
	std::list<uVit> local_selected_synthetic_uVit;
	orsa::Double    local_selected_synthetic_total_prob = orsa::zero();
	
	std::list<uVit>::const_iterator it = synthetic_uVit.begin();
	while (it != synthetic_uVit.end()) {
	  if (((*it).u*localTP.u) > cos_effectiveHalfFOV) {  
	    const orsa::Double prob = detectionProbability((*it).V);
	    // if (prob > 0) {
	    local_selected_synthetic_uVit.push_back(*it);
	    local_selected_synthetic_total_prob += prob;
	    // }
	  }
	  ++it;
	}
	
	// this should be refined and verified: better more NEOs or bigger probability?
	//
	if (local_selected_synthetic_total_prob > selected_synthetic_total_prob) {
	  if (boinc_is_standalone()) {
	    ORSA_DEBUG("better field: %3i synthetic NEOs   p.tot: %.3Ff   synthetic_uVit.size(): %i",
		       local_selected_synthetic_uVit.size(),
		       local_selected_synthetic_total_prob.get_mpf_t(),
		       synthetic_uVit.size());
	  }	
	  
	  selected_u                    = localTP.u;
	  selected_synthetic_uVit       = local_selected_synthetic_uVit;
	  selected_synthetic_total_prob = local_selected_synthetic_total_prob;
	}
      }
    }
    
    ++count;
    ++synthetic_uVit_it;
  }
  
  if (selected_u.isSet()) {
    
    localTP.epoch = t;
    localTP.u     = selected_u.getRef();
    
    alltp.push_back(localTP);
    
    lastObservationTime = localTP.epoch;
    
    writeTP(localTP,
	    sunPosition,
	    obsPosition);
    
    return true;
  }
  
  {
    // since SIEVE failed to find a TP (most likely because not enough synthetic NEOs)
    // we revert to ALLSKY to find just a TP...
    ORSA_DEBUG("reverting to allsky...");
    
    // temporary vars
    double sx, sy, sz;
    
    for (unsigned int m=0; m<1000000; ++m) {
      
      rnd->gsl_ran_dir_3d(&sx,&sy,&sz);
      localTP.u = orsa::Vector(sx,sy,sz);    
      localTP.u.normalize();
      
      // const orsa::Double moonDistanceAngle = orsa::acos(((moonPosition-obsPosition).normalized())*localTP.u);
      //
      const orsa::Double cos_moonDistanceAngle = ((moonPosition-obsPosition).normalized()) * localTP.u;
      
      // const orsa::Double zenithDistanceAngle = orsa::acos(obsNormal*localTP.u);
      // const orsa::Double elevation = orsa::halfpi() - orsa::acos(obsNormal*localTP.u);
      //
      const orsa::Double cos_zenithDistanceAngle = obsNormal*localTP.u;
      
      /* 
	 if ( ((moonDistanceAngle*orsa::radToDeg()) > 45.0) &&
	 ((elevation*orsa::radToDeg()) > 45.0) ) {
      */
      //
      if ( (cos_zenithDistanceAngle > cos_maxZenithDistanceAngle) &&
	   (cos_moonDistanceAngle   < cos_minMoonDistanceAngle) ) {
	
	bool tooClose=false;
	//
	{
	  TelescopeList::const_iterator tl_it = telescopeList.begin();
	  while (tl_it != telescopeList.end()) {
	    if ((*tl_it)->limitingMagnitude+1.05 <= limitingMagnitude) {
	      /* 
		 ORSA_DEBUG("NOT using external TPs, this Vlim: %Ff   other Vlim: %Ff",
		 limitingMagnitude.get_mpf_t(),
		 (*tl_it)->limitingMagnitude.get_mpf_t());
	      */
	      ++tl_it;
	      continue;
	    } 
	    /* 
	       ORSA_DEBUG("using external TPs, this Vlim: %Ff   other Vlim: %Ff",
	       limitingMagnitude.get_mpf_t(),
	       (*tl_it)->limitingMagnitude.get_mpf_t());
	    */
	    tpList::const_iterator tp_it = (*tl_it)->getTPList().begin();
	    while (tp_it != (*tl_it)->getTPList().end()) {
	      if (((*tp_it).u*localTP.u) > cos_effectiveHalfFOV) {
		tooClose=true;
	      }
	      if (tooClose) break;
	      ++tp_it;
	    }
	    if (tooClose) break;
	    ++tl_it;
	  }
	}
	
	if (!tooClose) {
	  
	  alltp.push_back(localTP);
	  
	  // ORSA_DEBUG("setting observation time for %s",name.c_str());
	  lastObservationTime = localTP.epoch;
	  
	  writeTP(localTP,
		  sunPosition,
		  obsPosition);
	  
	  return true;
	}
      }
    }
  }
  
  ORSA_DEBUG("telescope [%s] did not find any good TP...",
	     name.c_str());
  
  {
    ORSA_DEBUG("no good TP found, aborting.");
    boinc_finish(0);
  }
  
  return false;
}



bool AllSkyTelescope::sampleTP(TelescopePointing   & localTP,
			       const orsa::Time    & t,
			       const TelescopeList & telescopeList,
			       const NEOList       & ,
			       const orsa::Double  & ,
			       const orsa::Vector  &  sunPosition,
			       const orsa::Vector  & moonPosition,
			       const orsa::Vector  &  obsPosition,
			       const orsa::Vector  &  obsNormal,
			       const bool            ) {
  
  if (!isRunning(t)) {
    return false;
  }
  
  if (lastObservationTime.isSet()) {
    if ((t-lastObservationTime.getRef()) < effectiveDutyCycle) {
      // ORSA_DEBUG("--MARK-- name: %s [early exit]",name.c_str());
      return false;
    } 
  } else {
    // ORSA_DEBUG("--MARK-- name: %s [--observing--]",name.c_str());
  }	
  
  localTP.epoch = t;
  
  // const orsa::Vector u_opposition = -(sunPosition-obsPosition).normalized();
  // const orsa::Vector u_zenith     = (obsNormal).normalized();
  // const orsa::Vector u_northPole  = (equatorialToEcliptic*orsa::Vector(0,0,1)).normalized();
  // const orsa::Vector u_ortho      = (orsa::externalProduct(u_northPole,u_zenith)).normalized();
  
  // temporary vars
  double sx, sy, sz;
  
  for (unsigned int m=0; m<1000000; ++m) {
    
    rnd->gsl_ran_dir_3d(&sx,&sy,&sz);
    localTP.u = orsa::Vector(sx,sy,sz);    
    localTP.u.normalize();
    
    // const orsa::Double moonDistanceAngle = orsa::acos(((moonPosition-obsPosition).normalized())*localTP.u);
    //
    const orsa::Double cos_moonDistanceAngle = ((moonPosition-obsPosition).normalized()) * localTP.u;
    
    // const orsa::Double zenithDistanceAngle = orsa::acos(obsNormal*localTP.u);
    // const orsa::Double elevation = orsa::halfpi() - orsa::acos(obsNormal*localTP.u);
    //
    const orsa::Double cos_zenithDistanceAngle = obsNormal*localTP.u;
    
    /* 
       if ( ((moonDistanceAngle*orsa::radToDeg()) > 45.0) &&
       ((elevation*orsa::radToDeg()) > 45.0) ) {
    */
    //
    if ( (cos_zenithDistanceAngle > cos_maxZenithDistanceAngle) &&
	 (cos_moonDistanceAngle   < cos_minMoonDistanceAngle) ) {
      
      bool tooClose=false;
      //
      {
	TelescopeList::const_iterator tl_it = telescopeList.begin();
	while (tl_it != telescopeList.end()) {
	  if ((*tl_it)->limitingMagnitude+1.05 <= limitingMagnitude) {
	    /* 
	       ORSA_DEBUG("NOT using external TPs, this Vlim: %Ff   other Vlim: %Ff",
	       limitingMagnitude.get_mpf_t(),
	       (*tl_it)->limitingMagnitude.get_mpf_t());
	    */
	    ++tl_it;
	    continue;
	  } 
	  /* 
	     ORSA_DEBUG("using external TPs, this Vlim: %Ff   other Vlim: %Ff",
	     limitingMagnitude.get_mpf_t(),
	     (*tl_it)->limitingMagnitude.get_mpf_t());
	  */
	  tpList::const_iterator tp_it = (*tl_it)->getTPList().begin();
	  while (tp_it != (*tl_it)->getTPList().end()) {
	    if (((*tp_it).u*localTP.u) > cos_effectiveHalfFOV) {
	      tooClose=true;
	    }
	    if (tooClose) break;
	    ++tp_it;
	  }
	  if (tooClose) break;
	  ++tl_it;
	}
      }
      
      if (!tooClose) {
	
	alltp.push_back(localTP);
	
	// ORSA_DEBUG("setting observation time for %s",name.c_str());
	lastObservationTime = localTP.epoch;
	
	writeTP(localTP,
		sunPosition,
		obsPosition);
	
	return true;
      }
    }
  }
  
  ORSA_DEBUG("telescope [%s] did not find any good TP...",
	     name.c_str());
  
  {
    ORSA_DEBUG("no good TP found, aborting.");
    boinc_finish(0);
  }
  
  return false;
}


bool readTelescope(TelescopeList     & tl,
		   const std::string & fileName) {
  
#ifdef RUN_IN_BOINC_CLIENT
  std::string resolvedFileName;
  boinc_resolve_filename_s(fileName.c_str(),resolvedFileName);
  FILE * fp = boinc_fopen(resolvedFileName.c_str(),"r");
#else 
  FILE * fp = fopen(fileName.c_str(),"r");
#endif
  
  if (!fp) {
    return false;
  }
  
  char line[1024];
  
  int  mode;
  int  start_DAY, stop_DAY;
  int  randomSeed;
  orsa::Double limitingMagnitude;
  orsa::Double FOV_DEG;
  orsa::Double maxZenithDistanceAngle_DEG;
  orsa::Double minMoonDistanceAngle_DEG;
  orsa::Double minMoonPhase_DEG;
  int  recycleTime_DAY;
  int  dutyCycle_SEC;
  int  dutyCycleMultiplicity;
  char obscode[1024];
  char name[1024];
  
  // unique names
  std::map< std::string, orsa::Cache<bool> > nameMap;
  
  while (fgets(line,1024,fp) != 0) {
    if (strlen(line) >= 1) {
      if (line[0] == '#') {
	continue;
      }
      if (14 == gmp_sscanf(line,"%d %d %d %d %Ff %Ff %Ff %Ff %Ff %d %d %d %s %s",
			   &mode,
			   &start_DAY,
			   &stop_DAY,
			   &randomSeed,
			   limitingMagnitude.get_mpf_t(),
			   FOV_DEG.get_mpf_t(),
			   maxZenithDistanceAngle_DEG.get_mpf_t(),
			   minMoonDistanceAngle_DEG.get_mpf_t(),
			   minMoonPhase_DEG.get_mpf_t(),
			   &recycleTime_DAY,
			   &dutyCycle_SEC,
			   &dutyCycleMultiplicity,
			   obscode,
			   name)) {
	
	// unique names
	if (nameMap[name].isSet()) {
	  ORSA_ERROR("name [%s] used more than once", name);
	  return false;
	}
	//
	nameMap[name]=true;
	
	/* modes:
	   1 = opposition
	   2 = sieve
	   3 = allsky
	*/
	//
	switch (mode) {
	case 1:
	  tl.push_back(new OppositionTelescope(orsaSolarSystem::J2000() + orsa::Time(start_DAY,0,0,0,0),
					       orsaSolarSystem::J2000() + orsa::Time( stop_DAY,0,0,0,0),
					       randomSeed,
					       limitingMagnitude,
					       FOV_DEG,
					       maxZenithDistanceAngle_DEG,
					       minMoonDistanceAngle_DEG,
					       minMoonPhase_DEG,
					       recycleTime_DAY,
					       dutyCycle_SEC,
					       dutyCycleMultiplicity,
					       obscode,
					       name));
	  break;
	case 2:
	  tl.push_back(new SieveTelescope(orsaSolarSystem::J2000() + orsa::Time(start_DAY,0,0,0,0),
					  orsaSolarSystem::J2000() + orsa::Time( stop_DAY,0,0,0,0),
					  randomSeed,
					  limitingMagnitude,
					  FOV_DEG,
					  maxZenithDistanceAngle_DEG,
					  minMoonDistanceAngle_DEG,
					  minMoonPhase_DEG,
					  recycleTime_DAY,
					  dutyCycle_SEC,
					  dutyCycleMultiplicity,
					  obscode,
					  name));
	  break;
	case 3:
	  tl.push_back(new AllSkyTelescope(orsaSolarSystem::J2000() + orsa::Time(start_DAY,0,0,0,0),
					   orsaSolarSystem::J2000() + orsa::Time( stop_DAY,0,0,0,0),
					   randomSeed,
					   limitingMagnitude,
					   FOV_DEG,
					   maxZenithDistanceAngle_DEG,
					   minMoonDistanceAngle_DEG,
					   minMoonPhase_DEG,
					   recycleTime_DAY,
					   dutyCycle_SEC,
					   dutyCycleMultiplicity,
					   obscode,
					   name));
	  break;
	default:
	  ORSA_ERROR("unknown mode %i",mode);
	}
	
	gmp_fprintf(stderr,
		    "------------------------------------------------\n"
		    "start_DAY.................: %10i\n"
		    "stop_DAY..................: %10i\n"
		    "randomSeed................: %10i\n"
		    "limitingMagnitude.........: %20.9Ff\n"
		    "FOV_DEG...................: %20.9Ff\n"
		    "maxZenithDistanceAngle_DEG: %20.9Ff\n"
		    "minMoonDistanceAngle_DEG..: %20.9Ff\n"
		    "minMoonPhase_DEG..........: %20.9Ff\n"
		    "recycleTime_DAY...........: %10i\n"
		    "dutyCycle_SEC.............: %10i\n"
		    "dutyCycleMultiplicity.....: %10i\n"
		    "obscode...................: %s\n"
		    "name......................: %s\n"
		    ,
		    start_DAY,
		    stop_DAY,
		    randomSeed,
		    limitingMagnitude.get_mpf_t(),
		    FOV_DEG.get_mpf_t(),
		    maxZenithDistanceAngle_DEG.get_mpf_t(),
		    minMoonDistanceAngle_DEG.get_mpf_t(),
		    minMoonPhase_DEG.get_mpf_t(),
		    recycleTime_DAY,
		    dutyCycle_SEC,
		    dutyCycleMultiplicity,
		    obscode,
		    name);
      }
    }
  }
  
  fclose(fp);
  
  return true;
}


