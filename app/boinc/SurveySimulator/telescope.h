#ifndef _SS_TELESCOPE_H_
#define _SS_TELESCOPE_H_

#include <orsa/datetime.h>
#include <orsa/util.h>

#include <orsaSolarSystem/datetime.h>
#include <orsaSolarSystem/observatory.h>

#ifdef RUN_IN_BOINC_CLIENT
#include <boinc_api.h>
#include <filesys.h>
#endif

#include "boinc_util.h"
#include "SurveySimulator.h"


class TelescopePointing {
 public:
  orsa::Vector u;
 public:
  orsa::Time   epoch;
 public:
  bool operator < (const TelescopePointing & rhs) const {
    return (epoch < rhs.epoch);
  }
};


typedef std::list<TelescopePointing> tpList;


class Telescope;


typedef std::list< osg::ref_ptr<Telescope> > TelescopeList;


class Telescope : public osg::Referenced {
 public:
  Telescope(const orsa::Time   & start,
	    const orsa::Time   & stop,
	    const double & limitingMagnitude_in,
	    const double & FOV_DEG_in,
	    const double & maxZenithDistanceAngle_DEG_in,
	    const double & minMoonDistanceAngle_DEG_in,
	    const double & minMoonPhase_DEG_in,
	    const int            recycleTime_DAY_in,
	    const int            dutyCycle_SEC_in,
	    const int            dutyCycleMultiplicity_in,
	    const std::string  & obscode_in,
	    const std::string  & name_in);
 protected:
  virtual ~Telescope() { 
    close(fd_tp);
  }
 protected:
  int fd_tp;

 public:  
  const tpList & getTPList() const {
    return alltp;
  }
 protected:
  tpList alltp;
  
 public:
  std::string logFileName() const {
    char filename[1024];
    snprintf(filename,1024,"tp.%s.%s.log",obscode.c_str(),name.c_str());
    return filename;
  }
  
 public:
  std::string rngFileName() const {
    char filename[1024];
    snprintf(filename,1024,"sky.RNG.%s.%s.bin",obscode.c_str(),name.c_str());
    return filename;
  }
  
 public:
  virtual bool sampleTP(TelescopePointing   & tp,
			const orsa::Time    & t,
			const TelescopeList & telescopeList,
			const NEOList       & syntheticNEO,
			const double  & detectionProbabilityThreshold,
			const orsa::Vector  &  sunPosition,
			const orsa::Vector  & moonPosition,
			const orsa::Vector  &  obsPosition,
			const orsa::Vector  &  obsNormal,
			const bool            cacheON) = 0;
 public:
  // call this when sampleTP cannot be called, to keep track of tentatives or past TPs
  virtual void markObservationTime(const orsa::Time & t) {
    lastObservationTime = t;
  }
  
 protected:
  void writeTP(const TelescopePointing & tp,
	       const orsa::Vector      & sunPosition,
	       const orsa::Vector      & obsPosition) const;
  
 public:
  void applyRecycleTime(const orsa::Time & t) {
    alltp.tpList::sort();
    tpList::iterator it = alltp.begin();
    while (it != alltp.end()) {
      if ((t - (*it).epoch) > recycleTime) {
	it = alltp.erase(it);
      } else {
	++it;
      }
    }
  }
  
 public:  
  virtual void writeCheckpoint() const { }
  
 public:  
  virtual void readCheckpoint(const orsa::Time & t) {
    
    // ORSA_DEBUG("------------ [%s] ----------",name.c_str());
    
    char line[1024];
    char goodLine[1024];
    off_t lastGoodFilePosition;
    double ux, uy, uz;
    
    const mpz_class muSec = t.getMuSec();
    
    lastGoodFilePosition=lseek(fd_tp, 0, SEEK_SET);
    while (read(fd_tp,line,1024) > 0) {
      
      const std::string lineString(line);
      const size_t newLinePos = lineString.find_first_of('\n'); 
      strncpy(goodLine,line,newLinePos);
      goodLine[newLinePos] = '\0';
      
      // ORSA_DEBUG("goodLine: [%s]   newLinePos: %i",goodLine,newLinePos);
      
      if (gmp_sscanf(goodLine,"%*f %Zi %lf %lf %lf",
		     muSec.get_mpz_t(),
		     &ux,
		     &uy,
		     &uz) != 4) {
	break;
      }
      
      // ORSA_DEBUG("read time %Zi",muSec.get_mpz_t());
      
      if (muSec <= t.getMuSec()) {
	TelescopePointing localTP;
	//
	localTP.epoch = orsa::Time(muSec);
	localTP.u = orsa::Vector(ux,uy,uz).normalized();
	//
	alltp.push_back(localTP);
	
	markObservationTime(orsa::Time(muSec));
	
	if (newLinePos > 0) {
	  lastGoodFilePosition = lseek(fd_tp, 
				       lastGoodFilePosition+newLinePos+lineSkip, 
				       SEEK_SET);
	} else {
	  // ORSA_DEBUG("break...");
	  break;
	}	  
      } else {
	break;
      }
    }
    lseek(fd_tp, 
	  lastGoodFilePosition, 
	  SEEK_SET);
  }
  
 public:
  const orsaSolarSystem::Observatory & getObservatory() const;
 protected:
  mutable orsa::Cache<orsaSolarSystem::Observatory> observatoryCache;
  
 public:
  orsa::Time nextObservationTime (const orsa::Time & tMin) {
    
    if (tMin < tStart) {
      // ORSA_DEBUG("--MARK-- name: %s",name.c_str());
      return tStart;
    }
    
    if (tMin > tStop) {
      ORSA_ERROR("should not call this method if tMin >= tStop");
      return tMin;
    }
    
    if (!lastObservationTime.isSet()) {
      // ORSA_DEBUG("--MARK-- name: %s",name.c_str());
      return tMin;
    }
    
    if (lastObservationTime.getRef() == tMin) {
      // ORSA_DEBUG("--MARK-- name: %s",name.c_str());
      return (lastObservationTime.getRef()+effectiveDutyCycle);
    }
    
    const mpz_class muDelta = (tMin - lastObservationTime.getRef()).getMuSec();
    
    const mpz_class muDuty  = effectiveDutyCycle.getMuSec();
    
    const mpz_class k = muDelta/muDuty;
    
    orsa::Time retVal = lastObservationTime.getRef() + k*effectiveDutyCycle;
    if (retVal < tMin) {
      // ORSA_DEBUG("--MARK-- name: %s",name.c_str());
      retVal += effectiveDutyCycle;
    }
    
    // check
    if (retVal < tMin) {
      ORSA_ERROR("problems");
    }
    
    return retVal;
  }
 protected:
  orsa::Cache<orsa::Time> lastObservationTime;
  
 public:
  bool isRunning(const orsa::Time & t) const {
    if (t < tStart) return false;
    if (t > tStop)  return false;
    return true;
  }
  
 public:
  // NOT virtual: all telescopes should use the same function; added flexibility should be managed with more parameters
  double detectionProbability(const double & apparentMagnitude) {
    return detectionProbability(apparentMagnitude,
				limitingMagnitude);
  }
 public:
  static double detectionProbability(const double & apparentMagnitude,
					   const double & limitingMagnitude) {
    
    double retVal = 0.5 + (limitingMagnitude - apparentMagnitude) * 0.45;
    
    if (retVal > 1) retVal = 1;
    if (retVal < 0) retVal = 0;
    
    return retVal;
  }
  
 public:
  const orsa::Time tStart;
  const orsa::Time tStop;
 public:
  const double limitingMagnitude;
  const double FOV_DEG;
  const double maxZenithDistanceAngle_DEG;
  const double minMoonDistanceAngle_DEG;
  const double minMoonPhase_DEG;
 public:
  const int recycleTime_DAY;
  const int dutyCycle_SEC;
  const int dutyCycleMultiplicity;
 public:
  const std::string obscode;
 public:
  const std::string name;  
  
  // derived
 public:
  const orsa::Time   effectiveDutyCycle;
  const orsa::Time   recycleTime;
  const double FOV;
  // const double cos_FOV;
  const double     effectiveHalfFOV;
  const double cos_effectiveHalfFOV;
  const double maxZenithDistanceAngle;
  const double cos_maxZenithDistanceAngle;
  const double cos_minMoonDistanceAngle;
  
 protected:
  const orsa::Matrix eclipticToEquatorial;
  const orsa::Matrix equatorialToEcliptic;
};


class SingleModeTelescope : public Telescope {
 public:
  SingleModeTelescope(const orsa::Time         & start,
		      const orsa::Time         & stop,
		      const int                  randomSeed,
		      const double       & limitingMagnitude,
		      const double       & FOV_DEG,
		      const double       & maxZenithDistanceAngle_DEG,
		      const double       & minMoonDistanceAngle_DEG,
		      const double       & minMoonPhase_DEG,
		      const int                  recycleTime_DAY,
		      const int                  dutyCycle_SEC,
		      const int                  dutyCycleMultiplicity,
		      const std::string        & obscode,
		      const std::string        & name) : 
    Telescope(start,
	      stop,
	      limitingMagnitude,
	      FOV_DEG,
	      maxZenithDistanceAngle_DEG,
	      minMoonDistanceAngle_DEG,
	      minMoonPhase_DEG,
	      recycleTime_DAY,
	      dutyCycle_SEC,
	      dutyCycleMultiplicity,
	      obscode,
	      name) { 
    rnd = new orsa::RNG(randomSeed);
  }
    
 public:  
  void writeCheckpoint() const { 
    
    FILE * fp_rnd = fopen(rngFileName().c_str(),"wb");
    if (fp_rnd != 0) {
      rnd->gsl_rng_fwrite(fp_rnd);
      fclose(fp_rnd);
    } else {
      ORSA_ERROR("cannot open file %s",
		 rngFileName().c_str());  
      boinc_finish(1);
    }
    
    Telescope::writeCheckpoint();
  }
  
 public:  
  void readCheckpoint(const orsa::Time & t) {
    
    // ORSA_DEBUG("------------ [%s] ----------",name.c_str());
    
    FILE * fp_rnd = fopen(rngFileName().c_str(),"rb");
    if (fp_rnd != 0) {
      rnd->gsl_rng_fread(fp_rnd);
      fclose(fp_rnd);
    } else {
      ORSA_ERROR("cannot open file %s",
		 rngFileName().c_str());  
      boinc_finish(1);
    }
    
    Telescope::readCheckpoint(t);
  }
  
 protected:
  osg::ref_ptr<orsa::RNG> rnd;
};


class OppositionTelescope : public SingleModeTelescope {
 public:
  OppositionTelescope(const orsa::Time         & start,
		      const orsa::Time         & stop,
		      const int                  randomSeed,
		      const double       & limitingMagnitude,
		      const double       & FOV_DEG,
		      const double       & maxZenithDistanceAngle_DEG,
		      const double       & minMoonDistanceAngle_DEG,
		      const double       & minMoonPhase_DEG,
		      const int                  recycleTime_DAY,
		      const int                  dutyCycle_SEC,
		      const int                  dutyCycleMultiplicity,
		      const std::string        & obscode,
		      const std::string        & name) : 
    SingleModeTelescope(start,
			stop,
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
			name) { }
 public:
  bool sampleTP(TelescopePointing   & tp,
		const orsa::Time    & t,
		const TelescopeList & telescopeList,
		const NEOList       & syntheticNEO,
		const double  & detectionProbabilityThreshold,
		const orsa::Vector  &  sunPosition,
		const orsa::Vector  & moonPosition,
		const orsa::Vector  &  obsPosition,
		const orsa::Vector  &  obsNormal,
		const bool            cacheON);
};


class SieveTelescope : public SingleModeTelescope {
 public:
  SieveTelescope(const orsa::Time         & start,
		 const orsa::Time         & stop,
		 const int                  randomSeed,
		 const double       & limitingMagnitude,
		 const double       & FOV_DEG,
		 const double       & maxZenithDistanceAngle_DEG,
		 const double       & minMoonDistanceAngle_DEG,
		 const double       & minMoonPhase_DEG,
		 const int                  recycleTime_DAY,
		 const int                  dutyCycle_SEC,
		 const int                  dutyCycleMultiplicity,
		 const std::string        & obscode,
		 const std::string        & name) : 
    SingleModeTelescope(start,
			stop,
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
			name) { }
 public:
  bool sampleTP(TelescopePointing   & tp,
		const orsa::Time    & t,
		const TelescopeList & telescopeList,
		const NEOList       & syntheticNEO,
		const double  & detectionProbabilityThreshold,
		const orsa::Vector  &  sunPosition,
		const orsa::Vector  & moonPosition,
		const orsa::Vector  &  obsPosition,
		const orsa::Vector  &  obsNormal,
		const bool            cacheON);
};


class AllSkyTelescope : public SingleModeTelescope {
 public:
  AllSkyTelescope(const orsa::Time         & start,
		  const orsa::Time         & stop,
		  const int                  randomSeed,
		  const double       & limitingMagnitude,
		  const double       & FOV_DEG,
		  const double       & maxZenithDistanceAngle_DEG,
		  const double       & minMoonDistanceAngle_DEG,
		  const double       & minMoonPhase_DEG,
		  const int                  recycleTime_DAY,
		  const int                  dutyCycle_SEC,
		  const int                  dutyCycleMultiplicity,
		  const std::string        & obscode,
		  const std::string        & name) : 
    SingleModeTelescope(start,
			stop,
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
			name) { }
 public:
  bool sampleTP(TelescopePointing   & tp,
		const orsa::Time    & t,
		const TelescopeList & telescopeList,
		const NEOList       & syntheticNEO,
		const double  & detectionProbabilityThreshold,
		const orsa::Vector  &  sunPosition,
		const orsa::Vector  & moonPosition,
		const orsa::Vector  &  obsPosition,
		const orsa::Vector  &  obsNormal,
		const bool            cacheON);
};


bool readTelescope(TelescopeList     & tl,
		   const std::string & fileName);


#endif // _SS_TELESCOPE_H_
