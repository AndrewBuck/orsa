#include "parameter.h"

// BOINC API
#include <boinc_api.h>
#include <filesys.h>

bool ParameterFile::read_D(ParameterD & pD, FILE * fp) {
  char line[1024];
  bool skipLine;
  do {
    if (!fgets(line,1024,fp)) { 
      ORSA_DEBUG("problems...");
      return false; 
    }
    skipLine = false;
    if (strlen(line) >= 1) {
      if (line[0] == '#') {
	skipLine = true;
      }
    }
  } while (skipLine);
  orsa::Double tmp_D;
  gmp_sscanf(line,"%Ff",tmp_D.get_mpf_t());
  pD = tmp_D;
  // ORSA_DEBUG("%Ff",tmp_D.get_mpf_t());
  return true;
}

bool ParameterFile::read_I(ParameterI & pI, FILE * fp) {
  char line[1024];
  bool skipLine;
  do {
    if (!fgets(line,1024,fp)) { 
      ORSA_DEBUG("problems...");
      return false; 
    }
    skipLine = false;
    if (strlen(line) >= 1) {
      if (line[0] == '#') {
	skipLine = true;
      }
    }
  } while (skipLine);
  int tmp_I;
  gmp_sscanf(line,"%d",&tmp_I);
  pI = tmp_I;
  return true;
}

bool ParameterFile::read_S(ParameterS & pS, FILE * fp) {
  char line[1024];
  bool skipLine;
  do {
    if (!fgets(line,1024,fp)) { 
      ORSA_DEBUG("problems...");
      return false; 
    }
    skipLine = false;
    if (strlen(line) >= 1) {
      if (line[0] == '#') {
	skipLine = true;
      }
    }
  } while (skipLine);
  char tmp_s[1024];
  gmp_sscanf(line,"%s",tmp_s);
  pS = tmp_s;
  return true;
}

bool ParameterFile::mainRead(const std::string & fileName) {
  
  FILE * fp = boinc_fopen(fileName.c_str(),"r");
  
  if (!fp) {
    return false;
  }
  
  // read_I(numRealNEO,fp);
  // read_I(numSyntheticNEO,fp);
  // read_I(randomSeedRealNEO,fp);
  // read_I(randomSeedSyntheticNEO,fp);
  // read_I(randomSeedSky,fp);
  read_D(orbit_a_AU_min,fp);
  read_D(orbit_a_AU_max,fp);
  read_D(orbit_e_min,fp);
  read_D(orbit_e_max,fp);
  read_D(orbit_i_DEG_min,fp);
  read_D(orbit_i_DEG_max,fp);
  read_D(Hmax,fp);
  read_D(Ha,fp);
  // read_D(limitingMagnitude,fp);
  // read_D(FOV_DEG,fp);
  // read_D(maxZenithDistanceAngle_DEG,fp);
  // read_D(minMoonDistanceAngle_DEG,fp);
  // read_D(minMoonPhaseAngle_DEG,fp);
  read_D(detectionProbabilityThreshold,fp);
  // read_I(recycleTime_DAY,fp);
  // read_I(dutyCycle_SEC,fp);
  // read_I(dutyCycleMultiplicity,fp);
  read_I(cacheNEOVector,fp);
  
  fclose(fp);
  
  return true;
}

void ParameterFile::print() const {
  
  gmp_fprintf(stderr,
	      // "numRealNEO...................: %10i\n"
	      // "numSyntheticNEO..............: %10i\n"
	      // "randomSeedRealNEO............: %10i\n"
	      // "randomSeedSyntheticNEO.......: %10i\n"
	      // "randomSeedSky................: %10i\n"
	      "orbit_a_AU_min...............: %20.9Ff\n"
	      "orbit_a_AU_max...............: %20.9Ff\n"
	      "orbit_e_min..................: %20.9Ff\n"
	      "orbit_e_max..................: %20.9Ff\n"
	      "orbit_i_DEG_min..............: %20.9Ff\n"
	      "orbit_i_DEG_max..............: %20.9Ff\n"
	      "Hmax.........................: %20.9Ff\n"
	      "Ha...........................: %20.9Ff\n"
	      // "limitingMagnitude............: %20.9Ff\n"
	      // "FOV_DEG......................: %20.9Ff\n"
	      // "maxZenithDistanceAngle_DEG...: %20.9Ff\n"
	      // "minMoonDistanceAngle_DEG.....: %20.9Ff\n"
	      // "minMoonPhaseAngle_DEG........: %20.9Ff\n"
	      "detectionProbabilityThreshold: %20.9Ff\n"
	      // "recycleTime_DAY..............: %10i\n"
	      // "dutyCycle_SEC................: %10i\n"
	      // "dutyCycleMultiplicity........: %10i\n"
	      "cacheNEOVector...............: %10i\n"
	      ,
	      // numRealNEO.getRef(),
	      // numSyntheticNEO.getRef(),
	      // randomSeedRealNEO.getRef(),
	      // randomSeedSyntheticNEO.getRef(),
	      // randomSeedSky.getRef(),
	      orbit_a_AU_min.getRef().get_mpf_t(),
	      orbit_a_AU_max.getRef().get_mpf_t(),
	      orbit_e_min.getRef().get_mpf_t(),
	      orbit_e_max.getRef().get_mpf_t(),
	      orbit_i_DEG_min.getRef().get_mpf_t(),
	      orbit_i_DEG_max.getRef().get_mpf_t(),
	      Hmax.getRef().get_mpf_t(),
	      Ha.getRef().get_mpf_t(),
	      // limitingMagnitude.getRef().get_mpf_t(),
	      // FOV_DEG.getRef().get_mpf_t(),
	      // maxZenithDistanceAngle_DEG.getRef().get_mpf_t(),
	      // minMoonDistanceAngle_DEG.getRef().get_mpf_t(),
	      // minMoonPhaseAngle_DEG.getRef().get_mpf_t(),
	      detectionProbabilityThreshold.getRef().get_mpf_t(),
	      // recycleTime_DAY.getRef(),
	      // dutyCycle_SEC.getRef(),
	      // dutyCycleMultiplicity.getRef(),
	      cacheNEOVector.getRef());
}
