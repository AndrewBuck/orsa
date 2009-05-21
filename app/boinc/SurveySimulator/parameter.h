#ifndef _DUST_PARAMETER_H_
#define _DUST_PARAMETER_H_

#include <string>

#include <orsa/cache.h>
#include <orsa/double.h>

#include <osg/Referenced>
#include <osg/ref_ptr>

typedef orsa::Cache<double>      ParameterD;
typedef orsa::Cache<int>         ParameterI;
typedef orsa::Cache<std::string> ParameterS;

class ParameterFile : public osg::Referenced {
 public:
  ParameterFile() : osg::Referenced() { }
 protected:    
  virtual ~ParameterFile() { }
  
 public:
  bool mainRead(const std::string & fileName);
  
 public:
  void print() const;
  
 private:
  static bool read_D  (ParameterD  &, FILE *);
  static bool read_I  (ParameterI  &, FILE *);
  static bool read_S  (ParameterS  &, FILE *);
  
 public:
  // ParameterI numRealNEO;
  // ParameterI numSyntheticNEO;
  // ParameterI randomSeedRealNEO;
  // ParameterI randomSeedSyntheticNEO;
  // ParameterI randomSeedSky;
  ParameterD orbit_a_AU_min;
  ParameterD orbit_a_AU_max;
  ParameterD orbit_e_min;
  ParameterD orbit_e_max;
  ParameterD orbit_i_DEG_min;
  ParameterD orbit_i_DEG_max;
  ParameterD Hmax;
  ParameterD Ha;
  // ParameterD limitingMagnitude;
  // ParameterD FOV_DEG;
  // ParameterD maxZenithDistanceAngle_DEG;
  // ParameterD minMoonDistanceAngle_DEG;
  // ParameterD minMoonPhaseAngle_DEG;
  ParameterD detectionProbabilityThreshold;
  // ParameterI recycleTime_DAY;
  // ParameterI dutyCycle_SEC;
  // ParameterI dutyCycleMultiplicity;
  ParameterI cacheNEOVector;
};

#endif // _DUST_PARAMETER_H_
