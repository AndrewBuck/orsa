#ifndef SURVEY_REVIEW_SKY_COVERAGE
#define SURVEY_REVIEW_SKY_COVERAGE

#include <orsa/angle.h>
#include <orsa/cache.h>
#include <orsa/double.h>
#include <orsa/datetime.h>
#include <orsa/vector.h>

#include <orsaInputOutput/file.h>

#include <orsaSolarSystem/orbit.h>

#include <osg/Referenced>

#include <string>
#include <list>

// single night, single observatory
class SkyCoverage : public osg::Referenced {
 public:
  SkyCoverage();
 protected:
  virtual ~SkyCoverage();
  
  //! Unit Vectors are in Ecliptic coordinates
  
 public:
  static orsa::Vector unitVector(const orsa::Angle & ra,
				 const orsa::Angle & dec);  
 public:
  static void normalize(orsa::Angle & ra,
			orsa::Angle & dec);  
 public:
  static void normalize(double & ra,
			double & dec);  
 public:
  //! translates the observatory code used in sky-coverage file to a standard MPC observatory code
  static std::string alias(const std::string &);
  
 public:
  void reset();
  
 public:
  bool setField(const double & x1,
		const double & y1,
		const double & x2,
		const double & y2,
		const double & x3,
		const double & y3,
		const double & x4,
		const double & y4,
		const double & V);
  
 public:
  bool get(const orsa::Vector & u,
	   double & V,
	   const bool verbose=false) const;
  
 public:
  bool get(const orsa::Angle & ra,
	   const orsa::Angle & dec,
	   double & V,
	   const bool verbose=false) const {
    return get(unitVector(ra,dec),V,verbose);
  }
  
 public:
  // "fast" version, only to know if the given direction is covered
  bool fastGet(const orsa::Vector & u) const;
  bool fastGet(const orsa::Angle & ra,
	       const orsa::Angle & dec) const {
    return fastGet(unitVector(ra,dec));
  }
  
 public:
  //! smallest distance (angle) between the u direction and any of the fields centers
  double minDistance(const orsa::Vector & u,
		     const bool verbose=false) const;
  
 public:
  class SkyCoverageElement {
  public:
    orsa::Vector u_centerField;
    orsa::Vector u_RA;
    orsa::Vector u_DEC;
    double halfFieldSize_RA;
    double halfFieldSize_DEC;
    double minScalarProduct;
    double limitingMagnitude;
  };
  
 protected:
  std::list<SkyCoverageElement> data;
  
 public:
  orsa::Cache<orsa::Time>  epoch;               // local midnight epoch
  // orsa::Cache<orsa::Time>  halfObservingPeriod; // 1/2 of the duration of the observing night
  // orsa::Cache<orsa::Time>  nightStart, nightStop;  // epoch -/+ 1/2 observing period
  orsa::Cache<std::string> obscode;
  
  inline double eta_V(const double & V,
		      const double & V_limit) const {
    return eta_V(V,V_limit,eta0.getRef(),c.getRef(),V0.getRef(),w.getRef());
  }
  
 public:
  // V = apparent magnitude
  static inline double eta_V(const double & V,
			     const double & V_limit,
			     const double & eta0_V,
			     const double & c_V,
			     const double & V0,
			     const double & w_V) {
    if (V<V0) return eta0_V;
    double retVal = (eta0_V-c_V*orsa::square(V-V0))/(1+exp((V-V_limit)/w_V));
    // if (retVal < orsa::epsilon()) retVal=0.0;
    // if (retVal > 1.0) retVal=1.0;
    return retVal;
  }
  
 public:
  // U = apparent velocity
  static inline double eta_U(const double & U,
			     const double & U_limit_slow,
			     const double & w_U_slow,
			     const double & U0) {
    if (U>U0) return 1.0;
    double retVal = 1.0/(1+exp((U_limit_slow-U)/w_U_slow));
    // if (retVal < orsa::epsilon()) retVal=0.0;
    // if (retVal > 1.0) retVal=1.0;
    return retVal;
  }
  
 public:
  // coefficients for efficiency as function of apparent magnitude
  // should include a function of the galactic latitude <-> background stars number density
  orsa::Cache<double> eta0, c, V0, w;
  
 public:
  // return all the filename up to the first dot
  static std::string basename(const std::string & filename) {
    size_t found_dot = std::string(filename).find(".",0);
    if (found_dot == std::string::npos) {
      ORSA_DEBUG("not regular filename: %s",filename.c_str());
      exit(0);
    }
    std::string s;
    s.assign(filename,0,found_dot);
    return s;
  }
};

#endif // SURVEY_REVIEW_SKY_COVERAGE
