#ifndef SURVEY_REVIEW_SKY_COVERAGE
#define SURVEY_REVIEW_SKY_COVERAGE

#include <orsa/angle.h>
#include <orsa/cache.h>
#include <orsa/double.h>
#include <orsa/datetime.h>
#include <orsa/vector.h>

#include <orsaInputOutput/file.h>
#include <orsaInputOutput/MPC_asteroid.h>
#include <orsaInputOutput/MPC_observations.h>
#include <orsaInputOutput/MPC_obscode.h>

#include <orsaSolarSystem/observatory.h>
#include <orsaSolarSystem/orbit.h>

#include <osg/Referenced>

#include <string>
#include <list>
#include <libgen.h>

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
  orsa::Cache<orsa::Time>  epoch;   // local midnight epoch
  orsa::Cache<std::string> obscode;
  
 public:
  double eta(const double & V,
	     const double & U) const;
  /* inline double eta(const double & V,
     const double & U) const {
     return eta(V,
     V_limit.getRef(),
     eta0_V.getRef(),
     V0.getRef(),
     c_V.getRef(),
     w_V.getRef(),
     U,
     U_limit.getRef(),
     w_U.getRef(),
     beta.getRef());
     }
  */
  
 public:
  static double eta(const double & V,
		    const double & V_limit,
		    const double & eta0_V,
		    const double & V0,
		    const double & c_V,
		    const double & w_V,
		    const double & U,
		    const double & U_limit,
		    const double & w_U,
		    const double & beta);
  /* static inline double eta(const double & V,
     const double & V_limit,
     const double & eta0_V,
     const double & V0,
     const double & c_V,
     const double & w_V,
     const double & U,
     const double & U_limit,
     const double & w_U,
     const double & beta) {
     double retVal;
     if (V<V0) {
     retVal = eta0_V;
     } else {
     retVal = 
     (eta0_V-c_V*orsa::square(V-V0)) / 
     (1.0+exp( cos(beta)*(V-V_limit)/w_V + sin(beta)*(U_limit-U)/w_U)) / 
     (1.0+exp(-sin(beta)*(V-V_limit)/w_V + cos(beta)*(U_limit-U)/w_U));
     }
     if (retVal < 0.0) retVal=0.0;
     if (retVal > 1.0) retVal=1.0;
     return retVal;
     }
  */
  
 public:
  // coefficients for efficiency as function of apparent magnitude V
  orsa::Cache<double> V_limit, eta0_V, V0, c_V, w_V;
  // coefficients for efficiency as function of apparent velocity U
  orsa::Cache<double> U_limit, w_U;
  // mixing angle
  orsa::Cache<double> beta;
 public:
  // return filename, stripping path and suffix (after first dot)
  static std::string basename(const std::string & filename);
 public:
  static bool processFilename(const std::string & filename_in,
			      orsaInputOutput::MPCObsCodeFile * obsCodeFile,
			      std::string & obsCode,
			      orsa::Time & epoch,
			      int & year,
			      int & dayOfYear);    
};

#endif // SURVEY_REVIEW_SKY_COVERAGE
