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
  // total area of sky covered, as declared in file, no check for overlaps is performed
  double totalDegSq() const;
  
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

class SkyCoverageFile : 
public orsaInputOutput::InputFile < orsaInputOutput::PlainFile, osg::ref_ptr<SkyCoverage> > {
 public:
  SkyCoverageFile() : 
    orsaInputOutput::InputFile < orsaInputOutput::PlainFile, osg::ref_ptr<SkyCoverage> > () {
    _data = new SkyCoverage;
  }
 protected:
  ~SkyCoverageFile() { }
 public:
  void setFileName (const std::string & filename) {
    // first, read obscode and epoch from filename
    osg::ref_ptr<orsaInputOutput::MPCObsCodeFile> obsCodeFile = new orsaInputOutput::MPCObsCodeFile;
    obsCodeFile->setFileName("obscode.dat");
    obsCodeFile->read();
    std::string obsCode;
    orsa::Time epoch;
    int year;
    int dayOfYear;
    if (!SkyCoverage::processFilename(filename,
				      obsCodeFile.get(),
				      obsCode,
				      epoch,
				      year,
				      dayOfYear)) {
      ORSA_DEBUG("problems...");
    }
    _data->obscode=obsCode;
    _data->epoch=epoch;
    // now, regular call
    orsaInputOutput::InputFile < orsaInputOutput::PlainFile, osg::ref_ptr<SkyCoverage> >::setFileName(filename);
  }
 public:
  bool processLine(const char * line) {
    double x1,x2,x3,x4; // ra
    double y1,y2,y3,y4; // dec
    double V;           // limiting mag 
    if (9 == gmp_sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
			&x1,&y1,
			&x2,&y2,
			&x3,&y3,
			&x4,&y4,
			&V)) {
      SkyCoverage::normalize(x1,y1);
      SkyCoverage::normalize(x2,y2);
      SkyCoverage::normalize(x3,y3);
      SkyCoverage::normalize(x4,y4);
      //
      _data->setField(x1,y1,x2,y2,x3,y3,x4,y4,V);
      return true;
    } else {
      return false;
    }
  }
 public:
  bool goodLine(const char * /* line */ ) { 
    return true;
  }
};

#endif // SURVEY_REVIEW_SKY_COVERAGE
