#ifndef SURVEY_REVIEW_SKY_COVERAGE
#define SURVEY_REVIEW_SKY_COVERAGE

#include <orsa/angle.h>
#include <orsa/cache.h>
#include <orsa/double.h>
#include <orsa/datetime.h>
#include <orsa/vector.h>

#include <orsaInputOutput/file.h>

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
  bool get(const orsa::Vector u,
	   double & V) const;
  
 public:
  bool get(const orsa::Angle & ra,
	   const orsa::Angle & dec,
	   double & V) const {
    return get(unitVector(ra,dec),V);
  }
  
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
};

#endif // SURVEY_REVIEW_SKY_COVERAGE
