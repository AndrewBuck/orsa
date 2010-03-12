#include "skycoverage.h"

#include <orsa/double.h>

#include <orsaSolarSystem/obleq.h>

SkyCoverage::SkyCoverage() : osg::Referenced() { }

SkyCoverage::~SkyCoverage() { }

orsa::Vector SkyCoverage::unitVector(const orsa::Angle & ra,
				     const orsa::Angle & dec) {
  double s_ra, c_ra;
  sincos(ra.getRad(),
	 &s_ra,
	 &c_ra);
  double s_dec, c_dec;
  sincos(dec.getRad(),
	 &s_dec, 
	 &c_dec);
  return (orsaSolarSystem::equatorialToEcliptic() * 
	  orsa::Vector(c_dec*c_ra,
		       c_dec*s_ra,
		       s_dec));
}

void SkyCoverage::normalize(orsa::Angle & ra,
			    orsa::Angle & dec) {
  double  ra_rad =  ra.getRad();
  double dec_rad = dec.getRad();
  //
  SkyCoverage::normalize(ra_rad,dec_rad);
  //
  ra.setRad(ra_rad);
  dec.setRad(dec_rad);
}

void SkyCoverage::normalize(double & ra_rad,
			    double & dec_rad) {
  if (dec_rad > orsa::halfpi()) {
    dec_rad = orsa::pi() - dec_rad;
    ra_rad += orsa::pi();
  }
  if (dec_rad < -orsa::halfpi()) {
    dec_rad = -dec_rad - orsa::pi();
    ra_rad += orsa::pi();
  }
  if (fabs(dec_rad) > orsa::halfpi()) {
    ORSA_DEBUG("still problems after normalization...");
  }
  // ra after dec because depends on it, if close to halfpi
  ra_rad = fmod(fmod(ra_rad,orsa::twopi())+orsa::twopi(),orsa::twopi());
}

std::string SkyCoverage::alias(const std::string & fileCode) {
  if (fileCode=="CSS")        return "703";
  if (fileCode=="LINEAR")     return "704";
  if (fileCode=="SPACEWATCH") return "291";
  return fileCode;
} 

void SkyCoverage::reset() {
  data.clear();
}

bool SkyCoverage::setField(const double & x1,
			   const double & y1,
			   const double & x2,
			   const double & y2,
			   const double & x3,
			   const double & y3,
			   const double & x4,
			   const double & y4,
			   const double & V) {
  
  // check if data is in standard order and format [should check more]
  if ( (y1 != y2) || 
       (y3 != y4) ) {
    ORSA_DEBUG("data not in standard format");
    return false;
  }
  
  const orsa::Vector u1 = unitVector(x1,y1);
  const orsa::Vector u2 = unitVector(x2,y2);
  const orsa::Vector u3 = unitVector(x3,y3);
  const orsa::Vector u4 = unitVector(x4,y4);
  
  const double field_RA  = acos(u1*u2);
  const double field_DEC = acos(u2*u3);
  
  ORSA_DEBUG("field: %f x %f [deg^2] [RAxDEC]",
	     orsa::radToDeg()*field_RA,
	     orsa::radToDeg()*field_DEC);
  
  SkyCoverageElement e;
  
  // in Ecliptic coords...
  const orsa::Vector northEquatorialPole = orsaSolarSystem::equatorialToEcliptic()*orsa::Vector(0,0,1);
  
  e.u_centerField = (u1+u2+u3+u4).normalized();
  e.u_RA  = orsa::externalProduct(northEquatorialPole,e.u_centerField).normalized();
  e.u_DEC = orsa::externalProduct(e.u_centerField,e.u_RA).normalized();
  e.halfFieldSize_RA  = 0.5*field_RA;
  e.halfFieldSize_DEC = 0.5*field_DEC;
  e.minScalarProduct  = std::min(std::min(e.u_centerField*u1,
					  e.u_centerField*u2),
				 std::min(e.u_centerField*u3,
					  e.u_centerField*u4));
  e.limitingMagnitude = V;
  
  data.push_back(e);
  
  // ORSA_DEBUG("data size: %i   minScalarProduct: %f",data.size(),e.minScalarProduct);
  
  /* 
     {
     // some testing before leaving
     ORSA_DEBUG("consinstency checks");
     }
  */
  
  return true;
}

bool SkyCoverage::get(const orsa::Vector & u,
		      double & V,
		      const bool verbose) const {
  orsa::Cache<double> local_V;
  std::list<SkyCoverageElement>::const_iterator it = data.begin();
  while (it != data.end()) {
    if (verbose) {
      ORSA_DEBUG("acos(u*(*it).u_centerField): %7.3f [deg]",orsa::radToDeg()*acos(u*(*it).u_centerField));
      // orsa::print(u);
      // orsa::print((*it).u_centerField);
    }
    if (u*(*it).u_centerField > (*it).minScalarProduct) {
      const double delta_RA  = fabs(asin(u*(*it).u_RA));
      // if (verbose) ORSA_DEBUG("delta_RA: %f [deg]",orsa::radToDeg()*delta_RA);
      if (delta_RA < (*it).halfFieldSize_RA) {
	const double delta_DEC = fabs(asin(u*(*it).u_DEC));
	// if (verbose) ORSA_DEBUG("delta_DEC: %f [deg]",orsa::radToDeg()*delta_DEC);
	if (delta_DEC < (*it).halfFieldSize_DEC) {
	  // if (verbose)  ORSA_DEBUG("found in one field, V: %f",(*it).limitingMagnitude);
	  if (local_V.isSet()) {
	    local_V = std::max(local_V.getRef(),
			       (*it).limitingMagnitude);
	  } else {
	    local_V = (*it).limitingMagnitude;
	  }
	}
      }
    }
    ++it;
  }
  if (local_V.isSet()) {
    V = local_V.getRef();
    return true;
  } else {
    return false;
  }
}

bool SkyCoverage::fastGet(const orsa::Vector & u) const {
  orsa::Cache<double> local_V;
  std::list<SkyCoverageElement>::const_iterator it = data.begin();
  while (it != data.end()) {
    if (u*(*it).u_centerField > (*it).minScalarProduct) {
      if (fabs(asin(u*(*it).u_RA)) < (*it).halfFieldSize_RA) {
       	if (fabs(asin(u*(*it).u_DEC)) < (*it).halfFieldSize_DEC) {
	  return true;
	}
      }
    }
    ++it;
  }
  return false;
}

double SkyCoverage::minDistance(const orsa::Vector & u,
				const bool /* verbose */ ) const {
  double minArc=orsa::pi();
  std::list<SkyCoverageElement>::const_iterator it = data.begin();
  while (it != data.end()) {
    minArc=std::min(minArc,acos(u*(*it).u_centerField));
    ++it;
  }
  return minArc;
}

