#ifndef _ORSA_PAUL_MOMENT_H_
#define _ORSA_PAUL_MOMENT_H_

#include <vector>
#include <string>

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/double.h>
// #include <orsa/shape.h>
#include <orsa/vector.h>
// #include <orsa/massDistribution.h>
// #include <orsa/matrix.h>

namespace orsa {
  
  class PaulMoment : public osg::Referenced {
  public:
    PaulMoment(const unsigned int order);
  protected:
    virtual ~PaulMoment() { } 
    
  public:
    double M (const int i, const int j, const int k) const;
    double M_uncertainty (const int i, const int j, const int k) const;
  public:
    void setM (const double & val, const int i, const int j, const int k);
    void setM_uncertainty (const double & val, const int i, const int j, const int k);
  protected:
    std::vector< std::vector< std::vector<double> > > _M, _M_uncertainty;
  public:
    const unsigned int order;
  };
  
  // utility, just printing out values for now
  void convert(std::vector< std::vector<double> > & C,
	       std::vector< std::vector<double> > & S,
	       std::vector< std::vector<double> > & norm_C,
	       std::vector< std::vector<double> > & norm_S,
	       std::vector<double> & J,
	       const PaulMoment * const pm,
	       const double     & R0,
	       const bool         verbose=false);
  
  // tries to find a set of PaulMoments that matches the normalized C and S
  // the order of pm must be <= the order of norm_C and norm_S
  bool solve(PaulMoment * pm,
	     const std::vector< std::vector<double> > & norm_C,
	     const std::vector< std::vector<double> > & norm_S,
	     const double     & R0);
  
  // a,b,c are the 3 semiaxes along x,y,z respectively
  void EllipsoidExpansion(PaulMoment   * pm,
			  const double & a,
			  const double & b,
			  const double & c);
  
}; // namespace orsa

#endif // _ORSA_PAUL_MOMENT_H_
