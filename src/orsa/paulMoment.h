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
    // PaulMoment();
  public:
    PaulMoment(const unsigned int order);
  protected:
    virtual ~PaulMoment() { } 
    
  public:
    const double M (const int i, const int j, const int k) const;
    const double M_uncertainty (const int i, const int j, const int k) const;
  public:
    void setM (const double & val, const int i, const int j, const int k);
    void setM_uncertainty (const double & val, const int i, const int j, const int k);
  protected:
    std::vector< std::vector< std::vector<double> > > _M, _M_uncertainty;
  public:
    const unsigned int order;
  };
  
  // utility, just printing out values for now
  void convert(const PaulMoment * const pm,
	       const double     & R0);
  
}; // namespace orsa

#endif // _ORSA_PAUL_MOMENT_H_
