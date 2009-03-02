#ifndef _ORSA_COVARIANCE_
#define _ORSA_COVARIANCE_

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/cache.h>
#include <orsa/double.h>

#include <vector>

namespace orsa {
  
  // TO BE COMPLETED...
  
  class Covariance : public osg::Referenced {
  public:	
    Covariance();
  protected:
    virtual ~Covariance();
    
  private:
    std::vector<std::vector<orsa::Cache<orsa::Double> > > _covm;
  };
  
} // namespace orsa

#endif // _ORSA_COVARIANCE_
