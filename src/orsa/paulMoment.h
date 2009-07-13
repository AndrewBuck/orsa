#ifndef _ORSA_PAUL_MOMENT_H_
#define _ORSA_PAUL_MOMENT_H_

#include <vector>
#include <string>

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/double.h>
#include <orsa/shape.h>
#include <orsa/vector.h>
#include <orsa/massDistribution.h>
#include <orsa/matrix.h>

namespace orsa {
  
  class PaulMoment : public osg::Referenced {
  public:
    PaulMoment();
  public:
    PaulMoment(const int order);
  protected:
    virtual ~PaulMoment() { } 
    
  public:
    bool computeUsingShape(orsa::Vector & centerOfMass,
			   orsa::Matrix & principalAxis,
			   const orsa::MassDistribution * massDistribution,
			   const unsigned int sample_points,
			   const int random_seed);
    
  private:
    void _randomVectorInside(std::vector<bool> &,
			     const unsigned int sample_points,
			     const int random_seed);
  private:
    void _centerOfMass(Vector & center_of_mass,
		       Vector & center_of_mass_uncertainty,
		       const std::vector<bool> & in,
		       const unsigned int sample_points,
		       const int random_seed);
  private:    
    void _inertiaMoment(Matrix & inertia_moment,
			Matrix & inertia_moment_uncertainty,
			const Vector & center_of_mass,
			const Vector & center_of_mass_uncertainty,
			const std::vector<bool> & in,
			const unsigned int sample_points,
			const int random_seed);
  private:
    void _moment(std::vector< std::vector< std::vector< double > > > & M,
		 std::vector< std::vector< std::vector< double > > > & M_uncertainty,
		 const orsa::Vector                                        & center_of_mass,
		 const orsa::Vector                                        & center_of_mass_uncertainty,
		 const std::vector<bool>                                   & in,
		 const unsigned int                                          sample_points,
		 const int                                                   random_seed);
    
  public:
    const double M (const int i,
		    const int j, 
		    const int k) const {
      if (i<0) { ORSA_ERROR("negative index: i=%i",i); return 0; }
      if (j<0) { ORSA_ERROR("negative index: j=%i",j); return 0; }
      if (k<0) { ORSA_ERROR("negative index: k=%i",k); return 0; }
      //
      if (i+j+k > _order.getRef()) { ORSA_ERROR("index out of bound"); return 0; }
      //
      return _M[i][j][k];
    }
    
  public:
    void setM (const double & val,
	       const int i, 
	       const int j, 
	       const int k) {
      if (i<0) { ORSA_ERROR("negative index: i=%i",i); return; }
      if (j<0) { ORSA_ERROR("negative index: j=%i",j); return; }
      if (k<0) { ORSA_ERROR("negative index: k=%i",k); return; }
      //
      if (i+j+k > _order.getRef()) { ORSA_ERROR("index out of bound"); return; }
      //
      _M[i][j][k] = val;
    }
    
  protected:
    std::vector< std::vector< std::vector<double> > > _M, _M_uncertainty;
    
  protected:
    orsa::Cache<unsigned int> _sample_points;
    
    /* public:
       bool setOrder(const int newOrder) {
       if (newOrder <= order()) {
       _order = newOrder;
       return true;
       } else {
       return false;
       }
       }
    */
  public:
    const int order() const {
      return _order.getRef();
    }
  protected:
    orsa::Cache<int> _order;
  };
  
  
  
  
  
  
  
  
  
  // utility, just printing out values for now
  void convert(const PaulMoment * const pm,
	       const double     & R0);
  
}; // namespace orsa

#endif // _ORSA_PAUL_MOMENT_H_
