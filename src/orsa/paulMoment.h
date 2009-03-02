#ifndef _ORSA_PAUL_MOMENT_H_
#define _ORSA_PAUL_MOMENT_H_

#include <vector>
#include <string>

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/double.h>
#include <orsa/shape.h>
#include <orsa/vector.h>
#include <orsa/matrix.h>
#include <orsa/multipole.h>

namespace orsa {
  
  // taken from orsa::Multipole
  
  class PaulMoment : public osg::Referenced {
  public:
    PaulMoment();
  public:
    PaulMoment(const unsigned int order);
  protected:
    virtual ~PaulMoment();
    
  public:
    bool setShape(const orsa::Shape * s) {
      _shape = s;
      return true;
    }
  public:
    bool setShape(const orsa::Shape * s, const unsigned int order) {
      _shape = s;
      _order.set(order);
      return true;
    }
  public:
    const orsa::Shape * getShape() const {
      return _shape.get();
    }
  protected:
    osg::ref_ptr<const orsa::Shape> _shape;
    
  public:
    bool setMassDistribution(const orsa::MassDistribution * md) {
      if (md == 0) return false;
      _massDistribution = md;
      return true;
    }
  public:
    const orsa::MassDistribution * getMassDistribution() const {
      return _massDistribution.get();
    }
  protected:
    osg::ref_ptr<const orsa::MassDistribution> _massDistribution;
    
  public:
    bool computeUsingShape(const unsigned int sample_points,
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
    
    /* 
       public:
       bool writeToFile(const std::string & filename) const;
       bool readFromFile(const std::string & filename);
       
       public:
       bool readFromAltFile(const std::string & filename);
       
       public:
       bool readFromBisFile(const std::string & filename);
       
       protected:
       bool _write_in_vector_to_file(const std::vector<bool> & in,
       const std::string & filename) const;
       bool _read_in_vector_from_file(std::vector<bool> & in,
       const std::string & filename);
    */
    
  public:
    const double M (const unsigned int i,
		    const unsigned int j, 
		    const unsigned int k) const {
      if (i+j+k <= _order.getRef()) {
	return _M[i][j][k];
      } else {
	return 0;
      }
    }
  public:
    void setM (const double & val,
	       const unsigned int i, 
	       const unsigned int j, 
	       const unsigned int k) {
      if (i+j+k <= _order.getRef()) {
	_M[i][j][k] = val;
      } else {
        ORSA_ERROR("...");
      }
    }
  protected:
    std::vector< std::vector< std::vector<double> > > _M, _M_uncertainty;
    
  protected:
    orsa::Cache<unsigned int> _sample_points;
    
  public:
    void setCenterOfMass(const orsa::Vector & v) {
      _center_of_mass = v;
    }
  public:
    const orsa::Vector & getCenterOfMass() const {
      return _center_of_mass.getRef();
    }
  protected:
    orsa::Cache<orsa::Vector> _center_of_mass, _center_of_mass_uncertainty;
    
  public:
    void setInertiaMoment(const orsa::Matrix & I) {
      _inertia_moment = I;
    }
  public:
    const orsa::Matrix & getInertiaMoment() const {
      return _inertia_moment.getRef();
    }
  protected:
    orsa::Cache<orsa::Matrix> _inertia_moment, _inertia_moment_uncertainty;
    
  public:
    // the order can only be reduced...
    bool setOrder(const unsigned int newOrder) {
      if (newOrder <= order()) {
	_order = newOrder;
	return true;
      } else {
	return false;
      }
    }
  public:
    const unsigned int order() const {
      return _order.getRef();
    }
  protected:
    orsa::Cache<unsigned int> _order;
  };
  
}; // namespace orsa

#endif // _ORSA_PAUL_MOMENT_H_
