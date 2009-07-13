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

namespace orsa {
  
  // taken from orsa::Multipole
  
  class MassDistribution : public osg::Referenced {
  public:
    MassDistribution() : osg::Referenced(true) { }
  protected:
    virtual ~MassDistribution() { }
  public:
    virtual double density(const Vector &) const = 0;
  };
  
  class UniformMassDistribution : public MassDistribution {
  public:
    UniformMassDistribution() : MassDistribution() { }
  protected:
    ~UniformMassDistribution() { }
  public:
    double density(const Vector &) const {
      return 1;
    }
  };
  
  class SphericalCorePlusMantleMassDistribution : public MassDistribution {
  public:
    SphericalCorePlusMantleMassDistribution(const double & coreRadius,
					    const double & coreDensity,
					    const double & mantleDensity) :
      MassDistribution(),
      R2(coreRadius*coreRadius),
      dC(coreDensity),
      dM(mantleDensity) { }
  protected:
    ~SphericalCorePlusMantleMassDistribution() { }
  public:
    double density(const Vector & v) const {
      if (v.lengthSquared() > R2) {
	return dM;
      } else {
	return dC;
      }
    }
  protected:
    const double R2, dC, dM;
  };
  
  class PaulMoment : public osg::Referenced {
  public:
    PaulMoment();
  public:
    PaulMoment(const int order);
  protected:
    virtual ~PaulMoment();
    
  public:
    bool setShape(const orsa::Shape * s) {
      _shape = s;
      return true;
    }
  public:
    bool setShape(const orsa::Shape * s, const int order) {
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
    bool setOrder(const int newOrder) {
#warning "this should trigger some set_dirty method..."
      if (newOrder <= order()) {
	_order = newOrder;
	return true;
      } else {
	return false;
      }
    }
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
