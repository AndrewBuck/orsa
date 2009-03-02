#ifndef _ORSA_MULTIPOLE_
#define _ORSA_MULTIPOLE_

#include <vector>
#include <string>

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/double.h>
#include <orsa/shape.h>
#include <orsa/vector.h>
#include <orsa/matrix.h>

namespace orsa {
  
  class MassDistribution : public osg::Referenced {
  public:
    MassDistribution() : osg::Referenced() { }
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
  
  class Multipole : public osg::Referenced {
  public:
    Multipole();
  public:
    Multipole(const unsigned int order);
  protected:
    virtual ~Multipole();
    
  public:
    bool setShape(const orsa::Shape * s) {
      _shape = s;
      // return compute();
      return true;
    }
  public:
    bool setShape(const orsa::Shape * s, const unsigned int order) {
      _shape = s;
      _order.set(order);
      // return compute();
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
      // return compute();
      return true;
    }
  public:
    /* 
       bool setMassDistribution(const orsa::MassDistribution * md, const unsigned int order) {
       if (md == 0) return false;
       _massDistribution = md;
       _order.set(order);
       // return compute();
       return true;
       }
    */
  public:
    const orsa::MassDistribution * getMassDistribution() const {
      return _massDistribution.get();
    }
  protected:
    osg::ref_ptr<const orsa::MassDistribution> _massDistribution;
    
  public:
    bool computeUsingShape(const unsigned int sample_points,
			   const int random_seed);
  public:
    bool computeUsingShape(const unsigned int sample_points,
			   const int random_seed,
			   const double & R0);
  public:
    void setR0(const double & R0) {
      _R0.set(R0);
    }
  public:
    const double & getR0() const {
      return _R0.getRef();
    }
  protected:
    orsa::Cache<double> _R0;
    
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
    void _multipole(std::vector<std::vector<double> > & C,
		    std::vector<std::vector<double> > & C_uncertainty,
		    std::vector<std::vector<double> > & S, 
		    std::vector<std::vector<double> > & S_uncertainty, 
		    const Vector & center_of_mass,
		    const Vector & center_of_mass_uncertainty,
		    const std::vector<bool> & in,
		    const unsigned int sample_points,
		    const int random_seed);
    
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
    
  public:
    const double & C (const unsigned int l, const unsigned int m) const {
      // should check for l and m range
      return _C[l][m];
    }
  public:
    const double & S (const unsigned int l, const unsigned int m) const {
      // should check for l and m range
      return _S[l][m];
    }
  public:
    void setC (const double & val,
	       const unsigned int l, const unsigned int m) {
      _C[l][m] = val;
    }
  public:
    void setS (const double & val,
	       const unsigned int l, const unsigned int m) {
      _S[l][m] = val;
    }
  protected:
    std::vector<std::vector<double> > _C, _C_uncertainty;
    std::vector<std::vector<double> > _S, _S_uncertainty;  

  protected:
    orsa::Cache<unsigned int> _sample_points;
    
  public:
    void setCenterOfMass(const orsa::Vector & v) {
      _center_of_mass = v;
    }
  public:
    const orsa::Vector & getCenterOfMass() const {
      return _center_of_mass;
    }
  protected:
    orsa::Vector _center_of_mass, _center_of_mass_uncertainty;
    
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
    /* 
       const orsa::Cache<unsigned int> & order() const {
       return _order;
       } 
    */
    //
    const unsigned int order() const {
      return _order.getRef();
    }
  protected:
    orsa::Cache<unsigned int> _order;
  };
  
}; // namespace orsa

#endif // _ORSA_MULTIPOLE_
