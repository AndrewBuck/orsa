#ifndef _ORSA_MASS_DISTRIBUTION_H_
#define _ORSA_MASS_DISTRIBUTION_H_

#include <osg/Referenced>
#include <osg/ref_ptr>

namespace orsa {
  
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
  
  // add parameter for center of mass position
  // all in "shape" coordinates
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
   
}; // namespace orsa

#endif // _ORSA_MASS_DISTRIBUTION_H_
