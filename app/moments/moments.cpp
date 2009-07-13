// #include "vesta.h"
// #include "kleopatrashape.h"

//#include <orsa/multipole.h>
#include <orsa/unit.h>
#include <orsa/util.h>
#include <orsa/paulMoment.h>
#include <orsa/print.h>

using namespace std;
using namespace orsa;

class MassConcentration {
public:
  orsa::Vector position;
  double       radius;
  double       mass;
};

typedef std::list<MassConcentration> MassConcentrations;

void print(const MassConcentrations & mcs) {
  ORSA_DEBUG("MassConcentrations size: %i",mcs.size());
  
  orsa::Vector cm(0,0,0);
  double m=0;
  MassConcentrations::const_iterator it = mcs.begin();
  while (it != mcs.end()) {
    cm += (*it).mass*(*it).position;
    m  += (*it).mass;
    //
    /* ORSA_DEBUG("--MCS--");
       orsa::print((*it).position);
       orsa::print((*it).radius);
       orsa::print((*it).mass);
    */
    //
    ++it;
  }
  ORSA_DEBUG("cm: ...");
  orsa::print(cm/m);
}

class TestMassDistribution : public orsa::MassDistribution {
public:
  TestMassDistribution(const MassConcentrations & mcs_in) : MassDistribution(), mcs(mcs_in) { }
protected:
  ~TestMassDistribution() { }
public:
  double density(const Vector & v) const {
    double rho=0;
    MassConcentrations::const_iterator it = mcs.begin();
    while (it != mcs.end()) {
      if ((v-(*it).position).lengthSquared() < orsa::square((*it).radius)) {
	rho += (*it).mass / (4*orsa::pi()*orsa::cube((*it).radius)/3);
      }
      ++it;
    }
    if (rho < 0) {
      ORSA_ERROR("rho: %g",rho);
      orsa::print(v);
      exit(1);
    } 
    return rho;
  }
public:
  const MassConcentrations mcs;
};

void test_body(const double R,
	       const double totalMass,
	       const MassConcentrations & mcs) {
  
  print(mcs);
  
  // osg::ref_ptr<PaulMoment> pm = new PaulMoment(4);
  
  osg::ref_ptr<EllipsoidShape> vestaShapeEllipsoid = new EllipsoidShape(R,R,R);
  
  osg::ref_ptr<Shape> s = vestaShapeEllipsoid.get();
  
  const double   volume = 4*orsa::pi()*R*R*R/3; 
  
  const double meanDensity = totalMass/volume;
  //
  ORSA_DEBUG("volume: %g [km^3]",orsa::FromUnits(volume,orsa::Unit::KM,-3));
  //
  ORSA_DEBUG("total mass: %g [kg]",orsa::FromUnits(volume*meanDensity,orsa::Unit::KG,-1));
  //
  ORSA_DEBUG("mean density: %g g/cm^3",FromUnits(FromUnits(meanDensity,Unit::GRAM,-1),Unit::CM,3)); 
  
  // pm->setMassDistribution(new TestMassDistribution(mcs));
  osg::ref_ptr<orsa::MassDistribution> massDistribution = new TestMassDistribution(mcs);
  
  cerr << "shape bounding radius: " << FromUnits(s->boundingRadius(),Unit::KM,-1) << " KM" << endl;
  
  if (1) {
    const orsa::Box box = s->boundingBox();
    //
    cerr << "shape bounding box: [x] " << FromUnits(box.getXMax()-box.getXMin(),Unit::KM,-1) << " KM" << endl;
    cerr << "shape bounding box: [y] " << FromUnits(box.getYMax()-box.getYMin(),Unit::KM,-1) << " KM" << endl;
    cerr << "shape bounding box: [z] " << FromUnits(box.getZMax()-box.getZMin(),Unit::KM,-1) << " KM" << endl;
  }
  
  // pm->setShape(s.get());
  osg::ref_ptr<orsa::Shape> shape = s.get();
  
  const unsigned int _n_points = 100000000;
  //
  const unsigned int _random_seed = 85719;
  
  // pm->computeUsingShape(_n_points,_random_seed);
  // const orsa::Matrix I = pm->getInertiaMoment();
  
  const unsigned int order = 2;
  const unsigned int N = 10000;
  const int randomSeed = 95231;
  //
  double volumeAgain;
  orsa::Vector centerOfMass;
  orsa::Matrix shapeToLocal;
  orsa::Matrix localToShape;
  orsa::Matrix inertiaMatrix;
  orsa::PaulMoment * paulMoment;
  orsa::bodyInertialComputations(volumeAgain,
				 centerOfMass,
				 shapeToLocal,
				 localToShape,
				 inertiaMatrix,
				 &paulMoment,
				 order,
				 shape.get(),
				 massDistribution.get(),
				 N,
				 randomSeed);
  
  // print out...
  orsa::convert(paulMoment, R);
  
}

int main() {
  
  const double R = orsa::FromUnits(1000,orsa::Unit::KM);
  
  const double totalMass = orsa::FromUnits(1e22,orsa::Unit::KG);
  
  const unsigned int numMassConcentrations = 10;
  
  osg::ref_ptr<orsa::RNG> rng = new orsa::RNG(521031);
  
  MassConcentrations mcs;
  //
  for (unsigned int k=0; k<numMassConcentrations; ++k) {
    MassConcentration mc;
    //
    do {
      mc.position = 
	orsa::Vector(R*(2*rng->gsl_rng_uniform()-1),
		     R*(2*rng->gsl_rng_uniform()-1),
		     R*(2*rng->gsl_rng_uniform()-1));
    } while (mc.position.lengthSquared() > R*R);
    const double maxRadius = R-mc.position.length();
    mc.radius = maxRadius*rng->gsl_rng_uniform();
    mc.mass = totalMass/numMassConcentrations;
    //
    mcs.push_back(mc);
  }
  
  test_body(R,totalMass,mcs);
  
  return 0;
}
