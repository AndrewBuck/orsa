#include "vesta.h"

#include <orsa/multipole.h>

int main() {
  
  osg::ref_ptr<Multipole> m = new Multipole(8);
  
  // we must assume a mean density for Vesta if we choose to use a Core/Mantle mass distribution
  const double meanDensity = FromUnits(FromUnits(3.8,Unit::GRAM),Unit::CM,-3);
  
  osg::ref_ptr<VestaShape> vestaShapeThomas = new VestaShape;
  if (!vestaShapeThomas->read("vesta_thomas.dat")) {
    ORSA_ERROR("problems encountered while reading shape file...");
    exit(0);
  }	
  
  /* values +/- 12 km OLD!
     osg::ref_ptr<EllipsoidShape> vestaShapeEllipsoid = new EllipsoidShape(FromUnits(280,Unit::KM),
     FromUnits(272,Unit::KM),
     FromUnits(227,Unit::KM));
  */
  //
  // values +/- 5 km
  osg::ref_ptr<EllipsoidShape> vestaShapeEllipsoid = new EllipsoidShape(FromUnits(289,Unit::KM),
									FromUnits(280,Unit::KM),
									FromUnits(229,Unit::KM));
  
  // CHOOSE (uncomment) shape here, and uncomment also the relative volume value
  //
  osg::ref_ptr<Shape> s = vestaShapeThomas.get();
  const double   volume = FromUnits(7.87e7,Unit::KM,3); 
  //
  /* 
     osg::ref_ptr<Shape> s = vestaShapeEllipsoid.get();
     const Double   volume = FromUnits(7.80e7,Unit::KM,3);
  */
  
  // CHOOSE (uncomment) mass distribution here (default = uniform);
  //
  /* 
     {
     const Double   coreDensity = FromUnits(FromUnits(5,Unit::GRAM),Unit::CM,-3);
     const Double mantleDensity = FromUnits(FromUnits(3,Unit::GRAM),Unit::CM,-3);
     //
     const Double    coreRadius = cbrt((3.0/(4.0*pi()))*volume*(meanDensity-mantleDensity)/(coreDensity-mantleDensity));
     //
     ORSA_DEBUG("core radius: %Ff km",FromUnits(coreRadius,Unit::KM,-1).get_mpf_t());
     //
     m->setMassDistribution(new SphericalCorePlusMantleMassDistribution(coreRadius,
     coreDensity,
     mantleDensity));
     }
  */
  
  cerr << "shape bounding radius: " << FromUnits(s->boundingRadius(),Unit::KM,-1) << " KM" << endl;
  
  if (1) {
    const orsa::Box box = s->boundingBox();
    //
    cerr << "shape bounding box: [x] " << FromUnits(box.getXMax()-box.getXMin(),Unit::KM,-1) << " KM" << endl;
    cerr << "shape bounding box: [y] " << FromUnits(box.getYMax()-box.getYMin(),Unit::KM,-1) << " KM" << endl;
    cerr << "shape bounding box: [z] " << FromUnits(box.getZMax()-box.getZMin(),Unit::KM,-1) << " KM" << endl;
  }
  
  // osg::ref_ptr<Multipole> m = new Multipole(4);
  //
  m->setShape(s.get());
  //
  // m->setMassDistribution(md.get());
  //
  {
    // check for saved multipole file
    
    // const unsigned int _n_points = 100000000;
    // const unsigned int _n_points = 20000000;
    // const unsigned int _n_points = 10000000;
    // const unsigned int _n_points = 1000000;
    const unsigned int _n_points = 100000000;
    // const unsigned int _n_points = 40000;
    // const unsigned int _n_points = 10000;
    //
    const unsigned int _random_seed = 85719;
    //
    char filename[1024];
    sprintf(filename,"mp.dat");
    //
    if (!(m->readFromFile(filename))) {
      m->computeUsingShape(_n_points,_random_seed,FromUnits(258.0,Unit::KM));
      m->writeToFile(filename);
    }
  }
  
  return 0;
}
