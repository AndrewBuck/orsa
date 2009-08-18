#include "DAWN.h"

#include "multiminPhase.h"

#include <orsa/integrator_radau.h>
#include <orsa/orbit.h>
#include <orsa/paulMoment.h>
#include <orsa/print.h>
#include <orsa/util.h>

#include <orsaSolarSystem/attitude.h>
#include <orsaSolarSystem/data.h>
#include <orsaSolarSystem/datetime.h>

#include <orsaSPICE/spice.h>
#include <orsaSPICE/spiceBodyTranslationalCallback.h>
#include <orsaSPICE/spiceBodyRotationalCallback.h>

#include "vesta.h"

using namespace orsa;
using namespace orsaSolarSystem;
using namespace orsaSPICE;

orsa::BodyGroup * run() {
  
  const orsa::Time t0 = gregorTime(2012,
				   1,
				   1,
				   0,
				   0,
				   0,
				   0);
  
  orsaSPICE::SPICE::instance()->loadKernel("de405.bsp");
  orsaSPICE::SPICE::instance()->loadKernel("vesta-2003-2013.bsp");
  //
  orsaSPICE::SPICE::instance()->setDefaultObserver("SSB");
  
  const double vestaMass = 1.35e-10*orsaSolarSystem::Data::MSun();
  ORSA_DEBUG("vestaMass = %g x MSun = %.12e kg",
	     vestaMass / orsaSolarSystem::Data::MSun(),
	     orsa::FromUnits(vestaMass,orsa::Unit::KG,-1));
  
  const double vestaPeriod = FromUnits(5.342128799,Unit::HOUR);
  
  const double vestaPoleEclipticLatitude  = degToRad()*59.2;
  
  const double vestaPoleEclipticLongitude = degToRad()*319.5;
  
  const orsa::Time orbitEpoch = orsaSolarSystem::gregorTime(2012,
							    1,
							    1,
							    0,
							    0,
							    0,
							    0);
  const double orbitRadius      = orsa::FromUnits(950.0,orsa::Unit::KM);
  const double orbitInclination = orsa::degToRad()*90.0;
  const double orbitPhase       = orsa::degToRad()*60.0;
  
  orsa::BodyGroup * bg = new orsa::BodyGroup;
  
  osg::ref_ptr<Body> sun = new Body;
  {
    sun->setName("SUN");
    sun->isLightSource = true;
    SpiceBodyTranslationalCallback * sbtc = new SpiceBodyTranslationalCallback(sun->getName());
    orsa::IBPS ibps;
    ibps.inertial = new PointLikeConstantInertialBodyProperty(orsaSolarSystem::Data::MSun());
    ibps.translational = sbtc;
    sun->setInitialConditions(ibps);
    bg->addBody(sun.get());
  }
  
  // add some planets
  if (1) {
    
    {
      Body * mercury = new Body;
      //
      mercury->setName("MERCURY BARYCENTER");
      // mercury->setMass(FromUnits(0.33022e24,Unit::KG));
      //
      SpiceBodyTranslationalCallback * sbtc = new SpiceBodyTranslationalCallback(mercury->getName());
      // sbtc->setBodyName(mercury->getName());
      orsa::IBPS ibps;
      ibps.inertial = new PointLikeConstantInertialBodyProperty(orsaSolarSystem::Data::MMercury());
      ibps.translational = sbtc;
      mercury->setInitialConditions(ibps);
      //
      bg->addBody(mercury);
    }

    {
      Body * venus = new Body;
      //
      venus->setName("VENUS BARYCENTER");
      // venus->setMass(FromUnits(4.8690e24,Unit::KG));
      //
      SpiceBodyTranslationalCallback * sbtc = new SpiceBodyTranslationalCallback(venus->getName());
      // sbtc->setBodyName(venus->getName());
      orsa::IBPS ibps;
      ibps.inertial = new PointLikeConstantInertialBodyProperty(orsaSolarSystem::Data::MVenus());
      ibps.translational = sbtc;
      venus->setInitialConditions(ibps);
      //
      bg->addBody(venus);
    }
    
    {
      Body * earth = new Body;
      //
      earth->setName("EARTH BARYCENTER");
      // earth->setMass(FromUnits(5.9742e24,Unit::KG));
      //
      SpiceBodyTranslationalCallback * sbtc = new SpiceBodyTranslationalCallback(earth->getName());
      // sbtc->setBodyName(earth->getName());
      orsa::IBPS ibps;
      ibps.inertial = new PointLikeConstantInertialBodyProperty(orsaSolarSystem::Data::MEarthMoon());
      ibps.translational = sbtc;
      earth->setInitialConditions(ibps);
      //
      bg->addBody(earth);
    }

    {
      Body * mars = new Body;
      //
      mars->setName("MARS BARYCENTER");
      // mars->setMass(FromUnits(0.64191e24,Unit::KG));
      //
      SpiceBodyTranslationalCallback * sbtc = new SpiceBodyTranslationalCallback(mars->getName());
      // sbtc->setBodyName(mars->getName());
      orsa::IBPS ibps;
      ibps.inertial = new PointLikeConstantInertialBodyProperty(orsaSolarSystem::Data::MMars());
      ibps.translational = sbtc;
      mars->setInitialConditions(ibps);
      //
      bg->addBody(mars);
    }
    
    {
      Body * jupiter = new Body;
      //
      jupiter->setName("JUPITER BARYCENTER");
      // jupiter->setMass(FromUnits(1898.8e24,Unit::KG));
      //
      SpiceBodyTranslationalCallback * sbtc = new SpiceBodyTranslationalCallback(jupiter->getName());
      // sbtc->setBodyName(jupiter->getName());
      orsa::IBPS ibps;
      ibps.translational = sbtc;
      ibps.inertial = new PointLikeConstantInertialBodyProperty(orsaSolarSystem::Data::MJupiter());
      jupiter->setInitialConditions(ibps);
      //
      bg->addBody(jupiter);
    }

    {
      Body * saturn = new Body;
      //
      saturn->setName("SATURN BARYCENTER");
      // saturn->setMass(FromUnits(568.50e24,Unit::KG));
      //
      SpiceBodyTranslationalCallback * sbtc = new SpiceBodyTranslationalCallback(saturn->getName());
      // sbtc->setBodyName(saturn->getName());
      orsa::IBPS ibps;
      ibps.translational = sbtc;
      ibps.inertial = new PointLikeConstantInertialBodyProperty(orsaSolarSystem::Data::MSaturn());
      saturn->setInitialConditions(ibps);
      //
      bg->addBody(saturn);
    }

    {
      Body * uranus = new Body;
      //
      uranus->setName("URANUS BARYCENTER");
      // uranus->setMass(FromUnits(86.625e24,Unit::KG));
      //
      SpiceBodyTranslationalCallback * sbtc = new SpiceBodyTranslationalCallback(uranus->getName());
      // sbtc->setBodyName(uranus->getName());
      orsa::IBPS ibps;
      ibps.translational = sbtc;
      ibps.inertial = new PointLikeConstantInertialBodyProperty(FromUnits(86.625e24,Unit::KG));
      uranus->setInitialConditions(ibps);
      //
      bg->addBody(uranus);
    }
    
    {
      Body * neptune = new Body;
      //
      neptune->setName("NEPTUNE BARYCENTER");
      // neptune->setMass(FromUnits(102.78e24,Unit::KG));
      //
      SpiceBodyTranslationalCallback * sbtc = new SpiceBodyTranslationalCallback(neptune->getName());
      // sbtc->setBodyName(neptune->getName());
      orsa::IBPS ibps;
      ibps.translational = sbtc;
      ibps.inertial = new PointLikeConstantInertialBodyProperty(orsaSolarSystem::Data::MNeptune());
      neptune->setInitialConditions(ibps);
      //
      bg->addBody(neptune);
    }
    
  }
  
  osg::ref_ptr<Body> vesta = new Body;
  {
    vesta->setName("VESTA");
    
    osg::ref_ptr<SpiceBodyTranslationalCallback> sbtc = new SpiceBodyTranslationalCallback(vesta->getName());
    //
    // sbtc->setBodyName(vesta->getName());
    //
    orsa::IBPS ibps;
    
    osg::ref_ptr<orsa::Shape> shape;
    osg::ref_ptr<VestaShape> vestaShapeThomas = new VestaShape;
    if (!vestaShapeThomas->read("vesta_thomas.dat")) {
      ORSA_ERROR("problems encountered while reading shape file...");
    }
    shape = vestaShapeThomas.get();
    
    // mass dist.
    osg::ref_ptr<orsa::MassDistribution> massDistribution;
    
    // UPDATE this when changing shape...
    double volume = FromUnits(7.875e7,Unit::KM,3); 
    
    double   coreDensity =
      FromUnits(FromUnits(7.9,Unit::GRAM),Unit::CM,-3);
    
    double mantleDensity =
      FromUnits(FromUnits(3.12,Unit::GRAM),Unit::CM,-3);
    
    orsa::Vector coreCenter = 
      orsa::Vector(orsa::FromUnits(-0.071 + 50.0,orsa::Unit::KM),
		   orsa::FromUnits(-1.067,orsa::Unit::KM),
		   orsa::FromUnits( 8.487,orsa::Unit::KM));
    
    const double totalMass = vestaMass;
    
    const double meanDensity = totalMass/volume;
    
    double    coreRadius =
      cbrt((3.0/(4.0*pi()))*volume*(meanDensity-mantleDensity)/(coreDensity-mantleDensity));
    
    if (1) {
      
      // CORE + MANTLE
      
      ORSA_DEBUG("core radius: %f km",FromUnits(coreRadius,Unit::KM,-1));
      
      // Choose one:
      massDistribution = new orsa::SphericalCorePlusMantleMassDistribution(coreCenter,
									   coreRadius,
									   coreDensity,
									   mantleDensity);
    } else if (0) {
      
      double coreVolume = 4*orsa::pi()/3.0*orsa::cube(coreRadius);
      
      const double fragmentRadius = orsa::FromUnits(10.0,orsa::Unit::KM);
      const double fragmentVolume = 4*orsa::pi()/3.0*orsa::cube(fragmentRadius);
      
      coreVolume -= 8*fragmentVolume;
      
      // updated core radius
      coreRadius =
	cbrt((3.0/(4.0*pi()))*coreVolume);
      
#warning COMPLETE HERE!
      // use MultipleSphericalFragmentsPlusMantleMassDistribution
      
    } else {
      
      // Uniform mass distribution
      
      coreDensity = mantleDensity = meanDensity;
      
      coreCenter = orsa::Vector(0,0,0);
      
      coreRadius = 0.0;
      
      massDistribution = new orsa::UniformMassDistribution;
      
    }
    
    const unsigned int order = 4;
    const unsigned int N = 10000;
    const int randomSeed = 95231;
    //
#warning fix the problem of computing volume after it is actually needed...
    // double volume;
    orsa::Vector centerOfMass;
    orsa::Matrix shapeToLocal;
    orsa::Matrix localToShape;
    orsa::Matrix inertiaMatrix;
    orsa::PaulMoment * paulMoment;
    orsa::bodyInertialComputations(volume,
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
    
    ORSA_DEBUG("volume: %.12e",volume);
    
    orsa::print(shapeToLocal);
    orsa::print(localToShape);
    
    if (1) {
      // print out...
      std::vector< std::vector<double> > C, S, norm_C, norm_S;
      std::vector<double> J;
      orsa::convert(C, S, norm_C, norm_S, J,
		    paulMoment, 
		    FromUnits(300,orsa::Unit::KM));
      
      ORSA_DEBUG("$\\rho_{m}$ & $%9.3f$ \\\\",orsa::FromUnits(orsa::FromUnits(mantleDensity,orsa::Unit::GRAM,-1),orsa::Unit::CM,3));
      ORSA_DEBUG("$\\rho_{c}$ & $%9.3f$ \\\\",orsa::FromUnits(orsa::FromUnits(  coreDensity,orsa::Unit::GRAM,-1),orsa::Unit::CM,3));
      ORSA_DEBUG("$R_{c}$    & $%9.3f$ \\\\",orsa::FromUnits(coreRadius,orsa::Unit::KM,-1));
      ORSA_DEBUG("\%\\hline");
      ORSA_DEBUG("$x_{c}$    & $%+9.3f$ \\\\",orsa::FromUnits(coreCenter.getX(),orsa::Unit::KM,-1));
      ORSA_DEBUG("$y_{c}$    & $%+9.3f$ \\\\",orsa::FromUnits(coreCenter.getY(),orsa::Unit::KM,-1));
      ORSA_DEBUG("$z_{c}$    & $%+9.3f$ \\\\",orsa::FromUnits(coreCenter.getZ(),orsa::Unit::KM,-1));
      ORSA_DEBUG("\\hline");
      ORSA_DEBUG("$x_{0}$    & $%+9.3f$ \\\\",orsa::FromUnits(centerOfMass.getX(),orsa::Unit::KM,-1));
      ORSA_DEBUG("$y_{0}$    & $%+9.3f$ \\\\",orsa::FromUnits(centerOfMass.getY(),orsa::Unit::KM,-1));
      ORSA_DEBUG("$z_{0}$    & $%+9.3f$ \\\\",orsa::FromUnits(centerOfMass.getZ(),orsa::Unit::KM,-1));
      ORSA_DEBUG("\%\\hline");
      for (unsigned int l=2; l<=order; ++l) {
	// J_l is minus C_l0, where C_l0 is not normalized
	ORSA_DEBUG("$J_{%i}$    & $%+9.6f$ \\\\",l,-C[l][0]);
      }
      ORSA_DEBUG("\%\\hline");
      for (unsigned int l=2; l<=order; ++l) {
	for (unsigned int m=0; m<=l; ++m) {
	  // LaTeX Tabular style
	  ORSA_DEBUG("$C_{%i%i}$   & $%+9.6f$ \\\\",l,m,norm_C[l][m]);
	  if (m!=0) {
	    ORSA_DEBUG("$S_{%i%i}$   & $%+9.6f$ \\\\",l,m,norm_S[l][m]);
	  }
	}
      }
    }
    
    ibps.inertial = new ConstantInertialBodyProperty(vestaMass,
						     shape.get(),
						     centerOfMass,
						     shapeToLocal,
						     localToShape,
						     inertiaMatrix,
						     paulMoment);
    
    ibps.translational = sbtc.get();
    
    ibps.rotational = new orsaSolarSystem::ConstantZRotationEcliptic_RotationalBodyProperty(orsaSolarSystem::J2000(),
											    292.0*degToRad(),
											    twopi()/vestaPeriod,
											    vestaPoleEclipticLongitude,
											    vestaPoleEclipticLatitude);
    vesta->setInitialConditions(ibps);
    
    bg->addBody(vesta.get());
  }
  
  // ORSA_DEBUG("before DAWN...");

  osg::ref_ptr<Body> dawn = new Body;
  {

    double alpha = 0;
    //
    // computations for the phase
    {
      orsa::Vector rVesta, vVesta;
      orsa::Vector rSun, vSun;
      if (bg->getInterpolatedPosVel(rVesta,
				    vVesta,
				    vesta.get(),
				    orbitEpoch) &&
	  bg->getInterpolatedPosVel(rSun,
				    vSun,
				    sun.get(),
				    orbitEpoch)) {
	
	const orsa::Matrix g2l = orsa::globalToLocal(vesta.get(),bg,orbitEpoch);
	
	const Vector uVesta2Sun_local = (g2l*(rSun-rVesta).normalized()).normalized();
	
	{
	  // let's just call this in any case...
	  
	  const double i_z = cos(orbitInclination);
	  osg::ref_ptr<MultiminPhase> mmp = new MultiminPhase;
	  const double mmpAlpha = mmp->getAlpha(fmod(fmod(orbitPhase,twopi())+twopi(),twopi()),
						uVesta2Sun_local,
						orsa::Vector(0,
							     -sqrt(1-i_z*i_z),
							     i_z));
	  alpha = mmpAlpha;
	}
	
      }
    }
    
    /*
       ORSA_DEBUG("final alpha: %f",
       double(radToDeg()*alpha));
    */

    dawn->setName("DAWN");
    // dawn->setMass(0);
    //
    // osg::ref_ptr<orsa::BodyInitialConditions> dawn_bic = new orsa::BodyInitialConditions;
    IBPS ibps;
    //
    orsa::Orbit orbit;
    //
    orbit.mu = orsa::Unit::G() * vestaMass;
    orbit.a  = orbitRadius;
    orbit.e  = 0;
    orbit.i  = orbitInclination;
    orbit.omega_node       = alpha;
    orbit.omega_pericenter = 0;
    orbit.M                = 0;
    //
    orsa::Vector rOrbit, vOrbit;
    orbit.relativePosVel(rOrbit,vOrbit);
    
    ORSA_DEBUG("rOrbit: %.20e",rOrbit.length());
    
    // ORSA_DEBUG("rotate for attitude! (and phase?)");
    {
      // osg::ref_ptr<orsa::Attitude> vesta_attitude = new BodyAttitude(vesta.get(),bg.get());
      
      // const Matrix l2g = vesta->getAttitude()->localToGlobal(orbitEpoch.getRef());
      //
      // const Matrix l2g = vesta_attitude->localToGlobal(orbitEpoch.getRef());
      
      const orsa::Matrix l2g = orsa::localToGlobal(vesta.get(),bg,orbitEpoch);
      
      // ORSA_DEBUG("l2g.getM11(): %e",l2g.getM11());
      
      rOrbit = l2g * rOrbit;
      vOrbit = l2g * vOrbit;

      ORSA_DEBUG("rOrbit: %.20e",rOrbit.length());
    }
    
    {
      orsa::Vector rVesta, vVesta;
      if (bg->getInterpolatedPosVel(rVesta,
				    vVesta,
				    vesta.get(),
				    orbitEpoch)) {
	rOrbit += rVesta;
	vOrbit += vVesta;
      } else {
	ORSA_DEBUG("problems!");
      }
    }
    
    /* 
       dawn_bic->setTime(orbitEpoch.getRef());
       dawn_bic->setPosition(rOrbit);
       dawn_bic->setVelocity(vOrbit);
    */
    //
    /* 
       dawn_bic->time     = orbitEpoch.getRef();
       dawn_bic->position = rOrbit;
       dawn_bic->velocity = vOrbit;
    */
    //
    // dummy, positive mass, needed by SolarRadiationPressure (that is implemented as a propulsion...)
    const double dawn_mass = orsa::FromUnits(1.0,orsa::Unit::KG);
    //
    ibps.time = orbitEpoch;
    //
    ibps.inertial = new PointLikeConstantInertialBodyProperty(dawn_mass);
    //
    ibps.translational = new DynamicTranslationalBodyProperty;
    ibps.translational->setPosition(rOrbit);
    ibps.translational->setVelocity(vOrbit);
    
    // dawn->setInitialConditions(dawn_bic.get());
    dawn->setInitialConditions(ibps);
    
    /* {
       osg::ref_ptr<PaulMoment> dummy_pm = new PaulMoment(0);
       //
       dummy_pm->setM(1,0,0,0);
       // dummy_pm->setCenterOfMass(orsa::Vector(0,0,0));
       // dummy_pm->setInertiaMoment(orsa::Matrix::identity());
       //
       dawn->setPaulMoment(dummy_pm.get());
       }
    */
    
    // SRP
    ORSA_DEBUG("**** SRP ****");
    const double solarRadiationPressure_B = 25.0; // MKS, kg/m^2
    dawn->propulsion = new SolarRadiationPressure(dawn_mass,
						  solarRadiationPressure_B,
						  bg,
						  sun.get(),
						  dawn.get());
    
    // test
    // ORSA_DEBUG("========= DAWN time: %.6f",ibps.time.getRef().get_d());
    
    bg->addBody(dawn.get());
    
    // test
    /* ORSA_DEBUG("========= DAWN time: %.6f",
       dawn->getInitialConditions().time.getRef().get_d());
    */
  }
  
  
  const double integrationAccuracy = 1.0e-6;
  
  const orsa::Time duration = orsa::Time(0,1,0,0,0);
  const orsa::Time samplingPeriod = orsa::Time(0,0,0,10,0);  
  
  orsa::Time dummyTime(0);
  osg::ref_ptr<CustomIntegrator> radau = new CustomIntegrator(t0,true);
  radau->_accuracy = integrationAccuracy;
  radau->keepOnlyLastStep = true;
  // first call to output function
  radau->singleStepDone(bg,t0,dummyTime,dummyTime);
  const bool goodIntegration = radau->integrate(bg,
						t0,
						t0+duration,
						samplingPeriod); 
  // last call to output function
  radau->singleStepDone(bg,t0+duration,dummyTime,dummyTime);
  
  if (goodIntegration) {
    return bg;
  } else {
    return 0;
  }  
}

