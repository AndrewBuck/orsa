#include "NEAR.h"

#include <orsa/integrator_radau.h>
#include <orsa/orbit.h>
#include <orsa/paulMoment.h>
#include <orsa/print.h>
#include <orsa/util.h>

#include <orsaSPICE/spice.h>
#include <orsaSPICE/spiceBodyTranslationalCallback.h>
#include <orsaSPICE/spiceBodyRotationalCallback.h>

#include "shape.h"

using namespace orsa;
using namespace orsaSolarSystem;
using namespace orsaSPICE;

orsa::BodyGroup * run() {
  
  const orsa::Time t0 = gregorTime(2000,
				   7,
				   14,
				   12,
				   0,
				   0,
				   0);
  
  orsaSPICE::SPICE::instance()->loadKernel("de405.bsp");
  orsaSPICE::SPICE::instance()->loadKernel("ASTEPH.BSP");
  orsaSPICE::SPICE::instance()->loadKernel("ASTATT.BPC");
  orsaSPICE::SPICE::instance()->loadKernel("SCEPH.BSP");
  //
  orsaSPICE::SPICE::instance()->setDefaultObserver("SSB");
  
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
  
  osg::ref_ptr<Body> eros = new Body;
  {
    eros->setName("EROS");
    
    const double erosMass = orsa::FromUnits(6.69535497580e15,orsa::Unit::KG);
    
    orsa::IBPS ibps;
    
    osg::ref_ptr<ErosShape> erosShape = new ErosShape;
    if (!erosShape->read("eros001708.tab")) {
      ORSA_ERROR("problems encountered while reading shape file...");
    }
    osg::ref_ptr<orsa::Shape> shape = erosShape.get();
    
    if (0) {
      // compute coeffcients
      
      osg::ref_ptr<orsa::MassDistribution> massDistribution = new orsa::UniformMassDistribution;
      
      const unsigned int order = 16;
      const unsigned int N = 100000;
      const int randomSeed = 95231;
      //
      double volume;
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
      
      orsa::print(shapeToLocal);
      orsa::print(localToShape);
      
      // print out...
      std::vector< std::vector<double> > C, S, norm_C, norm_S;
      std::vector<double> J;
      const double R0 = FromUnits(16,orsa::Unit::KM);
      orsa::convert(C, S, norm_C, norm_S, J,
		    paulMoment, 
		    R0,
		    true);
      
      if (0) {
	// some tests...
	
	solve(paulMoment, 
	      norm_C, 
	      norm_S,
	      R0);
	
	// compare with original one...
	orsa::convert(C, S, norm_C, norm_S, J,
		      paulMoment, 
		      R0,
		      true);
	
      }
   
      ibps.inertial = new ConstantInertialBodyProperty(erosMass,
						       shape.get(),
						       centerOfMass,
						       shapeToLocal,
						       localToShape,
						       inertiaMatrix,
						       paulMoment);
      
    } else {
      
      // read from file
      
      unsigned int order;
      double eros_mu;
      double R0;
      std::vector< std::vector<double> > norm_C, norm_S;
      processGravityFile(order,
			 eros_mu,
			 R0,
			 norm_C, 
			 norm_S,
			 "GRE_004.dat");
      // print out
      ORSA_DEBUG("---- data from file ----");
      orsa::print(order);
      orsa::print(eros_mu);
      orsa::print(R0);
      for (unsigned int l=0; l<=order; ++l) {
	for (unsigned int m=0; m<=l; ++m) {
	  ORSA_DEBUG("norm_C[%02i][%02i] = %+e",l,m,norm_C[l][m]);
	  if (m!=0) {
	    ORSA_DEBUG("norm_S[%02i][%02i] = %+e",l,m,norm_S[l][m]);
	  }
	}
      }      
      ORSA_DEBUG("------------------------");
      
      // correct eros_mu, to get a better matching...
      // eros_mu = 446441.76205928; /* MKS */
      // eros_mu = 446450.77751572; /* MKS */
      
      // reduce order
      const unsigned int max_order=6;
      if (order > max_order) {
	order = max_order;
      }
      osg::ref_ptr<orsa::PaulMoment> paulMoment = 
	new orsa::PaulMoment(order);
      
      // find a paulmoment matching the spherical harmonics in the file
      solve(paulMoment.get(), 
	    norm_C, 
	    norm_S,
	    R0);
      
      // and check back if all matches
      {
	// new vars, overwriting outer ones...
	std::vector< std::vector<double> > C, S, norm_C, norm_S;
	std::vector<double> J; 
	orsa::convert(C, S, norm_C, norm_S, J,
		      paulMoment, 
		      R0,
		      true);
      }
      
      // placeholders
      //
      const orsa::Vector centerOfMass = orsa::Vector(0,0,0);
      //
      const orsa::Matrix shapeToLocal = orsa::Matrix::identity();
      const orsa::Matrix localToShape = orsa::Matrix::identity();
      //  
      const orsa::Matrix inertiaMatrix = orsa::Matrix::identity();
      
      ibps.inertial = new ConstantInertialBodyProperty(eros_mu/orsa::Unit::G(),
						       shape.get(),
						       centerOfMass,
						       shapeToLocal,
						       localToShape,
						       inertiaMatrix,
						       paulMoment.get());
    }
    
    ibps.translational = new SpiceBodyTranslationalCallback(eros->getName());
    
    ibps.rotational = new SpiceBodyRotationalCallback("IAU_EROS");
    
    eros->setInitialConditions(ibps);
    
    bg->addBody(eros.get());
  }
  
  osg::ref_ptr<Body> near = new Body;
  {  
    near->setName("NEAR");
    
    orsa::IBPS ibps;
    
    ibps.inertial = new PointLikeConstantInertialBodyProperty(0);
    
    ibps.translational = new SpiceBodyTranslationalCallback(near->getName());
    
    near->setInitialConditions(ibps);
    
    bg->addBody(near.get());
  }
  
  osg::ref_ptr<Body> clone = new Body;
  {
    clone->setName("CLONE");
    
    // dummy, positive mass, needed by SolarRadiationPressure (that is implemented as a propulstion...)
    const double CLONE_mass = orsa::FromUnits(1,orsa::Unit::KG);
    
    orsa::IBPS ibps;
    
    ibps.time = t0;
    
    ibps.inertial = new PointLikeConstantInertialBodyProperty(CLONE_mass);
    
    orsa::Vector r,v;
    bg->getInterpolatedPosVel(r,v,near.get(),t0);
    //
    ibps.translational = new orsa::DynamicTranslationalBodyProperty;
    // 
    ibps.translational->setPosition(r);
    ibps.translational->setVelocity(v);
    
    clone->setInitialConditions(ibps);
    
    // propulsion
    /* clone->propulsion = new SolarRadiationPressure(CLONE_mass,
       bg,
       sun.get(),
       clone.get());
    */
    
    bg->addBody(clone.get());
  }
  
  const double integrationAccuracy = 1.0e-6;
  
  if (0) {
    // multimin to find initial conditions for CLONE that better match NEAR in the long term
    const orsa::Time fitDuration = orsa::Time(0,0,10,0,0);
    orsa::Vector r,v;
    if (!betterInitialConditions(r,
				 v,
				 t0,
				 fitDuration,
				 integrationAccuracy,
				 near.get(),
				 clone.get(),
				 bg)) {
      ORSA_DEBUG("problems...");
      exit(0);
    }
    
    orsa::IBPS ibps = clone->getInitialConditions();
    ibps.translational->setPosition(r);
    ibps.translational->setVelocity(v);
    clone->setInitialConditions(ibps);
  }
  
  if (1) {
    // multimin to find an Eros mass that produces better match between NEAR and CLONE
    const orsa::Time fitDuration = orsa::Time(0,1,0,0,0);
    double erosMass;
    if (!betterErosMass(erosMass,
			t0,
			fitDuration,
			integrationAccuracy,
			eros.get(),
			near.get(),
			clone.get(),
			bg)) {
      ORSA_DEBUG("problems...");
      exit(0);
    }
    
    orsa::IBPS ibps = eros->getInitialConditions();
    osg::ref_ptr<orsa::ConstantInertialBodyProperty> inertial = new
      orsa::ConstantInertialBodyProperty(erosMass,
					 ibps.inertial->originalShape(),
					 ibps.inertial->centerOfMass(),
					 ibps.inertial->shapeToLocal(),
					 ibps.inertial->localToShape(),
					 ibps.inertial->inertiaMatrix(),
					 ibps.inertial->paulMoment());
    ibps.inertial = inertial.get();
    eros->setInitialConditions(ibps);
  }
  
  const orsa::Time duration = orsa::Time(2,0,0,0,0);
  const orsa::Time samplingPeriod = orsa::Time(0,0,0,10,0);  
  
  osg::ref_ptr<CustomIntegrator> radau = new CustomIntegrator(t0,true);
  radau->_accuracy = integrationAccuracy;
  radau->keepOnlyLastStep = true;
  // first call to output function
  // radau->singleStepDone(bg,t0,orsa::Time(0),orsa::Time(0));
  const bool goodIntegration = radau->integrate(bg,
						t0,
						t0+duration,
						samplingPeriod);
  
  if (0) {
    // test: estimate solar radiation pressure
    const orsa::Time dt = orsa::Time(0,0,0,1,0);
    
    // t0 already set
    const orsa::Time t1 = t0+dt;
    
    /* orsa::Vector sun_r1; 
       bg->getInterpolatedPosition(sun_r1,sun.get(),t1);
    */
    
    orsa::Vector NEAR_r1;
    bg->getInterpolatedPosition(NEAR_r1,near.get(),t1);
    
    orsa::Vector CLONE_r1;
    bg->getInterpolatedPosition(CLONE_r1,clone.get(),t1);
    
    const orsa::Vector dr = NEAR_r1 - CLONE_r1;
    
    // delta-acceleration
    /* const orsa::Vector da = ((2*dr)/orsa::square(dt.get_d()));
       ORSA_DEBUG("---acc---");
       orsa::print(da);
       orsa::print(da.normalized()*(NEAR_r1-sun_r1).normalized());
    */
  }
  
  if (0) {
    
    // this is done in the custom integrator now...
    
    // distance
    orsa::Time t=t0;
    while (t <= t0+duration) {
      
      orsa::Vector sun_r; 
      bg->getInterpolatedPosition(sun_r,sun.get(),t);
      
      orsa::Vector EROS_r;
      bg->getInterpolatedPosition(EROS_r,eros.get(),t);
      
      orsa::Vector NEAR_r;
      bg->getInterpolatedPosition(NEAR_r,near.get(),t);
      
      orsa::Vector CLONE_r;
      bg->getInterpolatedPosition(CLONE_r,clone.get(),t);
      
      ORSA_DEBUG("DD: %g %g %g %g %g",
		 (t-t0).get_d(),
		 (CLONE_r-NEAR_r).length(),
		 (EROS_r-NEAR_r).length(),
		 (CLONE_r-NEAR_r).normalized()*(EROS_r-NEAR_r).normalized(),
		 (sun_r-NEAR_r).length());
      
      orsa::Orbit orbit;
      //
      orbit.compute(near.get(),eros.get(),bg,t);
      ORSA_DEBUG("NEAR-ORBIT: %f %e %e %e",(t-t0).get_d(),orbit.a,orbit.e,orbit.i*orsa::radToDeg());
      //
      orbit.compute(clone.get(),eros.get(),bg,t);
      ORSA_DEBUG("CLONE-ORBIT: %f %e %e %e",(t-t0).get_d(),orbit.a,orbit.e,orbit.i*orsa::radToDeg());
      
      t += samplingPeriod;
    }
  }
  
  if (goodIntegration) {
    return bg;
  } else {
    return 0;
  }  
}


bool processGravityFile(unsigned int & order,
			double       & mu, /* G*m */
			double       & R0,
			std::vector< std::vector<double> > & norm_C, 
			std::vector< std::vector<double> > & norm_S,
			const std::string & fileName) {
  FILE * fp = fopen(fileName.c_str(),"r");
  if (fp==0) {
    ORSA_DEBUG("cannot open file [%s]",fileName.c_str());
    return false;
  }
  
  // units: km, kg, s
  
  char line[1024];
  unsigned int row=0;
  unsigned int l=0;
  unsigned int m=0;
  bool next_is_C=true;
  while (fgets(line,1024,fp) != 0) {
    if (row==0) {
      sscanf(line,"%d",&order);
      norm_C.resize(order+1);
      norm_S.resize(order+1);
      for (unsigned int l=0; l<=order; ++l) {
	norm_C[l].resize(l+1);
  	norm_S[l].resize(l+1);
      }
    } else if (row==1) {
      /* skip */
    } else if (row==2) {
      sscanf(line,"%lf %lf",&mu,&R0);
      mu = orsa::FromUnits(orsa::FromUnits(mu,orsa::Unit::KM,3),orsa::Unit::SECOND,-2);
      R0 = orsa::FromUnits(R0,orsa::Unit::KM);
    } else if (row==3) {
      // skip inertia matrix
    } else if (row==4) {
      // skip inertia matrix
    } else {
      double coeff[3];
      sscanf(line,"%lf %lf %lf",
	     &coeff[0],
	     &coeff[1],
	     &coeff[2]);
      for (unsigned int col=0; col<3; ++col) {
       	bool done=false;
	do {
	  if (next_is_C) {
	    norm_C[l][m]=coeff[col];
	    next_is_C=false;   
	    done=true;
	  } else {
	    if (m!=0) {
	      norm_S[l][m]=coeff[col];
	      next_is_C=true;
	      done=true;
	    } else {
	      next_is_C=true;
	    }
	  }	
	  if (next_is_C) {
	    if (m==l) {
	      ++l;
	      m=0;
	    } else {
	      ++m;
	    }
	  }
	  if (l>order) {
	    ORSA_DEBUG("problems...");
	    return false;
	  }
	} while (!done);
      }
    }
    ++row;
  }
  fclose(fp);
  return true;
}
