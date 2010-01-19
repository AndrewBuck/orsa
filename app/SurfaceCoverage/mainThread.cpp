#include "mainThread.h"

#include "multiminPhase.h"

#include <orsa/orbit.h>
#include <orsa/paulMoment.h>
#include <orsa/print.h>
#include <orsa/util.h>

#include <orsaSolarSystem/attitude.h>
#include <orsaSolarSystem/data.h>
#include <orsaSolarSystem/datetime.h>

#include <orsaSPICE/spice.h>
#include <orsaSPICE/spiceBodyTranslationalCallback.h>

#include <errno.h>

#include <SpiceUsr.h>

#include <gsl/gsl_rng.h>

#include "vesta.h"

using namespace orsa;
using namespace orsaSolarSystem;
using namespace orsaSPICE;

static void geo(double & latitude,
		double & longitude,
		const orsa::Vector & v) {
  const orsa::Vector u = v.normalized();
  longitude = orsa::radToDeg() * (atan2(u.getY(),u.getX()));  
  latitude  = orsa::radToDeg() * (orsa::halfpi() - acos(u.getZ()));
  
  if (longitude >  180) longitude -= 360;
  if (longitude < -180) longitude += 360;
}

static void geoSelect(double & tmpLatitude, 
		      double & tmpLongitude,
		      const unsigned int vertex,
		      const std::vector<orsa::Vector> & vertexVector,
		      const double faceLongitude) {
  geo(tmpLatitude, tmpLongitude, vertexVector[vertex]);
  if (fabs(tmpLongitude+360-faceLongitude) < fabs(tmpLongitude-faceLongitude)) { 
    tmpLongitude += 360;
  }	
  if (fabs(tmpLongitude-360-faceLongitude) < fabs(tmpLongitude-faceLongitude)) { 
    tmpLongitude -= 360;
  }
}

/*****/

MainThread::MainThread() :
  QThread() {
  // customIntegrator = new CustomIntegrator(this);
  
  SPICE::instance()->setDefaultObserver("SSB");
  //
  SPICE::instance()->loadKernel("de405.bsp");
  // SPICE::instance()->loadKernel("vesta_1900_2100.bsp");
  SPICE::instance()->loadKernel("vesta-2003-2013.bsp");
}

void MainThread::run() {
  
  aborted = false;
  
  // clear output data...
  /* faceVector_out.clear();
     vertexVector_out.clear();
     faceArea_out.clear();
  */
  //
  vertexCoverage_out.clear();
  allFaceAtOnce_out.clear();

  // const Time samplingPeriod(0,0,30,0,0);
  // const Time samplingPeriod(0,0,5,0,0);
  
  const bool accurateSPICE = true;
  
  // osg::ref_ptr<BodyGroup> bg = new BodyGroup;
  bg = new BodyGroup;

  osg::ref_ptr<Body> sun = new Body;
  {
    sun->setName("SUN");
    //
    // sun->setMass(FromUnits(1,Unit::MSUN));
    // sun->setMass(0);
    //
    sun->isLightSource = true;
    //
    if (accurateSPICE) {
      SpiceBodyTranslationalCallback * sbtc =
	new SpiceBodyTranslationalCallback(sun->getName());
      // sbpvc->setBodyName(sun->getName());
      orsa::IBPS ibps;
      ibps.inertial = new PointLikeConstantInertialBodyProperty(orsaSolarSystem::Data::MSun());
      ibps.translational = sbtc;
      sun->setInitialConditions(ibps);
    } else {
      /* 
	 SpiceBodyInterpolatedPosVelCallback * sbipvc =
	 new SpiceBodyInterpolatedPosVelCallback(sun->getName(),
	 orbitEpoch.getRef(),
	 orbitEpoch.getRef()+runDuration.getRef(),
	 samplingPeriod);
	 sun->setBodyPosVelCallback(sbipvc);
      */
    }
    //
    bg->addBody(sun.get());
  }
  
  osg::ref_ptr<Body> vesta = new Body;
  osg::ref_ptr<VestaShape> vestaShape = new VestaShape;
  {
    vesta->setName("VESTA");
    // vesta->setMass(vestaMass.getRef());
    //
    if (accurateSPICE) {
      
      osg::ref_ptr<SpiceBodyTranslationalCallback> sbtc = 
	new SpiceBodyTranslationalCallback(vesta->getName());
      
      orsa::IBPS ibps;
      ibps.inertial = new PointLikeConstantInertialBodyProperty(vestaMass.getRef());
      ibps.translational = sbtc.get();
      vesta->setInitialConditions(ibps);
    
    } else {
      /* 
	 SpiceBodyInterpolatedPosVelCallback * sbipvc =
	 new SpiceBodyInterpolatedPosVelCallback(vesta->getName(),
	 orbitEpoch.getRef(),
	 orbitEpoch.getRef()+runDuration.getRef(),
	 samplingPeriod);
	 vesta->setBodyPosVelCallback(sbipvc);
      */
    }
    
    {
      IBPS ibps = vesta->getInitialConditions();
      
      // ORSA_DEBUG("vestaPeriod: %g",vestaPeriod.getRef());
      
      {
	if (!vestaShape->read("vesta_thomas.dat")) {
	  ORSA_ERROR("problems encountered while reading shape file...");
	}
      }
      
      ibps.inertial = new ConstantInertialBodyProperty(vestaMass.getRef(),
						       vestaShape.get(),
						       orsa::Vector(0,0,0),
						       orsa::Matrix::identity(),
						       orsa::Matrix::identity(),
						       orsa::Matrix::identity(),
						       0);
      
      ibps.rotational = new orsaSolarSystem::ConstantZRotationEcliptic_RotationalBodyProperty(J2000(),
											      292.0*degToRad(),
											      twopi()/vestaPeriod.getRef(),
											      vestaPoleEclipticLongitude.getRef(),
											      vestaPoleEclipticLatitude.getRef());
      
      vesta->setInitialConditions(ibps);
    }
    
    bg->addBody(vesta.get());
  }
  
  // ORSA_DEBUG("before DAWN...");
  
  osg::ref_ptr<Body> dawn = new Body;
  orsa::Orbit dawnOrbit;
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
				    orbitEpoch.getRef()) &&
	  bg->getInterpolatedPosVel(rSun,
				    vSun,
				    sun.get(),
				    orbitEpoch.getRef())) {
	
	/* osg::ref_ptr<orsa::Attitude> vesta_attitude = new BodyAttitude(vesta.get(),
	   bg.get());
	*/
	
	// const Matrix g2l = vesta->getAttitude()->globalToLocal(orbitEpoch.getRef());
	
	// const Matrix g2l = vesta_attitude->globalToLocal(orbitEpoch.getRef());
	
	// const Vector uVesta2Sun_local = (g2l*(rSun-rVesta).normalized()).normalized();

	const orsa::Matrix g2l = orsa::globalToLocal(vesta.get(),bg.get(),orbitEpoch.getRef());
	
	const Vector uVesta2Sun_local = (g2l*(rSun-rVesta).normalized()).normalized();
	
	{
	  // let's just call this in any case...

	  const double i_z = cos(orbitInclination.getRef());
	  osg::ref_ptr<MultiminPhase> mmp = new MultiminPhase;
	  const double mmpAlpha = mmp->getAlpha(fmod(fmod(orbitPhase.getRef(),twopi())+twopi(),twopi()),
						uVesta2Sun_local,
						orsa::Vector(0,
							     -sqrt(1-i_z*i_z),
							     i_z));
	  alpha = mmpAlpha;
	}

	/* 
	   const double uSun_z = uVesta2Sun_local.getZ();
	   {
	   // let's just call this in any case...
	   const double i_z = cos(orbitInclination.getRef());
	   osg::ref_ptr<MultiminPhase> mmp = new MultiminPhase;
	   const double mmpAlpha = mmp->getAlpha(fmod(fmod(orbitPhase.getRef(),twopi())+twopi(),twopi()),
	   uVesta2Sun_local,
	   orsa::Vector(0,
	   -sqrt(1-i_z*i_z),
	   i_z));
	   alpha = mmpAlpha;
	   }
	*/
	
      }
      
    }
    
    dawn->setName("DAWN");
    // dawn->setMass(0);
    //
    // osg::ref_ptr<orsa::BodyInitialConditions> dawn_bic = new orsa::BodyInitialConditions;
    IBPS ibps;
    //
    orsa::Orbit orbit;
    //
    orbit.mu = orsa::Unit::G() * vestaMass.getRef();
    orbit.a  = orbitRadius.getRef();
    orbit.e  = 0;
    orbit.i  = orbitInclination.getRef();
    orbit.omega_node       = alpha;
    orbit.omega_pericenter = 0;
    orbit.M                = 0;
    //
    orsa::Vector rOrbit, vOrbit;
    orbit.relativePosVel(rOrbit,vOrbit);
    
    // ORSA_DEBUG("rOrbit: %.20e",rOrbit.length());
    
    {
      /* osg::ref_ptr<orsa::Attitude> vesta_attitude = new BodyAttitude(vesta.get(),
	 bg.get());
      */
      
      //  const Matrix l2g = vesta_attitude->localToGlobal(orbitEpoch.getRef());
      const orsa::Matrix l2g = orsa::localToGlobal(vesta.get(),bg.get(),orbitEpoch.getRef());
      
      // print(l2g);
      
      // ORSA_DEBUG("l2g.getM11(): %e",l2g.getM11());
      
      rOrbit = l2g * rOrbit;
      vOrbit = l2g * vOrbit;
      
      // ORSA_DEBUG("rOrbit: %.20e",rOrbit.length());
    }
    
    // global orbit wrt Vesta
    dawnOrbit.mu = orbit.mu;
    dawnOrbit.compute(rOrbit,vOrbit,orbit.mu);
    
    {
      orsa::Vector rVesta, vVesta;
      if (bg->getInterpolatedPosVel(rVesta,
				    vVesta,
				    vesta.get(),
				    orbitEpoch.getRef())) {
	rOrbit += rVesta;
	vOrbit += vVesta;
      } else {
	ORSA_DEBUG("problems!");
      }
    }
    
    ibps.time = orbitEpoch.getRef();
    //  
    ibps.inertial = new PointLikeConstantInertialBodyProperty(0);
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
       dummy_pm->setCenterOfMass(orsa::Vector(0,0,0));
       dummy_pm->setInertiaMoment(orsa::Matrix::identity());
       //
       dawn->setPaulMoment(dummy_pm.get());
       }
    */
    
    bg->addBody(dawn.get());
    
  }
  
  // OLD
  /* 
     emit progress(0);
     //
     const bool goodIntegration = customIntegrator->integrate(bg.get(),
     orbitEpoch.getRef(),
     orbitEpoch.getRef()+runDuration.getRef(),
     samplingPeriod);
     
     if (goodIntegration) {
     emit progress(100);
     }
  */
  
  {
    emit progress(0);
    
    /* osg::ref_ptr<orsa::Attitude> vesta_attitude = new BodyAttitude(vesta.get(),
       bg.get());
    */
    
    const double halfParallelFOV   = cameraParallelFOV.getRef()   / 2;
    const double halfOrthogonalFOV = cameraOrthogonalFOV.getRef() / 2;
    
    const double deltaMax = 
      sqrt(int_pow(halfParallelFOV,2) +
	   int_pow(halfOrthogonalFOV,2));
    
    const bool includeShadows = false;
    
    orsa::TriShape::AngleVector i, e, delta; 
    
    const orsa::TriShape::FaceVector faceVector = vestaShape->getFaceVector();
    
    // prepare output
    // one decisional trigger for all
    if (faceVector_out.size() != vestaShape->getFaceVector().size()) {
      faceVector_out   = vestaShape->getFaceVector();
      vertexVector_out = vestaShape->getVertexVector();
      faceArea_out.resize(vestaShape->getFaceVector().size());
      for (unsigned int f=0; f<faceVector.size(); ++f) {
	faceArea_out[f] = vestaShape->_getFaceArea(f);
      }
      
      // lat, lon
      latitude.resize(faceVector_out.size());
      longitude.resize(faceVector_out.size());
      for (unsigned int f=0; f<faceVector_out.size(); ++f) {
	latitude[f].resize(vertexVector_out.size());
	longitude[f].resize(vertexVector_out.size());
	
	const orsa::TriShape::TriIndex tri = faceVector_out[f];
	
	const orsa::Vector faceCenter =  
	  vertexVector_out[tri.i()] +
	  vertexVector_out[tri.j()] +
	  vertexVector_out[tri.k()];
	
	double faceLatitude, faceLongitude;
	geo(faceLatitude,
	    faceLongitude,
	    faceCenter);
	
	geoSelect(latitude[f][tri.i()], longitude[f][tri.i()], tri.i(), vertexVector_out, faceLongitude);
	geoSelect(latitude[f][tri.j()], longitude[f][tri.j()], tri.j(), vertexVector_out, faceLongitude);
	geoSelect(latitude[f][tri.k()], longitude[f][tri.k()], tri.k(), vertexVector_out, faceLongitude);
	
      }
    }
    vertexCoverage_out.resize(vestaShape->getVertexVector().size());
    allFaceAtOnce_out.resize(vestaShape->getFaceVector().size());
    
    const double min_i = cameraMinimumI.getRef();
    const double max_i = cameraMaximumI.getRef();
    const double max_e = cameraMaximumE.getRef();
    
    const orsa::Time tStop = orbitEpoch.getRef()+runDuration.getRef();
    const orsa::Time dt = cameraInterval.getRef();
    orsa::Time t = orbitEpoch.getRef();
    while (t <= tStop) {
      if (aborted) break;
      
      orsa::Vector rSun,   vSun;
      orsa::Vector rVesta, vVesta;
      if (bg->getInterpolatedPosVel(rSun,   vSun,   sun.get(), t) && 
	  bg->getInterpolatedPosVel(rVesta, vVesta, vesta.get(), t) ) {
	
        orsa::Orbit & localOrbit = dawnOrbit;
	
	const double orbitPeriod = localOrbit.period();
	
	const double original_M  = localOrbit.M;
	
	localOrbit.M += orsa::twopi() * (t - orbitEpoch.getRef()).get_d() / orbitPeriod;
	//
	orsa::Vector relPosDawn, relVelDawn; // relative to Vesta in terms of position and orientation
	//
	localOrbit.relativePosVel(relPosDawn,relVelDawn);
       	//
	localOrbit.M = original_M;
	
	const orsa::Vector rDawn = relPosDawn + rVesta;
	
	// const orsa::Matrix g2l = vesta_attitude->globalToLocal(t);
	const orsa::Matrix g2l = orsa::globalToLocal(vesta.get(),bg.get(),t);
	// print(g2l);
	
	// print(relPosDawn);
	relPosDawn = g2l * relPosDawn;
	relVelDawn = g2l * relVelDawn;
	// print(relPosDawn);
	
	const orsa::Vector relDawnToSun = g2l * (rSun-rDawn);
	
	// print(rSun - rVesta);
	const orsa::Vector relLightSource = g2l * (rSun - rVesta);
	// print(relLightSource);
	
	double phase;
	if (vestaShape->vertexIlluminationAngles(relLightSource,
						 relPosDawn,
						 phase,
						 i, 
						 e,
						 delta,
						 deltaMax,
						 includeShadows)) {
	  
	  // detailed check: within rectangular FOV?
	  const orsa::Vector u_FOV_center     = (-relPosDawn).normalized(); // Z component
	  const orsa::Vector u_tmp            = relDawnToSun.normalized();  // X component
	  const orsa::Vector u_FOV_orthogonal = externalProduct(u_FOV_center,u_tmp).normalized(); // Y direction
	  const orsa::Vector u_FOV_parallel   = externalProduct(u_FOV_center,u_FOV_orthogonal).normalized(); // redefine X direction
	  
	  std::vector<bool> localVertexCoverage;
	  localVertexCoverage.resize(vertexVector_out.size());
	  
	  for (unsigned int k=0; k<vertexVector_out.size(); ++k) {
	    
	    if (delta[k] > deltaMax) continue;
	    
	    const orsa::Vector & vv = vertexVector_out[k];
	    
	    // FOV components 
	    const orsa::Vector u_vv = (vv-relPosDawn).normalized();
	    //
	    const double angle = acos(u_vv*u_FOV_center);
	    const orsa::Vector u_projected = (u_vv - (u_vv*u_FOV_center)*u_FOV_center).normalized();
	    //
	    const double angle_orthogonal = delta[k] * u_projected * u_FOV_orthogonal;
	    const double angle_parallel   = delta[k] * u_projected * u_FOV_parallel;
	    
	    const bool inFOV = 
	      (fabs(angle_orthogonal) < halfOrthogonalFOV) && 
	      (fabs(angle_parallel)   < halfParallelFOV);
	    
	    if ( inFOV && 
		 (i[k] >= min_i) && 
		 (i[k] <= max_i) && 
		 (e[k] <= max_e) ) {
	      vertexCoverage_out[k].push_back(t);
	      localVertexCoverage[k] = true;
	    } else {
	      localVertexCoverage[k] = false;
	    }
	    
	    /* 
	       ORSA_DEBUG("delta: %f   deltaMax: %f   halfOrthogonalFOV: %f   halfParallelFOV: %f   parallel: %+f orthogonal: %+f   inFOV: %i",
	       delta[k], 
	       deltaMax,
	       halfOrthogonalFOV,
	       halfParallelFOV,
	       angle_parallel,
	       angle_orthogonal,
	       inFOV);
	    */
	    
	  }
	  
	  for (unsigned int f=0; f<faceVector.size(); ++f) {
	    const orsa::TriShape::TriIndex tri = faceVector[f];
	    if ( localVertexCoverage[tri.i()] && 
		 localVertexCoverage[tri.j()] && 
		 localVertexCoverage[tri.k()] ) {
	      allFaceAtOnce_out[f].push_back(t);
	    }
	  }
	  
	}
	
      } else {
	ORSA_DEBUG("problems!");
      }
      
      {
	if (runDuration.getRef().getMuSec() != 0) {
	  const int p = mpz_class((mpz_class("100")*(t-orbitEpoch.getRef()).getMuSec()) / runDuration.getRef().getMuSec()).get_ui();
	  emit progress(p);      
	}
      }
      
      t += dt;
    }
    
    if (!aborted) {
      emit progress(100);
    }
    
  }
  
}

void MainThread::abort() {
  // customIntegrator->abort();
  aborted = true;
}

void MainThread::reset() {

  vestaMass.reset();
  // vestaMassDistribution.reset();
  // vestaShapeModel.reset();
  vestaPeriod.reset();
  vestaPoleEclipticLatitude.reset();
  vestaPoleEclipticLongitude.reset();

  orbitEpoch.reset();
  orbitRadius.reset();
  orbitInclination.reset();
  orbitPhase.reset();
  
  cameraParallelFOV.reset();
  cameraOrthogonalFOV.reset();
  cameraInterval.reset();
  
  runDuration.reset();
  outputSPICEFile.reset();
  outputASCIIFile.reset();

}
