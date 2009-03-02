#include "mainThread.h"

#include "multiminPhase.h"

#include <orsa/orbit.h>
#include <orsa/paulMoment.h>
#include <orsa/print.h>
#include <orsa/util.h>

#include <orsaSolarSystem/attitude.h>
#include <orsaSolarSystem/datetime.h>

#include <orsaSPICE/spice.h>
#include <orsaSPICE/spiceBodyPosVelCallback.h>
#include <orsaSPICE/spiceBodyInterpolatedPosVelCallback.h>

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
  longitude = orsa::Double(orsa::radToDeg() * (orsa::atan2(u.getY(),u.getX()))).get_d();  
  latitude  = orsa::Double(orsa::radToDeg() * (orsa::halfpi() - orsa::acos(u.getZ()))).get_d();
  
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
    sun->setMass(FromUnits(one(),Unit::MSUN));
    // sun->setMass(zero());
    //
    sun->isLightSource = true;
    //
    if (accurateSPICE) {
      SpiceBodyPosVelCallback * sbpvc = new SpiceBodyPosVelCallback(sun->getName());
      // sbpvc->setBodyName(sun->getName());
      orsa::IBPS ibps;
      ibps.translational = sbpvc;
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
    vesta->setMass(vestaMass.getRef());
    //
    if (accurateSPICE) {
      
      osg::ref_ptr<SpiceBodyPosVelCallback> sbpvc = new SpiceBodyPosVelCallback(vesta->getName());
      
      orsa::IBPS ibps;
      ibps.translational = sbpvc.get();
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
      
      // ORSA_DEBUG("vestaPeriod: %Fg",vestaPeriod.getRef().get_mpf_t());
      
      ibps.rotational = new orsaSolarSystem::ConstantZRotationEcliptic_RotationalBodyProperty(J2000(),
											      292.0*degToRad(),
											      twopi()/vestaPeriod.getRef(),
											      vestaPoleEclipticLongitude.getRef(),
											      vestaPoleEclipticLatitude.getRef());
      
      vesta->setInitialConditions(ibps);
    }
    
    {
      if (!vestaShape->read("vesta_thomas.dat")) {
	ORSA_ERROR("problems encountered while reading shape file...");
      }
      vesta->setShape(vestaShape.get());
    }
    
    bg->addBody(vesta.get());
  }
  
  // ORSA_DEBUG("before DAWN...");
  
  osg::ref_ptr<Body> dawn = new Body;
  orsa::Orbit dawnOrbit;
  {
    
    Double alpha = zero();
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
	
	osg::ref_ptr<orsa::Attitude> vesta_attitude = new BodyAttitude(vesta.get(),
								       bg.get());
	
	// const Matrix g2l = vesta->getAttitude()->globalToLocal(orbitEpoch.getRef());
	
	const Matrix g2l = vesta_attitude->globalToLocal(orbitEpoch.getRef());
	
	const Vector uVesta2Sun_local = (g2l*(rSun-rVesta).normalized()).normalized();
	
	const Double uSun_z = uVesta2Sun_local.getZ();
	
	{
	  // let's just call this in any case...

	  const Double i_z = cos(orbitInclination.getRef());
	  osg::ref_ptr<MultiminPhase> mmp = new MultiminPhase;
	  const Double mmpAlpha = mmp->getAlpha(fmod(fmod(orbitPhase.getRef(),twopi())+twopi(),twopi()),
						uVesta2Sun_local,
						orsa::Vector(zero(),
							     -sqrt(one()-i_z*i_z),
							     i_z));
	  alpha = mmpAlpha;
	}
	
      }
      
    }
    
    /* 
       ORSA_DEBUG("final alpha: %Ff",
       Double(radToDeg()*alpha).get_mpf_t());
    */
    
    dawn->setName("DAWN");
    dawn->setMass(zero());
    //
    // osg::ref_ptr<orsa::BodyInitialConditions> dawn_bic = new orsa::BodyInitialConditions;
    IBPS ibps;
    //
    orsa::Orbit orbit;
    //
    orbit.mu = vesta->getMu();
    orbit.a  = orbitRadius.getRef();
    orbit.e  = zero();
    orbit.i  = orbitInclination.getRef();
    orbit.omega_node       = alpha;
    orbit.omega_pericenter = zero();
    orbit.M                = zero();
    //
    orsa::Vector rOrbit, vOrbit;
    orbit.relativePosVel(rOrbit,vOrbit);
    
    // ORSA_DEBUG("rOrbit: %.20Fe",rOrbit.length().get_mpf_t());
    
    {
      osg::ref_ptr<orsa::Attitude> vesta_attitude = new BodyAttitude(vesta.get(),
								     bg.get());
      
      const Matrix l2g = vesta_attitude->localToGlobal(orbitEpoch.getRef());
      
      // print(l2g);
      
      // ORSA_DEBUG("l2g.getM11(): %Fe",l2g.getM11().get_mpf_t());
      
      rOrbit = l2g * rOrbit;
      vOrbit = l2g * vOrbit;
      
      // ORSA_DEBUG("rOrbit: %.20Fe",rOrbit.length().get_mpf_t());
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
    ibps.translational = new DynamicTranslationalBodyProperty;
    ibps.translational->setPosition(rOrbit);
    ibps.translational->setVelocity(vOrbit);
    
    // dawn->setInitialConditions(dawn_bic.get());
    dawn->setInitialConditions(ibps);
    
    {
      osg::ref_ptr<PaulMoment> dummy_pm = new PaulMoment(0);
      //
      dummy_pm->setM(one(),0,0,0);
      dummy_pm->setCenterOfMass(orsa::Vector(0,0,0));
      dummy_pm->setInertiaMoment(orsa::Matrix::identity());
      //
      dawn->setPaulMoment(dummy_pm.get());
    }
    
    // test
    // ORSA_DEBUG("========= DAWN time: %.6Ff",ibps.time.getRef().asDouble().get_mpf_t());
    
    bg->addBody(dawn.get());
    
    // test
    /* 
       ORSA_DEBUG("========= DAWN time: %.6Ff",
       dawn->getInitialConditions().time.getRef().asDouble().get_mpf_t());
    */
  
  }
  
  // test
  /* 
     ORSA_DEBUG("========= DAWN time: %.6Ff",
     dawn->getInitialConditions().time.getRef().asDouble().get_mpf_t());
  */
  

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
    
    osg::ref_ptr<orsa::Attitude> vesta_attitude = new BodyAttitude(vesta.get(),
								   bg.get());
    
    const orsa::Double halfParallelFOV   = cameraParallelFOV.getRef()   / 2;
    const orsa::Double halfOrthogonalFOV = cameraOrthogonalFOV.getRef() / 2;
    
    const orsa::Double deltaMax = 
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
    
    const orsa::Double min_i = cameraMinimumI.getRef();
    const orsa::Double max_i = cameraMaximumI.getRef();
    const orsa::Double max_e = cameraMaximumE.getRef();
    
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
	
	const orsa::Double orbitPeriod = localOrbit.period();
	
	const orsa::Double original_M  = localOrbit.M;
	
	localOrbit.M += orsa::twopi() * (t - orbitEpoch.getRef()).asDouble() / orbitPeriod;
	//
	orsa::Vector relPosDawn, relVelDawn; // relative to Vesta in terms of position and orientation
	//
	localOrbit.relativePosVel(relPosDawn,relVelDawn);
       	//
	localOrbit.M = original_M;
	
	const orsa::Vector rDawn = relPosDawn + rVesta;
	
	const orsa::Matrix g2l = vesta_attitude->globalToLocal(t);
	// print(g2l);
	
	// print(relPosDawn);
	relPosDawn = g2l * relPosDawn;
	relVelDawn = g2l * relVelDawn;
	// print(relPosDawn);
	
	const orsa::Vector relDawnToSun = g2l * (rSun-rDawn);
	
	// print(rSun - rVesta);
	const orsa::Vector relLightSource = g2l * (rSun - rVesta);
	// print(relLightSource);
	
	orsa::Double phase;
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
	    const orsa::Double angle = acos(u_vv*u_FOV_center);
	    const orsa::Vector u_projected = (u_vv - (u_vv*u_FOV_center)*u_FOV_center).normalized();
	    //
	    const orsa::Double angle_orthogonal = delta[k] * u_projected * u_FOV_orthogonal;
	    const orsa::Double angle_parallel   = delta[k] * u_projected * u_FOV_parallel;
	    
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
	       ORSA_DEBUG("delta: %Ff   deltaMax: %Ff   halfOrthogonalFOV: %Ff   halfParallelFOV: %Ff   parallel: %+Ff orthogonal: %+Ff   inFOV: %i",
	       delta[k].get_mpf_t(), 
	       deltaMax.get_mpf_t(),
	       halfOrthogonalFOV.get_mpf_t(),
	       halfParallelFOV.get_mpf_t(),
	       angle_parallel.get_mpf_t(),
	       angle_orthogonal.get_mpf_t(),
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
