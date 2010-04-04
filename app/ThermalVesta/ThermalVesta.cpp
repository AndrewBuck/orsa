#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>

#include <orsa/orbit.h>
#include <orsa/util.h>
#include <orsaSolarSystem/attitude.h>
#include <orsaSolarSystem/data.h>
#include <orsaSolarSystem/print.h>

#include <orsaSPICE/spice.h>
#include <orsaSPICE/spiceBodyTranslationalCallback.h>
#include <orsaSPICE/spiceBodyRotationalCallback.h>

#include "vesta.h"

// All thermal equations and notation from:
// Spencer, Lebofsky, Sykes 1989.
// "Systematic Biases in Radiometric Diameter Determinations"
// Icarus 78, 337-254.

using namespace orsa;
using namespace orsaSolarSystem;
using namespace orsaSPICE;

const double bondAlbedo = 0.20;
const double emissivity = 0.90;

const double sigma = 5.67040e-8; // Stefan-Boltzmann
const double solar = 1370.0; // W m^-2 at 1 AU

const double rotationPeriod = 5.342128799*3600.0; // s

double omega() { return orsa::twopi()/rotationPeriod; } // s^-1

// k
const double thermalConductivity = 1.71e-4; // W m^-1 K^-1

// rho
const double density = 1750.0; // kg/m^3

// c
const double specificHeatCapacity = 750.0; // J kg^-1 K^-1

// Gamma
double thermalInertia() { return sqrt(thermalConductivity*density*specificHeatCapacity); } // Joule m^-2 s^-1/2 K^-1

// l_s
double skinDepth() { return sqrt(thermalConductivity/(density*specificHeatCapacity*omega())); } // m

// sub-solar Temperature, K
double Tss(const double hdist) {
  // return pow((1-bondAlbedo)*(solar/(hdist_AU*hdist_AU))/(emissivity*sigma),0.25);
  return pow((1-bondAlbedo)*(solar/pow(orsa::FromUnits(hdist,orsa::Unit::AU,-1),2))/(emissivity*sigma),0.25);
}

double theta(const double hdist) { return thermalInertia()*sqrt(omega())/(emissivity*sigma*pow(Tss(hdist),3)); }

class Slice {
public:
  double T;
};

typedef std::vector<Slice> ThermalData;
typedef std::vector<ThermalData>  History;

class Profile {
public:
  Profile(const int numSlices_in) :
    numSlices(numSlices_in),
    ls(skinDepth()) {
    data.resize(numSlices+1);
  }
public:
  void step(const double dx,
	    const double dt,
	    const double Fs) {
    
    // const double dX = dx/skinDepth();
    const double dX = dx/ls;
    
    old_data = data;
    
    // surface
    data[0].T =
      old_data[0].T
      + 2*dt*omega()/(dX*dX)*(old_data[1].T-old_data[0].T) 
      - 2*dt*sqrt(omega())/(thermalInertia()*dX)*(emissivity*sigma*pow(old_data[0].T,4)-(1-bondAlbedo)*Fs);
    
    // explicitly skip k==0 and k==numSlices
    for (int k=1; k<numSlices; ++k) {
      data[k].T = 
	old_data[k].T
	+ (dt*omega())/(dX*dX)*(old_data[k+1].T - 2*old_data[k].T + old_data[k-1].T);
    }
    
  }
public:
  // const double depth;
  const int numSlices;
  const double ls; // skinDepth
public:
  ThermalData data;
protected:
  ThermalData old_data;
};

/* class TimePeriod {
   public: 
   TimePeriod() { }
   TimePeriod(const orsa::Time & t1, const orsa::Time & t2) :
   tStart(std::min(t1,t2)),
   tStop(std::max(t1,t2)) { }
   public:
   orsa::Cache<orsa::Time> tStart, tStop;
   };
*/

int main(int argc,
	 char ** argv) {
  
  orsa::Debug::instance()->initTimer();
  
  if (argc != 3) {
    ORSA_DEBUG("Usage: %s <lat> <lon>",argv[0]);
    exit(0);
  }
  
  const double lat = atof(argv[1])*orsa::degToRad();
  const double lon = atof(argv[2])*orsa::degToRad();
  
  ORSA_DEBUG("lat: %+g [deg]   lon: %g [deg]",lat*orsa::radToDeg(),lon*orsa::radToDeg());
  
  ORSA_DEBUG("process ID: %i",getpid());
  
  orsaSPICE::SPICE::instance()->loadKernel("de405.bsp");
  orsaSPICE::SPICE::instance()->loadKernel("vesta-2003-2013.bsp");
  //
  orsaSPICE::SPICE::instance()->setDefaultObserver("SSB");
  
  const double vestaPeriod = FromUnits(rotationPeriod,Unit::SECOND);
  
  const double vestaPoleEclipticLatitude  = degToRad()*59.2;
  
  const double vestaPoleEclipticLongitude = degToRad()*319.5;
  
  osg::ref_ptr<orsa::BodyGroup> bg = new orsa::BodyGroup;
  
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
  
  osg::ref_ptr<Body> earth = new Body;
  {
    earth->setName("EARTH");
    SpiceBodyTranslationalCallback * sbtc = new SpiceBodyTranslationalCallback(earth->getName());
    orsa::IBPS ibps;
    ibps.inertial = new PointLikeConstantInertialBodyProperty(orsaSolarSystem::Data::MEarth());
    ibps.translational = sbtc;
    earth->setInitialConditions(ibps);
    bg->addBody(earth.get());
  }
  
  osg::ref_ptr<Body> vesta = new Body;
  {
    vesta->setName("VESTA");
    
    osg::ref_ptr<SpiceBodyTranslationalCallback> sbtc = new SpiceBodyTranslationalCallback(vesta->getName());
    
    orsa::IBPS ibps;
    
    osg::ref_ptr<orsa::Shape> shape;
    
    {
      osg::ref_ptr<VestaShape> vestaShapeThomas = new VestaShape;
      if (!vestaShapeThomas->read("vesta_thomas.dat")) {
	ORSA_ERROR("problems encountered while reading shape file...");
      }
      shape = vestaShapeThomas.get();
    }
    
    ibps.inertial = 
      new ConstantInertialBodyProperty(0.0,
				       shape.get(),
				       orsa::Vector(0,0,0),
				       orsa::Matrix::identity(),
				       orsa::Matrix::identity(),
				       orsa::Matrix::identity(),
				       0);
    
    ibps.translational = sbtc.get();
    
    ibps.rotational = new orsaSolarSystem::ConstantZRotationEcliptic_RotationalBodyProperty(orsaSolarSystem::J2000(),
											    292.0*degToRad(),
											    twopi()/vestaPeriod,
											    vestaPoleEclipticLongitude,
											    vestaPoleEclipticLatitude);
    vesta->setInitialConditions(ibps);
    
    bg->addBody(vesta.get());
  }
  
  if (0) {
    // Vesta's obliquity = angle between orbit pole and body pole
    const orsa::Time t = orsaSolarSystem::gregorTime(2010,1,1);
    osg::ref_ptr<orsa::Body> b = vesta.get();
    //
    orsa::Vector r,v;
    bg->getInterpolatedPosVel(r,v,b.get(),t);
    const orsa::Vector orbitPole = (orsa::externalProduct(r,v)).normalized();
    const orsa::Vector bodyPole = orsa::localToGlobal(b,bg,t)*orsa::Vector(0,0,1);
    const double obliquity = acos(orbitPole*bodyPole);
    ORSA_DEBUG("%s obliquity: %g [deg]",
	       b->getName().c_str(),
	       obliquity*orsa::radToDeg());
  }
  
  // normal element, using lat,lon
  orsa::Vector u_local_latlon = orsa::Vector(cos(lat)*cos(lon),
					     cos(lat)*sin(lon),
					     sin(lat));
  // orsa::print(u_local_latlon);
  
  orsa::Vector intersectionPoint;
  orsa::Vector normal_local;
  if (!vesta->getInitialConditions().inertial->originalShape()->rayIntersection(intersectionPoint,
										normal_local,
										u_local_latlon*orsa::FromUnits(1e4,orsa::Unit::KM), // apoint along the line, well outside the body
										-u_local_latlon,
										false)) {
    ORSA_DEBUG("problems...");
    exit(0);
  }	

  /* ORSA_DEBUG("intersectionPoint:");
     orsa::print(intersectionPoint);
     ORSA_DEBUG("normal_local:");
     orsa::print(normal_local);
  */
  //
  ORSA_DEBUG("u_local_latlon*normal_local: %f",u_local_latlon*normal_local);
  
  /* std::vector<TimePeriod> printPeriod;
     printPeriod.push_back(TimePeriod(orsaSolarSystem::gregorTime(2011,8,14.5),
     orsaSolarSystem::gregorTime(2011,8,16.5)));
  */
  
  // analysed period
  /* 
     const orsa::Time t0 = orsaSolarSystem::gregorTime(2011,8,15.0);
     // const orsa::Time tF = orsaSolarSystem::gregorTime(2011,8,16.0);
     const orsa::Time tF = t0+orsa::Time(orsa::FromUnits(vestaPeriod,orsa::Unit::MICROSECOND,-1));
  */
  // analysed period
  const orsa::Time t0 = orsaSolarSystem::gregorTime(2011,8,15.0);
  //
  orsa::Orbit orbit;
  orbit.compute(vesta.get(),sun.get(),bg.get(),t0);
  double orbitPeriod = orbit.period();
  ORSA_DEBUG("orbitPeriod: %g [year]",orsa::FromUnits(orbitPeriod,orsa::Unit::YEAR,-1));
  //
  const orsa::Time tF = t0+orsa::Time(orsa::FromUnits(orbitPeriod,orsa::Unit::MICROSECOND,-1));
  
  const double Q = orbit.a*(1.0+orbit.e);
  const double q = orbit.a*(1.0-orbit.e);
  
  const double ls = skinDepth();
  ORSA_DEBUG("skin depth: %g [m]",skinDepth());
  
  const int numSlicesPerSkinDepth = 3;
  const double dx = ls/numSlicesPerSkinDepth; // m
  
  ORSA_DEBUG("big-Theta at q: %g   at Q: %G",
	     theta(q),theta(Q));
  
  // steps in one thermal cycle = one rotation period
  /* const int numTimeStepsPerRotation = std::max(1000*numSlicesPerSkinDepth,
     (int)(10000*numSlicesPerSkinDepth/theta(q)));
     const double dt_s = rotationPeriod/numTimeStepsPerRotation; // s
  */
  //
  const int numTimeStepsPerRotation = std::max(10*numSlicesPerSkinDepth,
					       (int)(100*numSlicesPerSkinDepth/theta(q)));
  const double dt_s = rotationPeriod/numTimeStepsPerRotation; // s
  const orsa::Time dt = orsa::Time(dt_s*1000000);
  ORSA_DEBUG("dt: %f [s]",orsa::FromUnits(dt.get_d(),orsa::Unit::SECOND,-1));  
  
  // const int totalRotations = ceil(orbitPeriod/rotationPeriod); // sidereal rotations
  const int totalRotations = ceil((tF-t0).get_d()/rotationPeriod); // sidereal rotations
  ORSA_DEBUG("totalRotations: %i",totalRotations);
  
  const int numSkinDepths = 300;
  
  // enought slices to reach a depth with constant temperature
  Profile profile(numSkinDepths*numSlicesPerSkinDepth);  
  
  History history;
  // const unsigned int history_skip = numTimeStepsPerRotation/10000; // save memory, same results, affects print-out density of results too
  const unsigned int history_skip = 10; // save memory, similar results, affects print-out density of results too
  ORSA_DEBUG("history_skip: %i",history_skip);
  
  // initial temperature at max depth
  double deep_T = 0.0; // K
  for (unsigned int k=0; k<profile.data.size(); ++k) {
    profile.data[k].T = deep_T;
  }
  
  bool converged=false;
  while (!converged) {
    
    ORSA_DEBUG("deep_T: %g [K]",deep_T);
    
    ORSA_DEBUG("dt: %f [s]   dx: %g [m]",dt.get_d(),dx);
    
    bool stable=false;
    History old_history;
    unsigned int cycles=0;
    // const unsigned int stable_max_try=32;
    while (!stable) {
      old_history = history;
      history.clear();
      for (int p=0; p<=numTimeStepsPerRotation*totalRotations; ++p) {
	const orsa::Time t = t0+p*dt;
	// orsaSolarSystem::print(t);
	orsa::Vector rVesta;
	bg->getInterpolatedPosition(rVesta,vesta.get(),t);
	orsa::Vector rSun;
	bg->getInterpolatedPosition(rSun,sun.get(),t);
	const orsa::Vector dr = (rVesta-rSun);
	const orsa::Vector u_surface_to_sun = (-dr).normalized();
	const double hdist = dr.length();
	const double scaled_solar = solar/pow(orsa::FromUnits(hdist,orsa::Unit::AU,-1),2);
	const orsa::Matrix l2g = orsa::localToGlobal(vesta.get(),bg.get(),t);
	const orsa::Vector u_global_normal_surface_element = l2g*normal_local;
	const double Fs = scaled_solar*std::max(0.0,u_surface_to_sun*u_global_normal_surface_element);
	// ORSA_DEBUG("Fs: %g   scaled_solar: %g   p: %i   dt: %g",Fs,scaled_solar,p,dt.get_d());
	profile.step(dx,dt.get_d(),Fs);
	if (p%history_skip==0) history.push_back(profile.data);	
      }
      const double stable_eps = 1.0e-4;
      stable=true;
      if (old_history.size()==0) {
        stable=false;
      } else {
        for (unsigned int k=0; k<history[0].size(); ++k) {
          if (!stable) break;
          for (unsigned int j=0; j<history.size(); ++j) {
            if (!stable) break;
            if (fabs((history[j][k].T-old_history[j][k].T)/old_history[j][k].T) > stable_eps) {
	      ORSA_DEBUG("stable search did not converge yet, cycle # %03i   delta = %.9f > %g at bin: %03i/%03i time: %05i/%05i",
			 cycles,
			 fabs((history[j][k].T-old_history[j][k].T)/old_history[j][k].T),
			 stable_eps,
			 k,
			 history[0].size(),
			 j,
			 history.size());
              stable=false;
              break;
            }   
          }
        }
      }
      ++cycles;
      // if (cycles >= stable_max_try) break;
      
      // when debugging
      // #warning REMOVE!
      // if (cycles==10) exit(0);
    }
    //
    if (stable) ORSA_DEBUG("stable after %i cycles",cycles);
    
    // determine constant value of deep temperature
    std::vector<double> average;
    average.resize(history[0].size());
    for (unsigned int k=0; k<history[0].size(); ++k) {
      average[k]=0;
      for (unsigned int j=0; j<history.size(); ++j) {
	average[k] += history[j][k].T;
      }
      average[k] /= history.size();
      // ORSA_DEBUG("average[%03i] = %.6f",k,average[k]);
    }
    
    if (stable) {
      const double converged_eps = 1.0e-2;
      converged=true;
      for (unsigned int k=0; k<average.size(); ++k) {
	if (fabs((average[k]-average[0])/average[0]) > converged_eps) {
	  ORSA_DEBUG("not converged: %g > %g, k=%i",fabs((average[k]-average[0])/average[0]),converged_eps,k);
	  converged=false;
	  break;
	}
      }
    } else {
      converged=false;
    }
    
    // physical output
    if (1) {
      static unsigned int fileCounter=0;
      ++fileCounter;
      char filename[1024];
      sprintf(filename,"temperature.orbital.%+.2f.%.2f.v%03i.dat",lat*orsa::radToDeg(),lon*orsa::radToDeg(),fileCounter);
      ORSA_DEBUG("writing output file %s",filename);
      FILE * fp = fopen(filename,"w");
      for (unsigned int j=0; j<history.size(); ++j) {
       	const orsa::Time t = t0+j*history_skip*dt;
	orsa::Vector rVesta;
	bg->getInterpolatedPosition(rVesta,vesta.get(),t);
	orsa::Vector rSun;
	bg->getInterpolatedPosition(rSun,sun.get(),t);
	const orsa::Vector dr = (rVesta-rSun);
	const orsa::Vector u_surface_to_sun = (-dr).normalized();
	const double hdist = dr.length();
	const double scaled_solar = solar/pow(orsa::FromUnits(hdist,orsa::Unit::AU,-1),2);
	const orsa::Matrix l2g = orsa::localToGlobal(vesta.get(),bg.get(),t);
	const orsa::Matrix g2l = orsa::globalToLocal(vesta.get(),bg.get(),t);
	const orsa::Vector u_global_normal_surface_element = l2g*normal_local;
	// b2s vector to compute sub-solar longitude
	orsa::Vector b2s_local = g2l*(rSun-rVesta);
       	const double subSolarLongitude = atan2(b2s_local.getY(),b2s_local.getX());
	const double Fs = scaled_solar*std::max(0.0,u_surface_to_sun*u_global_normal_surface_element);
	fprintf(fp,"%10.5f %7.3f %10.6f %10.6f %10.6f %10.3f %.2f %.2f\n",
		orsaSolarSystem::timeToJulian(t), // j*history_skip*dt.get_d()/rotationPeriod,
		orsa::radToDeg()*fmod((lon+orsa::pi()-subSolarLongitude)+2*orsa::twopi(),orsa::twopi()), // centered at noon
		history[j][0].T,
		history[j][history[j].size()-2].T, // -2 because -1 is never changed by thermal algo
		orsa::FromUnits(hdist,orsa::Unit::AU,-1),
		Fs,
		lat*orsa::radToDeg(),
		lon*orsa::radToDeg());
      }
      fclose(fp);
      //
      char cmd[1024];
      sprintf(cmd,"cp %s temperature.orbital.%+.2f.%.2f.latest.dat",filename,lat*orsa::radToDeg(),lon*orsa::radToDeg());
      ORSA_DEBUG("executing: [%s]",cmd);
      int retval = system(cmd);
      if (retval != 0) ORSA_DEBUG("problems with the system call...");
    }
    
    // correction
    if (!converged) {
      const double new_deep_T = average[0];
      if (fabs((new_deep_T-deep_T)/deep_T) < 1e-9) {
	ORSA_DEBUG("temperature correction not converging, try to change integration steps, exiting");
	exit(0);
      }
      ORSA_DEBUG("new initial temperature: %.3f",new_deep_T);
      for (unsigned int l=0; l<profile.data.size(); ++l) {
	profile.data[l].T = new_deep_T;
      }
      deep_T = new_deep_T;
    } else {
      ORSA_DEBUG("converged");
    }
    
  }
  
  
  return 0;
}

