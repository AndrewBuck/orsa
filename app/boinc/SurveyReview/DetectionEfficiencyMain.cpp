#include "SurveyReview.h"
#include "skycoverage.h"
#include "eta.h"

#include <orsa/debug.h>
#include <orsa/multimin.h>
#include <orsa/print.h>
#include <orsa/statistic.h>
#include <orsa/util.h>

#include <orsaSolarSystem/datetime.h>
#include <orsaSolarSystem/observatory.h>
#include <orsaSolarSystem/orbit.h>

#include <orsaSPICE/spice.h>
#include <orsaSPICE/spiceBodyRotationalCallback.h>
#include <orsaSPICE/spiceBodyTranslationalCallback.h>

#include <orsaInputOutput/MPC_asteroid.h>
#include <orsaInputOutput/MPC_observations.h>

// CSPICE prototypes and definitions.      
#include <SpiceUsr.h>

class CustomMPCObservationsFile : public orsaInputOutput::MPCObservationsFile {
public:
  CustomMPCObservationsFile() :
    orsaInputOutput::MPCObservationsFile() {
    processedLines=0;
  }
public:
  bool processLine(const char * line) {
    
    const bool retVal = orsaInputOutput::MPCObservationsFile::processLine(line);
    /* if (retVal) {
       ORSA_DEBUG("ACCEPTED: [%s]",line);
       } else {
       ORSA_DEBUG("REJECTED: [%s]",line);
       }
    */
    ++processedLines;
    if ((processedLines>0) && (processedLines%100000==0)) {
      ORSA_DEBUG("lines processed: %i   selected: %i",processedLines,_data.size());
    }
    return retVal;
  }
protected:
  unsigned int processedLines;
};

class CustomMPCAsteroidFile : public orsaInputOutput::MPCAsteroidFile {
public:
  // absolute (ecliptic) positions at skyCoverage->epoch
  CustomMPCAsteroidFile(const orsa::Vector sunPosition_in,
			const orsa::Vector obsPosition_in) : 
    orsaInputOutput::MPCAsteroidFile(),
    sunPosition(sunPosition_in),
    obsPosition(obsPosition_in) {
    processedLines=0;
    observed=0;
  }
public:
  bool processLine(const char * line) {
    
    const bool retVal = orsaInputOutput::MPCAsteroidFile::processLine(line);
    
    if (!retVal) return retVal;
    
    // skip numbered orbits
    /* 
       if (_data[_data.size()-1].number.isSet()) {
       // remove from data
       _data.pop_back();
       return false;
       }
    */
    
    // a copy
    orsaSolarSystem::OrbitWithEpoch orbit = _data[_data.size()-1].orbit.getRef();
    
    const double orbitPeriod = orbit.period();
    
    const double original_M  = orbit.M;
    //
    orbit.M = original_M + fmod(orsa::twopi() * (skyCoverage->epoch.getRef()-orbit.epoch.getRef()).get_d() / orbitPeriod, orsa::twopi());
    orsa::Vector r;
    orbit.relativePosition(r);
    const orsa::Vector orbitPosition = r + sunPosition;
    // restore, important!
    orbit.M = original_M;
    
    // all in Ecliptic coords, as usual
    orsa::Vector dr = (orbitPosition - obsPosition).normalized();
    
    if (0) {
      // test:
      // if the object was observed,
      // check how far from the fields it is
      
      // does the orbit corresponds to an observed object?
      if (obsFile.get()) {
	
	bool present=false;
	orsaSolarSystem::OpticalObservation * obs;
	for (unsigned int k=0; k<obsFile->_data.size(); ++k) {
	  obs = dynamic_cast<orsaSolarSystem::OpticalObservation *> (obsFile->_data[k].get());
	  if (obs) {
	    if (!obs->mag.isSet()) {
	      // ORSA_DEBUG("mag not set, skipping");
	      continue;
	    }
	    if (obs->designation.isSet() && _data[_data.size()-1].designation.isSet()) {
	      if (obs->designation.getRef() == _data[_data.size()-1].designation.getRef()) {
		present=true;
		break;
	      }
	    }
	    if (obs->number.isSet() && _data[_data.size()-1].number.isSet()) {
	      if (obs->number.getRef() == _data[_data.size()-1].number.getRef()) {
		present=true;
		break;
	      }
	    }
	  }	  
	}
	if (present) {
	  const double minArc = skyCoverage->minDistance(dr.normalized());
	  
	  const orsa::Angle ra = fmod(atan2(dr.getY(),dr.getX())+orsa::twopi(),orsa::twopi());
	  const orsa::Angle dec = asin(dr.getZ()/dr.length());
	  
	  if (minArc > 0.0) {
	    if (obs->number.isSet()) {
	      ORSA_DEBUG("object (%i) present, min distance: %.2f [deg]",
			 obs->number.getRef(),
			 orsa::radToDeg()*minArc);
	    } else if (obs->designation.isSet()) {
	      ORSA_DEBUG("object [%s] present, min distance: %.2f [deg]",
			 obs->designation.getRef().c_str(),
			 orsa::radToDeg()*minArc);
	    } else {
	      ORSA_DEBUG("NONAME??");
	      exit(0);
	    }
	    /* orsa::print(orbit);
	       orsa::print(ra);
	       orsa::print(dec);
	    */
	  }
	}
      }
    }
    
    const orsa::Vector dr_nonint = dr.normalized();
    
    const double minArc = skyCoverage->minDistance(dr.normalized());
    if (0 && (minArc < 10.0*orsa::degToRad())) {
      // close enough, integrate
      
      osg::ref_ptr<orsa::BodyGroup> bg = new orsa::BodyGroup;
      
      osg::ref_ptr<orsa::Body> sun   = SPICEBody("SUN",orsaSolarSystem::Data::MSun());
      bg->addBody(sun.get());
      
      bg->addBody(SPICEBody("MERCURY BARYCENTER",orsaSolarSystem::Data::MMercury()));
      
      bg->addBody(SPICEBody("VENUS BARYCENTER",orsaSolarSystem::Data::MVenus()));
      
      bg->addBody(SPICEBody("EARTH BARYCENTER",orsaSolarSystem::Data::MEarthMoon()));
      
      bg->addBody(SPICEBody("MARS BARYCENTER",orsaSolarSystem::Data::MMars()));
      
      bg->addBody(SPICEBody("JUPITER BARYCENTER",orsaSolarSystem::Data::MJupiter()));
      
      bg->addBody(SPICEBody("SATURN BARYCENTER",orsaSolarSystem::Data::MSaturn()));
      
      bg->addBody(SPICEBody("URANUS BARYCENTER",orsaSolarSystem::Data::MUranus()));
      
      bg->addBody(SPICEBody("NEPTUNE BARYCENTER",orsaSolarSystem::Data::MNeptune()));
      
      orsa::Vector rSun, vSun;
      if (!bg->getInterpolatedPosVel(rSun,
				     vSun,
				     sun.get(),
				     orbit.epoch.getRef())) { 
	ORSA_DEBUG("problems");
      }
      
      // ORSA_DEBUG("orbit_a: %f [AU]",orsa::FromUnits(par->get("orbit_a"),orsa::Unit::AU,-1));
      
      orsa::Vector rOrbit, vOrbit;
      if (!orbit.relativePosVel(rOrbit,vOrbit)) {
	ORSA_DEBUG("problems");
      }
      rOrbit += rSun;
      vOrbit += vSun;
      
      osg::ref_ptr<orsa::Body> b = new orsa::Body;
      {
	b->setName("b");
	orsa::IBPS ibps;
	ibps.time = orbit.epoch.getRef();
	ibps.inertial = new orsa::PointLikeConstantInertialBodyProperty(0);
	ibps.translational = new orsa::DynamicTranslationalBodyProperty;
	ibps.translational->setPosition(rOrbit);
	ibps.translational->setVelocity(vOrbit);
	b->setInitialConditions(ibps);
	//
	bg->addBody(b.get());
      }
      
      osg::ref_ptr<orsa::IntegratorRadau> radau = new orsa::IntegratorRadau;
      radau->_accuracy = 1.0e-3;
      
      radau->integrate(bg.get(),
		       orbit.epoch.getRef(),
		       skyCoverage->epoch.getRef(),
		       orsa::Time(0,0,10,0,0));
      
      orsa::Vector bodyPosition, bodyVelocity;
      if (!bg->getInterpolatedPosVel(bodyPosition,
				     bodyVelocity,
				     b.get(),
				     skyCoverage->epoch.getRef())) { 
	ORSA_DEBUG("problems");
      }
      
      // replace dr
      dr = bodyPosition-obsPosition;
      // const double lightTimeDelay = dr.length()/orsa::Unit::c();
      // dr -= lightTimeDelay*(bodyVelocity);
      //
      // ORSA_DEBUG("dr: [AU]");
      // orsa::print(dr/orsa::FromUnits(1,orsa::Unit::AU));
      //
      // keep it in ecliptic coords, as usual
      // dr = orsaSolarSystem::eclipticToEquatorial()*dr;
      //
      // ORSA_DEBUG("dr.length(): %f [AU]",orsa::FromUnits(dr.length(),orsa::Unit::AU,-1)());
      //
      // const double ra_orbit  = fmod(atan2(dr.getY(),dr.getX())+orsa::twopi(),orsa::twopi());
      // const double dec_orbit = asin(dr.getZ()/dr.length());
      
      const orsa::Vector dr_integrated = dr.normalized();
      
      const double delta_arc = acos(dr_integrated.normalized()*dr_nonint.normalized());
      
      const double integrated_minArc = skyCoverage->minDistance(dr.normalized());
      
      ORSA_DEBUG("integrated minArc: %.2f [deg]   non-integrated: %.2f [deg]   distance: %.2f [deg] = %.2f [arcsec]",
		 orsa::radToDeg()*integrated_minArc,
		 orsa::radToDeg()*minArc,
		 orsa::radToDeg()*delta_arc,
		 orsa::radToArcsec()*delta_arc);
    }
    
    
    if (0) {
      const double minArc = skyCoverage->minDistance(dr.normalized());
      const bool inField  = skyCoverage->fastGet(dr.normalized());
      const orsaInputOutput::MPCAsteroidDataElement & orb = _data[_data.size()-1];
      if (orb.number.isSet()) {
	ORSA_DEBUG("object (%i) present, min distance: %.2f [deg]   in field: %i",
		   orb.number.getRef(),
		   orsa::radToDeg()*minArc,
		   inField);
      } else if (orb.designation.isSet()) {
	ORSA_DEBUG("object [%s] present, min distance: %.2f [deg]   in field: %i",
		   orb.designation.getRef().c_str(),
		   orsa::radToDeg()*minArc,
		   inField);
      }
      
    }    
    
    
    if (skyCoverage->fastGet(dr.normalized())) {
      
      // all good, keep the object
      
      if (0) {
	// plot
	
	// local dr, rotated, for plotting purposes only	
	
	// ORSA_DEBUG("---ROT---");
	// orsa::print(dr);
	const orsa::Vector dr_equatorial = orsaSolarSystem::eclipticToEquatorial()*dr;
	// orsa::print(dr);
	
	const orsa::Angle ra  = fmod(atan2(dr_equatorial.getY(),dr_equatorial.getX())+orsa::twopi(),orsa::twopi());
	const orsa::Angle dec = asin(dr_equatorial.getZ()/dr_equatorial.length());
	
	{
	  const orsaInputOutput::MPCAsteroidDataElement & orb = _data[_data.size()-1];
	  if (orb.number.isSet()) {
	    ORSA_DEBUG("[SKY-orb] (%i) %.6f %.6f",
		       orb.number.getRef(),ra.getRad(),dec.getRad());
	  } else if (orb.designation.isSet()) {
	    ORSA_DEBUG("[SKY-orb] [%s] %.6f %.6f",
		       orb.designation.getRef().c_str(),ra.getRad(),dec.getRad());
	  }
	}
      }
      
      // does the orbit corresponds to an observed object?
      if (obsFile.get()) {
	
	orsaSolarSystem::OpticalObservation * obs;
	for (unsigned int k=0; k<obsFile->_data.size(); ++k) {
	  obs = dynamic_cast<orsaSolarSystem::OpticalObservation *> (obsFile->_data[k].get());
	  if (obs) {
	    if (!obs->mag.isSet()) {
	      // ORSA_DEBUG("mag not set, skipping");
	      continue;
	    }
	    if (obs->designation.isSet() && _data[_data.size()-1].designation.isSet()) {
	      if (obs->designation.getRef() == _data[_data.size()-1].designation.getRef()) {
		++observed;
		break;
	      }
	    }
	    if (obs->number.isSet() && _data[_data.size()-1].number.isSet()) {
	      if (obs->number.getRef() == _data[_data.size()-1].number.getRef()) {
		++observed;
		break;
	      }
	    }
	  }	  
	}
      }
      
    } else {
      // object not in skyCoverage, remove from data
      _data.pop_back();
      return false;
    }
    
    ++processedLines;
    if ((processedLines>0) && (processedLines%100000==0)) {
      ORSA_DEBUG("lines processed: %i   selected: %i   observed: %i",processedLines,_data.size(),observed);
    }
    return retVal;
  }
public:
  osg::ref_ptr<SkyCoverage> skyCoverage;
  const orsa::Vector sunPosition;
  const orsa::Vector obsPosition;
public:
  osg::ref_ptr<CustomMPCObservationsFile> obsFile;
protected:
  unsigned int processedLines;
public:
  unsigned int observed;
};

int main(int argc, char ** argv) {
  
  orsa::Debug::instance()->initTimer();
  
  if (argc != 2) {
    printf("Usage: %s <sky_coverage_file>\n",argv[0]);
    exit(0);
  }
  
  ORSA_DEBUG("process ID: %i",getpid());
  
  const std::string basename = SkyCoverage::basename(argv[1]);
  
  orsaSPICE::SPICE::instance()->loadKernel("de405.bsp");
  
  const orsa::Time dt(0,0,10,0,0); // used to compute apparent velocity
  
  // extract observatory and date from input file name
  size_t found_underscore = std::string(argv[1]).find("_",0);
  size_t found_dot        = std::string(argv[1]).find(".",0);
  if ((found_underscore == std::string::npos) || (found_dot == std::string::npos)) {
    ORSA_DEBUG("not regular filename: %s",argv[1]);
    exit(0);
  }
  // ORSA_DEBUG("found: %i",found);
  std::string obsCode, compactDate;
  obsCode.assign(argv[1],0,found_underscore);
  compactDate.assign(argv[1],found_underscore+1,found_dot-found_underscore-1);
  //
  ORSA_DEBUG("    obsCode: [%s]",obsCode.c_str());
  ORSA_DEBUG("compactDate: [%s]",compactDate.c_str());
  
  // translate file obscode to MPC standard obscode
  obsCode = SkyCoverage::alias(obsCode);
  ORSA_DEBUG("MPC obsCode: [%s]",obsCode.c_str());
  
  osg::ref_ptr<orsaInputOutput::MPCObsCodeFile> obsCodeFile = new orsaInputOutput::MPCObsCodeFile;
  obsCodeFile->setFileName("obscode.dat");
  obsCodeFile->read();
  
  osg::ref_ptr<orsaSolarSystem::StandardObservatoryPositionCallback> obsPosCB =
    new orsaSolarSystem::StandardObservatoryPositionCallback(obsCodeFile.get());
  
  const orsaSolarSystem::Observatory & observatory = 
    obsCodeFile->_data.observatory[obsCode];
  
  // local midnight epoch
  orsa::Time epoch;
  int _y;
  if (strlen(compactDate.c_str())==7) {
    // seven character date format
    std::string s_year,s_dayOfYear;
    s_year.assign(compactDate,0,4);
    s_dayOfYear.assign(compactDate,4,3);
    _y = atoi(s_year.c_str());
    epoch = orsaSolarSystem::gregorTime(_y,
					1,
					atoi(s_dayOfYear.c_str())+1.0-observatory.lon.getRef()/orsa::twopi());
    orsa::print(epoch);
  } else {
    ORSA_DEBUG("problems...");
    exit(0);
  }
  const int year = _y;
  
  osg::ref_ptr<SkyCoverage> skyCoverage = new SkyCoverage;
  //
  skyCoverage->obscode = obsCode;
  skyCoverage->epoch   = epoch;
  //
  {
    
    FILE * fp;
    char line[1024];
    
    {
      fp = fopen(argv[1],"r");
      if (fp == 0) {
	ORSA_DEBUG("cannot open field file [%s]",argv[1]);
	exit(0);
      }
      double x1,x2,x3,x4; // ra
      double y1,y2,y3,y4; // dec
      double V;           // limiting mag
      //
      while (fgets(line,1024,fp)) {
	if (9 == gmp_sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
			    &x1,&y1,
			    &x2,&y2,
			    &x3,&y3,
			    &x4,&y4,
			    &V)) {
	  SkyCoverage::normalize(x1,y1);
	  SkyCoverage::normalize(x2,y2);
	  SkyCoverage::normalize(x3,y3);
	  SkyCoverage::normalize(x4,y4);
	  //
	  skyCoverage->setField(x1,y1,x2,y2,x3,y3,x4,y4,V);
	} 
      }
      fclose(fp);
    }
    
  }
  
  if (0) {
    // test skycoverage
    const int randomSeed=35092;
    osg::ref_ptr<orsa::RNG> rnd = new orsa::RNG(randomSeed); 
    for (unsigned int j=0; j<10000000; ++j) {
      const orsa::Angle ra  = rnd->gsl_rng_uniform()*orsa::twopi();
      const orsa::Angle dec = (2*rnd->gsl_rng_uniform()-1)*orsa::halfpi();
      if (skyCoverage->fastGet(ra,dec)) {
	ORSA_DEBUG("[SKY-rnd] (RND) %.6f %.6f",
		   ra.getRad(),dec.getRad());
      }
    }
  }
  

  osg::ref_ptr<orsa::BodyGroup> bg = new orsa::BodyGroup;
  
  // SUN
  osg::ref_ptr<orsa::Body> sun   = SPICEBody("SUN",orsaSolarSystem::Data::MSun());
  bg->addBody(sun.get());
  
  // EARTH
  /* osg::ref_ptr<orsa::Body> earth = SPICEBody("EARTH",orsaSolarSystem::Data::MEarth());
     bg->addBody(earth.get());
  */
  
  // MOON
  /* osg::ref_ptr<orsa::Body> moon  = SPICEBody("MOON",orsaSolarSystem::Data::MMoon());
     bg->addBody(moon.get());
  */
  
  orsa::Vector r;
  bg->getInterpolatedPosition(r,sun.get(),epoch);
  const orsa::Vector sunPosition = r;
  bg->getInterpolatedPosition(r,sun.get(),epoch+dt);
  const orsa::Vector sunPosition_dt = r;
  obsPosCB->getPosition(r,obsCode,epoch);
  const orsa::Vector obsPosition = r;
  obsPosCB->getPosition(r,obsCode,epoch+dt);
  const orsa::Vector obsPosition_dt = r;
  
  osg::ref_ptr<CustomMPCObservationsFile> obsFile =
    new CustomMPCObservationsFile;
  {
    obsFile->select_startEpoch = epoch - orsa::Time(0,12,0,0,0);
    obsFile->select_stopEpoch  = epoch + orsa::Time(0,12,0,0,0);
    ORSA_DEBUG("select start/stop:");
    orsa::print(obsFile->select_startEpoch.getRef());
    orsa::print(obsFile->select_stopEpoch.getRef());
    obsFile->select_obsCode = obsCode;
    /* obsFile->setFileName("mpn.arc.gz");
       obsFile->read();
       obsFile->setFileName("mpu.arc.gz");
       obsFile->read();
    */
    char filename[1024];
    sprintf(filename,"%i.obs.gz",year);
    obsFile->setFileName(filename);
    obsFile->read();
    //
    ORSA_DEBUG("total selected observations: %i",obsFile->_data.size());
  }
  
  {
    unsigned int inField=0,inFieldCandidates=0;
    orsaSolarSystem::OpticalObservation * obs;
    for (unsigned int k=0; k<obsFile->_data.size(); ++k) {
      obs = dynamic_cast<orsaSolarSystem::OpticalObservation *> (obsFile->_data[k].get());
      if (obs) {
	if (!obs->mag.isSet()) {
	  // ORSA_DEBUG("mag not set, skipping");
	  continue;
	}
	++inFieldCandidates;
	if (skyCoverage->fastGet(obs->ra.getRef(),
				 obs->dec.getRef())) {	
	  ++inField;
	  /* if (obs->number.isSet()) {
	     ORSA_DEBUG("V field: %.1f   V obs: %.1f   obj: (%i)",V,obs->mag.getRef(),obs->number.getRef());
	     } else if (obs->designation.isSet()) {
	     ORSA_DEBUG("V field: %.1f   V obs: %.1f   obj: [%s]",V,obs->mag.getRef(),obs->designation.getRef().c_str());
	     } else {
	     ORSA_DEBUG("NONAME??");
	     exit(0);
	     }
	  */
	} else {
	  // ORSA_DEBUG("this one out:");
	  // orsa::print(obs->epoch.getRef());
	  // orsa::print(obs->ra.getRef());
	  // orsa::print(obs->dec.getRef());
	}
      }
    }
    // debug
    ORSA_DEBUG("inField success: %i/%i",inField,inFieldCandidates);
    // this ratio should be smaller than 1.0, because the field declared
    // in the sky coverage files is smaller than the real data CCD field.
  }
  
  osg::ref_ptr<CustomMPCAsteroidFile> orbitFile = 
    new CustomMPCAsteroidFile(sunPosition,obsPosition);
  orbitFile->skyCoverage = skyCoverage.get();
  orbitFile->obsFile = obsFile.get();
  orbitFile->setFileName("MPCORB.DAT");
  orbitFile->read();
  orbitFile->setFileName("NEA.DAT");
  orbitFile->read();
  ORSA_DEBUG("selected orbits: %i   observed: %i",
	     orbitFile->_data.size(),
	     orbitFile->observed);
  
  {
    // dump lists
    if (0) {
      ORSA_DEBUG("--DUMP-OBS---");
      orsaSolarSystem::OpticalObservation * obs;
      for (unsigned int k=0; k<obsFile->_data.size(); ++k) {
	obs = dynamic_cast<orsaSolarSystem::OpticalObservation *> (obsFile->_data[k].get());
	if (obs) {
	  if (obs->number.isSet()) {
	    ORSA_DEBUG("[SKY-obs] (%i) %.6f %.6f",obs->number.getRef(),obs->ra.getRef().getRad(),obs->dec.getRef().getRad());
	  } else if (obs->designation.isSet()) {
	    ORSA_DEBUG("[SKY-obs] [%s] %.6f %.6f",obs->designation.getRef().c_str(),obs->ra.getRef().getRad(),obs->dec.getRef().getRad());
	  }
	}
      }
    }
    if (0) { 
      ORSA_DEBUG("--DUMP-ORB---");
      for (unsigned int k=0; k<orbitFile->_data.size(); ++k) {
	const orsaInputOutput::MPCAsteroidDataElement & orb = orbitFile->_data[k];
	if (orb.number.isSet()) {
	  ORSA_DEBUG("(%i)",orb.number.getRef());
	} else if (orb.designation.isSet()) {
	  ORSA_DEBUG("[%s]",orb.designation.getRef().c_str());
	}
      }
    }  
  }
  
  // start to work at efficiency
  std::vector<EfficiencyData> etaData;
  for (unsigned int korb=0; korb<orbitFile->_data.size(); ++korb) {
    orsaSolarSystem::OrbitWithEpoch orbit = orbitFile->_data[korb].orbit.getRef();
    const double orbitPeriod = orbit.period();
    const double original_M  = orbit.M;
    //
    orbit.M = original_M + fmod(orsa::twopi() * (skyCoverage->epoch.getRef()-orbit.epoch.getRef()).get_d() / orbitPeriod, orsa::twopi());
    orsa::Vector r;
    orbit.relativePosition(r);
    const orsa::Vector orbitPosition = r + sunPosition;
    // not at t+dt
    orbit.M = original_M + fmod(orsa::twopi() * (skyCoverage->epoch.getRef()+dt-orbit.epoch.getRef()).get_d() / orbitPeriod, orsa::twopi());
    orsa::Vector r_dt;
    orbit.relativePosition(r_dt);
    const orsa::Vector orbitPosition_dt = r_dt + sunPosition_dt;
    // restore, important!
    orbit.M = original_M;
    // all in Ecliptic coords, as usual
    const orsa::Vector orb2obs = obsPosition - orbitPosition;
    const orsa::Vector orb2sun = sunPosition - orbitPosition;
    const orsa::Vector orb2obs_dt = obsPosition_dt - orbitPosition_dt;
    const orsa::Vector orb2sun_dt = sunPosition_dt - orbitPosition_dt;
    const double phaseAngle = acos(orb2sun.normalized()*orb2obs.normalized());
    //
    bool observed=false;
    orsaSolarSystem::OpticalObservation * obs;
    for (unsigned int kobs=0; kobs<obsFile->_data.size(); ++kobs) {
      obs = dynamic_cast<orsaSolarSystem::OpticalObservation *> (obsFile->_data[kobs].get());
      if (obs) {
	if (!obs->mag.isSet()) {
	  // ORSA_DEBUG("mag not set, skipping");
	  continue;
	}
	if (obs->designation.isSet() && orbitFile->_data[korb].designation.isSet()) {
	  if (obs->designation.getRef() == orbitFile->_data[korb].designation.getRef()) {
	    observed=true;
	    break;
	  }
	}
	if (obs->number.isSet() && orbitFile->_data[korb].number.isSet()) {
	  if (obs->number.getRef() == orbitFile->_data[korb].number.getRef()) {
	    observed=true;
	    break;
	  }
	}
      }	  
    }
    EfficiencyData ed;
    ed.H = orbitFile->_data[korb].H.getRef();
    if (orbitFile->_data[korb].number.isSet()) {
      ed.number = orbitFile->_data[korb].number.getRef();
    }
    if (orbitFile->_data[korb].designation.isSet()) {
      ed.designation = orbitFile->_data[korb].designation.getRef();
    }
    ed.V = apparentMagnitude(orbitFile->_data[korb].H.getRef(),
			     phaseAngle,
			     orb2obs.length(),
			     orb2sun.length()); 
    ed.apparentVelocity = acos(orb2obs_dt.normalized()*orb2obs.normalized())/dt.get_d();
    ed.observed = observed;
    etaData.push_back(ed);
  }
  //
  {
    char filename[1024];
    sprintf(filename,"%s.allEta.dat",basename.c_str());
    FILE * fp_allEta = fopen(filename,"w");
    for (unsigned int k=0; k<etaData.size(); ++k) {
      const EfficiencyData & ed = etaData[k];
      char id[1024];
      if (ed.number.isSet()) {
	sprintf(id,"%i",ed.number.getRef());
      } else if (ed.designation.isSet()) {
	sprintf(id,"%s",ed.designation.getRef().c_str());
      }
      fprintf(fp_allEta,
	      "%5.2f %7s %5.2f %10.6f %1i\n",
	      ed.H.getRef(),
	      id,
	      ed.V.getRef(),
	      orsa::FromUnits(ed.apparentVelocity.getRef()*orsa::radToArcsec(),orsa::Unit::HOUR), // arcsec/hour
	      ed.observed.getRef());
    }
    fclose(fp_allEta);
  }
  
  return 0;
}
