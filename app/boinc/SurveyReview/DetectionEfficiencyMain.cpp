#include "SurveyReview.h"
#include "grain.h"
#include "skycoverage.h"

#include <orsa/debug.h>
#include <orsa/print.h>
#include <orsa/statistic.h>

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
    if (_data[_data.size()-1].number.isSet()) {
      // remove from data
      _data.pop_back();
      return false;
    }
    
    // a copy
    orsaSolarSystem::OrbitWithEpoch orbit = _data[_data.size()-1].orbit.getRef();
    
    const double orbitPeriod = orbit.period();
    
    const double original_M  = orbit.M;
    //
    orbit.M = original_M + fmod(orsa::twopi() * (skyCoverage->epoch.getRef()-orbit.epoch.getRef()).get_d() / orbitPeriod, orsa::twopi());
    orsa::Vector r;
    orbit.relativePosition(r);
    const orsa::Vector orbitPosition = r + sunPosition;
    //
    orbit.M = original_M;
    
    const orsa::Vector dr = orsaSolarSystem::eclipticToEquatorial()*(orbitPosition - obsPosition).normalized();
    
#warning lightTimeDelay... dr -= lightTimeDelay*(refBodyVelocity+relativeVelocity);
    
    // orsa::print(dr_nightStart);
    // orsa::print(dr_nightStop);
    
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
    

    {
      // debug
      if (_data[_data.size()-1].designation.isSet()) {
	if (_data[_data.size()-1].designation.getRef() == "K06VA1J") {
	  ORSA_DEBUG("---SPECIAL-CASE---");
	  const orsa::Angle ra = fmod(atan2(dr.getY(),dr.getX())+orsa::twopi(),orsa::twopi());
	  const orsa::Angle dec = asin(dr.getZ()/dr.length());
	  orsa::print(orbit);
	  orsa::print(ra);
	  orsa::print(dec);
	  const double minArc = skyCoverage->minDistance(dr.normalized(),true);
	  double V;
	  const bool inField = skyCoverage->get(dr.normalized(),V,true);
	  ORSA_DEBUG("object [%s] present, min distance: %.2f [deg] = min distance: %.0f [arcsec]  in-field: %i",
		     _data[_data.size()-1].designation.getRef().c_str(),
		     orsa::radToDeg()*minArc,
		     orsa::radToArcsec()*minArc,
		     inField);
	}
      }
    }
    
    double V; // apparent magnitude of field
    if (skyCoverage->get(dr.normalized(),V)) {
      // all good, keep the object
      
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
  
  orsaSPICE::SPICE::instance()->loadKernel("de405.bsp");
  
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
  if (strlen(compactDate.c_str())==7) {
    // seven character date format
    std::string year,dayOfYear;
    year.assign(compactDate,0,4);
    dayOfYear.assign(compactDate,4,3);
#warning ADD FRACTION OF DAY FOR LOCAL MIDNIGHT EXACT EPOCH
    ORSA_DEBUG("ADD FRACTION OF DAY FOR LOCAL MIDNIGHT EXACT EPOCH");
    //
    epoch = orsaSolarSystem::gregorTime(atoi(year.c_str()),
					1,
					atoi(dayOfYear.c_str())+1.0-observatory.lon.getRef()/orsa::twopi());
    orsa::print(epoch);
  }
  
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
  
  osg::ref_ptr<orsa::BodyGroup> bg = new orsa::BodyGroup;
  
  // SUN
  osg::ref_ptr<orsa::Body> sun   = SPICEBody("SUN",orsaSolarSystem::Data::MSun());
  bg->addBody(sun.get());
  
  // EARTH
  osg::ref_ptr<orsa::Body> earth = SPICEBody("EARTH",orsaSolarSystem::Data::MEarth());
  bg->addBody(earth.get());
  
  // MOON
  osg::ref_ptr<orsa::Body> moon  = SPICEBody("MOON",orsaSolarSystem::Data::MMoon());
  bg->addBody(moon.get());
  
  orsa::Vector r;
  bg->getInterpolatedPosition(r,sun.get(),epoch);
  const orsa::Vector sunPosition = r;
  obsPosCB->getPosition(r,obsCode,epoch);
  const orsa::Vector obsPosition = r;
  
  osg::ref_ptr<CustomMPCObservationsFile> obsFile =
    new CustomMPCObservationsFile;
#warning CHECK HERE...
  obsFile->select_startEpoch = epoch - orsa::Time(0,12,0,0,0);
  obsFile->select_stopEpoch  = epoch + orsa::Time(0,12,0,0,0);
  ORSA_DEBUG("select start/stop:");
  orsa::print(obsFile->select_startEpoch.getRef());
  orsa::print(obsFile->select_stopEpoch.getRef());
  obsFile->select_obsCode = obsCode;
  obsFile->setFileName("mpu.arc.gz");
  // obsFile->setFileName("mpu.2007010.obs");
  obsFile->read();
  ORSA_DEBUG("selected observations: %i",obsFile->_data.size());
  
  {
    unsigned int inField=0,inFieldCandidates=0;
    double V;
    orsaSolarSystem::OpticalObservation * obs;
    for (unsigned int k=0; k<obsFile->_data.size(); ++k) {
      obs = dynamic_cast<orsaSolarSystem::OpticalObservation *> (obsFile->_data[k].get());
      if (obs) {
	if (!obs->mag.isSet()) {
	  // ORSA_DEBUG("mag not set, skipping");
	  continue;
	}
	++inFieldCandidates;
	if (skyCoverage->get(obs->ra.getRef(),
			     obs->dec.getRef(),
			     V)) {	
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
  // orbitFile->setFileName("NEA.DAT");
  orbitFile->read();
  ORSA_DEBUG("selected orbits: %i   observed: %i",
	     orbitFile->_data.size(),
	     orbitFile->observed);
  
  if (0) {
    // dump lists
    ORSA_DEBUG("--DUMP-OBS---");
    {
      orsaSolarSystem::OpticalObservation * obs;
      for (unsigned int k=0; k<obsFile->_data.size(); ++k) {
	obs = dynamic_cast<orsaSolarSystem::OpticalObservation *> (obsFile->_data[k].get());
	if (obs) {
	  if (obs->number.isSet()) {
	    ORSA_DEBUG("(%i)",obs->number.getRef());
	  } else if (obs->designation.isSet()) {
	    ORSA_DEBUG("[%s]",obs->designation.getRef().c_str());
	  }
	}
      }
    }
    ORSA_DEBUG("--DUMP-ORB---");
    {
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
  
#warning COMPLETE WRITING OF obs.dat FILE
  // write obs.dat file
  {
    FILE * fp = fopen("obs.dat","w");
    
    fclose(fp);
  }

  return 0;
}
