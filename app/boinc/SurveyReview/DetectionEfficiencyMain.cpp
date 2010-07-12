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
#include <orsaSolarSystem/obleq.h>
#include <orsaSolarSystem/orbit.h>
#include <orsaSolarSystem/print.h>

#include <orsaUtil/observatory.h>

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
            // ORSA_DEBUG("lines processed: %i   selected: %i",processedLines,_data.size());
        }
        return retVal;
    }
public:
    orsaSolarSystem::OpticalObservation * getOpticalObservationNearEpoch(const orsa::Time & epoch) const {
        if (_data.size()==0) return 0;
        orsaSolarSystem::OpticalObservation * obs_near_epoch = 0;
        orsaSolarSystem::OpticalObservation * obs;
        for (unsigned int k=0; k<_data.size(); ++k) {
            obs = dynamic_cast<orsaSolarSystem::OpticalObservation *> (_data[k].get());
            if (obs) {
                if (obs_near_epoch != 0) {
                    if (fabs((obs->epoch.getRef()-epoch).get_d()) < fabs((obs_near_epoch->epoch.getRef()-epoch).get_d())) {
                        obs_near_epoch = obs;
                    }
                } else {
                    obs_near_epoch = obs;
                }
                if (obs_near_epoch->epoch.getRef() == epoch) {
                    break;
                }
            }
        }
        /* if (obs_near_epoch != 0) {
           ORSA_DEBUG("final epoch offset: %.3f [hours]",orsa::FromUnits((obs_near_epoch->epoch.getRef()-epoch).get_d(),orsa::Unit::HOUR,-1));
           }
        */
        return obs_near_epoch;
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
    // absolute (ecliptic) positions at skyCoverage->epoch
    /* CustomMPCAsteroidFile() : orsaInputOutput::MPCAsteroidFile() {
       processedLines=0;
       observed=0;
       }
    */
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
                        /* 
                           if (!obs->mag.isSet()) {
                           // ORSA_DEBUG("mag not set, skipping");
                           continue;
                           }
                        */
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
                        /* 
                           if (!obs->mag.isSet()) {
                           // ORSA_DEBUG("mag not set, skipping");
                           continue;
                           }
                        */
                        if (obs->designation.isSet() && _data[_data.size()-1].designation.isSet()) {
                            // ORSA_DEBUG("obs: [%s]   orb: [%s]",obs->designation.getRef().c_str(),_data[_data.size()-1].designation.getRef().c_str());
                            if (obs->designation.getRef() == _data[_data.size()-1].designation.getRef()) {
                                ++observed;
                                break;
                            }
                        }
                        if (obs->number.isSet() && _data[_data.size()-1].number.isSet()) {
                            // ORSA_DEBUG("obs: [%i]   orb: [%i]",obs->number.getRef(),_data[_data.size()-1].number.getRef());
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
    osg::ref_ptr<orsaUtil::StandardObservatoryPositionCallback> obsPosCB;
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
    
    // change random seed...
    // osg::ref_ptr<orsa::RNG> rnd = new orsa::RNG(getpid());
    osg::ref_ptr<orsa::RNG> rnd = new orsa::RNG(85719);
    
    const std::string basename = SkyCoverage::basename(argv[1]);
    
    char allEtaFilename[1024];
    sprintf(allEtaFilename,"%s.allEta.dat",basename.c_str());
    
    orsaSPICE::SPICE::instance()->loadKernel("de405.bsp");
    
    osg::ref_ptr<orsaInputOutput::MPCObsCodeFile> obsCodeFile = new orsaInputOutput::MPCObsCodeFile;
    obsCodeFile->setFileName("obscode.dat");
    obsCodeFile->read();
  
    osg::ref_ptr<orsaUtil::StandardObservatoryPositionCallback> obsPosCB =
        new orsaUtil::StandardObservatoryPositionCallback(obsCodeFile.get());
  
    std::string obsCode;
    orsa::Time epoch;
    int year;
    int dayOfYear;
    if (!SkyCoverage::processFilename(argv[1],
                                      obsCodeFile.get(),
                                      obsCode,
                                      epoch,
                                      year,
                                      dayOfYear)) {
        ORSA_DEBUG("problems...");
        exit(0);
    }
  
    osg::ref_ptr<SkyCoverageFile> skyCoverageFile = new SkyCoverageFile;
    skyCoverageFile->setFileName(argv[1]);
    skyCoverageFile->read();
    
    osg::ref_ptr<SkyCoverage> skyCoverage = skyCoverageFile->_data;
    
    if (0) {
        // test skycoverage
        // const int randomSeed=35092;
        // osg::ref_ptr<orsa::RNG> rnd = new orsa::RNG(randomSeed); 
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
    osg::ref_ptr<orsa::Body> earth = SPICEBody("EARTH",orsaSolarSystem::Data::MEarth());
    bg->addBody(earth.get());
    
    // MOON
    osg::ref_ptr<orsa::Body> moon  = SPICEBody("MOON",orsaSolarSystem::Data::MMoon());
    bg->addBody(moon.get());
    
    // const orsaSolarSystem::Observatory observatory = obsPosCB->getObservatory(obsCode);
    // const double observatoryLatitude = observatory.latitude();
  
    // earth north pole
    const orsa::Vector northPole = (orsaSolarSystem::equatorialToEcliptic()*orsa::Vector(0,0,1)).normalized();
  
    osg::ref_ptr<CustomMPCObservationsFile> obsFile =
        new CustomMPCObservationsFile;
    {
        obsFile->select_startEpoch = epoch - orsa::Time(0,12,0,0,0);
        obsFile->select_stopEpoch  = epoch + orsa::Time(0,12,0,0,0);
        /* ORSA_DEBUG("select start/stop:");
           orsaSolarSystem::print(obsFile->select_startEpoch.getRef());
           orsaSolarSystem::print(obsFile->select_stopEpoch.getRef());
        */
        //
        obsFile->select_obsCode = obsCode;
        //
        /* obsFile->setFileName("NumObs.txt.gz");
           obsFile->read();
           obsFile->setFileName("UnnObs.txt.gz");
           obsFile->read();
        */
        char filename[1024];
        sprintf(filename,"%i.obs.gz",year);
        obsFile->setFileName(filename);
        obsFile->read();
        //
        ORSA_DEBUG("total selected observations: %i",obsFile->_data.size());
    }

    if (obsFile->_data.size()==0) {
        // write dummy file and exit
        FILE * fp = fopen(allEtaFilename,"w");
        fclose(fp);
        exit(0);
    }
    
    osg::ref_ptr<orsaSolarSystem::OpticalObservation> obs_near_epoch =
        obsFile->getOpticalObservationNearEpoch(epoch);
    epoch = obs_near_epoch->epoch.getRef();
    /* {
       orsaSolarSystem::OpticalObservation * obs;
       for (unsigned int k=0; k<obsFile->_data.size(); ++k) {
       obs = dynamic_cast<orsaSolarSystem::OpticalObservation *> (obsFile->_data[k].get());
       if (obs) {
       if (obs_near_epoch.get()) {
       if (fabs((obs->epoch.getRef()-epoch).get_d()) < fabs((obs_near_epoch->epoch.getRef()-epoch).get_d())) {
       obs_near_epoch = obs;
       }
       } else {
       obs_near_epoch = obs;
       }
       }
       }
       }
    */
    
    orsa::Vector r;
    
    bg->getInterpolatedPosition(r,sun.get(),epoch);
    const orsa::Vector sunPosition = r;
    // bg->getInterpolatedPosition(r,sun.get(),epoch+dt);
    // const orsa::Vector sunPosition_dt = r;
    // bg->getInterpolatedPosition(r,moon.get(),epoch);
    // const orsa::Vector moonPosition = r;
    // bg->getInterpolatedPosition(r,moon.get(),epoch+dt);
    // const orsa::Vector moonPosition_dt = r;
    obsPosCB->getPosition(r,obs_near_epoch.get());
    const orsa::Vector obsPosition = r;
    // obsPosCB->getPosition(r,obsCode,epoch+dt);
    // const orsa::Vector obsPosition_dt = r;
    
    
    {
        unsigned int inField=0,inFieldCandidates=0;
        orsaSolarSystem::OpticalObservation * obs;
        for (unsigned int k=0; k<obsFile->_data.size(); ++k) {
            obs = dynamic_cast<orsaSolarSystem::OpticalObservation *> (obsFile->_data[k].get());
            if (obs) {
                /* 
                   if (!obs->mag.isSet()) {
                   // ORSA_DEBUG("mag not set, skipping");
                   continue;
                   }
                */
                ++inFieldCandidates;
                if (skyCoverage->fastGet(obs->ra.getRef(),
                                         obs->dec.getRef())) {	
                    ++inField;
                    skyCoverage->insertFieldTime(obs->epoch.getRef(),
                                                 obs->ra.getRef(),
                                                 obs->dec.getRef());
                } 
            }
        }
        // debug
        ORSA_DEBUG("inField success: %i/%i",inField,inFieldCandidates);
        // this ratio should be smaller than 1.0, because the field declared
        // in the sky coverage files is smaller than the real data CCD field.
    }
    
    {
        char filename[1024];
        sprintf(filename,"%s.fieldTime.dat",basename.c_str());
        ORSA_DEBUG("writing file: [%s]",filename);
        skyCoverage->writeFieldTimeFile(filename);
    }
    
    osg::ref_ptr<CustomMPCAsteroidFile> orbitFile = new CustomMPCAsteroidFile(sunPosition,obsPosition);
    // osg::ref_ptr<CustomMPCAsteroidFile> orbitFile = new CustomMPCAsteroidFile;
    orbitFile->skyCoverage = skyCoverage.get();
    orbitFile->obsPosCB = obsPosCB.get();
    orbitFile->obsFile = obsFile.get();
    orbitFile->setFileName("MPCORB.DAT.gz");
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
  
    // estimate of total time actively spent observing during the night
    // measured as multiple of dt, if there are observations recorded within t and t+dt
    orsa::Time activeTime(0);
    {
        orsa::Time           t = epoch - orsa::Time(0,12,0,0,0);
        const orsa::Time tStop = epoch + orsa::Time(0,12,0,0,0);
        const orsa::Time    dt = orsa::Time(0,0,15,0,0); // small for accurate estimate od activeTime, but not too small
        while (t < tStop) {
            for (unsigned int kobs=0; kobs<obsFile->_data.size(); ++kobs) {
                const orsaSolarSystem::Observation * obs = obsFile->_data[kobs].get();
                if (obs) {
                    if (obs->epoch.isSet()) {
                        if ( (obs->epoch.getRef()>=t) && 
                             (obs->epoch.getRef()<=t+dt) ) {
                            activeTime += dt;
                            break;
                        }
                    }
                }
            }
            t += dt;
        }
    }
    //
    ORSA_DEBUG("activeTime: %g [hour]",orsa::FromUnits(activeTime.get_d(),orsa::Unit::HOUR,-1));
  
    // start to work at efficiency
    std::vector<EfficiencyData> etaData;
    for (unsigned int korb=0; korb<orbitFile->_data.size(); ++korb) {
    
        bool observed=false;
        bool discovered=false;
        {
            // osg::ref_ptr< orsa::Statistic<double> > epochStat_JD = new orsa::Statistic<double>;
            std::vector<orsa::Time> epochVec;
            orsaSolarSystem::OpticalObservation * obs;
            for (unsigned int kobs=0; kobs<obsFile->_data.size(); ++kobs) {
                obs = dynamic_cast<orsaSolarSystem::OpticalObservation *> (obsFile->_data[kobs].get());
                if (obs) {
                    /* 
                       if (!obs->mag.isSet()) {
                       // ORSA_DEBUG("mag not set, skipping");
                       continue;
                       }
                    */
                    if (obs->designation.isSet() && orbitFile->_data[korb].designation.isSet()) {
                        if (obs->designation.getRef() == orbitFile->_data[korb].designation.getRef()) {
                            observed=true;
                            if (obs->discovery.isSet()) {
                                if (obs->discovery.getRef()) {
                                    discovered=true;
                                }
                            }  
                            // epoch=obs->epoch.getRef();
                            // epochStat_JD->insert(orsaSolarSystem::timeToJulian(obs->epoch.getRef()));
                            epochVec.push_back(obs->epoch.getRef());
                            // break; // no break, because it can skip the discovery asterisk and also because we want to include all relevant epochs
                        }
                    }
                    if (obs->number.isSet() && orbitFile->_data[korb].number.isSet()) {
                        if (obs->number.getRef() == orbitFile->_data[korb].number.getRef()) {
                            observed=true;
                            if (obs->discovery.isSet()) {
                                if (obs->discovery.getRef()) {
                                    discovered=true;
                                }
                            }  
                            // epoch=obs->epoch.getRef();
                            // epochStat_JD->insert(orsaSolarSystem::timeToJulian(obs->epoch.getRef()));
                            epochVec.push_back(obs->epoch.getRef());
                            // break; // no break, because it can skip the discovery asterisk and also because we want to include all relevant epochs
                        }
                    }
                }	  
            }
            
            // pick one, better than averaging to avoid time clustering
            if (epochVec.size()>0) {
                epoch = epochVec[rnd->gsl_rng_uniform_int(epochVec.size())];
            }
        }
        
        bool epochFromField=false;
        if (!observed) {
      
            // try to retrieve epoch from field
      
            orsaSolarSystem::OrbitWithEpoch orbit = orbitFile->_data[korb].orbit.getRef();
            const double orbitPeriod = orbit.period();
            const double original_M  = orbit.M;
            //
            orsa::Vector r;
            //
            orbit.M = original_M + fmod(orsa::twopi() * (epoch-orbit.epoch.getRef()).get_d() / orbitPeriod, orsa::twopi());
            orbit.relativePosition(r);
            const orsa::Vector orbitPosition = r + sunPosition;
            // restore, important!
            orbit.M = original_M;
            orsa::Vector dr = (orbitPosition - obsPosition).normalized();
            if (skyCoverage->pickFieldTime(epoch,dr,rnd.get())) {
                epochFromField=true;
            } else {
                // ORSA_DEBUG("problems retrieving epoch from field");
            }
        }
    
        orsa::Vector r;
        
        osg::ref_ptr<orsaSolarSystem::OpticalObservation> obs_near_epoch =
            obsFile->getOpticalObservationNearEpoch(epoch);
        
        osg::ref_ptr<orsaSolarSystem::OpticalObservation> obs_near_epoch_dt =
            obsFile->getOpticalObservationNearEpoch(epoch+orsa::Time(0,1,0,0,0));
        
        epoch = obs_near_epoch->epoch.getRef();
        
        const orsa::Time epoch_dt = obs_near_epoch_dt->epoch.getRef();

        if (epoch == epoch_dt) {
            ORSA_DEBUG("problems: null dt");
            continue;
        }
        
        // replace vectors for current epoch
        bg->getInterpolatedPosition(r,sun.get(),epoch);
        const orsa::Vector sunPosition = r;
        bg->getInterpolatedPosition(r,sun.get(),epoch_dt);
        const orsa::Vector sunPosition_dt = r;
        bg->getInterpolatedPosition(r,earth.get(),epoch);
        const orsa::Vector earthPosition = r;
        bg->getInterpolatedPosition(r,earth.get(),epoch_dt);
        const orsa::Vector earthPosition_dt = r;
        bg->getInterpolatedPosition(r,moon.get(),epoch);
        const orsa::Vector moonPosition = r;
        bg->getInterpolatedPosition(r,moon.get(),epoch_dt);
        const orsa::Vector moonPosition_dt = r;
        obsPosCB->getPosition(r,obs_near_epoch);
        const orsa::Vector obsPosition = r;
        obsPosCB->getPosition(r,obs_near_epoch_dt);
        const orsa::Vector obsPosition_dt = r;
        
        orsaSolarSystem::OrbitWithEpoch orbit = orbitFile->_data[korb].orbit.getRef();
        const double orbitPeriod = orbit.period();
        const double original_M  = orbit.M;
        //
        orbit.M = original_M + fmod(orsa::twopi() * (epoch-orbit.epoch.getRef()).get_d() / orbitPeriod, orsa::twopi());
        orbit.relativePosition(r);
        const orsa::Vector orbitPosition = r + sunPosition;
        // now at t+dt
        orbit.M = original_M + fmod(orsa::twopi() * (epoch_dt-orbit.epoch.getRef()).get_d() / orbitPeriod, orsa::twopi());
        orsa::Vector r_dt;
        orbit.relativePosition(r_dt);
        const orsa::Vector orbitPosition_dt = r_dt + sunPosition_dt;
        // restore, important!
        orbit.M = original_M;
        // all in Ecliptic coords, as usual
        const orsa::Vector orb2obs  = obsPosition - orbitPosition;
        const orsa::Vector orb2sun  = sunPosition - orbitPosition;
        const orsa::Vector obs2sun  = sunPosition - obsPosition;
        const orsa::Vector obs2moon = moonPosition - obsPosition;
        const orsa::Vector obs2orb  = -orb2obs;
        const orsa::Vector moon2sun = sunPosition - moonPosition;
        const orsa::Vector moon2obs = -obs2moon;
        const orsa::Vector orb2obs_dt = obsPosition_dt - orbitPosition_dt;
        const orsa::Vector orb2sun_dt = sunPosition_dt - orbitPosition_dt;
        const orsa::Vector obs2orb_dt = -orb2obs_dt;
        //
        const double phaseAngle = acos(orb2sun.normalized()*orb2obs.normalized());
        const double solarElongation = acos(obs2sun.normalized()*obs2orb.normalized());
        const double lunarElongation = acos(obs2moon.normalized()*obs2orb.normalized());
        const double lunarPhase = acos(moon2sun.normalized()*moon2obs.normalized());
        //
        // airmass & azimuth
        const orsa::Vector     zenith = (obsPosition - earthPosition).normalized();
        const orsa::Vector localEast  = orsa::externalProduct(northPole,zenith).normalized();
        const orsa::Vector localNorth = orsa::externalProduct(zenith,localEast).normalized();
        const double obs2orb_zenith     =     zenith*obs2orb.normalized();
        const double obs2orb_localEast  =  localEast*obs2orb.normalized();
        const double obs2orb_localNorth = localNorth*obs2orb.normalized();
        const double zenithAngle = acos(obs2orb_zenith);
        const double airMass = ((observed||epochFromField)&&(zenithAngle<orsa::halfpi())?(1.0/cos(zenithAngle)):10.0);
        const double azimuth = fmod(orsa::twopi()+atan2(obs2orb_localEast,obs2orb_localNorth),orsa::twopi());
        // lunar altitude
        const double solarAltitude = ((observed||epochFromField)?(orsa::halfpi()-acos(zenith*obs2sun.normalized())):-orsa::halfpi());
        const double lunarAltitude = ((observed||epochFromField)?(orsa::halfpi()-acos(zenith*obs2moon.normalized())):-orsa::halfpi());
        // galactic latitude
        const orsa::Vector obs2orb_Equatorial = orsaSolarSystem::eclipticToEquatorial()*obs2orb;
        const orsa::Vector dr_equatorial = obs2orb_Equatorial.normalized();
        const double  ra = fmod(atan2(dr_equatorial.getY(),dr_equatorial.getX())+orsa::twopi(),orsa::twopi());
        const double dec = asin(dr_equatorial.getZ()/dr_equatorial.length());
        double l,b;
        orsaSolarSystem::equatorialToGalactic(l,b,ra,dec);
        // format longitude between -180 and +180 deg
        l = fmod(l+2*orsa::twopi(),orsa::twopi());
        if (l > orsa::pi()) l -= orsa::twopi();
        const double galacticLongitude = l;
        const double galacticLatitude  = b;
        // ecliptic coordinates
        const orsa::Vector dr = obs2orb.normalized();
        const double phi      = fmod(atan2(dr.getY(),dr.getX())+orsa::twopi(),orsa::twopi());
        const double theta    = asin(dr.getZ()/dr.length());
        const orsa::Vector dr_sun = obs2sun.normalized();
        const double phi_sun      = fmod(atan2(dr_sun.getY(),dr_sun.getX())+orsa::twopi(),orsa::twopi());
        const double theta_sun    = asin(dr_sun.getZ()/dr_sun.length());
        const double tmp_eclipticLongitude = fmod(phi-phi_sun+orsa::twopi(),orsa::twopi());
        const double eclipticLongitude = (tmp_eclipticLongitude>orsa::pi()) ? (tmp_eclipticLongitude-orsa::twopi()) : (tmp_eclipticLongitude);
        const double eclipticLatitude  = theta-theta_sun;
        
        // ORSA_DEBUG("ra: %g  dec: %g",ra*orsa::radToDeg()/15.0,dec*orsa::radToDeg());
        
        EfficiencyData ed;
        ed.H = orbitFile->_data[korb].H.getRef();
        if (orbitFile->_data[korb].number.isSet()) {
            ed.number = orbitFile->_data[korb].number.getRef();
        }
        if (orbitFile->_data[korb].designation.isSet()) {
            ed.designation = orbitFile->_data[korb].designation.getRef();
        }
        ed.V = apparentMagnitude(orbitFile->_data[korb].H.getRef(),
                                 orbitFile->_data[korb].G.getRef(),
                                 phaseAngle,
                                 orb2obs.length(),
                                 orb2sun.length()); 
        ed.apparentVelocity = acos(obs2orb_dt.normalized()*obs2orb.normalized())/fabs((epoch_dt-epoch).get_d());
        ed.solarElongation = solarElongation;
        ed.lunarElongation = lunarElongation;
        ed.solarAltitude = solarAltitude;
        ed.lunarAltitude = lunarAltitude;
        ed.lunarPhase = lunarPhase;
        ed.airMass = airMass;
        ed.azimuth = azimuth;
        ed.galacticLongitude = galacticLongitude;
        ed.galacticLatitude  = galacticLatitude;
        ed.eclipticLongitude = eclipticLongitude;
        ed.eclipticLatitude  = eclipticLatitude;
        ed.activeTime = activeTime.get_d();
        ed.epochFromField = epochFromField;
        ed.observed = observed;
        ed.discovered = discovered;
        etaData.push_back(ed);
    }
    //    
    {
        writeEfficiencyDataFile(etaData,allEtaFilename);
    }
    
    return 0;
}
