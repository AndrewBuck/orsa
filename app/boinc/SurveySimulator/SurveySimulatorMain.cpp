#include "parameter.h"
#include "telescope.h"
#include "SurveySimulator.h"

#include <list>

#include <fcntl.h>

#include <orsa/body.h>
#include <orsa/bodygroup.h>
#include <orsa/double.h>
#include <orsa/matrix.h>
#include <orsa/unit.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <orsaSPICE/spice.h>
#include <orsaSPICE/spiceBodyAttitudeCallback.h>
#include <orsaSPICE/spiceBodyPosVelCallback.h>

#include <orsaSolarSystem/datetime.h>
#include <orsaSolarSystem/gmst.h>
#include <orsaSolarSystem/obleq.h>
#include <orsaSolarSystem/observatory.h>

// BOINC API
#include <boinc_api.h>
#include <filesys.h>

// CSPICE prototypes and definitions.      
#include "SpiceUsr.h"

void writeCheckpoint(const orsa::Time    & t,
		     const TelescopeList & telescopeList) {
  
  if (boinc_is_standalone()) {
    ORSA_DEBUG("muSec: %Zi   JD: %.5Ff",
	       t.getMuSec().get_mpz_t(),
	       orsaSolarSystem::timeToJulian(t).get_mpf_t());
  }
  
  {
    FILE * fp_checkpoint = boinc_fopen("checkpoint.dat","w");
    if (fp_checkpoint) {
      gmp_fprintf(fp_checkpoint,
		  "%Zi"
		  MODEOL,
		  t.getMuSec().get_mpz_t());
      fclose(fp_checkpoint);
    }
  }
  
  { 
    TelescopeList::const_iterator tl_ckpt = telescopeList.begin();
    while (tl_ckpt != telescopeList.end()) {
      (*tl_ckpt)->writeCheckpoint();
      ++tl_ckpt;
    }
  }
}


int main() {
  
  // tests
  /* 
     {
     osg::ref_ptr<orsa::RNG> rnd = new orsa::RNG(62342);
     for (unsigned int l=0; l<100000; ++l) {
     std::cout << rnd->gsl_rng_uniform() << std::endl;
     }
     exit(0);
     }
  */
  //
  /* {
     osg::ref_ptr<Telescope> telescope = new OppositionTelescope(orsaSolarSystem::J2000(),
     orsaSolarSystem::J2000(),
     1233,
     20.0,
     2.0,
     60,
     45,
     45,
     10,
     200,
     3,
     "703",
     "test");
     orsa::Double V = 10.0;
     while (V <= 25.0) {
     std::cout << V << " " <<  telescope->detectionProbability(V) << std::endl;
     V += 0.02134;
     }
     exit(0);
     }
  */
  
  orsa::Debug::instance()->initTimer();
  
  // mpf_set_default_prec(64);
  
  boinc_init();
  
  ORSA_DEBUG("process ID: %i",getpid());
  
  bool writeNEOLogFile=true;
  //
  /* 
     if (boinc_is_standalone()) {
     // only for testing
     writeNEOLogFile=false;
     } else {
     }
  */
  
  // test
  /* 
     if (boinc_is_standalone()) {
     ORSA_DEBUG("--- STANDALONE ---");
     } else {
     ORSA_DEBUG("--- NOT STANDALONE ---");  
     }
  */
  
  std::string resolvedFileName;
  
  // spice error file (should resolve filename?)
  SpiceChar spiceError[1024];
  sprintf(spiceError,"cspice.error");
  ::errdev_c("SET", 4096, spiceError);
  
  osg::ref_ptr<ParameterFile> parFile = new ParameterFile;
  //
  boinc_resolve_filename_s("param.conf",resolvedFileName);
  parFile->mainRead(resolvedFileName);
  
  // check values
  parFile->print();
  
  const orsa::Time genEpoch = orsaSolarSystem::J2000();
  
  TelescopeList tmpTL;
  boinc_resolve_filename_s("telescope.conf",resolvedFileName);
  if (!readTelescope(tmpTL,resolvedFileName)) {
    ORSA_ERROR("problems while reading telescope configuration file");
    boinc_finish(1);
  }
  //
  const TelescopeList telescopeList = tmpTL;
  
  /* 
     {
     // debug
     TelescopeList::const_iterator tl_dbg = telescopeList.begin();
     while (tl_dbg != telescopeList.end()) {
     ORSA_DEBUG("telescope: [%s]   obscode: [%s]",
     (*tl_dbg)->name.c_str(),
     (*tl_dbg)->obscode.c_str());
     ++tl_dbg;
     }
     }
  */
  
  orsa::Time endSimulation = genEpoch;
  { 
    TelescopeList::const_iterator tl_epoch = telescopeList.begin();
    while (tl_epoch != telescopeList.end()) {
      endSimulation = std::max((*tl_epoch)->tStop,endSimulation);
      ++tl_epoch;
    }
  }
  
  // const unsigned int Nreal      = parFile->numRealNEO.getRef();
  // const unsigned int Nsynthetic = parFile->numSyntheticNEO.getRef();
  
  // const int      realRandomSeed = parFile->randomSeedRealNEO.getRef();
  // const int syntheticRandomSeed = parFile->randomSeedSyntheticNEO.getRef();
  
  const orsa::Double detectionProbabilityThreshold = parFile->detectionProbabilityThreshold.getRef();
  
  // tpList tp;
  
  // re-observe same sky field, if possible, after this (long) time
  // const orsa::Time recycleTime = orsa::Time(14,0,0,0,0);
  //
  // const orsa::Time recycleTime = orsa::Time(parFile->recycleTime_DAY.getRef(),0,0,0,0);
  
  //
  // const orsa::Time dutyCycle = orsa::Time(0,0,0,parFile->dutyCycle_SEC.getRef(),0);
  // const int        dutyCycleMultiplicity = parFile->dutyCycleMultiplicity.getRef();
  
  // const orsa::Time effectiveDutyCycle = dutyCycle * dutyCycleMultiplicity;
  
  const bool cacheNEOVector = parFile->cacheNEOVector.getRef();
  
  if (cacheNEOVector) {
    ORSA_DEBUG("cacheNEOVector ENABLED");
  } else {
    ORSA_DEBUG("cacheNEOVector DISABLED");
  }
  
  // SPICE...
  orsaSPICE::SPICE::instance()->setDefaultObserver("SSB");
  {
    boinc_resolve_filename_s("de405.bsp",resolvedFileName);
    FILE * fp_dummy = boinc_fopen(resolvedFileName.c_str(),"r"); // little trick, is it helping?
    orsaSPICE::SPICE::instance()->loadKernel(resolvedFileName);
    fclose(fp_dummy);
  }
  
  // typedef std::list< osg::ref_ptr<NEO> > NEOList;
  //
  NEOList realNEO, syntheticNEO;
  //
  std::map<int, orsa::Cache<unsigned int> > generatedMap;
  //
  { 
    boinc_resolve_filename_s("realNEO.gen",resolvedFileName);
    FILE * fp_gen =  boinc_fopen(resolvedFileName.c_str(),"r");
    if (!fp_gen) {
      ORSA_ERROR("cannot open file");
      return false;
    } 
    char line[1024];
    unsigned int num;
    int randomSeed;
    while (fgets(line,1024,fp_gen) != 0) {
      if (2 == gmp_sscanf(line,"%d %d",
			  &num,
			  &randomSeed)) {
	if (generatedMap[randomSeed].isSet()) {
	  if (generatedMap[randomSeed].getRef() == num) {
	    ORSA_DEBUG("duplicate entry: already generated %i realNEOs using randomSeed: %i",
		       generatedMap[randomSeed].getRef(),
		       randomSeed);
	  } else {
	    ORSA_DEBUG("CONFLICT entry: already generated %i realNEOs using randomSeed: %i [now requested %i]",
		       generatedMap[randomSeed].getRef(),
		       randomSeed,
		       num);
	  }
	  // in either case, ignore this entry
	  continue;
	}
	generatedMap[randomSeed] = num;
	ORSA_DEBUG("generating %i realNEOs, randomSeed: %i",
		   num,
		   randomSeed);
	osg::ref_ptr<RealNEO::NEOFactory> realNEOFactory = new RealNEO::NEOFactory(parFile->orbit_a_AU_min.getRef(),
										   parFile->orbit_a_AU_max.getRef(),
										   parFile->orbit_e_min.getRef(),
										   parFile->orbit_e_max.getRef(),
										   parFile->orbit_i_DEG_min.getRef(),
										   parFile->orbit_i_DEG_max.getRef(),
										   randomSeed);
	for (unsigned int k=0; k<num; ++k) {
	  // ORSA_DEBUG("real, k: %i",k);
	  realNEO.push_back(realNEOFactory->sampleNEO(genEpoch,
						      parFile->Hmax.getRef(),
						      parFile->Ha.getRef(),
						      detectionProbabilityThreshold));
	}
      }
    }
    fclose(fp_gen);
  }
  //
  { 
    boinc_resolve_filename_s("syntheticNEO.gen",resolvedFileName);
    FILE * fp_gen = boinc_fopen(resolvedFileName.c_str(),"r");
    if (!fp_gen) {
      ORSA_ERROR("cannot open file");
      return false;
    } 
    char line[1024];
    unsigned int num;
    int randomSeed;
    while (fgets(line,1024,fp_gen) != 0) {
      if (2 == gmp_sscanf(line,"%d %d",
			  &num,
			  &randomSeed)) {
	if (generatedMap[randomSeed].isSet()) {
	  if (generatedMap[randomSeed].getRef() == num) {
	    ORSA_DEBUG("duplicate entry: already generated %i realNEOs, randomSeed: %i",
		       generatedMap[randomSeed].getRef(),
		       randomSeed);
	  } else {
	    ORSA_DEBUG("CONFLICT entry: already generated %i realNEOs, randomSeed: %i",
		       generatedMap[randomSeed].getRef(),
		       randomSeed);
	  }
	}
	generatedMap[randomSeed] = num;
	ORSA_DEBUG("generating %i syntheticNEOs, randomSeed: %i",
		   num,
		   randomSeed);
	osg::ref_ptr<SyntheticNEO::NEOFactory> syntheticNEOFactory = new SyntheticNEO::NEOFactory(parFile->orbit_a_AU_min.getRef(),
												  parFile->orbit_a_AU_max.getRef(),
												  parFile->orbit_e_min.getRef(),
												  parFile->orbit_e_max.getRef(),
												  parFile->orbit_i_DEG_min.getRef(),
												  parFile->orbit_i_DEG_max.getRef(),
												  randomSeed);
	for (unsigned int k=0; k<num; ++k) {
	  // ORSA_DEBUG("synthetic, k: %i",k);
	  syntheticNEO.push_back(syntheticNEOFactory->sampleNEO(genEpoch));
	}
      }
    }
    fclose(fp_gen);
  }
  //
  realNEO.sort(NEOListSortPredicate);
  syntheticNEO.sort(NEOListSortPredicate);
  //
  // lookup tables
  // random seed -> id -> NEO
  typedef std::map < int, std::map < unsigned int, osg::ref_ptr<NEO> > > NEOLookupMap;
  //
  NEOLookupMap      realNEOLookupMap;
  NEOLookupMap syntheticNEOLookupMap;
  //
  {
    NEOList::const_iterator it = realNEO.begin();
    while (it != realNEO.end()) {
      realNEOLookupMap[(*it)->randomSeed][(*it)->id] = (*it).get();
      ++it;
    }
  }
  //
  {
    NEOList::const_iterator it = syntheticNEO.begin();
    while (it != syntheticNEO.end()) {
      syntheticNEOLookupMap[(*it)->randomSeed][(*it)->id] = (*it).get();
      ++it;
    }
  }
  //
  {
    boinc_resolve_filename_s("NEO.real.log",resolvedFileName);
    FILE * fp_log = boinc_fopen(resolvedFileName.c_str(),"r");
    if (fp_log) {
      mpz_class    muSec;
      int          randomSeed;
      unsigned int id;
      orsa::Double p;
      char telescopeName[1024];
      char line[1024];
      osg::ref_ptr<NEO> neo;
      while (fgets(line,1024,fp_log) != 0) {
	if (5 == gmp_sscanf(line,"%Zi %d %d %Ff %s",
			    muSec.get_mpz_t(),
			    &randomSeed,
			    &id,
			    p.get_mpf_t(),
			    telescopeName)) {
	  /* 
	     if (boinc_is_standalone()) {
	     ORSA_DEBUG("log entry for real NEO, randomSeed: %i   id: %i",
	     randomSeed,
	     id);
	     }
	  */
	  bool matchFound=false;
	  if (neo.get()) {
	    if ((neo->randomSeed == randomSeed) && 
		(neo->id == id)) {
	      // neo pointer is already the correct one
	      // ORSA_DEBUG("recycling...");
	    } else {
	      neo = realNEOLookupMap[randomSeed][id].get();
	    }
	  } else {
	    neo = realNEOLookupMap[randomSeed][id].get();
	  }
	  if (neo.get()) {
	    RealNEO::LogEntry * re = new RealNEO::LogEntry;
	    //
	    re->p = p;
	    //
	    re->telescopeName = telescopeName;
	    //
	    // neo->logInsert(orsa::Time(muSec),re,writeNEOLogFile);
	    neo->logInsert(orsa::Time(muSec),re,false); // don't need to write at this time...
	    matchFound=true;
	    // ORSA_DEBUG("match found");
	  }
	  
	  /* 
	     NEOList::iterator it = realNEO.begin();
	     while (it != realNEO.end()) {
	     if (((*it)->id == id) && 
	     ((*it)->randomSeed == randomSeed)) {
	     RealNEO::LogEntry * re = new RealNEO::LogEntry;
	     //
	     re->p = p;
	     //
	     re->telescopeName = telescopeName;
	     //
	     (*it)->logInsert(orsa::Time(muSec),re);
	     matchFound=true;
	     ORSA_DEBUG("match found");
	     break;
	     }	    
	     ++it;
	     }
	  */
	  if (!matchFound) {
	    ORSA_ERROR("no matching real NEO found, randomSeed: %i   id: %i",
		       randomSeed,
		       id);
	  }
	}
      }
      fclose(fp_log);
    }
  }
  //
  {
    boinc_resolve_filename_s("NEO.synthetic.log",resolvedFileName);
    FILE * fp_log =  boinc_fopen(resolvedFileName.c_str(),"r");
    if (fp_log) {
      mpz_class    muSec;
      int          randomSeed;
      unsigned int id;
      orsa::Double neo2obs, neo2sun, phaseAngle;
      orsa::Double limitingMagnitude;
      char telescopeName[1024];
      char line[1024];
      osg::ref_ptr<NEO> neo;
      while (fgets(line,1024,fp_log) != 0) {
	if (8 == gmp_sscanf(line,"%Zi %d %d %Ff %Ff %Ff %Ff %s",
			    muSec.get_mpz_t(),
			    &randomSeed,
			    &id,
			    neo2obs.get_mpf_t(),
			    neo2sun.get_mpf_t(),
			    phaseAngle.get_mpf_t(),
			    limitingMagnitude.get_mpf_t(),
			    telescopeName)) {
	  /* 
	     if (boinc_is_standalone()) {
	     ORSA_DEBUG("log entry for synthetic NEO, randomSeed: %i   id: %i",
	     randomSeed,
	     id);
	     }
	  */
	  neo2obs     = FromUnits(neo2obs,orsa::Unit::AU);
	  neo2sun     = FromUnits(neo2sun,orsa::Unit::AU);
	  phaseAngle *= orsa::degToRad();
	  bool matchFound=false;
	  if (neo.get()) {
	    if ((neo->randomSeed == randomSeed) && 
		(neo->id == id)) {
	      // neo pointer is already the correct one
	      // ORSA_DEBUG("recycling...");
	    } else {
	      // first: call logPurify on previous neo
	      {
		SyntheticNEO * syntheticNEO = dynamic_cast<SyntheticNEO *> (neo.get());
		if (syntheticNEO) {
		  syntheticNEO->logPurify(detectionProbabilityThreshold);
		} else {
		  ORSA_DEBUG("problems");
		}
	      }
	      // then, update pointer
	      neo = syntheticNEOLookupMap[randomSeed][id].get();
	    }
	  } else {
	    neo = syntheticNEOLookupMap[randomSeed][id].get();
	  }
	  if (neo.get()) {
	    SyntheticNEO::LogEntry * se = new SyntheticNEO::LogEntry;
	    //
	    se->neo2obs    = neo2obs;
	    se->neo2sun    = neo2sun;
	    se->phaseAngle = phaseAngle;
	    //
	    se->limitingMagnitude = limitingMagnitude;
	    //
	    se->telescopeName = telescopeName;
	    //
	    // neo->logInsert(orsa::Time(muSec),se,writeNEOLogFile);
	    neo->logInsert(orsa::Time(muSec),se,false); // don't need to write at this time...
	    matchFound=true;
	    // ORSA_DEBUG("match found");
	  }	    
	  if (!matchFound) {
	    ORSA_ERROR("no matching synthetic NEO found, randomSeed: %i   id: %i",
		       randomSeed,
		       id);
	    // debug
	    /* 
	       NEOList::iterator it = syntheticNEO.begin();
	       while (it != syntheticNEO.end()) {
	       ORSA_DEBUG("id: %6i   randomSeed: %i",(*it)->id,(*it)->randomSeed);
	       ++it;
	       }
	    */
	  }
	} else {
	  // one last logPurify call
	  SyntheticNEO * syntheticNEO = dynamic_cast<SyntheticNEO *> (neo.get());
	  if (syntheticNEO) {
	    syntheticNEO->logPurify(detectionProbabilityThreshold);
	  } else {
	    ORSA_DEBUG("problems");
	  }
	}	
      }
      fclose(fp_log);
    }
  }
  
  // dump on files
  if (boinc_is_standalone()) {
    dumpSampledNEOs(     realNEO,     "realNEO.dat",detectionProbabilityThreshold);
    dumpSampledNEOs(syntheticNEO,"syntheticNEO.dat",detectionProbabilityThreshold);
  }
  
  osg::ref_ptr<orsa::BodyGroup> bg = new orsa::BodyGroup;
  
  // SUN
  osg::ref_ptr<orsa::Body> sun   = SPICEBody("SUN",
					     FromUnits(orsa::one(),orsa::Unit::MSUN));
  bg->addBody(sun.get());
  
  // EARTH
  osg::ref_ptr<orsa::Body> earth = SPICEBody("EARTH",
					     FromUnits(orsa::one(),orsa::Unit::MEARTH));
  bg->addBody(earth.get());
  
  // MOON
  osg::ref_ptr<orsa::Body> moon  = SPICEBody("MOON",orsa::zero());
  bg->addBody(moon.get());
  
  orsa::Matrix eclipticToEquatorial = orsa::Matrix::identity();
  eclipticToEquatorial.rotX(+orsaSolarSystem::obleqJ2000());
  
  orsa::Matrix equatorialToEcliptic = orsa::Matrix::identity(); 
  equatorialToEcliptic.rotX(-orsaSolarSystem::obleqJ2000());
  
  orsa::Time t = genEpoch;
  const orsa::Time tStop = endSimulation;
  //
  bool checkpointUsed = false;
  // 
  {
    // check if a checkpoint is available
    
    FILE * fp_checkpoint = boinc_fopen("checkpoint.dat","r");
    if (fp_checkpoint != 0) {
      
      ORSA_DEBUG("restoring from latest checkpoint data...");
      
      mpz_class muSec;
      gmp_fscanf(fp_checkpoint,
		 "%Zi",
		 muSec.get_mpz_t());
      fclose(fp_checkpoint);
      
      // ORSA_DEBUG("read time %Zi",muSec.get_mpz_t());
      
      // t = orsa::Time(muSec) + dt;
      // add dt later
      t = orsa::Time(muSec);
      
      /* 
	 FILE * fp_sky_rnd = fopen(sky_rnd_state_filename,"rb");
	 if (fp_sky_rnd != 0) {
	 sky_rnd->gsl_rng_fread(fp_sky_rnd);
	 fclose(fp_sky_rnd);
	 } else {
	 ORSA_ERROR("cannot open file %s resolved as %s",
	 sky_rnd_state_filename,
	 resolvedFileName.c_str());
	 boinc_finish(1);
	 }
      */
      //
      {
	TelescopeList::const_iterator tl_ckpt = telescopeList.begin();
	while (tl_ckpt != telescopeList.end()) {
	  /* 
	     ORSA_DEBUG("telescope: [%s]   obscode: [%s]   telescopeList.size(): %i",
	     (*tl_ckpt)->name.c_str(),
	     (*tl_ckpt)->obscode.c_str(),
	     telescopeList.size());
	  */
	  (*tl_ckpt)->readCheckpoint(t);
	  ++tl_ckpt;
	}
      }
      
      /* 
	 {
	 // debug
	 TelescopeList::const_iterator tl_dbg = telescopeList.begin();
	 while (tl_dbg != telescopeList.end()) {
	 ORSA_DEBUG("telescope: [%s]   obscode: [%s]",
	 (*tl_dbg)->name.c_str(),
	 (*tl_dbg)->obscode.c_str());
	 ++tl_dbg;
	 }
	 }
      */
      
      /* 
	 char line[1024];
	 char goodLine[1024];
	 off_t lastGoodFilePosition;
	 unsigned int NEO_id;
      */
      
      /* 
	 lastGoodFilePosition=lseek(fd_real, 0, SEEK_SET);
	 while (read(fd_real,line,1024) > 0) {
	 
	 const std::string lineString(line);
	 const size_t newLinePos = lineString.find_first_of('\n'); 
	 strncpy(goodLine,line,newLinePos);
	 goodLine[newLinePos] = '\0';
	 
	 // ORSA_DEBUG("goodLine: [%s]   newLinePos: %i",goodLine,newLinePos);
	 
	 if (gmp_sscanf(goodLine,"%*Ff %Zi %i",
	 muSec.get_mpz_t(),
	 &NEO_id) != 2) {
	 break;
	 }
	 
	 if (muSec <= t.getMuSec()) {
	 
	 NEOList::iterator it = realNEO.begin();
	 
	 while (it != realNEO.end()) {
	 if ((*it)->id == NEO_id) {
	 it = realNEO.erase(it);
	 break;
	 }
	 ++it;
	 }
	 
	 if (newLinePos > 0) {
	 lastGoodFilePosition = lseek(fd_real, 
	 lastGoodFilePosition+newLinePos+lineSkip, 
	 SEEK_SET);
	 } else {
	 break;
	 }
	 } else {
	 break;
	 }
	 }
	 lseek(fd_real, 
	 lastGoodFilePosition, 
	 SEEK_SET);
      */
      
      /* 
	 lastGoodFilePosition=lseek(fd_synthetic, 0, SEEK_SET);
	 while (read(fd_synthetic,line,1024) > 0) {
	 const std::string lineString(line);
	 const size_t newLinePos = lineString.find_first_of('\n'); 
	 strncpy(goodLine,line,newLinePos);
	 goodLine[newLinePos] = '\0';
	 // ORSA_DEBUG("goodLine: [%s]   newLinePos: %i",goodLine,newLinePos);
	 if (gmp_sscanf(goodLine,"%*Ff %Zi %i",
	 muSec.get_mpz_t(),
	 &NEO_id) != 2) {
	 break;
	 }
	 if (muSec <= t.getMuSec()) {
	 NEOList::iterator it = syntheticNEO.begin();
	 while (it != syntheticNEO.end()) {
	 if ((*it)->id == NEO_id) {
	 it = syntheticNEO.erase(it);
	 break;
	 }
	 ++it;
	 }
	 if (newLinePos > 0) {
	 lastGoodFilePosition = lseek(fd_synthetic, 
	 lastGoodFilePosition+newLinePos+lineSkip, 
	 SEEK_SET);
	 } else {
	 break;
	 }
	 } else {
	 break;
	 } 
	 }
	 lseek(fd_synthetic, 
	 lastGoodFilePosition, 
	 SEEK_SET);    
      */
      
      /* 
	 lastGoodFilePosition=lseek(fd_tp, 0, SEEK_SET);
	 while (read(fd_tp,line,1024) > 0) {
	 
	 const std::string lineString(line);
	 const size_t newLinePos = lineString.find_first_of('\n'); 
	 strncpy(goodLine,line,newLinePos);
	 goodLine[newLinePos] = '\0';
	 
	 // ORSA_DEBUG("goodLine: [%s]   newLinePos: %i",goodLine,newLinePos);
	 
	 if (gmp_sscanf(goodLine,"%*Ff %Zi %Ff %Ff %Ff",
	 muSec.get_mpz_t(),
	 ux.get_mpf_t(),
	 uy.get_mpf_t(),
	 uz.get_mpf_t()) != 4) {
	 break;
	 }
	 
	 // ORSA_DEBUG("read time %Zi",muSec.get_mpz_t());
	 
	 if (muSec <= t.getMuSec()) {
	 TelescopePointing localTP;
	 //
	 localTP.epoch = orsa::Time(muSec);
	 localTP.u = orsa::Vector(ux,uy,uz).normalized();
	 //
	 tp.push_back(localTP);
	 
	 if (newLinePos > 0) {
	 lastGoodFilePosition = lseek(fd_tp, 
	 lastGoodFilePosition+newLinePos+lineSkip, 
	 SEEK_SET);
	 } else {
	 // ORSA_DEBUG("break...");
	 break;
	 }	  
	 } else {
	 break;
	 }
	 }
	 lseek(fd_tp, 
	 lastGoodFilePosition, 
	 SEEK_SET);
      */
      
      // important! [no longer needed]
      // t += dt;
      
      checkpointUsed = true;
    }
    
  }

  if (0 && boinc_is_standalone()) {
    // real NEOs discovery history
    char line[1024];
    NEOList::const_iterator neo_it = realNEO.begin();
    while (neo_it != realNEO.end()) {
      RealNEO * realNEO = dynamic_cast<RealNEO *> ((*neo_it).get());
      if (realNEO) {
	const RealNEO::detectionLog & log = realNEO->getLog();
	RealNEO::detectionLog::const_iterator log_it = log.begin();
	while (log_it != log.end()) {
	  if ((*neo_it)->detectionProbabilityFromLog((*log_it).first) > detectionProbabilityThreshold) {
	    gmp_snprintf(line,1024,"%Zi %Ff %10i %10i %11.9Ff %s"
			 MODEOL,
			 (*log_it).first.getMuSec().get_mpz_t(),
			 orsaSolarSystem::timeToJulian((*log_it).first).get_mpf_t(),
			 realNEO->randomSeed,
			 realNEO->id,
			 (*neo_it)->detectionProbabilityFromLog((*log_it).first).get_mpf_t(),
			 (*log_it).second->telescopeName.getRef().c_str());
	    std::cout << line;
	    break;
	  }
	  
	  ++log_it;
	}
      }
      ++neo_it;
    }
    
    exit(0);
  }
  
  if (0 && boinc_is_standalone()) {
    // this section is for data analysis only
    
    t = orsaSolarSystem::J2000() + orsa::Time(365*5,0,0,0,0);
    
    // geocentric...
    const orsa::Vector obsPos(0,0,0);

    // temporary variables
    orsa::Vector r, v;
    
    bg->getInterpolatedPosVel(r,v,earth.get(),t);
    
    // no rotation needed for geocentric...
    const orsa::Vector obsPosition = r + equatorialToEcliptic*(obsPos);
    
    bg->getInterpolatedPosVel(r,v,sun.get(),t);
    const orsa::Vector sunPosition = r;
    const orsa::Vector sunVelocity = v;
    
    // bg->getInterpolatedPosition(r,moon.get(),t);
    // const orsa::Vector moonPosition = r;
    
    orsa::Vector u_obs2neo;
    orsa::Double V; // apparent magnitude
    orsa::Vector neo2obs;
    orsa::Vector neo2sun;
    orsa::Double phaseAngle; 
    
    orsa::Double ra,dec;
    
    vecToSky(ra,
	     dec,
	     (sunPosition-obsPosition).normalized());
    const orsa::Double sun_ra  = ra;
    const orsa::Double sun_dec = dec;
    
    // used for both H and V
    const unsigned int fp_mag_min = 15;
    const unsigned int fp_mag_max = 22;
    
    {
      char filename[1024];
      FILE * fp_V[fp_mag_max];
      FILE * fp_H[fp_mag_max];
      
      for (unsigned int k=fp_mag_min; k<=fp_mag_max; ++k) {
	gmp_sprintf(filename,"NEO.PP.real.V%0i.%04.Ffd.dat",
		    k,
		    orsa::Double(orsaSolarSystem::timeToJulian(t)-orsaSolarSystem::timeToJulian(orsaSolarSystem::J2000())).get_mpf_t());
	fp_V[k] = fopen(filename,"w");
   	
	gmp_sprintf(filename,"NEO.PP.real.H%0i.%04.Ffd.dat",
		    k,
		    orsa::Double(orsaSolarSystem::timeToJulian(t)-orsaSolarSystem::timeToJulian(orsaSolarSystem::J2000())).get_mpf_t());
	fp_H[k] = fopen(filename,"w");
      }
      //
      gmp_sprintf(filename,"NEO.PP.real.DETECTED.%04.Ffd.dat",
		  orsa::Double(orsaSolarSystem::timeToJulian(t)-orsaSolarSystem::timeToJulian(orsaSolarSystem::J2000())).get_mpf_t());
      FILE * fp_detected = fopen(filename,"w");
      
      NEOList::const_iterator it = realNEO.begin();
      while (it != realNEO.end()) {
	
	simpleNEOVector(u_obs2neo,
			V,
			neo2obs,
			neo2sun,
			phaseAngle,
			(*it).get(),
			t,
			sunPosition,
			obsPosition,
			orsa::Vector(0,0,1),
			orsa::pi(),
			-1,
			detectionProbabilityThreshold,
			cacheNEOVector);
	
	vecToSky(ra,
		 dec,
		 u_obs2neo);
	
	if ((*it)->detectionProbabilityFromLog(t) > detectionProbabilityThreshold) {
	  gmp_fprintf(fp_detected,
		      "%10i %10i %8.5Ff %+09.5Ff %8.3Ff %8.3Ff %12.9Ff %12.10Ff %12.8Ff"
		      MODEOL,
		      (*it)->randomSeed,
		      (*it)->id,
		      orsa::fmod(24+orsa::radToDeg()*(ra-sun_ra),24).get_mpf_t(),
		      orsa::Double(orsa::radToDeg()*dec).get_mpf_t(),
		      V.get_mpf_t(),
		      (*it)->getH(t,detectionProbabilityThreshold).get_mpf_t(),
		      orsa::FromUnits((*it)->orbit.a,orsa::Unit::AU,-1).get_mpf_t(),
		      (*it)->orbit.e.get_mpf_t(),
		      orsa::Double(orsa::radToDeg()*(*it)->orbit.i).get_mpf_t());
	  fflush(fp_detected);
	} else {
	  for (unsigned int k=fp_mag_min; k<=fp_mag_max; ++k) {
	    if (V <= k) {
	      gmp_fprintf(fp_V[k],
			  "%10i %10i %8.5Ff %+09.5Ff %8.3Ff %8.3Ff %12.9Ff %12.10Ff %12.8Ff"
			  MODEOL,
			  (*it)->randomSeed,
			  (*it)->id,
			  orsa::fmod(24+orsa::radToDeg()*(ra-sun_ra),24).get_mpf_t(),
			  orsa::Double(orsa::radToDeg()*dec).get_mpf_t(),
			  V.get_mpf_t(),
			  (*it)->getH(t,detectionProbabilityThreshold).get_mpf_t(),
			  orsa::FromUnits((*it)->orbit.a,orsa::Unit::AU,-1).get_mpf_t(),
			  (*it)->orbit.e.get_mpf_t(),
			  orsa::Double(orsa::radToDeg()*(*it)->orbit.i).get_mpf_t());
	      fflush(fp_V[k]);
	    }
	    if ((*it)->getH(t,detectionProbabilityThreshold) <= k) {
	      gmp_fprintf(fp_H[k],
			  "%10i %10i %8.5Ff %+09.5Ff %8.3Ff %8.3Ff %12.9Ff %12.10Ff %12.8Ff"
			  MODEOL,
			  (*it)->randomSeed,
			  (*it)->id,
			  orsa::fmod(24+orsa::radToDeg()*(ra-sun_ra),24).get_mpf_t(),
			  orsa::Double(orsa::radToDeg()*dec).get_mpf_t(),
			  V.get_mpf_t(),
			  (*it)->getH(t,detectionProbabilityThreshold).get_mpf_t(),
			  orsa::FromUnits((*it)->orbit.a,orsa::Unit::AU,-1).get_mpf_t(),
			  (*it)->orbit.e.get_mpf_t(),
			  orsa::Double(orsa::radToDeg()*(*it)->orbit.i).get_mpf_t());
	      fflush(fp_H[k]);
	    }
	  }
	}
	
	++it;
      }
      
      for (unsigned int k=fp_mag_min; k<=fp_mag_max; ++k) {
	fclose(fp_V[k]);
 	fclose(fp_H[k]);
      }
      //
      fclose(fp_detected);
      
    }
    
    {
      char filename[1024];
      FILE * fp_V[fp_mag_max];
      FILE * fp_H[fp_mag_max];
      
      for (unsigned int k=fp_mag_min; k<=fp_mag_max; ++k) {
	gmp_sprintf(filename,"NEO.PP.synthetic.V%0i.%04.Ffd.dat",
		    k,
		    orsa::Double(orsaSolarSystem::timeToJulian(t)-orsaSolarSystem::timeToJulian(orsaSolarSystem::J2000())).get_mpf_t());
	fp_V[k] = fopen(filename,"w");
	
	gmp_sprintf(filename,"NEO.PP.synthetic.H%0i.%04.Ffd.dat",
		    k,
		    orsa::Double(orsaSolarSystem::timeToJulian(t)-orsaSolarSystem::timeToJulian(orsaSolarSystem::J2000())).get_mpf_t());
	fp_H[k] = fopen(filename,"w");
      }
      //
      gmp_sprintf(filename,"NEO.PP.synthetic.UNOBSERVED.%04.Ffd.dat",
		  orsa::Double(orsaSolarSystem::timeToJulian(t)-orsaSolarSystem::timeToJulian(orsaSolarSystem::J2000())).get_mpf_t());
      FILE * fp_unobserved = fopen(filename,"w");
      
      NEOList::const_iterator it = syntheticNEO.begin();
      while (it != syntheticNEO.end()) {
	
	simpleNEOVector(u_obs2neo,
			V,
			neo2obs,
			neo2sun,
			phaseAngle,
			(*it).get(),
			t,
			sunPosition,
			obsPosition,
			orsa::Vector(0,0,1),
			orsa::pi(),
			-1,
			detectionProbabilityThreshold,
			cacheNEOVector);
	
	vecToSky(ra,
		 dec,
		 u_obs2neo);
	
	// debug
	/* 
	   {
	   SyntheticNEO * syntheticNEO = dynamic_cast<SyntheticNEO *> ((*it).get());
	   SyntheticNEO::detectionLog::const_iterator log_it = syntheticNEO->getLog().begin();
	   while (log_it != syntheticNEO->getLog().end()) {
	   ORSA_DEBUG("randomSeed: %10i   id: %10i   t: %Ff   H: %Ff",
	   (*it)->randomSeed,
	   (*it)->id,
	   orsaSolarSystem::timeToJulian((*log_it).first).get_mpf_t(),
	   (*it)->getH((*log_it).first,detectionProbabilityThreshold).get_mpf_t());
	   ++log_it;
	   }
	   }
	*/
	
	if ((*it)->getH(t,detectionProbabilityThreshold) < -10) {
	  // never observed synthetic NEOs
	  gmp_fprintf(fp_unobserved,
		      "%10i %10i %8.5Ff %+09.5Ff %8.3Ff %8.3Ff %12.9Ff %12.10Ff %12.8Ff"
		      MODEOL,
		      (*it)->randomSeed,
		      (*it)->id,
		      orsa::fmod(24+orsa::radToDeg()*(ra-sun_ra),24).get_mpf_t(),
		      orsa::Double(orsa::radToDeg()*dec).get_mpf_t(),
		      V.get_mpf_t(),
		      (*it)->getH(t,detectionProbabilityThreshold).get_mpf_t(),
		      orsa::FromUnits((*it)->orbit.a,orsa::Unit::AU,-1).get_mpf_t(),
		      (*it)->orbit.e.get_mpf_t(),
		      orsa::Double(orsa::radToDeg()*(*it)->orbit.i).get_mpf_t());
	  fflush(fp_unobserved);
	} else {
	  for (unsigned int k=fp_mag_min; k<=fp_mag_max; ++k) {
	    if (V <= k) {
	      gmp_fprintf(fp_V[k],
			  "%10i %10i %8.5Ff %+09.5Ff %8.3Ff %8.3Ff %12.9Ff %12.10Ff %12.8Ff"
			  MODEOL,
			  (*it)->randomSeed,
			  (*it)->id,
			  orsa::fmod(24+orsa::radToDeg()*(ra-sun_ra),24).get_mpf_t(),
			  orsa::Double(orsa::radToDeg()*dec).get_mpf_t(),
			  V.get_mpf_t(),
			  (*it)->getH(t,detectionProbabilityThreshold).get_mpf_t(),
			  orsa::FromUnits((*it)->orbit.a,orsa::Unit::AU,-1).get_mpf_t(),
			  (*it)->orbit.e.get_mpf_t(),
			  orsa::Double(orsa::radToDeg()*(*it)->orbit.i).get_mpf_t());
	      fflush(fp_V[k]);
	    }	   
	    if ((*it)->getH(t,detectionProbabilityThreshold) <= k) {
	      gmp_fprintf(fp_H[k],
			  "%10i %10i %8.5Ff %+09.5Ff %8.3Ff %8.3Ff %12.9Ff %12.10Ff %12.8Ff"
			  MODEOL,
			  (*it)->randomSeed,
			  (*it)->id,
			  orsa::fmod(24+orsa::radToDeg()*(ra-sun_ra),24).get_mpf_t(),
			  orsa::Double(orsa::radToDeg()*dec).get_mpf_t(),
			  V.get_mpf_t(),
			  (*it)->getH(t,detectionProbabilityThreshold).get_mpf_t(),
			  orsa::FromUnits((*it)->orbit.a,orsa::Unit::AU,-1).get_mpf_t(),
			  (*it)->orbit.e.get_mpf_t(),
			  orsa::Double(orsa::radToDeg()*(*it)->orbit.i).get_mpf_t());
	      fflush(fp_H[k]);
	    }
	  }
	}
	
	++it;
      }
      
      for (unsigned int k=fp_mag_min; k<=fp_mag_max; ++k) {
	fclose(fp_V[k]);
    	fclose(fp_H[k]);
      }
      //
      fclose(fp_unobserved);
      
    }
    
    exit(0);
  }
  
  while (t <= tStop) {
    
    {
      orsa::Cache<orsa::Time> goodTime = tStop;
      TelescopeList::const_iterator tl_obs = telescopeList.begin();
      /* 
	 if (tl_obs != telescopeList.end()) {
	 // init with first value
	 goodTime = (*tl_obs)->nextObservationTime(t);
	 ++tl_obs;
	 }
      */
      while (tl_obs != telescopeList.end()) {
	if (t <= (*tl_obs)->tStop) {
	  goodTime = std::min((*tl_obs)->nextObservationTime(t), 
			      goodTime.getRef());
	}
	++tl_obs;	    
      }
      if (goodTime.isSet()) {
	// ORSA_DEBUG("delta-time: %Zi [sec]",mpz_class((goodTime.getRef()-t).getMuSec()/mpz_class(1000000)).get_mpz_t());
	t = goodTime.getRef();
      } else {
	ORSA_ERROR("problems...");
      }	
    }
    
    if (t >= tStop) break;
    
    /* 
       ORSA_DEBUG("------------------------------------------------ t: %Ff",
       orsaSolarSystem::timeToJulian(t).get_mpf_t());
    */
    
    /* 
       ORSA_DEBUG("realNEO.size(): %i   syntheticNEO.size(): %i",
       realNEO.size(),
       syntheticNEO.size());
    */
    
    // ORSA_DEBUG("tp.size() = %i",tp.size());
    //
    /* 
       {
       tp.tpList::sort();
       tpList::iterator it = tp.begin();
       while (it != tp.end()) {
       if ((t - (*it).epoch) > recycleTime) {
       it = tp.erase(it);
       } else {
       ++it;
       }
       }
       }
    */
    //
    {
      TelescopeList::const_iterator tl_rcl = telescopeList.begin();
      while (tl_rcl != telescopeList.end()) {
	(*tl_rcl)->applyRecycleTime(t);
	++tl_rcl;	    
      }
    }
    
    TelescopeList::const_iterator tl_main = telescopeList.begin();
    while (tl_main != telescopeList.end()) {
      
      if (!(*tl_main)->isRunning(t)) {
	++tl_main;
	continue;
      }
      
      const orsaSolarSystem::Observatory observatory = (*tl_main)->getObservatory();
      
      // const orsa::Double FOV = (*tl_main)->FOV_DEG * orsa::degToRad();
      //
      // const orsa::Double halfFOV          = 0.5*FOV;
      // const orsa::Double effectiveHalfFOV = 1.15*halfFOV;
      //
      // const orsa::Double cos_FOV              = orsa::cos(FOV);
      // const orsa::Double cos_effectiveHalfFOV = orsa::cos(effectiveHalfFOV);
      // const orsa::Double sin_effectiveHalfFOV = orsa::sin(effectiveHalfFOV);
      
      // const orsa::Double     maxZenithDistanceAngle = parFile->maxZenithDistanceAngle_DEG.getRef() * orsa::degToRad();
      // const orsa::Double       minMoonDistanceAngle = parFile->minMoonDistanceAngle_DEG.getRef()   * orsa::degToRad();
      // const orsa::Double          minMoonPhaseAngle = parFile->minMoonPhaseAngle_DEG.getRef()      * orsa::degToRad();
      //
      // const orsa::Double cos_maxZenithDistanceAngle = orsa::cos(maxZenithDistanceAngle);
      // const orsa::Double sin_maxZenithDistanceAngle = orsa::sin(maxZenithDistanceAngle);
      // const orsa::Double   cos_minMoonDistanceAngle = orsa::cos(  minMoonDistanceAngle);
      
      
      //
      // ORSA_DEBUG("tp.size() = %i",tp.size());
      
      // temporary variables
      orsa::Vector r, v;
      
      /* 
	 const orsa::Vector obsPos(FromUnits(observatory.pxy.getRef(),orsa::Unit::ER)*orsa::cos(observatory.lon.getRef()),
	 FromUnits(observatory.pxy.getRef(),orsa::Unit::ER)*orsa::sin(observatory.lon.getRef()),
	 FromUnits(observatory.pz.getRef() ,orsa::Unit::ER));
      */
      //
      const orsa::Vector obsPos(observatory.pxy.getRef()*orsa::cos(observatory.lon.getRef()),
				observatory.pxy.getRef()*orsa::sin(observatory.lon.getRef()),
				observatory.pz.getRef());
      
      orsa::Matrix earthRotation = orsa::Matrix::identity();
      earthRotation.rotZ(orsaSolarSystem::gmst(t));
      
      bg->getInterpolatedPosVel(r,v,earth.get(),t);
      //
      const orsa::Vector obsPosition = r + equatorialToEcliptic*(earthRotation*obsPos);
      const orsa::Vector obsVelocity = v; // approximated: Earth rotation contribute not kept into account yet
      //
      /**** IMPORTANT: FIX obsVelocity ****/
      
      const orsa::Vector obsNormal = (equatorialToEcliptic*(earthRotation*obsPos.normalized())).normalized();
      
      /* 
	 ORSA_DEBUG("t: %6.3Ff   obsNormal: %+Ff %+Ff %+Ff   gmst: %Ff [deg]",
	 orsaSolarSystem::timeToJulian(t).get_mpf_t(),
	 obsNormal.getX().get_mpf_t(),
	 obsNormal.getY().get_mpf_t(),
	 obsNormal.getZ().get_mpf_t(),
	 orsa::Double((orsaSolarSystem::gmst(t)-orsaSolarSystem::gmst(genEpoch))*orsa::radToDeg()).get_mpf_t());
      */
      
      bg->getInterpolatedPosVel(r,v,sun.get(),t);
      const orsa::Vector sunPosition = r;
      const orsa::Vector sunVelocity = v;
      
      bg->getInterpolatedPosition(r,moon.get(),t);
      const orsa::Vector moonPosition = r;
      
      orsa::Double sunElevation = orsa::halfpi() - orsa::acos(obsNormal*((sunPosition-obsPosition).normalized()));
      
      orsa::Double moonElevation = orsa::halfpi() - orsa::acos(obsNormal*((moonPosition-obsPosition).normalized()));
      
      const orsa::Double moonPhaseAngle = orsa::acos(((obsPosition-moonPosition).normalized())*((sunPosition-moonPosition).normalized()));
      
      /* 
	 if ( (sunElevation*orsa::radToDeg() > (-18.0)) || 
	 (moonPhaseAngle*orsa::radToDeg() < 45.0) ) {
	 t += dt;
	 continue;
	 }	
      */
      
      const orsa::Double minMoonPhaseAngle = (*tl_main)->minMoonPhase_DEG * orsa::degToRad();
      
      if ( (sunElevation*orsa::radToDeg() > (-18.0)) || 
	   (moonPhaseAngle < minMoonPhaseAngle) ) {
	// ORSA_DEBUG("leaving [%s] for angles test",(*tl_main)->name.c_str());
	(*tl_main)->markObservationTime(t);
	++tl_main;
	continue;
      }	
      
      TelescopePointing localTP;
      
      const bool foundGoodTP = (*tl_main)->sampleTP(localTP,
						    t, 
						    telescopeList,
						    syntheticNEO,
						    detectionProbabilityThreshold,
						    sunPosition,
						    moonPosition,
						    obsPosition,
						    obsNormal,
						    cacheNEOVector);
      
      if (foundGoodTP) {
	
       	orsa::Vector u_obs2neo;
	orsa::Double V; // apparent magnitude
	orsa::Vector neo2obs;
	orsa::Vector neo2sun;
	orsa::Double phaseAngle;
	
	typedef std::list<NEOList *> NEOSuperList;
	NEOSuperList NEO_sl;
	//
	NEO_sl.push_back(&realNEO);
	NEO_sl.push_back(&syntheticNEO);
	
	NEOSuperList::const_iterator NEO_sl_it = NEO_sl.begin();
	while (NEO_sl_it != NEO_sl.end()) {
	  
	  NEOList & localNEOList = *(*NEO_sl_it);
	  
	  NEOList::const_iterator it = localNEOList.begin();
	  while (it != localNEOList.end()) {
	    
	    simpleNEOVector(u_obs2neo,
			    V,
			    neo2obs,
			    neo2sun,
			    phaseAngle,
			    (*it).get(),
			    t,
			    sunPosition,
			    obsPosition,
			    localTP.u,
			    (*tl_main)->effectiveHalfFOV,
			    (*tl_main)->cos_effectiveHalfFOV,
			    detectionProbabilityThreshold,
			    cacheNEOVector);
	    
	    /* 
	       ORSA_DEBUG("NEO id: %i   scalar product: %Ff",
	       (*it)->id,
	       orsa::Double(u_obs2neo*localTP.u).get_mpf_t());
	    */
	    
	    if ((u_obs2neo*localTP.u) > (*tl_main)->cos_effectiveHalfFOV) {  
	      
	      RealNEO      *      realNEO = dynamic_cast<     RealNEO *> ((*it).get());
	      SyntheticNEO * syntheticNEO = dynamic_cast<SyntheticNEO *> ((*it).get());
	      
	      if (realNEO) {
		
		if (boinc_is_standalone()) {
		  ORSA_DEBUG("real NEO in FOV: id: %10i   randomSeed: %10i   V: %5.2Ff   H: %5.2Ff   t: %.5Ff [%s]",
			     (*it)->id,
			     (*it)->randomSeed,
			     V.get_mpf_t(),
			     (*it)->getH(t,detectionProbabilityThreshold).get_mpf_t(),
			     orsaSolarSystem::timeToJulian(t).get_mpf_t(),
			     (*tl_main)->name.c_str());
		}
		
		const orsa::Double prob = (*tl_main)->detectionProbability(V);
		//
		if ((prob < orsa::zero()) || (prob > orsa::one())) {
		  ORSA_DEBUG("...probability problems...");
		}
		
		RealNEO::LogEntry * re = new RealNEO::LogEntry;
		//
		re->p = prob;
		//
		re->telescopeName = (*tl_main)->name;
		//
		(*it)->logInsert(t,re,writeNEOLogFile);
		
		if (prob > orsa::zero()) {
		  if (boinc_is_standalone()) {
		    ORSA_DEBUG("real NEO potentially observed: id: %10i   randomSeed: %10i   V: %5.2Ff   H: %5.2Ff   p: %4.2Ff   t: %.5Ff",
			       (*it)->id,
			       (*it)->randomSeed,
			       V.get_mpf_t(),
			       (*it)->getH(t,detectionProbabilityThreshold).get_mpf_t(),
			       prob.get_mpf_t(),
			       orsaSolarSystem::timeToJulian(t).get_mpf_t());
		  }
		}
		
		/* 
		   if ((*it)->detectionProbabilityFromLog(t) > detectionProbabilityThreshold) {
		   it = localNEOList.erase(it);
		   continue;
		   }
		*/
		
	      } else if (syntheticNEO) {
		
		if (boinc_is_standalone()) {
		  ORSA_DEBUG("synthetic NEO in FOV: id: %10i   randomSeed: %10i   neo2obs: %5.3Ff AU   neo2sun: %5.3Ff AU   phaseAngle: %7.3Ff DEG    t: %.5Ff [%s]",
			     (*it)->id,
			     (*it)->randomSeed,
			     orsa::FromUnits(neo2obs.length(),orsa::Unit::AU,-1).get_mpf_t(),
			     orsa::FromUnits(neo2sun.length(),orsa::Unit::AU,-1).get_mpf_t(),
			     orsa::Double(orsa::radToDeg()*phaseAngle).get_mpf_t(),
			     orsaSolarSystem::timeToJulian(t).get_mpf_t(),
			     (*tl_main)->name.c_str());
		}
		
		SyntheticNEO::LogEntry * se = new SyntheticNEO::LogEntry;
		//
		se->neo2obs    = neo2obs.length();
		se->neo2sun    = neo2sun.length();
		se->phaseAngle = phaseAngle;
		//
		se->limitingMagnitude = (*tl_main)->limitingMagnitude;
		//
		se->telescopeName = (*tl_main)->name;
		//
		syntheticNEO->logInsert(t,se,writeNEOLogFile);
		
	      } else {
		ORSA_DEBUG("case not handled");
	      }
	      
	    }
	    
	    ++it;
	  }
	  
	  /* 
	     if (num_in_FOV > 0) {
	     ORSA_DEBUG("FOV stats -- NEOs in FOV: %i   V range: %5.2Ff to %5.2Ff",
	     num_in_FOV,
	     V_min_in_FOV.get_mpf_t(),
	     V_max_in_FOV.get_mpf_t());
	     }
	  */
	  
	  ++NEO_sl_it;
	}
	
      }
      
      ++tl_main;
    }
    
    /* 
       ORSA_DEBUG("t: %6.3Ff   sunElevation: %6.2Ff [deg]   moonElevation: %6.2Ff [deg]   moonPhase: %6.2Ff [deg]",
       orsaSolarSystem::timeToJulian(t).get_mpf_t(),
       orsa::Double(sunElevation*orsa::radToDeg()).get_mpf_t(),
       orsa::Double(moonElevation*orsa::radToDeg()).get_mpf_t(),
       orsa::Double(moonPhaseAngle*orsa::radToDeg()).get_mpf_t());
    */
    
    {
      const orsa::Double fractionCompleted = (t-genEpoch).asDouble() / (tStop-genEpoch).asDouble();
      boinc_fraction_done(fractionCompleted.get_d());
    }
    
    // write checkpoint
    if (boinc_is_standalone()) {
      const  unsigned int ckpt_freq  = 256;
      static unsigned int ckpt_count = 0;
      // ORSA_DEBUG("ckpt_count: %i", ckpt_count);
      if ((++ckpt_count % ckpt_freq) == 0) {
	boinc_begin_critical_section();
	writeCheckpoint(t,
			telescopeList);
	boinc_end_critical_section();
      }
    } else if (boinc_time_to_checkpoint()) {
      writeCheckpoint(t,
		      telescopeList);
      boinc_checkpoint_completed();
    }
    
  }
  
  if (checkpointUsed) {
    ORSA_DEBUG(" -- this run was started from a checkpoint -- ");
  }
  
  boinc_resolve_filename_s("NEO.real.dat",resolvedFileName);
  int fd_real = open(resolvedFileName.c_str(), open_flag, open_mode);
  
  boinc_resolve_filename_s("NEO.synthetic.dat",resolvedFileName);
  int fd_synthetic = open(resolvedFileName.c_str(), open_flag, open_mode);
  
  {
    realNEO.sort(NEOListSortPredicate);
    char line[1024];
    NEOList::const_iterator neo_it = realNEO.begin();
    while (neo_it != realNEO.end()) {
      RealNEO * realNEO = dynamic_cast<RealNEO *> ((*neo_it).get());
      if (realNEO) {
	const RealNEO::detectionLog & log = realNEO->getLog();
	RealNEO::detectionLog::const_iterator log_it = log.begin();
	while (log_it != log.end()) {
	  gmp_snprintf(line,1024,"%Zi %10i %10i %11.9Ff %s"
		       MODEOL,
		       (*log_it).first.getMuSec().get_mpz_t(),
		       realNEO->randomSeed,
		       realNEO->id,
		       (*log_it).second->p.getRef().get_mpf_t(),
		       (*log_it).second->telescopeName.getRef().c_str());
	  write(fd_real,line,strlen(line));
	  ++log_it;
	}
      }
      ++neo_it;
    }
  }
  
  {
    syntheticNEO.sort(NEOListSortPredicate);
    char line[1024];
    NEOList::const_iterator neo_it = syntheticNEO.begin();
    while (neo_it != syntheticNEO.end()) {
      SyntheticNEO * syntheticNEO = dynamic_cast<SyntheticNEO *> ((*neo_it).get());
      if (syntheticNEO) {
	const SyntheticNEO::detectionLog & log = syntheticNEO->getLog();
	SyntheticNEO::detectionLog::const_iterator log_it = log.begin();
	while (log_it != log.end()) {
	  gmp_snprintf(line,1024,"%Zi %10i %10i %12.9Ff %12.9Ff %13.9Ff %6.3Ff %s"
		       MODEOL,
		       (*log_it).first.getMuSec().get_mpz_t(),
		       syntheticNEO->randomSeed,
		       syntheticNEO->id,
		       FromUnits((*log_it).second->neo2obs.getRef(),orsa::Unit::AU,-1).get_mpf_t(),
		       FromUnits((*log_it).second->neo2sun.getRef(),orsa::Unit::AU,-1).get_mpf_t(),
		       orsa::Double(orsa::radToDeg()*(*log_it).second->phaseAngle.getRef()).get_mpf_t(),
		       (*log_it).second->limitingMagnitude.getRef().get_mpf_t(),
		       (*log_it).second->telescopeName.getRef().c_str());
	  write(fd_synthetic,line,strlen(line));
	  ++log_it;
	}
      }
      ++neo_it;
    }
  }
  
  close(fd_real);
  close(fd_synthetic);
  
  boinc_finish(0);
  
  return 0;
}
