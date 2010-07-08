#include <orsa/util.h>
#include <orsaInputOutput/MPC_observations.h>
#include <orsaSolarSystem/datetime.h>

class CustomMPCObservationsFile : public orsaInputOutput::MPCObservationsFile {
public:
  CustomMPCObservationsFile() :
    orsaInputOutput::MPCObservationsFile() { }
public:
  // stripped-down code
  bool processLine(const char * line) {
    
    std::string s_epoch;
    
    s_epoch.assign(line,15,17);
    orsa::Time epoch;
    {
      int y, m; 
      double d;
      gmp_sscanf(s_epoch.c_str(),"%d %d %lf",&y,&m,&d);
      epoch = orsaSolarSystem::FromTimeScale(orsaSolarSystem::gregorTime(y,m,d),
					     orsaSolarSystem::TS_UTC);
      // orsa::print(epoch);
      if (select_startEpoch.isSet()) {
	if (epoch < select_startEpoch.getRef()) {
	  // ORSA_DEBUG("--OUT--");
	  return false;
	}
      }
      if (select_stopEpoch.isSet()) {
	if (epoch > select_stopEpoch.getRef()) {
	  // ORSA_DEBUG("--OUT--");
	  return false;
	}
      }
    }
    
    // output good line
    printf("%s\n",line);
    
    return true;
  }
public:
  bool processLines(const char *,
		    const char *) {
    return false;
  }
};

int main(int argc, char ** argv) {
  
  orsa::Debug::instance()->initTimer();
  
  if (argc != 2) {
    printf("Usage: %s <year>\n",argv[0]);
    exit(0);
  }
  
  const int year = atoi(argv[1]);
  
  ORSA_DEBUG("process ID: %i",getpid());
  
  osg::ref_ptr<CustomMPCObservationsFile> obsFile =
    new CustomMPCObservationsFile;
  // a couple days before and after
  obsFile->select_startEpoch = orsaSolarSystem::gregorTime(year-1,12,29);
  obsFile->select_stopEpoch  = orsaSolarSystem::gregorTime(year+1, 1, 3);
  ORSA_DEBUG("select start/stop:");
  orsa::print(obsFile->select_startEpoch.getRef());
  orsa::print(obsFile->select_stopEpoch.getRef());
  obsFile->setFileName("NumObs.txt.gz");
  obsFile->read();
  obsFile->setFileName("UnnObs.txt.gz");
  obsFile->read();
  
  return 0;
}
