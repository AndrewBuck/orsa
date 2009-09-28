#include "DAWN.h"

int main(int argc, char ** argv) {
  
  if ((argc != 4) && (argc !=5)) {
    ORSA_DEBUG("Usage: %s <orbitRadius_km> <SCENARIO> <duration_days> [thrust_mN]",argv[0]);
    exit(0);
  }
  
  // QApplication app(argc, argv);
  
  // orsaQt::Debug::instance()->initTimer();
  //
  orsa::Debug::instance()->initTimer();
  
  ORSA_DEBUG("process ID: %i",getpid());
  
  const double orbitRadius = orsa::FromUnits(atof(argv[1]),orsa::Unit::KM);
  
  const SCENARIO scenario = s2S(argv[2]);
  
  const orsa::Time duration(atoi(argv[3]),0,0,0,0);
  
  double thrust_mN=0.0;
  if (argc==5) {
    thrust_mN = atof(argv[4]);
  }
  
  osg::ref_ptr<orsa::BodyGroup> bg =
    run(orbitRadius,
	scenario,
	duration,
	thrust_mN);
  
  if (!bg.get()) {
    exit(0);
  }	
  
  return 0;
}
