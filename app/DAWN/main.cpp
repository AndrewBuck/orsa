#include "DAWN.h"

int main(int argc, char ** argv) {
  
    if ((argc != 4) && (argc !=5) && (argc !=6) && (argc !=7)) {
        ORSA_DEBUG("Usage: %s <orbitRadius_km> <SCENARIO> <duration_days> [degree] [phase_DEG] [thrust_mN]",argv[0]);
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
  
    unsigned int degree=8;
    if (argc>=5) {
        degree = atoi(argv[4]);
    }
    ORSA_DEBUG("degree: %i",degree);
  
    double phase_DEG=0.0;
    if (argc>=6) {
        phase_DEG = atof(argv[5]);
        ORSA_DEBUG("phase: %g [deg]",phase_DEG);
    }
  
    double thrust_mN=0.0;
    if (argc>=7) {
        thrust_mN = atof(argv[6]);
        ORSA_DEBUG("thrust: %g mN",thrust_mN);
    }
  
    osg::ref_ptr<orsa::BodyGroup> bg =
        run(orbitRadius,
            scenario,
            duration,
            degree,
            phase_DEG,
            thrust_mN);
  
    if (!bg.get()) {
        exit(0);
    }	
  
    return 0;
}
