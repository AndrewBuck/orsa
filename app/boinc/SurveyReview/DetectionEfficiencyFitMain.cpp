#include "skycoverage.h"
#include "eta.h"

int main(int argc, char ** argv) {
  
  orsa::Debug::instance()->initTimer();
  
  if (argc != 2) {
    printf("Usage: %s <efficiency-data-file>\n",argv[0]);
    exit(0);
  }
  
  EfficiencyMultifit::DataStorage data;
  {
    FILE * fp = fopen(argv[1],"r");
    if (!fp) {
      ORSA_DEBUG("cannot open file [%s]",argv[1]);
      exit(0);
    }
    char line[1024];
    double V, apparentVelocity, eta, sigmaEta;
    while (fgets(line,1024,fp)) {
      sscanf(line,"%lf %*s %lf %*s %lf %lf",
	     &V,
	     &apparentVelocity,
	     &eta,
	     &sigmaEta);
      // apparentVelocity is in arcsec/hour
      apparentVelocity = orsa::FromUnits(apparentVelocity*orsa::arcsecToRad(),orsa::Unit::HOUR,-1);
      //
      ORSA_DEBUG("reading: %g %g %g %g",
		 V,apparentVelocity,eta,sigmaEta);
      //
      EfficiencyMultifit::DataElement el;
      el.V=V;
      el.U=apparentVelocity;
      el.eta=eta;
      el.sigmaEta=sigmaEta;
      data.push_back(el);
    }
    fclose(fp);
  }
  
  const double V0=17.0;
  
  const double U0=orsa::FromUnits(100.0*orsa::arcsecToRad(),orsa::Unit::HOUR); // arcsec/hour
  
  osg::ref_ptr<EfficiencyMultifit> etaFit= new EfficiencyMultifit;
  const bool success = etaFit->fit(data,V0,U0);
  if (!success) { 
    ORSA_DEBUG("problems??");
  }
  osg::ref_ptr<const orsa::MultifitParameters> parFinal = etaFit->getMultifitParameters();
  // save final parameters
  
  exit(0);
}
