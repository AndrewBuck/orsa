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
  const double  eta0_V = parFinal->get("eta0_V");
  const double     c_V = parFinal->get("c_V");
  const double V_limit = parFinal->get("V_limit");
  const double     w_V = parFinal->get("w_V");
  const double U_limit = parFinal->get("U_limit");
  const double     w_U = parFinal->get("w_U");
  
  {
    // output for testing
    FILE * fpv = fopen("v.fit.dat","w");
    FILE * fpu = fopen("u.fit.dat","w");
    for (unsigned int k=0; k<data.size(); ++k) {
      const double eta_V = SkyCoverage::eta_V(data[k].V.getRef(),
					      V_limit,
					      eta0_V,
					      c_V,
					      V0,
					      w_V);
      const double eta_U = SkyCoverage::eta_U(data[k].U.getRef(),
					      U_limit,
					      1.0, // eta0_U,
					      0.0, // c_U,
					      U0,
					      w_U);
      if (eta_V<1e-2) continue;
      if (eta_U<1e-2) continue;
      fprintf(fpv,"%g %g %g %g\n",
	      data[k].V.getRef(),
	      eta_V,
	      data[k].eta.getRef()/eta_U,
	      data[k].sigmaEta.getRef()/eta_U);
      fprintf(fpu,"%g %g %g %g\n",
	      orsa::FromUnits(data[k].U.getRef()*orsa::radToArcsec(),orsa::Unit::HOUR),
	      eta_U,
	      data[k].eta.getRef()/eta_V,
	      data[k].sigmaEta.getRef()/eta_V);
    }
    fclose(fpv);
    fclose(fpu);
  }
  
  
  exit(0);
}
