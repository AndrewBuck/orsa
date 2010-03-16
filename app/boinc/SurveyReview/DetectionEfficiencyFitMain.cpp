#include <orsa/statistic.h>
#include "skycoverage.h"
#include "eta.h"

class EfficiencyStatistics : public orsa::WeightedStatistic<double> {
public:
  EfficiencyStatistics(const double & aux) : 
    orsa::WeightedStatistic<double>(),
    center(aux) { }
public:
  const double center; // V or U
  orsa::Cache<double> fit; // eta_V or eta_U
};

// for sorting
class EfficiencyStatistics_ptr_cmp {
public:
  bool operator() (EfficiencyStatistics * lhs,
		   EfficiencyStatistics * rhs) const {
    return (lhs->center < rhs->center);
  }
};

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
#warning keep vars in sync with fitting code
  
  {
    // output for testing
    typedef std::list< osg::ref_ptr<EfficiencyStatistics> > AllStat;
    AllStat stat_V, stat_U;  
    //
    for (unsigned int k=0; k<data.size(); ++k) {
      
      osg::ref_ptr<EfficiencyStatistics> sV;
      {
	AllStat::const_iterator it = stat_V.begin();
	while (it != stat_V.end()) {
	  if ((*it).get() != 0) {
	    if ((*it)->center == data[k].V.getRef()) {
	      sV = (*it).get();
	    }   
	  }
	  ++it;
	}
	if (sV.get() == 0) {
	  sV = new EfficiencyStatistics(data[k].V.getRef());
	  stat_V.push_back(sV);
	}
      }
      
      osg::ref_ptr<EfficiencyStatistics> sU;
      {
	AllStat::const_iterator it = stat_U.begin();
	while (it != stat_U.end()) {
	  if ((*it).get() != 0) {
	    if ((*it)->center == data[k].U.getRef()) {
	      sU = (*it).get();
	    }   
	  }
	  ++it;
	}
	if (sU.get() == 0) {
	  sU = new EfficiencyStatistics(data[k].U.getRef());
	  stat_U.push_back(sU);
	}
      }
      
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
      
      if (data[k].sigmaEta.getRef() > 0) {
	if (eta_U != 0) {
	  sV->insert(data[k].eta.getRef()/eta_U,
		     orsa::square(eta_U/data[k].sigmaEta.getRef()));
	}
	if (eta_V != 0) {
	  sU->insert(data[k].eta.getRef()/eta_V,
		     orsa::square(eta_V/data[k].sigmaEta.getRef()));
	}
      }
      
      if (!sV->fit.isSet()) sV->fit = eta_V;
      if (!sU->fit.isSet()) sU->fit = eta_U;
    }
    
    // sort lists
    stat_V.sort(EfficiencyStatistics_ptr_cmp());
    stat_U.sort(EfficiencyStatistics_ptr_cmp());
    
    {
      FILE * fpv = fopen("v.fit.dat","w");
      AllStat::const_iterator it = stat_V.begin();
      while (it != stat_V.end()) {
	if (!(*it)->fit.isSet()) {
	  ORSA_DEBUG("problems: fit value not set");
	  ++it;
	  continue;
	}
	fprintf(fpv,"%g %g %g %g\n",
		(*it)->center,
		(*it)->fit.getRef(),
		(*it)->average(),
		(*it)->averageError());
	++it;
      }
      fclose(fpv);
    }
    
    {
      FILE * fpu = fopen("u.fit.dat","w");
      AllStat::const_iterator it = stat_U.begin();
      while (it != stat_U.end()) {
	if (!(*it)->fit.isSet()) {
	  ORSA_DEBUG("problems: fit value not set");
	  ++it;
	  continue;
	}
	fprintf(fpu,"%g %g %g %g\n",
		orsa::FromUnits((*it)->center*orsa::radToArcsec(),orsa::Unit::HOUR), // arcsec/hour
		(*it)->fit.getRef(),
		(*it)->average(),
		(*it)->averageError());
	++it;
      }
      fclose(fpu);
    }
    
  }
  
  exit(0);
}
