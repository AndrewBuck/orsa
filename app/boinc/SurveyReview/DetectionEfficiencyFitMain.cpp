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
    printf("Usage: %s <allEta-file>\n",argv[0]);
    exit(0);
  }
  
  ORSA_DEBUG("process ID: %i",getpid());
  
  const std::string basename = SkyCoverage::basename(argv[1]);
  
  std::vector<EfficiencyData> etaData;
  {
    FILE * fp = fopen(argv[1],"r");
    if (!fp) {
      ORSA_DEBUG("cannot open file [%s]",argv[1]);
      exit(0);
    }
    char line[1024];
    double H, V, apparentVelocity;
    double solarElongation, lunarElongation;
    double lunarPhase;
    double minAirMass;
    char id[1024];
    int observed;
    int discovered;
    while (fgets(line,1024,fp)) {
      sscanf(line,"%lf %s %lf %lf %lf %lf %lf %lf %i %i",
	     &H,
	     id,
	     &V,
	     &apparentVelocity,
	     &solarElongation,
	     &lunarElongation,
	     &lunarPhase,
	     &minAirMass,
	     &observed,
	     &discovered);
      // convert
      apparentVelocity = orsa::FromUnits(apparentVelocity*orsa::arcsecToRad(),orsa::Unit::HOUR,-1);
      solarElongation = orsa::degToRad()*solarElongation;
      lunarElongation = orsa::degToRad()*lunarElongation;
      lunarPhase      = orsa::degToRad()*lunarPhase;
      //
      EfficiencyData ed;
      ed.H=H;
      if (atoi(id) != 0) {
	ed.number = atoi(id);
      } else {
	ed.designation = id;
      }
      ed.V = V;
      ed.apparentVelocity = apparentVelocity;
      ed.solarElongation  = solarElongation;
      ed.lunarElongation  = lunarElongation;
      ed.lunarPhase       = lunarPhase;
      ed.minAirMass       = minAirMass;
      ed.observed   = (observed==1);
      ed.discovered = (discovered==1);
      //
      etaData.push_back(ed);
    }
    fclose(fp);
  }
  
  const double V0=16.0;
  EfficiencyMultifit::DataStorage data;
  const bool write_fp_eta=false;
  {
    FILE * fp_eta;
    if (write_fp_eta) {
      char filename[1024];
      sprintf(filename,"%s.eta.dat",basename.c_str());
      fp_eta = fopen(filename,"w"); 
    }
    
    double V=V0;
    const double dV=0.2;
    while (V<=24.0) {
      
      double apparentVelocity=orsa::FromUnits(0.1*orsa::arcsecToRad(),orsa::Unit::HOUR,-1);
      const double apparentVelocityFactor=1.2;
      while (apparentVelocity<orsa::FromUnits(300.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1)) {
	
	double solarElongation=orsa::pi();
	const double deltaSolarElongation=-10.0*orsa::degToRad(); // negative, as we go from 180 deg to 0
	while (solarElongation>=0) {
	  
	  double lunarElongation=orsa::pi();
	  const double deltaLunarElongation=-10.0*orsa::degToRad(); // negative, as we go from 180 deg to 0
	  while (lunarElongation>=0) {
	    
	    // airmass
	    double AM=1.0;
	    const double AMFactor=1.1;
	    while (AM<10.0) {
	      
	      unsigned int Nobs=0,Ndsc=0,Ntot=0; // dsc=discovered
	      for (unsigned int k=0; k<etaData.size(); ++k) {
		if ((etaData[k].V.getRef()>=V) && 
		    (etaData[k].V.getRef()< V+dV) && 
		    (etaData[k].apparentVelocity.getRef()>=apparentVelocity) &&
		    (etaData[k].apparentVelocity.getRef()< apparentVelocityFactor*apparentVelocity) &&
		    (etaData[k].solarElongation.getRef()<=solarElongation) && 
		    (etaData[k].solarElongation.getRef()> solarElongation+deltaSolarElongation) &&
		    (etaData[k].lunarElongation.getRef()<=lunarElongation) && 
		    (etaData[k].lunarElongation.getRef()> lunarElongation+deltaLunarElongation) && 
		    (etaData[k].minAirMass.getRef()>=AM) &&
		    (etaData[k].minAirMass.getRef()< AMFactor*AM)) {
		  ++Ntot;
		  if (etaData[k].observed.getRef()) {
		    ++Nobs;
		  }
		  if (etaData[k].discovered.getRef()) {
		    ++Ndsc;
		  }	
		}
	      }
	      
	      // write point only if Ntot != 0
	      if (Ntot>=2) {
		const double      eta = (double)Nobs/(double)Ntot;
		const double sigmaEta = (double)(sqrt(Nobs+1))/(double)(Ntot); // Poisson counting statistics; using Nobs+1 instead of Nobs to have positive sigma even when Nobs=0
		
		if (write_fp_eta) {
		  fprintf(fp_eta,
			  "%5.2f %4.2f %6.2f %4.2f %6.2f %6.2f %6.2f %6.2f %6.3f %6.3f %5.3f %5.3f %5i %5i %5i\n",
			  V+0.5*dV, 
			  dV, 
			  orsa::FromUnits(apparentVelocity*0.5*(1.0+apparentVelocityFactor)*orsa::radToArcsec(),orsa::Unit::HOUR),
			  apparentVelocityFactor,
			  solarElongation+0.5*deltaSolarElongation,
			  deltaSolarElongation,
			  lunarElongation+0.5*deltaLunarElongation,
			  deltaLunarElongation,
			  AM*0.5*(1.0+AMFactor),
			  AMFactor,
			  eta,
			  sigmaEta,
			  Nobs,
			  Ndsc,
			  Ntot);
		}
		
		EfficiencyMultifit::DataElement el;
		//
		el.V=V+0.5*dV;
		el.U=apparentVelocity*0.5*(1.0+apparentVelocityFactor);		 
		el.SE=solarElongation+0.5*deltaSolarElongation;	 
		el.LE=lunarElongation+0.5*deltaLunarElongation;
		// add lunar phase here
		el.AM=AM;
		el.eta=eta;
		el.sigmaEta=sigmaEta;
		el.Nobs=Nobs;
		el.Ndsc=Ndsc;
		el.Ntot=Ntot;
		//
		data.push_back(el);
	      }
	      
	      AM *= AMFactor;
	    }
	    
	    lunarElongation += deltaLunarElongation;
	  }
	  
	  solarElongation += deltaSolarElongation;
	}
	
	apparentVelocity *= apparentVelocityFactor;
      }
      
      V += dV;
    }
    
    if (write_fp_eta) {
      fclose(fp_eta);
    }
    
  }    
  
  osg::ref_ptr<orsaInputOutput::MPCObsCodeFile> obsCodeFile = new orsaInputOutput::MPCObsCodeFile;
  obsCodeFile->setFileName("obscode.dat");
  obsCodeFile->read();
  
  char fitFilename[1024];
  sprintf(fitFilename,"%s.fit.dat",basename.c_str());
  
  unsigned int non_zero_n   = 0;
  for (unsigned int k=0; k<data.size(); ++k) {
    if (data[k].eta.getRef()>0.0) {
      ++non_zero_n;
    }	
  }
  ORSA_DEBUG("non_zero_n: %i",non_zero_n);
  
  // multifit parameters  
  osg::ref_ptr<orsa::MultifitParameters> par = new orsa::MultifitParameters;
  //
  const double initial_U_limit = orsa::FromUnits( 40.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1); // arcsec/hour
  const double initial_w_U     = orsa::FromUnits( 10.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1); // arcsec/hour
  //    
  // V pars
  par->insert("eta0_V",  0.90, 0.0001);
  par->insert("c_V",     0.01, 0.00001);
  par->insert("V_limit",19.00, 0.001);
  par->insert("w_V",     0.50, 0.0001); 
  // U pars
  par->insert("U_limit",    initial_U_limit, 0.01*initial_U_limit);
  par->insert("w_U",        initial_w_U,     0.01*initial_w_U); 
  // mixing 
  par->insert("beta",  0.0, 0.00001);
  // ranges
  par->setRange("eta0_V",0.0,1.0);
  // fixed?
  // par->setFixed("U_limit",true);
  // par->setFixed("w_U",true);
  // par->setFixed("beta",true);
  
  if (non_zero_n < par->sizeNotFixed()) {
    ORSA_DEBUG("not enought data for fit, exiting");
    // just "touch" the output file
    FILE * fp = fopen(fitFilename,"w");
    fclose(fp);
    exit(0);
  }
  
  osg::ref_ptr<EfficiencyMultifit> etaFit= new EfficiencyMultifit;
  const bool success = etaFit->fit(par.get(),
				   data,
				   V0,
				   fitFilename,
				   basename,
				   obsCodeFile.get());
  if (!success) { 
    ORSA_DEBUG("fitting did not converge, exiting");
    // just "touch" the output file
    FILE * fp = fopen(fitFilename,"w");
    fclose(fp);
    exit(0);
  }
  osg::ref_ptr<const orsa::MultifitParameters> parFinal = etaFit->getMultifitParameters();
  // save final parameters
  const double  eta0_V = parFinal->get("eta0_V");
  const double     c_V = parFinal->get("c_V");
  const double V_limit = parFinal->get("V_limit");
  const double     w_V = parFinal->get("w_V");
  //
  const double U_limit = parFinal->get("U_limit");
  const double     w_U = parFinal->get("w_U");
  //
  const double    beta = parFinal->get("beta");
  
  if (1) {
    // output for testing
    typedef std::list< osg::ref_ptr<EfficiencyStatistics> > AllStat;
    AllStat stat_V, stat_U; 
    AllStat stat_SE; // solarElongation
    AllStat stat_LE; // lunarElongation
    // AllStat stat_LP; // lunarPhase
    AllStat stat_AM; // minAirMass
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
      
      // solarElongation
      osg::ref_ptr<EfficiencyStatistics> sSE;
      {
	AllStat::const_iterator it = stat_SE.begin();
	while (it != stat_SE.end()) {
	  if ((*it).get() != 0) {
	    if ((*it)->center == data[k].SE.getRef()) {
	      sSE = (*it).get();
	    }   
	  }
	  ++it;
	}
	if (sSE.get() == 0) {
	  sSE = new EfficiencyStatistics(data[k].SE.getRef());
	  stat_SE.push_back(sSE);
	}
      }
      
      // lunarElongation
      osg::ref_ptr<EfficiencyStatistics> sLE;
      {
	AllStat::const_iterator it = stat_LE.begin();
	while (it != stat_LE.end()) {
	  if ((*it).get() != 0) {
	    if ((*it)->center == data[k].LE.getRef()) {
	      sLE = (*it).get();
	    }   
	  }
	  ++it;
	}
	if (sLE.get() == 0) {
	  sLE = new EfficiencyStatistics(data[k].LE.getRef());
	  stat_LE.push_back(sLE);
	}
      }
      
      // minAirMass
      osg::ref_ptr<EfficiencyStatistics> sAM;
      {
	AllStat::const_iterator it = stat_AM.begin();
	while (it != stat_AM.end()) {
	  if ((*it).get() != 0) {
	    if ((*it)->center == data[k].AM.getRef()) {
	      sAM = (*it).get();
	    }   
	  }
	  ++it;
	}
	if (sAM.get() == 0) {
	  sAM = new EfficiencyStatistics(data[k].AM.getRef());
	  stat_AM.push_back(sAM);
	}
      }
      
      const double eta = SkyCoverage::eta(data[k].V.getRef(),
					  V_limit,
					  eta0_V,
					  V0,
					  c_V,
					  w_V,
					  data[k].U.getRef(),
					  U_limit,
					  w_U,
					  beta);
      
      const double nominal_eta_V = SkyCoverage::nominal_eta_V(data[k].V.getRef(),
							      V_limit,
							      eta0_V,
							      V0,
							      c_V,
							      w_V);
      
      const double nominal_eta_U = SkyCoverage::nominal_eta_U(data[k].U.getRef(),
							      U_limit,
							      w_U);
      
      // at this point, small sigmas can become very large because divided by a very small eta's
      if (data[k].sigmaEta.getRef() > 0) {
	if (eta != 0) {
	  if (nominal_eta_V != 0) {
	    sV->insert(data[k].eta.getRef()*nominal_eta_V/eta,
		       orsa::square(eta/(nominal_eta_V*data[k].sigmaEta.getRef())));
	  }
	  if (nominal_eta_U != 0) {
	    sU->insert(data[k].eta.getRef()*nominal_eta_U/eta,
		       orsa::square(eta/(nominal_eta_U*data[k].sigmaEta.getRef())));
	  }
	  sSE->insert(data[k].eta.getRef()/eta,
		      orsa::square(eta/(data[k].sigmaEta.getRef())));
	  sLE->insert(data[k].eta.getRef()/eta,
		      orsa::square(eta/(data[k].sigmaEta.getRef())));
	  sAM->insert(data[k].eta.getRef()/eta,
		      orsa::square(eta/(data[k].sigmaEta.getRef())));
	}
	if ( !sV->fit.isSet())  sV->fit = nominal_eta_V;
	if ( !sU->fit.isSet())  sU->fit = nominal_eta_U;
	if (!sSE->fit.isSet()) sSE->fit = 1.0;
 	if (!sLE->fit.isSet()) sLE->fit = 1.0;
   	if (!sAM->fit.isSet()) sAM->fit = 1.0;
      }
    }
    
    // sort lists
    stat_V.sort(EfficiencyStatistics_ptr_cmp());
    stat_U.sort(EfficiencyStatistics_ptr_cmp());
    stat_SE.sort(EfficiencyStatistics_ptr_cmp());
    stat_LE.sort(EfficiencyStatistics_ptr_cmp());
    stat_AM.sort(EfficiencyStatistics_ptr_cmp());
    
    {
      FILE * fpv = fopen("v.fit.dat","w");
      AllStat::const_iterator it = stat_V.begin();
      while (it != stat_V.end()) {
	if (!(*it)->fit.isSet()) {
	  ORSA_DEBUG("problems: fit value not set");
	  ++it;
	  continue;
	}
	if ((*it)->averageError()==0.0) {
	  // ORSA_DEBUG("--SKIPPING--");
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
	if ((*it)->averageError()==0.0) {
	  // ORSA_DEBUG("--SKIPPING--");
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
    
    {
      FILE * fpse = fopen("se.fit.dat","w");
      AllStat::const_iterator it = stat_SE.begin();
      while (it != stat_SE.end()) {
	if (!(*it)->fit.isSet()) {
	  ORSA_DEBUG("problems: fit value not set");
	  ++it;
	  continue;
	}
	if ((*it)->averageError()==0.0) {
	  // ORSA_DEBUG("--SKIPPING--");
	  ++it;
	  continue;
	}
      	fprintf(fpse,"%g %g %g %g\n",
		orsa::radToDeg()*(*it)->center,
		(*it)->fit.getRef(),
		(*it)->average(),
		(*it)->averageError());
	++it;
      }
      fclose(fpse);
    }
    
    {
      FILE * fple = fopen("le.fit.dat","w");
      AllStat::const_iterator it = stat_LE.begin();
      while (it != stat_LE.end()) {
	if (!(*it)->fit.isSet()) {
	  ORSA_DEBUG("problems: fit value not set");
	  ++it;
	  continue;
	}
	if ((*it)->averageError()==0.0) {
	  // ORSA_DEBUG("--SKIPPING--");
	  ++it;
	  continue;
	}
      	fprintf(fple,"%g %g %g %g\n",
		orsa::radToDeg()*(*it)->center,
		(*it)->fit.getRef(),
		(*it)->average(),
		(*it)->averageError());
	++it;
      }
      fclose(fple);
    }
    
    {
      FILE * fpam = fopen("am.fit.dat","w");
      AllStat::const_iterator it = stat_AM.begin();
      while (it != stat_AM.end()) {
	if (!(*it)->fit.isSet()) {
	  ORSA_DEBUG("problems: fit value not set");
	  ++it;
	  continue;
	}
	if ((*it)->averageError()==0.0) {
	  // ORSA_DEBUG("--SKIPPING--");
	  ++it;
	  continue;
	}
      	fprintf(fpam,"%g %g %g %g\n",
		orsa::radToDeg()*(*it)->center,
		(*it)->fit.getRef(),
		(*it)->average(),
		(*it)->averageError());
	++it;
      }
      fclose(fpam);
    }
    
  }
  
  exit(0);
}
