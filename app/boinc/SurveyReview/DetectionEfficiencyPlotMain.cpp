#include "skycoverage.h"
#include "eta.h"
#include "fit.h"

int main(int argc, char ** argv) {
  
  orsa::Debug::instance()->initTimer();
  
  if (argc != 3) {
    printf("Usage: %s <allEta-file> <fit-file>\n",argv[0]);
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
    double solarAltitude, lunarAltitude;
    double lunarPhase;
    double airMass, azimuth;
    double galacticLongitude, galacticLatitude;
    double activeTime;
    char id[1024];
    int epochFromField;
    int observed;
    int discovered;
    while (fgets(line,1024,fp)) {
      sscanf(line,"%lf %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %i %i %i",
	     &H,
	     id,
	     &V,
	     &apparentVelocity,
	     &solarElongation,
	     &lunarElongation,
	     &solarAltitude,
	     &lunarAltitude,
	     &lunarPhase,
	     &airMass,
	     &azimuth,
	     &galacticLongitude,
	     &galacticLatitude,
	     &activeTime,
	     &epochFromField,
	     &observed,
	     &discovered);
      // convert
      apparentVelocity  = orsa::FromUnits(apparentVelocity*orsa::arcsecToRad(),orsa::Unit::HOUR,-1);
      solarElongation   = orsa::degToRad()*solarElongation;
      lunarElongation   = orsa::degToRad()*lunarElongation;
      solarAltitude     = orsa::degToRad()*solarAltitude;
      lunarAltitude     = orsa::degToRad()*lunarAltitude;
      lunarPhase        = orsa::degToRad()*lunarPhase;
      azimuth           = orsa::degToRad()*azimuth;
      galacticLongitude = orsa::degToRad()*galacticLongitude;
      galacticLatitude  = orsa::degToRad()*galacticLatitude;
      activeTime        = orsa::FromUnits(activeTime,orsa::Unit::HOUR);
      //
      EfficiencyData ed;
      ed.H=H;
      if (atoi(id) != 0) {
	ed.number = atoi(id);
      } else {
	ed.designation = id;
      }
      ed.V = V;
      ed.apparentVelocity  = apparentVelocity;
      ed.solarElongation   = solarElongation;
      ed.lunarElongation   = lunarElongation;
      ed.solarAltitude     = solarAltitude;
      ed.lunarAltitude     = lunarAltitude;
      ed.lunarPhase        = lunarPhase;
      ed.airMass           = airMass;
      ed.azimuth           = azimuth;
      ed.galacticLongitude = galacticLongitude;
      ed.galacticLatitude  = galacticLatitude;
      ed.activeTime        = activeTime;
      ed.epochFromField    = (epochFromField==1);
      ed.observed          = (observed==1);
      ed.discovered        = (discovered==1);
      //
      etaData.push_back(ed);
    }
    fclose(fp);
  }
  
  ORSA_DEBUG("etaData.size(): %i",etaData.size());
  
  // alternative method using CountStats
  std::vector< osg::ref_ptr<CountStats::Var> > varDefinition;
  //
  // [0] apparent magnitude
  osg::ref_ptr<CountStats::LinearVar> var_V = new CountStats::LinearVar(start_V,stop_V,step_V);
  varDefinition.push_back(var_V.get());
  
  // [1] apparent velocity
  osg::ref_ptr<CountStats::LinearVar> var_U = new CountStats::LinearVar(start_U,stop_U,step_U);
  varDefinition.push_back(var_U.get());
  
  // [2] solar elongation
  osg::ref_ptr<CountStats::LinearVar> var_SE = new CountStats::LinearVar(start_SE,stop_SE,step_SE);
  varDefinition.push_back(var_SE.get());
  
  // [3] lunar elongation
  osg::ref_ptr<CountStats::LinearVar> var_LE = new CountStats::LinearVar(start_LE,stop_LE,step_LE);
  varDefinition.push_back(var_LE.get());
  
  // [4] airmass
  osg::ref_ptr<CountStats::LinearVar> var_AM = new CountStats::LinearVar(start_AM,stop_AM,step_AM);
  varDefinition.push_back(var_AM.get());
  
  // [5] galactic longitude
  osg::ref_ptr<CountStats::LinearVar> var_GL = new CountStats::LinearVar(start_GL,stop_GL,step_GL);
  varDefinition.push_back(var_GL.get());
  
  // [6] galactic latitude
  osg::ref_ptr<CountStats::LinearVar> var_GB = new CountStats::LinearVar(start_GB,stop_GB,step_GB);
  varDefinition.push_back(var_GB.get());
  
  osg::ref_ptr<CountStats> countStats = 
    new CountStats(varDefinition);
  
  std::vector<double> xVector;
  xVector.resize(varDefinition.size());
  for (unsigned int k=0; k<etaData.size(); ++k) {
    // keep vars aligned with varDefinition content
    xVector[0] = etaData[k].V.getRef();
    xVector[1] = etaData[k].apparentVelocity.getRef();
    xVector[2] = etaData[k].solarElongation.getRef();
    xVector[3] = etaData[k].lunarElongation.getRef();
    xVector[4] = etaData[k].airMass.getRef();
    xVector[5] = etaData[k].galacticLongitude.getRef();
    xVector[6] = etaData[k].galacticLatitude.getRef();
    countStats->insert(xVector,
		       etaData[k].observed.getRef(),
		       etaData[k].discovered.getRef());
  }
  
  EfficiencyMultifit::DataStorage data;
  
  // for (unsigned int index=0; index<countStats->size(); ++index) {
  {
    const CountStats::DataType & csData = countStats->getData();
    CountStats::DataType::const_iterator it = csData.begin();
    while (it != csData.end()) {
      const CountStatsElement * cs = (*it).second.get();
      if (cs==0) {
	++it;
	continue;
      }
      if (cs->Ntot!=0) {
	const double      eta = (double)(cs->Nobs)/(double)(cs->Ntot);
	const double sigmaEta = (double)(sqrt(cs->Nobs+1))/(double)(cs->Ntot); // Poisson counting statistics; using Nobs+1 instead of Nobs to have positive sigma even when Nobs=0
	
	const std::vector<double> xVector = countStats->binCenterVector((*it).first); 
	
	EfficiencyMultifit::DataElement el;
	//
	el.V =xVector[0];
	el.U =xVector[1];		 
	el.SE=xVector[2];	 
	el.LE=xVector[3];
	// add lunar phase here
	el.AM=xVector[4];
	el.GL=xVector[5];
	el.GB=xVector[6];
	//
	el.eta=eta;
	el.sigmaEta=sigmaEta;
	//
	el.Nobs=cs->Nobs;
	el.Ndsc=cs->Ndsc;
	el.Ntot=cs->Ntot;
	//
	data.push_back(el);
      }
      ++it;
    }
  }
  
  osg::ref_ptr<orsaInputOutput::MPCObsCodeFile> obsCodeFile = new orsaInputOutput::MPCObsCodeFile;
  obsCodeFile->setFileName("obscode.dat");
  obsCodeFile->read();
  
  char plotFilename[1024];
  sprintf(plotFilename,"%s.fit.png",basename.c_str());
  
  // read fit file
  double JD, year;
  double V_limit, eta0_V, c_V, w_V;
  double U_limit, w_U;
  double beta;
  double GB_limit, w_GB;
  double Gmix;
  double chisq_dof;
  unsigned int Nobs, Ndsc, Ntot;
  double degSq;
  double V0;
  char jobID[1024];
  //
  {
    FILE * fp = fopen(argv[2],"r");
    if (!fp) {
      ORSA_DEBUG("cannot open file [%s]",argv[1]);
      exit(0);
    }
    
    char line[1024];
    while (fgets(line,1024,fp)) {
      if (line[0]=='#') continue; // comment
      if (19 == sscanf(line,
		       "%lf %lf %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %i %i %i %lf %lf %s",
		       &JD,
		       &year,
		       &V_limit,
		       &eta0_V,
		       &c_V,
		       &w_V,
		       &U_limit,
		       &w_U,
		       &beta,
		       &GB_limit,
		       &w_GB,
		       &Gmix,
		       &chisq_dof,
		       &Nobs,
		       &Ndsc,
		       &Ntot,
		       &degSq,
		       &V0,
		       jobID)) {
	// conversion
	U_limit   = orsa::FromUnits(U_limit*orsa::arcsecToRad(),orsa::Unit::HOUR,-1);
	w_U       = orsa::FromUnits(    w_U*orsa::arcsecToRad(),orsa::Unit::HOUR,-1);
	beta     *= orsa::degToRad();
	GB_limit *= orsa::degToRad();
	w_GB     *= orsa::degToRad();
	
	break;
      }
    }
    
    fclose(fp);
  }
  
  if (1) {
    // output for plots
    Histo<CountStats::LinearVar> histo_V (var_V.get());
    Histo<CountStats::LinearVar> histo_U (var_U.get());
    Histo<CountStats::LinearVar> histo_SE(var_SE.get());
    Histo<CountStats::LinearVar> histo_LE(var_LE.get());
    Histo<CountStats::LinearVar> histo_AM(var_AM.get());
    Histo<CountStats::LinearVar> histo_GL(var_GL.get());
    Histo<CountStats::LinearVar> histo_GB(var_GB.get());
    
    for (unsigned int k=0; k<data.size(); ++k) {
      
      /* osg::ref_ptr<EfficiencyStatistics> sV;
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
      */
      
      /* osg::ref_ptr<EfficiencyStatistics> sU;
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
      */
      
      // solarElongation
      /* osg::ref_ptr<EfficiencyStatistics> sSE;
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
      */
      
      // lunarElongation
      /* osg::ref_ptr<EfficiencyStatistics> sLE;
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
      */
      
      // airMass
      /* osg::ref_ptr<EfficiencyStatistics> sAM;
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
      */
      
      // galacticLatitude
      /* osg::ref_ptr<EfficiencyStatistics> sGL;
	 {
	 AllStat::const_iterator it = stat_GL.begin();
	 while (it != stat_GL.end()) {
	 if ((*it).get() != 0) {
	 if ((*it)->center == data[k].GL.getRef()) {
	 sGL = (*it).get();
	 }   
	 }
	 ++it;
	 }
	 if (sGL.get() == 0) {
	 sGL = new EfficiencyStatistics(data[k].GL.getRef());
	 stat_GL.push_back(sGL);
	 }
	 }
      */
      
      const double eta = SkyCoverage::eta(data[k].V.getRef(),
					  V_limit,
					  eta0_V,
					  V0,
					  c_V,
					  w_V,
					  data[k].U.getRef(),
					  U_limit,
					  w_U,
					  beta,
					  data[k].GL.getRef(),
					  data[k].GB.getRef(),
					  GB_limit,
					  w_GB,
					  Gmix);
      
      const double nominal_eta_V = SkyCoverage::nominal_eta_V(data[k].V.getRef(),
							      V_limit,
							      eta0_V,
							      V0,
							      c_V,
							      w_V);
      
      const double nominal_eta_U = SkyCoverage::nominal_eta_U(data[k].U.getRef(),
							      U_limit,
							      w_U);
      
      const double nominal_eta_GB = SkyCoverage::nominal_eta_GB(data[k].GL.getRef(),
								data[k].GB.getRef(),
								GB_limit,
								w_GB,
								Gmix);
      
      // at this point, small sigmas can become very large because divided by a very small eta's
      /* if (data[k].sigmaEta.getRef() > 0) {
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
	 if (nominal_eta_GL != 0) {
	 sGL->insert(data[k].eta.getRef()*nominal_eta_GL/eta,
	 orsa::square(eta/(nominal_eta_GL*data[k].sigmaEta.getRef())));
	 }
	 }
	 if ( !sV->fit.isSet())  sV->fit = nominal_eta_V;
	 if ( !sU->fit.isSet())  sU->fit = nominal_eta_U;
	 if (!sSE->fit.isSet()) sSE->fit = 1.0;
	 if (!sLE->fit.isSet()) sLE->fit = 1.0;
	 if (!sAM->fit.isSet()) sAM->fit = 1.0;
	 if (!sGL->fit.isSet()) sGL->fit = nominal_eta_GL;
	 } 
      */
      //
      if (data[k].sigmaEta.getRef() > 0) {
	if (eta != 0) {
	  if (nominal_eta_V != 0) {
	    histo_V.insert(data[k].V.getRef(),
			   data[k].eta.getRef()*nominal_eta_V/eta,
			   orsa::square(eta/(nominal_eta_V*data[k].sigmaEta.getRef())));
	  }
	  if (nominal_eta_U != 0) {
	    histo_U.insert(data[k].U.getRef(),
			   data[k].eta.getRef()*nominal_eta_U/eta,
			   orsa::square(eta/(nominal_eta_U*data[k].sigmaEta.getRef())));
	  }
	  histo_SE.insert(data[k].SE.getRef(),
			  data[k].eta.getRef()/eta,
			  orsa::square(eta/(data[k].sigmaEta.getRef())));
	  histo_LE.insert(data[k].LE.getRef(),
			  data[k].eta.getRef()/eta,
			  orsa::square(eta/(data[k].sigmaEta.getRef())));
	  histo_AM.insert(data[k].AM.getRef(),
			  data[k].eta.getRef()/eta,
			  orsa::square(eta/(data[k].sigmaEta.getRef())));
	  histo_GL.insert(data[k].GL.getRef(),
			  data[k].eta.getRef()/eta,
			  orsa::square(eta/(data[k].sigmaEta.getRef())));
	  if (nominal_eta_GB != 0) {
	    histo_GB.insert(data[k].GB.getRef(),
			    data[k].eta.getRef()*nominal_eta_GB/eta,
			    orsa::square(eta/(nominal_eta_GB*data[k].sigmaEta.getRef())));
	  }
	}
      }
    }
    
    // sort lists
    /* stat_V.sort(EfficiencyStatistics_ptr_cmp());
       stat_U.sort(EfficiencyStatistics_ptr_cmp());
       stat_SE.sort(EfficiencyStatistics_ptr_cmp());
       stat_LE.sort(EfficiencyStatistics_ptr_cmp());
       stat_AM.sort(EfficiencyStatistics_ptr_cmp());
       stat_GB.sort(EfficiencyStatistics_ptr_cmp());
    */
    
    /* {
       FILE * fpv = fopen("v.fit.dat","w");
       AllStat::const_iterator it = stat_V.begin();
       while (it != stat_V.end()) {
       if (!(*it)->fit.isSet()) {
       ORSA_DEBUG("problems: fit value not set");
       ++it;
       continue;
       }
       if ((*it)->standardDeviation()==0.0) {
       // ORSA_DEBUG("--SKIPPING--");
       ++it;
       continue;
       }
       fprintf(fpv,"%g %g %g %g\n",
       (*it)->center,
       (*it)->fit.getRef(),
       (*it)->average(),
       (*it)->standardDeviation());
       ++it;
       }
       fclose(fpv);
       }
    */
    
    {
      char filename[1024];
      sprintf(filename,"%s.histo.V.dat",basename.c_str());
      FILE * fp = fopen("V.fit.dat","w");
      Histo<CountStats::LinearVar>::HistoDataType::const_iterator it = histo_V.getData().begin();
      while (it != histo_V.getData().end()) {
	if ((*it)->standardDeviation()==0.0) {
	  // ORSA_DEBUG("--SKIPPING--");
	  ++it;
	  continue;
	}
	fprintf(fp,"%g %g %g\n",
		(*it)->center,
		(*it)->average(),
		(*it)->standardDeviation());
	++it;
      }
      fclose(fp);
    }
    
    {
      char filename[1024];
      sprintf(filename,"%s.histo.U.dat",basename.c_str());
      FILE * fp = fopen("U.fit.dat","w");
      Histo<CountStats::LinearVar>::HistoDataType::const_iterator it = histo_U.getData().begin();
      while (it != histo_U.getData().end()) {
	if ((*it)->standardDeviation()==0.0) {
	  // ORSA_DEBUG("--SKIPPING--");
	  ++it;
	  continue;
	}
	fprintf(fp,"%g %g %g\n",
		orsa::FromUnits((*it)->center*orsa::radToArcsec(),orsa::Unit::HOUR), // arcsec/hour
 		(*it)->average(),
		(*it)->standardDeviation());
	++it;
      }
      fclose(fp);
    }
    
    {
      char filename[1024];
      sprintf(filename,"%s.histo.GB.dat",basename.c_str());
      FILE * fp = fopen("GB.fit.dat","w");
      Histo<CountStats::LinearVar>::HistoDataType::const_iterator it = histo_GB.getData().begin();
      while (it != histo_GB.getData().end()) {
	if ((*it)->standardDeviation()==0.0) {
	  // ORSA_DEBUG("--SKIPPING--");
	  ++it;
	  continue;
	}
	fprintf(fp,"%g %g %g\n",
		orsa::radToDeg()*(*it)->center,
		(*it)->average(),
		(*it)->standardDeviation());
	++it;
      }
      fclose(fp);
    }
    
    /* 
       {
       FILE * fpu = fopen("u.fit.dat","w");
       AllStat::const_iterator it = stat_U.begin();
       while (it != stat_U.end()) {
       if (!(*it)->fit.isSet()) {
       ORSA_DEBUG("problems: fit value not set");
       ++it;
       continue;
       }
       if ((*it)->standardDeviation()==0.0) {
       // ORSA_DEBUG("--SKIPPING--");
       ++it;
       continue;
       }
       fprintf(fpu,"%g %g %g %g\n",
       orsa::FromUnits((*it)->center*orsa::radToArcsec(),orsa::Unit::HOUR), // arcsec/hour
       (*it)->fit.getRef(),
       (*it)->average(),
       (*it)->standardDeviation());
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
       if ((*it)->standardDeviation()==0.0) {
       // ORSA_DEBUG("--SKIPPING--");
       ++it;
       continue;
       }
       fprintf(fpse,"%g %g %g %g\n",
       orsa::radToDeg()*(*it)->center,
       (*it)->fit.getRef(),
       (*it)->average(),
       (*it)->standardDeviation());
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
       if ((*it)->standardDeviation()==0.0) {
       // ORSA_DEBUG("--SKIPPING--");
       ++it;
       continue;
       }
       fprintf(fple,"%g %g %g %g\n",
       orsa::radToDeg()*(*it)->center,
       (*it)->fit.getRef(),
       (*it)->average(),
       (*it)->standardDeviation());
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
       if ((*it)->standardDeviation()==0.0) {
       // ORSA_DEBUG("--SKIPPING--");
       ++it;
       continue;
       }
       fprintf(fpam,"%g %g %g %g\n",
       (*it)->center,
       (*it)->fit.getRef(),
       (*it)->average(),
       (*it)->standardDeviation());
       ++it;
       }
       fclose(fpam);
       }
       
       {
       FILE * fpgl = fopen("gl.fit.dat","w");
       AllStat::const_iterator it = stat_GB.begin();
       while (it != stat_GB.end()) {
       if (!(*it)->fit.isSet()) {
       ORSA_DEBUG("problems: fit value not set");
       ++it;
       continue;
       }
       if ((*it)->standardDeviation()==0.0) {
       // ORSA_DEBUG("--SKIPPING--");
       ++it;
       continue;
       }
       fprintf(fpgl,"%g %g %g %g\n",
       orsa::radToDeg()*(*it)->center,
       (*it)->fit.getRef(),
       (*it)->average(),
       (*it)->standardDeviation());
       ++it;
       }
       fclose(fpgl);
       }
    */
    
  }
  
  exit(0);
}
