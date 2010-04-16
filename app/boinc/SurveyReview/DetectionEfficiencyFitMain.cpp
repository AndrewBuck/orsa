#include "skycoverage.h"
#include "eta.h"
#include "fit.h"

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
  
  ORSA_DEBUG("etaData.size(): %i",etaData.size())
  
  // minimum apparent magnitude
  // const double V0=start_V;
  
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
    const bool goodInsert = countStats->insert(xVector,
					       etaData[k].observed.getRef(),
					       etaData[k].discovered.getRef());
    if (!goodInsert) {
      // ORSA_DEBUG("problems inserting this one: V: %g",etaData[k].V.getRef());
    }
  }
  
  EfficiencyMultifit::DataStorage data;
  
  // for (unsigned int index=0; index<countStats->size(); ++index) {
  {
    const CountStats::DataType & csData = countStats->getData();
    CountStats::DataType::const_iterator it = csData.begin();
    while (it != csData.end()) {
      
      const CountStatsElement * cs = (*it).second.get();
      if (cs==0) continue;
      // write point only if Ntot != 0
      if (cs->Ntot>0) {
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
  
  // write eta.dat file
  if (0) {
    char filename[1024];
    sprintf(filename,"%s.eta.dat",basename.c_str());
    FILE * fp_eta = fopen(filename,"w"); 
    for (unsigned int k=0; k<data.size(); ++k) {
      const EfficiencyMultifit::DataElement & el = data[k];
      // if (el.Ntot.getRef()>2) {
      fprintf(fp_eta,
	      "%5.2f %6.2f %6.2f %6.2f %6.3f %5.3f %5.3f %5i %5i %5i\n",
	      el.V.getRef(),
	      orsa::FromUnits(el.U.getRef()*orsa::radToArcsec(),orsa::Unit::HOUR),
	      el.SE.getRef()*orsa::radToDeg(),
	      el.LE.getRef()*orsa::radToDeg(),
	      el.AM.getRef(),
	      el.eta.getRef(),
	      el.sigmaEta.getRef(),
	      el.Nobs.getRef(),
	      el.Ndsc.getRef(),
	      el.Ntot.getRef());
      // }
    }
    fclose(fp_eta);
  }
  
  osg::ref_ptr<orsaInputOutput::MPCObsCodeFile> obsCodeFile = new orsaInputOutput::MPCObsCodeFile;
  obsCodeFile->setFileName("obscode.dat");
  obsCodeFile->read();
  
  char fitFilename[1024];
  sprintf(fitFilename,"%s.fit.dat",basename.c_str());
  
  unsigned int non_zero_n = 0;
  for (unsigned int k=0; k<data.size(); ++k) {
    if (data[k].eta.getRef()>0.0) {
      ++non_zero_n;
    }	
  }
  ORSA_DEBUG("non_zero_n: %i",non_zero_n);
  
  osg::ref_ptr<EfficiencyMultifit> etaFit= new EfficiencyMultifit;
  //
  etaFit->maxIter=100;
  etaFit->epsabs=1.0e-6;
  etaFit->epsrel=1.0e-6;
  
  // multifit parameters  
  osg::ref_ptr<orsa::MultifitParameters> par = new orsa::MultifitParameters;
  //
  const double initial_U_limit = orsa::FromUnits(5.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1); // arcsec/hour
  const double initial_w_U     = orsa::FromUnits(1.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1); // arcsec/hour
  //    
  // V pars
  par->insert("V_limit",19.00, 0.000001);
  par->insert("eta0_V",  0.90, 0.000001);
  par->insert("c_V",     0.00, 0.000001);
  par->insert("w_V",     0.50, 0.000001); 
  // U pars
  par->insert("U_limit",    initial_U_limit, 0.000001*initial_U_limit);
  par->insert("w_U",        initial_w_U,     0.000001*initial_w_U); 
  // mixing angle
  par->insert("beta",     0.0, 0.000001);
  // Galactic latitude
  par->insert("GB_limit",10.0*orsa::degToRad(),0.000001*orsa::degToRad());
  par->insert("w_GB",     1.0*orsa::degToRad(),0.000001*orsa::degToRad());
  // Galactic longitude-latitude mixing
  par->insert("Gmix",     0.0, 0.000001);
  
  osg::ref_ptr<orsa::MultifitParameters> lastGoodPar = new orsa::MultifitParameters;
  orsa::Cache<double> V0;
  bool success;
  
  // some parameters are not used now, fixed to a "null" value
  par->setFixed("beta",true);
  par->setFixed("Gmix",true);
  
  // first, try to fit apparent magnitude only
  par->setFixed("U_limit",true);
  par->setFixed("w_U",true);
  par->setFixed("GB_limit",true);
  par->setFixed("w_GB",true);
  
  if (non_zero_n < par->sizeNotFixed()) {
    ORSA_DEBUG("not enought data for fit, exiting");
    // just create an empty output file
    FILE * fp = fopen(fitFilename,"w");
    fclose(fp);
    exit(0);
  }
  
  {
    V0 = 14.0;
    while (V0.getRef() <= 20.0) {
      ORSA_DEBUG("V0: %g",V0.getRef());
      success = etaFit->fit(par.get(),
			    data,
			    V0.getRef(),
			    fitFilename,
			    basename,
			    obsCodeFile.get(),
			    false);
      if ( (success || (etaFit->getIter() == etaFit->maxIter)) &&
	   (par->get("c_V") > 0.0) ) {
	break;
      } else {
	V0 += 1.0;
      }
    }
    if (success || (etaFit->getIter() == etaFit->maxIter)) {
      // deep copy
      (*lastGoodPar.get()) = (*par.get());
    } else {
      ORSA_DEBUG("fitting did not converge, exiting");
      // just "touch" the output file
      FILE * fp = fopen(fitFilename,"w");
      fclose(fp);
      exit(0);
    }
  }
  
  // now try including apparent velocity
  par->setFixed("U_limit",false);
  par->setFixed("w_U",false);
  
  if (non_zero_n < par->sizeNotFixed()) {
    // restore, fit saving file, end exit
    (*par.get()) = (*lastGoodPar.get());
    etaFit->fit(par.get(),
		data,
		V0.getRef(),
		fitFilename,
		basename,
		obsCodeFile.get(),
		true);
    exit(0);
  }
  
  {
    V0 = 14.0;
    while (V0.getRef() <= 20.0) {
      ORSA_DEBUG("V0: %g",V0.getRef());
      success = etaFit->fit(par.get(),
			    data,
			    V0.getRef(),
			    fitFilename,
			    basename,
			    obsCodeFile.get(),
			    false);
      if ( (success || (etaFit->getIter() == etaFit->maxIter)) &&
	   (par->get("c_V") > 0.0) ) {
	break;
      } else {
	V0 += 1.0;
      }
    }
    if (success || (etaFit->getIter() == etaFit->maxIter)) {
      // deep copy
      (*lastGoodPar.get()) = (*par.get());
    } else {
      etaFit->fit(par.get(),
		  data,
		  V0.getRef(),
		  fitFilename,
		  basename,
		  obsCodeFile.get(),
		  true);
      exit(0);
    }
  }
  
  // now try including galactic latitude
  par->setFixed("GB_limit",false);
  par->setFixed("w_GB",false);
  
  if (non_zero_n < par->sizeNotFixed()) {
    // restore, fit saving file, end exit
    (*par.get()) =  (*lastGoodPar.get());
    etaFit->fit(par.get(),
		data,
		V0.getRef(),
		fitFilename,
		basename,
		obsCodeFile.get(),
		true);
    exit(0);
  }
  
  {
    V0 = 14.0;
    while (V0.getRef() <= 20.0) {
      ORSA_DEBUG("V0: %g",V0.getRef());
      success = etaFit->fit(par.get(),
			    data,
			    V0.getRef(),
			    fitFilename,
			    basename,
			    obsCodeFile.get(),
			    false);
      if ( (success || (etaFit->getIter() == etaFit->maxIter)) &&
	   (par->get("c_V") > 0.0) ) {
	break;
      } else {
	V0 += 1.0;
      }
    }
    if (success || (etaFit->getIter() == etaFit->maxIter)) {
      // deep copy
      (*lastGoodPar.get()) = (*par.get());
    } else {
      etaFit->fit(par.get(),
		  data,
		  V0.getRef(),
		  fitFilename,
		  basename,
		  obsCodeFile.get(),
		  true);
      exit(0);
    }
  }
  
  // now enforce range limits, call fabs() where needed
  par->setRange("eta0_V",0.0,1.0);
  par->setRangeMin("c_V",0.0);
  par->set("U_limit",fabs(par->get("U_limit")));
  par->set("GB_limit",fabs(par->get("GB_limit")));
  
  {
    success = etaFit->fit(par.get(),
			  data,
			  V0.getRef(),
			  fitFilename,
			  basename,
			  obsCodeFile.get(),
			  false);
    
    if (success || (etaFit->getIter() == etaFit->maxIter)) {
      // deep copy
      (*lastGoodPar.get()) = (*par.get());
    } else {
      etaFit->fit(par.get(),
		  data,
		  V0.getRef(),
		  fitFilename,
		  basename,
		  obsCodeFile.get(),
		  true);
      exit(0);
    }
  }
  
  // finally, review fit values and save if good
  bool goodFit=true;
  if (par->get("w_V")  < 0.0) goodFit=false;
  if (par->get("w_U")  < 0.0) goodFit=false;
  if (par->get("w_GB") < 0.0) goodFit=false;
  
  // save fit file
  if (goodFit) {
    (*par.get()) = (*lastGoodPar.get());
    etaFit->fit(par.get(),
		data,
		V0.getRef(),
		fitFilename,
		basename,
		obsCodeFile.get(),
		true);
    exit(0);
  } 
  
  // if we're here, goodFit is false, so we need to set the unfittable parameters
  // so that the part of the function becomes unity
  if (par->get("w_V") < 0.0) {
    // just "touch" the output file
    FILE * fp = fopen(fitFilename,"w");
    fclose(fp);
    exit(0);
  }
  if (par->get("w_U") < 0.0) {
    par->set("U_limit",orsa::FromUnits(0.000*orsa::arcsecToRad(),orsa::Unit::HOUR,-1));
    par->set("w_U",    orsa::FromUnits(0.001*orsa::arcsecToRad(),orsa::Unit::HOUR,-1));
    par->setFixed("U_limit",false);
    par->setFixed("w_U",false);
  }
  if (par->get("w_GB") < 0.0) {
    par->set("GB_limit",0.000*orsa::degToRad());
    par->set("w_GB",    0.001*orsa::degToRad());
    par->setFixed("GB_limit",false);
    par->setFixed("w_GB",false);
  }
  
  // last run and then save fit file
  {
    etaFit->fit(par.get(),
		data,
		V0.getRef(),
		fitFilename,
		basename,
		obsCodeFile.get(),
		true);
    exit(0);
  } 
  
#warning all code below needs to be updated and moved into the plotting program
  
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
  //
  const double GB_limit = parFinal->get("GB_limit");
  const double     w_GB = parFinal->get("w_GB");
  const double     Gmix = parFinal->get("Gmix");
  
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
					  V0.getRef(),
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
							      V0.getRef(),
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
