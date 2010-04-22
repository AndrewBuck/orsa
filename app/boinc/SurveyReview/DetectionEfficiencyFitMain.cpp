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
  
    // [ ] solar elongation
    // osg::ref_ptr<CountStats::LinearVar> var_SE = new CountStats::LinearVar(start_SE,stop_SE,step_SE);
    // varDefinition.push_back(var_SE.get());
  
    // [ ] lunar elongation
    // osg::ref_ptr<CountStats::LinearVar> var_LE = new CountStats::LinearVar(start_LE,stop_LE,step_LE);
    // varDefinition.push_back(var_LE.get());
  
    // [ ] airmass
    // osg::ref_ptr<CountStats::LinearVar> var_AM = new CountStats::LinearVar(start_AM,stop_AM,step_AM);
    // varDefinition.push_back(var_AM.get());
  
    // [ ] galactic longitude
    /* osg::ref_ptr<CountStats::LinearVar> var_GL = new CountStats::LinearVar(start_GL,stop_GL,step_GL);
       varDefinition.push_back(var_GL.get());
    */
    
    // [ ] galactic latitude
    /* osg::ref_ptr<CountStats::LinearVar> var_GB = new CountStats::LinearVar(start_GB,stop_GB,step_GB);
       varDefinition.push_back(var_GB.get());
    */
    
    osg::ref_ptr<CountStats> countStats = 
        new CountStats(varDefinition);
  
    std::vector<double> xVector;
    xVector.resize(varDefinition.size());
    for (unsigned int k=0; k<etaData.size(); ++k) {
        // keep vars aligned with varDefinition content
        xVector[0] = etaData[k].V.getRef();
        xVector[1] = etaData[k].apparentVelocity.getRef();
        // xVector[ ] = etaData[k].solarElongation.getRef();
        // xVector[ ] = etaData[k].lunarElongation.getRef();
        // xVector[ ] = etaData[k].airMass.getRef();
        // xVector[ ] = etaData[k].galacticLongitude.getRef();
        // xVector[ ] = etaData[k].galacticLatitude.getRef();
        const bool goodInsert = countStats->insert(xVector,
                                                   etaData[k].observed.getRef(),
                                                   etaData[k].discovered.getRef());
        if (0) if (!goodInsert) {
            if (etaData[k].number.isSet()) {
                ORSA_DEBUG("problems inserting: [%i] %.1f %.1f %.1f %.1f %.1f %.1f %.1f",
                           etaData[k].number.getRef(),
                           etaData[k].V.getRef(),
                           orsa::FromUnits(etaData[k].apparentVelocity.getRef()*orsa::radToArcsec(),orsa::Unit::HOUR),
                           etaData[k].solarElongation.getRef()*orsa::radToDeg(),
                           etaData[k].lunarElongation.getRef()*orsa::radToDeg(),
                           etaData[k].airMass.getRef(),
                           etaData[k].galacticLongitude.getRef()*orsa::radToDeg(),
                           etaData[k].galacticLatitude.getRef()*orsa::radToDeg());
            } else if (etaData[k].designation.isSet()) {
                ORSA_DEBUG("problems inserting: [%s] %.1f %.1f %.1f %.1f %.1f %.1f %.1f",
                           etaData[k].designation.getRef().c_str(),
                           etaData[k].V.getRef(),
                           orsa::FromUnits(etaData[k].apparentVelocity.getRef()*orsa::radToArcsec(),orsa::Unit::HOUR),
                           etaData[k].solarElongation.getRef()*orsa::radToDeg(),
                           etaData[k].lunarElongation.getRef()*orsa::radToDeg(),
                           etaData[k].airMass.getRef(),
                           etaData[k].galacticLongitude.getRef()*orsa::radToDeg(),
                           etaData[k].galacticLatitude.getRef()*orsa::radToDeg());
            } else {
                ORSA_DEBUG("problems inserting: [anonymous] %.1f %.1f %.1f %.1f %.1f %.1f %.1f",
                           etaData[k].V.getRef(),
                           orsa::FromUnits(etaData[k].apparentVelocity.getRef()*orsa::radToArcsec(),orsa::Unit::HOUR),
                           etaData[k].solarElongation.getRef()*orsa::radToDeg(),
                           etaData[k].lunarElongation.getRef()*orsa::radToDeg(),
                           etaData[k].airMass.getRef(),
                           etaData[k].galacticLongitude.getRef()*orsa::radToDeg(),
                           etaData[k].galacticLatitude.getRef()*orsa::radToDeg());
            }
        }
    
        if (0) if (goodInsert) {
            if (etaData[k].number.isSet()) {
                gmp_printf("inserting: [%i] %.1f %.1f %.1f %.1f %.1f %.1f %.1f\n",
                           etaData[k].number.getRef(),
                           etaData[k].V.getRef(),
                           orsa::FromUnits(etaData[k].apparentVelocity.getRef()*orsa::radToArcsec(),orsa::Unit::HOUR),
                           etaData[k].solarElongation.getRef()*orsa::radToDeg(),
                           etaData[k].lunarElongation.getRef()*orsa::radToDeg(),
                           etaData[k].airMass.getRef(),
                           etaData[k].galacticLongitude.getRef()*orsa::radToDeg(),
                           etaData[k].galacticLatitude.getRef()*orsa::radToDeg());
            } else if (etaData[k].designation.isSet()) {
                gmp_printf("inserting: [%s] %.1f %.1f %.1f %.1f %.1f %.1f %.1f\n",
                           etaData[k].designation.getRef().c_str(),
                           etaData[k].V.getRef(),
                           orsa::FromUnits(etaData[k].apparentVelocity.getRef()*orsa::radToArcsec(),orsa::Unit::HOUR),
                           etaData[k].solarElongation.getRef()*orsa::radToDeg(),
                           etaData[k].lunarElongation.getRef()*orsa::radToDeg(),
                           etaData[k].airMass.getRef(),
                           etaData[k].galacticLongitude.getRef()*orsa::radToDeg(),
                           etaData[k].galacticLatitude.getRef()*orsa::radToDeg());
            } else {
                gmp_printf("inserting: [anonymous] %.1f %.1f %.1f %.1f %.1f %.1f %.1f\n",
                           etaData[k].V.getRef(),
                           orsa::FromUnits(etaData[k].apparentVelocity.getRef()*orsa::radToArcsec(),orsa::Unit::HOUR),
                           etaData[k].solarElongation.getRef()*orsa::radToDeg(),
                           etaData[k].lunarElongation.getRef()*orsa::radToDeg(),
                           etaData[k].airMass.getRef(),
                           etaData[k].galacticLongitude.getRef()*orsa::radToDeg(),
                           etaData[k].galacticLatitude.getRef()*orsa::radToDeg());
            }
        }
    
        // specific cases
        if (0) if ( (etaData[k].V.getRef()>=10.0) && 
                    (etaData[k].V.getRef()< 10.1) )	 {
            if (etaData[k].number.isSet()) {
                gmp_printf("special: [%i] %.1f %.1f %.1f %.1f %.1f %.1f %.1f\n",
                           etaData[k].number.getRef(),
                           etaData[k].V.getRef(),
                           orsa::FromUnits(etaData[k].apparentVelocity.getRef()*orsa::radToArcsec(),orsa::Unit::HOUR),
                           etaData[k].solarElongation.getRef()*orsa::radToDeg(),
                           etaData[k].lunarElongation.getRef()*orsa::radToDeg(),
                           etaData[k].airMass.getRef(),
                           etaData[k].galacticLongitude.getRef()*orsa::radToDeg(),
                           etaData[k].galacticLatitude.getRef()*orsa::radToDeg());
            } else if (etaData[k].designation.isSet()) {
                gmp_printf("special: [%s] %.1f %.1f %.1f %.1f %.1f %.1f %.1f\n",
                           etaData[k].designation.getRef().c_str(),
                           etaData[k].V.getRef(),
                           orsa::FromUnits(etaData[k].apparentVelocity.getRef()*orsa::radToArcsec(),orsa::Unit::HOUR),
                           etaData[k].solarElongation.getRef()*orsa::radToDeg(),
                           etaData[k].lunarElongation.getRef()*orsa::radToDeg(),
                           etaData[k].airMass.getRef(),
                           etaData[k].galacticLongitude.getRef()*orsa::radToDeg(),
                           etaData[k].galacticLatitude.getRef()*orsa::radToDeg());
            } else {
                gmp_printf("special: [anonymous] %.1f %.1f %.1f %.1f %.1f %.1f %.1f\n",
                           etaData[k].V.getRef(),
                           orsa::FromUnits(etaData[k].apparentVelocity.getRef()*orsa::radToArcsec(),orsa::Unit::HOUR),
                           etaData[k].solarElongation.getRef()*orsa::radToDeg(),
                           etaData[k].lunarElongation.getRef()*orsa::radToDeg(),
                           etaData[k].airMass.getRef(),
                           etaData[k].galacticLongitude.getRef()*orsa::radToDeg(),
                           etaData[k].galacticLatitude.getRef()*orsa::radToDeg());
            }
        }
    
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
                // el.SE=xVector[ ];	 
                // el.LE=xVector[ ];
                // add lunar phase here
                // el.AM=xVector[ ];
                // el.GL=xVector[ ];
                // el.GB=xVector[ ];
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
  
    // debug
    /* if (1) {
       for (unsigned int k=0; k<data.size(); ++k) {
       const EfficiencyMultifit::DataElement & el = data[k];
       // if (el.Ntot.getRef()>2) {
       gmp_printf("%5.2f %6.2f %6.2f %6.2f %6.3f %5.3f %5.3f %5i %5i %5i\n",
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
       }
       }
    */
    
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
    etaFit->epsabs=1.0e-9;
    etaFit->epsrel=1.0e-9;
  
    // multifit parameters  
    osg::ref_ptr<orsa::MultifitParameters> par = new orsa::MultifitParameters;
    //
    const double initial_U_limit = orsa::FromUnits(20.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1); // arcsec/hour
    const double initial_w_U     = orsa::FromUnits( 2.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1); // arcsec/hour
    //    
    // V pars
    par->insert("V_limit",19.20, 0.000001);
    par->insert("eta0_V",  1.00, 0.000001);
    par->insert("c_V",     0.01, 0.000001);
    par->insert("w_V",     1.00, 0.000001); 
    // U pars
    par->insert("U_limit",    initial_U_limit, 0.000001*initial_U_limit);
    par->insert("w_U",        initial_w_U,     0.000001*initial_w_U); 
    // mixing angle
    // par->insert("beta",     0.0, 0.000001);
    // Galactic latitude
    // par->insert("GB_limit",0.0*orsa::degToRad(),0.000001*orsa::degToRad());
    // par->insert("w_GB",    0.001*orsa::degToRad(),0.000001*orsa::degToRad());
    // Galactic longitude-latitude mixing
    // par->insert("Gmix",     0.0, 0.000001);
    
    // try to enforce this from beginning
    // par->setRangeMin("c_V",0.0);
    /* par->setFixed("U_limit",true);
       par->setFixed("w_U",true);
       par->setFixed("V_limit",true);
       par->setFixed("eta0_V",true);
       par->setFixed("c_V",true);
    */
    
    osg::ref_ptr<orsa::MultifitParameters> lastGoodPar = new orsa::MultifitParameters;
    orsa::Cache<double> V0;
    bool success;

    if (non_zero_n < par->sizeNotFixed()) {
        ORSA_DEBUG("not enought data for fit, exiting");
        // just create an empty output file
        FILE * fp = fopen(fitFilename,"w");
        fclose(fp);
        exit(0);
    }
    
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
        V0 = 15.0;
        while (V0.getRef() <= 21.0) {
            // reset initial values
            par->set("V_limit",19.00);
            par->set("eta0_V",  0.90);
            par->set("c_V",     0.00);
            par->set("w_V",     1.00);
            par->set("U_limit",    initial_U_limit);
            par->set("w_U",        initial_w_U);
            
            ORSA_DEBUG("V0: %g",V0.getRef());
            success = etaFit->fit(par.get(),
                                  data,
                                  V0.getRef(),
                                  fitFilename,
                                  basename,
                                  obsCodeFile.get(),
                                  false);
            if ( (success || (etaFit->getIter() == etaFit->maxIter)) &&
                 (par->get("c_V") >= 0.0) ) {
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
    par->setRange("V_limit",10.00, 30.00);
    par->setRange("eta0_V",0.0,1.0);
    par->setRangeMin("c_V",0.0);
    // par->setRangeMin("w_V",0.0);
    par->set("U_limit",fabs(par->get("U_limit")));
    // par->setRangeMin("w_U",0.0);
    
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
    // if (par->get("w_GB") < 0.0) goodFit=false;
    
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
    /* if (par->get("w_GB") < 0.0) {
       par->set("GB_limit",0.000*orsa::degToRad());
       par->set("w_GB",    0.001*orsa::degToRad());
       par->setFixed("GB_limit",false);
       par->setFixed("w_GB",false);
       }
    */
    
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
  
    exit(0);
}
