#include "skycoverage.h"
#include "eta.h"
#include "fit.h"

int main(int argc, char ** argv) {
    
    orsa::Debug::instance()->initTimer();
    
    if (argc == 1) {
        printf("Usage: %s <allEta-file(s)>\n",argv[0]);
        exit(0);
    }
    
    ORSA_DEBUG("pid: %i",getpid());

    unsigned int numFiles=argc-1;
    
    std::vector< std::string > basename;
    basename.resize(numFiles);
    std::vector< std::vector <EfficiencyData> > etaData;
    etaData.resize(numFiles);
    
    for (unsigned int fileID=0; fileID<numFiles; ++fileID) {
        
        basename[fileID] = SkyCoverage::basename(argv[fileID+1]);
        
        readEfficiencyDataFile(etaData[fileID],
                               argv[fileID+1]);
        
        ORSA_DEBUG("etaData[%06i].size(): %i",fileID,etaData[fileID].size());
    }
    
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
    
    // [2] airmass
    osg::ref_ptr<CountStats::LinearVar> var_AM = new CountStats::LinearVar(start_AM,stop_AM,step_AM);
    varDefinition.push_back(var_AM.get());
    
    // [3] galactic latitude
    osg::ref_ptr<CountStats::LinearVar> var_GB = new CountStats::LinearVar(start_GB,stop_GB,step_GB);
    varDefinition.push_back(var_GB.get());
    
    // [4] galactic longitude
    osg::ref_ptr<CountStats::LinearVar> var_GL = new CountStats::LinearVar(start_GL,stop_GL,step_GL);
    varDefinition.push_back(var_GL.get());
    
    // [5] ecliptic latitude
    osg::ref_ptr<CountStats::LinearVar> var_EB = new CountStats::LinearVar(start_EB,stop_EB,step_EB);
    varDefinition.push_back(var_EB.get());
    
    // [6] ecliptic longitude
    osg::ref_ptr<CountStats::LinearVar> var_EL = new CountStats::LinearVar(start_EL,stop_EL,step_EL);
    varDefinition.push_back(var_EL.get());
    
    std::vector< osg::ref_ptr<CountStats> > countStats;
    countStats.resize(numFiles);
    // allocated in loop below
    
    for (unsigned int fileID=0; fileID<numFiles; ++fileID) {

        countStats[fileID] = new CountStats(varDefinition);

        std::vector<double> xVector;
        xVector.resize(varDefinition.size());
        for (unsigned int k=0; k<etaData[fileID].size(); ++k) {
            // keep vars aligned with varDefinition content
            xVector[0] = etaData[fileID][k].V.getRef();
            xVector[1] = etaData[fileID][k].apparentVelocity.getRef();
            // xVector[ ] = etaData[fileID][k].solarElongation.getRef();
            // xVector[ ] = etaData[fileID][k].lunarElongation.getRef();
            xVector[2] = etaData[fileID][k].airMass.getRef();
            xVector[3] = etaData[fileID][k].galacticLatitude.getRef();
            xVector[4] = etaData[fileID][k].galacticLongitude.getRef();
            xVector[5] = etaData[fileID][k].eclipticLatitude.getRef();
            xVector[6] = etaData[fileID][k].eclipticLongitude.getRef();
            const bool goodInsert = countStats[fileID]->insert(xVector,
                                                               etaData[fileID][k].observed.getRef(),
                                                               etaData[fileID][k].discovered.getRef());
        }
    }
    
    std::vector<EfficiencyMultifit::DataStorage> data;
    data.resize(numFiles);
    
    for (unsigned int fileID=0; fileID<numFiles; ++fileID) {

        const CountStats::DataType & csData = countStats[fileID]->getData();
        CountStats::DataType::const_iterator it = csData.begin();
        while (it != csData.end()) {
            const CountStatsElement * cs = (*it).second.get();
            if (cs==0) {
                ++it;
                continue;
            }
            if (cs->Ntot!=0) {
                // const double      eta = (double)(cs->Nobs)/(double)(cs->Ntot);
                // const double sigmaEta = (double)(sqrt(cs->Nobs+1))/(double)(cs->Ntot); // Poisson counting statistics; using Nobs+1 instead of Nobs to have positive sigma even when Nobs=0
                // const double sigmaEta = (double)(sqrt(cs->Ntot))/(double)(cs->Ntot);
                // #warning fix sigmaEta!
                // 
                // p = probability of binomial distribution
                // const double      eta = p*cs->Ntot;
                const double      eta = (double)(cs->Nobs)/(double)(cs->Ntot);
                const double        p = (double)(cs->Nobs+1)/(double)(cs->Ntot+2);
                const double sigmaEta = sqrt(p*(1-p)/cs->Ntot);
                
                const std::vector<double> xVector = countStats[fileID]->binCenterVector((*it).first); 
                
                EfficiencyMultifit::DataElement el;
                //
                el.V =xVector[0];
                el.U =xVector[1];		 
                // el.SE=xVector[ ];	 
                // el.LE=xVector[ ];
                // add lunar phase here
                el.AM=xVector[2];
                el.GB=xVector[3];
                el.GL=xVector[4];
                el.EB=xVector[5];
                el.EL=xVector[6];
                //
                el.eta=eta;
                el.sigmaEta=sigmaEta;
                //
                el.Nobs=cs->Nobs;
                el.Ndsc=cs->Ndsc;
                el.Ntot=cs->Ntot;
                //
                data[fileID].push_back(el);
                // ORSA_DEBUG("set el.AM to %g",el.AM.getRef());
            }
            ++it;
        }
    }
    
    {
        // debug
        ORSA_DEBUG("data.size(): %i",data.size());
        for (unsigned int l=0; l<data.size(); ++l) {
            ORSA_DEBUG("data[%i].size(): %i",l,data[l].size());
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

    std::vector<std::string> fitFilename;
    fitFilename.resize(numFiles);
    for (unsigned int fileID=0; fileID<numFiles; ++fileID) {
        char tmp_str[1024];
        sprintf(tmp_str,"%s.fit.dat",basename[fileID].c_str());
        fitFilename[fileID] = tmp_str;
    }
    
    /* std::vector<unsigned int> non_zero_n;
       non_zero_n.resize(argc);
       for (unsigned int fileID=0; fileID<numFiles; ++fileID) {
       non_zero_n[fileID]=0;
       for (unsigned int k=0; k<data.size(); ++k) {
       if (data[fileID][k].eta.getRef()>0.0) {
       ++non_zero_n[fileID];
       }	
       }
       ORSA_DEBUG("non_zero_n[%06i]: %i",fileID,non_zero_n[fileID]);
       }
    */
    
    osg::ref_ptr<EfficiencyMultifit> etaFit= new EfficiencyMultifit;
    //
    etaFit->maxIter=16;
    etaFit->epsabs=1.0e-3;
    etaFit->epsrel=1.0e-3;
    
    // multifit parameters  
    osg::ref_ptr<orsa::MultifitParameters> par = new orsa::MultifitParameters;
    //
    {
        const double initial_U_limit = orsa::FromUnits( 8.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1); // arcsec/hour
        const double initial_w_U     = orsa::FromUnits( 1.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1); // arcsec/hour
        //
        char varName[1024];
        for (unsigned int fileID=0; fileID<numFiles; ++fileID) {
            // V pars
            sprintf(varName,"V_limit_%06i",fileID); par->insert(varName, 19.00, 0.0000001);
            sprintf(varName, "eta0_V_%06i",fileID); par->insert(varName,  1.00, 0.0000001);
            sprintf(varName,    "c_V_%06i",fileID); par->insert(varName,  0.00, 0.0000001);
            sprintf(varName,    "w_V_%06i",fileID); par->insert(varName,  1.00, 0.0000001);
            // U pars
            sprintf(varName,"U_limit_%06i",fileID); par->insert(varName, initial_U_limit, 0.0000001*initial_U_limit);
            sprintf(varName,    "w_U_%06i",fileID); par->insert(varName,     initial_w_U, 0.0000001*initial_w_U); 
            
            // hard limits
            // sprintf(varName, "eta0_V_%06i",fileID); par->setRange(varName, 0.00, 1.00);
            // sprintf(varName,    "c_V_%06i",fileID); par->setRangeMin(varName, 0.00);
            
            /* 
               {
               // test: neutral V,U,
               sprintf(varName,"V_limit_%06i",fileID); par->set(varName, 30.00); par->setFixed(varName,true);
               sprintf(varName, "eta0_V_%06i",fileID); par->set(varName,  1.00); par->setFixed(varName,true);
               sprintf(varName,    "c_V_%06i",fileID); par->set(varName,  0.00); par->setFixed(varName,true);
               sprintf(varName,    "w_V_%06i",fileID); par->set(varName,  0.01); par->setFixed(varName,true);
               //
               sprintf(varName,    "w_U_%06i",fileID); par->set(varName, orsa::FromUnits( 0.01*orsa::arcsecToRad(),orsa::Unit::HOUR,-1)); par->setFixed(varName,true);
               sprintf(varName,"U_limit_%06i",fileID); par->set(varName, orsa::FromUnits( 0.01*orsa::arcsecToRad(),orsa::Unit::HOUR,-1)); par->setFixed(varName,true);
               }
            */
        }
        
        // AM
        par->insert("peak_AM", 1.00, 0.0000001);
        par->insert("scale_AM",1.00, 0.0000001);
        par->insert("shape_AM",0.10, 0.0000001);
        // GB
        par->insert("drop_GB",   0.20, 0.0000001);
        par->insert("scale_GB", 50.0*orsa::degToRad(), 0.0000001*orsa::degToRad()); // rad
        par->insert("center_GB", 0.0*orsa::degToRad(), 0.0000001*orsa::degToRad()); // rad
        // GL
        par->insert("scale_GL",100.0*orsa::degToRad(), 0.0000001*orsa::degToRad()); // rad
        par->insert("shape_GL",  0.1, 0.0000001);
        // EB
        par->insert("drop_EB",   0.20, 0.0000001);
        par->insert("scale_EB", 50.0*orsa::degToRad(), 0.0000001*orsa::degToRad()); // rad
        par->insert("center_EB", 0.0*orsa::degToRad(), 0.0000001*orsa::degToRad()); // rad
        // EL
        par->insert("scale_EL",100.0*orsa::degToRad(), 0.0000001*orsa::degToRad()); // rad
        par->insert("shape_EL",  0.1, 0.0000001);
        
        // hard limits
        par->setRange("drop_GB",
                      0.01,
                      1.00);
        par->setRange("scale_GB",
                      05.0*orsa::degToRad(),
                      50.0*orsa::degToRad());
        par->setRange("center_GB",
                      -7.0*orsa::degToRad(),
                      +7.0*orsa::degToRad());
        par->setFixed("center_GB",true);
        // par->setFixed("shape_GL",true);
        
        // hard limits
        par->setRange("drop_EB",
                      0.01,
                      1.00);
        par->setRange("scale_EB",
                      05.0*orsa::degToRad(),
                      50.0*orsa::degToRad());
        par->setRange("center_EB",
                      -7.0*orsa::degToRad(),
                      +7.0*orsa::degToRad());
        par->setFixed("center_EB",true);
        // par->setFixed("shape_GL",true);
        
        /* 
           {
           // test: neutral AM
           par->set("peak_AM", 100.00); par->setFixed("peak_AM", true);
           par->set("scale_AM",100.00); par->setFixed("scale_AM",true);
           par->set("shape_AM",  0.00); par->setFixed("shape_AM",true);
           }
        */
        
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
        //
        // par->setFixed("c_AM",true);
        //
        /* {
           par->setFixed("peak_AM",true);
           par->set("scale_AM",1.0e3);
           par->setFixed("scale_AM",true);
           par->setFixed("shape_AM",true);
           //
           par->setFixed("drop_GB",true);
           par->setFixed("scale_GB",true);
           }
        */        
    }
    
    osg::ref_ptr<orsa::MultifitParameters> lastGoodPar = new orsa::MultifitParameters;
    orsa::Cache<double> V0;
    bool success;
    
    {
        V0 = 16.0;
        /* while (V0.getRef() <= 21.0) {
        // reset initial values
        par->set("V_limit",19.00);
        par->set("eta0_V",  0.99);
        par->set("c_V",     0.001);
        par->set("w_V",     0.3);
        par->set("U_limit",    initial_U_limit);
        par->set("w_U",        initial_w_U);
        
        ORSA_DEBUG("V0: %g",V0.getRef());
        */
        success = etaFit->fit(par.get(),
                              data,
                              V0.getRef(),
                              fitFilename,
                              basename,
                              obsCodeFile.get(),
                              false); 
        /* if ( (success || (etaFit->getIter() == etaFit->maxIter)) &&
           (par->get("c_V") >= 0.0) ) {
           break;
           } else {
           V0 += 1.0;
           }
           }
        */
        
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
    
    {
        // now enforce range limits, and call fabs() where needed
        char varName[1024];
        for (unsigned int fileID=0; fileID<numFiles; ++fileID) {
            sprintf(varName,"V_limit_%06i",fileID); par->setRange(   varName, 10.00, 30.00);
            sprintf(varName, "eta0_V_%06i",fileID); par->setRange(   varName,  0.00,  1.00);
            sprintf(varName,    "c_V_%06i",fileID); par->setRangeMin(varName,  0.00);
            // w_V
            //
            sprintf(varName,"U_limit_%06i",fileID); par->set(varName,fabs(par->get(varName)));
            // w_U
        }
        par->setRangeMin("peak_AM",1.0);
        par->set("shape_AM",fabs(par->get("shape_AM")));
        //
        par->setRange("drop_GB",0.0,1.0);
        par->set("scale_GB",fabs(par->get("scale_GB")));
        par->setRangeMin("scale_GB",0.0);
        //
        par->set("scale_GL",fabs(par->get("scale_GL")));
        par->set("shape_GL",fabs(par->get("shape_GL")));
    }
    
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
    {
        char varName[1024];
        for (unsigned int fileID=0; fileID<numFiles; ++fileID) {
            sprintf(varName,"w_V_%06i",fileID); if (par->get(varName) < 0.0) goodFit=false;
            sprintf(varName,"w_U_%06i",fileID); if (par->get(varName) < 0.0) goodFit=false;
        }
    }
    
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
    
    {
        char varName[1024];
        for (unsigned int fileID=0; fileID<numFiles; ++fileID) {
            
            sprintf(varName,"w_V_%06i",fileID); 
            if (par->get(varName) < 0.0) {
                // just write an empty file and exit
                FILE * fp = fopen(fitFilename[fileID].c_str(),"w");
                fclose(fp);
                exit(0);
            }
            
            sprintf(varName,"w_U_%06i",fileID);
            if (par->get(varName) < 0.0) {
                // just write an empty file and exit
                FILE * fp = fopen(fitFilename[fileID].c_str(),"w");
                fclose(fp);
                exit(0);
            }
        }
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
