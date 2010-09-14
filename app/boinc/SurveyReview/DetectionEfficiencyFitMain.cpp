#include "skycoverage.h"
#include "eta.h"
#include "fit.h"

enum FitMode {
    UV, // U,V only
    ALL // all vars
};

FitMode getFitMode(const std::string & mode) {
    if (mode=="UV") {
        return UV;
    } else if (mode=="ALL") {
        return ALL;
    } else {
        exit(0);
    }
}

int main(int argc, char ** argv) {
    
    orsa::Debug::instance()->initTimer();
    
    if (argc < 3) {
        printf("Usage: %s <UV|ALL> <allEta-file(s)>\n",argv[0]);
        exit(0);
    }
    
    ORSA_DEBUG("pid: %i",getpid());
    
    const FitMode mode=getFitMode(argv[1]);
    
    const unsigned int numFiles=argc-2;
    
    std::vector< std::string > basename;
    basename.resize(numFiles);
    std::vector< std::vector <EfficiencyData> > etaData;
    etaData.resize(numFiles);
    
    for (unsigned int fileID=0; fileID<numFiles; ++fileID) {
        
        basename[fileID] = SkyCoverage::basename(argv[fileID+2]);
        
        readEfficiencyDataFile(etaData[fileID],
                               argv[fileID+2]);
        
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
    
    // [2] airmass
    osg::ref_ptr<CountStats::LinearVar> var_AM = new CountStats::LinearVar(start_AM,stop_AM,step_AM);
    varDefinition.push_back(var_AM.get());
    
    // [3] galactic latitude
    osg::ref_ptr<CountStats::LinearVar> var_GB = new CountStats::LinearVar(start_GB,stop_GB,step_GB);
    varDefinition.push_back(var_GB.get());
    
    // [4] galactic longitude
    osg::ref_ptr<CountStats::LinearVar> var_GL = new CountStats::LinearVar(start_GL,stop_GL,step_GL);
    varDefinition.push_back(var_GL.get());
    
    // [ ] solar altitude
    /* osg::ref_ptr<CountStats::LinearVar> var_SA = new CountStats::LinearVar(start_SA,stop_SA,step_SA);
       varDefinition.push_back(var_SA.get());
    */
    
    // [ ] lunar altitude
    /* osg::ref_ptr<CountStats::LinearVar> var_LA = new CountStats::LinearVar(start_LA,stop_LA,step_LA);
       varDefinition.push_back(var_LA.get());
    */
    
    // [ ] lunar phase
    /* osg::ref_ptr<CountStats::LinearVar> var_LP = new CountStats::LinearVar(start_LP,stop_LP,step_LP);
       varDefinition.push_back(var_LP.get());
    */
    
    // [ ] lunar illumination
    /* osg::ref_ptr<CountStats::LinearVar> var_LI = new CountStats::LinearVar(start_LI,stop_LI,step_LI);
       varDefinition.push_back(var_LI.get());
    */
    
    std::vector< osg::ref_ptr<CountStats> > countStats;
    countStats.resize(numFiles);
    // allocated in loop below
    
    for (unsigned int fileID=0; fileID<numFiles; ++fileID) {
        
        countStats[fileID] = new CountStats(varDefinition);
        
        if (etaData[fileID].size()==0) continue;
        
        std::vector<double> xVector;
        xVector.resize(varDefinition.size());
        unsigned int excluded=0;
        for (unsigned int k=0; k<etaData[fileID].size(); ++k) {
            // keep vars aligned with varDefinition content
            xVector[0] = etaData[fileID][k].V.getRef();
            xVector[1] = etaData[fileID][k].apparentVelocity.getRef();
            xVector[2] = etaData[fileID][k].airMass.getRef();
            xVector[3] = etaData[fileID][k].galacticLatitude.getRef();
            xVector[4] = etaData[fileID][k].galacticLongitude.getRef();
            // xVector[ ] = etaData[fileID][k].solarAltitude.getRef();
            // xVector[ ] = etaData[fileID][k].lunarAltitude.getRef();
            // xVector[ ] = LP2LI(etaData[fileID][k].lunarPhase.getRef());
            const bool goodInsert = countStats[fileID]->insert(xVector,
                                                               etaData[fileID][k].observed.getRef(),
                                                               etaData[fileID][k].discovered.getRef());
            if (!goodInsert) {
                // KEEP THIS!
                ++excluded;
                //
                if (0) {
                    ORSA_DEBUG("excluded: V=%4.1f U=%5.1f AM=%5.2f GB=%+3.0f GL=%+4.0f obs: %i",
                               etaData[fileID][k].V.getRef(),
                               orsa::FromUnits(etaData[fileID][k].apparentVelocity.getRef()*orsa::radToArcsec(),orsa::Unit::HOUR),
                               etaData[fileID][k].airMass.getRef(),
                               orsa::radToDeg()*etaData[fileID][k].galacticLatitude.getRef(),
                               orsa::radToDeg()*etaData[fileID][k].galacticLongitude.getRef(),
                               etaData[fileID][k].observed.getRef());
                }
            }
        }
        ORSA_DEBUG("etaData[%06i]: excluded %i out of %i [%.2f\%]",fileID,excluded,etaData[fileID].size(),(100.0*excluded)/etaData[fileID].size());
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
                el.AM=xVector[2];
                el.GB=xVector[3];
                el.GL=xVector[4];
                // el.SA=xVector[ ];
                // el.LA=xVector[ ];
                // el.LI=xVector[ ];
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
    if (!obsCodeFile->read()) {
        exit(0);
    }
    
    std::vector<std::string> fitFilename;
    fitFilename.resize(numFiles);
    for (unsigned int fileID=0; fileID<numFiles; ++fileID) {
        char tmp_str[1024];
        sprintf(tmp_str,"%s.fit.dat",basename[fileID].c_str());
        fitFilename[fileID] = tmp_str;
    }
    
    osg::ref_ptr<EfficiencyMultifit> etaFit= new EfficiencyMultifit;
    //
    etaFit->maxIter=16;
    etaFit->epsabs=1.0e-3;
    etaFit->epsrel=1.0e-3;
    
    // multifit parameters  
    osg::ref_ptr<orsa::MultifitParameters> par = new orsa::MultifitParameters;
    //
    {
        const double initial_U_limit = orsa::FromUnits(10.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1); // arcsec/hour
        const double initial_w_U     = orsa::FromUnits( 1.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1); // arcsec/hour
        //
        char varName[1024];
        for (unsigned int fileID=0; fileID<numFiles; ++fileID) {
            // V pars
            sprintf(varName,"V_limit_%06i",fileID); par->insert(varName, 19.00, 1.0e-6);
            sprintf(varName, "eta0_V_%06i",fileID); par->insert(varName,  1.00, 1.0e-6);
            sprintf(varName,    "c_V_%06i",fileID); par->insert(varName,  0.00, 1.0e-6);
            sprintf(varName,    "w_V_%06i",fileID); par->insert(varName,  1.00, 1.0e-6);
            // U pars
            sprintf(varName,"U_limit_%06i",fileID); par->insert(varName, initial_U_limit, 1.0e-6*initial_U_limit);
            sprintf(varName,    "w_U_%06i",fileID); par->insert(varName,     initial_w_U, 1.0e-6*initial_w_U); 
            
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
        
        // set AM, GB, GL fits to "neutral" values (close to unity);
        
        // AM
        par->insert("peak_AM", 1.00, 1.0e-6);
        par->insert("scale_AM",100.00, 1.0e-6);
        par->insert("shape_AM",0.10, 1.0e-6);
        // GB
        par->insert("drop_GB", 0.01, 1.0e-6);
        par->insert("scale_GB", 10.0*orsa::degToRad(), 1.0e-6*orsa::degToRad()); // rad
        // GL
        par->insert("scale_GL",100.0*orsa::degToRad(), 1.0e-6*orsa::degToRad()); // rad
        par->insert("shape_GL",  0.1, 1.0e-6);
        
        // hard limits
        /* par->setRange("drop_GB",
           0.001,
           1.00);
           par->setRange("scale_GB",
           01.0*orsa::degToRad(),
           75.0*orsa::degToRad());
        */
    }
    
    {
        // use base fit file, if present, to overwrite starting values   
        const std::string filename = "common.fit.dat";
        FitFileDataElement e;
        if (readFitFile(e,filename)) {
            ORSA_DEBUG("reading fit values from existing file [%s]",filename.c_str());
            char varName[1024];
            for (unsigned int fileID=0; fileID<numFiles; ++fileID) {
                sprintf(varName,"V_limit_%06i",fileID); par->set(varName, e.V_limit);
                sprintf(varName, "eta0_V_%06i",fileID); par->set(varName, e.eta0_V);
                sprintf(varName,    "c_V_%06i",fileID); par->set(varName, e.c_V);
                sprintf(varName,    "w_V_%06i",fileID); par->set(varName, e.w_V);
                // U pars
                sprintf(varName,"U_limit_%06i",fileID); par->set(varName, e.U_limit);
                sprintf(varName,    "w_U_%06i",fileID); par->set(varName, e.w_U); 
            }
            // AM
            par->set("peak_AM",  e.peak_AM);
            par->set("scale_AM", e.scale_AM);
            par->set("shape_AM", e.shape_AM);
            // GB
            par->set("drop_GB",  e.drop_GB);
            par->set("scale_GB", e.scale_GB); // rad
            // GL
            par->set("scale_GL", e.scale_GL); // rad
            par->set("shape_GL", e.shape_GL);
        }
    }
    
    
    {
        // use fit files, if present, to overwrite starting values   
        for (unsigned int fileID=0; fileID<numFiles; ++fileID) {
            FitFileDataElement e;
            if (readFitFile(e,fitFilename[fileID])) {
                ORSA_DEBUG("reading fit values from existing file [%s]",fitFilename[fileID].c_str());
                char varName[1024];
                sprintf(varName,"V_limit_%06i",fileID); par->set(varName, e.V_limit);
                sprintf(varName, "eta0_V_%06i",fileID); par->set(varName, e.eta0_V);
                sprintf(varName,    "c_V_%06i",fileID); par->set(varName, e.c_V);
                sprintf(varName,    "w_V_%06i",fileID); par->set(varName, e.w_V);
                // U pars
                sprintf(varName,"U_limit_%06i",fileID); par->set(varName, e.U_limit);
                sprintf(varName,    "w_U_%06i",fileID); par->set(varName, e.w_U); 
                // AM
                par->set("peak_AM",  e.peak_AM);
                par->set("scale_AM", e.scale_AM);
                par->set("shape_AM", e.shape_AM);
                // GB
                par->set("drop_GB",  e.drop_GB);
                par->set("scale_GB", e.scale_GB); // rad
                // GL
                par->set("scale_GL", e.scale_GL); // rad
                par->set("shape_GL", e.shape_GL);
            }
        }
    }
    
    // switch vars ON and OFF depending on "mode"
    if (mode==UV) {
        //AM
        par->setFixed("peak_AM",true);
        par->setFixed("scale_AM",true);
        par->setFixed("shape_AM",true);
        // GB
        par->setFixed("drop_GB",true);
        par->setFixed("scale_GB",true);
        // GL
        par->setFixed("scale_GL",true);
        par->setFixed("shape_GL",true);
    }
    
    osg::ref_ptr<orsa::MultifitParameters> lastGoodPar = new orsa::MultifitParameters;
    orsa::Cache<double> V0;
    bool success;
    
    {
        V0 = 16.0;
        success = etaFit->fit(par.get(),
                              data,
                              V0.getRef(),
                              fitFilename,
                              basename,
                              obsCodeFile.get(),
                              true); // false 
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
            // sprintf(varName,"V_limit_%06i",fileID); par->setRange(   varName, 10.00, 30.00);
            // sprintf(varName, "eta0_V_%06i",fileID); par->setRange(   varName,  0.00,  1.00);
            // sprintf(varName,    "c_V_%06i",fileID); par->setRangeMin(varName,  0.00);
            // w_V
            //
            sprintf(varName,"U_limit_%06i",fileID); par->set(varName,fabs(par->get(varName)));
            // w_U
        }
        // par->setRangeMin("peak_AM",1.0);
        par->set("shape_AM",fabs(par->get("shape_AM")));
        //
        // par->setRange("drop_GB",0.0,1.0);
        par->set("scale_GB",fabs(par->get("scale_GB")));
        // par->setRangeMin("scale_GB",0.0);
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
                              true); // false
        
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
