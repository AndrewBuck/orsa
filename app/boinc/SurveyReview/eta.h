#ifndef __ETA_H_
#define __ETA_H_

#include <orsa/multifit.h>
#include <orsaInputOutput/MPC_obscode.h>

class EfficiencyData {
public:
    orsa::Cache<double> H;   
    orsa::Cache<unsigned int> number;
    orsa::Cache<std::string> designation;
    orsa::Cache<double> V;
    orsa::Cache<double> apparentVelocity;
    orsa::Cache<double> solarElongation;
    orsa::Cache<double> lunarElongation;
    orsa::Cache<double> solarAltitude;
    orsa::Cache<double> lunarAltitude;
    orsa::Cache<double> lunarPhase;
    orsa::Cache<double> airMass;
    orsa::Cache<double> azimuth;
    orsa::Cache<double> galacticLongitude;
    orsa::Cache<double> galacticLatitude;
    orsa::Cache<double> eclipticLongitude; // Sun at longitude 0
    orsa::Cache<double> eclipticLatitude; 
    orsa::Cache<double> activeTime;
    orsa::Cache<bool>   epochFromField;
    orsa::Cache<bool>   observed;
    orsa::Cache<bool>   discovered;
};

void writeEfficiencyDataFile(const std::vector<EfficiencyData> & etaData,
                             const std::string & filename) {
    ORSA_DEBUG("writing file: [%s]",filename.c_str());
    FILE * fp_allEta = fopen(filename.c_str(),"w");
    for (unsigned int k=0; k<etaData.size(); ++k) {
        const EfficiencyData & ed = etaData[k];
        char id[1024];
        if (ed.number.isSet()) {
            sprintf(id,"%i",ed.number.getRef());
        } else if (ed.designation.isSet()) {
            sprintf(id,"%s",ed.designation.getRef().c_str());
        }
        fprintf(fp_allEta,
                "%5.2f %7s %5.2f %7.2f %6.2f %6.2f %+7.2f %+7.2f %6.2f %6.3f %5.1f %+8.3f %+7.3f %+8.3f %+7.3f %5.2f %1i %1i %1i\n",
                ed.H.getRef(),
                id,
                ed.V.getRef(),
                orsa::FromUnits(ed.apparentVelocity.getRef()*orsa::radToArcsec(),orsa::Unit::HOUR), // arcsec/hour
                orsa::radToDeg()*ed.solarElongation.getRef(),
                orsa::radToDeg()*ed.lunarElongation.getRef(),
                orsa::radToDeg()*ed.solarAltitude.getRef(),
                orsa::radToDeg()*ed.lunarAltitude.getRef(),
                orsa::radToDeg()*ed.lunarPhase.getRef(),
                ed.airMass.getRef(),
                orsa::radToDeg()*ed.azimuth.getRef(),
                orsa::radToDeg()*ed.galacticLongitude.getRef(),
                orsa::radToDeg()*ed.galacticLatitude.getRef(),
                orsa::radToDeg()*ed.eclipticLongitude.getRef(),
                orsa::radToDeg()*ed.eclipticLatitude.getRef(),
                orsa::FromUnits(ed.activeTime.getRef(),orsa::Unit::HOUR,-1), 
                ed.epochFromField.getRef(),
                ed.observed.getRef(),
                ed.discovered.getRef());
    }
    fclose(fp_allEta);
}

void readEfficiencyDataFile(std::vector<EfficiencyData> & etaData,
                            const std::string & filename) {
    FILE * fp = fopen(filename.c_str(),"r");
    if (!fp) {
        ORSA_DEBUG("cannot open file [%s]",filename.c_str());
        exit(0);
    }
    char line[1024];
    double H, V, apparentVelocity;
    double solarElongation, lunarElongation;
    double solarAltitude, lunarAltitude;
    double lunarPhase;
    double airMass, azimuth;
    double galacticLongitude, galacticLatitude;
    double eclipticLongitude, eclipticLatitude;
    double activeTime;
    char id[1024];
    int epochFromField;
    int observed;
    int discovered;
    while (fgets(line,1024,fp)) {
        sscanf(line,"%lf %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %i %i %i",
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
               &eclipticLongitude,
               &eclipticLatitude,
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
        eclipticLongitude = orsa::degToRad()*eclipticLongitude;
        eclipticLatitude  = orsa::degToRad()*eclipticLatitude;
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
        ed.eclipticLongitude = eclipticLongitude;
        ed.eclipticLatitude  = eclipticLatitude;
        ed.activeTime        = activeTime;
        ed.epochFromField    = (epochFromField==1);
        ed.observed          = (observed==1);
        ed.discovered        = (discovered==1);
        //
        etaData.push_back(ed);
    }
    fclose(fp);
}

class FitFileDataElement {
public:
    double JD, year;
    double V_limit, eta0_V, c_V, w_V;
    double U_limit, w_U;
    double peak_AM, scale_AM, shape_AM;
    double drop_GB, scale_GB;
    double scale_GL, shape_GL;
    // double peak_SA, scale_SA, shape_SA;
    // double LA_LI_limit_const, LA_LI_limit_linear, LA_LI_w_const, LA_LI_w_linear;
    double chisq_dof;
    unsigned int Nobs, Ndsc, Ntot;
    double degSq;
    double V0;
    char jobID[1024];
};
typedef std::vector<FitFileDataElement> FitFileData;

bool readFitFile(FitFileDataElement & e,
                 const std::string  & filename) {
    FILE * fp = fopen(filename.c_str(),"r");
    if (fp) {
        ORSA_DEBUG("reading file: [%s]",filename.c_str());
    } else {
        ORSA_DEBUG("problems reading file: [%s]",filename.c_str());
        return false;
    }
    char line[1024];
    while (fgets(line,1024,fp)) {
        if (line[0]=='#') continue; // comment
        // UPDATE THIS NUMBER
        if (22 == sscanf(line,
                         "%lf %lf %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %*s %lf %i %i %i %lf %lf %s",
                         &e.JD,
                         &e.year,
                         &e.V_limit,
                         &e.eta0_V,
                         &e.c_V,
                         &e.w_V,
                         &e.U_limit,
                         &e.w_U,
                         &e.peak_AM,
                         &e.scale_AM,
                         &e.shape_AM,
                         &e.drop_GB,
                         &e.scale_GB,
                         &e.scale_GL,
                         &e.shape_GL,
                         &e.chisq_dof,
                         &e.Nobs,
                         &e.Ndsc,
                         &e.Ntot,
                         &e.degSq,
                         &e.V0,
                         e.jobID)) {
            // conversion
            e.U_limit   = orsa::FromUnits(e.U_limit*orsa::arcsecToRad(),orsa::Unit::HOUR,-1);
            e.w_U       = orsa::FromUnits(    e.w_U*orsa::arcsecToRad(),orsa::Unit::HOUR,-1);
            e.scale_GB  = orsa::degToRad()*e.scale_GB;
            e.scale_GL  = orsa::degToRad()*e.scale_GL;
            break;
        } else {
            ORSA_DEBUG("problems reading fit file: [%s]",filename.c_str());
            // leaving
            fclose(fp);
            return false;
        }
    }
    fclose(fp);
    return true;
}

class EfficiencyMultifit : public orsa::Multifit {
public:
    class DataElement {
    public:
        orsa::Cache<double> V, U, AM, GB, GL, SA, LA, LI; // fitting on these
        orsa::Cache<double> SE, LE, EL, EB, AZ; // for plotting only
        orsa::Cache<double> eta, sigmaEta;
        orsa::Cache<unsigned int> Nobs, Ndsc, Ntot;
        orsa::Cache<unsigned int> fileID;
    };
    typedef std::vector<DataElement> DataStorage;
public: 
    bool fit(orsa::MultifitParameters * par,
             const std::vector<DataStorage> & data_in,
             const double      & V0_in,
             const std::vector<std::string> & outputFile_in,
             const std::vector<std::string> & jobID_in,
             orsaInputOutput::MPCObsCodeFile * obsCodeFile_in,
             const bool writeFile_in) {
        
        // local copies
        data = data_in;
        V0 = V0_in;
        //
        outputFile  = outputFile_in;
        jobID       = jobID_in;
        obsCodeFile = obsCodeFile_in;
        writeFile   = writeFile_in;
        
        // basic checks
        if (data.size()==0) return false;
        
        const unsigned int numFiles = data.size();
        
        setMultifitParameters(par);
        
        osg::ref_ptr<orsa::MultifitData> fitData = new orsa::MultifitData;
        //
        fitData->insertVariable("V");
        fitData->insertVariable("U");
        fitData->insertVariable("AM");
        fitData->insertVariable("GB");
        fitData->insertVariable("GL");
        /* fitData->insertVariable("SA");
           fitData->insertVariable("LA");
           fitData->insertVariable("LI"); */
        fitData->insertVariable("fileID");
        //
        {
            unsigned int row=0;
            for (unsigned int fileID=0; fileID<numFiles; ++fileID) {
                for (unsigned int l=0; l<data[fileID].size(); ++l) {
                    fitData->insertD("V", row,data[fileID][l].V.getRef());
                    fitData->insertD("U", row,data[fileID][l].U.getRef());
                    fitData->insertD("AM",row,data[fileID][l].AM.getRef());
                    fitData->insertD("GB",row,data[fileID][l].GB.getRef());
                    fitData->insertD("GL",row,data[fileID][l].GL.getRef());
                    /* fitData->insertD("SA",row,data[fileID][l].SA.getRef());
                       fitData->insertD("LA",row,data[fileID][l].LA.getRef());
                       fitData->insertD("LI",row,data[fileID][l].LI.getRef()); */
                    //
                    fitData->insertZ("fileID",row,fileID);
                    //
                    fitData->insertF(row,data[fileID][l].eta.getRef());
                    fitData->insertSigma(row,data[fileID][l].sigmaEta.getRef());
                    //
                    ++row;
                }
            }
            
        }
        setMultifitData(fitData.get());
        
        // if (!run_nmsimplex(4096,0.001)) {
        if (!run()) {
            ORSA_WARNING("the fit did not converge.");
            return false;
        } else {
            return true;
        }
    }
protected:
    double fun(const orsa::MultifitParameters * par, 
               const orsa::MultifitData * data,
               const unsigned int p, // par index
               const int          d, // delta
               const unsigned int row) const {    
        
        osg::ref_ptr<orsa::MultifitParameters> localPar = new orsa::MultifitParameters;
        
        // deep copy
        (*localPar) = (*par);
        
        // modify
        // localPar->set(p,par->get(p)+d*par->getDelta(p));
        {
            unsigned int gslIndex=0; 
            for (unsigned int k=0; k<_par->totalSize(); ++k) {
                if (!_par->isFixed(k)) {
                    if (p==gslIndex) {
                        localPar->set(k,par->get(k)+d*par->getDelta(k));
                        break;
                    }
                    ++gslIndex;
                }
            }
        }
        
        const unsigned int fileID = data->getZ("fileID",row).get_ui();
        
        // ORSA_DEBUG("fileID: %i   p: %i   d: %i",fileID,p,d);
        
        // V
        char varName_V_limit[1024]; sprintf(varName_V_limit,"V_limit_%06i",fileID);
        char  varName_eta0_V[1024]; sprintf(varName_eta0_V,  "eta0_V_%06i",fileID);
        char     varName_c_V[1024]; sprintf(varName_c_V,        "c_V_%06i",fileID);
        char     varName_w_V[1024]; sprintf(varName_w_V,        "w_V_%06i",fileID);
        // U
        char varName_U_limit[1024]; sprintf(varName_U_limit,"U_limit_%06i",fileID);
        char     varName_w_U[1024]; sprintf(varName_w_U,        "w_U_%06i",fileID);
        
        const double eta = SkyCoverage::eta(data->getD("V",row),
                                            localPar->get(varName_V_limit),
                                            localPar->get(varName_eta0_V),
                                            V0,
                                            localPar->get(varName_c_V),
                                            localPar->get(varName_w_V),
                                            data->getD("U",row),
                                            localPar->get(varName_U_limit),
                                            localPar->get(varName_w_U),
                                            data->getD("AM",row),
                                            localPar->get("peak_AM"),
                                            localPar->get("scale_AM"),
                                            localPar->get("shape_AM"),
                                            data->getD("GB",row),
                                            localPar->get("drop_GB"),
                                            localPar->get("scale_GB"),
                                            data->getD("GL",row),
                                            localPar->get("scale_GL"),
                                            localPar->get("shape_GL"));
        return eta;
    }
protected: 
    void singleIterationDone(const gsl_multifit_fdfsolver * s) const { 
        
        {
            unsigned int gslIndex=0;
            for (unsigned int k=0; k<_par->totalSize(); ++k) {
                if (!_par->isFixed(k)) {
                    _par->set(k, gsl_vector_get(s->x,gslIndex));
                    ++gslIndex;
                }
            }
        }
        
        double c = 1.0;
        //
        const unsigned int dof = _data->size() - _par->sizeNotFixed();
        if (dof > 0) {
            const double chi = gsl_blas_dnrm2(s->f);
            c = GSL_MAX_DBL(1.0, chi / sqrt(dof)); 
            ORSA_DEBUG("chisq/dof = %g",  chi*chi/dof);
            // gmp_fprintf(fp,"chisq/dof = %g\n",  chi*chi/dof);
        }
        //
        const double factor = c;
        
        gsl_matrix * covar = gsl_matrix_alloc(_par->sizeNotFixed(),_par->sizeNotFixed());
        
        gsl_multifit_covar(s->J, 0.0, covar);
        
        // call abort() if an uncertainty gets too large, or if its zero
        if (0) {
            unsigned int gslIndex=0;
            for (unsigned int p=0; p<_par->totalSize(); ++p) {
                if (!_par->isFixed(p)) {
                    if  (factor*sqrt(gsl_matrix_get(covar,gslIndex,gslIndex)) > fabs(1.0e9*(1.0e-3+_par->get(p)))) {
                        ORSA_DEBUG("uncertainty of parameter [%s] is too large, aborting",
                                   _par->name(p).c_str());
                        abort();
                    }
                    if  (factor*sqrt(gsl_matrix_get(covar,gslIndex,gslIndex)) == 0.0) {
                        ORSA_DEBUG("uncertainty of parameter [%s] is zero, aborting",
                                   _par->name(p).c_str());
                        abort();
                    }
                    ++gslIndex;
                }
            }
        }
        
        {
            unsigned int gslIndex=0;
            for (unsigned int p=0; p<_par->totalSize(); ++p) {
                // first, specific cases where unit conversions or factors are needed
                if ( (_par->name(p).find("U_limit") != std::string::npos) || 
                     (_par->name(p).find("w_U") != std::string::npos) ) {
                    ORSA_DEBUG("%20s: %g +/- %g [arcsec/hour]",
                               _par->name(p).c_str(),
                               orsa::FromUnits(_par->get(p)*orsa::radToArcsec(),orsa::Unit::HOUR),
                               _par->isFixed(p)?0.0:orsa::FromUnits(factor*sqrt(gsl_matrix_get(covar,gslIndex,gslIndex))*orsa::radToArcsec(),orsa::Unit::HOUR));
                } else if ( (_par->name(p).find("scale_GB")  != std::string::npos) ||
                            (_par->name(p).find("scale_GL")  != std::string::npos) ) {
                    ORSA_DEBUG("%20s: %g +/- %g [deg]",
                               _par->name(p).c_str(),
                               orsa::radToDeg()*_par->get(p),
                               _par->isFixed(p)?0.0:orsa::radToDeg()*factor*sqrt(gsl_matrix_get(covar,gslIndex,gslIndex)));
                } else {
                    // generic one
                    ORSA_DEBUG("%20s: %g +/- %g",
                               _par->name(p).c_str(),
                               _par->get(p),
                               _par->isFixed(p)?0.0:factor*sqrt(gsl_matrix_get(covar,gslIndex,gslIndex)));
                }
                if (!_par->isFixed(p)) {
                    ++gslIndex;
                }
            }
        }
        
        gsl_matrix_free(covar);
        
        if (writeFile) writeOutputFiles(s);
    }
protected:
    void runCompleted(const bool /* success */, const gsl_multifit_fdfsolver * s) const {
        if (writeFile) writeOutputFiles(s);
    }
protected:
    void writeOutputFiles(const gsl_multifit_fdfsolver * s) const {
        
        const unsigned int numFiles = data.size();

        // a vector of strings: "_000000", "_000001"...
        std::vector<std::string> fileID_str_vec;
        fileID_str_vec.resize(numFiles);
        {
            char fileID_str[1024];
            for (unsigned int fileID=0; fileID<numFiles; ++fileID) {
                sprintf(fileID_str,"_%06i",fileID);
                fileID_str_vec[fileID] = fileID_str;
            }
        }
        
        for (unsigned int fileID=0; fileID<numFiles; ++fileID) {

            std::string obsCode;
            orsa::Time epoch;
            int year;
            int dayOfYear;
            if (!SkyCoverage::processFilename(jobID[fileID].c_str(),
                                              obsCodeFile,
                                              obsCode,
                                              epoch,
                                              year,
                                              dayOfYear)) {
                ORSA_DEBUG("problems...");
                // exit(0);
            }
            
            double degSq=0;
            {
                osg::ref_ptr<SkyCoverageFile> skyCoverageFile = new SkyCoverageFile;
                char filename[1024];
                sprintf(filename,"%s.txt",jobID[fileID].c_str());
                skyCoverageFile->setFileName(filename);
                skyCoverageFile->read();
                degSq=skyCoverageFile->_data->totalDegSq();
            }
            
            unsigned int sNobs=0, sNdsc=0, sNtot=0;
            for (unsigned int row=0; row<data[fileID].size(); ++row) {
                sNobs += data[fileID][row].Nobs.getRef();
                sNdsc += data[fileID][row].Ndsc.getRef();
                sNtot += data[fileID][row].Ntot.getRef();
            }
            
            {
                unsigned int gslIndex=0;
                for (unsigned int k=0; k<_par->totalSize(); ++k) {
                    if (!_par->isFixed(k)) {
                        _par->set(k, gsl_vector_get(s->x,gslIndex));
                        ++gslIndex;
                    }
                }
            }
            
            double c = 1.0;
            //
            const unsigned int dof = _data->size() - _par->sizeNotFixed();
            const double chi = gsl_blas_dnrm2(s->f);
            if (dof > 0) {
                c = GSL_MAX_DBL(1.0, chi / sqrt(dof)); 
                // ORSA_DEBUG("chisq/dof = %g",  chi*chi/dof);
                // gmp_fprintf(fp,"chisq/dof = %g\n",  chi*chi/dof);
            } 
            //
            const double factor = c;
        
            gsl_matrix * covar = gsl_matrix_alloc(_par->sizeNotFixed(),_par->sizeNotFixed());
        
            gsl_multifit_covar(s->J, 0.0, covar);
        
            ORSA_DEBUG("writing file: [%s]",outputFile[fileID].c_str());
            FILE * fp = fopen(outputFile[fileID].c_str(),"w");
            
            // NOTE: all spaces in printed strings are measured to keep good alignment
            
            // first: var names
            fprintf(fp,"# ");
            fprintf(fp,"       JD ");
            fprintf(fp,"       year ");
            for (unsigned int p=0; p<_par->totalSize(); ++p) {
                const std::string::size_type pos = _par->name(p).find(fileID_str_vec[fileID]);
                if (pos != std::string::npos) {
                    // variable with correct fileID
                    fprintf(fp,"%11s +/- sigma ",_par->name(p).substr(0,pos).c_str());
                } else {
                    // check if it is a variable without any fileID
                    bool match=false;
                    for (unsigned int s=0; s<numFiles; ++s) {
                        if (_par->name(p).find(fileID_str_vec[s]) != std::string::npos) {
                            match=true;
                            break;
                        }
                    }
                    if (!match) {
                        fprintf(fp,"%11s +/- sigma ",_par->name(p).substr(0,pos).c_str());
                    }
                }
            }
            fprintf(fp," chisq/dof ");
            fprintf(fp,"   Nobs    Ndsc    Ntot ");
            fprintf(fp,"    degSq ");
            fprintf(fp,"   V0 ");
            fprintf(fp,"                   jobID ");
            fprintf(fp,"\n");

            // then: the values
            fprintf(fp,"%.3f ",orsaSolarSystem::timeToJulian(epoch));
            fprintf(fp,"%.6f ",orsaSolarSystem::fractionalYear(epoch));
            {
                unsigned int gslIndex=0;
                for (unsigned int p=0; p<_par->totalSize(); ++p) {
                    bool doPrint=false;
                    {
                        const std::string::size_type pos = _par->name(p).find(fileID_str_vec[fileID]);
                        if (pos != std::string::npos) {
                            // variable with correct fileID
                            doPrint=true;
                        } else {
                            // check if it is a variable without any fileID
                            bool match=false;
                            for (unsigned int s=0; s<numFiles; ++s) {
                                if (_par->name(p).find(fileID_str_vec[s]) != std::string::npos) {
                                    match=true;
                                    break;
                                }
                            }
                            if (!match) {
                                doPrint=true;
                            }
                        }
                    }
                    if (doPrint) {                        
                        // first, specific cases where unit conversions or factors are needed
                        if ((_par->name(p).find("U_limit") != std::string::npos) || 
                            (_par->name(p).find("w_U") != std::string::npos)) {
                            fprintf(fp,
                                    "%+.3e %+.3e ",
                                    orsa::FromUnits(_par->get(p)*orsa::radToArcsec(),orsa::Unit::HOUR),
                                    _par->isFixed(p)?0.0:orsa::FromUnits(factor*sqrt(gsl_matrix_get(covar,gslIndex,gslIndex))*orsa::radToArcsec(),orsa::Unit::HOUR));
                        } else if ( (_par->name(p).find("scale_GB")  != std::string::npos) ||
                                    (_par->name(p).find("scale_GL")  != std::string::npos) ) {
                            fprintf(fp,
                                    "%+.3e %+.3e ",
                                    orsa::radToDeg()*_par->get(p),
                                    _par->isFixed(p)?0.0:orsa::radToDeg()*factor*sqrt(gsl_matrix_get(covar,gslIndex,gslIndex)));
                        } else {
                            // generic one
                            fprintf(fp,
                                    "%+.3e %+.3e ",
                                    _par->get(p),
                                    _par->isFixed(p)?0.0:factor*sqrt(gsl_matrix_get(covar,gslIndex,gslIndex)));
                        }
                    }
                    if (!_par->isFixed(p)) {
                        ++gslIndex;
                    }
                } 
            }
            fprintf(fp,"%+.3e ",chi*chi/(dof>0?dof:1.0));
            fprintf(fp,"%7i %7i %7i ",sNobs,sNdsc,sNtot);
            fprintf(fp,"%.3e ",degSq);
            fprintf(fp,"%.2f ",V0);
            fprintf(fp,"%24s ",jobID[fileID].c_str());
            fprintf(fp,"\n");
        
            fclose(fp);
        
            gsl_matrix_free(covar);
        }        
    }
protected:
    void readOutputFiles() {
        
        const unsigned int numFiles = data.size();
        
        jobID.resize(numFiles);
        
        char varName[1024];
        
        for (unsigned int fileID=0; fileID<numFiles; ++fileID) {
            
            FitFileDataElement e;
            if (!readFitFile(e,outputFile[fileID].c_str())) {
                ORSA_DEBUG("problems...");
                return;
            }
            
            // V
            sprintf(varName,"V_limit_%06i",fileID); _par->set(varName,e.V_limit);
            sprintf(varName, "eta0_V_%06i",fileID); _par->set(varName,e.eta0_V);
            sprintf(varName,    "c_V_%06i",fileID); _par->set(varName,e.c_V);
            sprintf(varName,    "w_V_%06i",fileID); _par->set(varName,e.w_V);
            // U
            sprintf(varName,"U_limit_%06i",fileID); _par->set(varName,e.U_limit);
            sprintf(varName,    "w_U_%06i",fileID); _par->set(varName,e.w_U);
            // AM
            _par->set( "peak_AM",e.peak_AM);
            _par->set("scale_AM",e.scale_AM);
            _par->set("shape_AM",e.shape_AM);
            // GB
            _par->set( "drop_GB",e.drop_GB);
            _par->set("scale_GB",e.scale_GB);
            // GL
            _par->set("scale_GL",e.scale_GL);
            _par->set("shape_GL",e.shape_GL);
            // V0
            V0 = e.V0;
            // jobID
            jobID[fileID]=e.jobID;
        }           
    }    
protected:
    std::vector<DataStorage> data;
    double V0;
    std::vector<std::string> outputFile;
    std::vector<std::string> jobID;
    orsaInputOutput::MPCObsCodeFile * obsCodeFile;
    bool writeFile;
};

#endif // __ETA_H_

