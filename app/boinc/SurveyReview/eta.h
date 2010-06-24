#ifndef __ETA_H_
#define __ETA_H_

#include <orsa/multifit.h>

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
    orsa::Cache<double> activeTime;
    orsa::Cache<bool>   epochFromField;
    orsa::Cache<bool>   observed;
    orsa::Cache<bool>   discovered;
};

class EfficiencyMultifit : public orsa::Multifit {
public:
    class DataElement {
    public:
        orsa::Cache<double> V, U, AM, eta, sigmaEta;
        orsa::Cache<double> SE, LE, LP, GL, GB, AZ, LA, SA;
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
                                            localPar->get("scale_GB"));
        
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
                if ((_par->name(p).find("U_limit") != std::string::npos) || 
                    (_par->name(p).find("w_U") != std::string::npos)) {
                    ORSA_DEBUG("%20s: %g +/- %g [arcsec/hour]",
                               _par->name(p).c_str(),
                               orsa::FromUnits(_par->get(p)*orsa::radToArcsec(),orsa::Unit::HOUR),
                               _par->isFixed(p)?0.0:orsa::FromUnits(factor*sqrt(gsl_matrix_get(covar,gslIndex,gslIndex))*orsa::radToArcsec(),orsa::Unit::HOUR));
                } else if (_par->name(p).find("scale_GB") != std::string::npos) {
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
        
        // if (writeFile) writeOutputFiles(s);
    }
protected:
    void runCompleted(const bool /* success */, const gsl_multifit_fdfsolver * s) const {
        if (writeFile) writeOutputFiles(s);
    }
protected:
    void writeOutputFiles(const gsl_multifit_fdfsolver * s) const {
        
        const unsigned int numFiles = data.size();
        
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
                exit(0);
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
            
#warning need to re-align var names with values
            
            // first var names
            fprintf(fp,"# ");
            fprintf(fp,"       JD ");
            fprintf(fp,"       year ");
            for (unsigned int p=0; p<_par->totalSize(); ++p) {
                fprintf(fp,"%11s +/- sigma ",_par->name(p).c_str());
            }
            fprintf(fp," chisq/dof ");
            fprintf(fp," Nobs  Ndsc  Ntot ");
            fprintf(fp,"    degSq ");
            fprintf(fp,"   V0 ");
            fprintf(fp,"      jobID ");
            fprintf(fp,"\n");
        
            fprintf(fp,"%.3f ",orsaSolarSystem::timeToJulian(epoch));
            fprintf(fp,"%.6f ",orsaSolarSystem::fractionalYear(epoch));
            {
                unsigned int gslIndex=0;
                for (unsigned int p=0; p<_par->totalSize(); ++p) {
                    // first, specific cases where unit conversions or factors are needed
                    if ((_par->name(p).find("U_limit") != std::string::npos) || 
                        (_par->name(p).find("w_U") != std::string::npos)) {
                        fprintf(fp,
                                "%+.3e %+.3e ",
                                orsa::FromUnits(_par->get(p)*orsa::radToArcsec(),orsa::Unit::HOUR),
                                _par->isFixed(p)?0.0:orsa::FromUnits(factor*sqrt(gsl_matrix_get(covar,gslIndex,gslIndex))*orsa::radToArcsec(),orsa::Unit::HOUR));
                    } else if (_par->name(p).find("scale_GB") != std::string::npos) {
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
                    if (!_par->isFixed(p)) {
                        ++gslIndex;
                    }
                } 
            }
            fprintf(fp,"%+.3e ",chi*chi/(dof>0?dof:1.0));
            fprintf(fp,"%5i %5i %5i ",sNobs,sNdsc,sNtot);
            fprintf(fp,"%.3e ",degSq);
            fprintf(fp,"%.2f ",V0);
            fprintf(fp,"%s ",jobID[fileID].c_str());
            fprintf(fp,"\n");
        
            fclose(fp);
        
            gsl_matrix_free(covar);
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

