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
  orsa::Cache<double> lunarPhase;
  orsa::Cache<double> minAirMass;
  orsa::Cache<bool>   observed;
  orsa::Cache<bool>   discovered;
};

class EfficiencyMultifit : public orsa::Multifit {
 public:
  class DataElement {
  public:
    // V = apparent magnitude
    // U = apparent velocity
    orsa::Cache<double> V, U, eta, sigmaEta;
    orsa::Cache<unsigned int> Nobs, Ndsc, Ntot;
  };
  typedef std::vector<DataElement> DataStorage;
 public:
  // V = apparent magnitude
  // U = apparent velocity
  bool fit(orsa::MultifitParameters * par,
	   const DataStorage & data_in,
	   const double      & V0_in,
	   const std::string & outputFile_in,
	   const std::string & jobID_in,
	   orsaInputOutput::MPCObsCodeFile * obsCodeFile_in) {
    
    // local copies
    data = data_in;
    V0 = V0_in;
    //
    outputFile  = outputFile_in;
    jobID       = jobID_in;
    obsCodeFile = obsCodeFile_in;
    
    // basic checks
    if (data.size()==0) return false;
    
    /* 
       osg::ref_ptr<orsa::MultifitParameters> par = new orsa::MultifitParameters;
       //
       const double initial_U_limit = orsa::FromUnits( 10.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1); // arcsec/hour
       const double initial_w_U     = orsa::FromUnits(  1.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1); // arcsec/hour
       //    
       // V pars
       par->insert("eta0_V",  0.90, 0.0001);
       par->insert("c_V",     0.01, 0.00001);
       par->insert("V_limit",19.00, 0.001);
       par->insert("w_V",     0.50, 0.0001); 
       // U pars
       par->insert("U_limit",    initial_U_limit, 0.001*initial_U_limit);
       par->insert("w_U",        initial_w_U,     0.001*initial_w_U); 
       // mixing 
       par->insert("beta",  0.0, 0.00001);
       // ranges
       par->setRange("eta0_V",0.0,1.0);
    */
    //
    setMultifitParameters(par);
    
    osg::ref_ptr<orsa::MultifitData> fitData = new orsa::MultifitData;
    //
    fitData->insertVariable("V");
    fitData->insertVariable("U");
    //
    for (unsigned int row=0; row<data.size(); ++row) {
      fitData->insertD("V",row,data[row].V.getRef());
      fitData->insertD("U",row,data[row].U.getRef());
      fitData->insertF(row,data[row].eta.getRef());
      fitData->insertSigma(row,data[row].sigmaEta.getRef());
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
    localPar->set(p,par->get(p)+d*par->getDelta(p));
    
    const double eta = SkyCoverage::eta(data->getD("V",row),
					localPar->get("V_limit"),
					localPar->get("eta0_V"),
					V0,
					localPar->get("c_V"),
					localPar->get("w_V"),
					data->getD("U",row),
					localPar->get("U_limit"),
					localPar->get("w_U"),
					localPar->get("beta"));
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
    
    // call abort() if an uncertainty gets too large
    {
      unsigned int gslIndex=0;
      for (unsigned int p=0; p<_par->totalSize(); ++p) {
	if (!_par->isFixed(p)) {
	  if  (factor*sqrt(gsl_matrix_get(covar,gslIndex,gslIndex)) > fabs(1.0e9*_par->get(p))) {
	    ORSA_DEBUG("uncertainty of parameter [%s] is too large, aborting",
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
	if ((_par->name(p) == "U_limit") || 
	    (_par->name(p) == "w_U")) {
	  ORSA_DEBUG("%s: %g +/- %g [arcsec/hour]",
		     _par->name(p).c_str(),
		     orsa::FromUnits(_par->get(p)*orsa::radToArcsec(),orsa::Unit::HOUR),
		     _par->isFixed(p)?0.0:orsa::FromUnits(factor*sqrt(gsl_matrix_get(covar,gslIndex,gslIndex))*orsa::radToArcsec(),orsa::Unit::HOUR));
	} else if (_par->name(p) == "beta") {
	  ORSA_DEBUG("%s: %g +/- %g [deg]",
		     _par->name(p).c_str(),
		     orsa::radToDeg()*_par->get(p),
		     _par->isFixed(p)?0.0:orsa::radToDeg()*factor*sqrt(gsl_matrix_get(covar,gslIndex,gslIndex)));
	} else {
	  // generic one
	  ORSA_DEBUG("%s: %g +/- %g",
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
  }
 protected:
  void success(const gsl_multifit_fdfsolver * s) const {
    
    std::string obsCode;
    orsa::Time epoch;
    int year;
    int dayOfYear;
    if (!SkyCoverage::processFilename(jobID.c_str(),
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
      sprintf(filename,"%s.txt",jobID.c_str());
      skyCoverageFile->setFileName(filename);
      skyCoverageFile->read();
      degSq=skyCoverageFile->_data->totalDegSq();
    }
    
    unsigned int sNobs=0, sNdsc=0, sNtot=0;
    for (unsigned int row=0; row<data.size(); ++row) {
      sNobs += data[row].Nobs.getRef();
      sNdsc += data[row].Ndsc.getRef();
      sNtot += data[row].Ntot.getRef();
    }
    
    {
      unsigned int gslIndex=0;
      for (unsigned int k=0; k<_par->totalSize(); ++k) {
	if (!_par->isFixed(k)) {
	  _par->set(k, gsl_vector_get(s->x,k));
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
    
    ORSA_DEBUG("writing file: [%s]",outputFile.c_str());
    FILE * fp = fopen(outputFile.c_str(),"w");
    
    // NOTE: all spaces in printed strings are measured to keep good alignment
    
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
	if ((_par->name(p) == "U_limit") || 
	    (_par->name(p) == "w_U")) {
	  fprintf(fp,
		  "%+.3e %+.3e ",
		  orsa::FromUnits(_par->get(p)*orsa::radToArcsec(),orsa::Unit::HOUR),
		  _par->isFixed(p)?0.0:orsa::FromUnits(factor*sqrt(gsl_matrix_get(covar,gslIndex,gslIndex))*orsa::radToArcsec(),orsa::Unit::HOUR));
	} else if (_par->name(p) == "beta") {
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
    fprintf(fp,"%s ",jobID.c_str());
    fprintf(fp,"\n");
    
    fclose(fp);
    
    gsl_matrix_free(covar);  
 }
 protected:
  DataStorage data;
  double V0;
  std::string outputFile;
  std::string jobID;
  orsaInputOutput::MPCObsCodeFile * obsCodeFile;
};

#endif // __ETA_H_

