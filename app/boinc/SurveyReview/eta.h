#ifndef __ETA_H_
#define __ETA_H_

#include <orsa/multifit.h>

class EfficiencyMultifit : public orsa::Multifit {
 public:
  class DataElement {
  public:
    // V = apparent magnitude
    // U = apparent velocity
    orsa::Cache<double> V, U, eta, sigmaEta;
  };
  typedef std::vector<DataElement> DataStorage;
 public:
  // V = apparent magnitude
  // U = apparent velocity
  bool fit(const DataStorage & data_in,
	   const double      & V0_in,
	   const double      & U0_in) {
    
    // local copies
    data = data_in;
    V0 = V0_in;
    U0 = U0_in;
    
    // basic checks
    if (data.size()==0) return false;
    
    osg::ref_ptr<orsa::MultifitParameters> par = new orsa::MultifitParameters;
#warning improve initial guess!
    //
    const double initial_U_limit = orsa::FromUnits(10.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1); // arcsec/hour
    const double initial_w_U     = orsa::FromUnits( 1.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1); // arcsec/hour
    //    
    // V pars
    par->insert("eta0_V",     0.95,  0.001);
    par->insert("c_V",        0.001, 0.00001);
    par->insert("V_limit",   19.00,  0.001);
    par->insert("w_V",        0.50,  0.001); 
    // U pars
    // par->insert("eta0_U",     0.95,  0.001);
    // par->insert("c_U",        1.0e-9,1.0e-12);
    par->insert("U_limit",    initial_U_limit, initial_U_limit*0.001);
    par->insert("w_U",        initial_w_U,     initial_w_U*0.001); 
    // ranges
    // par->setRange("eta0_V",0.0,1.0);
    // excluding eta0_U for the moment...
    // par->setRange("eta0_U",0.0,1.0);
    // par->setRangeMin("c",0.0); // don't set this one, as often c becomes negative while iterating
    // par->setRangeMin("V_limit",V0);
    // par->setRangeMin("U_limit",0.0);
    // par->setRangeMax("U_limit",U0);
    //
    setMultifitParameters(par.get());
    
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
    (*(localPar.get())) = (*par);
    
    // modify
    localPar->set(p,par->get(p)+d*par->getDelta(p));
    
    const double eta_V = SkyCoverage::eta_V(data->getD("V",row),
					    localPar->get("V_limit"),
					    localPar->get("eta0_V"),
					    localPar->get("c_V"),
					    V0,
					    localPar->get("w_V"));
    const double eta_U = SkyCoverage::eta_U(data->getD("U",row),
					    localPar->get("U_limit"),
					    1.0, // localPar->get("eta0_U"),
					    0.0, // localPar->get("c_U"),
					    U0,
					    localPar->get("w_U"));
    
    return (eta_V*eta_U);
  }
 protected: 
  void singleIterationDone(const gsl_multifit_fdfsolver * s) const { 
    
    for (unsigned int k=0; k<_par->size(); ++k) {
      _par->set(k, gsl_vector_get(s->x,k));
    }
    
    double c = 1.0;
    //
    const unsigned int dof = _data->size() - _par->size();
    if (dof > 0) {
      const double chi = gsl_blas_dnrm2(s->f);
      c = GSL_MAX_DBL(1.0, chi / sqrt(dof)); 
      ORSA_DEBUG("chisq/dof = %g",  chi*chi/dof);
      // gmp_fprintf(fp,"chisq/dof = %g\n",  chi*chi/dof);
    }
    //
    const double factor = c;
    
    gsl_matrix * covar = gsl_matrix_alloc(_par->size(),_par->size());
    
    gsl_multifit_covar(s->J, 0.0, covar);
    
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
    
    for (unsigned int p=0; p<_par->size(); ++p) {
      // first, specific cases where unit conversions or factors are needed
      if ((_par->name(p) == "U_limit") || 
	  (_par->name(p) == "w_U")) {
	ORSA_DEBUG("%s: %g +/- %g [arcsec/hour]",
		   _par->name(p).c_str(),
		   orsa::FromUnits(_par->get(p)*orsa::radToArcsec(),orsa::Unit::HOUR),
		   orsa::FromUnits(factor*ERR(p)*orsa::radToArcsec(),orsa::Unit::HOUR));
      } else {
	// generic one
	ORSA_DEBUG("%s: %g +/- %g",
		   _par->name(p).c_str(),
		   _par->get(p),
		   factor*ERR(p));
      }
    }
    
    gsl_matrix_free(covar);  
  }
 protected:
  void success(const gsl_multifit_fdfsolver * s) const { 
#warning finish here!
    // implement this similarly to singleIterationDone above to save final parameters and corresponding uncertainties
  }
 protected:
  DataStorage data;
  double V0;
  double U0;
};

#endif // __ETA_H_

