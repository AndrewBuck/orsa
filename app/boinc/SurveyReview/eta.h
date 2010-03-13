#ifndef __ETA_H_
#define __ETA_H_

#include <orsa/multimin.h>

class EfficiencyMultimin : public orsa::Multimin {
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
    data = data_in;
    V0 = V0_in;
    U0 = U0_in;
    osg::ref_ptr<orsa::MultiminParameters> par = new orsa::MultiminParameters;
#warning improve initial guess!
    //
    const double initial_U_limit = orsa::FromUnits(10.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1); // arcsec/hour
    const double initial_w_U     = orsa::FromUnits( 1.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1); // arcsec/hour
    //    
    // V pars
    par->insert("eta0_V",     0.95,  0.01);
    par->insert("c_V",        0.001, 0.00001);
    par->insert("V_limit",   19.00,  0.01);
    par->insert("w_V",        0.50,  0.01); 
    // U pars
    par->insert("eta0_U",     0.95,  0.01);
    par->insert("c_U",        1.0e-9,1.0e-12);
    par->insert("U_limit",    initial_U_limit, initial_U_limit*0.01);
    par->insert("w_U",        initial_w_U,     initial_w_U*0.01); 
    // ranges
    par->setRange("eta0_V",0.0,1.0);
    par->setRange("eta0_U",0.0,1.0);
    // par->setRangeMin("c",0.0); // don't set this one, as often c becomes negative while iterating
    par->setRangeMin("V_limit",V0);
    par->setRangeMin("U_limit",0.0);
    par->setRangeMax("U_limit",U0);
    //
    setMultiminParameters(par.get());
    //
    if (!run_nmsimplex(4096,0.001)) {
      // if (!run_conjugate_fr(1024,0.01,0.1,0.1)) {
      ORSA_WARNING("the fit did not converge.");
      return false;
    } else {
      return true;
    }
  }
 protected:
  double fun(const orsa::MultiminParameters * par) const {
    if (data.size()==0) return 0.0;
    double retVal=0.0;
    for (unsigned int k=0; k<data.size(); ++k) {
      const double eta_V = SkyCoverage::eta_V(data[k].V.getRef(),
					      par->get("V_limit"),
					      par->get("eta0_V"),
					      par->get("c_V"),
					      V0,
					      par->get("w_V"));
      const double eta_U = SkyCoverage::eta_U(data[k].U.getRef(),
					      par->get("U_limit"),
					      par->get("eta0_U"),
					      par->get("c_U"),
					      U0,
					      par->get("w_U"));
      retVal += orsa::square((data[k].eta.getRef()-(eta_V*eta_U))/data[k].sigmaEta.getRef());
      
      /* ORSA_DEBUG("data.eta: %g   computed: %g   (V: %g, U: %g)",
	 data[k].eta.getRef(),
	 eta_V*eta_U,
	 eta_V,
	 eta_U);      
      */
    }
    
    ORSA_DEBUG("eta_V: %g   eta_U: %g   V_limit: %g   U_limit: %g [arcsec/hour]  c_V: %g   c_U: %g   w_V: %g   w_U: %g [arcsec/hour]   retVal: %g   retVal/Npoints: %g",
	       par->get("eta0_V"),
	       par->get("eta0_U"),
	       par->get("V_limit"),
	       orsa::FromUnits(par->get("U_limit")*orsa::radToArcsec(),orsa::Unit::HOUR),
	       par->get("c_V"),
	       par->get("c_U"),
	       par->get("w_V"),
	       orsa::FromUnits(par->get("w_U")*orsa::radToArcsec(),orsa::Unit::HOUR),
	       retVal,
	       retVal/data.size());
    return retVal;
  }
 protected:
  DataStorage data;
  double V0;
  double U0;
};

#endif // __ETA_H_
