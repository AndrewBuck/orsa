#include <orsa/paulMoment.h>

#include <orsa/box.h>
#include <orsa/double.h>
#include <orsa/legendre.h>
#include <orsa/multimin.h>
#include <orsa/statistic.h>
#include <orsa/unit.h>

#include <vector>

#include <gsl/gsl_rng.h>

#include <iostream>

using namespace orsa;

PaulMoment::PaulMoment(const unsigned int n) : osg::Referenced(true), order(n) {
  
    const unsigned int order_plus_one = order+1;
  
    {
        _M.resize(order_plus_one);
        _M_uncertainty.resize(order_plus_one);
        for (unsigned int i=0; i<order_plus_one; ++i) {
            _M[i].resize(order_plus_one-i);
            _M_uncertainty[i].resize(order_plus_one-i);
            for (unsigned int j=0; j<order_plus_one-i; ++j) {
                _M[i][j].resize(order_plus_one-i-j);
                _M_uncertainty[i][j].resize(order_plus_one-i-j);
            }
        }
    }
}

double PaulMoment::M (const int i,
                      const int j, 
                      const int k) const {
    if ( (i<0) || (j<0) || (k<0) || (i+j+k > (int)order) ) { 
        // ORSA_ERROR("index out of range [i:%i;j:%i;k:%i]",i,j,k);
        return 0; 
    } else {
        return _M[i][j][k];
    }
}

double PaulMoment::M_uncertainty (const int i,
                                  const int j, 
                                  const int k) const {
    if ( (i<0) || (j<0) || (k<0) || (i+j+k > (int)order) ) { 
        // ORSA_ERROR("index out of range [i:%i;j:%i;k:%i]",i,j,k);
        return 0;
    } else {
        return _M_uncertainty[i][j][k];
    }
}

void PaulMoment::setM (const double & val,
                       const int i, 
                       const int j, 
                       const int k) {
    if ( (i<0) || (j<0) || (k<0) || (i+j+k > (int)order) ) {
        ORSA_ERROR("index out of range [i:%i;j:%i;k:%i]",i,j,k);
        return;
    } else {
        _M[i][j][k] = val;
    }
}

void PaulMoment::setM_uncertainty (const double & val,
                                   const int i, 
                                   const int j, 
                                   const int k) {
    if ( (i<0) || (j<0) || (k<0) || (i+j+k > (int)order) ) { 
        ORSA_ERROR("index out of range [i:%i;j:%i;k:%i]",i,j,k);
        return;
    } else {
        _M_uncertainty[i][j][k] = val;
    }
}


/***/

/* 
   double normalization(const unsigned int l,
   const unsigned int m) {
   return sqrt(orsa::factorial(l+m).get_d() / ((2-orsa::kronecker(m,0))*(2*l+1)*orsa::factorial(l-m).get_d()));
   }
*/
//
static double normalization(const unsigned int l,
                            const unsigned int m) {
    // Cross checked with NEAR-Eros papers
    return sqrt( ((2-orsa::kronecker(m,0))*orsa::factorial(l-m).get_d()) / 
                 ((2*l+1)*orsa::factorial(l+m).get_d()) );
}

void orsa::convert(std::vector< std::vector<double> > & C,
                   std::vector< std::vector<double> > & S,
                   std::vector< std::vector<double> > & norm_C,
                   std::vector< std::vector<double> > & norm_S,
                   std::vector<double> & J,
                   const PaulMoment * const pm,
                   const double     & R0,
                   const bool         verbose) {
  
    if (verbose) {
        ORSA_DEBUG("R0: %f [km]",orsa::FromUnits(R0,orsa::Unit::KM,-1));
    }
  
    const unsigned int order = pm->order;
  
    // resize vectors
    C.resize(order+1);
    S.resize(order+1);
    norm_C.resize(order+1);
    norm_S.resize(order+1);
    //
    for (unsigned int l=0; l<=order; ++l) {
        C[l].resize(l+1);
        S[l].resize(l+1);
        norm_C[l].resize(l+1);
        norm_S[l].resize(l+1);
    }
    //
    J.resize(order+1);
  
    // C_lm coefficients
    //
    for (int l=0; l<=(int)order; ++l) {
        for (int m=0; m<=l; ++m) {
      
            double pq_factor=0;
            double pq_factor_uncertainty=0;
            //
            // integer division in limits
            for (int p=0;p<=(l/2);++p) {
                for (int q=0;q<=(m/2);++q) {
	  
                    double nu_factor=0;
                    double nu_factor_uncertainty=0;
                    //
                    for (int nu_x=0; nu_x<=p; ++nu_x) {
                        for (int nu_y=0; nu_y<=(p-nu_x); ++nu_y) {
	      
                            const int M_i = m-2*q+2*nu_x;
                            const int M_j = 2*q+2*nu_y;
                            const int M_k = l-m-2*nu_x-2*nu_y;
	      
                            if (M_i+M_j+M_k!=l) {
                                ORSA_DEBUG("WARNING!!! ***********************");
                            }
	      
                            if ( (M_i>=0) && 
                                 (M_j>=0) && 
                                 (M_k>=0) && 
                                 (M_i+M_j+M_k==l) ) {
		
                                // ORSA_DEBUG("requesting M[%i][%i][%i]...   l: %i",M_i, M_j, M_k, l);
		
                                const double nu_factor_base =
                                    (orsa::factorial(p).get_d() /(orsa::factorial(nu_x).get_d()*orsa::factorial(nu_y).get_d()*orsa::factorial(p-nu_x-nu_y).get_d()));
		
                                nu_factor += nu_factor_base * pm->M(M_i,M_j,M_k);
                                nu_factor_uncertainty += nu_factor_base * pm->M_uncertainty(M_i,M_j,M_k);
                            }
                        }
                    }
                    // need this because using M(i,j,k) instead of N(i,j,k)
                    nu_factor /= int_pow(R0,l);
                    nu_factor_uncertainty /= int_pow(R0,l);
	  
                    const double pq_factor_base = 
                        orsa::power_sign(p+q) *
                        orsa::binomial(l,p).get_d() *
                        orsa::binomial(2*l-2*p,l).get_d() *
                        orsa::binomial(m,2*q).get_d() *
                        orsa::pochhammer(l-m-2*p+1,m);
	  
                    pq_factor += 
                        pq_factor_base *
                        nu_factor;
	  
                    pq_factor_uncertainty += 
                        pq_factor_base *
                        nu_factor_uncertainty;
                }
            }
      
            pq_factor /= int_pow(2,l);
            pq_factor_uncertainty /= int_pow(2,l);
      
            const double C_lm = pq_factor;
            const double C_lm_uncertainty = fabs(pq_factor_uncertainty);
            //      
            const double norm_C_lm = C_lm*normalization(l,m);
            const double norm_C_lm_uncertainty = fabs(C_lm_uncertainty*normalization(l,m));
      
            if (verbose) {
                ORSA_DEBUG("     C[%i][%i] = %+16.12f +/- %16.12f",
                           l,m,     C_lm, C_lm_uncertainty);
                ORSA_DEBUG("norm_C[%i][%i] = %+16.12f +/- %16.12f   norm: %f",
                           l,m,norm_C_lm,norm_C_lm_uncertainty,normalization(l,m));
            }
      
            C[l][m]      = C_lm;
            norm_C[l][m] = norm_C_lm;
      
            // J_l is minus C_l0, where C_l0 is not normalized
            if (l>=2) {
                if (m==0) {
                    const double J_l = -C_lm;
                    const double J_l_uncertainty = -C_lm_uncertainty;
                    if (verbose) {
                        ORSA_DEBUG("J_%i = %+16.12f +/- %16.12f",l,J_l,J_l_uncertainty);
                    }
                    //
                    J[l] = J_l;
                }	
            } else {
                if (m==0) {
                    J[l] = 0.0;
                }
            }
      
        }
    }
  
    // S_lm coefficients
    //
    for (int l=0; l<=(int)order; ++l) {
        for (int m=1; m<=l; ++m) {
      
            double pq_factor=0;
            double pq_factor_uncertainty=0;
            //
            // integer division in limits
            for (int p=0;p<=(l/2);++p) {
                for (int q=0;q<=((m-1)/2);++q) {
	  
                    double nu_factor=0;
                    double nu_factor_uncertainty=0;
                    //
                    for (int nu_x=0; nu_x<=p; ++nu_x) {
                        for (int nu_y=0; nu_y<=(p-nu_x); ++nu_y) {
	      
                            const int M_i = m-2*q-1+2*nu_x;
                            const int M_j = 2*q+1+2*nu_y;
                            const int M_k = l-m-2*nu_x-2*nu_y;
	      
                            if (M_i+M_j+M_k!=l) {
                                ORSA_DEBUG("WARNING!!!");
                            }
	      
                            if ( (M_i>=0) && 
                                 (M_j>=0) && 
                                 (M_k>=0) && 
                                 (M_i+M_j+M_k==l) ) {
		
                                // ORSA_DEBUG("requesting M[%i][%i][%i]...   l: %i",M_i, M_j, M_k, l);
		
                                const double nu_factor_base =
                                    (orsa::factorial(p).get_d() /(orsa::factorial(nu_x).get_d()*orsa::factorial(nu_y).get_d()*orsa::factorial(p-nu_x-nu_y).get_d()));
		
                                nu_factor += nu_factor_base * pm->M(M_i,M_j,M_k);
                                nu_factor_uncertainty += nu_factor_base * pm->M_uncertainty(M_i,M_j,M_k);
                            }
                        }
                    }
                    // need this because using M(i,j,k) instead of N(i,j,k)
                    nu_factor /= int_pow(R0,l);
                    nu_factor_uncertainty /= int_pow(R0,l);
	  
                    const double pq_factor_base = 
                        orsa::power_sign(p+q) *
                        orsa::binomial(l,p).get_d() *
                        orsa::binomial(2*l-2*p,l).get_d() *
                        orsa::binomial(m,2*q+1).get_d() *
                        orsa::pochhammer(l-m-2*p+1,m);
	  
                    pq_factor += 
                        pq_factor_base *
                        nu_factor;
	  
                    pq_factor_uncertainty += 
                        pq_factor_base *
                        nu_factor_uncertainty;
                }
            }
            //
            pq_factor /= int_pow(2,l);
            pq_factor_uncertainty /= int_pow(2,l);
            //
            const double S_lm = pq_factor;
            const double S_lm_uncertainty = fabs(pq_factor_uncertainty);
            //      
            const double norm_S_lm = S_lm*normalization(l,m);
            const double norm_S_lm_uncertainty = fabs(S_lm_uncertainty*normalization(l,m));
      
            if (verbose) {
                ORSA_DEBUG("     S[%i][%i] = %+16.12f +/- %16.12f",
                           l,m,     S_lm, S_lm_uncertainty);
                ORSA_DEBUG("norm_S[%i][%i] = %+16.12f +/- %16.12f   norm: %f",
                           l,m,norm_S_lm,norm_S_lm_uncertainty,normalization(l,m));
            }
      
            S[l][m]      = S_lm;
            norm_S[l][m] = norm_S_lm;
        }
    }
  
}

/*** Use Multimin code to solve the inverse problem ***/

class PaulMomentsSolveMultimin : public orsa::Multimin {
public:
    PaulMomentsSolveMultimin(const unsigned int focusOrder_in, 
                             const std::vector< std::vector<double> > & norm_C_in,
                             const std::vector< std::vector<double> > & norm_S_in,
                             const double & R0_in) : 
        orsa::Multimin(), 
        focusOrder(focusOrder_in),
        norm_C(norm_C_in),
        norm_S(norm_S_in),
        R0(R0_in) { }
protected:
    const unsigned int focusOrder;
    const std::vector< std::vector<double> > & norm_C;
    const std::vector< std::vector<double> > & norm_S;
    const double R0;
public:
    double fun(const orsa::MultiminParameters * par) const {
        osg::ref_ptr<PaulMoment> local_pm = new PaulMoment(focusOrder);
        char parName[1024];
        for (unsigned int i=0; i<=focusOrder; ++i) {
            for (unsigned int j=0; j<=focusOrder; ++j) {
                for (unsigned int k=0; k<=focusOrder; ++k) {
                    if (i+j+k==focusOrder) {
                        sprintf(parName,"M_%i_%i_%i",i,j,k);
                        local_pm->setM(par->get(parName),i,j,k);
                    }
                }
            }
        }
    
        // new C,S values
        std::vector< std::vector<double> > local_C, local_S, local_norm_C, local_norm_S;
        std::vector<double> local_J;
        convert(local_C, local_S, local_norm_C, local_norm_S, local_J,
                local_pm.get(),
                R0);
    
        double retVal=0.0;
        {
            const unsigned int l=focusOrder;
            for (unsigned int m=0; m<=l; ++m) {
                retVal += orsa::square(local_norm_C[l][m]-norm_C[l][m]);
                if (m!=0) {
                    retVal += orsa::square(local_norm_S[l][m]-norm_S[l][m]);	  
                }
            }
        }
    
        return retVal;
    }
};

bool orsa::solve(PaulMoment * pm,
                 const std::vector< std::vector<double> > & norm_C,
                 const std::vector< std::vector<double> > & norm_S,
                 const double     & R0) {
  
    const unsigned int order = pm->order;
  
    for (unsigned int focusOrder=0; focusOrder<=order; ++focusOrder) {
    
        osg::ref_ptr<orsa::MultiminParameters> par = new orsa::MultiminParameters;
        //
        {
            char parName[1024];
            // for (unsigned int sum=0; sum<=order; ++sum) {
            for (unsigned int i=0; i<=focusOrder; ++i) {
                for (unsigned int j=0; j<=focusOrder; ++j) {
                    for (unsigned int k=0; k<=focusOrder; ++k) {
                        if (i+j+k==focusOrder) {
                            sprintf(parName,"M_%i_%i_%i",i,j,k);
                            // par->insert(parName,0.01*int_pow(R0,i+j+k),1e16*(int_pow(R0,i+j+k)));
                            par->insert(parName,0.01*int_pow(R0,focusOrder),int_pow(R0,focusOrder));
                        }
                    }
                }
            }
            // }
        }
    
        osg::ref_ptr<PaulMomentsSolveMultimin> multimin = 
            new PaulMomentsSolveMultimin(focusOrder,norm_C,norm_S,R0);
        //
        multimin->setMultiminParameters(par.get());
        //
        multimin->run_nmsimplex(1024*1024,1.0e-15*int_pow(R0,focusOrder));
        // multimin->run_conjugate_fr(1024*1024,1.0,1e-6,1.0e-15*int_pow(R0,focusOrder));
    
        // save results
        {
            char parName[1024];
            for (unsigned int i=0; i<=focusOrder; ++i) {
                for (unsigned int j=0; j<=focusOrder; ++j) {
                    for (unsigned int k=0; k<=focusOrder; ++k) {
                        if (i+j+k==focusOrder) { 
                            sprintf(parName,"M_%i_%i_%i",i,j,k);
                            pm->setM(par->get(parName),i,j,k);
                            pm->setM_uncertainty(0,i,j,k);
                        }
                    }
                }
            }
        }
    
    }
  
    return true;
}

static double EllipsoidExpansion_product_utility(const unsigned int n) {
    double product = 1;
    for (unsigned int k=1; k<=n; ++k) {
        product *= (2*k-1);
    }
    return product;
}

void orsa::EllipsoidExpansion(PaulMoment   * pm,
                              const double & a,
                              const double & b,
                              const double & c) {
  
    const unsigned int order = pm->order;
  
    for (unsigned int focusOrder=0; focusOrder<=order; ++focusOrder) {
        for (unsigned int i=0; i<=focusOrder; ++i) {
            for (unsigned int j=0; j<=focusOrder; ++j) {
                for (unsigned int k=0; k<=focusOrder; ++k) {
                    if (i+j+k==focusOrder) {
                        if (i%2==1) continue;
                        if (j%2==1) continue;
                        if (k%2==1) continue;
                        const double factor_i   = EllipsoidExpansion_product_utility(i/2);
                        const double factor_j   = EllipsoidExpansion_product_utility(j/2);
                        const double factor_k   = EllipsoidExpansion_product_utility(k/2);
                        const double factor_ijk = EllipsoidExpansion_product_utility((i+j+k)/2+2);
                        const double factor     = 3*(factor_i*factor_j*factor_k)/factor_ijk;
                        const double M_ijk      = factor*orsa::int_pow(a,i)*orsa::int_pow(b,j)*orsa::int_pow(c,k);
                        pm->setM(M_ijk,i,j,k);
                        ORSA_DEBUG("ijk: %i %i %i   factor: %g   M: %g",i,j,k,factor,M_ijk); 
                    }
                }
            }
        }
    }
  
}
