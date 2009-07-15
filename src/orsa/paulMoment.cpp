#include <orsa/paulMoment.h>

#include <orsa/box.h>
#include <orsa/double.h>
#include <orsa/legendre.h>
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

const double PaulMoment::M (const int i,
			    const int j, 
			    const int k) const {
  if ( (i<0) || (j<0) || (k<0) || (i+j+k > (int)order) ) { ORSA_ERROR("index out of range"); return 0; }
  return _M[i][j][k];
}

const double PaulMoment::M_uncertainty (const int i,
					const int j, 
					const int k) const {
  if ( (i<0) || (j<0) || (k<0) || (i+j+k > (int)order) ) { ORSA_ERROR("index out of range"); return 0; }
  return _M_uncertainty[i][j][k];
}

void PaulMoment::setM (const double & val,
		       const int i, 
		       const int j, 
		       const int k) {
  if ( (i<0) || (j<0) || (k<0) || (i+j+k > (int)order) ) { ORSA_ERROR("index out of range"); return; }
  _M[i][j][k] = val;
}

void PaulMoment::setM_uncertainty (const double & val,
				   const int i, 
				   const int j, 
				   const int k) {
  if ( (i<0) || (j<0) || (k<0) || (i+j+k > (int)order) ) { ORSA_ERROR("index out of range"); return; }
  _M_uncertainty[i][j][k] = val;
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

void orsa::convert(const PaulMoment * const pm,
		   const double     & R0) {
  
  ORSA_DEBUG("R0: %f [km]",orsa::FromUnits(R0,orsa::Unit::KM,-1));
  
  const unsigned int order = pm->order;
  
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
      
      ORSA_DEBUG("     C[%i][%i] = %+16.12f +/- %16.12f",
		 l,m,     C_lm, C_lm_uncertainty);
      ORSA_DEBUG("norm_C[%i][%i] = %+16.12f +/- %16.12f   norm: %f",
		 l,m,norm_C_lm,norm_C_lm_uncertainty,normalization(l,m));
      
      if ((l>=2) && (m==0)) {
	// J_l is minus C_l0, not-normalized
	const double J_l = -C_lm;
	const double J_l_uncertainty = -C_lm_uncertainty;
	ORSA_DEBUG("J_%i = %+16.12f +/- %16.12f",l,J_l,J_l_uncertainty);
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
      
      ORSA_DEBUG("     S[%i][%i] = %+16.12f +/- %16.12f",
		 l,m,     S_lm, S_lm_uncertainty);
      ORSA_DEBUG("norm_S[%i][%i] = %+16.12f +/- %16.12f   norm: %f",
		 l,m,norm_S_lm,norm_S_lm_uncertainty,normalization(l,m));
    }
  }
  
}
