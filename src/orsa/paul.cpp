#include <orsa/paul.h>

#include <orsa/paulMoment.h>
#include <orsa/util.h>

#include <algorithm>
#include <vector>

using namespace orsa;

orsa::Double Paul::C_lmn(const int l,
			 const int m,
			 const int n) {
  if ( (l==0) &&
       (m==0) &&
       (n==0) ) { 
    ORSA_ERROR("singular C value...");
    return zero();
  }
  
  const orsa::Double retVal = 
    orsa::one() / 
    orsa::Double(3.0 - 
		 orsa::kronecker(l,0) - 
		 orsa::kronecker(m,0) - 
		 orsa::kronecker(n,0) );
  
  /* 
     ORSA_DEBUG("C(%i,%i,%i) = %Fg",
     l,m,n,
     retVal.get_mpf_t());
  */
  
  return retVal;
}

//

Paul::t_lmnLMN * Paul::t_lmnLMN::_instance = 0;

Paul::t_lmnLMN * Paul::t_lmnLMN::instance() {
  if (_instance == 0) {
    _instance = new Paul::t_lmnLMN;
  }
  return _instance;
}

Paul::t_lmnLMN::t_lmnLMN() {
  
}

Paul::t_lmnLMN::~t_lmnLMN() { 
  _instance = 0;
}

orsa::Double Paul::t_lmnLMN::get(const int l,
				 const int m,
				 const int n,
				 const int L,
				 const int M,
				 const int N) const {
  
  /* 
     ORSA_DEBUG("called get(%i,%i,%i,%i,%i,%i)",
     l,m,n,L,M,N);
  */
  
  if ( (l==0) &&
       (m==0) &&
       (n==0) ) {
    if ( (L==0) &&
	 (M==0) &&
	 (N==0) ) {
      return orsa::one();
    } else {
      return orsa::zero();
    }	
  }
  
  if ( ((l+L)%2) ||
       ((m+M)%2) ||
       ((n+N)%2) ) {
    return zero();
  }
  
  if ( (l<0) ||
       (m<0) ||
       (n<0) ||
       (L<0) ||
       (M<0) ||
       (N<0) ) {
    return zero();
  }
  
  if ( (L>l) || 
       (M>m) || 
       (N>n) ) {
    return zero();
  }
  
  /* 
     int ls,ms,ns;
     sort(ls,ms,ns,l,m,n);
     
     int Ls,Ms,Ns;
     sort(Ls,Ms,Ns,L,M,N);
     
     return trueGet(ls,ms,ns,Ls,Ms,Ns);
  */
  
  return trueGet(l,m,n,L,M,N);
}

orsa::Double Paul::t_lmnLMN::trueGet(const int l,
				     const int m,
				     const int n,
				     const int L,
				     const int M,
				     const int N) const {
  
  /* 
     ORSA_DEBUG("called trueGet(%i,%i,%i,%i,%i,%i)",
     l,m,n,L,M,N);
  */
  
  /* 
     if ( (l==0) &&
     (m==0) &&
     (n==0) &&
     (L==0) &&
     (M==0) &&
     (N==0) ) { 
     // ORSA_DEBUG("--MARK-- 0,0,0");
     return orsa::one();
     }
     
     if ( (l==0) &&
     (m==0) &&
     (n==0) ) {
     return orsa::zero();
     }	
  */
  
  /* 
     if ( (l==0) &&
     (m==0) &&
     (n==0) ) {
     if ( (L==0) &&
     (M==0) &&
     (N==0) ) {
     return orsa::one();
     } else {
     return orsa::zero();
     }	
     }
  */
  
  /* 
     if ( ((l+L)%2) ||
     ((m+M)%2) ||
     ((n+N)%2) ) {
     return zero();
     }
  */
  
  /* 
     if ( (l<0) ||
     (m<0) ||
     (n<0) ||
     (L<0) ||
     (M<0) ||
     (N<0) ) {
     return zero();
     }
  */
  
  /* 
     if ( (L>l) || 
     (M>m) || 
     (N>n) ) {
     return zero();
     }   
  */
  
  // prepare the data container
  resize(l+m+n);
  
  if (0) {
    ORSA_DEBUG("s: %i",
	       _data.size());
    ORSA_DEBUG("s: %i",
	       _data[l].size());
    ORSA_DEBUG("s: %i",
	       _data[l][m].size());
    ORSA_DEBUG("s: %i",
	       _data[l][m][n].size());
    ORSA_DEBUG("s: %i",
	       _data[l][m][n][L].size());
    ORSA_DEBUG("s: %i",
	       _data[l][m][n][L][M].size());
  }
  
  // size checks
  /* 
     if ((int)_data.size() < (1+l)) {
     _data.resize(1+l);
     }	
     for (unsigned int _l=0; _l<_data.size(); ++_l) {
     if ((int)_data[_l].size() < (1+m)) {
     _data[_l].resize(1+m);
     }
     for (unsigned int _m=0; _m<_data[_l].size(); ++_m) {
     if ((int)_data[_l][_m].size() < (1+n)) {
     _data[_l][_m].resize(1+n);
     }
     for (unsigned int _n=0; _n<_data[_l][_m].size(); ++_n) {
     if ((int)_data[_l][_m][_n].size() < (1+L)) {
     _data[_l][_m][_n].resize(1+L);
     }
     for (unsigned int _L=0; _L<_data[_l][_m][_n].size(); ++_L) {
     if ((int)_data[_l][_m][_n][_L].size() < (1+M)) {
     _data[_l][_m][_n][_L].resize(1+M);
     }
     for (unsigned int _M=0; _M<_data[_l][_m][_n][_L].size(); ++_M) {
     if ((int)_data[_l][_m][_n][_L][_M].size() < (1+N)) {
     _data[_l][_m][_n][_L][_M].resize(1+N);
     }
     }
     }
     }
     }
     }
  */
  //  
  /* 
     if ((int)_data.size() < (1+l)) {
     _data.resize(1+l);
     }	
     if ((int)_data[l].size() < (1+m)) {
     _data[l].resize(1+m);
     }
     if ((int)_data[l][m].size() < (1+n)) {
     _data[l][m].resize(1+n);
     }
     if ((int)_data[l][m][n].size() < (1+L)) {
     _data[l][m][n].resize(1+L);
     }
     if ((int)_data[l][m][n][L].size() < (1+M)) {
     _data[l][m][n][L].resize(1+M);
     }
     if ((int)_data[l][m][n][L][M].size() < (1+N)) {
     _data[l][m][n][L][M].resize(1+N);
     }
  */
  
  // general rule
  if (!(_data[l][m][n][L][M][N].isSet())) {
    // ORSA_DEBUG("computing...");
    const orsa::Double _C_lmn = C_lmn(l,m,n);
    _data[l][m][n][L][M][N] =
      (orsa::kronecker(l,0)-1)*( (2*l-_C_lmn)*get(l-1,m,n,L-1,M,N) +
				 (l-1)*(l-_C_lmn)*get(l-2,m,n,L,M,N) ) +
      (orsa::kronecker(m,0)-1)*( (2*m-_C_lmn)*get(l,m-1,n,L,M-1,N) +
				 (m-1)*(m-_C_lmn)*get(l,m-2,n,L,M,N) ) +
      (orsa::kronecker(n,0)-1)*( (2*n-_C_lmn)*get(l,m,n-1,L,M,N-1) +
				 (n-1)*(n-_C_lmn)*get(l,m,n-2,L,M,N) );      
  }
  //
  /* 
     ORSA_DEBUG("get(%i,%i,%i,%i,%i,%i) = %Fg",
     l,m,n,L,M,N,
     _data[l][m][n][L][M][N].getRef().get_mpf_t());
  */
  //
  return _data[l][m][n][L][M][N].getRef();
}

void Paul::t_lmnLMN::resize(const size_t order) const {
  static size_t oldOrder = 0;
  // ORSA_DEBUG("order: %i   oldOrder: %i",order,oldOrder);
  if (order > oldOrder) {
    oldOrder = order;
    const size_t orderPlusOne = 1+order;
    // ORSA_DEBUG("--MARK--");
    _data.resize(std::max(_data.size(),orderPlusOne));
    for (unsigned int _l=0; _l<_data.size(); ++_l) {
      // ORSA_DEBUG("--MARK--");
      _data[_l].resize(std::max(_data[_l].size(),orderPlusOne-_l));
      for (unsigned int _m=0; _m<_data[_l].size(); ++_m) {
	// ORSA_DEBUG("--MARK--");
	_data[_l][_m].resize(std::max(_data[_l][_m].size(),orderPlusOne-_l-_m));
	for (unsigned int _n=0; _n<_data[_l][_m].size(); ++_n) {
	  // ORSA_DEBUG("--MARK--");
	  _data[_l][_m][_n].resize(std::max(_data[_l][_m][_n].size(),orderPlusOne));
	  for (unsigned int _L=0; _L<_data[_l][_m][_n].size(); ++_L) {
	    // ORSA_DEBUG("--MARK--");
	    _data[_l][_m][_n][_L].resize(std::max(_data[_l][_m][_n][_L].size(),orderPlusOne-_L));
	    for (unsigned int _M=0; _M<_data[_l][_m][_n][_L].size(); ++_M) {
	      // ORSA_DEBUG("--MARK--");
	      _data[_l][_m][_n][_L][_M].resize(std::max(_data[_l][_m][_n][_L][_M].size(),orderPlusOne-_L-_M));
	    }
	  }
	}
      }
    }
  }
}

void Paul::t_lmnLMN::sort(int & l_out,
			  int & m_out,
			  int & n_out,
			  const int l_in,
			  const int m_in,
			  const int n_in) const {
  // this is not efficient, should use some kind of sort3() algorithm...
  std::vector<int> v;
  v.push_back(l_in);
  v.push_back(m_in);
  v.push_back(n_in);
  std::sort(v.begin(),v.end());
  l_out = v[0];
  m_out = v[1];
  n_out = v[2];
  
  /* 
     ORSA_DEBUG("in: (%i,%i,%i)   out: (%i,%i,%i)",
     l_in,m_in,n_in,
     l_out,m_out,n_out);
  */
}

orsa::Double Paul::gravitationalPotential(const orsa::PaulMoment * M1,
					  const orsa::Matrix     & A1_g2l,
					  const orsa::PaulMoment * M2,
					  const orsa::Matrix     & A2_g2l,
					  const orsa::Vector     & R) {
  
  // const orsa::Double oneOverR = one()/R.length();
  //  
  IntPowCache oneOverR_PC(one()/R.length());
  
  const orsa::Matrix & R1 = A1_g2l;
  const orsa::Matrix & R2 = A2_g2l;
  //
  const orsa::Double l11 = R1.getM11();
  const orsa::Double l21 = R1.getM12();
  const orsa::Double l31 = R1.getM13();
  //
  const orsa::Double m11 = R1.getM21();
  const orsa::Double m21 = R1.getM22();
  const orsa::Double m31 = R1.getM23();
  //
  const orsa::Double n11 = R1.getM31();
  const orsa::Double n21 = R1.getM32();
  const orsa::Double n31 = R1.getM33();
  //
  const orsa::Double l12 = R2.getM11();
  const orsa::Double l22 = R2.getM12();
  const orsa::Double l32 = R2.getM13();
  //
  const orsa::Double m12 = R2.getM21();
  const orsa::Double m22 = R2.getM22();
  const orsa::Double m32 = R2.getM23();
  //
  const orsa::Double n12 = R2.getM31();
  const orsa::Double n22 = R2.getM32();
  const orsa::Double n32 = R2.getM33();
  //
  const orsa::Double lx = l11*l12+l21*l22+l31*l32;
  const orsa::Double mx = m11*l12+m21*l22+m31*l32;
  const orsa::Double nx = n11*l12+n21*l22+n31*l32;
  //
  const orsa::Double ly = l11*m12+l21*m22+l31*m32;
  const orsa::Double my = m11*m12+m21*m22+m31*m32;
  const orsa::Double ny = n11*m12+n21*m22+n31*m32;
  //
  const orsa::Double lz = l11*n12+l21*n22+l31*n32;
  const orsa::Double mz = m11*n12+m21*n22+m31*n32;
  const orsa::Double nz = n11*n12+n21*n22+n31*n32;
  //
  IntPowCache lx_PC(lx);
  IntPowCache mx_PC(mx);
  IntPowCache nx_PC(nx);
  IntPowCache ly_PC(ly);
  IntPowCache my_PC(my);
  IntPowCache ny_PC(ny);
  IntPowCache lz_PC(lz);
  IntPowCache mz_PC(mz);
  IntPowCache nz_PC(nz);
  //
  if (0) {
    ORSA_DEBUG("lx: %Fg",lx.get_mpf_t());
    ORSA_DEBUG("ly: %Fg",ly.get_mpf_t());
    ORSA_DEBUG("lz: %Fg",lz.get_mpf_t());
    ORSA_DEBUG("mx: %Fg",mx.get_mpf_t());
    ORSA_DEBUG("my: %Fg",my.get_mpf_t());
    ORSA_DEBUG("mz: %Fg",mz.get_mpf_t());
    ORSA_DEBUG("nx: %Fg",nx.get_mpf_t());
    ORSA_DEBUG("ny: %Fg",ny.get_mpf_t());
    ORSA_DEBUG("nz: %Fg",nz.get_mpf_t());
  }
  
  const orsa::Double & csi  = R.getX();
  const orsa::Double & eta  = R.getY();
  const orsa::Double & zeta = R.getZ();
  //
  /* 
     const orsa::Double csi1  = l11*csi + l21*eta + l31*zeta;
     const orsa::Double eta1  = m11*csi + m21*eta + m31*zeta;
     const orsa::Double zeta1 = n11*csi + n21*eta + n31*zeta;
  */
  //
  IntPowCache csi1_PC( l11*csi + l21*eta + l31*zeta);
  IntPowCache eta1_PC( m11*csi + m21*eta + m31*zeta);
  IntPowCache zeta1_PC(n11*csi + n21*eta + n31*zeta);
  
  orsa::Double outerSum = orsa::zero();
  
  for (unsigned int i1=0; i1<=M1->order(); ++i1) {
    for (unsigned int j1=0; j1<=M1->order()-i1; ++j1) {
      for (unsigned int k1=0; k1<=M1->order()-i1-j1; ++k1) {
	
	if (M1->M(i1,j1,k1) == zero()) continue;
	
	for (unsigned int i2=0; i2<=M2->order(); ++i2) {
	  for (unsigned int j2=0; j2<=M2->order()-i2; ++j2) {
	    for (unsigned int k2=0; k2<=M2->order()-i2-j2; ++k2) {
	      
	      if (M2->M(i2,j2,k2) == zero()) continue;
	      
	      orsa::Double innerSum = orsa::zero();
	      
	      for (unsigned int i3=0; i3<=i2; ++i3) {
		for (unsigned int j3=0; j3<=j2; ++j3) {
		  for (unsigned int k3=0; k3<=k2; ++k3) {
		    
		    for (unsigned int i4=0; i4<=i3; ++i4) {
		      for (unsigned int j4=0; j4<=j3; ++j4) {
			for (unsigned int k4=0; k4<=k3; ++k4) {
			  
			  const unsigned int i5 = i1 + i4 + j4 + k4;
			  const unsigned int j5 = j1 + i3 - i4 + j3 - j4 + k3 - k4;
			  const unsigned int k5 = k1 + i2 - i3 + j2 - j3 + k2 - k3;
			  
			  /* 
			     ORSA_DEBUG("i5: %i",i5);
			     ORSA_DEBUG("j5: %i",j5);
			     ORSA_DEBUG("k5: %i",k5);
			  */
			  
			  orsa::Double fiveSum = orsa::zero();
			  
			  for (unsigned int L=0; L<=i5; ++L) {
			    for (unsigned int M=0; M<=j5; ++M) {
			      for (unsigned int N=0; N<=k5; ++N) {
				fiveSum +=
				  orsa::Paul::t_lmnLMN::instance()->get(i5,j5,k5,L,M,N) *
				  // orsa::int_pow(csi1, L) * 
				  csi1_PC.get(L) *
				  // orsa::int_pow(eta1, M) * 
				  eta1_PC.get(M) *
				  // orsa::int_pow(zeta1,N) * 
				  zeta1_PC.get(N) *
				  // orsa::int_pow(oneOverR,L+M+N);
				  oneOverR_PC.get(L+M+N);
			      }
			    }
			  }
			  
			  innerSum +=
			    fiveSum *
			    orsa::binomial(i2,i3) * orsa::binomial(i3,i4) *
			    orsa::binomial(j2,j3) * orsa::binomial(j3,j4) *
			    orsa::binomial(k2,k3) * orsa::binomial(k3,k4) *
			    // orsa::int_pow(lx,i4) * orsa::int_pow(mx,i3-i4) * orsa::int_pow(nx,i2-i3) *  
			    lx_PC.get(i4) * mx_PC.get(i3-i4) * nx_PC.get(i2-i3) *
			    // orsa::int_pow(ly,j4) * orsa::int_pow(my,j3-j4) * orsa::int_pow(ny,j2-j3) *  
			    ly_PC.get(j4) * my_PC.get(j3-j4) * ny_PC.get(j2-j3) *
			    // orsa::int_pow(lz,k4) * orsa::int_pow(mz,k3-k4) * orsa::int_pow(nz,k2-k3);  
			    lz_PC.get(k4) * mz_PC.get(k3-k4) * nz_PC.get(k2-k3);
			}
		      }
		    }
		    
		  }
		}
	      }
	      
	      outerSum += 
		innerSum * 
		orsa::power_sign(i1+j1+k1) *
		M1->M(i1,j1,k1) * 
		M2->M(i2,j2,k2) *
		// orsa::int_pow(oneOverR,i1+j1+k1+i2+j2+k2+1) / 
		oneOverR_PC.get(i1+j1+k1+i2+j2+k2+1) / 
		orsa::Double( factorial(i1) * factorial(j1) * factorial(k1) *
			      factorial(i2) * factorial(j2) * factorial(k2) );
	      
	      if (0) {
		// debug
		
		const orsa::Double term = 
		  innerSum * 
		  orsa::power_sign(i1+j1+k1) *
		  M1->M(i1,j1,k1) * 
		  M2->M(i2,j2,k2) *
		  // orsa::int_pow(oneOverR,i1+j1+k1+i2+j2+k2+1) / 
		  oneOverR_PC.get(i1+j1+k1+i2+j2+k2+1) / 
		  orsa::Double( factorial(i1) * factorial(j1) * factorial(k1) *
				factorial(i2) * factorial(j2) * factorial(k2) );
		
		if (fabs(term*1.0e6) > fabs(outerSum)) {
		  ORSA_DEBUG("[%02i][%02i][%02i][%02i][%02i][%02i]   outerSum: %14.6Fe   term: %14.6Fe",
			     i1,j1,k1,i2,j2,k2,
			     outerSum.get_mpf_t(),
			     term.get_mpf_t());
		}
		
	      }
	      
	    }
	  }
	}
	
      }
    }
  }
  
  const orsa::Double U = outerSum;
  
  /* 
     ORSA_DEBUG("Paul U: %.20Fe",
     U.get_mpf_t());
  */
  
  return U;
}


orsa::Vector Paul::gravitationalForce(const orsa::PaulMoment * M1,
				      const orsa::Matrix     & A1_g2l,
				      const orsa::PaulMoment * M2,
				      const orsa::Matrix     & A2_g2l,
				      const orsa::Vector     & R) {
  
  /* 
     ORSA_DEBUG("M1: %x",M1);
     ORSA_DEBUG("M2: %x",M2);
     ORSA_DEBUG("A1: %x",A1);
     ORSA_DEBUG("A2: %x",A2);
  */
  
  // const orsa::Double oneOverR = one()/R.length();
  //  
  IntPowCache oneOverR_PC(one()/R.length());
  
  const orsa::Matrix & R1 = A1_g2l;
  const orsa::Matrix & R2 = A2_g2l;
  //
  const orsa::Double l11 = R1.getM11();
  const orsa::Double l21 = R1.getM12();
  const orsa::Double l31 = R1.getM13();
  //
  const orsa::Double m11 = R1.getM21();
  const orsa::Double m21 = R1.getM22();
  const orsa::Double m31 = R1.getM23();
  //
  const orsa::Double n11 = R1.getM31();
  const orsa::Double n21 = R1.getM32();
  const orsa::Double n31 = R1.getM33();
  //
  const orsa::Double l12 = R2.getM11();
  const orsa::Double l22 = R2.getM12();
  const orsa::Double l32 = R2.getM13();
  //
  const orsa::Double m12 = R2.getM21();
  const orsa::Double m22 = R2.getM22();
  const orsa::Double m32 = R2.getM23();
  //
  const orsa::Double n12 = R2.getM31();
  const orsa::Double n22 = R2.getM32();
  const orsa::Double n32 = R2.getM33();
  //
  const orsa::Double lx = l11*l12+l21*l22+l31*l32;
  const orsa::Double mx = m11*l12+m21*l22+m31*l32;
  const orsa::Double nx = n11*l12+n21*l22+n31*l32;
  //
  const orsa::Double ly = l11*m12+l21*m22+l31*m32;
  const orsa::Double my = m11*m12+m21*m22+m31*m32;
  const orsa::Double ny = n11*m12+n21*m22+n31*m32;
  //
  const orsa::Double lz = l11*n12+l21*n22+l31*n32;
  const orsa::Double mz = m11*n12+m21*n22+m31*n32;
  const orsa::Double nz = n11*n12+n21*n22+n31*n32;
  //
  IntPowCache lx_PC(lx);
  IntPowCache mx_PC(mx);
  IntPowCache nx_PC(nx);
  IntPowCache ly_PC(ly);
  IntPowCache my_PC(my);
  IntPowCache ny_PC(ny);
  IntPowCache lz_PC(lz);
  IntPowCache mz_PC(mz);
  IntPowCache nz_PC(nz);
  //
  if (0) {
    ORSA_DEBUG("lx: %Fg",lx.get_mpf_t());
    ORSA_DEBUG("ly: %Fg",ly.get_mpf_t());
    ORSA_DEBUG("lz: %Fg",lz.get_mpf_t());
    ORSA_DEBUG("mx: %Fg",mx.get_mpf_t());
    ORSA_DEBUG("my: %Fg",my.get_mpf_t());
    ORSA_DEBUG("mz: %Fg",mz.get_mpf_t());
    ORSA_DEBUG("nx: %Fg",nx.get_mpf_t());
    ORSA_DEBUG("ny: %Fg",ny.get_mpf_t());
    ORSA_DEBUG("nz: %Fg",nz.get_mpf_t());
  }
  
  const orsa::Double & csi  = R.getX();
  const orsa::Double & eta  = R.getY();
  const orsa::Double & zeta = R.getZ();
  //
  const orsa::Double csi1  = l11*csi + l21*eta + l31*zeta;
  const orsa::Double eta1  = m11*csi + m21*eta + m31*zeta;
  const orsa::Double zeta1 = n11*csi + n21*eta + n31*zeta;
  //
  IntPowCache csi1_PC(csi1);
  IntPowCache eta1_PC(eta1);
  IntPowCache zeta1_PC(zeta1);
  
  // orsa::Double outerSum = orsa::zero();
  //
  orsa::Double Fx = orsa::zero();
  orsa::Double Fy = orsa::zero();
  orsa::Double Fz = orsa::zero();
  
  for (unsigned int i1=0; i1<=M1->order(); ++i1) {
    for (unsigned int j1=0; j1<=M1->order()-i1; ++j1) {
      for (unsigned int k1=0; k1<=M1->order()-i1-j1; ++k1) {
	
	if (M1->M(i1,j1,k1) == zero()) continue;
	
	for (unsigned int i2=0; i2<=M2->order(); ++i2) {
	  for (unsigned int j2=0; j2<=M2->order()-i2; ++j2) {
	    for (unsigned int k2=0; k2<=M2->order()-i2-j2; ++k2) {
	      
	      if (M2->M(i2,j2,k2) == zero()) continue;
	      
	      /* 
		 ORSA_DEBUG("i1: %i   j1: %i   k1: %i   i2: %i   j2: %i   k2: %i",
		 i1,j1,k1,i2,j2,k2);
	      */
	      
	      // orsa::Double innerSum = orsa::zero();
	      //
	      orsa::Double innerFx = orsa::zero();
	      orsa::Double innerFy = orsa::zero();
	      orsa::Double innerFz = orsa::zero();
	      
	      for (unsigned int i3=0; i3<=i2; ++i3) {
		for (unsigned int j3=0; j3<=j2; ++j3) {
		  for (unsigned int k3=0; k3<=k2; ++k3) {
		    
		    for (unsigned int i4=0; i4<=i3; ++i4) {
		      for (unsigned int j4=0; j4<=j3; ++j4) {
			for (unsigned int k4=0; k4<=k3; ++k4) {
			  
			  /* 
			     ORSA_DEBUG("i3: %i   j3: %i   k3: %i   i4: %i   j4: %i   k4: %i",
			     i3,j3,k3,i4,j4,k4);
			  */
			  
			  const unsigned int i5 = i1 + i4 + j4 + k4;
			  const unsigned int j5 = j1 + i3 - i4 + j3 - j4 + k3 - k4;
			  const unsigned int k5 = k1 + i2 - i3 + j2 - j3 + k2 - k3;
			  
			  /* 
			     ORSA_DEBUG("i5: %i",i5);
			     ORSA_DEBUG("j5: %i",j5);
			     ORSA_DEBUG("k5: %i",k5);
			  */
			  
			  /* 
			     ORSA_DEBUG("i5: %i   j5: %i   k5: %i",
			     i5,j5,k5);
			  */
			  
			  // orsa::Double fiveSum = orsa::zero();
			  //
			  orsa::Double fiveSumFx = orsa::zero();
			  orsa::Double fiveSumFy = orsa::zero();
			  orsa::Double fiveSumFz = orsa::zero();
			  
		          for (unsigned int L=0; L<=i5; ++L) {
			    for (unsigned int M=0; M<=j5; ++M) {
			      for (unsigned int N=0; N<=k5; ++N) {
				
				const orsa::Double local_t = orsa::Paul::t_lmnLMN::instance()->get(i5,j5,k5,L,M,N);
				
				if (local_t == zero()) {
				  continue;
				}
				
				const orsa::Double tR = local_t * oneOverR_PC.get(L+M+N);
				
				if (L != 0) {
				  fiveSumFx +=
				    tR *
				    L * csi1_PC.get(L-1) *
				    eta1_PC.get(M) *
				    zeta1_PC.get(N);				  
				}
				//
				if (M != 0) {
				  fiveSumFy +=
				    tR *
				    csi1_PC.get(L) *
				    M * eta1_PC.get(M-1) *
				    zeta1_PC.get(N);
				}
				//
				if (N != 0) {
				  fiveSumFz +=
				    tR *
				    csi1_PC.get(L) *
				    eta1_PC.get(M) *
				    N * zeta1_PC.get(N-1);
				}
			      }
			    }
			  }
			  
			  for (unsigned int L=0; L<=i5; ++L) {
			    for (unsigned int M=0; M<=j5; ++M) {
			      for (unsigned int N=0; N<=k5; ++N) {
				
				const orsa::Double local_t = orsa::Paul::t_lmnLMN::instance()->get(i5,j5,k5,L,M,N);
				
				if (local_t == zero()) {
				  continue;
				}
				
				const orsa::Double ntR = (i5+j5+k5+L+M+N+1) * local_t * oneOverR_PC.get(L+M+N+2);
				
				fiveSumFx -=
				  ntR *
				  csi1_PC.get(L+1) *
				  eta1_PC.get(M) *
				  zeta1_PC.get(N);
				//
				fiveSumFy -=
				  ntR *
				  csi1_PC.get(L) *
				  eta1_PC.get(M+1) *
				  zeta1_PC.get(N);
				//
				fiveSumFz -=
				  ntR *
				  csi1_PC.get(L) *
				  eta1_PC.get(M) *
				  zeta1_PC.get(N+1);
			      }
			    }
			  }
			  
			  const orsa::Double commonInnerFactor = 
			    orsa::binomial(i2,i3) * orsa::binomial(i3,i4) *
			    orsa::binomial(j2,j3) * orsa::binomial(j3,j4) *
			    orsa::binomial(k2,k3) * orsa::binomial(k3,k4) *
			    lx_PC.get(i4) * mx_PC.get(i3-i4) * nx_PC.get(i2-i3) *
			    ly_PC.get(j4) * my_PC.get(j3-j4) * ny_PC.get(j2-j3) *
			    lz_PC.get(k4) * mz_PC.get(k3-k4) * nz_PC.get(k2-k3);
			  
			  /* 
			     ORSA_DEBUG("commonInnerFactor: %+12.6Ff",
			     commonInnerFactor.get_mpf_t());
			  */
			  
			  innerFx += fiveSumFx * commonInnerFactor;
			  innerFy += fiveSumFy * commonInnerFactor;
			  innerFz += fiveSumFz * commonInnerFactor;
			  
			}
		      }
		    }
		    
		  }
		}
	      }
	      
	      /* 
		 outerSum += 
		 innerSum * 
		 orsa::power_sign(i1+j1+k1) *
		 M1->M(i1,j1,k1) * 
		 M2->M(i2,j2,k2) *
		 oneOverR_PC.get(i1+j1+k1+i2+j2+k2+1) / 
		 orsa::Double( factorial(i1) * factorial(j1) * factorial(k1) *
		 factorial(i2) * factorial(j2) * factorial(k2) );
	      */
	      
	      const orsa::Double commonFactor = 
		orsa::power_sign(i1+j1+k1) *
		M1->M(i1,j1,k1) * 
		M2->M(i2,j2,k2) *
		oneOverR_PC.get(i1+j1+k1+i2+j2+k2+1) / 
		orsa::Double( factorial(i1) * factorial(j1) * factorial(k1) *
			      factorial(i2) * factorial(j2) * factorial(k2) );
	      
	      /* 
		 ORSA_DEBUG("commonFactor: %+18.12Fe",
		 commonFactor.get_mpf_t());
	      */
	      
	      Fx += innerFx * commonFactor;
	      Fy += innerFy * commonFactor;
	      Fz += innerFz * commonFactor;
	      
	      if (0) {
		ORSA_DEBUG("[%02i][%02i][%02i][%02i][%02i][%02i]   F: %20.12Fe %20.12Fe %20.12Fe",
			   i1,j1,k1,i2,j2,k2,
			   Fx.get_mpf_t(),
			   Fy.get_mpf_t(),
			   Fz.get_mpf_t());
	      }
	      
	    }
	  }
	}
	
      }
    }
  }
  
  /* 
     ORSA_DEBUG("Paul Fx: %.20Fe",Fx.get_mpf_t());
     ORSA_DEBUG("Paul Fy: %.20Fe",Fy.get_mpf_t());
     ORSA_DEBUG("Paul Fz: %.20Fe",Fz.get_mpf_t());
  */
  
  // rotate back, capital X,Y,Z
  const orsa::Double FX = l11*Fx + m11*Fy + n11*Fz;
  const orsa::Double FY = l21*Fx + m21*Fy + n21*Fz;
  const orsa::Double FZ = l31*Fx + m31*Fy + n31*Fz;
  
  const orsa::Vector F(-FX,-FY,-FZ);
  
  /* 
     ORSA_DEBUG("Paul F.X: %.20Fe",F.getX().get_mpf_t());
     ORSA_DEBUG("Paul F.Y: %.20Fe",F.getY().get_mpf_t());
     ORSA_DEBUG("Paul F.Z: %.20Fe",F.getZ().get_mpf_t());
  */
  //
  // ORSA_DEBUG("Paul Force: %.20Fe",F.length().get_mpf_t());
  
  return F;
}

/* 
orsa::Vector Paul::gravitationalTorque(const orsa::PaulMoment * M1,
				       const orsa::Attitude   * A1,
				       const orsa::PaulMoment * M2,
				       const orsa::Attitude   * A2,
				       const orsa::Vector     & R,
				       const orsa::Time       & t) {
  
  // const orsa::Double oneOverR = one()/R.length();
  //  
  IntPowCache oneOverR_PC(one()/R.length());
  
  const orsa::Matrix R1 = A1->globalToLocal(t);
  const orsa::Matrix R2 = A2->globalToLocal(t);
  //
  const orsa::Double l11 = R1.getM11();
  const orsa::Double l21 = R1.getM12();
  const orsa::Double l31 = R1.getM13();
  //
  const orsa::Double m11 = R1.getM21();
  const orsa::Double m21 = R1.getM22();
  const orsa::Double m31 = R1.getM23();
  //
  const orsa::Double n11 = R1.getM31();
  const orsa::Double n21 = R1.getM32();
  const orsa::Double n31 = R1.getM33();
  //
  const orsa::Double l12 = R2.getM11();
  const orsa::Double l22 = R2.getM12();
  const orsa::Double l32 = R2.getM13();
  //
  const orsa::Double m12 = R2.getM21();
  const orsa::Double m22 = R2.getM22();
  const orsa::Double m32 = R2.getM23();
  //
  const orsa::Double n12 = R2.getM31();
  const orsa::Double n22 = R2.getM32();
  const orsa::Double n32 = R2.getM33();
  //
  const orsa::Double lx = l11*l12+l21*l22+l31*l32;
  const orsa::Double mx = m11*l12+m21*l22+m31*l32;
  const orsa::Double nx = n11*l12+n21*l22+n31*l32;
  //
  const orsa::Double ly = l11*m12+l21*m22+l31*m32;
  const orsa::Double my = m11*m12+m21*m22+m31*m32;
  const orsa::Double ny = n11*m12+n21*m22+n31*m32;
  //
  const orsa::Double lz = l11*n12+l21*n22+l31*n32;
  const orsa::Double mz = m11*n12+m21*n22+m31*n32;
  const orsa::Double nz = n11*n12+n21*n22+n31*n32;
  //
  IntPowCache lx_PC(lx);
  IntPowCache mx_PC(mx);
  IntPowCache nx_PC(nx);
  IntPowCache ly_PC(ly);
  IntPowCache my_PC(my);
  IntPowCache ny_PC(ny);
  IntPowCache lz_PC(lz);
  IntPowCache mz_PC(mz);
  IntPowCache nz_PC(nz);
  //
  if (1) {
    ORSA_DEBUG("lx: %Fg",lx.get_mpf_t());
    ORSA_DEBUG("ly: %Fg",ly.get_mpf_t());
    ORSA_DEBUG("lz: %Fg",lz.get_mpf_t());
    ORSA_DEBUG("mx: %Fg",mx.get_mpf_t());
    ORSA_DEBUG("my: %Fg",my.get_mpf_t());
    ORSA_DEBUG("mz: %Fg",mz.get_mpf_t());
    ORSA_DEBUG("nx: %Fg",nx.get_mpf_t());
    ORSA_DEBUG("ny: %Fg",ny.get_mpf_t());
    ORSA_DEBUG("nz: %Fg",nz.get_mpf_t());
  }
  
  const orsa::Double & csi  = R.getX();
  const orsa::Double & eta  = R.getY();
  const orsa::Double & zeta = R.getZ();
  //
  const orsa::Double csi1  = l11*csi + l21*eta + l31*zeta;
  const orsa::Double eta1  = m11*csi + m21*eta + m31*zeta;
  const orsa::Double zeta1 = n11*csi + n21*eta + n31*zeta;
  //
  IntPowCache csi1_PC(csi1);
  IntPowCache eta1_PC(eta1);
  IntPowCache zeta1_PC(zeta1);
  
  orsa::Double Fx = orsa::zero();
  orsa::Double Fy = orsa::zero();
  orsa::Double Fz = orsa::zero();
  
  for (unsigned int i1=0; i1<=M1->order(); ++i1) {
    for (unsigned int j1=0; j1<=M1->order()-i1; ++j1) {
      for (unsigned int k1=0; k1<=M1->order()-i1-j1; ++k1) {
	
	if (M1->M(i1,j1,k1) == zero()) continue;
	
	for (unsigned int i2=0; i2<=M2->order(); ++i2) {
	  for (unsigned int j2=0; j2<=M2->order()-i2; ++j2) {
	    for (unsigned int k2=0; k2<=M2->order()-i2-j2; ++k2) {
	      
	      if (M2->M(i2,j2,k2) == zero()) continue;
	      
	      orsa::Double innerFx = orsa::zero();
	      orsa::Double innerFy = orsa::zero();
	      orsa::Double innerFz = orsa::zero();
	      
	      for (unsigned int i3=0; i3<=i2; ++i3) {
		for (unsigned int j3=0; j3<=j2; ++j3) {
		  for (unsigned int k3=0; k3<=k2; ++k3) {
		    
		    for (unsigned int i4=0; i4<=i3; ++i4) {
		      for (unsigned int j4=0; j4<=j3; ++j4) {
			for (unsigned int k4=0; k4<=k3; ++k4) {
			  
			  const unsigned int i5 = i1 + i4 + j4 + k4;
			  const unsigned int j5 = j1 + i3 - i4 + j3 - j4 + k3 - k4;
			  const unsigned int k5 = k1 + i2 - i3 + j2 - j3 + k2 - k3;
			  
			  orsa::Double fiveSumFx = orsa::zero();
			  orsa::Double fiveSumFy = orsa::zero();
			  orsa::Double fiveSumFz = orsa::zero();
			  
		          for (unsigned int L=0; L<=i5; ++L) {
			    for (unsigned int M=0; M<=j5; ++M) {
			      for (unsigned int N=0; N<=k5; ++N) {
				if (L != 0) {
				  fiveSumFx +=
				    orsa::Paul::t_lmnLMN::instance()->get(i5,j5,k5,L,M,N) *
				    L * csi1_PC.get(L-1) *
				    eta1_PC.get(M) *
				    zeta1_PC.get(N) *
				    oneOverR_PC.get(L+M+N);
				}
				//
				if (M != 0) {
				  fiveSumFy +=
				    orsa::Paul::t_lmnLMN::instance()->get(i5,j5,k5,L,M,N) *
				    csi1_PC.get(L) *
				    M * eta1_PC.get(M-1) *
				    zeta1_PC.get(N) *
				    oneOverR_PC.get(L+M+N);
				}
				//
				if (N != 0) {
				  fiveSumFz +=
				    orsa::Paul::t_lmnLMN::instance()->get(i5,j5,k5,L,M,N) *
				    csi1_PC.get(L) *
				    eta1_PC.get(M) *
				    N * zeta1_PC.get(N-1) *
				    oneOverR_PC.get(L+M+N);
				}
			      }
			    }
			  }
			  
		          for (unsigned int L=0; L<=i5; ++L) {
			    for (unsigned int M=0; M<=j5; ++M) {
			      for (unsigned int N=0; N<=k5; ++N) {
				fiveSumFx -=
				  (i5+j5+k5+L+M+N+1) *
				  orsa::Paul::t_lmnLMN::instance()->get(i5,j5,k5,L,M,N) *
				  csi1_PC.get(L+1) *
				  eta1_PC.get(M) *
				  zeta1_PC.get(N) *
				  oneOverR_PC.get(L+M+N+2);
				//
				fiveSumFy -=
				  (i5+j5+k5+L+M+N+1) *
				  orsa::Paul::t_lmnLMN::instance()->get(i5,j5,k5,L,M,N) *
				  csi1_PC.get(L) *
				  eta1_PC.get(M+1) *
				  zeta1_PC.get(N) *
				  oneOverR_PC.get(L+M+N+2);
				//
				fiveSumFz -=
				  (i5+j5+k5+L+M+N+1) *
				  orsa::Paul::t_lmnLMN::instance()->get(i5,j5,k5,L,M,N) *
				  csi1_PC.get(L) *
				  eta1_PC.get(M) *
				  zeta1_PC.get(N+1) *
				  oneOverR_PC.get(L+M+N+2);
			      }
			    }
			  }
			  
			  const orsa::Double commonInnerFactor = 
			    orsa::binomial(i2,i3) * orsa::binomial(i3,i4) *
			    orsa::binomial(j2,j3) * orsa::binomial(j3,j4) *
			    orsa::binomial(k2,k3) * orsa::binomial(k3,k4) *
			    lx_PC.get(i4) * mx_PC.get(i3-i4) * nx_PC.get(i2-i3) *
			    ly_PC.get(j4) * my_PC.get(j3-j4) * ny_PC.get(j2-j3) *
			    lz_PC.get(k4) * mz_PC.get(k3-k4) * nz_PC.get(k2-k3);
			  
			  innerFx += fiveSumFx * commonInnerFactor;
			  innerFy += fiveSumFy * commonInnerFactor;
			  innerFz += fiveSumFz * commonInnerFactor;
			  
			}
		      }
		    }
		    
		  }
		}
	      }
	      
	      const orsa::Double commonFactor = 
		orsa::power_sign(i1+j1+k1) *
		M1->M(i1,j1,k1) * 
		M2->M(i2,j2,k2) *
		oneOverR_PC.get(i1+j1+k1+i2+j2+k2+1) / 
		orsa::Double( factorial(i1) * factorial(j1) * factorial(k1) *
			      factorial(i2) * factorial(j2) * factorial(k2) );
	      
	      Fx += innerFx * commonFactor;
	      Fy += innerFy * commonFactor;
	      Fz += innerFz * commonFactor;
	      
	      if (1) {
		ORSA_DEBUG("[%02i][%02i][%02i][%02i][%02i][%02i]   F: %20.12Fe %20.12Fe %20.12Fe",
			   i1,j1,k1,i2,j2,k2,
			   Fx.get_mpf_t(),
			   Fy.get_mpf_t(),
			   Fz.get_mpf_t());
	      }
	      
	    }
	  }
	}
	
      }
    }
  }
  
  ORSA_DEBUG("Paul Fx: %.20Fe",Fx.get_mpf_t());
  ORSA_DEBUG("Paul Fy: %.20Fe",Fy.get_mpf_t());
  ORSA_DEBUG("Paul Fz: %.20Fe",Fz.get_mpf_t());
  
  // rotate back, capital X,Y,Z
  const orsa::Double FX = l11*Fx + m11*Fy + n11*Fz;
  const orsa::Double FY = l21*Fx + m21*Fy + n21*Fz;
  const orsa::Double FZ = l31*Fx + m31*Fy + n31*Fz;
  
  const orsa::Vector F(-FX,-FY,-FZ);
  
  ORSA_DEBUG("Paul F.X: %.20Fe",F.getX().get_mpf_t());
  ORSA_DEBUG("Paul F.Y: %.20Fe",F.getY().get_mpf_t());
  ORSA_DEBUG("Paul F.Z: %.20Fe",F.getZ().get_mpf_t());
  
  return F;
}
*/

orsa::Vector Paul::gravitationalTorque(const orsa::PaulMoment * M1,
				       const orsa::Matrix     & A1_g2l,
				       const orsa::PaulMoment * M2,
				       const orsa::Matrix     & A2_g2l,
				       const orsa::Vector     & R) {
  
  // const orsa::Double oneOverR = one()/R.length();
  //  
  IntPowCache oneOverR_PC(one()/R.length());
  
  const orsa::Matrix & R1 = A1_g2l;
  const orsa::Matrix & R2 = A2_g2l;
  //
  const orsa::Double l11 = R1.getM11();
  const orsa::Double l21 = R1.getM12();
  const orsa::Double l31 = R1.getM13();
  //
  const orsa::Double m11 = R1.getM21();
  const orsa::Double m21 = R1.getM22();
  const orsa::Double m31 = R1.getM23();
  //
  const orsa::Double n11 = R1.getM31();
  const orsa::Double n21 = R1.getM32();
  const orsa::Double n31 = R1.getM33();
  //
  const orsa::Double l12 = R2.getM11();
  const orsa::Double l22 = R2.getM12();
  const orsa::Double l32 = R2.getM13();
  //
  const orsa::Double m12 = R2.getM21();
  const orsa::Double m22 = R2.getM22();
  const orsa::Double m32 = R2.getM23();
  //
  const orsa::Double n12 = R2.getM31();
  const orsa::Double n22 = R2.getM32();
  const orsa::Double n32 = R2.getM33();
  //
  const orsa::Double lx = l11*l12+l21*l22+l31*l32;
  const orsa::Double mx = m11*l12+m21*l22+m31*l32;
  const orsa::Double nx = n11*l12+n21*l22+n31*l32;
  //
  const orsa::Double ly = l11*m12+l21*m22+l31*m32;
  const orsa::Double my = m11*m12+m21*m22+m31*m32;
  const orsa::Double ny = n11*m12+n21*m22+n31*m32;
  //
  const orsa::Double lz = l11*n12+l21*n22+l31*n32;
  const orsa::Double mz = m11*n12+m21*n22+m31*n32;
  const orsa::Double nz = n11*n12+n21*n22+n31*n32;
  //
  IntPowCache lx_PC(lx);
  IntPowCache mx_PC(mx);
  IntPowCache nx_PC(nx);
  IntPowCache ly_PC(ly);
  IntPowCache my_PC(my);
  IntPowCache ny_PC(ny);
  IntPowCache lz_PC(lz);
  IntPowCache mz_PC(mz);
  IntPowCache nz_PC(nz);
  //
  if (0) {
    ORSA_DEBUG("lx: %Fg",lx.get_mpf_t());
    ORSA_DEBUG("ly: %Fg",ly.get_mpf_t());
    ORSA_DEBUG("lz: %Fg",lz.get_mpf_t());
    ORSA_DEBUG("mx: %Fg",mx.get_mpf_t());
    ORSA_DEBUG("my: %Fg",my.get_mpf_t());
    ORSA_DEBUG("mz: %Fg",mz.get_mpf_t());
    ORSA_DEBUG("nx: %Fg",nx.get_mpf_t());
    ORSA_DEBUG("ny: %Fg",ny.get_mpf_t());
    ORSA_DEBUG("nz: %Fg",nz.get_mpf_t());
  }
  
  const orsa::Double & csi  = R.getX();
  const orsa::Double & eta  = R.getY();
  const orsa::Double & zeta = R.getZ();
  //
  const orsa::Double csi1  = l11*csi + l21*eta + l31*zeta;
  const orsa::Double eta1  = m11*csi + m21*eta + m31*zeta;
  const orsa::Double zeta1 = n11*csi + n21*eta + n31*zeta;
  //
  IntPowCache csi1_PC(csi1);
  IntPowCache eta1_PC(eta1);
  IntPowCache zeta1_PC(zeta1);
  
  orsa::Double Tx = orsa::zero();
  orsa::Double Ty = orsa::zero();
  orsa::Double Tz = orsa::zero();
  
  for (unsigned int i1=0; i1<=M1->order(); ++i1) {
    for (unsigned int j1=0; j1<=M1->order()-i1; ++j1) {
      for (unsigned int k1=0; k1<=M1->order()-i1-j1; ++k1) {
	
	if ( (M1->M(i1+1,j1,k1) == zero()) &&
	     (M1->M(i1,j1+1,k1) == zero()) &&
	     (M1->M(i1,j1,k1+1) == zero()) ) continue;
	
	for (unsigned int i2=0; i2<=M2->order(); ++i2) {
	  for (unsigned int j2=0; j2<=M2->order()-i2; ++j2) {
	    for (unsigned int k2=0; k2<=M2->order()-i2-j2; ++k2) {
	      
	      if (M2->M(i2,j2,k2) == zero()) continue;
	      
	      orsa::Double innerFx = orsa::zero();
	      orsa::Double innerFy = orsa::zero();
	      orsa::Double innerFz = orsa::zero();
	      
	      for (unsigned int i3=0; i3<=i2; ++i3) {
		for (unsigned int j3=0; j3<=j2; ++j3) {
		  for (unsigned int k3=0; k3<=k2; ++k3) {
		    
		    for (unsigned int i4=0; i4<=i3; ++i4) {
		      for (unsigned int j4=0; j4<=j3; ++j4) {
			for (unsigned int k4=0; k4<=k3; ++k4) {
			  
			  const unsigned int i5 = i1 + i4 + j4 + k4;
			  const unsigned int j5 = j1 + i3 - i4 + j3 - j4 + k3 - k4;
			  const unsigned int k5 = k1 + i2 - i3 + j2 - j3 + k2 - k3;
			  
			  /* 
			     ORSA_DEBUG("i5: %i",i5);
			     ORSA_DEBUG("j5: %i",j5);
			     ORSA_DEBUG("k5: %i",k5);
			  */
			  
			  orsa::Double fiveSumFx = orsa::zero();
			  orsa::Double fiveSumFy = orsa::zero();
			  orsa::Double fiveSumFz = orsa::zero();
			  
		          for (unsigned int L=0; L<=i5; ++L) {
			    for (unsigned int M=0; M<=j5; ++M) {
			      for (unsigned int N=0; N<=k5; ++N) {
				if (L != 0) {
				  fiveSumFx +=
				    orsa::Paul::t_lmnLMN::instance()->get(i5,j5,k5,L,M,N) *
				    L * csi1_PC.get(L-1) *
				    eta1_PC.get(M) *
				    zeta1_PC.get(N) *
				    oneOverR_PC.get(L+M+N);
				}
				//
				if (M != 0) {
				  fiveSumFy +=
				    orsa::Paul::t_lmnLMN::instance()->get(i5,j5,k5,L,M,N) *
				    csi1_PC.get(L) *
				    M * eta1_PC.get(M-1) *
				    zeta1_PC.get(N) *
				    oneOverR_PC.get(L+M+N);
				}
				//
				if (N != 0) {
				  fiveSumFz +=
				    orsa::Paul::t_lmnLMN::instance()->get(i5,j5,k5,L,M,N) *
				    csi1_PC.get(L) *
				    eta1_PC.get(M) *
				    N * zeta1_PC.get(N-1) *
				    oneOverR_PC.get(L+M+N);
				}
			      }
			    }
			  }
			  
		          for (unsigned int L=0; L<=i5; ++L) {
			    for (unsigned int M=0; M<=j5; ++M) {
			      for (unsigned int N=0; N<=k5; ++N) {
				fiveSumFx -=
				  (i5+j5+k5+L+M+N+1) *
				  orsa::Paul::t_lmnLMN::instance()->get(i5,j5,k5,L,M,N) *
				  csi1_PC.get(L+1) *
				  eta1_PC.get(M) *
				  zeta1_PC.get(N) *
				  oneOverR_PC.get(L+M+N+2);
				//
				fiveSumFy -=
				  (i5+j5+k5+L+M+N+1) *
				  orsa::Paul::t_lmnLMN::instance()->get(i5,j5,k5,L,M,N) *
				  csi1_PC.get(L) *
				  eta1_PC.get(M+1) *
				  zeta1_PC.get(N) *
				  oneOverR_PC.get(L+M+N+2);
				//
				fiveSumFz -=
				  (i5+j5+k5+L+M+N+1) *
				  orsa::Paul::t_lmnLMN::instance()->get(i5,j5,k5,L,M,N) *
				  csi1_PC.get(L) *
				  eta1_PC.get(M) *
				  zeta1_PC.get(N+1) *
				  oneOverR_PC.get(L+M+N+2);
			      }
			    }
			  }
			  
			  const orsa::Double commonInnerFactor = 
			    orsa::binomial(i2,i3) * orsa::binomial(i3,i4) *
			    orsa::binomial(j2,j3) * orsa::binomial(j3,j4) *
			    orsa::binomial(k2,k3) * orsa::binomial(k3,k4) *
			    lx_PC.get(i4) * mx_PC.get(i3-i4) * nx_PC.get(i2-i3) *
			    ly_PC.get(j4) * my_PC.get(j3-j4) * ny_PC.get(j2-j3) *
			    lz_PC.get(k4) * mz_PC.get(k3-k4) * nz_PC.get(k2-k3);
			  
			  innerFx += fiveSumFx * commonInnerFactor;
			  innerFy += fiveSumFy * commonInnerFactor;
			  innerFz += fiveSumFz * commonInnerFactor;
			  
			}
		      }
		    }
		    
		  }
		}
	      }
	      
	      /* 
		 const orsa::Double commonFactor = 
		 orsa::power_sign(i1+j1+k1) *
		 M1->M(i1,j1,k1) * 
		 M2->M(i2,j2,k2) *
		 oneOverR_PC.get(i1+j1+k1+i2+j2+k2+1) / 
		 orsa::Double( factorial(i1) * factorial(j1) * factorial(k1) *
		 factorial(i2) * factorial(j2) * factorial(k2) );
		 
		 Fx += innerFx * commonFactor;
		 Fy += innerFy * commonFactor;
		 Fz += innerFz * commonFactor;
	      */
	      
	      const orsa::Double commonFactor = 
		orsa::power_sign(i1+j1+k1) *
		M2->M(i2,j2,k2) *
		oneOverR_PC.get(i1+j1+k1+i2+j2+k2+1) / 
		orsa::Double( factorial(i1) * factorial(j1) * factorial(k1) *
			      factorial(i2) * factorial(j2) * factorial(k2) );
	      
	      /* 
		 ORSA_DEBUG("commonFactor: %Fe",commonFactor.get_mpf_t());
		 
		 ORSA_DEBUG("innerFx: %Fe",innerFx.get_mpf_t());
		 ORSA_DEBUG("innerFy: %Fe",innerFy.get_mpf_t());
		 ORSA_DEBUG("innerFz: %Fe",innerFz.get_mpf_t());
		 
		 ORSA_DEBUG("M1->M(i1+1,j1,k1): %Fe",M1->M(i1+1,j1,k1).get_mpf_t());
		 ORSA_DEBUG("M1->M(i1,j1+1,k1): %Fe",M1->M(i1,j1+1,k1).get_mpf_t());
		 ORSA_DEBUG("M1->M(i1,j1,k1+1): %Fe",M1->M(i1,j1,k1+1).get_mpf_t());
	      */
	      
	      Tx += commonFactor * (M1->M(i1,j1+1,k1)*innerFz - M1->M(i1,j1,k1+1)*innerFy);
	      Ty += commonFactor * (M1->M(i1,j1,k1+1)*innerFx - M1->M(i1+1,j1,k1)*innerFz);
	      Tz += commonFactor * (M1->M(i1+1,j1,k1)*innerFy - M1->M(i1,j1+1,k1)*innerFx);
	      
	      if (0) {
		ORSA_DEBUG("[%02i][%02i][%02i][%02i][%02i][%02i]   T: %20.12Fe %20.12Fe %20.12Fe",
			   i1,j1,k1,i2,j2,k2,
			   Tx.get_mpf_t(),
			   Ty.get_mpf_t(),
			   Tz.get_mpf_t());
	      }
	      
	    }
	  }
	}
	
      }
    }
  }
  
  /* 
     ORSA_DEBUG("Paul Tx: %.20Fe",Tx.get_mpf_t());
     ORSA_DEBUG("Paul Ty: %.20Fe",Ty.get_mpf_t());
     ORSA_DEBUG("Paul Tz: %.20Fe",Tz.get_mpf_t());
  */
  
  // rotate back, capital X,Y,Z
  const orsa::Double TX = l11*Tx + m11*Ty + n11*Tz;
  const orsa::Double TY = l21*Tx + m21*Ty + n21*Tz;
  const orsa::Double TZ = l31*Tx + m31*Ty + n31*Tz;
  
  const orsa::Vector T(-TX,-TY,-TZ);
  
  ORSA_DEBUG("Paul T.X: %.20Fe",T.getX().get_mpf_t());
  ORSA_DEBUG("Paul T.Y: %.20Fe",T.getY().get_mpf_t());
  ORSA_DEBUG("Paul T.Z: %.20Fe",T.getZ().get_mpf_t());
  
  return T;
}
