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

PaulMoment::PaulMoment() : osg::Referenced(true) {
  // _shape = 0;
  // _massDistribution = new UniformMassDistribution;
}

PaulMoment::PaulMoment(const int order) : osg::Referenced(true) {
  // _shape = 0;
  // _massDistribution = new UniformMassDistribution;
  _order = order;
  
  const unsigned int _order_plus_one = _order.getRef()+1;
  //
  {
    _M.resize(_order_plus_one);
    _M_uncertainty.resize(_order_plus_one);
    for (unsigned int i=0; i<_order_plus_one; ++i) {
      _M[i].resize(_order_plus_one-i);
      _M_uncertainty[i].resize(_order_plus_one-i);
      for (unsigned int j=0; j<_order_plus_one-i; ++j) {
	_M[i][j].resize(_order_plus_one-i-j);
	_M_uncertainty[i][j].resize(_order_plus_one-i-j);
      }
    }
  }
}

static Vector __randomVectorUtil(gsl_rng   * rnd,
				 const Box & boundingBox) {
  return Vector(boundingBox.getXMin()+(boundingBox.getXMax()-boundingBox.getXMin())*gsl_rng_uniform(rnd),
		boundingBox.getYMin()+(boundingBox.getYMax()-boundingBox.getYMin())*gsl_rng_uniform(rnd),
		boundingBox.getZMin()+(boundingBox.getZMax()-boundingBox.getZMin())*gsl_rng_uniform(rnd));
}

void PaulMoment::_randomVectorInside(std::vector<bool>  & in,
				     const unsigned int   sample_points,
				     const int            random_seed) {
  
  /* 
     ORSA_DEBUG("called _randomVectorInside(...): %i sample points, random seed: %i",
     sample_points,random_seed);
  */
  
  // GSL rng init
  gsl_rng * rnd = gsl_rng_alloc(gsl_rng_gfsr4);
  gsl_rng_set(rnd,random_seed);
  
  /* 
     in.clear();
     in.reserve(sample_points);
     const double boundingRadius = _shape->boundingRadius();
     unsigned int _iter = 0;
     while (_iter < sample_points) {
     in.push_back(_shape->isInside(__randomVectorUtil(rnd,boundingRadius)));
     ++_iter;
     }
  */
  
  in.clear();
  in.reserve(sample_points);
  const Box boundingBox = _shape->boundingBox();
  unsigned int _iter = 0;
  while (_iter < sample_points) {
    in.push_back(_shape->isInside(__randomVectorUtil(rnd,boundingBox)));
    ++_iter;
  }
  
  // GSL rng clean
  gsl_rng_free(rnd);
  
  if (0) {
    // debug output
    unsigned int _inside = 0;
    for (unsigned int j=0; j<sample_points; ++j) {
      if (in[j]) ++_inside;
    }
    ORSA_DEBUG("inside: %i (%5.2f\%)",_inside,double(100.0*double(_inside)/double(sample_points)));
  }
  
  // ORSA_DEBUG("done with _randomVectorInside(...).");
}

void PaulMoment::_centerOfMass(Vector                  & center_of_mass,
			       Vector                  & center_of_mass_uncertainty,
			       const std::vector<bool> & in,
			       const unsigned int        sample_points,
			       const int                 random_seed) {
  
  /* 
     ORSA_DEBUG("called _centerOfMass(...): %i sample points, random seed: %i",
     sample_points,random_seed);
  */
  
  // GSL rng init
  gsl_rng * rnd = gsl_rng_alloc(gsl_rng_gfsr4);
  gsl_rng_set(rnd,random_seed);
  
  osg::ref_ptr<orsa::Statistic<double> > _stat_M   = new orsa::Statistic<double>;
  //
  osg::ref_ptr<orsa::Statistic<double> > _stat_CMx = new orsa::Statistic<double>;
  osg::ref_ptr<orsa::Statistic<double> > _stat_CMy = new orsa::Statistic<double>;
  osg::ref_ptr<orsa::Statistic<double> > _stat_CMz = new orsa::Statistic<double>;
  //
  // const double _R_max = _shape->boundingRadius();
  // const double boundingRadius = _shape->boundingRadius();
  const Box boundingBox = _shape->boundingBox();
  unsigned int _iter = 0;
  while (_iter < sample_points) {
    const Vector & _v = __randomVectorUtil(rnd,boundingBox);
    if (in[_iter]) {
      const double density = _massDistribution->density(_v);
      //
      if (density > 0) {
	_stat_M->insert(density);
	//
	_stat_CMx->insert(density*_v.getX());
	_stat_CMy->insert(density*_v.getY());
	_stat_CMz->insert(density*_v.getZ());
      }
    }
    ++_iter;
  }
  
  /* 
     center_of_mass.set(_stat_CMx->average(),
     _stat_CMy->average(),
     _stat_CMz->average());
     
     center_of_mass_uncertainty.set(_stat_CMx->averageError(),
     _stat_CMy->averageError(),
     _stat_CMz->averageError());
  */
  
  center_of_mass.set(_stat_CMx->sum()/_stat_M->sum(),
		     _stat_CMy->sum()/_stat_M->sum(),
		     _stat_CMz->sum()/_stat_M->sum());
  
  // ORSA_DEBUG("IMPORTANT: complete density propagation to center_of_mass_uncertainty");
  center_of_mass_uncertainty.set(_stat_CMx->averageError(),
				 _stat_CMy->averageError(),
				 _stat_CMz->averageError());
  
  // GSL rng clean
  gsl_rng_free(rnd);
  
  if (1) {
    // debug output
    ORSA_DEBUG("cm.x: %14.6e +/- %14.6e",
	       center_of_mass.getX(),
	       center_of_mass_uncertainty.getX());
    ORSA_DEBUG("cm.y: %14.6e +/- %14.6e",
	       center_of_mass.getY(),
	       center_of_mass_uncertainty.getY());
    ORSA_DEBUG("cm.z: %14.6e +/- %14.6e",
	       center_of_mass.getZ(),
	       center_of_mass_uncertainty.getZ());
  }
  
  // ORSA_DEBUG("done with _centerOfMass(...).");
}

void PaulMoment::_inertiaMoment(Matrix                  & inertia_moment,
				Matrix                  & inertia_moment_uncertainty,
				const Vector            & center_of_mass,
				const Vector            & /* center_of_mass_uncertainty */,
				const std::vector<bool> & in,
				const unsigned int        sample_points,
				const int                 random_seed) {
  
  /* 
     ORSA_DEBUG("called _inertiaMoment(...): %i sample points, random seed: %i",
     sample_points,random_seed);
  */
  
  // GSL rng init
  gsl_rng * rnd = gsl_rng_alloc(gsl_rng_gfsr4);
  gsl_rng_set(rnd,random_seed);
  
  osg::ref_ptr<orsa::Statistic<double> > _stat_Ixx = new orsa::Statistic<double>;
  osg::ref_ptr<orsa::Statistic<double> > _stat_Iyy = new orsa::Statistic<double>;
  osg::ref_ptr<orsa::Statistic<double> > _stat_Izz = new orsa::Statistic<double>;
  osg::ref_ptr<orsa::Statistic<double> > _stat_Ixy = new orsa::Statistic<double>;
  osg::ref_ptr<orsa::Statistic<double> > _stat_Ixz = new orsa::Statistic<double>;
  osg::ref_ptr<orsa::Statistic<double> > _stat_Iyz = new orsa::Statistic<double>;
  //
  // const double _R_max = _shape->boundingRadius();
  // const double boundingRadius = _shape->boundingRadius();
  const Box boundingBox = _shape->boundingBox();
  unsigned int _iter = 0;
  while (_iter < sample_points) {
    // correct _v for the center-of-mass position
    // const Vector & _v = __randomVectorUtil(rnd,boundingRadius) - center_of_mass;
    const Vector & _v = __randomVectorUtil(rnd,boundingBox) - center_of_mass;
    if (in[_iter]) {
      const double density = _massDistribution->density(_v);
      //
      if (density > 0) {
	//
	_stat_Ixx->insert(density*(_v.getY()*_v.getY()+_v.getZ()*_v.getZ()));
	_stat_Iyy->insert(density*(_v.getX()*_v.getX()+_v.getZ()*_v.getZ()));
	_stat_Izz->insert(density*(_v.getX()*_v.getX()+_v.getY()*_v.getY()));
	//
	_stat_Ixy->insert(density*(-_v.getX()*_v.getY()));
	_stat_Ixz->insert(density*(-_v.getX()*_v.getZ()));
	_stat_Iyz->insert(density*(-_v.getY()*_v.getZ()));
      }
    }
    ++_iter;
  }
  
  inertia_moment.set(_stat_Ixx->average(),
		     _stat_Ixy->average(),
		     _stat_Ixz->average(),
		     _stat_Ixy->average(),
		     _stat_Iyy->average(),
		     _stat_Iyz->average(),
		     _stat_Ixz->average(),
		     _stat_Iyz->average(),
		     _stat_Izz->average());
  
  inertia_moment_uncertainty.set(_stat_Ixx->averageError(),
				 _stat_Ixy->averageError(),
				 _stat_Ixz->averageError(),
				 _stat_Ixy->averageError(),
				 _stat_Iyy->averageError(),
				 _stat_Iyz->averageError(),
				 _stat_Ixz->averageError(),
				 _stat_Iyz->averageError(),
				 _stat_Izz->averageError());
  
  // GSL rng clean
  gsl_rng_free(rnd);
  
  if (1) {
    // debug output    
    ORSA_DEBUG("Ixx:  %14.6e +/- %14.6e",
	       double(_stat_Ixx->average()),
	       double(_stat_Ixx->averageError()));
    ORSA_DEBUG("Iyy:  %14.6e +/- %14.6e",
	       double(_stat_Iyy->average()),
	       double(_stat_Iyy->averageError()));
    ORSA_DEBUG("Izz:  %14.6e +/- %14.6e",
	       double(_stat_Izz->average()),
	       double(_stat_Izz->averageError()));
    //
    ORSA_DEBUG("Ixy:  %14.6e +/- %14.6e",
	       double(_stat_Ixy->average()),
	       double(_stat_Ixy->averageError()));
    ORSA_DEBUG("Ixz:  %14.6e +/- %14.6e",
	       double(_stat_Ixz->average()),
	       double(_stat_Ixz->averageError()));
    ORSA_DEBUG("Iyz:  %14.6e +/- %14.6e",
	       double(_stat_Iyz->average()),
	       double(_stat_Iyz->averageError()));
  }
  
  // ORSA_DEBUG("done with _inertiaMoment(...).");
}

void PaulMoment::_moment(std::vector< std::vector< std::vector<double> > > & Mo,
			 std::vector< std::vector< std::vector<double> > > & Mo_uncertainty,
			 const Vector                                            & center_of_mass,
			 const Vector                                            & /* center_of_mass_uncertainty */,
			 const std::vector<bool>                                 & in,
			 const unsigned int                                        sample_points,
			 const int                                                 random_seed) {
  
  /* 
     ORSA_DEBUG("%i sample points, random seed: %i",
     sample_points,random_seed);
  */
  
  // GSL rng init
  gsl_rng * rnd = gsl_rng_alloc(gsl_rng_gfsr4);
  gsl_rng_set(rnd,random_seed);
  
  const unsigned int _order_plus_one = _order.getRef()+1;
  //
  osg::ref_ptr<orsa::Statistic<double> > _stat_M = new orsa::Statistic<double>;
  //
  std::vector< std::vector< std::vector< osg::ref_ptr< orsa::Statistic<double> > > > > _stat_Mo(_order_plus_one);
  //
  // std::vector< std::vector< std::vector< osg::ref_ptr< orsa::Statistic<double> > > > > _stat_Vo(_order_plus_one);
  //
  _stat_Mo.resize(_order_plus_one);
  // _stat_Vo.resize(_order_plus_one);
  for (unsigned int i=0; i<_order_plus_one; ++i) {
    _stat_Mo[i].resize(_order_plus_one-i);
    // _stat_Vo[i].resize(_order_plus_one-i);
    for (unsigned int j=0; j<_order_plus_one-i; ++j) {
      _stat_Mo[i][j].resize(_order_plus_one-i-j);
      // _stat_Vo[i][j].resize(_order_plus_one-i-j);
      for (unsigned int k=0; k<_order_plus_one-i-j; ++k) {
	_stat_Mo[i][j][k] = new orsa::Statistic<double>;
     	// _stat_Vo[i][j][k] = new orsa::Statistic<double>;
      }
    }
  }
  
  const Box boundingBox = _shape->boundingBox();
  unsigned int _iter = 0;
  while (_iter < sample_points) {
    // correct _v for the center-of-mass position
    const Vector & _v = __randomVectorUtil(rnd,boundingBox) - center_of_mass;
    //
    /* 
       for (unsigned int i=0; i<_order_plus_one; ++i) {
       for (unsigned int j=0; j<_order_plus_one-i; ++j) {
       for (unsigned int k=0; k<_order_plus_one-i-j; ++k) {
       _stat_Vo[i][j][k]->insert(int_pow(_v.getX(),i+1)*
       int_pow(_v.getY(),j+1)*
       int_pow(_v.getZ(),k+1)/double((i+1)*(j+1)*(k+1)));
       }
       }
       }
    */
    //
    if (in[_iter]) {
      const double density = _massDistribution->density(_v + center_of_mass);
      //
      if (density > 0) {
	//
	// ORSA_DEBUG("density: %f",density());
	_stat_M->insert(density);
	//
	for (unsigned int i=0; i<_order_plus_one; ++i) {
	  for (unsigned int j=0; j<_order_plus_one-i; ++j) {
	    for (unsigned int k=0; k<_order_plus_one-i-j; ++k) {
	      _stat_Mo[i][j][k]->insert(density*
					int_pow(_v.getX(),i)*
					int_pow(_v.getY(),j)*
					int_pow(_v.getZ(),k));
	    }
	  }
	}
      }
    }
    ++_iter;
  }
  
  // ORSA_DEBUG("IMPORTANT: complete the density propagation...");
  
  Mo.resize(_order_plus_one);
  Mo_uncertainty.resize(_order_plus_one);
  for (unsigned int i=0; i<_order_plus_one; ++i) {
    Mo[i].resize(_order_plus_one-i);
    Mo_uncertainty[i].resize(_order_plus_one-i);
    for (unsigned int j=0; j<_order_plus_one-i; ++j) {
      Mo[i][j].resize(_order_plus_one-i-j);
      Mo_uncertainty[i][j].resize(_order_plus_one-i-j);
    }
  }
  //
  // const double _test_norm = 1/(Legendre::norm(0,0)*_stat_C[0][0]->average());
  // const double _test_norm = _stat_M->sum()/(Legendre::norm(0,0)*_stat_C[0][0]->sum());
  // const double _test_norm = _stat_M->sum() / _stat_Mo[0][0][0]->sum();
  //
  /* 
     ORSA_DEBUG("_test_norm: %g",
     _test_norm());
  */
  //
  /* 
     for (unsigned int l=0; l<=_order.getRef(); ++l) {
     for (unsigned int m=0; m<=l; ++m) {
     C[l][m]             = _test_norm*Legendre::norm(l,m)*_stat_C[l][m]->average();
     C_uncertainty[l][m] = _test_norm*Legendre::norm(l,m)*_stat_C[l][m]->averageError();
     S[l][m]             = _test_norm*Legendre::norm(l,m)*_stat_S[l][m]->average();
     S_uncertainty[l][m] = _test_norm*Legendre::norm(l,m)*_stat_S[l][m]->averageError();
     }
     }
  */
  
  /* 
     for (unsigned int i=0; i<_order_plus_one; ++i) {
     for (unsigned int j=0; j<_order_plus_one-i; ++j) {
     for (unsigned int k=0; k<_order_plus_one-i-j; ++k) {
     Mo[i][j][k]             = _test_norm * _stat_Mo[i][j][k]->sum()/_stat_M->sum();
     Mo_uncertainty[i][j][k] = _test_norm * _stat_Mo[i][j][k]->averageError();
     }
     }
     }
  */
  
  for (unsigned int i=0; i<_order_plus_one; ++i) {
    for (unsigned int j=0; j<_order_plus_one-i; ++j) {
      for (unsigned int k=0; k<_order_plus_one-i-j; ++k) {
	/* 
	   Mo[i][j][k] =
	   _stat_Mo[i][j][k]->sum() / _stat_M->sum() /
	   ( int_pow(boundingBox.getXMax() - 
	   boundingBox.getXMin(),
	   i) *
	   int_pow(boundingBox.getYMax() - 
	   boundingBox.getYMin(),
	   j) *
	   int_pow(boundingBox.getZMax() - 
	   boundingBox.getZMin(),
	   k) );
	*/
	//
	/* 
	   Mo[i][j][k] =
	   _stat_Mo[i][j][k]->sum() / _stat_M->sum() /
	   int_pow( int_pow(boundingBox.getXMax() - 
	   boundingBox.getXMin(),
	   i+1) *
	   int_pow(boundingBox.getYMax() - 
	   boundingBox.getYMin(),
	   j+1) *
	   int_pow(boundingBox.getZMax() - 
	   boundingBox.getZMin(),
	   k+1),
	   1);
	*/
	//
	// Mo[i][j][k] = _stat_Mo[i][j][k]->average();
	//
	Mo[i][j][k] = _stat_Mo[i][j][k]->sum() / _stat_M->sum();
	//
	// Mo[i][j][k] = _stat_Mo[i][j][k]->average();
	//
	/* 
	   Mo_uncertainty[i][j][k] = 
	   _stat_Mo[i][j][k]->averageError() /
	   ( int_pow(boundingBox.getXMax() - 
	   boundingBox.getXMin(),
	   i) *
	   int_pow(boundingBox.getYMax() - 
	   boundingBox.getYMin(),
	   j) *
	   int_pow(boundingBox.getZMax() - 
	   boundingBox.getZMin(),
	   k) );
	*/
	//
	/* 
	   Mo_uncertainty[i][j][k] = 
	   _stat_Mo[i][j][k]->averageError() /
	   ( int_pow(boundingBox.getXMax() - 
	   boundingBox.getXMin(),
	   i) *
	   int_pow(boundingBox.getYMax() - 
	   boundingBox.getYMin(),
	   j) *
	   int_pow(boundingBox.getZMax() - 
	   boundingBox.getZMin(),
	   k) );
	*/
	//
	/* 
	   Mo_uncertainty[i][j][k] = 
	   _stat_Mo[i][j][k]->averageError() /
	   int_pow( int_pow(boundingBox.getXMax() - 
	   boundingBox.getXMin(),
	   i+1) *
	   int_pow(boundingBox.getYMax() - 
	   boundingBox.getYMin(),
	   j+1) *
	   int_pow(boundingBox.getZMax() - 
	   boundingBox.getZMin(),
	   k+1),
	   1);
	*/
	//
	Mo_uncertainty[i][j][k] = _stat_Mo[i][j][k]->averageError();
	//
	/* 
	   Mo_uncertainty[i][j][k] = 
	   _stat_Mo[i][j][k]->averageError();
	*/
      }
    }
  }
  // one more time, so that Mo[0][0][0] == 1.0 (by definition...)
  if (0) {
    ORSA_DEBUG("renormalizing...");
    const double localNorm = Mo[0][0][0];
    for (unsigned int i=0; i<_order_plus_one; ++i) {
      for (unsigned int j=0; j<_order_plus_one-i; ++j) {
	for (unsigned int k=0; k<_order_plus_one-i-j; ++k) {
	  Mo[i][j][k]             /= localNorm;
	  Mo_uncertainty[i][j][k] /= localNorm;
	}
      }
    }
  }
  
  if (0) {
    // remove "not significant" values
    for (unsigned int i=0; i<_order_plus_one; ++i) {
      for (unsigned int j=0; j<_order_plus_one-i; ++j) {
	for (unsigned int k=0; k<_order_plus_one-i-j; ++k) {
	  
	  if (3*Mo_uncertainty[i][j][k] > fabs(Mo[i][j][k])) {
	    Mo[i][j][k] = 0;
	  }	
	  
	}
      }
    }
  }
  
  // GSL rng clean
  gsl_rng_free(rnd);
  
  if (1) {
    // debug output
    
    for (unsigned int localOrder=0; localOrder<_order_plus_one; ++localOrder) {
      
      for (unsigned int i=0; i<_order_plus_one; ++i) {
	for (unsigned int j=0; j<_order_plus_one-i; ++j) {
	  for (unsigned int k=0; k<_order_plus_one-i-j; ++k) {
	    
	    if (i+j+k == localOrder) {
	      // if (Mo[i][j][k] > 0) {
	      ORSA_DEBUG("M[%02i][%02i][%02i]: %14.6e +/- %14.6e",
			 i,j,k,
			 Mo[i][j][k],
			 Mo_uncertainty[i][j][k]);
	      // }
	    }
	    
	  }
	}
      }
      
    }
    
  }
  
  
  if (1) {
    // debug output
    
    // MORE debugging, with N_{ijk}
    
    // volume
    /* 
       unsigned int _inside = 0;
       for (unsigned int j=0; j<sample_points; ++j) {
       if (in[j]) ++_inside;
       }
       const double volume = _shape->boundingBox().volume()*double(_inside)/double(sample_points);
    */
    
    const double r0 = _shape->boundingRadius();
    //
    ORSA_DEBUG("r0: %g [km]",
	       FromUnits(r0,Unit::KM,-1));    
    
    for (unsigned int localOrder=0; localOrder<_order_plus_one; ++localOrder) {
      
      for (unsigned int i=0; i<_order_plus_one; ++i) {
	for (unsigned int j=0; j<_order_plus_one-i; ++j) {
	  for (unsigned int k=0; k<_order_plus_one-i-j; ++k) {
	    
	    if (i+j+k == localOrder) {
	      // if (Mo[i][j][k] > 0) {
	      /* 
		 ORSA_DEBUG("M[%02i][%02i][%02i]: %14.6e +/- %14.6e",
		 i,j,k,
		 Mo[i][j][k](),
		 Mo_uncertainty[i][j][k]());
	      */
	      //
	      ORSA_DEBUG("N[%02i][%02i][%02i]: %14.6e +/- %14.6e",
			 i,j,k,
			 double(Mo[i][j][k]/(int_pow(r0,i+j+k))),
			 double(Mo_uncertainty[i][j][k]/(int_pow(r0,i+j+k))));
	      // }
	    }
	  }
	}
      }
    }
  }
}

bool PaulMoment::computeUsingShape(const unsigned int sample_points,
				   const int random_seed) {
  if (_shape==0) {
    ORSA_ERROR("unset Shape");
    return false;
  }
  
  if (!_order.isSet()) {
    ORSA_ERROR("unset order");
    return false;
  }
  
  /*
   * Three phases: 
   * 1) center-of-mass computation
   * 2) inertia moments computation
   * 3) multipole computation
   *
   * at each phase, the random number generation is reset,
   * so that we always use the same points at each phase.
   */
  
  _sample_points.set(sample_points);
  
  // vector of bool, for generated points inside Shape
  std::vector<bool> _in;
  _randomVectorInside(_in,sample_points,random_seed);
  // complete this version, but we need the name of the shape model...
  /* 
     std::vector<bool> _in;
     {
     char filename[1024];
     gmp_sprintf(filename,"eros_in_vector_%f_%i_%i.dat",
     _R0.get(),
     sample_points,
     random_seed);
     if (!_read_in_vector_from_file(_in,filename)) {
     _randomVectorInside(_in,sample_points,random_seed);
     _write_in_vector_to_file(_in,filename);
     }
     }
  */
  
  if (1) {
    // volume
    unsigned int _inside = 0;
    for (unsigned int j=0; j<sample_points; ++j) {
      if (_in[j]) ++_inside;
    }
    ORSA_DEBUG("volume: %f +/- %f  km^3",
	       FromUnits(_shape->boundingBox().volume()*double(_inside)/double(sample_points),Unit::KM,-3),
	       FromUnits(_shape->boundingBox().volume()*sqrt(double(_inside))/double(sample_points),Unit::KM,-3));
  }
  
  {
    // center-of-mass
    orsa::Vector tmp_center_of_mass, tmp_center_of_mass_uncertainty;
    _centerOfMass(tmp_center_of_mass,
		  tmp_center_of_mass_uncertainty,
		  _in,
		  sample_points,
		  random_seed);
    _center_of_mass             = tmp_center_of_mass;
    _center_of_mass_uncertainty = tmp_center_of_mass_uncertainty;
  }
  
  {
    // inertia moment
    orsa::Matrix tmp_inertia_moment, tmp_inertia_moment_uncertainty;
    _inertiaMoment(tmp_inertia_moment,
		   tmp_inertia_moment_uncertainty,
		   _center_of_mass.getRef(),
		   _center_of_mass_uncertainty.getRef(),
		   _in,
		   sample_points,
		   random_seed);
    _inertia_moment             = tmp_inertia_moment;
    _inertia_moment_uncertainty = tmp_inertia_moment_uncertainty;
  }
  
  {
    _moment(_M,
	    _M_uncertainty,
	    _center_of_mass.getRef(),
	    _center_of_mass_uncertainty.getRef(),
	    _in,
	    sample_points,
	    random_seed);
  }
  
  return true;
}

/* 
   double normalization(const unsigned int l,
   const unsigned int m) {
   return sqrt(orsa::factorial(l+m).get_d() / ((2-orsa::kronecker(m,0))*(2*l+1)*orsa::factorial(l-m).get_d()));
   }
*/
//
double normalization(const unsigned int l,
		     const unsigned int m) {
  // Cross checked with NEAR-Eros papers
  return sqrt( ((2-orsa::kronecker(m,0))*orsa::factorial(l-m).get_d()) / 
	       ((2*l+1)*orsa::factorial(l+m).get_d()) );
}

void orsa::convert(const PaulMoment * const pm,
		   const double     & R0) {
  
#warning add uncertainty estimates...
  
  ORSA_DEBUG("R0: %f [km]",orsa::FromUnits(R0,orsa::Unit::KM,-1));
  
  const int order = pm->order();
  
  // C_lm coefficients
  //
  for (int l=0; l<=order; ++l) {
    for (int m=0; m<=l; ++m) {
      
      double pq_factor=0;
      //
      // integer division in limits
      for (int p=0;p<=(l/2);++p) {
	for (int q=0;q<=(m/2);++q) {
	  
	  double nu_factor=0;
	  //
	  for (int nu_x=0; nu_x<=p; ++nu_x) {
	    for (int nu_y=0; nu_y<=(p-nu_x); ++nu_y) {
	      
	      const int M_i = m-2*q+2*nu_x;
	      const int M_j = 2*q+2*nu_y;
	      const int M_k = l-m-2*nu_x-2*nu_y;
	      
	      if (M_i+M_j+M_k!=l) {
		ORSA_DEBUG("WARNING!!!");
	      }
	      
	      if ( (M_i>=0) && 
		   (M_j>=0) && 
		   (M_k>=0) && 
		   (M_i+M_j+M_k==l) ) {
		// ORSA_DEBUG("requesting M[%i][%i][%i]...",m-2*q+2*nu_x,2*q+2*nu_y,l-m-2*nu_x-2*nu_y);
		//
		nu_factor += 
		  (orsa::factorial(p).get_d() /(orsa::factorial(nu_x).get_d()*orsa::factorial(nu_y).get_d()*orsa::factorial(p-nu_x-nu_y).get_d())) *
		  pm->M(M_i,M_j,M_k);
	      }
	    }
	  }
	  // need this because using M(i,j,k) instead of N(i,j,k)
	  nu_factor /= int_pow(R0,l);
	  
	  pq_factor += 
	    orsa::power_sign(p+q) *
	    orsa::binomial(l,p).get_d() *
	    orsa::binomial(2*l-2*p,l).get_d() *
	    orsa::binomial(m,2*q).get_d() *
	    orsa::pochhammer(l-m-2*p+1,m) * 
	    nu_factor;
	}
      }
      //
      pq_factor /= int_pow(2,l);
      //
      const double C_lm = pq_factor;
      //      
      const double norm_C_lm = C_lm*normalization(l,m);
      
      ORSA_DEBUG("     C[%i][%i] = %+16.12f",
		 l,m,     C_lm);
      ORSA_DEBUG("norm_C[%i][%i] = %+16.12f   norm: %f",
		 l,m,norm_C_lm,normalization(l,m));
      
      if ((l>=2) && (m==0)) {
	// J_l is minus C_l0, not-normalized
	const double J_l = -C_lm;
	ORSA_DEBUG("J_%i = %+16.12f",l,J_l);
      }	
      
    }
  }
  
  // S_lm coefficients
  //
  for (int l=0; l<=order; ++l) {
    for (int m=1; m<=l; ++m) {
      
      double pq_factor=0;
      //
      // integer division in limits
      for (int p=0;p<=(l/2);++p) {
	for (int q=0;q<=((m-1)/2);++q) {
	  
	  double nu_factor=0;
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
		// ORSA_DEBUG("requesting M[%i][%i][%i]...",m-2*q+2*nu_x,2*q+2*nu_y,l-m-2*nu_x-2*nu_y);
		//
		nu_factor += 
		  (orsa::factorial(p).get_d() /(orsa::factorial(nu_x).get_d()*orsa::factorial(nu_y).get_d()*orsa::factorial(p-nu_x-nu_y).get_d())) *
		  pm->M(M_i,M_j,M_k);
	      }
	    }
	  }
	  // need this because using M(i,j,k) instead of N(i,j,k)
	  nu_factor /= int_pow(R0,l);
	  
	  pq_factor += 
	    orsa::power_sign(p+q) *
	    orsa::binomial(l,p).get_d() *
	    orsa::binomial(2*l-2*p,l).get_d() *
	    orsa::binomial(m,2*q+1).get_d() *
	    orsa::pochhammer(l-m-2*p+1,m) * 
	    nu_factor;
	}
      }
      //
      pq_factor /= int_pow(2,l);
      //
      const double S_lm = pq_factor;
      //      
      const double norm_S_lm = S_lm*normalization(l,m);
      
      ORSA_DEBUG("     S[%i][%i] = %+16.12f",
		 l,m,     S_lm);
      ORSA_DEBUG("norm_S[%i][%i] = %+16.12f   norm: %f",
		 l,m,norm_S_lm,normalization(l,m));
      
    }
  }
  
}
