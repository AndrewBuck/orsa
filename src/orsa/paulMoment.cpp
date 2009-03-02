#include <orsa/paulMoment.h>

#include <orsa/box.h>
#include <orsa/double.h>
#include <orsa/multipole.h>
#include <orsa/legendre.h>
#include <orsa/statistic.h>
#include <orsa/unit.h>

#include <vector>

#include <gsl/gsl_rng.h>

#include <iostream>

using namespace orsa;

PaulMoment::PaulMoment() : osg::Referenced(true) {
  _shape = 0;
  _massDistribution = new UniformMassDistribution;
}

PaulMoment::PaulMoment(const unsigned int order) : osg::Referenced(true) {
  _shape = 0;
  _massDistribution = new UniformMassDistribution;
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

PaulMoment::~PaulMoment() { 
  
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
    ORSA_DEBUG("inside: %i (%F5.2f\%)",_inside,double(100.0*double(_inside)/double(sample_points)));
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
      _stat_M->insert(density);
      //
      _stat_CMx->insert(density*_v.getX());
      _stat_CMy->insert(density*_v.getY());
      _stat_CMz->insert(density*_v.getZ());
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
  
  ORSA_DEBUG("IMPORTANT: complete density propagation to center_of_mass_uncertainty");
  center_of_mass_uncertainty.set(_stat_CMx->averageError(),
				 _stat_CMy->averageError(),
				 _stat_CMz->averageError());
  
  // GSL rng clean
  gsl_rng_free(rnd);
  
  if (0) {
    // debug output
    ORSA_DEBUG("cm.x: %14.6Fe +/- %14.6Fe",
	       center_of_mass.getX(),
	       center_of_mass_uncertainty.getX());
    ORSA_DEBUG("cm.y: %14.6Fe +/- %14.6Fe",
	       center_of_mass.getY(),
	       center_of_mass_uncertainty.getY());
    ORSA_DEBUG("cm.z: %14.6Fe +/- %14.6Fe",
	       center_of_mass.getZ(),
	       center_of_mass_uncertainty.getZ());
  }
  
  // ORSA_DEBUG("done with _centerOfMass(...).");
}

void PaulMoment::_inertiaMoment(Matrix                  & inertia_moment,
				Matrix                  & inertia_moment_uncertainty,
				const Vector            & center_of_mass,
				const Vector            & center_of_mass_uncertainty,
				const std::vector<bool> & in,
				const unsigned int        sample_points,
				const int                 random_seed) {
  
  /* 
     ORSA_DEBUG("called _inertiaMoment(...): %i sample points, random seed: %i",
     sample_points,random_seed);
  */
  
  ORSA_DEBUG("IMPORTANT: complete density propagation in _inertiaMoment computation");
  
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
      _stat_Ixx->insert(_v.getY()*_v.getY()+_v.getZ()*_v.getZ());
      _stat_Iyy->insert(_v.getX()*_v.getX()+_v.getZ()*_v.getZ());
      _stat_Izz->insert(_v.getX()*_v.getX()+_v.getY()*_v.getY());
      //
      _stat_Ixy->insert(-_v.getX()*_v.getY());
      _stat_Ixz->insert(-_v.getX()*_v.getZ());
      _stat_Iyz->insert(-_v.getY()*_v.getZ());
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
  
  if (0) {
    // debug output    
    ORSA_DEBUG("Ixx:  %14.6Fe +/- %14.6Fe",
	       double(_stat_Ixx->average()),
	       double(_stat_Ixx->averageError()));
    ORSA_DEBUG("Iyy:  %14.6Fe +/- %14.6Fe",
	       double(_stat_Iyy->average()),
	       double(_stat_Iyy->averageError()));
    ORSA_DEBUG("Izz:  %14.6Fe +/- %14.6Fe",
	       double(_stat_Izz->average()),
	       double(_stat_Izz->averageError()));
    //
    ORSA_DEBUG("Ixy:  %14.6Fe +/- %14.6Fe",
	       double(_stat_Ixy->average()),
	       double(_stat_Ixy->averageError()));
    ORSA_DEBUG("Ixz:  %14.6Fe +/- %14.6Fe",
	       double(_stat_Ixz->average()),
	       double(_stat_Ixz->averageError()));
    ORSA_DEBUG("Iyz:  %14.6Fe +/- %14.6Fe",
	       double(_stat_Iyz->average()),
	       double(_stat_Iyz->averageError()));
  }
  
  // ORSA_DEBUG("done with _inertiaMoment(...).");
}

void PaulMoment::_moment(std::vector< std::vector< std::vector<double> > > & Mo,
			 std::vector< std::vector< std::vector<double> > > & Mo_uncertainty,
			 const Vector                                            & center_of_mass,
			 const Vector                                            & center_of_mass_uncertainty,
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
    ++_iter;
  }
  
  ORSA_DEBUG("IMPORTANT: complete the density propagation...");
  
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
     ORSA_DEBUG("_test_norm: %Fg",
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
  
  if (1) {
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
  
  if (0) {
    // debug output
    
    for (unsigned int localOrder=0; localOrder<_order_plus_one; ++localOrder) {
      
      for (unsigned int i=0; i<_order_plus_one; ++i) {
	for (unsigned int j=0; j<_order_plus_one-i; ++j) {
	  for (unsigned int k=0; k<_order_plus_one-i-j; ++k) {
	    
	    if (i+j+k == localOrder) {
	      if (Mo[i][j][k] > 0) {
		ORSA_DEBUG("M[%02i][%02i][%02i]: %14.6Fe +/- %14.6Fe",
			   i,j,k,
			   Mo[i][j][k],
			   Mo_uncertainty[i][j][k]);
	      }
	    }
	    
	  }
	}
      }
      
    }
    
  }
  
  
  if (0) {
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
    ORSA_DEBUG("r0: %Fg [km]",
	       FromUnits(r0,Unit::KM,-1));    
    
    for (unsigned int localOrder=0; localOrder<_order_plus_one; ++localOrder) {
      
      for (unsigned int i=0; i<_order_plus_one; ++i) {
	for (unsigned int j=0; j<_order_plus_one-i; ++j) {
	  for (unsigned int k=0; k<_order_plus_one-i-j; ++k) {
	    
	    if (i+j+k == localOrder) {
	      // if (Mo[i][j][k] > 0) {
	      /* 
		 ORSA_DEBUG("M[%02i][%02i][%02i]: %14.6Fe +/- %14.6Fe",
		 i,j,k,
		 Mo[i][j][k](),
		 Mo_uncertainty[i][j][k]());
	      */
	      //
	      ORSA_DEBUG("N[%02i][%02i][%02i]: %14.6Fe +/- %14.6Fe",
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
   bool PaulMoment::writeToFile(const std::string & filename) const {
   // ORSA_DEBUG("saving multipole to file [%s]",filename.c_str());
   FILE * fp = fopen(filename.c_str(),"w");
   if (!fp) {
   ORSA_WARNING("cannot open file %s",filename.c_str());
   return false;
   }
   gmp_fprintf(fp,"%18.10Fe\n",FromUnits(_R0.getRef(),Unit::M,-1));
   gmp_fprintf(fp,"%i\n",_order.getRef());
   gmp_fprintf(fp,"%i\n",_sample_points.getRef());
   //
   gmp_fprintf(fp,"%18.10Fe %18.10Fe\n",
   FromUnits(_center_of_mass.getX(),Unit::M,-1),
   FromUnits(_center_of_mass_uncertainty.getX(),Unit::M,-1));
   gmp_fprintf(fp,"%18.10Fe %18.10Fe\n",
   FromUnits(_center_of_mass.getY(),Unit::M,-1),
   FromUnits(_center_of_mass_uncertainty.getY(),Unit::M,-1));
   gmp_fprintf(fp,"%18.10Fe %18.10Fe\n",
   FromUnits(_center_of_mass.getZ(),Unit::M,-1),
   FromUnits(_center_of_mass_uncertainty.getZ(),Unit::M,-1));
   //
   for (unsigned int l=0; l<_C.size(); ++l) {
   for (unsigned int m=0; m<=l; ++m) {
   gmp_fprintf(fp,
   "%18.10Fe %18.10Fe\n",
   _C[l][m](),
   _C_uncertainty[l][m]());
   gmp_fprintf(fp,
   "%18.10Fe %18.10Fe\n",
   _S[l][m](),
   _S_uncertainty[l][m]());
   }
   }
   fclose(fp);
   return true;
   }
*/

/* 
   bool PaulMoment::readFromFile(const std::string & filename) {
   // ORSA_DEBUG("trying to read multipole from file [%s]",filename.c_str());
   FILE * fp = fopen(filename.c_str(),"r");
   if (!fp) {
   ORSA_WARNING("cannot open file %s",filename.c_str());
   return false;
   }
   {
   double R0;
   gmp_fscanf(fp,"%Fe",
   R0()); // '&' not required for GMP vars
   _R0.set(FromUnits(R0,Unit::M));
   }
   {
   unsigned int _n_order;
   gmp_fscanf(fp,"%i",&_n_order);
   _order.set(_n_order);
   // ORSA_DEBUG("read order: %i",_n_order);
   }
   //
   {
   unsigned int _n_sample_points;
   gmp_fscanf(fp,"%i",&_n_sample_points);
   _sample_points.set(_n_sample_points);
   // ORSA_DEBUG("read sample points: %i",_n_sample_points);
   }
   //
   {
   double x, dx, y, dy, z, dz;
   //
   gmp_fscanf(fp,"%Fe %Fe",
   FromUnits(x,Unit::M),
   FromUnits(dx,Unit::M)); // '&' not required for GMP vars
   gmp_fscanf(fp,"%Fe %Fe",
   FromUnits(y,Unit::M),
   FromUnits(dy,Unit::M)); // '&' not required for GMP vars
   gmp_fscanf(fp,"%Fe %Fe",
   FromUnits(z,Unit::M),
   FromUnits(dz,Unit::M)); // '&' not required for GMP vars
   //
   _center_of_mass.set(x,y,z);
   _center_of_mass_uncertainty.set(dx,dy,dz);
   //
   // ORSA_DEBUG("CM...");
   }
   //
   const unsigned int _order_plus_one = _order.getRef()+1;
   //
   {
   _C.resize(_order_plus_one);
   _C_uncertainty.resize(_order_plus_one);
   _S.resize(_order_plus_one);
   _S_uncertainty.resize(_order_plus_one);
   for (unsigned int l=0; l<_order_plus_one; ++l) {
   _C[l].resize(l+1);
   _C_uncertainty[l].resize(l+1);
   _S[l].resize(l+1);
   _S_uncertainty[l].resize(l+1);
   }
   }
   //
   {
   double _c, _dc, _s, _ds;
   for (unsigned int l=0; l<_order_plus_one; ++l) {
   for (unsigned int m=0; m<=l; ++m) {
   gmp_fscanf(fp,"%Fe %Fe",
   _c(),
   _dc()); // '&' not required for GMP vars
   _C[l][m] = _c;
   _C_uncertainty[l][m] = _dc;
   //
   gmp_fscanf(fp,"%Fe %Fe",
   _s(),
   _ds()); // '&' not required for GMP vars
   _S[l][m] = _s;
   _S_uncertainty[l][m] = _ds;
   //
   }
   }
   }
   fclose(fp);
   return true;
   }
*/

/* 
   bool PaulMoment::_write_in_vector_to_file(const std::vector<bool> & in,
   const std::string & filename) const {
   // ORSA_DEBUG("saving in vector to file [%s]",filename.c_str());
   FILE * fp = fopen(filename.c_str(),"w");
   if (!fp) {
   ORSA_WARNING("cannot open file %s",filename.c_str());
   return false;
   }
   gmp_fprintf(fp,"%i\n",in.size());
   for (unsigned int k=0; k<in.size(); ++k) {
   gmp_fprintf(fp,"%i\n",in[k]);
   }	
   fclose(fp);
   return true;
   }
*/

/* 
   bool PaulMoment::_read_in_vector_from_file(std::vector<bool> & in,
   const std::string & filename) {
   // ORSA_DEBUG("trying to read in vector from file [%s]",filename.c_str());
   FILE * fp = fopen(filename.c_str(),"r");
   if (!fp) {
   ORSA_WARNING("cannot open file %s",filename.c_str());
   return false;
   }
   in.clear();
   {
   unsigned int size;
   gmp_fscanf(fp,"%i",&size);
   in.resize(size);
   ORSA_DEBUG("read size: %i",size);
   }
   {
   unsigned int b;
   for (unsigned int k=0; k<in.size(); ++k) {
   gmp_fscanf(fp,"%i",&b);
   in[k] = b;
   }
   }
   return true;
   }	
*/

/* 
   bool PaulMoment::readFromAltFile(const std::string & filename) {
   // ORSA_DEBUG("trying to read multipole from alternative file [%s]",filename.c_str());
   FILE * fp = fopen(filename.c_str(),"r");
   if (!fp) {
   ORSA_WARNING("cannot open file %s",filename.c_str());
   return false;
   }
   
   _sample_points.set(0);
   
   char line[1024];
   
   // read multipole order
   if (fgets(line,1024,fp) != NULL) {
   unsigned int file_order;
   gmp_sscanf(line,"%i",&file_order);
   _order.set(file_order);
   // ORSA_DEBUG("read order: %i",file_order);
   } else {
   ORSA_ERROR("problems while reading file [%s]",filename.c_str());
   return false;
   }
   
   const unsigned int _order_plus_one = _order.getRef()+1;
   //
   {
   _C.resize(_order_plus_one);
   _C_uncertainty.resize(_order_plus_one);
   _S.resize(_order_plus_one);
   _S_uncertainty.resize(_order_plus_one);
   for (unsigned int l=0; l<_order_plus_one; ++l) {
   _C[l].resize(l+1);
   _C_uncertainty[l].resize(l+1);
   _S[l].resize(l+1);
   _S_uncertainty[l].resize(l+1);
   }
   }
   
   // don't use next line
   if (fgets(line,1024,fp) == NULL) {
   ORSA_ERROR("problems while reading file [%s]",filename.c_str());
   return false;
   }
   
   // read R0
   if (fgets(line,1024,fp) != NULL) {
   double R0;
   gmp_sscanf(line,"%*s %Fe",
   R0()); // '&' not required for GMP vars
   _R0.set(FromUnits(R0,Unit::KM));
   // ORSA_DEBUG("read R0: %Fg [km]",FromUnits(_R0.getRef(),Unit::KM,-1));
   } else {
   ORSA_ERROR("problems while reading file [%s]",filename.c_str());
   return false;
   }
   
   // don't use next line
   if (fgets(line,1024,fp) == NULL) {
   ORSA_ERROR("problems while reading file [%s]",filename.c_str());
   return false;
   }
   
   // don't use next line
   if (fgets(line,1024,fp) == NULL) {
   ORSA_ERROR("problems while reading file [%s]",filename.c_str());
   return false;
   }
   
   double file_entry;
   for (unsigned int l=0; l<_order_plus_one; ++l) {
   for (unsigned int m=0; m<=l; ++m) {
   // C
   gmp_fscanf(fp,"%Fg",file_entry());
   // ORSA_DEBUG("file_entry: %Fg",file_entry());      
   _C[l][m]             = file_entry;
   _C_uncertainty[l][m] = 0;
   //
   // ORSA_DEBUG("C[%02i][%02i] = %16.12f",l,m,_C[l][m]());      
   
   // S
   if (m==0) {
   _S[l][m]             = 0;
   _S_uncertainty[l][m] = 0;
   } else {
   gmp_fscanf(fp,"%Fg",file_entry());
   // ORSA_DEBUG("file_entry: %Fg",file_entry());     
   //
   _S[l][m]             = file_entry;
   _S_uncertainty[l][m] = 0;
   }
   //
   // ORSA_DEBUG("S[%02i][%02i] = %16.12f",l,m,_S[l][m]());      
   }
   }
   fclose(fp);
   return true;
   }
*/

/* 
   bool PaulMoment::readFromBisFile(const std::string & filename) {
   
   // not well tested, deduced by only one file...
   
   // ORSA_DEBUG("trying to read multipole from alternative file [%s]",filename.c_str());
   FILE * fp = fopen(filename.c_str(),"r");
   if (!fp) {
   ORSA_WARNING("cannot open file %s",filename.c_str());
   return false;
   }
   
   _sample_points.set(0);
   
   char line[1024];
   
   // skip first two lines, that look like:
   // &odpNamelist
   // OBNORM = .TRUE.,
   // 
   fgets(line,1024,fp);
   fgets(line,1024,fp);
   
   // read R0
   if (fgets(line,1024,fp) != NULL) {
   double R0;
   gmp_sscanf(line,"OBR = %f,",
   R0()); // '&' not required for GMP vars
   _R0.set(FromUnits(R0,Unit::KM));
   //
   ORSA_DEBUG("read R0: %Fg [km]",FromUnits(_R0.getRef(),Unit::KM,-1));
   } else {
   ORSA_ERROR("problems while reading file [%s]",filename.c_str());
   return false;
   }
   
   // read multipole order
   // 4 lines: skip one, read one, skip two
   fgets(line,1024,fp);
   if (fgets(line,1024,fp) != NULL) {
   unsigned int index;
   unsigned int file_order;
   gmp_sscanf(line,"OBF(%i) = %d,",&index,&file_order);
   _order.set(file_order);
   ORSA_DEBUG("read order: %i",file_order);
   } else {
   ORSA_ERROR("problems while reading file [%s]",filename.c_str());
   return false;
   }
   //
   fgets(line,1024,fp); 
   fgets(line,1024,fp);
   
   // GM
   if (fgets(line,1024,fp) != NULL) {
   unsigned int index;
   double GM;
   gmp_sscanf(line,"GM(%d) =   %f,",&index,GM());
   GM = FromUnits(FromUnits(GM,Unit::KM,3),Unit::SECOND,-2);
   //
   ORSA_DEBUG("M: %Fe kg",FromUnits(GM/Unit::instance()->getG(),Unit::KG,-1));
   }
   
   const unsigned int _order_plus_one = _order.getRef()+1;
   //
   {
   _C.resize(_order_plus_one);
   _C_uncertainty.resize(_order_plus_one);
   _S.resize(_order_plus_one);
   _S_uncertainty.resize(_order_plus_one);
   for (unsigned int l=0; l<_order_plus_one; ++l) {
   _C[l].resize(l+1);
   _C_uncertainty[l].resize(l+1);
   _S[l].resize(l+1);
   _S_uncertainty[l].resize(l+1);
   }
   }
   //
   {
   // set all to zero...
   // important, because not all the values are in the file;
   for (unsigned int l=0; l<_order_plus_one; ++l) {
   for (unsigned int m=0; m<=l; ++m) {
   _C[l][m] = 0;
   _S[l][m] = 0;
   }
   }
   // set C[0][0] to one
   _C[0][0] = 1;
   }  
   
   while (1) {
   unsigned int l;
   char valChar[1024];
   if (fgets(line,1024,fp) != NULL) {
   // ORSA_DEBUG("line: %s",line);
   if (gmp_sscanf(line,"OBAJ(%d) = %s",&l,valChar)) {
   for (unsigned int k=0;k<strlen(valChar);++k) {
   if (valChar[k] == 'D') valChar[k] = 'e';
   }
   // qremove trailing comma...
   if (strlen(valChar) > 0) valChar[strlen(valChar)-1] = '\0';
   //
   // ORSA_DEBUG("valChar: %s",valChar);
   //
   _C[l][0] = -double(valChar);
   // ORSA_DEBUG("C[%i][0]: %Fe",l,_C[l][0]());
   }
   }
   if (l==_order.getRef()) {
   // ORSA_DEBUG("out...");
   break;
   }
   }
   
   while (1) {
   unsigned int lC, mC, lS, mS;
   char valCharC[1024], valCharS[1024];
   if (fgets(line,1024,fp) != NULL) {
   if (strlen(line) < 20) break;
   if (gmp_sscanf(line,"OBAC(%d,%d) = %s OBAS(%d,%d) = %s",
   &lC, &mC, valCharC,
   &lS, &mS, valCharS)) {
   
   {
   for (unsigned int k=0;k<strlen(valCharC);++k) {
   if (valCharC[k] == 'D') valCharC[k] = 'e';
   }
   // qremove trailing comma...
   if (strlen(valCharC) > 0) valCharC[strlen(valCharC)-1] = '\0';
   //
   _C[lC][mC] = double(valCharC);
   // ORSA_DEBUG("C[%i][%i]: %Fe",lC,mC,_C[lC][mC]());
   }
   //
   {
   for (unsigned int k=0;k<strlen(valCharS);++k) {
   if (valCharS[k] == 'D') valCharS[k] = 'e';
   }
   // qremove trailing comma...
   if (strlen(valCharS) > 0) valCharS[strlen(valCharS)-1] = '\0';
   //
   _S[lS][mS] = double(valCharS);
   // ORSA_DEBUG("S[%i][%i]: %Fe",lS,mS,_S[lS][mS]());
   }
   }
   }
   if (lC==mC==_order.getRef()) break;
   }
   
   fclose(fp);
   return true;
   }
*/
