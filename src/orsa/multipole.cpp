#include <orsa/box.h>
#include <orsa/double.h>
#include <orsa/multipole.h>
#include <orsa/legendre.h>
#include <orsa/statistic.h>
#include <orsa/unit.h>

#include <vector>

#include <gsl/gsl_rng.h>

#include <iostream>

// debug
// #include <iostream>

using namespace orsa;

Multipole::Multipole() : osg::Referenced(true) {
  _shape = 0;
}

Multipole::Multipole(const unsigned int order) : osg::Referenced(true) {
  _shape = 0;
  _massDistribution = new UniformMassDistribution;
  _order.set(order);
  
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
}

Multipole::~Multipole() { 
  
}

/* 
   static Vector __randomVectorUtil(gsl_rng * rnd,
   const double & boundingRadius) {
   return Vector(boundingRadius*(two()*gsl_rng_uniform(rnd)-1),
   boundingRadius*(two()*gsl_rng_uniform(rnd)-1),
   boundingRadius*(two()*gsl_rng_uniform(rnd)-1));
   }
*/

static Vector __randomVectorUtil(gsl_rng * rnd,
				 const Box & boundingBox) {
  return Vector(boundingBox.getXMin()+(boundingBox.getXMax()-boundingBox.getXMin())*gsl_rng_uniform(rnd),
		boundingBox.getYMin()+(boundingBox.getYMax()-boundingBox.getYMin())*gsl_rng_uniform(rnd),
		boundingBox.getZMin()+(boundingBox.getZMax()-boundingBox.getZMin())*gsl_rng_uniform(rnd));
}

void Multipole::_randomVectorInside(std::vector<bool> & in,
				    const unsigned int sample_points,
				    const int random_seed) {
  
  ORSA_DEBUG("called _randomVectorInside(...): %i sample points, random seed: %i",
	     sample_points,random_seed);
  
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
  
  if (1) {
    // debug output
    unsigned int _inside = 0;
    for (unsigned int j=0; j<sample_points; ++j) {
      if (in[j]) ++_inside;
    }
    ORSA_DEBUG("inside: %i (%F5.2f\%)",_inside,(100.0*_inside)/sample_points);
  }
  
  ORSA_DEBUG("done with _randomVectorInside(...).");
}

void Multipole::_centerOfMass(Vector & center_of_mass,
			      Vector & center_of_mass_uncertainty,
			      const std::vector<bool> & in,
			      const unsigned int sample_points,
			      const int random_seed) {
  
  ORSA_DEBUG("called _centerOfMass(...): %i sample points, random seed: %i",
	     sample_points,random_seed);
  
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
    // const Vector & _v = __randomVectorUtil(rnd,boundingRadius);
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
  
  if (1) {
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
  
  ORSA_DEBUG("done with _centerOfMass(...).");
}

void Multipole::_inertiaMoment(Matrix & inertia_moment,
			       Matrix & inertia_moment_uncertainty,
			       const Vector & center_of_mass,
			       const Vector & center_of_mass_uncertainty,
			       const std::vector<bool> & in,
			       const unsigned int sample_points,
			       const int random_seed) {
  
  ORSA_DEBUG("called _inertiaMoment(...): %i sample points, random seed: %i",
	     sample_points,random_seed);
  
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
  
  if (1) {
    // debug output    
    ORSA_DEBUG("Ixx:  %14.6e +/- %14.6e",
	       _stat_Ixx->average(),
	       _stat_Ixx->averageError());
    ORSA_DEBUG("Iyy:  %14.6e +/- %14.6e",
	       _stat_Iyy->average(),
	       _stat_Iyy->averageError());
    ORSA_DEBUG("Izz:  %14.6e +/- %14.6e",
	       _stat_Izz->average(),
	       _stat_Izz->averageError());
    //
    ORSA_DEBUG("Ixy:  %14.6e +/- %14.6e",
	       _stat_Ixy->average(),
	       _stat_Ixy->averageError());
    ORSA_DEBUG("Ixz:  %14.6e +/- %14.6e",
	       _stat_Ixz->average(),
	       _stat_Ixz->averageError());
    ORSA_DEBUG("Iyz:  %14.6e +/- %14.6e",
	       _stat_Iyz->average(),
	       _stat_Iyz->averageError());
  }
  
  ORSA_DEBUG("done with _inertiaMoment(...).");
}

void Multipole::_multipole(std::vector<std::vector<double> > & C,
			   std::vector<std::vector<double> > & C_uncertainty,
			   std::vector<std::vector<double> > & S, 
			   std::vector<std::vector<double> > & S_uncertainty, 
			   const Vector & center_of_mass,
			   const Vector & center_of_mass_uncertainty,
			   const std::vector<bool> & in,
			   const unsigned int sample_points,
			   const int random_seed) {
  
  ORSA_DEBUG("called _multipole(...): %i sample points, random seed: %i",
	     sample_points,random_seed);
  
  // GSL rng init
  gsl_rng * rnd = gsl_rng_alloc(gsl_rng_gfsr4);
  gsl_rng_set(rnd,random_seed);
  
  const unsigned int _order_plus_one = _order.getRef()+1;
  //
  osg::ref_ptr<orsa::Statistic<double> > _stat_M = new orsa::Statistic<double>;
  //
  std::vector<std::vector<osg::ref_ptr<orsa::Statistic<double> > > > _stat_C(_order_plus_one);
  std::vector<std::vector<osg::ref_ptr<orsa::Statistic<double> > > > _stat_S(_order_plus_one);
  //
  for (unsigned int p=0; p<(_order_plus_one); ++p) {
    _stat_C[p].resize(p+1);
    _stat_S[p].resize(p+1);
    //
    for (unsigned int q=0; q<=p; ++q) {
      _stat_C[p][q] = new orsa::Statistic<double>;
      _stat_S[p][q] = new orsa::Statistic<double>;
    }
  }
  //
  // const double _R_max = _shape->boundingRadius();
  //
  std::vector<double> _R(_order_plus_one);
  _R[0] = 1;
  for (unsigned int p=1; p<_R.size(); ++p) {
    _R[p] = _R[p-1]*_R0.getRef();
  }
  //    
  std::vector<double> _r(_order_plus_one);
  //
  std::vector<double> _vec_c_m_phi(_order_plus_one), _vec_s_m_phi(_order_plus_one);
  //
  // const double boundingRadius = _shape->boundingRadius();
  const Box boundingBox = _shape->boundingBox();
  unsigned int _iter = 0;
  while (_iter < sample_points) {
    // correct _v for the center-of-mass position
    // const Vector & _v = __randomVectorUtil(rnd,boundingRadius) - center_of_mass;
    const Vector & _v = __randomVectorUtil(rnd,boundingBox) - center_of_mass;
    if (in[_iter]) {
      const Vector _unit_v = _v.normalized();
      const double _c_theta = _unit_v.getZ();
      const double _phi = atan2(_unit_v.getY(),
				_unit_v.getX());
      for (unsigned int m=0; m<_order_plus_one; ++m) {
	_vec_c_m_phi[m] = cos(m*_phi);
	_vec_s_m_phi[m] = sin(m*_phi);
	//
	// ORSA_DEBUG("m: %i  vec_s_m_phi[%i] =  %Fg",m,m,_vec_s_m_phi[m]);
      }
      //
      Legendre _L(_c_theta,_order.getRef());
      //
      const double _v_length = _v.length();
      //
      _r[0] = 1;
      for (unsigned int p=1; p<_order_plus_one; ++p) {
	_r[p] = _r[p-1]*_v_length;
      }
      //
      const double density = _massDistribution->density(_v + center_of_mass);
      _stat_M->insert(density);
      //
      for (unsigned int l=0; l<=_order.getRef(); ++l) {
	for (unsigned int m=0; m<=l; ++m) {
	  _stat_C[l][m]->insert(density*(_r[l]/_R[l])*_L.P(l,m)*_vec_c_m_phi[m]);
	  _stat_S[l][m]->insert(density*(_r[l]/_R[l])*_L.P(l,m)*_vec_s_m_phi[m]);
	}
      }
    }
    ++_iter;
  }
  
  ORSA_DEBUG("IMPORTANT: complete density propagation to C_uncertainty and S_uncertainty");
  
  C.resize(_order_plus_one);
  C_uncertainty.resize(_order_plus_one);
  S.resize(_order_plus_one);
  S_uncertainty.resize(_order_plus_one);
  for (unsigned int p=0; p<_order_plus_one; ++p) {
    C[p].resize(p+1);
    C_uncertainty[p].resize(p+1);
    S[p].resize(p+1);
    S_uncertainty[p].resize(p+1);
  }
  //
  // const double _test_norm = 1/(Legendre::norm(0,0)*_stat_C[0][0]->average());
  const double _test_norm = _stat_M->sum()/(Legendre::norm(0,0)*_stat_C[0][0]->sum());
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
  //
  for (unsigned int l=0; l<=_order.getRef(); ++l) {
    for (unsigned int m=0; m<=l; ++m) {
      C[l][m]             = _test_norm*Legendre::norm(l,m)*_stat_C[l][m]->sum()/_stat_M->sum();
      C_uncertainty[l][m] = _test_norm*Legendre::norm(l,m)*_stat_C[l][m]->averageError();
      S[l][m]             = _test_norm*Legendre::norm(l,m)*_stat_S[l][m]->sum()/_stat_M->sum();
      S_uncertainty[l][m] = _test_norm*Legendre::norm(l,m)*_stat_S[l][m]->averageError();
    }
  }
  
  // GSL rng clean
  gsl_rng_free(rnd);
  
  if (1) {
    // debug output
    const double _test_norm = 1/(Legendre::norm(0,0)*_stat_C[0][0]->average());
    ORSA_DEBUG("test_norm: %Fg",_test_norm);
    for (unsigned int l=0; l<=_order.getRef(); ++l) {
      for (unsigned int m=0; m<=l; ++m) {
	ORSA_DEBUG("C[%02i][%02i]: %14.6e +/- %14.6e   rln: %Fg",
		   l,m,
		   _test_norm*Legendre::norm(l,m)*_stat_C[l][m]->average(),
		   _test_norm*Legendre::norm(l,m)*_stat_C[l][m]->averageError(),
		   Legendre::norm(l,m));
	ORSA_DEBUG("S[%02i][%02i]: %14.6e +/- %14.6e   rln: %Fg",
		   l,m,
		   _test_norm*Legendre::norm(l,m)*_stat_S[l][m]->average(),
		   _test_norm*Legendre::norm(l,m)*_stat_S[l][m]->averageError(),
		   Legendre::norm(l,m));
	//
	// std::cerr << "norm: " << Legendre::norm(l,m) << std::endl;
      }
    }
  }
  
  ORSA_DEBUG("done with _multipole(...).");
}

// bool Multipole::compute() {
//
bool Multipole::computeUsingShape(const unsigned int sample_points,
				  const int random_seed) {
  if (_shape==0) {
    ORSA_ERROR("unset Shape");
    return false;
  }
  
  return (computeUsingShape(sample_points,
			    random_seed,
			    _shape->boundingRadius()));
}

bool Multipole::computeUsingShape(const unsigned int sample_points,
				  const int random_seed,
				  const double & R0) {
  if (_shape==0) {
    ORSA_ERROR("unset Shape");
    return false;
  }
  
  if (!_order.isSet()) {
    ORSA_ERROR("unset order");
    return false;
  }
  
  /* 
     if (R0 >= _shape->boundingRadius()) {
     _R0.set(R0);
     } else {
     ORSA_ERROR("R0=%f KM is smaler than the boundingRadius()=%f KM",
     FromUnits(R0,Unit::KM,-1), 
     FromUnits(_shape->boundingRadius(),Unit::KM,-1));
     return false;
     }
  */
  //
  if (R0 > 0) {
    ORSA_DEBUG("setting R0=%Fg",R0);
    _R0.set(R0);
    ORSA_DEBUG("_R0=%Fg",_R0.getRef());
  } else {
    ORSA_ERROR("cannot set negative R0, using Shape::boundingRadius()");
    _R0.set(_shape->boundingRadius());
    ORSA_DEBUG("_R0=%Fg",_R0.getRef());
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
  
  // center-of-mass
  _centerOfMass(_center_of_mass,
		_center_of_mass_uncertainty,
		_in,
		sample_points,
		random_seed);
  
  // inertia moment
  /* 
     Matrix _inertia_moment, _inertia_moment_uncertainty;
     _inertiaMoment(_inertia_moment,
     _inertia_moment_uncertainty,
     _center_of_mass,
     _center_of_mass_uncertainty,
     _in,
     sample_points,
     random_seed);
  */
  
  // multipole
  _multipole(_C,
	     _C_uncertainty,
	     _S, 
	     _S_uncertainty, 
	     _center_of_mass,
	     _center_of_mass_uncertainty,
	     _in,
	     sample_points,
	     random_seed);
  
  return true;
}

/* 
   bool Multipole::writeToFile(const std::string & filename) const {
   ORSA_DEBUG("saving multipole to file [%s]",filename.c_str());
   FILE * fp = fopen(filename.c_str(),"w");
   if (!fp) {
   ORSA_WARNING("cannot open file %s",filename.c_str());
   return false;
   }
   fprintf(fp,"%i\n",_order.getRef());
   for (unsigned int l=0; l<_C.size(); ++l) {
   for (unsigned int m=0; m<=l; ++m) {
   // gmp_fprintf(fp,"%16.8e %16.8e\n",
   gmp_fprintf(fp,"%18.10e %18.10e\n",
   _C[l][m],
   _S[l][m]);
   }
   }
   fclose(fp);
   return true;
   }
*/

/* 
   bool Multipole::readFromFile(const std::string & filename) {
   ORSA_DEBUG("reading multipole from file [%s]",filename.c_str());
   FILE * fp = fopen(filename.c_str(),"r");
   if (!fp) {
   ORSA_WARNING("cannot open file %s",filename.c_str());
   return false;
   }
   unsigned int _n_order;
   gmp_fscanf(fp,"%i",&_n_order);
   _order.set(_n_order);
   ORSA_DEBUG("read order: %i",_n_order);
   //
   _C.resize(_n_order+1);
   _S.resize(_n_order+1);
   for (unsigned int l=0; l<_n_order+1; ++l) {
   _C[l].resize(l+1);
   _S[l].resize(l+1);
   }
   //
   double _c, _s;
   for (unsigned int l=0; l<_n_order+1; ++l) {
   for (unsigned int m=0; m<=l; ++m) {
   // gmp_fprintf(fp,"%16.8e %16.8e\n",
   gmp_fscanf(fp,"%e %e",
   _c(),
   _s()); // '&' not required for GMP vars
   _C[l][m] = _c;
   _S[l][m] = _s;
   //
   ORSA_DEBUG("read C[%i][%i]: %18.10e   S[%i][%i]: %18.10e",
   l,m,_c(),
   l,m,_s());
   }
   }
   fclose(fp);
   return true;
   }
*/

bool Multipole::writeToFile(const std::string & filename) const {
  // ORSA_DEBUG("saving multipole to file [%s]",filename.c_str());
  FILE * fp = fopen(filename.c_str(),"w");
  if (!fp) {
    ORSA_WARNING("cannot open file %s",filename.c_str());
    return false;
  }
  gmp_fprintf(fp,"%18.10e\n",FromUnits(_R0.getRef(),Unit::M,-1));
  gmp_fprintf(fp,"%i\n",_order.getRef());
  gmp_fprintf(fp,"%i\n",_sample_points.getRef());
  //
  gmp_fprintf(fp,"%18.10e %18.10e\n",
	      FromUnits(_center_of_mass.getX(),Unit::M,-1),
	      FromUnits(_center_of_mass_uncertainty.getX(),Unit::M,-1));
  gmp_fprintf(fp,"%18.10e %18.10e\n",
	      FromUnits(_center_of_mass.getY(),Unit::M,-1),
	      FromUnits(_center_of_mass_uncertainty.getY(),Unit::M,-1));
  gmp_fprintf(fp,"%18.10e %18.10e\n",
	      FromUnits(_center_of_mass.getZ(),Unit::M,-1),
	      FromUnits(_center_of_mass_uncertainty.getZ(),Unit::M,-1));
  //
  for (unsigned int l=0; l<_C.size(); ++l) {
    for (unsigned int m=0; m<=l; ++m) {
      gmp_fprintf(fp,
		  "%18.10e %18.10e\n",
		  _C[l][m],
		  _C_uncertainty[l][m]);
      gmp_fprintf(fp,
		  "%18.10e %18.10e\n",
		  _S[l][m],
		  _S_uncertainty[l][m]);
      //
      /* 
	 gmp_fprintf(fp,
	 "%18.10e %18.10e  C %i %i %18.10e %18.10e %Fg\n",
	 _C[l][m],
	 _C_uncertainty[l][m],
	 l,
	 m,
	 double(_C[l][m]*Legendre::norm(l,m)),
	 double(_C[l][m]/Legendre::norm(l,m)),
	 double(Legendre::norm(l,m)));
	 gmp_fprintf(fp,
	 "%18.10e %18.10e  S %i %i %18.10e %18.10e %Fg\n",
	 _S[l][m],
	 _S_uncertainty[l][m],
	 l,
	 m,
	 double(_S[l][m]*Legendre::norm(l,m)),
	 double(_S[l][m]/Legendre::norm(l,m)),
	 double(Legendre::norm(l,m)));
      */
    }
  }
  fclose(fp);
  return true;
}

bool Multipole::readFromFile(const std::string & filename) {
  // ORSA_DEBUG("trying to read multipole from file [%s]",filename.c_str());
  FILE * fp = fopen(filename.c_str(),"r");
  if (!fp) {
    ORSA_WARNING("cannot open file %s",filename.c_str());
    return false;
  }
  {
    double R0;
    gmp_fscanf(fp,"%lf",
	       &R0); // '&' not required for GMP vars
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
    gmp_fscanf(fp,"%lf %lf",
	       &x,
	       &dx);
    gmp_fscanf(fp,"%lf %lf",
	       &y,
	       &dy);
    gmp_fscanf(fp,"%lf %lf",
	       &z,
	       &dz);
    //
    x  = FromUnits(x,Unit::M);
    dx = FromUnits(x,Unit::M);
    y  = FromUnits(x,Unit::M);
    dy = FromUnits(x,Unit::M);
    z  = FromUnits(x,Unit::M);
    dz = FromUnits(x,Unit::M);
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
	gmp_fscanf(fp,"%lf %lf",
		   &_c,
		   &_dc); // '&' not required for GMP vars
	_C[l][m] = _c;
	_C_uncertainty[l][m] = _dc;
	//
	gmp_fscanf(fp,"%lf %lf",
		   &_s,
		   &_ds); // '&' not required for GMP vars
	_S[l][m] = _s;
	_S_uncertainty[l][m] = _ds;
	//
	/* 
	   ORSA_DEBUG("read C[%i][%i]: %18.10e   S[%i][%i]: %18.10e",
	   l,m,_c(),
	   l,m,_s());
	*/
      }
    }
  }
  
  fclose(fp);
  return true;
}

bool Multipole::_write_in_vector_to_file(const std::vector<bool> & in,
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

bool Multipole::_read_in_vector_from_file(std::vector<bool> & in,
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

bool Multipole::readFromAltFile(const std::string & filename) {
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
    gmp_sscanf(line,"%*s %lf",
	       &R0); // '&' not required for GMP vars
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
      gmp_fscanf(fp,"%lf",&file_entry);
      // ORSA_DEBUG("file_entry: %Fg",file_entry());      
      _C[l][m]             = file_entry;
      _C_uncertainty[l][m] = 0;
      //
      // ORSA_DEBUG("C[%02i][%02i] = %16.12f",l,m,_C[l][m]);      
      
      // S
      if (m==0) {
	_S[l][m]             = 0;
	_S_uncertainty[l][m] = 0;
      } else {
	gmp_fscanf(fp,"%lf",&file_entry);
	// ORSA_DEBUG("file_entry: %Fg",file_entry());     
	//
	_S[l][m]             = file_entry;
	_S_uncertainty[l][m] = 0;
      }
      //
      // ORSA_DEBUG("S[%02i][%02i] = %16.12f",l,m,_S[l][m]);      
    }
  }
  
  fclose(fp);
  return true;
}

bool Multipole::readFromBisFile(const std::string & filename) {
  
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
    gmp_sscanf(line,"OBR = %lf,",
	       &R0); // '&' not required for GMP vars
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
    ORSA_DEBUG("read order: %i",&file_order);
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
    gmp_sscanf(line,"GM(%d) =   %lf,",&index,&GM);
    GM = FromUnits(FromUnits(GM,Unit::KM,3),Unit::SECOND,-2);
    //
    ORSA_DEBUG("M: %e kg",FromUnits(GM/Unit::G(),Unit::KG,-1));
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
	_C[l][0] = -atof(valChar);
	// ORSA_DEBUG("C[%i][0]: %e",l,_C[l][0]);
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
	  _C[lC][mC] = atof(valCharC);
	  // ORSA_DEBUG("C[%i][%i]: %e",lC,mC,_C[lC][mC]);
	}
	//
	{
	  for (unsigned int k=0;k<strlen(valCharS);++k) {
	    if (valCharS[k] == 'D') valCharS[k] = 'e';
	  }
	  // qremove trailing comma...
	  if (strlen(valCharS) > 0) valCharS[strlen(valCharS)-1] = '\0';
	  //
	  _S[lS][mS] = atof(valCharS);
	  // ORSA_DEBUG("S[%i][%i]: %e",lS,mS,_S[lS][mS]);
	}
      }
    }
    if (lC==mC==_order.getRef()) break;
  }
  
  fclose(fp);
  return true;
}
