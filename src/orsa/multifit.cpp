#include <orsa/multifit.h>

using namespace orsa;

// MultifitGlue

MultifitGlue * MultifitGlue::_instance = 0;

MultifitGlue * MultifitGlue::instance() {
  if (_instance == 0) {
    _instance = new MultifitGlue;
  }
  return _instance;
}

MultifitGlue::MultifitGlue() {
  _mf = 0;
}

MultifitGlue::~MultifitGlue() {
  _instance = 0;
}

// gsl functions

int orsa::multifit_global_f_gsl (const gsl_vector * v, 
				 void * dataPoints, 
				 gsl_vector * f) {
  return MultifitGlue::instance()->_mf->f_gsl(v,dataPoints,f);
}

int orsa::multifit_global_df_gsl (const gsl_vector * v, 
				  void * dataPoints, 
				  gsl_matrix * J) {
  return MultifitGlue::instance()->_mf->df_gsl(v,dataPoints,J);
}

int orsa::multifit_global_fdf_gsl (const gsl_vector * v, 
				   void * dataPoints, 
				   gsl_vector * f, 
				   gsl_matrix * J) {
  return MultifitGlue::instance()->_mf->fdf_gsl(v,dataPoints,f,J);
}

// MultifitParameters

MultifitParameters::MultifitParameters() : Referenced(true) { } 

MultifitParameters::~MultifitParameters() { }

bool MultifitParameters::insert(const  std::string & name,
				const double & initialValue,
				const double & delta) {
  dataType::const_iterator it = _data.begin();
  while (it != _data.end()) {
    if ((*it).name == name) {
      ORSA_ERROR("trying to insert duplicate variable");
      return false;
    } 
    ++it;
  }
  {
    NVD tmp;
    //
    tmp.name  = name;
    tmp.value = initialValue;
    tmp.delta = delta;
    //
    _data.push_back(tmp);
  }
  return true;
}

void MultifitParameters::clear() {
  _data.clear();
}

unsigned int MultifitParameters::index(const std::string & name) const {
  for (unsigned int k=0; k<_data.size(); ++k) {
    if (_data[k].name == name) {
      return k;
    }	
  }
  ORSA_ERROR("index not found: name=%s",name.c_str());
  return 0;
}

std::string MultifitParameters::name(const unsigned int index) const {
  if (_data.size() <= index) {
    ORSA_ERROR("name not found");
    return "";
  }
  return _data[index].name;
}

bool MultifitParameters::set(const std::string & name,
			     const double      & value) {
  return set(MultifitParameters::index(name),value);
}

bool MultifitParameters::set(const unsigned int   index,
			     const double       & value) {
  _data[index].value = value;
  setInRange(index);
  return true;
}

bool MultifitParameters::setRange(const std::string  & name,
				  const double       & min,
				  const double       & max) {
  return setRange(MultifitParameters::index(name),
		  min,
		  max); 
}

bool MultifitParameters::setRange(const unsigned int   index,
				  const double       & min,
				  const double       & max) {
  if (min > max) {
    ORSA_DEBUG("problems: min > max");
    return false;
  }
  _data[index].min = min;
  _data[index].max = max;
  setInRange(index);
  return true;
}

bool MultifitParameters::setRangeMin(const std::string & name,
				     const double      & min) {
  return setRangeMin(MultifitParameters::index(name),min);
}

bool MultifitParameters::setRangeMax(const std::string & name,
				     const double      & max) {
  return setRangeMax(MultifitParameters::index(name),max);
}

bool MultifitParameters::setRangeMin(const unsigned int   index,
				     const double       & min) {
  _data[index].min = min;
  return true;
}

bool MultifitParameters::setRangeMax(const unsigned int   index,
				     const double       & max) {
  _data[index].max = max;
  return true;
}

void MultifitParameters::setInRange(const std::string & name) {
  setInRange(MultifitParameters::index(name));
}

void MultifitParameters::setInRange(const unsigned int index) {
  if (_data[index].min.isSet()) {
    if (_data[index].value < _data[index].min.getRef()) {
      _data[index].value = _data[index].min.getRef();
    }
  }
  if (_data[index].max.isSet()) {
    if (_data[index].value > _data[index].max.getRef()) {
      _data[index].value = _data[index].max.getRef();
    }
  }
}

double MultifitParameters::get(const std::string & name) const {
  return get(MultifitParameters::index(name));
}

double MultifitParameters::get(const unsigned int index) const {
  return _data[index].value;
}

/* bool MultifitParameters::haveRange(const std::string & name) const {
   return haveRange(MultifitParameters::index(name));
   }
   
   bool MultifitParameters::haveRange(const unsigned int index) const {
   return (_data[index].min.isSet() || _data[index].max.isSet());
   }
*/

const orsa::Cache<double> & MultifitParameters::getRangeMin(const std::string & name) const {
  return getRangeMin(MultifitParameters::index(name));
}

const orsa::Cache<double> & MultifitParameters::getRangeMax(const std::string & name) const {
  return getRangeMax(MultifitParameters::index(name));
}

const orsa::Cache<double> & MultifitParameters::getRangeMin(const unsigned int index) const {
  return (_data[index].min);
}

const orsa::Cache<double> & MultifitParameters::getRangeMax(const unsigned int index) const {
  return (_data[index].max);
}

bool MultifitParameters::setDelta(const std::string  & name,
				  const double & delta) {
  _data[MultifitParameters::index(name)].delta = delta;
  return true;
}

bool MultifitParameters::setDelta(const unsigned int   index,
				  const double & delta) {
  _data[index].delta = delta;
  return true;
}

double MultifitParameters::getDelta(const std::string & name) const {
  return _data[MultifitParameters::index(name)].delta;
}

double MultifitParameters::getDelta(const unsigned int index) const {
  return _data[index].delta;
}

bool MultifitParameters::writeToFile(const std::string & fileName) const {
  FILE * fp = fopen(fileName.c_str(),"w");
  if (fp == 0) {
    ORSA_ERROR("cannot open file %s",fileName.c_str());
    return false;
  }
  for (unsigned int k=0; k<size(); ++k) {
    gmp_fprintf(fp,"%16s %+22.16e %+22.16e\n",
		name(k).c_str(),
		get(k),
		getDelta(k));
  } 
  fclose(fp);
  return true;
}

bool MultifitParameters::readFromFile(const std::string & fileName) {
  clear();
  FILE * fp = fopen(fileName.c_str(),"r");
  if (fp == 0) {
    ORSA_ERROR("cannot open file %s",fileName.c_str());
    return false;
  }
  char name[1024];
  double value, delta;
  while (gmp_fscanf(fp,"%s %lf %lf\n",
		    name,
		    &value,
		    &delta) == 3) {
    insert(name,value,delta);
  }	
  fclose(fp);
  return true;
}

bool orsa::operator == (const MultifitParameters & p1, 
			const MultifitParameters & p2) {
  
  if (0) {
    // debug output 
    for (unsigned int k=0; k<p1.size(); ++k) {
      ORSA_DEBUG("p1[%02i] = [%20s] = %18.8g",
		 k,
		 p1.name(k).c_str(),
		 p1.get(k));
    }
    for (unsigned int k=0; k<p2.size(); ++k) {
      ORSA_DEBUG("p2[%02i] = [%20s] = %18.8g",
		 k,
		 p2.name(k).c_str(),
		 p2.get(k));
    }
  }
  
  if (p1.size() != p2.size()) {
    return false;
  }
  
  for (unsigned int k=0; k<p1.size(); ++k) {
    if (p1.name(k) != p2.name(k)) {
      return false;
    } 
  }
  
  /* 
     for (unsigned int k=0; k<p1.size(); ++k) {
     if (p1.get(k) != p2.get(k)) {
     return false;
     } 
     }
  */
  
  for (unsigned int k=0; k<p1.size(); ++k) {
    if (p1.get(k) != p2.get(k)) {
      return false;
    } else {
      /* 
	 ORSA_DEBUG("EQUAL: %f == %f",
	 p1.get(k),
	 p2.get(k));
      */
    }
  }
  
  return true;
}

bool orsa::operator != (const MultifitParameters & p1, 
			const MultifitParameters & p2) {
  return (!(p1 == p2));
}

// MultifitData

MultifitData::MultifitData() : Referenced(true) { } 

MultifitData::~MultifitData() { }

bool MultifitData::insertVariable(const std::string & name) {
  std::vector<NZDV>::const_iterator it = _data.var.begin();
  while (it != _data.var.end()) {
    if ((*it).name.getRef() == name) {
      ORSA_ERROR("trying to insert duplicate variable");
      return false;
    }
    ++it;
  }
  {
    NZDV tmp;
    //
    tmp.name.set(name);
    //
    _data.var.push_back(tmp);
  }
  return true;
}

bool MultifitData::insertZ(const std::string  & name,
			   const unsigned int   row,
			   const mpz_class    & value) {
  return insertZ(MultifitData::index(name),row,value);
}

bool MultifitData::insertZ(const unsigned int   index,
			   const unsigned int   row,
			   const mpz_class    & value) {
  if (_data.var.size() <= index) {
    ORSA_ERROR("index bigger than container size");
    return false;
  }
  if (_data.var[index].z.size() <= row) {
    _data.var[index].z.resize(row+1);
  }
  _data.var[index].z[row].set(value);
  return true;
}

bool MultifitData::insertD(const std::string  & name,
			   const unsigned int   row,
			   const double & value) {
  return insertD(MultifitData::index(name),row,value);
}

bool MultifitData::insertD(const unsigned int   index,
			   const unsigned int   row,
			   const double & value) {
  if (_data.var.size() <= index) {
    ORSA_ERROR("index bigger than container size");
    return false;
  }
  if (_data.var[index].d.size() <= row) {
    _data.var[index].d.resize(row+1);
  }
  _data.var[index].d[row].set(value);
  return true;
}

bool MultifitData::insertV(const std::string  & name,
			   const unsigned int   row,
			   const orsa::Vector & value) {
  return insertV(MultifitData::index(name),row,value);
}

bool MultifitData::insertV(const unsigned int   index,
			   const unsigned int   row,
			   const orsa::Vector & value) {
  if (_data.var.size() <= index) {
    ORSA_ERROR("index bigger than container size");
    return false;
  }
  if (_data.var[index].v.size() <= row) {
    _data.var[index].v.resize(row+1);
  }
  _data.var[index].v[row].set(value);
  return true;
}

bool MultifitData::insertF(const unsigned int   row,
			   const double & value) {
  if (_data.f.size() <= row) {
    _data.f.resize(row+1);
  }
  _data.f[row].set(value);
  return true;
}

bool MultifitData::insertSigma(const unsigned int   row,
			       const double & value) {
  if (_data.sigma.size() <= row) {
    _data.sigma.resize(row+1);
  }
  _data.sigma[row].set(value);
  return true;
}

unsigned int MultifitData::index(const std::string & name) const {
  // ORSA_DEBUG("called name: %s",name.c_str());
  for (unsigned int k=0; k<_data.var.size(); ++k) {
    // ORSA_DEBUG("tested name: %s",_data.var[k].name.getRef().c_str());
    if (_data.var[k].name.getRef() == name) {
      // ORSA_DEBUG("found: returning %i",k);
      return k;
    }	
  }
  ORSA_ERROR("index not found");
  return 0;
}

std::string MultifitData::name(const unsigned int index) const {
  if (_data.var.size() <= index) {
    ORSA_ERROR("name not found: %i",index);
    return "";
  }
  return _data.var[index].name.getRef();
}

mpz_class MultifitData::getZ(const std::string & name,
			     const unsigned int  row) const {
  return getZ(MultifitData::index(name),row);
}

mpz_class MultifitData::getZ(const unsigned int index,
			     const unsigned int row) const {
  if (_data.var[index].z.size() <= row) {
    ORSA_ERROR("row too big, var: %s row=%i",name(index).c_str(),row);
    return 0;
  }
  return _data.var[index].z[row].getRef();
}

double MultifitData::getD(const std::string & name,
			  const unsigned int  row) const {
  return getD(MultifitData::index(name),row);
}

double MultifitData::getD(const unsigned int index,
			  const unsigned int row) const {
  if (_data.var[index].d.size() <= row) {
    ORSA_ERROR("row too big");
    return 0;
  }
  return _data.var[index].d[row].getRef();
}

Vector MultifitData::getV(const std::string & name,
			  const unsigned int  row) const {
  return getV(MultifitData::index(name),row);
}

Vector MultifitData::getV(const unsigned int index,
			  const unsigned int row) const {
  if (_data.var[index].v.size() <= row) {
    ORSA_ERROR("row too big");
    return Vector();
  }
  return _data.var[index].v[row].getRef();
}

double MultifitData::getF(const unsigned int row) const {
  // ORSA_DEBUG("v: %i",_data.var.size());
  // ORSA_DEBUG("f: %i",_data.f.size());
  // ORSA_DEBUG("s: %i",_data.sigma.size());
  if (_data.f.size() <= row) {
    ORSA_ERROR("row too big");
    return 0;
  }
  return _data.f[row].getRef();
}

double MultifitData::getSigma(const unsigned int row) const {
  if (_data.sigma.size() <= row) {
    ORSA_ERROR("row too big");
    return 0;
  }
  return _data.sigma[row].getRef();
}

unsigned int MultifitData::size() const {
  // a bit naive, assumes that the data is not "corrupted"...
  return _data.f.size();
}

unsigned int MultifitData::vars() const {
  return _data.var.size();
}

// Multifit

Multifit::Multifit() : Referenced(true) { }

Multifit::~Multifit() { }

double Multifit::__fun__(const orsa::MultifitParameters * par, 
			 const orsa::MultifitData       * data,
			 const unsigned int p, // par index
			 const int          d, // delta
			 const unsigned int row) const {
  
  double effectiveFun = fun(par,data,p,d,row);
  
  bool outsideRange=false;
  //
  if (par->getRangeMin(p).isSet()) {
    if (par->get(p) < par->getRangeMin(p).getRef()) {
      outsideRange=true;
    }    
  }
  // 
  if (par->getRangeMax(p).isSet()) {
    if (par->get(p) > par->getRangeMax(p).getRef()) {
        outsideRange=true;
    }    
  }
  // 
  if (outsideRange) {
    effectiveFun += 1.0e3 * (1.0 + fabs(effectiveFun));
  }
  
  return effectiveFun;
}

int Multifit::f_gsl (const gsl_vector * parameters, 
		     void * dataPoints, 
		     gsl_vector * f) {
  
  // ORSA_DEBUG("f...");
  
  for(unsigned int k=0; k<_par->size(); ++k) {
    _par->set(k,gsl_vector_get(parameters,k));
  }
  
  const orsa::MultifitData * data = (orsa::MultifitData *) dataPoints;
  
  // get advantage of parallel architectures
  // 'fun' must use its own cache to store all the values
  computeAllFunctionCalls(_par.get(), data, MODE_F);
  
  for (unsigned int j=0; j<data->size(); ++j) {
    const double fj = (__fun__(_par.get(),data,0,0,j) - data->getF(j))/(data->getSigma(j));
    // ORSA_DEBUG("f[%02i] = %10.3f", j, fj);
    gsl_vector_set (f, j, fj);
  }
  
  return GSL_SUCCESS;
}

double Multifit::_diff_two_points_ (const double & y_m,
				    const double & y_p) {
  return 0.5*(y_p-y_m);
}

double Multifit::_diff_five_points_ (const double & y_mm,
				     const double & y_m,
				     const double & y_p,
				     const double & y_pp) {
  
  // static const double coeff = 1.0/24.0;
  // five points rule
  // diff_ra  = coeff*(2.0*d_ra_mm -16.0*d_ra_m +16.0*d_ra_p -2.0*d_ra_pp);
  
  static const double coeff = 1.0/24.0;
  return (coeff*(2*(y_mm-y_pp) + 16.0*(y_p-y_m)));
}

int Multifit::df_gsl (const gsl_vector * v, 
		      void             * dataPoints, 
		      gsl_matrix       * J) {
  
  // ORSA_DEBUG("df...");
  
  for(unsigned int k=0; k<_par->size(); ++k) {
    _par->set(k,gsl_vector_get(v,k));
  }
  
  const orsa::MultifitData * data = (orsa::MultifitData *) dataPoints;
  
  // get advantage of parallel architectures
  // 'fun' must use its own cache to store all the values
  computeAllFunctionCalls(_par.get(), data, MODE_DF);
  
  for (unsigned int k=0; k<_par->size(); ++k) {
    for (unsigned int j=0; j<data->size(); ++j) {
      gsl_matrix_set (J, j, k, _diff_two_points_(__fun__(_par.get(),data,k,-1,j),
						 __fun__(_par.get(),data,k,+1,j)) / 
		      (data->getSigma(j) * _par->getDelta(k)));  
    }
  }
  
  return GSL_SUCCESS;
}

int Multifit::fdf_gsl (const gsl_vector * v, 
		       void * dataPoints, 
		       gsl_vector * f, 
		       gsl_matrix * J) {
  
  // ORSA_DEBUG("fdf...");
  
  for(unsigned int k=0; k<_par->size(); ++k) {
    _par->set(k,gsl_vector_get(v,k));
  }
  
  const orsa::MultifitData * data = (orsa::MultifitData *) dataPoints;
  
  // get advantage of parallel architectures
  // 'fun' must use its own cache to store all the values
  computeAllFunctionCalls(_par.get(), data, MODE_FDF);
  
  f_gsl( v, dataPoints, f);
  
  df_gsl(v, dataPoints, J);
  
  return GSL_SUCCESS;
}

bool Multifit::run() {
  
  // initial checks
  if (_data.get() == 0) {
    ORSA_ERROR("no data");
    return false;
  }	
  //
  if (_par.get() == 0) {
    ORSA_ERROR("no parameters");
    return false;
  }
  //
  ORSA_DEBUG("n: %i  p: %i",_data->size(),_par->size());
  //
  if (_data->size() == 0) {
    ORSA_ERROR("no data");
    return false;
  }
  //
  if (_par->size() == 0) {
    ORSA_ERROR("no parameters");
    return false;
  }
  
  gsl_multifit_fdfsolver * s = gsl_multifit_fdfsolver_alloc(gsl_multifit_fdfsolver_lmsder, 
							    _data->size(),
							    _par->size());
  
  MultifitGlue::instance()->_mf = this;
  
  gsl_multifit_function_fdf mf;
  // using the glue...
  mf.f   = &multifit_global_f_gsl;
  mf.df  = &multifit_global_df_gsl;
  mf.fdf = &multifit_global_fdf_gsl;
  mf.n   = _data->size();
  mf.p   = _par->size();
  mf.params = _data.get(); // sorry for the name confusion here
  
  gsl_vector * x = gsl_vector_alloc(_par->size());
  //
  for (unsigned int k=0; k<_par->size(); ++k) {
    gsl_vector_set(x, k, _par->get(k));
  }
  
  gsl_multifit_fdfsolver_set(s,&mf,x);
  
  unsigned int iter = 0;
  const unsigned int local_max_iter = 100000;
  int it_status;
  int cv_status;
  do {
    
    ++iter;
    
    it_status = gsl_multifit_fdfsolver_iterate(s);
    //
    ORSA_DEBUG("itaration status = %s",gsl_strerror(it_status));
    
    cv_status = gsl_multifit_test_delta(s->dx, s->x, 1.0e-12, 1.0e-12);
    //
    /* 
       {
       gsl_vector * g = gsl_vector_alloc(_par->size());
       gsl_multifit_gradient(s->J, s->f, g);
       //
       cv_status = gsl_multifit_test_gradient(g, 1.0e-3); 
       //
       gsl_vector_free(g);
       }
    */
    //
    ORSA_DEBUG("convergence status = %s", gsl_strerror(cv_status));
    
    ORSA_DEBUG("iter: %i",iter);
    
    if (0) {
      // debug only
      ORSA_DEBUG("iter: %i",iter);
      for (unsigned int k=0; k<_par->size(); ++k) {
	ORSA_DEBUG("par[%i] = \"%s\" = %+18.8f",
		   k,
		   _par->name(k).c_str(),
		   _par->get(k));
      } 
      
      double c = 1.0;
      //
      if (mf.n > mf.p) {
	double chi = gsl_blas_dnrm2(s->f);
	double dof = mf.n - mf.p;
	c = GSL_MAX_DBL(1.0, chi / sqrt(dof)); 
	ORSA_DEBUG("chisq/dof = %g",  chi*chi/dof);
      }
    }
    
    {
      // stop if chisq/dof < epsilon
      if (mf.n > mf.p) {
	const double chi = gsl_blas_dnrm2(s->f);
	const double dof = mf.n - mf.p;
	if (chi*chi/dof < orsa::epsilon()) {
	  ORSA_DEBUG("chisq/dof too small, interrupting...");
	  break;
	}
      }
    }

    {
      // stop if a parameter or its estimated uncertainty is not finite
      gsl_matrix * covar = gsl_matrix_alloc(_par->size(),_par->size());
      gsl_multifit_covar(s->J, 0.0, covar);
      bool break_main_iteration=false;
      for (unsigned int k=0; k<_par->size(); ++k) {
	if (!finite(gsl_vector_get(s->x, k))) {
	  ORSA_DEBUG("interrupting because parameter [%s] is not finite.",
		     _par->name(k).c_str());
	  break_main_iteration=true;
	  break;
	}
	if (!finite(gsl_matrix_get(covar,k,k))) {
	  ORSA_DEBUG("interrupting because uncertainty of parameter [%s] is not finite.",
		     _par->name(k).c_str());
	  break_main_iteration=true;
	  break;
	}
      }
      gsl_matrix_free(covar);	
      if (break_main_iteration) break; // this trick avoids leaking memory, i.e. leaving before freeing covar memory
    }
    
    if (logFile.isSet()) {
      
      FILE * fp = fopen(logFile.getRef().c_str(),"a");
      if (fp == 0) {
	ORSA_ERROR("cannot open file %s",logFile.getRef().c_str());
	return false;
      }
      
      gmp_fprintf(fp,"iter: %i\n",iter);
      
      double c = 1.0;
      //
      if (mf.n > mf.p) {
	double chi = gsl_blas_dnrm2(s->f);
	double dof = mf.n - mf.p;
	c = GSL_MAX_DBL(1.0, chi / sqrt(dof)); 
	// ORSA_DEBUG("chisq/dof = %g",  chi*chi/dof);
	gmp_fprintf(fp,"chisq/dof = %g\n",  chi*chi/dof);
      }
      
      gsl_matrix * covar = gsl_matrix_alloc(_par->size(),_par->size());
      gsl_multifit_covar(s->J, 0.0, covar);
      
      // #define FIT(i) gsl_vector_get(s->x, i)
      // #define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
      
      ORSA_DEBUG("appending to file %s",logFile.getRef().c_str());
      
      for (unsigned int k=0; k<_par->size(); ++k) {
	_par->set(k, gsl_vector_get(s->x,k));
      }
      
      for (unsigned int k=0; k<_par->size(); ++k) {
	gmp_fprintf(fp,"par[%02i] = [%32s] = %18.8g +/- %18.8g (%g * %g)\n",
		    k,
		    _par->name(k).c_str(),
		    _par->get(k),
		    c*sqrt(gsl_matrix_get(covar,k,k)),
		    c,
		    sqrt(gsl_matrix_get(covar,k,k)));
      }
      
      {
	gmp_fprintf(fp,"dx:");
	for (unsigned int k=0; k<_par->size(); ++k) {
	  gmp_fprintf(fp," %+16.12e",
		      gsl_vector_get(s->dx,k));
	}
	gmp_fprintf(fp,"\n");
	
	// normalized dx
	double dl = 0;
	for (unsigned int k=0; k<_par->size(); ++k) {
	  double dk = gsl_vector_get(s->dx,k);
	  dl += dk*dk;
	}
	const double l = sqrt(dl);
	if (l>0) {
	  gmp_fprintf(fp,"nx:");
	  for (unsigned int k=0; k<_par->size(); ++k) {
	    gmp_fprintf(fp," %+8.6f",
			gsl_vector_get(s->dx,k)/l);
	  }
	  gmp_fprintf(fp,"\n");
	}
      }
      
      gsl_matrix_free(covar); 
      fclose(fp);
    }
    
    for (unsigned int k=0; k<_par->size(); ++k) {
      _par->set(k, gsl_vector_get(s->x,k));
    }
    
    /* 
       if (0) {
       // update par delta...
       double c = 1.0;
       //
       if (mf.n > mf.p) {
       double chi = gsl_blas_dnrm2(s->f);
       double dof = mf.n - mf.p;
       c = GSL_MAX_DBL(1.0, chi / sqrt(dof)); 
       }
       gsl_matrix * covar = gsl_matrix_alloc(_par->size(),_par->size());
       gsl_multifit_covar(s->J, 0.0, covar);
       #define FIT(i) gsl_vector_get(s->x, i)
       #define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
       const double deltaFactor = 0.001;
       for (unsigned int k=0; k<_par->size(); ++k) {
       _par->setDelta(k, deltaFactor*c*ERR(k));
       }
       gsl_matrix_free(covar); 
       }
    */
    
    singleIterationDone(s); 
    singleIterationDone(_par.get());  
    
  } while (((cv_status == GSL_CONTINUE) || (it_status == GSL_CONTINUE)) && (iter < local_max_iter));
  
  if (cv_status == GSL_SUCCESS) {
    success(s);
    success(_par.get());
  }	
  
  for (unsigned int k=0; k<_par->size(); ++k) {
    _par->set(k, gsl_vector_get(s->x,k));
  }
  
  // covariance matrix here...
  /* 
     {
     gsl_matrix * covar = gsl_matrix_alloc(_par->size(),_par->size());
     gsl_multifit_covar(s->J, 0.0, covar);
     
     #define FIT(i) gsl_vector_get(s->x, i)
     #define ERR(i) sqrt(gsl_matrix_get(covar,i,i))
     
     for (unsigned int k=0; k<_par->size(); ++k) {
     ORSA_DEBUG("par[%i] = \"%s\" = %+18.8f +/- %+18.8f",
     k,
     _par->name(k).c_str(),
     FIT(k),
     ERR(k));
     }
     
     gsl_matrix_free(covar); 
     }
  */
  
  gsl_multifit_fdfsolver_free (s);
  gsl_vector_free (x);
  
  if (MultifitGlue::instance()->_mf != this) {
    ORSA_ERROR("mf should never change during the run... use a lock?");
    return false;
  }
  
  return true;
}

void Multifit::setMultifitParameters(orsa::MultifitParameters * mp) {
  _par = mp;
}	

const orsa::MultifitParameters * Multifit::getMultifitParameters() const {
  return _par.get();
}	

void Multifit::setMultifitData(orsa::MultifitData * md) {
  _data = md;
}	

const orsa::MultifitData * Multifit::getMultifitData() const {
  return _data.get();
}	

