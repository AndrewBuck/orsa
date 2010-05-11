#include <orsa/multimin.h>

using namespace orsa;

// MultiminGlue

MultiminGlue * MultiminGlue::_instance = 0;

MultiminGlue * MultiminGlue::instance() {
    if (_instance == 0) {
        _instance = new MultiminGlue;
    }
    return _instance;
}

MultiminGlue::MultiminGlue() {
    _mm = 0;
}

MultiminGlue::~MultiminGlue() {
    _instance = 0;
}

// gsl functions

double orsa::multimin_global_f_gsl (const gsl_vector * v, 
                                    void *) {
    return MultiminGlue::instance()->_mm->f_gsl(v,0);
}

void orsa::multimin_global_df_gsl(const gsl_vector * v, 
                                  void *,
                                  gsl_vector * df) {
    return MultiminGlue::instance()->_mm->df_gsl(v,0,df);
}

void orsa::multimin_global_fdf_gsl(const gsl_vector * v, 
                                   void *,
                                   double * f, 
                                   gsl_vector * df) {
    return MultiminGlue::instance()->_mm->fdf_gsl(v,0,f,df);
}

// MultiminParameters

MultiminParameters::MultiminParameters() : Referenced(true) { } 

MultiminParameters::~MultiminParameters() { }

unsigned int MultiminParameters::insert(const double & initialValue,
                                        const double & step) {
    {
        NVS tmp;
        //
        tmp.name  = ""; // empty name
        tmp.value = initialValue;
        tmp.step  = step;
        //
        _data.push_back(tmp);
    }
    return (_data.size()-1);
}

unsigned int MultiminParameters::insert(const  std::string & name,
                                        const double & initialValue,
                                        const double & step) {
    dataType::const_iterator it = _data.begin();
    while (it != _data.end()) {
        if ((*it).name == name) {
            ORSA_ERROR("trying to insert duplicate variable: [%s] already inserted",name.c_str());
            // return false;
            return index(name);
        } 
        ++it;
    }
    {
        NVS tmp;
        //
        tmp.name  = name;
        tmp.value = initialValue;
        tmp.step  = step;
        //
        _data.push_back(tmp);
    }
    // return true;
    return index(name);
}

void MultiminParameters::clear() {
    _data.clear();
}

unsigned int MultiminParameters::index(const std::string & name) const {
    for (unsigned int k=0; k<_data.size(); ++k) {
        if (_data[k].name == name) {
            return k;
        }	
    }
    ORSA_ERROR("index not found: name=%s",name.c_str());
    return 0;
}

std::string MultiminParameters::name(const unsigned int index) const {
    if (_data.size() <= index) {
        ORSA_ERROR("name not found");
        return "";
    }
    return _data[index].name;
}

bool MultiminParameters::set(const std::string  & name,
                             const double & value) {
    return set(MultiminParameters::index(name),value);
}

bool MultiminParameters::set(const unsigned int   index,
                             const double & value) {
    _data[index].value = value;
    setInRange(index);
    return true;
}

bool MultiminParameters::setRange(const std::string  & name,
                                  const double & min,
                                  const double & max) {
    return setRange(MultiminParameters::index(name),
                    min,
                    max); 
}

bool MultiminParameters::setRange(const unsigned int   index,
                                  const double & min,
                                  const double & max) {
    _data[index].min = min;
    _data[index].max = max;
    setInRange(index);
    return true;
}

bool MultiminParameters::setRangeMin(const std::string & name,
                                     const double      & min) {
    return setRangeMin(MultiminParameters::index(name),min);
}

bool MultiminParameters::setRangeMax(const std::string & name,
                                     const double      & max) {
    return setRangeMax(MultiminParameters::index(name),max);
}

bool MultiminParameters::setRangeMin(const unsigned int   index,
                                     const double       & min) {
    _data[index].min = min;
    return true;
}

bool MultiminParameters::setRangeMax(const unsigned int   index,
                                     const double       & max) {
    _data[index].max = max;
    return true;
}

void MultiminParameters::setInRange(const std::string & name) {
    setInRange(MultiminParameters::index(name));
}

void MultiminParameters::setInRange(const unsigned int index) {
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

const double & MultiminParameters::get(const std::string & name) const {
    return get(MultiminParameters::index(name));
}

const double & MultiminParameters::get(const unsigned int index) const {
    return _data[index].value;
}

/* 
   bool MultiminParameters::haveRange(const std::string & name) const {
   return haveRange(MultiminParameters::index(name));
   }
*/

/* 
   bool MultiminParameters::haveRange(const unsigned int index) const {
   return (_data[index].min.isSet() _data[index].max.isSet());
   }
*/

const orsa::Cache<double> & MultiminParameters::getRangeMin(const std::string & name) const {
    return getRangeMin(MultiminParameters::index(name));
}

const orsa::Cache<double> & MultiminParameters::getRangeMax(const std::string & name) const {
    return getRangeMax(MultiminParameters::index(name));
}

const orsa::Cache<double> & MultiminParameters::getRangeMin(const unsigned int index) const {
    return (_data[index].min);
}

const orsa::Cache<double> & MultiminParameters::getRangeMax(const unsigned int index) const {
    return (_data[index].max);
}

const double & MultiminParameters::getStep(const std::string & name) const {
    return getStep(MultiminParameters::index(name));
}

const double & MultiminParameters::getStep(const unsigned int index) const {
    return _data[index].step;
}

bool MultiminParameters::writeToFile(const std::string & fileName) const {
    FILE * fp = fopen(fileName.c_str(),"w");
    if (fp == 0) {
        ORSA_ERROR("cannot open file %s",fileName.c_str());
        return false;
    }
    for (unsigned int k=0; k<size(); ++k) {
        gmp_fprintf(fp,"%16s %+22.16e %+22.16e\n",
                    name(k).c_str(),
                    get(k),
                    getStep(k));
    } 
    fclose(fp);
    return true;
}

bool MultiminParameters::readFromFile(const std::string & fileName) {
    clear();
    FILE * fp = fopen(fileName.c_str(),"r");
    if (fp == 0) {
        ORSA_ERROR("cannot open file %s",fileName.c_str());
        return false;
    }
    char name[1024];
    double value;
    double step;
    while (gmp_fscanf(fp,"%s %lf %lf\n",
                      name,
                      &value,
                      &step) == 3) {
        insert(name,value,step);
    }	
    fclose(fp);
    return true;
}

bool orsa::operator == (const MultiminParameters & p1, 
                        const MultiminParameters & p2) {
  
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

bool orsa::operator != (const MultiminParameters & p1, 
                        const MultiminParameters & p2) {
    return (!(p1 == p2));
}

// Multimin

Multimin::Multimin() : Referenced(true) { }

Multimin::~Multimin() { }

double Multimin::__fun__(const gsl_vector * parameters,
                         const orsa::MultiminParameters * par,
                         const bool verbose) const {
    double effectiveFun = fun(par);
    if (parametersOutsideRange(parameters, par)) {
        effectiveFun += 1.0e3 * (1.0 + fabs(effectiveFun));
    }
    if (verbose) {
        ORSA_DEBUG("effectiveFun: %g",effectiveFun);
    }
    //
    return effectiveFun;
}

double Multimin::f_gsl(const gsl_vector * parameters, 
                       void *) const {
  
    for(unsigned int k=0; k<_par->size(); ++k) {
        _par->set(k,gsl_vector_get(parameters,k));
    }
  
    // get advantage of parallel architectures
    // 'fun' must use its own cache to store all the values
    computeAllFunctionCalls(_par.get(), MODE_F);
  
    return __fun__(parameters,_par.get(),false);
}

void Multimin::df_gsl(const gsl_vector * parameters, 
                      void *,
                      gsl_vector * df) const {
  
    for(unsigned int k=0; k<_par->size(); ++k) {
        _par->set(k,gsl_vector_get(parameters,k));
    }
  
    osg::ref_ptr<orsa::MultiminParameters> par_mm = new orsa::MultiminParameters;
    osg::ref_ptr<orsa::MultiminParameters> par_m  = new orsa::MultiminParameters;
    osg::ref_ptr<orsa::MultiminParameters> par_p  = new orsa::MultiminParameters;
    osg::ref_ptr<orsa::MultiminParameters> par_pp = new orsa::MultiminParameters;
  
    // get advantage of parallel architectures
    // 'fun' must use its own cache to store all the values
    computeAllFunctionCalls(_par.get(), MODE_DF);
  
    for (unsigned int k=0; k<_par->size(); ++k) {
    
        // ORSA_DEBUG("df[%03i]...",k);
    
        // deep copies of the reference _par
        (*(par_mm.get())) = (*(_par.get()));
        (*(par_m.get()))  = (*(_par.get()));
        (*(par_p.get()))  = (*(_par.get()));
        (*(par_pp.get())) = (*(_par.get()));
    
        par_mm->set(k, _par->get(k) - 2*_par->getStep(k));
        par_m-> set(k, _par->get(k) -   _par->getStep(k));
        par_p-> set(k, _par->get(k) +   _par->getStep(k));
        par_pp->set(k, _par->get(k) + 2*_par->getStep(k));
    
        /* 
           gsl_vector_set (df, k, double(_diff_five_points_(__fun__(parameters,par_mm.get()),
           __fun__(parameters,par_m.get()),
           __fun__(parameters,par_p.get()),
           __fun__(parameters,par_pp.get())) / 
           _par->getStep(k) 
           ).get_d());  
        */
        //
        gsl_vector_set (df, k, _diff_two_points_(__fun__(parameters,par_m.get()),
                                                 __fun__(parameters,par_p.get())) / 
                        _par->getStep(k));  
        //
        /* 
           const double df2 = _diff_two_points_(__fun__(parameters,par_m.get()),
           __fun__(parameters,par_p.get()))/_par->getStep(k);
        */
        //
        // ORSA_DEBUG("df[%02i] = %+20.14f",k,gsl_vector_get(df,k));
        /*     
               ORSA_DEBUG("par[%02i] = %+020.14f   df5[%02i] = %+020.14f   df2[%02i] = %+020.14f",
               k,_par->get(k),
               k,gsl_vector_get(df,k),
               k,df2());
        */
    }
}

void Multimin::fdf_gsl(const gsl_vector * parameters, 
                       void *,
                       double * f, 
                       gsl_vector * df) const {
  
    for(unsigned int k=0; k<_par->size(); ++k) {
        _par->set(k,gsl_vector_get(parameters,k));
    }
  
    (*f) = f_gsl(parameters,0);
  
    df_gsl(parameters,0,df);
}

double Multimin::_diff_two_points_ (const double & y_m,
                                    const double & y_p) {
    /* 
       ORSA_DEBUG("d2p: y_m=%f   y_p=%f   retVal: %f",
       y_m(),
       y_p(),
       double(0.5*(y_p-y_m)));
    */
  
    return 0.5*(y_p-y_m);
}

double Multimin::_diff_five_points_ (const double & y_mm,
                                     const double & y_m,
                                     const double & y_p,
                                     const double & y_pp) {
  
    // static const double coeff = 1.0/24.0;
    // five points rule
    // diff_ra  = coeff*(2.0*d_ra_mm -16.0*d_ra_m +16.0*d_ra_p -2.0*d_ra_pp);
  
    static const double coeff = 1.0/24.0;
    return (coeff*(2*(y_mm-y_pp) +
                   16*(y_p-y_m)));
}

bool Multimin::run_nmsimplex(const unsigned int maxIter,
                             const double       epsAbs) const {
  
    if (_par.get() == 0) {
        ORSA_ERROR("no parameters");
        return false;
    }
  
    gsl_multimin_fminimizer * s = 
        gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex, 
                                      _par->size());
  
    MultiminGlue::instance()->_mm = this;
  
    gsl_multimin_function mf;
  
    mf.f      = &multimin_global_f_gsl;
    mf.n      = _par->size();
    mf.params = 0; // it's a (void *)
  
    gsl_vector * x = gsl_vector_alloc(_par->size());
    //
    for (unsigned int k=0; k<_par->size(); ++k) {
        gsl_vector_set(x, k, _par->get(k));
    }
  
    gsl_vector * step = gsl_vector_alloc (_par->size());
    //
    for (unsigned int k=0; k<_par->size(); ++k) {
        gsl_vector_set(step, k, _par->getStep(k));
    }
  
    gsl_multimin_fminimizer_set(s,&mf,x,step);
  
    unsigned int iter = 0;
    const unsigned int local_max_iter = maxIter;
    int it_status;
    int cv_status;
    do {
        ++iter;
    
        it_status = gsl_multimin_fminimizer_iterate (s);
        //
        // ORSA_DEBUG("itaration status = %s",gsl_strerror(it_status));
    
        // if (it_status == GSL_ENOPROG) break;
    
        // ORSA_DEBUG("s->size: %g",s->size);
    
        cv_status = gsl_multimin_test_size(s->size, epsAbs);
        //
        /* 
           {
           gsl_vector * g = gsl_vector_alloc(_par->size());
           gsl_multimin_gradient(s->J, s->f, g);
           //
           cv_status = gsl_multimin_test_gradient(g, 1.0e-3); 
           //
           gsl_vector_free(g);
           }
        */
        //
        // ORSA_DEBUG("convergence status = %s", gsl_strerror(cv_status));
    
        if (0) {
            // debug only
            ORSA_DEBUG("iter: %i",iter);
            for (unsigned int k=0; k<_par->size(); ++k) {
                ORSA_DEBUG("par[%03i] = \"%s\" = %+24.15e",
                           k,
                           _par->name(k).c_str(),
                           _par->get(k));
            }
        }
    
        if (0) {
            if (logFile.isSet()) {
	
                FILE * fp = fopen(logFile.getRef().c_str(),"a");
                if (fp == 0) {
                    ORSA_ERROR("cannot open file %s",logFile.getRef().c_str());
                    return false;
                }
	
                fclose(fp);
            }
        }
    
        { 
            for (unsigned int k=0; k<_par->size(); ++k) {
                _par->set(k, gsl_vector_get(s->x,k));
            }
      
            singleIterationDone(_par.get());
        }
    
    } while (((cv_status == GSL_CONTINUE) || (it_status == GSL_CONTINUE)) && (iter < local_max_iter));
  
    if (cv_status == GSL_SUCCESS) {
        success(_par.get());
    }	
  
    // ORSA_DEBUG("total iterations: %i",iter);
  
    for (unsigned int k=0; k<_par->size(); ++k) {
        _par->set(k, gsl_vector_get(s->x,k));
    }
  
    gsl_multimin_fminimizer_free (s);
    gsl_vector_free (x);
    gsl_vector_free (step);
  
    if (MultiminGlue::instance()->_mm != this) {
        ORSA_ERROR("mf should never change during the run... use a lock?");
        return false;
    }
  
    if ((cv_status == GSL_CONTINUE) || (it_status == GSL_CONTINUE)) {
        return false;
    } else {
        return true;
    }
}

bool Multimin::run_conjugate_fr(const unsigned int maxIter,
                                const double       initialStepSize,
                                const double       tollerance,
                                const double       epsAbs) const {
  
    if (_par.get() == 0) {
        ORSA_ERROR("no parameters");
        return false;
    }

    /* 
       gsl_multimin_fdfminimizer * s = 
       gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr, 
       _par->size());
    */
    //
    gsl_multimin_fdfminimizer * s = 
        gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs, 
                                        _par->size());
  
    MultiminGlue::instance()->_mm = this;
  
    gsl_multimin_function_fdf mf;
  
    mf.f      = &multimin_global_f_gsl;
    mf.df     = &multimin_global_df_gsl;
    mf.fdf    = &multimin_global_fdf_gsl;
    mf.n      = _par->size();
    mf.params = 0; // it's a (void *)
  
    gsl_vector * x = gsl_vector_alloc(_par->size());
    //
    for (unsigned int k=0; k<_par->size(); ++k) {
        gsl_vector_set(x, k, _par->get(k));
    }
  
    gsl_multimin_fdfminimizer_set(s,&mf,x,initialStepSize,tollerance);
  
    unsigned int iter = 0;
    const unsigned int local_max_iter = maxIter;
    int it_status;
    int cv_status;
    do {
        ++iter;
    
        ORSA_DEBUG("iter: %i",iter);
    
        if (0) {
            // debug only
            ORSA_DEBUG("iter: %i",iter);
            for (unsigned int k=0; k<_par->size(); ++k) {
                ORSA_DEBUG("par[%03i] = \"%s\" = %+24.16e",
                           k,
                           _par->name(k).c_str(),
                           _par->get(k));
            }
        }
    
        {
            //debug
            for (unsigned int k=0; k<1; ++k) {
                ORSA_DEBUG("gradient[%i] = %e",k,gsl_vector_get(s->gradient,k));
            }
        }
    
        it_status = gsl_multimin_fdfminimizer_iterate(s);
        //
        ORSA_DEBUG("itaration status = %s",gsl_strerror(it_status));
    
        // if (it_status == GSL_ENOPROG) break;
    
        // cv_status = gsl_multimin_test_size(s->size, epsAbs);
        //
        cv_status = gsl_multimin_test_gradient(s->gradient, epsAbs); 
    
        ORSA_DEBUG("convergence status = %s", gsl_strerror(cv_status));
    
        if (0) {
            if (logFile.isSet()) {
	
                FILE * fp = fopen(logFile.getRef().c_str(),"a");
                if (fp == 0) {
                    ORSA_ERROR("cannot open file %s",logFile.getRef().c_str());
                    return false;
                }
	
                fclose(fp);
            }
        }
    
        { 
            for (unsigned int k=0; k<_par->size(); ++k) {
                _par->set(k, gsl_vector_get(s->x,k));
            }
      
            singleIterationDone(_par.get());
        }
    
    } while (((cv_status == GSL_CONTINUE) || (it_status == GSL_CONTINUE)) && (iter < local_max_iter));
  
    if (cv_status == GSL_SUCCESS) {
        success(_par.get());
    }	
  
    for (unsigned int k=0; k<_par->size(); ++k) {
        _par->set(k, gsl_vector_get(s->x,k));
    }
  
    gsl_multimin_fdfminimizer_free (s);
    gsl_vector_free (x);
  
    if (MultiminGlue::instance()->_mm != this) {
        ORSA_ERROR("mf should never change during the run... use a lock?");
        return false;
    }
  
    return true;
}

void Multimin::setMultiminParameters(orsa::MultiminParameters * mp) {
    _par = mp;
}	

const orsa::MultiminParameters * Multimin::getMultiminParameters() const {
    return _par.get();
}	
