#ifndef _ORSA_MULTIFIT_
#define _ORSA_MULTIFIT_

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>

#include <orsa/cache.h>
#include <orsa/double.h>
#include <orsa/vector.h>

#include <string>
#include <vector>

namespace orsa {
  
  class Multifit;
  
  class MultifitGlue {
  public:
    static MultifitGlue * instance();
  protected:
    MultifitGlue();
  public:
    virtual ~MultifitGlue();
    
  public:
    osg::ref_ptr<Multifit> _mf;
    
  protected:
    static MultifitGlue * _instance;
  };
  
  class MultifitParameters : public osg::Referenced {
  public:
    MultifitParameters();
  protected:
    virtual ~MultifitParameters();
  public:
    bool insert(const  std::string & name,
		const double & initialValue,
		const double & delta);
  public:
    void clear();
  public:
    unsigned int index(const std::string & name) const;
    std::string   name(const unsigned int  index) const;
  public:
    bool set(const std::string  & name,
	     const double & value);
    bool set(const unsigned int   index,
	     const double & value);
  public:
    bool setDelta(const std::string  & name,
		  const double & value);
    bool setDelta(const unsigned int   index,
		  const double & value);
  public:
    double get(const std::string & name) const;
    double get(const unsigned int index) const;
  public:
    double getDelta(const std::string & name) const;
    double getDelta(const unsigned int index) const;
  public:
    unsigned int size() const {
      return _data.size();
    }
  private:
    class NVD {
    public:
      std::string  name;
      double value;
      double delta;
    };
  private: 
    typedef std::vector<NVD> dataType;
    dataType _data;
    
  public:
    bool writeToFile(const std::string & fileName) const;
  public:
    bool readFromFile(const std::string & fileName);
    
    /* 
       public:
       double covariance(const std::string & name1,
       const std::string & name2) const;
       public:
       double covariance(const unsigned int index1,
       const unsigned int index2) const;
    */
  };
  
  bool operator == (const MultifitParameters &, 
		    const MultifitParameters &);
  
  bool operator != (const MultifitParameters &, 
		    const MultifitParameters &);
  
  class MultifitData : public osg::Referenced {
  public:
    MultifitData();
  protected:
    virtual ~MultifitData();
  public:
    bool insertVariable(const std::string & name);
  public:
    bool insertZ(const std::string  & name,
		 const unsigned int   row,
		 const mpz_class & value);
    bool insertZ(const unsigned int   index,
		 const unsigned int   row,
		 const mpz_class & value);
  public:
    bool insertD(const std::string  & name,
		 const unsigned int   row,
		 const double & value);
    bool insertD(const unsigned int   index,
		 const unsigned int   row,
		 const double & value);
  public:
    bool insertV(const std::string  & name,
		 const unsigned int   row,
		 const orsa::Vector & value);
    bool insertV(const unsigned int   index,
		 const unsigned int   row,
		 const orsa::Vector & value);
  public:
    bool insertF(const unsigned int   row,
		 const double & value);
  public:
    bool insertSigma(const unsigned int   row,
		     const double & value);
  public:
    unsigned int index(const std::string & name) const;
    std::string   name(const unsigned int  index) const;
  public:
    mpz_class getZ(const std::string & name,
		   const unsigned int  row) const;
    mpz_class getZ(const unsigned int  index,
		   const unsigned int  row) const;
  public:
    double getD(const std::string & name,
		const unsigned int  row) const;
    double getD(const unsigned int  index,
		const unsigned int  row) const;
  public:
    Vector getV(const std::string & name,
		const unsigned int  row) const;
    Vector getV(const unsigned int  index,
		const unsigned int  row) const;
  public:
    double getF(const unsigned int row) const;
  public:
    double getSigma(const unsigned int row) const;
  public:
    unsigned int size() const;
    unsigned int vars() const;
    /* 
       private:
       class NV {
       public:
       orsa::Cache<std::string>                name;
       std::vector<orsa::Cache<double> > value;
       };
    */
  private:
    /* class NZD {
       public:
       orsa::Cache<std::string>                name;
       std::vector<orsa::Cache<mpz_class> >    z;
       std::vector<orsa::Cache<double> > d;
       };
    */
  private:
    class NZDV {
    public:
      orsa::Cache<std::string>                name;
      std::vector<orsa::Cache<mpz_class> >    z;
      std::vector<orsa::Cache<double> > d;
      std::vector<orsa::Cache<orsa::Vector> > v;
    };
  private:
    class VFS {
    public:
      std::vector<MultifitData::NZDV>         var;
      std::vector<orsa::Cache<double> > f;
      std::vector<orsa::Cache<double> > sigma;
    };
  private: 
    MultifitData::VFS _data;
  public:
    void clear() {
      _data.var.clear();
      _data.f.clear();
      _data.sigma.clear();
    }
  };
  
  // trick, using the glue...
  int multifit_global_f_gsl (const gsl_vector * v, 
			     void * dataPoints, 
			     gsl_vector * f);
  int multifit_global_df_gsl (const gsl_vector * v, 
			      void * dataPoints, 
			      gsl_matrix * J);
  int multifit_global_fdf_gsl (const gsl_vector * v, 
			       void * dataPoints, 
			       gsl_vector * f, 
			       gsl_matrix * J);
  
  class Multifit : public osg::Referenced {
  public:
    Multifit();
    
  protected:
    virtual ~Multifit();
    
    /* protected:
       virtual double fun(const orsa::MultifitParameters *, 
       const orsa::MultifitData *,
       const unsigned int row) const = 0;
    */
    //
  protected:
    virtual double fun(const orsa::MultifitParameters *, 
			     const orsa::MultifitData *,
			     const unsigned int p, // par index
			     const int          d, // delta
			     const unsigned int row) const = 0;
    
  protected:
    enum computeAllCallsMode { 
      MODE_F   = 0x1,
      MODE_DF  = 0x2,
      MODE_FDF = MODE_F | MODE_DF
    };
    
  protected:
    virtual void computeAllFunctionCalls(const orsa::MultifitParameters *, 
					 const orsa::MultifitData *,
					 const computeAllCallsMode) const = 0;
    
  public:
    virtual int f_gsl (const gsl_vector * parameters, 
		       void * dataPoints, 
		       gsl_vector * f);
    
  public:
    static double _diff_two_points_(const double & y_m,
				    const double & y_p);
    
  public:
    static double _diff_five_points_(const double & y_mm,
				     const double & y_m,
				     const double & y_p,
				     const double & y_pp);
    
  public:
    virtual int df_gsl (const gsl_vector * v, 
			void * dataPoints, 
			gsl_matrix * J);
    
  public:
    virtual int fdf_gsl (const gsl_vector * v, 
			 void * dataPoints, 
			 gsl_vector * f, 
			 gsl_matrix * J);
    
  public:
    bool run();
    
  protected:
    virtual void singleIterationDone(const gsl_multifit_fdfsolver *) const { }
    
  protected:
    // called only if converged with GSL_SUCCESS
    virtual void success(const gsl_multifit_fdfsolver *) const { }
    
  public:
    void setMultifitParameters(orsa::MultifitParameters *);
  public:
    const orsa::MultifitParameters * getMultifitParameters() const;
    
  public:
    void setMultifitData(orsa::MultifitData *);
  public:
    const orsa::MultifitData * getMultifitData() const;
    
  private:
    orsa::Cache<std::string> logFile;
  public:
    virtual void setLogFile(const std::string & lf) { 
      logFile.set(lf);
      // erase it
      FILE * fp = fopen(logFile.getRef().c_str(),"w");
      if (fp == 0) {
	ORSA_ERROR("cannot open file %s",logFile.getRef().c_str());
      }
      fclose(fp);
    }
    virtual const std::string & getLogFile() const { return logFile.getRef(); }
    
  protected:
    osg::ref_ptr<orsa::MultifitParameters> _par;
    osg::ref_ptr<orsa::MultifitData>       _data;
  };
  
}; // namespace orsa

#endif // _ORSA_MULTIFIT_
