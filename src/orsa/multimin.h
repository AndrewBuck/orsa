#ifndef _ORSA_MULTIMIN_
#define _ORSA_MULTIMIN_

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <gsl/gsl_multimin.h>

#include <orsa/cache.h>
#include <orsa/double.h>

#include <string>
#include <vector>

namespace orsa {
  
    class Multimin;
  
    class MultiminGlue {
    public:
        static MultiminGlue * instance();
    protected:
        MultiminGlue();
    public:
        virtual ~MultiminGlue();
    
    public:
        osg::ref_ptr<const Multimin> _mm;
    
    protected:
        static MultiminGlue * _instance;
    };
  
    class MultiminParameters : public osg::Referenced {
    public:
        MultiminParameters();
    protected:
        virtual ~MultiminParameters();
    
    public:
        //! returns the index
        unsigned int insert(const double & initialValue,
                            const double & step);
    public:
        //! returns the index
        unsigned int insert(const  std::string & name,
                            const double & initialValue,
                            const double & step);
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
        bool setRange(const std::string  & name,
                      const double & min,
                      const double & max);
        bool setRange(const unsigned int index,
                      const double & min,
                      const double & max);
    public:
        bool setRangeMin(const std::string & name,
                         const double      & min);
        bool setRangeMax(const std::string & name,
                         const double      & max);
        bool setRangeMin(const unsigned int   index,
                         const double       & min);
        bool setRangeMax(const unsigned int   index,
                         const double       & max);
    public:
        void setInRange(const std::string & name);
        void setInRange(const unsigned int index);
    public:
        const double & get(const std::string & name) const;
        const double & get(const unsigned int index) const;
    public:
        // bool haveRange(const std::string & name) const;
        // bool haveRange(const unsigned int index) const;
    public:
        const orsa::Cache<double> & getRangeMin(const std::string & name) const;
        const orsa::Cache<double> & getRangeMax(const std::string & name) const;
        const orsa::Cache<double> & getRangeMin(const unsigned int index) const;
        const orsa::Cache<double> & getRangeMax(const unsigned int index) const;
    public:
        const double & getStep(const std::string & name) const;
        const double & getStep(const unsigned int index) const;
    public:
        unsigned int size() const {
            return _data.size();
        }
    private:
        class NVS {
        public:
            std::string name;
            double      value;
            double      step;
            //
            orsa::Cache<double> min, max;
        };
    private: 
        typedef std::vector<NVS> dataType;
        dataType _data;
    
    public:
        bool writeToFile(const std::string & fileName) const;
    public:
        bool readFromFile(const std::string & fileName);
    };
  
    bool operator == (const MultiminParameters &, 
                      const MultiminParameters &);
  
    bool operator != (const MultiminParameters &, 
                      const MultiminParameters &);
  
    // trick, using the glue...
    double multimin_global_f_gsl(const gsl_vector * v, 
                                 void *);
  
    void multimin_global_df_gsl(const gsl_vector * v, 
                                void *,
                                gsl_vector * df);
  
    void multimin_global_fdf_gsl(const gsl_vector * v, 
                                 void *,
                                 double * f, 
                                 gsl_vector * df);
  
    class Multimin : public osg::Referenced {
    public:
        Multimin();
    
    protected:
        virtual ~Multimin();
    
    public:
        virtual double fun(const orsa::MultiminParameters *) const = 0;
    
    protected:
        double __fun__(const gsl_vector *,
                       const orsa::MultiminParameters *,
                       const bool verbose = false) const;
    
    protected:
        enum computeAllCallsMode { 
            MODE_F   = 0x1,
            MODE_DF  = 0x2,
            MODE_FDF = MODE_F | MODE_DF
        };
    
    protected:
        virtual void computeAllFunctionCalls(const orsa::MultiminParameters *, 
                                             const computeAllCallsMode) const { }; 
    
    public:
        virtual double f_gsl(const gsl_vector * parameters, 
                             void *) const;
    
        virtual void df_gsl(const gsl_vector * parameters, 
                            void *,
                            gsl_vector * df) const;
    
        virtual void fdf_gsl(const gsl_vector * parameters, 
                             void *,
                             double * f, 
                             gsl_vector * df) const;
    
    protected:
        static bool parametersOutsideRange(const gsl_vector * parameters,
                                           const orsa::MultiminParameters * par) {
            for (unsigned int k=0; k<par->size(); ++k) {
	
                const double xVal = gsl_vector_get(parameters,k);
	
                if (par->getRangeMin(k).isSet()) {
                    if (xVal < par->getRangeMin(k).getRef()) {
                        return true;
                    }
                }
	
                if (par->getRangeMax(k).isSet()) {
                    if (xVal > par->getRangeMax(k).getRef()) {
                        return true;
                    }
                }
	
            }
      
            return false;
        }
    
    public:
        static double _diff_two_points_(const double & y_m,
                                        const double & y_p);
    
    public:
        static double _diff_five_points_(const double & y_mm,
                                         const double & y_m,
                                         const double & y_p,
                                         const double & y_pp);
    
    public:
        bool run_nmsimplex(const unsigned int maxIter = 4096,
                           const double       epsAbs  = 1.0e-3) const;
    
    public:
        bool run_conjugate_fr(const unsigned int maxIter         = 4096,
                              const double       initialStepSize = 0.01,
                              const double       tollerance      = 1.0e-3,
                              const double       epsAbs          = 1.0e-3) const;
    
    protected:
        virtual void singleIterationDone(const orsa::MultiminParameters *) const { }
    
    protected:
        // called only if converged with GSL_SUCCESS
        virtual void success(const orsa::MultiminParameters *) const { }
    
    public:
        void setMultiminParameters(orsa::MultiminParameters *);
    public:
        const orsa::MultiminParameters * getMultiminParameters() const;
    
    private:
        orsa::Cache<std::string> logFile;
    public:
        virtual void setLogFile(const std::string & lf) { logFile.set(lf); }
        virtual const std::string & getLogFile() const { return logFile.getRef(); }
    
    protected:
        osg::ref_ptr<orsa::MultiminParameters> _par;
    };
  
}; // namespace orsa

#endif // _ORSA_MULTIMIN_
