#ifndef _ORSA_UTIL_
#define _ORSA_UTIL_

#include <string>

#include <orsa/body.h>
#include <orsa/bodygroup.h>
#include <orsa/cache.h>
#include <orsa/datetime.h>
#include <orsa/double.h>
#include <orsa/massDistribution.h>
#include <orsa/matrix.h>
#include <orsa/quaternion.h>
#include <orsa/shape.h>

#include <map>

#include <QHash>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace orsa {
  
    std::string & removeLeadingAndTrailingSpaces(std::string & s);
  
    std::string & removeAllSpaces(std::string & s);
  
    std::string & removeLeadingPlusSign(std::string & s);
  
    bool FrenetSerret(const orsa::Body * b,
                      orsa::BodyGroup  * bg,
                      const orsa::Time & t,
                      const orsa::Time & dt,
                      orsa::Vector & T,
                      orsa::Vector & N,
                      orsa::Vector & B);
  
    class IntPowCache {
    public:
        IntPowCache(const double & x) : _x(x) { }
    public:
        virtual ~IntPowCache() { }
    public:
        /* 
           const double get(const int p) const {
           if (!_data[p].isSet()) {
           // this could be more efficient...
           _data[p] = int_pow(_x,p);
           }
           return _data[p].getRef();
           }
        */
        //
        inline double get(const int p) const {      
            if (_data[p].isSet()) {
                return _data[p].getRef();
            } else {
                if (p > 0) {
                    _data[p] = get(p-1) * _x;
                } else if (p < 0) {
                    _data[p] = get(p+1) / _x;
                } else {
                    _data[p] = 1;
                }
                return _data[p].getRef();
            }
        }
    protected:
        const double _x;
    protected:
        // mutable QHash <int, orsa::Cache<double> > _data;
        mutable std::map <int, orsa::Cache<double> > _data;
    };
  
    bool eulerAnglesToMatrix(orsa::Matrix & m,
                             const double & psi,
                             const double & theta,
                             const double & phi);
  
    bool matrixToEulerAngles(double       & psi,
                             double       & theta,
                             double       & phi,
                             const orsa::Matrix & m);
  
    // wrapper
    /* 
       inline bool eulerAnglesToMatrix(orsa::Matrix                 & m,
       const RotationalBodyProperty * rotational) {
       if (rotational) {
       return eulerAnglesToMatrix(m, 
       rotational->getPsi(),
       rotational->getTheta(),
       rotational->getPhi());
       } else {
       // setting m should not be needed, but doesn't hurt
       m = orsa::Matrix::identity();
       return false;
       }
       }
    */
  
    orsa::Matrix QuaternionToMatrix (const orsa::Quaternion &);
  
    orsa::Quaternion MatrixToQuaternion (const orsa::Matrix &);
  
    orsa::Matrix localToGlobal(const orsa::Body       * b,
                               const orsa::BodyGroup  * bg,
                               const orsa::Time       & t);
  
    orsa::Matrix globalToLocal(const orsa::Body       * b,
                               const orsa::BodyGroup  * bg,
                               const orsa::Time       & t);
  
    // p = albedo, H = absolute magnitude
    double asteroidDiameter(const double & p, 
                            const double & H);
  
    /* 
       class ConstantZRotation : public PrecomputedRotationalBodyProperty {
       public:	
       ConstantZRotation(const orsa::Time   & t0,
       const double & phi0,
       const double & omega) : 
       PrecomputedRotationalBodyProperty(),
       _t0(t0), 
       _phi0(phi0),
       _omega(omega) { }
       public:    
       bool update(const orsa::Time & t);
       public:
       bool get(double & phi,
       double & theta,
       double & psi,
       double & phiDot,
       double & thetaDot,
       double & psiDot) const {
       if (!matrixToEulerAngles(psi,
       theta,
       phi,
       _m.getRef())) {
       ORSA_DEBUG("problems...");
       }
     
       // ORSA_DEBUG("implement dot section...");
       #warning "implement dot section..."
     
       return true;
       }
       public:
       double getPhi() const { 
       double phi,theta,psi,phiDot,thetaDot,psiDot;
       if (!get(phi, theta, psi, phiDot, thetaDot, psiDot)) {
       ORSA_DEBUG("problems...");
       }
       return phi;
       }
       double getTheta() const { 
       double phi,theta,psi,phiDot,thetaDot,psiDot;
       if (!get(phi, theta, psi, phiDot, thetaDot, psiDot)) {
       ORSA_DEBUG("problems...");
       }
       return theta;
       }
       double getPsi() const { 
       double phi,theta,psi,phiDot,thetaDot,psiDot;
       if (!get(phi, theta, psi, phiDot, thetaDot, psiDot)) {
       ORSA_DEBUG("problems...");
       }
       return psi;
       }
       double getPhiDot()   const { ORSA_DEBUG("implement dot section..."); return 0; }
       double getThetaDot() const { ORSA_DEBUG("implement dot section..."); return 0; }
       double getPsiDot()   const { ORSA_DEBUG("implement dot section..."); return 0; }
       public:
       RotationalBodyProperty * clone() const {
       return new ConstantZRotation(*this);
       }	
       private:
       const orsa::Time   _t0;
       const double _phi0;
       const double _omega;
       private:
       orsa::Cache<orsa::Matrix> _m;
       };
    */
  
    void principalAxis(orsa::Matrix & genericToPrincipal,
                       orsa::Matrix & principalInertiaMatrix,
                       const orsa::Matrix & inertiaMatrix);
  
    // RNG class has two main purposes: 
    // - automatic memory handling; 
    // - state saving/restoring for checkpointing
    //
    // checkpoints: gsl_rng_fwrite and gsl_rng_fread (WARNING: binary files not portable across different architectures)
    // http://www.gnu.org/software/gsl/manual/html_node/Reading-and-writing-random-number-generator-state.html
    //
    // TODO: make this class thread safe?
    // 
    class RNG : public osg::Referenced  {
    public:
        RNG(int randomSeed_in) : 
            osg::Referenced(true),
            randomSeed(randomSeed_in) {
            commonInit();
        }
    protected:
        void commonInit() {
            rnd = ::gsl_rng_alloc(gsl_rng_gfsr4);
            ::gsl_rng_set(rnd,randomSeed);
        }
    protected:
        virtual ~RNG() {
            ::gsl_rng_free(rnd); 
        }
    public:
        int gsl_rng_fwrite(FILE * stream) const {
            return ::gsl_rng_fwrite(stream,rnd);
        }
    public:
        int gsl_rng_fread(FILE * stream) const {
            return ::gsl_rng_fread(stream,rnd);
        }
    public:
        double gsl_rng_uniform() const {
            return ::gsl_rng_uniform(rnd);
        }
    public:
        double gsl_rng_uniform_pos() const {
            return ::gsl_rng_uniform_pos(rnd);
        }
    public:
        unsigned long int gsl_rng_uniform_int(unsigned long int n) const {
            return ::gsl_rng_uniform_int(rnd,n);
        }
    public:
        void gsl_ran_dir_2d(double * x, double * y) const {
            ::gsl_ran_dir_2d(rnd,x,y);
        }
    public:
        void gsl_ran_dir_3d(double * x, double * y, double * z) const {
            ::gsl_ran_dir_3d(rnd,x,y,z);
        }
    public:
        double gsl_ran_gaussian(double sigma) const {
            return ::gsl_ran_gaussian(rnd,sigma);
        }
    public:
        double gsl_ran_laplace(double a) const {
            return ::gsl_ran_laplace(rnd,a);
        }    
    public:
        const int randomSeed;
    protected:  
        gsl_rng * rnd;
    };
  
    /***/
  
    class RandomPointsInShape : public osg::Referenced {
    public:
        RandomPointsInShape(const orsa::Shape * shape,
                            const unsigned int N,
                            const int randomSeed);
    protected:
        virtual ~RandomPointsInShape() { } 
    public:
        osg::ref_ptr<const orsa::Shape> shape;
    public:
        const unsigned int size;
        const int randomSeed;
    public:
        // if get(v) fails, you used all vectors
        // in order to use again all vectors, call reset()
        bool get(orsa::Vector & v) const;
    public:
        void reset() const { 
            rng = new orsa::RNG(randomSeed);
            counter = 0; 
        }
    protected:
        mutable unsigned int counter;
        mutable osg::ref_ptr<orsa::RNG> rng;
    protected:
        static orsa::Vector __randomVectorUtil(const orsa::RNG * rng,
                                               const Box & boundingBox);
    protected:
        std::vector<bool> in;
    };
  
    double volume(const orsa::RandomPointsInShape * randomPointsInShape);
  
    orsa::Vector centerOfMass(const orsa::RandomPointsInShape * randomPointsInShape,
                              const orsa::MassDistribution * massDistribution);
  
    // InertiaMatrix about the centerOfMass, along the PrincipalAxis
    /* void principalAxisAndInertiaMatrix(orsa::Matrix & principalAxis,
       orsa::Matrix & inertiaMatrix,
       const orsa::Vector & centerOfMass,
       const orsa::RandomPointsInShape * randomPointsInShape,
       const orsa::MassDistribution * massDistribution);
    */
    //
    // InertiaMatrix about the centerOfMass
    // local here is the principal axis system
    void diagonalizedInertiaMatrix(orsa::Matrix & shapeToLocal,
                                   orsa::Matrix & localToShape,
                                   orsa::Matrix & inertiaMatrix,
                                   const orsa::Vector & centerOfMass,
                                   const orsa::RandomPointsInShape * randomPointsInShape,
                                   const orsa::MassDistribution * massDistribution);
  
    orsa::PaulMoment * computePaulMoment(const unsigned int order,
                                         const orsa::Matrix & shapeToLocal,
                                         const orsa::Matrix & localToShape,
                                         const orsa::Vector & centerOfMass,
                                         const orsa::RandomPointsInShape * randomPointsInShape,
                                         const orsa::MassDistribution * massDistribution);
  
    // utility, to perform all the computations above...
    void bodyInertialComputations(double & volume,
                                  orsa::Vector & centerOfMass,
                                  orsa::Matrix & shapeToLocal,
                                  orsa::Matrix & localToShape,
                                  orsa::Matrix & inertiaMatrix,
                                  orsa::PaulMoment * * paulMoment,
                                  const unsigned int order,
                                  const orsa::Shape * shape,
                                  const orsa::MassDistribution * massDistribution,
                                  const unsigned int N,
                                  const int randomSeed);
  
} // namespace orsa

#endif // _ORSA_UTIL_
