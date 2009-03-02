#ifndef _ORSA_UTIL_
#define _ORSA_UTIL_

#include <string>

#include <orsa/body.h>
#include <orsa/bodygroup.h>
#include <orsa/cache.h>
#include <orsa/datetime.h>
#include <orsa/double.h>
#include <orsa/matrix.h>
#include <orsa/quaternion.h>

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
    IntPowCache(const orsa::Double & x) : _x(x) { }
  public:
    virtual ~IntPowCache() { }
  public:
    /* 
       const orsa::Double get(const int p) const {
       if (!_data[p].isSet()) {
       // this could be more efficient...
       _data[p] = int_pow(_x,p);
       }
       return _data[p].getRef();
       }
    */
    //
    inline const orsa::Double get(const int p) const {      
      if (_data[p].isSet()) {
	return _data[p].getRef();
      } else {
	if (p > 0) {
	  _data[p] = get(p-1) * _x;
	} else if (p < 0) {
	  _data[p] = get(p+1) / _x;
	} else {
	  _data[p] = orsa::one();
	}
	return _data[p].getRef();
      }
    }
  protected:
    const orsa::Double _x;
  protected:
    // mutable QHash <int, orsa::Cache<orsa::Double> > _data;
    mutable std::map <int, orsa::Cache<orsa::Double> > _data;
  };
  
  bool eulerAnglesToMatrix(orsa::Matrix       & m,
			   const orsa::Double & psi,
			   const orsa::Double & theta,
			   const orsa::Double & phi);
  
  bool matrixToEulerAngles(orsa::Double       & psi,
			   orsa::Double       & theta,
			   orsa::Double       & phi,
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
  orsa::Double asteroidDiameter(const orsa::Double & p, 
				const orsa::Double & H);
  
  /* 
     class ConstantZRotation : public PrecomputedRotationalBodyProperty {
     public:	
     ConstantZRotation(const orsa::Time   & t0,
     const orsa::Double & phi0,
     const orsa::Double & omega) : 
     PrecomputedRotationalBodyProperty(),
     _t0(t0), 
     _phi0(phi0),
     _omega(omega) { }
     public:    
     bool update(const orsa::Time & t);
     public:
     bool get(orsa::Double & phi,
     orsa::Double & theta,
     orsa::Double & psi,
     orsa::Double & phiDot,
     orsa::Double & thetaDot,
     orsa::Double & psiDot) const {
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
     orsa::Double getPhi() const { 
     orsa::Double phi,theta,psi,phiDot,thetaDot,psiDot;
     if (!get(phi, theta, psi, phiDot, thetaDot, psiDot)) {
     ORSA_DEBUG("problems...");
     }
     return phi;
     }
     orsa::Double getTheta() const { 
     orsa::Double phi,theta,psi,phiDot,thetaDot,psiDot;
     if (!get(phi, theta, psi, phiDot, thetaDot, psiDot)) {
     ORSA_DEBUG("problems...");
     }
     return theta;
     }
     orsa::Double getPsi() const { 
     orsa::Double phi,theta,psi,phiDot,thetaDot,psiDot;
     if (!get(phi, theta, psi, phiDot, thetaDot, psiDot)) {
     ORSA_DEBUG("problems...");
     }
     return psi;
     }
     orsa::Double getPhiDot()   const { ORSA_DEBUG("implement dot section..."); return zero(); }
     orsa::Double getThetaDot() const { ORSA_DEBUG("implement dot section..."); return zero(); }
     orsa::Double getPsiDot()   const { ORSA_DEBUG("implement dot section..."); return zero(); }
     public:
     RotationalBodyProperty * clone() const {
     return new ConstantZRotation(*this);
     }	
     private:
     const orsa::Time   _t0;
     const orsa::Double _phi0;
     const orsa::Double _omega;
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
    RNG(int random_seed) : 
      osg::Referenced(),
      randomSeed(random_seed) {
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
    double gsl_ran_laplace(double a) const {
      return ::gsl_ran_laplace(rnd,a);
    }    
  public:
    const int randomSeed;
  protected:  
    gsl_rng * rnd;
  };
  
} // namespace orsa

#endif // _ORSA_UTIL_
