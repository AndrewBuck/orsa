#ifndef _ORSA_ORBIT_
#define _ORSA_ORBIT_

#include <orsaTBB/malloc.h>

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/datetime.h>
#include <orsa/interval.h>
#include <orsa/vector.h>

namespace orsa {
  
  class Orbit;
  class Body;
  class BodyGroup;  
  
  class Orbit {
  public:   
    Double a,e,i,omega_node,omega_pericenter,M;
  public:
    Double mu;
  public:
    Double period() const {
      // ORSA_DEBUG("a: %Fg   mu: %Fg",a.get_mpf_t(),mu.get_mpf_t());
      return (sqrt(Double("4")*pi()*pi()*fabs(a*a*a)/mu));
    }
  public:
    static Double eccentricAnomaly(const Double & e, const Double & M);
  public:
    Double eccentricAnomaly() const;
  public:
    bool compute(const orsa::Body * b, 
		 const orsa::Body * ref_b, 
		 orsa::BodyGroup  * bg, 
		 const Time       & t);
  public:
    bool compute(const orsa::Vector & relative_position,
		 const orsa::Vector & relative_velocity,
		 const orsa::Double & mu);
  public:
    bool relativePosVel(Vector & relativePosition, Vector & relativeVelocity) const;
    
  public:
    inline bool relativePosition(Vector & relativePosition) const {
      orsa::Vector tmp_v;
      return relativePosVel(relativePosition, tmp_v);
    }
    
    // try to speed up relativePosVel
    // comment this define to disable...
#define _ORBIT_RPV_SPEEDUP_
    
#ifdef _ORBIT_RPV_SPEEDUP_
    /* 
       public:
       Orbit() { }
       public:
       Orbit(const Orbit & o) {
       a                = o.a;
       e                = o.e;
       i                = o.i;
       omega_node       = o.omega_node;
       omega_pericenter = o.omega_pericenter;
       M                = o.M;
       mu               = o.mu;
       //
       _cached_omega_pericenter = o._cached_omega_pericenter;
       _sp                      = o._sp;
       _cp                      = o._cp;
       
       ORSA_DEBUG("copy-op... _cached_omega_pericenter = %Ff   o._cached_omega_pericenter = %Ff",
       _cached_omega_pericenter.getRef().get_mpf_t(),
       o._cached_omega_pericenter.getRef().get_mpf_t());
       }
    */
  protected:
    mutable orsa::Cache<orsa::Double> _cached_omega_pericenter;
    mutable orsa::Double              _sp, _cp;
  protected:
    mutable orsa::Cache<orsa::Double> _cached_omega_node;
    mutable orsa::Double              _so, _co;
  protected:
    mutable orsa::Cache<orsa::Double> _cached_i;
    mutable orsa::Double              _si, _ci;
  protected:    
    mutable orsa::Cache<orsa::Double> _cached_e;
    mutable orsa::Double              _sqe;
  protected:    
    mutable orsa::Cache<orsa::Double> _cached_mu;
    mutable orsa::Cache<orsa::Double> _cached_a;
    mutable orsa::Double              _sqgma;
#endif // _ORBIT_RPV_SPEEDUP_
  };
  
  bool MOID(orsa::Double       & moid,
	    orsa::Double       & M1,
	    orsa::Double       & M2,
	    const orsa::Orbit  & o1,
	    const orsa::Orbit  & o2,
	    const orsa::Double & epsAbs = 1e-6);
  
  Double HillRadius(const Double & a,
		    const Double & m,
		    const Double & M);
  
  //! Hill-radius-based parent body determination
  const Body * HillParentBody(const Body * b,
			      BodyGroup  * bg,
			      const Time & t);
  
  //! The closest body among all the bodies in bg with mass > b::mass...
  const Body * simpleParentBody(const Body * b,
				BodyGroup  * bg,
				const Time & t);
  
  // mainly used for graphical purposes, as 3D visualization
  class OrbitProxy : public osg::Referenced {
  public:
    OrbitProxy(const orsa::Body   * b, 
	       orsa::BodyGroup    * bg,
	       const orsa::Double & accuracy = orsa::Double("0.01"),
	       const orsa::Time   & maxPeriod = orsa::Time(1,0,0,0,0));
    
  protected:
    virtual ~OrbitProxy();
    
  protected:
    class BOT {
    public:
      osg::ref_ptr<const orsa::Body> b; // orbit's reference body 
      orsa::Orbit                    o;
      orsa::Time                     t;
    public:
      inline bool operator == (const BOT & rhs) const {
	return (t == rhs.t);
      }
    public:
      inline bool operator != (const BOT & rhs) const {
	return (t != rhs.t);
      }
    public:
      inline bool operator < (const BOT & rhs) const {
	return (t < rhs.t);
      }
    public:
      inline bool operator > (const BOT & rhs) const {
	return (t > rhs.t);
      }
    public:
      inline bool operator <= (const BOT & rhs) const {
	return (t <= rhs.t);
      }
    public:
      inline bool operator >= (const BOT & rhs) const {
	return (t >= rhs.t);
      }
    };
  protected:
    osg::ref_ptr< orsa::Interval<BOT> > _interval;
    
  public:
    bool getOrbit(orsa::Orbit      & orbit,
		  const orsa::Time & t);
    
  protected:
    virtual bool insertBOT(BOT              & bot,
			   const orsa::Time & t);
    
  protected:
    virtual bool interpolatedBOT(BOT              & bot,
				 const BOT        & bot1,
				 const BOT        & bot2,
				 const orsa::Time & t) const;
    
  protected:
    orsa::Double delta(const BOT & bot1,
		       const BOT & bot2);
    
  protected:   
    osg::ref_ptr<const orsa::Body> _b;
    osg::ref_ptr<orsa::BodyGroup>  _bg;
  protected:
    const orsa::Double             _accuracy;
    const orsa::Time               _maxPeriod;
  };
  
  class OrbitSolution : public orsa::Orbit {
  public:
    orsa::Cache<orsa::Time> epoch;
  public:
    // to be added: covariance matrix
    // for RMS residuals handling, another class can be used....
  };
  
} // namespace orsa

#endif // _ORSA_ORBIT_
