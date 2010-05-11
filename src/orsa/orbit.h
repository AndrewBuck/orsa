#ifndef _ORSA_ORBIT_
#define _ORSA_ORBIT_

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
        double a,e,i,omega_node,omega_pericenter,M;
    public:
        double mu;
    public:
        double period() const {
            return (sqrt(4*pisquared()*fabs(a*a*a)/mu));
        }
    public:
        static double eccentricAnomaly(const double & e, const double & M);
    public:
        double eccentricAnomaly() const;
    public:
        bool compute(const orsa::Body * b, 
                     const orsa::Body * ref_b, 
                     orsa::BodyGroup  * bg, 
                     const Time       & t);
    public:
        bool compute(const orsa::Vector & relative_position,
                     const orsa::Vector & relative_velocity,
                     const double & mu);
    public:
        bool relativePosVel(Vector & relativePosition, Vector & relativeVelocity) const;
    
    public:
        inline bool relativePosition(Vector & relativePosition) const {
            orsa::Vector tmp_v;
            return relativePosVel(relativePosition, tmp_v);
        }
    
        // try to speed up relativePosVel
        // comment this define to disable...
        // #define _ORBIT_RPV_SPEEDUP_
    
#ifdef _ORBIT_RPV_SPEEDUP_
    protected:
        mutable orsa::Cache<double> _cached_omega_pericenter;
        mutable double              _sp, _cp;
    protected:
        mutable orsa::Cache<double> _cached_omega_node;
        mutable double              _so, _co;
    protected:
        mutable orsa::Cache<double> _cached_i;
        mutable double              _si, _ci;
    protected:    
        mutable orsa::Cache<double> _cached_e;
        mutable double              _sqe;
    protected:    
        mutable orsa::Cache<double> _cached_mu;
        mutable orsa::Cache<double> _cached_a;
        mutable double              _sqgma;
#endif // _ORBIT_RPV_SPEEDUP_
    };
  
    bool MOID(double             & moid,
              double             & M1,
              double             & M2,
              const orsa::Orbit  & o1,
              const orsa::Orbit  & o2,
              const int          & randomSeed,
              const unsigned int & numPoints = 16,
              const double       & epsAbs    = 1e-6);
    
    double HillRadius(const double & a,
                      const double & m,
                      const double & M);
    
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
                   const double       & accuracy = 0.01,
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
        double delta(const BOT & bot1,
                     const BOT & bot2);
    
    protected:   
        osg::ref_ptr<const orsa::Body> _b;
        osg::ref_ptr<orsa::BodyGroup>  _bg;
    protected:
        const double     _accuracy;
        const orsa::Time _maxPeriod;
    };
  
    class OrbitSolution : public orsa::Orbit {
    public:
        orsa::Cache<orsa::Time> epoch;
    public:
        // to be added: covariance matrix
        // for RMS residuals handling, another class can be used....
    };
  
    class EquinoctialOrbit {
    public:
        /* p = a*(1-e^2)
         * f = e*cos(node+peri)
         * g = e*sin(node+peri)
         * h = tan(i/2)*cos(node)
         * k = tan(i/2)*sin(node)
         * L = node+peri+M
         */
        double p,f,g,h,k,L;
    public:
        double mu;
    public:
        void set(const orsa::Orbit & orbit);
    public:
        void get(orsa::Orbit & orbit) const;
    public:
        double period() const {
            orsa::Orbit orbit;
            get(orbit);
            return orbit.period();
        }
    };
  
} // namespace orsa

#endif // _ORSA_ORBIT_
