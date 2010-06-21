#ifndef SURVEY_REVIEW_H
#define SURVEY_REVIEW_H

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/cache.h>
#include <orsa/double.h>
#include <orsa/orbit.h>
#include <orsa/util.h>

#include <orsaInputOutput/MPC_obscode.h>

#include <orsaSolarSystem/data.h>
#include <orsaSolarSystem/observatory.h>
#include <orsaSolarSystem/gmst.h>
#include <orsaSolarSystem/obleq.h>

#include <orsaSPICE/spice.h>
#include <orsaSPICE/spiceBodyTranslationalCallback.h>

// magnitude function
// alpha = solar phase angle = angle Sun-Asteroid-Observer
// G = slope parameter (G ~= 0.15)
inline double P (const double & alpha, 
                 const double & G) {
    // ORSA_DEBUG("P:   alpha = %f",alpha.get_mpf_t());
    const double phi_1 = exp(-3.33*pow(tan(0.5*alpha),0.63));
    const double phi_2 = exp(-1.87*pow(tan(0.5*alpha),1.22));
    /* 
       ORSA_DEBUG("P = %f   alpha: %f   p1: %f   p2: %f",
       -2.5*log10((1.0-G)*phi_1+G*phi_2),
       alpha.get_mpf_t(),
       phi_1,
       phi_2);
    */
    return (-2.5*log10((1.0-G)*phi_1+G*phi_2));
}

/**** function interpolation, inspired from OrbitProxy ****/

// to be moved into ORSA library when tested and working

/* 
   template <class X, class Y> class FunctionProxyEntry {
   public:
   virtual ~FunctionProxyEntry() { }
   public:
   virtual double delta(const FunctionProxyEntry * e1,
   const FunctionProxyEntry * e2) const = 0;
   public:	
   X x;
   Y y;
   public:
   inline bool operator == (const FunctionProxyEntry & rhs) const { return (x == rhs.x); }
   inline bool operator != (const FunctionProxyEntry & rhs) const { return (x != rhs.x); }
   inline bool operator <  (const FunctionProxyEntry & rhs) const { return (x <  rhs.x); }
   inline bool operator >  (const FunctionProxyEntry & rhs) const { return (x >  rhs.x); }
   inline bool operator <= (const FunctionProxyEntry & rhs) const { return (x <= rhs.x); }
   inline bool operator >= (const FunctionProxyEntry & rhs) const { return (x >= rhs.x); }
   public:
   inline virtual bool interpolatedEntry(FunctionProxyEntry       * e0,      
   const X                  & x0,
   const FunctionProxyEntry * e1,
   const FunctionProxyEntry * e2) const {
   const X beta1 = (e2->x-x0) / (e2->x-e1->x);
   const X beta2 = (x0-e1->x) / (e2->x-e1->x);
   
   e0->x = x0;
   e0->y = beta1*e1->y + beta2*e2->y;
   
   return true;
   }
   };
*/

/* 
   template <class X, class Y, class E> class FunctionProxy : public osg::Referenced {
   public:
   FunctionProxy(const double & accuracy_in) :
   osg::Referenced(),
   accuracy(accuracy_in) {
   if (accuracy <= 0) {
   ORSA_ERROR("non-positive accuracy");
   }
   entryInterval = new orsa::Interval<E>;
   entryInterval->enableDataStoring();
   }
   protected:
   ~FunctionProxy() { }
   public:
   inline virtual bool get(Y & y,
   const X & x) const {
   E e;
   e.x = x;
   if (!entryInterval->size()) {
   if (insert(e,x)) {
   y = e.y;
   return true;
   } else {
   return false;
   }
   }
   if ( (x < entryInterval->min().x) ||
   (x > entryInterval->max().x) )  {
   if (insert(e,x)) {
   y = e.y;
   return true;
   } else {
   return false;
   }
   }
   E eMin;
   E eMax;
   if (!entryInterval->getSubInterval(e,eMin,eMax)) {
   return false;
   }
   if (eMin.x == eMax.x) {
   y = eMin.y;
   return true;
   }
   if (e.delta(&eMin,&eMax) > accuracy) {
   if (insert(e,x)) {
   y = e.y;
   return true;
   } else {
   return false;
   }
   } else {
   if (e.interpolatedEntry(&e,      
   x,
   &eMin,
   &eMax)) {
   y = e.y;
   return true;
   } else {
   return false;
   }
   }
   }
   protected:	
   virtual Y function(const X &) const = 0;
   protected:
   // computes new entry and inserts it into interval
   virtual bool insert(E & e,
   const X & x) const {
   e.x = x;
   e.y = function(x);
   ORSA_DEBUG("-- inserted new entry, size: %i --",entryInterval->size()+1);
   return entryInterval->insert(e);
   }
   protected:
   // virtual osg::ref_ptr<E> createEntry() const { return new E; }
   protected:
   const double accuracy;
   protected:
   // mutable osg::ref_ptr< orsa::Interval< osg::ref_ptr<E> > > entryInterval;
   mutable osg::ref_ptr< orsa::Interval<E> > entryInterval;
   };
*/

/**** FunctionProxy: using it for the P(phase,G) function ****/

/* class PhaseComponentProxyEntry : public FunctionProxyEntry < double, double > {
   public:
   PhaseComponentProxyEntry() : FunctionProxyEntry<double,double>() { }
   public:
   double delta(const FunctionProxyEntry<double,double> * e1,
   const FunctionProxyEntry<double,double> * e2) const {
   const PhaseComponentProxyEntry * p1 = dynamic_cast<const PhaseComponentProxyEntry *> (e1);
   const PhaseComponentProxyEntry * p2 = dynamic_cast<const PhaseComponentProxyEntry *> (e2);
   const double d = fabs((p2->y-p1->y)/(std::min(fabs(p1->y),fabs(p2->y))+orsa::epsilon()));
   return d;
   }
   };
*/

/* 
   class PhaseComponentProxy : public FunctionProxy <double,double,PhaseComponentProxyEntry> { 
   public:
   PhaseComponentProxy(const double & accuracy) : 
   FunctionProxy<double,double,PhaseComponentProxyEntry>(accuracy) { }
   protected:
   double function(const double & x) const {
   return P(x);
   }	
   };
*/

// globaly accessible proxy
// extern osg::ref_ptr<PhaseComponentProxy> phaseComponentProxy;

/**** FunctionProxy: using it for the log10(x) function ****/

/* class Log10ProxyEntry : public FunctionProxyEntry < double, double > {
   public:
   Log10ProxyEntry() : FunctionProxyEntry<double,double>() { }
   public:
   double delta(const FunctionProxyEntry<double,double> * e1,
   const FunctionProxyEntry<double,double> * e2) const {
   const Log10ProxyEntry * p1 = dynamic_cast<const Log10ProxyEntry *> (e1);
   const Log10ProxyEntry * p2 = dynamic_cast<const Log10ProxyEntry *> (e2);
   const double d = fabs((p2->y-p1->y)/(std::min(fabs(p1->y),fabs(p2->y))+orsa::epsilon()));
   return d;
   }
   };
*/

/* 
   class Log10Proxy : public FunctionProxy <double,double,Log10ProxyEntry> { 
   public:
   Log10Proxy(const double & accuracy) : 
   FunctionProxy<double,double,Log10ProxyEntry>(accuracy) { }
   protected:
   double function(const double & x) const {
   return log10(x);
   }	
   };
*/

// globaly accessible proxy
// extern osg::ref_ptr<Log10Proxy> log10Proxy;

/*********/

// choose here if you want to use the functionProxy or not

/* 
   inline double apparentMagnitude(const double & H,
   const double & phaseAngle,
   const double & neo2obs,
   const double & neo2sun) {
   double proxyP;
   if (!phaseComponentProxy->get(proxyP,phaseAngle)) {
   ORSA_DEBUG("problems");
   }
   double proxyLog10;
   if (!log10Proxy->get(proxyLog10,
   FromUnits(neo2obs,orsa::Unit::AU,-1)*FromUnits(neo2sun,orsa::Unit::AU,-1))) {
   ORSA_DEBUG("problems");
   }
   const double V = H + proxyP + 5*proxyLog10;
   return V;
   }
*/
//
inline double apparentMagnitude(const double & H,
                                const double & G,
                                const double & phaseAngle,
                                const double & neo2obs,
                                const double & neo2sun) {
  
    const double V = H + P(phaseAngle,G) + 
        5*log10(FromUnits(neo2obs,orsa::Unit::AU,-1)*FromUnits(neo2sun,orsa::Unit::AU,-1));
  
    return V;
}

/****/

orsa::Body * SPICEBody (const std::string  & bodyName,
                        const double & bodyMass);

class OrbitID : public orsa::Orbit, public osg::Referenced {
public:  
    OrbitID(const mpz_class id_in,
            const orsa::Orbit & earthOrbit_in) :
        // const          int random_seed) : 
        orsa::Orbit(),
        osg::Referenced(),
        id(id_in),
        earthOrbit(earthOrbit_in) { }
protected:
    virtual ~OrbitID() { }
public:
    bool isNEO()    const;
    bool isIEO()    const;
    bool isAten()   const;
    bool isApollo() const;
    bool isAmor()   const;
    bool isPHO()    const;
public:
    const mpz_class id;
    const orsa::Orbit & earthOrbit;
    // double H;
public:
    static const double NEO_max_q, ONE_AU, EARTH_q, EARTH_Q;
};

class OrbitFactory : public osg::Referenced {
public:
    OrbitFactory(const double & a_AU_min_in,
                 const double & a_AU_max_in,
                 const double & e_min_in,
                 const double & e_max_in,
                 const double & i_DEG_min_in,
                 const double & i_DEG_max_in,
                 const double & node_DEG_min_in,
                 const double & node_DEG_max_in,
                 const double & peri_DEG_min_in,
                 const double & peri_DEG_max_in,
                 const double & M_DEG_min_in,
                 const double & M_DEG_max_in,
                 // const double & H_min_in,
                 // const double & H_max_in,
                 const orsa::RNG * rnd_in,
                 const orsa::Orbit & earthOrbit_in) :
        osg::Referenced(),
        a_AU_min(a_AU_min_in),
        a_AU_max(a_AU_max_in),
        e_min(e_min_in),
        e_max(e_max_in),
        i_DEG_min(i_DEG_min_in),
        i_DEG_max(i_DEG_max_in),
        node_DEG_min(node_DEG_min_in),
        node_DEG_max(node_DEG_max_in),
        peri_DEG_min(peri_DEG_min_in),
        peri_DEG_max(peri_DEG_max_in),
        M_DEG_min(M_DEG_min_in),
        M_DEG_max(M_DEG_max_in),
        // H_min(H_min_in),
        // H_max(H_max_in),
        rnd(rnd_in),
        earthOrbit(earthOrbit_in),
        GMSun(orsaSolarSystem::Data::GMSun()) {
        idCounter = 0;
    
        // debug
        /* ORSA_DEBUG("new factory object: %g-%g %g-%g %g-%g %g-%g %g-%g %g-%g %g-%g",
           a_AU_min,
           a_AU_max,
           e_min,
           e_max,
           i_DEG_min,
           i_DEG_max,
           node_DEG_min,
           node_DEG_max,
           peri_DEG_min,
           peri_DEG_max,
           M_DEG_min,
           M_DEG_max,
           H_min,
           H_max);
        */
    }
protected:
    virtual ~OrbitFactory() { }
  
public:
    virtual OrbitID * sample() const;
  
protected:
    const double a_AU_min;
    const double a_AU_max;
    const double e_min;
    const double e_max;
    const double i_DEG_min;
    const double i_DEG_max;
    const double node_DEG_min;
    const double node_DEG_max;
    const double peri_DEG_min;
    const double peri_DEG_max;
    const double M_DEG_min;
    const double M_DEG_max;
    // const double H_min;
    // const double H_max;
    
protected:
    osg::ref_ptr<const orsa::RNG> rnd;
  
protected:
    const orsa::Orbit & earthOrbit;
  
private:
    const double GMSun;
  
private:
    mutable mpz_class idCounter;
};

#endif // SURVEY_REVIEW_H
