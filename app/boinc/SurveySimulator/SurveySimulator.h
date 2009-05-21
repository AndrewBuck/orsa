#ifndef _SURVEY_SIMULATOR_H_
#define _SURVEY_SIMULATOR_H_

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/datetime.h>
#include <orsa/double.h>
#include <orsa/orbit.h>
#include <orsa/util.h>

#include <map>
#include <list>

#include "boinc_util.h"


class r_Cache {
 public:
  orsa::Cache<orsa::Time> epoch;
 public:
  orsa::Cache<orsa::Vector> r; // absolute position at epochj
 public:
  orsa::Cache<double> vMax; // maximum orbital velocity
 public:
  void update_vMax(const orsa::Orbit & orbit) {
    vMax = sqrt((orbit.mu/orbit.a)*(1+orbit.e)/(1-orbit.e));
  }
};

// magnitude function
// alpha = solar phase angle = angle Sun-Asteroid-Observer
// G = slope parameter (G ~= 0.15)
inline double P (const double & alpha, 
		       const double & G = 0.15) {
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
    // ORSA_DEBUG("-- inserted new entry, size: %i --",entryInterval->size()+1);
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

/**** FunctionProxy: using it for the P(phase,G) function ****/

class PhaseComponentProxyEntry : public FunctionProxyEntry < double, double > {
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

class PhaseComponentProxy : public FunctionProxy <double,double,PhaseComponentProxyEntry> { 
 public:
  PhaseComponentProxy(const double & accuracy) : 
    FunctionProxy<double,double,PhaseComponentProxyEntry>(accuracy) { }
 protected:
  double function(const double & x) const {
    return P(x);
  }	
};

// globaly accessible proxy
extern osg::ref_ptr<PhaseComponentProxy> phaseComponentProxy;

/**** FunctionProxy: using it for the log10(x) function ****/

class Log10ProxyEntry : public FunctionProxyEntry < double, double > {
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

class Log10Proxy : public FunctionProxy <double,double,Log10ProxyEntry> { 
 public:
  Log10Proxy(const double & accuracy) : 
    FunctionProxy<double,double,Log10ProxyEntry>(accuracy) { }
 protected:
  double function(const double & x) const {
    return log10(x);
  }	
};

// globaly accessible proxy
extern osg::ref_ptr<Log10Proxy> log10Proxy;

/*********/

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
  
  // debug
  /* 
     const double nominalLog10 = log10(FromUnits(neo2obs,orsa::Unit::AU,-1)*FromUnits(neo2sun,orsa::Unit::AU,-1));
     ORSA_DEBUG("nominal: %f   proxy: %f   diff: %f",
     nominalLog10,
     proxyLog10,
     fabs(proxyLog10-nominalLog10));
  */
  
  const double V = H + proxyP + 5*proxyLog10;
  
  return V;
}

/* 
   class DiscreteDistributionElement {
   public:
   DiscreteDistributionElement() { }
   public:
   DiscreteDistributionElement(const double & x,
   const double & val) :
   _x(x),
   _val(val) { }
   public:
   orsa::Cache<double> _x;   // left (initial) position
   orsa::Cache<double> _val;
   public:
   bool operator < (const DiscreteDistributionElement & rhs) const {
   return (_x.getRef() < rhs._x.getRef());
   }
   };
*/


class DiscreteDistributionElement {
 public:
  DiscreteDistributionElement() { }
 public:
  DiscreteDistributionElement(const double x,
			      const double val) :
    _x(x),
    _val(val) { }
 public:
  orsa::Cache<double> _x;   // left (initial) position
  orsa::Cache<double> _val;
 public:
  bool operator < (const DiscreteDistributionElement & rhs) const {
    return (_x.getRef() < rhs._x.getRef());
  }
};


// after loading it, you must sort it
typedef std::list<DiscreteDistributionElement> DiscreteDistribution;


class DiscreteDistributionUtility : public osg::Referenced {
 public:  
  // the list must be sorted...
  DiscreteDistributionUtility(const DiscreteDistribution & dd) :
    osg::Referenced(),
    _dd(dd) {
    
    if (dd.size() == 0) {
      ORSA_DEBUG("problems...");
    }
    
    totalVal = 0;
    //
    DiscreteDistribution::const_iterator itLeft  = _dd.begin();
    DiscreteDistribution::const_iterator itRight = _dd.begin(); 
    ++itRight;
    //
    xMin = (*itLeft)._x.getRef();
    //
    while (itRight != _dd.end()) {
      
      const double xLeft  =  (*itLeft)._x.getRef();
      const double xRight = (*itRight)._x.getRef();
      const double val    =  (*itLeft)._val.getRef();
      
      /* 
	 ORSA_DEBUG("considering range: %g to %g  value: %g",
	 xLeft,
	 xRight,
	 val);
      */
      
      if (xLeft > xRight) {
	ORSA_ERROR("using unsorted list");
      }
      
      if (val < 0) {
	ORSA_ERROR("list should contain only positive values or zero");
      }
      
      totalVal += val;
      
      ++itLeft;
      ++itRight;
    }
    //
    xMax = (*itLeft)._x.getRef();
    
    // last value should be zero
    if ((*itLeft)._val.getRef() != 0) {
      ORSA_ERROR("the last value of the list must be zero");
    }
    
    if (xMin.getRef() == xMax.getRef()) {
      ORSA_DEBUG("dubiously singular distribution");
    }
    // ORSA_DEBUG("totalVal: %g",totalVal);
  }
 protected:
  virtual ~DiscreteDistributionUtility() { }
 public:
  double sample(orsa::RNG * rnd) const {
    
    if (0) {
      // debug
      DiscreteDistribution::const_iterator it = _dd.begin();
      while (it != _dd.end()) {
	ORSA_DEBUG("dd: _x=%f _val=%f",(*it)._x.getRef(),(*it)._val.getRef());
	++it;
      }
    }
    
    const double midVal = totalVal.getRef() * rnd->gsl_rng_uniform();
    
    double runningVal = 0;
    //
    DiscreteDistribution::const_iterator itLeft  = _dd.begin();
    DiscreteDistribution::const_iterator itRight = _dd.begin(); 
    ++itRight;
    //
    while (itRight != _dd.end()) {
      
      const double xLeft  =  (*itLeft)._x.getRef();
      const double xRight = (*itRight)._x.getRef();
      const double val    =  (*itLeft)._val.getRef();
      
      /* 
	 ORSA_DEBUG("xLeft: %f xRight: %f val: %f   runningVal: %f runningVal+val: %f   midVal: %f",
	 xLeft,
	 xRight,
	 val,
	 runningVal,
	 runningVal+val,
	 midVal);
      */
      
      if ((   val >  0) &&
	  (midVal >= runningVal) && 
	  (midVal < (runningVal + val))) {
	// ORSA_DEBUG("SELECTED, returning: %f", (xLeft + (xRight - xLeft)*((midVal - runningVal)/val)));
	return (xLeft + (xRight - xLeft)*((midVal - runningVal)/val));
      } else {
	runningVal += val;
      }
      ++itLeft;
      ++itRight;
    }
    ORSA_ERROR("nothing found...");
    return 0;
  }
 public:
  double sample(orsa::RNG    * rnd,
		const double   rangeMin,
		const double   rangeMax) const {
    
    /* 
       ORSA_DEBUG("rangeMin: %f   rangeMax: %f   xMin: %f   xMax: %f",
       rangeMin,
       rangeMax,
       xMin,
       xMax);
    */
    
    if ( (rangeMin < xMin.getRef()) || 
	 (rangeMax > xMax.getRef()) ) {
      ORSA_ERROR("problem: out of boundaries");
      return 0;
    }
    
    if (rangeMin == rangeMax) {
      return rangeMin;
    }
    
    if (rangeMin > rangeMax) {
      ORSA_DEBUG("inverting arguments...");
      return sample(rnd, rangeMax, rangeMin);
    }
    
    if (sampleRangeMin.isSet() &&
	sampleRangeMax.isSet()) {
      if ( (sampleRangeMin.getRef() == rangeMin) &&
	   (sampleRangeMax.getRef() == rangeMax) ) {
	if (sampleRangeDDU.get()) {
	  // ORSA_DEBUG("using cache DDU");
	  return sampleRangeDDU->sample(rnd);
	}
      }
    }
    
    DiscreteDistribution newDD;
    //
    DiscreteDistribution::const_iterator itLeft  = _dd.begin();
    DiscreteDistribution::const_iterator itRight = _dd.begin(); 
    ++itRight;
    //
    while (itRight != _dd.end()) {
      const double xLeft  =  (*itLeft)._x.getRef();
      const double xRight = (*itRight)._x.getRef();
      const double val    =  (*itLeft)._val.getRef();
      
      /* 
	 ORSA_DEBUG("xLeft: %f   xRight: %f   rangeMin: %f   rangeMax: %f",
	 xLeft,
	 xRight,
	 rangeMin,
	 rangeMax);
      */
      
      /* 
	 if (xRight > rangeMin) {
	 ORSA_DEBUG("regular push back...");
	 if (xRight != xLeft) {
	 newDD.push_back(DiscreteDistributionElement(std::max(xLeft,rangeMin),
	 val*(std::min(xRight,rangeMax)-std::max(xLeft,rangeMin))/(xRight-xLeft)));
	 } else {
	 newDD.push_back(DiscreteDistributionElement(std::max(xLeft,rangeMin),
	 val));
	 }
	 if (xRight >= rangeMax) {
	 // last one
	 ORSA_DEBUG("final push back...");
	 newDD.push_back(DiscreteDistributionElement(rangeMax, 0));
	 break;
	 }
	 }	
      */
      //
      if ( (xLeft  < rangeMax) && 
	   (xRight > rangeMin) ) {
	// ORSA_DEBUG("regular push back...");
	if (xRight != xLeft) {
	  newDD.push_back(DiscreteDistributionElement(std::max(xLeft,rangeMin),
						      val*(std::min(xRight,rangeMax)-std::max(xLeft,rangeMin))/(xRight-xLeft)));
	} else {
	  newDD.push_back(DiscreteDistributionElement(std::max(xLeft,rangeMin),
						      val));
	}
	if (xRight >= rangeMax) {
	  // last one
	  // ORSA_DEBUG("final push back...");
	  newDD.push_back(DiscreteDistributionElement(rangeMax, 0));
	  break;
	}
      }	
      
      ++itLeft;
      ++itRight;
    }
    
    // cache for future use
    sampleRangeMin = rangeMin;
    sampleRangeMax = rangeMax;
    sampleRangeDDU = new DiscreteDistributionUtility(newDD);
    
    return sampleRangeDDU->sample(rnd);
  }
 private:
  const DiscreteDistribution _dd;
 private:
  orsa::Cache<double> totalVal;
  orsa::Cache<double> xMin;
  orsa::Cache<double> xMax;
 private:
  mutable orsa::Cache<double> sampleRangeMin, sampleRangeMax;
  mutable osg::ref_ptr<DiscreteDistributionUtility> sampleRangeDDU;
};


class NEO;


class NEOFactory : public osg::Referenced {
 public:
  NEOFactory(const double & a_AU_min_in,
	     const double & a_AU_max_in,
	     const double & e_min_in,
	     const double & e_max_in,
	     const double & i_DEG_min_in,
	     const double & i_DEG_max_in,
	     const int            randomSeed) :
    osg::Referenced(),
    a_AU_min(a_AU_min_in),
    a_AU_max(a_AU_max_in),
    e_min(e_min_in),
    e_max(e_max_in),
    i_DEG_min(i_DEG_min_in),
    i_DEG_max(i_DEG_max_in) {
    
#warning "should read each distribution from external file"
    
    idCounter = 0;
    
    rnd = new orsa::RNG(randomSeed);
    
    // semi-major axis distribution
    DiscreteDistribution dd_a_AU;
    //
    dd_a_AU.push_back(DiscreteDistributionElement(0.0,0.000));
    dd_a_AU.push_back(DiscreteDistributionElement(0.2,0.000));
    dd_a_AU.push_back(DiscreteDistributionElement(0.4,0.003));
    dd_a_AU.push_back(DiscreteDistributionElement(0.6,0.019));
    dd_a_AU.push_back(DiscreteDistributionElement(0.8,0.042));
    dd_a_AU.push_back(DiscreteDistributionElement(1.0,0.062));
    dd_a_AU.push_back(DiscreteDistributionElement(1.2,0.085));
    dd_a_AU.push_back(DiscreteDistributionElement(1.4,0.095));
    dd_a_AU.push_back(DiscreteDistributionElement(1.6,0.107));
    dd_a_AU.push_back(DiscreteDistributionElement(1.8,0.088));
    dd_a_AU.push_back(DiscreteDistributionElement(2.0,0.112));
    dd_a_AU.push_back(DiscreteDistributionElement(2.2,0.128));
    dd_a_AU.push_back(DiscreteDistributionElement(2.4,0.108));
    dd_a_AU.push_back(DiscreteDistributionElement(2.6,0.063));
    dd_a_AU.push_back(DiscreteDistributionElement(2.8,0.050));
    dd_a_AU.push_back(DiscreteDistributionElement(3.0,0.022));
    dd_a_AU.push_back(DiscreteDistributionElement(3.2,0.009));
    dd_a_AU.push_back(DiscreteDistributionElement(3.4,0.008));
    dd_a_AU.push_back(DiscreteDistributionElement(3.6,0.005));
    dd_a_AU.push_back(DiscreteDistributionElement(3.8,0.004));
    dd_a_AU.push_back(DiscreteDistributionElement(4.0,0.002));
    dd_a_AU.push_back(DiscreteDistributionElement(4.2,0.001));
    dd_a_AU.push_back(DiscreteDistributionElement(4.4,0.000));
    //
    dd_a_AU.sort();
    //
    ddu_a_AU = new DiscreteDistributionUtility(dd_a_AU);
    
    // eccentricity distribution
    DiscreteDistribution dd_e;
    //
    dd_e.push_back(DiscreteDistributionElement(0.0,0.005));
    dd_e.push_back(DiscreteDistributionElement(0.1,0.020));
    dd_e.push_back(DiscreteDistributionElement(0.2,0.050));
    dd_e.push_back(DiscreteDistributionElement(0.3,0.084));
    dd_e.push_back(DiscreteDistributionElement(0.4,0.168));
    dd_e.push_back(DiscreteDistributionElement(0.5,0.205));
    dd_e.push_back(DiscreteDistributionElement(0.6,0.190));
    dd_e.push_back(DiscreteDistributionElement(0.7,0.148));
    dd_e.push_back(DiscreteDistributionElement(0.8,0.100));
    dd_e.push_back(DiscreteDistributionElement(0.9,0.030));
    dd_e.push_back(DiscreteDistributionElement(1.0,0.000));
    //
    dd_e.sort();
    //
    ddu_e = new DiscreteDistributionUtility(dd_e);
    
    // inclination distribution
    DiscreteDistribution dd_i_DEG;
    //
    dd_i_DEG.push_back(DiscreteDistributionElement(00.0,0.080));
    dd_i_DEG.push_back(DiscreteDistributionElement(05.0,0.180));
    dd_i_DEG.push_back(DiscreteDistributionElement(10.0,0.162));
    dd_i_DEG.push_back(DiscreteDistributionElement(15.0,0.128));
    dd_i_DEG.push_back(DiscreteDistributionElement(20.0,0.095));
    dd_i_DEG.push_back(DiscreteDistributionElement(25.0,0.076));
    dd_i_DEG.push_back(DiscreteDistributionElement(30.0,0.060));
    dd_i_DEG.push_back(DiscreteDistributionElement(35.0,0.055));
    dd_i_DEG.push_back(DiscreteDistributionElement(40.0,0.040));
    dd_i_DEG.push_back(DiscreteDistributionElement(45.0,0.035));
    dd_i_DEG.push_back(DiscreteDistributionElement(50.0,0.030));
    dd_i_DEG.push_back(DiscreteDistributionElement(55.0,0.023));
    dd_i_DEG.push_back(DiscreteDistributionElement(60.0,0.015));
    dd_i_DEG.push_back(DiscreteDistributionElement(65.0,0.010));
    dd_i_DEG.push_back(DiscreteDistributionElement(70.0,0.000));
    //
    dd_i_DEG.sort();
    //
    ddu_i_DEG = new DiscreteDistributionUtility(dd_i_DEG);
  }
 protected:
  virtual ~NEOFactory() { }
  
 protected:
  virtual osg::ref_ptr<NEO> sampleNEO(const orsa::Time & orbitEpoch) const;
  
 protected:
  const double a_AU_min;
  const double a_AU_max;
  const double e_min;
  const double e_max;
  const double i_DEG_min;
  const double i_DEG_max;
  
 protected:
  osg::ref_ptr<orsa::RNG> rnd;
  
 protected:
  osg::ref_ptr<DiscreteDistributionUtility> ddu_a_AU, ddu_e, ddu_i_DEG;
  
 public:	
  virtual osg::ref_ptr<NEO> createNEO(const unsigned int id,
				      const int randomSeed) const = 0;
 private:
  mutable unsigned int idCounter;
};


// note: NEO.id is not unique, but the combination NEO.id and NEO.randomSeed is unique

class NEO : public osg::Referenced {
 public:  
  NEO(const unsigned int id_in,
      const int random_seed) : 
    osg::Referenced(), 
    id(id_in),
    randomSeed(random_seed) { }
 protected:
  virtual ~NEO() { }
  
 public:
  const unsigned int id;
 public:
  // the original random seed used when this NEO was generated, only for tracking purposes
  const int randomSeed;
  
 public:
  // mutable only because I need to change M on this instead of on a copy,
  // to take advantage of the speedups enabled by _ORBIT_RPV_SPEEDUP_
  mutable orsa::Orbit orbit;
 public:  
  orsa::Time orbitEpoch;
  
 public:
  virtual bool setH(const double & H) = 0;
 public:
  virtual double getH(const orsa::Time   & t,
			    const double & detectionProbabilityThreshold) const = 0;
 public:
  mutable orsa::Cache<r_Cache> r_cache;
  
 public:
  class LogEntry : public osg::Referenced {
  public:
    LogEntry() : osg::Referenced() { }
  protected:
    virtual ~LogEntry() { }
  };
 public:  
  // writeFile should always be true for checkpointing to work, change only for testing or post-processing purposes
  virtual void logInsert(const orsa::Time   & t, 
			 NEO::LogEntry      * e,
			 const bool           writeFile = true) const = 0;
 public:
  std::string logFileName() const {
    char filename[1024];
    snprintf(filename,1024,"NEO.%i.%i.log",randomSeed,id);
    return filename;
  }
 public:
  virtual double detectionProbabilityFromLog (const orsa::Time & t) const = 0;
};


class RealNEO : public NEO {
  
 public:
  class LogEntry : 
  public NEO::LogEntry {
  public:
    //! detection probability at each observation epoch;
    orsa::Cache<double> p;
  public:
    orsa::Cache<std::string> telescopeName;
  }; 
  
 public:
  class NEOFactory : 
  public ::NEOFactory {
  public:
    NEOFactory(const double & a_AU_min,
	       const double & a_AU_max,
	       const double & e_min,
	       const double & e_max,
	       const double & i_DEG_min,
	       const double & i_DEG_max,
	       const int            randomSeed) :
      ::NEOFactory(a_AU_min,
		   a_AU_max,
		   e_min,
		   e_max,
		   i_DEG_min,
		   i_DEG_max,
		   randomSeed) { }
  public:
    osg::ref_ptr<NEO> sampleNEO(const orsa::Time   & orbitEpoch,
				const double & Hmax,
				const double & Ha,
				const double & detectionProbabilityThreshold) const {
      
      while (1) {
	
	osg::ref_ptr<NEO> neo = ::NEOFactory::sampleNEO(orbitEpoch);
	
	const double q = neo->orbit.a*(1-neo->orbit.e);
	const double Q = neo->orbit.a*(1+neo->orbit.e);
	
	neo->setH(Hmax - fabs(rnd->gsl_ran_laplace(Ha)));
	
	// approximate minimum apparent magnitude
	// all in AU already // P(0) = 0 
	/* 
	   const double V_min = 
	   neo->getH(orbitEpoch,detectionProbabilityThreshold) + 
	   5*log10(0.3*1.3); // all in AU already // P(0) = 0 
	*/
	//
	/* const double V_min = 
	   apparentMagnitude(neo->getH(orbitEpoch,detectionProbabilityThreshold),
	   0,
	   orsa::FromUnits(0.3,orsa::Unit::AU),
	   orsa::FromUnits(1.3,orsa::Unit::AU));
	*/
	
	if ( (q < FromUnits(1.3,  orsa::Unit::AU)) &&
	     (Q > FromUnits(0.983,orsa::Unit::AU)) ) {
	  // done!
	  return neo;
	} else {
	  // ORSA_DEBUG("NEO discarded...");
	}	
      }
    }
  public:
    osg::ref_ptr<NEO> createNEO(const unsigned int id,
				const int randomSeed) const {
      return new RealNEO(id,randomSeed);
    }
  };
  
 public:
  RealNEO(const unsigned int id,
	  const int randomSeed) : NEO(id,randomSeed) { 
    trustLogFile=false;
    FILE * fp = fopen(logFileName().c_str(),"r");
    if (fp) {
      mpz_class z;
      double p;
      char telescopeName[1024];
      char line[1024];
      while (fgets(line,1024,fp)) {
	gmp_sscanf(line,"%Zi %lf %s",
		   z.get_mpz_t(),
		   &p,
		   telescopeName);
	const orsa::Time t = orsa::Time(z);
	if (p > 0) {
	  // ORSA_DEBUG("id: %i   t: %Zi   p: %f",id,t.getMuSec().get_mpz_t(),p);
	  RealNEO::LogEntry * e = new RealNEO::LogEntry;
	  //
	  e->p = p;
	  //
	  e->telescopeName = telescopeName;
	  //
	  log[t] = e;
	}
      }
      fclose(fp);
    }
  }
    
 protected:
  orsa::Cache<double> _H; //! absolute magnitude 
 public:
  bool setH(const double & H) {
    _H = H;
    return true;
  }
 public:
  double getH(const orsa::Time   &,
		    const double &) const {
    return _H.getRef();
  }
  
 public:
  typedef std::map<orsa::Time, osg::ref_ptr<RealNEO::LogEntry> > detectionLog;
 protected:
  mutable detectionLog log;
 public:
  const detectionLog & getLog() const {
    return log;
  }
  
 protected:
  mutable bool trustLogFile;
 public:
  void logInsert(const orsa::Time   & t, 
		 NEO::LogEntry      * e,
		 const bool           writeFile = true) const {
    if (!e) return;
    RealNEO::LogEntry * re = dynamic_cast<RealNEO::LogEntry *>(e);
    if (!re) return;
    if (re->p.getRef() > 0) {
      log[t] = re;
      if (writeFile) {
       	if (trustLogFile) {
	  FILE * fp = fopen(logFileName().c_str(),"a");
	  if (fp) {
	    boinc_begin_critical_section();
	    gmp_fprintf(fp,"%Zi %.9f %s"
			MODEOL,
			t.getMuSec().get_mpz_t(),
			re->p.getRef(),
			re->telescopeName.getRef().c_str());
	    fflush(fp);
	    boinc_end_critical_section();
	    fclose(fp);
	  }
	} else {
	  FILE * fp = fopen(logFileName().c_str(),"w");
	  if (fp) {
	    detectionLog::const_iterator it = log.begin();
	    boinc_begin_critical_section();
	    while (it != log.end()) {
	      gmp_fprintf(fp,"%Zi %.9f %s"
			  MODEOL,
			  (*it).first.getMuSec().get_mpz_t(),
			  (*it).second->p.getRef(),
			  (*it).second->telescopeName.getRef().c_str());
	      ++it;
	    }
	    fflush(fp);
	    boinc_end_critical_section();
	    fclose(fp);
	    trustLogFile=true;
	  }
	}
      }
    }
  }
 public:
  double detectionProbabilityFromLog (const orsa::Time & t) const {
    double p = 1;
    detectionLog::const_iterator it = log.begin();
    while (it != log.end()) {
      if ((*it).first <= t) {
	p *= (1 - (*it).second->p.getRef());
      }
      ++it;
    }
    return (1 - p);
  }
  
};


class SyntheticNEO : public NEO {
  
 public:
  class LogEntry : 
  public NEO::LogEntry {
  public:
    orsa::Cache<double> neo2obs;
    orsa::Cache<double> neo2sun;
    orsa::Cache<double> phaseAngle;
  public:
    // data needed to reconstruct the probability curve for the telescope detecting this NEO
    orsa::Cache<double> limitingMagnitude;
  public:
    orsa::Cache<std::string> telescopeName;
  };
  
 public:
  class NEOFactory : 
  public ::NEOFactory {
  public:
    NEOFactory(const double & a_AU_min,
	       const double & a_AU_max,
	       const double & e_min,
	       const double & e_max,
	       const double & i_DEG_min,
	       const double & i_DEG_max,
	       const int            randomSeed) :
      ::NEOFactory(a_AU_min,
		   a_AU_max,
		   e_min,
		   e_max,
		   i_DEG_min,
		   i_DEG_max,
		   randomSeed) { }
  public:
    osg::ref_ptr<NEO> sampleNEO(const orsa::Time & orbitEpoch) const {
      while (1) {
	osg::ref_ptr<NEO> neo = ::NEOFactory::sampleNEO(orbitEpoch);
	const double q = neo->orbit.a*(1-neo->orbit.e);
	const double Q = neo->orbit.a*(1+neo->orbit.e);
	if ( (q < FromUnits(1.3,  orsa::Unit::AU)) &&
	     (Q > FromUnits(0.983,orsa::Unit::AU)) ) {
	  // done
	  return neo;
	}
      }      
    }
  public:
    osg::ref_ptr<NEO> createNEO(const unsigned int id,
				const int randomSeed) const {
      return new SyntheticNEO(id,randomSeed);
    }
  };
  
 public:
  SyntheticNEO(const unsigned int id,
	       const int randomSeed) : NEO(id,randomSeed) { 
    detectionInterval = new orsa::Interval<detectionIntervalEntry>;
    detectionInterval->enableDataStoring();
    trustLogFile=false;
    FILE * fp = fopen(logFileName().c_str(),"r");
    if (fp) {
      mpz_class z;
      double neo2obs, neo2sun, phaseAngle;
      double limitingMagnitude;
      char telescopeName[1024];
      char line[1024];
      while (fgets(line,1024,fp)) {
	gmp_sscanf(line,"%Zi %lf %lf %lf %lf %s",
		   z.get_mpz_t(),
		   &neo2obs,
		   &neo2sun,
		   &phaseAngle,
		   &limitingMagnitude,
		   telescopeName);
	neo2obs     = FromUnits(neo2obs,orsa::Unit::AU);
	neo2sun     = FromUnits(neo2sun,orsa::Unit::AU);
	phaseAngle *= orsa::degToRad();
	
	const orsa::Time t = orsa::Time(z);
	
	SyntheticNEO::LogEntry * e = new SyntheticNEO::LogEntry;
	//
	e->neo2obs    = neo2obs;
	e->neo2sun    = neo2sun;
	e->phaseAngle = phaseAngle;
	//
	e->limitingMagnitude = limitingMagnitude;
	//
	e->telescopeName = telescopeName;	
	//
	log[t] = e;	
	/* 
	   ORSA_DEBUG("id: %i   neo2obs: %5.3f AU   neo2sun: %5.3f AU   phaseAngle: %7.4f DEG   t: %Zi",
	   id,
	   orsa::FromUnits(neo2obs,orsa::Unit::AU,-1),
	   orsa::FromUnits(neo2sun,orsa::Unit::AU,-1),
	   orsa::radToDeg()*phaseAngle,
	   t.getMuSec().get_mpz_t());
	*/
      }
      fclose(fp);
    }
  }
    
 public:
  bool setH(const double &) {
    ORSA_ERROR("cannot set H in a SyntheticNEO");
    return false;
  }
 public:
  double getH(const orsa::Time   & t,
		    const double & detectionProbabilityThreshold) const;
  
 public:
  typedef std::map<orsa::Time, osg::ref_ptr<SyntheticNEO::LogEntry> > detectionLog;
 protected:
  mutable detectionLog log;
 public:
  const detectionLog & getLog() const {
    return log;
  }
  
 protected:
  class detectionIntervalEntry {
  public:	
    orsa::Cache<orsa::Time>   t;
    orsa::Cache<double> H;
  public:
    inline bool operator == (const detectionIntervalEntry & rhs) const { return (t.getRef() == rhs.t.getRef()); }
    inline bool operator != (const detectionIntervalEntry & rhs) const { return (t.getRef() != rhs.t.getRef()); }
    inline bool operator <  (const detectionIntervalEntry & rhs) const { return (t.getRef() <  rhs.t.getRef()); }
    inline bool operator >  (const detectionIntervalEntry & rhs) const { return (t.getRef() >  rhs.t.getRef()); }
    inline bool operator <= (const detectionIntervalEntry & rhs) const { return (t.getRef() <= rhs.t.getRef()); }
    inline bool operator >= (const detectionIntervalEntry & rhs) const { return (t.getRef() >= rhs.t.getRef()); }
  };
 protected:
  typedef orsa::Interval<detectionIntervalEntry> DetectionIntervalType;
 protected:
  osg::ref_ptr<DetectionIntervalType> detectionInterval;
  
 protected:
  mutable bool trustLogFile;  
 public:
  void logInsert(const orsa::Time   & t, 
		 NEO::LogEntry      * e,
		 const bool           writeFile = true) const;
 public:	
  // call logPurify after a series of logInsert calls, to remove all the log entries that do not contribute to H
  void logPurify(const double & detectionProbabilityThreshold) const;
 public:
  double detectionProbabilityFromLog (const orsa::Time &) const {
    ORSA_ERROR("you should not need to call this...");
    return 0;
  }
};


typedef std::list< osg::ref_ptr<NEO> > NEOList;


inline bool NEOListSortPredicate (const osg::ref_ptr<NEO> & lhs, const osg::ref_ptr<NEO> & rhs) {
  if (lhs->id == rhs->id) {
    return (lhs->randomSeed < rhs->randomSeed);
  } else {
    return (lhs->id < rhs->id);
  }
}
//
/* 
   inline bool NEOListSortPredicate (const osg::ref_ptr<NEO> & lhs, const osg::ref_ptr<NEO> & rhs) {
   return ((lhs->id) < (rhs->id));
   }
*/


//! pointing and apparent magnitude
class uVit {
 public:
  orsa::Vector u;
  double V;
 public:
  NEOList::const_iterator it;
};

void dumpSampledNEOs (const NEOList      & neoList,
		      const std::string  & fileName,
		      const double & detectionProbabilityThreshold);

osg::ref_ptr<orsa::Body> SPICEBody (const std::string  & bodyName,
				    const double & bodyMass);

void simpleNEOVector(orsa::Vector & u_obs2neo,
		     double & V,
		     orsa::Vector & neo2obs,
		     orsa::Vector & neo2sun,
		     double & phaseAngle,
		     const NEO * neo,
		     const orsa::Time   & epoch,
	    	     const orsa::Vector & sunPosition,
		     const orsa::Vector & obsPosition,
		     const orsa::Vector & tp_u,
		     const double & apertureAngle,
		     const double & cos_apertureAngle,
		     const double & detectionProbabilityThreshold,
		     const bool           cacheON);

// same refsys as u_obs2neo
inline void vecToSky(double & ra,
		     double & dec,
		     const orsa::Vector & u_obs2neo) {
  ra  = atan2(u_obs2neo.getY(),
	      u_obs2neo.getX());
  dec = orsa::halfpi() - acos(u_obs2neo.getZ());
}

#endif // _SURVEY_SIMULATOR_H_
