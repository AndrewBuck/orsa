#ifndef NICEMODEL_H
#define NICEMODEL_H

#include <orsaTBB/malloc.h>

#include <orsa/body.h>
#include <orsa/bodygroup.h>
#include <orsa/integrator_leapfrog.h>
// #include <orsa/integrator_radau.h>
#include <orsa/orbit.h>
#include <orsa/statistic.h>
#include <orsa/unit.h>
#include <orsa/util.h>

#include <QHash>

class QSIntegrator : public orsa::IntegratorLeapFrog {
  
 public:
  QSIntegrator(const int rs) : orsa::IntegratorLeapFrog(), randomSeed(rs) {
    rnd = new orsa::RNG(randomSeed);
    filesInitialized = false;
    
    stat = new orsa::Statistic<orsa::Double>;
  }
    
 protected:
  ~QSIntegrator() { }
  
 protected:
  const int randomSeed;
  
 public:
  mutable osg::ref_ptr<orsa::Body> sun;
  // mutable osg::ref_ptr<orsa::Body> planet;
  
 private:
  mutable bool filesInitialized;
  
 private:
  void checkFilesInitialization(orsa::BodyGroup * bg) const {
    
    // call this before any fprintf...
    
    if (!filesInitialized) { 
      orsa::BodyGroup::BodyList::const_iterator b_it = bg->getBodyList().begin();
      while (b_it != bg->getBodyList().end()) {
	
	if ((*b_it) == sun.get()) {
	  ++b_it;
	  continue;
	}
	
	// bodyOutputFile.insert((*b_it).get(),fopen(outputFileName((*b_it)->getName(),randomSeed).c_str(),"w"));
	// bodyPerturbationsOutputFile.insert((*b_it).get(),fopen(perturbationsOutputFileName((*b_it)->getName(),randomSeed).c_str(),"w"));
	
	bodyOutputFileName.insert((*b_it).get(),outputFileName((*b_it)->getName(),randomSeed));
	fclose(fopen(bodyOutputFileName[(*b_it).get()].c_str(),"w"));
	
	++b_it;
      }
      
      filesInitialized = true;
    }
    
  }
  
 private:
  mutable orsa::Cache<orsa::Time> lastPerturbTime;
  
 protected:
  osg::ref_ptr< orsa::Statistic<orsa::Double> > stat;
  
  /* 
     public:
     typedef QHash <
     const orsa::Body * , 
     std::FILE * > BodyFile;
     public:
     mutable BodyFile bodyOutputFile;
     mutable BodyFile bodyPerturbationsOutputFile;
  */
  //
 public:
  typedef QHash <
    const orsa::Body * , 
    std::string > BodyFileName;
 public:
  mutable BodyFileName bodyOutputFileName;
  // mutable BodyFile bodyPerturbationsOutputFile;
  
 private:
  std::string outputFileName(const std::string & bodyName, const int randomSeed) const {
    char filename[1024];
    snprintf(filename,1024,"%i.%s.dat",randomSeed,bodyName.c_str());
    return filename;
  }
  /* 
     private:
     std::string perturbationsOutputFileName(const std::string & bodyName, const int randomSeed) const {
     char filename[1024];
     snprintf(filename,1024,"%i.%s.ptb.dat",randomSeed,bodyName.c_str());
     return filename;
     }
  */
  
 public:
  void output(orsa::BodyGroup  * bg,
	      const orsa::Time & t) const {
    
    /* 
       orsa::Orbit orbitPlanet;
       orsa::Vector rp, vp;
       if (bg->getInterpolatedPosVel(rp,vp,planet.get(),t)) {
       if (!orbitPlanet.compute(planet.get(),sun.get(),bg,t)) {
       ORSA_DEBUG("problems...");
       }
       } else {
       ORSA_DEBUG("problems...");
       }
    */
    
    orsa::Vector r, v;
    orsa::BodyGroup::BodyList::const_iterator b_it = bg->getBodyList().begin();
    while (b_it != bg->getBodyList().end()) {
      
      if ((*b_it) == sun.get()) {
	++b_it;
	continue;
      }
      
      if (!(*b_it)->alive(t)) {
	++b_it;
	continue;
      }
      
      if (bg->getInterpolatedPosVel(r,v,(*b_it).get(),t)) {
	
	orsa::Orbit orbit;
	if (!orbit.compute((*b_it).get(),sun.get(),bg,t)) {
	  ORSA_DEBUG("problems");
	}
	
	/* 
	   orsa::Double criticalArgument = 
	   (      orbit.omega_node +       orbit.omega_pericenter +       orbit.M) - 
	   (orbitPlanet.omega_node + orbitPlanet.omega_pericenter + orbitPlanet.M);
	   criticalArgument = orsa::fmod(6*orsa::twopi() + criticalArgument, orsa::twopi());
	   if (criticalArgument > orsa::pi()) criticalArgument -= orsa::twopi(); 
	*/
	
	checkFilesInitialization(bg);
	FILE * fp = fopen(bodyOutputFileName[(*b_it).get()].c_str(),"a");
	gmp_fprintf(fp,
		    // "%14.3Ff %12.9Fe %12.9Fe %10.6Ff %10.6Ff %10.6Ff %10.6Ff %+12.6Ff %12.6Fe %12.6Fe\n",
		    "%14.3Ff %12.9Fe %12.9Fe %10.6Ff %10.6Ff %10.6Ff %10.6Ff\n",
		    FromUnits(t.asDouble(),orsa::Unit::YEAR,-1).get_mpf_t(),
		    orsa::FromUnits(orbit.a,orsa::Unit::AU,-1).get_mpf_t(),
		    orbit.e.get_mpf_t(),
		    orsa::Double(orsa::radToDeg()*orbit.i).get_mpf_t(),
		    orsa::Double(orsa::radToDeg()*orbit.omega_node).get_mpf_t(),  
		    orsa::Double(orsa::radToDeg()*orbit.omega_pericenter).get_mpf_t(),
		    orsa::Double(orsa::radToDeg()*orbit.M).get_mpf_t());
		    // orsa::Double(orsa::radToDeg()*criticalArgument).get_mpf_t(),
		    // orsa::FromUnits(orbit.a-orbitPlanet.a,orsa::Unit::AU,-1).get_mpf_t(),
		    // orsa::FromUnits((r-rp).length(),orsa::Unit::AU,-1).get_mpf_t());
	// fflush(fp);
	fclose(fp);
      } else {
	ORSA_DEBUG("problems");
	return;
      }
      
      ++b_it;
    }
    
  }
  
 protected:
  void cleanup(orsa::BodyGroup  * bg,
	       const orsa::Time & t) const {
    
    // remove all data other not at time = call_t+call_dt
    {
      // ORSA_DEBUG("remember, you're cleaning the BodyGroup intervals...");
      
      // const orsa::Time t = call_t + call_dt;
      
      orsa::BodyGroup::BodyList::const_iterator _b_it = bg->getBodyList().begin();
      while (_b_it != bg->getBodyList().end()) { 
	if (!((*_b_it)->getInitialConditions().translational->dynamic())) { 
	  ++_b_it;
	  continue;
	}
	orsa::BodyGroup::BodyInterval * _b_interval = bg->getBodyInterval((*_b_it).get());
	orsa::BodyGroup::BodyInterval::DataType & _b_interval_data = _b_interval->getData();
	orsa::BodyGroup::BodyInterval::DataType::iterator _b_interval_data_it = _b_interval_data.begin();
	//
#warning "can improve this: by copying the .end() element in a new interval, and wiping out all other elements..."
	//
	while (_b_interval_data_it != _b_interval_data.end()) {
	  if ((*_b_interval_data_it).time.getRef() != t) {
	    _b_interval_data_it = _b_interval_data.erase(_b_interval_data_it);
	  } else {
	    ++_b_interval_data_it;
	  }
	}     
	// IMPORTANT!
	_b_interval->update();
	
	++_b_it;
      }
    }
    
    // "ejected" bodies?
    {
      
      orsa::Vector rSun, vSun;
      bg->getInterpolatedPosVel(rSun,vSun,sun.get(),t);
      
      orsa::Vector r, v;
      orsa::BodyGroup::BodyList::const_iterator b_it = bg->getBodyList().begin();
      while (b_it != bg->getBodyList().end()) {
	
	if ((*b_it) == sun.get()) {
	  ++b_it;
	  continue;
	}
	
	/* 
	   if ((*b_it) == planet.get()) {
	   ++b_it;
	   continue;
	   }
	*/
	
	if (!(*b_it)->alive(t)) {
	  ++b_it;
	  continue;
	}
	
	if (bg->getInterpolatedPosVel(r,v,(*b_it).get(),t)) {
	  if ((r-rSun).length() > orsa::FromUnits(1000,orsa::Unit::AU)) {
	    ORSA_DEBUG("ejecting body [%s]",
		       (*b_it)->getName().c_str());
	    (*b_it)->deathTime = t;
	  }
	}
	
	++b_it;
      }
      
    }
    
  }
  
 public:
  osg::ref_ptr<orsa::RNG> rnd;
  
 protected:
  mutable orsa::Cache<bool> resetNeeded;
  
 protected:    
  bool step(orsa::BodyGroup  * bg,
	    const orsa::Time & start,
	    const orsa::Time & timestep,
	    orsa::Time       & next_timestep) {
    if (resetNeeded.isSet()) {
      if (resetNeeded.getRef()) {
	// ORSA_DEBUG("calling reset()...");
	reset();
	resetNeeded = false;
      }
    }
    return orsa::IntegratorLeapFrog::step(bg,start,timestep,next_timestep);
  }
  
 protected:
  void singleStepDone(orsa::BodyGroup  * bg,
		      const orsa::Time & call_t,
		      const orsa::Time & call_dt,
		      orsa::Time       & ) const {
    
    // output + cleanup
    {
      const unsigned int sample = 64; // 256;
      static unsigned int iter = 0;
      if (iter == sample) {
	output(bg, call_t+call_dt);
	cleanup(bg, call_t+call_dt);
	iter = 0;
      }
      ++iter;
    }
    
  }
  
};

#endif // NICEMODEL_H
