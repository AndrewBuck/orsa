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
#include <orsaSPICE/spiceBodyPosVelCallback.h>

orsa::Body * SPICEBody (const std::string  & bodyName,
			const double & bodyMass);

class OrbitID : public orsa::Orbit, public osg::Referenced {
 public:  
  OrbitID(const unsigned int id_in) :
    // const          int random_seed) : 
    orsa::Orbit(),
    osg::Referenced(),
    id(id_in),
    // randomSeed(random_seed),
    NEO_max_q(FromUnits(1.3,  orsa::Unit::AU)),
    ONE_AU(   FromUnits(1.0,  orsa::Unit::AU)),
    EARTH_q(  FromUnits(0.983,orsa::Unit::AU)),
    EARTH_Q(  FromUnits(1.017,orsa::Unit::AU))
      { }
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
  const unsigned int id;
  // const          int randomSeed;
 public:
  double H;
 protected:
  const double NEO_max_q, ONE_AU, EARTH_q, EARTH_Q;
};

class OrbitFactory : public osg::Referenced {
 public:
  OrbitFactory(const double & a_AU_min_in,
	       const double & a_AU_max_in,
	       const double & e_min_in,
	       const double & e_max_in,
	       const double & i_DEG_min_in,
	       const double & i_DEG_max_in,
	       const double & H_min_in,
	       const double & H_max_in,
	       const orsa::RNG * rnd_in) :
    osg::Referenced(),
    a_AU_min(a_AU_min_in),
    a_AU_max(a_AU_max_in),
    e_min(e_min_in),
    e_max(e_max_in),
    i_DEG_min(i_DEG_min_in),
    i_DEG_max(i_DEG_max_in),
    H_min(H_min_in),
    H_max(H_max_in),
    rnd(rnd_in),
    GMSun(orsaSolarSystem::Data::GMSun()) {
    idCounter = 0;

    // debug
    ORSA_DEBUG("new factory object: %g-%g %g-%g %g-%g %g-%g",
	       a_AU_min,
	       a_AU_max,
	       e_min,
	       e_max,
	       i_DEG_min,
	       i_DEG_max,
	       H_min,
	       H_max);
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
  const double H_min;
  const double H_max;
  
 protected:
  osg::ref_ptr<const orsa::RNG> rnd;
  
 private:
  const double GMSun;
  
 private:
  mutable unsigned int idCounter;
};


class StandardObservatoryPositionCallback : public orsaSolarSystem::ObservatoryPositionCallback {
public:
  StandardObservatoryPositionCallback(orsaInputOutput::MPCObsCodeFile * ocf) : 
    orsaSolarSystem::ObservatoryPositionCallback(),
    obsCodeFile(ocf) {
    
    bg = new orsa::BodyGroup;
    
    earth = new orsa::Body;
    //
    earth->setName("EARTH");
    orsaSPICE::SpiceBodyPosVelCallback * sbpvc = new orsaSPICE::SpiceBodyPosVelCallback(earth->getName());
    orsa::IBPS ibps;
    ibps.inertial = new orsa::ConstantMassBodyProperty(orsaSolarSystem::Data::MEarth());
    ibps.translational = sbpvc;
    earth->setInitialConditions(ibps);
    //
    bg->addBody(earth.get());
  }
protected:
  osg::ref_ptr<orsa::Body> earth;
  osg::ref_ptr<orsa::BodyGroup> bg;
  osg::ref_ptr<orsaInputOutput::MPCObsCodeFile> obsCodeFile;
  
public:
  bool getPosition(orsa::Vector      & position,
		   const std::string & obsCode,
		   const orsa::Time  & t) const {
    
    // ORSA_DEBUG("++ obsPos ++    obsCode: %s",obsCode.c_str());
    // orsa::print(t);
    
    orsa::Vector rEarth;
    if (!bg->getInterpolatedPosition(rEarth,earth.get(),t)) { ORSA_DEBUG("problems"); }
    
    // ORSA_DEBUG("earth position: [km]");
    // print(rEarth/orsa::FromUnits(1,orsa::Unit::KM));
    
    /* {
       osg::ref_ptr<orsa::Body> sun = new orsa::Body;
       //
       sun->setName("SUN");
       orsaSPICE::SpiceBodyPosVelCallback * sbpvc = new orsaSPICE::SpiceBodyPosVelCallback(sun->getName());
       orsa::IBPS ibps;
       ibps.inertial = new orsa::ConstantMassBodyProperty(FromUnits(1,orsa::Unit::MSUN));
       ibps.translational = sbpvc;
       sun->setInitialConditions(ibps);
       //
       osg::ref_ptr<orsa::BodyGroup> bg = new orsa::BodyGroup;
       bg->addBody(sun.get());
       orsa::Vector rSun;
       if (!bg->getInterpolatedPosition(rSun,sun.get(),t)) { ORSA_DEBUG("problems"); }
       ORSA_DEBUG("(earth-sun) position: [km]");
       print((rEarth-rSun)/orsa::FromUnits(1,orsa::Unit::KM));
       }
    */
    
    const orsaSolarSystem::Observatory & observatory = obsCodeFile->_data.observatory[obsCode];
    
    // orsa::print(observatory.lon.getRef());
    // orsa::print(observatory.pxy.getRef());
    // orsa::print(observatory.pz.getRef());
    
    double s, c;
    sincos(observatory.lon.getRef(),&s,&c);
    orsa::Vector obsPos(observatory.pxy.getRef()*c,
			observatory.pxy.getRef()*s,
			observatory.pz.getRef());
    // orsa::print(obsPos);
    orsa::Matrix m = orsa::Matrix::identity();
    m.rotZ(orsaSolarSystem::gmst(t));
    // orsa::print(m);
    obsPos = m*obsPos;
    // orsa::print(obsPos);
    obsPos = orsaSolarSystem::equatorialToEcliptic()*obsPos;
    // orsa::print(obsPos);
    
    // ORSA_DEBUG("obsPos relative to Earth [eclip]");
    // print(obsPos);
    
    position = rEarth + obsPos;
    
    // orsa::print(rEarth);
    // orsa::print(position);
    
    return true;
  }
};

#endif // SURVEY_REVIEW_H
