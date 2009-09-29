#ifndef _DAWN_H_
#define _DAWN_H_

#include <orsa/body.h>
#include <orsa/bodygroup.h>
#include <orsa/print.h>
#include <orsa/integrator.h>
#include <orsa/integrator_radau.h>
#include <orsa/multimin.h>
#include <orsa/util.h>

#include <algorithm>
#include <string>

enum SCENARIO {U,
	       C0,
	       CX10,
	       CX20,
	       CX30,
	       CX40,
	       CX50,
	       CZ,
	       C0F20,
	       EU};

inline SCENARIO s2S(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), ::toupper);
  //
  if (s=="U")     return U;
  if (s=="C0")    return C0;
  if (s=="CX10")  return CX10;
  if (s=="CX20")  return CX20;
  if (s=="CX30")  return CX30;
  if (s=="CX40")  return CX40;
  if (s=="CX50")  return CX50;
  if (s=="CZ")    return CZ;
  if (s=="C0F20") return C0F20; 
  if (s=="EU")    return EU; 
  //
  ORSA_DEBUG("problem: could not recognize scenario %s, using U",s.c_str());
  return U;
}

inline std::string S2s (const SCENARIO S) {
  switch (S) {
  case U:     return "U";     break;
  case C0:    return "C0";    break;
  case CX10:  return "CX10";  break;
  case CX20:  return "CX20";  break;
  case CX30:  return "CX30";  break;
  case CX40:  return "CX40";  break;
  case CX50:  return "CX50";  break;
  case CZ:    return "CZ";    break;
  case C0F20: return "C0F20"; break;
  case EU:    return "EU";    break;
  default:    return "undefined"; break;
  }
}

// for a constant mass body
class SRP_and_Engine : public orsa::Propulsion {
 public:
  SRP_and_Engine(const double     & bodyMass_in,
		 const double     & B_in,
		 const double     & thrust_mN_in,
		 orsa::BodyGroup  * bg_in,
		 const orsa::Body * sun_in,
		 const orsa::Body * asteroid_in,
		 const orsa::Body * body_in) :
    orsa::Propulsion(),
    bodyMass(bodyMass_in),
    B(B_in),
    thrust_mN(thrust_mN_in),
    bg(bg_in),
    sun(sun_in),
    asteroid(asteroid_in),
    body(body_in),
    newton(orsa::FromUnits(orsa::FromUnits(orsa::FromUnits(1,orsa::Unit::KG),orsa::Unit::METER),orsa::Unit::SECOND,-2)) {
    
  }
 protected:
  const double bodyMass;
  const double B;
  const double thrust_mN;
  osg::ref_ptr<orsa::BodyGroup>  bg;
  osg::ref_ptr<const orsa::Body> sun;
  osg::ref_ptr<const orsa::Body> asteroid;
  osg::ref_ptr<const orsa::Body> body;
 protected:
  const double newton;
 public:
  orsa::Vector getThrust(const orsa::Time & t) const {
    // compute it once only, use forever...
    // if (!thrust.isSet()) {
    {
      orsa::Vector rSun;
      if (!bg->getInterpolatedPosition(rSun,sun.get(),t)) {
	ORSA_DEBUG("problems...");
      }	
      
      orsa::Vector rAsteroid, vAsteroid;
      if (!bg->getInterpolatedPosVel(rAsteroid,vAsteroid,asteroid.get(),t)) {
	ORSA_DEBUG("problems...");
      }
      orsa::Vector rBody,vBody;
      if (!bg->getInterpolatedPosVel(rBody,vBody,body.get(),t)) {
	ORSA_DEBUG("problems...");
      }	
      
      orsa::Vector thrust(0,0,0);
      
      // first, Solar Radiation Pressure
      {
	// Scheeres (1999) formalism
	const orsa::Vector   s2b = (rBody-rSun);
    	const orsa::Vector u_s2b = s2b.normalized();
	const double G1  = orsa::FromUnits(orsa::FromUnits(orsa::FromUnits(orsa::FromUnits(1.0e8,orsa::Unit::KG),orsa::Unit::KM,3),orsa::Unit::SECOND,-2),orsa::Unit::METER,-2);
	const double R   = s2b.length();
	const double acc = G1/(B*R*R);
       	thrust += (bodyMass*acc*u_s2b);
      }
      
      // then, ION thrust
      /* {
	 const orsa::Vector u_dv = (vBody-vAsteroid).normalized();
	 const double newton = orsa::FromUnits(orsa::FromUnits(orsa::FromUnits(1,orsa::Unit::KG),orsa::Unit::METER),orsa::Unit::SECOND,-2);
	 thrust += (-0.040*newton*u_dv); // nominal: between 0.019 and 0.092 Newton
	 }
      */
      //
      /* {
	 const orsa::Vector   dr =  rBody-rAsteroid;
	 const orsa::Vector u_dr = dr.normalized();
	 const orsa::Vector u_dv = (vBody-vAsteroid).normalized();
	 const double newton = orsa::FromUnits(orsa::FromUnits(orsa::FromUnits(1,orsa::Unit::KG),orsa::Unit::METER),orsa::Unit::SECOND,-2);
	 thrust += (-0.001*thrust_mN*newton*u_dv); // nominal: between 0.019 and 0.092 Newton
	 }
      */
      //
      {
	const orsa::Vector   dr =  rBody-rAsteroid;
	const orsa::Vector u_dr = dr.normalized();
	const orsa::Vector u_dv = (vBody-vAsteroid).normalized();
	const double  target_dr = orsa::FromUnits(420,orsa::Unit::KM); // IMPORTANT: TARGET DISTANCE (little higher than that...)
	if (dr.length() > target_dr) {
	  thrust += (-0.001*thrust_mN*newton*u_dv); // nominal: between 0.019 and 0.092 Newton
	}
      }
      
      return thrust;
    }
    // return thrust.getRef();
  }
 protected:
  // mutable orsa::Cache<orsa::Vector> thrust;
 public:
  bool nextEventTime(orsa::Time      &,
		     const mpz_class &) const {
    return false; // always ON
  }
};

orsa::BodyGroup * run(const double orbitRadius,
		      const SCENARIO scenario,
		      const orsa::Time duration,
		      const double phase_DEG=0,
		      const double thrust_mN=0);

class CustomIntegrator : public orsa::IntegratorRadau {
 public:
  CustomIntegrator(const orsa::Time & t0_in, const bool _verbose=true) :
    orsa::IntegratorRadau(),
    t0(t0_in),
    verbose(_verbose) { }
 protected:
  const orsa::Time t0;
  const bool verbose;
  mutable orsa::Cache<orsa::Time> lastOutput;
 public:
  void singleStepDone(orsa::BodyGroup  * bg,
		      const orsa::Time & call_t,
		      const orsa::Time & call_dt,
		      orsa::Time       & next_dt) const {
    
    if (!verbose) return;
    
    const orsa::Time t = call_t+call_dt;
    
    const orsa::Time minOutputInterval = orsa::Time(0,0,5,0,0);
    bool print=false;
    if (lastOutput.isSet()) {
      if (t>lastOutput.getRef()) {
	if ((t-lastOutput.getRef())>minOutputInterval) {
	  print=true;
	}
      } else {
	if ((lastOutput.getRef()-t)>minOutputInterval) {
	  print=true;
	}
      }
    } else {
      print=true;
    }
    
    // one more chance
    if (next_dt == orsa::Time(0)) {
      print=true;
    }
    
    // force printing, while debugging
    // print=true;
    
    if (!print) return;
    
    lastOutput=t;
    
    osg::ref_ptr<const orsa::Body> sun   = bg->getBody("SUN");
    osg::ref_ptr<const orsa::Body> vesta = bg->getBody("VESTA");
    osg::ref_ptr<const orsa::Body> dawn  = bg->getBody("DAWN");
    
    osg::ref_ptr<const orsa::Shape> vesta_shape = vesta->getInitialConditions().inertial->localShape();
    
    orsa::Vector rDAWN,  vDAWN;
    orsa::Vector rVesta, vVesta;
    orsa::Vector rSun,   vSun;
    //
    // Vector dr, dv, dru, dvu, uPole;
    orsa::Orbit  orbit_rot, orbit_norot;
    orsa::Cache<double> precNode;
    orsa::Cache<orsa::Time> precNodeTime;
    double nodeDot;
    
    if (bg->getInterpolatedPosVel(rDAWN,
				  vDAWN,
				  dawn.get(),
				  t) &&
	bg->getInterpolatedPosVel(rVesta,
				  vVesta,
				  vesta.get(),
				  t) &&
	bg->getInterpolatedPosVel(rSun,
				  vSun,
				  sun.get(),
				  t)) {
      const orsa::Vector dr  = rDAWN-rVesta;
      const orsa::Vector dv  = vDAWN-vVesta;
      //
      const orsa::Vector dr_eclip = dr;
      //
      const orsa::Vector dru = dr.normalized();
      const orsa::Vector dvu = dv.normalized();
      // pole, in global coordinates, unitary vector
      const orsa::Vector uPole = externalProduct(dru,dvu).normalized();
      //
      const double orbitPhaseAngle = fabs(orsa::halfpi() - acos(uPole*(rSun-rVesta).normalized()));
      const double realPhaseAngle  = acos(dru*(rSun-rVesta).normalized());
      
      // fixed time, so no spin rotation is included (fixed body-equatorial frame)
      const orsa::Matrix g2l_norot = orsa::globalToLocal(vesta.get(),bg,t0);
      const orsa::Vector dr_norot  =  g2l_norot * dr;
      const orsa::Vector dv_norot  =  g2l_norot * dv;
      const orsa::Vector dru_norot = (g2l_norot * dru).normalized();
      const orsa::Vector dvu_norot = (g2l_norot * dvu).normalized();
      
      const orsa::Matrix g2l_rot = orsa::globalToLocal(vesta.get(),bg,t);
      const orsa::Vector dr_rot  =  g2l_rot * dr;
      const orsa::Vector dv_rot  =  g2l_rot * dv;
      const orsa::Vector dru_rot = (g2l_rot * dru).normalized();
      const orsa::Vector dvu_rot = (g2l_rot * dvu).normalized();
      
      double vestaMass_t;
      if (!bg->getInterpolatedMass(vestaMass_t,vesta.get(),t)) {
	ORSA_DEBUG("problems...");
      } 
      
      orsa::Vector intersectionPoint;
      orsa::Vector normal;
      if (!vesta_shape->rayIntersection(intersectionPoint,
					normal,
					dr_rot,
					(-dr_rot).normalized(),
					false)) {
	ORSA_DEBUG("problems...");
      } 
      
      // ORSA_DEBUG("remember: globalToLocal!");
      orbit_norot.compute(dr_norot,
			  dv_norot,
			  orsa::Unit::G() * vestaMass_t);
      orbit_rot.compute(dr_rot,
			dv_rot,
			orsa::Unit::G() * vestaMass_t);
      const double latitude  = orsa::halfpi() - acos(dru_rot.getZ());
      const double longitude = fmod(orsa::twopi() + atan2(dru_rot.getY(),dru_rot.getX()),orsa::twopi());
      //
      if (precNode.isSet()) {
	nodeDot = (orbit_norot.omega_node-precNode.getRef())/(t-precNodeTime.getRef()).get_d();
      } else {
	nodeDot = 0.0;
      }
      //
      precNode = orbit_norot.omega_node;
      precNodeTime = t;
      //
      ORSA_DEBUG(
		 /*    1      2      3     4      5      6      7      8      9     10     11     12     13     14      15     16     17     18     19 */
		 "%15.5f %14.3f %12.6f %8.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %12.6f %12.6f %+12.9e %12.6f %12.6f %12.6f\n",
		 orsaSolarSystem::timeToJulian(t),
		 FromUnits((t-t0).get_d(),orsa::Unit::SECOND,-1),
		 FromUnits(orbit_norot.a,orsa::Unit::KM,-1),
		 orbit_norot.e,
		 orsa::radToDeg()*orbit_norot.i,
		 orsa::radToDeg()*orbit_norot.omega_node,
		 orsa::radToDeg()*orbit_norot.omega_pericenter,  
		 orsa::radToDeg()*orbit_norot.M,
		 orsa::FromUnits(orbit_norot.period(),orsa::Unit::HOUR,-1),
		 orsa::radToDeg()*orbitPhaseAngle,
		 orsa::radToDeg()*realPhaseAngle,
		 orsa::radToDeg()*latitude,
		 orsa::radToDeg()*longitude,
		 FromUnits(dr.length(),orsa::Unit::KM,-1),
		 FromUnits(intersectionPoint.length(),orsa::Unit::KM,-1),
		 FromUnits(orsa::radToDeg()*nodeDot,orsa::Unit::HOUR),
		 FromUnits(dr_eclip.getX(),orsa::Unit::KM,-1),
		 FromUnits(dr_eclip.getY(),orsa::Unit::KM,-1),
		 FromUnits(dr_eclip.getZ(),orsa::Unit::KM,-1));
    }
  }
};

#endif // _DAWN_H_
