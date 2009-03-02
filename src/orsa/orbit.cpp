#include <orsa/orbit.h>

#include <orsa/bodygroup.h>
#include <orsa/debug.h>
#include <orsa/multimin.h>
#include <orsa/print.h>
#include <orsa/util.h>

using namespace orsa;

// Double Orbit::eccentricAnomaly(const Double & e, const Double & M) {
//
Double Orbit::eccentricAnomaly(const Double & e, const Double & M) {
  
  if (e >= one()) {
    ORSA_WARNING("static orsa::Orbit::eccentricAnomaly(e,M) called with eccentricity = %Fg (greater than 1.0); returning M.",e.get_mpf_t());
    //
    return M;
  }
  
  Double E = 0.0;
  if (e < Double("0.8")) {
    
    const Double sm = sin(M);
    const Double cm = cos(M);
    
    // begin with a guess accurate to order e^3 
    Double x = M+e*sm*(one()+e*(cm+e*(one()-Double("1.5")*sm*sm)));
    
    Double sx,cx;
    E = x;
    Double old_E;
    Double es,ec,f,fp,fpp,fppp,dx;
    
    unsigned int count = 0;
    const unsigned int max_count = 1024;
    do {
      
      sx = sin(x);
      cx = cos(x);
      
      es = e*sx;
      ec = e*cx;
      f = x - es  - M;
      fp = one() - ec; 
      fpp = es;
      fppp = ec; 
      dx = -f/fp;
      dx = -f/(fp + dx*fpp/two());
      dx = -f/(fp + dx*fpp/two() + dx*dx*fppp/Double("6.0"));
      //
      old_E = E;
      E = x + dx;
      ++count;
      // update x, ready for the next iteration
      x = E;
      
    } while ((fabs(E-old_E) > 10*(fabs(E)+fabs(M))*epsilon()) && (count < max_count));
    
    if (count >= max_count) {
      ORSA_ERROR("Orbit::eccentricAnomaly(...): max count reached");
      // ORSA_ERROR("Orbit::eccentricAnomaly(): max count reached, e = %Fg    E = %Fg   fabs(E-old_E) = %Fg   10*(fabs(E)+fabs(M))*epsilon() = %g",e,E,fabs(E-old_E),10*(fabs(E)+fabs(M))*std::numeric_limits<Double>::epsilon());
    }
    
  } else {
    
    Double m = fmod(Double("10")*twopi()+fmod(M,twopi()),twopi());
    bool iflag = false;
    if (m > pi()) {
      m = twopi() - m;
      iflag = true;
    }
    
    // Make a first guess that works well for e near 1.  
    // Double x = pow(Double("6.0")*m,one()/Double("3.0")) - m;
    Double x = cbrt(Double("6.0")*m) - m;
    E = x;
    Double old_E;
    Double sa,ca,esa,eca,f,fp,dx;
    
    unsigned int count = 0;
    const unsigned int max_count = 128;
    do {
      
      sa = sin(x+m);
      ca = cos(x+m);
      
      esa = e*sa;
      eca = e*ca;
      f = x - esa;
      fp = one() - eca;
      dx = -f/fp;
      dx = -f/(fp + 0.5*dx*esa);
      dx = -f/(fp + 0.5*dx*(esa+one()/Double("3.0")*eca*dx));
      x += dx;
      //
      old_E = E;
      E = x + m;
      ++count;
      
    } while ((fabs(E-old_E) > 10*(fabs(E)+fabs(M)+twopi())*epsilon()) && (count < max_count));
    
    if (iflag) {
      E = twopi() - E;
      old_E = twopi() - old_E;
    }
    
    if (count >= max_count) {
      ORSA_WARNING("Orbit::eccentricAnomaly(...): max count reached...");
      // ORSA_WARNING("Orbit::GetEccentricAnomaly(): max count reached, e = %g    E = %g   fabs(E-old_E) = %g   10*(fabs(E)+fabs(M))*std::numeric_limits<Double>::epsilon() = %g",e,E,fabs(E-old_E),10*(fabs(E)+fabs(M))*std::numeric_limits<Double>::epsilon());
    }
  }
  
  return E;
}

Double Orbit::eccentricAnomaly() const {
  return eccentricAnomaly(e,M);
}

bool Orbit::compute(const Body * b, const Body * ref_b, BodyGroup * bg, const Time & t) {
  
  // relative_position, relative_velocity...
  
  if (b == ref_b) {
    return false;
  }
  
  orsa::Vector rb,vb,rrb,vrb;  
  if (!(bg->getInterpolatedPosVel(rb,vb,b,t) &&
	bg->getInterpolatedPosVel(rrb,vrb,ref_b,t) ) ) {
    ORSA_DEBUG("orbit problems");
    return false;
  }
  //
  const Vector relative_position = rb - rrb;
  const Vector relative_velocity = vb - vrb;
  
  // ORSA_DEBUG("--AGAIN--");
  // orsa::print(t);
  // orsa::print(relative_position);
  // orsa::print(relative_velocity);
  
  // mu = b->getMu()+ref_b->getMu();
  // 
  orsa::Double m_b, m_ref_b;
  if (!bg->getInterpolatedMass(m_b,b,t)) {
    ORSA_DEBUG("problems...");
  } 
  if (!bg->getInterpolatedMass(m_ref_b,ref_b,t)) {
    ORSA_DEBUG("problems...");
  }
  //
  if (b->betaSun == ref_b) {
    ORSA_DEBUG("beta-orbit...");
    // mu = (one() - b->beta.getRef())*b->getMu() + ref_b->getMu();
    mu = orsa::Unit::instance()->getG() * 
      ((one()-b->beta.getRef())*m_b + m_ref_b);
  } else {
    // mu = b->getMu()+ref_b->getMu();
    mu = orsa::Unit::instance()->getG() * 
      (m_b + m_ref_b);
  }
  
  return (compute(relative_position,
		  relative_velocity,
		  mu));
}

bool Orbit::compute(const orsa::Vector & relative_position,
		    const orsa::Vector & relative_velocity,
		    const orsa::Double & mu_in) {
  
  ////////////////////////////////////////////////////
  // This alghoritm is taken from the swift package, 
  // ORBEL_XV2EL.F file written by M. Duncan.
  ////////////////////////////////////////////////////
  
  mu = mu_in;
  
  // const Double tiny = 1.0e-100; // about 4.0e-15
  
  // internals
  Double  face,cape,capf,tmpf,cw,sw,w,u;
  int ialpha = 0;
  
  //  ORSA_DEBUG("--MARK--");
  // orsa::print(relative_position);
  // orsa::print(relative_velocity);
  
  // angular momentum
  Vector h = externalProduct(relative_position,relative_velocity);
  
  Double h2 = h.lengthSquared();
  Double hh = h.length();
  
  // ORSA_DEBUG("hh: %Fg",hh.get_mpf_t());
  
  // inclination
  i = acos(h.getZ()/hh);
  // my test...
  // i = asin(sqrt(h.getX()*h.getX()+h.getY()*h.getY())/hh);
  
  // Compute longitude of ascending node omega_node and the argument of latitude u
  // Double fac = secure_sqrt(secure_pow(h.x,2)+secure_pow(h.y,2))/h2;
  Double fac = sqrt(h.getX()*h.getX()+h.getY()*h.getY())/h2;
  
  if (fac < epsilon()) {
    omega_node = zero();
    u = atan2(relative_position.getY(), relative_position.getX());
    if ( fabs(i-pi()) < epsilon()) u = -u;
  } else {
    omega_node = atan2(h.getX(),-h.getY());
    u = atan2(relative_position.getZ()/sin(i), relative_position.getX()*cos(omega_node)+relative_position.getY()*sin(omega_node));
  }
  
  if (omega_node < zero()) omega_node += twopi();
  if (u < zero()) u += twopi();
  
  //  Compute the radius r and velocity squared v2, and the dot
  //  product rdotv, the energy per unit mass energy 
  Double r  = relative_position.length();
  Double v2 = relative_velocity.lengthSquared();
  
  Double vdotr  = relative_position*relative_velocity;
  
  Double energy = v2/two() - mu/r;
  
  // Determine type of conic section and label it via ialpha
  if (fabs(energy*r/mu) < epsilon()) {
    ialpha = 0;
  } else {
    if (energy < zero()) ialpha = -1;
    if (energy > zero()) ialpha = +1;
  }
  
  // Depending on the conic type, determine the remaining elements 
  
  // ellipse 
  if (ialpha == -1) {
    
    a   = -mu/(two()*energy); 
    
    fac = one() - h2/(mu*a); 
    
    if (fac > epsilon()) {
      e = sqrt(fac);
      face = (a-r)/(a*e);
      
      if (face > one()) {
	cape = zero();
      } else {
	if (face > -one())
	  cape = acos(face);
	else
	  cape = pi();
      }
      
      if (vdotr < zero()) cape = twopi() - cape;
      cw = (cos(cape)-e)/(one()-e*cos(cape));                  
      sw = sqrt(one()-e*e)*sin(cape)/(one()-e*cos(cape));  
      w = atan2(sw,cw);
      if (w < zero()) w += twopi();
    } else {
      e = zero();
      w = u;
      cape = u;
    }
    
    M = cape - e*sin(cape);
    omega_pericenter = u - w;
    if (omega_pericenter < zero()) omega_pericenter += twopi();
    omega_pericenter = fmod(omega_pericenter,twopi());
  }
  
  // hyperbola
  if (ialpha == 1) {
    
    a   = mu/(two()*energy); 
    fac = h2/(mu*a); 
    
    if (fac > epsilon()) {
      
      e = sqrt(one()+fac);
      tmpf = (a+r)/(a*e);
      if (tmpf < one()) tmpf = one();
      
      capf = log(tmpf+sqrt(tmpf*tmpf-one()));
      
      if (vdotr < zero()) capf = -capf;
      
      cw = (e-cosh(capf))/(e*cosh(capf)-one()); 
      sw = sqrt(e*e-one())*sinh(capf)/(e*cosh(capf)-one());
      w  = atan2(sw,cw);
      if (w < zero()) w += twopi();
    } else {
      // we only get here if a hyperbola is essentially a parabola 
      // so we calculate e and w accordingly to avoid singularities
      e = one();
      tmpf = h2/(two()*mu); 
      w = acos(two()*tmpf/r - one());
      if (vdotr < zero()) w = twopi() - w;
      tmpf = (a+r)/(a*e);
      capf = log(tmpf+sqrt(tmpf*tmpf-one()));
    }
    
    M = e * sinh(capf) - capf; 
    omega_pericenter = u - w;
    if (omega_pericenter < zero()) omega_pericenter += twopi();
    omega_pericenter = fmod(omega_pericenter,twopi());
  }
  
  // parabola
  //  NOTE - in this case we use "a" to mean pericentric distance
  if (ialpha == 0) {
    
    a = Double("0.5")*h2/mu;
    e = one();
    w = acos(two()*a/r - one());
    if (vdotr < zero()) w = twopi() - w;
    tmpf = tan(w/two());
    
    M = tmpf*(one()+tmpf*tmpf/three());
    omega_pericenter = u - w;
    if (omega_pericenter < zero()) omega_pericenter += twopi();
    omega_pericenter = fmod(omega_pericenter,twopi());
  }
  
  return true;
}

bool Orbit::relativePosVel(Vector & relativePosition, Vector & relativeVelocity) const {
  
  /////////////////////////////////////////////////////
  // This alghoritm is taken from the swift package, 
  // ORBEL_EL2XV.F file written by M. Duncan.
  /////////////////////////////////////////////////////
  
  Double s,c;
  
#ifdef _ORBIT_RPV_SPEEDUP_
  
  // ORSA_DEBUG("-- speedup area --");
  
  if (_cached_omega_pericenter.isSet()) {
    if (_cached_omega_pericenter.getRef() != omega_pericenter) {
      _cached_omega_pericenter = omega_pericenter;
      sincos(omega_pericenter,_sp,_cp);
    } else {
      // ORSA_DEBUG("-- speedup --");
    }
  } else {
    _cached_omega_pericenter = omega_pericenter;
    sincos(omega_pericenter,_sp,_cp);
  } 
  //
  const Double sp = _sp;
  const Double cp = _cp;
  
  if (_cached_omega_node.isSet()) {
    if (_cached_omega_node.getRef() != omega_node) {
      _cached_omega_node = omega_node;
      sincos(omega_node,_so,_co);
    } else {
      // ORSA_DEBUG("-- speedup --");
    }
  } else {
    _cached_omega_node = omega_node;
    sincos(omega_node,_so,_co);
  } 
  //
  const Double so = _so;
  const Double co = _co;
  
  if (_cached_i.isSet()) {
    if (_cached_i.getRef() != i) {
      _cached_i = i;
      sincos(i,_si,_ci);
    } else {
      // ORSA_DEBUG("-- speedup --");
    }
  } else {
    _cached_i = i;
    sincos(i,_si,_ci);
  } 
  //
  const Double si = _si;
  const Double ci = _ci;
  
#else // _ORBIT_RPV_SPEEDUP_
  
  sincos(omega_pericenter,s,c);
  const Double sp = s;
  const Double cp = c;
  
  sincos(omega_node,s,c);
  const Double so = s;
  const Double co = c;
  
  sincos(i,s,c);
  const Double si = s;
  const Double ci = c;
  
#endif // _ORBIT_RPV_SPEEDUP_
  
  const Vector d1(cp*co - sp*so*ci,
		  cp*so + sp*co*ci,
		  sp*si);
  
  const Vector d2(-sp*co - cp*so*ci,
		  -sp*so + cp*co*ci,
		  cp*si);
  
  // Get the other quantities depending on orbit type
  
  // Double  cape,scap,ccap,sqe,sqgma,tmp;
  Double  cape,tmp;
  Double  xfac1,xfac2,vfac1,vfac2;
  
  Double  capf,shcap,chcap;
  Double  zpara;
  
  // ORSA_DEBUG("e: %Ff",e.get_mpf_t());
  
  if (e < one()) {
    
    cape = eccentricAnomaly();
    
    sincos(cape,s,c);
    const Double scap = s;
    const Double ccap = c;

#ifdef _ORBIT_RPV_SPEEDUP_
    
    if (_cached_e.isSet()) {
      if (_cached_e.getRef() != e) {
	_cached_e = e;
	_sqe = sqrt(one() - e*e);
      } else {
	// ORSA_DEBUG("-- speedup --");
      }  
    } else {
      _cached_e = e;
      _sqe = sqrt(one() - e*e);
    } 
    //
    const Double sqe = _sqe;
    
    if (_cached_mu.isSet() && _cached_a.isSet()) {
      if ((_cached_mu.getRef() != mu) || (_cached_a.getRef() != a)) {
	_cached_mu = mu;
	_cached_a  = a;
	_sqgma = sqrt(fabs(mu*a));
      } else {
	// ORSA_DEBUG("-- speedup --");
      }
    } else {
      _cached_mu = mu;
      _cached_a  = a;
      _sqgma = sqrt(fabs(mu*a));
    }
    //    
    const Double sqgma = _sqgma;
    
#else // _ORBIT_RPV_SPEEDUP_
    const Double sqe   = sqrt(one() - e*e);
    const Double sqgma = sqrt(fabs(mu*a));
#endif // _ORBIT_RPV_SPEEDUP_

    xfac1 = a*(ccap - e);
    xfac2 = a*sqe*scap;
    // ri = 1/r
    const Double ri = one()/(a*(one() - e*ccap));
    vfac1 = -ri * sqgma * scap;
    vfac2 =  ri * sqgma * ccap * sqe;
    
  } else if (e > 1.0) {
    
    Double x,shx,chx,esh,ech;
    Double f,fp,fpp,fppp,dx;
    
    // use the 'right' value for M -- NEEDED!
    Double local_M = M;
    if (fabs(local_M-twopi()) < fabs(local_M)) {
      local_M -= twopi();
    }
    
    // begin with a guess proposed by Danby	
    if (local_M < zero()) {
      tmp = -two()*local_M/e + Double("1.8");
      x = -log(tmp);
    } else {
      tmp = +two()*local_M/e + Double("1.8");
      x =  log(tmp);
    }
    
    capf = x;
    
    int count = 0;
    do {
      x = capf;
      shx = sinh(x);
      chx = cosh(x);
      esh = e*shx;
      ech = e*chx;
      f = esh-x-local_M;
      fp = ech - one(); 
      fpp = esh; 
      fppp = ech; 
      dx = -f/fp;
      dx = -f/(fp + dx*fpp/two());
      dx = -f/(fp + dx*fpp/two() + dx*dx*fppp/Double("6.0"));
      capf = x + dx;
      ++count;
      // } while ((fabs(dx) > 1.0e-14) && (count < 100));
    } while ((fabs(dx) > epsilon()) && (count < 256));
    
    shcap = sinh(capf);
    chcap = cosh(capf);
    
    const Double sqe   = sqrt(e*e-one());
    const Double sqgma = sqrt(fabs(mu*fabs(a)));
    xfac1 = a*(e-chcap);
    xfac2 = a*sqe*shcap;
    const Double ri = one()/(a*(e*chcap - one()));
    vfac1 = -ri * sqgma * shcap;
    vfac2 =  ri * sqgma * chcap * sqe;
    
  } else { // e = 1.0 within roundoff errors
    
    Double q = M;
    if (q < 1.0e-3) {
      zpara = q*(one() - (q*q/three())*(one()-q*q));
    } else {
      Double x = Double("0.5")*(three()*q+sqrt(Double("9.0")*(q*q)+Double("4.0")));
      // Double tmp = secure_pow(x,(1.0/3.0));
      Double tmp = cbrt(x);
      zpara = tmp - one()/tmp;
    }
    
    const Double sqgma = sqrt(fabs(two()*mu*a));
    xfac1 = a*(one() - zpara*zpara);
    xfac2 = two()*a*zpara;
    const Double ri = one()/(a*(one() + zpara*zpara));
    vfac1 = -ri * sqgma * zpara;
    vfac2 =  ri * sqgma;
    
  }
  
  relativePosition = d1*xfac1 + d2*xfac2;
  relativeVelocity = d1*vfac1 + d2*vfac2;
  
  return true;
}

/** MOID **/

class MOID_Multimin : public orsa::Multimin {
public:
  MOID_Multimin(const orsa::Orbit & o1_in,
		const orsa::Orbit & o2_in) :
    orsa::Multimin(),
    o1(o1_in),
    o2(o2_in) { }
public:
  orsa::Double fun(const orsa::MultiminParameters * par) const {
    o1.M = par->get("M1");
    o2.M = par->get("M2");
    orsa::Vector r1, r2;
    o1.relativePosition(r1);
    o2.relativePosition(r2);
    return (r2-r1).length();
  }
public:
  mutable orsa::Orbit o1, o2;
};

bool orsa::MOID(orsa::Double       & moid,
		orsa::Double       & M1,
		orsa::Double       & M2,
		const orsa::Orbit  & o1,
		const orsa::Orbit  & o2,
		const orsa::Double & epsAbs) {
  
  ORSA_DEBUG("called...");
  
  osg::ref_ptr<orsa::MultiminParameters> par = new MultiminParameters;
  //
  par->insert("M1",o1.M,1.0*orsa::degToRad());
  par->insert("M2",o2.M,1.0*orsa::degToRad());
  
  osg::ref_ptr<MOID_Multimin> multimin = new MOID_Multimin(o1,o2);
  //
  multimin->setMultiminParameters(par.get());
  
  bool found=false;
  
  // inital tentative solution, using original orbit::M values
  if (multimin->run_nmsimplex(128,
			      epsAbs.get_d())) {
    moid = multimin->fun(multimin->getMultiminParameters());
    M1   = multimin->getMultiminParameters()->get("M1");
    M2   = multimin->getMultiminParameters()->get("M2");
    found=true;
    
    ORSA_DEBUG("moid: %.16Ff [AU]",orsa::FromUnits(moid,orsa::Unit::AU,-1).get_mpf_t());
  }
  
  const unsigned int numTests=5;
  for (unsigned int k=0; k<numTests; ++k) {
    par->set("M1", k*orsa::twopi()/numTests);
    par->set("M2",-k*orsa::twopi()/numTests);
    if (multimin->run_nmsimplex(128,
				epsAbs.get_d())) {
      if ((found && (multimin->fun(multimin->getMultiminParameters()) < moid)) ||
	  (!found) ) {      
	moid = multimin->fun(multimin->getMultiminParameters());
	M1   = multimin->getMultiminParameters()->get("M1");
	M2   = multimin->getMultiminParameters()->get("M2");
	found=true;
	
	ORSA_DEBUG("better value! k: %i *******",k);
	ORSA_DEBUG("moid: %.16Ff [AU]",orsa::FromUnits(moid,orsa::Unit::AU,-1).get_mpf_t());
      }
    }
  }
  
  // #warning remove this in production, testing only...
  /* {
     osg::ref_ptr<orsa::RNG> rng = new orsa::RNG(997345);
     for (unsigned int k=0; k<1024; ++k) {
     par->set("M1",rng->gsl_rng_uniform()*orsa::twopi());
     par->set("M2",rng->gsl_rng_uniform()*orsa::twopi());
     const orsa::Double tmp_moid = multimin->fun(par.get());
     if (tmp_moid < moid) {
     ORSA_DEBUG("************ FOUND SMALLER MOID IN STRESS TEST ***********");
     ORSA_DEBUG("moid: %.16Ff [AU]",orsa::FromUnits(moid,orsa::Unit::AU,-1).get_mpf_t());
     exit(0);
     }	
     }
     }
  */
  
  ORSA_DEBUG("done.");
  
  return found;
}

/**********/

Double orsa::HillRadius(const Double & a,
			const Double & m,
			const Double & M) {
  ORSA_DEBUG("check this equation...");
  return (a*cbrt(m/(three()*M)));
}


const Body * orsa::HillParentBody(const Body * b,
				  BodyGroup * bg,
				  const Time & t) {
  
  ORSA_DEBUG("this function is still incomplete...");
  
  const Body * pb = 0;
  Double pb_d_H = -1;
  
  Vector p_b;
  if (!(bg->getInterpolatedPosition(p_b,b,t))) {
    ORSA_DEBUG("problems with BodyGroup::getInterpolatedPosition(...), b: [%s]",
	       b->getName().c_str());
    return (pb);
  }
  
  orsa::Double m_b, m_b_it;
  if (!bg->getInterpolatedMass(m_b,b,t)) {
    ORSA_DEBUG("problems...");
  }
  
  BodyGroup::BodyList::const_iterator _b_it = bg->getBodyList().begin();
  while (_b_it != bg->getBodyList().end()) {
    if ((*_b_it)==b) {
      ++_b_it;
      // ORSA_DEBUG("continue...");
      continue;
    }
    
    if (!bg->getInterpolatedMass(m_b_it,(*_b_it).get(),t)) {
      ORSA_DEBUG("problems...");
    }
    
    if ((m_b_it == zero()) && 
	(m_b    == zero())) {
      ++_b_it;
      // ORSA_DEBUG("continue...");
      continue;
    }
    
    if (m_b_it > m_b) {
      const Body * parent_b = orsa::HillParentBody((*_b_it).get(),bg,t);
      if (parent_b) {
	Vector p_b_parent, p_b_it;
	if (bg->getInterpolatedPosition(p_b_parent,parent_b,t) &&
	    bg->getInterpolatedPosition(p_b_it,(*_b_it).get(),t)) {
	  orsa::Double m_parent;
	  if (!bg->getInterpolatedMass(m_parent,parent_b,t)) {
	    ORSA_DEBUG("problems...");
	  }
	  const Double a = (p_b_parent-p_b_it).length();
	  const Double m = m_b_it;
	  const Double M = m_parent;
	  const Double r_H = orsa::HillRadius(a,m,M);
	  const Double d_H = (p_b-p_b_it).length()/r_H;
	  if ((d_H < pb_d_H) || 
	      (pb_d_H < zero())) {
	    pb     = (*_b_it).get();
	    pb_d_H = d_H;
	  }
	} else {
	  ORSA_DEBUG("problems with BodyGroup::getInterpolatedPosition(...)");
	  //
	  ++_b_it;
	  // ORSA_DEBUG("continue...");
	  continue;
	}
      } else {
	ORSA_DEBUG("code needed here...");
      }
    }
    
    ++_b_it;
  }
  
  return (pb);
}

const Body * orsa::simpleParentBody(const Body * b,
				    BodyGroup  * bg,
				    const Time & t) {
  
  const Body * pb = 0;
  Double pb_a = -1;
  
  Vector r_b, v_b;
  if (!(bg->getInterpolatedPosVel(r_b,v_b,b,t))) {
    ORSA_DEBUG("problems, b: [%s]",
	       b->getName().c_str());
    return (pb);
  }
  
  orsa::Double m_b, m_b_it;
  if (!bg->getInterpolatedMass(m_b,b,t)) {
    ORSA_DEBUG("problems...");
  }
  
  BodyGroup::BodyList::const_iterator _b_it = bg->getBodyList().begin();
  while (_b_it != bg->getBodyList().end()) {
    if ((*_b_it)==b) {
      ++_b_it;
      // ORSA_DEBUG("continue...");
      continue;
    }
    
    if (!bg->getInterpolatedMass(m_b_it,(*_b_it).get(),t)) {
      ORSA_DEBUG("problems...");
    }
    
    if (!(*_b_it)->alive(t)) {
      ++_b_it;
      // ORSA_DEBUG("continue...");
      continue;
    }
    
    if ((m_b_it == zero()) && 
	(m_b    == zero())) {
      ++_b_it;
      // ORSA_DEBUG("continue...");
      continue;
    }
    
    /* 
       if ((*_b_it)->getMass() > b->getMass()) {
       Vector p_b_it;
       if (bg->getInterpolatedPosition(p_b_it,(*_b_it).get(),t)) {
       const Double d2 = (p_b-p_b_it).lengthSquared();
       if ((d2 < pb_d2) || 
       (pb_d2 < zero())) {
       pb   = (*_b_it).get();
       pb_d2 = d2;
       }
       }
       }
    */
    //
    /* 
       if ((*_b_it)->getMass() > b->getMass()) {
       Vector p_b_it;
       if (bg->getInterpolatedPosition(p_b_it,(*_b_it).get(),t)) {
       const Double m2d2 = ((*_b_it)->getMass()*(*_b_it)->getMass())/(p_b-p_b_it).lengthSquared();
       if ((m2d2 > pb_m2d2) || 
       (pb_m2d2 < zero())) {
       pb   = (*_b_it).get();
       pb_m2d2 = m2d2;
       }
       }
       }
    */
    //
    if (m_b_it > m_b) {
      Vector r_b_it, v_b_it;
      if (bg->getInterpolatedPosVel(r_b_it,v_b_it,(*_b_it).get(),t)) {
	// const Double E = Double("0.5")*(v_b-v_b_it).lengthSquared() - (*_b_it)->getMu()/(r_b-r_b_it).length();
	const Double E = Double("0.5")*(v_b-v_b_it).lengthSquared() - 
	  orsa::Unit::instance()->getG() * m_b_it/(r_b-r_b_it).length();
	// const Double a = - (*_b_it)->getMu() / (2*E);
	const Double a = - orsa::Unit::instance()->getG() * m_b_it / (2*E);
	
	/* 
	   ORSA_DEBUG("---------------- b: %s  b_it: %s   a: %Fg    dr: %Fg  dv: %Fg  mu: %Fg",
	   b->getName().c_str(),
	   (*_b_it)->getName().c_str(),
	   a.get_mpf_t(),
	   (r_b-r_b_it).length().get_mpf_t(),
	   (v_b-v_b_it).length().get_mpf_t(),
	   (*_b_it)->getMu().get_mpf_t());
	*/
	
	if ((a > zero()) &&
	    ((a < pb_a) || 
	     (pb_a < zero()))) {
	  pb   = (*_b_it).get();
	  pb_a = a;
	}
      }
    } 
    
    ++_b_it;
  }
  
  return (pb);
}

////
OrbitProxy::OrbitProxy(const orsa::Body   * b, 
		       orsa::BodyGroup    * bg,
		       const orsa::Double & accuracy,
		       const orsa::Time   & maxPeriod) :
  osg::Referenced(),
  _b(b),
  _bg(bg),
  _accuracy(accuracy),
  _maxPeriod(maxPeriod) {
  if (b==0) ORSA_DEBUG("problems...");
  if (bg==0) ORSA_DEBUG("problems...");
  _interval = new orsa::Interval<BOT>;
  _interval->enableDataStoring();
}

OrbitProxy::~OrbitProxy() { }

bool OrbitProxy::insertBOT(OrbitProxy::BOT  & bot,
			   const orsa::Time & t) {
  
  bot.t = t;
  
  bot.b = orsa::simpleParentBody(_b.get(),
				 _bg.get(),
				 t);
  
  /* 
     if (bot.b.get()) {
     ORSA_DEBUG("b: [%s]   parent: [%s]",
     _b->getName().c_str(),
     bot.b->getName().c_str());
     } else {
     ORSA_DEBUG("-----NO_PARENT_BODY-----");
     }
  */
  
  if (bot.b.get()) {
    
    bot.o.compute(_b.get(),
		  bot.b.get(),
		  _bg.get(),
		  t);
    
    _interval->insert(bot);
    
    /* 
       ORSA_DEBUG("inserting, body: [%s] parent: [%s] size: %i   a: %Fg [km]   e: %Fg",
       _b->getName().c_str(),
       bot.b->getName().c_str(),
       _interval->size(),
       FromUnits(bot.o.a,Unit::KM,-1).get_mpf_t(),
       bot.o.e.get_mpf_t());
    */
    
    return true;
    
  } else {
   
    return false;
  }
}

bool OrbitProxy::interpolatedBOT(OrbitProxy::BOT       & bot,
				 const OrbitProxy::BOT & bot1,
				 const OrbitProxy::BOT & bot2,
				 const orsa::Time      & t) const {
  
  /* 
     ORSA_DEBUG("interpolating, body: [%s] size: %i",
     _b->getName().c_str(),
     _interval->size());
  */
  
  /* 
     if ( (t < std::min(bot1.t,bot2.t)) || 
     (t > std::max(bot1.t,bot2.t)) ||
     (bot1.b.get() != bot2.b.get()) ) {
     return false;
     }
  */

  // t can be outside the limits of bot1 and bot2:
  // that's the whole purpose of the OrbitProxy...
  if (bot1.b.get() != bot2.b.get()) {
    return false;
  }
  
  if ( (bot1.b.get()==0) || 
       (bot2.b.get()==0) ){
    return false;
  }
  
  bot.b = bot1.b;
  
  bot.t = t;
  
  const orsa::Double beta = 
    (t-bot1.t).asDouble() / 
    (bot2.t-bot1.t).asDouble();
  const orsa::Double oneMinusBeta = one()-beta;
  
  bot.o.a =
    oneMinusBeta*bot1.o.a + 
    beta*bot2.o.a;
  
  bot.o.e =
    oneMinusBeta*bot1.o.e + 
    beta*bot2.o.e;
  
  bot.o.i =
    oneMinusBeta*bot1.o.i + 
    beta*bot2.o.i;
  
  /* 
     bot.o.omega_node =
     oneMinusBeta*bot1.o.omega_node + 
     beta*bot2.o.omega_node;
  */
  //
  orsa::Double delta_omega_node = bot2.o.omega_node - bot1.o.omega_node;
  if (fabs(bot2.o.omega_node - bot1.o.omega_node + twopi()) < 
      fabs(delta_omega_node)) delta_omega_node = bot2.o.omega_node - bot1.o.omega_node + twopi();
  if (fabs(bot2.o.omega_node - bot1.o.omega_node - twopi()) < 
      fabs(delta_omega_node)) delta_omega_node = bot2.o.omega_node - bot1.o.omega_node - twopi();
  bot.o.omega_node =
    bot1.o.omega_node +
    beta*delta_omega_node;
  
  /* 
     bot.o.omega_pericenter =
     oneMinusBeta*bot1.o.omega_pericenter + 
     beta*bot2.o.omega_pericenter;
  */
  //
  orsa::Double delta_omega_pericenter = bot2.o.omega_pericenter - bot1.o.omega_pericenter;
  if (fabs(bot2.o.omega_pericenter - bot1.o.omega_pericenter + twopi()) < 
      fabs(delta_omega_pericenter)) delta_omega_pericenter = bot2.o.omega_pericenter - bot1.o.omega_pericenter + twopi();
  if (fabs(bot2.o.omega_pericenter - bot1.o.omega_pericenter - twopi()) < 
      fabs(delta_omega_pericenter)) delta_omega_pericenter = bot2.o.omega_pericenter - bot1.o.omega_pericenter - twopi();
  bot.o.omega_pericenter =
    bot1.o.omega_pericenter +
    beta*delta_omega_pericenter;
  
  /* 
     bot.o.M =
     oneMinusBeta*(bot1.o.M + twopi()*(t-bot1.t).asDouble()/bot1.o.period()) +
     beta*(bot2.o.M + twopi()*(t-bot2.t).asDouble()/bot2.o.period());
  */
  //
  orsa::Double delta_M = 
    (bot2.o.M-bot1.o.M) + 
    twopi()*((t-bot2.t).asDouble()/bot2.o.period()-(t-bot1.t).asDouble()/bot1.o.period());
  if (fabs(delta_M + twopi()) < fabs(delta_M)) delta_M += twopi();
  if (fabs(delta_M - twopi()) < fabs(delta_M)) delta_M -= twopi();
  bot.o.M = 
    (bot1.o.M + twopi()*(t-bot1.t).asDouble()/bot1.o.period()) +
    beta*delta_M;
  
  bot.o.mu =
    oneMinusBeta*bot1.o.mu + 
    beta*bot2.o.mu;
  
  // normalize
  bot.o.omega_node       = fmod(twopi()+fmod(bot.o.omega_node,
					     twopi()),twopi());
  bot.o.omega_pericenter = fmod(twopi()+fmod(bot.o.omega_pericenter,
					     twopi()),twopi());
  bot.o.M                = fmod(twopi()+fmod(bot.o.M,
					     twopi()),twopi());
  
  return true;
}

bool OrbitProxy::getOrbit(orsa::Orbit      & orbit,
			  const orsa::Time & t) {
  
  // ORSA_DEBUG("b: [%s]",_b->getName().c_str());
  
  OrbitProxy::BOT bot;
  bot.t = t;
  
  if (_interval->size() < 2) {
    if (OrbitProxy::insertBOT(bot,t)) {
      orbit = bot.o;
      return true;
    } else {
      return false;
    }
  }
  
  if (t < _interval->min().t) {
    const Interval<BOT>::DataType & data = _interval->getData();
    const OrbitProxy::BOT           min  = *(data.begin());
    const OrbitProxy::BOT           more = *(++(data.begin()));
    const orsa::Double              d    = delta(min, more);
    const orsa::Double              s    = fabs((min.t-t).asDouble()/(more.t-min.t).asDouble());
    if ( (d*s < _accuracy) && 
	 (interpolatedBOT(bot,
			  min,
			  more,
			  t)) ) {
      orbit = bot.o;
      return true;
    } else if (OrbitProxy::insertBOT(bot,t)) {
      orbit = bot.o;
      return true;
    } else {
      return false;
    }
  } else if (t > _interval->max().t) {
    const Interval<BOT>::DataType & data = _interval->getData();
    const OrbitProxy::BOT           max  = *(--data.end());
    const OrbitProxy::BOT           less = *(--(--data.end()));
    const orsa::Double              d    = delta(less, max);
    const orsa::Double              s    = fabs((t-max.t).asDouble()/(max.t-less.t).asDouble());
    if ( (d*s < _accuracy) &&
	 (interpolatedBOT(bot,
			  less,
			  max,
			  t)) ) {
      orbit = bot.o;
      return true;
    } else if (OrbitProxy::insertBOT(bot,t)) {
      orbit = bot.o;
      return true;
    } else {
      return false;
    }
  }
  
  /* 
     if ( (t < _interval->min().t) || 
     (t > _interval->max().t) ) {
     if (OrbitProxy::insertBOT(bot,t)) {
     orbit = bot.o;
     return true;
     } else {
     return false;
     }
     }
  */
  
  OrbitProxy::BOT botMin, botMax;
  const bool sub = _interval->getSubInterval(bot,botMin,botMax);
  
  if (!sub) {
    ORSA_ERROR("problems, b: [%s] size: %i",
	       _b->getName().c_str(),
	       _interval->size());
    return false;
  }
  
  if ((botMax.t-botMin.t) > _maxPeriod) {
    if (OrbitProxy::insertBOT(bot,t)) {
      orbit = bot.o;
      return true;
    } else {
      return false;
    }
  }  
  
  if (botMax.t == botMin.t) {
    orbit = botMin.o;
    return true;
  }
  
  const orsa::Double d = delta(botMin,botMax);
  if (d > _accuracy) {
    if (OrbitProxy::insertBOT(bot,t)) {
      orbit = bot.o;
      return true;
    } else {
      return false;
    }
  } else {
    if (interpolatedBOT(bot,
			botMin,
			botMax,
			t)) {
      orbit = bot.o;
      return true;
    } else {
      return false;
    }
  }  
}
  
orsa::Double OrbitProxy::delta(const BOT & bot1,
			       const BOT & bot2) {
  if (bot1.b.get() != bot2.b.get()) return one();
  Double d = 
    fabs((bot2.o.a-bot1.o.a)/(fabs(bot1.o.a)+epsilon())) +
    fabs((bot2.o.e-bot1.o.e)/(fabs(bot1.o.e)+epsilon())) +
    fabs((bot2.o.i-bot1.o.i)/(fabs(bot1.o.i)+epsilon())) +
    fabs((bot2.o.omega_node-bot1.o.omega_node)/(fabs(bot1.o.omega_node)+epsilon())) +
    fabs((bot2.o.omega_pericenter-bot1.o.omega_pericenter)/(fabs(bot1.o.omega_pericenter)+epsilon()));
  // ORSA_DEBUG("delta: %Ff",d.get_mpf_t());
  return d;
}
