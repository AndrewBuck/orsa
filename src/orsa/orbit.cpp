#include <orsa/orbit.h>

#include <orsa/bodygroup.h>
#include <orsa/debug.h>
#include <orsa/multimin.h>
#include <orsa/print.h>
#include <orsa/util.h>

using namespace orsa;

// double Orbit::eccentricAnomaly(const double & e, const double & M) {
//
double Orbit::eccentricAnomaly(const double & e, const double & M) {
  
    if (e >= 1) {
        ORSA_WARNING("static orsa::Orbit::eccentricAnomaly(e,M) called with eccentricity = %g (greater than 1.0); returning M.",e);
        //
        return M;
    }
  
    double E = 0.0;
    if (e < 0.8) {
    
        const double sm = sin(M);
        const double cm = cos(M);
    
        // begin with a guess accurate to order e^3 
        double x = M+e*sm*(1+e*(cm+e*(1-1.5*sm*sm)));
    
        double sx,cx;
        E = x;
        double old_E;
        double es,ec,f,fp,fpp,fppp,dx;
    
        unsigned int count = 0;
        const unsigned int max_count = 1024;
        do {
      
            sx = sin(x);
            cx = cos(x);
      
            es = e*sx;
            ec = e*cx;
            f = x - es  - M;
            fp = 1 - ec; 
            fpp = es;
            fppp = ec; 
            dx = -f/fp;
            dx = -f/(fp + dx*fpp/2);
            dx = -f/(fp + dx*fpp/2 + dx*dx*fppp/6);
            //
            old_E = E;
            E = x + dx;
            ++count;
            // update x, ready for the next iteration
            x = E;
      
        } while ((fabs(E-old_E) > 10*(fabs(E)+fabs(M))*epsilon()) && (count < max_count));
    
        if (count >= max_count) {
            ORSA_ERROR("Orbit::eccentricAnomaly(...): max count reached");
            // ORSA_ERROR("Orbit::eccentricAnomaly(): max count reached, e = %g    E = %g   fabs(E-old_E) = %g   10*(fabs(E)+fabs(M))*epsilon() = %g",e,E,fabs(E-old_E),10*(fabs(E)+fabs(M))*std::numeric_limits<double>::epsilon());
        }
    
    } else {
    
        double m = fmod(10*twopi()+fmod(M,twopi()),twopi());
        bool iflag = false;
        if (m > pi()) {
            m = twopi() - m;
            iflag = true;
        }
    
        // Make a first guess that works well for e near 1.  
        double x = cbrt(6*m) - m;
        E = x;
        double old_E;
        double sa,ca,esa,eca,f,fp,dx;
    
        unsigned int count = 0;
        const unsigned int max_count = 128;
        do {
      
            sa = sin(x+m);
            ca = cos(x+m);
      
            esa = e*sa;
            eca = e*ca;
            f = x - esa;
            fp = 1 - eca;
            dx = -f/fp;
            dx = -f/(fp + 0.5*dx*esa);
            dx = -f/(fp + 0.5*dx*(esa+eca*dx/3));
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
            // ORSA_WARNING("Orbit::GetEccentricAnomaly(): max count reached, e = %g    E = %g   fabs(E-old_E) = %g   10*(fabs(E)+fabs(M))*std::numeric_limits<double>::epsilon() = %g",e,E,fabs(E-old_E),10*(fabs(E)+fabs(M))*std::numeric_limits<double>::epsilon());
        }
    }
  
    return E;
}

double Orbit::eccentricAnomaly() const {
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
    double m_b, m_ref_b;
    if (!bg->getInterpolatedMass(m_b,b,t)) {
        ORSA_DEBUG("problems...");
    } 
    if (!bg->getInterpolatedMass(m_ref_b,ref_b,t)) {
        ORSA_DEBUG("problems...");
    }
    //
    if (b->betaSun == ref_b) {
        ORSA_DEBUG("beta-orbit...");
        // mu = (1 - b->beta.getRef())*b->getMu() + ref_b->getMu();
        mu = orsa::Unit::G() * 
            ((1-b->beta.getRef())*m_b + m_ref_b);
    } else {
        // mu = b->getMu()+ref_b->getMu();
        mu = orsa::Unit::G() * 
            (m_b + m_ref_b);
    }
  
    return (compute(relative_position,
                    relative_velocity,
                    mu));
}

bool Orbit::compute(const orsa::Vector & relative_position,
                    const orsa::Vector & relative_velocity,
                    const double & mu_in) {
  
    ////////////////////////////////////////////////////
    // This alghoritm is taken from the swift package, 
    // ORBEL_XV2EL.F file written by M. Duncan.
    ////////////////////////////////////////////////////
  
    mu = mu_in;
  
    // const double tiny = 1.0e-100; // about 4.0e-15
  
    // internals
    double  face,cape,capf,tmpf,cw,sw,w,u;
    int ialpha = 0;
  
    //  ORSA_DEBUG("--MARK--");
    // orsa::print(relative_position);
    // orsa::print(relative_velocity);
  
    // angular momentum
    Vector h = externalProduct(relative_position,relative_velocity);
  
    double h2 = h.lengthSquared();
    double hh = h.length();
  
    // ORSA_DEBUG("hh: %g",hh);
  
    // inclination
    i = acos(h.getZ()/hh);
    // my test...
    // i = asin(sqrt(h.getX()*h.getX()+h.getY()*h.getY())/hh);
  
    // Compute longitude of ascending node omega_node and the argument of latitude u
    // double fac = sqrt(h.getX()*h.getX()+h.getY()*h.getY())/hh;
    double fac = (h.getX()*h.getX()+h.getY()*h.getY())/h2;
  
    // ORSA_DEBUG("fac: %e   h2: %e   eps: %e",fac,h2,epsilon());
  
    if (fac < (epsilon()*epsilon())) {
        // ORSA_DEBUG("--MARK--");
        omega_node = 0;
        u = atan2(relative_position.getY(), relative_position.getX());
        if ( fabs(i-pi()) < epsilon()) u = -u;
    } else {  
        // ORSA_DEBUG("--MARK--");
        omega_node = atan2(h.getX(),-h.getY());
        u = atan2(relative_position.getZ()/sin(i), relative_position.getX()*cos(omega_node)+relative_position.getY()*sin(omega_node));
    }
  
    if (omega_node < 0) omega_node += twopi();
    if (u < 0) u += twopi();
  
    //  Compute the radius r and velocity squared v2, and the dot
    //  product rdotv, the energy per unit mass energy 
    double r  = relative_position.length();
    double v2 = relative_velocity.lengthSquared();
  
    double vdotr  = relative_position*relative_velocity;
  
    double energy = 0.5*v2 - mu/r;
  
    // Determine type of conic section and label it via ialpha
    if (fabs(energy*r/mu) < epsilon()) {
        ialpha = 0;
    } else {
        if (energy < 0) ialpha = -1;
        if (energy > 0) ialpha = +1;
    }
  
    // Depending on the conic type, determine the remaining elements 
  
    // ellipse 
    if (ialpha == -1) {
    
        a   = -mu/(2*energy); 
    
        fac = 1 - h2/(mu*a); 
    
        if (fac > epsilon()) {
            e = sqrt(fac);
            face = (a-r)/(a*e);
      
            if (face > 1) {
                cape = 0;
            } else {
                if (face > -1)
                    cape = acos(face);
                else
                    cape = pi();
            }
      
            if (vdotr < 0) cape = twopi() - cape;
            cw = (cos(cape)-e)/(1-e*cos(cape));                  
            sw = sqrt(1-e*e)*sin(cape)/(1-e*cos(cape));  
            w = atan2(sw,cw);
            if (w < 0) w += twopi();
        } else {
            e = 0;
            w = u;
            cape = u;
        }
    
        M = cape - e*sin(cape);
        omega_pericenter = u - w;
        if (omega_pericenter < 0) omega_pericenter += twopi();
        omega_pericenter = fmod(omega_pericenter,twopi());
    }
  
    // hyperbola
    if (ialpha == 1) {
    
        a   = mu/(2*energy); 
        fac = h2/(mu*a); 
    
        if (fac > epsilon()) {
      
            e = sqrt(1+fac);
            tmpf = (a+r)/(a*e);
            if (tmpf < 1) tmpf = 1;
      
            if (tmpf+sqrt(tmpf*tmpf-1)<0) { ORSA_ERROR("log of negative!"); }
            capf = log(tmpf+sqrt(tmpf*tmpf-1));
      
            if (vdotr < 0) capf = -capf;
      
            cw = (e-cosh(capf))/(e*cosh(capf)-1); 
            sw = sqrt(e*e-1)*sinh(capf)/(e*cosh(capf)-1);
            w  = atan2(sw,cw);
            if (w < 0) w += twopi();
        } else {
            // we only get here if a hyperbola is essentially a parabola 
            // so we calculate e and w accordingly to avoid singularities
            e = 1;
            tmpf = h2/(2*mu); 
            w = acos(2*tmpf/r - 1);
            if (vdotr < 0) w = twopi() - w;
            tmpf = (a+r)/(a*e);
            if (tmpf+sqrt(tmpf*tmpf-1)<0) { ORSA_ERROR("log of negative!"); }
            capf = log(tmpf+sqrt(tmpf*tmpf-1));
        }
    
        M = e * sinh(capf) - capf; 
        omega_pericenter = u - w;
        if (omega_pericenter < 0) omega_pericenter += twopi();
        omega_pericenter = fmod(omega_pericenter,twopi());
    }
  
    // parabola
    //  NOTE - in this case we use "a" to mean pericentric distance
    if (ialpha == 0) {
    
        a = 0.5*h2/mu;
        e = 1;
        w = acos(2*a/r - 1);
        if (vdotr < 0) w = twopi() - w;
        tmpf = tan(w/2);
    
        M = tmpf*(1+tmpf*tmpf/3);
        omega_pericenter = u - w;
        if (omega_pericenter < 0) omega_pericenter += twopi();
        omega_pericenter = fmod(omega_pericenter,twopi());
    }
  
    return true;
}

bool Orbit::relativePosVel(Vector & relativePosition, Vector & relativeVelocity) const {
  
    /////////////////////////////////////////////////////
    // This alghoritm is taken from the swift package, 
    // ORBEL_EL2XV.F file written by M. Duncan.
    /////////////////////////////////////////////////////
  
    double s,c;
  
#ifdef _ORBIT_RPV_SPEEDUP_
  
    // ORSA_DEBUG("-- speedup area --");
  
    if (_cached_omega_pericenter.isSet()) {
        if (_cached_omega_pericenter.getRef() != omega_pericenter) {
            _cached_omega_pericenter = omega_pericenter;
            orsa::sincos(omega_pericenter,&_sp,&_cp);
        } else {
            // ORSA_DEBUG("-- speedup --");
        }
    } else {
        _cached_omega_pericenter = omega_pericenter;
        orsa::sincos(omega_pericenter,&_sp,&_cp);
    } 
    //
    const double sp = _sp;
    const double cp = _cp;
  
    if (_cached_omega_node.isSet()) {
        if (_cached_omega_node.getRef() != omega_node) {
            _cached_omega_node = omega_node;
            orsa::sincos(omega_node,&_so,&_co);
        } else {
            // ORSA_DEBUG("-- speedup --");
        }
    } else {
        _cached_omega_node = omega_node;
        orsa::sincos(omega_node,&_so,&_co);
    } 
    //
    const double so = _so;
    const double co = _co;
  
    if (_cached_i.isSet()) {
        if (_cached_i.getRef() != i) {
            _cached_i = i;
            orsa::sincos(i,&_si,&_ci);
        } else {
            // ORSA_DEBUG("-- speedup --");
        }
    } else {
        _cached_i = i;
        orsa::sincos(i,&_si,&_ci);
    } 
    //
    const double si = _si;
    const double ci = _ci;
  
#else // _ORBIT_RPV_SPEEDUP_
  
    orsa::sincos(omega_pericenter,&s,&c);
    const double sp = s;
    const double cp = c;
  
    orsa::sincos(omega_node,&s,&c);
    const double so = s;
    const double co = c;
  
    orsa::sincos(i,&s,&c);
    const double si = s;
    const double ci = c;
  
#endif // _ORBIT_RPV_SPEEDUP_
  
    const Vector d1(cp*co - sp*so*ci,
                    cp*so + sp*co*ci,
                    sp*si);
  
    const Vector d2(-sp*co - cp*so*ci,
                    -sp*so + cp*co*ci,
                    cp*si);
  
    // Get the other quantities depending on orbit type
  
    // double  cape,scap,ccap,sqe,sqgma,tmp;
    double  cape,tmp;
    double  xfac1,xfac2,vfac1,vfac2;
  
    double  capf,shcap,chcap;
    double  zpara;
  
    if (e < 1) {
    
        cape = eccentricAnomaly();
    
        orsa::sincos(cape,&s,&c);
        const double scap = s;
        const double ccap = c;

#ifdef _ORBIT_RPV_SPEEDUP_
    
        if (_cached_e.isSet()) {
            if (_cached_e.getRef() != e) {
                _cached_e = e;
                _sqe = sqrt(1 - e*e);
            } else {
                // ORSA_DEBUG("-- speedup --");
            }  
        } else {
            _cached_e = e;
            _sqe = sqrt(1 - e*e);
        } 
        //
        const double sqe = _sqe;
    
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
        const double sqgma = _sqgma;
    
#else // _ORBIT_RPV_SPEEDUP_
        const double sqe   = sqrt(1 - e*e);
        const double sqgma = sqrt(fabs(mu*a));
#endif // _ORBIT_RPV_SPEEDUP_

        xfac1 = a*(ccap - e);
        xfac2 = a*sqe*scap;
        // ri = 1/r
        const double ri = 1/(a*(1 - e*ccap));
        vfac1 = -ri * sqgma * scap;
        vfac2 =  ri * sqgma * ccap * sqe;
    
    } else if (e > 1.0) {
    
        double x,shx,chx,esh,ech;
        double f,fp,fpp,fppp,dx;
    
        // use the 'right' value for M -- NEEDED!
        double local_M = M;
        if (fabs(local_M-twopi()) < fabs(local_M)) {
            local_M -= twopi();
        }
    
        // begin with a guess proposed by Danby	
        if (local_M < 0) {
            tmp = -2*local_M/e + 1.8;
            if (tmp<0) { ORSA_ERROR("log of negative!"); }
            x = -log(tmp);
        } else {
            tmp = +2*local_M/e + 1.8;
            if (tmp<0) { ORSA_ERROR("log of negative!"); }
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
            fp = ech - 1; 
            fpp = esh; 
            fppp = ech; 
            dx = -f/fp;
            dx = -f/(fp + dx*fpp/2);
            dx = -f/(fp + dx*fpp/2 + dx*dx*fppp/6);
            capf = x + dx;
            ++count;
            // } while ((fabs(dx) > 1.0e-14) && (count < 100));
        } while ((fabs(dx) > epsilon()) && (count < 256));
    
        shcap = sinh(capf);
        chcap = cosh(capf);
    
        const double sqe   = sqrt(e*e-1);
        const double sqgma = sqrt(fabs(mu*fabs(a)));
        xfac1 = a*(e-chcap);
        xfac2 = a*sqe*shcap;
        const double ri = 1/(a*(e*chcap - 1));
        vfac1 = -ri * sqgma * shcap;
        vfac2 =  ri * sqgma * chcap * sqe;
    
    } else { // e = 1.0 within roundoff errors
    
        double q = M;
        if (q < 1.0e-3) {
            zpara = q*(1 - (q*q/3)*(1-q*q));
        } else {
            double x = 0.5*(3*q+sqrt(9*(q*q)+4));
            double tmp = cbrt(x);
            zpara = tmp - 1/tmp;
        }
    
        const double sqgma = sqrt(fabs(2*mu*a));
        xfac1 = a*(1 - zpara*zpara);
        xfac2 = 2*a*zpara;
        const double ri = 1/(a*(1 + zpara*zpara));
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
    double fun(const orsa::MultiminParameters * par) const {
        o1.M = par->get("M1");
        o2.M = par->get("M2");
        o1.M = fmod(fmod(o1.M,orsa::twopi())+orsa::twopi(),orsa::twopi());
        o2.M = fmod(fmod(o2.M,orsa::twopi())+orsa::twopi(),orsa::twopi());
        orsa::Vector r1, r2;
        o1.relativePosition(r1);
        o2.relativePosition(r2);
        return (r2-r1).length();
    }
public:
    mutable orsa::Orbit o1, o2;
};

bool orsa::MOID(double             & moid,
                double             & M1,
                double             & M2,
                const orsa::Orbit  & o1,
                const orsa::Orbit  & o2,
                const int          & randomSeed,
                const unsigned int & numPoints,
                const double       & epsAbs) {
  
    osg::ref_ptr<orsa::MultiminParameters> par = new MultiminParameters;
    //
    par->insert("M1",o1.M,orsa::arcsecToRad());
    par->insert("M2",o2.M,orsa::arcsecToRad());
  
    osg::ref_ptr<MOID_Multimin> multimin = new MOID_Multimin(o1,o2);
    //
    multimin->setMultiminParameters(par.get());
  
    bool found=false;
  
    {
        osg::ref_ptr<orsa::RNG> rng = new orsa::RNG(randomSeed);  
        for (unsigned int k=0; k<numPoints; ++k) {
            par->set("M1",rng->gsl_rng_uniform()*orsa::twopi());
            par->set("M2",rng->gsl_rng_uniform()*orsa::twopi());
            if (multimin->run_nmsimplex(128,
                                        epsAbs)) {
                if ((found && (multimin->fun(multimin->getMultiminParameters()) < moid)) ||
                    (!found) ) {      
                    moid = multimin->fun(multimin->getMultiminParameters());
                    M1   = multimin->getMultiminParameters()->get("M1");
                    M2   = multimin->getMultiminParameters()->get("M2");
                    found=true;
	  
                    M1 = fmod(fmod(M1,orsa::twopi())+orsa::twopi(),orsa::twopi());
                    M2 = fmod(fmod(M2,orsa::twopi())+orsa::twopi(),orsa::twopi());
	  
                    /* ORSA_DEBUG("moid: %.16f [AU]   M1: %f  M2: %f   k: %i",
                       orsa::FromUnits(moid,orsa::Unit::AU,-1),
                       M1,
                       M2,
                       k);
                    */
                }
            }
        }
    }
  
    /* 
       #warning remove this in production, testing only...
       {
       static osg::ref_ptr<orsa::RNG> rng = new orsa::RNG(23523);
       for (unsigned int k=0; k<256; ++k) {
       par->set("M1",rng->gsl_rng_uniform()*orsa::twopi());
       par->set("M2",rng->gsl_rng_uniform()*orsa::twopi());
       if (multimin->run_nmsimplex(128,
       epsAbs)) {
     
       const double tmp_moid = 1.001 * multimin->fun(par.get());
       if (tmp_moid < moid) {
       ORSA_DEBUG("************ FOUND SMALLER MOID IN STRESS TEST ***********");
       ORSA_DEBUG("moid: %.16f [AU] M1: %f  M2: %f   (deeper k: %i)",
       orsa::FromUnits(tmp_moid,orsa::Unit::AU,-1),
       par->get("M1"), 
       par->get("M2"),
       k);
       exit(0);
       }	
       }
       }
       }
    */	
  
    return found;
}

/**********/

double orsa::HillRadius(const double & a,
                        const double & m,
                        const double & M) {
    // ORSA_DEBUG("check this equation...");
    return (a*cbrt(m/(3*M)));
}

const Body * orsa::HillParentBody(const Body * b,
                                  BodyGroup * bg,
                                  const Time & t) {
  
    ORSA_DEBUG("this function is still incomplete...");
  
    const Body * pb = 0;
    double pb_d_H = -1;
  
    Vector p_b;
    if (!(bg->getInterpolatedPosition(p_b,b,t))) {
        ORSA_DEBUG("problems with BodyGroup::getInterpolatedPosition(...), b: [%s]",
                   b->getName().c_str());
        return (pb);
    }
  
    double m_b, m_b_it;
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
    
        if ((m_b_it == 0) && 
            (m_b    == 0)) {
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
                    double m_parent;
                    if (!bg->getInterpolatedMass(m_parent,parent_b,t)) {
                        ORSA_DEBUG("problems...");
                    }
                    const double a = (p_b_parent-p_b_it).length();
                    const double m = m_b_it;
                    const double M = m_parent;
                    const double r_H = orsa::HillRadius(a,m,M);
                    const double d_H = (p_b-p_b_it).length()/r_H;
                    if ((d_H < pb_d_H) || 
                        (pb_d_H < 0)) {
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
    double pb_a = -1;
  
    Vector r_b, v_b;
    if (!(bg->getInterpolatedPosVel(r_b,v_b,b,t))) {
        ORSA_DEBUG("problems, b: [%s]",
                   b->getName().c_str());
        return (pb);
    }
  
    double m_b, m_b_it;
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
    
        if ((m_b_it == 0) && 
            (m_b    == 0)) {
            ++_b_it;
            // ORSA_DEBUG("continue...");
            continue;
        }
    
        /* 
           if ((*_b_it)->getMass() > b->getMass()) {
           Vector p_b_it;
           if (bg->getInterpolatedPosition(p_b_it,(*_b_it).get(),t)) {
           const double d2 = (p_b-p_b_it).lengthSquared();
           if ((d2 < pb_d2) || 
           (pb_d2 < 0)) {
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
           const double m2d2 = ((*_b_it)->getMass()*(*_b_it)->getMass())/(p_b-p_b_it).lengthSquared();
           if ((m2d2 > pb_m2d2) || 
           (pb_m2d2 < 0)) {
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
                // const double E = 0.5*(v_b-v_b_it).lengthSquared() - (*_b_it)->getMu()/(r_b-r_b_it).length();
                const double E = 0.5*(v_b-v_b_it).lengthSquared() - 
                    orsa::Unit::G() * m_b_it/(r_b-r_b_it).length();
                // const double a = - (*_b_it)->getMu() / (2*E);
                const double a = - orsa::Unit::G() * m_b_it / (2*E);
	
                /* 
                   ORSA_DEBUG("---------------- b: %s  b_it: %s   a: %g    dr: %g  dv: %g  mu: %g",
                   b->getName().c_str(),
                   (*_b_it)->getName().c_str(),
                   a(),
                   (r_b-r_b_it).length(),
                   (v_b-v_b_it).length(),
                   (*_b_it)->getMu());
                */
	
                if ((a > 0) &&
                    ((a < pb_a) || 
                     (pb_a < 0))) {
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
                       const double       & accuracy,
                       const orsa::Time   & maxPeriod) :
    osg::Referenced(true),
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
    
        _interval->insert(bot,false,false);
    
        /* 
           ORSA_DEBUG("inserting, body: [%s] parent: [%s] size: %i   a: %g [km]   e: %g",
           _b->getName().c_str(),
           bot.b->getName().c_str(),
           _interval->size(),
           FromUnits(bot.o.a,Unit::KM,-1)(),
           bot.o.e());
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
  
    const double beta = 
        (t-bot1.t).get_d() / 
        (bot2.t-bot1.t).get_d();
    const double oneMinusBeta = 1-beta;
  
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
    double delta_omega_node = bot2.o.omega_node - bot1.o.omega_node;
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
    double delta_omega_pericenter = bot2.o.omega_pericenter - bot1.o.omega_pericenter;
    if (fabs(bot2.o.omega_pericenter - bot1.o.omega_pericenter + twopi()) < 
        fabs(delta_omega_pericenter)) delta_omega_pericenter = bot2.o.omega_pericenter - bot1.o.omega_pericenter + twopi();
    if (fabs(bot2.o.omega_pericenter - bot1.o.omega_pericenter - twopi()) < 
        fabs(delta_omega_pericenter)) delta_omega_pericenter = bot2.o.omega_pericenter - bot1.o.omega_pericenter - twopi();
    bot.o.omega_pericenter =
        bot1.o.omega_pericenter +
        beta*delta_omega_pericenter;
  
    /* 
       bot.o.M =
       oneMinusBeta*(bot1.o.M + twopi()*(t-bot1.t).get_d()/bot1.o.period()) +
       beta*(bot2.o.M + twopi()*(t-bot2.t).get_d()/bot2.o.period());
    */
    //
    double delta_M = 
        (bot2.o.M-bot1.o.M) + 
        twopi()*((t-bot2.t).get_d()/bot2.o.period()-(t-bot1.t).get_d()/bot1.o.period());
    if (fabs(delta_M + twopi()) < fabs(delta_M)) delta_M += twopi();
    if (fabs(delta_M - twopi()) < fabs(delta_M)) delta_M -= twopi();
    bot.o.M = 
        (bot1.o.M + twopi()*(t-bot1.t).get_d()/bot1.o.period()) +
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
        const double              d    = delta(min, more);
        const double              s    = fabs((min.t-t).get_d()/(more.t-min.t).get_d());
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
        const double              d    = delta(less, max);
        const double              s    = fabs((t-max.t).get_d()/(max.t-less.t).get_d());
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
  
    const double d = delta(botMin,botMax);
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
  
double OrbitProxy::delta(const BOT & bot1,
                         const BOT & bot2) {
    if (bot1.b.get() != bot2.b.get()) return 1;
    double d = 
        fabs((bot2.o.a-bot1.o.a)/(fabs(bot1.o.a)+epsilon())) +
        fabs((bot2.o.e-bot1.o.e)/(fabs(bot1.o.e)+epsilon())) +
        fabs((bot2.o.i-bot1.o.i)/(fabs(bot1.o.i)+epsilon())) +
        fabs((bot2.o.omega_node-bot1.o.omega_node)/(fabs(bot1.o.omega_node)+epsilon())) +
        fabs((bot2.o.omega_pericenter-bot1.o.omega_pericenter)/(fabs(bot1.o.omega_pericenter)+epsilon()));
    // ORSA_DEBUG("delta: %f",d());
    return d;
}

// EquinoctialOrbit

void EquinoctialOrbit::set(const orsa::Orbit & orbit) {
    p = orbit.a*(1-orbit.e*orbit.e);
    f = orbit.e*cos(orbit.omega_node+orbit.omega_pericenter);
    g = orbit.e*sin(orbit.omega_node+orbit.omega_pericenter);
    h = tan(orbit.i/2)*cos(orbit.omega_node);
    k = tan(orbit.i/2)*sin(orbit.omega_node);
    L = orbit.omega_node+orbit.omega_pericenter+orbit.M;
    //
    mu = orbit.mu;
}

void EquinoctialOrbit::get(orsa::Orbit & orbit) const {
    // should check if e==1 -> "a" is undetermined
    orbit.e = sqrt(f*f+g*g);
    orbit.a = p/(1-orbit.e*orbit.e);
    orbit.i = 2*atan(sqrt(h*h+k*k));
    orbit.omega_node = atan2(k,h);
    orbit.omega_pericenter = atan2(g,f)-orbit.omega_node;
    orbit.M = L-orbit.omega_node-orbit.omega_pericenter;
    //
    orbit.mu = mu;
}
