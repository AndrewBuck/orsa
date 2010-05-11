#include <orsaUtil/orbit.h>

#include <orsaSolarSystem/obleq.h>

double orsaUtil::RMS(const std::vector<orsaUtil::Residual> & r) {
    // ORSA_DEBUG("r.size(): %i",r.size());
    if (r.size() == 0) return 0;
    double rms=0, d2ra, d2dec;
    for (unsigned int k=0; k<r.size(); ++k) {
    
        d2ra  = r[k].delta_ra.getRef();
        d2ra *= d2ra;
    
        d2dec  = r[k].delta_dec.getRef();
        d2dec *= d2dec;
    
        rms += d2ra + d2dec;
    }
    rms /= 2*r.size();
    rms = sqrt(rms);
    // ORSA_DEBUG("RMS: %f",rms);
    return rms;
}

// convert uncertainties from Equinoctial variables to regular Orbit variables
// uncertainty factor (due to non-unitary chi-squared) not included
void orsaUtil::ConvertEquinoctialUncertainties(double & sigma_a,
                                               double & sigma_e,
                                               double & sigma_i,
                                               double & sigma_omega_node,
                                               double & sigma_omega_pericenter,
                                               double & sigma_M,
                                               const orsa::EquinoctialOrbit & equinoctialOrbit,
                                               const gsl_matrix * equinoctialCovariance,
                                               const int index_p,
                                               const int index_f,
                                               const int index_g,
                                               const int index_h,
                                               const int index_k,
                                               const int index_L) {
  
    const orsa::EquinoctialOrbit & eo = equinoctialOrbit;
    const gsl_matrix * covar = equinoctialCovariance;
  
    const double e = sqrt(orsa::square(eo.f)+orsa::square(eo.g));
    const double one_minus_e2 = 1-orsa::square(e);
    //
    const double dadp = 1/one_minus_e2;
    const double dadf = 2*eo.p*eo.f/orsa::square(one_minus_e2);
    const double dadg = 2*eo.p*eo.g/orsa::square(one_minus_e2);
    //
    sigma_a = sqrt( orsa::square(dadp)*gsl_matrix_get(covar,index_p,index_p) +
                    orsa::square(dadf)*gsl_matrix_get(covar,index_f,index_f) +
                    orsa::square(dadg)*gsl_matrix_get(covar,index_g,index_g) +
                    2*dadp*dadf*gsl_matrix_get(covar,index_p,index_f) +
                    2*dadp*dadg*gsl_matrix_get(covar,index_p,index_g) +
                    2*dadf*dadg*gsl_matrix_get(covar,index_f,index_g) );
    
    const double dedf = eo.f/e;
    const double dedg = eo.g/e;
    //
    sigma_e = sqrt( orsa::square(dedf)*gsl_matrix_get(covar,index_f,index_f) +
                    orsa::square(dedg)*gsl_matrix_get(covar,index_g,index_g) +
                    2*dedf*dedg*gsl_matrix_get(covar,index_f,index_g) );
    
    const double h2_k2 = orsa::square(eo.h)+orsa::square(eo.k);
    const double one_plus_h2_k2 = 1+h2_k2;
    const double sqrt_h2_k2 = sqrt(h2_k2);
    //
    const double didh = 2*eo.h/(sqrt_h2_k2*one_plus_h2_k2);
    const double didk = 2*eo.k/(sqrt_h2_k2*one_plus_h2_k2);
    //
    sigma_i = sqrt( orsa::square(didh)*gsl_matrix_get(covar,index_h,index_h) + 
                    orsa::square(didk)*gsl_matrix_get(covar,index_k,index_k) + 
                    2*didh*didk*gsl_matrix_get(covar,index_h,index_k) );
    
    const double dOdh = -eo.k/h2_k2;
    const double dOdk =  eo.h/h2_k2;
    //
    sigma_omega_node = sqrt( orsa::square(dOdh)*gsl_matrix_get(covar,index_h,index_h) + 
                             orsa::square(dOdk)*gsl_matrix_get(covar,index_k,index_k) + 
                             2*dOdh*dOdk*gsl_matrix_get(covar,index_h,index_k) );
    
    const double f2_g2 = orsa::square(eo.f)+orsa::square(eo.g);
    //
    const double dwdf = -eo.g/f2_g2;
    const double dwdg =  eo.f/f2_g2;
    const double dwdh =  eo.k/h2_k2;
    const double dwdk = -eo.h/h2_k2;
    //
    sigma_omega_pericenter = sqrt( orsa::square(dwdf)*gsl_matrix_get(covar,index_f,index_f) + 
                                   orsa::square(dwdg)*gsl_matrix_get(covar,index_g,index_g) + 
                                   orsa::square(dwdh)*gsl_matrix_get(covar,index_h,index_h) + 
                                   orsa::square(dwdk)*gsl_matrix_get(covar,index_k,index_k) + 
                                   2*dwdf*dwdg*gsl_matrix_get(covar,index_f,index_g) + 
                                   2*dwdf*dwdh*gsl_matrix_get(covar,index_f,index_h) + 
                                   2*dwdf*dwdk*gsl_matrix_get(covar,index_f,index_k) + 
                                   2*dwdg*dwdh*gsl_matrix_get(covar,index_g,index_h) + 
                                   2*dwdg*dwdk*gsl_matrix_get(covar,index_g,index_k) + 
                                   2*dwdh*dwdk*gsl_matrix_get(covar,index_h,index_k) );
    
    const double dMdf =  eo.g/f2_g2;
    const double dMdg = -eo.f/f2_g2;
    const double dMdh =  eo.k/h2_k2;
    const double dMdk = -eo.h/h2_k2;
    const double dMdL =  1;
    //
    sigma_M = sqrt( orsa::square(dMdf)*gsl_matrix_get(covar,index_f,index_f) + 
                    orsa::square(dMdg)*gsl_matrix_get(covar,index_g,index_g) + 
                    orsa::square(dMdh)*gsl_matrix_get(covar,index_h,index_h) + 
                    orsa::square(dMdk)*gsl_matrix_get(covar,index_k,index_k) + 
                    orsa::square(dMdL)*gsl_matrix_get(covar,index_L,index_L) + 
                    2*dMdf*dMdg*gsl_matrix_get(covar,index_f,index_g) + 
                    2*dMdf*dMdh*gsl_matrix_get(covar,index_f,index_h) + 
                    2*dMdf*dMdk*gsl_matrix_get(covar,index_f,index_k) + 
                    2*dMdf*dMdL*gsl_matrix_get(covar,index_f,index_L) + 
                    2*dMdg*dMdh*gsl_matrix_get(covar,index_g,index_h) + 
                    2*dMdg*dMdk*gsl_matrix_get(covar,index_g,index_k) + 
                    2*dMdg*dMdL*gsl_matrix_get(covar,index_g,index_L) + 
                    2*dMdh*dMdk*gsl_matrix_get(covar,index_h,index_k) +
                    2*dMdh*dMdL*gsl_matrix_get(covar,index_h,index_L) +
                    2*dMdk*dMdL*gsl_matrix_get(covar,index_k,index_L) );
}
  
  


void orsaUtil::ComputeResidual(std::vector<orsaUtil::Residual>   & residual,
                               const orsaSolarSystem::OrbitWithEpoch    & orbit,
                               const orsaSolarSystem::OpticalObservationVector & obs,
                               const orsaSolarSystem::ObservatoryPositionCallback * obsPosCB, 
                               const orsa::Body * refBody,
                               orsa::BodyGroup  * bg) {
  
    orsa::Vector refBodyPosition, refBodyVelocity;
    orsa::Vector obsPosition;
    orsa::Vector relativePosition, relativeVelocity;
    orsa::Vector dr;
  
    residual.resize(obs.size());
  
    for (unsigned int k=0; k<obs.size(); ++k) {
    
        orsaSolarSystem::OpticalObservation * opticalObservation =  
            dynamic_cast<orsaSolarSystem::OpticalObservation *>(obs[k].get());
        if (!opticalObservation) {
            ORSA_DEBUG("observation is not optical");
        }
    
        if (!bg->getInterpolatedPosVel(refBodyPosition,refBodyVelocity,refBody,obs[k]->epoch.getRef())) { ORSA_DEBUG("problems"); }
    
        if (!obsPosCB->getPosition(obsPosition,obs[k].get())) { ORSA_DEBUG("problems"); }
    
        {
            orsa::Orbit localOrbit = orbit;
            localOrbit.M = fmod(orbit.M + orsa::twopi()*(obs[k]->epoch.getRef()-orbit.epoch.getRef()).get_d()/orbit.period(),orsa::twopi()); 
            if (!localOrbit.relativePosVel(relativePosition,relativeVelocity)) { ORSA_DEBUG("problems"); }
        }
    
        dr = refBodyPosition+relativePosition-obsPosition;
        const double lightTimeDelay = dr.length()/orsa::Unit::c();
        dr -= lightTimeDelay*(refBodyVelocity+relativeVelocity);
        dr = orsaSolarSystem::eclipticToEquatorial()*dr;
    
        const double  ra_orbit = fmod(atan2(dr.getY(),dr.getX())+orsa::twopi(),orsa::twopi());
        const double dec_orbit = asin(dr.getZ()/dr.length());
    
        double delta_ra = opticalObservation->ra.getRef().getRad()-ra_orbit;
        if (fabs(opticalObservation->ra.getRef().getRad()-ra_orbit+orsa::twopi()) < fabs(delta_ra)) delta_ra = opticalObservation->ra.getRef().getRad()-ra_orbit+orsa::twopi();
        if (fabs(opticalObservation->ra.getRef().getRad()-ra_orbit-orsa::twopi()) < fabs(delta_ra)) delta_ra = opticalObservation->ra.getRef().getRad()-ra_orbit-orsa::twopi();
        const double cos_dec = cos(0.5*(opticalObservation->dec.getRef().getRad()+dec_orbit));
        delta_ra *= cos_dec;
    
        double delta_dec = opticalObservation->dec.getRef().getRad()-dec_orbit;
    
        residual[k].delta_ra  = delta_ra;
        residual[k].delta_dec = delta_dec;
    }
  
}

void orsaUtil::OrbitMultifit::singleIterationDone(const gsl_multifit_fdfsolver * s) const {
  
    {
        unsigned int gslIndex=0;
        for (unsigned int k=0; k<_par->totalSize(); ++k) {
            if (!_par->isFixed(k)) {
                _par->set(k, gsl_vector_get(s->x,gslIndex));
                ++gslIndex;
            }
        }
    }
  
    double c = 1.0;
    //
    const unsigned int dof = _data->size() - _par->sizeNotFixed();
    if (dof > 0) {
        const double chi = gsl_blas_dnrm2(s->f);
        c = GSL_MAX_DBL(1.0, chi / sqrt(dof)); 
        // ORSA_DEBUG("chisq/dof = %g",  chi*chi/dof);
        // gmp_fprintf(fp,"chisq/dof = %g\n",  chi*chi/dof);
    }
    //
    const double factor = c;
  
    gsl_matrix * covar = gsl_matrix_alloc(_par->sizeNotFixed(),_par->sizeNotFixed());
  
    gsl_multifit_covar(s->J, 0.0, covar);
  
    orsa::EquinoctialOrbit equinoctialOrbit;
    //
    equinoctialOrbit.p = _par->get("equinoctialOrbit_p"); 
    equinoctialOrbit.f = _par->get("equinoctialOrbit_f"); 
    equinoctialOrbit.g = _par->get("equinoctialOrbit_g"); 
    equinoctialOrbit.h = _par->get("equinoctialOrbit_h"); 
    equinoctialOrbit.k = _par->get("equinoctialOrbit_k"); 
    equinoctialOrbit.L = _par->get("equinoctialOrbit_L"); 
  
    /* orsaSolarSystem::OrbitWithEpoch orbit;
       orbit.epoch = orbitEpoch.getRef();
       //
       orbit.a = _par->get("orbit_a");
       orbit.e = _par->get("orbit_e");
       orbit.i = _par->get("orbit_i");
       orbit.omega_node       = _par->get("orbit_omega_node");
       orbit.omega_pericenter = _par->get("orbit_omega_pericenter");
       orbit.M                = _par->get("orbit_M");
    */
    //
    // equinoctialOrbit.get(orbit);
  
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

    ORSA_DEBUG("p: %14.10f +/- %14.10f [AU]",
               orsa::FromUnits(equinoctialOrbit.p,orsa::Unit::AU,-1),
               orsa::FromUnits(factor*ERR(0),orsa::Unit::AU,-1));

    ORSA_DEBUG("f: %14.10f +/- %14.10f",
               equinoctialOrbit.f,
               factor*ERR(1));
  
    ORSA_DEBUG("g: %14.10f +/- %14.10f",
               equinoctialOrbit.g,
               factor*ERR(2));
  
    ORSA_DEBUG("h: %14.10f +/- %14.10f",
               equinoctialOrbit.h,
               factor*ERR(3));
  
    ORSA_DEBUG("k: %14.10f +/- %14.10f",
               equinoctialOrbit.k,
               factor*ERR(4));
  
    ORSA_DEBUG("L: %14.10f +/- %14.10f [deg]",
               orsa::radToDeg()*equinoctialOrbit.L,
               orsa::radToDeg()*factor*ERR(5));
  
    {
        orsa::Orbit orbit;
        equinoctialOrbit.get(orbit);
        // orsa::print(orbit);
    
        double sigma_a, sigma_e, sigma_i, sigma_omega_node, sigma_omega_pericenter, sigma_M;
    
        orsaUtil::ConvertEquinoctialUncertainties(sigma_a,
                                                  sigma_e,
                                                  sigma_i,
                                                  sigma_omega_node,
                                                  sigma_omega_pericenter,
                                                  sigma_M,
                                                  equinoctialOrbit,
                                                  covar,
                                                  0,1,2,3,4,5);
    
        ORSA_DEBUG("a: %14.10f +/- %14.10f [AU]",
                   orsa::FromUnits(orbit.a,orsa::Unit::AU,-1),
                   orsa::FromUnits(factor*sigma_a,orsa::Unit::AU,-1));
    
        ORSA_DEBUG("e: %14.10f +/- %14.10f",
                   orbit.e,
                   factor*sigma_e);
    
        ORSA_DEBUG("i: %14.10f +/- %14.10f [deg]",
                   orsa::radToDeg()*orbit.i,
                   orsa::radToDeg()*factor*sigma_i);
    
        ORSA_DEBUG("O: %14.10f +/- %14.10f [deg]",
                   orsa::radToDeg()*orbit.omega_node,
                   orsa::radToDeg()*factor*sigma_omega_node);
    
        ORSA_DEBUG("w: %14.10f +/- %14.10f [deg]",
                   orsa::radToDeg()*orbit.omega_pericenter,
                   orsa::radToDeg()*factor*sigma_omega_pericenter);
    
        ORSA_DEBUG("M: %14.10f +/- %14.10f [deg]",
                   orsa::radToDeg()*orbit.M,
                   orsa::radToDeg()*factor*sigma_M);
        
        {
            // approx moid
            double moid, M1, M2;
            orsa::Orbit eo; // earth orbit
            eo.mu = orsaSolarSystem::Data::GMSun();
            eo.a  = orsa::FromUnits(1,orsa::Unit::AU);
            eo.e  = 0.0;
            eo.i  = 0.0;
            eo.omega_node = 0.0;
            eo.omega_pericenter = 0.0;
            eo.M  = 0.0;
            if (orsa::MOID(moid,M1,M2,orbit,eo,134923)) {
                ORSA_DEBUG("approx. MOID: %.3f [AU]",orsa::FromUnits(moid,orsa::Unit::AU,-1));
            }
        }
        
    }
    
    gsl_matrix_free(covar);  
}
