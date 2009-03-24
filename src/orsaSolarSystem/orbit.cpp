#include <orsaSolarSystem/orbit.h>

#include <orsaSolarSystem/obleq.h>

double orsaSolarSystem::RMS(const std::vector<orsaSolarSystem::Residual> & r) {
  ORSA_DEBUG("r.size(): %i",r.size());
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

void orsaSolarSystem::ComputeResidual(std::vector<orsaSolarSystem::Residual>   & residual,
				      const orsaSolarSystem::OrbitWithEpoch    & orbit,
				      const orsaSolarSystem::ObservationVector & obs,
				      const orsaSolarSystem::ObservatoryPositionCallback * obsPosCB, 
				      const orsa::Body * refBody,
				      orsa::BodyGroup  * bg) {
  
  orsa::Vector refBodyPosition, refBodyVelocity;
  orsa::Vector obsPosition;
  orsa::Vector relativePosition, relativeVelocity;
  orsa::Vector dr;
  
  residual.resize(obs.size());
  
  for (unsigned int k=0; k<obs.size(); ++k) {
    
    if (!bg->getInterpolatedPosVel(refBodyPosition,refBodyVelocity,refBody,obs[k]->epoch.getRef())) { ORSA_DEBUG("problems"); }
    
    if (!obsPosCB->getPosition(obsPosition,obs[k]->obsCode.getRef(),obs[k]->epoch.getRef())) { ORSA_DEBUG("problems"); }
    
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
    
    double delta_ra = obs[k]->ra.getRef().getRad()-ra_orbit;
    if (fabs(obs[k]->ra.getRef().getRad()-ra_orbit+orsa::twopi()) < fabs(delta_ra)) delta_ra = obs[k]->ra.getRef().getRad()-ra_orbit+orsa::twopi();
    if (fabs(obs[k]->ra.getRef().getRad()-ra_orbit-orsa::twopi()) < fabs(delta_ra)) delta_ra = obs[k]->ra.getRef().getRad()-ra_orbit-orsa::twopi();
    const double cos_dec = cos(0.5*(obs[k]->dec.getRef().getRad()+dec_orbit));
    delta_ra *= cos_dec;
    
    double delta_dec = obs[k]->dec.getRef().getRad()-dec_orbit;
    
    residual[k].delta_ra  = delta_ra;
    residual[k].delta_dec = delta_dec;
  }
  
}

void orsaSolarSystem::OrbitMultifit::singleIterationDone(const gsl_multifit_fdfsolver * s) const {
  
  for (unsigned int k=0; k<_par->size(); ++k) {
    _par->set(k, gsl_vector_get(s->x,k));
  }
  
  double c = 1.0;
  //
  const unsigned int dof = _data->size() - _par->size();
  if (dof > 0) {
    const double chi = gsl_blas_dnrm2(s->f);
    c = GSL_MAX_DBL(1.0, chi / sqrt(dof)); 
    // ORSA_DEBUG("chisq/dof = %g",  chi*chi/dof);
    // gmp_fprintf(fp,"chisq/dof = %g\n",  chi*chi/dof);
  }
  //
  const double factor = c;
  
  gsl_matrix * covar = gsl_matrix_alloc(_par->size(),_par->size());
  
  gsl_multifit_covar(s->J, 0.0, covar);
  
  orsaSolarSystem::OrbitWithEpoch orbit;
  //
  orbit.epoch = orbitEpoch.getRef();
  //
  orbit.a = _par->get("orbit_a");
  orbit.e = _par->get("orbit_e");
  orbit.i = _par->get("orbit_i");
  orbit.omega_node       = _par->get("orbit_omega_node");
  orbit.omega_pericenter = _par->get("orbit_omega_pericenter");
  orbit.M                = _par->get("orbit_M");
  //
#warning remember to set orbit.mu
  
#define ERR(i) sqrt(gsl_matrix_get(covar,i,i))

  ORSA_DEBUG("a: %14.10f +/- %14.10f [AU]",
	     orsa::FromUnits(orbit.a,orsa::Unit::AU,-1),
	     orsa::FromUnits(factor*ERR(0),orsa::Unit::AU,-1));

  ORSA_DEBUG("e: %14.10f +/- %14.10f",
	     orbit.e,
	     factor*ERR(1));
  
  ORSA_DEBUG("i: %14.10f +/- %14.10f [deg]",
	     orsa::radToDeg()*orbit.i,
	     orsa::radToDeg()*factor*ERR(2));
  
  ORSA_DEBUG("O: %14.10f +/- %14.10f [deg]",
	     orsa::radToDeg()*orbit.omega_node,
	     orsa::radToDeg()*factor*ERR(3));
  
  ORSA_DEBUG("w: %14.10f +/- %14.10f [deg]",
	     orsa::radToDeg()*orbit.omega_pericenter,
	     orsa::radToDeg()*factor*ERR(4));
  
  ORSA_DEBUG("M: %14.10f +/- %14.10f [deg]",
	     orsa::radToDeg()*orbit.M,
	     orsa::radToDeg()*factor*ERR(5));
  
  gsl_matrix_free(covar);  
}
