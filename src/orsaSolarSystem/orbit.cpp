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
  return sqrt(rms);
}

void orsaSolarSystem::ComputeResidual(std::vector<orsaSolarSystem::Residual> & residual,
				      const orsaSolarSystem::OrbitWithEpoch & orbit,
				      const orsaSolarSystem::ObservationVector & obs,
				      const orsaSolarSystem::ObservatoryPositionCallback * obsPosCB, 
				      const orsa::Body * refBody,
				      orsa::BodyGroup * bg) {
  
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
    const double lightTimeDelay = dr.length()/orsa::Unit::instance()->getC();
    dr -= lightTimeDelay*(refBodyVelocity+relativeVelocity);
    dr = orsaSolarSystem::eclipticToEquatorial()*dr;
    
    const double  ra_orbit = fmod(atan2(dr.getY(),dr.getX())+orsa::twopi(),orsa::twopi());
    const double dec_orbit = asin(dr.getZ()/dr.length());
    
    double delta_ra = obs[k]->ra.getRef().getRad()-ra_orbit;
    if (fabs(obs[k]->ra.getRef().getRad()-ra_orbit+orsa::twopi()) < fabs(delta_ra)) delta_ra = obs[k]->ra.getRef().getRad()-ra_orbit+orsa::twopi();
    if (fabs(obs[k]->ra.getRef().getRad()-ra_orbit-orsa::twopi()) < fabs(delta_ra)) delta_ra = obs[k]->ra.getRef().getRad()-ra_orbit-orsa::twopi();
    const double cos_dec = orsa::cos(0.5*(obs[k]->dec.getRef().getRad()+dec_orbit));
    delta_ra *= cos_dec;
    
    double delta_dec = obs[k]->dec.getRef().getRad()-dec_orbit;
    
    residual[k].delta_ra  = delta_ra;
    residual[k].delta_dec = delta_dec;
  }
  
}

