#include <orsaSolarSystem/orbit.h>

#include <orsaSolarSystem/obleq.h>

orsa::Double orsaSolarSystem::RMS(const std::vector<orsaSolarSystem::Residual> & r) {
  ORSA_DEBUG("r.size(): %i",r.size());
  if (r.size() == 0) return 0;
  orsa::Double rms=0, d2ra, d2dec;
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
      localOrbit.M = orsa::fmod(orbit.M + orsa::twopi()*(obs[k]->epoch.getRef()-orbit.epoch.getRef()).asDouble()/orbit.period(),orsa::twopi()); 
      if (!localOrbit.relativePosVel(relativePosition,relativeVelocity)) { ORSA_DEBUG("problems"); }
    }
    
    dr = refBodyPosition+relativePosition-obsPosition;
    const orsa::Double lightTimeDelay = dr.length()/orsa::Unit::instance()->getC();
    dr -= lightTimeDelay*(refBodyVelocity+relativeVelocity);
    dr = orsaSolarSystem::eclipticToEquatorial()*dr;
    
    const orsa::Double  ra_orbit = orsa::fmod(orsa::atan2(dr.getY(),dr.getX())+orsa::twopi(),orsa::twopi());
    const orsa::Double dec_orbit = orsa::halfpi() - orsa::acos(dr.getZ()/dr.length());
    
    orsa::Double delta_ra = obs[k]->ra.getRef().getRad()-ra_orbit;
    if (orsa::fabs(obs[k]->ra.getRef().getRad()-ra_orbit+orsa::twopi()) < orsa::fabs(delta_ra)) delta_ra = obs[k]->ra.getRef().getRad()-ra_orbit+orsa::twopi();
    if (orsa::fabs(obs[k]->ra.getRef().getRad()-ra_orbit-orsa::twopi()) < orsa::fabs(delta_ra)) delta_ra = obs[k]->ra.getRef().getRad()-ra_orbit-orsa::twopi();
    const orsa::Double cos_dec = orsa::cos(0.5*(obs[k]->dec.getRef().getRad()+dec_orbit));
    delta_ra *= cos_dec;
    
    orsa::Double delta_dec = obs[k]->dec.getRef().getRad()-dec_orbit;
    
    residual[k].delta_ra  = delta_ra;
    residual[k].delta_dec = delta_dec;
  }
  
}

