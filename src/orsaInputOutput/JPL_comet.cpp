#include <orsaInputOutput/JPL_comet.h>

#include <orsa/util.h>

#include <orsaSolarSystem/datetime.h>

using namespace orsa;
using namespace orsaInputOutput;
using namespace orsaSolarSystem;

bool JPLCometFile::goodLine(const char * line) {
  
  if (strlen(line) < 113) return false;
  
  if (!isdigit(line[45])) return false;
  if (!isdigit(line[52])) return false;
  
  return true;
}

bool JPLCometFile::processLine(const char * line) {
  
  std::string s_designation;
  std::string s_epoch_MJD;
  std::string s_q, s_e, s_i, s_peri, s_node;
  std::string s_Tp_y, s_Tp_m, s_Tp_d;
  
  s_designation.assign(line,0,38);
  removeLeadingAndTrailingSpaces(s_designation);
  
  // ORSA_DEBUG("des: %s",s_designation.c_str());
  
  s_epoch_MJD.assign(line,39,7);
  
  s_q.assign(line,47,11);
  s_e.assign(line,59,10);
  s_i.assign(line,70,9);
  s_peri.assign(line,80,9);
  s_node.assign(line,90,9);
  
  s_Tp_y.assign(line,100,4);
  s_Tp_m.assign(line,104,2);
  s_Tp_d.assign(line,106,8);
  
  const Time Tp = orsaSolarSystem::gregorTime(atoi(s_Tp_y.c_str()),
					      atoi(s_Tp_m.c_str()),
					      Double(s_Tp_d));
  
  // timescale?
  const Time epoch = orsaSolarSystem::julianToTime(orsaSolarSystem::MJD2JD(Double(s_epoch_MJD)));
  
  Orbit orbit;
  orbit.mu = orsa::Unit::instance()->getG()*orsa::FromUnits(orsa::one(),orsa::Unit::MSUN); 
  orbit.e  = Double(s_e);
  //
  if (orbit.e > 0.99) {
    // non-periodic orbit, not included for the moment
    return false;
  }
  //
  orbit.a                = FromUnits(Double(s_q)/(one()-orbit.e),orsa::Unit::AU);
  orbit.i                = degToRad() * Double(s_i);
  orbit.omega_node       = degToRad() * Double(s_node);
  orbit.omega_pericenter = degToRad() * Double(s_peri);
  orbit.M                = twopi()*(epoch-Tp).asDouble()/orbit.period();
  
  JPLCometDataElement element;
  //
  element.orbit       = orbit;
  element.epoch       = epoch;
  element.designation = s_designation;
  
  _data.push_back(element);
  
  return true;
}
