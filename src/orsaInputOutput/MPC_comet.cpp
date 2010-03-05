#include <orsaInputOutput/MPC_comet.h>

#include <orsa/util.h>

#include <orsaSolarSystem/data.h>
#include <orsaSolarSystem/datetime.h>

using namespace orsa;
using namespace orsaInputOutput;
using namespace orsaSolarSystem;

bool MPCCometFile::goodLine(const char * line) {
  if (strlen(line) < 103) return false;
  return true;
}

bool MPCCometFile::processLine(const char * line) {
  
  std::string s_Tp_y, s_Tp_m, s_Tp_d;
  std::string s_q, s_e, s_peri, s_node, s_i;
  std::string s_epoch_y, s_epoch_m, s_epoch_d;
  std::string s_designation;
  
  s_Tp_y.assign(line,14,4);
  s_Tp_m.assign(line,19,2);
  s_Tp_d.assign(line,22,7);
  
  s_q.assign(line,30,9);
  s_e.assign(line,40,9);
  s_peri.assign(line,51,8);
  s_node.assign(line,61,8);
  s_i.assign(line,71,8);
  
  s_epoch_y.assign(line,81,4);
  s_epoch_m.assign(line,85,2);
  s_epoch_d.assign(line,87,2);
  
  s_designation.assign(line,102,strlen(line)-102); // the \n has already been removed  
  
  removeLeadingAndTrailingSpaces(s_epoch_y);
  
  if (strlen(s_epoch_y.c_str()) != 4) {
    // orbit with undefined epoch
    return false;
  }
  
  const Time Tp = orsaSolarSystem::gregorTime(atoi(s_Tp_y.c_str()),
					      atoi(s_Tp_m.c_str()),
					      atof(s_Tp_d.c_str()));
  
  const Time epoch = orsaSolarSystem::gregorTime(atoi(s_epoch_y.c_str()),
						 atoi(s_epoch_m.c_str()),
						 atof(s_epoch_d.c_str()));
  
  OrbitWithEpoch orbit;
  orbit.epoch = epoch;
  orbit.mu = orsaSolarSystem::Data::GMSun(); 
  orbit.e  = atof(s_e.c_str());
  //
  if (orbit.e > 0.99) {
    // non-periodic orbit, not included for the moment
    return false;
  }
  //
  orbit.a                = FromUnits(atof(s_q.c_str())/(1-orbit.e),orsa::Unit::AU);
  orbit.i                = degToRad() * atof(s_i.c_str());
  orbit.omega_node       = degToRad() * atof(s_node.c_str());
  orbit.omega_pericenter = degToRad() * atof(s_peri.c_str());
  orbit.M                = twopi()*(epoch-Tp).get_d()/orbit.period();
  
  MPCCometDataElement element;
  //
  element.orbit       = orbit;
  element.designation = s_designation;
  
  _data.push_back(element);
  
  return true;
}
