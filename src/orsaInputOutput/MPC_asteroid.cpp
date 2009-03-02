#include <orsaInputOutput/MPC_asteroid.h>

#include <orsaInputOutput/MPC.h>

#include <orsa/util.h>

#include <orsaSolarSystem/datetime.h>

using namespace orsa;
using namespace orsaInputOutput;
using namespace orsaSolarSystem;

bool MPCAsteroidFile::goodLine(const char * line) {
  if (strlen(line) < 202) return false;
  return true;
}

bool MPCAsteroidFile::processLine(const char * line) {
  
  std::string s_designation;
  std::string s_H;
  std::string s_epoch;
  std::string s_a, s_e, s_i, s_node, s_peri, s_M;
  
  s_designation.assign(line,0,7);
  removeLeadingAndTrailingSpaces(s_designation);
  
  s_H.assign(line,8,5);
  
  s_epoch.assign(line,20,5);
  
  s_M.assign(line,26,9);
  s_peri.assign(line,37,9);
  s_node.assign(line,48,9);
  s_i.assign(line,59,9);
  s_e.assign(line,70,9);
  s_a.assign(line,92,11);
  
  const Time epoch = MPC_packedToTime(s_epoch);
  
  Orbit orbit;
  orbit.mu = orsa::Unit::instance()->getG()*orsa::FromUnits(orsa::one(),orsa::Unit::MSUN); 
  orbit.e  = Double(s_e);
  //
  if (orbit.e > 0.99) {
    // non-periodic orbit, not included for the moment
    return false;
  }
  //
  orbit.a                = FromUnits(Double(s_a),orsa::Unit::AU);
  orbit.i                = degToRad() * Double(s_i);
  orbit.omega_node       = degToRad() * Double(s_node);
  orbit.omega_pericenter = degToRad() * Double(s_peri);
  orbit.M                = degToRad() * Double(s_M);
  
  MPCAsteroidDataElement element;
  //
  element.orbit       = orbit;
  element.epoch       = epoch;
  element.H           = Double(s_H);
  element.designation = s_designation;
  
  _data.push_back(element);
  
  return true;
}
