#include <orsaInputOutput/MPC_asteroid.h>

#include <orsaInputOutput/MPC.h>

#include <orsa/util.h>

#include <orsaSolarSystem/data.h>
#include <orsaSolarSystem/datetime.h>

using namespace orsa;
using namespace orsaInputOutput;
using namespace orsaSolarSystem;

// static members of MPCAsteroidFile
const double MPCAsteroidFile::NEO_max_q = FromUnits(1.3, orsa::Unit::AU);

bool MPCAsteroidFile::goodLine(const char * line) {
    if (strlen(line) < 202) return false;
    return true;
}

bool MPCAsteroidFile::processLine(const char * line) {
    
    std::string s_designation;
    std::string s_H;
    std::string s_G;
    std::string s_epoch;
    std::string s_a, s_e, s_i, s_node, s_peri, s_M;
    
    // first, fields affected by select_* code
    
    s_designation.assign(line,0,7);
    removeLeadingAndTrailingSpaces(s_designation);
    
    if (select_number.isSet()) {
        if (select_number.getRef() != MPC_packedNumber(s_designation)) {
            return false;
        }
    }
    if (select_designation.isSet()) {
        if (select_designation.getRef() != s_designation.c_str()) {
            return false;
        }
    }
    
    orsaSolarSystem::OrbitWithEpoch orbit;
    s_e.assign(line,70,9);
    orbit.e = atof(s_e.c_str());
    s_a.assign(line,92,11);
    orbit.a = FromUnits(atof(s_a.c_str()),orsa::Unit::AU);
    
    if (select_NEO.isSet()) {
        if (select_NEO.getRef()) {
            if (orbit.a*(1-orbit.e) > NEO_max_q) {
                return false;
            }
        }
    }
    
    // now, all the remaining fields
    
    s_H.assign(line,8,5);
    
    s_G.assign(line,14,5);
    
    s_epoch.assign(line,20,5);
  
    s_M.assign(line,26,9);
    s_peri.assign(line,37,9);
    s_node.assign(line,48,9);
    s_i.assign(line,59,9);
    // s_e.assign(line,70,9);
    // s_a.assign(line,92,11);
    
    // orsaSolarSystem::OrbitWithEpoch orbit;
    orbit.epoch = MPC_packedToTime(s_epoch.c_str());
    orbit.mu = orsaSolarSystem::Data::GMSun(); 
    // orbit.e  = atof(s_e.c_str());
    //
    if (orbit.e > 0.99) {
        // non-periodic orbit, not included for the moment
        return false;
    }
    //
    // orbit.a                = FromUnits(atof(s_a.c_str()),orsa::Unit::AU);
    orbit.i                = degToRad() * atof(s_i.c_str());
    orbit.omega_node       = degToRad() * atof(s_node.c_str());
    orbit.omega_pericenter = degToRad() * atof(s_peri.c_str());
    orbit.M                = degToRad() * atof(s_M.c_str());
  
    MPCAsteroidDataElement element;
    //
    element.orbit       = orbit;
    element.H           = atof(s_H.c_str());
    element.G           = atof(s_G.c_str());
    if (MPC_packedNumber(s_designation) != 0) {
        element.number = MPC_packedNumber(s_designation);
    }
    if (strlen(s_designation.c_str()) != 0) {
        element.designation = s_designation;
    }
  
    _data.push_back(element);
  
    return true;
}
