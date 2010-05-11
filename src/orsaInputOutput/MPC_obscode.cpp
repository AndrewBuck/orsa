#include <orsaInputOutput/MPC_obscode.h>

#include <orsa/unit.h>

#include <orsa/util.h>

#include <orsaSolarSystem/data.h>

using namespace orsa;
using namespace orsaInputOutput;
using namespace orsaSolarSystem;

bool MPCObsCodeFile::goodLine(const char * line) {
  
    if (strlen(line) < 30) return false;
  
    if (!isalnum(line[0])) return false;
    if (!isalnum(line[1])) return false;
    if (!isalnum(line[2])) return false;
    if (!isspace(line[3])) return false;
  
    return true;
}

bool MPCObsCodeFile::processLine(const char * line) {
  
    std::string s_obscode, s_longitude, s_pos_xy, s_pos_z, s_name;
  
    s_obscode.assign(line,0,3);
    s_longitude.assign(line,4,9);
    s_pos_xy.assign(line,13,8);
    s_pos_z.assign(line,21,9);
    // name.assign(line,30,strlen(line)-30-1); // the last -1 is set to avoid the '\n' character in the name
    s_name.assign(line,30,strlen(line)-30); // the \n has already been removed  
  
    removeLeadingAndTrailingSpaces(s_longitude);
    removeLeadingAndTrailingSpaces(s_pos_xy);
    removeLeadingAndTrailingSpaces(s_pos_z);
    removeLeadingAndTrailingSpaces(s_name);
  
    removeLeadingPlusSign(s_longitude);
    removeLeadingPlusSign(s_pos_xy);
    removeLeadingPlusSign(s_pos_z);
    removeLeadingPlusSign(s_name);
  
    // another test
    {
        const unsigned int size = strlen(s_longitude.c_str());
        for (unsigned int k=0; k<size; ++k) {
            if (isalpha(s_longitude[k])) {
                return false;
            }
        }
    }
  
    Observatory workObs;
  
    if (s_longitude.size()) workObs.lon = degToRad() * atof(s_longitude.c_str());
    if (s_pos_xy.size())    workObs.pxy = atof(s_pos_xy.c_str()) * orsaSolarSystem::Data::REarth();
    if (s_pos_z.size())     workObs.pz  = atof(s_pos_z.c_str())  * orsaSolarSystem::Data::REarth();
  
    workObs.name = s_name;
  
    _data.observatory[s_obscode] = workObs;
  
    _data.code.push_back(s_obscode);
  
    return true;
}
