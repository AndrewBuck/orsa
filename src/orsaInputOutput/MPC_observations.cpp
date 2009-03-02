#include <orsaInputOutput/MPC_observations.h>

// #include <orsa/print.h>
#include <orsa/util.h>

#include <orsaSolarSystem/datetime.h>

using namespace orsa;
using namespace orsaSolarSystem;
using namespace orsaInputOutput;

bool MPCObservationsFile::goodLine(const char * line) {
  
  if (strlen(line) < 80) return false;
  
  // the following line positions MUST be spaces: 19, 22, 34, 37, 47, 50;
  // the following CAN be non-spaces: 31, 43, 55;
  //
  if (!isspace(line[19])) return false;
  if (!isspace(line[22])) return false;
  if (!isspace(line[34])) return false;
  if (!isspace(line[37])) return false;
  if (!isspace(line[47])) return false;
  if (!isspace(line[50])) return false;
  
  return true;
}

bool MPCObservationsFile::processLine(const char * line) {
  
  // std::string number;
  std::string s_designation, s_discovery, s_note1, s_note2;
  std::string s_epoch, s_ra, s_decSign, s_dec;
  std::string s_mag, s_magCode;
  std::string s_obsCode;
  
  s_designation.assign(line,5,7); 
  // s_designation.assign(line,0,12); 
  removeLeadingAndTrailingSpaces(s_designation);
  
  s_discovery.assign(line,12,1); 
  removeAllSpaces(s_discovery);
  
  s_note1.assign(line,13,1);
  s_note2.assign(line,14,1);
  
  s_epoch.assign(line,15,17);
 
  s_ra.assign(line,32,12);

  s_decSign.assign(line,44,1);
  s_dec.assign(line,45,11);
  
  s_mag.assign(line,65,5);
  removeLeadingAndTrailingSpaces(s_mag);
  
  s_magCode.assign(line,70,1);
  
  s_obsCode.assign(line,77,3);
  removeLeadingAndTrailingSpaces(s_obsCode);
  
  osg::ref_ptr<Observation> workObs = new Observation;
  
  workObs->designation = s_designation;
  workObs->obsCode     = s_obsCode;
  workObs->discovery   = (strlen(s_discovery.c_str()) > 0);
  
  if (strlen(s_mag.c_str()) > 0) {
    workObs->mag     = Double(s_mag);
    workObs->magCode = s_magCode;
  }
  
  {
    int y, m; 
    Double d;
    gmp_sscanf(s_epoch.c_str(),"%d %d %Ff",&y,&m,&d);
    // ORSA_DEBUG("d: %Ff",d.get_mpf_t());
    // ORSA_DEBUG("remember: UTC!!");
    workObs->epoch = 
      orsaSolarSystem::FromTimeScale(orsaSolarSystem::gregorTime(y,m,d),
				     orsaSolarSystem::TS_UTC);
    // ORSA_DEBUG("%s",s_epoch.c_str());
    // ORSA_DEBUG("y: %i   m: %i   d: %Ff",y,m,d.get_mpf_t());
    // orsa::print(workObs->epoch.getRef());
  }
  
  {
    int h, m;
    Double s;
    gmp_sscanf(s_ra.c_str(),"%d %d %Ff",&h,&m,&s);
    Angle tmp; tmp.setHMS(h,m,s);
    workObs->ra = tmp;
    // ORSA_DEBUG("h: %i   m: %i   s: %Ff",h,m,s.get_mpf_t());
  }
  
  {
    const int sign = s_decSign == "-" ? -1 : +1;
    int d, p;
    Double s;
    gmp_sscanf(s_dec.c_str(),"%d %d %Ff",&d,&p,&s);
    Angle tmp; tmp.setDPS(d,p,s,sign);
    workObs->dec = tmp;
  }
  
  if ((s_designation != "") && 
      (s_obsCode != "") &&
      (strlen(s_obsCode.c_str())) == 3) {
    if ( (isalnum(s_obsCode[0])) &&
	 (isalnum(s_obsCode[1])) &&
	 (isalnum(s_obsCode[2]))) {
      _data.push_back(workObs);
    } else {
      return false;
    }
  } else {
    return false;
  }
  
  return true;
}
