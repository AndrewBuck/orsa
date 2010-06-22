#include <orsaInputOutput/RWO.h>

// #include <orsa/print.h>
#include <orsa/util.h>

#include <orsaSolarSystem/datetime.h>

using namespace orsa;
using namespace orsaSolarSystem;
using namespace orsaInputOutput;

bool RWOFile::goodLine(const char * line) {
  
    if (strlen(line) < 179) return false;
  
    if (!isdigit(line[17])) return false;
    if (!isspace(line[21])) return false;
  
    // ORSA_DEBUG("%s   [GOOD LINE]",line);
  
    return true;
}

bool RWOFile::processLine(const char * line) {
  
    // std::string number;
    std::string s_designation;
    std::string s_epoch, s_ra, s_decSign, s_dec;
    std::string s_mag;
    char        c_band;
    std::string s_obsCode;
  
    s_designation.assign(line,0,10); 
    removeLeadingAndTrailingSpaces(s_designation);
  
    // no discovery information in RWO files
  
    s_epoch.assign(line,17,21);
 
    s_ra.assign(line,50,12);

    s_decSign.assign(line,103,1);
    s_dec.assign(line,104,11);
  
    s_mag.assign(line,156,4);
    removeLeadingAndTrailingSpaces(s_mag);
  
    // s_magCode.assign(line,161,1);
    c_band = line[161];
    
    s_obsCode.assign(line,176,3);
    removeLeadingAndTrailingSpaces(s_obsCode);
  
    osg::ref_ptr<OpticalObservation> workObs = new OpticalObservation;
  
    workObs->designation = s_designation;
    workObs->obsCode     = s_obsCode;
  
    if (strlen(s_mag.c_str()) > 0) {
        workObs->mag  = atof(s_mag.c_str());
        workObs->band = c_band;
    }
  
    {
        int y, m; 
        double d;
        gmp_sscanf(s_epoch.c_str(),"%d %d %lf",&y,&m,&d);
        // ORSA_DEBUG("d: %f",d());
        // ORSA_DEBUG("remember: UTC!!");
        workObs->epoch = 
            orsaSolarSystem::FromTimeScale(orsaSolarSystem::gregorTime(y,m,d),
                                           orsaSolarSystem::TS_UTC);
        // ORSA_DEBUG("%s",s_epoch.c_str());
        // ORSA_DEBUG("y: %i   m: %i   d: %f",y,m,d());
        // orsa::print(workObs->epoch.getRef());
    }
  
    {
        int h, m;
        double s;
        gmp_sscanf(s_ra.c_str(),"%d %d %lf",&h,&m,&s);
        Angle tmp; tmp.setHMS(h,m,s);
        workObs->ra = tmp;
        // ORSA_DEBUG("h: %i   m: %i   s: %f",h,m,s());
    }
  
    {
        const int sign = s_decSign == "-" ? -1 : +1;
        int d, p;
        double s;
        gmp_sscanf(s_dec.c_str(),"%d %d %lf",&d,&p,&s);
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
