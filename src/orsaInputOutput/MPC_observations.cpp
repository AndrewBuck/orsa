#include <orsaInputOutput/MPC_observations.h>
#include <orsaInputOutput/MPC.h>

#include <orsa/util.h>

#include <orsaSolarSystem/datetime.h>

using namespace orsa;
using namespace orsaSolarSystem;
using namespace orsaInputOutput;

bool MPCObservationsFile::goodLine(const char * line) {
  
    if (line==0) return false;
  
    if (strlen(line) < 80) return false;
  
    if ( isspace(line[15])) return false;
    if ( isspace(line[16])) return false;
    if ( isspace(line[17])) return false;
    if ( isspace(line[18])) return false;
    if (!isspace(line[19])) return false;
    if (!isspace(line[22])) return false;
    if ( isspace(line[77])) return false;
    if ( isspace(line[78])) return false;
    if ( isspace(line[79])) return false;  
  
    if (!isdigit(line[17])) return false;
  
    // ORSA_DEBUG("good line: [%s]",line);
  
    return true;
}

bool MPCObservationsFile::processLine(const char * line) {
  
    // std::string number;
    std::string s_number, s_designation, s_discovery, s_note1, s_note2;
    std::string s_epoch, s_ra, s_decSign, s_dec;
    std::string s_mag;
    char        c_band;
    std::string s_obsCode;
  
    // first, fields affected by select_* code
  
    s_obsCode.assign(line,77,3);
    removeLeadingAndTrailingSpaces(s_obsCode);
    if (select_obsCode.isSet()) {
        if (s_obsCode != select_obsCode.getRef()) {
            // ORSA_DEBUG("--OUT--   select: [%s] read: [%s]",select_obsCode.getRef().c_str(),s_obsCode.c_str());
            return false;
        }
    }
  
    s_epoch.assign(line,15,17);
    orsa::Time epoch;
    {
        int y, m; 
        double d;
        gmp_sscanf(s_epoch.c_str(),"%d %d %lf",&y,&m,&d);
        epoch = orsaSolarSystem::FromTimeScale(orsaSolarSystem::gregorTime(y,m,d),
                                               orsaSolarSystem::TS_UTC);
        // orsa::print(epoch);
        if (select_startEpoch.isSet()) {
            if (epoch < select_startEpoch.getRef()) {
                // ORSA_DEBUG("--OUT--");
                return false;
            }
        }
        if (select_stopEpoch.isSet()) {
            if (epoch > select_stopEpoch.getRef()) {
                // ORSA_DEBUG("--OUT--");
                return false;
            }
        }
    }
  
    // now, all the remaining fields
  
    // it's a number in MPC format, can contain letters...
    s_number.assign(line,0,5);
    removeLeadingAndTrailingSpaces(s_number);
  
    s_designation.assign(line,5,7); 
    removeLeadingAndTrailingSpaces(s_designation);
  
    s_discovery.assign(line,12,1); 
    removeAllSpaces(s_discovery);
  
    s_note1.assign(line,13,1);
    s_note2.assign(line,14,1);
  
    if (!twoLinesCall) {
        if (s_note2[0] == 'R') return false; // RadarObservation,     requires two lines
        if (s_note2[0] == 'V') return false; // RovingObservation,    requires two lines
        if (s_note2[0] == 'S') return false; // SatelliteObservation, requires two lines
    } else {
        if ( (s_note2[0] != 'R') && 
             (s_note2[0] != 'V') &&
             (s_note2[0] != 'S') ) {
            return false;
        }    
    }
  
    // always check this
    if (s_note2[0] == 'r') return false; // is the second line, ignore
    if (s_note2[0] == 'v') return false; // is the second line, ignore
    if (s_note2[0] == 's') return false; // is the second line, ignore
  
    // moved earlier, for faster select_* code
    // s_epoch.assign(line,15,17);
  
    s_ra.assign(line,32,12);

    s_decSign.assign(line,44,1);
    s_dec.assign(line,45,11);
  
    s_mag.assign(line,65,5);
    removeLeadingAndTrailingSpaces(s_mag);
    
    // s_magCode.assign(line,70,1);
    c_band = line[70];
    
    // moved earlier, for faster select_* code
    // s_obsCode.assign(line,77,3);
    // removeLeadingAndTrailingSpaces(s_obsCode);
  
    // osg::ref_ptr<OpticalObservation> workObs = new OpticalObservation;
    //
    osg::ref_ptr<orsaSolarSystem::Observation> workObs;
    if (twoLinesCall) {
        switch (s_note2[0]) {
            case 'R': workObs = new orsaSolarSystem::RadarObservation;     break;
            case 'V': workObs = new orsaSolarSystem::RovingObservation;    break;
            case 'S': workObs = new orsaSolarSystem::SatelliteObservation; break;
            default:  workObs = new orsaSolarSystem::OpticalObservation;   break;
        }
    } else {
        workObs = new orsaSolarSystem::OpticalObservation;
    }
  
    if (MPC_packedNumber(s_number) != 0) {
        workObs->number = MPC_packedNumber(s_number);
    }
    if (strlen(s_designation.c_str()) != 0) {
        workObs->designation = s_designation;
    }
    workObs->obsCode     = s_obsCode;
    workObs->discovery   = (strlen(s_discovery.c_str()) > 0);
  
    // moved earlier, for faster select_* code
    /* {
       int y, m; 
       double d;
       gmp_sscanf(s_epoch.c_str(),"%d %d %lf",&y,&m,&d);
       workObs->epoch = 
       orsaSolarSystem::FromTimeScale(orsaSolarSystem::gregorTime(y,m,d),
       orsaSolarSystem::TS_UTC);
       // orsa::print(workObs->epoch.getRef());
       }
    */
    //
    workObs->epoch = epoch;
  
    {
        orsaSolarSystem::OpticalObservation * opticalObservation =  
            dynamic_cast<orsaSolarSystem::OpticalObservation *>(workObs.get());
    
        if (opticalObservation) {
      
            {
                int h, m;
                double s;
                gmp_sscanf(s_ra.c_str(),"%d %d %lf",&h,&m,&s);
                Angle tmp; tmp.setHMS(h,m,s);
                opticalObservation->ra = tmp;
                // orsa::print(workObs->ra.getRef());
            }
      
            {
                const int sign = s_decSign == "-" ? -1 : +1;
                int d, p;
                double s;
                gmp_sscanf(s_dec.c_str(),"%d %d %lf",&d,&p,&s);
                Angle tmp; tmp.setDPS(d,p,s,sign);
                opticalObservation->dec = tmp;
                // orsa::print(workObs->dec.getRef());
            }
      
            if (strlen(s_mag.c_str()) > 0) {
                opticalObservation->mag  = atof(s_mag.c_str());
                opticalObservation->band = c_band;
            }
        }
    }
  
    if ( ((MPC_packedNumber(s_number) != 0) || (s_designation != "")) && 
         (s_obsCode != "") &&
         (strlen(s_obsCode.c_str())) == 3) {
        if ( (isalnum(s_obsCode[0])) &&
             (isalnum(s_obsCode[1])) &&
             (isalnum(s_obsCode[2]))) {
            _data.push_back(workObs);
        } else {
            // ORSA_DEBUG("--OUT--");
            return false;
        }
    } else {
        // ORSA_DEBUG("--OUT--");
        return false;
    }
  
    return true;
}

bool MPCObservationsFile::processLines(const char * line1, const char * line2) {
    twoLinesCall=true;
    if (!processLine(line1)) {
        twoLinesCall=false;
        return false;
    }
    std::string s_number, s_designation, s_note1, s_note2;
    std::string s_epoch;
    std::string s_obsCode;
    // check that designation, epoch, and obscode are the same on both lines
    s_number.assign(line2,0,5); 
    removeLeadingAndTrailingSpaces(s_number);
    if (_data[_data.size()-1]->number.isSet()) {
        if (MPC_packedNumber(s_number) != _data[_data.size()-1]->number.getRef()) {
            // ORSA_DEBUG("found data inconsistency");
            twoLinesCall=false;
            return false;
        }
    }
    s_designation.assign(line2,5,7); 
    removeLeadingAndTrailingSpaces(s_designation);
    if (_data[_data.size()-1]->designation.isSet()) {
        if (s_designation != _data[_data.size()-1]->designation.getRef()) {
            // ORSA_DEBUG("found data inconsistency");
            twoLinesCall=false;
            return false;
        }
    }
    s_note1.assign(line2,13,1);
    s_note2.assign(line2,14,1);
    s_epoch.assign(line2,15,17);
    {
        int y, m; 
        double d;
        gmp_sscanf(s_epoch.c_str(),"%d %d %lf",&y,&m,&d);
        const orsa::Time epoch =
            orsaSolarSystem::FromTimeScale(orsaSolarSystem::gregorTime(y,m,d),
                                           orsaSolarSystem::TS_UTC);
        if (epoch != _data[_data.size()-1]->epoch.getRef()) {
            // ORSA_DEBUG("found data inconsistency:");
            /* orsa::print(epoch);
               orsa::print(_data[_data.size()-1]->epoch.getRef());
               ORSA_DEBUG("line1: [%s]",line1);
               ORSA_DEBUG("line2: [%s]",line2);
            */
            twoLinesCall=false;
            return false;
        }
    }
    s_obsCode.assign(line2,77,3);
    removeLeadingAndTrailingSpaces(s_obsCode);
    if (s_obsCode != _data[_data.size()-1]->obsCode.getRef()) {
        // ORSA_DEBUG("found data inconsistency");
        twoLinesCall=false;
        return false;
    }
    switch (s_note2[0]) {
        case 'r': 
            break;
        case 'v':
            break;
        case 's':
        {
            orsaSolarSystem::SatelliteObservation * satelliteObservation =  
                dynamic_cast<orsaSolarSystem::SatelliteObservation *>(_data[_data.size()-1].get());
            if (!satelliteObservation) {
                ORSA_DEBUG("observation is not a SatelliteObservation");
            }
            std::string s_12;
            std::string s_sign_x,s_sign_y,s_sign_z;
            std::string s_vx,s_vy,s_vz;
            //
            s_12.assign(line2,32,1);
            //
            s_sign_x.assign(line2,34,1);
            s_vx.assign(line2,35,10);
            //
            s_sign_y.assign(line2,46,1);
            s_vy.assign(line2,47,10);
            //
            s_sign_z.assign(line2,58,1);
            s_vz.assign(line2,59,10);
            //
            double vx,vy,vz;
            const int sx = s_sign_x == "-" ? -1 : +1;
            const int sy = s_sign_y == "-" ? -1 : +1;
            const int sz = s_sign_z == "-" ? -1 : +1;
            gmp_sscanf(s_vx.c_str(),"%lf",&vx);
            gmp_sscanf(s_vy.c_str(),"%lf",&vy);
            gmp_sscanf(s_vz.c_str(),"%lf",&vz);
            const int onetwo = atoi(s_12.c_str());
            if (onetwo==1) {
                vx = orsa::FromUnits(sx*vx,orsa::Unit::KM);
                vy = orsa::FromUnits(sy*vy,orsa::Unit::KM);
                vz = orsa::FromUnits(sz*vz,orsa::Unit::KM);
            } else if (onetwo==2) {
                vx = orsa::FromUnits(sx*vx,orsa::Unit::AU);
                vy = orsa::FromUnits(sy*vy,orsa::Unit::AU);
                vz = orsa::FromUnits(sz*vz,orsa::Unit::AU);
            } else {
                ORSA_DEBUG("problem: case not handled");
            }
            satelliteObservation->obsPos = orsa::Vector(vx,vy,vz);
        }
        break;
        default:
            ORSA_DEBUG("should not be here...");
            break;
    }
  
    twoLinesCall=false;
    return true;
}
