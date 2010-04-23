#include "skycoverage.h"

#include <orsa/double.h>

#include <orsaSolarSystem/obleq.h>

SkyCoverage::SkyCoverage() : osg::Referenced() { }

SkyCoverage::~SkyCoverage() { }

orsa::Vector SkyCoverage::unitVector(const orsa::Angle & ra,
                                     const orsa::Angle & dec) {
    double s_ra, c_ra;
    orsa::sincos(ra.getRad(),
                 &s_ra,
                 &c_ra);
    double s_dec, c_dec;
    orsa::sincos(dec.getRad(),
                 &s_dec, 
                 &c_dec);
    return (orsaSolarSystem::equatorialToEcliptic() * 
            orsa::Vector(c_dec*c_ra,
                         c_dec*s_ra,
                         s_dec));
}

void SkyCoverage::normalize(orsa::Angle & ra,
                            orsa::Angle & dec) {
    double  ra_rad =  ra.getRad();
    double dec_rad = dec.getRad();
    //
    SkyCoverage::normalize(ra_rad,dec_rad);
    //
    ra.setRad(ra_rad);
    dec.setRad(dec_rad);
}

void SkyCoverage::normalize(double & ra_rad,
                            double & dec_rad) {
    if (dec_rad > orsa::halfpi()) {
        dec_rad = orsa::pi() - dec_rad;
        ra_rad += orsa::pi();
    }
    if (dec_rad < -orsa::halfpi()) {
        dec_rad = -dec_rad - orsa::pi();
        ra_rad += orsa::pi();
    }
    if (fabs(dec_rad) > orsa::halfpi()) {
        ORSA_DEBUG("still problems after normalization...");
    }
    // ra after dec because depends on it, if close to halfpi
    ra_rad = fmod(fmod(ra_rad,orsa::twopi())+orsa::twopi(),orsa::twopi());
}

std::string SkyCoverage::alias(const std::string & fileCode) {
    if (fileCode=="CATALINA")   return "703";
    if (fileCode=="CSS")        return "703";
    if (fileCode=="LINEAR")     return "704";
    if (fileCode=="LONEOS")     return "699";
    if (fileCode=="SPACEWATCH") return "291";
    return fileCode;
} 

void SkyCoverage::reset() {
    data.clear();
}

bool SkyCoverage::setField(const double & x1,
                           const double & y1,
                           const double & x2,
                           const double & y2,
                           const double & x3,
                           const double & y3,
                           const double & x4,
                           const double & y4,
                           const double & V) {
  
    // check if data is in standard order and format [should check more]
    if ( (y1 != y2) || 
         (y3 != y4) ) {
        ORSA_DEBUG("data not in standard format");
        return false;
    }
  
    const orsa::Vector u1 = unitVector(x1,y1);
    const orsa::Vector u2 = unitVector(x2,y2);
    const orsa::Vector u3 = unitVector(x3,y3);
    const orsa::Vector u4 = unitVector(x4,y4);
  
    const double field_RA  = acos(u1*u2);
    const double field_DEC = acos(u2*u3);
  
    /* ORSA_DEBUG("field: %f x %f [deg^2] [RAxDEC]",
       orsa::radToDeg()*field_RA,
       orsa::radToDeg()*field_DEC);
    */
  
    SkyCoverageElement e;
  
    // in Ecliptic coords...
    const orsa::Vector northEquatorialPole = orsaSolarSystem::equatorialToEcliptic()*orsa::Vector(0,0,1);
  
    e.u_centerField = (u1+u2+u3+u4).normalized();
    e.u_RA  = orsa::externalProduct(northEquatorialPole,e.u_centerField).normalized();
    e.u_DEC = orsa::externalProduct(e.u_centerField,e.u_RA).normalized();
    e.halfFieldSize_RA  = 0.5*field_RA;
    e.halfFieldSize_DEC = 0.5*field_DEC;
    e.minScalarProduct  = std::min(std::min(e.u_centerField*u1,
                                            e.u_centerField*u2),
                                   std::min(e.u_centerField*u3,
                                            e.u_centerField*u4));
    e.limitingMagnitude = V;
  
    data.push_back(e);
  
    // ORSA_DEBUG("data size: %i   minScalarProduct: %f",data.size(),e.minScalarProduct);
  
    /* 
       {
       // some testing before leaving
       ORSA_DEBUG("consinstency checks");
       }
    */
  
    return true;
}

bool SkyCoverage::get(const orsa::Vector & u,
                      double & V,
                      const bool verbose) const {
    orsa::Cache<double> local_V;
    std::list<SkyCoverageElement>::const_iterator it = data.begin();
    while (it != data.end()) {
        if (verbose) {
            ORSA_DEBUG("acos(u*(*it).u_centerField): %7.3f [deg]",orsa::radToDeg()*acos(u*(*it).u_centerField));
            // orsa::print(u);
            // orsa::print((*it).u_centerField);
        }
        if (u*(*it).u_centerField > (*it).minScalarProduct) {
            const double delta_RA  = fabs(asin(u*(*it).u_RA));
            // if (verbose) ORSA_DEBUG("delta_RA: %f [deg]",orsa::radToDeg()*delta_RA);
            if (delta_RA < (*it).halfFieldSize_RA) {
                const double delta_DEC = fabs(asin(u*(*it).u_DEC));
                // if (verbose) ORSA_DEBUG("delta_DEC: %f [deg]",orsa::radToDeg()*delta_DEC);
                if (delta_DEC < (*it).halfFieldSize_DEC) {
                    // if (verbose)  ORSA_DEBUG("found in one field, V: %f",(*it).limitingMagnitude);
                    if (local_V.isSet()) {
                        local_V = std::max(local_V.getRef(),
                                           (*it).limitingMagnitude);
                    } else {
                        local_V = (*it).limitingMagnitude;
                    }
                }
            }
        }
        ++it;
    }
    if (local_V.isSet()) {
        V = local_V.getRef();
        return true;
    } else {
        return false;
    }
}

bool SkyCoverage::fastGet(const orsa::Vector & u) const {
    orsa::Cache<double> local_V;
    std::list<SkyCoverageElement>::const_iterator it = data.begin();
    while (it != data.end()) {
        if (u*(*it).u_centerField > (*it).minScalarProduct) {
            if (fabs(asin(u*(*it).u_RA)) < (*it).halfFieldSize_RA) {
                if (fabs(asin(u*(*it).u_DEC)) < (*it).halfFieldSize_DEC) {
                    return true;
                }
            }
        }
        ++it;
    }
    return false;
}

bool SkyCoverage::insertFieldTime(const orsa::Time & epoch,
                                  const orsa::Vector & u) {
    // sets time on ALL interested fields, not just the first one found
    bool setSome=false;
    std::list<SkyCoverageElement>::iterator it = data.begin();
    while (it != data.end()) {
        if (u*(*it).u_centerField > (*it).minScalarProduct) {
            if (fabs(asin(u*(*it).u_RA)) < (*it).halfFieldSize_RA) {
                if (fabs(asin(u*(*it).u_DEC)) < (*it).halfFieldSize_DEC) {
                    (*it).epochVec.push_back(epoch);
                    setSome=true;
                }
            }
        }
        ++it;
    }
    return setSome;
}

bool SkyCoverage::getFieldAverageTime(orsa::Time & epoch,
                                      const orsa::Vector & u) const {
    // average on all entries in all fields
    osg::ref_ptr< orsa::Statistic<double> > epochStat_JD = new orsa::Statistic<double>;
    std::list<SkyCoverageElement>::const_iterator it = data.begin();
    while (it != data.end()) {
        if (u*(*it).u_centerField > (*it).minScalarProduct) {
            if (fabs(asin(u*(*it).u_RA)) < (*it).halfFieldSize_RA) {
                if (fabs(asin(u*(*it).u_DEC)) < (*it).halfFieldSize_DEC) {
                    for (unsigned int k=0; k<(*it).epochVec.size(); ++k) {
                        epochStat_JD->insert(orsaSolarSystem::timeToJulian((*it).epochVec[k]));
                    }
                }
            }
        }
        ++it;
    }
    if (epochStat_JD->entries()>0) {
        epoch = orsaSolarSystem::julianToTime(epochStat_JD->average());
        return true;
    } else {
        return false;
    }
}

bool SkyCoverage::pickFieldTime(orsa::Time & epoch,
                                const orsa::Vector & u,
                                const orsa::RNG * rnd) const {
    // average on all entries in all fields
    // osg::ref_ptr< orsa::Statistic<double> > epochStat_JD = new orsa::Statistic<double>;
    std::vector<orsa::Time> epochVec;
    std::list<SkyCoverageElement>::const_iterator it = data.begin();
    while (it != data.end()) {
        if (u*(*it).u_centerField > (*it).minScalarProduct) {
            if (fabs(asin(u*(*it).u_RA)) < (*it).halfFieldSize_RA) {
                if (fabs(asin(u*(*it).u_DEC)) < (*it).halfFieldSize_DEC) {
                    for (unsigned int k=0; k<(*it).epochVec.size(); ++k) {
                        epochVec.push_back((*it).epochVec[k]);
                    }
                }
            }
        }
        ++it;
    }
    if (epochVec.size()>0) {
        epoch = epochVec[rnd->gsl_rng_uniform_int(epochVec.size())];
        return true;
    } else {
        return false;
    }
}

double SkyCoverage::totalDegSq() const {
    double area=0;
    std::list<SkyCoverageElement>::const_iterator it = data.begin();
    while (it != data.end()) {
        area += (*it).halfFieldSize_RA*(*it).halfFieldSize_DEC;
        ++it;
    }
    // factor four because dealing with half-sizes above
    area *= 4*orsa::square(orsa::radToDeg());
    return area;
}

double SkyCoverage::minDistance(const orsa::Vector & u,
                                const bool /* verbose */ ) const {
    double minArc=orsa::pi();
    std::list<SkyCoverageElement>::const_iterator it = data.begin();
    while (it != data.end()) {
        minArc=std::min(minArc,acos(u*(*it).u_centerField));
        ++it;
    }
    return minArc;
}

double SkyCoverage::eta(const double & V,
                        const double & U) const {
    return SkyCoverage::eta(V,
                            V_limit.getRef(),
                            eta0_V.getRef(),
                            V0.getRef(),
                            c_V.getRef(),
                            w_V.getRef(),
                            U,
                            U_limit.getRef(),
                            w_U.getRef());
}

double SkyCoverage::eta(const double & V,
                        const double & V_limit,
                        const double & eta0_V,
                        const double & V0,
                        const double & c_V,
                        const double & w_V,
                        const double & U,
                        const double & U_limit,
                        const double & w_U) {
    
/* const double & beta,
   const double & GL,
   const double & GB,
   const double & GB_limit,
   const double & w_GB,
   const double & Gmix) {
*/
    double retVal;
    if (V<V0) {
        retVal = eta0_V;
    } else {
        retVal = 
            (eta0_V-c_V*orsa::square(V-V0)) / 
            (1.0+exp((V-V_limit)/w_V)) / 
            (1.0+exp((fabs(U_limit)-U)/w_U));
        /* retVal = 
           (eta0_V-c_V*orsa::square(V-V0)) / 
           (1.0+exp( cos(beta)*(V-V_limit)/w_V + sin(beta)*(fabs(U_limit)-U)/w_U)) / 
           (1.0+exp(-sin(beta)*(V-V_limit)/w_V + cos(beta)*(fabs(U_limit)-U)/w_U)) /
           (1.0+exp((fabs(GB_limit)-fabs(GB)-Gmix*fabs(GL))/w_GB));
        */
    }
    if (retVal < 0.0) retVal=0.0;
    // if (retVal > 1.0) retVal=1.0;
    return retVal;
}

double SkyCoverage::nominal_eta_V(const double & V,
                                  const double & V_limit,
                                  const double & eta0_V,
                                  const double & V0,
                                  const double & c_V,
                                  const double & w_V) {
    double retVal;
    if (V<V0) {
        retVal = eta0_V;
    } else {
        retVal = 
            (eta0_V-c_V*orsa::square(V-V0)) / 
            (1.0+exp((V-V_limit)/w_V));
    }
    if (retVal < 0.0) retVal=0.0;
    // if (retVal > 1.0) retVal=1.0;
    return retVal;
}

double SkyCoverage::nominal_eta_U(const double & U,
                                  const double & U_limit,
                                  const double & w_U) {
    return (1.0/(1.0+exp((fabs(U_limit)-U)/w_U)));
}

/* double SkyCoverage::nominal_eta_GB(const double & GL,
   const double & GB,
   const double & GB_limit,
   const double & w_GB,
   const double & Gmix) {
   return 1.0/(1.0+exp((fabs(GB_limit)-fabs(GB)-Gmix*fabs(GL))/w_GB));
   }
*/

std::string SkyCoverage::basename(const std::string & filename) {
    const size_t found_last_slash = std::string(filename).find_last_of("//");
    const size_t found_dot = std::string(filename).find(".",(found_last_slash==std::string::npos?0:found_last_slash+1));
    // ORSA_DEBUG("[%s] -> last_slash: %i first dot after slash: %i",filename.c_str(),found_last_slash,found_dot);
    if (found_dot == std::string::npos) {
        ORSA_DEBUG("not regular filename: %s",filename.c_str());
        exit(0);
    }
    std::string s;
    if (found_last_slash!=std::string::npos) {
        s.assign(filename,found_last_slash+1,found_dot-found_last_slash-1);
    } else {
        s.assign(filename,0,found_dot);
    }
    // ORSA_DEBUG("returning: [%s]",s.c_str());
    return s;
}

bool SkyCoverage::processFilename(const std::string & filename_in,
                                  orsaInputOutput::MPCObsCodeFile * obsCodeFile,
                                  std::string & obsCode,
                                  orsa::Time & epoch,
                                  int & year,
                                  int & dayOfYear) {
  
    char filename[1024];
    sprintf(filename,"%s",filename_in.c_str());
  
    // extract observatory and date from input file name
    size_t found_underscore = std::string(::basename(filename)).find("_",0);
    size_t found_dot        = std::string(::basename(filename)).find(".",0);
    if (found_underscore == std::string::npos) {
        ORSA_DEBUG("no underscore found in filename: %s",::basename(filename));
        return false;
    }
    if (found_dot == std::string::npos) {
        // filename without dot, using full file size
        found_dot = strlen(filename);
    }
    // ORSA_DEBUG("found: %i",found);
    std::string compactDate;
    obsCode.assign(::basename(filename),0,found_underscore);
    compactDate.assign(::basename(filename),found_underscore+1,found_dot-found_underscore-1);
    //
    // ORSA_DEBUG("    obsCode: [%s]",obsCode.c_str());
    // ORSA_DEBUG("compactDate: [%s]",compactDate.c_str());
  
    // translate file obscode to MPC standard obscode
    obsCode = SkyCoverage::alias(obsCode);
    // ORSA_DEBUG("MPC obsCode: [%s]",obsCode.c_str());
  
    /* osg::ref_ptr<orsaInputOutput::MPCObsCodeFile> obsCodeFile = new orsaInputOutput::MPCObsCodeFile;
       obsCodeFile->setFileName("obscode.dat");
       obsCodeFile->read();
     
       osg::ref_ptr<orsaSolarSystem::StandardObservatoryPositionCallback> obsPosCB =
       new orsaSolarSystem::StandardObservatoryPositionCallback(obsCodeFile.get());
    */
  
    const orsaSolarSystem::Observatory & observatory = 
        obsCodeFile->_data.observatory[obsCode];
  
    // local midnight epoch
    if (strlen(compactDate.c_str())==7) {
        // seven character date format
        std::string s_year,s_dayOfYear;
        s_year.assign(compactDate,0,4);
        s_dayOfYear.assign(compactDate,4,3);
        year = atoi(s_year.c_str());
        dayOfYear = atoi(s_dayOfYear.c_str());
        epoch = orsaSolarSystem::gregorTime(year,
                                            1,
                                            dayOfYear+1.0-observatory.lon.getRef()/orsa::twopi());
        // orsa::print(epoch);
    } else {
        ORSA_DEBUG("problems...");
        return false;
    }   
    return true;
}


