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
    if (fileCode=="NEAT")       return "644";
    if (fileCode=="SPACEWATCH") return "291";
    if (fileCode=="WISE")       return "C51";
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
    
    /* Here, the only assumption is that the four points are numbered consecutively,
     * either CW or CCW, without cutting trought the diagonal, i.e.:
     *
     *  2------3    2------1
     *  |      |    |      |
     *  |      | OR |      |
     *  |      |    |      |
     *  1------4    3------4
     *
     */
    
    // temporary vectors
    orsa::Vector s1 = unitVector(x1,y1);
    orsa::Vector s2 = unitVector(x2,y2);
    orsa::Vector s3 = unitVector(x3,y3);
    orsa::Vector s4 = unitVector(x4,y4);
    
    // check assumption: diagonals must be longer than sides
    {
        // using scalar products instead of acos(...), so min and max are inverted
        const double min_side = std::min(std::min(s1*s2,s3*s4),std::min(s2*s3,s1*s4));
        const double max_diag = std::max(s1*s3,s2*s4);
        if (max_diag > min_side) {
            // ORSA_DEBUG("error: data points not numbered consecutively along field perimeter)");
            /* ORSA_DEBUG("s1 -- phi: %6.2f [deg]  theta: %+6.2f [deg]",atan2(s1.getY(),s1.getX())*orsa::radToDeg(),asin(s1.getZ())*orsa::radToDeg());
               ORSA_DEBUG("s2 -- phi: %6.2f [deg]  theta: %+6.2f [deg]",atan2(s2.getY(),s2.getX())*orsa::radToDeg(),asin(s2.getZ())*orsa::radToDeg());
               ORSA_DEBUG("s3 -- phi: %6.2f [deg]  theta: %+6.2f [deg]",atan2(s3.getY(),s3.getX())*orsa::radToDeg(),asin(s3.getZ())*orsa::radToDeg());
               ORSA_DEBUG("s4 -- phi: %6.2f [deg]  theta: %+6.2f [deg]",atan2(s4.getY(),s4.getX())*orsa::radToDeg(),asin(s4.getZ())*orsa::radToDeg());
            */
            return false;
        }
    }
    
    if (!obscode.isSet()) {
        ORSA_DEBUG("problems: obscode is not set...");
        return false;
    }
    if (obscode.getRef()=="C51") {
        // fix for C51, where the field width in the direction of the ecliptic longitude
        // is not scaled by cos(eclipticLatitude), so we do it here
        // NOTE: code still has some problems when the field includes one of the ecliptic poles
        const double phi1=atan2(s1.getY(),s1.getX());
        const double phi2=atan2(s2.getY(),s2.getX());
        const double phi3=atan2(s3.getY(),s3.getX());
        const double phi4=atan2(s4.getY(),s4.getX());
        //
        const double theta1=asin(s1.getZ());
        const double theta2=asin(s2.getZ());
        const double theta3=asin(s3.getZ());
        const double theta4=asin(s4.getZ());
        //
        double delta_phi_2_3=phi3-phi2;
        if (fabs(delta_phi_2_3+orsa::twopi()) < fabs(delta_phi_2_3)) delta_phi_2_3+=orsa::twopi();
        if (fabs(delta_phi_2_3-orsa::twopi()) < fabs(delta_phi_2_3)) delta_phi_2_3-=orsa::twopi();
        delta_phi_2_3/=cos(0.5*(theta2+theta3));
        const double new_phi3=0.5*(phi2+phi3+delta_phi_2_3);
        const double new_phi2=0.5*(phi2+phi3-delta_phi_2_3);
        //
        double delta_phi_4_1=phi1-phi4;
        if (fabs(delta_phi_4_1+orsa::twopi()) < fabs(delta_phi_4_1)) delta_phi_4_1+=orsa::twopi();
        if (fabs(delta_phi_4_1-orsa::twopi()) < fabs(delta_phi_4_1)) delta_phi_4_1-=orsa::twopi();
        delta_phi_4_1/=cos(0.5*(theta4+theta1));
        //
        const double new_phi1=0.5*(phi4+phi1+delta_phi_4_1);
        const double new_phi4=0.5*(phi4+phi1-delta_phi_4_1);
        //
        s1 = orsa::Vector(cos(theta1)*cos(new_phi1),cos(theta1)*sin(new_phi1),sin(theta1));
        s2 = orsa::Vector(cos(theta2)*cos(new_phi2),cos(theta2)*sin(new_phi2),sin(theta2));
        s3 = orsa::Vector(cos(theta3)*cos(new_phi3),cos(theta3)*sin(new_phi3),sin(theta3));
        s4 = orsa::Vector(cos(theta4)*cos(new_phi4),cos(theta4)*sin(new_phi4),sin(theta4));
    }
    
    const orsa::Vector u1 = s1;
    const orsa::Vector u2 = s2;
    const orsa::Vector u3 = s3;
    const orsa::Vector u4 = s4;
    
    /* ORSA_DEBUG("u1 -- phi: %6.2f [deg]  theta: %+6.2f [deg]",atan2(u1.getY(),u1.getX())*orsa::radToDeg(),asin(u1.getZ())*orsa::radToDeg());
       ORSA_DEBUG("u2 -- phi: %6.2f [deg]  theta: %+6.2f [deg]",atan2(u2.getY(),u2.getX())*orsa::radToDeg(),asin(u2.getZ())*orsa::radToDeg());
       ORSA_DEBUG("u3 -- phi: %6.2f [deg]  theta: %+6.2f [deg]",atan2(u3.getY(),u3.getX())*orsa::radToDeg(),asin(u3.getZ())*orsa::radToDeg());
       ORSA_DEBUG("u4 -- phi: %6.2f [deg]  theta: %+6.2f [deg]",atan2(u4.getY(),u4.getX())*orsa::radToDeg(),asin(u4.getZ())*orsa::radToDeg());
    */
    
    /* ORSA_DEBUG("acos(u1*u2): %.3f [deg]",acos(u1*u2)*orsa::radToDeg());
       ORSA_DEBUG("acos(u2*u3): %.3f [deg]",acos(u2*u3)*orsa::radToDeg());
       ORSA_DEBUG("acos(u3*u4): %.3f [deg]",acos(u3*u4)*orsa::radToDeg());
       ORSA_DEBUG("acos(u4*u1): %.3f [deg]",acos(u4*u1)*orsa::radToDeg());
    */
    
    SkyCoverageElement e;
    
    // reference axis in ecliptic coords
    // default: north equatorial pole, for all terrestrial observatories
    // except for C51 (WISE satellite) where it is north ecliptic pole
    /* const orsa::Vector zeta_axis =
       (obscode.getRef()=="C51") ?
       (orsa::Vector(0,0,1)) :
       (orsaSolarSystem::equatorialToEcliptic()*orsa::Vector(0,0,1));
       
       e.u_centerField = (u1+u2+u3+u4).normalized();
       e.u_X = orsa::externalProduct(zeta_axis,e.u_centerField).normalized();
       e.u_Y = orsa::externalProduct(e.u_centerField,e.u_X).normalized();
       
       if (fabs((u2-u1)*e.u_X) > fabs((u3-u2)*e.u_X)) {
       e.halfFieldSize_X = 0.25*(acos(u1*u2)+acos(u3*u4));
       e.halfFieldSize_Y = 0.25*(acos(u2*u3)+acos(u4*u1));
       } else {
       e.halfFieldSize_X = 0.25*(acos(u2*u3)+acos(u4*u1));
       e.halfFieldSize_Y = 0.25*(acos(u1*u2)+acos(u3*u4));
       }
    */
    //    
    // better and more generic, as it does not require to assume a reference axis (equatorial or ecliptic)
    e.u_centerField = (u1+u2+u3+u4).normalized();
    e.u_X = (u4-u1+u3-u2).normalized();
    e.u_Y = (u2-u1+u3-u4).normalized();
    e.halfFieldSize_X = 0.25*(acos(u1*u2)+acos(u3*u4));
    e.halfFieldSize_Y = 0.25*(acos(u2*u3)+acos(u4*u1));
    
    e.minScalarProduct  = std::min(std::min(e.u_centerField*u1,
                                            e.u_centerField*u2),
                                   std::min(e.u_centerField*u3,
                                            e.u_centerField*u4));
    e.limitingMagnitude = V;
    
    data.push_back(e);
    
    // ORSA_DEBUG("field: %.2f x %.2f [deg x deg]",2*e.halfFieldSize_X*orsa::radToDeg(),2*e.halfFieldSize_Y*orsa::radToDeg());
    
    // ORSA_DEBUG("data size: %i   minScalarProduct: %f",data.size(),e.minScalarProduct);
    
    return true;
}

bool SkyCoverage::get(const orsa::Vector & u,
                      double & V,
                      const bool verbose) const {
    orsa::Cache<double> local_V;
    DataType::const_iterator it = data.begin();
    while (it != data.end()) {
        if (verbose) {
            ORSA_DEBUG("acos(u*(*it).u_centerField): %7.3f [deg]",orsa::radToDeg()*acos(u*(*it).u_centerField));
            // orsa::print(u);
            // orsa::print((*it).u_centerField);
        }
        if (u*(*it).u_centerField > (*it).minScalarProduct) {
            const double delta_X  = fabs(asin(u*(*it).u_X));
            // if (verbose) ORSA_DEBUG("delta_X: %f [deg]",orsa::radToDeg()*delta_X);
            if (delta_X < (*it).halfFieldSize_X) {
                const double delta_Y = fabs(asin(u*(*it).u_Y));
                // if (verbose) ORSA_DEBUG("delta_Y: %f [deg]",orsa::radToDeg()*delta_Y);
                if (delta_Y < (*it).halfFieldSize_Y) {
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
    DataType::const_iterator it = data.begin();
    while (it != data.end()) {
        if (u*(*it).u_centerField > (*it).minScalarProduct) {
            if (fabs(asin(u*(*it).u_X)) < (*it).halfFieldSize_X) {
                if (fabs(asin(u*(*it).u_Y)) < (*it).halfFieldSize_Y) {
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
    DataType::iterator it = data.begin();
    while (it != data.end()) {
        if (u*(*it).u_centerField > (*it).minScalarProduct) {
            if (fabs(asin(u*(*it).u_X)) < (*it).halfFieldSize_X) {
                if (fabs(asin(u*(*it).u_Y)) < (*it).halfFieldSize_Y) {
                    bool unique=true;
                    for (unsigned int z=0; z<(*it).epochVec.size(); ++z) {
                        if (epoch == (*it).epochVec[z]) {
                            unique=false;
                            break;
                        }
                    }
                    if (unique) {
                        (*it).epochVec.push_back(epoch);
                        setSome=true;
                    }
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
    DataType::const_iterator it = data.begin();
    while (it != data.end()) {
        if (u*(*it).u_centerField > (*it).minScalarProduct) {
            if (fabs(asin(u*(*it).u_X)) < (*it).halfFieldSize_X) {
                if (fabs(asin(u*(*it).u_Y)) < (*it).halfFieldSize_Y) {
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

bool SkyCoverage::writeFieldTimeFile(const std::string & filename) const {
    FILE * fp = fopen(filename.c_str(),"w");
    if (!fp) {
        ORSA_DEBUG("cannot write file [%s]",filename.c_str());
        return false;
    }
    for (unsigned int k=0; k<data.size(); ++k) {
        for (unsigned int z=0; z<data[k].epochVec.size(); ++z) {
            fprintf(fp,"%3i %.5f\n",k,orsaSolarSystem::timeToJulian(data[k].epochVec[z]));
        }
    }
    fclose(fp);
    return true;
}

bool SkyCoverage::readFieldTimeFile(const std::string & filename) {
    FILE * fp = fopen(filename.c_str(),"r");
    if (!fp) {
        ORSA_DEBUG("cannot read file [%s]",filename.c_str());
        return false;
    }
    unsigned int goodEntries=0;
    char line[1024];
    unsigned int fieldID;
    double JD;
    while (fgets(line,1024,fp)) {
        if (2 == sscanf(line,"%i %lf",&fieldID,&JD)) {
            if (fieldID >= data.size()) {
                ORSA_DEBUG("fieldID outside range");
                continue;
            }
            data[fieldID].epochVec.push_back(orsaSolarSystem::julianToTime(JD));
            ++goodEntries;
        }
    }
    fclose(fp);
    // ORSA_DEBUG("read %i field times",goodEntries);
    return true;
}

bool SkyCoverage::pickFieldTime(orsa::Time & epoch,
                                const orsa::Vector & u,
                                const orsa::RNG * rnd) const {
    // average on all entries in all fields
    // osg::ref_ptr< orsa::Statistic<double> > epochStat_JD = new orsa::Statistic<double>;
    std::vector<orsa::Time> epochVec;
    DataType::const_iterator it = data.begin();
    while (it != data.end()) {
        if (u*(*it).u_centerField > (*it).minScalarProduct) {
            if (fabs(asin(u*(*it).u_X)) < (*it).halfFieldSize_X) {
                if (fabs(asin(u*(*it).u_Y)) < (*it).halfFieldSize_Y) {
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
    DataType::const_iterator it = data.begin();
    while (it != data.end()) {
        area += (*it).halfFieldSize_X*(*it).halfFieldSize_Y;
        ++it;
    }
    // factor four because dealing with half-sizes above
    area *= 4*orsa::square(orsa::radToDeg());
    return area;
}

double SkyCoverage::minDistance(const orsa::Vector & u,
                                const bool /* verbose */ ) const {
    double minArc=orsa::pi();
    DataType::const_iterator it = data.begin();
    while (it != data.end()) {
        minArc=std::min(minArc,acos(u*(*it).u_centerField));
        ++it;
    }
    return minArc;
}

/* 
   double SkyCoverage::eta(const double & V,
   const double & U,
   const double & AM,
   const double & GB,
   const double & GL,
   const double & SA,
   const double & LA,
   const double & LI) const {
   return SkyCoverage::eta(V,
   V_limit.getRef(),
   eta0_V.getRef(),
   V0.getRef(),
   c_V.getRef(),
   w_V.getRef(),
   U,
   U_limit.getRef(),
   w_U.getRef(),
   AM,
   peak_AM.getRef(),
   scale_AM.getRef(),
   shape_AM.getRef(),
   GB,
   drop_GB.getRef(),
   scale_GB.getRef(),
   center_GB.getRef(),
   GL,
   scale_GL.getRef(),
   shape_GL.getRef(),
   SA,
   peak_SA.getRef(),
   scale_SA.getRef(),
   shape_SA.getRef(),
   LA,
   LI,
   LA_LI_limit_const.getRef(),
   LA_LI_limit_linear.getRef(),
   LA_LI_w_const.getRef(),
   LA_LI_w_linear.getRef());
   }
*/

double SkyCoverage::eta(const double & V,
                        const double & U,
                        const double & AM,
                        const double & GB,
                        const double & GL) const {
    return SkyCoverage::eta(V,
                            V_limit.getRef(),
                            eta0_V.getRef(),
                            V0.getRef(),
                            c_V.getRef(),
                            w_V.getRef(),
                            U,
                            U_limit.getRef(),
                            w_U.getRef(),
                            AM,
                            peak_AM.getRef(),
                            scale_AM.getRef(),
                            shape_AM.getRef(),
                            GB,
                            drop_GB.getRef(),
                            scale_GB.getRef(),
                            GL,
                            scale_GL.getRef(),
                            shape_GL.getRef());
}

/* 
   double SkyCoverage::eta(const double & V,
   const double & V_limit,
   const double & eta0_V,
   const double & V0,
   const double & c_V,
   const double & w_V,
   const double & U,
   const double & U_limit,
   const double & w_U,
   const double & AM,
   const double & peak_AM,
   const double & scale_AM,
   const double & shape_AM,
   const double & GB,
   const double & drop_GB,
   const double & scale_GB,
   const double & center_GB,
   const double & GL,
   const double & scale_GL,
   const double & shape_GL,
   const double & SA,
   const double & peak_SA,
   const double & scale_SA,
   const double & shape_SA,
   const double & LA,
   const double & LI,
   const double & LA_LI_limit_const,
   const double & LA_LI_limit_linear,
   const double & LA_LI_w_const,
   const double & LA_LI_w_linear) {
   double retVal =
   nominal_eta_V(V,V_limit,eta0_V,V0,c_V,w_V) *
   nominal_eta_U(U,U_limit,w_U) *
   nominal_eta_AM(AM,peak_AM,scale_AM,shape_AM) *
   nominal_eta_GB_GL(GB,drop_GB,scale_GB,center_GB,GL,scale_GL,shape_GL) *
   nominal_eta_SA(SA,peak_SA,scale_SA,shape_SA) *
   nominal_eta_LA_LI(LA,LI,LA_LI_limit_const,LA_LI_limit_linear,LA_LI_w_const,LA_LI_w_linear);
   if (retVal < 0.0) retVal=0.0;
   if (retVal > 1.0) retVal=1.0;
   return retVal;
   }
*/

double SkyCoverage::eta(const double & V,
                        const double & V_limit,
                        const double & eta0_V,
                        const double & V0,
                        const double & c_V,
                        const double & w_V,
                        const double & U,
                        const double & U_limit,
                        const double & w_U,
                        const double & AM,
                        const double & peak_AM,
                        const double & scale_AM,
                        const double & shape_AM,
                        const double & GB,
                        const double & drop_GB,
                        const double & scale_GB,
                        const double & GL,
                        const double & scale_GL,
                        const double & shape_GL) {
    double retVal =
        nominal_eta_V(V,V_limit,eta0_V,V0,c_V,w_V) *
        nominal_eta_U(U,U_limit,w_U) *
        nominal_eta_AM(AM,peak_AM,scale_AM,shape_AM) *
        nominal_eta_GB_GL(GB,drop_GB,scale_GB,GL,scale_GL,shape_GL);
    if (retVal < 0.0) retVal=0.0;
    if (retVal > 1.0) retVal=1.0;
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
    if (retVal > 1.0) retVal=1.0;
    return retVal;
}

double SkyCoverage::nominal_eta_U(const double & U,
                                  const double & U_limit,
                                  const double & w_U) {
    double retVal = (1.0/(1.0+exp((fabs(U_limit)-U)/w_U)));
    if (retVal < 0.0) retVal=0.0;
    if (retVal > 1.0) retVal=1.0;
    return retVal;
}

double SkyCoverage::nominal_eta_AM(const double & AM,
                                   const double & peak_AM,
                                   const double & scale_AM,
                                   const double & shape_AM) {
    double retVal;
    if (AM<peak_AM) {
        retVal = 1.0;
    } else {
        retVal = 1.0+fabs(shape_AM)-sqrt(orsa::square((AM-peak_AM)/scale_AM)+orsa::square(shape_AM));
    }
    if (retVal < 0.0) retVal=0.0;
    if (retVal > 1.0) retVal=1.0;
    return retVal;
}

double SkyCoverage::nominal_eta_GB_GL(const double & GB,
                                      const double & drop_GB,
                                      const double & scale_GB,
                                      const double & GL,
                                      const double & scale_GL,
                                      const double & shape_GL) {
    double retVal = 1.0-drop_GB*(1+fabs(shape_GL)-sqrt(orsa::square(shape_GL)+orsa::square(GL/scale_GL)))/(1+orsa::square(GB/scale_GB));
    if (retVal < 0.0) retVal=0.0;
    if (retVal > 1.0) retVal=1.0;
    return retVal;
}

/* double SkyCoverage::nominal_eta_SA(const double & SA,
   const double & peak_SA,
   const double & scale_SA,
   const double & shape_SA) {
   double retVal;
   if (SA<peak_SA) {
   retVal = 1.0;
   } else {
   retVal = 1.0+fabs(shape_SA)-sqrt(orsa::square((SA-peak_SA)/scale_SA)+orsa::square(shape_SA));
   }
   if (retVal < 0.0) retVal=0.0;
   if (retVal > 1.0) retVal=1.0;
   return retVal;
   }
*/
#warning TEST!! REWRITE THIS CORRECTLY or REMOVE IT
double SkyCoverage::nominal_eta_SA(const double & ,
                                   const double & ,
                                   const double & ,
                                   const double & ) {
    return 1.0;
}

/* 
   #warning TEST!! REWRITE THIS CORRECTLY
   double SkyCoverage::nominal_eta_SA(const double & SA,
   const double & peak_SA,
   const double & scale_SA,
   const double & shape_SA) {
   double retVal = (1.0/(1.0+exp((SA-peak_SA)/scale_SA)));
   if (retVal < 0.0) retVal=0.0;
   if (retVal > 1.0) retVal=1.0;
   return retVal;
   }
*/

double SkyCoverage::nominal_eta_LA_LI(const double & LA,
                                      const double & LI,
                                      const double & LA_LI_limit_const,
                                      const double & LA_LI_limit_linear,
                                      const double & LA_LI_w_const,
                                      const double & LA_LI_w_linear) {
    const double LA_limit = LA_LI_limit_const + (1.0 - LI) * LA_LI_limit_linear;
    const double LA_w     = LA_LI_w_const     + (1.0 - LI) * LA_LI_w_linear;
    double retVal = 1.0 / (1.0+exp((LA-LA_limit)/LA_w));
    if (retVal < 0.0) retVal=0.0;
    if (retVal > 1.0) retVal=1.0;
    return retVal;
}

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
        // ORSA_DEBUG("no underscore found in filename: %s",::basename(filename));
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
        if (!observatory.moving()) {
            epoch = orsaSolarSystem::gregorTime(year,
                                                1,
                                                dayOfYear+1.0-observatory.lon.getRef()/orsa::twopi());
        } else {
            // tested on C51/WISE data
            epoch = orsaSolarSystem::gregorTime(year,
                                                1,
                                                dayOfYear+0.5);
        }
        // orsa::print(epoch);
    } else {
        // ORSA_DEBUG("problems...");
        return false;
    }   
    return true;
}
