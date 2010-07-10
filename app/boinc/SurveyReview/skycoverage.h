#ifndef SURVEY_REVIEW_SKY_COVERAGE
#define SURVEY_REVIEW_SKY_COVERAGE

#include <orsa/angle.h>
#include <orsa/cache.h>
#include <orsa/double.h>
#include <orsa/datetime.h>
#include <orsa/statistic.h>
#include <orsa/util.h>
#include <orsa/vector.h>

#include <orsaInputOutput/file.h>
#include <orsaInputOutput/MPC_asteroid.h>
#include <orsaInputOutput/MPC_observations.h>
#include <orsaInputOutput/MPC_obscode.h>

#include <orsaSolarSystem/observatory.h>
#include <orsaSolarSystem/orbit.h>

#include <osg/Referenced>

#include <string>
#include <list>
#include <libgen.h>

// single night, single observatory
class SkyCoverage : public osg::Referenced {
public:
    SkyCoverage();
protected:
    virtual ~SkyCoverage();
  
    //! Unit Vectors are in Ecliptic coordinates
  
public:
    static orsa::Vector unitVector(const orsa::Angle & ra,
                                   const orsa::Angle & dec);  
public:
    static void normalize(orsa::Angle & ra,
                          orsa::Angle & dec);  
public:
    static void normalize(double & ra,
                          double & dec);  
public:
    //! translates the observatory code used in sky-coverage file to a standard MPC observatory code
    static std::string alias(const std::string &);
  
public:
    void reset();
  
public:
    bool setField(const double & x1,
                  const double & y1,
                  const double & x2,
                  const double & y2,
                  const double & x3,
                  const double & y3,
                  const double & x4,
                  const double & y4,
                  const double & V);
  
public:
    bool get(const orsa::Vector & u,
             double & V,
             const bool verbose=false) const;
  
public:
    bool get(const orsa::Angle & ra,
             const orsa::Angle & dec,
             double & V,
             const bool verbose=false) const {
        return get(unitVector(ra,dec),V,verbose);
    }
  
public:
    // "fast" version, only to know if the given direction is covered
    bool fastGet(const orsa::Vector & u) const;
    bool fastGet(const orsa::Angle & ra,
                 const orsa::Angle & dec) const {
        return fastGet(unitVector(ra,dec));
    }
  
public:
    // add the time of one observation to the interested field
    bool insertFieldTime(const orsa::Time & epoch,
                         const orsa::Vector & u);
    bool insertFieldTime(const orsa::Time & epoch,
                         const orsa::Angle & ra,
                         const orsa::Angle & dec) {
        return insertFieldTime(epoch,unitVector(ra,dec));
    }
public:  
    // estimate the field time as the average of all the inserted times
    bool getFieldAverageTime(orsa::Time & epoch,
                             const orsa::Vector & u) const;
    bool getFieldAverageTime(orsa::Time & epoch,
                             const orsa::Angle & ra,
                             const orsa::Angle & dec) const {
        return getFieldAverageTime(epoch,unitVector(ra,dec));
    }
public:
    // picks from all the inserted field times
    bool pickFieldTime(orsa::Time & epoch,
                       const orsa::Vector & u,
                       const orsa::RNG * rnd) const;
    bool pickFieldTime(orsa::Time & epoch,
                       const orsa::Angle & ra,
                       const orsa::Angle & dec,
                       const orsa::RNG * rnd) const {
        return pickFieldTime(epoch,unitVector(ra,dec),rnd);
    }
    
public:
    bool writeFieldTimeFile(const std::string & filename) const;
    bool  readFieldTimeFile(const std::string & filename);
    
public:
    //! smallest distance (angle) between the u direction and any of the fields centers
    double minDistance(const orsa::Vector & u,
                       const bool verbose=false) const;
  
public:
    // total area of sky covered, as declared in file, no check for overlaps is performed
    double totalDegSq() const;
    
public:
    class SkyCoverageElement {
    public:
        SkyCoverageElement() {
            // epochStat_JD = new orsa::Statistic<double>;
        }
    public:
        orsa::Vector u_centerField;
        orsa::Vector u_X; // old R.A.
        orsa::Vector u_Y; // old Dec.
        double halfFieldSize_X; // old R.A.
        double halfFieldSize_Y; // old Dec.
        double minScalarProduct;
        double limitingMagnitude;
    public:
        // per-field epoch, to be preferred to global one
        // osg::ref_ptr< orsa::Statistic<double> > epochStat_JD;
        std::vector<orsa::Time> epochVec;
    };
    
protected:
    typedef std::vector<SkyCoverageElement> DataType;
    DataType data;
    
public:
    orsa::Cache<orsa::Time>  epoch;   // local midnight epoch
    orsa::Cache<std::string> obscode;
  
public:
    double eta(const double & V,
               const double & U,
               const double & AM,
               const double & GB,
               const double & GL,
               const double & EB,
               const double & EL) const;
    
public:
    static double eta(const double & V,
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
                      const double & EB,
                      const double & drop_EB,
                      const double & scale_EB,
                      const double & center_EB,
                      const double & EL,
                      const double & scale_EL,
                      const double & shape_EL);
public:
    // nominal eta values, mostly for plotting purposes, no mixing angle
    static double nominal_eta_V(const double & V,
                                const double & V_limit,
                                const double & eta0_V,
                                const double & V0,
                                const double & c_V,
                                const double & w_V);
    static double nominal_eta_U(const double & U,
                                const double & U_limit,
                                const double & w_U);
    static double nominal_eta_AM(const double & AM,
                                 const double & peak_AM,
                                 const double & scale_AM,
                                 const double & shape_AM);
    static double nominal_eta_GB_GL(const double & GB,
                                    const double & drop_GB,
                                    const double & scale_GB,
                                    const double & center_GB,
                                    const double & GL,
                                    const double & scale_GL,
                                    const double & shape_GL);
public:
    // coefficients for efficiency as function of apparent magnitude V
    orsa::Cache<double> V_limit, eta0_V, V0, c_V, w_V;
    // coefficients for efficiency as function of apparent velocity U
    orsa::Cache<double> U_limit, w_U;
    //
    orsa::Cache<double> peak_AM, scale_AM, shape_AM;
    //
    orsa::Cache<double> drop_GB, scale_GB, center_GB;
    orsa::Cache<double> scale_GL, shape_GL;
    //
    orsa::Cache<double> drop_EB, scale_EB, center_EB;
    orsa::Cache<double> scale_EL, shape_EL;
    
public:
    // return filename, stripping path and suffix (after first dot)
    static std::string basename(const std::string & filename);
public:
    static bool processFilename(const std::string & filename_in,
                                orsaInputOutput::MPCObsCodeFile * obsCodeFile,
                                std::string & obsCode,
                                orsa::Time & epoch,
                                int & year,
                                int & dayOfYear);    
};

class SkyCoverageFile : 
    public orsaInputOutput::InputFile < orsaInputOutput::PlainFile, osg::ref_ptr<SkyCoverage> > {
public:
    SkyCoverageFile() : 
        orsaInputOutput::InputFile < orsaInputOutput::PlainFile, osg::ref_ptr<SkyCoverage> > () {
        _data = new SkyCoverage;
    }
protected:
    ~SkyCoverageFile() { }
public:
    void setFileName (const std::string & filename) {
        // first, read obscode and epoch from filename
        osg::ref_ptr<orsaInputOutput::MPCObsCodeFile> obsCodeFile = new orsaInputOutput::MPCObsCodeFile;
        obsCodeFile->setFileName("obscode.dat");
        obsCodeFile->read();
        std::string obsCode;
        orsa::Time epoch;
        int year;
        int dayOfYear;
        if (SkyCoverage::processFilename(filename,
                                         obsCodeFile.get(),
                                         obsCode,
                                         epoch,
                                         year,
                                         dayOfYear)) {
            _data->obscode=obsCode;
            _data->epoch=epoch;
        } else {
            // ORSA_DEBUG("problems...");
        }
        // now, regular call
        orsaInputOutput::InputFile < orsaInputOutput::PlainFile, osg::ref_ptr<SkyCoverage> >::setFileName(filename);
    }
public:
    bool processLine(const char * line) {
        double x1,x2,x3,x4; // ra
        double y1,y2,y3,y4; // dec
        double V;           // limiting mag 
        if (9 == gmp_sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",
                            &x1,&y1,
                            &x2,&y2,
                            &x3,&y3,
                            &x4,&y4,
                            &V)) {
            SkyCoverage::normalize(x1,y1);
            SkyCoverage::normalize(x2,y2);
            SkyCoverage::normalize(x3,y3);
            SkyCoverage::normalize(x4,y4);
            //
            _data->setField(x1,y1,x2,y2,x3,y3,x4,y4,V);
            return true;
        } else {
            return false;
        }
    }
public:
    bool goodLine(const char * /* line */ ) { 
        return true;
    }
};

#endif // SURVEY_REVIEW_SKY_COVERAGE
