#ifndef _ORSA_SOLAR_SYSTEM_OBSERVATION_
#define _ORSA_SOLAR_SYSTEM_OBSERVATION_

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/angle.h>
#include <orsa/cache.h>
#include <orsa/datetime.h>
#include <orsa/double.h>
#include <orsa/vector.h>

#include <string>

namespace orsaSolarSystem {
  
    class Observation : public osg::Referenced {
    public:
        Observation() : Referenced(true) { }
    protected:
        virtual ~Observation() { }
    public:
        orsa::Cache<unsigned int> number;
        orsa::Cache<std::string>  designation, obsCode;
        orsa::Cache<bool>         discovery;
        orsa::Cache<orsa::Time>   epoch;
    };
  
    class OpticalObservation : public orsaSolarSystem::Observation {
    public:
        OpticalObservation() : orsaSolarSystem::Observation()  { }
    protected:
        virtual ~OpticalObservation() { }
    public:
        orsa::Cache<orsa::Angle>  ra, dec;
        orsa::Cache<orsa::Angle>  sigma_ra, sigma_dec;
        orsa::Cache<double>       mag;
        orsa::Cache<char>         band; // or filter: V,R,...
    };
    
    double MPC_band_correction(const char band);
    
    class SatelliteObservation : public orsaSolarSystem::OpticalObservation {
    public:
        SatelliteObservation() : orsaSolarSystem::OpticalObservation() { }
    protected:
        virtual ~SatelliteObservation() { }
    public:
        orsa::Cache<orsa::Vector> obsPos; // Equatorial coordinates
    };
  
    class RovingObservation : public orsaSolarSystem::OpticalObservation {
    public:
        RovingObservation() : orsaSolarSystem::OpticalObservation() { }
    protected:
        virtual ~RovingObservation() { }
        // TO BE CONTINUED...
    };
  
    class RadarObservation : public orsaSolarSystem::Observation {
    public:
        RadarObservation() : orsaSolarSystem::Observation() { }
    protected:
        virtual ~RadarObservation() { }
        // TO BE CONTINUED...
    };
  
    typedef std::vector< osg::ref_ptr<orsaSolarSystem::Observation> > ObservationVector;
  
    typedef std::vector< osg::ref_ptr<orsaSolarSystem::OpticalObservation> > OpticalObservationVector;
  
} // namespace orsaSolarSystem

#endif // _ORSA_SOLAR_SYSTEM_OBSERVATION_
