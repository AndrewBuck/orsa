#ifndef _ORSA_INPUT_OUTPUT_JPL_ASTEROID_
#define _ORSA_INPUT_OUTPUT_JPL_ASTEROID_

#include <orsaInputOutput/file.h>

#include <orsa/cache.h>
#include <orsa/orbit.h>

#include <vector>
#include <string>

namespace orsaInputOutput {
  
    class JPLAsteroidDataElement {
    public:
        orsa::Orbit  orbit;
        orsa::Time   epoch;
        double H;
        orsa::Cache<mpz_class> number;
        std::string  designation;
    };
  
    typedef std::vector<JPLAsteroidDataElement> JPLAsteroidData;
  
    class JPLNumberedAsteroidFile : 
        public orsaInputOutput::InputFile <
        orsaInputOutput::CompressedFile,
        orsaInputOutput::JPLAsteroidData
        > {
    
    public:
        JPLNumberedAsteroidFile() : 
            orsaInputOutput::InputFile <
        orsaInputOutput::CompressedFile,
        orsaInputOutput::JPLAsteroidData
        > () { }
      
    protected:
        ~JPLNumberedAsteroidFile() { }
    
    public:
        bool goodLine(const char * line);
    
    public:
        bool processLine(const char * line);
    };
  
    class JPLUnnumberedAsteroidFile : 
        public orsaInputOutput::InputFile <
        orsaInputOutput::CompressedFile,
        orsaInputOutput::JPLAsteroidData
        > {
    
    public:
        JPLUnnumberedAsteroidFile() : 
            orsaInputOutput::InputFile <
        orsaInputOutput::CompressedFile,
        orsaInputOutput::JPLAsteroidData
        > () { }
      
    protected:
        ~JPLUnnumberedAsteroidFile() { }
    
    public:
        bool goodLine(const char * line);
    
    public:
        bool processLine(const char * line);
    };
  
}; // orsaInputOutput

#endif // _ORSA_INPUT_OUTPUT_JPL_ASTEROID_
