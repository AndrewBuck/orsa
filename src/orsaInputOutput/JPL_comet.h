#ifndef _ORSA_INPUT_OUTPUT_JPL_COMET_
#define _ORSA_INPUT_OUTPUT_JPL_COMET_

#include <orsaInputOutput/file.h>

#include <orsa/orbit.h>

#include <list>
// #include <map>
#include <string>

namespace orsaInputOutput {
  
  class JPLCometDataElement {
  public:
    orsa::Orbit orbit;
    orsa::Time  epoch;
    std::string designation;
  };
  
  typedef std::list<JPLCometDataElement> JPLCometData;
  
  class JPLCometFile : 
  public orsaInputOutput::InputFile <
    orsaInputOutput::CompressedFile,
    orsaInputOutput::JPLCometData
    > {
    
  public:
    JPLCometFile() : 
      orsaInputOutput::InputFile <
      orsaInputOutput::CompressedFile,
      orsaInputOutput::JPLCometData
      > () { }
      
  protected:
    ~JPLCometFile() { }
    
  public:
    bool goodLine(const char * line);
    
  public:
    bool processLine(const char * line);
  };
  
}; // orsaInputOutput

#endif // _ORSA_INPUT_OUTPUT_JPL_COMET_
