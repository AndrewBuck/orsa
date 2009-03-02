#ifndef _ORSA_INPUT_OUTPUT_MPC_COMET_
#define _ORSA_INPUT_OUTPUT_MPC_COMET_

#include <orsaInputOutput/file.h>

#include <orsa/orbit.h>

#include <list>
// #include <map>
#include <string>

namespace orsaInputOutput {
  
  class MPCCometDataElement {
  public:
    orsa::Orbit orbit;
    orsa::Time  epoch;
    std::string designation;
  };
  
  typedef std::list<MPCCometDataElement> MPCCometData;
  
  class MPCCometFile : 
  public orsaInputOutput::InputFile <
    orsaInputOutput::CompressedFile,
    orsaInputOutput::MPCCometData
    > {
    
  public:
    MPCCometFile() : 
      orsaInputOutput::InputFile <
      orsaInputOutput::CompressedFile,
      orsaInputOutput::MPCCometData
      > () { }
      
  protected:
    ~MPCCometFile() { }
    
  public:
    bool goodLine(const char * line);
    
  public:
    bool processLine(const char * line);
  };
  
}; // orsaInputOutput

#endif // _ORSA_INPUT_OUTPUT_MPC_COMET_
