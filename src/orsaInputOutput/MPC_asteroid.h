#ifndef _ORSA_INPUT_OUTPUT_MPC_ASTEROID_
#define _ORSA_INPUT_OUTPUT_MPC_ASTEROID_

#include <orsaInputOutput/file.h>

#include <orsa/orbit.h>

#include <list>
// #include <map>
#include <string>

namespace orsaInputOutput {
  
  class MPCAsteroidDataElement {
  public:
    orsa::Orbit orbit;
    orsa::Time  epoch;
    orsa::Double H;
    std::string designation;
  };
  
  typedef std::list<MPCAsteroidDataElement> MPCAsteroidData;
  
  class MPCAsteroidFile : 
  public orsaInputOutput::InputFile <
    orsaInputOutput::CompressedFile,
    orsaInputOutput::MPCAsteroidData
    > {
    
  public:
    MPCAsteroidFile() : 
      orsaInputOutput::InputFile <
      orsaInputOutput::CompressedFile,
      orsaInputOutput::MPCAsteroidData
      > () { }
      
  protected:
    ~MPCAsteroidFile() { }
    
  public:
    bool goodLine(const char * line);
    
  public:
    bool processLine(const char * line);
  };
  
}; // orsaInputOutput

#endif // _ORSA_INPUT_OUTPUT_MPC_ASTEROID_
