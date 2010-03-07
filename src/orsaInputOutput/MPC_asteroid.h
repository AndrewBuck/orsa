#ifndef _ORSA_INPUT_OUTPUT_MPC_ASTEROID_
#define _ORSA_INPUT_OUTPUT_MPC_ASTEROID_

#include <orsaInputOutput/file.h>

#include <orsaSolarSystem/orbit.h>

#include <vector>
#include <string>

namespace orsaInputOutput {
  
  class MPCAsteroidDataElement {
  public:
    orsa::Cache<orsaSolarSystem::OrbitWithEpoch> orbit;
    orsa::Cache<double> H;
    orsa::Cache<unsigned int> number;
    orsa::Cache<std::string> designation;
  };
  
  typedef std::vector<MPCAsteroidDataElement> MPCAsteroidData;
  
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
