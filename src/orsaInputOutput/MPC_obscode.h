#ifndef _ORSA_INPUT_OUTPUT_MPC_OBSCODE_
#define _ORSA_INPUT_OUTPUT_MPC_OBSCODE_

#include <orsaInputOutput/file.h>

#include <orsaSolarSystem/observatory.h>

#include <orsa/double.h>

#include <list>
#include <map>
#include <string>

namespace orsaInputOutput {
  
  class MPCObsCodeData {
  public:
    // we could use a hash-map
    std::map<std::string,orsaSolarSystem::Observatory> observatory;
  public:
    std::list<std::string> code;
  };
  
  class MPCObsCodeFile : 
  public orsaInputOutput::InputFile <
    orsaInputOutput::CompressedFile,
    orsaInputOutput::MPCObsCodeData
    > {
    
  public:
    MPCObsCodeFile() : 
      orsaInputOutput::InputFile <
      orsaInputOutput::CompressedFile,
      orsaInputOutput::MPCObsCodeData
      > () { }
      
  protected:
    ~MPCObsCodeFile() { }
    
  public:
    bool goodLine(const char * line);
    
  public:
    bool processLine(const char * line);
  };
  
}; // orsaInputOutput

#endif // _ORSA_INPUT_OUTPUT_MPC_OBSCODE_
