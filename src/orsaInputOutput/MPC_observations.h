#ifndef _ORSA_INPUT_OUTPUT_MPC_OBSERVATIONS_
#define _ORSA_INPUT_OUTPUT_MPC_OBSERVATIONS_

#include <orsaInputOutput/file.h>

#include <orsaSolarSystem/observation.h>

#include <vector>

namespace orsaInputOutput {
  
  class MPCObservationsFile : 
  public orsaInputOutput::InputFile <
    orsaInputOutput::CompressedFile,
    std::vector< osg::ref_ptr<orsaSolarSystem::Observation> > 
    > {
    
  public:    
    MPCObservationsFile() : 
      InputFile <
      orsaInputOutput::CompressedFile,
      std::vector< osg::ref_ptr<orsaSolarSystem::Observation> > 
      > () {
      twoLinesCall=false;
    }
      
  protected:
    ~MPCObservationsFile() { }
    
  protected:
    bool twoLinesCall;
    
  public:
    bool goodLine(const char * line);
    
  public:
    bool processLine(const char * line);
    
  public:
    bool processLines(const char * lineAbove, const char * line);
  };
  
}; // orsaInputOutput

#endif // _ORSA_INPUT_OUTPUT_MPC_OBSERVATIONS_
