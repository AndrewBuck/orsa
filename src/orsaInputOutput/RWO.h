#ifndef _ORSA_INPUT_OUTPUT_RWO_
#define _ORSA_INPUT_OUTPUT_RWO_

#include <orsaInputOutput/file.h>

#include <orsaSolarSystem/observation.h>

#include <vector>

namespace orsaInputOutput {
  
  typedef std::vector< osg::ref_ptr<orsaSolarSystem::Observation> > RWOData;
  
  class RWOFile : 
  public orsaInputOutput::InputFile <
    orsaInputOutput::CompressedFile,
    RWOData
    > {
    
  public:    
    RWOFile() : 
      InputFile <
      orsaInputOutput::CompressedFile,
      RWOData 
      > () { }
      
  protected:
    ~RWOFile() { }
    
  public:
    bool goodLine(const char * line);
    
  public:
    bool processLine(const char * line);
  };
  
}; // orsaInputOutput

#endif // _ORSA_INPUT_OUTPUT_RWO_
