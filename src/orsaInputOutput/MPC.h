#ifndef _ORSA_INPUT_OUTPUT_MPC_
#define _ORSA_INPUT_OUTPUT_MPC_

#include <orsa/datetime.h>

#include <string>

namespace orsaInputOutput {
  
    int MPC_charToInt(const char c);
  
    orsa::Time MPC_packedToTime(const std::string & packedEpoch);
  
    unsigned int MPC_packedNumber(const std::string & packedNumber);
  
}; // orsaInputOutput

#endif // _ORSA_INPUT_OUTPUT_MPC_
