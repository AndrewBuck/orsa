#include <orsaInputOutput/MPC.h>

using namespace orsa;
using namespace orsaInputOutput;

#include <orsaSolarSystem/datetime.h>

int orsaInputOutput::MPC_charToInt(const char c) {
    int d;
    const int ch = (int)c;
    switch (ch) {
        case '0': d=0;  break;
        case '1': d=1;  break;
        case '2': d=2;  break;
        case '3': d=3;  break;
        case '4': d=4;  break;
        case '5': d=5;  break;
        case '6': d=6;  break;
        case '7': d=7;  break;
        case '8': d=8;  break;
        case '9': d=9;  break;
        case 'A': d=10; break; 
        case 'B': d=11; break;
        case 'C': d=12; break;
        case 'D': d=13; break;
        case 'E': d=14; break;
        case 'F': d=15; break;
        case 'G': d=16; break;
        case 'H': d=17; break;
        case 'I': d=18; break;
        case 'J': d=19; break;
        case 'K': d=20; break;
        case 'L': d=21; break;
        case 'M': d=22; break;
        case 'N': d=23; break;
        case 'O': d=24; break;
        case 'P': d=25; break;
        case 'Q': d=26; break;
        case 'R': d=27; break;
        case 'S': d=28; break;
        case 'T': d=29; break;
        case 'U': d=30; break;
        case 'V': d=31; break;
        default:
            ORSA_DEBUG("case not handled, c: [%c]",c);
            d=0;
    }
    return d;
}

orsa::Time orsaInputOutput::MPC_packedToTime(const std::string & packedEpoch) {
    const char * s = packedEpoch.c_str();
    if (strlen(s) != 5) {
        ORSA_DEBUG("problems");
        return orsa::Time();
    }
    int y, m, d;
    y  = 100*MPC_charToInt(s[0]);
    y +=  10*MPC_charToInt(s[1]);
    y +=     MPC_charToInt(s[2]);
    m  = MPC_charToInt(s[3]);
    d  = MPC_charToInt(s[4]);
    const orsa::Time t = 
        orsaSolarSystem::FromTimeScale(orsaSolarSystem::gregorTime(y,m,d,0,0,0,0), 
                                       orsaSolarSystem::TS_TDT);
    return t;
}

unsigned int orsaInputOutput::MPC_packedNumber(const std::string & packedNumber) {
    if (strlen(packedNumber.c_str()) == 0) {
        return 0;
    }
    if (strlen(packedNumber.c_str()) != 5) {
        // ORSA_DEBUG("MPC packed number use only 5 digits (arg: [%s])",packedNumber.c_str());
        return 0;
    }
    // if last digit is a 'P', issue an error, as it is probably a comet number
    if (packedNumber[strlen(packedNumber.c_str())-1] == 'P') {
        // ORSA_DEBUG("[%s] interpreted as comet number, returning zero.",packedNumber.c_str());
        return 0;
    }
    unsigned int num=0;
    for (unsigned int k=0; k<strlen(packedNumber.c_str()); ++k) {
        if (!isalnum(packedNumber[k])) {
            ORSA_DEBUG("cannot handle this character: [%s]",packedNumber[k]);
        }
        num *= 10;
        num += MPC_charToInt(packedNumber[k]);
    }
    // ORSA_DEBUG("[%s] -> %i",packedNumber.c_str(),num);
    return num;
}

