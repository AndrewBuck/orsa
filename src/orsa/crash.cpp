#include <orsa/crash.h>

#include <orsa/debug.h>

void orsa::crash() {
    ORSA_ERROR("time to die!");
    double * q = 0; q[2000000000] = 0; // voluntary segfault, useful for debugging purposes
}
