#ifndef _ORSA_BOINC_UTIL_H_
#define _ORSA_BOINC_UTIL_H_

// BOINC API
#include <boinc_api.h>
#include <filesys.h>

#include <fcntl.h>

#ifdef _WIN32
#define fsync(fd) _commit(fd) 
#endif

#ifdef _WIN32
#define MODEOL "\n"
const unsigned int lineSkip = 2;
#else 
#define MODEOL "\r\n"
const unsigned int lineSkip = 1;
#endif

const int open_flag = O_RDWR | O_CREAT;

#ifdef _WIN32
const int open_mode = S_IWUSR | S_IRUSR;
#else
const int open_mode = S_IWUSR | S_IRUSR | S_IRGRP | S_IROTH;
#endif

#endif // _ORSA_BOINC_UTIL_H_
