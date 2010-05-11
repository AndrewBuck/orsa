#include <orsa/debug.h>

#include <cstdio>
#include <gmpxx.h>

using namespace orsa;

#ifdef __ADVANCED_TIMER__
static const clockid_t debugClockID = CLOCK_MONOTONIC;
#endif

Debug * Debug::_instance = 0;

Debug::Debug() {
#ifdef __ADVANCED_TIMER__
    tp_init = false;
#endif
}

Debug::~Debug() {
    _instance = 0;
}

void Debug::initTimer() {
#ifdef __ADVANCED_TIMER__
    const int retval = clock_gettime(debugClockID, &tp_zero);
    if (retval == 0) {
        tp_init = true;
    } else {
        ORSA_DEBUG("problems in orsa::Debug::initTimer()...");
    }
#else
    // ORSA_DEBUG("timer not available");
#endif
}

const char * Debug::__head__(const DebugType type) {
    switch(type) {
        case DT_DEBUG:   return "debug:";   break;
        case DT_WARNING: return "warning:"; break;
        case DT_ERROR:   return "error:";   break;
    }
    return "";
}

Debug * Debug::instance() {
    if (_instance == 0) {
        _instance = new Debug;
    }
    return _instance;
}

void Debug::create() {
    mutex.lock();
    queue.push_back(orsa::Debug::DebugMessage());
    mutex.unlock();
}

void Debug::flush() {
    mutex.lock();
    bool queueReady = true;
    QueueType::const_iterator it = queue.begin();
    while (it != queue.end()) {
        if ( (!(*it).headIsSet()) ||
             (!(*it).bodyIsSet()) ) {
            queueReady = false;
            break;
        }
        ++it;
    }
    if (queueReady) {
        char message[2048];
        QueueType::const_iterator it = queue.begin();
        while (it != queue.end()) {
            gmp_sprintf(message,
                        "%s%s",
                        (*it).getHead().c_str(),
                        (*it).getBody().c_str());
            printOut(message);
            ++it;
        }
        queue.clear();
    }
    mutex.unlock();
}

/* // useful when debugging under win32
   void Debug::printOut(const std::string & s) {
   FILE * fp = fopen("debug.txt","a");
   gmp_fprintf(fp,"%s\n",s.c_str());
   fclose(fp);
   }
*/

void Debug::printOut(const std::string & s) {
    gmp_fprintf(stderr,"%s\n",s.c_str());
}

/*
  void Debug::head(const DebugType type, const char * file, const int line) {
  gmp_fprintf(stderr, "ORSA[%s:%i] %s ", file, line, __head__(type));
  }
*/

// with local time!
/*
  void Debug::head(const DebugType type, const char * file, const int line) {
  const time_t tt_now = time(0);
  struct tm * tm_struct = localtime(&tt_now);
  gmp_fprintf(stderr,
  "ORSA[%02i:%02i:%02i][%s:%i] %s ",
  tm_struct->tm_hour,
  tm_struct->tm_min,
  tm_struct->tm_sec,
  file,
  line,
  __head__(type));
  }
*/

// with local time & total time elapsed
/*
  void Debug::head(const DebugType type, const char * file, const int line) {

  const time_t tt_now = time(0);
  struct tm * tm_struct = localtime(&tt_now);

  if (tp_init) {

  struct timespec tp;

  clock_gettime(debugClockID, &tp);

  if (tp.tv_nsec > tp_zero.tv_nsec) {
  tp.tv_sec  -= tp_zero.tv_sec;
  tp.tv_nsec -= tp_zero.tv_nsec;
  } else {
  tp.tv_sec  -= (tp_zero.tv_sec + 1);
  tp.tv_nsec += (1000000000 - tp_zero.tv_nsec);
  }

  // should be nanosec resolution, but we just print the microseconds

  gmp_fprintf(stderr,
  "ORSA[%02i:%02i:%02i][%06i.%06i][%s:%i] %s ",
  tm_struct->tm_hour,
  tm_struct->tm_min,
  tm_struct->tm_sec,
  tp.tv_sec,
  tp.tv_nsec/1000,
  file,
  line,
  __head__(type));

  } else {

  gmp_fprintf(stderr,
  "ORSA[%02i:%02i:%02i][%s:%i] %s ",
  tm_struct->tm_hour,
  tm_struct->tm_min,
  tm_struct->tm_sec,
  file,
  line,
  __head__(type));
  }

  }
*/

// with local time & total time elapsed
void Debug::head(const DebugType type, const char * file, const int line, const char * function) {

    const time_t tt_now = time(0);
    struct tm * tm_struct = localtime(&tt_now);

    char headString[1024];

#ifdef __ADVANCED_TIMER__

    if (tp_init) {

        struct timespec tp;

        clock_gettime(debugClockID, &tp);

        /*
          tp.tv_sec -= tp_zero.tv_sec;
          if (tp.tv_nsec > tp_zero.tv_nsec) {
          tp.tv_nsec -= tp_zero.tv_nsec;
          } else {
          tp.tv_nsec += (1000000000 - tp_zero.tv_nsec);
          --tp.tv_sec;
          }
        */

        if (tp.tv_nsec > tp_zero.tv_nsec) {
            tp.tv_sec  -= tp_zero.tv_sec;
            tp.tv_nsec -= tp_zero.tv_nsec;
        } else {
            tp.tv_sec  -= (tp_zero.tv_sec + 1);
            tp.tv_nsec += (1000000000 - tp_zero.tv_nsec);
        }

        // should be nanosec resolution, but we just print the microseconds

        gmp_sprintf(headString,
                    "ORSA[%02i:%02i:%02i][%06i.%06i][%s:%i|%s] %s ",
                    tm_struct->tm_hour,
                    tm_struct->tm_min,
                    tm_struct->tm_sec,
                    tp.tv_sec,
                    tp.tv_nsec/1000,
                    file,
                    line,
                    function,
                    __head__(type));

    } else {

#endif

        gmp_sprintf(headString,
                    "ORSA[%02i:%02i:%02i][%s:%i|%s] %s ",
                    tm_struct->tm_hour,
                    tm_struct->tm_min,
                    tm_struct->tm_sec,
                    file,
                    line,
                    function,
                    __head__(type));

#ifdef __ADVANCED_TIMER__
    }
#endif

    mutex.lock();
    QueueType::iterator it = queue.begin();
    while (it != queue.end()) {
        if (!(*it).headIsSet()) {
            (*it).setHead(headString);
            break;
        }
        ++it;
    }
    mutex.unlock();
}

void Debug::trace(const char * fmt, ...) {
    va_list list;
    va_start(list, fmt);
    vtrace(fmt, list);
}

/*
  void Debug::vtrace(const char * fmt, std::va_list list) {
  gmp_vfprintf(stderr, fmt, list);
  gmp_fprintf(stderr, "\n");
  }
*/

void Debug::vtrace(const char * fmt, std::va_list list) {
    char traceString[1024];
    gmp_vsprintf(traceString,
                 fmt,
                 list);

    /*
      gmp_sprintf(traceString,
      "%s\n",
      traceString);
    */

    mutex.lock();
    QueueType::iterator it = queue.end();
    while (it != queue.begin()) {
        --it;
        if (!(*it).bodyIsSet()) {
            (*it).setBody(traceString);
            break;
        }
    }
    mutex.unlock();
}
