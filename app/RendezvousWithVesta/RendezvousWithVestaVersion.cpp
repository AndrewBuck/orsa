#include "RendezvousWithVestaVersion.h"

#include "BuildCounter.h"

// #include <iostream>

const QString programVersion() {
  QString s;
  s.sprintf("Version 1.0 Build %s. Compiled %s.",
	    __BUILD_COUNTER__,
	    __DATE__);
  return s;
}
