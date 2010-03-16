#include <orsa/debug.h>
#include <stdlib.h>
#include "skycoverage.h"

int main(int argc, char ** argv) {
  
  orsa::Debug::instance()->initTimer();
  
  if (argc == 1) {
    printf("Usage: %s <Sky-Coverage-File(s)>\n",argv[0]);
    exit(0);
  }
  
  char line[4096];
  char patch[1024];
  
  // first, the targets
  sprintf(line,"all:");
  for (int arg=1; arg<argc; ++arg) {
    sprintf(patch," %s.allEta.dat",SkyCoverage::basename(argv[arg]).c_str());
    strcat(line,patch);
  }
  strcat(line,"\n");
  printf("%s",line);
  
  // now the rules
  printf("%%.allEta.dat: %%.txt\n");
  printf("\t./DetectionEfficiency $*.txt\n");
  
  return 0;
}
