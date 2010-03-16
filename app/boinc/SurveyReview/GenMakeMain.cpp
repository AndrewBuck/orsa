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
    sprintf(patch," %s.fit.dat",SkyCoverage::basename(argv[arg]).c_str());
    strcat(line,patch);
  }
  strcat(line,"\n");
  printf("%s",line);
  
  printf("\n");
  
  // now the rules
  printf("%%.allEta.dat: %%.txt\n");
  printf("\t./DetectionEfficiency $*.txt\n");
  printf("\n");
  printf("%%.fit.dat: %%.allEta.dat\n");
  printf("\t./DetectionEfficiencyFit $*.allEta.dat\n");
  printf("\n");
  printf("%%.eta.dat: %%.allEta.dat\n");
  printf("\t./DetectionEfficiencyFit $*.allEta.dat\n");
  printf("\n");
  
  return 0;
}
