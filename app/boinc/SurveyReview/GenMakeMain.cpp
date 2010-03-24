#include <orsa/debug.h>
#include <stdlib.h>
#include <libgen.h>
#include "skycoverage.h"

int main(int argc, char ** argv) {
  
  orsa::Debug::instance()->initTimer();
  
  if (argc == 1) {
    printf("Usage: %s <Sky-Coverage-File(s)>\n",argv[0]);
    exit(0);
  }
  
  // before all, write this rule
  // with this, all intermediate files are kept instead of deleted automatically
  // this is needed to keep the *.allEta.dat files around
  printf(".SECONDARY:\n");
  printf("\n");
  
  // first, the targets
  printf("all:");
  for (int arg=1; arg<argc; ++arg) {
    printf(" %s/%s.fit.dat",
	   dirname(argv[arg]),
	   SkyCoverage::basename(argv[arg]).c_str());
  }
  printf("\n"); 
  printf("\n");
  
  // now the rules
  printf("%%.allEta.dat: %%.txt\n");
  printf("\t./DetectionEfficiency $*.txt\n");
  printf("\n");
  //
  printf("%%.fit.dat: %%.allEta.dat\n");
  printf("\t./DetectionEfficiencyFit $*.allEta.dat\n");
  printf("\n");
  //
  /* printf("%%.eta.dat: %%.allEta.dat\n");
     printf("\t./DetectionEfficiencyFit $*.allEta.dat\n");
     printf("\n");
  */
  //
  
  return 0;
}
