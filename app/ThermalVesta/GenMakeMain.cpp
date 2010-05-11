#include <orsa/debug.h>
#include <cstdlib>
#include <libgen.h>
#include <cstdio>

int main(int argc, char ** argv) {
  
  orsa::Debug::instance()->initTimer();
  
  if (argc < 5) {
    printf("Usage: %s <JD-start> <JD-stop> <JD-step> <ThermalVesta-output-file(s)>\n",argv[0]);
    exit(0);
  }
  
  const double JD_start = atof(argv[1]);
  const double JD_stop  = atof(argv[2]);
  const double JD_step  = atof(argv[3]);
  
  // before all, write this rule
  // with this, all intermediate files are kept instead of deleted automatically
  // this is needed to keep the *.allEta.dat files around
  printf(".SECONDARY:\n");
  printf("\n");
  
  // first, the targets
  printf("all:");
  double JD=JD_start;
  while (JD <= JD_stop) {
    printf(" vesta.thermal_%.2f.gif",JD);
    JD += JD_step;
  }
  printf("\n"); 
  printf("\n");
  
  printf("DATA =");
  for (int arg=4; arg<argc; ++arg) {
    printf(" %s",argv[arg]);
  }
  printf("\n"); 
  printf("\n");
  
  // now the rules
  printf("vesta.thermal_%%.gif:\n");
  printf("\t./density_plot $* $(DATA)\n");
  printf("\n");
  
  return 0;
}
