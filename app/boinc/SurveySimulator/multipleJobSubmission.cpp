#include <stdlib.h>

#include <iostream>
#include <string>

using namespace std;

int main (int argc, char ** argv) {
  
  if (argc != 7) {
    cout << "Usage: " << argv[0] << " baseName numJobs realNEORandomSeed numRealNEO syntheticNEORandomSeed numSyntheticNEO" << endl;
    exit(0);
  }
  
  const std::string baseName               = argv[1];
  const int         numJobs                = atoi(argv[2]);
  const int         realNEORandomSeed      = atoi(argv[3]);
  const int         numRealNEO             = atoi(argv[4]);
  const int         syntheticNEORandomSeed = atoi(argv[5]);
  const int         numSyntheticNEO        = atoi(argv[6]);
  
  if (numJobs < 1) {
    cout << "numJobs must be positive" << endl;
    exit(0);
  }
  
  if (numJobs > 9999) {
    cout << "numJobs is too big, are you sure of what you're doing?" << endl;
    exit(0);
  }
  
  if ( (numJobs > 1) &&
       (numSyntheticNEO == 0) ) {
    cout << "generating identical jobs, are you sure of what you're doing?" << endl;
    exit(0);
  }
  
  printf("Submission of %i job(s)\n"
	 "real.....:  seed: %10i  num: %i\n"
	 "synthetic:  seed: %10i  num: %i\n"
	 "proceed? [y/n] ",
	 numJobs,
	 realNEORandomSeed,
	 numRealNEO,
	 syntheticNEORandomSeed,
	 numSyntheticNEO);
  
  const std::string yes = "y";
  
  std::string yesno;
  
  cin >> yesno;
  
  // cout << "yesno: " << yesno << endl;
  
  if (yesno == yes) {
    {
      int waitSec = 5;
      while (waitSec > 0) {
	cout << "Starting in " << waitSec << " seconds..." << endl;
	sleep(1);
	--waitSec;
      }
    }
    char cmd[1024];
    FILE * fp;
    for (int k=0; k<numJobs; ++k) {
      fp = fopen("realNEO.gen","w");
      fprintf(fp,"%i %i\n",numRealNEO,realNEORandomSeed);
      fclose(fp);
      fp = fopen("syntheticNEO.gen","w");
      fprintf(fp,"%i %i\n",numSyntheticNEO,syntheticNEORandomSeed+k);
      fclose(fp);
      snprintf(cmd,1024,"./SurveySimulatorWorkGenerator %s",baseName.c_str());
      cout << "executing: " << cmd << endl;
      system(cmd);
    }
  }
  
  return 0;
}
