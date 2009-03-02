// BOINC API
#include <boinc_api.h>
#include <filesys.h>

#include <gsl/gsl_rng.h>

#include <string>
#include <iostream>

int main() {
  
  boinc_init();
  
  std::string resolvedFileName;
  
  // boinc_resolve_filename_s("rng.input.dat", resolvedFileName);
  // FILE * fp_in  = boinc_fopen(resolvedFileName.c_str(), "r");
  
  int randomSeed     = 85719;
  unsigned int nIter = 10000;
  // fscanf(fp_in,"%d %d",&randomSeed,&nIter);
  // fclose(fp_in);
  
  boinc_resolve_filename_s("rng.output.dat", resolvedFileName);
  FILE * fp_out = boinc_fopen(resolvedFileName.c_str(), "w");
  
  // init random generator
  gsl_rng * rnd = gsl_rng_alloc(gsl_rng_gfsr4);
  gsl_rng_set(rnd,randomSeed);
  
  for (unsigned int j=0; j<nIter; ++j) {
    fprintf(fp_out,"%li\n",gsl_rng_get(rnd));
  }
  
  for (unsigned int j=0; j<nIter; ++j) {
    fprintf(fp_out,"%.20f\n",gsl_rng_uniform(rnd));
  }
  
  // rng clean
  gsl_rng_free(rnd);
  
  fclose(fp_out);
  
  boinc_finish(0);
  
  return 0;
}
