#include <orsa/debug.h>
#include <orsa/unit.h>

#include <stdio.h>
#include <stdlib.h> 

#include <vector>

class DataElement {
public:
  double t;
  double a;
  double e;
  double i;
  double node;
  double peri;
  double M;
};

class FourierElement {
public:
  double t;
  double re;
  double im;
};

int main(int argc, char **argv) {
  
  if (argc != 2) {
    ORSA_DEBUG("Usage: %s <DAWN-extracted-output-file>");
    exit(0);
  }
  
  FILE * fp = fopen(argv[1],"r");
  
  if (!fp) {
    ORSA_DEBUG("cannot open input file: [%s]",argv[1]);
    exit(0);
  }
  
  char line[1024];
  
  std::vector<DataElement> data;
  
  while (fgets(line,1024,fp)) {
    DataElement row;
    sscanf(line,"%*s %*s %*s %lf %lf %lf %lf %lf %lf %lf",
	   &row.t,
	   &row.a,
	   &row.e,
	   &row.i,
	   &row.node,
	   &row.peri,
	   &row.M);
    // convert units
    row.t = orsa::FromUnits(row.t,orsa::Unit::SECOND);
    row.a = orsa::FromUnits(row.a,orsa::Unit::KM);
    // row.e is fine
    row.i    *= orsa::degToRad();
    row.node *= orsa::degToRad();
    row.peri *= orsa::degToRad();
    row.M    *= orsa::degToRad();
    
    data.push_back(row);
  }
  
  const double     T = data[data.size()-1].t - data[0].t;
  
  /* ORSA_DEBUG("rows: %i   T: %g [h]",
     data.size(),
     orsa::FromUnits(T,orsa::Unit::HOUR,-1));
  */
  
  std::vector<FourierElement> FourierData;
  for (unsigned int k=0; k<data.size(); ++k) {
    
    FourierElement fe;
    
    fe.t = data[k].t;
    //
    fe.re = data[k].e*cos(data[k].node+data[k].peri);
    fe.im = data[k].e*sin(data[k].node+data[k].peri);
    
    FourierData.push_back(fe);
  }
  
#warning remember to apply some windowing...

#warning check normalizations!!
  
  double freq=0;
  const double freqMax  = data.size()/T;
  const double freqStep = freqMax/(32*data.size());
  const double norm     = 1.0/(orsa::twopi()*sqrt(data.size()));
  while (freq <= freqMax) {
    double re = 0;
    double im = 0;
    for (unsigned int k=0; k<FourierData.size(); ++k) {
      const FourierElement & z = FourierData[k];
      const double c = cos(-orsa::twopi()*freq*z.t);
      const double s = sin(-orsa::twopi()*freq*z.t);
      re += z.re*c-z.im*s;
      im += z.re*s+z.im*c;
    }
    /* ORSA_DEBUG("%g %g",
       3600.0*freq,
       norm*sqrt(re*re+im*im));
    */
    ORSA_DEBUG("%g %g",
	       1/(3600.0*freq), // period
	       norm*sqrt(re*re+im*im));
    freq += freqStep;
  }
  
  fclose(fp);
  
  return 0;
}

