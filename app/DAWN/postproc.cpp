#include <orsa/debug.h>
#include <orsa/unit.h>

#include <stdio.h>
#include <stdlib.h> 

#include <vector>

#include <fftw3.h>

double  norm( const fftw_complex z) {
  return sqrt( z[0] * z[0] + z[1] * z[1] );
}

double  norm_sq( const fftw_complex z) {
  return ( z[0] * z[0] + z[1] * z[1] );
}

void phi(fftw_complex & result, 
	 double omega, 
	 const fftw_complex * in, 
	 const unsigned int size,
	 const bool Hanning=false) {
  
  result[0] = 0;
  result[1] = 0;
  
  double arg,c,s; 
  double window_factor;
  for(unsigned int k=0; k<size; ++k) {  
    arg = orsa::twopi()*k*omega;
    c = cos(arg);
    s = sin(arg);
    if (Hanning) {
      window_factor = (1-cos(k*orsa::twopi()/size));
    } else {
      window_factor = 1;
    }	
    result[0] += window_factor * (in[k][0] * c + in[k][1] * s);
    result[1] -= window_factor * (in[k][0] * s - in[k][1] * c);            
  }
  
  result[0] /= size;
  result[1] /= size;
}

double phi_amp(double omega, fftw_complex * in, const int size, const bool Hanning=false) {
  
  fftw_complex phiPlus;
  phi(phiPlus,  omega,in,size,Hanning);
  fftw_complex phiMinus;
  phi(phiMinus,-omega,in,size,Hanning);
  
  return sqrt( norm_sq(phiPlus) + norm_sq(phiMinus) );
}

class FFTDataStructure {
public:
  double time,amplitude,phase;
};

class FFTDataStream : public std::vector<FFTDataStructure> {
public:
  double timestep;
};

class FFTPowerSpectrumBaseElement {
public:
  double frequency;
  double power;
};

class FFTPowerSpectrum : public std::vector<FFTPowerSpectrumBaseElement> {
};

void FFT(FFTPowerSpectrum & FFTps,
       	 const FFTDataStream & data_stream,
	 const bool Hanning=false,
	 const bool HiRes=false) {
  
  const unsigned int size = data_stream.size();
  
  const double timestep = data_stream.timestep;
  
  if (HiRes) {
    
    fftw_complex * in = (fftw_complex*)malloc(size*sizeof(fftw_complex));
    
    for (unsigned int j=0; j<size; j++) {
      in[j][0] = data_stream[j].amplitude * cos(data_stream[j].phase);
      in[j][1] = data_stream[j].amplitude * sin(data_stream[j].phase);
    }
    
    // high resoultion spectrum
    int resample = 16;
    FFTPowerSpectrumBaseElement tps;
    FFTps.resize(0);
    for (unsigned int l=0;l<(size/2+1)*resample;l++) {
      tps.frequency = ((double)l/(double)(size*resample))/timestep;
      
      // Hanning or flat (NONE) windowing...
      // tps.power = phi_amp(((double)l/(double)(size*resample)),in,size);
      tps.power = phi_amp(((double)l/(double)(size*resample)),in,size,Hanning);
      
      FFTps.push_back(tps);
    }
    
    free(in); 
    
  } else {
    
#warning todo: replace "size" with "N"
    
    fftw_complex *in, *out;
    
    in   = (fftw_complex*)malloc(size*sizeof(fftw_complex));
    out  = (fftw_complex*)malloc(size*sizeof(fftw_complex));
    
    for (unsigned int j=0; j<size; j++) {
      in[j][0] = data_stream[j].amplitude * cos(data_stream[j].phase);
      in[j][1] = data_stream[j].amplitude * sin(data_stream[j].phase);
    }
    
    fftw_plan plan = fftw_plan_dft_1d(size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    
    fftw_execute(plan);
    
    fftw_destroy_plan(plan);
    
    
    std::vector<double> psd;
    psd.resize(size/2+1);
    psd[0] = norm(out[0])/size;
    for (unsigned int k=1; k<(size+1)/2; ++k)
      psd[k] = sqrt(norm_sq(out[k])+norm_sq(out[size-k]))/size;
    if (size % 2 == 0) // size is even
      psd[size/2] = norm(out[size/2])/size;  // Nyquist freq.
  
    
    FFTPowerSpectrumBaseElement tps;
    FFTps.resize(0);
    for (unsigned int l=0;l<psd.size();l++) {
      tps.frequency = ((double)l/(double)size)/timestep;
      tps.power     = psd[l];
      
      FFTps.push_back(tps);
    }
    
    free(in); free(out);
  }
  
}

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

/* class FourierElement {
   public:
   double t;
   double re;
   double im;
   };
*/

int main(int argc, char **argv) {
  
  if (argc != 2) {
    ORSA_DEBUG("Usage: %s <DAWN-extracted-output-file>");
    exit(0);
  }
  
  std::vector<DataElement> data;
  
  
  {
    // read file
    
    FILE * fp = fopen(argv[1],"r");
    
    if (!fp) {
      ORSA_DEBUG("cannot open input file: [%s]",argv[1]);
      exit(0);
    }
    
    char line[1024];
    
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
    
    fclose(fp);
  }
  
  const double T = data[data.size()-1].t - data[0].t;
  //
  const double samplingPeriod    = T / (data.size()-1);
  const double samplingFrequency = 1 / samplingPeriod;
  
  /* ORSA_DEBUG("rows: %i   T: %g [h]",
     data.size(),
     orsa::FromUnits(T,orsa::Unit::HOUR,-1));
  */
  
  
  FFTDataStream dataStream;
  //
  dataStream.timestep = samplingPeriod;
  
  for (unsigned int k=0; k<data.size(); ++k) {
    
    FFTDataStructure row;
    
    row.time      = data[k].t;
    row.amplitude = data[k].e;
    row.phase     = data[k].node+data[k].peri;
    
    dataStream.push_back(row);
  }
  
  FFTPowerSpectrum FFTps;
  const bool Hanning = false;
  const bool HiRes   = true;
  FFT(FFTps,
      dataStream,
      Hanning,
      HiRes);
  
  for (unsigned int k=0; k<FFTps.size(); ++k) {
    ORSA_DEBUG("%g %g",
	       FFTps[k].frequency,
	       FFTps[k].power);
  }
  
  return 0;
}

