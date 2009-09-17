#include <orsa/debug.h>
#include <orsa/unit.h>

#include <stdio.h>
#include <stdlib.h>
 
#include <algorithm>
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
	 const bool HiRes=false,
	 const int resample=4) {
  
  const unsigned int size = data_stream.size();
  
  const double timestep = data_stream.timestep;
  
  if (HiRes) {
    
    fftw_complex * in = (fftw_complex*)malloc(size*sizeof(fftw_complex));
    
    for (unsigned int j=0; j<size; j++) {
      in[j][0] = data_stream[j].amplitude * cos(data_stream[j].phase);
      in[j][1] = data_stream[j].amplitude * sin(data_stream[j].phase);
    }
    
    // high resoultion spectrum
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
  //
  double a;
  double e;
  double i;
  double node;
  double peri;
  double M;
  //
  double R;
};


int main(int argc, char **argv) {
  
  if (argc == 1) {
    ORSA_DEBUG("Usage: %s <DAWN-extracted-output-file(s)>");
    exit(0);
  }
  
  int fileID = 1;
  
  while (fileID < argc) {
    
    std::vector<DataElement> data;
    
    char filename[1024];
    
    {
      // read file
      
      FILE * fp = fopen(argv[fileID],"r");
      
      if (!fp) {
	ORSA_DEBUG("cannot open input file: [%s]",argv[fileID]);
	exit(0);
      }
      
      ORSA_DEBUG("working on file [%s]",argv[fileID]);
      
      char line[1024];
      
      while (fgets(line,1024,fp)) {
	DataElement row;
	//           1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16
	sscanf(line,"%*s %*s %*s %lf %lf %lf %lf %lf %lf %lf %*s %*s %*s %*s %*s %lf",
	       &row.t,
	       &row.a,
	       &row.e,
	       &row.i,
	       &row.node,
	       &row.peri,
	       &row.M,
	       &row.R);
	// convert units
	row.t = orsa::FromUnits(row.t,orsa::Unit::SECOND);
	//
	row.a = orsa::FromUnits(row.a,orsa::Unit::KM);
	// row.e is fine
	row.i    *= orsa::degToRad();
	row.node *= orsa::degToRad();
	row.peri *= orsa::degToRad();
	row.M    *= orsa::degToRad();
	//
	row.R = orsa::FromUnits(row.R,orsa::Unit::KM);
	
	data.push_back(row);
      }
      
      fclose(fp);
    }
    
    if (data.size()==0) {
      ++fileID;
      continue;
    }
    

    /*** R ***/
    
    std::vector<double> data_R;
    double mean_R=0;
    for (unsigned int k=0; k<data.size(); ++k) {
      data_R.push_back(data[k].R);
      mean_R += data[k].R;
    }
    mean_R /= data.size();
    std::sort(data_R.begin(),data_R.end());
    const double median_R = data_R[data_R.size()/2];
    //
    const double range_R_50 = data_R[3*data_R.size()/4]   - data_R[data_R.size()/4];
    const double range_R_90 = data_R[19*data_R.size()/20] - data_R[data_R.size()/20];
    const double range_R    = data_R[data_R.size()-1]     - data_R[0];
    
    sprintf(filename,"%s_R.dat",argv[fileID]);
    FILE * fp_R = fopen(filename,"w");
    {
      fprintf(fp_R,
	      "%g %g %g %g %g %g\n",
	      orsa::FromUnits(data[0].a,orsa::Unit::KM,-1),
	      orsa::FromUnits(mean_R,orsa::Unit::KM,-1),
	      orsa::FromUnits(median_R,orsa::Unit::KM,-1),
	      orsa::FromUnits(range_R_50,orsa::Unit::KM,-1),
	      orsa::FromUnits(range_R_90,orsa::Unit::KM,-1),
	      orsa::FromUnits(range_R,orsa::Unit::KM,-1));
    }
    fclose(fp_R);
    
    
    

    /*** FFT ***/
    
    const double T = data[data.size()-1].t - data[0].t;
    //
    const double samplingPeriod    = T / (data.size()-1);
    // const double samplingFrequency = 1 / samplingPeriod;
    
    FFTDataStream dataStream;
    //
    dataStream.timestep = samplingPeriod;
    
    FFTPowerSpectrum FFTps;
    const bool Hanning = true;
    const bool HiRes   = true;
    
    // test
    /* {
       dataStream.clear();
       dataStream.timestep = 1.0;
       for (unsigned int k=0; k<1000; ++k) {
       FFTDataStructure row;
       row.time      = dataStream.timestep*k;
       row.amplitude = 40.0;
       row.phase     = row.time*orsa::twopi()*0.28;
       dataStream.push_back(row);
       }  
       }
    */
    
    /*** FG ***/
    
    {
      dataStream.clear();
      for (unsigned int k=0; k<data.size(); ++k) {
	FFTDataStructure row;
	row.time      = data[k].t;
	row.amplitude = data[k].e;
	row.phase     = data[k].node+data[k].peri;
	dataStream.push_back(row);
      }
    }
    
    FFT(FFTps,
	dataStream,
	Hanning,
	HiRes);
    
    sprintf(filename,"%s_FG_a0.dat",argv[fileID]);
    FILE * fp_FG = fopen(filename,"w");
    for (unsigned int k=0; k<FFTps.size(); ++k) {
      fprintf(fp_FG,
	      "%g %g %g\n",
	      FFTps[k].frequency*3600.0,
	      orsa::FromUnits(data[0].a,orsa::Unit::KM,-1),
	      FFTps[k].power);
    }
    fclose(fp_FG);
    
    /*** HK ***/
    
    if (0) {
      
      {
	dataStream.clear();
	for (unsigned int k=0; k<data.size(); ++k) {
	  FFTDataStructure row;
	  row.time      = data[k].t;
	  row.amplitude = tan(data[k].i/2);
	  row.phase     = data[k].node;
	  dataStream.push_back(row);
	}
      }
      
      FFT(FFTps,
	  dataStream,
	  Hanning,
	  HiRes);
      
      sprintf(filename,"%s_HK_a0.dat",argv[fileID]);
      FILE * fp_HK = fopen(filename,"w");
      for (unsigned int k=0; k<FFTps.size(); ++k) {
	fprintf(fp_HK,
		"%g %g %g\n",
		FFTps[k].frequency*3600.0,
		orsa::FromUnits(data[0].a,orsa::Unit::KM,-1),
		FFTps[k].power);
      }
      fclose(fp_HK);
      
    }
    
    /*** PL ***/
    
    if (0) {
      
      {
	dataStream.clear();
	for (unsigned int k=0; k<data.size(); ++k) {
	  FFTDataStructure row;
	  row.time      = data[k].t;
	  row.amplitude = orsa::FromUnits(data[k].a,orsa::Unit::KM,-1)*(1-orsa::square(data[k].e));
	  row.phase     = data[k].node+data[k].peri+data[k].M;
	  dataStream.push_back(row);
	}
      }
      
      FFT(FFTps,
	  dataStream,
	  Hanning,
	  HiRes);
      
      sprintf(filename,"%s_PL_a0.dat",argv[fileID]);
      FILE * fp_PL = fopen(filename,"w");
      for (unsigned int k=0; k<FFTps.size(); ++k) {
	fprintf(fp_PL,
		"%g %g %g\n",
		FFTps[k].frequency*3600.0,
		orsa::FromUnits(data[0].a,orsa::Unit::KM,-1),
		FFTps[k].power);
      }
      fclose(fp_PL);
      
    }
    
    
    ++fileID;
  }
  
  return 0;
}

