#include <orsa/debug.h>
#include <orsa/orbit.h>
#include <orsa/print.h>
#include <orsa/unit.h>

#include <orsaSolarSystem/data.h>

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
  //
  double rot; // reconstructed rotation angle of the asteroid
  //
  std::string line;
};


/* 
   class BasicOrbitProxy : public osg::Referenced {
   public:
   BasicOrbitProxy() : osg::Referenced() { 
   data = new orsa::Interval<OT>;
   data->enableDataStoring();
   }
   protected:
   virtual ~BasicOrbitProxy() { }
   protected:
   class OT {
   public:
   orsa::Orbit o;
   orsa::Time  t;
   public:
   inline bool operator == (const OT & rhs) const {
   return (t == rhs.t);
   }
   public:
   inline bool operator != (const OT & rhs) const {
   return (t != rhs.t);
   }
   public:
   inline bool operator < (const OT & rhs) const {
   return (t < rhs.t);
   }
   public:
   inline bool operator > (const OT & rhs) const {
   return (t > rhs.t);
   }
   public:
   inline bool operator <= (const OT & rhs) const {
   return (t <= rhs.t);
   }
   public:
   inline bool operator >= (const OT & rhs) const {
   return (t >= rhs.t);
   }
   };
   protected:
   osg::ref_ptr< orsa::Interval<OT> > data;
   public:
   bool getOrbit(orsa::Orbit      & orbit,
   const orsa::Time & t) const {
   OT ot,otMin,otMax;
   ot.t = t;
   if (!data->getSubInterval(ot,otMin,otMax)) {
   ORSA_DEBUG("problems...");
   return false;
   }
   if (otMin.t==otMax.t) {
   orbit=otMin.o;
   return true;
   }
   //orsa::print(otMin.t);
   // orsa::print(otMax.t);
   // orsa::print(otMin.o);
   // orsa::print(otMax.o);
   if ((t<otMin.t) || (t>otMax.t)) {
   ORSA_DEBUG("problems...");
   return false;  
   }
   const double beta = 
   (t-otMin.t).get_d() / 
   (otMax.t-otMin.t).get_d();
   const double oneMinusBeta = 1-beta;
   // orsa::print(beta);
   // orsa::print(oneMinusBeta);
   ot.o.a =
   oneMinusBeta*otMin.o.a + 
   beta*otMax.o.a;
   ot.o.e =
   oneMinusBeta*otMin.o.e + 
   beta*otMax.o.e;
   ot.o.i =
   oneMinusBeta*otMin.o.i + 
   beta*otMax.o.i;
   double delta_omega_node = otMax.o.omega_node - otMin.o.omega_node;
   if (fabs(otMax.o.omega_node - otMin.o.omega_node + orsa::twopi()) < 
   fabs(delta_omega_node)) delta_omega_node = otMax.o.omega_node - otMin.o.omega_node + orsa::twopi();
   if (fabs(otMax.o.omega_node - otMin.o.omega_node - orsa::twopi()) < 
   fabs(delta_omega_node)) delta_omega_node = otMax.o.omega_node - otMin.o.omega_node - orsa::twopi();
   ot.o.omega_node =
   otMin.o.omega_node +
   beta*delta_omega_node;
   double delta_omega_pericenter = otMax.o.omega_pericenter - otMin.o.omega_pericenter;
   if (fabs(otMax.o.omega_pericenter - otMin.o.omega_pericenter + orsa::twopi()) < 
   fabs(delta_omega_pericenter)) delta_omega_pericenter = otMax.o.omega_pericenter - otMin.o.omega_pericenter + orsa::twopi();
   if (fabs(otMax.o.omega_pericenter - otMin.o.omega_pericenter - orsa::twopi()) < 
   fabs(delta_omega_pericenter)) delta_omega_pericenter = otMax.o.omega_pericenter - otMin.o.omega_pericenter - orsa::twopi();
   ot.o.omega_pericenter =
   otMin.o.omega_pericenter +
   beta*delta_omega_pericenter;
   double delta_M = 
   (otMax.o.M-otMin.o.M) + 
   orsa::twopi()*((t-otMax.t).get_d()/otMax.o.period()-(t-otMin.t).get_d()/otMin.o.period());
   if (fabs(delta_M + orsa::twopi()) < fabs(delta_M)) delta_M += orsa::twopi();
   if (fabs(delta_M - orsa::twopi()) < fabs(delta_M)) delta_M -= orsa::twopi();
   ot.o.M = 
   (otMin.o.M + orsa::twopi()*(t-otMin.t).get_d()/otMin.o.period()) +
   beta*delta_M;
   ot.o.mu =
   oneMinusBeta*otMin.o.mu + 
   beta*otMax.o.mu;
   // normalize
   ot.o.omega_node       = fmod(orsa::twopi()+fmod(ot.o.omega_node,
   orsa::twopi()),orsa::twopi());
   ot.o.omega_pericenter = fmod(orsa::twopi()+fmod(ot.o.omega_pericenter,
   orsa::twopi()),orsa::twopi());
   ot.o.M                = fmod(orsa::twopi()+fmod(ot.o.M,
   orsa::twopi()),orsa::twopi());
   orbit = ot.o;
   // orsa::print(ot.o);
   // orsa::print(orbit);
   return true;
   }
   public:
   bool insertOrbit(const orsa::Orbit & orbit,
   const orsa::Time  & t) {
   OT ot;
   ot.o = orbit;
   ot.t = t;
   return data->insert(ot);
   }
   };
*/

// this version uses the equinoctial elements, easier to interpolate
class BasicOrbitProxy : public osg::Referenced {
public:
  BasicOrbitProxy() : osg::Referenced() { 
    data = new orsa::Interval<EOT>;
    data->enableDataStoring();
  }
protected:
  virtual ~BasicOrbitProxy() { }
protected:
  class EOT {
  public:
    orsa::EquinoctialOrbit eo;
    orsa::Time             t;
  public:
    inline bool operator == (const EOT & rhs) const {
      return (t == rhs.t);
    }
  public:
    inline bool operator != (const EOT & rhs) const {
      return (t != rhs.t);
    }
  public:
    inline bool operator < (const EOT & rhs) const {
      return (t < rhs.t);
    }
  public:
    inline bool operator > (const EOT & rhs) const {
      return (t > rhs.t);
    }
  public:
    inline bool operator <= (const EOT & rhs) const {
      return (t <= rhs.t);
    }
  public:
    inline bool operator >= (const EOT & rhs) const {
      return (t >= rhs.t);
    }
  };
protected:
  osg::ref_ptr< orsa::Interval<EOT> > data;
public:
  bool getOrbit(orsa::Orbit      & orbit,
		const orsa::Time & t) const {
    EOT eot,eotMin,eotMax;
    eot.t = t;
    if (!data->getSubInterval(eot,eotMin,eotMax)) {
      // ORSA_DEBUG("problems...");
      return false;
    }
    if ((t<eotMin.t) || (t>eotMax.t)) {
      // ORSA_DEBUG("problems...");
      return false;  
    }
    if (eotMin.t==eotMax.t) {
      eotMin.eo.get(orbit);
      return true;
    }
    const double beta = 
      (t-eotMin.t).get_d() / 
      (eotMax.t-eotMin.t).get_d();
    const double oneMinusBeta = 1-beta;
    
    // orsa::print(beta);
    // orsa::print(oneMinusBeta);
    
    eot.eo.p = 
      oneMinusBeta*eotMin.eo.p + 
      beta*eotMax.eo.p;
    
    eot.eo.f = 
      oneMinusBeta*eotMin.eo.f + 
      beta*eotMax.eo.f;
    
    eot.eo.g = 
      oneMinusBeta*eotMin.eo.g + 
      beta*eotMax.eo.g;
    
    eot.eo.h = 
      oneMinusBeta*eotMin.eo.h + 
      beta*eotMax.eo.h;
    
    eot.eo.k = 
      oneMinusBeta*eotMin.eo.k + 
      beta*eotMax.eo.k;
    
    // more care for the L angle
    double delta_L = 
      (eotMax.eo.L-eotMin.eo.L) + 
      orsa::twopi()*((t-eotMax.t).get_d()/eotMax.eo.period()-(t-eotMin.t).get_d()/eotMin.eo.period());
    if (fabs(delta_L + orsa::twopi()) < fabs(delta_L)) delta_L += orsa::twopi();
    if (fabs(delta_L - orsa::twopi()) < fabs(delta_L)) delta_L -= orsa::twopi();
    eot.eo.L = 
      (eotMin.eo.L + orsa::twopi()*(t-eotMin.t).get_d()/eotMin.eo.period()) +
      beta*delta_L;
    
    eot.eo.mu =
      oneMinusBeta*eotMin.eo.mu + 
      beta*eotMax.eo.mu;
    
    // normalize
    eot.eo.L = fmod(orsa::twopi()+fmod(eot.eo.L,
				       orsa::twopi()),orsa::twopi());
    eot.eo.get(orbit);
    return true;
  }
public:
  bool insertOrbit(const orsa::Orbit & orbit,
		   const orsa::Time  & t) {
    EOT eot;
    eot.eo.set(orbit);
    eot.t  = t;
    return data->insert(eot);
  }
};


int main(int argc, char **argv) {
  
  if (argc != 3) {
    ORSA_DEBUG("Usage: %s <DAWN-extracted-output-file-1> <DAWN-extracted-output-file-2>",argv[0]);
    exit(0);
  }
  
  // keep this in sync with main simulation code 
  const double vestaPeriod = FromUnits(5.342128799,orsa::Unit::HOUR);
  const double vestaMass   = 1.35e-10*orsaSolarSystem::Data::MSun();
  const orsa::Time t0 = orsaSolarSystem::gregorTime(2012,
						    1,
						    1,
						    0,
						    0,
						    0,
						    0);
  
  int fileID = 1;
  
  std::vector< std::vector<DataElement> > data;
  
  while (fileID < argc) {
    
    data.resize(fileID);
    
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
	//
	// asteroid rotation angle
	// keep this in sync with main simulation code
	row.rot = fmod(292.0*orsa::degToRad() + (row.t-t0.get_d())*(orsa::twopi()/vestaPeriod),120*orsa::twopi());
	//
	line[strlen(line)-1] = '\0'; // remove trailing \n
	row.line = line;
	
	data[fileID-1].push_back(row);
      }
      
      fclose(fp);
    }
    
    ++fileID;
  }
  
  /*** distance ***/
  
  if (1) {
    
    orsa::Orbit orbit;
    orbit.mu = vestaMass*orsa::Unit::G(); // set once for all
    
    osg::ref_ptr<BasicOrbitProxy> basicOrbitProxy1 = new BasicOrbitProxy;
    for (unsigned int k=0; k<data[0].size(); ++k) {
      orbit.a = data[0][k].a;
      orbit.e = data[0][k].e;
      orbit.i = data[0][k].i;
      orbit.omega_node       = data[0][k].node;
      orbit.omega_pericenter = data[0][k].peri;
      orbit.M                = data[0][k].M;
      //
      basicOrbitProxy1->insertOrbit(orbit,orsa::Time(orsa::FromUnits(data[0][k].t,orsa::Unit::MICROSECOND,-1)));
    }
    
    osg::ref_ptr<BasicOrbitProxy> basicOrbitProxy2 = new BasicOrbitProxy;
    for (unsigned int k=0; k<data[1].size(); ++k) {
      orbit.a = data[1][k].a;
      orbit.e = data[1][k].e;
      orbit.i = data[1][k].i;
      orbit.omega_node       = data[1][k].node;
      orbit.omega_pericenter = data[1][k].peri;
      orbit.M                = data[1][k].M;
      //
      basicOrbitProxy2->insertOrbit(orbit,orsa::Time(orsa::FromUnits(data[1][k].t,orsa::Unit::MICROSECOND,-1)));
    }
    
    const orsa::Time dt(0,0,10,0,0);
    orsa::Time t = dt; // skip first point at 0
    orsa::Orbit o1, o2;
    orsa::Vector r1, r2, v1, v2;
    while (1) {
      if (!basicOrbitProxy1->getOrbit(o1,t)) break;
      if (!basicOrbitProxy2->getOrbit(o2,t)) break;
      if (!o1.relativePosVel(r1,v1)) break;
      if (!o2.relativePosVel(r2,v2)) break;
      printf("%8.5f %12.3f %12.6f\n",
	     orsa::FromUnits(t.get_d(),orsa::Unit::DAY,-1),
	     orsa::FromUnits((r2-r1).length(),orsa::Unit::METER,-1),
	     orsa::FromUnits(orsa::FromUnits((v2-v1).length(),orsa::Unit::METER,-1),orsa::Unit::SECOND));
      t += dt;
    }
    
  }
  
  return 0;
}

