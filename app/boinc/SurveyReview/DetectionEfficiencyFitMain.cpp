#include <orsa/statistic.h>
#include "skycoverage.h"
#include "eta.h"

class EfficiencyStatistics : public orsa::WeightedStatistic<double> {
public:
  EfficiencyStatistics(const double & aux) : 
    orsa::WeightedStatistic<double>(),
    center(aux) { }
public:
  const double center; // V or U
  orsa::Cache<double> fit; // eta_V or eta_U
};

// for sorting
class EfficiencyStatistics_ptr_cmp {
public:
  bool operator() (EfficiencyStatistics * lhs,
		   EfficiencyStatistics * rhs) const {
    return (lhs->center < rhs->center);
  }
};

class CountStatsElement : public osg::Referenced {
public:
  CountStatsElement() : osg::Referenced() {
    Nobs=0;
    Ndsc=0;
    Ntot=0;
  }
protected:
  ~CountStatsElement() { }
public:
  unsigned int Nobs, Ndsc, Ntot;  
};

class CountStats : public osg::Referenced {  
  // first, the classes to handle linear and logarithmic variables
public:
  class Var : public osg::Referenced {
  public:
    Var() : osg::Referenced() { }
  protected:
    ~Var() { }
  public:
    virtual size_t size() const = 0;
    virtual size_t bin(const double x) const = 0;
    virtual double binStart(const size_t bin) const = 0;
    virtual double binStop(const size_t bin) const = 0;
    double binCenter(const size_t bin) const {
      return (0.5*(binStart(bin)+binStop(bin)));
    }
  };
public:	
  class LinearVar : public Var {
  public:
    LinearVar(const double startValue,
	      const double stopValue,
	      const double incrValue) : 
      Var(),
      start(startValue),
      stop(stopValue),
      incr(incrValue) {
    }
  public:
    size_t size() const {
      return (size_t)(1+(stop-start)/incr);
    }
    size_t bin(const double x) const {
      return (size_t)((x-start)/incr);
    }
    double binStart(const size_t bin) const {
      return (start+bin*incr);
    }
    double binStop(const size_t bin) const {
      return (start+(bin+1)*incr);
    }
  protected:
    const double start, stop, incr;
  };
public:	
  class LogarithmicVar : public Var {
  public:
    LogarithmicVar(const double startValue,
		   const double stopValue,
		   const double factorValue) : 
      Var(),
      start(startValue),
      stop(stopValue),
      factor(factorValue) {
    }
  public:
    size_t size() const {
      return (size_t)(1+log(stop/start)/log(factor));
    }
    size_t bin(const double x) const {
      return (size_t)(log(x/start)/log(factor));
    }
    double binStart(const size_t bin) const {
      return (start*orsa::int_pow(factor,bin));
    }
    double binStop(const size_t bin) const {
      return (start*orsa::int_pow(factor,bin+1));
    }
  protected:
    const double start, stop, factor;
  };
  
public:
  CountStats(const std::vector< osg::ref_ptr<Var> > & varDefinition) :
    osg::Referenced(),
    var(varDefinition) {
    size_t totalSize=1;
    for (unsigned int k=0; k<varDefinition.size(); ++k) {
      totalSize *= varDefinition[k]->size();
    }  
    ORSA_DEBUG("totalSize: %i",totalSize);
    data.resize(totalSize);
  }
protected:
  ~CountStats() { }
  
public:
  bool insert(const std::vector<double> & xVector,
	      const bool obs, 
	      const bool dsc) {
    if (xVector.size() != var.size()) {
      ORSA_DEBUG("dimension mismatch");
      return false;
    }
    std::vector<size_t> binVector;
    if (!bin(binVector,xVector)) {
      return false;
    }   
    
    
    if (0) { 
      // test only
      ORSA_DEBUG("TEST: ---");
      ORSA_DEBUG("original binVector:");
      for (unsigned int k=0; k<binVector.size(); ++k) {
	ORSA_DEBUG("binVector[%i] = %i",k,binVector[k]);
      } 
      const size_t idx = index(binVector);
      ORSA_DEBUG("index: %i",idx);
      bin(idx);
      ORSA_DEBUG("reconstructed binVector:");
      for (unsigned int k=0; k<binVector.size(); ++k) {
	ORSA_DEBUG("binVector[%i] = %i",k,binVector[k]);
      } 
    }
    
    const size_t idx = index(binVector);
    if (data[idx].get()==0) {
      data[idx] = new CountStatsElement;
    }
    data[idx]->Ntot++;
    if (obs) data[idx]->Nobs++;
    if (dsc) data[idx]->Ndsc++;
    return true; 
  }
public:
  bool bin(std::vector<size_t> & binVector,
	   const std::vector<double> & xVector) const {
    binVector.resize(xVector.size());
    for (unsigned int k=0; k<xVector.size(); ++k) {
      binVector[k] = var[k]->bin(xVector[k]);
      if (binVector[k] > var[k]->size()) {
	// ORSA_DEBUG("bin larger than size, k: %i   bin: %i   size: %i",k,binVector[k],var[k]->size());
	return false;
      }	
    }
    return true;
  }
public:
  // from binVector to index
  size_t index(const std::vector<size_t> & binVector) const {
    size_t idx=0;
    size_t mul=1;
    for (unsigned int k=0; k<var.size(); ++k) {
      idx += binVector[k] * mul;
      mul *= var[k]->size();
    }
    // ORSA_DEBUG("index: %i",idx);
    return idx;
  }
public:
  // from index to binVector
  std::vector<size_t> bin(size_t index) const {
    size_t mul=data.size();
    std::vector<size_t> binVector;
    binVector.resize(var.size());
    {
      unsigned int k=var.size();
      do {
	--k;
	mul /= var[k]->size();
	binVector[k] = index/mul;
	index -= binVector[k]*mul;
      } while (k>0);
    }
    return binVector;
  }
public:
  // vector of the center of each bin
  std::vector<double> binCenterVector(const size_t index) const {
    std::vector<size_t> binVector = bin(index);
    std::vector<double> xVector;
    xVector.resize(var.size());
    for (unsigned int k=0; k<var.size(); ++k) {
      xVector[k] = var[k]->binCenter(binVector[k]);
    }
    return xVector;
  }
public:
  size_t size() const { return data.size(); }
public:
  const CountStatsElement * stats(const size_t index) const {
    return data[index].get();
  }
protected:
  const std::vector< osg::ref_ptr<Var> > & var;
protected:
  std::vector< osg::ref_ptr<CountStatsElement> > data;
};

#warning use this class?
class LinearHisto_ES {
public:
  LinearHisto_ES(const double startValue,
		 const double stopValue,
		 const double incrValue) :
    start(startValue),
    stop(stopValue),
    incr(incrValue) {
    histo.resize((size_t)(1+(stop-start)/incr));  
    for (unsigned int k=0; k< histo.size(); ++k) {
      histo[k] = new EfficiencyStatistics(start+0.5*incr*k);
    }
  }
public:
  bool insert(const double x,
	      const double val,
	      const double sigma) {
    const size_t histo_bin = bin(x);
    if (histo_bin>=histo.size()) {
      ORSA_DEBUG("problems...  bin: %i size: %i",histo_bin,histo.size());
      return false;
    }
    histo[histo_bin]->insert(val,sigma);
    return true;
  }
public:
  size_t bin(const double x) const {
    return (size_t)((x-start)/incr);
  }
public:
  const std::vector< osg::ref_ptr< EfficiencyStatistics > > & getHisto() const {
    return histo;
  }
public:
  const double start, stop, incr;
protected:
  std::vector< osg::ref_ptr< EfficiencyStatistics > > histo;
};

int main(int argc, char ** argv) {
  
  orsa::Debug::instance()->initTimer();
  
  if (argc != 2) {
    printf("Usage: %s <allEta-file>\n",argv[0]);
    exit(0);
  }
  
  ORSA_DEBUG("process ID: %i",getpid());
  
  const std::string basename = SkyCoverage::basename(argv[1]);
  
  std::vector<EfficiencyData> etaData;
  {
    FILE * fp = fopen(argv[1],"r");
    if (!fp) {
      ORSA_DEBUG("cannot open file [%s]",argv[1]);
      exit(0);
    }
    char line[1024];
    double H, V, apparentVelocity;
    double solarElongation, lunarElongation;
    double lunarPhase;
    double airMass;
    double galacticLatitude;
    char id[1024];
    int observed;
    int discovered;
    while (fgets(line,1024,fp)) {
      sscanf(line,"%lf %s %lf %lf %lf %lf %lf %lf %lf %i %i",
	     &H,
	     id,
	     &V,
	     &apparentVelocity,
	     &solarElongation,
	     &lunarElongation,
	     &lunarPhase,
	     &airMass,
	     &galacticLatitude,
	     &observed,
	     &discovered);
      // convert
      apparentVelocity = orsa::FromUnits(apparentVelocity*orsa::arcsecToRad(),orsa::Unit::HOUR,-1);
      solarElongation  = orsa::degToRad()*solarElongation;
      lunarElongation  = orsa::degToRad()*lunarElongation;
      lunarPhase       = orsa::degToRad()*lunarPhase;
      galacticLatitude = orsa::degToRad()*galacticLatitude;
      //
      EfficiencyData ed;
      ed.H=H;
      if (atoi(id) != 0) {
	ed.number = atoi(id);
      } else {
	ed.designation = id;
      }
      ed.V = V;
      ed.apparentVelocity = apparentVelocity;
      ed.solarElongation  = solarElongation;
      ed.lunarElongation  = lunarElongation;
      ed.lunarPhase       = lunarPhase;
      ed.airMass       = airMass;
      ed.galacticLatitude = galacticLatitude;
      ed.observed   = (observed==1);
      ed.discovered = (discovered==1);
      //
      etaData.push_back(ed);
    }
    fclose(fp);
  }
  
  // minimum apparent magnitude
  const double V0=18.0;
  
  // alternative method using CountStats
  std::vector< osg::ref_ptr<CountStats::Var> > varDefinition;
  //
  // apparent magnitude
  varDefinition.push_back(new CountStats::LinearVar(V0,
						    24,
						    0.2));
  // apparent velocity
  /* varDefinition.push_back(new CountStats::LogarithmicVar(orsa::FromUnits(  0.1*orsa::arcsecToRad(),orsa::Unit::HOUR,-1),
     orsa::FromUnits(300.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1),
     1.2));
  */
  //
  varDefinition.push_back(new CountStats::LinearVar(orsa::FromUnits( 0.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1),
						    orsa::FromUnits(50.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1),
						    orsa::FromUnits( 2.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1)));
  // solar elongation
  varDefinition.push_back(new CountStats::LinearVar(0.0,
						    orsa::pi(),
						    10.0*orsa::degToRad()));
  // lunar elongation
  varDefinition.push_back(new CountStats::LinearVar(0.0,
						    orsa::pi(),
						    10.0*orsa::degToRad()));
  // airmass
  /* varDefinition.push_back(new CountStats::LogarithmicVar(1.0,
     3.0,
     1.03));
  */
  //
  varDefinition.push_back(new CountStats::LinearVar(1.0,
						    2.0,
						    0.05));
  // galactic latitude
  varDefinition.push_back(new CountStats::LinearVar(-orsa::halfpi(),
						    orsa::halfpi(),
						    5.0*orsa::degToRad()));
  
  osg::ref_ptr<CountStats> countStats = 
    new CountStats(varDefinition);
  
  std::vector<double> xVector;
  xVector.resize(varDefinition.size());
  for (unsigned int k=0; k<etaData.size(); ++k) {
    // keep vars aligned with varDefinition content
    xVector[0] = etaData[k].V.getRef();
    xVector[1] = etaData[k].apparentVelocity.getRef();
    xVector[2] = etaData[k].solarElongation.getRef();
    xVector[3] = etaData[k].lunarElongation.getRef();
    xVector[4] = etaData[k].airMass.getRef();
    xVector[5] = etaData[k].galacticLatitude.getRef();
    const bool goodInsert = countStats->insert(xVector,
					       etaData[k].observed.getRef(),
					       etaData[k].discovered.getRef());
    if (!goodInsert) {
      // ORSA_DEBUG("problems with this one: V: %g",etaData[k].V.getRef())
    }
  }
  
  EfficiencyMultifit::DataStorage data;
  
  for (unsigned int index=0; index<countStats->size(); ++index) {
    const CountStatsElement * cs = countStats->stats(index);
    if (cs==0) continue;
    // write point only if Ntot != 0
    if (cs->Ntot>=2) {
      const double      eta = (double)(cs->Nobs)/(double)(cs->Ntot);
      const double sigmaEta = (double)(sqrt(cs->Nobs+1))/(double)(cs->Ntot); // Poisson counting statistics; using Nobs+1 instead of Nobs to have positive sigma even when Nobs=0
      
      const std::vector<double> xVector = countStats->binCenterVector(index); 
      
      EfficiencyMultifit::DataElement el;
      //
      el.V=xVector[0];
      el.U=xVector[1];		 
      el.SE=xVector[2];	 
      el.LE=xVector[3];
      // add lunar phase here
      el.AM=xVector[4];
      el.GL=xVector[5];
      //
      el.eta=eta;
      el.sigmaEta=sigmaEta;
      //
      el.Nobs=cs->Nobs;
      el.Ndsc=cs->Ndsc;
      el.Ntot=cs->Ntot;
      //
      data.push_back(el);
    }
  }
  
  // write eta.dat file
  if (1) {
    char filename[1024];
    sprintf(filename,"%s.eta.dat",basename.c_str());
    FILE * fp_eta = fopen(filename,"w"); 
    for (unsigned int k=0; k<data.size(); ++k) {
      const EfficiencyMultifit::DataElement & el = data[k];
      if (el.Ntot.getRef()>2) {
	fprintf(fp_eta,
		"%5.2f %6.2f %6.2f %6.2f %6.3f %5.3f %5.3f %5i %5i %5i\n",
		el.V.getRef(),
		orsa::FromUnits(el.U.getRef()*orsa::radToArcsec(),orsa::Unit::HOUR),
		el.SE.getRef()*orsa::radToDeg(),
		el.LE.getRef()*orsa::radToDeg(),
		el.AM.getRef(),
		el.eta.getRef(),
		el.sigmaEta.getRef(),
		el.Nobs.getRef(),
		el.Ndsc.getRef(),
		el.Ntot.getRef());
      }
    }
    fclose(fp_eta);
  }
  
  osg::ref_ptr<orsaInputOutput::MPCObsCodeFile> obsCodeFile = new orsaInputOutput::MPCObsCodeFile;
  obsCodeFile->setFileName("obscode.dat");
  obsCodeFile->read();
  
  char fitFilename[1024];
  sprintf(fitFilename,"%s.fit.dat",basename.c_str());
  
  unsigned int non_zero_n = 0;
  for (unsigned int k=0; k<data.size(); ++k) {
    if (data[k].eta.getRef()>0.0) {
      ++non_zero_n;
    }	
  }
  ORSA_DEBUG("non_zero_n: %i",non_zero_n);
  
  // multifit parameters  
  osg::ref_ptr<orsa::MultifitParameters> par = new orsa::MultifitParameters;
  //
  const double initial_U_limit = orsa::FromUnits( 40.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1); // arcsec/hour
  const double initial_w_U     = orsa::FromUnits( 10.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1); // arcsec/hour
  //    
  // V pars
  par->insert("eta0_V",  0.90, 0.0001);
  par->insert("c_V",     0.01, 0.00001);
  par->insert("V_limit",19.00, 0.001);
  par->insert("w_V",     0.50, 0.0001); 
  // U pars
  par->insert("U_limit",    initial_U_limit, 0.01*initial_U_limit);
  par->insert("w_U",        initial_w_U,     0.01*initial_w_U); 
  // mixing angle
  par->insert("beta",  0.0, 0.00001);
  // Galactic latitude
  par->insert("GL_limit",30.0*orsa::degToRad(),0.01*orsa::degToRad());
  // par->insert("c_GL",     0.000,               0.000001);
  par->insert("w_GL",     1.0*orsa::degToRad(),0.001*orsa::degToRad());
  // ranges
  par->setRange("eta0_V",0.0,1.0);
  // par->setRangeMin("c_GL",0.0);
  // fixed?
  // par->setFixed("U_limit",true);
  // par->setFixed("w_U",true);
  // par->setFixed("beta",true);
  // par->setFixed("c_GL",true);
  
  if (non_zero_n < par->sizeNotFixed()) {
    ORSA_DEBUG("not enought data for fit, exiting");
    // just create an empty output file
    FILE * fp = fopen(fitFilename,"w");
    fclose(fp);
    exit(0);
  }
  
  osg::ref_ptr<EfficiencyMultifit> etaFit= new EfficiencyMultifit;
  const bool success = etaFit->fit(par.get(),
				   data,
				   V0,
				   fitFilename,
				   basename,
				   obsCodeFile.get());
  if (!success) { 
    ORSA_DEBUG("fitting did not converge, exiting");
    // just "touch" the output file
    FILE * fp = fopen(fitFilename,"w");
    fclose(fp);
    exit(0);
  }
  osg::ref_ptr<const orsa::MultifitParameters> parFinal = etaFit->getMultifitParameters();
  // save final parameters
  const double  eta0_V = parFinal->get("eta0_V");
  const double     c_V = parFinal->get("c_V");
  const double V_limit = parFinal->get("V_limit");
  const double     w_V = parFinal->get("w_V");
  //
  const double U_limit = parFinal->get("U_limit");
  const double     w_U = parFinal->get("w_U");
  //
  const double    beta = parFinal->get("beta");
  //
  const double GL_limit = parFinal->get("GL_limit");
  // const double     c_GL = parFinal->get("c_GL"); 
  const double     w_GL = parFinal->get("w_GL");
  
  if (1) {
    // output for testing
    typedef std::list< osg::ref_ptr<EfficiencyStatistics> > AllStat;
    AllStat stat_V, stat_U; 
    AllStat stat_SE; // solarElongation
    AllStat stat_LE; // lunarElongation
    // AllStat stat_LP; // lunarPhase
    AllStat stat_AM; // airMass
    AllStat stat_GL; // galacticLatitude
    //
    for (unsigned int k=0; k<data.size(); ++k) {
      
      osg::ref_ptr<EfficiencyStatistics> sV;
      {
	AllStat::const_iterator it = stat_V.begin();
	while (it != stat_V.end()) {
	  if ((*it).get() != 0) {
	    if ((*it)->center == data[k].V.getRef()) {
	      sV = (*it).get();
	    }   
	  }
	  ++it;
	}
	if (sV.get() == 0) {
	  sV = new EfficiencyStatistics(data[k].V.getRef());
	  stat_V.push_back(sV);
	}
      }
      
      osg::ref_ptr<EfficiencyStatistics> sU;
      {
	AllStat::const_iterator it = stat_U.begin();
	while (it != stat_U.end()) {
	  if ((*it).get() != 0) {
	    if ((*it)->center == data[k].U.getRef()) {
	      sU = (*it).get();
	    }   
	  }
	  ++it;
	}
	if (sU.get() == 0) {
	  sU = new EfficiencyStatistics(data[k].U.getRef());
	  stat_U.push_back(sU);
	}
      }
      
      // solarElongation
      osg::ref_ptr<EfficiencyStatistics> sSE;
      {
	AllStat::const_iterator it = stat_SE.begin();
	while (it != stat_SE.end()) {
	  if ((*it).get() != 0) {
	    if ((*it)->center == data[k].SE.getRef()) {
	      sSE = (*it).get();
	    }   
	  }
	  ++it;
	}
	if (sSE.get() == 0) {
	  sSE = new EfficiencyStatistics(data[k].SE.getRef());
	  stat_SE.push_back(sSE);
	}
      }
      
      // lunarElongation
      osg::ref_ptr<EfficiencyStatistics> sLE;
      {
	AllStat::const_iterator it = stat_LE.begin();
	while (it != stat_LE.end()) {
	  if ((*it).get() != 0) {
	    if ((*it)->center == data[k].LE.getRef()) {
	      sLE = (*it).get();
	    }   
	  }
	  ++it;
	}
	if (sLE.get() == 0) {
	  sLE = new EfficiencyStatistics(data[k].LE.getRef());
	  stat_LE.push_back(sLE);
	}
      }
      
      // airMass
      osg::ref_ptr<EfficiencyStatistics> sAM;
      {
	AllStat::const_iterator it = stat_AM.begin();
	while (it != stat_AM.end()) {
	  if ((*it).get() != 0) {
	    if ((*it)->center == data[k].AM.getRef()) {
	      sAM = (*it).get();
	    }   
	  }
	  ++it;
	}
	if (sAM.get() == 0) {
	  sAM = new EfficiencyStatistics(data[k].AM.getRef());
	  stat_AM.push_back(sAM);
	}
      }
      
      // galacticLatitude
      osg::ref_ptr<EfficiencyStatistics> sGL;
      {
	AllStat::const_iterator it = stat_GL.begin();
	while (it != stat_GL.end()) {
	  if ((*it).get() != 0) {
	    if ((*it)->center == data[k].GL.getRef()) {
	      sGL = (*it).get();
	    }   
	  }
	  ++it;
	}
	if (sGL.get() == 0) {
	  sGL = new EfficiencyStatistics(data[k].GL.getRef());
	  stat_GL.push_back(sGL);
	}
      }
      
      const double eta = SkyCoverage::eta(data[k].V.getRef(),
					  V_limit,
					  eta0_V,
					  V0,
					  c_V,
					  w_V,
					  data[k].U.getRef(),
					  U_limit,
					  w_U,
					  beta,
					  data[k].GL.getRef(),
					  GL_limit,
					  // c_GL,
					  w_GL);
      
      const double nominal_eta_V = SkyCoverage::nominal_eta_V(data[k].V.getRef(),
							      V_limit,
							      eta0_V,
							      V0,
							      c_V,
							      w_V);
      
      const double nominal_eta_U = SkyCoverage::nominal_eta_U(data[k].U.getRef(),
							      U_limit,
							      w_U);
     
      const double nominal_eta_GL = SkyCoverage::nominal_eta_GL(data[k].GL.getRef(),
								GL_limit,
								// c_GL,
								w_GL);
      
      // at this point, small sigmas can become very large because divided by a very small eta's
      if (data[k].sigmaEta.getRef() > 0) {
	if (eta != 0) {
	  if (nominal_eta_V != 0) {
	    sV->insert(data[k].eta.getRef()*nominal_eta_V/eta,
		       orsa::square(eta/(nominal_eta_V*data[k].sigmaEta.getRef())));
	  }
	  if (nominal_eta_U != 0) {
	    sU->insert(data[k].eta.getRef()*nominal_eta_U/eta,
		       orsa::square(eta/(nominal_eta_U*data[k].sigmaEta.getRef())));
	  }
	  sSE->insert(data[k].eta.getRef()/eta,
		      orsa::square(eta/(data[k].sigmaEta.getRef())));
	  sLE->insert(data[k].eta.getRef()/eta,
		      orsa::square(eta/(data[k].sigmaEta.getRef())));
	  sAM->insert(data[k].eta.getRef()/eta,
		      orsa::square(eta/(data[k].sigmaEta.getRef())));
	  if (nominal_eta_GL != 0) {
	    sGL->insert(data[k].eta.getRef()*nominal_eta_GL/eta,
			orsa::square(eta/(nominal_eta_GL*data[k].sigmaEta.getRef())));
	  }
	}
	if ( !sV->fit.isSet())  sV->fit = nominal_eta_V;
	if ( !sU->fit.isSet())  sU->fit = nominal_eta_U;
	if (!sSE->fit.isSet()) sSE->fit = 1.0;
 	if (!sLE->fit.isSet()) sLE->fit = 1.0;
   	if (!sAM->fit.isSet()) sAM->fit = 1.0;
   	if (!sGL->fit.isSet()) sGL->fit = nominal_eta_GL;
      }
    }
    
    // sort lists
    stat_V.sort(EfficiencyStatistics_ptr_cmp());
    stat_U.sort(EfficiencyStatistics_ptr_cmp());
    stat_SE.sort(EfficiencyStatistics_ptr_cmp());
    stat_LE.sort(EfficiencyStatistics_ptr_cmp());
    stat_AM.sort(EfficiencyStatistics_ptr_cmp());
    stat_GL.sort(EfficiencyStatistics_ptr_cmp());
    
    {
      FILE * fpv = fopen("v.fit.dat","w");
      AllStat::const_iterator it = stat_V.begin();
      while (it != stat_V.end()) {
	if (!(*it)->fit.isSet()) {
	  ORSA_DEBUG("problems: fit value not set");
	  ++it;
	  continue;
	}
	if ((*it)->standardDeviation()==0.0) {
	  // ORSA_DEBUG("--SKIPPING--");
	  ++it;
	  continue;
	}
       	fprintf(fpv,"%g %g %g %g\n",
		(*it)->center,
		(*it)->fit.getRef(),
		(*it)->average(),
		(*it)->standardDeviation());
	++it;
      }
      fclose(fpv);
    }
    
    {
      FILE * fpu = fopen("u.fit.dat","w");
      AllStat::const_iterator it = stat_U.begin();
      while (it != stat_U.end()) {
	if (!(*it)->fit.isSet()) {
	  ORSA_DEBUG("problems: fit value not set");
	  ++it;
	  continue;
	}
	if ((*it)->standardDeviation()==0.0) {
	  // ORSA_DEBUG("--SKIPPING--");
	  ++it;
	  continue;
	}
      	fprintf(fpu,"%g %g %g %g\n",
		orsa::FromUnits((*it)->center*orsa::radToArcsec(),orsa::Unit::HOUR), // arcsec/hour
		(*it)->fit.getRef(),
		(*it)->average(),
		(*it)->standardDeviation());
	++it;
      }
      fclose(fpu);
    }
    
    {
      FILE * fpse = fopen("se.fit.dat","w");
      AllStat::const_iterator it = stat_SE.begin();
      while (it != stat_SE.end()) {
	if (!(*it)->fit.isSet()) {
	  ORSA_DEBUG("problems: fit value not set");
	  ++it;
	  continue;
	}
	if ((*it)->standardDeviation()==0.0) {
	  // ORSA_DEBUG("--SKIPPING--");
	  ++it;
	  continue;
	}
      	fprintf(fpse,"%g %g %g %g\n",
		orsa::radToDeg()*(*it)->center,
		(*it)->fit.getRef(),
		(*it)->average(),
		(*it)->standardDeviation());
	++it;
      }
      fclose(fpse);
    }
    
    {
      FILE * fple = fopen("le.fit.dat","w");
      AllStat::const_iterator it = stat_LE.begin();
      while (it != stat_LE.end()) {
	if (!(*it)->fit.isSet()) {
	  ORSA_DEBUG("problems: fit value not set");
	  ++it;
	  continue;
	}
	if ((*it)->standardDeviation()==0.0) {
	  // ORSA_DEBUG("--SKIPPING--");
	  ++it;
	  continue;
	}
      	fprintf(fple,"%g %g %g %g\n",
		orsa::radToDeg()*(*it)->center,
		(*it)->fit.getRef(),
		(*it)->average(),
		(*it)->standardDeviation());
	++it;
      }
      fclose(fple);
    }
    
    {
      FILE * fpam = fopen("am.fit.dat","w");
      AllStat::const_iterator it = stat_AM.begin();
      while (it != stat_AM.end()) {
	if (!(*it)->fit.isSet()) {
	  ORSA_DEBUG("problems: fit value not set");
	  ++it;
	  continue;
	}
	if ((*it)->standardDeviation()==0.0) {
	  // ORSA_DEBUG("--SKIPPING--");
	  ++it;
	  continue;
	}
      	fprintf(fpam,"%g %g %g %g\n",
		(*it)->center,
		(*it)->fit.getRef(),
		(*it)->average(),
		(*it)->standardDeviation());
	++it;
      }
      fclose(fpam);
    }
    
    {
      FILE * fpgl = fopen("gl.fit.dat","w");
      AllStat::const_iterator it = stat_GL.begin();
      while (it != stat_GL.end()) {
	if (!(*it)->fit.isSet()) {
	  ORSA_DEBUG("problems: fit value not set");
	  ++it;
	  continue;
	}
	if ((*it)->standardDeviation()==0.0) {
	  // ORSA_DEBUG("--SKIPPING--");
	  ++it;
	  continue;
	}
      	fprintf(fpgl,"%g %g %g %g\n",
		orsa::radToDeg()*(*it)->center,
		(*it)->fit.getRef(),
		(*it)->average(),
		(*it)->standardDeviation());
	++it;
      }
      fclose(fpgl);
    }
    
  }
  
  exit(0);
}
