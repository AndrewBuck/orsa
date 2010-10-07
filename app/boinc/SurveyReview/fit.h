#ifndef __FIT_H__
#define __FIT_H__

#include <map>
#include <orsa/statistic.h>
#include <orsa/unit.h>

// only when debugging...
#include <orsa/crash.h>

// CountStats::LinearVar ranges
const double start_V  =   12.00;
const double  stop_V  =   24.00;
const double  step_V  =    0.25;
const double start_U  = orsa::FromUnits(  0.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1);
const double  stop_U  = orsa::FromUnits(600.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1);
const double  step_U  = orsa::FromUnits(  1.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1);
const double start_SE =    0.0*orsa::degToRad();
const double  stop_SE =  180.0*orsa::degToRad();
const double  step_SE =   10.0*orsa::degToRad();
const double start_LE =    0.0*orsa::degToRad();
const double  stop_LE =  180.0*orsa::degToRad();
const double  step_LE =    5.0*orsa::degToRad();
const double start_AM =    1.00;
const double  stop_AM =    3.00;
const double  step_AM =    0.05;
const double start_GL = -180.0*orsa::degToRad();
const double  stop_GL =  180.0*orsa::degToRad();
const double  step_GL =   10.0*orsa::degToRad();
const double start_GB =  -90.0*orsa::degToRad();
const double  stop_GB =   90.0*orsa::degToRad();
const double  step_GB =    5.0*orsa::degToRad();
/* const double start_EL = -180.0*orsa::degToRad();
   const double  stop_EL =  180.0*orsa::degToRad();
   const double  step_EL =   10.0*orsa::degToRad();
   const double start_EB =  -90.0*orsa::degToRad();
   const double  stop_EB =   90.0*orsa::degToRad();
   const double  step_EB =    5.0*orsa::degToRad(); */
const double start_AZ =    0.0*orsa::degToRad();
const double  stop_AZ =  360.0*orsa::degToRad();
const double  step_AZ =   20.0*orsa::degToRad();
const double start_LA =  -90.0*orsa::degToRad();
const double  stop_LA =   90.0*orsa::degToRad();
const double  step_LA =    5.0*orsa::degToRad();
const double start_SA =  -90.0*orsa::degToRad();
const double  stop_SA =    0.0*orsa::degToRad();
const double  step_SA =    5.0*orsa::degToRad();
/* const double start_LP =    0.0*orsa::degToRad();
   const double  stop_LP =  180.0*orsa::degToRad();
   const double  step_LP =    5.0*orsa::degToRad(); */
const double start_LI =       0.0;
const double  stop_LI =       1.0;
const double  step_LI =       0.05;

class EfficiencyStatistics : public orsa::WeightedStatistic<double> {
public:
    EfficiencyStatistics(const double & Center) : 
        orsa::WeightedStatistic<double>(),
        center(Center) { }
public:
    const double center;
};

typedef EfficiencyStatistics EfficiencyStatistics1D;

class EfficiencyStatistics2D : public orsa::WeightedStatistic<double> {
public:
    EfficiencyStatistics2D(const double & CenterX,
                           const double & CenterY) : 
        orsa::WeightedStatistic<double>(),
        centerX(CenterX), centerY(CenterY) { }
public:
    const double centerX, centerY;
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

template <typename T> class BinStats : public osg::Referenced {  
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
            return (size_t)ceil((stop-start)/incr);
        }
        size_t bin(const double x) const {
            const size_t retVal = (size_t)((x-start)/incr);
            // ORSA_DEBUG("x: %g\tbin: %i\tstart: %g\tincr: %g",x,retVal,start,incr);
            // the correct check follows; think twice before changing it!
            if (retVal>=size()) return ((size_t)-1);
            return retVal;
        }
        double binStart(const size_t bin) const {
            return (start+bin*incr);
        }
        double binStop(const size_t bin) const {
            return (start+(bin+1)*incr);
        }
    public:
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
            return (size_t)ceil(log(stop/start)/log(factor));
        }
        size_t bin(const double x) const {
            const size_t retVal = (size_t)(log(x/start)/log(factor));
            // the correct check follows; think twice before changing it!
            if (retVal>=size()) return ((size_t)-1);
            return retVal;
        }
        double binStart(const size_t bin) const {
            return (start*orsa::int_pow(factor,bin));
        }
        double binStop(const size_t bin) const {
            return (start*orsa::int_pow(factor,bin+1));
        }
    public:
        const double start, stop, factor;
    };
    
public:
    BinStats(const std::vector< osg::ref_ptr<Var> > & varDefinition) :
        osg::Referenced(),
        var(varDefinition) {
        totalSize=1;
        for (unsigned int k=0; k<varDefinition.size(); ++k) {
            totalSize *= varDefinition[k]->size();
        }  
        // ORSA_DEBUG("totalSize: %Zi",totalSize.get_mpz_t());
    }
protected:
    ~BinStats() { }    
public:
    bool bin(std::vector<size_t> & binVector,
             const std::vector<double> & xVector) const {
        if (xVector.size() != var.size()) {
            ORSA_DEBUG("dimension mismatch");
            return false;
        }
        binVector.resize(var.size());
        for (unsigned int k=0; k<var.size(); ++k) {
            binVector[k] = var[k]->bin(xVector[k]);
            if (binVector[k] == ((size_t)-1)) {
                // out of boundaries
                return false;
            }	
        }
        return true;
    }
public:
// from binVector to index
    mpz_class index(const std::vector<size_t> & binVector) const {
        mpz_class idx=0;
        mpz_class mul=1;
        for (unsigned int k=0; k<var.size(); ++k) {
            idx += binVector[k] * mul;
            mul *= var[k]->size();
        }
        // ORSA_DEBUG("index: %i",idx);
        return idx;
    }
public:
// from index to binVector
    std::vector<size_t> bin(mpz_class index) const {
        mpz_class mul=totalSize;
        std::vector<size_t> binVector;
        binVector.resize(var.size());
        {
            unsigned int k=var.size();
            do {
                --k;
                mul /= var[k]->size();
                binVector[k] = mpz_class(index/mul).get_si();
                index -= binVector[k]*mul;
            } while (k>0);
        }
        return binVector;
    }
public:
// vector of the center of each bin
    std::vector<double> binCenterVector(const mpz_class index) const {
        std::vector<size_t> binVector = bin(index);
        std::vector<double> xVector;
        xVector.resize(var.size());
        for (unsigned int k=0; k<var.size(); ++k) {
            xVector[k] = var[k]->binCenter(binVector[k]);
        }
        return xVector;
    }
public:
    mpz_class size() const { return totalSize; }
public:
    const T * stats(const mpz_class index) {
        typename DataType::iterator it = data.find(index);
        if (it != data.end()) {
            return (*it).second.get();
        } else {
            return 0;
        }
    }
protected:
    const std::vector< osg::ref_ptr<Var> > & var;
public:
    typedef typename std::map< mpz_class, osg::ref_ptr<T> > DataType;
protected:
    DataType data;
public:
    const DataType & getData() const {
        return data;
    }
protected:
    mpz_class totalSize;
};

class CountStats : public BinStats<CountStatsElement> {
public:
    CountStats(const std::vector< osg::ref_ptr<Var> > & varDefinition) :
        BinStats<CountStatsElement>(varDefinition) { }
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
        const mpz_class idx = index(binVector);
        if (data[idx].get()==0) {
            // lazy allocation
            data[idx] = new CountStatsElement;
        }
        
        data[idx]->Ntot++;
        if (obs) data[idx]->Nobs++;
        if (dsc) data[idx]->Ndsc++;
        
        return true; 
    }
};

/* 
   template <class T> class Histo : public T {
   public:
   Histo(const T * orig) : T(*orig) {
   histo.resize(T::size());
   for (unsigned int k=0; k<histo.size(); ++k) {
   histo[k] = new EfficiencyStatistics(T::binCenter(k));
   }
   }
   public:
   bool insert(const double x,
   const double val,
   const double sigma) {
   const size_t histo_bin = T::bin(x);
   if (histo_bin == ((size_t)-1)) {
   // out of boundaries
   return false;
   } else {
   histo[histo_bin]->insert(val,sigma);
   return true;
   }
   }
   public:
   typedef std::vector< osg::ref_ptr< EfficiencyStatistics > > HistoDataType;
   protected:
   HistoDataType histo;
   public:
   const HistoDataType & getData() const { return histo; }
   };
*/

class Histo {
public:
    Histo(const CountStats::Var * Var) : var(Var) {
        histo.resize(var->size());
        for (unsigned int k=0; k<histo.size(); ++k) {
            histo[k] = new EfficiencyStatistics(var->binCenter(k));
        }
    }
public:
    bool insert(const double x,
                const double val,
                const double sigma) {
        const size_t histo_bin = var->bin(x);
        if (histo_bin == ((size_t)-1)) {
            // out of boundaries
            return false;
        } else {
            histo[histo_bin]->insert(val,sigma);
            return true;
        }
    }
public:
    typedef std::vector< osg::ref_ptr< EfficiencyStatistics > > HistoDataType;
protected:
    HistoDataType histo;
public:
    const HistoDataType & getData() const { return histo; }
protected:
    const CountStats::Var * var;  
};

typedef Histo Histo1D;

class Histo2D {
public:
    Histo2D(const CountStats::Var * VarX,
            const CountStats::Var * VarY) : varx(VarX), vary(VarY) {
        histo.resize(varx->size());
        for (unsigned int j=0; j<varx->size(); ++j) {
            histo[j].resize(vary->size());
            for (unsigned int k=0; k<vary->size(); ++k) {
                histo[j][k] = new EfficiencyStatistics2D(varx->binCenter(j),vary->binCenter(k));
            }
        }
    }
public:
    bool insert(const double x,
                const double y,
                const double val,
                const double sigma) {
        const size_t histo_binx = varx->bin(x);
        if (histo_binx == ((size_t)-1)) {
            // out of boundaries
            return false;
        }
        const size_t histo_biny = vary->bin(y);
        if (histo_biny == ((size_t)-1)) {
            // out of boundaries
            return false;
        }
        histo[histo_binx][histo_biny]->insert(val,sigma);
        return true;
    }
public:
    typedef std::vector< std::vector< osg::ref_ptr< EfficiencyStatistics2D > > > HistoDataType;
protected:
    HistoDataType histo;
public:
    const HistoDataType & getData() const { return histo; }
protected:
    const CountStats::Var * varx;
    const CountStats::Var * vary;
};

#endif // __FIT_H__
