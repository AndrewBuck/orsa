#ifndef __FIT_H__
#define __FIT_H__

#include <map>
#include <orsa/statistic.h>

// CountStats::LinearVar ranges
const double start_V  = 10.0;
const double  stop_V  = 24.0;
const double  step_V  =  0.1;
const double start_U  = orsa::FromUnits(  0.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1);
const double  stop_U  = orsa::FromUnits(100.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1);
const double  step_U  = orsa::FromUnits(  1.0*orsa::arcsecToRad(),orsa::Unit::HOUR,-1);
const double start_SE = 0.0*orsa::degToRad();
const double  stop_SE = orsa::pi();
const double  step_SE = 2.0*orsa::degToRad();
const double start_LE = 0.0*orsa::degToRad();
const double  stop_LE = orsa::pi();
const double  step_LE = 2.0*orsa::degToRad();
const double start_AM = 1.00;
const double  stop_AM = 3.00;
const double  step_AM = 0.02;
const double start_GL = -orsa::pi();
const double  stop_GL =  orsa::pi();
const double  step_GL = 2.0*orsa::degToRad();
const double start_GB = -orsa::halfpi();
const double  stop_GB =  orsa::halfpi();
const double  step_GB = 2.0*orsa::degToRad();
const double start_AZ = 0.0*orsa::degToRad();
const double  stop_AZ = orsa::twopi();
const double  step_AZ = 5.0*orsa::degToRad();
const double start_LA = -orsa::halfpi();
const double  stop_LA =  orsa::halfpi();
const double  step_LA = 2.0*orsa::degToRad();
const double start_SA = -orsa::halfpi();
const double  stop_SA =  orsa::halfpi();
const double  step_SA = 2.0*orsa::degToRad();

class EfficiencyStatistics : public orsa::WeightedStatistic<double> {
public:
    EfficiencyStatistics(const double & aux) : 
        orsa::WeightedStatistic<double>(),
        center(aux) { }
public:
    const double center;
    // orsa::Cache<double> fit;
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
            return (size_t)ceil((stop-start)/incr);
        }
        size_t bin(const double x) const {
            if (x<start) return -1;
            if (x>stop)  return -1;
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
            return (size_t)ceil(log(stop/start)/log(factor));
        }
        size_t bin(const double x) const {
            if (x<start) return -1;
            if (x>stop)  return -1;
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
        totalSize=1;
        for (unsigned int k=0; k<varDefinition.size(); ++k) {
            totalSize *= varDefinition[k]->size();
        }  
        ORSA_DEBUG("totalSize: %Zi",totalSize.get_mpz_t());
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
            const mpz_class idx = index(binVector);
            ORSA_DEBUG("index: %Zi",idx.get_mpz_t());
            bin(idx);
            ORSA_DEBUG("reconstructed binVector:");
            for (unsigned int k=0; k<binVector.size(); ++k) {
                ORSA_DEBUG("binVector[%i] = %i",k,binVector[k]);
            } 
        }
        
        const mpz_class idx = index(binVector);
        if (data[idx].get()==0) {
            // lazy allocation
            data[idx] = new CountStatsElement;
            
            /* { // testing only
               static mpz_class numAllocated=0;
               ++numAllocated;
               if (numAllocated%1000==0) ORSA_DEBUG("numAllocated: %6Zi",numAllocated.get_mpz_t());
               }
            */
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
            if (binVector[k] == (size_t)-1) {
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
    const CountStatsElement * stats(const mpz_class index) {
        DataType::iterator it = data.find(index);
        if (it != data.end()) {
            return (*it).second.get();
        } else {
            return 0;
        }
    }
protected:
    const std::vector< osg::ref_ptr<Var> > & var;
public:
    typedef std::map< mpz_class, osg::ref_ptr<CountStatsElement> > DataType;
protected:
    DataType data;
public:
    const DataType & getData() const {
        return data;
    }
protected:
    mpz_class totalSize;
};

template <class T> class Histo : public T {
public:
    // dx can be a factor if T is log
    /* Histo(const double x1,
       const double x2,
       const double dx) : T(x1,x2,dx) {
       histo.resize(T::size());
       for (unsigned int k=0; k< histo.size(); ++k) {
       histo[k] = new EfficiencyStatistics(T::binCenter(k));
       }
       }
    */
public:
    Histo(const T * orig) : T(*orig) {
        histo.resize(T::size());
        for (unsigned int k=0; k< histo.size(); ++k) {
            histo[k] = new EfficiencyStatistics(T::binCenter(k));
        }
    }
public:
    bool insert(const double x,
                const double val,
                const double sigma) {
        const size_t histo_bin = T::bin(x);
        if (histo_bin == (size_t)-1) {
            // out of boundaries
            return false;
        }
        histo[histo_bin]->insert(val,sigma);
        return true;
    }
public:
    typedef std::vector< osg::ref_ptr< EfficiencyStatistics > > HistoDataType;
protected:
    HistoDataType histo;
public:
    const HistoDataType & getData() const { return histo; }
};

#endif // __FIT_H__
