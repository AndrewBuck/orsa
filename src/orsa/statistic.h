#ifndef _ORSA_STATISTIC_
#define _ORSA_STATISTIC_

#include <vector>

#include <osg/Referenced>
#include <osg/ref_ptr>

#include <orsa/cache.h>

namespace orsa {
  
  // add min() and max() to the Statistic class...
  
  template <class T> class Statistic : public osg::Referenced {
  public:
    Statistic() : Referenced(true) {
      reset();
    }
  protected:
    ~Statistic() { }
    
  public:
    void reset() {
      _s = _s2 = 0;
      _n = 0;
    }
  public:
    void insert(const T & val) {
      _s  += val;
      _s2 += val*val;
      ++_n;
    }
  public:
    T sum() const {
      return _s;
    }
  public:
    T average() const {
      // return (_s/(T)_n);
      return (_s/_n.get_d());
    }
  public:
    T variance() const {
      if (_n > 1) {
	const T _a = average();
	// return (((T)_n/(T)(_n-1))*(_s2/(T)_n - _a*_a));
   	// fabs is needed to avoid sqrt of negative
	// values for functions calling variance(),
	// expecially when the inserted values are almost constant,
	// so that the difference (_s2/_n - _a*_a) can be negative
	// at the limit of the precision/roundoff errors...
	return fabs(mpz_class(_n/(_n-1)).get_d()*(_s2/_n.get_d() - _a*_a));
      } else {
	return 0;
      }
    }
  public:
    T standardDeviation() const {
      return (sqrt(variance()));
    }
  public:
    //! Error of the average, i.e. average() +/- averageError()
    T averageError() const {
      // return sqrt(variance()/(T)_n);
      return sqrt(variance()/_n.get_d());
    }
  public:
    const mpz_class & entries() const {
      return _n;
    }
  protected:
    T _s, _s2;
    mpz_class _n;
  };
  
  template <class T> class RunningStatistic : public osg::Referenced {
  public:
    RunningStatistic() : Referenced(true) {
      _v.resize(0);
      _loop_index = 0;
      _dirty = true;
      _stat = new Statistic<T>;
    }
  public:
    RunningStatistic(const size_t s) : Referenced(true) {
      _v.resize(s);
      _loop_index = 0;
      _dirty = true;
      _stat = new Statistic<T>;
    }
  protected:
    ~RunningStatistic() { }
    
  public:
    void setLength(const size_t s) {
      _v.resize(s);
      _loop_index = _loop_index%_v.size();
      _dirty = true;
    }
  public:
    const size_t getLength() const {
      return _v.size();
    } 
    
  public:
    void reset() {
      const size_t old_size = _v.size();
      _v.clear();
      _v.resize(old_size);
      _loop_index = 0;
      _dirty = true;
    }
  public:
    void insert(const T & val) {
      _v[_loop_index].set(val);
      _loop_index = ((_loop_index+1)%_v.size());
      _dirty = true;
    }
    
  public:
    const size_t size() const {
      _checkDirty();
      return (size_t)(_stat->entries().get_ui());
    }  
  public:
    const bool isFull() const {
      _checkDirty();
      return (size()==_v.size());
    } 
    
  public:
    T sum() const {
      _checkDirty();
      return _stat->sum();
    }
  public:
    T average() const {
      _checkDirty();
      return _stat->average();
    }
  public:
    T variance() const {
      _checkDirty();
      return _stat->variance();
    }
  public:
    T standardDeviation() const {
      _checkDirty();
      return _stat->standardDeviation();
    }
  public:
    //! Error of the average, i.e. average() +/- averageError()
    T averageError() const {
      _checkDirty();
      return _stat->averageError();
    }
    
  protected:
    void _checkDirty() const {
      if (_dirty) {
	_stat->reset();
	for (unsigned int k=0; k<_v.size(); ++k) {
	  if (_v[k].isSet()) {
	    _stat->insert(_v[k].getRef());
	  }
	}
      }
      _dirty = false;
    }
    
  protected:
    std::vector< orsa::Cache<T> > _v;
    size_t                        _loop_index;
  protected:
    mutable bool                          _dirty;
    mutable osg::ref_ptr< Statistic<T> >  _stat;
  };
  
} // namespace orsa

#endif // _ORSA_STATISTIC_
