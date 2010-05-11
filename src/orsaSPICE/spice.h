#ifndef _ORSA_SPICE_SPICE_
#define _ORSA_SPICE_SPICE_

#include <QMutex>

#include <string>

#include <orsa/cache.h>
#include <orsa/datetime.h>
#include <orsa/matrix.h>
#include <orsa/vector.h>

namespace orsaSPICE {
  
    class SPICE {
    public:
        static SPICE * instance() {
            if (_instance == 0) {
                _instance = new SPICE;
            }
            return _instance;
        }
    protected:
        SPICE() { 
            _observer.set("SSB");
            // _global.set("J2000");
            _global.set("ECLIPJ2000");  
        }
    public:
        virtual ~SPICE() {
            _instance = 0;
        }
    protected:
        static SPICE * _instance;
    
    public:    
        void loadKernel(const std::string & filename);
        void unloadKernel(const std::string & filename);
    
    public:
        void getPosVel(const std::string & target,
                       const orsa::Time  & ephemerisTime,
                       const std::string & observer,
                       orsa::Vector      & relativePosition,
                       orsa::Vector      & relativeVelocity);
    public:
        void getPosVel(const std::string & target,
                       const orsa::Time  & ephemerisTime,
                       orsa::Vector      & relativePosition,
                       orsa::Vector      & relativeVelocity);
    public:
        void setDefaultObserver(const std::string & observer) {
            _observer.set(observer);
        }
    public:
        const std::string & getDefaultObserver() const {
            return _observer.getRef();
        }
    protected:
        orsa::Cache<std::string> _observer;
    
    public:
        // seconds since J2000
        static double SPICETime(const orsa::Time & t);
    
    public:
        orsa::Matrix localToGlobal(const std::string & local,
                                   const orsa::Time  & ephemerisTime);
    public:
        orsa::Matrix globalToLocal(const std::string & local,
                                   const orsa::Time  & ephemerisTime);
    protected:
        orsa::Cache<std::string> _global;
    
        // public mutex calls, to sync with spice calls placed in code outside this class..
        // for the future, all the spice calls should be inside this class only.
    public:
        void lock() {
            mutex.lock();
        }
    public:
        void unlock() {
            mutex.unlock();
        }
    private:
        // a mutex that provides access serialization to the SPICE library calls
        QMutex mutex;
    };
  
}; // namespace orsaSPICE

#endif // _ORSA_SPICE_SPICE_
