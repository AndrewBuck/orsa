#ifndef _ORSA_INTEGRATOR_
#define _ORSA_INTEGRATOR_

#include <orsaTBB/malloc.h>

#include <osg/Referenced>
#include <osg/ref_ptr>

// #include <orsa/bodygroup.h>
#include <orsa/datetime.h>
#include <orsa/interval.h>
#include <orsa/vector.h>

// #include <map>

namespace orsa {
  
    class Body;
    class BodyGroup;  
  
    class Integrator : public osg::Referenced {
        /* 
           public:
           enum IntegratorType {
           ITG_RADAU = 1,
           ITG_RK4 = 2,  
           ITG_LEAPFROG=3,
           //
           ITG_DEFAULT=ITG_LEAPFROG
           };
        */
    
    public:
        Integrator();
        /* 
           public:
           Integrator(const IntegratorType);
        */
    protected:
        virtual ~Integrator() { };
    
        /* 
           public:
           IntegratorType getType() const { return _type; };
           protected:
           const IntegratorType _type;
        */
    
    public:
        virtual bool canHandleVelocityDependantInteractions() const = 0;
    
    public:
        //! used by variable timestep integrators
        virtual bool lastCallRejected() const = 0;
    
    public:
        // reset all the internal variables
        virtual void reset() = 0;
    
        /* 
         * We should implement more than one integration approach, i.e:
         * - integrate over a period, and keep everything in memory (or on a file),
         *   so that you can refine the integration inside the period (splines & Co.);
         * - integrate over a period, given a timestep, and save all the timesteps 
         *   on file without saving anything in memory;
         * - integrate to a given time, saving the final condition only;
         * - ...
         *
         */
    
    public:
        //! Integrate all the bodies in the BodyGroup over the given period
        virtual bool integrate(orsa::BodyGroup  * bg,
                               const orsa::Time & start,
                               const orsa::Time & stop,
                               const orsa::Time & samplig_period);
    
    protected:
        /* 
           virtual bool step(orsa::BodyGroup  * bg,
           const orsa::Time & start,
           const orsa::Time & timestep) const = 0;	
        */
    
    protected:    
        virtual bool step(orsa::BodyGroup  * bg,
                          const orsa::Time & start,
                          const orsa::Time & timestep,
                          orsa::Time       & next_timestep) = 0;	
    
    protected:
        virtual void singleStepDone(orsa::BodyGroup  *,
                                    const orsa::Time &,
                                    const orsa::Time &,
                                    orsa::Time       &) const { }
    
    public:
        virtual void abort() const {
            doAbort = true;
        }
    
    private:
        mutable bool doAbort;
    
    public:
        orsa::Cache<unsigned int> progressiveCleaningSteps;
        orsa::Cache<bool>         keepOnlyLastStep;
    };
  
} // namespace orsa

#endif // _ORSA_INTEGRATOR_
