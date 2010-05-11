#ifndef _ORSA_INTEGRATOR_LEAPFROG_
#define _ORSA_INTEGRATOR_LEAPFROG_

#include <orsa/integrator.h>

namespace orsa {
  
    class IntegratorLeapFrog : public Integrator {
    public:
        IntegratorLeapFrog();
    protected:
        ~IntegratorLeapFrog();
    public:
        bool step(orsa::BodyGroup  * bg,
                  const orsa::Time & start,
                  const orsa::Time & timestep,
                  orsa::Time       & next_timestep);	
    public:
        bool canHandleVelocityDependantInteractions() const { 
            ORSA_DEBUG("check this!");
            return false;
        }
    
    public:
        bool lastCallRejected() const { return false; }
    
    public:
        void reset() { }
    
    };
  
} // namespace orsa

#endif // _ORSA_INTEGRATOR_LEAPFROG_
