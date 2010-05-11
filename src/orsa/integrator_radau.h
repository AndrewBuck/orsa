#ifndef _ORSA_INTEGRATOR_RADAU_
#define _ORSA_INTEGRATOR_RADAU_

// #include <map>
#include <vector>

#include <QHash>

#include <orsa/integrator.h>

#include <orsa/interaction.h>
#include <orsa/quaternion.h>

namespace orsa {
  
    class IntegratorRadau : public Integrator {
    public:
        IntegratorRadau();
    protected:
        ~IntegratorRadau();
    public:
        bool step(orsa::BodyGroup  * bg,
                  const orsa::Time & start,
                  const orsa::Time & timestep,
                  orsa::Time       & next_timestep);
    
    public:
        bool canHandleVelocityDependantInteractions() const { return true; }
    
    public:	
        bool lastCallRejected() const {
            /* if (_lastCallRejected.get()) {
               ORSA_DEBUG("LAST CALL REJECTED");
               }
            */
            //
            return _lastCallRejected.get(); 
        }
    private:
        orsa::Cache<bool> _lastCallRejected;
    
    public:
        void reset();
    
    private:
        void _init();
        void _body_mass_or_number_changed(orsa::BodyGroup  * bg,
                                          const orsa::Time & t);
    
    public:
        orsa::Cache<double> _accuracy;
    
    private:
        // Gauss-Radau spacings for substeps within a sequence, for the 15th order 
        // integrator. The sum of the H values should be 3.733333333333333 
        // const vector<double> h;
        // vector<double> h;
        double h[8];
    
        // Constant coefficients used in series expansions for X and V
        //  XC: 1/2,  1/6,  1/12, 1/20, 1/30, 1/42, 1/56, 1/72
        //  VC: 1/2,  1/3,  1/4,  1/5,  1/6,  1/7,  1/8
        // vector<double> xc,vc;
        double xc[8],vc[7];
    
        // vector<double> r,c,d,s;
        // double r[28],c[21],d[21],s[9];
        //
        double r[28],c[21],d[21];
    
        // double g[7][3000],b[7][3000],e[7][3000];
        // std::vector< std::vector<double> > g,b,e;
        // std::vector< std::vector<double> > g,b,e;
        // std::vector< std::map<const orsa::Body *, double> > g,b,e;
        //
        // std::vector< std::map<const orsa::Body *, orsa::Vector> > g,b,e;
        std::vector< QHash<const orsa::Body *, orsa::Vector> > g,b,e;
        //
        // rotational
        /* 
           std::vector< QHash<const orsa::Body *, double> > gPhi, gTheta, gPsi;
           std::vector< QHash<const orsa::Body *, double> > bPhi, bTheta, bPsi;
           std::vector< QHash<const orsa::Body *, double> > ePhi, eTheta, ePsi;
        */
        //
        std::vector< QHash<const orsa::Body *, orsa::Quaternion> > gQ, bQ, eQ; 
    
        // unsigned int nv, niter;
    
        // unsigned int niter;
    
        // double x[3000],v[3000],a[3000],x1[3000],v1[3000],a1[3000];
        // std::vector<double> x,v,a,x1,v1,a1;
        // std::map<const orsa::Body *, orsa::Vector> x,v,a,x1,v1,a1;
        QHash<const orsa::Body *, orsa::Vector> x,v,a,x1,v1,a1;
    
        // rotational components
        /* 
           QHash<const orsa::Body *, double> phi,       theta,       psi;
           QHash<const orsa::Body *, double> phiDot,    thetaDot,    psiDot;
           QHash<const orsa::Body *, double> phiDotDot, thetaDotDot, psiDotDot;
           //
           QHash<const orsa::Body *, double> phi1,       theta1,       psi1;
           QHash<const orsa::Body *, double> phi1Dot,    theta1Dot,    psi1Dot;
           QHash<const orsa::Body *, double> phi1DotDot, theta1DotDot, psi1DotDot;
        */
        //
        QHash<const orsa::Body *, orsa::Quaternion> Q, QDot, QDotDot;
        //
        QHash<const orsa::Body *, orsa::Quaternion> Q1, Q1Dot, Q1DotDot;
    
        // std::vector<double> mass;
        // std::map<const orsa::Body *, double> mass;
        QHash<const orsa::Body *, double> mass;
    
        // std::vector<Vector> acc;
        // std::map<const orsa::Body *, orsa::Vector> acc;
        // QHash<const orsa::Body *, orsa::Vector> acc;
        // orsa::Interaction::VectorHash acc;
        orsa::Interaction::InteractionVector acc;
        //
        // torque
        // orsa::Interaction::VectorHash torque;
        orsa::Interaction::InteractionVector torque;
    
        unsigned int size;
    };
  
} // namespace orsa

#endif // _ORSA_INTEGRATOR_RADAU_
