#include <orsaSolarSystem/gauss.h>

#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>

#include <orsa/bodygroup.h>
#include <orsa/print.h>

#include <orsaSolarSystem/obleq.h>

struct poly_8_params {
    double coeff_6, coeff_3, coeff_0;
};

double poly_8     (double x, void *params);
double poly_8_df  (double x, void *params);
void   poly_8_fdf (double x, void *params, double *y, double *dy);

double poly_8 (double x, void *params) {
    struct poly_8_params *p = (struct poly_8_params *) params;
    return (pow(x,8) - p->coeff_6*pow(x,6) - p->coeff_3*pow(x,3) - p->coeff_0);
}

double poly_8_df (double x, void *params) {
    struct poly_8_params *p = (struct poly_8_params *) params;
    return (8*pow(x,7) - 6*p->coeff_6*pow(x,5) - 3*p->coeff_3*pow(x,2));
}

void poly_8_fdf (double x, void *params, double *y, double *dy) {
    struct poly_8_params *p = (struct poly_8_params *) params;
    *y  = (pow(x,8) - p->coeff_6*pow(x,6) - p->coeff_3*pow(x,3) - p->coeff_0);
    *dy = (8*pow(x,7) - 6*p->coeff_6*pow(x,5) - 3*p->coeff_3*pow(x,2));
  
    // debug
    /* 
       ORSA_DEBUG("poly_8_fdf:   x: %g   y: %g   dy: %g   coeff_6: %g   coeff_3: %g   coeff_0: %g",
       x,y,dy,p->coeff_6,p->coeff_6,p->coeff_6);
    */
}

class poly_8_solution {
public:
    double value, error;
};

void poly_8_gsl_solve(poly_8_params & params, 
                      std::vector<poly_8_solution> & solutions,
                      const double & minRange,
                      const double & maxRange,
                      const unsigned int steps) {
  
    // ORSA_DEBUG("inside poly_8...");
  
    // const double x_start  = FromUnits(100,orsa::Unit::KM).get_d();
    // const double x_start  = FromUnits(100,orsa::Unit::KM);
    // const double x_incr   = FromUnits(0.2,AU);
    // smaller increment needed for Earth's artificial satellites...
    /* 
       const double x_incr   = FromUnits(0.01,AU);
       const int    max_iter = 1500;
    */
    //
    // const double x_incr   = FromUnits(50,orsa::Unit::KM);
    // const int    max_iter = 1000;
  

    const double x_start = minRange;
    const double x_incr  = (maxRange-minRange)/steps;
    const int    max_iter = steps;

    const double nominal_relative_accuracy = 1.0e-5;
  
    solutions.clear();
  
    poly_8_solution tmp_solution;
  
    double x, x0; // this value is not used
    int    gsl_status;
    // const  gsl_root_fdfsolver_type *T;
  
    //   gsl_function_fdf FDF;
  
    /* 
       cerr << " poly_8_gsl_solve(): params.coeff_6 = " << params->coeff_6 << endl;
       cerr << " poly_8_gsl_solve(): params.coeff_3 = " << params->coeff_3 << endl;
       cerr << " poly_8_gsl_solve(): params.coeff_0 = " << params->coeff_0 << endl;
    */
  
    gsl_function_fdf FDF;
    //
    FDF.f      = &poly_8;
    FDF.df     = &poly_8_df;
    FDF.fdf    = &poly_8_fdf;
    FDF.params = &params;
  
    // T = gsl_root_fdfsolver_steffenson;
    // s = gsl_root_fdfsolver_alloc (T);
    gsl_root_fdfsolver * s = gsl_root_fdfsolver_alloc(gsl_root_fdfsolver_steffenson);
  
    int iter = -1;
    while (iter<max_iter) {
    
        ++iter;
        x = x_start+iter*x_incr;
        // std::cerr << "x: " << x << std::endl;
        gsl_root_fdfsolver_set (s, &FDF, x);
    
        int iter_gsl=0;
        const int max_iter_gsl = 10;
        do {
            ++iter_gsl;
            gsl_status = gsl_root_fdfsolver_iterate(s);
            x0 = x;
            x = gsl_root_fdfsolver_root(s);
            gsl_status = gsl_root_test_delta (x, x0, nominal_relative_accuracy, 0);
            // if (gsl_status == GSL_SUCCESS) printf ("Converged:\n");
            // printf ("%5d %10.7f %10.7f\n",iter_gsl, x, x - x0);
        } while ((gsl_status == GSL_CONTINUE) && (iter_gsl < max_iter_gsl));
    
        if (gsl_status == GSL_SUCCESS) {
            tmp_solution.value = x;
            // gsl_root_fdfsolver_iterate(s); tmp_solution.error = fabs(x-gsl_root_fdfsolver_root(s)); // GSL doc: ...the error can be estimated more accurately by taking the difference between the current iterate and next iterate rather than the previous iterate.
            tmp_solution.error = nominal_relative_accuracy;
      
            unsigned int k=0;
            bool duplicate=false;
            while (k<solutions.size()) {
                if (fabs(solutions[k].value-tmp_solution.value) < (solutions[k].error+tmp_solution.error)) {
                    duplicate= true;
                    break;
                }
                ++k;
            }
      
            if (!duplicate) {
                solutions.push_back(tmp_solution);
            }
        }
    }
  
    gsl_root_fdfsolver_free(s);
  
    /* 
       if (1) { 
       const double AU = orsa::FromUnits(1,orsa::Unit::AU).get_d();
       std::cerr << "solutions found: " << solutions.size() << std::endl;
       unsigned int k=0;
       while (k<solutions.size()) {
       std::cerr << k << ": "  << solutions[k].value << " accuracy: " << solutions[k].error << std::endl;
       //
       ORSA_DEBUG("solution[%i]: %f +/- %f [AU]",k,
       solutions[k].value/AU,
       solutions[k].error/AU);
       ++k;
       }
       }
    */
  
    // ORSA_DEBUG("leaving poly_8...");
}

void orsaSolarSystem::GaussMethod(std::vector<orsaSolarSystem::OrbitWithEpoch> & preliminaryOrbitVector,
                                  const orsaSolarSystem::ObservationVector & obs,
                                  const double & minRange,
                                  const double & maxRange,
                                  const unsigned int steps,
                                  const orsaSolarSystem::ObservatoryPositionCallback * obsPosCB,
                                  const orsa::Body * refBody,
                                  orsa::BodyGroup * bg) {
  
    if (obs.size() != 3) {
        ORSA_ERROR("exactly three observations needed, passed %i",obs.size());
        return;
    }
  
    orsaSolarSystem::OrbitWithEpoch orbit;
  
    // std::vector< osg::ref_ptr<orsaSolarSystem::Observation> > obs(obs_in);
  
    // test observations:
    // if same observatory, must have different epoch
    {
        for (unsigned int p=0; p<3; ++p) {
            for (unsigned int q=0; q<p; ++q) {
                if (obs[p]->obsCode.getRef() == obs[q]->obsCode.getRef()) {
                    if (obs[p]->epoch.getRef() == obs[q]->epoch.getRef()) {
                        ORSA_ERROR("observations not unique");
                        orsa::print(obs[p]->epoch.getRef());
                        orsa::print(obs[q]->epoch.getRef());
                        return;
                    }
                }
            }
        }	
    }
  
    // this does not do what we want... anyway, do we really need to sort?
    // sort(obs.begin(),obs.end());
  
    double tau[3];
    double sqrtGM[3];
    //
    {
        double mass;
        for (unsigned int k=0; k<3; ++k) {
            if (!bg->getInterpolatedMass(mass,refBody,obs[k]->epoch.getRef())) { ORSA_DEBUG("problems..."); }	
            sqrtGM[k] = sqrt(orsa::Unit::G()*mass);
            // ORSA_DEBUG("sqrtGM[%i]: %f",k,sqrtGM[k]());
        }
    
        tau[2] = 
            sqrtGM[1]*obs[1]->epoch.getRef().get_d() -
            sqrtGM[0]*obs[0]->epoch.getRef().get_d();
    
        tau[0] = 
            sqrtGM[2]*obs[2]->epoch.getRef().get_d() -
            sqrtGM[1]*obs[1]->epoch.getRef().get_d();
    
        tau[1] = 
            sqrtGM[2]*obs[2]->epoch.getRef().get_d() -
            sqrtGM[0]*obs[0]->epoch.getRef().get_d();
    }
  
    // debug
    /* 
       for (unsigned int k=0; k<3; ++k) {     
       // ORSA_DEBUG("epoch[%i]: %f",k,obs[k]->epoch.getRef().get_d());
       ORSA_DEBUG("tau[%i]: %f",k,tau[k]());
       }
    */
  
    orsa::Vector refBodyPosition[3];
    orsa::Vector refBodyVelocity[3];
    for (unsigned int k=0; k<3; ++k) {
        if (!bg->getInterpolatedPosVel(refBodyPosition[k],refBodyVelocity[k],refBody,obs[k]->epoch.getRef())) {
            ORSA_DEBUG("problems");
        }
    }	    
  
    orsa::Vector obsPosition[3];
    for (unsigned int k=0; k<3; ++k) {       
        if (!obsPosCB->getPosition(obsPosition[k],obs[k].get())) { ORSA_DEBUG("problems"); }
    }
  
    orsa::Vector R[3]; 
    for (unsigned int k=0; k<3; ++k) { 
        R[k] = obsPosition[k] - refBodyPosition[k];
    }
  
    orsa::Vector u_rho[3];
    {
        double c_ra,  s_ra;
        double c_dec, s_dec;
        for (unsigned int k=0; k<3; ++k) { 
            orsaSolarSystem::OpticalObservation * opticalObservation =  
                dynamic_cast<orsaSolarSystem::OpticalObservation *>(obs[k].get());
            if (!opticalObservation) {
                ORSA_DEBUG("observation is not optical");
            }
            orsa::sincos(opticalObservation->ra.getRef().getRad(), 
                         &s_ra, 
                         &c_ra);
            orsa::sincos(opticalObservation->dec.getRef().getRad(),
                         &s_dec,
                         &c_dec);
            u_rho[k].set(c_dec*c_ra,
                         c_dec*s_ra,
                         s_dec);
            u_rho[k] = orsaSolarSystem::equatorialToEcliptic()*u_rho[k];
        }
    }
  
    const orsa::Vector f = orsa::externalProduct(u_rho[0],u_rho[2]).normalized();
  
    const double rho_1_f = u_rho[1]*f;
  
    const double R_0_f = R[0]*f;
    const double R_1_f = R[1]*f;
    const double R_2_f = R[2]*f;
  
    const double A = (tau[0]/tau[1]*R_0_f + tau[2]/tau[1]*R_2_f - R_1_f)/rho_1_f;
    const double B = (tau[0]/tau[1]*(tau[1]*tau[1]-tau[0]*tau[0])*R_0_f + tau[2]/tau[1]*(tau[1]*tau[1]-tau[2]*tau[2])*R_2_f)/rho_1_f/6.0;
  
    const double Xl_Ym_Zn = R[1]*u_rho[1];
  
    poly_8_params params;
    params.coeff_6 = R[1].lengthSquared() + A*A + 2*A*Xl_Ym_Zn;
    params.coeff_3 = 2*A*B + 2*B*Xl_Ym_Zn;
    params.coeff_0 = B*B;
    std::vector<poly_8_solution> solutions;
    poly_8_gsl_solve(params,solutions,minRange,maxRange,steps);
  
    // ORSA_DEBUG("solutions: %i",solutions.size());
  
    if (solutions.size() > 0) {
    
        orsa::Vector rho[3];
        orsa::Vector r[3];
        orsa::Vector v; 
        double c[3];
    
        double tmp_length;
        double tmp_value;
    
        for (unsigned int p=0; p<solutions.size(); ++p) {
      
            // ORSA_DEBUG("solutions[%i] value: %f  error: %f",p,solutions[p].value,solutions[p].error);
      
            // rho[1] = u_rho[1]*(A+(B/secure_pow(solutions[p].value,3)));
            // check
            // tmp_length = A + (B/secure_pow(solutions[p].value,3));
            tmp_value = solutions[p].value;
            //
            if (tmp_value == 0.0) {
                // cerr << "out..." << endl;
                continue;
            }
            //
            tmp_length = A + (B/(tmp_value*tmp_value*tmp_value));
            // cerr << "tmp_length: " << tmp_length << endl;
      
            /* 
               ORSA_DEBUG("tmp_value: %Fg   tmp_length: %Fg   A: %f   B: %f",
               tmp_value(),
               tmp_length(),
               A(),
               B());
            */
      
            //
            if (tmp_length <= 0.0) {
                // cerr << "out..." << endl;
                continue;
            }
            //
            rho[1] = u_rho[1]*tmp_length;
      
            r[1] = R[1] + rho[1];
      
            // standard relation
            for (unsigned int j=0; j<3; ++j) { 
                // c[j] = tau[j]/tau[1]*(1+(secure_pow(tau[1],2)-secure_pow(tau[j],2))/(6*secure_pow(r[1].Length(),3)));
                c[j] = tau[j]/tau[1]*(1+(tau[1]*tau[1]-tau[j]*tau[j])/(6*orsa::int_pow(r[1].length(),3)));
                // printf("OLD c[%i] = %g\n",j,c[j]);
            }
      
            {
	
                const orsa::Vector v_k = rho[1] - (c[0]*R[0] + c[2]*R[2] - R[1]);
                const double   k = v_k.length();
                // const Vector u_k = v_k/v_k.Length();
                const orsa::Vector u_k = v_k.normalized();
	
                const double s02 = u_rho[0]*u_rho[2];
                const double s0k = u_rho[0]*u_k;
                const double s2k = u_rho[2]*u_k;
	
                // rho[0] = u_rho[0]*(k*(s0k-s02*s2k)/(1-secure_pow(s02,2)))/c[0];
                // tmp_length = (k*(s0k-s02*s2k)/(1-secure_pow(s02,2)))/c[0];
                tmp_length = (k*(s0k-s02*s2k)/(1-s02*s02))/c[0];
                if (tmp_length <= 0.0) {
                    // cerr << "out..." << endl;
                    continue;
                }
                //
                rho[0] = u_rho[0]*tmp_length;
	
                // rho[2] = u_rho[2]*(k*(s2k-s02*s0k)/(1-secure_pow(s02,2)))/c[2];
                // tmp_length = (k*(s2k-s02*s0k)/(1-secure_pow(s02,2)))/c[2];
                tmp_length = (k*(s2k-s02*s0k)/(1-s02*s02))/c[2];
                if (tmp_length <= 0.0) {
                    // cerr << "out..." << endl;
                    continue;
                }
                //
                rho[2] = u_rho[2]*tmp_length;
	
                r[0] = rho[0] + R[0];
                r[2] = rho[2] + R[2];
	
            }
      
            // ORSA_DEBUG("tic...");
      
      
            // try a simpler rule
            // v = (r[1]-r[0])/(FromUnits(obs[0].date.GetJulian()-obs[1].date.GetJulian(),DAY));
            /* v = ( (r[1]-r[0])/(FromUnits(obs[1].date.GetJulian()-obs[0].date.GetJulian(),DAY)) + 
               (r[2]-r[1])/(FromUnits(obs[2].date.GetJulian()-obs[1].date.GetJulian(),DAY)) ) / 2.0;
            */
            // v = velocity at the epoch of the first obs, obs[0]
            // Vector v = (r[1]-r[0])/(FromUnits(obs[0].date.GetJulian()-obs[1].date.GetJulian(),DAY));
            // orsa::Vector v = (r[1]-r[0])/(FromUnits(obs[1].date.GetJulian()-obs[0].date.GetJulian(),DAY));
            //
            orsa::Vector v = (r[1]-r[0]) / (obs[1]->epoch.getRef()-obs[0]->epoch.getRef()).get_d();
      
            // light-time correction [to be checked!]
            r[0] += (refBodyVelocity[0]+v)*(r[0]-R[0]).length()/orsa::Unit::c();
      
            /* 
               ORSA_DEBUG("r and v");
               orsa::print(r[0]);
               orsa::print(v);
               ORSA_DEBUG("r: %f [AU]   v: %f [km/s]",
               orsa::FromUnits(r[0].length(),orsa::Unit::AU,-1),
               orsa::FromUnits(orsa::FromUnits(v.length(),orsa::Unit::KM,-1),orsa::Unit::SECOND));
            */
      
            // orbit.ref_body = Body("Sun",GetMSun(),Vector(0,0,0),Vector(0,0,0));
            // orbit.Compute(r[1],v,GetG()*GetMSun());
            // orbit.compute(r[0],v,ref_jpl_planet,obs[0].date);
            //
            orbit.compute(r[0],v,sqrtGM[0]*sqrtGM[0]);
            orbit.epoch = obs[0]->epoch.getRef();
      
            // ORSA_DEBUG("orbit.mu: %f",orbit.mu);
      
            /*
              ORSA_DEBUG("tentative orbit: a=%f [au]   e=%f   i=%f [deg]",
              orsa::FromUnits(orbit.a,orsa::Unit::AU,-1),
              orbit.e(),
              orbit.i*orsa::radToDeg());
            */
      
            // #warning "check limit on dr..."
            // if ((orbit.e < 1.0) && 
            // ((r[0]-R[0]).length() > orsa::FromUnits(0.1,orsa::Unit::AU))) {
      
            // 1.1 to include comets...
            if (orbit.e < 1.1) {
	
                preliminaryOrbitVector.push_back(orbit);
	
                /* 
                   {
                   // test
                   orbit.computeRMS(obs,obsPosCB,refBody,bg);
                   //
                   ORSA_DEBUG("a: %f [AU]   e: %f   i: %f [deg]   rms: %f   ***** [TRIPLET ONLY]",
                   orsa::FromUnits(orbit.a,orsa::Unit::AU,-1),
                   orbit.e(),
                   orbit.i*orsa::radToDeg(),
                   orbit.rms.getRef());
                   }
                */
	
            }
        }
    }
}
