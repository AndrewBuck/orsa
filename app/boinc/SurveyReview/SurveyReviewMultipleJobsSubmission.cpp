#include <stdlib.h>

#include <cstdio>
#include <iostream>
#include <string>

#include <gmpxx.h>
#include <cmath>
#include <unistd.h>

#include "SurveyReview.h"
#include "grain.h"

#include <gsl/gsl_rng.h>

#include <orsa/unit.h>

using namespace std;

class SkipData {
public:
    int z_a_min, z_a_max;
    int z_e_min, z_e_max;
    int z_i_min, z_i_max;    
};

int main (int argc, char ** argv) {
    
    orsa::Debug::instance()->initTimer();
    
    if ( (argc != 11) &&
         (argc != 12) ) {
        // only part of the parameters are arguments
        cout << "Usage: " << argv[0] << " baseName samplesPerBin a_AU_min a_AU_max e_min e_max i_DEG_min i_DEG_max H_min H_max [skip_file]" << endl;
        exit(0);
    }
    
    // use pid as seed for rng
    const int pid = getpid();
    // cout << "process ID: " << pid << endl;
    gsl_rng * rnd = gsl_rng_alloc(gsl_rng_gfsr4);
    gsl_rng_set(rnd,pid);
    
    const std::string baseName = argv[1];
    const unsigned int samplesPerBin = atoi(argv[2]);
    const int z_a_min = lrint(atof(argv[3])/grain_a_AU);
    const int z_a_max = lrint(atof(argv[4])/grain_a_AU);
    const int z_e_min = lrint(atof(argv[5])/grain_e);
    const int z_e_max = lrint(atof(argv[6])/grain_e);
    const int z_i_min = lrint(atof(argv[7])/grain_i_DEG);
    const int z_i_max = lrint(atof(argv[8])/grain_i_DEG);
    const int z_H_min = lrint(atof(argv[9])/grain_H);
    const int z_H_max = lrint(atof(argv[10])/grain_H);

    // parse skip_file
    std::list<SkipData> skipList;
    if (argc == 12) {
        FILE * fp = fopen(argv[11],"r");
        if (fp) {
            char line[1024];
            SkipData S;
            while (fgets(line,1024,fp)) {
                if (6 == sscanf(line,
                                "%i %i %i %i %i %i",
                                &S.z_a_min,
                                &S.z_a_max,
                                &S.z_e_min,
                                &S.z_e_max,
                                &S.z_i_min,
                                &S.z_i_max)) {
                    skipList.push_back(S);
                }
            }            
            fclose(fp);
        }
    }
    
    const int z_a_delta = lrint(0.05/grain_a_AU);
    const int z_e_delta = lrint(0.05/grain_e);
    const int z_i_delta = lrint(5.00/grain_i_DEG);
    const int z_H_delta = lrint(1.00/grain_H);
    
    // the other angles all go from 0 to 360 deg in steps of 30 deg, 12 iter factor each
    
    // other args, constant for now
    // const double JD = 2455198.0; // epoch of orbits
    const double JD = 2455650; // epoch of orbits

    // approx. number of jobs, does not take into account regions with no NEOs or regions to skip
    /* mpz_class numJobs = 1;
       numJobs *= (z_a_max-z_a_min)/z_a_delta;
       numJobs *= (z_e_max-z_e_min)/z_e_delta;
       numJobs *= (z_i_max-z_i_min)/z_i_delta;
    */
    //
    // first count the real number of jobs
    // keep this code in sync with the code below
    mpz_class numJobs=0;
    for (int z_a=z_a_min; z_a<z_a_max; z_a+=z_a_delta) {
        for (int z_e=z_e_min; z_e<z_e_max; z_e+=z_e_delta) {
            for (int z_i=z_i_min; z_i<z_i_max; z_i+=z_i_delta) {
                
                {
                    // quick check if NEO
                    // minimum perihelion: q = a_min*(1-e_max), with min and max of this specific interval
                    const double q_min = orsa::FromUnits(grain_a_AU*z_a*(1.0-grain_e*(z_e+z_e_delta)),orsa::Unit::AU);
                    if (q_min > OrbitID::NEO_max_q) {
                        // ORSA_DEBUG("skipping, no NEOs in this interval");
                        continue;
                    }
                }                                
                
                // check skip list
                {
                    bool matching=false;
                    std::list<SkipData>::const_iterator it = skipList.begin();
                    while (it != skipList.end()) {
                        if ( (z_a==(*it).z_a_min) && (z_a+z_a_delta==(*it).z_a_max) &&
                             (z_e==(*it).z_e_min) && (z_e+z_e_delta==(*it).z_e_max) &&
                             (z_i==(*it).z_i_min) && (z_i+z_i_delta==(*it).z_i_max) ) {
                            // ORSA_DEBUG("skipping, matching entry in skip list");
                            matching=true;
                            break;
                        }
                        ++it;
                    }
                    if (matching) continue;
                }
                
                ++numJobs;                    
            }
        }
    }
    
    // should check directly for bins with 0 NEOs
    
    if (numJobs < 1) {
        cout << "numJobs must be greater than 0" << endl;
        exit(0);
    }
    
    if (numJobs > 99999) {
        cout << "numJobs is too big, are you sure of what you're doing?" << endl;
        exit(0);
    }

    bool batch=false;
    {
        FILE * fp = fopen(".stage.lock","r");
        if (fp!=0) {
            fclose(fp);
            gmp_printf("Submitting %Zi job(s), assuming batch call...\n",numJobs.get_mpz_t());
            batch=true;
        }
    }
    
    std::string yesno;
    const std::string yes = "y";
    if (!batch) {
        
        
        gmp_printf("Submitting %Zi job(s), proceed? [y/n] ",numJobs.get_mpz_t());
        cin >> yesno;
        // cout << "yesno: " << yesno << endl;
    }
    
    if (batch || (yesno == yes)) {
        
        {
            int waitSec = 5;
            while (waitSec > 0) {
                printf("\rStarting in %i seconds...",waitSec);
                fflush(stdout);
                sleep(1);
                --waitSec;
            }
            cout << endl;
        }
        
        char cmd[1024];
        FILE * fp;
        for (int z_a=z_a_min; z_a<z_a_max; z_a+=z_a_delta) {
            for (int z_e=z_e_min; z_e<z_e_max; z_e+=z_e_delta) {
                for (int z_i=z_i_min; z_i<z_i_max; z_i+=z_i_delta) {
                    
                    {
                        // quick check if NEO
                        // minimum perihelion: q = a_min*(1-e_max), with min and max of this specific interval
                        const double q_min = orsa::FromUnits(grain_a_AU*z_a*(1.0-grain_e*(z_e+z_e_delta)),orsa::Unit::AU);
                        if (q_min > OrbitID::NEO_max_q) {
                            ORSA_DEBUG("skipping, no NEOs in this interval");
                            continue;
                        }
                    }                                

                    // check skip list
                    {
                        bool matching=false;
                        std::list<SkipData>::const_iterator it = skipList.begin();
                        while (it != skipList.end()) {
                            if ( (z_a==(*it).z_a_min) && (z_a+z_a_delta==(*it).z_a_max) &&
                                 (z_e==(*it).z_e_min) && (z_e+z_e_delta==(*it).z_e_max) &&
                                 (z_i==(*it).z_i_min) && (z_i+z_i_delta==(*it).z_i_max) ) {
                                ORSA_DEBUG("skipping, matching entry in skip list");
                                matching=true;
                                break;
                            }
                            ++it;
                        }
                        if (matching) continue;
                    }
                    
                    fp = fopen("grid.dat","w");
                    fprintf(fp,"%.2f %i   %i %i %i   %i %i %i   %i %i %i   %i %i %i   %i %i %i   %i %i %i   %i %i %i\n",
                            JD,
                            samplesPerBin,
                            z_a,z_a+z_a_delta,z_a_delta,
                            z_e,z_e+z_e_delta,z_e_delta,
                            z_i,z_i+z_i_delta,z_i_delta,
                            0,360000,30000,
                            0,360000,30000,
                            0,360000,30000,
                            z_H_min,z_H_max,z_H_delta); // all range in H
                    fclose(fp);
                    //
                    fp = fopen("randomSeed.dat","w");
                    fprintf(fp,"%li\n",gsl_rng_uniform_int(rnd,1000000000));
                    fclose(fp);
                    //
                    /* snprintf(cmd,1024,"./SurveyReviewWorkGenerator %s_%i_%i-%i_%i-%i_%i-%i_%i-%i",
                       baseName.c_str(),
                       samplesPerBin,
                       z_a,z_a+z_a_delta,
                       z_e,z_e+z_e_delta,
                       z_i,z_i+z_i_delta,
                       z_H_min,z_H_max);
                    */
                    //
                    snprintf(cmd,1024,"./SurveyReviewWorkGenerator %s",
                             baseName.c_str());
                    cout << "executing: [" << cmd << "]" << endl;
                    if (system(cmd)) {
                        ORSA_DEBUG("problems with system call: [%s]",cmd);
                        exit(1); 
                    }
                    usleep(10000);
                }
            }
        }
    }
    
    gsl_rng_free(rnd); 
    
    return 0;
}
