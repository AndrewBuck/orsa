#include "SurveyReview.h"
#include "grain.h"
#include "skycoverage.h"
#include "eta.h"

#include <orsa/debug.h>
#include <orsa/print.h>
#include <orsa/statistic.h>

#include <orsaSolarSystem/datetime.h>
#include <orsaSolarSystem/observatory.h>
#include <orsaSolarSystem/print.h>

#include <orsaUtil/observatory.h>

#include <orsaSPICE/spice.h>
#include <orsaSPICE/spiceBodyRotationalCallback.h>
#include <orsaSPICE/spiceBodyTranslationalCallback.h>

// CSPICE prototypes and definitions.      
#include <SpiceUsr.h>

// BOINC API
#include <boinc_api.h>
#include <filesys.h>

// SQLite3
#include "sqlite3.h"

int main() {
    
    orsa::Debug::instance()->initTimer();
    
    boinc_init();
    
    ORSA_DEBUG("process ID: %i",getpid());
    
    /* 
       ORSA_DEBUG("NEO_max_q: %g = %g [AU]",
       OrbitID::NEO_max_q,
       orsa::FromUnits(OrbitID::NEO_max_q,orsa::Unit::AU,-1));
    */
    
    std::string resolvedFileName;
  
    // needed to work with SQLite database
    sqlite3     * db;
    char        * zErr;
    int           rc;
    std::string   sql;
  
    {
        // open database
        boinc_resolve_filename_s("results.db",resolvedFileName);
        rc = sqlite3_open(resolvedFileName.c_str(),&db);
        //
        if (rc) {
            fprintf(stderr,"Can't open db: %s\n",sqlite3_errmsg(db));
            sqlite3_close(db);
            boinc_finish(0); 
        }
    }
  
    {
        // get list of tables, to see if some results have been obtained already
        // if no table is present, create it
        char **result;
        int nrows, ncols;
        sql = "SELECT name FROM sqlite_master";
        rc = sqlite3_get_table(db,sql.c_str(),&result,&nrows,&ncols,&zErr);
        //
        if (rc != SQLITE_OK) {
            if (zErr != NULL) {
                fprintf(stderr,"SQL error: %s\n",zErr);
                sqlite3_free(zErr);
                boinc_finish(0); 
            }
        }
        // ORSA_DEBUG("nrows: %i  ncols: %i",nrows, ncols);
        //
        /* for (int i=0; i<nrows; ++i) {
           for (int j=0; j<ncols; ++j) {
           // i=0 is the header
           const int index = (i+1)*ncols+j;
           ORSA_DEBUG("result[%i] = %s",index, result[index]);
           }
           }
        */
        //
        bool createTable=true;
        //
        if (nrows==2) {
            createTable=false;
        }
        //
        sqlite3_free_table(result);
    
        if (createTable) {
            // create results table
            sql = "CREATE TABLE grid(z_a_min INTEGER, z_a_max INTEGER, z_e_min INTEGER, z_e_max INTEGER, z_i_min INTEGER, z_i_max INTEGER, z_node_min INTEGER, z_node_max INTEGER, z_peri_min INTEGER, z_peri_max INTEGER, z_M_min INTEGER, z_M_max INTEGER, z_H INTEGER, N_NEO INTEGER, N_PHO INTEGER, NEO_in_field INTEGER, PHO_in_field INTEGER, eta_NEO REAL, sigma_eta_NEO REAL, eta_PHO REAL, sigma_eta_PHO REAL)";
            rc = sqlite3_exec(db,sql.c_str(),NULL,NULL,&zErr);
            //
            if (rc != SQLITE_OK) {
                if (zErr != NULL) {
                    fprintf(stderr,"SQL error: %s\n",zErr);
                    sqlite3_free(zErr);
                    boinc_finish(0); 
                }
            }
      
            // create random number state table
            sql = "CREATE TABLE rng(binary_state BLOB)";
            rc = sqlite3_exec(db,sql.c_str(),NULL,NULL,&zErr);
            //
            if (rc != SQLITE_OK) {
                if (zErr != NULL) {
                    fprintf(stderr,"SQL error: %s\n",zErr);
                    sqlite3_free(zErr);
                    boinc_finish(0); 
                }
            }
        }
    }
  
    osg::ref_ptr<orsaInputOutput::MPCObsCodeFile> obsCodeFile = new orsaInputOutput::MPCObsCodeFile;
    boinc_resolve_filename_s("obscode.dat",resolvedFileName);
    obsCodeFile->setFileName(resolvedFileName);
    if (!obsCodeFile->read()) {
        ORSA_DEBUG("cannot open file [%s]",resolvedFileName.c_str());
        boinc_finish(0); 
    }
    
    osg::ref_ptr<orsaUtil::StandardObservatoryPositionCallback> obsPosCB =
        new orsaUtil::StandardObservatoryPositionCallback(obsCodeFile.get());
  
    {
        // spice error file (should resolve filename?)
        SpiceChar spiceError[1024];
        sprintf(spiceError,"cspice.error");
        ::errdev_c("SET", 4096, spiceError);
    
        orsaSPICE::SPICE::instance()->setDefaultObserver("SSB");
    
        boinc_resolve_filename_s("de405.bsp",resolvedFileName);
		FILE * fp = boinc_fopen(resolvedFileName.c_str(),"r"); 
        if (!fp) {
            ORSA_DEBUG("cannot open file [%s]",resolvedFileName.c_str());
            boinc_finish(0); 
        }
        orsaSPICE::SPICE::instance()->loadKernel(resolvedFileName);
    }
    
    osg::ref_ptr<SkyCoverageFile> skyCoverageFile = new SkyCoverageFile;
    //
    osg::ref_ptr<SkyCoverage> skyCoverage = skyCoverageFile->_data;
    
    {
        // read fit.dat file
        
        FitFileDataElement e;
        boinc_resolve_filename_s("fit.dat",resolvedFileName);
        if (!readFitFile(e,resolvedFileName.c_str())) {
            ORSA_DEBUG("cannot open file [%s]",resolvedFileName.c_str());
            boinc_finish(0); 
        }
        
        std::string obsCode;
        orsa::Time epoch;
        int int_year;
        int int_dayOfYear;
        if (!SkyCoverage::processFilename(e.jobID,
                                          obsCodeFile.get(),
                                          obsCode,
                                          epoch,
                                          int_year,
                                          int_dayOfYear)) {
            ORSA_DEBUG("problems...");
            boinc_finish(0); 
        }
        
        // update skyCoverage
        skyCoverage->obscode  = obsCode;
        skyCoverage->epoch    = epoch;
        //
        skyCoverage->V_limit  = e.V_limit;
        skyCoverage->eta0_V   = e.eta0_V;
        skyCoverage->c_V      = e.c_V;
        skyCoverage->w_V      = e.w_V;
        //
        skyCoverage->U_limit  = e.U_limit;
        skyCoverage->w_U      = e.w_U;
        //
        skyCoverage->peak_AM  = e.peak_AM;
        skyCoverage->scale_AM = e.scale_AM;
        skyCoverage->shape_AM = e.shape_AM;
        //
        skyCoverage->drop_GB   = e.drop_GB;
        skyCoverage->scale_GB  = e.scale_GB;
        //
        skyCoverage->scale_GL = e.scale_GL;
        skyCoverage->shape_GL = e.shape_GL;
        //
        skyCoverage->V0       = e.V0;
    }
    
    boinc_resolve_filename_s("field.dat",resolvedFileName);
    skyCoverageFile->setFileName(resolvedFileName);
    if (!skyCoverageFile->read()) {
        ORSA_DEBUG("cannot open file [%s]",resolvedFileName.c_str());
        boinc_finish(0); 
    }
    //
    boinc_resolve_filename_s("fieldTime.dat",resolvedFileName);
    if (!skyCoverageFile->_data->readFieldTimeFile(resolvedFileName)) {
        ORSA_DEBUG("cannot open file [%s]",resolvedFileName.c_str());
        boinc_finish(0); 
    }
    
    // orsa::Vector sunPosition, earthPosition, moonPosition;
    // orsa::Vector sunVelocity, earthVelocity, moonVelocity;
    //
  
    osg::ref_ptr<orsa::BodyGroup> bg = new orsa::BodyGroup;
  
    // SUN
    osg::ref_ptr<orsa::Body> sun   = SPICEBody("SUN",orsaSolarSystem::Data::MSun());
    bg->addBody(sun.get());
  
    // EARTH
    osg::ref_ptr<orsa::Body> earth = SPICEBody("EARTH",orsaSolarSystem::Data::MEarth());
    bg->addBody(earth.get());
  
    // MOON
    osg::ref_ptr<orsa::Body> moon  = SPICEBody("MOON",orsaSolarSystem::Data::MMoon());
    bg->addBody(moon.get());
    
    int rs;
    //
    {
        boinc_resolve_filename_s("randomSeed.dat",resolvedFileName);
        FILE * fp = boinc_fopen(resolvedFileName.c_str(),"r"); 
        if (fp) { 
            gmp_fscanf(fp,"%i",&rs);
            fclose(fp);
        } else {
            ORSA_DEBUG("cannot open file [%s]",resolvedFileName.c_str());
            boinc_finish(0); 
        }
    }
    //
    const int randomSeed = rs;
    //
    ORSA_DEBUG("randomSeed: %i",randomSeed);
    //
    osg::ref_ptr<orsa::RNG> rnd = new orsa::RNG(randomSeed);
    //
    {
    
        // check if a state is available
        char **result;
        int nrows, ncols;
        sql = "select binary_state from rng";
        rc = sqlite3_get_table(db,sql.c_str(),&result,&nrows,&ncols,&zErr);
        //
        if (nrows==1) {
            ORSA_DEBUG("restoring rnd state...");
            // read last rng state
            std::string sqlselect = "select binary_state from rng";
            sqlite3_stmt * selectstmt;
            rc=sqlite3_prepare(db, sqlselect.c_str(), strlen(sqlselect.c_str()), &selectstmt, NULL);
            sqlite3_step(selectstmt);
            FILE * fp_rnd = fopen("rng.bin", "wb");
            fwrite (sqlite3_column_blob(selectstmt, 0), sqlite3_column_bytes(selectstmt, 0), 1, fp_rnd);
            fclose (fp_rnd);
            sqlite3_finalize(selectstmt);
            //
            fp_rnd = fopen("rng.bin", "rb");
            rnd->gsl_rng_fread(fp_rnd);
            fclose (fp_rnd);
        }
    }
  
    // osg::ref_ptr<OrbitFactory> orbitFactory;
    //
    double JD;
    unsigned int numGen;
    //
    int z_a_min, z_a_max, z_a_delta;
    int z_e_min, z_e_max, z_e_delta;
    int z_i_min, z_i_max, z_i_delta;
    int z_node_min, z_node_max, z_node_delta;
    int z_peri_min, z_peri_max, z_peri_delta;
    int z_M_min, z_M_max, z_M_delta;
    int z_H_min, z_H_max, z_H_delta;
    //
    {
        boinc_resolve_filename_s("grid.dat",resolvedFileName);
        FILE * fp = boinc_fopen(resolvedFileName.c_str(),"r"); 
        if (!fp) {
            ORSA_DEBUG("cannot open file [%s]",resolvedFileName.c_str());
            boinc_finish(0);  
        }
        bool done=false;
        char line[1024];
        while (fgets(line,1024,fp)) {
            if (strlen(line) > 0) {
                if (line[0] == '#') {
                    continue;
                }
            }
            if (23 == gmp_sscanf(line,
                                 "%lf %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i",
                                 &JD,
                                 &numGen,
                                 &z_a_min, &z_a_max, &z_a_delta,
                                 &z_e_min, &z_e_max, &z_e_delta,
                                 &z_i_min, &z_i_max, &z_i_delta,
                                 &z_node_min, &z_node_max, &z_node_delta,
                                 &z_peri_min, &z_peri_max, &z_peri_delta,
                                 &z_M_min, &z_M_max, &z_M_delta,
                                 &z_H_min, &z_H_max, &z_H_delta)) {
                done=true;
            }
        }
        fclose(fp);
        if (!done) {
            ORSA_DEBUG("cannot parse file");
            boinc_finish(0);   
        }
    }
    //
    const orsa::Time orbitEpoch = orsaSolarSystem::julianToTime(JD);
    // orsaSolarSystem::print(orbitEpoch);
    
    orsa::Orbit earthOrbit;
    // earthOrbit.compute(earth.get(),sun.get(),bg.get(),skyCoverage->epoch.getRef());
    earthOrbit.compute(earth.get(),sun.get(),bg.get(),orbitEpoch);
    
    const orsa::Time apparentMotion_dt_T = orsa::Time(0,0,1,0,0);
    //
    const double apparentMotion_dt = apparentMotion_dt_T.get_d();
    
    // osg::ref_ptr<OrbitID> orbit;
    
    orsa::Vector observerPosition_epoch;
    orsa::Vector observerPosition_epoch_plus_dt;
    orsa::Vector sunPosition_epoch;
    orsa::Vector sunPosition_epoch_plus_dt;
    
    obsPosCB->getPosition(observerPosition_epoch,
                          skyCoverage->obscode.getRef(),
                          skyCoverage->epoch.getRef());
    
    obsPosCB->getPosition(observerPosition_epoch_plus_dt,
                          skyCoverage->obscode.getRef(),
                          skyCoverage->epoch.getRef()+apparentMotion_dt_T);
    
    bg->getInterpolatedPosition(sunPosition_epoch,
                                sun.get(),
                                skyCoverage->epoch.getRef());
    
    bg->getInterpolatedPosition(sunPosition_epoch_plus_dt,
                                sun.get(),
                                skyCoverage->epoch.getRef()+apparentMotion_dt_T);
    
    // computed at skyCoverage->epoch
    const orsa::Vector observerPosition_sk_epoch         = observerPosition_epoch;
    const orsa::Vector observerPosition_sk_epoch_plus_dt = observerPosition_epoch_plus_dt;
    const orsa::Vector sunPosition_sk_epoch              = sunPosition_epoch;
    const orsa::Vector sunPosition_sk_epoch_plus_dt      = sunPosition_epoch_plus_dt;
    
    // double V_nightStart, V_nightStop;
    //
    double V_field;
    
    // estimate total number of iterations, to report progress
    mpz_class znum = 1;
    znum *= (z_a_max-z_a_min)/z_a_delta;
    znum *= (z_e_max-z_e_min)/z_e_delta;
    znum *= (z_i_max-z_i_min)/z_i_delta;
    znum *= (z_node_max-z_node_min)/z_node_delta;
    znum *= (z_peri_max-z_peri_min)/z_peri_delta;
    znum *= (z_M_max-z_M_min)/z_M_delta;
    // znum *= (z_H_max-z_H_min)/z_H_delta;
    const mpz_class totalIterations = znum;
    ORSA_DEBUG("maximum expected iterations: %Zi",totalIterations.get_mpz_t());
    
    // H not included in totalIterations; each iteration runs all H values
    const int size_H = 1 + (z_H_max-z_H_min)/z_H_delta;

    // defined outside the loop for better performance, but cleared at every loop iteration
    std::vector< osg::ref_ptr< orsa::Statistic<double> > > eta_NEO;
    std::vector< osg::ref_ptr< orsa::Statistic<double> > > eta_PHO;
    {
        eta_NEO.resize(size_H);
        eta_PHO.resize(size_H);
        for (int h=0; h<size_H; ++h) {
            eta_NEO[h] = new orsa::Statistic<double>;
            eta_PHO[h] = new orsa::Statistic<double>;
        }
    }

    {
        // begin first transaction, and wait for boinc_time_to_checkpoint() to commit it and start a new one
        sql = "begin";
        rc = sqlite3_exec(db,sql.c_str(),NULL,NULL,&zErr);
    }
    
    bool firstIter=true;
    for (int z_a=z_a_min; z_a<z_a_max; z_a+=z_a_delta) {
        for (int z_e=z_e_min; z_e<z_e_max; z_e+=z_e_delta) {
            for (int z_i=z_i_min; z_i<z_i_max; z_i+=z_i_delta) {
                for (int z_node=z_node_min; z_node<z_node_max; z_node+=z_node_delta) {
                    for (int z_peri=z_peri_min; z_peri<z_peri_max; z_peri+=z_peri_delta) {
                        for (int z_M=z_M_min; z_M<z_M_max; z_M+=z_M_delta) {

                            // for (int z_H=z_H_min; z_H<z_H_max; z_H+=z_H_delta) {
                            
                            {
                                // report progress
                                mpz_class idx=0;
                                mpz_class mul=totalIterations;
                                //
                                mul /= (z_a_max-z_a_min) / z_a_delta;
                                idx += mul*(z_a-z_a_min) / z_a_delta;
                                //
                                mul /= (z_e_max-z_e_min) / z_e_delta;
                                idx += mul*(z_e-z_e_min) / z_e_delta;
                                //
                                mul /= (z_i_max-z_i_min) / z_i_delta;
                                idx += mul*(z_i-z_i_min) / z_i_delta;
                                //
                                mul /= (z_node_max-z_node_min) / z_node_delta;
                                idx += mul*(z_node-z_node_min) / z_node_delta;
                                //
                                mul /= (z_peri_max-z_peri_min) / z_peri_delta;
                                idx += mul*(z_peri-z_peri_min) / z_peri_delta;
                                //
                                mul /= (z_M_max-z_M_min) / z_M_delta;
                                idx += mul*(z_M-z_M_min) / z_M_delta;
                                //
                                /* mul /= (z_H_max-z_H_min) / z_H_delta;
                                   idx += mul*(z_H-z_H_min) / z_H_delta;
                                */
                                //
                                const double fractionDone = idx.get_d()/totalIterations.get_d();
                                //
                                boinc_fraction_done(fractionDone);
                                if (boinc_is_standalone()) {
                                    ORSA_DEBUG("fraction completed: %5.2f%%",100*fractionDone);
                                    /* ORSA_DEBUG("fraction completed: %5.2f%%    idx: %Zi mul: %Zi totalIterations: %Zi",
                                       100*fractionDone,
                                       idx.get_mpz_t(),
                                       mul.get_mpz_t(),
                                       totalIterations.get_mpz_t());
                                    */
                                }
                            }
                            
                            if (firstIter) {
                                // if grid table is already populated,
                                // then this is a restore, and we skip directly to the last row
                                char **result;
                                int nrows, ncols;
                                char sql_line[1024];
                                // select last inserted row
                                sprintf(sql_line,"SELECT * FROM grid ORDER BY rowid desc LIMIT 1");
                                rc = sqlite3_get_table(db,sql_line,&result,&nrows,&ncols,&zErr);
                                //
                                if (rc != SQLITE_OK) {
                                    if (zErr != NULL) {
                                        fprintf(stderr,"SQL error: %s\n",zErr);
                                        sqlite3_free(zErr);
                                        boinc_finish(0); 
                                    }
                                }
                                // ORSA_DEBUG("nrows: %i  ncols: %i",nrows, ncols);
                                //
                                /* for (int i=0; i<nrows; ++i) {
                                   for (int j=0; j<ncols; ++j) {
                                   // i=0 is the header
                                   const int index = (i+1)*ncols+j;
                                   ORSA_DEBUG("result[%i] = %s = %s",j,result[j],result[index]);
                                   }
                                   }
                                */
                                //
                                if (nrows==1) {
                                    ORSA_DEBUG("skipping to last computed row");
                                    for (int j=0; j<ncols; ++j) {
                                        // first row is the header
                                        if (std::string(result[j])=="z_a_min")    { z_a    = atoi(result[ncols+j]); continue; }
                                        if (std::string(result[j])=="z_e_min")    { z_e    = atoi(result[ncols+j]); continue; }
                                        if (std::string(result[j])=="z_i_min")    { z_i    = atoi(result[ncols+j]); continue; }
                                        if (std::string(result[j])=="z_node_min") { z_node = atoi(result[ncols+j]); continue; }
                                        if (std::string(result[j])=="z_peri_min") { z_peri = atoi(result[ncols+j]); continue; }
                                        if (std::string(result[j])=="z_M_min")    { z_M    = atoi(result[ncols+j]); continue; }
                                        // if (std::string(result[j])=="z_H_min")    { z_H    = atoi(result[ncols+j]); continue; }
                                    }
                                }
                                sqlite3_free_table(result);
                                firstIter=false;
                            }
                            
                            {
                                // quick check if NEO
                                // minimum perihelion: q = a_min*(1-e_max), with min and max of this specific interval
                                const double q_min = orsa::FromUnits(grain_a_AU*z_a*(1.0-grain_e*(z_e+z_e_delta)),orsa::Unit::AU);
                                if (q_min > OrbitID::NEO_max_q) {
                                    // ORSA_DEBUG("skipping, no NEOs in this interval");
                                    continue;
                                }
                            }                                
                            
                            {
                                // check if already computed this one
                                char **result;
                                int nrows, ncols;
                                char sql_line[1024];
                                sprintf(sql_line,
                                        "SELECT * FROM grid WHERE z_a_min=%i and z_a_max=%i and z_e_min=%i and z_e_max=%i and z_i_min=%i and z_i_max=%i and z_node_min=%i and z_node_max=%i and z_peri_min=%i and z_peri_max=%i and z_M_min=%i and z_M_max=%i",
                                        z_a,z_a+z_a_delta,
                                        z_e,z_e+z_e_delta,
                                        z_i,z_i+z_i_delta,
                                        z_node,z_node+z_node_delta,
                                        z_peri,z_peri+z_peri_delta,
                                        z_M,z_M+z_M_delta);
                                // z_H,z_H+z_H_delta);
                                rc = sqlite3_get_table(db,sql_line,&result,&nrows,&ncols,&zErr);
                                //
                                if (rc != SQLITE_OK) {
                                    if (zErr != NULL) {
                                        fprintf(stderr,"SQL error: %s\n",zErr);
                                        sqlite3_free(zErr);
                                        boinc_finish(0); 
                                    }
                                }
                                // ORSA_DEBUG("nrows: %i  ncols: %i",nrows, ncols);
                                //
                                /* for (int i=0; i<nrows; ++i) {
                                   for (int j=0; j<ncols; ++j) {
                                   // i=0 is the header
                                   const int index = (i+1)*ncols+j;
                                   ORSA_DEBUG("result[%i] = %s",index, result[index]);
                                   }
                                   } 
                                */
                                //
                                bool skip=false;
                                //
                                if (nrows==0) {
                                    // nothing, but must keep this case!
                                } else if (nrows==size_H) {
                                    // ORSA_DEBUG("skipping value already computed...");
                                    /* ORSA_DEBUG("skipping: (%i,%i,%i,%i,%i,%i,%i)",
                                       z_a,
                                       z_e,
                                       z_i,
                                       z_node,
                                       z_peri,
                                       z_M,
                                       z_H);
                                    */
                                    skip=true;
                                } else { // if (nrows>size_H) {
                                    ORSA_ERROR("database corrupted, only (size_H=%i) entries per grid element are admitted",size_H);
                                    boinc_finish(0); 
                                }
                                //
                                sqlite3_free_table(result);
                                //
                                if (skip) {
                                    // ORSA_DEBUG("skipping...");
                                    continue;
                                }
                            }
                            
                            osg::ref_ptr<OrbitFactory> orbitFactory =
                                new OrbitFactory(grain_a_AU* z_a,
                                                 grain_a_AU*(z_a+z_a_delta),
                                                 grain_e* z_e,
                                                 grain_e*(z_e+z_e_delta),
                                                 grain_i_DEG* z_i,
                                                 grain_i_DEG*(z_i+z_i_delta),
                                                 grain_node_DEG* z_node,
                                                 grain_node_DEG*(z_node+z_node_delta),
                                                 grain_peri_DEG* z_peri,
                                                 grain_peri_DEG*(z_peri+z_peri_delta),
                                                 grain_M_DEG* z_M,
                                                 grain_M_DEG*(z_M+z_M_delta),
                                                 // grain_H* z_H,
                                                 // grain_H*(z_H+z_H_delta),
                                                 rnd.get(),
                                                 earthOrbit);
                            
                            unsigned int N_NEO=0;
                            unsigned int N_PHO=0;
                            unsigned int NEO_in_field=0;
                            unsigned int PHO_in_field=0;
                            // 
                            // now defined outside the loop
                            //
                            // std::vector< osg::ref_ptr< orsa::Statistic<double> > > eta_NEO;
                            // std::vector< osg::ref_ptr< orsa::Statistic<double> > > eta_PHO;
                            {
                                // eta_NEO.resize(size_H);
                                // eta_PHO.resize(size_H);
                                for (int h=0; h<size_H; ++h) {
                                    eta_NEO[h]->reset();
                                    eta_PHO[h]->reset();
                                }
                            }
                            
                            for (unsigned int j=0; j<numGen; ++j) {
                                
                                // if (j%10==0) ORSA_DEBUG("j: %i",j);
                                
                                osg::ref_ptr<OrbitID> orbit = orbitFactory->sample();
								
                                if (!orbit->isNEO()) {
                                    continue;
                                }
                                ++N_NEO;
                                const bool isPHO = orbit->isPHO();
                                if (isPHO) {
                                    ++N_PHO;
                                }
                                
                                // ++count;
                                
                                const double orbitPeriod = orbit->period();
                                
                                orsa::Vector r;
                                
                                const double original_M  = orbit->M;
                                //
                                orbit->M = original_M + fmod(orsa::twopi() * (skyCoverage->epoch.getRef()-orbitEpoch).get_d() / orbitPeriod, orsa::twopi());
                                orbit->relativePosition(r);
                                orsa::Vector orbitPosition_epoch = r + sunPosition_sk_epoch;
                                //
                                /* orbit->M = original_M + fmod(orsa::twopi() * (skyCoverage->epoch.getRef()+apparentMotion_dt_T-orbitEpoch).get_d() / orbitPeriod, orsa::twopi());
                                   orbit->relativePosition(r);
                                   orsa::Vector orbitPosition_epoch_plus_dt = r + sunPosition_sk_epoch_plus_dt;
                                */
                                //
                                orbit->M = original_M;
                                
                                orsa::Vector dr_epoch          = (orbitPosition_epoch         - observerPosition_sk_epoch);
                                // orsa::Vector dr_epoch_plus_dt  = (orbitPosition_epoch_plus_dt - observerPosition_sk_epoch_plus_dt);
                                
                                if (skyCoverage->get(dr_epoch.normalized(),V_field)) {
                                    
                                    // we already checked, it's a NEO
                                    ++NEO_in_field;
                                    if (isPHO) {
                                        ++PHO_in_field;
                                    }
                                    
                                    // retrieve one accurate observation epoch from the field, and recompute the orbit position
                                    orsa::Time epoch = skyCoverage->epoch.getRef();
                                    bool epochFromField = skyCoverage->pickFieldTime(epoch,dr_epoch.normalized(),rnd.get());
                                    
                                    obsPosCB->getPosition(observerPosition_epoch,
                                                          skyCoverage->obscode.getRef(),
                                                          epoch);
                                    
                                    obsPosCB->getPosition(observerPosition_epoch_plus_dt,
                                                          skyCoverage->obscode.getRef(),
                                                          epoch+apparentMotion_dt_T);
                                    
                                    bg->getInterpolatedPosition(sunPosition_epoch,
                                                                sun.get(),
                                                                epoch);
                                    
                                    bg->getInterpolatedPosition(sunPosition_epoch_plus_dt,
                                                                sun.get(),
                                                                epoch+apparentMotion_dt_T);
                                    
                                    orsa::Vector earthPosition_epoch;
                                    bg->getInterpolatedPosition(earthPosition_epoch,
                                                                earth.get(),
                                                                epoch);
                                    
                                    orsa::Vector moonPosition_epoch;
                                    bg->getInterpolatedPosition(moonPosition_epoch,
                                                                moon.get(),
                                                                epoch);
                                    
                                    // earth north pole
                                    const orsa::Vector northPole = (orsaSolarSystem::equatorialToEcliptic()*orsa::Vector(0,0,1)).normalized();
                                    
                                    orbit->M = original_M + fmod(orsa::twopi() * (epoch-orbitEpoch).get_d() / orbitPeriod, orsa::twopi());
                                    orbit->relativePosition(r);
                                    orbitPosition_epoch = r + sunPosition_epoch;
                                    //
                                    orbit->M = original_M + fmod(orsa::twopi() * (epoch+apparentMotion_dt_T-orbitEpoch).get_d() / orbitPeriod, orsa::twopi());
                                    orbit->relativePosition(r);
                                    orsa::Vector orbitPosition_epoch_plus_dt = r + sunPosition_epoch_plus_dt;
                                    //
                                    orbit->M = original_M;
                                    
                                    dr_epoch          = (orbitPosition_epoch         - observerPosition_epoch);
                                    orsa::Vector dr_epoch_plus_dt  = (orbitPosition_epoch_plus_dt - observerPosition_epoch_plus_dt);
                                    
                                    // ++inField;
                                    
                                    const orsa::Vector orb2obs    = observerPosition_epoch - orbitPosition_epoch;
                                    const orsa::Vector obs2orb    = -orb2obs;
                                    const orsa::Vector orb2sun    = sunPosition_epoch      - orbitPosition_epoch;
                                    const orsa::Vector obs2sun    = sunPosition_epoch      - observerPosition_epoch;
                                    const orsa::Vector obs2moon   = moonPosition_epoch     - observerPosition_epoch;
                                    const double       phaseAngle = acos((orb2obs.normalized())*(orb2sun.normalized()));
                                    
                                    const orsa::Vector moon2obs = observerPosition_epoch - moonPosition_epoch;
                                    const orsa::Vector moon2sun =      sunPosition_epoch - moonPosition_epoch;
                                    
                                    // apparent magnitude V moved later, in H loop
                                    
                                    // apparent velocity
                                    const double U = acos(dr_epoch_plus_dt.normalized()*dr_epoch.normalized())/apparentMotion_dt;
                                    
                                    // airmass
                                    //
                                    // aliases
                                    const orsa::Vector & earthPosition = earthPosition_epoch;
                                    const orsa::Vector &   obsPosition = observerPosition_epoch;
                                    const orsa::Vector     zenith = (obsPosition - earthPosition).normalized();
                                    // const orsa::Vector localEast  = orsa::externalProduct(northPole,zenith).normalized();
                                    //  const orsa::Vector localNorth = orsa::externalProduct(zenith,localEast).normalized();
                                    const double obs2orb_zenith     =     zenith*obs2orb.normalized();
                                    // const double obs2orb_localEast  =  localEast*obs2orb.normalized();
                                    // const double obs2orb_localNorth = localNorth*obs2orb.normalized();
                                    const double zenithAngle = acos(obs2orb_zenith);
                                    // const double airMass = ((observed||epochFromField)&&(zenithAngle<orsa::halfpi())?(1.0/cos(zenithAngle)):-1.0);
                                    // const double azimuth = fmod(orsa::twopi()+atan2(obs2orb_localEast,obs2orb_localNorth),orsa::twopi());
                                    const double AM = ((epochFromField)&&(zenithAngle<orsa::halfpi())?(1.0/cos(zenithAngle)):100.0);
                                    
                                    // const double solarAltitude = orsa::halfpi()-acos(zenith*obs2sun.normalized());
                                    // const double lunarAltitude = orsa::halfpi()-acos(zenith*obs2moon.normalized());
                                    // const double lunarElongation = acos(obs2moon.normalized()*obs2orb.normalized());
                                    // const double lunarPhase = acos(moon2obs.normalized()*moon2sun.normalized());
                                    //
                                    // const double SA = solarAltitude;
                                    // const double LA = lunarAltitude;
                                    // const double LE = lunarElongation;
                                    // const double LP = lunarPhase;
                                    
                                    // galactic latitude
                                    const orsa::Vector obs2orb_Equatorial = orsaSolarSystem::eclipticToEquatorial()*obs2orb;
                                    const orsa::Vector dr_equatorial = obs2orb_Equatorial.normalized();
                                    const double  ra = fmod(atan2(dr_equatorial.getY(),dr_equatorial.getX())+orsa::twopi(),orsa::twopi());
                                    const double dec = asin(dr_equatorial.getZ()/dr_equatorial.length());
                                    double l,b;
                                    orsaSolarSystem::equatorialToGalactic(l,b,ra,dec);
                                    // format longitude between -180 and +180 deg
                                    l = fmod(l+2*orsa::twopi(),orsa::twopi());
                                    if (l > orsa::pi()) l -= orsa::twopi();
                                    // const double galacticLongitude = l;
                                    // const double galacticLatitude  = b;
                                    const double GB = b;
                                    const double GL = l;
                                    // ecliptic coordinates
                                    /* const orsa::Vector dr = obs2orb.normalized();
                                       const double phi      = fmod(atan2(dr.getY(),dr.getX())+orsa::twopi(),orsa::twopi());
                                       const double theta    = asin(dr.getZ()/dr.length());
                                       const orsa::Vector dr_sun = obs2sun.normalized();
                                       const double phi_sun      = fmod(atan2(dr_sun.getY(),dr_sun.getX())+orsa::twopi(),orsa::twopi());
                                       const double theta_sun    = asin(dr_sun.getZ()/dr_sun.length());
                                       const double tmp_eclipticLongitude = fmod(phi-phi_sun+orsa::twopi(),orsa::twopi());
                                       const double eclipticLongitude = (tmp_eclipticLongitude>orsa::pi()) ? (tmp_eclipticLongitude-orsa::twopi()) : (tmp_eclipticLongitude);
                                       const double eclipticLatitude  = theta-theta_sun;
                                       const double EL = eclipticLongitude;
                                       const double EB = eclipticLatitude;
                                    */
                                    
                                    // total: size_H iterations
                                    for (int z_H=z_H_min; z_H<=z_H_max; z_H+=z_H_delta) {
                                        
                                        const double H = z_H*grain_H;
                                        
                                        const double G = 0.15;
                                        
                                        // apparent magnitude
                                        const double V = apparentMagnitude(H,
                                                                           G,
                                                                           phaseAngle,
                                                                           orb2obs.length(),
                                                                           orb2sun.length());

                                        // detection efficiency
                                        // const double eta = skyCoverage->eta(V,U,AM,GB,GL,EB,EL,SA,LA,LE);
                                        // const double eta = skyCoverage->eta(V,U,AM,GB,GL,SA,LA,LP2LI(LP));
                                        const double eta = skyCoverage->eta(V,U,AM,GB,GL);
                                        
                                        if (boinc_is_standalone()) {
                                            ORSA_DEBUG("a: %f [AU] e: %f i: %f [deg] H: %f V: %f eta: %e",
                                                       orsa::FromUnits(orbit->a,orsa::Unit::AU,-1),
                                                       orbit->e,
                                                       orsa::radToDeg()*orbit->i,
                                                       H,
                                                       V,
                                                       eta);
                                            
                                            // modified for best airmass
                                            /* const double eta = skyCoverage->eta(V,U,1.0,GB,GL);
                                               ORSA_DEBUG("RA: %6.3f  DEC: %+7.3f  eta: %.3f  V: %.2f  U: %6.2f  GB: %+6.2g  GL: %6.2f  orb2obs: %.2f  orb2sun: %.2f",
                                               orsa::radToDeg()*ra/15.0,
                                               orsa::radToDeg()*dec,
                                               eta,
                                               V,
                                               orsa::FromUnits(U*orsa::radToArcsec(),orsa::Unit::HOUR),
                                               GB*orsa::radToDeg(),
                                               GL*orsa::radToDeg(),
                                               orsa::FromUnits(orb2obs.length(),orsa::Unit::AU,-1),
                                               orsa::FromUnits(orb2sun.length(),orsa::Unit::AU,-1));
                                            */
                                        }
                                        
                                        const unsigned int pos = (z_H-z_H_min)/z_H_delta;
                                        
                                        eta_NEO[pos]->insert(eta);
                                        //
                                        if (isPHO) {
                                            eta_PHO[pos]->insert(eta);
                                        }
                                    }
                                    
                                } else {
                                    
                                    // no inserts if out-of-field
                                    // the average then represents only those in field
                                    
                                }
                            }
                            
                            {
                                // all in a transaction
                                /* MOVED OUT OF LOOP
                                   sql = "begin";
                                   rc = sqlite3_exec(db,sql.c_str(),NULL,NULL,&zErr);
                                */
                                
                                // total: size_H iterations
                                for (int z_H=z_H_min; z_H<=z_H_max; z_H+=z_H_delta) {
                                    
                                    const unsigned int pos = (z_H-z_H_min)/z_H_delta;
                                    
                                    // save values in db
                                    double       good_eta_NEO = (eta_NEO[pos]->entries() > 0) ? eta_NEO[pos]->average()           : 0;
                                    double good_sigma_eta_NEO = (eta_NEO[pos]->entries() > 1) ? eta_NEO[pos]->standardDeviation() : 0;
                                    double       good_eta_PHO = (eta_PHO[pos]->entries() > 0) ? eta_PHO[pos]->average()           : 0;
                                    double good_sigma_eta_PHO = (eta_PHO[pos]->entries() > 1) ? eta_PHO[pos]->standardDeviation() : 0;
                                    //
                                    // ORSA_DEBUG("eta_NEO: %f +/- %f",eta_NEO->average(),eta_NEO->averageError());
                                    //
#warning should use _bind_ sqlite commands for better floating point accuracy
                                    //
                                    char sql_line[1024];
                                    
                                    sprintf(sql_line,
                                            "INSERT INTO grid VALUES(%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%g,%g,%g,%g)",
                                            z_a,z_a+z_a_delta,
                                            z_e,z_e+z_e_delta,
                                            z_i,z_i+z_i_delta,
                                            z_node,z_node+z_node_delta,
                                            z_peri,z_peri+z_peri_delta,
                                            z_M,z_M+z_M_delta,
                                            z_H, // z_H,z_H+z_H_delta,
                                            N_NEO,
                                            N_PHO,
                                            NEO_in_field,
                                            PHO_in_field,
                                            good_eta_NEO,
                                            good_sigma_eta_NEO,
                                            good_eta_PHO,
                                            good_sigma_eta_PHO);
                                    // ORSA_DEBUG("executing: [%s]",sql_line);
                                    // ORSA_DEBUG("eta_NEO->entries(): %Zi",eta_NEO->entries().get_mpz_t());
                                    // ORSA_DEBUG("eta_PHO->entries(): %Zi",eta_PHO->entries().get_mpz_t());
                                    // do {
                                    rc = sqlite3_exec(db,sql_line,NULL,NULL,&zErr);
                                    // if (rc==SQLITE_BUSY) {
                                    // ORSA_DEBUG("database busy, retrying...");
                                    // }
                                    // } while (rc==SQLITE_BUSY);
                                    //
                                    if (rc != SQLITE_OK) {
                                        if (zErr != NULL) {
                                            fprintf(stderr,"SQL error: %s\n",zErr);
                                            sqlite3_free(zErr);
                                            boinc_finish(0); 
                                        }
                                    }
                                }
                                
                                // save rng state on database...
                                {
                                    // first delete old entry
                                    rc = sqlite3_exec(db,"delete from rng",NULL,NULL,&zErr);
                                    //
                                    if (rc != SQLITE_OK) {
                                        if (zErr != NULL) {
                                            fprintf(stderr,"SQL error: %s\n",zErr);
                                            sqlite3_free(zErr);
                                            boinc_finish(0); 
                                        }		
                                    }
                                    
                                    FILE * fp_rnd = fopen("rng.bin","wb");
                                    if (fp_rnd != 0) {
                                        rnd->gsl_rng_fwrite(fp_rnd);
                                        fclose(fp_rnd);
                                    } else {
                                        ORSA_ERROR("cannot open file rng.bin");  
                                        boinc_finish(1);
                                    }
                                    //
                                    fp_rnd = fopen("rng.bin","rb");
                                    // get the size f1Size of the input file
                                    fseek(fp_rnd, 0, SEEK_END);
                                    const unsigned int fileSize=ftell(fp_rnd);
                                    fseek(fp_rnd, 0, SEEK_SET);
                                    //
                                    char * copyBuffer = (char*)malloc(fileSize+1);
                                    //
                                    if (fileSize != fread(copyBuffer, sizeof(char), fileSize, fp_rnd)) {
                                        free (copyBuffer);
                                        boinc_finish(0); 
                                    }
                                    //
                                    fclose(fp_rnd);
                                    //
                                    sqlite3_stmt * insertstmt;
                                    std::string sqlinsert = "insert into rng (binary_state) values (?);";
                                    rc = sqlite3_prepare(db, sqlinsert.c_str(), strlen(sqlinsert.c_str()), &insertstmt, NULL);
                                    sqlite3_bind_blob(insertstmt, 1, (const void*)copyBuffer, fileSize, SQLITE_STATIC);
                                    sqlite3_step(insertstmt);
                                    sqlite3_finalize(insertstmt);
                                    free (copyBuffer);
                                }
                                
                                if (boinc_time_to_checkpoint()) {
                                    
                                    // commit open transaction, and then start a new one
                                    if (boinc_is_standalone()) {
                                        ORSA_DEBUG("checkpointing...");
                                    }
                                    
                                    sql = "commit";
                                    do {
                                        rc = sqlite3_exec(db,sql.c_str(),NULL,NULL,&zErr);
                                        if (rc==SQLITE_BUSY) {
                                            ORSA_DEBUG("database busy, retrying...");
                                        }
                                    } while (rc==SQLITE_BUSY);
                                    
                                    sql = "begin";
                                    rc = sqlite3_exec(db,sql.c_str(),NULL,NULL,&zErr); 
                                    
                                    boinc_checkpoint_completed();
                                }
                            }
                            
                            // six // seven "for" iterations...
                            
                            // }
                        }
                    }
                }
            }
        }
    }
    
    {
        // finally close last transaction
        sql = "commit";
        do {
            rc = sqlite3_exec(db,sql.c_str(),NULL,NULL,&zErr);
            if (rc==SQLITE_BUSY) {
                ORSA_DEBUG("database busy, retrying...");
            }
        } while (rc==SQLITE_BUSY);
    }
    
    // 100% done
    boinc_fraction_done(1.0);
    
    // cleanup
    do {
        rc = sqlite3_exec(db,"delete from rng",NULL,NULL,&zErr);
    } while (rc==SQLITE_BUSY);
  
    // vacuum
    do {
        rc = sqlite3_exec(db,"vacuum",NULL,NULL,&zErr);
    } while (rc==SQLITE_BUSY);
    
    // close db
    sqlite3_close(db);
    
    ORSA_DEBUG("done.");
    boinc_finish(0);
    
    return 0;
}

