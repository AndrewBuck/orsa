#include "grain.h"

#include <orsa/debug.h>
#include <orsa/print.h>
#include <orsa/statistic.h>

#include <orsaSolarSystem/datetime.h>
#include <orsaSolarSystem/print.h>

// SQLite3
#include "sqlite3.h"

int main(int argc, char ** argv) {
    
    if (argc == 1) {
        ORSA_DEBUG("Usage: %s <sqlite-result-db(s)>",argv[0]);
        exit(0);
    }
    
    orsa::Debug::instance()->initTimer();
    
    ORSA_DEBUG("process ID: %i",getpid());
    
    // needed to work with SQLite database
    sqlite3     * db;
    char        * zErr;
    int           rc;
    // std::string   sql;

    for (int fileID=1; fileID<argc; ++fileID) {
        
        ORSA_DEBUG("inspecting db file [%s]",argv[fileID]);
        
        // open db
        // rc = sqlite3_open(argv[1],&db);
        rc = sqlite3_open_v2(argv[fileID],&db,SQLITE_OPEN_READONLY,NULL);
        if (rc) {
            ORSA_DEBUG("Can't open db: %s\n",sqlite3_errmsg(db));
            sqlite3_close(db);
            exit(0);
        }
        
        // find columns outside main loop
        int col_z_H=-1;
        int col_N_NEO=-1;
        int col_N_PHO=-1;
        int col_NEO_in_field=-1;
        int col_PHO_in_field=-1;
        int col_eta_NEO=-1;
        int col_sigma_eta_NEO=-1;
        int col_eta_PHO=-1;
        int col_sigma_eta_PHO=-1;
        // just to get the columns
        {
            char * * sql_result;
            int nrows, ncols;
            char sql_line[1024];
            sprintf(sql_line,"SELECT * FROM grid LIMIT 1");
            rc = sqlite3_get_table(db,sql_line,&sql_result,&nrows,&ncols,&zErr);
            if (rc != SQLITE_OK) {
                if (zErr != NULL) {
                    ORSA_DEBUG("SQL error: %s\n",zErr);
                    sqlite3_free(zErr);
                    exit(0); 
                }
            }
            
            if (ncols==0) {
                ORSA_DEBUG("skipping empty db");
                continue;
            }
            
            for (int col=0; col<ncols; ++col) {
                if (sql_result[col] == std::string("z_H")) {
                    col_z_H=col;
                } else if (sql_result[col] == std::string("N_NEO")) {
                    col_N_NEO=col;
                } else if (sql_result[col] == std::string("N_PHO")) {
                    col_N_PHO=col;
                } else if (sql_result[col] == std::string("NEO_in_field")) {
                    col_NEO_in_field=col;
                } else if (sql_result[col] == std::string("PHO_in_field")) {
                    col_PHO_in_field=col;
                } else if (sql_result[col] == std::string("eta_NEO")) {
                    col_eta_NEO=col;
                } else if (sql_result[col] == std::string("sigma_eta_NEO")) {
                    col_sigma_eta_NEO=col;
                } else if (sql_result[col] == std::string("eta_PHO")) {
                    col_eta_PHO=col;
                } else if (sql_result[col] == std::string("sigma_eta_PHO")) {
                    col_sigma_eta_PHO=col;
                }
            }
        
            sqlite3_free_table(sql_result);
        
            if ( (col_z_H==-1) ||
                 (col_N_NEO==-1) ||
                 (col_N_PHO==-1) ||
                 (col_NEO_in_field==-1) ||
                 (col_PHO_in_field==-1) ||
                 (col_eta_NEO==-1) ||
                 (col_sigma_eta_NEO==-1) ||
                 (col_eta_PHO==-1) ||
                 (col_sigma_eta_PHO==-1) ) {
                ORSA_DEBUG("problem: did not find some columns");
                exit(0); 
            }
        }


        {
            // simple version for now, assuming {a,e,i} range is constant in all entries (single bin)
        
            char **sql_result;
            int nrows, ncols;
            char sql_line[1024];
            sprintf(sql_line,"SELECT * FROM grid");
            rc = sqlite3_get_table(db,sql_line,&sql_result,&nrows,&ncols,&zErr);
            if (rc != SQLITE_OK) {
                if (zErr != NULL) {
                    ORSA_DEBUG("SQL error: %s\n",zErr);
                    sqlite3_free(zErr);
                    exit(0); 
                }
            }
        
            // read z_a, z_e, z_i from first row only (assume are constant in all rows)
#warning should read column names for z_* as well, to make sure code works even if column order changes
            int z_a_min = atoi(sql_result[ncols+0]);
            int z_a_max = atoi(sql_result[ncols+1]);
            int z_e_min = atoi(sql_result[ncols+2]);
            int z_e_max = atoi(sql_result[ncols+3]);
            int z_i_min = atoi(sql_result[ncols+4]);
            int z_i_max = atoi(sql_result[ncols+5]);

            std::vector<int> z_H_vec;
            //
            std::vector< osg::ref_ptr< orsa::WeightedStatistic<double> > > NEO_eta_field_ws_vec;
            std::vector< osg::ref_ptr< orsa::WeightedStatistic<double> > > NEO_eta_detect_ws_vec;
            //
            std::vector< osg::ref_ptr< orsa::WeightedStatistic<double> > > PHO_eta_field_ws_vec;
            std::vector< osg::ref_ptr< orsa::WeightedStatistic<double> > > PHO_eta_detect_ws_vec;
        
            for (int row=1; row<=nrows; ++row) {
                const int z_H              = atoi(sql_result[row*ncols+col_z_H]);
                const int N_NEO            = atoi(sql_result[row*ncols+col_N_NEO]);
                const int N_PHO            = atoi(sql_result[row*ncols+col_N_PHO]);
                const int NEO_in_field     = atoi(sql_result[row*ncols+col_NEO_in_field]);
                const int PHO_in_field     = atoi(sql_result[row*ncols+col_PHO_in_field]);
                const double eta_NEO       = atof(sql_result[row*ncols+col_eta_NEO]);
                const double sigma_eta_NEO = atof(sql_result[row*ncols+col_sigma_eta_NEO]);
                const double eta_PHO       = atof(sql_result[row*ncols+col_eta_PHO]);
                const double sigma_eta_PHO = atof(sql_result[row*ncols+col_sigma_eta_PHO]);
                //
                orsa::Cache<unsigned int> index;
                for (unsigned int k=0; k<z_H_vec.size(); ++k) {
                    if (z_H_vec[k]==z_H) {
                        index = k;
                        break;
                    }
                }
                if (!index.isSet()) {
                    // add one entry
                    z_H_vec.push_back(z_H);
                    index = z_H_vec.size()-1;
                    //
                    NEO_eta_field_ws_vec.resize(z_H_vec.size());
                    NEO_eta_field_ws_vec[index.get()] = new orsa::WeightedStatistic<double>;
                    //
                    NEO_eta_detect_ws_vec.resize(z_H_vec.size());
                    NEO_eta_detect_ws_vec[index.get()] = new orsa::WeightedStatistic<double>;
                    //
                    PHO_eta_field_ws_vec.resize(z_H_vec.size());
                    PHO_eta_field_ws_vec[index.get()] = new orsa::WeightedStatistic<double>;
                    //
                    PHO_eta_detect_ws_vec.resize(z_H_vec.size());
                    PHO_eta_detect_ws_vec[index.get()] = new orsa::WeightedStatistic<double>;
                }
                //
                if (N_NEO>0) {
                    const double NEO_eta_field       = (double)(NEO_in_field)/(double)(N_NEO);
                    const double p                   = (double)(NEO_in_field+1)/(double)(N_NEO+2);
                    const double NEO_sigma_eta_field = sqrt(p*(1-p)/N_NEO);
                    if (NEO_sigma_eta_field>0) NEO_eta_field_ws_vec[index.get()]->insert(NEO_eta_field,orsa::square(1.0/NEO_sigma_eta_field));
                    //
                    const double NEO_eta_detect       = eta_NEO;
                    const double NEO_sigma_eta_detect = sigma_eta_NEO;
                    if (NEO_sigma_eta_detect>0) NEO_eta_detect_ws_vec[index.get()]->insert(NEO_eta_detect,orsa::square(1.0/NEO_sigma_eta_detect));
                }
                //
                if (N_PHO>0) {
                    const double PHO_eta_field       = (double)(PHO_in_field)/(double)(N_PHO);
                    const double p                   = (double)(PHO_in_field+1)/(double)(N_PHO+2);
                    const double PHO_sigma_eta_field = sqrt(p*(1-p)/N_PHO);
                    if (PHO_sigma_eta_field>0) PHO_eta_field_ws_vec[index.get()]->insert(PHO_eta_field,orsa::square(1.0/PHO_sigma_eta_field));
                    //
                    const double PHO_eta_detect       = eta_PHO;
                    const double PHO_sigma_eta_detect = sigma_eta_PHO;
                    if (PHO_sigma_eta_detect>0) PHO_eta_detect_ws_vec[index.get()]->insert(PHO_eta_detect,orsa::square(1.0/PHO_sigma_eta_detect));
                }
            }

            if (0) {

                // verbose output
            
                ORSA_DEBUG("a: %g-%g [AU]   e: %g-%g   i: %g-%g [deg]",
                           z_a_min*grain_a_AU,
                           z_a_max*grain_a_AU,
                           z_e_min*grain_e,
                           z_e_max*grain_e,
                           z_i_min*grain_i_DEG,
                           z_i_max*grain_i_DEG);
            
                for (unsigned int index=0; index<z_H_vec.size(); ++index) {
                    if (NEO_eta_detect_ws_vec[index]->entries()==0) continue;
                    const double NEO_eta_obs =
                        NEO_eta_field_ws_vec[index]->average() *
                        NEO_eta_detect_ws_vec[index]->average();
                    const double NEO_sigma_eta_obs =
                        sqrt(orsa::square(NEO_eta_field_ws_vec[index]->average()*NEO_eta_detect_ws_vec[index]->standardDeviation()) +
                             orsa::square(NEO_eta_detect_ws_vec[index]->average()*NEO_eta_field_ws_vec[index]->standardDeviation()));
                    ORSA_DEBUG("NEO at H=%g -- in field: %.2e +/- %.2e (entries: %4Zi)   detect: %.2e +/- %.2e (entries: %4Zi)   obs: %.2e +/- %.2e",
                               z_H_vec[index]*grain_H,
                               NEO_eta_field_ws_vec[index]->average(),
                               NEO_eta_field_ws_vec[index]->standardDeviation(),
                               NEO_eta_field_ws_vec[index]->entries().get_mpz_t(),
                               NEO_eta_detect_ws_vec[index]->average(),
                               NEO_eta_detect_ws_vec[index]->standardDeviation(),
                               NEO_eta_detect_ws_vec[index]->entries().get_mpz_t(),
                               NEO_eta_obs,
                               NEO_sigma_eta_obs);
                }
            
                for (unsigned int index=0; index<z_H_vec.size(); ++index) {
                    if (PHO_eta_detect_ws_vec[index]->entries()==0) continue;
                    const double PHO_eta_obs =
                        PHO_eta_field_ws_vec[index]->average() *
                        PHO_eta_detect_ws_vec[index]->average();
                    const double PHO_sigma_eta_obs =
                        sqrt(orsa::square(PHO_eta_field_ws_vec[index]->average()*PHO_eta_detect_ws_vec[index]->standardDeviation()) +
                             orsa::square(PHO_eta_detect_ws_vec[index]->average()*PHO_eta_field_ws_vec[index]->standardDeviation()));
                    ORSA_DEBUG("PHO at H=%g -- in field: %.2e +/- %.2e (entries: %4Zi)   detect: %.2e +/- %.2e (entries: %4Zi)   obs: %.2e +/- %.2e",
                               z_H_vec[index]*grain_H,
                               PHO_eta_field_ws_vec[index]->average(),
                               PHO_eta_field_ws_vec[index]->standardDeviation(),
                               PHO_eta_field_ws_vec[index]->entries().get_mpz_t(),
                               PHO_eta_detect_ws_vec[index]->average(),
                               PHO_eta_detect_ws_vec[index]->standardDeviation(),
                               PHO_eta_detect_ws_vec[index]->entries().get_mpz_t(),
                               PHO_eta_obs,
                               PHO_sigma_eta_obs);
                }

            }

            // now, real output to file
            
            for (unsigned int index=0; index<z_H_vec.size(); ++index) {
                if (NEO_eta_detect_ws_vec[index]->entries()==0) continue;
                const double NEO_eta_obs =
                    NEO_eta_field_ws_vec[index]->average() *
                    NEO_eta_detect_ws_vec[index]->average();
                const double NEO_sigma_eta_obs =
                    sqrt(orsa::square(NEO_eta_field_ws_vec[index]->average()*NEO_eta_detect_ws_vec[index]->standardDeviation()) +
                         orsa::square(NEO_eta_detect_ws_vec[index]->average()*NEO_eta_field_ws_vec[index]->standardDeviation()));
                gmp_fprintf(stdout,
                            "NEO %5i %5i %4i %4i %6i %6i %3i %.2e %.2e %4Zi %.2e %.2e %4Zi %.2e %.2e\n",
                            z_a_min,
                            z_a_max,
                            z_e_min,
                            z_e_max,
                            z_i_min,
                            z_i_max,                           
                            z_H_vec[index],
                            NEO_eta_field_ws_vec[index]->average(),
                            NEO_eta_field_ws_vec[index]->standardDeviation(),
                            NEO_eta_field_ws_vec[index]->entries().get_mpz_t(),
                            NEO_eta_detect_ws_vec[index]->average(),
                            NEO_eta_detect_ws_vec[index]->standardDeviation(),
                            NEO_eta_detect_ws_vec[index]->entries().get_mpz_t(),
                            NEO_eta_obs,
                            NEO_sigma_eta_obs);
            }
            
            for (unsigned int index=0; index<z_H_vec.size(); ++index) {
                if (PHO_eta_detect_ws_vec[index]->entries()==0) continue;
                const double PHO_eta_obs =
                    PHO_eta_field_ws_vec[index]->average() *
                    PHO_eta_detect_ws_vec[index]->average();
                const double PHO_sigma_eta_obs =
                    sqrt(orsa::square(PHO_eta_field_ws_vec[index]->average()*PHO_eta_detect_ws_vec[index]->standardDeviation()) +
                         orsa::square(PHO_eta_detect_ws_vec[index]->average()*PHO_eta_field_ws_vec[index]->standardDeviation()));
                gmp_fprintf(stdout,
                            "PHO %5i %5i %4i %4i %6i %6i %3i %.2e %.2e %4Zi %.2e %.2e %4Zi %.2e %.2e\n",
                            z_a_min,
                            z_a_max,
                            z_e_min,
                            z_e_max,
                            z_i_min,
                            z_i_max,                           
                            z_H_vec[index],
                            PHO_eta_field_ws_vec[index]->average(),
                            PHO_eta_field_ws_vec[index]->standardDeviation(),
                            PHO_eta_field_ws_vec[index]->entries().get_mpz_t(),
                            PHO_eta_detect_ws_vec[index]->average(),
                            PHO_eta_detect_ws_vec[index]->standardDeviation(),
                            PHO_eta_detect_ws_vec[index]->entries().get_mpz_t(),
                            PHO_eta_obs,
                            PHO_sigma_eta_obs);
            }
            
            sqlite3_free_table(sql_result);
        }    
    
        // close db
        sqlite3_close(db);
    
    } // loop on files
    
    return 0;
}

