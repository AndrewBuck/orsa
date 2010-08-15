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
            continue;
        }
        
        {
            // integrity_check
            
            int nrows, ncols;
            char * * sql_result;
            rc = sqlite3_get_table(db,"pragma integrity_check",&sql_result,&nrows,&ncols,&zErr);
            //
            if (rc != SQLITE_OK) {
                if (zErr != NULL) {
                    ORSA_DEBUG("SQL error: %s\n",zErr); 
                    sqlite3_free(zErr);
                    sqlite3_close(db);
                    continue;
                }
            }
            //
            bool integrity_check_passed=false;
            if (nrows==1) {
                if (ncols==1) {
                    if (sql_result[1]==std::string("ok")) {
                        integrity_check_passed=true;
                    }
                }
            }
            if (!integrity_check_passed) {
                ORSA_DEBUG("SQLite problem: integrity_check failed\n"); 
                sqlite3_close(db);
                continue; 
            }    
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
                    sqlite3_close(db);
                    continue; 
                }
            }
            
            if (ncols==0) {
                ORSA_DEBUG("skipping empty db");
                sqlite3_close(db);
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
                sqlite3_close(db);
                continue;
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
                    sqlite3_close(db);
                    continue;
                }
            }
            
            // read z_a, z_e, z_i from first row only (assume are constant in all rows)
            // if NOT constant, move inside loop and change "nrows" to "row*nrows"
#warning should read column names for z_* as well, to make sure code works even if column order changes
            const int z_a_min    = atoi(sql_result[ncols+0]);
            const int z_a_max    = atoi(sql_result[ncols+1]);
            const int z_e_min    = atoi(sql_result[ncols+2]);
            const int z_e_max    = atoi(sql_result[ncols+3]);
            const int z_i_min    = atoi(sql_result[ncols+4]);
            const int z_i_max    = atoi(sql_result[ncols+5]);
            
            // extend by one more dimension: L, taking 12 values (0 to 360 step 30 deg)
            const unsigned int z_L_size = 12;
            
            std::vector<int> z_H_vec;
            //
            std::vector< std::vector< osg::ref_ptr< orsa::WeightedStatistic<double> > > > NEO_eta_field_ws_vec;
            std::vector< std::vector< osg::ref_ptr< orsa::WeightedStatistic<double> > > > NEO_eta_detect_ws_vec;
            //
            std::vector< std::vector< osg::ref_ptr< orsa::WeightedStatistic<double> > > > PHO_eta_field_ws_vec;
            std::vector< std::vector< osg::ref_ptr< orsa::WeightedStatistic<double> > > > PHO_eta_detect_ws_vec;
            
            for (int row=1; row<=nrows; ++row) {
                const int z_node_min = atoi(sql_result[row*ncols+6]);
                const int z_node_max = atoi(sql_result[row*ncols+7]);
                const int z_peri_min = atoi(sql_result[row*ncols+8]);
                const int z_peri_max = atoi(sql_result[row*ncols+9]);
                const int z_M_min    = atoi(sql_result[row*ncols+10]);
                const int z_M_max    = atoi(sql_result[row*ncols+11]);
                const int z_H              = atoi(sql_result[row*ncols+col_z_H]);
                const int N_NEO            = atoi(sql_result[row*ncols+col_N_NEO]);
                const int N_PHO            = atoi(sql_result[row*ncols+col_N_PHO]);
                const int NEO_in_field     = atoi(sql_result[row*ncols+col_NEO_in_field]);
                const int PHO_in_field     = atoi(sql_result[row*ncols+col_PHO_in_field]);
                const double eta_NEO       = atof(sql_result[row*ncols+col_eta_NEO]);
                // const double sigma_eta_NEO = atof(sql_result[row*ncols+col_sigma_eta_NEO]);
                const double eta_PHO       = atof(sql_result[row*ncols+col_eta_PHO]);
                // const double sigma_eta_PHO = atof(sql_result[row*ncols+col_sigma_eta_PHO]);
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
                    NEO_eta_field_ws_vec[index.get()].resize(z_L_size);
                    for (unsigned int l=0; l<z_L_size; ++l) {
                        NEO_eta_field_ws_vec[index.get()][l] = new orsa::WeightedStatistic<double>;
                    }
                    //
                    NEO_eta_detect_ws_vec.resize(z_H_vec.size());
                    NEO_eta_detect_ws_vec[index.get()].resize(z_L_size);
                    for (unsigned int l=0; l<z_L_size; ++l) {
                        NEO_eta_detect_ws_vec[index.get()][l] = new orsa::WeightedStatistic<double>;
                    }
                    //
                    PHO_eta_field_ws_vec.resize(z_H_vec.size());
                    PHO_eta_field_ws_vec[index.get()].resize(z_L_size);
                    for (unsigned int l=0; l<z_L_size; ++l) {
                        PHO_eta_field_ws_vec[index.get()][l] = new orsa::WeightedStatistic<double>;
                    }
                    //
                    PHO_eta_detect_ws_vec.resize(z_H_vec.size());
                    PHO_eta_detect_ws_vec[index.get()].resize(z_L_size);
                    for (unsigned int l=0; l<z_L_size; ++l) {
                        PHO_eta_detect_ws_vec[index.get()][l] = new orsa::WeightedStatistic<double>;
                    }
                }
                //
                std::vector<int> z_L_vec;
                {
                    // assume step is the same for node,peri,M
                    // #warning this is grain-size dependent!
                    z_L_vec.push_back(((z_node_min+z_peri_min+z_M_min)%360000)/30000);
                    z_L_vec.push_back(((z_node_max+z_peri_min+z_M_min)%360000)/30000);
                    z_L_vec.push_back(((z_node_max+z_peri_max+z_M_min)%360000)/30000);
                }
                
                for (unsigned int l=0; l<z_L_vec.size(); ++l) {
                    const int index_L = z_L_vec[l];
                    if (N_NEO>0) {
                        const double NEO_eta_field       = (double)(NEO_in_field)/(double)(N_NEO);
                        /* const double p                   = (double)(NEO_in_field+1)/(double)(N_NEO+2);
                           const double NEO_sigma_eta_field = sqrt(p*(1-p)/N_NEO);
                           if (NEO_sigma_eta_field>0) NEO_eta_field_ws_vec[index.get()]->insert(NEO_eta_field,orsa::square(1.0/NEO_sigma_eta_field));
                        */
                        // try to use weight=N_NEO
                        NEO_eta_field_ws_vec[index.get()][index_L]->insert(NEO_eta_field,N_NEO);
                        //
                        const double NEO_eta_detect       = eta_NEO;
                        /* 
                           const double NEO_sigma_eta_detect = sigma_eta_NEO;
                           if (NEO_sigma_eta_detect>0) NEO_eta_detect_ws_vec[index.get()]->insert(NEO_eta_detect,orsa::square(1.0/NEO_sigma_eta_detect));
                        */
                        // try to use weight=N_NEO
                        NEO_eta_detect_ws_vec[index.get()][index_L]->insert(NEO_eta_detect,N_NEO);
                    }
                    //
                    if (N_PHO>0) {
                        const double PHO_eta_field       = (double)(PHO_in_field)/(double)(N_PHO);
                        /* const double p                   = (double)(PHO_in_field+1)/(double)(N_PHO+2);
                           const double PHO_sigma_eta_field = sqrt(p*(1-p)/N_PHO);
                           if (PHO_sigma_eta_field>0) PHO_eta_field_ws_vec[index.get()]->insert(PHO_eta_field,orsa::square(1.0/PHO_sigma_eta_field));
                        */
                        // try to use weight=N_PHO
                        PHO_eta_field_ws_vec[index.get()][index_L]->insert(PHO_eta_field,N_PHO);
                        //
                        const double PHO_eta_detect       = eta_PHO;
                        /* const double PHO_sigma_eta_detect = sigma_eta_PHO;
                           if (PHO_sigma_eta_detect>0) PHO_eta_detect_ws_vec[index.get()]->insert(PHO_eta_detect,orsa::square(1.0/PHO_sigma_eta_detect));
                        */
                        // try to use weight=N_PHO
                        PHO_eta_detect_ws_vec[index.get()][index_L]->insert(PHO_eta_detect,N_PHO);
                    }
                }
            }
            
            // now, real output to file
            
            for (unsigned int index=0; index<z_H_vec.size(); ++index) {
                for (unsigned int index_L=0; index_L<z_L_size; ++index_L) {
                    // ORSA_DEBUG("index: %i  index_L: %i  entries: %Zi",index,index_L,NEO_eta_detect_ws_vec[index][index_L]->entries().get_mpz_t());
                    if (NEO_eta_detect_ws_vec[index][index_L]->entries()==0) continue; 
                    const double NEO_eta_obs =
                        NEO_eta_field_ws_vec[index][index_L]->average() *
                        NEO_eta_detect_ws_vec[index][index_L]->average();
                    const double NEO_sigma_eta_obs =
                        sqrt(orsa::square(NEO_eta_field_ws_vec[index][index_L]->average()*NEO_eta_detect_ws_vec[index][index_L]->standardDeviation()) +
                             orsa::square(NEO_eta_detect_ws_vec[index][index_L]->average()*NEO_eta_field_ws_vec[index][index_L]->standardDeviation()));
                    gmp_fprintf(stdout,
                                "NEO %5i %5i %4i %4i %6i %6i %6i %6i %3i %.2e %.2e %4Zi %.2e %.2e %4Zi %.2e %.2e\n",
                                z_a_min,
                                z_a_max,
                                z_e_min,
                                z_e_max,
                                z_i_min,
                                z_i_max,
                                30000*index_L,
                                30000*(1+index_L),
                                z_H_vec[index],
                                NEO_eta_field_ws_vec[index][index_L]->average(),
                                NEO_eta_field_ws_vec[index][index_L]->standardDeviation(),
                                NEO_eta_field_ws_vec[index][index_L]->entries().get_mpz_t(),
                                NEO_eta_detect_ws_vec[index][index_L]->average(),
                                NEO_eta_detect_ws_vec[index][index_L]->standardDeviation(),
                                NEO_eta_detect_ws_vec[index][index_L]->entries().get_mpz_t(),
                                NEO_eta_obs,
                                NEO_sigma_eta_obs);
                }
            }
            
            for (unsigned int index=0; index<z_H_vec.size(); ++index) {
                for (unsigned int index_L=0; index_L<z_L_size; ++index_L) {
                    // ORSA_DEBUG("index: %i  index_L: %i  entries: %Zi",index,index_L,NEO_eta_detect_ws_vec[index][index_L]->entries().get_mpz_t());
                    if (PHO_eta_detect_ws_vec[index][index_L]->entries()==0) continue;
                    const double PHO_eta_obs =
                        PHO_eta_field_ws_vec[index][index_L]->average() *
                        PHO_eta_detect_ws_vec[index][index_L]->average();
                    const double PHO_sigma_eta_obs =
                        sqrt(orsa::square(PHO_eta_field_ws_vec[index][index_L]->average()*PHO_eta_detect_ws_vec[index][index_L]->standardDeviation()) +
                             orsa::square(PHO_eta_detect_ws_vec[index][index_L]->average()*PHO_eta_field_ws_vec[index][index_L]->standardDeviation()));
                    gmp_fprintf(stdout,
                                "PHO %5i %5i %4i %4i %6i %6i %6i %6i %3i %.2e %.2e %4Zi %.2e %.2e %4Zi %.2e %.2e\n",
                                z_a_min,
                                z_a_max,
                                z_e_min,
                                z_e_max,
                                z_i_min,
                                z_i_max,
                                30000*index_L,
                                30000*(1+index_L),
                                z_H_vec[index],
                                PHO_eta_field_ws_vec[index][index_L]->average(),
                                PHO_eta_field_ws_vec[index][index_L]->standardDeviation(),
                                PHO_eta_field_ws_vec[index][index_L]->entries().get_mpz_t(),
                                PHO_eta_detect_ws_vec[index][index_L]->average(),
                                PHO_eta_detect_ws_vec[index][index_L]->standardDeviation(),
                                PHO_eta_detect_ws_vec[index][index_L]->entries().get_mpz_t(),
                                PHO_eta_obs,
                                PHO_sigma_eta_obs);
                }
            }
            
            sqlite3_free_table(sql_result);
        }
    
        // close db
        sqlite3_close(db);
    
    } // loop on files
    
    return 0;
}

