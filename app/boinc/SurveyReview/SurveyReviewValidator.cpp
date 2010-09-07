#include <string>
#include <vector>
#include <cmath>

#include <gmpxx.h>

#include "boinc_db.h"
#include "error_numbers.h"
#include "parse.h"
#include "util.h"
#include "filesys.h"

#include "sched_util.h"
#include "sched_config.h"
#include "sched_msgs.h"
#include "validator.h"
#include "validate_util.h"

// SQLite3
#include "sqlite3.h"

// same as assimilator, should share in a common header...
int get_input_file_infos(const WORKUNIT & wu, 
                         std::vector<FILE_INFO>& fis) {
    char tag[256], path[1024];
    bool is_tag;
    MIOFILE mf;
    std::string name;
    mf.init_buf_read(wu.xml_doc);
    XML_PARSER xp(&mf);
    fis.clear();
    while (!xp.get(tag, sizeof(tag), is_tag)) {
        if (!is_tag) continue;
        if (!strcmp(tag, "file_ref")) {
            FILE_INFO fi;
            int retval = fi.parse(xp);
            if (retval) return retval;
            dir_hier_path(fi.name.c_str(), config.download_dir, config.uldl_dir_fanout, path);
            fi.path = path;
            fis.push_back(fi);
        }
    }
    return 0;
}

int init_result(RESULT & result,
                void * & data) {
    
    FILE_INFO fi;
    int retval;
    
    retval = get_output_file_path(result, fi.path);
    if (retval) return retval;
    
    DB_WORKUNIT wu;
    wu.lookup_id(g_wup->id);

    log_messages.printf(MSG_DEBUG, 
                        "[WORKUNIT#%d %s] opening SQLite db [%s]\n", 
                        wu.id,  
                        wu.name,  
                        fi.path.c_str());
    
    sqlite3     * db;
    int           rc;
    char        * zErr;
    
    // rc = sqlite3_open(fi.path.c_str(),&db);
    rc = sqlite3_open_v2(fi.path.c_str(),&db,SQLITE_OPEN_READONLY,NULL);
    //
    if (rc) {
        fprintf(stderr,"Can't open db: %s\n",sqlite3_errmsg(db));
        sqlite3_close(db);
        return rc;
    }

    int nrows, ncols;
    char * * sql_result;
    rc = sqlite3_get_table(db,"pragma integrity_check",&sql_result,&nrows,&ncols,&zErr);
    //
    if (rc != SQLITE_OK) {
        if (zErr != NULL) {
            log_messages.printf(MSG_CRITICAL, 
                                "[WORKUNIT#%d %s] SQL error: %s\n", 
                                wu.id,  
                                wu.name,
                                zErr); 
            sqlite3_close(db);
            return rc;
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
        log_messages.printf(MSG_CRITICAL, 
                            "[WORKUNIT#%d %s] SQLite problem: integrity_check failed\n", 
                            wu.id,  
                            wu.name); 
        sqlite3_close(db);
        return -1; 
    }
    
    {
        // copy "grid" table into TEMP table "memgrid", stored in memory, for better performance
        
        rc = sqlite3_exec(db,"pragma temp_store=memory",NULL,NULL,&zErr);
        //
        if (rc != SQLITE_OK) {
            if (zErr != NULL) {
                log_messages.printf(MSG_CRITICAL, 
                                    "[WORKUNIT#%d %s] SQL error: %s\n", 
                                    wu.id,  
                                    wu.name,
                                    zErr); 
                sqlite3_close(db);
                return rc;
            }
        }
        
        rc = sqlite3_exec(db,"create temp table memgrid as select * from grid",NULL,NULL,&zErr);
        //
        if (rc != SQLITE_OK) {
            if (zErr != NULL) {
                log_messages.printf(MSG_CRITICAL, 
                                    "[WORKUNIT#%d %s] SQL error: %s\n", 
                                    wu.id,  
                                    wu.name,
                                    zErr); 
                sqlite3_close(db);
                return rc;
            }
        }
        
    }
    
    data = (void *) db;
    
    return 0;
}

int compare_results(RESULT & /* r1 */,
                    void * data1,
                    const RESULT & /* r2 */,
                    void * data2,
                    bool & match) {
    
    // set match to false while initializing, in case something goes wrong
    match = false;
    
    sqlite3 * db1 = (sqlite3 *) data1;
    sqlite3 * db2 = (sqlite3 *) data2;

    DB_WORKUNIT wu;
    wu.lookup_id(g_wup->id);
    
    // read grid.dat input file
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
        bool done=false;
        std::vector<FILE_INFO> fis;
        get_input_file_infos(wu,fis);
        std::vector<FILE_INFO>::const_iterator it = fis.begin();
        while (it != fis.end()) {
            // test if file exists
            FILE * fp = fopen((*it).path.c_str(),"r"); 
            if (fp) {
                log_messages.printf(MSG_DEBUG, 
                                    "[WORKUNIT#%d %s] trying to read file %s\n", 
                                    wu.id,  
                                    wu.name,  
                                    (*it).path.c_str()); 
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

                        log_messages.printf(MSG_DEBUG, 
                                            "[WORKUNIT#%d %s] successfully read grid.dat file\n", 
                                            wu.id,  
                                            wu.name);
                        
                        /* log_messages.printf(MSG_DEBUG, 
                           "[WORKUNIT#%d %s] %.5f %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n", 
                           wu.id,  
                           wu.name,
                           JD,
                           numGen,
                           z_a_min, z_a_max, z_a_delta,
                           z_e_min, z_e_max, z_e_delta,
                           z_i_min, z_i_max, z_i_delta,
                           z_node_min, z_node_max, z_node_delta,
                           z_peri_min, z_peri_max, z_peri_delta,
                           z_M_min, z_M_max, z_M_delta,
                           z_H_min, z_H_max, z_H_delta);
                        */

                        done=true;
                        break;
                    }
                }
                fclose(fp);
            } else { 
                log_messages.printf(MSG_CRITICAL, 
                                    "[WORKUNIT#%d %s] Couldn't open %s\n", 
                                    wu.id,  
                                    wu.name,  
                                    (*it).path.c_str()); 
                return 1; 
            }
            if (done) break;
            ++it;  
        }
        if (!done) {
            log_messages.printf(MSG_CRITICAL, 
                                "[WORKUNIT#%d %s] Couldn't open grid.dat input file\n", 
                                wu.id,  
                                wu.name); 
        }
    }
    
    // H not included in totalIterations; each iteration runs all H values
    // const int size_H = 1 + (z_H_max-z_H_min)/z_H_delta;
    
    match=true;
    char sql_line[1024];
    int rc1, rc2;
    char * zErr;
    char * * sql_result1;
    char * * sql_result2;
    int nrows1, ncols1;
    int nrows2, ncols2;
    // 
    int col_N_NEO=-1;
    int col_N_PHO=-1;
    int col_NEO_in_field=-1;
    int col_PHO_in_field=-1;
    int col_eta_NEO=-1;
    int col_sigma_eta_NEO=-1;
    int col_eta_PHO=-1;
    int col_sigma_eta_PHO=-1;
    //
    {
        // first perform basic checks outside loop
        
        // just to get the columns
        // sprintf(sql_line,"SELECT * FROM grid LIMIT 1");
        sprintf(sql_line,"SELECT * FROM memgrid LIMIT 1");
        
        rc1 = sqlite3_get_table(db1,sql_line,&sql_result1,&nrows1,&ncols1,&zErr);
        //
        if (rc1 != SQLITE_OK) {
            if (zErr != NULL) {
                log_messages.printf(MSG_CRITICAL, 
                                    "[WORKUNIT#%d %s] SQL error: %s\n", 
                                    wu.id,  
                                    wu.name,
                                    zErr); 
                sqlite3_free(zErr);
                match=false;
            }
        }
        
        rc2 = sqlite3_get_table(db2,sql_line,&sql_result2,&nrows2,&ncols2,&zErr);
        //
        if (rc2 != SQLITE_OK) {
            if (zErr != NULL) {
                log_messages.printf(MSG_CRITICAL, 
                                    "[WORKUNIT#%d %s] SQL error: %s\n", 
                                    wu.id,  
                                    wu.name,
                                    zErr); 
                sqlite3_free(zErr);
                match=false;
            }
        }

        if ( (nrows1==0) &&
             (nrows2==0) ) {
            // empty and matching, so cleanup and return            
            sqlite3_free_table(sql_result1);
            sqlite3_free_table(sql_result2);
            return 0;
        }
        
        const int ncols = ncols1;
        for (int col=0; col<ncols; ++col) {
            if (sql_result1[col] == std::string("N_NEO")) {
                col_N_NEO=col;
            } else if (sql_result1[col] == std::string("N_PHO")) {
                col_N_PHO=col;
            } else if (sql_result1[col] == std::string("NEO_in_field")) {
                col_NEO_in_field=col;
            } else if (sql_result1[col] == std::string("PHO_in_field")) {
                col_PHO_in_field=col;
            } else if (sql_result1[col] == std::string("eta_NEO")) {
                col_eta_NEO=col;
            } else if (sql_result1[col] == std::string("sigma_eta_NEO")) {
                col_sigma_eta_NEO=col;
            } else if (sql_result1[col] == std::string("eta_PHO")) {
                col_eta_PHO=col;
            } else if (sql_result1[col] == std::string("sigma_eta_PHO")) {
                col_sigma_eta_PHO=col;
            }
        }
        
        sqlite3_free_table(sql_result1);
        sqlite3_free_table(sql_result2);

        if ( (col_N_NEO==-1) ||
             (col_N_PHO==-1) ||
             (col_NEO_in_field==-1) ||
             (col_PHO_in_field==-1) ||
             (col_eta_NEO==-1) ||
             (col_sigma_eta_NEO==-1) ||
             (col_eta_PHO==-1) ||
             (col_sigma_eta_PHO==-1) ) {
            log_messages.printf(MSG_CRITICAL, 
                                "[WORKUNIT#%d %s] could not find some columns\n", 
                                wu.id,  
                                wu.name); 
            match=false;
        } else {
            /* log_messages.printf(MSG_DEBUG, 
               "[WORKUNIT#%d %s] columns: %i %i %i %i %i %i %i %i\n",
               wu.id,  
               wu.name,
               col_N_NEO,
               col_N_PHO,
               col_NEO_in_field,
               col_PHO_in_field,
               col_eta_NEO,
               col_sigma_eta_NEO,
               col_eta_PHO,
               col_sigma_eta_PHO);
            */
        }
        
        // quick exit
        if (!match) return 0;
    }
    //
    /* for (int z_a=z_a_min; z_a<z_a_max; z_a+=z_a_delta) {
       for (int z_e=z_e_min; z_e<z_e_max; z_e+=z_e_delta) {
       for (int z_i=z_i_min; z_i<z_i_max; z_i+=z_i_delta) {
       for (int z_node=z_node_min; z_node<z_node_max; z_node+=z_node_delta) {
       for (int z_peri=z_peri_min; z_peri<z_peri_max; z_peri+=z_peri_delta) {
       for (int z_M=z_M_min; z_M<z_M_max; z_M+=z_M_delta) {
    */
    
    // different loop for H, checking BOTH extremes
    // for (int z_H=z_H_min; z_H<=z_H_max; z_H+=z_H_delta) {
    
    /* sprintf(sql_line,
       "SELECT * FROM grid WHERE z_a_min=%i and z_a_max=%i and z_e_min=%i and z_e_max=%i and z_i_min=%i and z_i_max=%i and z_node_min=%i and z_node_max=%i and z_peri_min=%i and z_peri_max=%i and z_M_min=%i and z_M_max=%i and z_H=%i",
       z_a,z_a+z_a_delta,
       z_e,z_e+z_e_delta,
       z_i,z_i+z_i_delta,
       z_node,z_node+z_node_delta,
       z_peri,z_peri+z_peri_delta,
       z_M,z_M+z_M_delta,
       z_H);
    */
    //
    /* sprintf(sql_line,
       "SELECT * FROM grid WHERE z_a_min=%i and z_a_max=%i and z_e_min=%i and z_e_max=%i and z_i_min=%i and z_i_max=%i and z_node_min=%i and z_node_max=%i and z_peri_min=%i and z_peri_max=%i and z_M_min=%i and z_M_max=%i order by z_H",
       z_a,z_a+z_a_delta,
       z_e,z_e+z_e_delta,
       z_i,z_i+z_i_delta,
       z_node,z_node+z_node_delta,
       z_peri,z_peri+z_peri_delta,
       z_M,z_M+z_M_delta);
    */
    //
    /* sprintf(sql_line,
       "SELECT * FROM memgrid WHERE z_a_min=%i and z_a_max=%i and z_e_min=%i and z_e_max=%i and z_i_min=%i and z_i_max=%i and z_node_min=%i and z_node_max=%i and z_peri_min=%i and z_peri_max=%i and z_M_min=%i and z_M_max=%i order by z_H",
       z_a,z_a+z_a_delta,
       z_e,z_e+z_e_delta,
       z_i,z_i+z_i_delta,
       z_node,z_node+z_node_delta,
       z_peri,z_peri+z_peri_delta,
       z_M,z_M+z_M_delta);
    */
    //
    sprintf(sql_line,
            "SELECT * FROM memgrid order by z_a_min,z_a_max,z_e_min,z_e_max,z_i_min,z_i_max,z_node_min,z_node_max,z_peri_min,z_peri_max,z_M_min,z_M_max,z_H");
    
    rc1 = sqlite3_get_table(db1,sql_line,&sql_result1,&nrows1,&ncols1,&zErr);
    //
    if (rc1 != SQLITE_OK) {
        if (zErr != NULL) {
            log_messages.printf(MSG_CRITICAL, 
                                "[WORKUNIT#%d %s] SQL error: %s\n", 
                                wu.id,  
                                wu.name,
                                zErr); 
            sqlite3_free(zErr);
            match=false;
            sqlite3_free_table(sql_result1);
            // sqlite3_free_table(sql_result2);
            return 1;
        }
    }
    // 0 is admisible if bin is not NEO
    /* if (nrows1 > size_H) {
       log_messages.printf(MSG_CRITICAL, 
       "[WORKUNIT#%d %s] corrupted result\n", 
       wu.id,  
       wu.name); 
       match=false;
       sqlite3_free_table(sql_result1);
       // sqlite3_free_table(sql_result2);
       return 1;
       }
    */
    
    rc2 = sqlite3_get_table(db2,sql_line,&sql_result2,&nrows2,&ncols2,&zErr);
    //
    if (rc2 != SQLITE_OK) {
        if (zErr != NULL) {
            log_messages.printf(MSG_CRITICAL, 
                                "[WORKUNIT#%d %s] SQL error: %s\n", 
                                wu.id,  
                                wu.name,
                                zErr); 
            sqlite3_free(zErr);
            match=false;
            sqlite3_free_table(sql_result1);
            sqlite3_free_table(sql_result2);
            return 1;
        }
    }
    // 0 is admissible if bin is not NEO
    /* if (nrows2 > size_H) {
       log_messages.printf(MSG_CRITICAL, 
       "[WORKUNIT#%d %s] corrupted result\n", 
       wu.id,  
       wu.name); 
       match=false;
       sqlite3_free_table(sql_result1);
       sqlite3_free_table(sql_result2);
       return 1;
       }
    */

    // checked already
    /* if ( (nrows1==0) &&
       (nrows2==0) ) {
       continue;
       }
    */
    
    // comparison
    if ( (nrows1 != nrows2) ||
         (ncols1 != ncols2) ) {
        log_messages.printf(MSG_CRITICAL, 
                            "[WORKUNIT#%d %s] corrupted result\n", 
                            wu.id,  
                            wu.name); 
        match=false;
        sqlite3_free_table(sql_result1);
        sqlite3_free_table(sql_result2);
        return 1;
    }
    // ncols1 and ncols2 are equal now
    const int nrows = nrows1;
    const int ncols = ncols1;
    /* if ( (nrows!=last_nrows) ||
       (ncols!=last_ncols) ) {
       // print only if changed...
       log_messages.printf(MSG_DEBUG,"[WORKUNIT#%d %s] z_a_min=%i and z_a_max=%i and z_e_min=%i and z_e_max=%i and z_i_min=%i and z_i_max=%i and z_node_min=%i and z_node_max=%i and z_peri_min=%i and z_peri_max=%i and z_M_min=%i and z_M_max=%i\n",
       wu.id,wu.name,
       z_a,z_a+z_a_delta,
       z_e,z_e+z_e_delta,
       z_i,z_i+z_i_delta,
       z_node,z_node+z_node_delta,
       z_peri,z_peri+z_peri_delta,
       z_M,z_M+z_M_delta);
       log_messages.printf(MSG_DEBUG,"[WORKUNIT#%d %s] nrows=%i ncols=%i\n",wu.id,wu.name,nrows,ncols);
       last_nrows=nrows;
       last_ncols=ncols;
       }
    */
    
    // for (int col=0; col<ncols; ++col) {
    for (int row=1; row<=nrows; ++row) {
        {
            // if (sql_result1[col] == std::string("N_NEO")) {
            // log_messages.printf(MSG_DEBUG,"[WORKUNIT#%d %s] testing N_NEO\n",wu.id,wu.name);
            int col = col_N_NEO;
            int N_NEO1 = atoi(sql_result1[row*ncols+col]);
            int N_NEO2 = atoi(sql_result2[row*ncols+col]);
            if (N_NEO1 != N_NEO2) {
                match=false;
                log_messages.printf(MSG_DEBUG,"[WORKUNIT#%d %s] not maching: N_NEO1=%i N_NEO2=%i\n",wu.id,wu.name,N_NEO1,N_NEO2);
            }
        }
        {
            // } else if (sql_result1[col] == std::string("N_PHO")) {
            // log_messages.printf(MSG_DEBUG,"[WORKUNIT#%d %s] testing N_PHO\n",wu.id,wu.name);
            int col = col_N_PHO;
            int N_PHO1 = atoi(sql_result1[row*ncols+col]);
            int N_PHO2 = atoi(sql_result2[row*ncols+col]);
            if (N_PHO1 != N_PHO2) {
                match=false;
                log_messages.printf(MSG_DEBUG,"[WORKUNIT#%d %s] not maching: N_PHO1=%i N_PHO2=%i\n",wu.id,wu.name,N_PHO1,N_PHO2);
            }
        }
        {
            // } else if (sql_result1[col] == std::string("NEO_in_field")) {
            // log_messages.printf(MSG_DEBUG,"[WORKUNIT#%d %s] testing NEO_in_field\n",wu.id,wu.name);
            int col = col_NEO_in_field;
            int NEO_in_field1 = atoi(sql_result1[row*ncols+col]);
            int NEO_in_field2 = atoi(sql_result2[row*ncols+col]);
            if (NEO_in_field1 != NEO_in_field2) {
                match=false;
                log_messages.printf(MSG_DEBUG,"[WORKUNIT#%d %s] not maching: NEO_in_field1=%i NEO_in_field2=%i\n",wu.id,wu.name,NEO_in_field1,NEO_in_field2);
            }
        }
        {
            // } else if (sql_result1[col] == std::string("PHO_in_field")) {
            // log_messages.printf(MSG_DEBUG,"[WORKUNIT#%d %s] testing PHO_in_field\n",wu.id,wu.name);
            int col = col_PHO_in_field;
            int PHO_in_field1 = atoi(sql_result1[row*ncols+col]);
            int PHO_in_field2 = atoi(sql_result2[row*ncols+col]);
            if (PHO_in_field1 != PHO_in_field2) {
                match=false;
                log_messages.printf(MSG_DEBUG,"[WORKUNIT#%d %s] not maching: PHO_in_field1=%i PHO_in_field2=%i\n",wu.id,wu.name,PHO_in_field1,PHO_in_field2);
            }
        }
        {
            // } else if (sql_result1[col] == std::string("eta_NEO")) {
            // log_messages.printf(MSG_DEBUG,"[WORKUNIT#%d %s] testing eta_NEO\n",wu.id,wu.name);
            int col = col_eta_NEO;
            double eta_NEO1 = atof(sql_result1[row*ncols+col]);
            double eta_NEO2 = atof(sql_result2[row*ncols+col]);
            if (eta_NEO2!=eta_NEO1) { // this takes care of the {0,0} case
                if (fabs((eta_NEO2-eta_NEO1)/(fabs(eta_NEO1)+fabs(eta_NEO2))) > 0.01) {
                    match=false;
                    log_messages.printf(MSG_DEBUG,"[WORKUNIT#%d %s] not maching: eta_NEO1=%g eta_NEO2=%g\n",wu.id,wu.name,eta_NEO1,eta_NEO2);
                }
            }
        }
        {
            // } else if (sql_result1[col] == std::string("sigma_eta_NEO")) {
            // log_messages.printf(MSG_DEBUG,"[WORKUNIT#%d %s] testing sigma_eta_NEO\n",wu.id,wu.name);
            int col = col_sigma_eta_NEO;
            double sigma_eta_NEO1 = atof(sql_result1[row*ncols+col]);
            double sigma_eta_NEO2 = atof(sql_result2[row*ncols+col]);
            if (sigma_eta_NEO2!=sigma_eta_NEO1) { // this takes care of the {0,0} case
                if (fabs((sigma_eta_NEO2-sigma_eta_NEO1)/(fabs(sigma_eta_NEO1)+fabs(sigma_eta_NEO2))) > 0.01) {
                    match=false;
                    log_messages.printf(MSG_DEBUG,"[WORKUNIT#%d %s] not maching: sigma_eta_NEO1=%g sigma_eta_NEO2=%g\n",wu.id,wu.name,sigma_eta_NEO1,sigma_eta_NEO2);
                }
            }
        }
        {
            // } else if (sql_result1[col] == std::string("eta_PHO")) {
            // log_messages.printf(MSG_DEBUG,"[WORKUNIT#%d %s] testing eta_PHO\n",wu.id,wu.name);
            int col = col_eta_PHO;
            double eta_PHO1 = atof(sql_result1[row*ncols+col]);
            double eta_PHO2 = atof(sql_result2[row*ncols+col]);
            if (eta_PHO2!=eta_PHO1) { // this takes care of the {0,0} case
                if (fabs((eta_PHO2-eta_PHO1)/(fabs(eta_PHO1)+fabs(eta_PHO2))) > 0.01) {
                    match=false;
                    log_messages.printf(MSG_DEBUG,"[WORKUNIT#%d %s] not maching: eta_PHO1=%g eta_PHO2=%g\n",wu.id,wu.name,eta_PHO1,eta_PHO2);
                }
            }
        }
        {
            // } else if (sql_result1[col] == std::string("sigma_eta_PHO")) {
            // log_messages.printf(MSG_DEBUG,"[WORKUNIT#%d %s] testing sigma_eta_PHO\n",wu.id,wu.name);
            int col = col_sigma_eta_PHO;
            double sigma_eta_PHO1 = atof(sql_result1[row*ncols+col]);
            double sigma_eta_PHO2 = atof(sql_result2[row*ncols+col]);
            if (sigma_eta_PHO2!=sigma_eta_PHO1) { // this takes care of the {0,0} case
                if (fabs((sigma_eta_PHO2-sigma_eta_PHO1)/(fabs(sigma_eta_PHO1)+fabs(sigma_eta_PHO2))) > 0.01) {
                    match=false;
                    log_messages.printf(MSG_DEBUG,"[WORKUNIT#%d %s] not maching: sigma_eta_PHO1=%g sigma_eta_PHO2=%g\n",wu.id,wu.name,sigma_eta_PHO1,sigma_eta_PHO2);
                }
            }
        }
    }
    
    sqlite3_free_table(sql_result1);
    sqlite3_free_table(sql_result2);
    
    // quick exit
    if (!match) return 0;
    // }

/* }
   }
   }
   }
   }
   }
*/
    return 0;
}

int cleanup_result(RESULT const & /* result */,
                   void * data) {
    sqlite3_close((sqlite3 *)data);
    return 0;
}

double compute_granted_credit(WORKUNIT & /* wu */, vector<RESULT> & /* results */) {
    // return median_mean_credit(wu, results);
    return 25.0;
}
