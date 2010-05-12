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

    sqlite3     * db;
    int           rc;
    
    rc = sqlite3_open(fi.path.c_str(),&db);
    //
    if (rc) {
        fprintf(stderr,"Can't open db: %s\n",sqlite3_errmsg(db));
        sqlite3_close(db);
        return rc;
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
    
    match=true;
    char sql_line[1024];
    int rc1, rc2;
    char * zErr;
    char * * sql_result1;
    char * * sql_result2;
    int nrows1, ncols1;
    int nrows2, ncols2;
    std::string sql;
    for (int z_a=z_a_min; z_a<z_a_max; z_a+=z_a_delta) {
        for (int z_e=z_e_min; z_e<z_e_max; z_e+=z_e_delta) {
            for (int z_i=z_i_min; z_i<z_i_max; z_i+=z_i_delta) {
                for (int z_node=z_node_min; z_node<z_node_max; z_node+=z_node_delta) {
                    for (int z_peri=z_peri_min; z_peri<z_peri_max; z_peri+=z_peri_delta) {
                        for (int z_M=z_M_min; z_M<z_M_max; z_M+=z_M_delta) {
                            // different loop for H, checking BOTH extremes
                            for (int z_H=z_H_min; z_H<=z_H_max; z_H+=z_H_delta) {

                                sprintf(sql_line,
                                        "SELECT * FROM grid WHERE z_a_min=%i and z_a_max=%i and z_e_min=%i and z_e_max=%i and z_i_min=%i and z_i_max=%i and z_node_min=%i and z_node_max=%i and z_peri_min=%i and z_peri_max=%i and z_M_min=%i and z_M_max=%i and z_H=%i",
                                        z_a,z_a+z_a_delta,
                                        z_e,z_e+z_e_delta,
                                        z_i,z_i+z_i_delta,
                                        z_node,z_node+z_node_delta,
                                        z_peri,z_peri+z_peri_delta,
                                        z_M,z_M+z_M_delta,
                                        z_H);

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
                                        return rc1;
                                    }
                                }
                                //
                                if (nrows1 != 1) {
                                    log_messages.printf(MSG_CRITICAL, 
                                                        "[WORKUNIT#%d %s] corrupted result\n", 
                                                        wu.id,  
                                                        wu.name); 
                                    match=false;
                                    return 1;
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
                                        return rc2;
                                    }
                                }
                                //
                                if (nrows2 != 1) {
                                    log_messages.printf(MSG_CRITICAL, 
                                                        "[WORKUNIT#%d %s] corrupted result\n", 
                                                        wu.id,  
                                                        wu.name); 
                                    match=false;
                                    return 1;
                                }
                                
                                // comparison
                                if (ncols1 != ncols2) {
                                    log_messages.printf(MSG_CRITICAL, 
                                                        "[WORKUNIT#%d %s] corrupted result\n", 
                                                        wu.id,  
                                                        wu.name); 
                                    match=false;
                                    return 1;
                                }
                                const int ncols = ncols1;
                                for (int col=0; col<ncols; ++col) {
                                    if (sql_result1[col] == std::string("N_NEO")) {
                                        int N_NEO1 = atoi(sql_result1[ncols+col]);
                                        int N_NEO2 = atoi(sql_result2[ncols+col]);
                                        if (N_NEO1 != N_NEO2) match=false;
                                    } else if (sql_result1[col] == std::string("N_PHO")) {
                                        int N_PHO1 = atoi(sql_result1[ncols+col]);
                                        int N_PHO2 = atoi(sql_result2[ncols+col]);
                                        if (N_PHO1 != N_PHO2) match=false;
                                    } else if (sql_result1[col] == std::string("NEO_in_field")) {
                                        int NEO_in_field1 = atoi(sql_result1[ncols+col]);
                                        int NEO_in_field2 = atoi(sql_result2[ncols+col]);
                                        if (NEO_in_field1 != NEO_in_field2) match=false;
                                    } else if (sql_result1[col] == std::string("PHO_in_field")) {
                                        int PHO_in_field1 = atoi(sql_result1[ncols+col]);
                                        int PHO_in_field2 = atoi(sql_result2[ncols+col]);
                                        if (PHO_in_field1 != PHO_in_field2) match=false;
                                    } else if (sql_result1[col] == std::string("eta_NEO")) {
                                        double eta_NEO1 = atof(sql_result1[ncols+col]);
                                        double eta_NEO2 = atof(sql_result2[ncols+col]);
                                        if (eta_NEO2!=eta_NEO1) { // this takes care of the {0,0} case
                                            if (fabs((eta_NEO2-eta_NEO1)/(fabs(eta_NEO1)+fabs(eta_NEO2))) > 0.01) match=false;
                                        }
                                    } else if (sql_result1[col] == std::string("sigma_eta_NEO")) {
                                        double sigma_eta_NEO1 = atof(sql_result1[ncols+col]);
                                        double sigma_eta_NEO2 = atof(sql_result2[ncols+col]);
                                        if (sigma_eta_NEO2!=sigma_eta_NEO1) { // this takes care of the {0,0} case
                                            if (fabs((sigma_eta_NEO2-sigma_eta_NEO1)/(fabs(sigma_eta_NEO1)+fabs(sigma_eta_NEO2))) > 0.01) match=false;
                                        }
                                    } else if (sql_result1[col] == std::string("eta_PHO")) {
                                        double eta_PHO1 = atof(sql_result1[ncols+col]);
                                        double eta_PHO2 = atof(sql_result2[ncols+col]);
                                        if (eta_PHO2!=eta_PHO1) { // this takes care of the {0,0} case
                                            if (fabs((eta_PHO2-eta_PHO1)/(fabs(eta_PHO1)+fabs(eta_PHO2))) > 0.01) match=false;
                                        }
                                    } else if (sql_result1[col] == std::string("sigma_eta_PHO")) {
                                        double sigma_eta_PHO1 = atof(sql_result1[ncols+col]);
                                        double sigma_eta_PHO2 = atof(sql_result2[ncols+col]);
                                        if (sigma_eta_PHO2!=sigma_eta_PHO1) { // this takes care of the {0,0} case
                                            if (fabs((sigma_eta_PHO2-sigma_eta_PHO1)/(fabs(sigma_eta_PHO1)+fabs(sigma_eta_PHO2))) > 0.01) match=false;
                                        }
                                    }
                                }
                                
                                sqlite3_free_table(sql_result1);
                                sqlite3_free_table(sql_result2);
                                
                                // quick exit
                                if (!match) return 0;
                            }
                        }
                    }
                }
            }
        }
    }
    
    return 0;
}

int cleanup_result(RESULT const & /* result */,
                   void * data) {
    sqlite3_close((sqlite3 *)data);
    return 0;
}

double compute_granted_credit(WORKUNIT& wu, vector<RESULT>& results) {
    return median_mean_credit(wu, results);
}
