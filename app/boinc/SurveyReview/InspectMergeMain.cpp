#include "grain.h"

#include <orsa/debug.h>
#include <orsa/print.h>
#include <orsa/statistic.h>

#include <orsaSolarSystem/datetime.h>
#include <orsaSolarSystem/print.h>

#include "grain.h"
#include "fit.h"

// SQLite3
#include "sqlite3.h"

class PlotStatsElement : public orsa::WeightedStatistic<double> {
    
};

class PlotStats : public BinStats<PlotStatsElement> {
public:
    PlotStats(const std::vector< osg::ref_ptr<Var> > & varDefinition) :
        BinStats<PlotStatsElement>(varDefinition) { }
public:
    bool insert(const std::vector<double> & xVector,
                const double & val,
                const double & sigma) {
        if (xVector.size() != var.size()) {
            ORSA_DEBUG("dimension mismatch");
            return false;
        }
        std::vector<size_t> binVector;
        if (!bin(binVector,xVector)) {
            return false;
        }
        const mpz_class idx = index(binVector);
        if (data[idx].get()==0) {
            // lazy allocation
            data[idx] = new PlotStatsElement;
        }
        data[idx]->insert(val,orsa::square(1.0/sigma));
        // ORSA_DEBUG("xV: %g %g bV: %i %i",xVector[0],xVector[1],binVector[0],binVector[1]);
        return true;         
    }
};

int main(int argc, char ** argv) {
    
    if (argc != 2) {
        ORSA_DEBUG("Usage: %s <sqlite-merged-db>",argv[0]);
        exit(0);
    }
    
    orsa::Debug::instance()->initTimer();
    
    ORSA_DEBUG("process ID: %i",getpid());
    
    // needed to work with SQLite database
    sqlite3     * db;
    char        * zErr;
    int           rc;
    std::string   sql;
    
    {
        // open database
        // rc = sqlite3_open(argv[1],&db);
        rc = sqlite3_open_v2(argv[1],&db,SQLITE_OPEN_READONLY,NULL);
        //
        if (rc) {
            fprintf(stderr,"Can't open db: %s\n",sqlite3_errmsg(db));
            sqlite3_close(db);
            exit(0);
        }
    }
    
    {
        // integrity_check
        
        int nrows, ncols;
        char * * sql_result;
        do {
            rc = sqlite3_get_table(db,"pragma integrity_check",&sql_result,&nrows,&ncols,&zErr);
            if (rc==SQLITE_BUSY) {
                ORSA_DEBUG("database busy, retrying...");
                usleep(100000);
            }
        } while (rc==SQLITE_BUSY);
        if (rc != SQLITE_OK) {
            if (zErr != NULL) {
                ORSA_DEBUG("SQL error: %s\n",zErr); 
                sqlite3_free(zErr);
                sqlite3_close(db);
                exit(0);
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
            exit(0); 
        }
        sqlite3_free_table(sql_result);
    }
    
    {
        // now the real work
        
        char **sql_result;
        int nrows, ncols;
        char sql_line[1024];
        sprintf(sql_line,"SELECT * FROM grid");
        do {
            rc = sqlite3_get_table(db,sql_line,&sql_result,&nrows,&ncols,&zErr);
            if (rc==SQLITE_BUSY) {
                ORSA_DEBUG("database busy, retrying...");
                usleep(100000);
            }
        } while (rc==SQLITE_BUSY);
        if (rc != SQLITE_OK) {
            if (zErr != NULL) {
                ORSA_DEBUG("SQL error: %s\n",zErr);
                sqlite3_free(zErr);
                sqlite3_close(db);
                exit(0);
            }
        }

        // a,e
        const double x_step = 0.05;
        const double x_min  = 1.90;
        const double x_max  = 2.30;
        //
        const double y_step =  0.05;
        const double y_min  =  0.00;
        const double y_max  =  1.00;
        
        std::vector< osg::ref_ptr<PlotStats::Var> > varDefinition;
        //
        osg::ref_ptr<PlotStats::LinearVar> var_x = new PlotStats::LinearVar(x_min,x_max,x_step);
        varDefinition.push_back(var_x.get());
        //
        osg::ref_ptr<PlotStats::LinearVar> var_y = new PlotStats::LinearVar(y_min,y_max,y_step);
        varDefinition.push_back(var_y.get());
        
        osg::ref_ptr<PlotStats> plotStats = 
            new PlotStats(varDefinition);

        std::vector<double> xVector;
        xVector.resize(varDefinition.size());
        
        for (int row=1; row<=nrows; ++row) {
            
            if (row%10000==0) ORSA_DEBUG("progress: %i/%i (%6.3f\%)",row,nrows,100*(double)row/(double)nrows);
            
            const int z_a_min          = atoi(sql_result[row*ncols+0]);
            const int z_a_max          = atoi(sql_result[row*ncols+1]);
            const int z_e_min          = atoi(sql_result[row*ncols+2]);
            const int z_e_max          = atoi(sql_result[row*ncols+3]);
            /* const int z_i_min          = atoi(sql_result[row*ncols+4]);
               const int z_i_max          = atoi(sql_result[row*ncols+5]);
               const int z_node_min       = atoi(sql_result[row*ncols+6]);
               const int z_node_max       = atoi(sql_result[row*ncols+7]);
               const int z_peri_min       = atoi(sql_result[row*ncols+8]);
               const int z_peri_max       = atoi(sql_result[row*ncols+9]);
               const int z_M_min          = atoi(sql_result[row*ncols+10]);
               const int z_M_max          = atoi(sql_result[row*ncols+11]);
            */
            //
            const int z_H              = atoi(sql_result[row*ncols+12]);
            //
            /* const int N_NEO            = atoi(sql_result[row*ncols+13]);
               const int N_PHO            = atoi(sql_result[row*ncols+14]);
               const int NEO_in_field     = atoi(sql_result[row*ncols+15]);
               const int PHO_in_field     = atoi(sql_result[row*ncols+16]);
            */
            const double eta_NEO       = atof(sql_result[row*ncols+17]);
            const double sigma_eta_NEO = atof(sql_result[row*ncols+18]);
            const double eta_PHO       = atof(sql_result[row*ncols+19]);
            const double sigma_eta_PHO = atof(sql_result[row*ncols+20]);


            // local center bin
            const double center_a    = 0.5*(z_a_max+z_a_min)*grain_a_AU;
            const double center_e    = 0.5*(z_e_max+z_e_min)*grain_e;
            // const double center_i    = 0.5*(z_i_max+z_i_min)*grain_i_DEG;
            // assuming grain size for node=peri=M=L
            // const double center_L = 0.5*(z_L_max+z_L_min)*grain_node_DEG;

#warning FIX z_H HERE
            if (z_H == 180) {
                // a,e
                xVector[0] = center_a;
                xVector[1] = center_e;
                
#warning FIX eta HERE
                plotStats->insert(xVector, eta_NEO, 1.0);
            }
        }
        
        for (unsigned j=0; j<var_x->size(); ++j) {
            for (unsigned k=0; k<var_y->size(); ++k) {
                // const unsigned int mesh_id = j*var_y->size()+k;
                xVector[0] = x_min+x_step*(j+0.5);
                xVector[1] = y_min+y_step*(k+0.5);
                std::vector<size_t> binVector;
                if (plotStats->bin(binVector,xVector)) {
                    const PlotStatsElement * e =  plotStats->stats(plotStats->index(binVector));
                    if (e) {
                        // printf("%g %g %g\n",xVector[0],xVector[1],e->average());
                        gmp_printf("%g %g %8.6f %8.6f %6Zi\n",
                                   xVector[0],
                                   xVector[1],
                                   e->average(),
                                   e->average()*e->entries().get_d()/(18*12*12*12), // number of sub-bins in each a,e bin -- to account for missing zero entries
                                   e->entries().get_mpz_t());
#warning the two should be the same for a complete analysis, because over several years, each a,e bin has been observed at least once!
                    }
                }
            }
        }
        
        sqlite3_free_table(sql_result);
    }
    
    sqlite3_close(db);
    
    return 0;
}

