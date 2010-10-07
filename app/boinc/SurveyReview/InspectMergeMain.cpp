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

// Nsub is the number of sub-bins in each x,y bin -- to account for missing zero entries
// i.e. if writing a,e then Nsub=Ni*Nnode*Nperi*Nm=18*12*12*12...
void writeOutputFile(const std::string & filename,
                     PlotStats * plotStats,
                     const PlotStats::LinearVar * var_x,
                     const PlotStats::LinearVar * var_y,
                     const unsigned int Nsub) {
    
    FILE * fp = fopen(filename.c_str(),"w");    
    
    std::vector<double> xVector;
    xVector.resize(2);
    
    for (unsigned j=0; j<var_x->size(); ++j) {
        for (unsigned k=0; k<var_y->size(); ++k) {
            
            xVector[0] = var_x->start+var_x->incr*(j+0.5);
            xVector[1] = var_y->start+var_y->incr*(k+0.5);
            
            std::vector<size_t> binVector;
            if (plotStats->bin(binVector,xVector)) {
                const PlotStatsElement * e = plotStats->stats(plotStats->index(binVector));
                if (e) {
                    gmp_fprintf(fp,
                                "%8g %8g %8.6f %8.6f %6Zi\n",
                                xVector[0],
                                xVector[1],
                                e->average(),
                                e->average()*e->entries().get_d()/Nsub, // divide by Nsub to account for missing zero entries
                                e->entries().get_mpz_t());
#warning the two should be the same for a complete analysis, because over several years, each a,e bin has been observed at least once!
                }
            }
        }
    }
    
    fclose(fp);
}

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


        // vars
        const double a_min  = 1.90;
        const double a_max  = 2.30;
        const double a_step = 0.05;
        osg::ref_ptr<PlotStats::LinearVar> var_a = new PlotStats::LinearVar(a_min,a_max,a_step);
        //
        const double e_min  = 0.00;
        const double e_max  = 1.00;
        const double e_step = 0.05;
        osg::ref_ptr<PlotStats::LinearVar> var_e = new PlotStats::LinearVar(e_min,e_max,e_step);
        //
        const double i_min  =  0.00;
        const double i_max  = 90.00;
        const double i_step =  5.00;
        osg::ref_ptr<PlotStats::LinearVar> var_i = new PlotStats::LinearVar(i_min,i_max,i_step);
        //
        const double L_min  =   0.00;
        const double L_max  = 360.00;
        const double L_step =  30.00;
        osg::ref_ptr<PlotStats::LinearVar> var_L = new PlotStats::LinearVar(L_min,L_max,L_step);
        
        // a,e
        std::vector< osg::ref_ptr<PlotStats::Var> > varDefinition_ae;
        varDefinition_ae.push_back(var_a.get());
        varDefinition_ae.push_back(var_e.get());
        //
        osg::ref_ptr<PlotStats> plotStats_ae_NEO_H18 = new PlotStats(varDefinition_ae);
        osg::ref_ptr<PlotStats> plotStats_ae_PHO_H18 = new PlotStats(varDefinition_ae);
        // a,i
        std::vector< osg::ref_ptr<PlotStats::Var> > varDefinition_ai;
        varDefinition_ai.push_back(var_a.get());
        varDefinition_ai.push_back(var_i.get());
        //
        osg::ref_ptr<PlotStats> plotStats_ai_NEO_H18 = new PlotStats(varDefinition_ai);
        osg::ref_ptr<PlotStats> plotStats_ai_PHO_H18 = new PlotStats(varDefinition_ai);
        // a,L
        std::vector< osg::ref_ptr<PlotStats::Var> > varDefinition_aL;
        varDefinition_aL.push_back(var_a.get());
        varDefinition_aL.push_back(var_L.get());
        osg::ref_ptr<PlotStats> plotStats_aL_NEO_H18 = new PlotStats(varDefinition_aL);
        osg::ref_ptr<PlotStats> plotStats_aL_PHO_H18 = new PlotStats(varDefinition_aL);
        
        // size=2 should work for most cases
        std::vector<double> xVector; xVector.resize(2);
        // xVector.resize(varDefinition.size());
        
        for (int row=1; row<=nrows; ++row) {
            
            if (row%10000==0) ORSA_DEBUG("progress: %i/%i (%6.3f\%)",row,nrows,100*(double)row/(double)nrows);
            
            const int z_a_min          = atoi(sql_result[row*ncols+0]);
            const int z_a_max          = atoi(sql_result[row*ncols+1]);
            const int z_e_min          = atoi(sql_result[row*ncols+2]);
            const int z_e_max          = atoi(sql_result[row*ncols+3]);
            const int z_i_min          = atoi(sql_result[row*ncols+4]);
            const int z_i_max          = atoi(sql_result[row*ncols+5]);
            const int z_node_min       = atoi(sql_result[row*ncols+6]);
            const int z_node_max       = atoi(sql_result[row*ncols+7]);
            const int z_peri_min       = atoi(sql_result[row*ncols+8]);
            const int z_peri_max       = atoi(sql_result[row*ncols+9]);
            const int z_M_min          = atoi(sql_result[row*ncols+10]);
            const int z_M_max          = atoi(sql_result[row*ncols+11]);
            //
            const int z_H              = atoi(sql_result[row*ncols+12]);
            //
            const int N_NEO            = atoi(sql_result[row*ncols+13]);
            const int N_PHO            = atoi(sql_result[row*ncols+14]);
            const int NEO_in_field     = atoi(sql_result[row*ncols+15]);
            const int PHO_in_field     = atoi(sql_result[row*ncols+16]);
            //
            const double eta_NEO       = atof(sql_result[row*ncols+17]);
            const double sigma_eta_NEO = atof(sql_result[row*ncols+18]);
            const double eta_PHO       = atof(sql_result[row*ncols+19]);
            const double sigma_eta_PHO = atof(sql_result[row*ncols+20]);
            
            // local center bin
            const double center_a = 0.5*(z_a_max+z_a_min)*grain_a_AU;
            const double center_e = 0.5*(z_e_max+z_e_min)*grain_e;
            const double center_i = 0.5*(z_i_max+z_i_min)*grain_i_DEG;
            
            if (z_H == 180) {
                xVector[0] = center_a;
                xVector[1] = center_e;
                plotStats_ae_NEO_H18->insert(xVector, eta_NEO, 1.0);
                plotStats_ae_PHO_H18->insert(xVector, eta_PHO, 1.0);
                //
                xVector[0] = center_a;
                xVector[1] = center_i;
                plotStats_ai_NEO_H18->insert(xVector, eta_NEO, 1.0);
                plotStats_ai_PHO_H18->insert(xVector, eta_PHO, 1.0);
            }
            
            // L
            std::vector<int> z_L_vec;
            {
                // assume step is the same for node,peri,M
                // #warning this is grain-size dependent!
                z_L_vec.push_back(((z_node_min+z_peri_min+z_M_min)%360000)/30000);
                z_L_vec.push_back(((z_node_max+z_peri_min+z_M_min)%360000)/30000);
                z_L_vec.push_back(((z_node_max+z_peri_max+z_M_min)%360000)/30000);
            }
            //
            for (unsigned int l=0; l<z_L_vec.size(); ++l) {
                const int index_L = z_L_vec[l]; // 0-11
                
                // assuming grain size for node=peri=M=L
                const double center_L = (0.5+index_L)*30.0;
                
                if (z_H == 180) {
                    xVector[0] = center_a;
                    xVector[1] = center_L;
                    plotStats_aL_NEO_H18->insert(xVector, eta_NEO, 1.0);
                    plotStats_aL_PHO_H18->insert(xVector, eta_PHO, 1.0);
                }
            }
            
            
        }
        
        // numbers for Nsub
#warning these can change at each run
        // no sub_a since it is always plotted
        const unsigned int sub_e = 20;
        const unsigned int sub_i = 18;
        const unsigned int sub_node = 12;
        const unsigned int sub_peri = 12;
        const unsigned int sub_M    = 12;
        
        // time to write output files
        char filename[1024];
        //
        sprintf(filename,"%s_inspect_ae_NEO_H18.dat",argv[1]);
        writeOutputFile(filename, plotStats_ae_NEO_H18, var_a, var_e, sub_i*sub_node*sub_peri*sub_M);
        //
        sprintf(filename,"%s_inspect_ae_PHO_H18.dat",argv[1]);
        writeOutputFile(filename, plotStats_ae_PHO_H18, var_a, var_e, sub_i*sub_node*sub_peri*sub_M);
        //
        sprintf(filename,"%s_inspect_ai_NEO_H18.dat",argv[1]);
        writeOutputFile(filename, plotStats_ai_NEO_H18, var_a, var_i, sub_e*sub_node*sub_peri*sub_M);
        //
        sprintf(filename,"%s_inspect_ai_PHO_H18.dat",argv[1]);
        writeOutputFile(filename, plotStats_ai_PHO_H18, var_a, var_i, sub_e*sub_node*sub_peri*sub_M);
        //
        sprintf(filename,"%s_inspect_aL_NEO_H18.dat",argv[1]);
        writeOutputFile(filename, plotStats_aL_NEO_H18, var_a, var_L, sub_e*sub_i*sub_node*sub_peri*3);
        //
        sprintf(filename,"%s_inspect_aL_PHO_H18.dat",argv[1]);
        writeOutputFile(filename, plotStats_aL_PHO_H18, var_a, var_L, sub_e*sub_i*sub_node*sub_peri*3);
        
        
        sqlite3_free_table(sql_result);
    }
    
    sqlite3_close(db);
    
    return 0;
}

