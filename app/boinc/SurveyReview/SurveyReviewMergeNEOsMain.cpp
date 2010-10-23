#include "grain.h"

#include <orsa/debug.h>
#include <orsa/print.h>
#include <orsa/statistic.h>

#include <orsaSolarSystem/datetime.h>
#include <orsaSolarSystem/print.h>

#include <orsa/util.h>
#include <orsaInputOutput/MPC_asteroid.h>

#include "grain.h"
#include "fit.h"
#include "skycoverage.h"
#include "eta.h"

// SQLite3
#include "sqlite3.h"

int main(int argc, char ** argv) {
    
    if (argc < 3) {
        ORSA_DEBUG("Usage: %s <sqlite-merged-db> <allEta-file(s)>",argv[0]);
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

    // read allEta files and extract ID (number or designation) of all observed objects (not only NEOs)
    std::list<unsigned int>      number_observed;
    std::list<std::string>  designation_observed;
    {
        const unsigned int numFiles=argc-2;
        std::vector <EfficiencyData> etaData;
        for (unsigned int fileID=0; fileID<numFiles; ++fileID) {
            etaData.clear();
            readEfficiencyDataFile(etaData,argv[fileID+2]);
            ORSA_DEBUG("file: [%s]   etaData.size(): %i",argv[fileID+2],etaData.size());
            for (unsigned int k=0; k<etaData.size(); ++k) {
                if (etaData[k].number.isSet()) {
                    number_observed.push_back(etaData[k].number.getRef());
                } else if (etaData[k].designation.isSet()) {
                    designation_observed.push_back(etaData[k].designation.getRef());
                } else {
                    ORSA_DEBUG("problem: neither number or designation is set");
                }
            }
        }
        number_observed.sort();
        number_observed.unique();
        designation_observed.sort();
        designation_observed.unique();
    }

    osg::ref_ptr<orsaInputOutput::MPCAsteroidFile> orbitFile =
        new orsaInputOutput::MPCAsteroidFile;
    orbitFile->select_NEO = true;
    orbitFile->setFileName("MPCORB.DAT.gz");
    orbitFile->read();
    ORSA_DEBUG("selected orbits: %i",orbitFile->_data.size());
    
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
        //
        // loop here...
        //
        sqlite3_free_table(sql_result);
    }
    
    sqlite3_close(db);
    
    return 0;
}

