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
    
    if (0) {
        // debug output
        {
            std::list<unsigned int>::const_iterator it = number_observed.begin();
            while (it != number_observed.end()) {
                ORSA_DEBUG("observed: %6i",(*it));
                ++it;
            }
        }
        {
            std::list<std::string>::const_iterator it = designation_observed.begin();
            while (it != designation_observed.end()) {
                ORSA_DEBUG("observed: %10s",(*it).c_str());
                ++it;
            }
        }
    }
    
    osg::ref_ptr<orsaInputOutput::MPCAsteroidFile> orbitFile =
        new orsaInputOutput::MPCAsteroidFile;
    orbitFile->select_NEO = true;
    orbitFile->setFileName("MPCORB.DAT.gz");
    orbitFile->read();
    ORSA_DEBUG("selected orbits: %i",orbitFile->_data.size());
    
    {
        // now the real work
        
        
        
#warning remember to propagate the mean anomaly
        
        
        char **sql_result;
        int nrows, ncols;
        char sql_line[1024];

        orsa::Cache<orsaInputOutput::MPCAsteroidDataElement> found;
        orsaInputOutput::MPCAsteroidData::const_iterator it_orb =
            orbitFile->_data.begin();
        while (it_orb != orbitFile->_data.end()) {
            if ((*it_orb).number.isSet()) {
                std::list<unsigned int>::const_iterator it_num = number_observed.begin();
                while (it_num != number_observed.end()) {
                    if ((*it_num) == (*it_orb).number.getRef()) {
                        found=(*it_orb);
                        break;
                    }
                    ++it_num;
                }
            } else if ((*it_orb).designation.isSet()) {
                std::list<std::string>::const_iterator it_des = designation_observed.begin();
                while (it_des != designation_observed.end()) {
                    if ((*it_des) == (*it_orb).designation.getRef()) {
                        found=(*it_orb);
                        break;
                    }
                    ++it_des;
                }
            } else {
                ORSA_DEBUG("neither one is set?");
            }
            
            if (found.isSet()) break;
            ++it_orb;
        }
        
        if (found.isSet()) {

#warning should have these delta values somewhere else, once and for all
            const int z_a_delta = lrint(0.05/grain_a_AU);
            const int z_e_delta = lrint(0.05/grain_e);
            const int z_i_delta = lrint(5.00/grain_i_DEG);
            const int z_H_delta = lrint(1.00/grain_H);
            
            const int z_a_min = (  lrint(orsa::FromUnits(found.getRef().orbit.getRef().a,orsa::Unit::AU,-1)/grain_a_AU)/z_a_delta)*z_a_delta;
            const int z_a_max = (1+lrint(orsa::FromUnits(found.getRef().orbit.getRef().a,orsa::Unit::AU,-1)/grain_a_AU)/z_a_delta)*z_a_delta;
            const int z_e_min = (  lrint(found.getRef().orbit.getRef().e/grain_e)/z_e_delta)*z_e_delta;
            const int z_e_max = (1+lrint(found.getRef().orbit.getRef().e/grain_e)/z_e_delta)*z_e_delta;
            const int z_i_min = (  lrint(found.getRef().orbit.getRef().i*orsa::radToDeg()/grain_i_DEG)/z_i_delta)*z_i_delta;
            const int z_i_max = (1+lrint(found.getRef().orbit.getRef().i*orsa::radToDeg()/grain_i_DEG)/z_i_delta)*z_i_delta;
            
#warning check this "+1"
            const int z_H     = (1+lrint(found.getRef().H.getRef()/grain_H)/z_H_delta)*z_H_delta;
            
            
            if (found.getRef().number.isSet()) {
                ORSA_DEBUG("found: [%i] z_a: [%i,%i] z_e: [%i,%i] z_i: [%i,%i] z_H: %i",
                           found.getRef().number.getRef(),
                           z_a_min,
                           z_a_max,
                           z_e_min,
                           z_e_max,
                           z_i_min,
                           z_i_max,
                           z_H);
            } else if (found.getRef().designation.isSet()) {
                ORSA_DEBUG("found: [%s]",found.getRef().designation.getRef().c_str());
            } else {
                ORSA_DEBUG("neither one is set?");
            }
            
            /* sprintf(sql_line,"SELECT * FROM grid");
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
            */
        }
        
    }
    
    sqlite3_close(db);
    
    return 0;
}

