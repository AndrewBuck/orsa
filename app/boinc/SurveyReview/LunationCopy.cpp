#include "skycoverage.h"
#include "eta.h"
#include "fit.h"
#include <dirent.h>

void run(char * cmd) {
    ORSA_DEBUG("system call: [%s]",cmd);
    system(cmd);
}

int main(int argc, char ** argv) {
    
    if (argc < 3) {
        printf("Usage: %s <lunation-dir> <sky_coverage_file(s)>\n",argv[0]);
        exit(0);
    }
    
    const std::string baseDir = argv[1];
    
    {
        DIR * dir = opendir(baseDir.c_str());
        if (dir==0) {
            int errnum = errno;
            ORSA_DEBUG("problem with [%s]: %s",baseDir.c_str(),strerror(errnum));
            exit(0);
        } else {
            closedir(dir);
        }
    }
    
    // first, read obscode and epoch from filename
    osg::ref_ptr<orsaInputOutput::MPCObsCodeFile> obsCodeFile = new orsaInputOutput::MPCObsCodeFile;
    obsCodeFile->setFileName("obscode.dat");
    obsCodeFile->read();
    
    char cmd[1024];
    
    for (int arg=2; arg<argc; ++arg) {

        FILE * fp_sky = fopen(argv[arg],"r");
        if (fp_sky!=0) {
            fclose(fp_sky);
        } else {
            ORSA_DEBUG("cannot open file [%s], skipping...",argv[arg]);
            continue;
        }
     
        // longBasename includes the whole path
        const size_t found_last_slash = std::string(argv[arg]).find_last_of("//");
        const size_t found_dot = std::string(argv[arg]).find(".",(found_last_slash==std::string::npos?0:found_last_slash+1));
        std::string longBasename;
        longBasename.assign(argv[arg],0,found_dot);
        char allEtaFilename[1024];
        sprintf(allEtaFilename,"%s.allEta.dat",longBasename.c_str());
        FILE * fp_allEta = fopen(allEtaFilename,"r");
        if (fp_allEta!=0) {
            fclose(fp_allEta);
        } else {
            ORSA_DEBUG("cannot open file [%s], skipping...",allEtaFilename);
            continue;
        }
        
        std::string obsCode;
        orsa::Time epoch;
        int year;
        int dayOfYear;
        //
        SkyCoverage::processFilename(argv[arg],
                                     obsCodeFile.get(),
                                     obsCode,
                                     epoch,
                                     year,
                                     dayOfYear);
        orsa::Time begin, end;
        int lunationID;
        orsaSolarSystem::ApproximatedLunation(epoch,
                                              begin,
                                              end,
                                              lunationID);
        char lunationDir[1024];
        sprintf(lunationDir,"%s/L%03i",baseDir.c_str(),lunationID);
        sprintf(cmd,"mkdir -p %s",lunationDir);
        run(cmd);
        sprintf(cmd,"cp -a %s %s %s",argv[arg],allEtaFilename,lunationDir);
        run(cmd);
    }
    
    return 0;
}
