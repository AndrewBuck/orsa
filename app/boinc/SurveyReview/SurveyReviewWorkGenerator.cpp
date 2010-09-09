#include "SurveyReviewWorkGenerator.h"

#include "boinc_db.h"
#include "backend_lib.h"
#include "lib/util.h"

#include <orsa/datetime.h>
#include <orsaSolarSystem/datetime.h>

#include <osg/ref_ptr>

bool writeInputTemplate(const std::string  & fileName,
                        const unsigned int   min_quorum,
                        const unsigned int   target_nresults,
                        const orsa::Time   & delay_bound,
                        const double & rsc_fpops_est,
                        const double & rsc_fpops_bound,
                        const double & rsc_memory_bound,
                        const double & rsc_disk_bound) {
    FILE * fp = fopen(fileName.c_str(),"w");
    if (!fp) {
        ORSA_ERROR("problems");
        return false;
    }
    gmp_fprintf(fp,
                "<file_info>\n"
                "  <number>0</number>\n"
                "</file_info>\n"
                "<file_info>\n"
                "  <number>1</number>\n"
                "</file_info>\n"
                "<file_info>\n"
                "  <number>2</number>\n"
                "</file_info>\n"
                "<file_info>\n"
                "  <number>3</number>\n"
                "</file_info>\n"
                "<file_info>\n"
                "  <number>4</number>\n"
                "</file_info>\n"
                "<file_info>\n"
                "  <number>5</number>\n"
                "</file_info>\n"
                "<workunit>\n"
                "  <file_ref>\n"
                "    <file_number>0</file_number>\n"
                "    <open_name>field.dat</open_name>\n"
                "  </file_ref>\n"
                "  <file_ref>\n"
                "    <file_number>1</file_number>\n"
                "    <open_name>fieldTime.dat</open_name>\n"
                "  </file_ref>\n"
                "  <file_ref>\n"
                "    <file_number>2</file_number>\n"
                "    <open_name>fit.dat</open_name>\n"
                "  </file_ref>\n"
                "  <file_ref>\n"
                "    <file_number>3</file_number>\n"
                "    <open_name>grid.dat</open_name>\n"
                "  </file_ref>\n"
                "  <file_ref>\n"
                "    <file_number>4</file_number>\n"
                "    <open_name>obscode.dat</open_name>\n"
                "  </file_ref>\n"
                "  <file_ref>\n"
                "    <file_number>5</file_number>\n"
                "    <open_name>randomSeed.dat</open_name>\n"
                "  </file_ref>\n"
                "  <min_quorum>%i</min_quorum>\n"
                "  <target_nresults>%i</target_nresults>\n"
                "  <max_error_results>6</max_error_results>\n"
                "  <max_total_results>12</max_total_results>\n"
                "  <max_success_results>6</max_success_results>\n"
                "  <delay_bound>%f</delay_bound>\n"
                "  <rsc_fpops_est>%f</rsc_fpops_est>\n"
                "  <rsc_fpops_bound>%f</rsc_fpops_bound>\n"
                "  <rsc_memory_bound>%f</rsc_memory_bound>\n"
                "  <rsc_disk_bound>%f</rsc_disk_bound>\n"
                "</workunit>\n",
                min_quorum,
                target_nresults,
                orsa::FromUnits(delay_bound.get_d(),orsa::Unit::SECOND,-1),
                rsc_fpops_est,
                rsc_fpops_bound,
                rsc_memory_bound,
                rsc_disk_bound);
    fclose(fp);
    return true;
}

bool writeOutputTemplate(const std::string          & fileName,
                         const OutputTemplateVector & outputTemplateVector) {
    FILE * fp = fopen(fileName.c_str(),"w");
    if (!fp) {
        ORSA_ERROR("problems");
        return false;
    }
    {
        unsigned int uploadCounter=0;
        for (unsigned int k=0; k<outputTemplateVector.size(); ++k) {
            if (outputTemplateVector[k].upload.getRef()) {
                gmp_fprintf(fp,
                            "<file_info>\n"
                            "  <name><OUTFILE_%i/></name>\n"
                            "  <generated_locally/>\n"
                            "  <upload_when_present/>\n"
                            "  <max_nbytes>%i</max_nbytes>\n"
                            "  <url><UPLOAD_URL/></url>\n"
                            "</file_info>\n",
                            uploadCounter,
                            outputTemplateVector[k].fileSize.getRef());
                ++uploadCounter;
            }
        }
    }
    gmp_fprintf(fp,"<result>\n");
    {
        unsigned int uploadCounter=0;
        for (unsigned int k=0; k<outputTemplateVector.size(); ++k) {
            if (outputTemplateVector[k].upload.getRef()) {
                gmp_fprintf(fp,
                            "  <file_ref>\n"
                            "    <file_name><OUTFILE_%i/></file_name>\n"
                            "    <open_name>%s</open_name>\n",
                            uploadCounter,
                            outputTemplateVector[k].fileName.getRef().c_str());
                if (outputTemplateVector[k].optional.getRef()) {
                    gmp_fprintf(fp,
                                "    <optional/>\n");
                }
                gmp_fprintf(fp,
                            "  </file_ref>\n");
                ++uploadCounter;
            }
        }
    }
    gmp_fprintf(fp,"</result>\n");
    fclose(fp);
    return true;
}

int main(int argc, char ** argv) {
    
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " baseName" << std::endl;
        exit(0);
    }
    
    const std::string baseName = argv[1];
    
    OutputTemplateVector outputTemplateVector;
    //
    outputTemplateVector.push_back(OutputTemplateEntry("results.db", 100*1024*1024, true, true));
    
    // tentative values
    const double flops_est = 2e12;
    
    const unsigned int min_quorum       = 2;
    const unsigned int target_nresults  = 2;
    const double rsc_fpops_est          =    flops_est;
    const double rsc_fpops_bound        = 32*flops_est;
    const orsa::Time   delay_bound      = orsa::Time(3,0,0,0,0);
    const double rsc_memory_bound       = 134217728;
    //
    double rsc_disk_bound = 0;
    for (unsigned int k=0; k<outputTemplateVector.size(); ++k) {
        rsc_disk_bound += outputTemplateVector[k].fileSize.getRef();
    }
    
    /* const std::string inTemplateName = 
       inputTemplateFileName(baseName,
       min_quorum,
       target_nresults,
       delay_bound,
       rsc_fpops_est,
       rsc_fpops_bound,
       rsc_memory_bound,
       rsc_disk_bound);    
    */
    //
    char inTemplateName[1024];
    sprintf(inTemplateName,"templates/wu.SR.input.%s.xml",baseName.c_str());
    
    if (!writeInputTemplate(inTemplateName,
                            min_quorum,
                            target_nresults,
                            delay_bound,
                            rsc_fpops_est,
                            rsc_fpops_bound,
                            rsc_memory_bound,
                            rsc_disk_bound)) {
        exit(0);
    }
    
    char outTemplateName[1024];
    sprintf(outTemplateName,"templates/wu.SR.output.%s.xml",baseName.c_str());
    
    if (!writeOutputTemplate(outTemplateName,
                             outputTemplateVector)) {
        exit(0);
    }
  
    const orsa::Time timeNow = orsaSolarSystem::now();
    int y,m,d,H,M,S,ms;
    orsaSolarSystem::gregorDay(timeNow,
                               y,
                               m,
                               d,
                               H,
                               M,
                               S,
                               ms);
    
    // unique!
    char wuName[1024];
    //
    {
        gmp_snprintf(wuName,
                     1024,
                     "%s_%ld_%i",
                     baseName.c_str(),
                     time(0),
                     getpid());
    }
    //
    ORSA_DEBUG("WU name: [%s]",wuName);
    
    int retVal;
  
    SCHED_CONFIG config;
  
    config.parse_file();
  
    retVal = boinc_db.open(config.db_name, config.db_host, config.db_user, config.db_passwd);
    //
    if (retVal != 0) {
        ORSA_DEBUG("problems with boinc_db.open(...)");
        exit(1);
    }
  
    DB_APP app;
    //
    retVal = app.lookup("where name='SurveyReview'");
    //
    if (retVal != 0) {
        ORSA_DEBUG("problems with app.lookup(...)");
        exit(1);
    }
  
    char cmd[1024];
    char path[1024];
  
    char copyFieldDat[1024];
    //
    gmp_snprintf(copyFieldDat,1024,"field.dat.%s",wuName);
    config.download_path(copyFieldDat,path);
    gmp_snprintf(cmd,1024,"cp -f field.dat %s",path);
    if (system(cmd)) {
        ORSA_DEBUG("problems with system call: [%s]",cmd);
        exit(1);
    }   
    
    char copyFieldTimeDat[1024];
    //
    gmp_snprintf(copyFieldTimeDat,1024,"fieldTime.dat.%s",wuName);
    config.download_path(copyFieldTimeDat,path);
    gmp_snprintf(cmd,1024,"cp -f fieldTime.dat %s",path);
    if (system(cmd)) {
        ORSA_DEBUG("problems with system call: [%s]",cmd);
        exit(1);
    }   
    
    char copyFitDat[1024];
    //
    gmp_snprintf(copyFitDat,1024,"fit.dat.%s",wuName);
    config.download_path(copyFitDat,path);
    gmp_snprintf(cmd,1024,"cp -f fit.dat %s",path);
    if (system(cmd)) {
        ORSA_DEBUG("problems with system call: [%s]",cmd);
        exit(1);
    }   
    
    char copyGridDat[1024];
    //
    gmp_snprintf(copyGridDat,1024,"grid.dat.%s",wuName);
    config.download_path(copyGridDat,path);
    gmp_snprintf(cmd,1024,"cp -f grid.dat %s",path);
    if (system(cmd)) {
        ORSA_DEBUG("problems with system call: [%s]",cmd);
        exit(1);
    }   
    
    char copyObscodeDat[1024];
    //
    gmp_snprintf(copyObscodeDat,1024,"obscode.dat.%s",wuName);
    config.download_path(copyObscodeDat,path);
    gmp_snprintf(cmd,1024,"cp -f obscode.dat %s",path);
    if (system(cmd)) {
        ORSA_DEBUG("problems with system call: [%s]",cmd);
        exit(1);
    }   
    
    char copyRandomSeedDat[1024];
    //
    gmp_snprintf(copyRandomSeedDat,1024,"randomSeed.dat.%s",wuName);
    config.download_path(copyRandomSeedDat,path);
    gmp_snprintf(cmd,1024,"cp -f randomSeed.dat %s",path);
    if (system(cmd)) {
        ORSA_DEBUG("problems with system call: [%s]",cmd);
        exit(1);
    }   
    
    const char * infiles[] = 
        { copyFieldDat,
          copyFieldTimeDat,
          copyFitDat,
          copyGridDat,
          copyObscodeDat,
          copyRandomSeedDat 
        };
    
    char * wu_template;
    read_file_malloc(inTemplateName, wu_template);
    
    DB_WORKUNIT wu;
    wu.clear(); // zeroes all fields
    wu.appid = app.id;
    strcpy(wu.name, wuName);
  
    create_work(wu,
                wu_template,
                outTemplateName,
                outTemplateName,
                infiles,
                6, // infiles size
                config);
}
