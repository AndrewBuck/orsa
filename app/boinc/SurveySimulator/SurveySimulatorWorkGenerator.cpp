#include "SurveySimulatorWorkGenerator.h"

#include "boinc_db.h"
#include "backend_lib.h"
#include "lib/util.h"

#include "parameter.h"
#include "telescope.h"

#include <osg/ref_ptr>

std::string inputTemplateFileName(const std::string  & baseName,
				  const unsigned int   ,
				  const unsigned int   ,
				  const orsa::Time   & delay_bound,
				  const orsa::Double & ,
				  const orsa::Double & ,
				  const orsa::Double & rsc_memory_bound,
				  const orsa::Double & rsc_disk_bound) {
  char fileName[1024];
  gmp_snprintf(fileName,
	       1024,
	       "templates/wu.SS.input.%s_%.Ff_%.Ff_%.Ff.xml",  
	       baseName.c_str(),
	       orsa::Double(delay_bound.getMuSec()/1000000.0).get_mpf_t(),
	       rsc_memory_bound.get_mpf_t(),
	       rsc_disk_bound.get_mpf_t());
  return fileName;
}

bool writeInputTemplate(const std::string  & fileName,
			const unsigned int   min_quorum,
			const unsigned int   target_nresults,
			const orsa::Time   & delay_bound,
			const orsa::Double & rsc_fpops_est,
			const orsa::Double & rsc_fpops_bound,
			const orsa::Double & rsc_memory_bound,
			const orsa::Double & rsc_disk_bound) {
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
	      "<workunit>\n"
	      "  <file_ref>\n"
	      "    <file_number>0</file_number>\n"
	      "    <open_name>param.conf</open_name>\n"
	      "  </file_ref>\n"
	      "  <file_ref>\n"
	      "    <file_number>1</file_number>\n"
	      "    <open_name>telescope.conf</open_name>\n"
	      "  </file_ref>\n"
	      "  <file_ref>\n"
	      "    <file_number>2</file_number>\n"
	      "    <open_name>obscode.dat</open_name>\n"
	      "  </file_ref>\n"
	      "  <file_ref>\n"
	      "    <file_number>3</file_number>\n"
	      "    <open_name>realNEO.gen</open_name>\n"
	      "  </file_ref>\n"
	      "  <file_ref>\n"
	      "    <file_number>4</file_number>\n"
	      "    <open_name>syntheticNEO.gen</open_name>\n"
	      "  </file_ref>\n"
	      "  <min_quorum>%i</min_quorum>\n"
	      "  <target_nresults>%i</target_nresults>\n"
	      "  <max_error_results>5</max_error_results>\n"
	      "  <max_total_results>10</max_total_results>\n"
	      "  <max_success_results>10</max_success_results>\n"
	      "  <delay_bound>%Ff</delay_bound>\n"
	      "  <rsc_fpops_est>%Ff</rsc_fpops_est>\n"
	      "  <rsc_fpops_bound>%Ff</rsc_fpops_bound>\n"
	      "  <rsc_memory_bound>%Ff</rsc_memory_bound>\n"
	      "  <rsc_disk_bound>%Ff</rsc_disk_bound>\n"
	      "</workunit>\n",
	      min_quorum,
	      target_nresults,
	      orsa::Double(delay_bound.getMuSec()/mpz_class("1000000")).get_mpf_t(),
	      rsc_fpops_est.get_mpf_t(),
	      rsc_fpops_bound.get_mpf_t(),
	      rsc_memory_bound.get_mpf_t(),
	      rsc_disk_bound.get_mpf_t());
  fclose(fp);
  return true;
}

std::string outputTemplateFileName(const std::string  & baseName,
				   const unsigned int   numRealNEO,
				   const unsigned int   numSyntheticNEO,
				   const mpz_class    & tp_lines) {
  char fileName[1024];
  gmp_snprintf(fileName,
	       1024,
	       "templates/wu.SS.output.%s.%i_%i_%Zi.xml",  
	       baseName.c_str(),
	       numRealNEO,
	       numSyntheticNEO,
	       tp_lines.get_mpz_t());
  return fileName;
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
		    "  <max_nbytes>%Zi</max_nbytes>\n"
		    "  <url><UPLOAD_URL/></url>\n"
		    "</file_info>\n",
		    uploadCounter,
		    outputTemplateVector[k].fileSize.getRef().get_mpz_t());
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
  
  // keep this in sync with other source files in this dir!
  const orsa::Time genEpoch = orsaSolarSystem::J2000();
  
  osg::ref_ptr<ParameterFile> parFile = new ParameterFile;
  parFile->mainRead("param.conf");
  
  TelescopeList telescopeList;
  //
  if (!readTelescope(telescopeList,"telescope.conf")) {
    ORSA_ERROR("cannot read telescope configuration file");
    exit(0);
  }
  
  orsa::Time endSimulation = genEpoch;
  { 
    TelescopeList::const_iterator it = telescopeList.begin();
    while (it != telescopeList.end()) {
      endSimulation = std::max((*it)->tStop,endSimulation);
      ++it;
    }
  }
  
  
  FILE         * fp_gen;
  char           line[1024];
  unsigned int   num;
  int            randomSeed;
  
  int numRealNEO=0;
  //
  fp_gen = boinc_fopen("realNEO.gen","r");
  if (!fp_gen) {
    ORSA_ERROR("cannot open file");
    return 0;
  } 
  while (fgets(line,1024,fp_gen) != 0) {
    if (2 == gmp_sscanf(line,"%d %d",
			&num,
			&randomSeed)) {
      numRealNEO += num;
    }
  }
  fclose(fp_gen);
  
  int numSyntheticNEO=0;
  orsa::Cache<int> firstSyntheticRandomSeed;
  //
  fp_gen = boinc_fopen("syntheticNEO.gen","r");
  if (!fp_gen) {
    ORSA_ERROR("cannot open file");
    return 0;
  } 
  while (fgets(line,1024,fp_gen) != 0) {
    if (2 == gmp_sscanf(line,"%d %d",
			&num,
			&randomSeed)) {
      numSyntheticNEO += num;
      if (!firstSyntheticRandomSeed.isSet()) {
	firstSyntheticRandomSeed = randomSeed;
      }
    }
  }
  fclose(fp_gen);
  
  OutputTemplateVector outputTemplateVector;
  //
  mpz_class    total_tp_lines   = 0;
  orsa::Double total_FOV_DEG_SQ = 0;
  // visits = observing time / recycle time
  mpz_class    total_visits = 0;
  {
    TelescopeList::const_iterator it = telescopeList.begin();
    while (it != telescopeList.end()) {
      
      const mpz_class local_tp_lines = ((*it)->tStop - (*it)->tStart).getMuSec() / (*it)->effectiveDutyCycle.getMuSec();
      //
      total_tp_lines += 1 + local_tp_lines;
      
      total_FOV_DEG_SQ += orsa::int_pow((*it)->FOV_DEG,2);
      
      total_visits += 1 + ((*it)->tStop - (*it)->tStart).getMuSec() / (*it)->recycleTime.getMuSec();
      
      outputTemplateVector.push_back(OutputTemplateEntry((*it)->logFileName(), mpz_class("1024")*local_tp_lines, false, true));
      
      ++it;
    }
  }
  //
  // ORSA_DEBUG("total_visits: %Zi",total_visits.get_mpz_t());
  //
  outputTemplateVector.push_back(OutputTemplateEntry("NEO.real.dat",      mpz_class("1024")*total_visits*numRealNEO,      true, true));
  outputTemplateVector.push_back(OutputTemplateEntry("NEO.synthetic.dat", mpz_class("1024")*total_visits*numSyntheticNEO, true, true));
  
  // tentative value
  const orsa::Double flops_est = 
    1e12 * 
    (orsa::Double(std::max(512,numRealNEO))/1024.0) * 
    (orsa::Double(std::max(512,numSyntheticNEO))/1024.0) * 
    (orsa::Double(total_tp_lines)/3600.0) *
    (total_FOV_DEG_SQ/4.0);
  
  const unsigned int min_quorum       = 2;
  const unsigned int target_nresults  = 2;
  const orsa::Double rsc_fpops_est    =    flops_est;
  const orsa::Double rsc_fpops_bound  = 32*flops_est;
  const orsa::Time   delay_bound      = std::max(orsa::Time(0,0,0,flops_est/1e8,0),
						 orsa::Time(1,0,0,0,0)); 
  const orsa::Double rsc_memory_bound = mpz_class("134217728");
  //
  orsa::Double rsc_disk_bound = 0;
  for (unsigned int k=0; k<outputTemplateVector.size(); ++k) {
    rsc_disk_bound += outputTemplateVector[k].fileSize.getRef();
  }
  
  const std::string inTemplateName = 
    inputTemplateFileName(baseName,
			  min_quorum,
			  target_nresults,
			  delay_bound,
			  rsc_fpops_est,
			  rsc_fpops_bound,
			  rsc_memory_bound,
			  rsc_disk_bound);    
  
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
  
  const std::string outTemplateName = 
    outputTemplateFileName(baseName,
			   numRealNEO,
			   numSyntheticNEO,
			   total_tp_lines);
  
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
    
  char wuName[1024];
  if (firstSyntheticRandomSeed.isSet()) {
    gmp_snprintf(wuName,
		 1024,
		 "%s_%i_%i_%i.%02i.%02i_%02i.%02i.%02i",
		 baseName.c_str(),
		 firstSyntheticRandomSeed.getRef(),
		 getpid(),
		 y,m,d,H,M,S);
  } else {
    gmp_snprintf(wuName,
		 1024,
		 "%s_%i_%i.%02i.%02i_%02i.%02i.%02i",
		 baseName.c_str(),
		 getpid(),
		 y,m,d,H,M,S);
  }
  
  ORSA_DEBUG("[%s]",wuName);
  
  int retVal;
  
  SCHED_CONFIG config;
  
  config.parse_file();
  
  retVal = boinc_db.open(config.db_name, config.db_host, config.db_user, config.db_passwd);
  //
  if (retVal != 0) {
    ORSA_DEBUG("problems with boinc_db.open(...)");
    exit(0);
  }
  
  DB_APP app;
  //
  retVal = app.lookup("where name='SurveySimulator'");
  //
  if (retVal != 0) {
    ORSA_DEBUG("problems with app.lookup(...)");
    exit(0);
  }
  
  char cmd[1024];
  char path[1024];
  
  char copyParamConf[1024];
  //
  gmp_snprintf(copyParamConf,1024,"param.conf.%s",wuName);
  config.download_path(copyParamConf,path);
  gmp_snprintf(cmd,1024,"cp -f param.conf %s",path);
  system(cmd);
  
  char copyTelescopeConf[1024];
  //
  gmp_snprintf(copyTelescopeConf,1024,"telescope.conf.%s",wuName);
  config.download_path(copyTelescopeConf,path);
  gmp_snprintf(cmd,1024,"cp -f telescope.conf %s",path);
  system(cmd);
  
  char copyObscodeDat[1024];
  //
  gmp_snprintf(copyObscodeDat,1024,"obscode.dat.%s",wuName);
  config.download_path(copyObscodeDat,path);
  gmp_snprintf(cmd,1024,"cp -f obscode.dat %s",path);
  system(cmd);
  
  char copyRealNEOGen[1024];
  //
  gmp_snprintf(copyRealNEOGen,1024,"realNEO.gen.%s",wuName);
  config.download_path(copyRealNEOGen,path);
  gmp_snprintf(cmd,1024,"cp -f realNEO.gen %s",path);
  system(cmd);
  
  char copySyntheticNEOGen[1024];
  //
  gmp_snprintf(copySyntheticNEOGen,1024,"syntheticNEO.gen.%s",wuName);
  config.download_path(copySyntheticNEOGen,path);
  gmp_snprintf(cmd,1024,"cp -f syntheticNEO.gen %s",path);
  system(cmd);
  
  const char * infiles[] = 
    { copyParamConf,
      copyTelescopeConf,
      copyObscodeDat,
      copyRealNEOGen,
      copySyntheticNEOGen 
    };
  
  char * wu_template;
  read_file_malloc(inTemplateName.c_str(), wu_template);
  
  DB_WORKUNIT wu;
  wu.clear(); // zeroes all fields
  wu.appid = app.id;
  strcpy(wu.name, wuName);
  
  create_work(wu,
	      wu_template,
	      outTemplateName.c_str(),
	      outTemplateName.c_str(),
	      infiles,
	      5,
	      config);
}
