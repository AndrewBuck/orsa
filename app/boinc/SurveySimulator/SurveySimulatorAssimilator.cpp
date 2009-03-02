#include <boinc_db.h>
#include <validate_util.h>

#include <stdlib.h>

#include "error_numbers.h"
#include "parse.h"
#include "util.h"
#include "filesys.h"

#include "sched_util.h"
#include "sched_config.h"
#include "sched_msgs.h"
#include "validator.h"
#include "validate_util.h"

#include <iostream>

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

static void print_cmd(const char   * cmd,
		      const RESULT & canonical_result) {
  log_messages.printf(MSG_DEBUG, 
		      "[RESULT#%d %s] executing command: %s\n", 
		      canonical_result.id,  
		      canonical_result.name,
		      cmd); 
}

int assimilate_handler(WORKUNIT            & wu, 
		       std::vector<RESULT> & , 
		       RESULT              & canonical_result) {
  
  if (wu.canonical_resultid==0) {
    log_messages.printf(MSG_CRITICAL,  
			"[WU#%d %s] no canonical result, error_mask: %i\n",  
			wu.id,   
			wu.name, 
			wu.error_mask);
    return 0;
  }
  
  if (wu.error_mask!=0) {
    log_messages.printf(MSG_CRITICAL,   
                        "[WU#%d %s] error_mask: %i\n",   
                        wu.id,    
                        wu.name,  
                        wu.error_mask); 
    return 0;
  }
  
  const char * output_dir_base = "/home/boinc/SurveySimulatorStorage";
  
  char cmd[1024];  
 
  char output_dir[1024];
  sprintf(output_dir,"%s/%s",
	  output_dir_base,
	  wu.name);
  sprintf(cmd,"mkdir -p %s",
	  output_dir);
  print_cmd(cmd,canonical_result);
  system(cmd);
  
  /* 
     log_messages.printf(MSG_DEBUG, 
     "[RESULT#%d %s] wu.xml_doc:\n%s", 
     canonical_result.id,  
     canonical_result.name,  
     wu.xml_doc); 
  */
  
  {
    std::vector<FILE_INFO> fis;
    get_input_file_infos(wu,fis);
    std::vector<FILE_INFO>::const_iterator it = fis.begin();
    while (it != fis.end()) {
      // test if file exists 
      FILE * fp = fopen((*it).path.c_str(),"r"); 
      if (fp) { 
        fclose(fp); 
        sprintf(cmd,"cp %s %s", 
                (*it).path.c_str(), 
                output_dir); 
        print_cmd(cmd,canonical_result); 
        system(cmd); 
      } else { 
	log_messages.printf(MSG_CRITICAL, 
			    "[RESULT#%d %s] Couldn't open %s\n", 
			    canonical_result.id,  
			    canonical_result.name,  
			    (*it).path.c_str()); 
	return 1; 
      }
      ++it;
    }
  }
  
  /* 
     log_messages.printf(MSG_DEBUG, 
     "[RESULT#%d %s] canonical_result.xml_doc_in:\n%s", 
     canonical_result.id,  
     canonical_result.name,  
     canonical_result.xml_doc_in); 
     log_messages.printf(MSG_DEBUG,  
     "[RESULT#%d %s] canonical_result.xml_doc_out:\n%s",  
     canonical_result.id,   
     canonical_result.name,   
     canonical_result.xml_doc_out);  
  */
  
  {
    std::vector<FILE_INFO> fis;
    get_output_file_infos(canonical_result,fis);
    std::vector<FILE_INFO>::const_iterator it = fis.begin();
    while (it != fis.end()) {
      /* 
	 log_messages.printf(MSG_DEBUG,  
	 "[RESULT#%d %s] name: %s  path: %s  optional: %i\n",  
	 canonical_result.id,   
	 canonical_result.name,   
	 (*it).name.c_str(),
	 (*it).path.c_str(),
	 (*it).optional);  
      */
      // test if file exists
      FILE * fp = fopen((*it).path.c_str(),"r");
      if (fp) {
	fclose(fp);
	sprintf(cmd,"cp %s %s",
		(*it).path.c_str(),
		output_dir);
	print_cmd(cmd,canonical_result);
	system(cmd);
      } else {
	if (!(*it).optional) {
	  log_messages.printf(MSG_CRITICAL,
			      "[RESULT#%d %s] Couldn't open %s\n",
			      canonical_result.id, 
			      canonical_result.name, 
			      (*it).path.c_str());
	  return 1;
	}
      }
      ++it;
    }
  }
  
  return 0;
}

