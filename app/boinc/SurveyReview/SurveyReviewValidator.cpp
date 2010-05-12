#include <string>
#include <vector>
#include <cmath>

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

struct DATA {
    int i;
    double x;
};

int init_result(RESULT & result , void * & data) {
    
    FILE* f;
    FILE_INFO fi;
    int i, n, retval;
    double x;

    retval = get_output_file_path(result, fi.path);
    if (retval) return retval;
    retval = try_fopen(fi.path.c_str(), f, "r");
    if (retval) return retval;
    n = fscanf(f, "%d %f", &i, &x);
    fclose(f);
    if (n != 2) return ERR_XML_PARSE;
    DATA* dp = new DATA;
    dp->i = i;
    dp->x = x;
    data = (void*) dp;
    return 0;
}

int compare_results(
    RESULT& r1, void* _data1, RESULT const& r2, void* _data2, bool& match
) {
    DATA* data1 = (DATA*)_data1;
    DATA* data2 = (DATA*)_data2;
    match = true;
    if (data1->i != data2->i) match = false;
    if (fabs(data1->x - data2->x) > 0.01) match = false;
    return 0;
}

int cleanup_result(RESULT const& r, void* data) {
    if (data) delete (DATA*) data;
    return 0;
}

double compute_granted_credit(WORKUNIT& wu, vector<RESULT>& results) {
    return median_mean_credit(wu, results);
}

