#include <orsa/debug.h>
#include <stdlib.h>
#include <libgen.h>
#include "skycoverage.h"

int main(int argc, char ** argv) {
    
    orsa::Debug::instance()->initTimer();
    
    if (argc == 1) {
        printf("Usage: %s <Sky-Coverage-File(s)>\n",argv[0]);
        exit(0);
    }
    
    // before all, write this rule
    // with this, all intermediate files are kept instead of deleted automatically
    // this is needed to keep the *.allEta.dat files around
    printf(".SECONDARY:\n");
    printf("\n");
    
    // first, the targets
    printf("all:");
    for (int arg=1; arg<argc; ++arg) {
        printf(" %s/%s.fit.dat %s/%s.fit.Vxx.dat %s/%s.fit.Uxx.dat %s/%s.fit.pdf %s/%s.fit.jpg",
               dirname(argv[arg]),
               SkyCoverage::basename(argv[arg]).c_str(),
               dirname(argv[arg]),
               SkyCoverage::basename(argv[arg]).c_str(),
               dirname(argv[arg]),
               SkyCoverage::basename(argv[arg]).c_str(),
               dirname(argv[arg]),
               SkyCoverage::basename(argv[arg]).c_str(),
               dirname(argv[arg]),
               SkyCoverage::basename(argv[arg]).c_str());
    }
    printf("\n"); 
    printf("\n");
    
    // now the rules
    printf("%%.allEta.dat: %%.txt\n");
    printf("\t./DetectionEfficiency $*.txt > $*.allEta.log 2>&1\n");
    printf("\n");
    //
    printf("%%.fit.dat: %%.allEta.dat\n");
    printf("\t./DetectionEfficiencyFit $*.allEta.dat > $*.fit.log 2>&1\n");
    printf("\n");
    //
    printf("%%.fit.pdf %%.fit.Vxx.dat %%.fit.Uxx.dat: %%.allEta.dat %%.fit.dat\n");
    printf("\t./DetectionEfficiencyPlot $*.allEta.dat $*.fit.dat > $*.fit.pdf.log 2>&1\n");
    printf("\n");
    //
    printf("%%.fit.jpg: %%.fit.pdf\n");
    printf("\tconvert -density 400 $*.fit.pdf $*.fit.jpg\n");
    printf("\n");
    
    return 0;
}
