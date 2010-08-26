#include <orsa/debug.h>
#include <cstdlib>
#include <libgen.h>
#include <cstdio>

int main(int argc, char ** argv) {
  
    orsa::Debug::instance()->initTimer();
  
    if (argc < 2) {
        printf("Usage: %s <GRaND-scale(g/cm^2)>\n",argv[0]);
        exit(0);
    }
    
    const char * GRaND_scale = argv[1];

    const double lon = +10.0;
    
    // include this rule to keep intermediate files
    {
        printf(".SECONDARY:\n");
        printf("\n");
    }
    
    // first, the targets
    {
        printf("all:");
        double lat=-90;
        while (lat <= +90) {
            printf(" temperature_orbital_%+05.2f_%g_%sgcm2_final.dat",
                   lat,
                   lon,
                   GRaND_scale);
            lat += 5;
        }
        printf("\n"); 
        printf("\n");
    }
    
    // now the rules
    {      
        double lat=-90;
        while (lat <= +90) {
            printf("temperature_orbital_%+05.2f_%g_%sgcm2_final.dat:\n"
                   "\t./ThermalVesta %g %g %s > temperature_orbital_%+05.2f_%g_%sgcm2.log 2>&1\n\n",
                   lat,
                   lon,
                   GRaND_scale,
                   lat,
                   lon,
                   GRaND_scale,
                   lat,
                   lon,
                   GRaND_scale);
            lat += 5;
        }
        printf("\n"); 
        printf("\n");
    }
    
    return 0;
}
