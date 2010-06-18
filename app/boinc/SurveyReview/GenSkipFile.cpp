#include <orsa/debug.h>
#include <orsa/double.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

int main(int argc, char ** argv) {

    orsa::Debug::instance()->initTimer();
    
    if (argc == 1) {
        printf("Usage: %s <grid.dat(s)>\n",argv[0]);
        exit(0);
    }
    
    // read grid.dat input file
    double JD;
    unsigned int numGen;
    //
    int z_a_min, z_a_max, z_a_delta;
    int z_e_min, z_e_max, z_e_delta;
    int z_i_min, z_i_max, z_i_delta;
    int z_node_min, z_node_max, z_node_delta;
    int z_peri_min, z_peri_max, z_peri_delta;
    int z_M_min, z_M_max, z_M_delta;
    int z_H_min, z_H_max, z_H_delta;   
    //
    for (int arg=1; arg<argc; ++arg) {
        FILE * fp = fopen(argv[arg],"r"); 
        if (fp) {
            char line[1024];
            while (fgets(line,1024,fp)) {
                if (strlen(line) > 0) {
                    if (line[0] == '#') {
                        continue;
                    }
                }
                if (23 == gmp_sscanf(line,
                                     "%lf %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i",
                                     &JD,
                                     &numGen,
                                     &z_a_min, &z_a_max, &z_a_delta,
                                     &z_e_min, &z_e_max, &z_e_delta,
                                     &z_i_min, &z_i_max, &z_i_delta,
                                     &z_node_min, &z_node_max, &z_node_delta,
                                     &z_peri_min, &z_peri_max, &z_peri_delta,
                                     &z_M_min, &z_M_max, &z_M_delta,
                                     &z_H_min, &z_H_max, &z_H_delta)) {
                    // main output
                    printf("%4i %4i %4i %4i %6i %6i\n",
                           z_a_min, z_a_max,
                           z_e_min, z_e_max,
                           z_i_min, z_i_max);
                    break;
                }
            }
            fclose(fp);
        }
    }
    
    return 0;
}
