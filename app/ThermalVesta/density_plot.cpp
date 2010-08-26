#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <orsaSolarSystem/datetime.h>

#include <iostream>

#include "dislin.h"

/* #define   NX  36
   #define   NY  36
*/

#define   NX  73
#define   NY  37

float mesh[NX][NY];
int   entries[NX][NY];

typedef struct {
    double x,y,z;
} data;

using namespace std;

int main(int argc, char **argv) { 
    
    const double round_T = 20.0; // K, used for contour levels
    
    FILE * input_file;
  
    if (argc < 3) {
        fprintf(stderr,"Usage: %s <JD> <input_file(s)>\n",argv[0]);
        exit(0);
    }
  
    // real requested epoch
    const double JD_input = atof(argv[1]);
  
    // orbit periodicity
    double tmp_JD = JD_input;
    if (tmp_JD < 2454832.50000) tmp_JD += 1325.51455;
    if (tmp_JD > 2456158.01455) tmp_JD -= 1325.51455;
    // lookup epoch in data files
    const double JD = tmp_JD;
  
    const orsa::Time t_JD = orsaSolarSystem::julianToTime(JD_input);
    int y,m,d;
    double fd;
    orsaSolarSystem::gregorDay(t_JD,y,m,d,fd);
  
    for (int nx=0; nx<NX; ++nx) {
        for (int ny=0; ny<NY; ++ny) {
            entries[nx][ny]=0;
            mesh[nx][ny]=0.0;
        }
    }
  
    double T_min = 1.0e10;
    double T_max = 0.0;
    
    for (int infile=2; infile<argc; ++infile) {
    
        if ( (input_file = fopen(argv[infile],"r")) == 0) {
            fprintf(stderr,"ERROR: can't open input file %s.\n",argv[infile]);
            exit(0);
        }
    
        double jd,x,y,z;
    
        // epoch selection
        const double T_rot = .2226; // days
        // const double t0 = 2455788.50000; // 2011 Aug 15
        // const double t0 = 2455880.50000; // 2011 Nov 15
        const double t0 = JD;
        
        // include data from several rotations, to avoid empty bins
        const double T_select = 25*T_rot;
        
        // while (fscanf(input_file,"%lf %lf %lf %*s %*s %*s %*s %lf %*s",&jd,&x,&z,&y) != EOF) { // Temperature
        while (fscanf(input_file,"%lf %lf %*s %*s %*s %*s %lf %lf %*s",&jd,&x,&z,&y) != EOF) { // GRaND coefficient
        
            // exclude some points
            if (jd < t0) continue;
            if (jd > t0+T_select) continue;
      
            // if strictly time ordered, can break after jd > t0+T_select
            if (jd > t0+T_select) break;
      
            const int i_x = (int)round((NX-1)*(x/360.0));
            const int i_y = (int)round((NY-1)*((y+90.0)/180.0));
      
            if (i_x < 0)  printf("ERROR!! i_x=%i\n",i_x);
            if (i_y < 0)  printf("ERROR!! i_y=%i\n",i_y);
            if (i_x >=NX) printf("ERROR!! i_x=%i\n",i_x);
            if (i_y >=NY) printf("ERROR!! i_y=%i\n",i_y);
      
            // printf("%i %i\n",i_x, i_y);
      
            mesh[i_x][i_y] += z;
            entries[i_x][i_y] += 1;
      
            if (z < T_min) T_min=z;
            if (z > T_max) T_max=z;
        }
    
        fclose(input_file);
    }
  
    for (int nx=0; nx<NX; ++nx) {
        for (int ny=0; ny<NY; ++ny) {
            if (entries[nx][ny]>0) mesh[nx][ny] /= entries[nx][ny];
        }
    }
  
    // correct temperature range
    T_max = ceil(T_max/round_T)*round_T;
    T_min = floor(T_min/round_T)*round_T;
  
    // change if needed
    /* T_min =   0.0;
       T_max = 260.0;
    */
    
    // good for printing: metafl("POST") + psfont("AvantGarde-Book")
  
    // metafl("POST");
    // metafl("PSCL");
    // metafl("PNG");
    metafl("GIF");
    // metafl("XWIN"); x11mod("STORE"); clrmod("FULL");
    // metafl("CONS");
    // metafl("TIFF");
  
    pagmod("LAND");
    // page(3600,1400);
    // pagfll(255);
    winsiz(800,600);
  
    // output file name
    char outfilename[1024];
    // keep extension in sync with metafl command
    sprintf(outfilename,"Vesta_%.2f.gif",JD_input);
    setfil(outfilename);
    
    // new files overwrite old ones
    filmod("DELETE");
  
    // background color
    // scrmod("REVERS");
  
    disini();
    // pagera();
    // hwfont();
    simplx();
    // duplx();
    // triplx();
    // helve();
    // psfont("AvantGarde-Book");
    color("fore");
  
    // paghdr("","",2,0);
  
    // select a color table
    setvlt("TEMP"); // TEMP,GREY,RGREY,VGA,RAIN,SPEC...
    // setvlt("RGREY");

    // titlin("3-D Colour Plot of the Function",2);
    // titlin("F(X,Y) = 2 * SIN(X) * SIN(Y)",4);
    // titlin("Saturn Trojans Predictor",4);
  
    // titlin("Vesta Surface Temperature",1);
    titlin("GRaND Modified Temperature at Vesta",1);
    char str_date[1024];
    sprintf(str_date,"%4i/%02i/%02i",y,m,d);
    titlin(str_date,3);
  
    // name("initial semi-major axis [AU]","x");
    // name("libration amplitude","x");
    // name("eccentricity","y");
    // name("Z-axis","z");

    // name("Longitude, Degrees Past Midnight","x");
    name("Longitude, East of Sub-Solar Meridian [deg]","x");
    name("Latitude [deg]","y");
    name("Temperature [K]","z");
   
    // name("Long.","x");
    // name("Lat.","y");
  
    intax();
    autres(NX,NY);
    // axspos(300,1850);
    // ax3len(2200,1400,1400);
    // digits(0,"X");
    digits(0,"X"); 
    // digits(2,"Y");
    digits(0,"Y");
    digits(0,"Z");
  
    ticks(1,"X");
    ticks(1,"Y");
    ticks(1,"Z");
  
    /* graf3(0,360,0,30,
       -90,90,-90,30,
       T_min,T_max,T_min,round_T);
    */
    //
    graf3(-180,180,-180,30,
          -90,90,-90,30,
          T_min,T_max,T_min,round_T);
  
    crvmat((float *)mesh,NX,NY,1,1);
  
    if (1) {
        // contour plot, level curves
        digits(0,"contour");
        labels("FLOAT","CONTUR");
        double T=T_min+round_T;
        while (T<T_max) {
            conmat((float *)mesh,NX,NY,T);
            T += round_T;
        }
    }
  
    // height(50); // chars height
    title();
    // mpaepl(3);
  
    // setvlt("VGA"); // TEMP,GREY,RGREY,VGA,RAIN,SPEC...
    // color("blue");
    // color("red");
  
    disfin();
  
    /* 
       char filename[256], cmd[256];
       sscanf(epoch_filename,"%s",filename);
       sprintf(cmd,"cp -f dislin.tif stp_%s.tif",filename);
       // sprintf(cmd,"cp -f dislin.psc stp_%s.eps",filename);
       system(cmd);
    */

    {
        // write GRaND output file
        if ((NX!=73) || (NY!=37)) {
            ORSA_DEBUG("problems...");
        } else {
            sprintf(outfilename,"Vesta_GRaND_%.2f.dat",JD_input);
            ORSA_DEBUG("writing file [%s]",outfilename);
            FILE * fp = fopen(outfilename,"w");
            for (unsigned int j=0; j<NX; ++j) {
                for (unsigned int k=0; k<NY; ++k) {
                    fprintf(fp,"%+4i %+3i %6.2f\n",j*5-180,k*5-90,mesh[j][k]);
                }
            }
            fclose(fp);
        }
    }
    
    return 0;
}
