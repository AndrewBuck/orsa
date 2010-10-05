#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <iostream>

#include "grain.h"
#include "fit.h"

#include "dislin.h"

#include <orsa/debug.h>
#include <orsa/print.h>
#include <orsa/statistic.h>

using namespace std;

class PlotStatsElement : public orsa::WeightedStatistic<double> {
    
};

class PlotStats : public BinStats<PlotStatsElement> {
public:
    PlotStats(const std::vector< osg::ref_ptr<Var> > & varDefinition) :
        BinStats<PlotStatsElement>(varDefinition) { }
public:
    bool insert(const std::vector<double> & xVector,
                const double & val,
                const double & sigma) {
        if (xVector.size() != var.size()) {
            ORSA_DEBUG("dimension mismatch");
            return false;
        }
        std::vector<size_t> binVector;
        if (!bin(binVector,xVector)) {
            return false;
        }   
        const mpz_class idx = index(binVector);
        if (data[idx].get()==0) {
            // lazy allocation
            data[idx] = new PlotStatsElement;
        }
        data[idx]->insert(val,orsa::square(1.0/sigma));
        return true;         
    }
};

int main(int argc, char **argv) { 
    
    // choose below if plotting "field" coverage, detection efficiency, or observation efficiency
    
    if (argc != 2) {
        ORSA_DEBUG("Usage: %s <inspect_merge_file>",argv[0]);
        exit(0);
    }
    
    FILE * fp;
    if ( (fp = fopen(argv[1],"r")) == 0) {
        ORSA_DEBUG("ERROR: can't open input file %s.",argv[1]);
        exit(0);
    }

    // a,e
    const double x_step = 0.05;
    const double x_min  = 1.90;
    const double x_max  = 2.30;
    //
    const double y_step =  0.05;
    const double y_min  =  0.30;
    const double y_max  =  1.00;
    // a,i
    /* const double x_step = 0.05;
       const double x_min  = 0.90;
       const double x_max  = 2.20;
       //
       const double y_step =  5.00;
       const double y_min  =  0.00;
       const double y_max  = 85.00;
       #warning check why points with inclination between 85 and 90 deg are not plotted
    */
    //
    // a,L (L=node+peri+M
    /* const double x_step = 0.05;
       const double x_min  = 0.90;
       const double x_max  = 2.20;
       //
       const double y_step =  30.0;
       const double y_min  =   0.0;
       const double y_max  = 360.0;
    */
    
    std::vector< osg::ref_ptr<PlotStats::Var> > varDefinition;
    //
    osg::ref_ptr<PlotStats::LinearVar> var_x = new PlotStats::LinearVar(x_min,x_max+2*x_step,x_step);
    varDefinition.push_back(var_x.get());
    //
    osg::ref_ptr<PlotStats::LinearVar> var_y = new PlotStats::LinearVar(y_min,y_max+2*y_step,y_step);
    varDefinition.push_back(var_y.get());
    
    osg::ref_ptr<PlotStats> plotStats = 
        new PlotStats(varDefinition);
    {
        char line[1024];
        std::vector<double> xVector;
        xVector.resize(varDefinition.size());
        double a, e, eta_obs;
        while (fgets(line,1024,fp)) {
            if (3 != sscanf(line,
                            "%lf %lf %*s %lf %*s",
                            &a,
                            &e,
                            &eta_obs)) {
                ORSA_DEBUG("problems with line [%s]",line);
                continue;
            } else {
                
                // a,e
                xVector[0] = a + x_step;
                xVector[1] = e + y_step;
                // a,i
                /* xVector[0] = x_step+center_a;
                   xVector[1] = y_step+center_i;
                */
                //
                // a,L
                /* xVector[0] = x_step+center_a;
                   xVector[1] = y_step+fmod(center_L,360.0);
                */
                
                // ORSA_DEBUG("xVector: %g %g",xVector[0],xVector[1]);
                
                // CHOOSE one insert here
                // plotStats->insert(xVector, eta_field,  1.0);
                    // plotStats->insert(xVector, eta_detect, 1.0);
                plotStats->insert(xVector, eta_obs, 1.0);
                
                // OLD
                // OLD // plotStats->insert(xVector, eta_obs,    sigma_eta_obs);
                // OLD // plotStats->insert(xVector, eta_field,  sigma_eta_field);
                // OLD // plotStats->insert(xVector, eta_detect, sigma_eta_detect);
            }
        }
    }
    
    const double empty_mesh_val=+1000;
    float * mesh;
    const size_t meshSize = plotStats->size().get_si();
    //
    {
        mesh = (float*)calloc(meshSize,sizeof(float));
        std::vector<double> xVector;
        xVector.resize(varDefinition.size());
        for (unsigned j=0; j<var_x->size(); ++j) {
            for (unsigned k=0; k<var_y->size(); ++k) {
                const unsigned int mesh_id = j*var_y->size()+k;
                xVector[0] = x_min+x_step*(j+0.5);
                xVector[1] = y_min+y_step*(k+0.5);
                std::vector<size_t> binVector;
                if (plotStats->bin(binVector,xVector)) {
                    const PlotStatsElement * e =  plotStats->stats(plotStats->index(binVector));
                    if (e) {
                        if (e->average() > 0) {
                            // mesh[mesh_id] = pow10(e->average());
                            mesh[mesh_id] = e->average();
                            // mesh[mesh_id] = log10(e->average());
                        } else {
                            // mesh[mesh_id] = 1e-20;
                            // mesh[mesh_id] = -1000;
                            mesh[mesh_id] = empty_mesh_val;
                         }
                    } else {
                        mesh[mesh_id] = empty_mesh_val;
                    }
                } else {
                    mesh[mesh_id] = empty_mesh_val;
                }
            }
        }
        
        // OLD
        /* const PlotStats::DataType & plotData = plotStats->getData();
           PlotStats::DataType::const_iterator it = plotData.begin();
           while (it != plotData.end()) {
           mesh[(*it).first.get_si()] = (*it).second->average();
           // ORSA_DEBUG("index: %i   val: %g",(*it).first.get_si(),(*it).second->average());
           ++it;
           }
        */
    }

    float mesh_min =  1.0e100;
    float mesh_max = -1.0e100;
    for (unsigned int k=0; k<meshSize; ++k) {
        // ORSA_DEBUG("mesh[%06i] = %g",k,mesh[k]);
        if (mesh[k]==empty_mesh_val) continue;
        if (mesh[k]<mesh_min) mesh_min=mesh[k];
        if (mesh[k]>mesh_max) mesh_max=mesh[k];
    }
    // ORSA_DEBUG("mesh_min: %g   mesh_max: %g",mesh_min,mesh_max);
    //
    // linear
    double mesh_step = 1.0;
    while ((mesh_max/mesh_step)<1.0) mesh_step /= 10.0;
    // correct mesh_max
    mesh_max = ceil(mesh_max/mesh_step)*mesh_step;
    // lower
    mesh_max  *= 0.1;
    mesh_step *= 0.1;
    //
    // logarithmic
    // mesh_max =  ceil(mesh_max);
    // mesh_min = floor(mesh_min);
    // const double mesh_step = 1.0;
    // ORSA_DEBUG("mesh_step: %g",mesh_step);
    // while ((mesh_max/mesh_step)<1.0) mesh_step /= 10.0;
    // correct mesh_max
    // mesh_max = ceil(mesh_max/mesh_step)*mesh_step;
    // lower
    // mesh_max  *= 0.1;
    // mesh_step *= 0.1;
    
    // good for printing: metafl("POST") + psfont("AvantGarde-Book")
    
    /*** DISLIN ***/
    page(2300,2300);
    pagmod("LAND");
    
    // output file name
    /* char plotFilename[1024];
       sprintf(plotFilename,"%s.fit.pdf",basename.c_str());
       setfil(plotFilename);
    */
    // new files overwrite old ones
    filmod("DELETE");
    
    // metafl("POST");
    metafl("PSCL");
    // metafl("PNG");
    // metafl("XWIN"); x11mod("STORE"); clrmod("FULL");
    // metafl("CONS");
    // metafl("TIFF");
    
    pagmod("LAND");
    // page(3600,1400);
    // pagfll(255);
    // winsiz(800,600);
    
    // background color
    scrmod("REVERS");
    
    disini();
    // pagera();
    hwfont();
    simplx();
    // triplx();
    // helve();
    // helves();
    // winfnt();
    // disalf();
    // psfont("AvantGarde-Book");
    // color("fore");
    
    penwid(1.0);
    height(36); // text height
    
    // paghdr("","",2,0);
    
    /* axspos(200,1300);
       axslen(1350,1050);
    */
    axspos( 300,1700);
    axslen( 900,1500);
    
    // select a color table
    setvlt("RAIN"); // TEMP,GREY,RGREY,VGA,RAIN,SPEC...
    // setvlt("GREY");
    
    hwmode("ON","LINE");
    
    texmod("ON"); // TeX text
    
    // NEOs
    // titlin("NEOs In-Field Probability for 703",4);
    // titlin("NEOs In-Field Probability for G96",4);
    // titlin("H=18 NEOs Detection Efficiency for 703",4);
    // titlin("H=18 NEOs Detection Efficiency for G96",4);
    // titlin("H=18 NEOs Observation Probability for 703",4);
    // titlin("H=18 NEOs Observation Probability for G96",4);
    // 
    // PHOs
    // titlin("PHOs In-Field Probability for 703",4);
    // titlin("PHOs In-Field Probability for G96",4);
    // titlin("H=18 PHOs Detection Efficiency for 703",4);
    // titlin("H=18 PHOs Detection Efficiency for G96",4);
    // titlin("H=18 PHOs Observation Probability for 703",4);
    // titlin("H=18 PHOs Observation Probability for G96",4);
    //
    // misc
    // titlin("H=20 PHOs Observation Probability for 703",4);
    // titlin("H=20 PHOs Observation Probability for G96",4);
    // titlin("H=22 PHOs Observation Probability for 703",4);
    // titlin("H=22 PHOs Observation Probability for G96",4);
    
    // titlin("3-D Colour Plot of the Function",2);
    // titlin("F(X,Y) = 2 * SIN(X) * SIN(Y)",4);
    // titlin("Saturn Trojans Predictor",4);

    titlin("NEOs H=18",4);
    
    // name("initial semi-major axis [AU]","x");
    // name("libration amplitude","x");
    // name("eccentricity","y");
    // name("Z-axis","z");

    // a,e
    name("Semi-Major Axis [AU]","x");
    name("Eccentricity","y");
    // a,i
    /* name("Semi-Major Axis [AU]","x");
       name("Inclination [deg]","y");
    */
    // a,L
    /* name("Semi-Major Axis [AU]","x");
       name("True Longitude [deg]","y");
    */
    //
    // name("Probability","z");
    // name("Detection Efficiency","z");
    // name("Probability","z");
    
    // name("Long.","x");
    // name("Lat.","y");
    
    // axsscl("log","z");
    // labels("float","y");  

    setgrf("NAME","NAME","NONE","NONE");
    
    intax();

    // a,e
    /* autres(var_x->size()+3,
       var_y->size()+3);
    */
    /* autres(var_x->size()+1,
       var_y->size()+2);
    */
    autres(var_x->size(),
           var_y->size());
    // a,i
    /* autres(var_x->size()+3,
       var_y->size()+3);
    */
    // a,L
    /* autres(var_x->size()+3,
       var_y->size()+1);
    */
    
    labtyp("vert","z"); // vertical labels for z axis

    // a,e
    //
    // NEOs
    // const double z_min=1e-4; const double z_max=1e-1;
    // const double z_min=1e-3; const double z_max=1e-1;
    // const double z_min=1e-5; const double z_max=1e-3;
    //
    // PHOs
    // const double z_min=1e-4; const double z_max=1e-1;
    // const double z_min=1e-3; const double z_max=1e0;
    // const double z_min=1e-6; const double z_max=1e-3;
    //
    // other... LINEAR 0->1
    const double z_min=0.0; const double z_max=0.05;
    //
    // a,i
    // const double z_min=1e-5; const double z_max=1e-3;
    //
    // a,L
    // const double z_min=1e-5; const double z_max=1e-3;
    //
    // 
    {
        // bound z
        for (unsigned int k=0; k<meshSize; ++k) {
            if (mesh[k]==empty_mesh_val) continue;
            if (mesh[k]<z_min) mesh[k]=z_min;
            if (mesh[k]>z_max) mesh[k]=z_max;
        }
    }
    // a,e
    /* {
       digits(1,"X"); 
       digits(1,"Y");
       digits(0,"Z");
       ticks(1,"X");
       ticks(1,"Y");
       ticks(5,"Z");
       axsscl("log","z");
       labels("log","z");
       frame(5); // frame thickness
       graf3(x_min-x_step/2,x_max+x_step/2,1.0,0.2,
       y_min-y_step/2,y_max+y_step/2,y_min,0.1,
       log10(z_min),log10(z_max),log10(z_min),1);
       crvmat(mesh,var_x->size(),var_y->size(),1,1);
       }
    */
    // a,e
    {
        digits(1,"X"); 
        digits(1,"Y");
        digits(2,"Z");
        ticks(2,"X");
        ticks(2,"Y");
        ticks(2,"Z");
        // axsscl("log","z");
        // labels("log","z");
        frame(5); // frame thickness
        /* graf3(x_min-x_step/2,x_max+x_step/2,1.0,0.2,
           y_min-y_step/2,y_max+y_step/2,y_min,0.1,
           log10(z_min),log10(z_max),log10(z_min),1);
        */
        graf3(x_min-0.5*x_step,x_max+0.5*x_step,x_min,0.1,
              y_min-0.5*y_step,y_max+0.5*y_step,y_min,0.1,
              z_min,z_max,z_min,0.01);
        crvmat(mesh,var_x->size(),var_y->size(),1,1);
    }
    // a,i
    /* {
       digits(1,"X"); 
       digits(0,"Y");
       digits(0,"Z");
       ticks(1,"X");
       ticks(1,"Y");
       ticks(5,"Z");
       axsscl("log","z");
       labels("log","z");
       frame(5); // frame thickness
       graf3(x_min-x_step/2,x_max+x_step/2,1.0,0.2,
       y_min-y_step/2,y_max+y_step/2,0.0,15.0,
       log10(z_min),log10(z_max),log10(z_min),1);
       crvmat(mesh,var_x->size(),var_y->size(),1,1);
       }
    */
    // a,L
    /* {
       digits(1,"X"); 
       digits(0,"Y");
       digits(0,"Z");
       ticks(1,"X");
       ticks(3,"Y");
       ticks(5,"Z");
       axsscl("log","z");
       labels("log","z");
       frame(5); // frame thickness
       graf3(x_min-x_step/2,x_max+x_step/2,1.0,0.2,
       y_min-y_step/2,y_max+y_step/2,y_min,90.0,
       log10(z_min),log10(z_max),log10(z_min),1);
       crvmat(mesh,var_x->size(),var_y->size(),1,1);
       }
    */
    
    
    /* 
       if (1) {
       // contour plot, level curves
       digits(1,"contour");
       labels("FLOAT","CONTUR");
       double T=eta_min+round_eta;
       while (T<eta_max) {
       conmat((float *)mesh,NX,NY,T);
       T += round_eta;
       }
       }
    */
    
    // title only
    vkytit(-50); // title closer to plot
    height(50); // text height
    title();
    
    disfin();
    
    fclose(fp);
    
    return 0;
}
