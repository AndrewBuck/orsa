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
    
    // choose H
    const int z_H_fix = 180;
    
    // choose NEO or PHO
    const std::string OBJ = "NEO"; 
    // const std::string OBJ = "PHO"; 
    
    // also choose below if plotting "field" coverage, detection efficiency, or observation efficiency
    
    if (argc != 2) {
        ORSA_DEBUG("Usage: %s <inspect_file>",argv[0]);
        exit(0);
    }

    FILE * fp;
    if ( (fp = fopen(argv[1],"r")) == 0) {
        ORSA_DEBUG("ERROR: can't open input file %s.",argv[1]);
        exit(0);
    }

    const double a_step = 0.05;
    const double a_min  = 0.90;
    const double a_max  = 2.20;
    //
    const double e_step =  0.05;
    const double e_min  =  0.00;
    const double e_max  =  1.00;
    //
    /* const double i_step =   5.0;
       const double i_min  =   0.0;
       const double i_max  = 150.0;
    */
    
    std::vector< osg::ref_ptr<PlotStats::Var> > varDefinition;
    //
    // [0] a
    osg::ref_ptr<PlotStats::LinearVar> var_a = new PlotStats::LinearVar(a_min+1*a_step,a_max+3*a_step,a_step);
    varDefinition.push_back(var_a.get());
    //
    // [1] e
    osg::ref_ptr<PlotStats::LinearVar> var_e = new PlotStats::LinearVar(e_min+1*e_step,e_max+3*e_step,e_step);
    varDefinition.push_back(var_e.get());
    //
    // [2]
    /* osg::ref_ptr<PlotStats::LinearVar> var_i = new PlotStats::LinearVar(i_min,i_max,i_step);
       varDefinition.push_back(var_i.get());
    */
    //
    osg::ref_ptr<PlotStats> plotStats = 
        new PlotStats(varDefinition);
    {
        char line[1024];
        std::vector<double> xVector;
        xVector.resize(varDefinition.size());
        char obj[1024];
        int z_a_min, z_a_max;
        int z_e_min, z_e_max;
        int z_i_min, z_i_max;
        int z_H;
        double eta_field, sigma_eta_field;
        int entries_field;
        double eta_detect, sigma_eta_detect;
        int entries_detect;
        double eta_obs, sigma_eta_obs;
        while (fgets(line,1024,fp)) {
            if (16 != sscanf(line,
                             "%s %i %i %i %i %i %i %i %lf %lf %i %lf %lf %i %lf %lf",
                             obj,
                             &z_a_min, &z_a_max,
                             &z_e_min, &z_e_max,
                             &z_i_min, &z_i_max,
                             &z_H,
                             &eta_field, &sigma_eta_field,
                             &entries_field,
                             &eta_detect, &sigma_eta_detect,
                             &entries_detect,
                             &eta_obs, &sigma_eta_obs)) {
                ORSA_DEBUG("problems with line [%s]",line);
                continue;
            } else {

                if ( (z_H == z_H_fix) &&
                     (obj == OBJ) &&
                     (entries_detect >= 1) ) {
                    // keep vars aligned with varDefinition content
                    xVector[0] = a_step+0.5*(z_a_max+z_a_min)*grain_a_AU;
                    xVector[1] = e_step+0.5*(z_e_max+z_e_min)*grain_e;
                    // xVector[2] = 0.5*(z_i_max+z_i_min)*grain_i_DEG;
                    // ORSA_DEBUG("xVector: %g %g",xVector[0],xVector[1]);
                    //
                    /* xVector[0] = z_a_min*grain_a_AU;
                       xVector[1] = z_e_min*grain_e;
                       // xVector[2] = z_i_min*grain_i_DEG;
                       */
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
    }
    
    const double empty_mesh_val=+1000;
    float * mesh;
    const size_t meshSize = plotStats->size().get_si();
    //
    {
        mesh = (float*)calloc(meshSize,sizeof(float));
        std::vector<double> xVector;
        xVector.resize(varDefinition.size());
        for (unsigned j=0; j<var_a->size(); ++j) {
            for (unsigned k=0; k<var_e->size(); ++k) {
                const unsigned int mesh_id = j*var_e->size()+k;
                xVector[0] = a_min+a_step*(j+0.5);
                xVector[1] = e_min+e_step*(k+0.5);
                std::vector<size_t> binVector;
                if (plotStats->bin(binVector,xVector)) {
                    // ORSA_DEBUG("mesh[%i]  totalSize: %i  j: %i k: %i bv[0]: %i bv[1]: %i",j*var_e->size()+k,plotStats->size().get_si(),j,k,binVector[0],binVector[1]);
                    const PlotStatsElement * e =  plotStats->stats(plotStats->index(binVector));
                    if (e) {
                        // mesh[mesh_id] = e->average();

                        if (e->average() > 0) {
                            mesh[mesh_id] = log10(e->average());
                        } else {
                            mesh[mesh_id] = -1000;
                        }
                        
                        // ORSA_DEBUG("a: %g   e: %g   j: %i   k: %i   mesh[%06i] = %g",xVector[0],xVector[1],j,k,mesh_id,mesh[mesh_id]);
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
    /* double mesh_step = 1.0;
       while ((mesh_max/mesh_step)<1.0) mesh_step /= 10.0;
       // correct mesh_max
       mesh_max = ceil(mesh_max/mesh_step)*mesh_step;
       // lower
       mesh_max  *= 0.1;
       mesh_step *= 0.1;
    */
    //
    // logarithmic
    mesh_max =  ceil(mesh_max);
    mesh_min = floor(mesh_min);
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
    page(2000,1500);
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
    // hwfont();
    simplx();
    // triplx();
    // helve();
    // psfont("AvantGarde-Book");
    // color("fore");
    
    penwid(0.2);
    height(36); // text height
    
    // paghdr("","",2,0);
    
    axspos(200,1300);
    axslen(1350,1050);
    
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
    titlin("H=18 NEOs Observation Probability for 703",4);
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
    // titlin("H=22 PHOs Observation Probability for 703",4);
    // titlin("H=22 PHOs Observation Probability for G96",4);
    
    // titlin("3-D Colour Plot of the Function",2);
    // titlin("F(X,Y) = 2 * SIN(X) * SIN(Y)",4);
    // titlin("Saturn Trojans Predictor",4);
    
    // name("initial semi-major axis [AU]","x");
    // name("libration amplitude","x");
    // name("eccentricity","y");
    // name("Z-axis","z");

    name("Semi-Major Axis [AU]","x");
    name("Eccentricity","y");
    //
    // name("Log$_{10}$ Probability","z");
    // name("Log$_{10}$ Detection Efficiency","z");
    name("Log$_{10}$ Probability","z");
    
    // name("Long.","x");
    // name("Lat.","y");
    
    // axsscl("log","z");
    // labels("float","y");  

    setgrf("NAME","NAME","NONE","NONE");
    
    intax();
    autres(var_a->size()+3,
           var_e->size()+3);
    // axspos(300,1850);
    // ax3len(2200,1400,1400);
    // digits(0,"X");
    digits(1,"X"); 
    // digits(2,"Y");
    digits(1,"Y");
    digits(0,"Z");
    
    ticks(1,"X");
    ticks(1,"Y");
    ticks(5,"Z");
    
    // lin
    /* graf3(a_min-a_step/2,a_max+a_step/2,a_min,0.1,
       e_min-e_step/2,e_max+e_step/2,e_min,0.1,
       0,mesh_max,0,mesh_step);
       crvmat(mesh,var_a->size(),var_e->size(),1,1);
    */
    // log
    /* graf3(a_min-a_step/2,a_max+a_step/2,a_min,0.1,
       e_min-e_step/2,e_max+e_step/2,e_min,0.1,
       mesh_min,mesh_max,mesh_min,mesh_step);
    */
    /* graf3(a_min,a_max,a_min,0.1,
       e_min,e_max,e_min,0.1,
       mesh_min,mesh_max,mesh_min,mesh_step);
       crvmat(mesh,var_a->size(),var_e->size(),1,1);
    */
    /* graf3(a_min-a_step/2,a_max+a_step/2,a_min,0.1,
       e_min-e_step/2,e_max+e_step/2,e_min,0.1,
       mesh_min,mesh_max,mesh_min,mesh_step);
       crvmat(mesh,var_a->size(),var_e->size(),1,1);
    */
    //
    // NEOs
    // const double z_min=-4; const double z_max=-1;
    // const double z_min=-3; const double z_max=-1;
    const double z_min=-5; const double z_max=-3;
    //
    // PHOs
    // const double z_min=-4; const double z_max=-1;
    // const double z_min=-3; const double z_max=0;
    // const double z_min=-6; const double z_max=-3;
    //
    // other...
    // const double z_min=-6; const double z_max=-3;
    // 
    {
        // bound z
        for (unsigned int k=0; k<meshSize; ++k) {
            if (mesh[k]==empty_mesh_val) continue;
            if (mesh[k]<z_min) mesh[k]=z_min;
            if (mesh[k]>z_max) mesh[k]=z_max;
        }
    }
    //
    graf3(a_min-a_step/2,a_max+a_step/2,a_min,0.1,
          e_min-e_step/2,e_max+e_step/2,e_min,0.1,
          z_min,z_max,z_min,1);
    crvmat(mesh,var_a->size(),var_e->size(),1,1);
    
    // ORSA_DEBUG("var_a->size(): %i",var_a->size());
    
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
    height(50); // text height
    title();
    
    disfin();
    
    fclose(fp);
    
    return 0;
}
