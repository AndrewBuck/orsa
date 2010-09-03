#include "skycoverage.h"
#include "eta.h"
#include "fit.h"
#include <vector>
#include "dislin.h"

double bigfun(double (*fun)(const double, const FitFileDataElement &),
              const FitFileData & ffd,
              const double & x) {
    double y;
    if (ffd.size() != 0) {
        y = 0;
        for (unsigned int ffde=0; ffde<ffd.size(); ++ffde) {
            y += fun(x,ffd[ffde]);
        }
        y /= ffd.size();
    } else {
        ORSA_DEBUG("empty ffd?");
        y = 0;
    }
    return y;
}

class XYE {
public:
    double x,y,e;
};

void PlotUtil(double (*fun)(const double, const FitFileDataElement &),
              const FitFileData & ffd,
              const Histo & histo,
              const double xfactor,
              const int xpos,
              const int ypos,
              const int xlen,
              const int ylen,
              const double x1,
              const double x2,
              const double dx,
              const int xdigits,
              const int xticks,
              const std::string & xlabel,
              const double y1,
              const double y2,
              const double dy,
              const int ydigits,
              const int yticks,
              const std::string & ylabel) {
    penwid(0.2);
    height(30); // text height
    axslen(xlen,ylen);
    axspos(xpos,ypos);
    digits(xdigits,"x");
    digits(ydigits,"y");
    ticks(xticks,"x");
    ticks(yticks,"y");
    namdis(25,"x"); // distance between axis name and labels
    namdis(25,"y");
    name(xlabel.c_str(),"x");
    name(ylabel.c_str(),"y");
    frame(0);
    setgrf("NAME","NAME","NONE","NONE");
    graf(x1,x2,x1,dx,y1,y2,y1,dy);
    // xcross();
    // xaxgit();
    // setgrf("NAME","NAME","NONE","NONE");
    {
        const int nPoints=1000;
        float xray[nPoints], yray[nPoints];
        for (int j=0; j<nPoints; ++j) {
            xray[j] = x1+(j*(x2-x1))/nPoints;
            yray[j] = std::max(0.0,bigfun(fun,ffd,xray[j]/xfactor));
        }
        penwid(0.3);
        dashl();
        curve(xray,yray,nPoints);
    }        
    {
        std::vector<XYE> xye_vec;
        // ORSA_DEBUG("histo.getData().size: %i",histo.getData().size());
        for (unsigned int j=0; j<histo.getData().size(); ++j) {
            XYE xye;
            xye.x=xfactor*histo.getData()[j]->center;
            xye.y=histo.getData()[j]->average();
            xye.e=histo.getData()[j]->standardDeviation();
            // ORSA_DEBUG("x: %g y: %g e: %g",xye.x,xye.y,xye.e);
            if ( (xye.x>=x1) &&
                 (xye.x<=x2) &&
                 (xye.y>=y1) &&
                 (xye.y<=y2) ) {
                xye_vec.push_back(xye);
            }
        }
        const int nPoints=xye_vec.size();
        // ORSA_DEBUG("nPoints: %i",nPoints);
        float xray[nPoints], yray[nPoints];  
        float eray[nPoints];
        for (int j=0; j<nPoints; ++j) {
            xray[j]=xye_vec[j].x;
            yray[j]=xye_vec[j].y;
            eray[j]=xye_vec[j].e;
            // ORSA_DEBUG("x: %g y: %g e: %g",xray[j],yray[j],eray[j]);
        }
        marker(1); // marker type
        hsymbl(5); // marker size
        // smaller pen for smaller lines
        penwid(0.2); 
        solid();
        errbar(xray,yray,eray,eray,nPoints);
    }
    endgrf();
}

double PlotUtil_fun_V(const double V,
                      const FitFileDataElement & e) {
    return SkyCoverage::nominal_eta_V(V,
                                      e.V_limit,
                                      e.eta0_V,
                                      e.V0,
                                      e.c_V,
                                      e.w_V);
}

double PlotUtil_fun_U(const double U,
                      const FitFileDataElement & e) {
    return SkyCoverage::nominal_eta_U(U,
                                      e.U_limit,
                                      e.w_U);
}

double PlotUtil_fun_AM(const double AM,
                       const FitFileDataElement & e) {
    return SkyCoverage::nominal_eta_AM(AM,
                                       e.peak_AM,
                                       e.scale_AM,
                                       e.shape_AM);
}

double PlotUtil_fun_GB(const double GB,
                       const FitFileDataElement & e) {
    return SkyCoverage::nominal_eta_GB_GL(GB,
                                          e.drop_GB,
                                          e.scale_GB,
                                          0.0,
                                          e.scale_GL,
                                          0.0);
}

double PlotUtil_fun_GL(const double GL,
                       const FitFileDataElement & e) {
    return SkyCoverage::nominal_eta_GB_GL(0.0,
                                          e.drop_GB,
                                          e.scale_GB,
                                          GL,
                                          e.scale_GL,
                                          e.shape_GL);
}

/* 
   double PlotUtil_fun_SA(const double SA,
   const FitFileDataElement & e) {
   return SkyCoverage::nominal_eta_SA(SA,
   e.peak_SA,
   e.scale_SA,
   e.shape_SA);
   }
*/

/* double PlotUtil_fun_LA(const double LA,
   const FitFileDataElement & e) {
   return SkyCoverage::nominal_eta_AM(LA,
   e.peak_LA,
   e.scale_LA,
   e.shape_LA);
   }
*/

/* double PlotUtil_fun_LE(const double LE,
   const FitFileDataElement & e) {
   return SkyCoverage::nominal_eta_LE(LE,
   e.peak_LE,
   e.scale_LE,
   e.shape_LE);
   }
*/

double PlotUtil_fun_one(const double,
                        const FitFileDataElement &) {
    return 1.0;
}

int main(int argc, char ** argv) {
    
    orsa::Debug::instance()->initTimer();
    
    if ((argc==1) || (argc%2==0)) {
        printf("Usage: %s <allEta-file-1> <fit-file-1> [eta-file-2 fit-file-2 ... eta-file-N fit-file-N]\n",argv[0]);
        exit(0);
    }
    
    ORSA_DEBUG("pid: %i",getpid());
    
    // this is HALF the number of files, because each data set is composed by two files
    const unsigned int numFiles=(argc-1)/2;
    
    std::vector< osg::ref_ptr<CountStats::Var> > varDefinition;
    //
    // [0] apparent magnitude
    osg::ref_ptr<CountStats::LinearVar> var_V = new CountStats::LinearVar(start_V,stop_V,step_V);
    varDefinition.push_back(var_V.get());
  
    // [1] apparent velocity
    osg::ref_ptr<CountStats::LinearVar> var_U = new CountStats::LinearVar(start_U,stop_U,step_U);
    varDefinition.push_back(var_U.get());
  
    // [ ] solar elongation
    /* osg::ref_ptr<CountStats::LinearVar> var_SE = new CountStats::LinearVar(start_SE,stop_SE,step_SE);
       varDefinition.push_back(var_SE.get());
    */
    
    // [ ] lunar elongation
    /* osg::ref_ptr<CountStats::LinearVar> var_LE = new CountStats::LinearVar(start_LE,stop_LE,step_LE);
       varDefinition.push_back(var_LE.get());
    */
    
    // [2] airmass
    osg::ref_ptr<CountStats::LinearVar> var_AM = new CountStats::LinearVar(start_AM,stop_AM,step_AM);
    varDefinition.push_back(var_AM.get());
  
    // [3] galactic longitude
    osg::ref_ptr<CountStats::LinearVar> var_GL = new CountStats::LinearVar(start_GL,stop_GL,step_GL);
    varDefinition.push_back(var_GL.get());
  
    // [4] galactic latitude
    osg::ref_ptr<CountStats::LinearVar> var_GB = new CountStats::LinearVar(start_GB,stop_GB,step_GB);
    varDefinition.push_back(var_GB.get());
    
    // [ ] ecliptic longitude
    /* osg::ref_ptr<CountStats::LinearVar> var_EL = new CountStats::LinearVar(start_EL,stop_EL,step_EL);
       varDefinition.push_back(var_EL.get());
    */
    
    // [ ] ecliptic latitude
    /* osg::ref_ptr<CountStats::LinearVar> var_EB = new CountStats::LinearVar(start_EB,stop_EB,step_EB);
       varDefinition.push_back(var_EB.get());
    */
    
    // [5] azimuth
    osg::ref_ptr<CountStats::LinearVar> var_AZ = new CountStats::LinearVar(start_AZ,stop_AZ,step_AZ);
    varDefinition.push_back(var_AZ.get());
  
    // [ ] lunar altitude
    /* osg::ref_ptr<CountStats::LinearVar> var_LA = new CountStats::LinearVar(start_LA,stop_LA,step_LA);
       varDefinition.push_back(var_LA.get());
    */
    
    // [6] solar altitude
    osg::ref_ptr<CountStats::LinearVar> var_SA = new CountStats::LinearVar(start_SA,stop_SA,step_SA);
    varDefinition.push_back(var_SA.get());
    
    // [ ] lunar phase
    /* osg::ref_ptr<CountStats::LinearVar> var_LP = new CountStats::LinearVar(start_LP,stop_LP,step_LP);
       varDefinition.push_back(var_LP.get());
    */
    
    // [7] lunar altitude
    osg::ref_ptr<CountStats::LinearVar> var_LA = new CountStats::LinearVar(start_LA,stop_LA,step_LA);
    varDefinition.push_back(var_LA.get());
    
    // [8] lunar illumination
    osg::ref_ptr<CountStats::LinearVar> var_LI = new CountStats::LinearVar(start_LI,stop_LI,step_LI);
    varDefinition.push_back(var_LI.get());

    // for data plots
    Histo histo_V (var_V.get());
    Histo histo_U (var_U.get());
    // Histo histo_SE(var_SE.get());
    // Histo histo_LE(var_LE.get());
    Histo histo_AM(var_AM.get());
    Histo histo_GB(var_GB.get());
    Histo histo_GL(var_GL.get());
    // Histo histo_EB(var_EB.get());
    // Histo histo_EL(var_EL.get());
    Histo histo_AZ(var_AZ.get());
    Histo histo_SA(var_SA.get());
    // Histo histo_LP(var_LP.get());
    Histo histo_LA(var_LA.get());
    Histo histo_LI(var_LI.get());
    
    // specialized, 2D Histo
    Histo2D histo2D_galactic(var_GB.get(),var_GL.get());
    Histo2D histo2D_galactic_eta(var_GB.get(),var_GL.get());
    // ecliptic
    /* Histo2D histo2D_ecliptic(var_EB.get(),var_EL.get());
       Histo2D histo2D_ecliptic_eta(var_EB.get(),var_EL.get());
    */
    //
    // lunar (choose ONE and update code everywhere)
    // Histo2D histo2D_lunar(var_LA.get(),var_LE.get());
    // Histo2D histo2D_lunar_eta(var_LA.get(),var_LE.get());
    //
    // Histo2D histo2D_lunar(var_LE.get(),var_LP.get());
    // Histo2D histo2D_lunar_eta(var_LE.get(),var_LP.get());
    //
    // Histo2D histo2D_lunar(var_LA.get(),var_LP.get());
    // Histo2D histo2D_lunar_eta(var_LA.get(),var_LP.get());
    // 
    Histo2D histo2D_lunar(var_LA.get(),var_LI.get());
    Histo2D histo2D_lunar_eta(var_LA.get(),var_LI.get());

    // EXTRA
    //
    // [ ] apparent velocity
    osg::ref_ptr<CountStats::LinearVar> var_U_more = new CountStats::LinearVar(start_U,stop_U,20*step_U);
    varDefinition.push_back(var_U_more.get());
    //
    Histo histo_U_more (var_U_more.get());
    
    // read fit files
    FitFileData ffd;
    ffd.resize(numFiles);
    for (unsigned int fileID=0; fileID<numFiles; ++fileID) {
        FitFileDataElement e;
        if (!readFitFile(e,argv[2*fileID+2])) {
            ORSA_DEBUG("problems...");
            exit(0);
        }
        ffd[fileID] = e;
    }
    
    std::vector<std::string> basename;
    basename.resize(numFiles);
    
    for (unsigned int fileID=0; fileID<numFiles; ++fileID) {
        
        basename[fileID] = SkyCoverage::basename(argv[2*fileID+1]);
        
        std::vector<EfficiencyData> etaData;
        
        readEfficiencyDataFile(etaData,argv[2*fileID+1]);
        
        ORSA_DEBUG("etaData.size(): %i",etaData.size());
        
        if (etaData.size()==0) continue;
        
        osg::ref_ptr<CountStats> countStats = 
            new CountStats(varDefinition);
        
        {
            std::vector<double> xVector;
            xVector.resize(varDefinition.size());
            unsigned int excluded=0;
            for (unsigned int k=0; k<etaData.size(); ++k) {
                // keep vars aligned with varDefinition content
                xVector[0]  = etaData[k].V.getRef();
                xVector[1]  = etaData[k].apparentVelocity.getRef();
                // xVector[ ]  = etaData[k].solarElongation.getRef();
                // xVector[ ]  = etaData[k].lunarElongation.getRef();
                xVector[2]  = etaData[k].airMass.getRef();
                xVector[3]  = etaData[k].galacticLongitude.getRef();
                xVector[4]  = etaData[k].galacticLatitude.getRef();
                // xVector[ ]  = etaData[k].eclipticLongitude.getRef();
                // xVector[ ]  = etaData[k].eclipticLatitude.getRef();
                xVector[5]  = etaData[k].azimuth.getRef();
                xVector[6]  = etaData[k].solarAltitude.getRef();
                // xVector[ ] = etaData[k].lunarPhase.getRef();
                xVector[7]  = etaData[k].lunarAltitude.getRef();
                xVector[8]  = LP2LI(etaData[k].lunarPhase.getRef());
                const bool goodInsert = countStats->insert(xVector,
                                                           etaData[k].observed.getRef(),
                                                           etaData[k].discovered.getRef());
                if (!goodInsert) {
                    // KEEP THIS!
                    ++excluded;
                    //
                    if (0) {
                        ORSA_DEBUG("excluded: V=%4.1f U=%5.1f AM=%5.2f GB=%+3.0f GL=%+4.0f obs: %i",
                                   etaData[k].V.getRef(),
                                   orsa::FromUnits(etaData[k].apparentVelocity.getRef()*orsa::radToArcsec(),orsa::Unit::HOUR),
                                   etaData[k].airMass.getRef(),
                                   orsa::radToDeg()*etaData[k].galacticLatitude.getRef(),
                                   orsa::radToDeg()*etaData[k].galacticLongitude.getRef(),
                                   etaData[k].observed.getRef());
                    }
                }
            }
            ORSA_DEBUG("etaData: excluded %i out of %i [%.2f\%]",excluded,etaData.size(),(100.0*excluded)/etaData.size());
        }
        
        EfficiencyMultifit::DataStorage data;
        
        {
            const CountStats::DataType & csData = countStats->getData();
            CountStats::DataType::const_iterator it = csData.begin();
            while (it != csData.end()) {
                const CountStatsElement * cs = (*it).second.get();
                if (cs==0) {
                    ++it;
                    continue;
                }
                if (cs->Ntot!=0) {
                    // const double      eta = (double)(cs->Nobs)/(double)(cs->Ntot);
                    // const double sigmaEta = (double)(sqrt(cs->Nobs+1))/(double)(cs->Ntot); // Poisson counting statistics; using Nobs+1 instead of Nobs to have positive sigma even when Nobs=0
                    // const double sigmaEta = (double)(sqrt(cs->Ntot))/(double)(cs->Ntot);
                    // #warning fix sigmaEta!
                    // 
                    // p = probability of binomial distribution
                    // const double      eta = p*cs->Ntot;
                    const double      eta = (double)(cs->Nobs)/(double)(cs->Ntot);
                    const double        p = (double)(cs->Nobs+1)/(double)(cs->Ntot+2);
                    const double sigmaEta = sqrt(p*(1-p)/cs->Ntot);
                
                    const std::vector<double> xVector = countStats->binCenterVector((*it).first); 
	
                    EfficiencyMultifit::DataElement el;
                    //
                    el.V =xVector[0];
                    el.U =xVector[1];		 
                    // el.SE=xVector[ ];	 
                    // el.LE=xVector[ ];
                    el.AM=xVector[2];
                    el.GL=xVector[3];
                    el.GB=xVector[4];
                    // el.EL=xVector[ ];
                    // el.EB=xVector[ ];
                    el.AZ=xVector[5];
                    el.SA=xVector[6];
                    // el.LP=xVector[ ];
                    el.LA=xVector[7];
                    el.LI=xVector[8];                        
                    //
                    el.eta=eta;
                    el.sigmaEta=sigmaEta;
                    //
                    el.Nobs=cs->Nobs;
                    el.Ndsc=cs->Ndsc;
                    el.Ntot=cs->Ntot;
                    //
                    data.push_back(el);
                }
                ++it;
            }
            
            ORSA_DEBUG("data.size(): %i",data.size());
        }
        
        {
            const FitFileDataElement & e = ffd[fileID];
            
            for (unsigned int k=0; k<data.size(); ++k) {
                
                /* const double eta = SkyCoverage::eta(data[k].V.getRef(),
                   e.V_limit,
                   e.eta0_V,
                   e.V0,
                   e.c_V,
                   e.w_V,
                   data[k].U.getRef(),
                   e.U_limit,
                   e.w_U,
                   data[k].AM.getRef(),
                   e.peak_AM,
                   e.scale_AM,
                   e.shape_AM,
                   data[k].GB.getRef(),
                   e.drop_GB,
                   e.scale_GB,
                   e.center_GB,
                   data[k].GL.getRef(),
                   e.scale_GL,
                   e.shape_GL,
                   data[k].SA.getRef(),
                   e.peak_SA,
                   e.scale_SA,
                   e.shape_SA,
                   data[k].LA.getRef(),
                   data[k].LI.getRef(),
                   e.LA_LI_limit_const,
                   e.LA_LI_limit_linear,
                   e.LA_LI_w_const,
                   e.LA_LI_w_linear);
                */
                //
                const double eta = SkyCoverage::eta(data[k].V.getRef(),
                                                    e.V_limit,
                                                    e.eta0_V,
                                                    e.V0,
                                                    e.c_V,
                                                    e.w_V,
                                                    data[k].U.getRef(),
                                                    e.U_limit,
                                                    e.w_U,
                                                    data[k].AM.getRef(),
                                                    e.peak_AM,
                                                    e.scale_AM,
                                                    e.shape_AM,
                                                    data[k].GB.getRef(),
                                                    e.drop_GB,
                                                    e.scale_GB,
                                                    data[k].GL.getRef(),
                                                    e.scale_GL,
                                                    e.shape_GL);
                
                const double nominal_eta_V = SkyCoverage::nominal_eta_V(data[k].V.getRef(),
                                                                        e.V_limit,
                                                                        e.eta0_V,
                                                                        e.V0,
                                                                        e.c_V,
                                                                        e.w_V);
                
                const double nominal_eta_U = SkyCoverage::nominal_eta_U(data[k].U.getRef(),
                                                                        e.U_limit,
                                                                        e.w_U);
                
                const double nominal_eta_AM = SkyCoverage::nominal_eta_AM(data[k].AM.getRef(),
                                                                          e.peak_AM,
                                                                          e.scale_AM,
                                                                          e.shape_AM);
                
                const double nominal_eta_GB = SkyCoverage::nominal_eta_GB_GL(data[k].GB.getRef(),
                                                                             e.drop_GB,
                                                                             e.scale_GB,
                                                                             0.0,
                                                                             e.scale_GL,
                                                                             0.0);
                
                const double nominal_eta_GL = SkyCoverage::nominal_eta_GB_GL(0.0,
                                                                             e.drop_GB,
                                                                             e.scale_GB,
                                                                             data[k].GL.getRef(),
                                                                             e.scale_GL,
                                                                             e.shape_GL);
                
                /* const double nominal_eta_SA = SkyCoverage::nominal_eta_SA(data[k].SA.getRef(),
                   e.peak_SA,
                   e.scale_SA,
                   e.shape_SA);
                */
                //
                const double nominal_eta_SA = 1.0;
                    
                // at this point, small sigmas can become very large because divided by a very small eta's
                if (data[k].sigmaEta.getRef() > 0) {
                    if (eta != 0) {
                        if (nominal_eta_V != 0) {
                            histo_V.insert(data[k].V.getRef(),
                                           data[k].eta.getRef()*nominal_eta_V/eta,
                                           orsa::square(eta/(nominal_eta_V*data[k].sigmaEta.getRef())));
                        }
                        if (nominal_eta_U != 0) {
                            histo_U.insert(data[k].U.getRef(),
                                           data[k].eta.getRef()*nominal_eta_U/eta,
                                           orsa::square(eta/(nominal_eta_U*data[k].sigmaEta.getRef())));
                            histo_U_more.insert(data[k].U.getRef(),
                                                data[k].eta.getRef()*nominal_eta_U/eta,
                                                orsa::square(eta/(nominal_eta_U*data[k].sigmaEta.getRef())));
                        }
                        /* histo_SE.insert(data[k].SE.getRef(),
                           data[k].eta.getRef()/eta,
                           orsa::square(eta/(data[k].sigmaEta.getRef())));
                        */
                        if (nominal_eta_AM != 0) {
                            histo_AM.insert(data[k].AM.getRef(),
                                            data[k].eta.getRef()*nominal_eta_AM/eta,
                                            orsa::square(eta/(nominal_eta_AM*data[k].sigmaEta.getRef())));
                        }
                        if (nominal_eta_GB != 0) {
                            histo_GB.insert(data[k].GB.getRef(),
                                            data[k].eta.getRef()*nominal_eta_GB/eta,
                                            orsa::square(eta/(nominal_eta_GB*data[k].sigmaEta.getRef())));
                        }
                        if (nominal_eta_GL != 0) {
                            histo_GL.insert(data[k].GL.getRef(),
                                            data[k].eta.getRef()*nominal_eta_GL/eta,
                                            orsa::square(eta/(nominal_eta_GL*data[k].sigmaEta.getRef())));
                        }
                        histo_AZ.insert(data[k].AZ.getRef(),
                                        data[k].eta.getRef()/eta,
                                        orsa::square(eta/(data[k].sigmaEta.getRef())));
                        if (nominal_eta_SA != 0) {
                            histo_SA.insert(data[k].SA.getRef(),
                                            data[k].eta.getRef()*nominal_eta_SA/eta,
                                            orsa::square(eta/(nominal_eta_SA*data[k].sigmaEta.getRef())));
                        }
                        {
                            const double nominal_eta_galactic =
                                SkyCoverage::nominal_eta_GB_GL(data[k].GB.getRef(),
                                                               e.drop_GB,
                                                               e.scale_GB,
                                                               data[k].GL.getRef(),
                                                               e.scale_GL,
                                                               e.shape_GL);  
                            if (nominal_eta_galactic != 0) {
                                histo2D_galactic.insert(data[k].GB.getRef(),
                                                        data[k].GL.getRef(),
                                                        data[k].eta.getRef()*nominal_eta_galactic/eta,
                                                        orsa::square(eta/(nominal_eta_galactic*data[k].sigmaEta.getRef())));
                                // this tracks the fitting function only, not the data!
                                histo2D_galactic_eta.insert(data[k].GB.getRef(),
                                                            data[k].GL.getRef(),
                                                            nominal_eta_galactic,
                                                            1.0); // OK to use constant sigma for fitting function?
                            }
                        }
                        
                        {
                            const double nominal_eta_lunar = 1.0;
                            /* SkyCoverage::nominal_eta_LA_LI(data[k].LA.getRef(),
                               data[k].LI.getRef(),
                               e.LA_LI_limit_const,
                               e.LA_LI_limit_linear,
                               e.LA_LI_w_const,
                               e.LA_LI_w_linear);
                            */
                            if (nominal_eta_lunar != 0) {
                                histo2D_lunar.insert(data[k].LA.getRef(),
                                                     data[k].LI.getRef(),
                                                     data[k].eta.getRef()*nominal_eta_lunar/eta,
                                                     orsa::square(eta/(nominal_eta_lunar*data[k].sigmaEta.getRef())));
                                // this tracks the fitting function only, not the data!
                                histo2D_lunar_eta.insert(data[k].LA.getRef(),
                                                         data[k].LI.getRef(),
                                                         nominal_eta_lunar,
                                                         1.0); // OK to use constant sigma for fitting function?
                            }
                        }
                    }
                }      
            }
        }
    }
    
    // now the histo_* containers are ready to be used for plotting
    
    for (unsigned int fileID=0; fileID<numFiles; ++fileID) {
        
        const FitFileDataElement & e = ffd[fileID];
        
        char filename[1024];
        sprintf(filename,"%s.fit.Vxx.dat",basename[fileID].c_str());
        ORSA_DEBUG("writing file [%s]",filename);
        FILE * fp = fopen(filename,"w");
        double V01=0.0;
        double V05=0.0;
        double V10=0.0;
        double V50=0.0;
        double V90=0.0;
        double V95=0.0;
        double V99=0.0;
        // 
        double V=24.0;
        double old_eta_V = PlotUtil_fun_V(V,e);
        while (V>=14.0) {
            V -= 0.01;
            const double eta_V = PlotUtil_fun_V(V,e);
            if (eta_V>=0.01 && old_eta_V<=0.01) V01=V;
            if (eta_V>=0.05 && old_eta_V<=0.05) V05=V;
            if (eta_V>=0.10 && old_eta_V<=0.10) V10=V;
            if (eta_V>=0.50 && old_eta_V<=0.50) V50=V;
            if (eta_V>=0.90 && old_eta_V<=0.90) V90=V;
            if (eta_V>=0.95 && old_eta_V<=0.95) V95=V;
            if (eta_V>=0.99 && old_eta_V<=0.99) V99=V;
            old_eta_V=eta_V;
        }
        // check for consinstency
        bool doPrint=true;
        if ( (e.eta0_V>0.99) && (V99>=V95) ) doPrint=false;
        if ( (e.eta0_V>0.95) && (V95>=V90) ) doPrint=false;
        if ( (e.eta0_V>0.90) && (V90>=V50) ) doPrint=false;
        if ( (e.eta0_V>0.50) && (V50>=V10) ) doPrint=false;
        if ( (e.eta0_V>0.10) && (V10>=V05) ) doPrint=false;
        if ( (e.eta0_V>0.05) && (V05>=V01) ) doPrint=false;
        if (doPrint) {
            fprintf(fp,"%.3f %.6f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n",
                    e.JD,
                    e.year,
                    V01,V05,V10,V50,V90,V95,V99);
        }
        fclose(fp);
    }
    
    for (unsigned int fileID=0; fileID<numFiles; ++fileID) {
        
        const FitFileDataElement & e = ffd[fileID];
        
        char filename[1024];
        sprintf(filename,"%s.fit.Uxx.dat",basename[fileID].c_str());
        ORSA_DEBUG("writing file [%s]",filename);
        FILE * fp = fopen(filename,"w");
        const double arcsecPerHour = orsa::FromUnits(orsa::arcsecToRad(),orsa::Unit::HOUR,-1);
        double U01=0.0;
        double U05=0.0;
        double U10=0.0;
        double U50=0.0;
        double U90=0.0;
        double U95=0.0;
        double U99=0.0;
        // 
        double U=0.0;
        // double old_eta_U = SkyCoverage::nominal_eta_U(U,U_limit,w_U);
        double old_eta_U = PlotUtil_fun_U(U,e);
        while (U<=100.0*arcsecPerHour) {
            U += 0.01*arcsecPerHour;
            const double eta_U = PlotUtil_fun_U(U,e);
            if (eta_U>=0.01 && old_eta_U<=0.01) U01=U;
            if (eta_U>=0.05 && old_eta_U<=0.05) U05=U;
            if (eta_U>=0.10 && old_eta_U<=0.10) U10=U;
            if (eta_U>=0.50 && old_eta_U<=0.50) U50=U;
            if (eta_U>=0.90 && old_eta_U<=0.90) U90=U;
            if (eta_U>=0.95 && old_eta_U<=0.95) U95=U;
            if (eta_U>=0.99 && old_eta_U<=0.99) U99=U;
            old_eta_U=eta_U;
        }
        // check for consinstency
        if ( (U01<U05) &&
             (U05<U10) &&
             (U10<U50) &&
             (U50<U90) &&
             (U90<U95) &&
             (U95<U99) ) { 
            fprintf(fp,"%.3f %.6f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n",
                    e.JD,
                    e.year,
                    U01/arcsecPerHour,
                    U05/arcsecPerHour,
                    U10/arcsecPerHour,
                    U50/arcsecPerHour,
                    U90/arcsecPerHour,
                    U95/arcsecPerHour,
                    U99/arcsecPerHour);
        }
        fclose(fp);
    }
    
    for (unsigned int fileID=0; fileID<numFiles; ++fileID) {
        char filename[1024];
        sprintf(filename,"%s.galactic.dat",basename[fileID].c_str());
        ORSA_DEBUG("writing file [%s]",filename);
        FILE * fp = fopen(filename,"w");
        for (unsigned int j=0; j<histo2D_galactic.getData().size(); ++j) {
            for (unsigned int k=0; k<histo2D_galactic.getData()[j].size(); ++k) {
                const EfficiencyStatistics2D * es = histo2D_galactic.getData()[j][k];
                if (es->entries()>0) {
                    gmp_fprintf(fp,
                                "%+5.1f %+6.1f %5.3f %5.3f %5.3f %Zi\n",
                                orsa::radToDeg()*es->centerX,
                                orsa::radToDeg()*es->centerY,
                                es->average(),
                                es->standardDeviation(),
                                histo2D_galactic_eta.getData()[j][k]->average(),
                                es->entries().get_mpz_t());
                }
            }
        }        
        fclose(fp);
    }

    // 1 deg by 1 deg mesh of fiitting function, for contour plots
    for (unsigned int fileID=0; fileID<numFiles; ++fileID) {
        char filename[1024];
        sprintf(filename,"%s.fit.galactic.contour.dat",basename[fileID].c_str());
        ORSA_DEBUG("writing file [%s]",filename);
        FILE * fp = fopen(filename,"w");
        const FitFileDataElement & e = ffd[fileID];
        for (int gl=-180; gl<=180; ++gl) {
            for (int gb=-90; gb<=90; ++gb) {
                const double nominal_eta_galactic =
                    SkyCoverage::nominal_eta_GB_GL(gb*orsa::degToRad(),
                                                   e.drop_GB,
                                                   e.scale_GB,
                                                   gl*orsa::degToRad(),
                                                   e.scale_GL,
                                                   e.shape_GL);
                gmp_fprintf(fp,
                            "%+3i %+2i %5.3f\n",
                            gl,
                            gb,
                            nominal_eta_galactic);
            }
        }
        fclose(fp);
    }
    
    for (unsigned int fileID=0; fileID<numFiles; ++fileID) {
        char filename[1024];
        sprintf(filename,"%s.lunar.dat",basename[fileID].c_str());
        ORSA_DEBUG("writing file [%s]",filename);
        FILE * fp = fopen(filename,"w");
        for (unsigned int j=0; j<histo2D_lunar.getData().size(); ++j) {
            for (unsigned int k=0; k<histo2D_lunar.getData()[j].size(); ++k) {
                const EfficiencyStatistics2D * es = histo2D_lunar.getData()[j][k];
                if (es->entries()>0) {
                    gmp_fprintf(fp,
                                "%+5.1f %+6.3f %5.3f %5.3f %5.3f %Zi\n",
                                orsa::radToDeg()*es->centerX,
                                es->centerY,
                                es->average(),
                                es->standardDeviation(),
                                histo2D_lunar_eta.getData()[j][k]->average(),
                                es->entries().get_mpz_t());
                }
            }
        }        
        fclose(fp);
    }
    
    // 1 deg by 1 deg mesh of fiitting function, for contour plots
    for (unsigned int fileID=0; fileID<numFiles; ++fileID) {
        char filename[1024];
        sprintf(filename,"%s.fit.lunar.contour.dat",basename[fileID].c_str());
        ORSA_DEBUG("writing file [%s]",filename);
        FILE * fp = fopen(filename,"w");
        // const FitFileDataElement & e = ffd[fileID];
        for (int li=0; li<=100; ++li) { 
            for (int la=-90; la<=90; ++la) {
                const double nominal_eta_lunar = 1.0;
                /* SkyCoverage::nominal_eta_LA_LI(la*orsa::degToRad(),
                   li*0.01,
                   e.LA_LI_limit_const,
                   e.LA_LI_limit_linear,
                   e.LA_LI_w_const,
                   e.LA_LI_w_linear);
                */
                gmp_fprintf(fp,
                            "%+3i %+2i %5.3f\n",
                            li,
                            la,
                            nominal_eta_lunar);
            }
        }
        fclose(fp);
    }
    
    std::string bigbasename;
    if (numFiles==1) {
        bigbasename=basename[0];
    } else {
        bigbasename="multi";
        for (unsigned int fileID=0; fileID<numFiles; ++fileID) {
            bigbasename += "_";
            bigbasename += basename[fileID];
        }
    }
    ORSA_DEBUG("bigbasename: [%s]",bigbasename.c_str());
    
    // #warning remove this exit() call to generate PDF!
    // exit(0);
        
    // dislin
    page(3000,4000);
    // page(2390,3092); // approx letter size, since setpag("USAP") does not seem to work properly
    // setpag("USAP"); // USAP = US letter portrait
    pagmod("PORT");
    //
    // metafl("PNG");
    // metafl("PS");
    metafl("PDF");
    
    // output file name
    char plotFilename[1024];
    // sprintf(plotFilename,"%s.fit.png",bigbasename.c_str());
    // sprintf(plotFilename,"%s.fit.ps",bigbasename.c_str());
    sprintf(plotFilename,"%s.fit.pdf",bigbasename.c_str());
    setfil(plotFilename);
    // new files overwrite old ones
    filmod("DELETE");
    
    disini();
    simplx();
    // helve();
    // psfont("AvantGarde-Book");
    // height(20); // text height
    color("fore");
    // paghdr("preceding","following",4,0); // page header
    pagera(); // border around page

    hwmode("ON","LINE");
    // hwmode("ON","SHADING");
    // linwid(1);    
    // barmod("VARIABLE","WIDTH");
    // barwth(4.0);
    texmod("ON"); // TeX text
        
    // header stuff
    {
#warning RE-INSERT this, and use longer header for multiple input files
        /* const orsa::Time epoch = orsaSolarSystem::julianToTime(JD);
           int y,m,d; double fd; orsaSolarSystem::gregorDay(epoch,y,m,d,fd);
           penwid(0.5);
           height(30); // text height
           char line[1024];
           sprintf(line,"[%s] %4i/%02i/%02i",jobID,y,m,d);
           messag(line,250,100);
           //
           orsaSolarSystem::gregorDay(orsaSolarSystem::now(),y,m,d,fd);
           sprintf(line,"Processed on %4i/%02i/%02i by P. Tricarico <tricaric@psi.edu>",y,m,d);
           messag(line,250,150);
        */
    }
    
    const double etaMin  = 0.0;
    const double etaMax  = 1.4;
    const double etaStep = 0.2;
    const int etaDigits  = 1;
    const int etaTicks   = 2;
    const std::string etaLabel="detection efficiency";
    // const std::string etaLabel="$\\eta$";
        
    // formatting
    const int px = 1100;
    const int py =  450;
    // spacing
    const int sx =  250;
    const int sy =  250;
    // plot size + spacing
    const int lx = px+sx;
    const int ly = py+sy;
    // initial position
    const int nx0 = 250;
    const int ny0 = 900;
        
    // first row of plots
    int nx=nx0;
    int ny=ny0;
        
    PlotUtil((&PlotUtil_fun_V),
             ffd,
             histo_V,
             1.0,
             nx,ny,px,py,
             14,24,1.0,0,1,"apparent magnitude",
             etaMin,etaMax,etaStep,etaDigits,etaTicks,etaLabel);
        
    ny += ly;
        
    PlotUtil((&PlotUtil_fun_AM),
             ffd,
             histo_AM,
             1.0,
             nx,ny,px,py,
             1.0,3.0,0.2,1,1,"airmass",
             etaMin,etaMax,etaStep,etaDigits,etaTicks,etaLabel);
        
    ny += ly;
        
    PlotUtil((&PlotUtil_fun_one),
             ffd,
             histo_AZ,
             orsa::radToDeg(),
             nx,ny,px,py,
             0,360,60,0,3,"azimuth [deg]",
             etaMin,etaMax,etaStep,etaDigits,etaTicks,etaLabel);
        
    ny += ly;

    /* 
       PlotUtil((&PlotUtil_fun_SA),
       ffd,
       histo_SA,
       orsa::radToDeg(),
       nx,ny,px,py,
       -90,90,30,0,3,"solar altitude [deg]",
       etaMin,etaMax,etaStep,etaDigits,etaTicks,etaLabel);
    */
    //
    PlotUtil((&PlotUtil_fun_one),
             ffd,
             histo_SA,
             orsa::radToDeg(),
             nx,ny,px,py,
             -90,90,30,0,3,"solar altitude [deg]",
             etaMin,etaMax,etaStep,etaDigits,etaTicks,etaLabel);
        
    ny += ly;
        
    /* PlotUtil((&PlotUtil_fun_one),
       ffd,
       histo_SE,
       orsa::radToDeg(),
       nx,ny,px,py,
       0,180,30,0,3,"solar elongation [deg]",
       etaMin,etaMax,etaStep,etaDigits,etaTicks,etaLabel);
    */
    
    // second row of plots
    nx += lx;
    ny  = ny0;
    
    PlotUtil((&PlotUtil_fun_U),
             ffd,
             histo_U,
             orsa::FromUnits(orsa::radToArcsec(),orsa::Unit::HOUR),
             nx,ny,px,py,
             0,100,10,0,1,"apparent velocity [arcsec/hour]",
             etaMin,etaMax,etaStep,etaDigits,etaTicks,etaLabel);
        
    ny += ly;
    
    PlotUtil((&PlotUtil_fun_U),
             ffd,
             histo_U_more,
             orsa::FromUnits(orsa::radToArcsec(),orsa::Unit::HOUR),
             nx,ny,px,py,
             0,600,100,0,2,"apparent velocity [arcsec/hour]",
             etaMin,etaMax,etaStep,etaDigits,etaTicks,etaLabel);
    
    ny += ly;
        
    PlotUtil((&PlotUtil_fun_GB),
             ffd,
             histo_GB,
             orsa::radToDeg(),
             nx,ny,px,py,
             -90,90,30,0,3,"galactic latitude [deg]",
             etaMin,etaMax,etaStep,etaDigits,etaTicks,etaLabel);
        
    ny += ly;
    
#warning re-introduce the lunar phase plot somewhere else...
    /* PlotUtil((&PlotUtil_fun_one),
       ffd,
       histo_LP,
       orsa::radToDeg(),
       nx,ny,px,py,
       0,180,30,0,3,"lunar phase [deg]",
       etaMin,etaMax,etaStep,etaDigits,etaTicks,etaLabel);
    */
    //
    /* PlotUtil((&PlotUtil_fun_one),
       ffd,
       histo_GL,
       orsa::radToDeg(),
       nx,ny,px,py,
       -180,180,60,0,3,"galactic longitude [deg]",
       etaMin,etaMax,etaStep,etaDigits,etaTicks,etaLabel);
    */
    //
    PlotUtil((&PlotUtil_fun_GL),
             ffd,
             histo_GL,
             orsa::radToDeg(),
             nx,ny,px,py,
             -180,180,60,0,3,"galactic longitude [deg]",
             etaMin,etaMax,etaStep,etaDigits,etaTicks,etaLabel);
        
    ny += ly;
        
    /* PlotUtil((&PlotUtil_fun_one),
       ffd,
       histo_LA,
       orsa::radToDeg(),
       nx,ny,px,py,
       -90,90,30,0,3,"lunar altitude [deg]",
       etaMin,etaMax,etaStep,etaDigits,etaTicks,etaLabel);
    */
    
    ny += ly;
        
    /* PlotUtil((&PlotUtil_fun_one),
       ffd,
       histo_LE,
       orsa::radToDeg(),
       nx,ny,px,py,
       0,180,30,0,3,"lunar elongation [deg]",
       etaMin,etaMax,etaStep,etaDigits,etaTicks,etaLabel);
    */
    
    disfin();
    
    exit(0);
}
