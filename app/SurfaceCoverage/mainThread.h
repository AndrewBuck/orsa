#ifndef __VESTA_MAIN_THREAD_H__
#define __VESTA_MAIN_THREAD_H__

#include <string>

#include <QThread>

#include <orsa/bodygroup.h>
#include <orsa/double.h>
#include <orsa/cache.h>
#include <orsa/integrator.h>
#include <orsa/integrator_leapfrog.h>
#include <orsa/integrator_radau.h>
#include <orsa/shape.h>

#include "SurfaceCoverage.h"

class MainThread : public QThread {
  
  Q_OBJECT;

 public:
  MainThread();
 public:
  virtual ~MainThread() { }

 protected:
  void run();
  
 public slots:
  void abort();
 protected:
  bool aborted;
  
 signals:
  void progress(int value) const;

 public slots:
  void reset();

 public:
  orsa::Cache<double> vestaMass;
  // orsa::Cache<ComboMassDistribution::MassDistributionType> vestaMassDistribution;
  // orsa::Cache<ComboShapeModel::ShapeModelType>             vestaShapeModel;
  orsa::Cache<double> vestaPeriod;
  orsa::Cache<double> vestaPoleEclipticLatitude;
  orsa::Cache<double> vestaPoleEclipticLongitude;
 public:
  orsa::Cache<orsa::Time>   orbitEpoch;
  orsa::Cache<double> orbitRadius;
  orsa::Cache<double> orbitInclination;
  orsa::Cache<double> orbitPhase;
 public:
  orsa::Cache<double> cameraParallelFOV;
  orsa::Cache<double> cameraOrthogonalFOV;
  orsa::Cache<orsa::Time>   cameraInterval;
  orsa::Cache<double> cameraMinimumI;
  orsa::Cache<double> cameraMaximumI;
  orsa::Cache<double> cameraMaximumE;
 public:
  orsa::Cache<orsa::Time>   runDuration;
  orsa::Cache<std::string>  outputSPICEFile;
  orsa::Cache<std::string>  outputASCIIFile;
  
 protected:
  // osg::ref_ptr<CustomIntegrator> customIntegrator;
  
 public:
  osg::ref_ptr<orsa::BodyGroup> bg;
  
  // output
 public:	
  std::vector< std::vector<double> > latitude, longitude;
  std::vector<double>    faceArea_out;
  orsa::TriShape::FaceVector   faceVector_out;
  orsa::TriShape::VertexVector vertexVector_out;
  std::vector< std::vector<orsa::Time> > vertexCoverage_out;
  // test
 public:
  // std::vector<unsigned int> allFaceAtOnce_out;
  std::vector< std::vector<orsa::Time> > allFaceAtOnce_out;
};

#endif // __VESTA_MAIN_THREAD_H__
