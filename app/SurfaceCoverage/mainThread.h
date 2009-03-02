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
  orsa::Cache<orsa::Double> vestaMass;
  // orsa::Cache<ComboMassDistribution::MassDistributionType> vestaMassDistribution;
  // orsa::Cache<ComboShapeModel::ShapeModelType>             vestaShapeModel;
  orsa::Cache<orsa::Double> vestaPeriod;
  orsa::Cache<orsa::Double> vestaPoleEclipticLatitude;
  orsa::Cache<orsa::Double> vestaPoleEclipticLongitude;
 public:
  orsa::Cache<orsa::Time>   orbitEpoch;
  orsa::Cache<orsa::Double> orbitRadius;
  orsa::Cache<orsa::Double> orbitInclination;
  orsa::Cache<orsa::Double> orbitPhase;
 public:
  orsa::Cache<orsa::Double> cameraParallelFOV;
  orsa::Cache<orsa::Double> cameraOrthogonalFOV;
  orsa::Cache<orsa::Time>   cameraInterval;
  orsa::Cache<orsa::Double> cameraMinimumI;
  orsa::Cache<orsa::Double> cameraMaximumI;
  orsa::Cache<orsa::Double> cameraMaximumE;
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
  std::vector<orsa::Double>    faceArea_out;
  orsa::TriShape::FaceVector   faceVector_out;
  orsa::TriShape::VertexVector vertexVector_out;
  std::vector< std::vector<orsa::Time> > vertexCoverage_out;
  // test
 public:
  // std::vector<unsigned int> allFaceAtOnce_out;
  std::vector< std::vector<orsa::Time> > allFaceAtOnce_out;
};

#endif // __VESTA_MAIN_THREAD_H__
