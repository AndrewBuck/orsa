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

#include "RendezvousWithVesta.h"

class MainThread;

// class CustomIntegrator : public QObject, public orsa::IntegratorLeapFrog {
class CustomIntegrator : public QObject, public orsa::IntegratorRadau {
  
  Q_OBJECT;

 public:
  CustomIntegrator(const MainThread * mt);

 signals:
  void progress(int value) const;

 public:
  const MainThread * mainThread;
 protected:
  void singleStepDone(orsa::BodyGroup  *,
		      const orsa::Time &,
		      const orsa::Time &,
		      orsa::Time       &) const;
};

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
  // protected:
  // bool aborted;

 signals:
  void progress(int value) const;

 public slots:
  void reset();

 public:
  orsa::Cache<double> vestaMass;
  orsa::Cache<ComboMassDistribution::MassDistributionType> vestaMassDistribution;
  orsa::Cache<ComboShapeModel::ShapeModelType>             vestaShapeModel;
  orsa::Cache<double> vestaPeriod;
  orsa::Cache<double> vestaPoleEclipticLatitude;
  orsa::Cache<double> vestaPoleEclipticLongitude;
 public:
  orsa::Cache<orsa::Time>   orbitEpoch;
  orsa::Cache<double> orbitRadius;
  orsa::Cache<double> orbitInclination;
  orsa::Cache<double> orbitPhase;
 public:
  orsa::Cache<orsa::Time>   runDuration;
  orsa::Cache<std::string>  outputSPICEFile;
  orsa::Cache<std::string>  outputASCIIFile;

 protected:
  osg::ref_ptr<CustomIntegrator> customIntegrator;

 public:
  osg::ref_ptr<orsa::BodyGroup> bg;

};

#endif // __VESTA_MAIN_THREAD_H__
