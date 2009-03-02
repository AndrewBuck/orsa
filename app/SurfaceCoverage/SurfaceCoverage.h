#ifndef __SURFACE_COVERAGE_H__
#define __SURFACE_COVERAGE_H__

#include <orsaTBB/malloc.h>

#include <QComboBox>
#include <QLineEdit> 
#include <QProgressBar>
#include <QPushButton>
#include <QStatusBar>
#include <QVariant>
#include <QWidget>

#include <osgGA/NodeTrackerManipulator>

#include <orsa/debug.h>

// #include "plot2D.h"

class MainThread;
class Plot2D;
class PlotTimeline;

/* 
   class ComboShapeModel : public QComboBox {
   public:
   ComboShapeModel(QWidget * parent = 0) : 
   QComboBox(parent) {
   insertItem(smt_ellipsoid,tr("ellipsoid"));
   insertItem(smt_thomas,   tr("Thomas et al. (1997)"));
   }
   public:
   enum ShapeModelType {
   smt_ellipsoid = 0,
   smt_thomas    = 1
   };
   };
*/

class SurfaceCoverageWidget : public QWidget {
  
  Q_OBJECT;
  
 public:
  SurfaceCoverageWidget();
  
 public:
  virtual ~SurfaceCoverageWidget() {
    ORSA_DEBUG("~RendezvousWithVestaWidget() called...");
  }
  
 public slots:
  void reset();
  
 public slots:
  virtual void run() const;
  
 public slots:
  void abort();
  
 public slots:
  void progressChanged(int) const;
  
 public slots:
  void runFinished() const;
  
 public slots:
  void showAboutDialog();
  
 public slots:
  void computePeriodsRatio(const QString &);
  
 protected:
  void setInputEnabled(const bool) const;
  
 protected:
  void closeEvent(QCloseEvent *); 
  
 protected:
  Plot2D * plot2D;
  
 protected:
  PlotTimeline * plotTimeline;
  
 protected:
  QLineEdit * lineVestaMass;
  // ComboShapeModel * comboShapeModel;
  QLineEdit * lineVestaPeriod;
  QLineEdit * lineVestaPoleEclipticLatitude;
  QLineEdit * lineVestaPoleEclipticLongitude; 
 protected:
  QLineEdit * lineOrbitEpoch;  
  QLineEdit * lineOrbitRadius;
  QLineEdit * lineOrbitPeriodRatio;
  QLineEdit * lineOrbitInclination;
  QLineEdit * lineOrbitPhase;
 protected:
  QLineEdit * lineCameraParallelFOV;
  QLineEdit * lineCameraOrthogonalFOV;
  QLineEdit * lineCameraInterval;
  QLineEdit * lineCameraMinimumI;
  QLineEdit * lineCameraMaximumI;
  QLineEdit * lineCameraMaximumE;
 protected:
  QLineEdit * lineRunDuration;
  // QLineEdit * lineRunSPICEFile;
  // QPushButton * buttonBrowseSPICE;
  // QLineEdit * lineRunASCIIFile;
  // QPushButton * buttonBrowseASCII;
  QPushButton * buttonRun;
  QProgressBar * progressRun;
  QPushButton * buttonAbort;
  
 protected:
  QPushButton * buttonAbout;
  
 protected:  
  QStatusBar * statusBar;
  
 protected slots:
  void lineEdit_textChanged(const QString &);
  
 protected:
  MainThread * mainThread;
};

// main purpose: disable the rotation after the release of the mouse button
class VestaNodeTrackerManipulator : public osgGA::NodeTrackerManipulator {
 public:
  bool handle(const osgGA::GUIEventAdapter & ea,
	      osgGA::GUIActionAdapter      & us) {
    if (ea.getEventType() == osgGA::GUIEventAdapter::RELEASE) {
      osg::ref_ptr<osgGA::GUIEventAdapter> mod_ea = new osgGA::GUIEventAdapter(ea);
      mod_ea->setButtonMask(99);
      return osgGA::NodeTrackerManipulator::handle((*mod_ea.get()),us);
    } else {
      return osgGA::NodeTrackerManipulator::handle(ea,us);
    }
  }
};

#endif // __SURFACE_COVERAGE_H__
