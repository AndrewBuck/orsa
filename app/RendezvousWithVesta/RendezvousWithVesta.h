#ifndef __RENDEZVOUS_WITH_VESTA_H__
#define __RENDEZVOUS_WITH_VESTA_H__

#include <QComboBox>
#include <QLineEdit> 
#include <QProgressBar>
#include <QPushButton>
#include <QStatusBar>
#include <QVariant>
#include <QWidget>

#include <osgGA/NodeTrackerManipulator>

#include <orsa/debug.h>

// #include "mainThread.h"
class MainThread;

class ComboMassDistribution : public QComboBox {
 public:
  ComboMassDistribution(QWidget * parent = 0) : 
    QComboBox(parent) {
    insertItem(mdt_uniform,tr("uniform density"));
    insertItem(mdt_core,   tr("core + mantle")); 
  }
 public:
  enum MassDistributionType {
    mdt_uniform = 0,
    mdt_core    = 1
  };
};


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


class RendezvousWithVestaWidget : public QWidget {
  
  Q_OBJECT;
  
 public:
  RendezvousWithVestaWidget();
  
 public:
  virtual ~RendezvousWithVestaWidget() {
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
  void browseSPICE();
  void browseASCII();
  
 public slots:
  void showAboutDialog();
  
 protected:
  void setInputEnabled(const bool) const;
  
 protected:
  void closeEvent(QCloseEvent *); 
  
 protected:
  QLineEdit * lineVestaMass;
  ComboMassDistribution * comboMassDistribution;
  ComboShapeModel * comboShapeModel;
  QLineEdit * lineVestaPeriod;
  QLineEdit * lineVestaPoleEclipticLatitude;
  QLineEdit * lineVestaPoleEclipticLongitude; 
 protected:
  QLineEdit * lineOrbitEpoch;  
  QLineEdit * lineOrbitRadius;
  QLineEdit * lineOrbitInclination;
  QLineEdit * lineOrbitPhase;
 protected:
  QLineEdit * lineRunDuration;
  QLineEdit * lineRunSPICEFile;
  QPushButton * buttonBrowseSPICE;
  QLineEdit * lineRunASCIIFile;
  QPushButton * buttonBrowseASCII;
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

#endif // __RENDEZVOUS_WITH_VESTA_H__
