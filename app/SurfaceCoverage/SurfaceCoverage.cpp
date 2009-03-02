#include "SurfaceCoverage.h"
#include "SurfaceCoverageVersion.h"

#include <orsa/orbit.h>
#include <orsa/unit.h>
#include <orsaSolarSystem/datetime.h>

#include <QComboBox> 
#include <QCoreApplication>
#include <QDoubleValidator>
#include <QFileDialog>
#include <QFrame>
#include <QGridLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QIntValidator>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QPalette>
#include <QProgressBar>
#include <QPushButton>
#include <QSizePolicy>
#include <QStatusBar>
#include <QTabWidget>
#include <QVariant>
#include <QVBoxLayout>
#include <QWidget>

#include "plot2D.h"

#include "plotTimeline.h"

#include "mainThread.h"

SurfaceCoverageWidget::SurfaceCoverageWidget() : QWidget(0) {
  
  QGridLayout * grid = new QGridLayout(this);
  
  QTabWidget * tabCoverage = new QTabWidget(this);
  
  plot2D = new Plot2D(tabCoverage);
  
  tabCoverage->addTab(plot2D, tr("2D Surface Coverage"));
  
  // const int tabCoverage3DIndex = tabCoverage->addTab(new QWidget(tabCoverage),tr("3D Surface Coverage"));
  
  QVBoxLayout * layoutParam = new QVBoxLayout;
  
  {
    // physical properties group
    
    QGroupBox * group = new QGroupBox(tr("Vesta properties"),this);
    layoutParam->addWidget(group);
    
    QGridLayout * grid = new QGridLayout(group);
    //
    // grid->setMargin(3);

    lineVestaMass                  = new QLineEdit(group);
    // comboMassDistribution          = new ComboMassDistribution(group);
    // comboShapeModel                = new ComboShapeModel(group);
    lineVestaPeriod                = new QLineEdit(group);
    lineVestaPoleEclipticLatitude  = new QLineEdit(group);
    lineVestaPoleEclipticLongitude = new QLineEdit(group);

    int localRow = 0;

    grid->addWidget(new QLabel(tr("mass [kg]"),group),localRow,0);
    grid->addWidget(lineVestaMass                          ,localRow,1);

    // ++localRow;
    // 
    // grid->addWidget(new QLabel(tr("mass distribution"),group),localRow,0);
    // grid->addWidget(comboMassDistribution                    ,localRow,1);
    
    // ++localRow;
    // 
    // grid->addWidget(new QLabel(tr("shape model"),group),localRow,0);
    // grid->addWidget(comboShapeModel                    ,localRow,1);
    
    ++localRow;

    grid->addWidget(new QLabel(tr("rotation period [hour]"),group),localRow,0);
    grid->addWidget(lineVestaPeriod                               ,localRow,1);

    ++localRow;

    grid->addWidget(new QLabel(tr("pole ecliptic latitude [deg]"),group),localRow,0);
    grid->addWidget(lineVestaPoleEclipticLatitude                       ,localRow,1);

    ++localRow;

    grid->addWidget(new QLabel(tr("pole ecliptic longitude [deg]"),group),localRow,0);
    grid->addWidget(lineVestaPoleEclipticLongitude                       ,localRow,1);

    ++localRow;

  }
  
  {
    // orbital parameters group
    
    QGroupBox * group = new QGroupBox(tr("Orbit"),this);
    layoutParam->addWidget(group);
    
    QGridLayout * grid = new QGridLayout(group);
    //
    // grid->setMargin(3);
    
    lineOrbitEpoch       = new QLineEdit(group);
    lineOrbitRadius      = new QLineEdit(group);
    lineOrbitPeriodRatio = new QLineEdit(group);
    lineOrbitInclination = new QLineEdit(group);
    lineOrbitPhase       = new QLineEdit(group);
    
    lineOrbitPeriodRatio->setReadOnly(true);
    lineOrbitPeriodRatio->setEnabled(false);
    
    int localRow = 0;
    
    grid->addWidget(new QLabel(tr("epoch [yyyy-mm-dd HH:MM:SS] [UTC]"),group),localRow,0);
    grid->addWidget(lineOrbitEpoch                                           ,localRow,1);
    
    ++localRow;
    
    grid->addWidget(new QLabel(tr("radius [km]"),group),localRow,0);
    grid->addWidget(lineOrbitRadius                    ,localRow,1);
    
    ++localRow;
    
    grid->addWidget(new QLabel(tr("periods ratio [orbital:rotational]"),group),localRow,0);
    grid->addWidget(lineOrbitPeriodRatio                                      ,localRow,1);
    
    ++localRow;
    
    grid->addWidget(new QLabel(tr("equatorial inclination [deg]"),group),localRow,0);
    grid->addWidget(lineOrbitInclination                                ,localRow,1);

    ++localRow;

    grid->addWidget(new QLabel(tr("phase angle [deg]"),group),localRow,0);
    grid->addWidget(lineOrbitPhase                           ,localRow,1);

    ++localRow;

  }
  
  {
    // camera group
    
    QGroupBox * group = new QGroupBox(tr("Camera"),this);
    layoutParam->addWidget(group);
    
    QGridLayout * grid = new QGridLayout(group);
    
    lineCameraParallelFOV   = new QLineEdit(group);
    lineCameraOrthogonalFOV = new QLineEdit(group);
    lineCameraInterval      = new QLineEdit(group);
    lineCameraMinimumI      = new QLineEdit(group);
    lineCameraMaximumI      = new QLineEdit(group);
    lineCameraMaximumE      = new QLineEdit(group);
    
    int localRow = 0;   
    
    grid->addWidget(new QLabel(tr("FOV (parallel) [deg]"),group),localRow,0);
    grid->addWidget(lineCameraParallelFOV                       ,localRow,1);
    
    ++localRow;

    grid->addWidget(new QLabel(tr("FOV (orthogonal) [deg]"),group),localRow,0);
    grid->addWidget(lineCameraOrthogonalFOV                       ,localRow,1);
    
    ++localRow;

    grid->addWidget(new QLabel(tr("Interval [s]"),group),localRow,0);
    grid->addWidget(lineCameraInterval                  ,localRow,1);
    
    ++localRow;
    
    grid->addWidget(new QLabel(tr("minimum light incidence angle [deg]"),group),localRow,0);
    grid->addWidget(lineCameraMinimumI                                         ,localRow,1);
    
    ++localRow;
    
    grid->addWidget(new QLabel(tr("maximum light incidence angle [deg]"),group),localRow,0);
    grid->addWidget(lineCameraMaximumI                                         ,localRow,1);
    
    ++localRow;
    
    grid->addWidget(new QLabel(tr("maximum surface inclination angle [deg]"),group),localRow,0);
    grid->addWidget(lineCameraMaximumE                                             ,localRow,1);
    
  }
  
  {
    // run group

    QGroupBox * group = new QGroupBox(tr("Simulation"),this);
    layoutParam->addWidget(group);

    QGridLayout * grid = new QGridLayout(group);
    //
    // grid->setMargin(3);

    lineRunDuration          = new QLineEdit(group);
    // lineRunSPICEFile         = new QLineEdit(group);
    // buttonBrowseSPICE        = new QPushButton("browse",group);
    // lineRunASCIIFile         = new QLineEdit(group);
    // buttonBrowseASCII        = new QPushButton("browse",group);
    buttonRun                = new QPushButton("run simulation",group);
    progressRun              = new QProgressBar(group);
    buttonAbort              = new QPushButton("abort",group);
    
    int localRow = 0;

    grid->addWidget(new QLabel(tr("duration [ddd:HH:MM:SS]"),group),localRow,0);
    grid->addWidget(lineRunDuration                                ,localRow,1);
    
    ++localRow;

    grid->addWidget(buttonRun  ,localRow,0);
    grid->addWidget(progressRun,localRow,1);
    grid->addWidget(buttonAbort,localRow,2);

    ++localRow;

  }
  
  layoutParam->addStretch();
  
  {
    // about group

    QGroupBox * group = new QGroupBox(tr("About"),this);
    layoutParam->addWidget(group);

    QGridLayout * grid = new QGridLayout(group);
    //
    // grid->setMargin(3);

    buttonAbout = new QPushButton("more",group);

    int localRow = 0;

    {
      QLabel * aboutLabel = new QLabel(tr("Copyright (c) 2008 Pasquale Tricarico <a href=\"mailto:tricaric@psi.edu\">tricaric@psi.edu</a>"),group);
      aboutLabel->setTextInteractionFlags(Qt::TextSelectableByMouse | Qt::TextSelectableByKeyboard | Qt::LinksAccessibleByMouse | Qt::LinksAccessibleByKeyboard);
      grid->addWidget(aboutLabel,localRow,0);
    }

    ++localRow;

    {
      QHBoxLayout * boxLayout = new QHBoxLayout;
      QLabel * versionLabel = new QLabel(programVersion(),group);
      versionLabel->setTextInteractionFlags(Qt::TextSelectableByMouse | Qt::TextSelectableByKeyboard | Qt::LinksAccessibleByMouse | Qt::LinksAccessibleByKeyboard);
      boxLayout->addWidget(versionLabel,100);
      boxLayout->addWidget(buttonAbout);
      grid->addLayout(boxLayout,localRow,0,1,1);
    }

    ++localRow;

  }


  {
    statusBar = new QStatusBar(this);
    statusBar->setSizeGripEnabled(false);
    statusBar->showMessage(tr("Ready"));
    grid->addWidget(statusBar,100,0);
  }
  

  {
    // validation
    
    QDoubleValidator * positiveValidator = new QDoubleValidator(this);
    positiveValidator->setBottom(0.0);
    
    QDoubleValidator * doubleValidator = new QDoubleValidator(this);
    
    QIntValidator * positiveIntValidator = new QIntValidator(this);
    positiveIntValidator->setBottom(0);
    
    lineVestaMass->setValidator(positiveValidator);
    lineVestaPeriod->setValidator(positiveValidator);
    lineVestaPoleEclipticLatitude->setValidator(doubleValidator);
    lineVestaPoleEclipticLongitude->setValidator(doubleValidator);
    
    lineOrbitEpoch->setInputMask("9999-99-99 99:99:99;_");
    //
    lineOrbitRadius->setValidator(positiveValidator);
    lineOrbitInclination->setValidator(positiveValidator);
    lineOrbitPhase->setValidator(doubleValidator);
    
    lineCameraParallelFOV->setValidator(positiveValidator);
    lineCameraOrthogonalFOV->setValidator(positiveValidator);
    lineCameraInterval->setValidator(positiveIntValidator);
    lineCameraMinimumI->setValidator(positiveIntValidator);
    lineCameraMaximumI->setValidator(positiveIntValidator);
    lineCameraMaximumE->setValidator(positiveIntValidator);
    
    lineRunDuration->setInputMask("999:99:99:99;_");
  }
  
  {
    // main thread
    mainThread = new MainThread;
  }
  
  {
    // connect signals and slots
    
    connect(buttonRun,
	    SIGNAL(clicked()),
	    this,
	    SLOT(run()));

    connect(buttonAbort,
	    SIGNAL(clicked()),
	    this,
	    SLOT(abort()));

    connect(buttonAbout,
	    SIGNAL(clicked()),
	    this,
	    SLOT(showAboutDialog()));

    connect(mainThread,
	    SIGNAL(progress(int)),
	    progressRun,
	    SLOT(setValue(int)));

    connect(mainThread,
	    SIGNAL(progress(int)),
	    this,
	    SLOT(progressChanged(int)));
    
    connect(mainThread,
	    SIGNAL(finished()),
	    this,
	    SLOT(runFinished()));
    
    connect(lineOrbitEpoch,
	    SIGNAL(textChanged(const QString &)),
	    this,
	    SLOT(lineEdit_textChanged(const QString &)));
    
    connect(lineRunDuration,
	    SIGNAL(textChanged(const QString &)),
	    this,
	    SLOT(lineEdit_textChanged(const QString &)));  
  }
  
  {
    // orbitPeriodsRatio signals
    
    connect(lineVestaMass,
	    SIGNAL(textChanged(const QString &)),
	    this,
	    SLOT(computePeriodsRatio(const QString &)));
 
    connect(lineVestaPeriod,
	    SIGNAL(textChanged(const QString &)),
	    this,
	    SLOT(computePeriodsRatio(const QString &)));
    
    connect(lineOrbitRadius,
	    SIGNAL(textChanged(const QString &)),
	    this,
	    SLOT(computePeriodsRatio(const QString &)));
    
  }
  
  QTabWidget * tabPlotTimeline = new QTabWidget(this);
  
  plotTimeline = new PlotTimeline(tabPlotTimeline);
  
  tabPlotTimeline->addTab(plotTimeline,tr("Timeline"));
  
  grid->addWidget(tabCoverage,0,0);
  grid->addLayout(layoutParam,0,1,2,1);
  grid->addWidget(tabPlotTimeline,1,0);
  
  {
    QSizePolicy sp = QSizePolicy(QSizePolicy::MinimumExpanding,QSizePolicy::MinimumExpanding);
    sp.setHorizontalStretch(10);
    sp.setVerticalStretch(3);
    tabCoverage->setSizePolicy(sp);
  }
  
  {
    QSizePolicy sp = QSizePolicy(QSizePolicy::Minimum,QSizePolicy::Minimum);
    sp.setHorizontalStretch(10);
    sp.setVerticalStretch(1);
    tabPlotTimeline->setSizePolicy(sp);
  }
  
  // input default values
  reset();
}

void SurfaceCoverageWidget::reset() {
  
  lineVestaMass->setText("2.71e20");
  // comboMassDistribution->setCurrentIndex((int)ComboMassDistribution::mdt_uniform);
  // comboShapeModel->setCurrentIndex((int)ComboShapeModel::smt_thomas);
  lineVestaPeriod->setText("5.342128799");
  lineVestaPoleEclipticLatitude->setText("59.2");
  lineVestaPoleEclipticLongitude->setText("319.5");
  
  lineOrbitEpoch->setText("2012-01-01 00:00:00");
  //
  lineOrbitRadius->setText("950.0");
  lineOrbitInclination->setText("85.0");
  lineOrbitPhase->setText("60.0");
  
  lineCameraParallelFOV->setText("5.4718");
  lineCameraOrthogonalFOV->setText("5.4718");
  lineCameraInterval->setText("300");
  lineCameraMinimumI->setText("0");
  lineCameraMaximumI->setText("75");
  lineCameraMaximumE->setText("55");
  
  lineRunDuration->setText("001:00:00:00");
  // lineRunSPICEFile->setText("test_orbit.bsp");
  // lineRunASCIIFile->setText("test_orbit.dat");

  progressRun->setEnabled(false);
  buttonAbort->setEnabled(false);
  
  setInputEnabled(true);
}

void SurfaceCoverageWidget::run() const {
  
  // ORSA_DEBUG("run called...");
  
  setInputEnabled(false);
  
  statusBar->showMessage(tr("Initializing Simulation..."));
  
  buttonRun->setEnabled(false);
  //
  progressRun->setEnabled(true);
  buttonAbort->setEnabled(true);
  
  buttonAbort->setFocus(Qt::OtherFocusReason);
  
  // reset before setting all the new values...
  mainThread->reset();

  {
    // set all variables

    mainThread->vestaMass                  = FromUnits(lineVestaMass->text().toDouble(),orsa::Unit::KG);
    // mainThread->vestaMassDistribution      = (ComboMassDistribution::MassDistributionType)comboMassDistribution->currentIndex();
    // mainThread->vestaShapeModel            = (ComboShapeModel::ShapeModelType)comboShapeModel->currentIndex();
    mainThread->vestaPeriod                = FromUnits(lineVestaPeriod->text().toDouble(),orsa::Unit::HOUR);
    mainThread->vestaPoleEclipticLatitude  = orsa::degToRad()*lineVestaPoleEclipticLatitude->text().toDouble();
    mainThread->vestaPoleEclipticLongitude = orsa::degToRad()*lineVestaPoleEclipticLongitude->text().toDouble();
    //
    {
      int y, m, d, H, M, S;
      sscanf(lineOrbitEpoch->text().toStdString().c_str(),
	     "%d-%d-%d %d:%d:%d",
	     &y,
	     &m,
	     &d,
	     &H,
	     &M,
	     &S);

      mainThread->orbitEpoch = orsaSolarSystem::gregorTime(y,
							   m,
							   d,
							   H,
							   M,
							   S,
							   0);
      /*
	ORSA_DEBUG("lineOrbitEpoch->text().toStdString().c_str(): [%s]",
	 lineOrbitEpoch->text().toStdString().c_str());
      */
    }
    //
    mainThread->orbitRadius      = orsa::FromUnits(lineOrbitRadius->text().toDouble(),orsa::Unit::KM);
    mainThread->orbitInclination = orsa::degToRad()*lineOrbitInclination->text().toDouble();
    mainThread->orbitPhase       = orsa::degToRad()*lineOrbitPhase->text().toDouble();
    //
    mainThread->cameraParallelFOV   = orsa::degToRad()*lineCameraParallelFOV->text().toDouble();
    mainThread->cameraOrthogonalFOV = orsa::degToRad()*lineCameraOrthogonalFOV->text().toDouble();
    //
    {
      int S;
      sscanf(lineCameraInterval->text().toStdString().c_str(),
	     "%d", &S);
      
      // ORSA_DEBUG("S: %i",S);
      
      mainThread->cameraInterval = orsa::Time(0,0,0,S,0);
    }
    //
    mainThread->cameraMinimumI = orsa::degToRad()*lineCameraMinimumI->text().toDouble();
    mainThread->cameraMaximumI = orsa::degToRad()*lineCameraMaximumI->text().toDouble();
    mainThread->cameraMaximumE = orsa::degToRad()*lineCameraMaximumE->text().toDouble();
    //
    {
      int d, H, M, S;
      sscanf(lineRunDuration->text().toStdString().c_str(),
	     "%d:%d:%d:%d",
	     &d,
	     &H,
	     &M,
	     &S);

      mainThread->runDuration = orsa::Time(d,H,M,S,0);

      /*
	 ORSA_DEBUG("lineRunDuration->text().toStdString().c_str(): [%s]",
	 lineRunDuration->text().toStdString().c_str());
	 ORSA_DEBUG("d: %i   H: %i   M: %i   S: %i",
	 d,H,M,S);
      */
    }
    //
    // mainThread->outputSPICEFile = lineRunSPICEFile->text().toStdString();
    // mainThread->outputASCIIFile = lineRunASCIIFile->text().toStdString();
    
  }
  
  // this returns instantly!
  mainThread->start();
}

void SurfaceCoverageWidget::abort() {

  const QMessageBox::StandardButton sb = QMessageBox::question(this,
							       tr("Abort Simulation?"),
							       tr("Do you want to abort this simulation?"),
							       QMessageBox::Yes | QMessageBox::No,
							       QMessageBox::No);

  if (sb == QMessageBox::Yes) {

    statusBar->showMessage(tr("Aborting Simulation..."));

    mainThread->abort();
  }

}

void SurfaceCoverageWidget::progressChanged(int value) const {
  if (value == 100) {
    statusBar->showMessage(tr("Writing Output Files..."));
  } else {
    statusBar->showMessage(tr("Simulation Running..."));
  }
}

void SurfaceCoverageWidget::runFinished() const {
  bool completed;
  if (progressRun->value() == 100) {
    statusBar->showMessage(tr("Simulation Completed"));
    completed = true;
  } else {
    statusBar->showMessage(tr("Simulation Aborted"));
    completed = false;
  }
  
  plot2D->setMainThread(mainThread);
  
  plotTimeline->setMainThread(mainThread);
  
  buttonRun->setEnabled(true);
  //
  progressRun->setEnabled(false);
  buttonAbort->setEnabled(false);
  
  setInputEnabled(true);
}

void SurfaceCoverageWidget::showAboutDialog() {
  
  const QString mainText =
    "This program is part of ORSA."
    "<br>"
    "<br>"
    "ORSA - Orbit Reconstruction, Simulation and Analysis"
    "<br>"
    "Copyright (C) 2002-2008 Pasquale Tricarico"
    "<br>"
    "<br>"
    "ORSA is free software and comes with ABSOLUTELY NO WARRANTY."
    "<br>"
    "See the GNU General Public License for more details."
    "<br>"
    "<br>"
    // "Website: <a href=\"http://orsa.sourceforge.net\">http://orsa.sourceforge.net</a>"
    "Website: <a href=\"http://orsa.sf.net\">http://orsa.sf.net</a>"
    "<br>";

  QMessageBox::about(this,
		     "About SurfaceCoverage",
		     mainText);
}

void SurfaceCoverageWidget::computePeriodsRatio(const QString &) {
  
  QString s;  
  
  s.sprintf("--.--- : 1");
  
  const orsa::Double vestaMass = FromUnits(lineVestaMass->text().toDouble(),orsa::Unit::KG);
  
  if (lineVestaPeriod->text().length() != 0) { 
    
    const orsa::Double vestaPeriod = FromUnits(lineVestaPeriod->text().toDouble(),orsa::Unit::HOUR);
    
    const orsa::Double orbitRadius = orsa::FromUnits(lineOrbitRadius->text().toDouble(),orsa::Unit::KM);
    
    orsa::Orbit orbit;
    // only 'mu' and 'a' needed
    orbit.mu = orsa::Unit::instance()->getG() * vestaMass;
    orbit.a  = orbitRadius;
    //
    if ((orbit.mu > orsa::zero()) && (orbit.a > orsa::zero())) {
      const orsa::Double orbitPeriod = orbit.period();
      s.sprintf("%.3f : 1",orsa::Double(orbitPeriod/vestaPeriod).get_d());
    }
    
  }
  
  lineOrbitPeriodRatio->setText(s);
}

void SurfaceCoverageWidget::setInputEnabled(const bool b) const {
  
  lineVestaMass->setEnabled(b);
  // comboMassDistribution->setEnabled(b);
  // comboShapeModel->setEnabled(b);
  lineVestaPeriod->setEnabled(b);
  lineVestaPoleEclipticLatitude->setEnabled(b);
  lineVestaPoleEclipticLongitude->setEnabled(b);

  lineOrbitEpoch->setEnabled(b);
  //
  lineOrbitRadius->setEnabled(b);
  // lineOrbitPeriodRatio->setEnabled(b);
  lineOrbitInclination->setEnabled(b);
  lineOrbitPhase->setEnabled(b);
  
  lineCameraParallelFOV->setEnabled(b);
  lineCameraOrthogonalFOV->setEnabled(b);
  lineCameraInterval->setEnabled(b);
  lineCameraMinimumI->setEnabled(b);  
  lineCameraMaximumI->setEnabled(b);  
  lineCameraMaximumE->setEnabled(b);
  
  lineRunDuration->setEnabled(b);
  // lineRunSPICEFile->setEnabled(b);
  // buttonBrowseSPICE->setEnabled(b);
  // lineRunASCIIFile->setEnabled(b);
  // buttonBrowseASCII->setEnabled(b);

  // buttonAbout->setEnabled(b);

}

void SurfaceCoverageWidget::closeEvent(QCloseEvent *) {
  QCoreApplication::quit();
}

// mainly useful with QLineEdit widgets using an InputMask
void SurfaceCoverageWidget::lineEdit_textChanged(const QString &) {
  
  bool pass = true;
  
  {
    // lineOrbitEpoch

    /*
       ORSA_DEBUG("RendezvousWithVestaWidget::lineOrbitEpoch_textChanged(): [%s]   length: %i",
       text.toStdString().c_str(),
       strlen(text.toStdString().c_str()));
    */

    QPalette pal = lineOrbitEpoch->palette();
    //
    if (strlen(lineOrbitEpoch->text().toStdString().c_str()) == 19) {
      pal.setColor(QPalette::Active,QPalette::Base,Qt::white);
      // buttonRun->setEnabled(true);
    } else {
      pal.setColor(QPalette::Active,QPalette::Base,Qt::yellow);
      pass = false;
      // buttonRun->setEnabled(false);
    }
    //
    lineOrbitEpoch->setPalette(pal);
  }

  {
    // lineRunDuration

    /*
       ORSA_DEBUG("RendezvousWithVestaWidget::lineRunDuration_textChanged(): [%s]   length: %i",
       text.toStdString().c_str(),
       strlen(text.toStdString().c_str()));
    */

    QPalette pal = lineRunDuration->palette();
    //
    if (strlen(lineRunDuration->text().toStdString().c_str()) == 12) {
      pal.setColor(QPalette::Active,QPalette::Base,Qt::white);
      // buttonRun->setEnabled(true);
    } else {
      pal.setColor(QPalette::Active,QPalette::Base,Qt::yellow);
      pass = false;
      // buttonRun->setEnabled(false);
    }
    //
    lineRunDuration->setPalette(pal);
  }

  if (pass) {
    buttonRun->setEnabled(true);
  } else {
    buttonRun->setEnabled(false);
  }
}
