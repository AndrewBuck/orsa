// this is to fix some random preprocessor problems
// #include <Qt>

#include <QComboBox>
#include <QDoubleValidator>
#include <QFileDialog>
#include <QFrame>
#include <QGridLayout>
#include <QGroupBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QPalette>
#include <QProgressBar>
#include <QPushButton>
#include <QStatusBar>
#include <QVariant>
#include <QVBoxLayout>

#include "vestaViz.h"

#include "RendezvousWithVesta.h"

#include "RendezvousWithVestaVersion.h"

#include "mainThread.h"

#include <orsaOSG/viz.h>
#include <orsaOSG/FindNamedNodeVisitor.h>
#include <orsaOSG/Track.h>
#include <orsaOSG/DepthPartitionNode.h>
#include <orsaOSG/DistanceAccumulator.h>

#include <orsa/debug.h>

#include <orsaSolarSystem/datetime.h>

#include <osgGA/KeySwitchMatrixManipulator>
#include <osgGA/NodeTrackerManipulator>

#include <osgUtil/Optimizer>

#include <osgViewer/Viewer>

#include "plot.h"

// explicitly load all required OSG plugins
// ONLY when the OSG library is statically linked
#include <osgDB/Registry>
USE_OSGPLUGIN(freetype);

using namespace std;
using namespace orsa;
using namespace orsaSolarSystem;

RendezvousWithVestaWidget::RendezvousWithVestaWidget() : QWidget(0) {

  QVBoxLayout * layout = new QVBoxLayout(this);
  //
  layout->setMargin(3);

  //

  {
    // physical properties group

    QGroupBox * group = new QGroupBox(tr("Vesta properties"),this);
    layout->addWidget(group);

    QGridLayout * grid = new QGridLayout(group);
    //
    grid->setMargin(3);

    lineVestaMass                  = new QLineEdit(group);
    comboMassDistribution          = new ComboMassDistribution(group);
    comboShapeModel                = new ComboShapeModel(group);
    lineVestaPeriod                = new QLineEdit(group);
    lineVestaPoleEclipticLatitude  = new QLineEdit(group);
    lineVestaPoleEclipticLongitude = new QLineEdit(group);

    int localRow = 0;

    grid->addWidget(new QLabel(tr("mass [kg]"),group),localRow,0);
    grid->addWidget(lineVestaMass                          ,localRow,1);

    ++localRow;

    grid->addWidget(new QLabel(tr("mass distribution"),group),localRow,0);
    grid->addWidget(comboMassDistribution                    ,localRow,1);

    ++localRow;

    grid->addWidget(new QLabel(tr("shape model"),group),localRow,0);
    grid->addWidget(comboShapeModel                    ,localRow,1);

    ++localRow;

    grid->addWidget(new QLabel(tr("rotation period [hour]"),group),localRow,0);
    grid->addWidget(lineVestaPeriod                               ,localRow,1);

    ++localRow;

    grid->addWidget(new QLabel(tr("pole ecliptic latitude [DEG]"),group),localRow,0);
    grid->addWidget(lineVestaPoleEclipticLatitude                       ,localRow,1);

    ++localRow;

    grid->addWidget(new QLabel(tr("pole ecliptic longitude [DEG]"),group),localRow,0);
    grid->addWidget(lineVestaPoleEclipticLongitude                       ,localRow,1);

    ++localRow;

  }

  {
    // orbital parameters group

    QGroupBox * group = new QGroupBox(tr("DAWN initial orbit"),this);
    layout->addWidget(group);

    QGridLayout * grid = new QGridLayout(group);
    //
    grid->setMargin(3);

    lineOrbitEpoch       = new QLineEdit(group);
    lineOrbitRadius      = new QLineEdit(group);
    lineOrbitInclination = new QLineEdit(group);
    lineOrbitPhase       = new QLineEdit(group);

    int localRow = 0;

    grid->addWidget(new QLabel(tr("epoch [yyyy-mm-dd HH:MM:SS] [UTC]"),group),localRow,0);
    grid->addWidget(lineOrbitEpoch                                           ,localRow,1);

    ++localRow;

    grid->addWidget(new QLabel(tr("radius [km]"),group),localRow,0);
    grid->addWidget(lineOrbitRadius                    ,localRow,1);

    ++localRow;

    grid->addWidget(new QLabel(tr("equatorial inclination [DEG]"),group),localRow,0);
    grid->addWidget(lineOrbitInclination                                ,localRow,1);

    ++localRow;

    grid->addWidget(new QLabel(tr("phase angle [DEG]"),group),localRow,0);
    grid->addWidget(lineOrbitPhase                           ,localRow,1);

    ++localRow;

  }

  {
    // run group

    QGroupBox * group = new QGroupBox(tr("Simulation"),this);
    layout->addWidget(group);

    QGridLayout * grid = new QGridLayout(group);
    //
    grid->setMargin(3);

    lineRunDuration          = new QLineEdit(group);
    lineRunSPICEFile         = new QLineEdit(group);
    buttonBrowseSPICE        = new QPushButton("browse",group);
    lineRunASCIIFile         = new QLineEdit(group);
    buttonBrowseASCII        = new QPushButton("browse",group);
    buttonRun                = new QPushButton("run simulation",group);
    progressRun              = new QProgressBar(group);
    buttonAbort              = new QPushButton("abort",group);

    int localRow = 0;

    grid->addWidget(new QLabel(tr("duration [ddd:HH:MM:SS]"),group),localRow,0);
    grid->addWidget(lineRunDuration                                ,localRow,1);

    ++localRow;

    grid->addWidget(new QLabel(tr("output file [SPICE kernel]"),group),localRow,0);
    grid->addWidget(lineRunSPICEFile                                  ,localRow,1);
    grid->addWidget(buttonBrowseSPICE                                 ,localRow,2);

    ++localRow;

    grid->addWidget(new QLabel(tr("output file [ASCII]"),group),localRow,0);
    grid->addWidget(lineRunASCIIFile                           ,localRow,1);
    grid->addWidget(buttonBrowseASCII                          ,localRow,2);

    ++localRow;

    grid->addWidget(buttonRun  ,localRow,0);
    grid->addWidget(progressRun,localRow,1);
    grid->addWidget(buttonAbort,localRow,2);

    ++localRow;

  }

  {
    // about group

    QGroupBox * group = new QGroupBox(tr("About"),this);
    layout->addWidget(group);

    QGridLayout * grid = new QGridLayout(group);
    //
    grid->setMargin(3);

    buttonAbout = new QPushButton("more",group);

    int localRow = 0;

    {
      // QLabel * aboutLabel = new QLabel(tr("Copyright (c) 2007 Pasquale Tricarico <tricaric@psi.edu>"),group);
      QLabel * aboutLabel = new QLabel(tr("Copyright (c) 2009 Pasquale Tricarico <a href=\"mailto:tricaric@psi.edu\">tricaric@psi.edu</a>"),group);
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
    layout->addWidget(statusBar);
  }

  {
    // validation

    QDoubleValidator * positiveValidator = new QDoubleValidator(this);
    positiveValidator->setBottom(0.0);

    QDoubleValidator * doubleValidator = new QDoubleValidator(this);

    lineVestaMass->setValidator(positiveValidator);
    lineVestaPeriod->setValidator(positiveValidator);
    lineVestaPoleEclipticLatitude->setValidator(doubleValidator);
    lineVestaPoleEclipticLongitude->setValidator(doubleValidator);

    lineOrbitEpoch->setInputMask("9999-99-99 99:99:99;_");
    //
    lineOrbitRadius->setValidator(positiveValidator);
    lineOrbitInclination->setValidator(positiveValidator);
    lineOrbitPhase->setValidator(doubleValidator);

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

    connect(buttonBrowseSPICE,
	    SIGNAL(clicked()),
	    this,
	    SLOT(browseSPICE()));

    connect(buttonBrowseASCII,
	    SIGNAL(clicked()),
	    this,
	    SLOT(browseASCII()));

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

    /*
       connect(buttonAbort,
       SIGNAL(clicked()),
       mainThread,
       SLOT(abort()));
    */

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

  // input default values
  reset();

  // size policy
  setMaximumSize(width(),minimumSizeHint().height());

  // slightly bigger size
  resize(minimumSizeHint().width()+32,height());
}


void RendezvousWithVestaWidget::reset() {

  lineVestaMass->setText("2.685041251213e20");
  comboMassDistribution->setCurrentIndex((int)ComboMassDistribution::mdt_uniform);
  comboShapeModel->setCurrentIndex((int)ComboShapeModel::smt_thomas);
  // lineVestaPeriod->setText("5.3424");
  lineVestaPeriod->setText("5.342128799");
  // lineVestaPoleEclipticLatitude->setText("62.0");
  // lineVestaPoleEclipticLongitude->setText("343.0");
  lineVestaPoleEclipticLatitude->setText("59.2");
  lineVestaPoleEclipticLongitude->setText("319.5");

  lineOrbitEpoch->setText("2012-01-01 00:00:00");
  //
  lineOrbitRadius->setText("950.0");
  lineOrbitInclination->setText("85.0");
  lineOrbitPhase->setText("60.0");
  
  lineRunDuration->setText("001:00:00:00");
  lineRunSPICEFile->setText("test_orbit.bsp");
  lineRunASCIIFile->setText("test_orbit.dat");
  
  progressRun->setEnabled(false);
  buttonAbort->setEnabled(false);

  setInputEnabled(true);
}

void RendezvousWithVestaWidget::run() const {

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

    mainThread->vestaMass                  = FromUnits(lineVestaMass->text().toDouble(),Unit::KG);
    mainThread->vestaMassDistribution      = (ComboMassDistribution::MassDistributionType)comboMassDistribution->currentIndex();
    mainThread->vestaShapeModel            = (ComboShapeModel::ShapeModelType)comboShapeModel->currentIndex();
    mainThread->vestaPeriod                = FromUnits(lineVestaPeriod->text().toDouble(),Unit::HOUR);
    mainThread->vestaPoleEclipticLatitude  = degToRad()*lineVestaPoleEclipticLatitude->text().toDouble();
    mainThread->vestaPoleEclipticLongitude = degToRad()*lineVestaPoleEclipticLongitude->text().toDouble();
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

      mainThread->orbitEpoch = gregorTime(y,
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
    mainThread->orbitRadius      = FromUnits(lineOrbitRadius->text().toDouble(),Unit::KM);
    mainThread->orbitInclination = degToRad()*lineOrbitInclination->text().toDouble();
    mainThread->orbitPhase       = degToRad()*lineOrbitPhase->text().toDouble();
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
    mainThread->outputSPICEFile = lineRunSPICEFile->text().toStdString();
    mainThread->outputASCIIFile = lineRunASCIIFile->text().toStdString();

  }

  // this returns instantly!
  mainThread->start();
}

void RendezvousWithVestaWidget::abort() {

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

void RendezvousWithVestaWidget::progressChanged(int value) const {
  if (value == 100) {
    statusBar->showMessage(tr("Writing Output Files..."));
  } else {
    statusBar->showMessage(tr("Simulation Running..."));
  }
}

void RendezvousWithVestaWidget::runFinished() const {

  bool completed;
  if (progressRun->value() == 100) {
    statusBar->showMessage(tr("Simulation Completed"));
    completed = true;
  } else {
    statusBar->showMessage(tr("Simulation Aborted"));
    completed = false;
  }

  if (1) if (completed) {

    // viz
    osg::ref_ptr<orsaOSG::Viz> viz = new orsaOSG::Viz(mainThread->bg.get(),
						      600.0);

    ViewerQT * viewerWindow = new ViewerQT;
    //
    // viewerWindow->getSceneView()->setComputeNearFarMode( osgUtil::CullVisitor::DO_NOT_COMPUTE_NEAR_FAR );
    //
    // viewerWindow->getCamera()->setCullingMode(osg::CullSettings::VIEW_FRUSTUM_SIDES_CULLING);
    // viewerWindow->getCamera()->setCullingMode(osg::CullSettings::VIEW_FRUSTUM_CULLING);
    //
    //
    // IMPORTANT!!
    // this ration shouls be close to the ratio between
    // the size of the smallest object in the scene and
    // the size of the biggest object (orbit?) in the scene.
    // viewerWindow->getSceneView()->setNearFarRatio(1.0e-6);
    //
    // viewerWindow->getSceneView()->setNearFarRatio(1.0e-8);
    //
    // viewerWindow->getCamera()->setNearFarRatio(1.0e-8);
    // viewerWindow->getCamera()->setNearFarRatio(1.0e-12);
    
    // enable multisampling
    osg::DisplaySettings::instance()->setNumMultiSamples(4);
    
    osg::Group * rootNode = viz->createRoot();
    //
    // very important!
    viz->_at->setCentralBody(mainThread->bg->getBody("VESTA"));
    // viz->_at->setCentralBody(mainThread->bg->getBody("DAWN"));

    // run optimization over the scene graph
    /*
       {
       osgUtil::Optimizer optimzer;
       ORSA_DEBUG("calling Optimizer...");
       optimzer.optimize(rootNode);
       ORSA_DEBUG("Optimization finished...");
       }
    */

    {
      // track...
      const orsa::Time dt(0,0,10,0,0);
      const bool groundTrack = true;
      viz->getBodyAttitudeTransform(mainThread->bg->getBody("VESTA"))->addChild(new orsaOSG::Track(mainThread->bg.get(),
												   mainThread->bg->getBody("DAWN"),
												   mainThread->bg->getBody("VESTA"),
												   dt,
												   viz->_at.get(),
												   groundTrack));
    }

    /*
       {
       // track...
       const orsa::Time dt(0,0,10,0,0);
       const bool groundTrack = false;
       viz->getBodyPositionTransform(mainThread->bg->getBody("VESTA"))->addChild(new orsaOSG::Track(mainThread->bg.get(),
       mainThread->bg->getBody("DAWN"),
       mainThread->bg->getBody("VESTA"),
       dt,
       viz->_at.get(),
       groundTrack));
       }
    */

    /*
       {
       // track...
       const orsa::Time dt(0,0,10,0,0);
       const bool groundTrack = false;
       rootNode->addChild(new orsaOSG::Track(mainThread->bg.get(),
       mainThread->bg->getBody("DAWN"),
       0,
       dt,
       viz->_at.get(),
       groundTrack));
       }
    */
    
    // viewerWindow->setSceneData(rootNode);
    //
    osg::ref_ptr<DepthPartitionNode> dpn = new DepthPartitionNode;
    dpn->addChild(rootNode);
    dpn->setActive(true);
    //
    viewerWindow->setSceneData(dpn.get());
    //
    // depth partion node only supports single window/single threaded at present.
    viewerWindow->setThreadingModel(osgViewer::Viewer::SingleThreaded);
    
    /*
      {
       double ortho_d = 2*FromUnits(lineOrbitRadius->text().toDouble(),Unit::KM).get_d();
       
       viewerWindow->getCamera()->setProjectionMatrixAsOrtho(ortho_d,
       ortho_d,
       ortho_d,
       ortho_d,
       ortho_d,
       3*ortho_d);

       const osg::Vec3 eye(-2*ortho_d,0,0);
       const osg::Vec3 center(0,0,0);
       const osg::Vec3 up(0,0,1);
       //
       viewerWindow->getCamera()->setViewMatrixAsLookAt(eye,
       center,
       up);
       }
    */

    /*
       osg::ref_ptr<osgGA::KeySwitchMatrixManipulator> keyswitchManipulator =
       new osgGA::KeySwitchMatrixManipulator;
       //
       viewerWindow->setCameraManipulator( keyswitchManipulator.get() );
    */

    //viewerWindow->setCameraManipulator(new osgGA::TrackballManipulator);
    //
    if (1) {
      orsaOSG::FindNamedNodeVisitor fnnv("VESTA");
      // orsaOSG::FindNamedNodeVisitor fnnv("DAWN");
      // orsaOSG::FindNamedNodeVisitor fnnv("SUN");
      //
      rootNode->accept(fnnv);

      if (!fnnv._foundNodes.empty()) {

	// set up the node tracker.
	// osgGA::NodeTrackerManipulator * tm = new osgGA::NodeTrackerManipulator;
	VestaNodeTrackerManipulator * tm = new VestaNodeTrackerManipulator;

	osgGA::NodeTrackerManipulator::TrackerMode   trackerMode =
	  osgGA::NodeTrackerManipulator::NODE_CENTER; // NODE_CENTER_AND_ROTATION

	osgGA::NodeTrackerManipulator::RotationMode rotationMode =
	  osgGA::NodeTrackerManipulator::TRACKBALL;   // TRACKBALL ELEVATION_AZIM

	tm->setTrackerMode(  trackerMode);
	tm->setRotationMode(rotationMode);
       	//
	tm->setTrackNode(fnnv._foundNodes.front().get());

	/*
	   tm->setHomePosition(osg::Vec3d(0.0,FromUnits(-1500.0,Unit::KM).get_d(),0.0),
	   osg::Vec3d(0.0,0.0,0.0),
	   osg::Vec3d(0.0,0.0,1.0));
	*/

	viewerWindow->setCameraManipulator(tm);

      }

    } else {

      viewerWindow->setCameraManipulator(new VestaNodeTrackerManipulator);

    }

    // test
    {
      /*
	 osg::Camera * camera = viewerWindow->getSceneView()->getCamera();
	 camera->setProjectionMatrixAsPerspective(45.0,
	 1.0,
	 FromUnits(1.0e2,Unit::KM).get_d(),
	 FromUnits(1.0e6,Unit::KM).get_d());
	 viewerWindow->getSceneView()->setCamera(camera);
      */
      //
      // viewerWindow->getSceneView()->setComputeNearFarMode( osgUtil::CullVisitor::DO_NOT_COMPUTE_NEAR_FAR );
      //
      /*
	 viewerWindow->getSceneView()->setProjectionMatrixAsPerspective(60.0,
	 1.4,
	 FromUnits(1.0e1,Unit::KM).get_d(),
	 FromUnits(1.0e5,Unit::KM).get_d());
      */
      //
      /*
	 viewerWindow->getSceneView()->setViewMatrixAsLookAt(osg::Vec3d(0.0,-FromUnits(2000.0,Unit::KM).get_d(),0.0),
	 osg::Vec3d(0.0,0.0,0.0),
	 osg::Vec3d(0.0,0.0,1.0));
      */
      //
      /*
	 viewerWindow->getSceneView()->setProjectionMatrixAsFrustum(FromUnits(1.0e3,Unit::KM).get_d(),
	 FromUnits(1.0e3,Unit::KM).get_d(),
	 FromUnits(1.0e3,Unit::KM).get_d(),
	 FromUnits(1.0e3,Unit::KM).get_d(),
	 FromUnits(1.0e1,Unit::KM).get_d(),
	 FromUnits(1.0e6,Unit::KM).get_d());
      */

      if (0) {
	double fovy;
	double aspectRatio;
	double zNear;
	double zFar;
	//
	viewerWindow->getCamera()->getProjectionMatrixAsPerspective(fovy, aspectRatio,
								    zNear, zFar);
	//
	ORSA_DEBUG("fovy.......: %f",fovy);
	ORSA_DEBUG("aspectRatio: %f",aspectRatio);
	ORSA_DEBUG("zNear......: %f",zNear);
	ORSA_DEBUG("zFar.......: %f",zFar);
      }
    }
    
    viewerWindow->show();
    
    if (1) {
      // plot

      VestaPlot * plot = new VestaPlot(mainThread->bg.get(),
				       viz->_at.get());

      plot->show();
    }

  }

  buttonRun->setEnabled(true);
  //
  progressRun->setEnabled(false);
  buttonAbort->setEnabled(false);

  setInputEnabled(true);
}

void RendezvousWithVestaWidget::browseSPICE() {

  QFileDialog dialog(this,"Select SPICE output file",".");
  dialog.setFileMode(QFileDialog::AnyFile);
  dialog.setAcceptMode(QFileDialog::AcceptSave);
  dialog.setConfirmOverwrite(true);

  QStringList filters;
  filters << "SPICE kernel (*.bsp)"
	  << "Any files (*)";
  dialog.setFilters(filters);

  if (dialog.exec()) {
    QStringList fileNames = dialog.selectedFiles();
    if (fileNames[0].length() != 0) {
      lineRunSPICEFile->setText(fileNames[0]);
    }
  }
}

void RendezvousWithVestaWidget::browseASCII() {

  QFileDialog dialog(this,"Select ASCII output file",".");
  dialog.setFileMode(QFileDialog::AnyFile);
  dialog.setAcceptMode(QFileDialog::AcceptSave);
  dialog.setConfirmOverwrite(true);

  QStringList filters;
  filters << "ASCII file (*.dat *.txt)"
	  << "Any files (*)";
  dialog.setFilters(filters);

  if (dialog.exec()) {
    QStringList fileNames = dialog.selectedFiles();
    if (fileNames[0].length() != 0) {
      lineRunASCIIFile->setText(fileNames[0]);
    }
  }
}

void RendezvousWithVestaWidget::showAboutDialog() {

  const QString mainText =
    "This program is part of ORSA."
    "<br>"
    "<br>"
    "ORSA - Orbit Reconstruction, Simulation and Analysis"
    "<br>"
    "Copyright (C) 2002-2009 Pasquale Tricarico"
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
		     "About RendezvousWithVesta",
		     mainText);
}

void RendezvousWithVestaWidget::setInputEnabled(const bool b) const {

  lineVestaMass->setEnabled(b);
  comboMassDistribution->setEnabled(b);
  comboShapeModel->setEnabled(b);
  lineVestaPeriod->setEnabled(b);
  lineVestaPoleEclipticLatitude->setEnabled(b);
  lineVestaPoleEclipticLongitude->setEnabled(b);

  lineOrbitEpoch->setEnabled(b);
  //
  lineOrbitRadius->setEnabled(b);
  lineOrbitInclination->setEnabled(b);
  lineOrbitPhase->setEnabled(b);

  lineRunDuration->setEnabled(b);
  lineRunSPICEFile->setEnabled(b);
  buttonBrowseSPICE->setEnabled(b);
  lineRunASCIIFile->setEnabled(b);
  buttonBrowseASCII->setEnabled(b);

  // buttonAbout->setEnabled(b);

}

void RendezvousWithVestaWidget::closeEvent(QCloseEvent *) {
  QCoreApplication::quit();
}

// mainly useful with QLineEdit widgets using an InputMask
void RendezvousWithVestaWidget::lineEdit_textChanged(const QString &) {

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
