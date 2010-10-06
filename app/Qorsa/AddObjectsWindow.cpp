#include <iostream>
using namespace std;

#include <QtCore/QtAlgorithms>

#include "MainWindow.h"
#include "AddObjectsWindow.h"

AddObjectsWindow::AddObjectsWindow(MainWindow *nSpawningWindow, QWidget *nParent)
	: QWidget(nParent)
{
	spawningWindow = nSpawningWindow;

	cout << "Adding objects...\n";
	QString tempString = "Add Objects to the Universe: " + spawningWindow->getWorkspaceName();
	setWindowTitle(tempString);

	resize(600, 800);

	// Create the three main group boxes for the add objects window and arrange them in a QVBoxLayout.
	addObjectsVBoxLayout = new QVBoxLayout();
	objectSourceGroupBox = new QGroupBox("Object Source");
	objectSelectionGroupBox = new QGroupBox("Objects Available for Selection");
	selectedObjectsGroupBox = new QGroupBox("Selected Objects to be Added to the Universe");
	addObjectsVBoxLayout->addWidget(objectSourceGroupBox);
	addObjectsVBoxLayout->addWidget(objectSelectionGroupBox);
	addObjectsVBoxLayout->addWidget(selectedObjectsGroupBox);

	// Create the widgets for the sub-display for each of the object sources.
	// The "Local File" tab
	localFileSourceGridLayout = new QGridLayout();
	localFileSourceWidget = new QWidget();
	localFileSourceFilenameLabel = new QLabel("Filename");
	localFileSourceLineEdit = new QLineEdit();
	localFileSourceBrowseButton = new QPushButton("Browse");
	localFileSourceGridLayout->addWidget(localFileSourceFilenameLabel, 0, 0);
	localFileSourceGridLayout->addWidget(localFileSourceLineEdit, 0, 1);
	localFileSourceGridLayout->addWidget(localFileSourceBrowseButton, 0, 2);
	localFileSourceWidget->setLayout(localFileSourceGridLayout);

	// The "MYSQL Query" tab
	mysqlDatabaseSourceGridLayout = new QGridLayout();
	mysqlDatabaseSourceWidget = new QWidget();
	mysqlDatabaseSourceHostnameLabel = new QLabel("Hostname");
	mysqlDatabaseSourcePortLabel = new QLabel("Port");
	mysqlDatabaseSourceUsernameLabel = new QLabel("Username");
	mysqlDatabaseSourcePasswordLabel = new QLabel("Password");
	mysqlDatabaseSourceDatabaseLabel = new QLabel("Database");
	mysqlDatabaseSourceQueryLabel = new QLabel("Query:");
	mysqlDatabaseSourceHostnameLineEdit = new QLineEdit("localhost");
	mysqlDatabaseSourcePortLineEdit = new QLineEdit("3306");
	mysqlDatabaseSourceUsernameLineEdit = new QLineEdit();
	mysqlDatabaseSourcePasswordLineEdit = new QLineEdit();
	mysqlDatabaseSourcePasswordLineEdit->setEchoMode(QLineEdit::Password);
	mysqlDatabaseSourceDatabaseLineEdit = new QLineEdit();
	mysqlDatabaseSourceQueryTextEdit = new QTextEdit();
	mysqlDatabaseSourceQueryTextEdit->setAcceptRichText(false);
	mysqlDatabaseSourceQuerySpacer = new QSpacerItem(1, 15);

	mysqlDatabaseSourceGridLayout->addWidget(mysqlDatabaseSourceHostnameLabel, 0, 0);
	mysqlDatabaseSourceGridLayout->addWidget(mysqlDatabaseSourcePortLabel, 0, 4);
	mysqlDatabaseSourceGridLayout->addWidget(mysqlDatabaseSourceUsernameLabel, 1, 0);
	mysqlDatabaseSourceGridLayout->addWidget(mysqlDatabaseSourcePasswordLabel, 1, 2);
	mysqlDatabaseSourceGridLayout->addWidget(mysqlDatabaseSourceDatabaseLabel, 1, 4);
	mysqlDatabaseSourceGridLayout->addItem(mysqlDatabaseSourceQuerySpacer, 2, 0);
	mysqlDatabaseSourceGridLayout->addWidget(mysqlDatabaseSourceQueryLabel, 3, 0);

	mysqlDatabaseSourceGridLayout->addWidget(mysqlDatabaseSourceHostnameLineEdit, 0, 1, 1, 3);
	mysqlDatabaseSourceGridLayout->addWidget(mysqlDatabaseSourcePortLineEdit, 0, 5);
	mysqlDatabaseSourceGridLayout->addWidget(mysqlDatabaseSourceUsernameLineEdit, 1, 1);
	mysqlDatabaseSourceGridLayout->addWidget(mysqlDatabaseSourcePasswordLineEdit, 1, 3);
	mysqlDatabaseSourceGridLayout->addWidget(mysqlDatabaseSourceDatabaseLineEdit, 1, 5);
	mysqlDatabaseSourceGridLayout->addWidget(mysqlDatabaseSourceQueryTextEdit, 4, 0, 1, -1);

	mysqlDatabaseSourceWidget->setLayout(mysqlDatabaseSourceGridLayout);

	// The "Custom Object" tab
	customObjectSourceGridLayout = new QGridLayout();
	customObjectSourceWidget = new QWidget();
	customObjectSourceVectorGroupBox = new QGroupBox("Vector State");
	customObjectSourceKeplerianGroupBox = new QGroupBox("Keplerian Elements");
	customObjectSourceVectorGroupBox->setCheckable(true);
	customObjectSourceKeplerianGroupBox->setCheckable(true);
	customObjectSourceVectorGroupBoxGridLayout = new QGridLayout();
	customObjectSourceKeplerianGroupBoxGridLayout = new QGridLayout();
	customObjectSourceVectorGroupBox->setLayout(customObjectSourceVectorGroupBoxGridLayout);
	customObjectSourceKeplerianGroupBox->setLayout(customObjectSourceKeplerianGroupBoxGridLayout);
	customObjectSourceVectorGroupBox->setChecked(true);
	customObjectSourceKeplerianGroupBox->setChecked(false);
	QObject::connect(customObjectSourceVectorGroupBox, SIGNAL(toggled(bool)), this, SLOT(customObjectSourceVectorGroupBoxToggled(bool)));
	QObject::connect(customObjectSourceKeplerianGroupBox, SIGNAL(toggled(bool)), this, SLOT(customObjectSourceKeplerianGroupBoxToggled(bool)));
	customObjectSourceGridLayout->addWidget(customObjectSourceVectorGroupBox, 2, 0, 1, 2);
	customObjectSourceGridLayout->addWidget(customObjectSourceKeplerianGroupBox, 2, 3, 1, 2);

	customObjectSourceNameLabel = new QLabel("Name");
	customObjectSourceMassLabel = new QLabel("Mass");

	customObjectSourceGridLayout->addWidget(customObjectSourceNameLabel, 0, 0);
	customObjectSourceGridLayout->addWidget(customObjectSourceMassLabel, 0, 3);

	customObjectSourceXLabel = new QLabel("X");
	customObjectSourceYLabel = new QLabel("Y");
	customObjectSourceZLabel = new QLabel("Z");
	customObjectSourceVXLabel = new QLabel("VX");
	customObjectSourceVYLabel = new QLabel("VY");
	customObjectSourceVZLabel = new QLabel("VZ");

	customObjectSourceVectorGroupBoxGridLayout->addWidget(customObjectSourceXLabel, 0, 0);
	customObjectSourceVectorGroupBoxGridLayout->addWidget(customObjectSourceYLabel, 1, 0);
	customObjectSourceVectorGroupBoxGridLayout->addWidget(customObjectSourceZLabel, 2, 0);
	customObjectSourceVectorGroupBoxGridLayout->addWidget(customObjectSourceVXLabel, 0, 3);
	customObjectSourceVectorGroupBoxGridLayout->addWidget(customObjectSourceVYLabel, 1, 3);
	customObjectSourceVectorGroupBoxGridLayout->addWidget(customObjectSourceVZLabel, 2, 3);

	customObjectSourceEccentricityLabel = new QLabel("Ecc.");
	customObjectSourceSMALabel = new QLabel("Semi Maj.\nAxis");
	customObjectSourceInclinationLabel = new QLabel("Incl.");
	customObjectSourceNodeLabel = new QLabel("Lon. Asc.\nNode");
	customObjectSourcePeriapsisLabel = new QLabel("Arg.\nPeriapsis");
	customObjectSourceAnomalyLabel = new QLabel("Mean\nAnomaly");

	customObjectSourceKeplerianGroupBoxGridLayout->addWidget(customObjectSourceSMALabel, 0, 0);
	customObjectSourceKeplerianGroupBoxGridLayout->addWidget(customObjectSourceEccentricityLabel, 0, 3);
	customObjectSourceKeplerianGroupBoxGridLayout->addWidget(customObjectSourceInclinationLabel, 1, 0);
	customObjectSourceKeplerianGroupBoxGridLayout->addWidget(customObjectSourceNodeLabel, 1, 3);
	customObjectSourceKeplerianGroupBoxGridLayout->addWidget(customObjectSourcePeriapsisLabel, 2, 0);
	customObjectSourceKeplerianGroupBoxGridLayout->addWidget(customObjectSourceAnomalyLabel, 2, 3);

	customObjectSourceNameLineEdit = new QLineEdit();
	customObjectSourceMassLineEdit = new QLineEdit();

	customObjectSourceGridLayout->addWidget(customObjectSourceNameLineEdit, 0, 1);
	customObjectSourceGridLayout->addWidget(customObjectSourceMassLineEdit, 0, 4);

	customObjectSourceXLineEdit = new QLineEdit();
	customObjectSourceYLineEdit = new QLineEdit();
	customObjectSourceZLineEdit = new QLineEdit();
	customObjectSourceVXLineEdit = new QLineEdit();
	customObjectSourceVYLineEdit = new QLineEdit();
	customObjectSourceVZLineEdit = new QLineEdit();

	customObjectSourceVectorGroupBoxGridLayout->addWidget(customObjectSourceXLineEdit, 0, 1);
	customObjectSourceVectorGroupBoxGridLayout->addWidget(customObjectSourceYLineEdit, 1, 1);
	customObjectSourceVectorGroupBoxGridLayout->addWidget(customObjectSourceZLineEdit, 2, 1);
	customObjectSourceVectorGroupBoxGridLayout->addWidget(customObjectSourceVXLineEdit, 0, 4);
	customObjectSourceVectorGroupBoxGridLayout->addWidget(customObjectSourceVYLineEdit, 1, 4);
	customObjectSourceVectorGroupBoxGridLayout->addWidget(customObjectSourceVZLineEdit, 2, 4);

	customObjectSourceEccentricityLineEdit = new QLineEdit();
	customObjectSourceSMALineEdit = new QLineEdit();
	customObjectSourceInclinationLineEdit = new QLineEdit();
	customObjectSourceNodeLineEdit = new QLineEdit();
	customObjectSourcePeriapsisLineEdit = new QLineEdit();
	customObjectSourceAnomalyLineEdit = new QLineEdit();

	customObjectSourceKeplerianGroupBoxGridLayout->addWidget(customObjectSourceSMALineEdit, 0, 1);
	customObjectSourceKeplerianGroupBoxGridLayout->addWidget(customObjectSourceEccentricityLineEdit, 0, 4);
	customObjectSourceKeplerianGroupBoxGridLayout->addWidget(customObjectSourceInclinationLineEdit, 1, 1);
	customObjectSourceKeplerianGroupBoxGridLayout->addWidget(customObjectSourceNodeLineEdit, 1, 4);
	customObjectSourceKeplerianGroupBoxGridLayout->addWidget(customObjectSourcePeriapsisLineEdit, 2, 1);
	customObjectSourceKeplerianGroupBoxGridLayout->addWidget(customObjectSourceAnomalyLineEdit, 2, 4);

	customObjectSourceInsertButton = new QPushButton("Add Object for Selection");
	customObjectSourceGridLayout->addWidget(customObjectSourceInsertButton, 3, 2, 1, 2);
	QObject::connect(customObjectSourceInsertButton, SIGNAL(released()), this, SLOT(customObjectSourceInsertButtonPressed()));

	customObjectSourceWidget->setLayout(customObjectSourceGridLayout);

	// Set up the tab widget at the top of the window to allow the selection of the input source.
	objectSourceDetailsTabWidget = new QTabWidget();
	objectSourceDetailsTabWidgetLayout = new QVBoxLayout();
	objectSourceDetailsTabWidget->addTab(localFileSourceWidget, "Local File");
	objectSourceDetailsTabWidget->addTab(mysqlDatabaseSourceWidget, "MYSQL Query");
	objectSourceDetailsTabWidget->addTab(customObjectSourceWidget, "Custom Object");
	objectSourceDetailsTabWidgetLayout->addWidget(objectSourceDetailsTabWidget);
	objectSourceGroupBox->setLayout(objectSourceDetailsTabWidgetLayout);

	// Fill in the object selection group box.
	objectSelectionVBoxLayout = new QVBoxLayout();
	objectSelectionTableView = new QTableView();
	objectSelectionTableModel = new BodyTableModel();
	objectSelectionTableView->setModel(objectSelectionTableModel);
	objectSelectionTableView->setShowGrid(true);
	objectSelectionTableView->setSelectionBehavior(QAbstractItemView::SelectRows);
	objectSelectionSelectButton = new QPushButton("Add Selected Objects to Resultset");
	objectSelectionVBoxLayout->addWidget(objectSelectionTableView);
	objectSelectionVBoxLayout->addWidget(objectSelectionSelectButton);
	objectSelectionGroupBox->setLayout(objectSelectionVBoxLayout);

	QObject::connect(objectSelectionSelectButton, SIGNAL(released()), this, SLOT(objectSelectionSelectButtonPressed()));

	// Fill in the selected objects group box.
	selectedObjectsVBoxLayout = new QVBoxLayout();
	selectedObjectsTableView = new QTableView();
	selectedObjectsHBoxLayout = new QHBoxLayout();
	selectedObjectsOkButton = new QPushButton("Ok");
	selectedObjectsCancelButton = new QPushButton("Cancel");
	selectedObjectsTableModel = new BodyTableModel();
	selectedObjectsTableView->setModel(selectedObjectsTableModel);
	selectedObjectsTableView->setShowGrid(true);
	selectedObjectsTableView->setSelectionBehavior(QAbstractItemView::SelectRows);
	selectedObjectsHBoxLayout->addWidget(selectedObjectsOkButton);
	selectedObjectsHBoxLayout->addWidget(selectedObjectsCancelButton);
	selectedObjectsVBoxLayout->addWidget(selectedObjectsTableView);
	selectedObjectsVBoxLayout->addLayout(selectedObjectsHBoxLayout);
	selectedObjectsGroupBox->setLayout(selectedObjectsVBoxLayout);

	// Connect the Ok and Cancel buttons to their respective handler functions.
	QObject::connect(selectedObjectsOkButton, SIGNAL(released()), this, SLOT(okButtonPressed()));
	QObject::connect(selectedObjectsCancelButton, SIGNAL(released()), this, SLOT(cancelButtonPressed()));

	setLayout(addObjectsVBoxLayout);
}

void AddObjectsWindow::closeEvent(QCloseEvent *event)
{
	emit addObjectsWindowClosed();
}

void AddObjectsWindow::okButtonPressed()
{
	// Delete the bodies we created that are still in the object selection list since they will not be moved into the main simulation.
	objectSelectionTableModel->clearAndEraseAllBodies();

	// Return the bodies in the selection window back to the main window and add them to the body group there.
	for(int i = 0; i < selectedObjectsTableModel->rowCount(); i++)
	{
		// Add the current body to the objects in the main window.
		spawningWindow->addBody(selectedObjectsTableModel->getBody(i));
	}

	// Now that all the bodies have been added to the main window we need to remove
	// the references from the table model (but do not delete them, just unlink them).
	selectedObjectsTableModel->clearAllBodies();

	//TODO:  Write the destructor to clean up all the qt widgets we created.
	close();
}

void AddObjectsWindow::cancelButtonPressed()
{
	close();
}

void AddObjectsWindow::customObjectSourceVectorGroupBoxToggled(bool newState)
{
	customObjectSourceKeplerianGroupBox->setChecked(!newState);
}

void AddObjectsWindow::customObjectSourceKeplerianGroupBoxToggled(bool newState)
{
	customObjectSourceVectorGroupBox->setChecked(!newState);
}

void AddObjectsWindow::customObjectSourceInsertButtonPressed()
{
	orsa::Orbit tempOrbit;
	orsa::Body *tempBody;
	QString tempString;
	orsa::Vector tempPosition, tempVelocity;

	QString tempName = customObjectSourceNameLineEdit->text();
	if(tempName.length() == 0)  return;

	tempString = customObjectSourceMassLineEdit->text();
	if(tempString.length() == 0)  return;
	double tempMass = tempString.toDouble();

	if(customObjectSourceVectorGroupBox->isChecked())
	{
		if(customObjectSourceXLineEdit->text().length() == 0 || \
			customObjectSourceYLineEdit->text().length() == 0 || \
			customObjectSourceZLineEdit->text().length() == 0 || \
			customObjectSourceVXLineEdit->text().length() == 0 || \
			customObjectSourceVYLineEdit->text().length() == 0 || \
			customObjectSourceVZLineEdit->text().length() == 0)
		{
			return;
		}

		tempPosition.set(customObjectSourceXLineEdit->text().toDouble(),
			customObjectSourceYLineEdit->text().toDouble(),
			customObjectSourceZLineEdit->text().toDouble());

		tempVelocity.set(customObjectSourceVXLineEdit->text().toDouble(),
			customObjectSourceVYLineEdit->text().toDouble(),
			customObjectSourceVZLineEdit->text().toDouble());

		cout << "Object with vector elements added:\n";
	}
	else if(customObjectSourceKeplerianGroupBox->isChecked())
	{
		tempString = customObjectSourceEccentricityLineEdit->text();
		if(tempString.length() == 0)  return;
		tempOrbit.e = tempString.toDouble();

		tempString = customObjectSourceSMALineEdit->text();
		if(tempString.length() == 0)  return;
		tempOrbit.a = tempString.toDouble();

		tempString = customObjectSourceInclinationLineEdit->text();
		if(tempString.length() == 0)  return;
		tempOrbit.i = tempString.toDouble();

		tempString = customObjectSourceNodeLineEdit->text();
		if(tempString.length() == 0)  return;
		tempOrbit.omega_node = tempString.toDouble();

		tempString = customObjectSourcePeriapsisLineEdit->text();
		if(tempString.length() == 0)  return;
		tempOrbit.omega_pericenter = tempString.toDouble();

		tempString = customObjectSourceAnomalyLineEdit->text();
		if(tempString.length() == 0)  return;
		tempOrbit.M = tempString.toDouble();

		tempOrbit.relativePosVel(tempPosition, tempVelocity);

		cout << "Object with keplerian elements added:\n";
	}

	// Create the body in memory and set its name.
	tempBody = new orsa::Body();
	tempBody->setName(tempName.toStdString());

	orsa::IBPS tempIBPS;

	// Set initialize and store the body's inertial parameters (its mass, shape, etc).
	orsa::PointLikeConstantInertialBodyProperty *tempInertial = new orsa::PointLikeConstantInertialBodyProperty(tempMass);
	tempIBPS.inertial = tempInertial;

	// Set initialize and store the body's translational parameters (position and velocity).
	orsa::DynamicTranslationalBodyProperty *tempTranslational = new orsa::DynamicTranslationalBodyProperty();
	tempTranslational->setPosition(tempPosition);
	tempTranslational->setVelocity(tempVelocity);
	tempIBPS.translational = tempTranslational;

	//TODO: Set initialize and store the body's rotational parameters.

	// Copy the IBPS into into the newloy created body.
	tempBody->setInitialConditions(tempIBPS);

	cout << "Created new body with ID " << tempBody->id() << " and name \'" << tempBody->getName() << "\'" << endl;

	//TODO:  Add the body to the object selection box.
	objectSelectionTableModel->addBody(tempBody);
}

void AddObjectsWindow::objectSelectionSelectButtonPressed()
{
	QModelIndexList selectedRows = objectSelectionTableView->selectionModel()->selectedRows();

	// Loop over the selectedRows and move them into the other table.  We sort the rows and then loop backwards over the
	// list so that when we remove items from this list it doesn't change the indexes of other elements yet to be moved over.
	qSort(selectedRows);
	for(int i = selectedRows.size()-1; i >= 0; i--)
	{
		// Add the body to the selectedObjectsTable.
		selectedObjectsTableModel->addBody(objectSelectionTableModel->getBody(selectedRows[i].row()));

		// Remove the body from this table.
		objectSelectionTableModel->removeBody(selectedRows[i].row());
	}
}

