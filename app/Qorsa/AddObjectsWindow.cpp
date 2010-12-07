#include <QtCore/QtAlgorithms>
#include <QtSql/QSqlDatabase>
#include <QtSql/QSqlRecord>

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
	mysqlDatabaseSourceUsernameLineEdit = new QLineEdit("qorsa");
	mysqlDatabaseSourcePasswordLineEdit = new QLineEdit("qorsa");
	mysqlDatabaseSourcePasswordLineEdit->setEchoMode(QLineEdit::Password);
	mysqlDatabaseSourceDatabaseLineEdit = new QLineEdit("qorsa");
	mysqlDatabaseSourceQueryTextEdit = new QTextEdit();
	mysqlDatabaseSourceQueryTextEdit->setPlainText("select * from objects\norder by id\nlimit 10\n;");
	mysqlDatabaseSourceQueryTextEdit->setAcceptRichText(false);
	mysqlDatabaseSourceQuerySpacer = new QSpacerItem(1, 15);
	mysqlDatabaseSourceExecuteQueryButton = new QPushButton("Execute Query");

	QObject::connect(mysqlDatabaseSourceExecuteQueryButton, SIGNAL(released()), this, SLOT(mysqlDatabaseSourceExecuteQueryButtonPressed()));

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
	mysqlDatabaseSourceGridLayout->addWidget(mysqlDatabaseSourceExecuteQueryButton, 5, 2, 1, 2);

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
	customObjectSourceMassLabel = new QLabel("Mass (kg)");

	customObjectSourceGridLayout->addWidget(customObjectSourceNameLabel, 0, 0);
	customObjectSourceGridLayout->addWidget(customObjectSourceMassLabel, 0, 3);

	customObjectSourceXLabel = new QLabel("X\n(m)");
	customObjectSourceYLabel = new QLabel("Y\n(m)");
	customObjectSourceZLabel = new QLabel("Z\n(m)");
	customObjectSourceVXLabel = new QLabel("VX\n(m/s)");
	customObjectSourceVYLabel = new QLabel("VY\n(m/s)");
	customObjectSourceVZLabel = new QLabel("VZ\n(m/s)");

	customObjectSourceVectorGroupBoxGridLayout->addWidget(customObjectSourceXLabel, 0, 0);
	customObjectSourceVectorGroupBoxGridLayout->addWidget(customObjectSourceYLabel, 1, 0);
	customObjectSourceVectorGroupBoxGridLayout->addWidget(customObjectSourceZLabel, 2, 0);
	customObjectSourceVectorGroupBoxGridLayout->addWidget(customObjectSourceVXLabel, 0, 3);
	customObjectSourceVectorGroupBoxGridLayout->addWidget(customObjectSourceVYLabel, 1, 3);
	customObjectSourceVectorGroupBoxGridLayout->addWidget(customObjectSourceVZLabel, 2, 3);

	customObjectSourceReferenceBodyLabel = new QLabel("Ref. Body");
	customObjectSourceEccentricityLabel = new QLabel("Ecc.");
	customObjectSourceSMALabel = new QLabel("Semi Maj.\nAxis (m)");
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
	customObjectSourceKeplerianGroupBoxGridLayout->addWidget(customObjectSourceReferenceBodyLabel, 3, 0);

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
	customObjectSourceReferenceBodyComboBox = new QComboBox();

	// Loop over the bodies in the main window body group and add any that have a name and a mass to the reference body combo box.
	const orsa::BodyGroup::BodyList *bl = &(spawningWindow->bodyGroup->getBodyList());
	for(unsigned int i = 0; i < bl->size(); i++)
	{
		if((*bl)[i]->getName().length() != 0 && (*bl)[i]->getInitialConditions().inertial->mass() > 0)
		{
			string tempName = (*bl)[i]->getName();
			//TODO: Should also add the body's unique integer ID as the "data" parameter here to avoid name conflicts.
			customObjectSourceReferenceBodyComboBox->addItem(QString(tempName.c_str()));
		}
	}

	customObjectSourceKeplerianGroupBoxGridLayout->addWidget(customObjectSourceSMALineEdit, 0, 1);
	customObjectSourceKeplerianGroupBoxGridLayout->addWidget(customObjectSourceEccentricityLineEdit, 0, 4);
	customObjectSourceKeplerianGroupBoxGridLayout->addWidget(customObjectSourceInclinationLineEdit, 1, 1);
	customObjectSourceKeplerianGroupBoxGridLayout->addWidget(customObjectSourceNodeLineEdit, 1, 4);
	customObjectSourceKeplerianGroupBoxGridLayout->addWidget(customObjectSourcePeriapsisLineEdit, 2, 1);
	customObjectSourceKeplerianGroupBoxGridLayout->addWidget(customObjectSourceAnomalyLineEdit, 2, 4);
	customObjectSourceKeplerianGroupBoxGridLayout->addWidget(customObjectSourceReferenceBodyComboBox, 3, 1, 1, 4);

	customObjectSourceEpochLabel = new QLabel("Epoch");
	customObjectSourceEpochDateEdit = new QDateEdit();
	customObjectSourceEpochTimeEdit = new QTimeEdit();
	customObjectSourceGridLayout->addWidget(customObjectSourceEpochLabel, 3, 0);
	customObjectSourceGridLayout->addWidget(customObjectSourceEpochDateEdit, 3, 1);
	customObjectSourceGridLayout->addWidget(customObjectSourceEpochTimeEdit, 3, 2);

	customObjectSourceInsertButton = new QPushButton("Add Object for Selection");
	customObjectSourceGridLayout->addWidget(customObjectSourceInsertButton, 4, 2, 1, 2);
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

double AddObjectsWindow::parseLengthUnit(double num, string unit)
{
	static bool lookupMapInitialized = false;
	static map<string,orsa::Unit::LengthUnit> lengthUnitMap;

	if(!lookupMapInitialized)
	{
		lookupMapInitialized = true;
		lengthUnitMap["Mpc"] = orsa::Unit::MPARSEC;
		lengthUnitMap["Kpc"] = orsa::Unit::KPARSEC;
		lengthUnitMap["pc"] = orsa::Unit::PARSEC;
		lengthUnitMap["ly"] = orsa::Unit::LY;
		lengthUnitMap["AU"] = orsa::Unit::AU;
		lengthUnitMap["au"] = orsa::Unit::AU;
		lengthUnitMap["km"] = orsa::Unit::KM;
		lengthUnitMap["m"] = orsa::Unit::M;
		lengthUnitMap["cm"] = orsa::Unit::CM;
	}

	map<string,orsa::Unit::LengthUnit>::iterator itr = lengthUnitMap.find(unit);
	if(itr != lengthUnitMap.end())
	{
		return orsa::Unit::FromUnits(num, itr->second);
	}
	else
	{
		cout << "\n\n\nERROR:  Unrecognized unit \'" << unit << "\' in unit conversion.\n\n\n";
		return -1;
	}
}

orsa::Time AddObjectsWindow::parseTimeUnit(double num, string unit)
{
	if(unit.compare("jd") == 0)
	{
		return orsaSolarSystem::julianToTime(num);
	}
	else if(unit.compare("mjd") == 0)
	{
		return orsaSolarSystem::julianToTime(orsaSolarSystem::MJD2JD(num));
	}

	cout << "\n\n\nERROR:  Unrecognized unit \'" << unit << "\' in unit conversion.\n\n\n";
	return orsa::Time(0);
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

void AddObjectsWindow::mysqlDatabaseSourceExecuteQueryButtonPressed()
{
	map<int, orsa::Body*> referenceBodyMap;
	vector< pair<int, orsa::Body*> > unresolvedBodies;
	stack< pair<QSqlQuery*, int> > queryStack;

	if(!(mysqlDatabaseSourceHostnameLineEdit->text().length() != 0 && \
			mysqlDatabaseSourcePortLineEdit->text().length() != 0 && \
			mysqlDatabaseSourceUsernameLineEdit->text().length() != 0 && \
			mysqlDatabaseSourcePasswordLineEdit->text().length() != 0 && \
			mysqlDatabaseSourceDatabaseLineEdit->text().length() != 0 && \
			mysqlDatabaseSourceQueryTextEdit->toPlainText().length() != 0))
	{
		cout << "Some of the required input fields were left empty, aborting query.\n";
		return;
	}

	//TODO: Add a QComboBox that lets you select the database driver which will be passed to addDatabase().
	QSqlDatabase db = QSqlDatabase::addDatabase("QMYSQL");
	db.setHostName(mysqlDatabaseSourceHostnameLineEdit->text());
	db.setPort(mysqlDatabaseSourcePortLineEdit->text().toInt());
	db.setDatabaseName(mysqlDatabaseSourceDatabaseLineEdit->text());
	db.setUserName(mysqlDatabaseSourceUsernameLineEdit->text());
	db.setPassword(mysqlDatabaseSourcePasswordLineEdit->text());
	if(!db.open())
	{
		cerr << "ERROR: Could not open the database connection.\n";
		return;
	}

	cout << "Database connection opened successfully.  Performing query...\n";

	QSqlQuery *query = new QSqlQuery();
	query->prepare(mysqlDatabaseSourceQueryTextEdit->toPlainText());
	if(!query->exec())
	{
		cerr << "ERROR: Failed to perform the query.\n";
		return;
	}

	cout << "Query completed successfully.  Returned " << query->size() << " rows.\n";

	while(query->next())
	{
		resolveSqlQueryResult(referenceBodyMap, unresolvedBodies, queryStack, query);
	}

	if(unresolvedBodies.size() > 0)
	{
		cout << "\n\n\nWARNING:  There were bodies in the query results which were not added to the object selection box since the reference parent bodies which could not be resolved.\n\n\n";
	}
}

void AddObjectsWindow::resolveSqlQueryResult(map<int, orsa::Body*> &referenceBodyMap, vector< pair<int, orsa::Body*> > &unresolvedBodies, stack< pair<QSqlQuery*, int> > &queryStack, QSqlQuery *query)
{
	orsa::Body *tempBody;
	orsa::Orbit tempOrbit;
	orsa::Vector tempPosition, tempVelocity;
	orsa::Time epochTime;

	int index;
	QString name;
	double mass, epoch;
	int bodyID, parentBodyID;

	index = query->record().indexOf("id");
	bodyID = query->value(index).toInt();

	// If this body has already been read in, skip over it.
	map<int, orsa::Body*>::iterator itr = referenceBodyMap.find(bodyID);
	if(itr != referenceBodyMap.end())
		return;

	index = query->record().indexOf("name");
	if(!query->isNull(index))
		name = query->value(index).toString();
	else
		name = "";

	index = query->record().indexOf("mass");
	if(!query->isNull(index))
		mass = query->value(index).toDouble();
	else
		mass = 0;

	index = query->record().indexOf("epoch");
	epoch = query->value(index).toDouble();
	index = query->record().indexOf("epochUnit");
	string epochUnit = query->value(index).toString().toStdString();
	epochTime = parseTimeUnit(epoch, epochUnit);

	index = query->record().indexOf("sma");
	double tempSMA = query->value(index).toDouble();
	index = query->record().indexOf("smaUnit");
	string tempSMAUnit = query->value(index).toString().toStdString();
	tempOrbit.a = parseLengthUnit(tempSMA, tempSMAUnit);

	index = query->record().indexOf("eccentricity");
	tempOrbit.e = query->value(index).toDouble();

	index = query->record().indexOf("inclination");
	tempOrbit.i = query->value(index).toDouble()*M_PI/180.0;

	index = query->record().indexOf("lan");
	tempOrbit.omega_node = (M_PI/180.0)*query->value(index).toDouble();

	index = query->record().indexOf("peri");
	tempOrbit.omega_pericenter = (M_PI/180.0)*query->value(index).toDouble();

	index = query->record().indexOf("anomaly");
	tempOrbit.M = (M_PI/180.0)*query->value(index).toDouble();

	index = query->record().indexOf("referenceObject");
	if(!query->isNull(index))
		parentBodyID = query->value(index).toInt();
	else
		parentBodyID = -1;

	bool referenceBodyResolved = false;

	if(parentBodyID < 0)
	{
		// The referenceObject field was NULL so this object should
		// use the user selected default as the reference body.
		cout << name.toStdString() << ":  Using the user defined object as the reference body.\n";
	}
	else if(parentBodyID == bodyID)
	{
		// Place this object at the origin and ignore the keplerian elements.
		cout << name.toStdString() << ":  Placing this object at the origin.\n";

		tempPosition = tempVelocity = orsa::Vector(0.0, 0.0, 0.0);
		referenceBodyResolved = true;
	}
	else
	{
		// Try to resolve the parent body this object references.
		cout << name.toStdString() << ":  Resolving the reference body for this object...  ";

		map<int, orsa::Body*>::iterator itr = referenceBodyMap.find(parentBodyID);
		if(itr != referenceBodyMap.end())
		{
			// Found the parent body, compute the pos and velocity vectors from it.
			cout << "Success.\n";
			tempOrbit.mu = orsa::Unit::G() * (*itr).second->getInitialConditions().inertial->mass();
			tempOrbit.relativePosVel(tempPosition, tempVelocity);
			tempPosition += (*itr).second->getInitialConditions().translational->position();
			tempVelocity += (*itr).second->getInitialConditions().translational->velocity();

			referenceBodyResolved = true;
		}
		else
		{
			// Did not find the parent body, save the current place in the query and look at
			// the rest of the result set to see if the parent body is somewhere later.
			cout << "Failure.  Checking the rest of the result set...  ";
			queryStack.push(pair<QSqlQuery*, int>(query, query->at()));
			bool lookaheadSuccessful = false;
			while(query->next())
			{
				index = query->record().indexOf("id");
				int tempID = query->value(index).toInt();

				if(tempID == parentBodyID)
				{
					cout << "Success.\n";
					resolveSqlQueryResult(referenceBodyMap, unresolvedBodies, queryStack, query);
					map<int, orsa::Body*>::iterator itr = referenceBodyMap.find(parentBodyID);
					tempOrbit.mu = orsa::Unit::G() * (*itr).second->getInitialConditions().inertial->mass();
					tempOrbit.relativePosVel(tempPosition, tempVelocity);
					tempPosition += (*itr).second->getInitialConditions().translational->position();
					tempVelocity += (*itr).second->getInitialConditions().translational->velocity();
					referenceBodyResolved = true;

					lookaheadSuccessful = true;
					break;
				}
			}

			// Check to see if looking ahead in the result set properly resolved the parent body,
			// if not we need to perform and additional query to retrieve the parent body.
			if(!lookaheadSuccessful)
			{
				cout << "Failure.\nPerforming additional query:  ";
				query = new QSqlQuery();
				QString tempQueryString = QString("select * from objects where id = %1;").arg(parentBodyID);
				cout << tempQueryString.toStdString() << endl;
				query->prepare(tempQueryString);
				if(!query->exec())
				{
					cerr << "ERROR: Failed to perform the query.\n";
					return;
				}

				cout << "Query completed successfully.  Returned " << query->size() << " rows.\n";
				query->next();
				resolveSqlQueryResult(referenceBodyMap, unresolvedBodies, queryStack, query);
				map<int, orsa::Body*>::iterator itr = referenceBodyMap.find(parentBodyID);
				tempOrbit.mu = orsa::Unit::G() * (*itr).second->getInitialConditions().inertial->mass();
				tempOrbit.relativePosVel(tempPosition, tempVelocity);
				tempPosition += (*itr).second->getInitialConditions().translational->position();
				tempVelocity += (*itr).second->getInitialConditions().translational->velocity();
				referenceBodyResolved = true;
			}

			// Restore the query to where it was before we ran ahead through the rest of the result set.
			query = queryStack.top().first;
			query->seek(queryStack.top().second);
			queryStack.pop();
		}
	}

	tempBody = createNewBody(name, mass, epochTime, tempPosition, tempVelocity);

	if(referenceBodyResolved)
	{
		referenceBodyMap[bodyID] = tempBody;

		// Add the body to the object selection box.
		objectSelectionTableModel->addBody(tempBody);
	}
	else
	{
		unresolvedBodies.push_back(pair<int,orsa::Body*>(parentBodyID, tempBody));
	}

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
	double tempMass = orsa::FromUnits(tempString.toDouble(), orsa::Unit::KG);

	if(customObjectSourceVectorGroupBox->isChecked())
	{
		if(customObjectSourceXLineEdit->text().length() == 0 || \
			customObjectSourceYLineEdit->text().length() == 0 || \
			customObjectSourceZLineEdit->text().length() == 0 || \
			customObjectSourceVXLineEdit->text().length() == 0 || \
			customObjectSourceVYLineEdit->text().length() == 0 || \
			customObjectSourceVZLineEdit->text().length() == 0)
		{
			cout << "Not entering object since at least one required field is left blank.\n";
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
		if(customObjectSourceEccentricityLineEdit->text().length() == 0 || \
			customObjectSourceSMALineEdit->text().length() == 0 || \
			customObjectSourceInclinationLineEdit->text().length() == 0 || \
			customObjectSourceNodeLineEdit->text().length() == 0 || \
			customObjectSourcePeriapsisLineEdit->text().length() == 0 || \
			customObjectSourceAnomalyLineEdit->text().length() == 0 || \
			customObjectSourceReferenceBodyComboBox->currentIndex() == -1)
		{
			cout << "Not entering object since at least one required field is left blank.\n";
			return;
		}

		tempString = customObjectSourceEccentricityLineEdit->text();
		tempOrbit.e = tempString.toDouble();

		tempString = customObjectSourceSMALineEdit->text();
		tempOrbit.a = orsa::FromUnits(tempString.toDouble(), orsa::Unit::M);

		tempString = customObjectSourceInclinationLineEdit->text();
		tempOrbit.i = tempString.toDouble();

		tempString = customObjectSourceNodeLineEdit->text();
		tempOrbit.omega_node = tempString.toDouble();

		tempString = customObjectSourcePeriapsisLineEdit->text();
		tempOrbit.omega_pericenter = tempString.toDouble();

		tempString = customObjectSourceAnomalyLineEdit->text();
		tempOrbit.M = tempString.toDouble();

		const orsa::Body *refBody = spawningWindow->bodyGroup->getBody(customObjectSourceReferenceBodyComboBox->currentText().toStdString());
		tempOrbit.mu = orsa::Unit::G() * refBody->getInitialConditions().inertial->mass();
		//FIXME: These position and velocity vectors should be checked to make sure they are the same epoch.
		orsa::Vector refPosition = refBody->getInitialConditions().translational->position();
		orsa::Vector refVelocity = refBody->getInitialConditions().translational->velocity();

		tempOrbit.relativePosVel(tempPosition, tempVelocity);
		tempPosition += refPosition;
		tempVelocity += refVelocity;

		cout << "Object with keplerian elements added:\n";
	}

	QDate tempDate = customObjectSourceEpochDateEdit->date();
	QTime tempTime = customObjectSourceEpochTimeEdit->time();
	orsa::Time epochTime = orsaSolarSystem::gregorTime(tempDate.year(), tempDate.month(), tempDate.day(), tempTime.hour(), tempTime.minute(), tempTime.second(), 0);

	tempBody = createNewBody(tempName, tempMass, epochTime, tempPosition, tempVelocity);

	// Add the body to the object selection box.
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

orsa::Body* AddObjectsWindow::createNewBody(QString name, double mass, orsa::Time epoch, orsa::Vector position, orsa::Vector velocity)
{
	orsa::Body *tempBody;

	// Create the body in memory and set its name.
	tempBody = new orsa::Body();
	tempBody->setName(name.toStdString());

	orsa::IBPS tempIBPS;
	tempIBPS.time = epoch;

	// Set initialize and store the body's inertial parameters (its mass, shape, etc).
	orsa::PointLikeConstantInertialBodyProperty *tempInertial = new orsa::PointLikeConstantInertialBodyProperty(mass);
	tempIBPS.inertial = tempInertial;

	// Set initialize and store the body's translational parameters (position and velocity).
	orsa::DynamicTranslationalBodyProperty *tempTranslational = new orsa::DynamicTranslationalBodyProperty();
	tempTranslational->setPosition(position);
	tempTranslational->setVelocity(velocity);
	tempIBPS.translational = tempTranslational;

	//TODO: Set initialize and store the body's rotational parameters.

	// Copy the IBPS into into the newloy created body.
	tempBody->setInitialConditions(tempIBPS);

	cout << "Created new body with ID " << tempBody->id() << " and name \'" << tempBody->getName() << "\'";
	cout << " and mass " << tempBody->getInitialConditions().inertial->mass() << " kg." << endl;

	return tempBody;
}

