#include <iostream>
using namespace std;

#include "Globals.h"
#include "MainWindow.h"
#include "AddObjectsWindow.h"
#include "NewIntegrationWindow.h"

MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
{
	// Create the workspace name and set it as the window title.
	workspaceName = "Untitled ";
	workspaceName += QString::number(untitledWorkspaceID);
	untitledWorkspaceID++;
	setWindowTitle("QOrsa: " + workspaceName);

	bodyGroup = new orsa::BodyGroup();
	bodyGroup->setName(workspaceName.toStdString());

	cout << "Creating a new workspace: " << workspaceName.toStdString() << endl;

	resize(800, 600);

	addObjectsWindow = NULL;
	newIntegrationWindow = NULL;

	// Create the menubar and menus at the top of the main window.
	menuBar = new QMenuBar(this);
	createMenus();
	setMenuBar(menuBar);

	// Create the Table View's for the objects on the left and the integrations on the right of the main window.
	objectsLabel = new QLabel("Objects in Universe");
	objectsTableView = new QTableView(this);
	objectsTableModel = new BodyTableModel();
	objectsTableView->setModel(objectsTableModel);
	objectsTableView->setShowGrid(true);
	objectsTableView->setSelectionBehavior(QAbstractItemView::SelectRows);

	integrationsLabel = new QLabel("Integrations Performed in this Universe");
	integrationsTableView = new IntegrationTableView(this);
	integrationsTableModel = new IntegrationTableModel();
	integrationsTableView->setModel(integrationsTableModel);
	integrationsTableView->setShowGrid(true);
	integrationsTableView->setSelectionBehavior(QAbstractItemView::SelectRows);

	// Create the layouts used by the main window's central widget and fill them with widgets.
	objectsVBoxLayout = new QVBoxLayout();
	objectsWidget = new QWidget();
	addRemoveObjectsHBoxLayout = new QHBoxLayout();
	addObjectsPushButton = new QPushButton("Add Objects");
	removeObjectsPushButton = new QPushButton("Remove Objects");
	removeObjectsPushButton->setEnabled(false);
	objectsWidget->setLayout(objectsVBoxLayout);
	objectsVBoxLayout->addWidget(objectsLabel);
	objectsVBoxLayout->addWidget(objectsTableView);
	addRemoveObjectsHBoxLayout->addWidget(addObjectsPushButton);
	addRemoveObjectsHBoxLayout->addWidget(removeObjectsPushButton);
	objectsVBoxLayout->addLayout(addRemoveObjectsHBoxLayout);

	integrationsVBoxLayout = new QVBoxLayout();
	integrationsWidget = new QWidget();
	newIntegrationPushButton = new QPushButton("New Integration");
	newIntegrationPushButton->setEnabled(false);
	integrationsWidget->setLayout(integrationsVBoxLayout);
	integrationsVBoxLayout->addWidget(integrationsLabel);  //TODO:  Make this and the one for objects a QGroupBox.
	integrationsVBoxLayout->addWidget(integrationsTableView);
	integrationsVBoxLayout->addWidget(newIntegrationPushButton);

	// Connect the add/remove objects to their respective button handlers.
	QObject::connect(addObjectsPushButton, SIGNAL(released()), this, SLOT(addObjects()));
	QObject::connect(removeObjectsPushButton, SIGNAL(released()), this, SLOT(removeObjects()));
	QObject::connect(newIntegrationPushButton, SIGNAL(released()), this, SLOT(newIntegration()));

	// Create the splitter widget which divides the main window into the objects and integrations tables. 
	mainWindowSplitter = new QSplitter();
	mainWindowSplitter->addWidget(objectsWidget);
	mainWindowSplitter->addWidget(integrationsWidget);

	// Initialize the dummy central widget for this window with the layout set up above.
	setCentralWidget(mainWindowSplitter);

	// Connect the X button in the window frame to the quit method so it behaves like File->quit.
	QObject::connect(this, SIGNAL(destroyed()), this, SLOT(quit()));
}

QString MainWindow::getWorkspaceName()
{
	return workspaceName;
}

void MainWindow::addBody(const orsa::Body *b)
{
		// Add the current body to the bodyList and the body table model.
		bodyGroup->addBody(b);
		objectsTableModel->addBody(b);

		// Enable the removeObject and newIntegration push buttons.
		removeObjectsPushButton->setEnabled(true);
		newIntegrationPushButton->setEnabled(true);
}

void MainWindow::performIntegration(orsa::Time startTime, orsa::Time endTime, orsa::Time timeStep, NewIntegrationWindow::IntegratorType integratorType)
{
	// Create an intance of an integrator of the type passed as an argument.
	orsa::Integrator *integrator = NULL;
	cout << "Integrator value is " << integratorType << "\n";
       	switch(integratorType)
	{
		case NewIntegrationWindow::LeapFrogIntegrator:
			cout << "Using Leap Frog integrator.\n";
			integrator = new orsa::IntegratorLeapFrog();
			break;

		case NewIntegrationWindow::RadauIntegrator:
			cout << "Using Radau integrator.\n";
			integrator = new orsa::IntegratorRadau();
			break;
	}

	cout << "Beginning integration...\n";

	// Actually carry out the integration.
	integrator->integrate(bodyGroup, startTime, endTime, timeStep);

	// Store the resulting integration into the results table and allocate a new body group to be used for the next integration.
	//TODO: Most of this would be done by something like a BodyGroup::clone() function.
	orsa::BodyGroup *tempBG = bodyGroup;
	integrationsTableModel->addIntegration(bodyGroup);
	bodyGroup = new orsa::BodyGroup();
	objectsTableModel->clearAllBodies();
	orsa::BodyGroup::BodyList bl = tempBG->getBodyList();
	for(unsigned int i = 0; i < bl.size(); i++)
	{
		orsa::Body *tempBody = new orsa::Body();

		//TODO: The Body::operator= should probably be implemented to do this.
		tempBody->setInitialConditions(bl[i]->getInitialConditions());
		tempBody->setName(bl[i]->getName());
		tempBody->propulsion = bl[i]->propulsion;
		tempBody->beta = bl[i]->beta;
		tempBody->betaSun = bl[i]->betaSun;
		tempBody->birthTime = bl[i]->birthTime;
		tempBody->deathTime = bl[i]->deathTime;
		tempBody->isLightSource = bl[i]->isLightSource;
		tempBody->nonInteractingGroup = bl[i]->nonInteractingGroup;

		objectsTableModel->addBody(tempBody);
		bodyGroup->addBody(tempBody);
	}

	/*
	cout << "\n\nPerforming integration:\nInitial conditions\nBody\tx\ty\tz\t\tvx\tvy\tvz\n";
	orsa::BodyGroup::BodyList bl = bodyGroup->getBodyList();
	for(unsigned int i = 0; i < bl.size(); i++)
	{
		cout << "\n\n";
		const orsa::Body *b = bl[i];
		orsa::BodyGroup::BodyInterval *bodyInterval = bodyGroup->getBodyInterval(b);

		deque<orsa::IBPS>::iterator itr = bodyInterval->getData().begin();
		while(itr != bodyInterval->getData().end())
		{
			cout << b->getName() << "\t";
			orsa::Vector position = (*itr).translational->position();
			orsa::Vector velocity = (*itr).translational->velocity();
			cout << position.getX() << "\t";
			cout << position.getY() << "\t";
			cout << position.getZ() << "\t\t";
			cout << velocity.getX() << "\t";
			cout << velocity.getY() << "\t";
			cout << velocity.getZ() << "\n";

			itr++;
		}
	}
	*/
}

void MainWindow::newWorkspace()
{
	MainWindow *tempWindow = new MainWindow();
	tempWindow->show();
}

void MainWindow::closeEvent(QCloseEvent *event)
{
	QString tempString = "Quit called on workspace: ";
	tempString += workspaceName;
	cout << tempString.toStdString() << endl;

	if(addObjectsWindow != NULL)
		addObjectsWindow->close();

	if(addObjectsWindow != NULL)
		newIntegrationWindow->close();

	QMainWindow::closeEvent(event);
}

void MainWindow::quit()
{
	close();
}

void MainWindow::about()
{
	cout << "About QOrsa...\n";
}

void MainWindow::addObjects()
{
	if(addObjectsWindow == NULL)
	{
		addObjectsWindow = new AddObjectsWindow(this, NULL);
		addObjectsWindow->show();
		QObject::connect(addObjectsWindow, SIGNAL(addObjectsWindowClosed()), this, SLOT(addObjectsWindowClosed()));
	}
}

void MainWindow::removeObjects()
{
	QModelIndexList selectedRows = objectsTableView->selectionModel()->selectedRows();

	// Loop over the selectedRows and remove.  We sort the rows and then loop backwards over the list so that
	// when we remove items from this list it doesn't change the indexes of other elements yet to be moved over.
	qSort(selectedRows);
	for(int i = selectedRows.size()-1; i >= 0; i--)
	{
		int index = selectedRows[i].row();

		const orsa::Body *tempBody = objectsTableModel->getBody(index);
		bodyGroup->removeBody(tempBody);
		objectsTableModel->removeBody(index);
		//TODO:  After this object has been removed from the bodyGroup it should be deleted to reclaim the memory.
	}

	if(objectsTableModel->rowCount() == 0)
	{
		removeObjectsPushButton->setEnabled(false);
		newIntegrationPushButton->setEnabled(false);
	}
}

void MainWindow::addObjectsWindowClosed()
{
	addObjectsWindow = NULL;
}

void MainWindow::newIntegration()
{
	if(newIntegrationWindow == NULL)
	{
		newIntegrationWindow = new NewIntegrationWindow(this, NULL);
		newIntegrationWindow->show();
		QObject::connect(newIntegrationWindow, SIGNAL(newIntegrationWindowClosed()), this, SLOT(newIntegrationWindowClosed()));
	}
}

void MainWindow::newIntegrationWindowClosed()
{
	newIntegrationWindow = NULL;
}

void MainWindow::createMenus()
{
	// Create each menu and populate it with the menu entries it should contain.
	fileMenu = new QMenu("File", this);
		fileMenu->addAction("New Workspace", this, SLOT(newWorkspace()));
		fileMenu->addAction("Clone Workspace");
		fileMenu->addAction("Open Workspace");
		fileMenu->addSeparator();
		fileMenu->addAction("Save Workspace");
		fileMenu->addAction("Save Workspace As");
		fileMenu->addSeparator();
		fileMenu->addAction("Close Workspace", this, SLOT(quit()));

	helpMenu = new QMenu("Help", this);
		helpMenu->addAction("QOrsa Help");
		helpMenu->addAction("About", this, SLOT(about()));

	menuBar->addMenu(fileMenu);
	menuBar->addMenu(helpMenu);
}

