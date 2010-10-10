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
	integrationsTableView = new QTableView(this);

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
	//newIntegrationPushButton->setEnabled(false);
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

void MainWindow::addBody(orsa::Body *b)
{
		// Add the current body to the bodyList and the body table model.
		bodyGroup->addBody(b);
		objectsTableModel->addBody(b);

		// Enable the removeObject and newIntegration push buttons.
		removeObjectsPushButton->setEnabled(true);
		newIntegrationPushButton->setEnabled(true);
}

void MainWindow::performIntegration(orsa::Time startTime, orsa::Time endTime, orsa::Time timeStep)
{
	cout << "\n\nPerforming integration:\nInitial conditions\nBody\tx\ty\tz\t\tvx\tvy\tvz\n";

	//TODO Make the integrator user selectable.
	orsa::Integrator *integrator = new orsa::IntegratorLeapFrog();
	integrator->integrate(bodyGroup, startTime, endTime, timeStep);

	orsa::BodyGroup::BodyList bl = bodyGroup->getBodyList();
	for(unsigned int i = 0; i < bl.size(); i++)
	{
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

		orsa::Body *tempBody = objectsTableModel->getBody(index);
		bodyGroup->removeBody(tempBody);
		objectsTableModel->removeBody(index);
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

