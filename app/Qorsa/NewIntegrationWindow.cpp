#include <iostream>
using namespace std;

#include <orsaSolarSystem/datetime.h>

#include "MainWindow.h"
#include "NewIntegrationWindow.h"

NewIntegrationWindow::NewIntegrationWindow(MainWindow *nSpawningWindow, QWidget *parent)
	: QWidget(parent)
{
	spawningWindow = nSpawningWindow;

	setWindowTitle("New Integration");

	// Create the widgets.
	newIntegrationWindowGridLayout = new QGridLayout();
	startTimeLabel = new QLabel("Start Time");
	endTimeLabel = new QLabel("End Time");
	timeStepLabel = new QLabel("Time Step (s)");
	integratorLabel = new QLabel("Integrator");
	startDateEdit = new QDateEdit();
	endDateEdit = new QDateEdit();
	startDateEdit->setCalendarPopup(true);
	endDateEdit->setCalendarPopup(true);
	startTimeEdit = new QTimeEdit();
	endTimeEdit = new QTimeEdit();
	timeStepLineEdit = new QLineEdit("1");
	okPushButton = new QPushButton("Ok");
	cancelPushButton = new QPushButton("Cancel");

	//TODO:  Move this into its own class to make a reusable widget to select an integrator.
	integratorComboBox = new QComboBox();
	integratorComboBox->addItem("Leap Frog", LeapFrogIntegrator);
	integratorComboBox->addItem("Radau", RadauIntegrator);

	// Add the widgets to the grid layout.
	newIntegrationWindowGridLayout->addWidget(startTimeLabel, 0, 0);
	newIntegrationWindowGridLayout->addWidget(endTimeLabel, 1, 0);
	newIntegrationWindowGridLayout->addWidget(timeStepLabel, 2, 0);
	newIntegrationWindowGridLayout->addWidget(startDateEdit, 0, 1);
	newIntegrationWindowGridLayout->addWidget(endDateEdit, 1, 1);
	newIntegrationWindowGridLayout->addWidget(startTimeEdit, 0, 2);
	newIntegrationWindowGridLayout->addWidget(endTimeEdit, 1, 2);
	newIntegrationWindowGridLayout->addWidget(timeStepLineEdit, 2, 1);
	newIntegrationWindowGridLayout->addWidget(integratorLabel, 3, 0);
	newIntegrationWindowGridLayout->addWidget(integratorComboBox, 3, 1, 1, 2);
	newIntegrationWindowGridLayout->addWidget(okPushButton, 4, 0);
	newIntegrationWindowGridLayout->addWidget(cancelPushButton, 4, 1);

	// Connect the Ok and Cancel buttons to their respective handler functions.
	QObject::connect(okPushButton, SIGNAL(released()), this, SLOT(okButtonPressed()));
	QObject::connect(cancelPushButton, SIGNAL(released()), this, SLOT(cancelButtonPressed()));

	setLayout(newIntegrationWindowGridLayout);
}

void NewIntegrationWindow::closeEvent(QCloseEvent *event)
{
	emit newIntegrationWindowClosed();
}

void NewIntegrationWindow::okButtonPressed()
{
	orsa::Time startTime, endTime;

	QDate tempDate;
	QTime tempTime;

	// Read in the starting time for the integration.
	tempDate = startDateEdit->date();
	tempTime = startTimeEdit->time();
	startTime = orsaSolarSystem::gregorTime(tempDate.year(), tempDate.month(), tempDate.day(), tempTime.hour(), tempTime.minute(), tempTime.second(), 0);

	// Read in the end time for the integration.
	tempDate = endDateEdit->date();
	tempTime = endTimeEdit->time();
	endTime = orsaSolarSystem::gregorTime(tempDate.year(), tempDate.month(), tempDate.day(), tempTime.hour(), tempTime.minute(), tempTime.second(), 0);

	// Read in the time step.
	orsa::Time timeStep;
	timeStep.set(0, 0, 0, timeStepLineEdit->text().toDouble(), 0);

	// Tell the main window to begin computing the integration.
	spawningWindow->performIntegration(startTime, endTime, timeStep, NewIntegrationWindow::IntegratorType(integratorComboBox->itemData(integratorComboBox->currentIndex()).toInt()));

	//TODO:  Write the destructor to clean up all the qt widgets we created.
	close();
}

void NewIntegrationWindow::cancelButtonPressed()
{
	close();
}

