#include <iostream>
using namespace std;

#include "AnalyzeIntegrationWindow.h"
#include "IntegrationTableModel.h"

AnalyzeIntegrationWindow::AnalyzeIntegrationWindow(IntegrationTableView *nSpawningWindow, int index, QWidget *parent)
	: QWidget(parent)
{
	spawningWindow = nSpawningWindow;
	bodyGroup = ((IntegrationTableModel*)spawningWindow->model())->getIntegration(index);

	resize(800, 600);

	mainGridLayout = new QGridLayout();

	tableGroupBox = new QGroupBox("Table");
	tableGroupBox->setCheckable(true);
	tableGroupBox->setChecked(false);
	mainGridLayout->addWidget(tableGroupBox, 0, 0);

	graphGroupBox = new QGroupBox("Graph");
	graphGroupBox->setCheckable(true);
	graphGroupBox->setChecked(false);
	graphGridLayout = new QGridLayout();
	graphTypeXLabel = new QLabel("Graph Type X");
	graphTypeXComboBox = new QComboBox();
	graphTypeXComboBox->addItem("X", BodyDataAccessor::XData);
	graphTypeXComboBox->addItem("Y", BodyDataAccessor::YData);
	graphTypeXComboBox->addItem("Z", BodyDataAccessor::ZData);
	graphTypeYLabel = new QLabel("Graph Type Y");
	graphTypeYComboBox = new QComboBox();
	graphTypeYComboBox->addItem("Y", BodyDataAccessor::YData);
	graphTypeYComboBox->addItem("Y", BodyDataAccessor::YData);
	graphTypeYComboBox->addItem("Z", BodyDataAccessor::ZData);
	graphTypeYComboBox->setCurrentIndex(1);
	graph = new QwtPlot();
	graphGridLayout->addWidget(graphTypeXLabel, 0, 0);
	graphGridLayout->addWidget(graphTypeXComboBox, 0, 1);
	graphGridLayout->addWidget(graphTypeYLabel, 0, 2);
	graphGridLayout->addWidget(graphTypeYComboBox, 0, 3);
	graphGridLayout->addWidget(graph, 3, 0, 1, -1);
	graphGroupBox->setLayout(graphGridLayout);
	mainGridLayout->addWidget(graphGroupBox, 0, 1);

	performAnalysisPushButton = new QPushButton("Perform Analysis on Selected Bodies");
	QObject::connect(performAnalysisPushButton, SIGNAL(released()), this, SLOT(performAnalysisButtonPressed()));
	mainGridLayout->addWidget(performAnalysisPushButton, 2, 0, 1, -1);

	objectSelectionTableView = new QTableView();
	objectSelectionTableModel = new BodyTableModel();
	objectSelectionTableView->setModel(objectSelectionTableModel);
	objectSelectionTableView->setSelectionBehavior(QAbstractItemView::SelectRows);
	objectSelectionTableView->setSelectionMode(QAbstractItemView::MultiSelection);
	orsa::BodyGroup::BodyList bl = bodyGroup->getBodyList();
	for(unsigned int i = 0; i < bl.size(); i++)
	{
		objectSelectionTableModel->addBody(bl[i]);
	}

	mainGridLayout->addWidget(objectSelectionTableView, 3, 0, 1, -1);

	okPushButton = new QPushButton("Ok");
	cancelPushButton = new QPushButton("Cancel");

	mainGridLayout->addWidget(okPushButton, 5, 0);
	mainGridLayout->addWidget(cancelPushButton, 5, 1);

	setLayout(mainGridLayout);

	QObject::connect(okPushButton, SIGNAL(released()), this, SLOT(okButtonPressed()));
	QObject::connect(cancelPushButton, SIGNAL(released()), this, SLOT(cancelButtonPressed()));
}

void AnalyzeIntegrationWindow::closeEvent(QCloseEvent *event)
{
	emit analyzeIntegrationWindowClosed();
}

void AnalyzeIntegrationWindow::okButtonPressed()
{
	//TODO:  Write the destructor to clean up all the qt widgets we created.
	close();
}

void AnalyzeIntegrationWindow::cancelButtonPressed()
{
	close();
}

void AnalyzeIntegrationWindow::performAnalysisButtonPressed()
{
	cout << "Performing analysis.\n";
	QModelIndexList selectedRows = objectSelectionTableView->selectionModel()->selectedRows();

	// If there are no objects selected and no output panes enabled, then there is nothing to do.
	if(selectedRows.size() == 0 || 
			(!tableGroupBox->isChecked() && !graphGroupBox->isChecked()))
	{
		cout << "Returning because no output panes are selected or no objects are selected.\n";
		return;
	}

	if(graphGroupBox->isChecked())
	{
		for(int i = 0; i < selectedRows.size(); i++)
		{
		}
	}
}

