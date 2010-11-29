#include <iostream>
using namespace std;

#include "AnalyzeIntegrationWindow.h"
#include "IntegrationTableModel.h"

AnalyzeIntegrationWindow::AnalyzeIntegrationWindow(IntegrationTableView *nSpawningWindow, int index, QWidget *parent)
	: QWidget(parent)
{
	spawningWindow = nSpawningWindow;
	bodyGroup = ((IntegrationTableModel*)spawningWindow->model())->getIntegration(index);

	resize(570, 768);

	mainGridLayout = new QGridLayout();

	tableGroupBox = new QGroupBox("Table");
	tableGroupBox->setCheckable(true);
	tableGroupBox->setChecked(false);
	mainGridLayout->addWidget(tableGroupBox, 0, 0);

	graphWidget = new AnalyzeIntegrationGraphSubwindow();
	//QObject::connect(graphReplotButton, SIGNAL(released()), this, SLOT(performAnalysisButtonPressed()));
	//mainGridLayout->addWidget(graphGroupBox, 0, 1);

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
		// Clear the old list of data accessors and create a list based on the current selection and axis settings.
		// Also create a plot curve for each item and attach it to the graph.
		graphDataAccessors.clear();

		// Clear the old list of graphPlotCurves and detach them from the graph.
		for(int i = 0; i < graphPlotCurves.size(); i++)
		{
			graphPlotCurves[i]->detach();
		}
		graphPlotCurves.clear();

		for(int i = 0; i < selectedRows.size(); i++)
		{
			graphDataAccessors.append(BodyDataAccessor(bodyGroup, objectSelectionTableModel->getBody(selectedRows[i].row())));
			graphDataAccessors[i].xDataType = (BodyDataAccessor::GraphType)(graphTypeXComboBox->itemData(graphTypeXComboBox->currentIndex()).toInt());
			graphDataAccessors[i].yDataType = (BodyDataAccessor::GraphType)(graphTypeYComboBox->itemData(graphTypeYComboBox->currentIndex()).toInt());

			graphPlotCurves.append(new QwtPlotCurve());
			graphPlotCurves[i]->setData(graphDataAccessors[i]);
			graphPlotCurves[i]->attach(graph);
		}

		graph->replot();
	}
}

/* ---------------------------------------------------------------------------- */

BodyDataAccessor::BodyDataAccessor(const orsa::BodyGroup *nBodyGroup, const orsa::Body *nBody)
{
	bodyGroup = nBodyGroup;
	body = nBody;
	bodyInterval = bodyGroup->getBodyInterval(body);

	xDataType = XData;
	yDataType = YData;
}

BodyDataAccessor::~BodyDataAccessor()
{
}

double BodyDataAccessor::getDataItem(Axis axis, size_t i) const
{
	GraphType dataType;

	switch(axis)
	{
		case XAxis:  dataType = xDataType;  break;
		case YAxis:  dataType = yDataType;  break;
	}

	switch(dataType)
	{
		case XData:  return bodyInterval->getData()[i].translational->position().getX();  break;
		case YData:  return bodyInterval->getData()[i].translational->position().getY();  break;
		case ZData:  return bodyInterval->getData()[i].translational->position().getZ();  break;

		case VXData:  return bodyInterval->getData()[i].translational->velocity().getX();  break;
		case VYData:  return bodyInterval->getData()[i].translational->velocity().getY();  break;
		case VZData:  return bodyInterval->getData()[i].translational->velocity().getZ();  break;

		case TData:   return bodyInterval->getData()[i].time.get().getMuSec().get_d();  break;
	}
}

QwtData* BodyDataAccessor::copy() const
{
	return new BodyDataAccessor(*this);
}

size_t BodyDataAccessor::size() const
{
	return bodyInterval->size();
}

double BodyDataAccessor::x(size_t i) const
{
	return getDataItem(XAxis, i);
}

double BodyDataAccessor::y(size_t i) const
{
	return getDataItem(YAxis, i);
}

AnalyzeIntegrationSubwindow::AnalyzeIntegrationSubwindow()
	: QWidget()
{
}

AnalyzeIntegrationGraphSubwindow::AnalyzeIntegrationGraphSubwindow()
	: AnalyzeIntegrationSubwindow()
{
	graphGroupBox = new QGroupBox("Graph", NULL);
	resize(960, 768);
	show();
	graphWidgetGridLayout = new QGridLayout();
	graphWidgetGridLayout->addWidget(graphGroupBox, 0, 1);
	setLayout(graphWidgetGridLayout);
	graphGroupBox->setCheckable(true);
	graphGroupBox->setChecked(false);
	graphGridLayout = new QGridLayout();
	graphReplotButton = new QPushButton("Replot");

	graphTypeXLabel = new QLabel("Horizontal\nAxis");
	//TODO:  Make a QWidget for this.
	graphTypeXComboBox = new QComboBox();
	graphTypeXComboBox->addItem("X", BodyDataAccessor::XData);
	graphTypeXComboBox->addItem("Y", BodyDataAccessor::YData);
	graphTypeXComboBox->addItem("Z", BodyDataAccessor::ZData);
	graphTypeXComboBox->addItem("VX", BodyDataAccessor::VXData);
	graphTypeXComboBox->addItem("VY", BodyDataAccessor::VYData);
	graphTypeXComboBox->addItem("VZ", BodyDataAccessor::VZData);
	graphTypeXComboBox->addItem("Time", BodyDataAccessor::TData);
	graphTypeYLabel = new QLabel("Vertical\nAxis");
	graphTypeYComboBox = new QComboBox();
	graphTypeYComboBox->addItem("X", BodyDataAccessor::XData);
	graphTypeYComboBox->addItem("Y", BodyDataAccessor::YData);
	graphTypeYComboBox->addItem("Z", BodyDataAccessor::ZData);
	graphTypeYComboBox->addItem("VX", BodyDataAccessor::VXData);
	graphTypeYComboBox->addItem("VY", BodyDataAccessor::VYData);
	graphTypeYComboBox->addItem("VZ", BodyDataAccessor::VZData);
	graphTypeYComboBox->addItem("Time", BodyDataAccessor::TData);
	graphTypeYComboBox->setCurrentIndex(1);
	graph = new QwtPlot();
	graphGridLayout->addWidget(graphReplotButton, 0, 0);
	graphGridLayout->addWidget(graphTypeXLabel, 0, 1);
	graphGridLayout->addWidget(graphTypeXComboBox, 0, 2);
	graphGridLayout->addWidget(graphTypeYLabel, 0, 3);
	graphGridLayout->addWidget(graphTypeYComboBox, 0, 4);
	graphGridLayout->addWidget(graph, 3, 0, 1, -1);
	graphGroupBox->setLayout(graphGridLayout);
}

