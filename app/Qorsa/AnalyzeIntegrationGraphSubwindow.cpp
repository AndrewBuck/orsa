#include "AnalyzeIntegrationGraphSubwindow.h"

AnalyzeIntegrationGraphSubwindow::AnalyzeIntegrationGraphSubwindow(AnalyzeIntegrationWindow *nSpawningWindow, QWidget *parent)
	: AnalyzeIntegrationSubwindow(nSpawningWindow, parent)
{
	graphGroupBox = new QGroupBox("Graph", NULL);
	resize(960, 768);
	//setTitle("Graph");
	graphWidgetGridLayout = new QGridLayout();
	graphWidgetGridLayout->addWidget(graphGroupBox, 0, 1);
	setLayout(graphWidgetGridLayout);
	graphGroupBox->setCheckable(true);
	graphGroupBox->setChecked(true);
	graphGridLayout = new QGridLayout();
	graphReplotButton = new QPushButton("Replot");
	QObject::connect(graphReplotButton, SIGNAL(released()), spawningWindow, SLOT(performAnalysisButtonPressed()));

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
	plotZoomer = new QwtPlotZoomer(graph->canvas());

	graphGridLayout->addWidget(graphReplotButton, 0, 0);
	graphGridLayout->addWidget(graphTypeXLabel, 0, 1);
	graphGridLayout->addWidget(graphTypeXComboBox, 0, 2);
	graphGridLayout->addWidget(graphTypeYLabel, 0, 3);
	graphGridLayout->addWidget(graphTypeYComboBox, 0, 4);
	graphGridLayout->addWidget(graph, 3, 0, 1, -1);
	graphGroupBox->setLayout(graphGridLayout);

	show();
}

void AnalyzeIntegrationGraphSubwindow::performAnalysis()
{
	QModelIndexList selectedRows = spawningWindow->objectSelectionTableView->selectionModel()->selectedRows();

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
			graphDataAccessors.append(BodyDataAccessor(spawningWindow->bodyGroup, spawningWindow->objectSelectionTableModel->getBody(spawningWindow->proxyModel->mapToSource(selectedRows[i]).row())));
			graphDataAccessors[i].xDataType = (BodyDataAccessor::GraphType)(graphTypeXComboBox->itemData(graphTypeXComboBox->currentIndex()).toInt());
			graphDataAccessors[i].yDataType = (BodyDataAccessor::GraphType)(graphTypeYComboBox->itemData(graphTypeYComboBox->currentIndex()).toInt());

			graphPlotCurves.append(new QwtPlotCurve());
			graphPlotCurves[i]->setData(graphDataAccessors[i]);
			graphPlotCurves[i]->attach(graph);
		}

		graph->setAxisAutoScale(QwtPlot::xBottom);
		graph->setAxisAutoScale(QwtPlot::yLeft);
		graph->replot();
		plotZoomer->setZoomBase();
	}
}


