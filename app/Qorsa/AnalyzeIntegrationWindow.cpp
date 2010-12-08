#include <iostream>
#include <list>
using namespace std;

#include "AnalyzeIntegrationWindow.h"
#include "AnalyzeIntegrationOppositionSubwindow.h"
#include "IntegrationTableModel.h"

AnalyzeIntegrationWindow::AnalyzeIntegrationWindow(IntegrationTableView *nSpawningWindow, int index, QWidget *parent)
	: QWidget(parent)
{
	spawningWindow = nSpawningWindow;
	bodyGroup = ((IntegrationTableModel*)spawningWindow->model())->getIntegration(index);

	resize(570, 768);

	mainGridLayout = new QGridLayout();

	spawnGraphWindowPushButton = new QPushButton("Graph");
	mainGridLayout->addWidget(spawnGraphWindowPushButton, 1, 0);
	QObject::connect(spawnGraphWindowPushButton, SIGNAL(released()), this, SLOT(spawnGraphWindow()));

	spawnEncounterWindowPushButton = new QPushButton("Encounter");
	mainGridLayout->addWidget(spawnEncounterWindowPushButton, 1, 1);
	QObject::connect(spawnEncounterWindowPushButton, SIGNAL(released()), this, SLOT(spawnEncounterWindow()));

	spawnOppositionWindowPushButton = new QPushButton("Opposition");
	mainGridLayout->addWidget(spawnOppositionWindowPushButton, 1, 2);
	QObject::connect(spawnOppositionWindowPushButton, SIGNAL(released()), this, SLOT(spawnOppositionWindow()));

	performAnalysisPushButton = new QPushButton("Perform Analysis on Selected Bodies");
	QObject::connect(performAnalysisPushButton, SIGNAL(released()), this, SLOT(performAnalysisButtonPressed()));
	mainGridLayout->addWidget(performAnalysisPushButton, 2, 0, 1, -1);

	objectSelectionTableView = new QTableView();
	objectSelectionTableModel = new BodyTableModel(bodyGroup);

	proxyModel = new QSortFilterProxyModel(this);
	proxyModel->setSourceModel(objectSelectionTableModel);
	proxyModel->setDynamicSortFilter(true);
	objectSelectionTableView->setModel(proxyModel);
	objectSelectionTableView->setSortingEnabled(true);

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
	// Close all of the subwindows that were opened as a result of this one.
	while(subwindows.size() > 0)
	{
		subwindows[0]->close();
	}

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

	// If there are no objects selected and no output panes enabled, then there is nothing to do.
	for(unsigned int i = 0; i < subwindows.size(); i++)
	{
		subwindows[i]->performAnalysis();
	}
}

void AnalyzeIntegrationWindow::spawnGraphWindow()
{
	addSubwindow(new AnalyzeIntegrationGraphSubwindow(this));
}

void AnalyzeIntegrationWindow::spawnEncounterWindow()
{
	addSubwindow(new AnalyzeIntegrationEncounterSubwindow(this));
}

void AnalyzeIntegrationWindow::spawnOppositionWindow()
{
	addSubwindow(new AnalyzeIntegrationOppositionSubwindow(this));
}

void AnalyzeIntegrationWindow::addSubwindow(AnalyzeIntegrationSubwindow *s)
{
	subwindows.push_back(s);
	QObject::connect(s, SIGNAL(subwindowClosed(AnalyzeIntegrationSubwindow*)), this, SLOT(subwindowClosedSlot(AnalyzeIntegrationSubwindow*)));
}

void AnalyzeIntegrationWindow::subwindowClosedSlot(AnalyzeIntegrationSubwindow *s)
{
	for(unsigned int i = 0; i < subwindows.size(); i++)
	{
		if(subwindows[i] == s)
		{
			subwindows.erase(subwindows.begin()+i);
			break;
		}
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

		// This really should never happen, this is just here to prevent the compiler from complaining.
		default:     dataType = xDataType;  break;
	}

	switch(dataType)
	{
		case Zero:   return 0;                                                            break;

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

AnalyzeIntegrationSubwindow::AnalyzeIntegrationSubwindow(AnalyzeIntegrationWindow *nSpawningWindow, QWidget *parent)
	: QWidget(parent)
{
	spawningWindow = nSpawningWindow;
}

void AnalyzeIntegrationSubwindow::performAnalysis()
{
}

void AnalyzeIntegrationSubwindow::closeEvent(QCloseEvent *event)
{
	emit subwindowClosed(this);
}

