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

	// If there are no objects selected and no output panes enabled, then there is nothing to do.
	for(unsigned int i = 0; i < subwindows.size(); i++)
	{
		subwindows[i]->performAnalysis();
	}
}

void AnalyzeIntegrationWindow::spawnGraphWindow()
{
	subwindows.push_back(new AnalyzeIntegrationGraphSubwindow(this));
}

void AnalyzeIntegrationWindow::spawnEncounterWindow()
{
	subwindows.push_back(new AnalyzeIntegrationEncounterSubwindow(this));
}

void AnalyzeIntegrationWindow::spawnOppositionWindow()
{
	subwindows.push_back(new AnalyzeIntegrationOppositionSubwindow(this));
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
			graphDataAccessors.append(BodyDataAccessor(spawningWindow->bodyGroup, spawningWindow->objectSelectionTableModel->getBody(selectedRows[i].row())));
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

AnalyzeIntegrationEncounterSubwindow::AnalyzeIntegrationEncounterSubwindow(AnalyzeIntegrationWindow *nSpawningWindow, QWidget *parent)
	: AnalyzeIntegrationSubwindow(nSpawningWindow, parent),
	resultSet(nSpawningWindow->bodyGroup)
{
	encounterGroupBox = new QGroupBox("Encounter", NULL);
	resize(960, 768);
	//setTitle("Graph");
	encounterWidgetGridLayout = new QGridLayout();
	encounterWidgetGridLayout->addWidget(encounterGroupBox, 0, 1);
	setLayout(encounterWidgetGridLayout);
	encounterGroupBox->setCheckable(true);
	encounterGroupBox->setChecked(true);
	encounterGridLayout = new QGridLayout();
	encounterRecalculateButton = new QPushButton("Recalculate");
	QObject::connect(encounterRecalculateButton, SIGNAL(released()), spawningWindow, SLOT(performAnalysisButtonPressed()));

	encounterGridLayout->addWidget(encounterRecalculateButton, 0, 0);

	numStepsLabel = new QLabel("Num Steps:");
	numResultsLabel = new QLabel("Num Results\nper Step:");
	numStepsLineEdit = new QLineEdit("10");
	numResultsLineEdit = new QLineEdit("50");

	encounterGridLayout->addWidget(numStepsLabel, 0, 1);
	encounterGridLayout->addWidget(numStepsLineEdit, 0, 2);
	encounterGridLayout->addWidget(numResultsLabel, 0, 3);
	encounterGridLayout->addWidget(numResultsLineEdit, 0, 4);

	graph = new QwtPlot();
	graph->setAxisTitle(QwtPlot::xBottom, "Relative Distance");
	graph->setAxisTitle(QwtPlot::yLeft, "Relative Velocity");
	plotZoomer = new QwtPlotZoomer(graph->canvas());

	resultsTableView = new QTableView();
	resultsTableViewModel = new ResultTableModel();
	resultsTableView->setModel(resultsTableViewModel);
	resultsTableView->setSortingEnabled(true);

	splitter = new QSplitter(Qt::Horizontal);
	splitter->addWidget(graph);
	splitter->addWidget(resultsTableView);

	encounterGridLayout->addWidget(splitter, 2, 0, 1, -1);

	encounterGroupBox->setLayout(encounterGridLayout);

	show();
}

void AnalyzeIntegrationEncounterSubwindow::performAnalysis()
{
	QModelIndexList selectedRows = spawningWindow->objectSelectionTableView->selectionModel()->selectedRows();

	unsigned int numSteps = atoi(numStepsLineEdit->text().toStdString().c_str());
	unsigned int numResults = atoi(numResultsLineEdit->text().toStdString().c_str());
	if(encounterGroupBox->isChecked())
	{
		resultSet.results.clear();
		orsa::BodyGroup::BodyList bodyList = spawningWindow->bodyGroup->getBodyList();

		// Loop over the body list and remove the objects not in the selectedRows.
		for(unsigned int i = 0; i < bodyList.size(); i++)
		{
			bool bodyFound = false;
			for(int j = 0; j < selectedRows.size(); j++)
			{
				if(bodyList[i] == spawningWindow->objectSelectionTableModel->getBody(selectedRows[j].row()))
				{
					bodyFound = true;
					break;
				}
			}

			if(!bodyFound)
			{
				bodyList.erase(bodyList.begin() + i);
				// Decrement i to cancel the effect of the auto increment done by the loop.
				i--;
			}
		}

		// Loop over the common time interval of the objects in the body group.
		//FIXME: This really should be the common interval of just the bodies in the subset we are investigating.
		orsa::Time startTime, stopTime;
		orsa::Time currentTime, timeStep;
		spawningWindow->bodyGroup->getCommonInterval(startTime, stopTime, false);
		currentTime = startTime;
		timeStep = (stopTime - startTime) / (double)numSteps;
		for(unsigned int step = 0; step < numSteps; step++)
		{
			currentTime = startTime + step*timeStep;

			std::cout << "\nStep " << step << endl;
			orsa::print(currentTime);

			// Loop over all pairs of objects and compute their separation distance.
			std::list<EncounterResult*> closeEncounters;
			//FIXME: Since this list contains pointers when it is cleared the objects pointed to should eventually be deleted.
			closeEncounters.clear();
			for(unsigned int i = 0; i < bodyList.size(); i++)
			{
				for(unsigned int j = i+1; j < bodyList.size(); j++)
				{
					// Compute the distance between bodies i and j.
					orsa::Vector r1, r2;
					spawningWindow->bodyGroup->getInterpolatedPosition(r1, bodyList[i], currentTime);
					spawningWindow->bodyGroup->getInterpolatedPosition(r2, bodyList[j], currentTime);

					double dx = r2.getX() - r1.getX();
					double dy = r2.getY() - r1.getY();
					double dz = r2.getZ() - r1.getZ();

					double distSquared = dx*dx + dy*dy + dz*dz;

					if(closeEncounters.size() == 0)
					{
						closeEncounters.push_back(new EncounterResult(spawningWindow->bodyGroup, bodyList[i], bodyList[j], distSquared, currentTime));
					}
					else
					{
						// Loop backwards over the list of closest encounters seen thus far this time step and try to find where this one fits in.
						std::list<EncounterResult*>::reverse_iterator itr = closeEncounters.rbegin();
						std::list<EncounterResult*>::reverse_iterator nextItr;
						std::list<EncounterResult*>::iterator forwardItr;
						while(itr != closeEncounters.rend())
						{
							nextItr = itr;
							nextItr++;
							if(nextItr == closeEncounters.rend())
							{
								if(distSquared < (*itr)->getDistanceSquared())
								{
									closeEncounters.push_front(new EncounterResult(spawningWindow->bodyGroup, bodyList[i], bodyList[j], distSquared, currentTime));
								}

								break;
							}

							if((*nextItr)->getDistanceSquared() <= distSquared && distSquared <= (*itr)->getDistanceSquared())
							{
								closeEncounters.insert(nextItr.base(), new EncounterResult(spawningWindow->bodyGroup, bodyList[i], bodyList[j], distSquared, currentTime));
								break;
							}

							itr++;
						}

						if(closeEncounters.size() > numResults)
							closeEncounters.pop_back();

					}
					
				}
			}

			// Now that we have computed the list of close encounters for this time step, print out the results.
			std::list<EncounterResult*>::iterator itr = closeEncounters.begin();
			std::cout << "\nDistance\tRel. Vel.\tBody 1 - Body 2";
			std::cout << "\n--------\t---------\t------   ------\n";
			while(itr != closeEncounters.end())
			{
				EncounterResult *enc = *itr;
				std::cout << "\n" << enc->getDistance() << "\t" << enc->getRelVel(spawningWindow->bodyGroup).length();
				std::cout << "\t" << enc->b1->getName() << " - " << enc->b2->getName();

				resultSet.results.push_back(*itr);

				// Add the result to the results table.
				resultsTableViewModel->addResult(*itr);

				itr++;
			}

			std::cout.flush();
		}

		// Attach the resultSet (which implements the QwtData interface) to the graph.
		graphPlotCurve.detach();
		graphPlotCurve.setData(resultSet);
		graphPlotCurve.setStyle(QwtPlotCurve::Dots);
		QPen pen = graphPlotCurve.pen();
		pen.setWidth(3);
		graphPlotCurve.setPen(pen);
		graphPlotCurve.attach(graph);
		graph->setAxisAutoScale(QwtPlot::xBottom);
		graph->setAxisAutoScale(QwtPlot::yLeft);
		graph->replot();
		plotZoomer->setZoomBase();
	}
}

AnalyzeIntegrationEncounterSubwindow::ResultTableModel::ResultTableModel(QObject *parent)
	: QAbstractTableModel(parent)
{
}

int AnalyzeIntegrationEncounterSubwindow::ResultTableModel::rowCount(const QModelIndex & parent) const
{
	return resultList.size();
}

int AnalyzeIntegrationEncounterSubwindow::ResultTableModel::columnCount(const QModelIndex & parent) const
{
	return 5;
}

QVariant AnalyzeIntegrationEncounterSubwindow::ResultTableModel::data(const QModelIndex & index, int role) const
{
	orsa::Time tempTime;
	stringstream tempSS;
	int y, m, d, H, M, S, ms;

	if (!index.isValid())
		return QVariant();

	if (index.row() >= resultList.size())
		return QVariant();

	if (role == Qt::DisplayRole)
	{
		switch(index.column())
		{
			case 0:
				return QString(resultList.at(index.row())->b1->getName().c_str());
				break;
			case 1:
				return QString(resultList.at(index.row())->b2->getName().c_str());
				break;

			case 2:
				tempTime = resultList.at(index.row())->time;
				orsaSolarSystem::gregorDay(tempTime, y, m, d, H, M, S, ms);
				tempSS << y << "-" << m << "-" << d << " " << H << ":" << M << ":" << S << "." << ms/1e6;
				return tempSS.str().c_str();
				break;

			case 3:
				return resultList.at(index.row())->getDistance();
				break;

			case 4:
				return resultList.at(index.row())->getRelVel(resultList.at(index.row())->bodyGroup).length();
				break;

			default:  return QString("---");                                                     break;
		}

		//QVariant tempVariant(QVariant::UserType);
		//tempVariant.setValue(bodyList.at(index.row()));
		//return tempVariant;
	}
	else
	{
		return QVariant();
	}
}

QVariant AnalyzeIntegrationEncounterSubwindow::ResultTableModel::headerData(int section, Qt::Orientation orientation, int role) const
{
	if(role != Qt::DisplayRole)
		return QVariant();

	// For the labels across the top of the QTableView.
	if(orientation == Qt::Horizontal)
	{
		switch(section)
		{
			case 0:   return QString("Body 1");                  break;
			case 1:   return QString("Body 2");                  break;
			case 2:   return QString("Date");                    break;
			case 3:   return QString("Rel Dist");                break;
			case 4:   return QString("Rel Vel.");                break;
			default:  return QString("Column %1").arg(section);  break;
		}
	}
	// For the labels across down the side of the QTableView.
	else if(orientation == Qt::Vertical)
	{
		return section+1;
	}

	return QVariant();
}

void AnalyzeIntegrationEncounterSubwindow::ResultTableModel::sort(int column, Qt::SortOrder order)
{
}

void AnalyzeIntegrationEncounterSubwindow::ResultTableModel::addResult(const EncounterResult *r)
{
	beginResetModel();

	//TODO: Check that this result has not already been added to the model.

	resultList.append(r);

	//TODO: Could this be optimized somehow?
	reset();
	endResetModel();
}

void AnalyzeIntegrationEncounterSubwindow::ResultTableModel::removeResult(int index)
{
	beginResetModel();

	resultList.erase(resultList.begin() + index);

	//TODO:  Could this be optimized somehow?
	reset();
	endResetModel();
}

void AnalyzeIntegrationEncounterSubwindow::ResultTableModel::clearAllResults()
{
	beginResetModel();
	resultList.clear();
	reset();
	endResetModel();
}

const AnalyzeIntegrationEncounterSubwindow::EncounterResult* AnalyzeIntegrationEncounterSubwindow::ResultTableModel::getResult(int index)
{
	return resultList[index];
}

