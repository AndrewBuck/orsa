#include <iostream>
#include <sstream>

#include <QtSql/QSqlDatabase>
#include <QtSql/QSqlQuery>
#include <QtSql/QSqlError>

#include "AnalyzeIntegrationEncounterSubwindow.h"

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
	resultsTableView->setSortingEnabled(true);
	resultsTableView->setSelectionBehavior(QAbstractItemView::SelectRows);
	resultsTableView->setSelectionMode(QAbstractItemView::MultiSelection);

	proxyModel = new QSortFilterProxyModel(this);
	proxyModel->setSourceModel(resultsTableViewModel);
	proxyModel->setDynamicSortFilter(true);
	resultsTableView->setModel(proxyModel);

	mysqlConnectionWidget = new MysqlConnectionWidget(false);
	QObject::connect(mysqlConnectionWidget, SIGNAL(executeQueryButtonPressed()), this, SLOT(insertResultsIntoDatabase()));

	lrSplitter = new QSplitter(Qt::Horizontal);
	lrSplitter->addWidget(graph);
	lrSplitter->addWidget(resultsTableView);

	tbSplitter = new QSplitter(Qt::Vertical);
	tbSplitter->addWidget(lrSplitter);
	tbSplitter->addWidget(mysqlConnectionWidget);

	encounterGridLayout->addWidget(tbSplitter, 2, 0, 1, -1);

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
				if(bodyList[i] == spawningWindow->objectSelectionTableModel->getBody(spawningWindow->proxyModel->mapToSource(selectedRows[j]).row()))
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

			std::cout << "\nStep " << step << std::endl;
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

void AnalyzeIntegrationEncounterSubwindow::insertResultsIntoDatabase()
{
	if(!mysqlConnectionWidget->isReady())
		return;

	//TODO: Add a QComboBox that lets you select the database driver which will be passed to addDatabase().
	QSqlDatabase db = QSqlDatabase::addDatabase("QMYSQL");
	db.setHostName(mysqlConnectionWidget->getHostname());
	db.setPort(mysqlConnectionWidget->getPort());
	db.setDatabaseName(mysqlConnectionWidget->getDatabase());
	db.setUserName(mysqlConnectionWidget->getUsername());
	db.setPassword(mysqlConnectionWidget->getPassword());
	if(!db.open())
	{
		std::cerr << "ERROR: Could not open the database connection.\n";
		return;
	}

	std::cout << "Database connection opened successfully.  Performing query...\n";

	QSqlQuery query;

	QModelIndexList selectedRows = resultsTableView->selectionModel()->selectedRows();
	for(int i = 0; i < selectedRows.size(); i++)
	{
		const EncounterResult *tempResult = resultsTableViewModel->getResult(proxyModel->mapToSource(selectedRows[i]).row());

		std::cout << "\nInserting result for (" << tempResult->b1->getName() << ", ";
		std::cout << tempResult->b2->getName() << ")...\n";

		std::stringstream tempSS;
		tempSS << "insert into gravInteractions (name1, name2, date, distance, velocity) values ('";
		tempSS << tempResult->b1->getName() << "', '" << tempResult->b2->getName() << "', ";
		tempSS << orsaSolarSystem::timeToJulian(tempResult->time);
		tempSS << tempResult->getDistance() << ", " << tempResult->getRelVel(spawningWindow->bodyGroup).length() << ");";

		query.prepare(tempSS.str().c_str());
		if(query.exec())
		{
			std::cout << "Query completed successfully.";
		}
		else
		{
			std::cerr << "ERROR: Failed to perform the query.  The query was:\n\n";
			std::cerr << tempSS.str();
			std::cerr << "\n\nThe response from the DB was:\n";
			std::cerr << query.lastError().databaseText().toStdString();
		}


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
	std::stringstream tempSS;
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

	return QVariant();
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

