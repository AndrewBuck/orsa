#include <iostream>
#include <sstream>
#include <utility>
#include <cmath>

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
	if(encounterGroupBox->isChecked())
	{
		orsa::Vector r1, r2;

		QModelIndexList selectedRows = spawningWindow->objectSelectionTableView->selectionModel()->selectedRows();

		unsigned int numSteps = atoi(numStepsLineEdit->text().toStdString().c_str());
		unsigned int numResults = atoi(numResultsLineEdit->text().toStdString().c_str());

		resultSet.results.clear();
		resultsTableViewModel->clearAllResults();
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

		// Create a lattice of cells which will contain the addresses of
		// objects based on their position to allow the lookups to proceed faster.
		const int numCells = 20;  //NOTE: This should be divisible by 2 to allow for +/- in each axis.
		double xSize = orsa::Unit::FromUnits(0.5, orsa::Unit::AU);
		double ySize = orsa::Unit::FromUnits(0.5, orsa::Unit::AU);
		double zSize = orsa::Unit::FromUnits(0.5, orsa::Unit::AU);
		double xStart = -0.5 * xSize * (double)numCells;
		double yStart = -0.5 * ySize * (double)numCells;
		double zStart = -0.5 * zSize * (double)numCells;
		double xMin, xMax, yMin, yMax, zMin, zMax;
		std::vector<const orsa::Body*> cells[numCells][numCells][numCells];
		std::vector<const orsa::Body*> overflowBodies;

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
			std::cout << "\nxSize: " << xSize << "\tySize: " << ySize << "\tzSize: " << zSize;
			std::cout.flush();

			// Clear the overflow list and the list in each cell.
			overflowBodies.clear();
			for(int i = 0; i < numCells; i++)
				for(int j = 0; j < numCells; j++)
					for(int k = 0; k < numCells; k++)
						cells[i][j][k].clear();

			// Coarse grain the orbiting bodies into cells.
			//TODO:  Adjust xStart and xSize (same for other dims too) to try to keep the boxes small but with few objects in the overflow list.
			bool firstLoop = true;
			for(unsigned int i = 0; i < bodyList.size(); i++)
			{
				spawningWindow->bodyGroup->getInterpolatedPosition(r1, bodyList[i], currentTime);

				double x, y, z;
				x = r1.getX();
				y = r1.getY();
				z = r1.getZ();

				int xIndex = floor( x/xSize ) + numCells/2;
				int yIndex = floor( y/ySize ) + numCells/2;
				int zIndex = floor( z/zSize ) + numCells/2;

				//std::cout << "\nx = " << x << "\txStart = " << xStart << "\txSize = " << xSize << "\txIndex = " << xIndex;

				if(firstLoop)
				{
					firstLoop = false;
					xMin = xMax = x;
					yMin = yMax = y;
					zMin = zMax = z;
				}
				else
				{
					if(x < xMin)  xMin = x;
					if(x > xMax)  xMax = x;
					if(y < yMin)  yMin = y;
					if(y > yMax)  yMax = y;
					if(z < zMin)  zMin = z;
					if(z > zMax)  zMax = z;
				}

				if(xIndex >= 0 && yIndex >= 0 && zIndex >= 0 && xIndex < numCells && yIndex < numCells && zIndex < numCells)
					cells[xIndex][yIndex][zIndex].push_back(bodyList[i]);
				else
					overflowBodies.push_back(bodyList[i]);

			}

			std::cout << "\nThere are " << overflowBodies.size() << " of " << bodyList.size() << " bodies in the overflow list.";
			std::cout.flush();

			// Create and fill up a list of all the pairs which should be checked for close interactions.
			std::vector< std::pair<const orsa::Body*, const orsa::Body*> > interactingPairs;
			overflowBodies.clear();
			for(int i = 0; i < numCells; i++)
			{
				for(int j = 0; j < numCells; j++)
				{
					for(int k = 0; k < numCells; k++)
					{
						// Add all the pairs in this cell.
						for(unsigned int l = 0; l < cells[i][j][k].size(); l++)
							for(unsigned int m = l+1; m < cells[i][j][k].size(); m++)
								interactingPairs.push_back(std::pair<const orsa::Body*,const orsa::Body*>(cells[i][j][k][l],cells[i][j][k][m]));

						/*
						// Add all the pairs with the objects in the overflow list.
						for(unsigned int l = 0; l < cells[i][j][k].size(); l++)
						{
							if(i+1 < numCells)
								for(unsigned int m = 0; m < cells[i+1][j][k].size(); m++)
									interactingPairs.push_back(std::pair<const orsa::Body*,const orsa::Body*>(cells[i][j][k][l],cells[i][j][k][m]));

							if(j+1 < numCells)
								for(unsigned int m = 0; m < cells[i][j+1][k].size(); m++)
									interactingPairs.push_back(std::pair<const orsa::Body*,const orsa::Body*>(cells[i][j][k][l],cells[i][j][k][m]));

							if(k+1 < numCells)
								for(unsigned int m = 0; m < cells[i][j][k+1].size(); m++)
									interactingPairs.push_back(std::pair<const orsa::Body*,const orsa::Body*>(cells[i][j][k][l],cells[i][j][k][m]));
						}


						// Add all the pairs with the cells 1 "higher" in each dimension.
						for(unsigned int l = 0; l < cells[i][j][k].size(); l++)
						{
							for(unsigned int m = 0; m < overflowBodies.size(); m++)
							{
									interactingPairs.push_back(std::pair<const orsa::Body*,const orsa::Body*>(cells[i][j][k][l],overflowBodies[m]));
							}
						}
						*/
					}
				}
			}

			std::cout << "\nBeginning pair evaluation.";
			std::cout.flush();
			// Loop over all pairs of objects and compute their separation distance.
			std::list<EncounterResult*> closeEncounters;
			//FIXME: Since this list contains pointers when it is cleared the objects pointed to should eventually be deleted.
			closeEncounters.clear();
			for(unsigned int i = 0; i < interactingPairs.size(); i++)
			{
				if(interactingPairs[i].first == interactingPairs[i].second)
					continue;

				// Compute the distance between bodies i and j.
				spawningWindow->bodyGroup->getInterpolatedPosition(r1, interactingPairs[i].first, currentTime);
				spawningWindow->bodyGroup->getInterpolatedPosition(r2, interactingPairs[i].second, currentTime);

				double dx = r2.getX() - r1.getX();
				double dy = r2.getY() - r1.getY();
				double dz = r2.getZ() - r1.getZ();

				double distSquared = dx*dx + dy*dy + dz*dz;

				if(closeEncounters.size() == 0)
				{
					closeEncounters.push_back(new EncounterResult(spawningWindow->bodyGroup, interactingPairs[i].first, interactingPairs[i].second, distSquared, currentTime));
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
								closeEncounters.push_front(new EncounterResult(spawningWindow->bodyGroup, interactingPairs[i].first, interactingPairs[i].second, distSquared, currentTime));
							}

							break;
						}

						if((*nextItr)->getDistanceSquared() <= distSquared && distSquared <= (*itr)->getDistanceSquared())
						{
							closeEncounters.insert(nextItr.base(), new EncounterResult(spawningWindow->bodyGroup, interactingPairs[i].first, interactingPairs[i].second, distSquared, currentTime));
							break;
						}

						itr++;
					}

					if(closeEncounters.size() > numResults)
						closeEncounters.pop_back();

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
		tempSS << orsaSolarSystem::timeToJulian(tempResult->time) << ", ";
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
	return 6;
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

			case 5:
				double r, v;
				r = resultList.at(index.row())->getRelVel(resultList.at(index.row())->bodyGroup).length();
				v = resultList.at(index.row())->getDistance();
				return log10(r*r*v);
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
			case 5:   return QString("Score");                break;
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

