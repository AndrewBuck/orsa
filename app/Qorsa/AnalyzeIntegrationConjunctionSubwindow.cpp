#include <iostream>

#include "AnalyzeIntegrationConjunctionSubwindow.h"

AnalyzeIntegrationConjunctionSubwindow::AnalyzeIntegrationConjunctionSubwindow(AnalyzeIntegrationWindow *nSpawningWindow, QWidget *parent)
	: AnalyzeIntegrationSubwindow(nSpawningWindow, parent)
{
	spawningWindow = nSpawningWindow;
	bodyGroup = spawningWindow->bodyGroup;

	resize(960, 768);

	mainGridLayout = new QGridLayout();

	mainGroupBox = new QGroupBox("Conjunction", NULL);
	mainGridLayout->addWidget(mainGroupBox, 0, 1);
	setLayout(mainGridLayout);
	mainGroupBox->setCheckable(true);
	mainGroupBox->setChecked(true);
	gridLayout = new QGridLayout();
	recalculateButton = new QPushButton("Recalculate");
	QObject::connect(recalculateButton, SIGNAL(released()), spawningWindow, SLOT(performAnalysisButtonPressed()));

	gridLayout->addWidget(recalculateButton, 0, 0);

	observingBodyLabel = new QLabel("Observing Body");
	observingBodyComboBox = new QComboBox();

	// Fill the observing body combo box with the selected objects in the analysis window.
	int indexOfEarth = -1;
	BodyTableModel *tempModel = spawningWindow->objectSelectionTableModel;
	for(int i = 0; i < tempModel->rowCount(); i++)
	{
		const orsa::Body *tempBody = tempModel->getBody(i);
		std::string tempString = tempBody->getName();
		if(tempString.length() != 0)
		{
			observingBodyComboBox->addItem(QString(tempString.c_str()), QVariant(tempBody->id()));

			if(tempString.compare("Earth") == 0)
			{
				indexOfEarth = i;
				break;
			}
		}

	}

	// If a body named "Earth" was found, make this the default object to observe from.
	if(indexOfEarth != -1)
		observingBodyComboBox->setCurrentIndex(indexOfEarth);

	gridLayout->addWidget(observingBodyLabel, 0, 1);
	gridLayout->addWidget(observingBodyComboBox, 0, 2);

	numTimeStepsLabel = new QLabel("Num time steps");
	numTimeStepsLineEdit = new QLineEdit("200");
	numResultsLabel = new QLabel("Resuts per step");
	numResultsLineEdit = new QLineEdit("5");

	gridLayout->addWidget(numTimeStepsLabel, 0, 3);
	gridLayout->addWidget(numTimeStepsLineEdit, 0, 4);
	gridLayout->addWidget(numResultsLabel, 0, 5);
	gridLayout->addWidget(numResultsLineEdit, 0, 6);

	mainGroupBox->setLayout(gridLayout);
	show();
}

void AnalyzeIntegrationConjunctionSubwindow::performAnalysis()
{
	if(mainGroupBox->isChecked())
	{
		std::cout << "\n\nPerforming analysis of conjunctions...\n\n";

		QModelIndexList selectedRows = spawningWindow->objectSelectionTableView->selectionModel()->selectedRows();

		// Parse the input boxes to read in the parameters for the calculation.
		unsigned int numSteps = atoi(numTimeStepsLineEdit->text().toStdString().c_str());
		unsigned int numResults = atoi(numResultsLineEdit->text().toStdString().c_str());

		// Find the body that we are making the observations from.
		orsa::Body::BodyID bodyID = observingBodyComboBox->itemData(observingBodyComboBox->currentIndex()).toUInt();
		const orsa::Body *observingBody = bodyGroup->getBody(bodyID);

		// Loop over the body list and remove the objects not in the selectedRows.
		orsa::BodyGroup::BodyList bodyList = spawningWindow->bodyGroup->getBodyList();
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

			orsa::Vector observerPosition;
			spawningWindow->bodyGroup->getInterpolatedPosition(observerPosition, observingBody, currentTime);

			// Loop over all pairs of objects and compute their separation distance.
			//std::list<EncounterResult*> closeEncounters;
			//FIXME: Since this list contains pointers when it is cleared the objects pointed to should eventually be deleted.
			//closeEncounters.clear();
			for(unsigned int i = 0; i < bodyList.size(); i++)
			{
				for(unsigned int j = i+1; j < bodyList.size(); j++)
				{
					// Compute the distance between bodies i and j.
					orsa::Vector r1, r2;
					spawningWindow->bodyGroup->getInterpolatedPosition(r1, bodyList[i], currentTime);
					spawningWindow->bodyGroup->getInterpolatedPosition(r2, bodyList[j], currentTime);
					r1 -= observerPosition;
					r2 -= observerPosition;

					double dotProduct = orsa::internalProduct(r1, r2);
					dotProduct /= r1.length() * r2.length();
					std::cout << "\ndotProd = " << dotProduct;

					/*
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
					*/
					
				}
			}

			/*
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
			*/
		}
	}
}

