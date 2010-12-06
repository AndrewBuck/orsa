#include <string>

#include "AnalyzeIntegrationOppositionSubwindow.h"
#include "IntegrationTableModel.h"

AnalyzeIntegrationOppositionSubwindow::AnalyzeIntegrationOppositionSubwindow(AnalyzeIntegrationWindow *nSpawningWindow, QWidget *parent)
	: AnalyzeIntegrationSubwindow(nSpawningWindow, parent)
{
	spawningWindow = nSpawningWindow;
	bodyGroup = spawningWindow->bodyGroup;

	resize(960, 768);

	mainGridLayout = new QGridLayout();

	groupBox = new QGroupBox("Graph", NULL);
	mainGridLayout->addWidget(groupBox, 0, 1);
	setLayout(mainGridLayout);
	groupBox->setCheckable(true);
	groupBox->setChecked(true);
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
				indexOfEarth = i;
		}

	}

	// If a body named "Earth" was found, make this the default object to observe from.
	if(indexOfEarth != -1)
		observingBodyComboBox->setCurrentIndex(indexOfEarth);

	gridLayout->addWidget(observingBodyLabel, 0, 1);
	gridLayout->addWidget(observingBodyComboBox, 0, 2);

	resultsTableView = new QTableView();
	resultsTableViewModel = new OppositionResultTableModel();
	resultsTableView->setModel(resultsTableViewModel);
	resultsTableView->setColumnWidth(0, 175);
	resultsTableView->setColumnWidth(1, 175);
	resultsTableView->setColumnWidth(2, 50);
	resultsTableView->setColumnWidth(3, 100);
	resultsTableView->setColumnWidth(4, 100);
	resultsTableView->setColumnWidth(5, 100);

	gridLayout->addWidget(resultsTableView, 1, 0, 1, -1);

	groupBox->setLayout(gridLayout);
	show();
}

void AnalyzeIntegrationOppositionSubwindow::performAnalysis()
{
	if(groupBox->isChecked())
	{
		// Find the body that we are making the observations from.
		orsa::Body::BodyID bodyID = observingBodyComboBox->itemData(observingBodyComboBox->currentIndex()).toUInt();
		const orsa::Body *observingBody = bodyGroup->getBody(bodyID);
		orsa::Vector r1, r2;

		// Loop over the selected bodies and calculate the opposition times for each body.
		QModelIndexList selectedRows = spawningWindow->objectSelectionTableView->selectionModel()->selectedRows();
		for(int i = 0; i < selectedRows.size(); i++)
		{
			const orsa::Body *tempBody = spawningWindow->objectSelectionTableModel->getBody(selectedRows[i].row());
		}
	}
}

AnalyzeIntegrationOppositionSubwindow::OppositionResultTableModel::OppositionResultTableModel(QObject *parent)
	: QAbstractTableModel(parent)
{
}

int AnalyzeIntegrationOppositionSubwindow::OppositionResultTableModel::rowCount(const QModelIndex & parent) const
{
	return resultList.size();
}

int AnalyzeIntegrationOppositionSubwindow::OppositionResultTableModel::columnCount(const QModelIndex & parent) const
{
	return 6;
}

QVariant AnalyzeIntegrationOppositionSubwindow::OppositionResultTableModel::data(const QModelIndex & index, int role) const
{
	if (!index.isValid())
		return QVariant();

	if (index.row() >= resultList.size())
		return QVariant();

	if (role == Qt::DisplayRole)
	{
		switch(index.column())
		{
			case 0:
				return QString(resultList.at(index.row())->body->getName().c_str());
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

QVariant AnalyzeIntegrationOppositionSubwindow::OppositionResultTableModel::headerData(int section, Qt::Orientation orientation, int role) const
{
	if(role != Qt::DisplayRole)
		return QVariant();

	// For the labels across the top of the QTableView.
	if(orientation == Qt::Horizontal)
	{
		switch(section)
		{
			case 0:   return QString("Name");                    break;
			case 1:   return QString("Approx Opposition Date");  break;
			case 2:   return QString("Mag");                     break;
			case 3:   return QString("App Motion");              break;
			case 4:   return QString("Rel Dist");                break;
			case 5:   return QString("Rel Vel");                 break;
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

void AnalyzeIntegrationOppositionSubwindow::OppositionResultTableModel::sort(int column, Qt::SortOrder order)
{
}

void AnalyzeIntegrationOppositionSubwindow::OppositionResultTableModel::addResult(const OppositionResult *r)
{
	beginResetModel();

	//TODO: Check that this result has not already been added to the model.

	resultList.append(r);

	//TODO: Could this be optimized somehow?
	reset();
	endResetModel();
}

void AnalyzeIntegrationOppositionSubwindow::OppositionResultTableModel::removeResult(int index)
{
	beginResetModel();

	resultList.erase(resultList.begin() + index);

	//TODO:  Could this be optimized somehow?
	reset();
	endResetModel();
}

void AnalyzeIntegrationOppositionSubwindow::OppositionResultTableModel::clearAllResults()
{
	beginResetModel();
	resultList.clear();
	reset();
	endResetModel();
}

