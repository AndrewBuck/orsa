#include <sstream>
using namespace std;

#include "IntegrationTableModel.h"

IntegrationTableModel::IntegrationTableModel(QObject *parent)
	: QAbstractTableModel(parent)
{
}

int IntegrationTableModel::rowCount(const QModelIndex & parent) const
{
	return integrationList.size();
}

int IntegrationTableModel::columnCount(const QModelIndex & parent) const
{
	return 5;
}

QVariant IntegrationTableModel::data(const QModelIndex & index, int role) const
{
	if (!index.isValid())
		return QVariant();

	if (index.row() >= integrationList.size())
		return QVariant();

	orsa::Time startTime, endTime;
	stringstream tempSS;
	int y, m, d, H, M, S, ms;
	if (role == Qt::DisplayRole)
	{
		switch(index.column())
		{
			case 0: return QString(integrationList.at(index.row())->getName().c_str());       break;
			case 1: return int(integrationList.at(index.row())->getBodyList().size());        break;

			case 2:
				integrationList.at(index.row())->getGlobalInterval(startTime, endTime);
				//TODO: This should be moved into the orsa library.
				orsaSolarSystem::gregorDay(startTime, y, m, d, H, M, S, ms);
				tempSS << y << "-" << m << "-" << d << " " << H << ":" << M << ":" << S << "." << ms/1e6;
				return tempSS.str().c_str();
				break;
			case 3:
				integrationList.at(index.row())->getGlobalInterval(startTime, endTime);
				//TODO: This should be moved into the orsa library.
				orsaSolarSystem::gregorDay(endTime, y, m, d, H, M, S, ms);
				tempSS << y << "-" << m << "-" << d << " " << H << ":" << M << ":" << S << "." << ms/1e6;
				return tempSS.str().c_str();
				break;

			default: return QString("---");                                                   break;
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

QVariant IntegrationTableModel::headerData(int section, Qt::Orientation orientation, int role) const
{
	if(role != Qt::DisplayRole)
		return QVariant();

	// For the labels across the top of the QTableView.
	if(orientation == Qt::Horizontal)
	{
		switch(section)
		{
			case 0:   return QString("Name");                    break;
			case 1:   return QString("Num Objs.");               break;
			case 2:   return QString("Start Time");              break;
			case 3:   return QString("End Time");                break;
			case 4:   return QString("Time Step");               break;
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

void IntegrationTableModel::addIntegration(orsa::BodyGroup *bg)
{
	beginResetModel();

	//TODO:  Check that this integration has not already been added to the model.

	integrationList.append(bg);

	//TODO:  Could this be optimized somehow?
	reset();
	endResetModel();
}

void IntegrationTableModel::removeIntegration(int index)
{
	beginResetModel();

	integrationList.erase(integrationList.begin() + index);

	//TODO:  Could this be optimized somehow?
	reset();
	endResetModel();
}

void IntegrationTableModel::clearAndEraseAllIntegrations()
{
	//FIXME:  This needs to actually erase the bodyGroups from memory as well.
	clearAllIntegrations();
}

void IntegrationTableModel::clearAllIntegrations()
{
	beginResetModel();
	integrationList.clear();
	reset();
	endResetModel();
}

orsa::BodyGroup* IntegrationTableModel::getIntegration(int index)
{
	return integrationList[index];
}

