#include "BodyTableModel.h"

BodyTableModel::BodyTableModel(QObject *parent)
	: QAbstractTableModel(parent)
{
}

int BodyTableModel::rowCount(const QModelIndex & parent) const
{
	return bodyList.size();
}

int BodyTableModel::columnCount(const QModelIndex & parent) const
{
	return 5;
}

QVariant BodyTableModel::data(const QModelIndex & index, int role) const
{
	if (!index.isValid())
		return QVariant();

	if (index.row() >= bodyList.size())
		return QVariant();

	if (role == Qt::DisplayRole)
	{
		switch(index.column())
		{
			case 0:   return QString(bodyList.at(index.row())->getName().c_str());               break;
			case 1:   return bodyList.at(index.row())->getInitialConditions().inertial->mass();  break;
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

QVariant BodyTableModel::headerData(int section, Qt::Orientation orientation, int role) const
{
	if(role != Qt::DisplayRole)
		return QVariant();

	// For the labels across the top of the QTableView.
	if(orientation == Qt::Horizontal)
	{
		switch(section)
		{
			case 0:   return QString("Name");                    break;
			case 1:   return QString("Mass");                    break;
			case 2:   return QString("S-maj axis");              break;
			case 3:   return QString("eccentricity");            break;
			case 4:   return QString("inclination");             break;
			default:  return QString("Column %1").arg(section);  break;
		}
	}
	// For the labels across down the side of the QTableView.
	else if(orientation == Qt::Vertical)
	{
		return section+1;
	}
}

void BodyTableModel::addBody(orsa::Body *b)
{
	beginResetModel();

	//TODO:  Check that this body has not already been added to the model.

	bodyList.append(b);

	//TODO:  Could this be optimized somehow?
	reset();
	endResetModel();
}

void BodyTableModel::removeBody(int index)
{
	beginResetModel();

	bodyList.erase(bodyList.begin() + index);

	//TODO:  Could this be optimized somehow?
	reset();
	endResetModel();
}

void BodyTableModel::clearAndEraseAllBodies()
{
	beginResetModel();

	//FIXME: This is commented out right now because the body destructor is protected, need to find a solution to this.
	// First loop over all the bodies in the list and delete them since the structure is no
	// longer needed and will be inaccessable once the pointer is removed from the list.
	//for(int i = 0; i < bodyList.size(); i++)
		//delete bodyList[i];

	// Now that all of the referenced objects have been deleted the list itself can be erased.
	bodyList.clear();

	reset();
	endResetModel();
}

void BodyTableModel::clearAllBodies()
{
	beginResetModel();
	bodyList.clear();
	reset();
	endResetModel();
}

orsa::Body* BodyTableModel::getBody(int index)
{
	return bodyList[index];
}

