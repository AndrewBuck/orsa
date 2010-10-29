#ifndef QORSAINTEGRATIONTABLEMODEL_H
#define QORSAINTEGRATIONTABLEMODEL_H

#include <QtCore/QAbstractTableModel>

#include <orsa/bodygroup.h>

class IntegrationTableModel : public QAbstractTableModel
{
	public:
		IntegrationTableModel(QObject *parent = NULL);

		// Functions to enable Read-only access to the model.
		int rowCount(const QModelIndex & parent = QModelIndex()) const;
		int columnCount(const QModelIndex & parent = QModelIndex()) const;
		QVariant data(const QModelIndex & index, int role = Qt::DisplayRole) const;
		QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const;

		// My functions
		void addIntegration(orsa::BodyGroup *bg);
		void removeIntegration(int index);
		void clearAndEraseAllIntegrations();
		void clearAllIntegrations();
		orsa::BodyGroup* getIntegration(int index);

		// Functions to enable Read-write access to the model.
		//bool setData(const QModelIndex & index, const QVariant & value, int role = Qt::EditRole);
		//Qt::ItemFlags flags(const QModelIndex & index) const;

	private:
		QList<orsa::BodyGroup*> integrationList;
};

#endif

