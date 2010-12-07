#ifndef QORSABODYTABLEMODEL_H
#define QORSABODYTABLEMODEL_H

#include <QtCore/QAbstractTableModel>

#include <orsa/orbit.h>
#include <orsa/body.h>
#include <orsa/bodygroup.h>

class BodyTableModel : public QAbstractTableModel
{
	public:
		BodyTableModel(orsa::BodyGroup *nBodyGroup, QObject *parent = NULL);

		// Functions to enable Read-only access to the model.
		int rowCount(const QModelIndex & parent = QModelIndex()) const;
		int columnCount(const QModelIndex & parent = QModelIndex()) const;
		QVariant data(const QModelIndex & index, int role = Qt::DisplayRole) const;
		QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const;

		// My functions
		void addBody(const orsa::Body *b);
		void removeBody(int index);
		void clearAndEraseAllBodies();
		void clearAllBodies();
		const orsa::Body* getBody(int index);
		orsa::Orbit getOrbitForBody(orsa::BodyGroup *bodyGroup, const orsa::Body *b, orsa::Time t) const;

		// Functions to enable Read-write access to the model.
		//bool setData(const QModelIndex & index, const QVariant & value, int role = Qt::EditRole);
		//Qt::ItemFlags flags(const QModelIndex & index) const;

	private:
		QList<const orsa::Body*> bodyList;
		orsa::BodyGroup *bodyGroup;
};

Q_DECLARE_METATYPE(orsa::Body*)

#endif

