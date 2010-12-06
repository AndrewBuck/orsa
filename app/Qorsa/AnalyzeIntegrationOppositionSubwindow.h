#ifndef ANALYZIINTEGRATIONOPPOSITIONSUBWINDOW_H
#define ANALYZIINTEGRATIONOPPOSITIONSUBWINDOW_H

#include "AnalyzeIntegrationWindow.h"

class AnalyzeIntegrationOppositionSubwindow : public AnalyzeIntegrationSubwindow
{
	Q_OBJECT

	public:
		AnalyzeIntegrationOppositionSubwindow(AnalyzeIntegrationWindow *nSpawningWindow, QWidget *parent = NULL);

		void performAnalysis();

		class OppositionResult
		{
			public:
				orsa::Body *body;
				orsa::Time oppositionTime;
				double alpha;
		};

		class OppositionResultTableModel : public QAbstractTableModel
		{
			public:
				OppositionResultTableModel(QObject *parent = NULL);

				// Functions to enable Read-only access to the model.
				int rowCount(const QModelIndex & parent = QModelIndex()) const;
				int columnCount(const QModelIndex & parent = QModelIndex()) const;
				QVariant data(const QModelIndex & index, int role = Qt::DisplayRole) const;
				QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const;

				void sort(int column, Qt::SortOrder order = Qt::AscendingOrder);

				// My functions
				void addResult(const OppositionResult *r);
				void removeResult(int index);
				void clearAllResults();
				const OppositionResult* getResult(int index);

				// Functions to enable Read-write access to the model.
				//bool setData(const QModelIndex & index, const QVariant & value, int role = Qt::EditRole);
				//Qt::ItemFlags flags(const QModelIndex & index) const;

			private:
				QList<const OppositionResult*> resultList;
		};


	private:
		QGridLayout *mainGridLayout;
		QGroupBox *mainGroupBox;
		QGroupBox *groupBox;
		QGridLayout *gridLayout;
		QPushButton *recalculateButton;

		QLabel *observingBodyLabel;
		QComboBox *observingBodyComboBox;

		QTableView *resultsTableView;
		OppositionResultTableModel *resultsTableViewModel;

		orsa::BodyGroup *bodyGroup;
};

#endif

