#ifndef ANALYZEINTEGRATIONENCOUNTERSUBWINDOW_H
#define ANALYZEINTEGRATIONENCOUNTERSUBWINDOW_H

#include "AnalyzeIntegrationWindow.h"

class AnalyzeIntegrationEncounterSubwindow : public AnalyzeIntegrationSubwindow
{
	Q_OBJECT

	public:
		AnalyzeIntegrationEncounterSubwindow(AnalyzeIntegrationWindow *nSpawningWindow, QWidget *parent = NULL);

		void performAnalysis();

		class EncounterResult
		{
			public:
				EncounterResult(orsa::BodyGroup *nBodyGroup, const orsa::Body *nb1, const orsa::Body *nb2, double nDistSquared, orsa::Time nTime)
				{ bodyGroup = nBodyGroup, b1 = nb1;   b2 = nb2;   distSquared = nDistSquared;   time = nTime; }

				const orsa::Body *b1, *b2;
				orsa::BodyGroup *bodyGroup;
				orsa::Time time;
				double getDistance() const {return sqrt(distSquared);}
				double getDistanceSquared() const {return distSquared;}
				orsa::Vector getRelVel(orsa::BodyGroup *bg) const
				{
					orsa::Vector v1, v2;
					bg->getInterpolatedVelocity(v1, b1, time);
					bg->getInterpolatedVelocity(v2, b2, time);
					return v2-v1;
				}

			private:
				double distSquared;
		};

		class EncounterResultSet : public QwtData
		{
			public:
				EncounterResultSet(orsa::BodyGroup *nBg) {bg = nBg;}

				// Functions inherited from QwtData.
				QwtData *copy() const {return new EncounterResultSet(*this);}
				size_t size() const {return results.size();}
				double x(size_t i) const {return results[i]->getDistance();}
				double y(size_t i) const {return results[i]->getRelVel(bg).length();}
				//TODO: Implement this function since the default implementation is slow.
				//QwtDoubleRect boundingRect() const;

				orsa::BodyGroup *bg;
				std::vector<EncounterResult*> results;
		};

		class ResultTableModel : public QAbstractTableModel
		{
			public:
				ResultTableModel(QObject *parent = NULL);

				// Functions to enable Read-only access to the model.
				int rowCount(const QModelIndex & parent = QModelIndex()) const;
				int columnCount(const QModelIndex & parent = QModelIndex()) const;
				QVariant data(const QModelIndex & index, int role = Qt::DisplayRole) const;
				QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const;

				//void sort(int column, Qt::SortOrder order = Qt::AscendingOrder);
				//static bool distanceComparitor(const EncounterResult *r1, const EncounterResult *r2);

				// My functions
				void addResult(const EncounterResult *r);
				void removeResult(int index);
				void clearAllResults();
				const EncounterResult* getResult(int index);

				// Functions to enable Read-write access to the model.
				//bool setData(const QModelIndex & index, const QVariant & value, int role = Qt::EditRole);
				//Qt::ItemFlags flags(const QModelIndex & index) const;

			private:
				QList<const EncounterResult*> resultList;
		};

	protected:
		QGroupBox *encounterGroupBox;
		QGridLayout *encounterWidgetGridLayout;
		QGridLayout *encounterGridLayout;
		QPushButton *encounterRecalculateButton;
		QList<BodyDataAccessor> dataAccessors;
		QLabel *numStepsLabel;
		QLabel *numResultsLabel;
		QLineEdit *numStepsLineEdit;
		QLineEdit *numResultsLineEdit;

		EncounterResultSet resultSet;
		QwtPlot *graph;
		QwtPlotCurve graphPlotCurve;
		QwtPlotZoomer *plotZoomer;

		QTableView *resultsTableView;
		ResultTableModel *resultsTableViewModel;
		QSortFilterProxyModel *proxyModel;

		QSplitter *splitter;

};

#endif

