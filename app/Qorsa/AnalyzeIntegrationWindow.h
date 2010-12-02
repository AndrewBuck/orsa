#ifndef QORSAANALYZEINTEGRATIONWINDOW_H
#define QORSAANALYZEINTEGRATIONWINDOW_H

#include <QtGui/QWidget>
#include <QtGui/QLabel>
#include <QtGui/QComboBox>
#include <QtGui/QGridLayout>
#include <QtGui/QGroupBox>
#include <QtGui/QTableView>
#include <QtGui/QPushButton>
#include <QtGui/QLineEdit>

#include <qwt-qt4/qwt_data.h>
#include <qwt-qt4/qwt_plot.h>
#include <qwt-qt4/qwt_plot_zoomer.h>
#include <qwt-qt4/qwt_plot_curve.h>

#include <orsa/bodygroup.h>

#include "BodyTableModel.h"
#include "IntegrationTableView.h"

class BodyDataAccessor;
class AnalyzeIntegrationSubwindow;
class AnalyzeIntegrationGraphSubwindow;

class AnalyzeIntegrationWindow : public QWidget
{
	Q_OBJECT

	public:
		AnalyzeIntegrationWindow(IntegrationTableView *nSpawningWindow, int index, QWidget *parent = NULL);

		// Public datamembers
		QTableView *objectSelectionTableView;
		BodyTableModel *objectSelectionTableModel;
		orsa::BodyGroup *bodyGroup;

	signals:
		void analyzeIntegrationWindowClosed();

	protected:
		// Protected Functions
		void closeEvent(QCloseEvent *event);

	private slots:
		void okButtonPressed();
		void cancelButtonPressed();
		void performAnalysisButtonPressed();
		void spawnGraphWindow();
		void spawnEncounterWindow();

	private:
		IntegrationTableView *spawningWindow;

		QGridLayout *mainGridLayout;

		QPushButton *spawnGraphWindowPushButton;
		QPushButton *spawnEncounterWindowPushButton;

		std::vector<AnalyzeIntegrationSubwindow*> subwindows;

		QPushButton *performAnalysisPushButton;

		QPushButton *okPushButton;
		QPushButton *cancelPushButton;
};

class BodyDataAccessor : public QwtData
{
	public:
		enum GraphType {Zero=0, XData=1, YData, ZData, 
			VXData, VYData, VZData,
			TData};

		enum Axis {XAxis=1, YAxis};

		BodyDataAccessor(const orsa::BodyGroup *nBodyGroup, const orsa::Body *nBody);
		~BodyDataAccessor();

		double getDataItem(Axis axis, size_t i) const;

		// Functions inherited from QwtData.
		QwtData *copy() const;
		size_t size() const;
		double x(size_t i) const;
		double y(size_t i) const;
		//TODO: Implement this function since the default implementation is slow.
		//QwtDoubleRect boundingRect() const;

		// Public variables.
		GraphType xDataType;
		GraphType yDataType;

	private:
		const orsa::BodyGroup *bodyGroup;
		const orsa::Body *body;
		const orsa::BodyGroup::BodyInterval *bodyInterval;
};

class AnalyzeIntegrationSubwindow : public QWidget
{
	Q_OBJECT

	public:
		AnalyzeIntegrationSubwindow(AnalyzeIntegrationWindow *nSpawningWindow, QWidget *parent = NULL);

		virtual void performAnalysis();

	protected:
		AnalyzeIntegrationWindow *spawningWindow;
};

class AnalyzeIntegrationGraphSubwindow : public AnalyzeIntegrationSubwindow
{
	Q_OBJECT

	public:
		AnalyzeIntegrationGraphSubwindow(AnalyzeIntegrationWindow *nSpawningWindow, QWidget *parent = NULL);

		void performAnalysis();

	protected:
		QGroupBox *graphGroupBox;
		QGridLayout *graphWidgetGridLayout;
		QGridLayout *graphGridLayout;
		QPushButton *graphReplotButton;
		QLabel *graphTypeXLabel;
		QComboBox *graphTypeXComboBox;
		QLabel *graphTypeYLabel;
		QComboBox *graphTypeYComboBox;
		QList<BodyDataAccessor> graphDataAccessors;
		QList<QwtPlotCurve*> graphPlotCurves;
		QwtPlot *graph;
		QwtPlotZoomer *plotZoomer;
};

class AnalyzeIntegrationEncounterSubwindow : public AnalyzeIntegrationSubwindow
{
	Q_OBJECT

	public:
		AnalyzeIntegrationEncounterSubwindow(AnalyzeIntegrationWindow *nSpawningWindow, QWidget *parent = NULL);

		void performAnalysis();

		class EncounterResult
		{
			public:
				EncounterResult(const orsa::Body *nb1, const orsa::Body *nb2, double nDistSquared, orsa::Time nTime)
				{ b1 = nb1;   b2 = nb2;   distSquared = nDistSquared;   time = nTime; }

				const orsa::Body *b1, *b2;
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
				double x(size_t i) const {return results[i].getDistance();}
				double y(size_t i) const {return results[i].getRelVel(bg).length();}
				//TODO: Implement this function since the default implementation is slow.
				//QwtDoubleRect boundingRect() const;

				orsa::BodyGroup *bg;
				std::vector<EncounterResult> results;
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
};

#endif

