#ifndef QORSAANALYZEINTEGRATIONWINDOW_H
#define QORSAANALYZEINTEGRATIONWINDOW_H

#include <QtGui/QWidget>
#include <QtGui/QLabel>
#include <QtGui/QComboBox>
#include <QtGui/QGridLayout>
#include <QtGui/QGroupBox>
#include <QtGui/QTableView>
#include <QtGui/QPushButton>

#include <qwt-qt4/qwt_data.h>
#include <qwt-qt4/qwt_plot.h>
#include <qwt-qt4/qwt_plot_curve.h>

#include <orsa/bodygroup.h>

#include "BodyTableModel.h"
#include "IntegrationTableView.h"

class BodyDataAccessor;

class AnalyzeIntegrationWindow : public QWidget
{
	Q_OBJECT

	public:
		AnalyzeIntegrationWindow(IntegrationTableView *nSpawningWindow, int index, QWidget *parent = NULL);

	signals:
		void analyzeIntegrationWindowClosed();

	protected:
		// Protected Functions
		void closeEvent(QCloseEvent *event);

	private slots:
		void okButtonPressed();
		void cancelButtonPressed();
		void performAnalysisButtonPressed();

	private:
		IntegrationTableView *spawningWindow;
		orsa::BodyGroup *bodyGroup;

		QGridLayout *mainGridLayout;

		QGroupBox *tableGroupBox;

		QGroupBox *graphGroupBox;
		QWidget *graphWidget;
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

		QPushButton *performAnalysisPushButton;

		QTableView *objectSelectionTableView;
		BodyTableModel *objectSelectionTableModel;

		QPushButton *okPushButton;
		QPushButton *cancelPushButton;
};

class BodyDataAccessor : public QwtData
{
	public:
		enum GraphType {XData=1, YData, ZData, 
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

#endif

