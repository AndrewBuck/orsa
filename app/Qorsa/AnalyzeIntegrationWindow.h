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

#include <orsa/bodygroup.h>

#include "BodyTableModel.h"
#include "IntegrationTableView.h"

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
		QGridLayout *graphGridLayout;
		QLabel *graphTypeXLabel;
		QComboBox *graphTypeXComboBox;
		QLabel *graphTypeYLabel;
		QComboBox *graphTypeYComboBox;
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

		BodyDataAccessor(orsa::BodyGroup *nBodyGroup, orsa::Body *nBody);

	private:
};

#endif

