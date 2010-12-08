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
#include <QtGui/QSplitter>
#include <QtGui/QSortFilterProxyModel>

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
		QSortFilterProxyModel *proxyModel;
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
		void spawnOppositionWindow();
		void spawnConjunctionWindow();
		void subwindowClosedSlot(AnalyzeIntegrationSubwindow *s);

	private:
		// Private functions
		void addSubwindow(AnalyzeIntegrationSubwindow *s);

		// Private datamembers
		IntegrationTableView *spawningWindow;

		QGridLayout *mainGridLayout;

		QPushButton *spawnGraphWindowPushButton;
		QPushButton *spawnEncounterWindowPushButton;
		QPushButton *spawnOppositionWindowPushButton;
		QPushButton *spawnConjunctionWindowPushButton;

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

	signals:
		void subwindowClosed(AnalyzeIntegrationSubwindow *s);

	protected:
		void closeEvent(QCloseEvent *event);

		AnalyzeIntegrationWindow *spawningWindow;
};

#include "AnalyzeIntegrationConjunctionSubwindow.h"
#include "AnalyzeIntegrationEncounterSubwindow.h"
#include "AnalyzeIntegrationGraphSubwindow.h"

#endif

