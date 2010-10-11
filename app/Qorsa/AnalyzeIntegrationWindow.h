#ifndef QORSAANALYZEINTEGRATIONWINDOW_H
#define QORSAANALYZEINTEGRATIONWINDOW_H

#include <QtGui/QWidget>
#include <QtGui/QGridLayout>
#include <QtGui/QGroupBox>
#include <QtGui/QTableView>
#include <QtGui/QPushButton>

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

	private:
		IntegrationTableView *spawningWindow;
		orsa::BodyGroup *bodyGroup;

		QGridLayout *mainGridLayout;

		QGroupBox *tableGroupBox;
		QGroupBox *graphGroupBox;

		QTableView *objectSelectionTableView;
		BodyTableModel *objectSelectionTableModel;

		QPushButton *okPushButton;
		QPushButton *cancelPushButton;
};

#endif

