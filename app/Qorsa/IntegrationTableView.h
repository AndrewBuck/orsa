#ifndef QORSAINTEGRATIONTABLEVIEW_H
#define QORSAINTEGRATIONTABLEVIEW_H

#include <QtGui/QWidget>
#include <QtGui/QContextMenuEvent>
#include <QtGui/QTableView>
#include <QtGui/QMenu>

class AnalyzeIntegrationWindow;

class IntegrationTableView : public QTableView
{
	Q_OBJECT

	public:
		IntegrationTableView(QWidget *parent = NULL);

		void contextMenuEvent(QContextMenuEvent *event);

	private slots:
		void analyzeIntegration();
		void analyzeIntegrationWindowClosed();

	private:
		QModelIndex contextIndex;
		AnalyzeIntegrationWindow* analyzeIntegrationWindow;
};

#include "AnalyzeIntegrationWindow.h"

#endif

