#include <iostream>
using namespace std;

#include "IntegrationTableView.h"

IntegrationTableView::IntegrationTableView(QWidget *parent)
	: QTableView(parent)
{
	analyzeIntegrationWindow = NULL;
}

void IntegrationTableView::contextMenuEvent(QContextMenuEvent *event)
{
	QModelIndex index = indexAt(event->pos());
	if (index.isValid())
	{
		contextIndex = index;
		QMenu menu(this);
		menu.addAction("Analyze", this, SLOT(analyzeIntegration()));
		menu.addAction("Remove");
		menu.exec(event->globalPos());
	}
}

AnalyzeIntegrationWindow* IntegrationTableView::getAnalyzeIntegrationWindow()
{
	return analyzeIntegrationWindow;
}

void IntegrationTableView::closeAnalyzeIntegrationWindow()
{
	analyzeIntegrationWindow->close();
}

void IntegrationTableView::analyzeIntegration()
{
	cout << "Analyzing Integration\n"; 
	if(analyzeIntegrationWindow == NULL)
	{
		analyzeIntegrationWindow = new AnalyzeIntegrationWindow(this, contextIndex.row());
		analyzeIntegrationWindow->show();
		QObject::connect(analyzeIntegrationWindow, SIGNAL(analyzeIntegrationWindowClosed()), this, SLOT(analyzeIntegrationWindowClosed()));
	}
}

void IntegrationTableView::analyzeIntegrationWindowClosed()
{
	analyzeIntegrationWindow = NULL;
}

