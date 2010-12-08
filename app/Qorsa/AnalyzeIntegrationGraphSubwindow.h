#ifndef ANALYZEINTEGRATIONGRAPHSUBWINDOW_H
#define ANALYZEINTEGRATIONGRAPHSUBWINDOW_H

#include "AnalyzeIntegrationWindow.h"

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

#endif

