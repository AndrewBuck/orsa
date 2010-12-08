#ifndef ANALYZEINTEGRATIONCONJUNCTIONSUBWINDOW_H
#define ANALYZEINTEGRATIONCONJUNCTIONSUBWINDOW_H

#include "AnalyzeIntegrationWindow.h"

class AnalyzeIntegrationConjunctionSubwindow : public AnalyzeIntegrationSubwindow
{
	Q_OBJECT

	public:
		AnalyzeIntegrationConjunctionSubwindow(AnalyzeIntegrationWindow *nSpawningWindow, QWidget *parent = NULL);

		void performAnalysis();

	protected:
		QGridLayout *mainGridLayout;
		QGroupBox *mainGroupBox;
		QGridLayout *gridLayout;
		QPushButton *recalculateButton;

		QLabel *observingBodyLabel;
		QComboBox *observingBodyComboBox;

		QLabel *numTimeStepsLabel;
		QLineEdit *numTimeStepsLineEdit;
		QLabel *numResultsLabel;
		QLineEdit *numResultsLineEdit;

		orsa::BodyGroup *bodyGroup;
};




#endif

