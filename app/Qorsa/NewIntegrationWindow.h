#ifndef QORSANEWINTEGRATIONWINDOW_H
#define QORSANEWINTEGRATIONWINDOW_H

#include <QtGui/QWidget>
#include <QtGui/QGridLayout>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QDateEdit>
#include <QtGui/QTimeEdit>
#include <QtGui/QComboBox>
#include <QtGui/QPushButton>

class MainWindow;

class NewIntegrationWindow : public QWidget
{
	Q_OBJECT

	public:
		enum IntegratorType {LeapFrogIntegrator=1, RadauIntegrator};
		NewIntegrationWindow(MainWindow *nSpawningWindow, QWidget *parent = NULL);

	signals:
		void newIntegrationWindowClosed();

	protected:
		// Protected Functions
		void closeEvent(QCloseEvent *event);

	private slots:
		void okButtonPressed();
		void cancelButtonPressed();

	private:
		MainWindow *spawningWindow;

		QGridLayout *newIntegrationWindowGridLayout;
		QLabel *startTimeLabel;
		QLabel *endTimeLabel;
		QLabel *timeStepLabel;
		QLabel *integratorLabel;
		QDateEdit *startDateEdit;
		QDateEdit *endDateEdit;
		QTimeEdit *startTimeEdit;
		QTimeEdit *endTimeEdit;
		QLineEdit *timeStepLineEdit;
		QComboBox *integratorComboBox;
		QPushButton *okPushButton;
		QPushButton *cancelPushButton;
};

#endif

