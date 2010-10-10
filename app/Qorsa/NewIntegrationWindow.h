#ifndef QORSANEWINTEGRATIONWINDOW_H
#define QORSANEWINTEGRATIONWINDOW_H

#include <QtGui/QWidget>
#include <QtGui/QGridLayout>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QDateEdit>
#include <QtGui/QTimeEdit>
#include <QtGui/QPushButton>

class MainWindow;

class NewIntegrationWindow : public QWidget
{
	Q_OBJECT

	public:
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
		QDateEdit *startDateEdit;
		QDateEdit *endDateEdit;
		QTimeEdit *startTimeEdit;
		QTimeEdit *endTimeEdit;
		QLineEdit *timeStepLineEdit;
		QPushButton *okPushButton;
		QPushButton *cancelPushButton;
};

#endif

