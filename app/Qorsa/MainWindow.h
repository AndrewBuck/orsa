#ifndef QORSAMAINWINDOW_H
#define QORSAMAINWINDOW_H

#include <QtGui/QWidget>
#include <QtGui/QMainWindow>
#include <QtGui/QMenuBar>
#include <QtGui/QPushButton>
#include <QtGui/QLabel>
#include <QtGui/QSplitter>
#include <QtGui/QVBoxLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QTableView>

#include <orsa/bodygroup.h>
#include <orsa/datetime.h>
#include <orsa/integrator.h>
#include <orsa/integrator_leapfrog.h>
#include <orsa/integrator_radau.h>

#include "BodyTableModel.h"
#include "IntegrationTableModel.h"
#include "NewIntegrationWindow.h"

class AddObjectsWindow;

class MainWindow : public QMainWindow
{
	Q_OBJECT

	public:
		// Public Functions
		MainWindow(QWidget *parent = NULL);
		QString getWorkspaceName();
		void addBody(orsa::Body *b);
		void performIntegration(orsa::Time startTime, orsa::Time endTime, orsa::Time timeStep, NewIntegrationWindow::IntegratorType integratorType);

		// Public variables
		//TODO:  Prtoect this via a semaphore lock.
		orsa::BodyGroup *bodyGroup;

	protected:
		// Protected Functions
		void closeEvent(QCloseEvent *event);

	private slots:
		// Private Slot Functions
		void newWorkspace();
		void quit();
		void about();
		void addObjects();
		void addObjectsWindowClosed();
		void removeObjects();
		void newIntegration();
		void newIntegrationWindowClosed();

	private:
		// Private Functions
		void createMenus();

		// Private Variables
		QMenuBar *menuBar;
		QMenu *fileMenu;
		QMenu *helpMenu;

		QWidget *objectsWidget;
		QVBoxLayout *objectsVBoxLayout;
		QHBoxLayout *addRemoveObjectsHBoxLayout;
		QLabel *objectsLabel;
		QTableView *objectsTableView;
		QPushButton *addObjectsPushButton;
		QPushButton *removeObjectsPushButton;
		BodyTableModel *objectsTableModel;

		QWidget *integrationsWidget;
		QVBoxLayout *integrationsVBoxLayout;
		QLabel *integrationsLabel;
		QTableView *integrationsTableView;
		IntegrationTableModel *integrationsTableModel;
		QPushButton *newIntegrationPushButton;

		QSplitter *mainWindowSplitter;

		AddObjectsWindow *addObjectsWindow;
		NewIntegrationWindow *newIntegrationWindow;

		QString workspaceName;
};

#endif

