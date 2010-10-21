#ifndef QORSAADDOBJECTSWINDOW_H
#define QORSAADDOBJECTSWINDOW_H

#include <iostream>
#include <map>
#include <vector>
#include <stack>
#include <utility>
using namespace std;

#include <QtGui/QWidget>
#include <QtGui/QVBoxLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QGridLayout>
#include <QtGui/QGroupBox>
#include <QtGui/QRadioButton>
#include <QtGui/QPushButton>
#include <QtGui/QButtonGroup>
#include <QtGui/QTableView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QTextEdit>
#include <QtGui/QTabWidget>
#include <QtGui/QSpacerItem>
#include <QtGui/QDateEdit>
#include <QtGui/QTimeEdit>
#include <QtGui/QComboBox>
#include <QtSql/QSqlQuery>

#include <orsa/orbit.h>
#include <orsa/body.h>

#include "BodyTableModel.h"

class MainWindow;

class AddObjectsWindow : public QWidget
{
	Q_OBJECT

	public:
		AddObjectsWindow(MainWindow *nSpawningWindow, QWidget *nParent = NULL);

	signals:
		void addObjectsWindowClosed();

	protected:
		// Protected Functions
		void closeEvent(QCloseEvent *event);

	private slots:
		void okButtonPressed();
		void cancelButtonPressed();
		void mysqlDatabaseSourceExecuteQueryButtonPressed();
		void customObjectSourceVectorGroupBoxToggled(bool newState);
		void customObjectSourceKeplerianGroupBoxToggled(bool newState);
		void customObjectSourceInsertButtonPressed();
		void objectSelectionSelectButtonPressed();

	private:
		enum objectSourceType {localFile = 1, mysqlDatabase};

		// Private functions
		orsa::Body* createNewBody(QString name, double mass, orsa::Time epoch, orsa::Vector position, orsa::Vector velocity);
		void resolveSqlQueryResult(map<int, orsa::Body*> &referenceBodyMap, vector< pair<int, orsa::Body*> > &unresolvedBodies, stack< pair<QSqlQuery*, int> > &queryStack, QSqlQuery *query);

		MainWindow *spawningWindow;

		QVBoxLayout *addObjectsVBoxLayout;

		// Widgets for the Object Source Box
		QGroupBox *objectSourceGroupBox;
		QTabWidget *objectSourceDetailsTabWidget;
		QVBoxLayout *objectSourceDetailsTabWidgetLayout;

		// Widgets for the individual Object Sources
		// The "Local File" tab
		QGridLayout *localFileSourceGridLayout;
		QWidget *localFileSourceWidget;
		QLabel *localFileSourceFilenameLabel; 
		QLineEdit *localFileSourceLineEdit; 
		QPushButton *localFileSourceBrowseButton;

		// The "MYSQL Query" tab
		QGridLayout *mysqlDatabaseSourceGridLayout;
		QWidget *mysqlDatabaseSourceWidget;
		QLabel *mysqlDatabaseSourceHostnameLabel;
		QLabel *mysqlDatabaseSourcePortLabel;
		QLabel *mysqlDatabaseSourceUsernameLabel;
		QLabel *mysqlDatabaseSourcePasswordLabel;
		QLabel *mysqlDatabaseSourceDatabaseLabel;
		QLabel *mysqlDatabaseSourceQueryLabel;
		QLineEdit *mysqlDatabaseSourceHostnameLineEdit;
		QLineEdit *mysqlDatabaseSourcePortLineEdit;
		QLineEdit *mysqlDatabaseSourceUsernameLineEdit;
		QLineEdit *mysqlDatabaseSourcePasswordLineEdit;
		QLineEdit *mysqlDatabaseSourceDatabaseLineEdit;
		QTextEdit *mysqlDatabaseSourceQueryTextEdit;
		QSpacerItem *mysqlDatabaseSourceQuerySpacer;
		QPushButton *mysqlDatabaseSourceExecuteQueryButton;

		// The "Custom Object" tab
		QGridLayout *customObjectSourceGridLayout;
		QWidget *customObjectSourceWidget;
		QGroupBox *customObjectSourceVectorGroupBox;
		QGroupBox *customObjectSourceKeplerianGroupBox;
		QGridLayout *customObjectSourceVectorGroupBoxGridLayout;
		QGridLayout *customObjectSourceKeplerianGroupBoxGridLayout;
		QLabel *customObjectSourceEpochLabel;
		QLabel *customObjectSourceNameLabel;
		QLabel *customObjectSourceMassLabel;
		QLabel *customObjectSourceXLabel;
		QLabel *customObjectSourceYLabel;
		QLabel *customObjectSourceZLabel;
		QLabel *customObjectSourceVXLabel;
		QLabel *customObjectSourceVYLabel;
		QLabel *customObjectSourceVZLabel;
		QLabel *customObjectSourceEccentricityLabel;
		QLabel *customObjectSourceSMALabel;
		QLabel *customObjectSourceInclinationLabel;
		QLabel *customObjectSourceNodeLabel;
		QLabel *customObjectSourcePeriapsisLabel;
		QLabel *customObjectSourceAnomalyLabel;
		QLabel *customObjectSourceReferenceBodyLabel;
		QDateEdit *customObjectSourceEpochDateEdit;
		QTimeEdit *customObjectSourceEpochTimeEdit;
		QLineEdit *customObjectSourceNameLineEdit;
		QLineEdit *customObjectSourceMassLineEdit;
		QLineEdit *customObjectSourceXLineEdit;
		QLineEdit *customObjectSourceYLineEdit;
		QLineEdit *customObjectSourceZLineEdit;
		QLineEdit *customObjectSourceVXLineEdit;
		QLineEdit *customObjectSourceVYLineEdit;
		QLineEdit *customObjectSourceVZLineEdit;
		QLineEdit *customObjectSourceEccentricityLineEdit;
		QLineEdit *customObjectSourceSMALineEdit;
		QLineEdit *customObjectSourceInclinationLineEdit;
		QLineEdit *customObjectSourceNodeLineEdit;
		QLineEdit *customObjectSourcePeriapsisLineEdit;
		QLineEdit *customObjectSourceAnomalyLineEdit;
		QComboBox *customObjectSourceReferenceBodyComboBox;
		QPushButton *customObjectSourceInsertButton;

		// Widgets for the Object Selection Box
		QGroupBox *objectSelectionGroupBox;
		QVBoxLayout *objectSelectionVBoxLayout;
		QTableView *objectSelectionTableView;
		QPushButton *objectSelectionSelectButton;
		BodyTableModel *objectSelectionTableModel;

		// Widgets for the Selected Objects Box
		QGroupBox *selectedObjectsGroupBox;
		QVBoxLayout *selectedObjectsVBoxLayout;
		QTableView *selectedObjectsTableView;
		QHBoxLayout *selectedObjectsHBoxLayout;
		QPushButton *selectedObjectsOkButton;
		QPushButton *selectedObjectsCancelButton;
		BodyTableModel *selectedObjectsTableModel;
};

#endif

