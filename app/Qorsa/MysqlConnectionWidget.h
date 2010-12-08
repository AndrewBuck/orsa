#ifndef MYSQLCONNECTIONWIDGET_H
#define MYSQLCONNECTIONWIDGET_H

#include <QtGui/QWidget>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QPushButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QTextEdit>
#include <QtGui/QGridLayout>

class MysqlConnectionWidget : public QWidget
{
	Q_OBJECT

	public:
		MysqlConnectionWidget(bool nHasQueryField = true, QWidget *parent = NULL);

		bool isReady();

		QString getHostname();
		QString getUsername();
		QString getPassword();
		QString getDatabase();
		QString getQuery();
		int getPort();

	signals:
		void executeQueryButtonPressed();

	private slots:
		void queryButtonHandler();

	private:
		bool hasQueryField;
		QGridLayout *mysqlDatabaseSourceGridLayout;
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
};

#endif

