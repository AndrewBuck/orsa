#include "MysqlConnectionWidget.h"

MysqlConnectionWidget::MysqlConnectionWidget(bool nHasQueryField, QWidget *parent)
	: QWidget(parent)
{
	hasQueryField = nHasQueryField;

	mysqlDatabaseSourceGridLayout = new QGridLayout();

	mysqlDatabaseSourceHostnameLabel = new QLabel("Hostname");
	mysqlDatabaseSourcePortLabel = new QLabel("Port");
	mysqlDatabaseSourceUsernameLabel = new QLabel("Username");
	mysqlDatabaseSourcePasswordLabel = new QLabel("Password");
	mysqlDatabaseSourceDatabaseLabel = new QLabel("Database");
	mysqlDatabaseSourceHostnameLineEdit = new QLineEdit("localhost");
	mysqlDatabaseSourcePortLineEdit = new QLineEdit("3306");
	mysqlDatabaseSourceUsernameLineEdit = new QLineEdit("qorsa");
	mysqlDatabaseSourcePasswordLineEdit = new QLineEdit("qorsa");
	mysqlDatabaseSourcePasswordLineEdit->setEchoMode(QLineEdit::Password);
	mysqlDatabaseSourceDatabaseLineEdit = new QLineEdit("qorsa");
	mysqlDatabaseSourceExecuteQueryButton = new QPushButton("Execute Query");
	QObject::connect(mysqlDatabaseSourceExecuteQueryButton, SIGNAL(released()), this, SLOT(queryButtonHandler()));

	mysqlDatabaseSourceGridLayout->addWidget(mysqlDatabaseSourceHostnameLabel, 0, 0);
	mysqlDatabaseSourceGridLayout->addWidget(mysqlDatabaseSourcePortLabel, 0, 4);
	mysqlDatabaseSourceGridLayout->addWidget(mysqlDatabaseSourceUsernameLabel, 1, 0);
	mysqlDatabaseSourceGridLayout->addWidget(mysqlDatabaseSourcePasswordLabel, 1, 2);
	mysqlDatabaseSourceGridLayout->addWidget(mysqlDatabaseSourceDatabaseLabel, 1, 4);

	mysqlDatabaseSourceGridLayout->addWidget(mysqlDatabaseSourceHostnameLineEdit, 0, 1, 1, 3);
	mysqlDatabaseSourceGridLayout->addWidget(mysqlDatabaseSourcePortLineEdit, 0, 5);
	mysqlDatabaseSourceGridLayout->addWidget(mysqlDatabaseSourceUsernameLineEdit, 1, 1);
	mysqlDatabaseSourceGridLayout->addWidget(mysqlDatabaseSourcePasswordLineEdit, 1, 3);
	mysqlDatabaseSourceGridLayout->addWidget(mysqlDatabaseSourceDatabaseLineEdit, 1, 5);
	mysqlDatabaseSourceGridLayout->addWidget(mysqlDatabaseSourceExecuteQueryButton, 5, 2, 1, 2);

	if(hasQueryField)
	{
		mysqlDatabaseSourceQueryLabel = new QLabel("Query:");
		mysqlDatabaseSourceQueryTextEdit = new QTextEdit();
		mysqlDatabaseSourceQueryTextEdit->setPlainText("select * from objects\norder by id\nlimit 10\n;");
		mysqlDatabaseSourceQueryTextEdit->setAcceptRichText(false);
		mysqlDatabaseSourceQuerySpacer = new QSpacerItem(1, 15);

		mysqlDatabaseSourceGridLayout->addItem(mysqlDatabaseSourceQuerySpacer, 2, 0);
		mysqlDatabaseSourceGridLayout->addWidget(mysqlDatabaseSourceQueryLabel, 3, 0);
		mysqlDatabaseSourceGridLayout->addWidget(mysqlDatabaseSourceQueryTextEdit, 4, 0, 1, -1);
	}

	setLayout(mysqlDatabaseSourceGridLayout);

}

void MysqlConnectionWidget::queryButtonHandler()
{
	emit  executeQueryButtonPressed();
}

bool MysqlConnectionWidget::isReady()
{
	return (mysqlDatabaseSourceHostnameLineEdit->text().length() != 0 && \
			mysqlDatabaseSourcePortLineEdit->text().length() != 0 && \
			mysqlDatabaseSourceUsernameLineEdit->text().length() != 0 && \
			mysqlDatabaseSourcePasswordLineEdit->text().length() != 0 && \
			mysqlDatabaseSourceDatabaseLineEdit->text().length() != 0 && \
			mysqlDatabaseSourceQueryTextEdit->toPlainText().length() != 0);
}

QString MysqlConnectionWidget::getHostname()
{
	return mysqlDatabaseSourceHostnameLineEdit->text();
}

QString MysqlConnectionWidget::getUsername()
{
	return mysqlDatabaseSourceUsernameLineEdit->text();
}

QString MysqlConnectionWidget::getPassword()
{
	return mysqlDatabaseSourcePasswordLineEdit->text();
}

QString MysqlConnectionWidget::getDatabase()
{
	return mysqlDatabaseSourceDatabaseLineEdit->text();
}

QString MysqlConnectionWidget::getQuery()
{
	if(hasQueryField)
		return mysqlDatabaseSourceQueryTextEdit->toPlainText();
	else
		return "";
}

int MysqlConnectionWidget::getPort()
{
	return mysqlDatabaseSourcePortLineEdit->text().toInt();
}

