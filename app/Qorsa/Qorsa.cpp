#include <iostream>
using namespace std;

#include <QtGui/QApplication>

#include "MainWindow.h"

unsigned int untitledWorkspaceID = 1;

int main(int argc, char **argv)
{
	cout << "Starting Qt - ORSA...\n";

	QApplication app(argc, argv);

	MainWindow *mw = new MainWindow();
	mw->show();

	return app.exec();
}

