#include "RendezvousWithVesta.h"

#include <QApplication>

#include <orsaQt/debug.h>

#include <iostream>

// #include <tbb/task_scheduler_init.h>

int main(int argc, char ** argv) {
  
  // TBB
  // tbb::task_scheduler_init init;
  
  QApplication app(argc, argv);
  
  orsaQt::Debug::instance()->initTimer();
  //
  // orsa::Debug::instance()->initTimer();
  
  QWidget * mainWidget = new RendezvousWithVestaWidget;
  
  mainWidget->show();
  
  app.connect(&app,
	      SIGNAL(lastWindowClosed()), 
	      &app, 
	      SLOT(quit()));
  
  return app.exec();
}
