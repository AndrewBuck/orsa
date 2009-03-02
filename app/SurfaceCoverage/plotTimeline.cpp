#include "plotTimeline.h"

PlotTimeline::PlotTimeline (QWidget * parent) : QwtPlot(parent) {
  
  mainThread = 0;
  
  setCanvasBackground(QColor(Qt::black));
  
  setAxisScale(yLeft, 0, 100, 20);
  
  setAxisTitle(xBottom, "time [h]");
  setAxisTitle(yLeft,   "coverage [%]");
  
}

void PlotTimeline::setMainThread(MainThread * mt) {
  
  if (!mt) return;
  
  const bool populate = (mainThread == 0);
  
  mainThread = mt;
  
  /* 
     setAxisScale(xBottom, 
     0, 
     mpz_class(mainThread->runDuration.getRef().getMuSec() / mpz_class("3600000000")).get_ui());
  */
  
  if (populate) {
    curve_yellow    = new TimelineCurve(mainThread,Qt::yellow);
    curve_yellow->attach(this);
    //
    curve_green     = new TimelineCurve(mainThread,Qt::green);
    curve_green->attach(this);
    //
    curve_darkGreen = new TimelineCurve(mainThread,Qt::darkGreen);
    curve_darkGreen->attach(this);
  }
  
  curve_yellow->updateData();
  curve_green->updateData();
  curve_darkGreen->updateData();
  
  replot();
}

