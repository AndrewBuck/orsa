#include "plot2D.h"

Plot2D::Plot2D (QWidget * parent) : QwtPlot(parent) {
  
  mainThread = 0;
  
  setCanvasBackground(QColor(Qt::black));
  
  setAxisScale(xBottom, -180, 180, 30);
  setAxisScale(yLeft,    -90,  90, 30);
  
  setAxisTitle(xBottom, "longitude");
  setAxisTitle(yLeft,   "latitude");
  
}

void Plot2D::setMainThread(MainThread * mt) {
  
  if (!mt) return;
  
  const bool populate = (mainThread == 0);
  
  mainThread = mt;
  
  // clear();
  
  if (populate) {
    for (unsigned int f=0; f<mainThread->faceVector_out.size(); ++f) {
      TriangleItem * item = new TriangleItem(mainThread,f);
      item->attach(this);
    }
  }
  
  replot();
}

