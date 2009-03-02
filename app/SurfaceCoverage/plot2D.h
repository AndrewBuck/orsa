#ifndef __SURFACE_COVERAGE_PLOT_2D_H__
#define __SURFACE_COVERAGE_PLOT_2D_H__

#include <qpainter.h>

#include <qwt_painter.h>
#include <qwt_plot.h>
#include <qwt_plot_item.h>

#include "mainThread.h"

#include "FaceColor.h"

class TriangleItem: public QwtPlotItem {
  
 public:
  TriangleItem(MainThread * mt,
	       const unsigned int f) :
    QwtPlotItem(),
    face(f) { 
    
    mainThread = mt;
    
    const orsa::TriShape::TriIndex tri = mainThread->faceVector_out[face];
    
    d_pen = QPen(Qt::gray);
    
    lat_i = mainThread->latitude[face][tri.i()];
    lat_j = mainThread->latitude[face][tri.j()];
    lat_k = mainThread->latitude[face][tri.k()];
    //
    lon_i = mainThread->longitude[face][tri.i()];
    lon_j = mainThread->longitude[face][tri.j()];
    lon_k = mainThread->longitude[face][tri.k()];
    //
    d_rect = QwtDoubleRect(std::min(lon_i,std::min(lon_j,lon_k)),
			   std::min(lat_i,std::min(lat_j,lat_k)),
			   std::max(lon_i,std::max(lon_j,lon_k)),
			   std::max(lat_i,std::max(lat_j,lat_k)));
  }
    
 public:
  virtual QwtDoubleRect boundingRect() const {
    return d_rect;
  }
  
 protected:
  void setBrush() const {
    d_brush = FaceColor(mainThread,
			face,
			mainThread->orbitEpoch.getRef() + mainThread->runDuration.getRef());
  }
  
 public:
  void draw(QPainter * painter,
	    const QwtScaleMap & xMap, 
	    const QwtScaleMap & yMap,
	    const QRect &) const {
    
    setBrush();
    
    QVector<QPoint> point;
    
    point.push_back(QPoint(xMap.transform(lon_i),yMap.transform(lat_i)));
    point.push_back(QPoint(xMap.transform(lon_j),yMap.transform(lat_j)));
    point.push_back(QPoint(xMap.transform(lon_k),yMap.transform(lat_k)));
    
    painter->setPen(d_pen);
    painter->setBrush(d_brush);
    
    QwtPainter::drawPolygon(painter, QwtPolygon(point));
  }
  
 private:
  mutable QPen   d_pen;
  mutable QBrush d_brush;
  QwtDoubleRect  d_rect;
  
 private: 
  MainThread * mainThread;
  
 private:
  const unsigned int face;
  
 private:
  double lat_i, lon_i;
  double lat_j, lon_j;
  double lat_k, lon_k;
};

class Plot2D : public QwtPlot {
  
 public:
  Plot2D (QWidget * parent = 0);
  
 public:
  void setMainThread(MainThread * mt);
  
 private: 
  MainThread * mainThread;
};

#endif // __SURFACE_COVERAGE_PLOT_2D_H__
