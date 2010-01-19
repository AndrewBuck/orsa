#ifndef __SURFACE_COVERAGE_PLOT_TIMELINE_H__
#define __SURFACE_COVERAGE_PLOT_TIMELINE_H__

#include <qwt_plot.h>
#include <qwt_plot_curve.h>

#include "mainThread.h"

#include "FaceColor.h"

class TimelineCurve : public QwtPlotCurve {
 public:
  TimelineCurve(MainThread * mt,
		Qt::GlobalColor c) :
    QwtPlotCurve(),
    color(c) {
    mainThread = mt;
    setBrush(color);
    setRenderHint(QwtPlotItem::RenderAntialiased, true);
    totalArea=0;
    for (unsigned int f=0; f<mainThread->faceVector_out.size(); ++f) {
      totalArea += mainThread->faceArea_out[f];
      // ORSA_DEBUG("area of face %i: %f",f,mainThread->faceArea_out[f]);
    }
    if (!(totalArea > 0)) {
      ORSA_DEBUG("problems...");
    }
  }
 public:
  void updateData() {
    QVector<double> xData;
    QVector<double> yData;
    const orsa::Time tStop = mainThread->orbitEpoch.getRef() + mainThread->runDuration.getRef();
    const orsa::Time dt = mainThread->cameraInterval.getRef();
    orsa::Time t = mainThread->orbitEpoch.getRef();
    double tmpArea;
    while (t <= tStop) {
      tmpArea=0;
      for (unsigned int face=0; face<mainThread->faceVector_out.size(); ++face) {
       	const Qt::GlobalColor faceColor = FaceColor(mainThread, face, t);
	if (CumulativeColorMatch(faceColor,color)) {
	  tmpArea += mainThread->faceArea_out[face];
	}
      }
      //
      /* 
	 ORSA_DEBUG("coverage: %f   t: %f",
	 100 * tmpArea / totalArea,
	 FromUnits((t-mainThread->orbitEpoch.getRef()),
	 orsa::Unit::HOUR,-1));
      */
      //
      xData.push_back(FromUnits((t-mainThread->orbitEpoch.getRef()).get_d(),
				orsa::Unit::HOUR,-1));
      yData.push_back(100 * tmpArea / totalArea); // percent
      //
      t += dt;
    }
    setData(QwtPointArrayData(xData,
			      yData));
  }
 protected: 
  MainThread * mainThread;
 protected:
  const Qt::GlobalColor color;
 protected:
  double totalArea;
};

class PlotTimeline : public QwtPlot {
  
 public:
  PlotTimeline (QWidget * parent = 0);
  
 public:
  void setMainThread(MainThread * mt);
  
 private:	
  TimelineCurve * curve_yellow;
  TimelineCurve * curve_green;
  TimelineCurve * curve_darkGreen;
  
 private: 
  MainThread * mainThread;
};

#endif // __SURFACE_COVERAGE_PLOT_TIMELINE_H__
