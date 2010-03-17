#ifndef __VESTA_PLOT_H__
#define __VESTA_PLOT_H__

#include <qwt_plot.h>
#include <qwt_plot_curve.h>
#include <qwt_plot_marker.h>

#include <orsa/bodygroup.h>
#include <orsa/interval.h>
#include <orsa/unit.h>

#include <orsaOSG/AnimationTime.h>

#include <QList>
#include <QMutex>
#include <QThread>
#include <QTime>
#include <QTimer>

class PlotData {
 public:
  orsa::Time   t;
  double distance, sphere_radius, vesta_profile;
  double q, a, Q;
 public:
  inline bool operator == (const PlotData & rhs) const {
    return (t == rhs.t);
  }
 public:
  inline bool operator != (const PlotData & rhs) const {
    return (t != rhs.t);
  }
 public:
  inline bool operator < (const PlotData & rhs) const {
    return (t < rhs.t);
  }
 public:
  inline bool operator > (const PlotData & rhs) const {
    return (t > rhs.t);
  }
 public:
  inline bool operator <= (const PlotData & rhs) const {
    return (t <= rhs.t);
  }
 public:
  inline bool operator >= (const PlotData & rhs) const {
    return (t >= rhs.t);
  }
};

class PlotFillThread : public QThread {
  
  Q_OBJECT;
  
 public:
  PlotFillThread(orsa::BodyGroup  * bg,
		 const orsa::Time & dt);
 public:
  virtual ~PlotFillThread() { }
  
 protected:
  void run();
  
 public slots:
  void abort();
 private:
  bool doAbort;
  
 protected:    
  osg::ref_ptr<orsa::BodyGroup> _bg;
  const orsa::Time              _dt;
  
 public:
  // std::vector<PlotData> data;
  // QList<PlotData> data;
  typedef orsa::Interval<PlotData> DataType;
  osg::ref_ptr<DataType> data;
  QMutex   dataMutex;
};

class VestaPlotCurve : public QwtPlotCurve {
 public:
  VestaPlotCurve(bool staticData) : 
    QwtPlotCurve(),
    _staticData(staticData) { }
 public:
  const bool _staticData;
};

class VestaPlotMarker : public QwtPlotMarker {
 public:
  VestaPlotMarker(bool staticData) : 
    QwtPlotMarker(),
    _staticData(staticData) { }
 public:
  const bool _staticData;
};

class VestaPlot : public QwtPlot {
  
  Q_OBJECT;
  
 public:
  VestaPlot(orsa::BodyGroup        * bg,
	    orsaOSG::AnimationTime * at,
	    QWidget                * parent = 0);
  
 public:
  ~VestaPlot();
  
 protected slots:
  void plotFill();
  
 protected:
  VestaPlotCurve * curve_distance; 
  VestaPlotCurve * curve_sphere_radius; 
  VestaPlotCurve * curve_vesta_profile; 
  
 protected:
  VestaPlotCurve * curve_q; 
  VestaPlotCurve * curve_a; 
  VestaPlotCurve * curve_Q; 
  
 protected:
  VestaPlotCurve  *  curve_altitude; 
  VestaPlotMarker * marker_altitude;
  
 protected slots:
  void moveAltitudeCurve(const orsa::Time & t);
 protected:
  mutable QTime moveTime;
  
 protected:
  PlotFillThread * plotFillThread;
  mutable QTimer   plotFillTimer;
  void drawItems(QPainter *, 
		 const QRect &,
		 const QwtScaleMap maps[axisCnt],
		 const QwtPlotPrintFilter &) const;
  // this below will work with Qwt 5.3
  /* void drawItems(QPainter *, 
     const QRectF &,
     const QwtScaleMap maps[axisCnt]) const;
  */
 protected:
  void closeEvent(QCloseEvent *);  
  
 protected:    
  osg::ref_ptr<orsa::BodyGroup> _bg;
 protected:
  osg::ref_ptr<orsaOSG::AnimationTime> _at;
 protected:
  mutable orsa::Cache<orsa::Time> _lastSimulationTime;
 protected:
  mutable QPixmap * _canvasCache;
  mutable bool      _canvasCacheDirty;
};

#endif // __VESTA_PLOT_H__
