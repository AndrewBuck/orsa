#include "plot.h"

#include <qwt_plot_canvas.h>

#include <QPainter>

// #include <orsa/attitude.h>
#include <orsa/debug.h>
#include <orsa/orbit.h>
#include <orsa/util.h>

PlotFillThread::PlotFillThread(orsa::BodyGroup  * bg,
			       const orsa::Time & dt) :
  QThread(),
  _bg(bg),
  _dt(dt) {
  data = new DataType;
  data->enableDataStoring();
}

void PlotFillThread::run() {
  
  dataMutex.lock();
  data->reset();
  dataMutex.unlock();
  
  doAbort = false;
  
  orsa::Time t_start, t_stop;
  _bg->getGlobalInterval(t_start,t_stop,false);
  
  osg::ref_ptr<const orsa::Body> vesta = _bg->getBody("VESTA");
  osg::ref_ptr<const orsa::Body> dawn  = _bg->getBody("DAWN");
  
  if (vesta.get() && dawn.get()) {
    
    osg::ref_ptr<const orsa::Shape> vesta_shape = vesta->getInitialConditions().inertial->localShape();
    
    /* 
       osg::ref_ptr<const orsa::Attitude> vesta_attitude = 
       new orsa::BodyAttitude(vesta.get(),
       _bg.get());
    */
    
    orsa::Vector rVesta, rDAWN;
    orsa::Vector intersectionPoint;
    orsa::Vector normal;
    
    PlotData tmpData;
    
    orsa::Time t = t_start;
    // strictly less than "t_stop+_dt"
    while (t < (t_stop+_dt)) {
      
      if (doAbort) {
	break;
      }
      
      // going out at the next step, when we'll have t==t_stop+_dt
      if (t > t_stop) {
	t = t_stop;
      }
      
      if (_bg->getInterpolatedPosition(rVesta,vesta.get(),t) &&
	  _bg->getInterpolatedPosition(rDAWN,  dawn.get(),t) ) {
	
	const orsa::Vector dr = rDAWN-rVesta;
	
	// xData[k] = FromUnits((t-t_start).get_d(),orsa::Unit::SECOND,-1);
	// tmpData.t = t-t_start;
	tmpData.t = t;
	
      	tmpData.distance = dr.length();
	
	if (vesta_shape.get()) {
	  
	  tmpData.sphere_radius = vesta_shape->boundingRadius();
	  
	  // ORSA_DEBUG("radius: %Ff",vesta_shape->boundingRadius().get_mpf_t());
	  
	  const orsa::Matrix g2l = orsa::globalToLocal(vesta.get(),_bg.get(),t);
	  
	  /* const orsa::Vector dr_local = (vesta_attitude.get()) ? 
	     (vesta_attitude->globalToLocal(t)*dr) : dr;
	  */
	  //
	  const orsa::Vector dr_local = g2l*dr;
	  
	  if (vesta_shape->rayIntersection(intersectionPoint,
					   normal,
					   dr_local,
					   (-dr_local).normalized(),
					   false)) {
	    tmpData.vesta_profile = intersectionPoint.length();
	  } else {
	    ORSA_DEBUG("problems...");
	    tmpData.vesta_profile = 0;
	  }
	  
	} else {
	  tmpData.sphere_radius = 0;
	  tmpData.vesta_profile = 0;
	}
	
	// q,a,Q
	{
	  orsa::Orbit orbit;
	  if (orbit.compute(dawn.get(),
			    vesta.get(),
			    _bg.get(),
			    t)) {
	    tmpData.q = orbit.a * (1-orbit.e);
	    tmpData.a = orbit.a;
	    tmpData.Q = orbit.a * (1+orbit.e);
	  } else {
	    tmpData.q = 0;
	    tmpData.a = 0;
	    tmpData.Q = 0;
	  }	      
	}
	
	dataMutex.lock();
	data->insert(tmpData,false,false);
	dataMutex.unlock();
	
      }
      
      t += _dt;
    }
    
  } else {
    
    ORSA_WARNING("bodies not found...");
    
  } 
  
}

void PlotFillThread::abort() {
  doAbort = true;
}

////

VestaPlot::~VestaPlot() {
  // ORSA_DEBUG("called ~VestaPlot()...");
  
  delete curve_distance;
  curve_distance = 0;
  
  delete curve_sphere_radius;
  curve_sphere_radius = 0;
  
  delete curve_vesta_profile;
  curve_vesta_profile = 0;
  
  delete curve_q;
  curve_q = 0;
  
  delete curve_a;
  curve_a = 0;
  
  delete curve_Q;
  curve_Q = 0;
  
  delete curve_altitude;
  curve_altitude = 0;
  
  delete marker_altitude;
  marker_altitude = 0;

  delete plotFillThread;
  plotFillThread = 0;
  
  if (_canvasCache) {
    delete _canvasCache;
    _canvasCache = 0;
  }
}

VestaPlot::VestaPlot(orsa::BodyGroup        * bg,
		     orsaOSG::AnimationTime * at,
		     QWidget                * parent) :
  QwtPlot(parent),
  _bg(bg),
  _at(at),
  _canvasCache(0),
  _canvasCacheDirty(true) {
  
  setAttribute(Qt::WA_DeleteOnClose,true);
  
  setAutoReplot(false);
  
  setMargin(3);
  
  setAxisTitle(xBottom,"time [s]");
  setAxisTitle(yLeft,  "[km]");
  
  curve_distance = new VestaPlotCurve(true);
  // curve_distance->setRenderHint(QwtPlotItem::RenderAntialiased);
  curve_distance->setYAxis(QwtPlot::yLeft);
  curve_distance->attach(this);
  //
  {
    QPen pen;
    pen.setWidth(2);
    pen.setStyle(Qt::SolidLine);
    curve_distance->setPen(pen);
  }  
  
  curve_sphere_radius = new VestaPlotCurve(true);
  // curve_sphere_radius->setRenderHint(QwtPlotItem::RenderAntialiased);
  curve_sphere_radius->setYAxis(QwtPlot::yLeft);
  curve_sphere_radius->attach(this);
  //
  {
    QPen pen;
    pen.setStyle(Qt::DashDotLine);
    curve_sphere_radius->setPen(pen);
  }  
  
  curve_vesta_profile = new VestaPlotCurve(true);
  // curve_vesta_profile->setRenderHint(QwtPlotItem::RenderAntialiased);
  curve_vesta_profile->setYAxis(QwtPlot::yLeft);
  curve_vesta_profile->attach(this);
  //
  {
    QColor c(Qt::gray);
    c.setAlpha(150);
    
    // curve_vesta_profile->setPen(c);
    curve_vesta_profile->setBrush(c);
  }
  
  curve_q = new VestaPlotCurve(true);
  // curve_q->setRenderHint(QwtPlotItem::RenderAntialiased);
  curve_q->setYAxis(QwtPlot::yLeft);
  curve_q->attach(this);
  //
  {
    QPen pen;
    pen.setStyle(Qt::DotLine);
    curve_q->setPen(pen);
  }  
  
  curve_a = new VestaPlotCurve(true);
  // curve_a->setRenderHint(QwtPlotItem::RenderAntialiased);
  curve_a->setYAxis(QwtPlot::yLeft);
  curve_a->attach(this);
  //
  {
    QPen pen;
    pen.setStyle(Qt::DashLine);
    curve_a->setPen(pen);
  }  
  
  curve_Q = new VestaPlotCurve(true);
  // curve_Q->setRenderHint(QwtPlotItem::RenderAntialiased);
  curve_Q->setYAxis(QwtPlot::yLeft);
  curve_Q->attach(this);
  //
  {
    QPen pen;
    pen.setStyle(Qt::DotLine);
    curve_Q->setPen(pen);
  }  
  
  curve_altitude = new VestaPlotCurve(false);
  // curve_altitude->setRenderHint(QwtPlotItem::RenderAntialiased);
  curve_altitude->setYAxis(QwtPlot::yLeft);
  curve_altitude->attach(this);
  //
  {
    QPen pen;
    pen.setStyle(Qt::DashLine);
    curve_altitude->setPen(pen);
  }  
  
  marker_altitude = new VestaPlotMarker(false);
  marker_altitude->setLineStyle(QwtPlotMarker::NoLine);
  // marker_altitude->setLineStyle(QwtPlotMarker::VLine);
  marker_altitude->setLabelAlignment(Qt::AlignRight | Qt::AlignVCenter);
  marker_altitude->attach(this);
  
  {
    const orsa::Time dt(0,0,15,0,0);
    plotFillThread = new PlotFillThread(_bg.get(),
					dt);
    //
    // this returns instantly!
    plotFillThread->start(QThread::LowestPriority);
  }
  
  plotFillTimer.start(500);
  
  connect(&plotFillTimer,
	  SIGNAL(timeout()),
	  this,
	  SLOT(plotFill()));
  
  connect(plotFillThread,
	  SIGNAL(finished()),
	  this,
	  SLOT(plotFill()));
  
  connect(plotFillThread,
	  SIGNAL(finished()),
	  &plotFillTimer,
	  SLOT(stop()));
  
  /* 
     connect(_at.get(),
     SIGNAL(simulationTimeChanged(const orsa::Time &)),
     this,
     SLOT(moveMarker(const orsa::Time &)));
  */
  
  connect(_at.get(),
	  SIGNAL(simulationTimeChanged(const orsa::Time &)),
	  this,
	  SLOT(moveAltitudeCurve(const orsa::Time &)));
  
  resize(512,256);
  
  // one call, to fill a bit...
  plotFill();
  
  moveTime.start();
}

/* 
   void VestaPlot::moveMarker(const orsa::Time & t) const {
   
   // ORSA_DEBUG("inside VestaPlot::moveMaker()...");
   
   if (_lastSimulationTime.isSet() &&
   _lastSimulationTime.getRef() == t) {
   // nothing to do here...
   return; 
   }
   
   orsa::Time t_start, t_stop;
   _bg->getGlobalInterval(t_start,t_stop,false);
   marker->setValue(FromUnits((t-t_start).get_d(),orsa::Unit::SECOND,-1), 
   0.0);
   
   _lastSimulationTime = t;
   }
*/

void VestaPlot::moveAltitudeCurve(const orsa::Time & t) {
  
  // don't update too often...
  if (moveTime.elapsed() < 40) return;
  
  // ORSA_DEBUG("inside VestaPlot::moveAltitudeCurve()...");
  
  if (_lastSimulationTime.isSet() &&
      _lastSimulationTime.getRef() == t) {
    // nothing to do here...
    return; 
  }
  
  orsa::Time t_start, t_stop;
  _bg->getGlobalInterval(t_start,t_stop,false);
  
  double xData[2];
  double yData[2];
  
  xData[0] = xData[1] = FromUnits((t-t_start).get_d(),orsa::Unit::SECOND,-1);
  
  plotFillThread->dataMutex.lock();
  //
  {
    const unsigned int vSize = plotFillThread->data->size();
    
    // ORSA_DEBUG("vSize: %i",vSize);
    
    if (vSize >= 2) {
      
      /* 
	 ORSA_DEBUG("t: %Ff   min().t: %Ff   max().t: %Ff   size: %i",
	 t.get_d().get_mpf_t(),
	 plotFillThread->data->min().t.get_d().get_mpf_t(),
	 plotFillThread->data->max().t.get_d().get_mpf_t(),
	 plotFillThread->data->size());
      */
      
      if (t < plotFillThread->data->min().t) {
	yData[0] = FromUnits(plotFillThread->data->min().vesta_profile,
			     orsa::Unit::KM,-1);
	yData[1] = FromUnits(plotFillThread->data->min().distance,
			     orsa::Unit::KM,-1);
      } else if (t > plotFillThread->data->max().t) {
       	yData[0] = FromUnits(plotFillThread->data->max().vesta_profile,
			     orsa::Unit::KM,-1);
	yData[1] = FromUnits(plotFillThread->data->max().distance,
			     orsa::Unit::KM,-1);
      } else {
	
      	PlotData pd, pdMin, pdMax;
	pd.t = t;
	//
	const bool sub = plotFillThread->data->getSubInterval(pd, pdMin, pdMax);
	//
	if (!sub) {
	  
	  yData[0] = yData[1] = 0.0;
	  
	} else if (pdMin.t != pdMax.t) {
	  
	  /* 
	     ORSA_DEBUG("t: %Ff   pdMin.t: %Ff   pdMax.t: %Ff   size: %i",
	     t.get_d().get_mpf_t(),
	     pdMin.t.get_d().get_mpf_t(),
	     pdMax.t.get_d().get_mpf_t(),
	     plotFillThread->data->size());
	  */
	  
	  const double yDistanceLeft  = FromUnits(pdMin.distance,
						  orsa::Unit::KM,-1);
	  const double yDistanceRight = FromUnits(pdMax.distance,
						  orsa::Unit::KM,-1);
	  
	  const double yProfileLeft   = FromUnits(pdMin.vesta_profile,
						  orsa::Unit::KM,-1);
	  const double yProfileRight  = FromUnits(pdMax.vesta_profile,
						  orsa::Unit::KM,-1);
	  
	  const double left_t_k = 
	    FromUnits((pdMin.t-t_start).get_d(),
		      orsa::Unit::SECOND,-1);
	  const double right_t_k = 
	    FromUnits((pdMax.t-t_start).get_d(),
		      orsa::Unit::SECOND,-1);
	  
	  const double beta = (xData[0]-left_t_k)/(right_t_k-left_t_k);
	  
	  // ORSA_DEBUG("beta: %f",beta);
	  
	  yData[0] = 
	    (1.0-beta)*yProfileLeft + 
	    beta*yProfileRight;
	  
	  yData[1] = 
	    (1.0-beta)*yDistanceLeft + 
	    beta*yDistanceRight;
	  
	} else {
	  yData[0] = FromUnits(pdMin.vesta_profile,
			       orsa::Unit::KM,-1);
	  yData[1] = FromUnits(pdMin.distance,
			       orsa::Unit::KM,-1);
	}
	
      }
      
    } else if (vSize == 1) {
      yData[0] = FromUnits(plotFillThread->data->min().vesta_profile,
			   orsa::Unit::KM,-1);
      yData[1] = FromUnits(plotFillThread->data->min().distance,
			   orsa::Unit::KM,-1);
    } else {
      // vSize == 0;
      yData[0] = yData[1] = 0.0;
    }
  }
  //
  plotFillThread->dataMutex.unlock();
  
  curve_altitude->setData(xData,
			  yData,
			  2);
  // this below will work with Qwt 5.3
  /* curve_altitude->setData(QwtPointArrayData(xData,
     yData,
     2));
  */
  
  {
    // altitude marker
    
    marker_altitude->setValue(xData[0], 
			      0.5*(yData[0]+yData[1]));
    
    QString label;
    label.sprintf("%.1f [km]",yData[1]-yData[0]);
    
    QwtText text(label);
    text.setFont(QFont("Helvetica", 10, QFont::Bold));
    
    marker_altitude->setLabel(text); 
  }
  
  _lastSimulationTime = t;
  
  moveTime.restart();
  
  replot();
}

void VestaPlot::plotFill() {
  
  // ORSA_DEBUG("inside VestaPlot::plotFill()...");
  
  orsa::Time t_start, t_stop;
  _bg->getGlobalInterval(t_start,t_stop,false);
  
  plotFillThread->dataMutex.lock();
  
  // const unsigned int vSize = plotFillThread->data->size();
  
  QVector<double> xData;
  //
  QVector<double> yData_distance;
  QVector<double> yData_sphere_radius;
  QVector<double> yData_vesta_profile;
  //
  QVector<double> yData_q;
  QVector<double> yData_a;
  QVector<double> yData_Q;
  
  /* 
     xData.resize(vSize);
     //
     yData_distance.resize(vSize);
     yData_sphere_radius.resize(vSize);
     yData_vesta_profile.resize(vSize);
     //
     yData_q.resize(vSize);
     yData_a.resize(vSize);
     yData_Q.resize(vSize);
  */
  
  /* 
     for (unsigned int k=0; k<vSize; ++k) {
     xData[k]               = FromUnits(((plotFillThread->data->getData())[k].t-t_start).get_d(),
     orsa::Unit::SECOND,-1);
     
     yData_distance[k]      = FromUnits((plotFillThread->data->getData())[k].distance,
     orsa::Unit::KM,-1);
     yData_sphere_radius[k] = FromUnits((plotFillThread->data->getData())[k].sphere_radius,
     orsa::Unit::KM,-1);
     yData_vesta_profile[k] = FromUnits((plotFillThread->data->getData())[k].vesta_profile,
     orsa::Unit::KM,-1);
     
     yData_q[k] = FromUnits((plotFillThread->data->getData())[k].q,
     orsa::Unit::KM,-1);
     yData_a[k] = FromUnits((plotFillThread->data->getData())[k].a,
     orsa::Unit::KM,-1);
     yData_Q[k] = FromUnits((plotFillThread->data->getData())[k].Q,
     orsa::Unit::KM,-1);
     }
  */
  
  PlotFillThread::DataType::DataType::const_iterator it = plotFillThread->data->getData().begin();
  while (it != plotFillThread->data->getData().end()) {
    const PlotData & pd = (*it);
     
    xData.push_back(FromUnits((pd.t-t_start).get_d(),
			      orsa::Unit::SECOND,-1));
    
    yData_distance.push_back(FromUnits(pd.distance,
				       orsa::Unit::KM,-1));
    yData_sphere_radius.push_back(FromUnits(pd.sphere_radius,
					    orsa::Unit::KM,-1));
    yData_vesta_profile.push_back(FromUnits(pd.vesta_profile,
					    orsa::Unit::KM,-1));
    
    yData_q.push_back(FromUnits(pd.q,
				orsa::Unit::KM,-1));
    yData_a.push_back(FromUnits(pd.a,
				orsa::Unit::KM,-1));
    yData_Q.push_back(FromUnits(pd.Q,
				orsa::Unit::KM,-1));
    ++it;
  }
  
  curve_distance->setData(xData,
			  yData_distance);
  curve_sphere_radius->setData(xData,
			       yData_sphere_radius);
  curve_vesta_profile->setData(xData,
			       yData_vesta_profile);
  
  curve_q->setData(xData,
		   yData_q);
  curve_a->setData(xData,
		   yData_a);
  curve_Q->setData(xData,
		   yData_Q);
  
  // this below will work with Qwt 5.3
  /* curve_distance->setData(QwtPointArrayData(xData,
     yData_distance));
     curve_sphere_radius->setData(QwtPointArrayData(xData,
     yData_sphere_radius));
     curve_vesta_profile->setData(QwtPointArrayData(xData,
     yData_vesta_profile));
     
     curve_q->setData(QwtPointArrayData(xData,
     yData_q));
     curve_a->setData(QwtPointArrayData(xData,
     yData_a));
     curve_Q->setData(QwtPointArrayData(xData,
     yData_Q));
  */
  
  plotFillThread->dataMutex.unlock();
  
  _canvasCacheDirty = true;
  replot();
}

void VestaPlot::closeEvent(QCloseEvent *) {
  if (plotFillThread) {
    plotFillThread->abort();
  }
}

void VestaPlot::drawItems(QPainter          * painter, 
			  const QRect       & rect,
			  const QwtScaleMap   maps[axisCnt],
			  const QwtPlotPrintFilter & filter) const {
  
  // this below will work with Qwt 5.3
  // const QRect rect = rectf.toRect();
  
  if (autoReplot()) {
    ORSA_ERROR("warning: this method does not work properly when autoReplot is enabled");
    return QwtPlot::drawItems(painter, 
			      rect, 
			      maps,
			      filter);
  }
  
  /* 
     if (_canvasCache) {
     // debug info
     std::cerr << " rect: " <<          rect.size().width() << "x" <<          rect.size().height() << " at (" <<          rect.x() << "," <<          rect.y() << ")" << std::endl;
     std::cerr << "cache: " << _canvasCache->size().width() << "x" << _canvasCache->size().height() << std::endl;
     } 
  */
  
  if (!_canvasCache) {
    _canvasCache      = new QPixmap(rect.size());
    _canvasCacheDirty = true;
  }
  
  const QwtPlotItemList & itmList = itemList();
  if ( (!_canvasCacheDirty) && 
       (rect.size() == _canvasCache->size()) ) {
    // using cache...
    // ORSA_DEBUG("using cache...");
    painter->drawPixmap(rect,(*_canvasCache));
    // painter->drawPixmap(rect,(*_canvasCache),rect); // this call offsets the moving items by 2 pixels, and the rect (x,y) = (2,2) so maybe the solution to this problem lies in modifying the QRect appropriately
  } else {
    // not using cache...
    // ORSA_DEBUG("not using cache...");
    for (QwtPlotItemIterator it = itmList.begin(); it != itmList.end(); ++it) {
      {
	VestaPlotCurve * curve = dynamic_cast< VestaPlotCurve * > (*it);
	if (curve) {
	  if (!curve->_staticData) {
	    curve->setVisible(false);
	  }	
	}
      }
      //
      {
	VestaPlotMarker * marker = dynamic_cast< VestaPlotMarker * > (*it);
	if (marker) {
	  if (!marker->_staticData) {
	    marker->setVisible(false);
	  }	
	}
      }
    }
    //
    QwtPlot::drawItems(painter, 
		       rect, 
		       maps,
		       filter);
    //
    for (QwtPlotItemIterator it = itmList.begin(); it != itmList.end(); ++it) {
      {
	VestaPlotCurve * curve = dynamic_cast< VestaPlotCurve * > (*it);
	if (curve) {
	  if (!curve->_staticData) {
	    curve->setVisible(true);
	  }	
	}
      }
      //
      {
	VestaPlotMarker * marker = dynamic_cast< VestaPlotMarker * > (*it);
	if (marker) {
	  if (!marker->_staticData) {
	    marker->setVisible(true);
	  }	
	}
      }
    }
    //
    (*_canvasCache)   = canvas()->paintCache()->copy();
    _canvasCacheDirty = false;
  }
  
  for (QwtPlotItemIterator it = itmList.begin(); it != itmList.end(); ++it) {
    {
      VestaPlotCurve * curve = dynamic_cast< VestaPlotCurve * > (*it);
      if (curve) {
	if (!curve->_staticData) {
	  curve->draw(painter,
		      maps[curve->xAxis()],
		      maps[curve->yAxis()],
		      rect);
	}	
      }
    }
    //
    {
      VestaPlotMarker * marker = dynamic_cast< VestaPlotMarker * > (*it);
      if (marker) {
	if (!marker->_staticData) {
	  marker->draw(painter,
		       maps[marker->xAxis()],
		       maps[marker->yAxis()],
		       rect);
	}	
      }
    }
  }
}
