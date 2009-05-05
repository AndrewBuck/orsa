#ifndef __VESTA_VIZ_H__
#define __VESTA_VIZ_H__

// From the file osgviewerQT.cpp, version of October 10 2007

#include <orsa/debug.h>
#include <orsa/unit.h>

#include <osgViewer/Viewer>
#include <osgViewer/CompositeViewer>
#include <osgViewer/ViewerEventHandlers>
#include <osgGA/TrackballManipulator>
#include <osgDB/ReadFile>


#include <QtCore/QString>
#include <QtCore/QTimer>
#include <QtGui/QKeyEvent>
#include <QtGui/QApplication>
#include <QtOpenGL/QGLWidget>

using Qt::WindowFlags;

#include <iostream>

class AdapterWidget : public QGLWidget
{
    public:

        AdapterWidget( QWidget * parent = 0, const char * name = 0, const QGLWidget * shareWidget = 0, WindowFlags f = 0 );

        virtual ~AdapterWidget() {}

        osgViewer::GraphicsWindow* getGraphicsWindow() { return _gw.get(); }
        const osgViewer::GraphicsWindow* getGraphicsWindow() const { return _gw.get(); }

    protected:

        void init();

        virtual void resizeGL( int width, int height );
        virtual void keyPressEvent( QKeyEvent* event );
        virtual void keyReleaseEvent( QKeyEvent* event );
        virtual void mousePressEvent( QMouseEvent* event );
        virtual void mouseReleaseEvent( QMouseEvent* event );
        virtual void mouseMoveEvent( QMouseEvent* event );

        osg::ref_ptr<osgViewer::GraphicsWindowEmbedded> _gw;
};

AdapterWidget::AdapterWidget( QWidget * parent, const char *, const QGLWidget * shareWidget, WindowFlags f) :
  QGLWidget(parent, shareWidget, f) {
  _gw = new osgViewer::GraphicsWindowEmbedded(0,0,width(),height());
}

void AdapterWidget::resizeGL( int width, int height )
{
    _gw->getEventQueue()->windowResize(0, 0, width, height );
    _gw->resized(0,0,width,height);
}

void AdapterWidget::keyPressEvent( QKeyEvent* event ) {
  _gw->getEventQueue()->keyPress( (osgGA::GUIEventAdapter::KeySymbol) *(event->text().toAscii().data() ) );
}

void AdapterWidget::keyReleaseEvent( QKeyEvent* event ) {
  _gw->getEventQueue()->keyRelease( (osgGA::GUIEventAdapter::KeySymbol) *(event->text().toAscii().data() ) );
}

void AdapterWidget::mousePressEvent( QMouseEvent* event )
{
    int button = 0;
    switch(event->button())
    {
        case(Qt::LeftButton): button = 1; break;
        case(Qt::MidButton): button = 2; break;
        case(Qt::RightButton): button = 3; break;
        case(Qt::NoButton): button = 0; break;
        default: button = 0; break;
    }
    _gw->getEventQueue()->mouseButtonPress(event->x(), event->y(), button);
}

void AdapterWidget::mouseReleaseEvent( QMouseEvent* event )
{
    int button = 0;
    switch(event->button())
    {
        case(Qt::LeftButton): button = 1; break;
        case(Qt::MidButton): button = 2; break;
        case(Qt::RightButton): button = 3; break;
        case(Qt::NoButton): button = 0; break;
        default: button = 0; break;
    }
    _gw->getEventQueue()->mouseButtonRelease(event->x(), event->y(), button);
}

void AdapterWidget::mouseMoveEvent( QMouseEvent* event )
{
    _gw->getEventQueue()->mouseMotion(event->x(), event->y());
}


class ViewerQT : public osgViewer::Viewer, public AdapterWidget
{
    public:

        ViewerQT(QWidget * parent = 0, const char * name = 0, const QGLWidget * shareWidget = 0, WindowFlags f = 0):
            AdapterWidget( parent, name, shareWidget, f )
	  {
	    
	    // usually true is a good value
	    bool use_multiple_cameras = false;
	    
	    if (use_multiple_cameras) {
	      
	      const double FOV = 45.0;
	    
	    /* const double nearClip = orsa::FromUnits(  1,orsa::Unit::METER);
	       const double  farClip = orsa::FromUnits(1e3,orsa::Unit::AU);
	    */
	    //
	    const double nearClip = orsa::FromUnits(1e3,orsa::Unit::KM);
	    const double  farClip = orsa::FromUnits(1e2,orsa::Unit::AU);
	    //
	    const double subRatio = 1e-5;
	    
	    // the main camera (getCamera()) controls the view, but does not render,
	    // only the slave cameras render
	    //
	    getCamera()->setComputeNearFarMode(osg::CullSettings::DO_NOT_COMPUTE_NEAR_FAR);
	    getCamera()->setCullingMode(osg::CullSettings::VIEW_FRUSTUM_CULLING | 
					osg::CullSettings::SMALL_FEATURE_CULLING);
	    getCamera()->setProjectionMatrixAsPerspective(FOV, 
							  static_cast<double>(width())/static_cast<double>(height()), 
							  nearClip, 
							  farClip);
	    
	    setThreadingModel(osgViewer::Viewer::SingleThreaded);
	    
	    connect(&_timer, SIGNAL(timeout()), this, SLOT(updateGL()));
	    _timer.start(40);
	    
	    
	    const osg::Matrixd invertedMainProjection = osg::Matrixd::inverse(getCamera()->getProjectionMatrix());
	    
	    // from FAR to NEAR...
	    double clipDistance = farClip;
	    unsigned int cameraID = 0;
	    while (clipDistance > nearClip) {
	      
	      osg::ref_ptr<osg::Camera> camera = new osg::Camera;
	      camera->setGraphicsContext(getGraphicsWindow());
	      // camera->setViewport(new osg::Viewport(0,0,400,400));
	      camera->setViewport(new osg::Viewport(0,0,width(),height()));
	      if (cameraID == 0) {
		camera->setClearColor(osg::Vec4(1.0f,0.65f,0.0f,1.0f));
		camera->setClearMask(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	      } else {
		camera->setClearMask(GL_DEPTH_BUFFER_BIT);
	      }
	      //
	      osg::Matrixd cameraProjection;
	      cameraProjection.makePerspective(FOV, 
					       static_cast<double>(width())/static_cast<double>(height()), 
					       clipDistance*subRatio, 
					       clipDistance);
	      //
	      osg::Matrixd relativeProjection = invertedMainProjection;
	      relativeProjection.postMult(cameraProjection);
	      //
	      addSlave(camera.get(),
		       relativeProjection, 
		       osg::Matrixd());
	      
	      // step
	      clipDistance *= subRatio;
	      ++cameraID;
	    }
	    
	    // NOTE: the HUD should be a child of the latest camera, otherwise is rendered N times (N=num.cameras)
	    
	    ORSA_DEBUG("total cameras: %i",cameraID);
	    
	    } else {

	      // use_multiple_cameras == false
	      
	      getCamera()->setViewport(new osg::Viewport(0,0,width(),height()));
	      getCamera()->setProjectionMatrixAsPerspective(30.0f, static_cast<double>(width())/static_cast<double>(height()), 1.0f, 10000.0f);
	      getCamera()->setGraphicsContext(getGraphicsWindow());
	      
	      setThreadingModel(osgViewer::Viewer::SingleThreaded);
	      
	      connect(&_timer, SIGNAL(timeout()), this, SLOT(updateGL()));
	      _timer.start(40);
	      
	    }
	}
	
        virtual void paintGL()
        {
          
	    const bool debug = false;
	    
	    if (debug) {
	      double fovy;
	      double aspectRatio;
	      double zNear; 
	      double zFar;
	      //
	      // getSceneView()->getProjectionMatrixAsPerspective(fovy, aspectRatio,
	      // zNear, zFar);
	      //
	      getCamera()->getProjectionMatrixAsPerspective(fovy, aspectRatio,
							    zNear, zFar);
	      
	      ORSA_DEBUG("fovy........: %e",fovy);
	      ORSA_DEBUG("aspectRatio.: %e",aspectRatio);
	      ORSA_DEBUG("zNear.......: %e",zNear);
	      ORSA_DEBUG("zFar........: %e",zFar);
	      ORSA_DEBUG("NearFarRatio: %e",zNear/zFar);
	    }
	    
	    if (debug) {
	      osg::Vec3 eye;
	      osg::Vec3 center;
	      osg::Vec3 up;
	      //
	      // getSceneView()->getViewMatrixAsLookAt(eye, center, up);
	      //
	      getCamera()->getViewMatrixAsLookAt(eye, center, up);
	      
	      ORSA_DEBUG("eye...: %e %e %e",   eye[0],   eye[1],   eye[2]);
	      ORSA_DEBUG("center: %e %e %e",center[0],center[1],center[2]);
	      ORSA_DEBUG("up....: %e %e %e",    up[0],    up[1],    up[2]);
	      //
	      ORSA_DEBUG("center-eye: %e",
			 sqrt((center[0]-eye[0])*(center[0]-eye[0])+
			      (center[1]-eye[1])*(center[1]-eye[1])+
			      (center[2]-eye[2])*(center[2]-eye[2])));
	    }
	    
	    frame();
        }
    
    protected:

        QTimer _timer;
};

#endif // __VESTA_VIZ_H__
