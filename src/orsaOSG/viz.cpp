#include <orsaOSG/viz.h>

#include <osg/BlendFunc>
#include <osg/Depth>
#include <osg/Geode>
#include <osg/Geometry>
#include <osg/Light>
#include <osg/LightModel>
#include <osg/LightSource>
#include <osg/LineWidth>
#include <osg/Material>
#include <osg/MatrixTransform>
#include <osg/Node>
#include <osg/PositionAttitudeTransform>
#include <osg/Point>
#include <osg/PolygonMode>
#include <osg/ShapeDrawable>
#include <osg/StateAttribute>
#include <osg/Switch>

#include <osgText/Font>
#include <osgText/Text>

#include <orsaOSG/bioCallback.h>
#include <orsaOSG/FindNamedNodeVisitor.h>
#include <orsaOSG/BodyTranslationCallback.h>
#include <orsaOSG/CenterOfMassPositionCallback.h>
#include <orsaOSG/BodyAttitudeCallback.h>
#include <orsaOSG/LightSourceUpdateCallback.h>
#include <orsaOSG/OrbitNodeCallback.h>
#include <orsaOSG/AnimationTime.h>

#include <orsaOSG/HUD.h>

#include <gsl/gsl_rng.h>

using namespace osg;
using namespace orsa;
using namespace orsaOSG;

// by using this callback, we prevent the "erroneous" oscillation of asymmetric Bodies rotating quickly
class BodySymmetricBoundingBoxCallback : public Drawable::ComputeBoundingBoxCallback {
public:
    BodySymmetricBoundingBoxCallback(const orsa::Body * b) :
        ComputeBoundingBoxCallback(),
        _b(b) { }	
public:
    BoundingBox computeBound(const osg::Drawable &) const  { 
        if (_b.get()) {
            if (_b->getInitialConditions().inertial->localShape()) {
                const orsa::Box sb = _b->getInitialConditions().inertial->localShape()->symmetricBoundingBox();
                // extend using center-of-mass components (could simply offset, but then would not be symmetric anymore)
                // no rotations here...
                const orsa::Vector cm  = _b->getInitialConditions().inertial->centerOfMass();
                const orsa::Box lb = 
                    orsa::Box(sb.getXMin()-fabs(cm.getX()),
                              sb.getXMax()+fabs(cm.getX()),
                              sb.getYMin()-fabs(cm.getY()),
                              sb.getYMax()+fabs(cm.getY()),
                              sb.getZMin()-fabs(cm.getZ()),
                              sb.getZMax()+fabs(cm.getZ()));		    
                return lb.getOSGBoundingBox();
            }
        }
        return osg::BoundingBox();
    }
protected:
    osg::ref_ptr<const orsa::Body> _b;
};

// test
/* 
   class OrbitComputeBoundingSphereCallback : public osg::Node::ComputeBoundingSphereCallback {
   BoundingSphere computeBound(const osg::Node &) const { 
   ORSA_DEBUG("returning zero size BoundingSphere...");
   return BoundingSphere(osg::Vec3(0,0,0),0); 
   }
   };
*/

Viz::Viz(orsa::BodyGroup * bg, 
         const double & timeMultiplier) :
    osg::Referenced(true),
    _bg(bg),
    _timeMultiplier(timeMultiplier) { }

osg::Group * Viz::createRoot() {
  
    // time
    _at = new AnimationTime(_bg.get());
    _at->setTimeMultiplier(_timeMultiplier);
  
    // root
    osg::Group * rootNode = new osg::Group;
    //
    rootNode->setName("root");
  
    /* osg::ClearNode * clearNode = new osg::ClearNode;
       clearNode->setClearColor(osg::Vec4(1.0f,0.65f,0.0f,1.0f));
       rootNode->addChild(clearNode);
    */
  
    osg::Group * lastLightSource = rootNode;
  
    // light sources
    {
        BodyGroup::BodyList::const_iterator _b_it = _bg->getBodyList().begin();
        while (_b_it != _bg->getBodyList().end()) {
            if ((*_b_it)->isLightSource.getRef()) {

                /* 
                   ORSA_DEBUG("body [%s] is a light source...",
                   (*_b_it)->getName().c_str());
                */
	
                osg::LightSource * sunLightSource = new osg::LightSource;
	
                /* 
                   osg::Light * sunLight2 = sunLightSource->getLight();
                   // sunLight2->setPosition( osg::Vec4d(1e10, 0.0, 0.0, 1.0));
                   // sunLight2->setPosition( osg::Vec4d(0.0, FromUnits(0.0,Unit::KM), 0.0, 1.0));
                   // sunLight2->setDirection( osg::Vec3d(0.0, -1.0, 0.0));
                   // sunLight2->setConstantAttenuation(0.3);
                   //
                   sunLight2->setAmbient( osg::Vec4( 0.f, 0.f, 0.f, 1.0f ) );
	   
                   sunLightSource->setLight( sunLight2 );
                   sunLightSource->setLocalStateSetModes( osg::StateAttribute::ON );
                   sunLightSource->getOrCreateStateSet()->setMode(GL_LIGHTING, osg::StateAttribute::ON);
                */
	
                /* 
                   osg::LightModel * lightModel = new osg::LightModel;
                   // lightModel->setAmbientIntensity(osg::Vec4(0.0f,0.0f,0.0f,1.0f));
                   // lightModel->setAmbientIntensity(osg::Vec4(1.0f,1.0f,1.0f,1.0f));
                   lightModel->setAmbientIntensity(osg::Vec4(0.,0.,0.,1.0f));
                   sunLightSource->getOrCreateStateSet()->setAttribute(lightModel);
                */
	
                sunLightSource->setUpdateCallback(new LightSourceUpdateCallback(_bg.get(),
                                                                                (*_b_it).get(),
                                                                                _at.get()));
	
                lastLightSource->addChild(sunLightSource);
	
                lastLightSource = sunLightSource;
	
            } 
      
            ++_b_it;
        }
    }
  
    // big global 3D cage
    if (0) {
    
        osg::Geode * cageGeode = new osg::Geode();
    
        cageGeode->setName("cage");
      
        osg::Geometry * cageGeometry = new osg::Geometry();
      
        // const double step = FromUnits(16.0,Unit::KM);
        const double step = FromUnits(1.0,Unit::AU);
        //
        const          int num_lines_1d_one_side = abs(5); // must be positive
        const unsigned int num_lines_1d          = 1+2*abs(num_lines_1d_one_side);
        const unsigned int num_lines             = 3*num_lines_1d*num_lines_1d;
        // 
        const double edge = num_lines_1d_one_side*step;
    
        osg::Vec3Array * vertices = new osg::Vec3Array(2*num_lines);
        unsigned int _vertices_index = 0;
      
        // X 
        for (int j = -num_lines_1d_one_side; j <= num_lines_1d_one_side; ++j) {
            for (int k = -num_lines_1d_one_side; k <= num_lines_1d_one_side; ++k) {
                (*vertices)[_vertices_index].set(-edge,j*step,k*step);
                ++_vertices_index;
                (*vertices)[_vertices_index].set( edge,j*step,k*step);
                ++_vertices_index;
            }
        }
      
        // Y 
        for (int i = -num_lines_1d_one_side; i <= num_lines_1d_one_side; ++i) {
            for (int k = -num_lines_1d_one_side; k <= num_lines_1d_one_side; ++k) {
                (*vertices)[_vertices_index].set(i*step,-edge,k*step);
                ++_vertices_index;
                (*vertices)[_vertices_index].set(i*step, edge,k*step);
                ++_vertices_index;
            }
        }
      
        // Z
        for (int i = -num_lines_1d_one_side; i <= num_lines_1d_one_side; ++i) {
            for (int j = -num_lines_1d_one_side; j <= num_lines_1d_one_side; ++j) {
                (*vertices)[_vertices_index].set(i*step,j*step,-edge);
                ++_vertices_index;
                (*vertices)[_vertices_index].set(i*step,j*step, edge);
                ++_vertices_index;
            }
        }
      
        // check
        if (_vertices_index != (2*num_lines)) {
            ORSA_ERROR("index problem...");
        }
      
        cageGeometry->setVertexArray(vertices);
             
        osg::Vec4Array * colors = new osg::Vec4Array;
        colors->push_back(osg::Vec4(0.0f,1.0f,0.0f,1.0f));
        cageGeometry->setColorArray(colors);
        cageGeometry->setColorBinding(osg::Geometry::BIND_OVERALL);
      
        cageGeometry->getOrCreateStateSet()->setMode(GL_LIGHTING,osg::StateAttribute::OFF);
      
        cageGeometry->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINES,0,2*num_lines));
      
        cageGeode->addDrawable(cageGeometry);
    
        lastLightSource->addChild(cageGeode);
    }
  
    // for the orbit...
    //
    // unit offset circle
    //
    osg::Geometry * unitOffsetCircle = new osg::Geometry();
    //
    {
        osg::Vec4Array * colours = new osg::Vec4Array(1);
        // (*colours)[0] = osg::Vec4d(1.0,1.0,1.0,1.0);
        (*colours)[0] = osg::Vec4d(1.0,1.0,1.0,0.3);
        // (*colours)[0] = osg::Vec4d(0.0,1.0,0.0,0.3);
        unitOffsetCircle->setColorArray(colours);
        unitOffsetCircle->setColorBinding(osg::Geometry::BIND_OVERALL);
    
        const unsigned int n_points = 1024;
        //
        osg::Vec3Array * coords = new osg::Vec3Array(n_points);
        const double dx = twopi()/n_points;
        double s,c;
        for (unsigned int j=0; j<n_points; ++j) {
            orsa::sincos(dx*j,&s,&c);
            (*coords)[j].set(osg::Vec3d(c-1.0,s,0.0));
            // sincos(dx*(j+1),&s,&c);
            // (*coords)[j].set(osg::Vec3d(c,s,0.0));
        }
        //
        unitOffsetCircle->setVertexArray(coords);
        //
        unitOffsetCircle->getOrCreateStateSet()->setMode(GL_LIGHTING,osg::StateAttribute::OFF);
        unitOffsetCircle->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP,0,n_points));
    }
    //
    /* const unsigned int numSegments = 64;
       osg::Geometry * unitOffsetCircleSegment[numSegments];
       //
       {
       osg::Vec4Array * colours = new osg::Vec4Array(1);
       // (*colours)[0] = osg::Vec4d(1.0,1.0,1.0,1.0);
       (*colours)[0] = osg::Vec4d(1.0,1.0,1.0,0.3);
       // (*colours)[0] = osg::Vec4d(0.0,1.0,0.0,0.3);
     
       const unsigned int numCirclePoints = 1024;
     
       const unsigned int pointsPerSegment = 
       numCirclePoints / numSegments;
     
       for (unsigned int s=0; s<numSegments; ++s) {
     
       unitOffsetCircleSegment[s] = new osg::Geometry();
     
       unitOffsetCircleSegment[s]->setColorArray(colours);
       unitOffsetCircleSegment[s]->setColorBinding(osg::Geometry::BIND_OVERALL);
     
       // const unsigned int n_points = 1024;
       //
       osg::Vec3Array * coords = new osg::Vec3Array(pointsPerSegment);
       {
       // const double dx = twopi()/numCirclePoints;
       const double dx = twopi()/(numCirclePoints-numSegments);
       double ds, dc;
       const double AU = FromUnits(1,orsa::Unit::AU);
       for (unsigned int j=0; j<pointsPerSegment; ++j) {
       // for (unsigned int j=s*pointsPerSegment; j<(s+1)*pointsPerSegment; ++j) {
       orsa::sincos(dx*(j+s*pointsPerSegment),&ds,&dc);
       (*coords)[j].set(osg::Vec3d(AU*(dc-1.0),AU*ds,0.0));
       }
       }
       //
       unitOffsetCircleSegment[s]->setVertexArray(coords);
       //
       unitOffsetCircleSegment[s]->getOrCreateStateSet()->setMode(GL_LIGHTING,osg::StateAttribute::OFF);
       // unitOffsetCircleSegment[s]->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_LOOP,0,pointsPerSegment));
       unitOffsetCircleSegment[s]->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_STRIP,0,pointsPerSegment));
       }
       }
    */
  
    // body by body, Group creation
    {
        BodyGroup::BodyList::const_iterator _b_it = _bg->getBodyList().begin();
        while (_b_it != _bg->getBodyList().end()) {
            osg::Switch * s = new osg::Switch;
            _bodySwitch[(*_b_it).get()] = s; 
            s->setUpdateCallback(new bioCallback(_bg.get(),
                                                 (*_b_it).get(),
                                                 _at.get()));
            ++_b_it;
        }
    }
  
    // body by body
    BodyGroup::BodyList::const_iterator _b_it = _bg->getBodyList().begin();
    while (_b_it != _bg->getBodyList().end()) {
    
        osg::Group * bodyGroup = new osg::Group;
        //
        bodyGroup->setName((*_b_it)->getName());
    
        // bodyGeode->getOrCreateStateSet()->setMode(GL_LIGHTING,osg::StateAttribute::OFF);
    
        /* 
           osg::Geometry * bodyGeometry = new osg::Geometry();
       
           // bodyGeometry->setUseDisplayList(false);
       
           osg::Vec3Array * coords = new osg::Vec3Array(1);
           (*coords)[0].set(osg::Vec3d(0.0f,0.0f,0.0f));
           bodyGeometry->setVertexArray(coords);
       
           osg::Vec4Array * colours = new osg::Vec4Array;
           bodyGeometry->setColorArray(colours);
           bodyGeometry->setColorBinding(osg::Geometry::BIND_OVERALL);
           colours->push_back(osg::Vec4(0.0f,1.0f,1.0f,1.0f));
       
           bodyGeometry->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::POINTS,0,2));
       
           bodyGeode->getOrCreateStateSet()->setMode(GL_LIGHTING,osg::StateAttribute::OFF);
       
           bodyGeode->addDrawable(bodyGeometry);
        */
    
        // see osgspacewarp example...
        /* 
           osg::DrawElementsUShort * points = new osg::DrawElementsUShort(GL_POINTS,1); // size
           bodyGeometry->addPrimitiveSet(points);
           (*points)[0] = 0;
        */
    
        /* 
           if (1) {
           osg::Point * point = new osg::Point;
           //
           point->setSize(1.0);
           //
           // bodyGeometry->getOrCreateStateSet()->setAttribute(point.get());
           bodyGeode->getOrCreateStateSet()->setAttribute(point);
           }
        */
    
        // label
        /* 
           { 
           osg::Vec4 characterSizeModeColor(1.0f,0.0f,0.5f,1.0f);
       
           osg::Vec3 center(0.0f,0.0f,0.0f);
       
           osgText::Text * text5 = new osgText::Text;
           text5->setColor(characterSizeModeColor);
           // text5->setFont("fonts/times.ttf");
           //text5->setCharacterSize(characterSize);
           text5->setCharacterSize(32.0f); // medium
           text5->setPosition(center);
           text5->setAxisAlignment(osgText::Text::SCREEN);
           text5->setCharacterSizeMode(osgText::Text::SCREEN_COORDS);
           text5->setText((*_b_it)->getName());
       
           bodyGeode->addDrawable(text5);
           }
        */
        //
        /* 
           if (0) {
       
           // double characterSize=FromUnits(.01,Unit::AU);
           double characterSize=FromUnits(30,Unit::KM);
       
           osg::Vec3 center(0.0,0.0,0.0);
       
           osgText::Text * text4 = new osgText::Text;
           text4->setFont("fonts/times.ttf");
           text4->setCharacterSize(characterSize);
           text4->setPosition(center);
           text4->setAxisAlignment(osgText::Text::SCREEN);
           text4->getOrCreateStateSet()->setMode(GL_LIGHTING,osg::StateAttribute::OFF);
           text4->setText((*_b_it)->getName());
       
           bodyGeode->addDrawable(text4);
           }
        */
    
        /* 
           if (0) {
           // sphere	
           osg::ShapeDrawable * shape = new osg::ShapeDrawable(new osg::Sphere(osg::Vec3d(0,0,0),
           FromUnits(1000.0,Unit::M)));
       
           shape->getOrCreateStateSet()->setMode(GL_LIGHTING,osg::StateAttribute::OFF);
       
           bodyGeode->addDrawable(shape);
           }
       
           if (0) {
           // cube
           osg::ShapeDrawable * shape = new osg::ShapeDrawable(new osg::Box(osg::Vec3d(0,0,0),
           FromUnits(.01,Unit::AU)));
       
           shape->getOrCreateStateSet()->setMode(GL_LIGHTING,osg::StateAttribute::OFF);
       
           bodyGeode->addDrawable(shape);
           }
        */
    
    
        if ((*_b_it)->getInitialConditions().inertial->localShape() != 0) {
      
            // bad style for the moment...
            switch ((*_b_it)->getInitialConditions().inertial->localShape()->getType()) {
                case orsa::Shape::SHAPE_TRI:
                {
	  
                    osg::Geode * bodyGeode = new osg::Geode;
	  
                    // bodyGeode->setName((*_b_it)->getName());
	  
                    // create a container that makes the body drawable
                    osg::Geometry * bodyGeometry = new osg::Geometry;
	  
                    bodyGeometry->setComputeBoundingBoxCallback(new BodySymmetricBoundingBoxCallback((*_b_it).get()));
	  
                    // should use osg::TriangleMesh ??
	  
                    // ORSA_DEBUG("SHAPE_TRI...");
                    osg::ref_ptr<const orsa::TriShape> _ts = 
                        (const orsa::TriShape *) (*_b_it)->getInitialConditions().inertial->localShape();
                    //
                    const TriShape::VertexVector _vertex = _ts->getVertexVector();
                    const TriShape::FaceVector   _face   = _ts->getFaceVector();
	  
                    // set the single colour so bind overall	  
                    osg::Vec4Array * colours = new osg::Vec4Array(1);
                    // 
                    // (*colours)[0] = osg::Vec4d(0.3,0.3,0.3,1);
                    (*colours)[0] = osg::Vec4d(0.6,0.6,0.6,1.0);
                    // (*colours)[0] = osg::Vec4d(0.0,0.69,0.94,0.5);
                    //
                    bodyGeometry->setColorArray(colours);
                    bodyGeometry->setColorBinding(osg::Geometry::BIND_OVERALL);
	  
                    osg::Vec3Array * coords = new osg::Vec3Array(_vertex.size());
                    // bodyGeometry->setVertexArray(coords);
                    //
                    for (unsigned int j=0; j<_vertex.size(); ++j) {
                        (*coords)[j].set(_vertex[j].getVec3d());
                    }
                    //
                    bodyGeometry->setVertexArray(coords);
	  
                    osg::Vec3Array * normals = new osg::Vec3Array(_vertex.size());
                    // bodyGeometry->setNormalArray(normals);
                    // bodyGeometry->setNormalBinding(osg::Geometry::BIND_PER_VERTEX);
                    //
                    for (unsigned int j=0; j<_vertex.size(); ++j) {
                        (*normals)[j].set(_ts->_getVertexNormal(j).getVec3d());
                    }
                    //
                    bodyGeometry->setNormalArray(normals);
                    bodyGeometry->setNormalBinding(osg::Geometry::BIND_PER_VERTEX);
	  
                    {
                        osg::DrawElementsUShort * elements = new osg::DrawElementsUShort(GL_TRIANGLES);
                        for (unsigned int j=0; j<_face.size(); ++j) {
                            const TriShape::TriIndex & _t = _face[j];
                            // osg::DrawElementsUShort * elements = new osg::DrawElementsUShort(GL_TRIANGLES);
                            //
                            elements->push_back(_t.i());
                            elements->push_back(_t.j());
                            elements->push_back(_t.k());
                            //
                            // bodyGeometry->addPrimitiveSet(elements);
                            // ORSA_DEBUG("loaded triangle %i",j);	
                        }
                        bodyGeometry->addPrimitiveSet(elements);
                    }
	  
                    // bodyGeometry->getOrCreateStateSet()->setMode(GL_LIGHTING,osg::StateAttribute::ON);
                    //
                    bodyGeode->addDrawable(bodyGeometry);
	  
                    bodyGroup->addChild(bodyGeode);
	  
                }	
                break;
                case orsa::Shape::SHAPE_ELLIPSOID:
                {
	  
                    double a,b,c;
	  
                    {
                        const orsa::EllipsoidShape * es = dynamic_cast<const orsa::EllipsoidShape * > ((*_b_it)->getInitialConditions().inertial->localShape());
                        if (!es) {
                            ORSA_ERROR("problems...");
                        }
                        es->getABC(a,b,c);
                    }
	  
                    PositionAttitudeTransform * pat = new PositionAttitudeTransform;
                    //
                    pat->setScale(osg::Vec3d(a,b,c));
	  
                    ORSA_DEBUG("make sure localShape is rotated and offset correctly! (is that even possible?)");
	  
                    osg::ShapeDrawable * shape = 
                        new osg::ShapeDrawable(new osg::Sphere(osg::Vec3d(0,0,0),1.0));
	  
                    shape->getOrCreateStateSet()->setMode(GL_LIGHTING,osg::StateAttribute::ON);
                    shape->getOrCreateStateSet()->setMode(GL_RESCALE_NORMAL,osg::StateAttribute::ON );
	  
                    osg::Geode * bodyGeode = new osg::Geode;
                    //
                    bodyGeode->addDrawable(shape);
	  
                    pat->addChild(bodyGeode);
	  
                    bodyGroup->addChild(pat);
                }
                break;
                default:
                    ORSA_WARNING("switch case not handled yet...");
                    break;
            }
      
      
            /* 
               if (0) {
               // add sample points inside shape
               const unsigned int n_points = 4000;
	 
               // GSL rng init
               gsl_rng * rnd = gsl_rng_alloc(gsl_rng_gfsr4);
               gsl_rng_set(rnd,85719); // random seed
               // const double _R_max = (*_b_it)->getShape()->boundingRadius();
               const orsa::Box boundingBox = (*_b_it)->getShape()->boundingBox();
	 
               for (unsigned int k=0; k<n_points; ++k) {
	 
               //
               const Vector v(boundingBox.getXMin()+(boundingBox.getXMax()-boundingBox.getXMin())*gsl_rng_uniform(rnd),
               boundingBox.getYMin()+(boundingBox.getYMax()-boundingBox.getYMin())*gsl_rng_uniform(rnd),
               boundingBox.getZMin()+(boundingBox.getZMax()-boundingBox.getZMin())*gsl_rng_uniform(rnd));
               //
               if ((*_b_it)->getShape()->isInside(v)) {
               osg::Geometry * bodyGeometry = new osg::Geometry();
               //
               osg::Vec4Array * colours = new osg::Vec4Array(1);
               (*colours)[0] = osg::Vec4d(0,1,1,1);
               bodyGeometry->setColorArray(colours);
               bodyGeometry->setColorBinding(osg::Geometry::BIND_OVERALL);
               //
               osg::Vec3Array * coords = new osg::Vec3Array(1);
               (*coords)[0].set(v.getVec3d());
               bodyGeometry->setVertexArray(coords);  
               //
               bodyGeometry->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::POINTS,0,coords->size()));
               //
               bodyGeometry->getOrCreateStateSet()->setMode(GL_LIGHTING,osg::StateAttribute::OFF);
               //
               bodyGeode->addDrawable(bodyGeometry);
               }
               }
	 
               // GSL rng clean
               gsl_rng_free(rnd);
               }
            */
      
            // bounding box 
            if (0) {
	
                // NOTE: plot the symmetric box as well??
	
                osg::Geometry * boxGeometry = new osg::Geometry();
	
                osg::Vec4Array * colours = new osg::Vec4Array(1);
                (*colours)[0] = osg::Vec4d(1,1,0,1.0);
                boxGeometry->setColorArray(colours);
                boxGeometry->setColorBinding(osg::Geometry::BIND_OVERALL);
	
                const orsa::Box & boundingBox = (*_b_it)->getInitialConditions().inertial->localShape()->boundingBox();
	
                osg::Vec3Array * coords = new osg::Vec3Array(24);
                // x_min face
                (*coords)[0].set(osg::Vec3d(boundingBox.getXMin(),
                                            boundingBox.getYMin(),
                                            boundingBox.getZMin()));
                (*coords)[1].set(osg::Vec3d(boundingBox.getXMin(),
                                            boundingBox.getYMax(),
                                            boundingBox.getZMin()));
                (*coords)[2].set(osg::Vec3d(boundingBox.getXMin(),
                                            boundingBox.getYMax(),
                                            boundingBox.getZMax()));
                (*coords)[3].set(osg::Vec3d(boundingBox.getXMin(),
                                            boundingBox.getYMin(),
                                            boundingBox.getZMax()));
                // x_max face
                (*coords)[4].set(osg::Vec3d(boundingBox.getXMax(),
                                            boundingBox.getYMin(),
                                            boundingBox.getZMin()));
                (*coords)[5].set(osg::Vec3d(boundingBox.getXMax(),
                                            boundingBox.getYMin(),
                                            boundingBox.getZMin()));
                (*coords)[6].set(osg::Vec3d(boundingBox.getXMax(),
                                            boundingBox.getYMin(),
                                            boundingBox.getZMin()));
                (*coords)[7].set(osg::Vec3d(boundingBox.getXMax(),
                                            boundingBox.getYMin(),
                                            boundingBox.getZMin()));
                // y_min face
                (*coords)[8].set(osg::Vec3d(boundingBox.getXMin(),
                                            boundingBox.getYMin(),
                                            boundingBox.getZMin()));
                (*coords)[9].set(osg::Vec3d(boundingBox.getXMax(),
                                            boundingBox.getYMin(),
                                            boundingBox.getZMin()));
                (*coords)[10].set(osg::Vec3d(boundingBox.getXMax(),
                                             boundingBox.getYMin(),
                                             boundingBox.getZMax()));
                (*coords)[11].set(osg::Vec3d(boundingBox.getXMin(),
                                             boundingBox.getYMin(),
                                             boundingBox.getZMax()));
                // y_max face
                (*coords)[12].set(osg::Vec3d(boundingBox.getXMin(),
                                             boundingBox.getYMax(),
                                             boundingBox.getZMin()));
                (*coords)[13].set(osg::Vec3d(boundingBox.getXMax(),
                                             boundingBox.getYMax(),
                                             boundingBox.getZMin()));
                (*coords)[14].set(osg::Vec3d(boundingBox.getXMax(),
                                             boundingBox.getYMax(),
                                             boundingBox.getZMax()));
                (*coords)[15].set(osg::Vec3d(boundingBox.getXMin(),
                                             boundingBox.getYMax(),
                                             boundingBox.getZMax()));
                // z_min face
                (*coords)[16].set(osg::Vec3d(boundingBox.getXMin(),
                                             boundingBox.getYMin(),
                                             boundingBox.getZMin()));
                (*coords)[17].set(osg::Vec3d(boundingBox.getXMax(),
                                             boundingBox.getYMin(),
                                             boundingBox.getZMin()));
                (*coords)[18].set(osg::Vec3d(boundingBox.getXMax(),
                                             boundingBox.getYMax(),
                                             boundingBox.getZMin()));
                (*coords)[19].set(osg::Vec3d(boundingBox.getXMin(),
                                             boundingBox.getYMax(),
                                             boundingBox.getZMin()));
                // z_max face
                (*coords)[20].set(osg::Vec3d(boundingBox.getXMin(),
                                             boundingBox.getYMin(),
                                             boundingBox.getZMax()));
                (*coords)[21].set(osg::Vec3d(boundingBox.getXMax(),
                                             boundingBox.getYMin(),
                                             boundingBox.getZMax()));
                (*coords)[22].set(osg::Vec3d(boundingBox.getXMax(),
                                             boundingBox.getYMax(),
                                             boundingBox.getZMax()));
                (*coords)[23].set(osg::Vec3d(boundingBox.getXMin(),
                                             boundingBox.getYMax(),
                                             boundingBox.getZMax()));
                //
                boxGeometry->setVertexArray(coords);
	
                boxGeometry->getOrCreateStateSet()->setMode(GL_LIGHTING,osg::StateAttribute::OFF);      
                boxGeometry->getOrCreateStateSet()->setMode(GL_CULL_FACE,osg::StateAttribute::OFF);      
                boxGeometry->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::QUADS,0,24));
	
                osg::PolygonMode * polymode = new osg::PolygonMode;
                polymode->setMode(osg::PolygonMode::FRONT_AND_BACK, // FRONT_AND_BACK,
                                  osg::PolygonMode::LINE);
                /* 
                   boxGeometry->getOrCreateStateSet()->setAttributeAndModes(polymode,
                   osg::StateAttribute::OVERRIDE | osg::StateAttribute::ON);
                */
                //
                boxGeometry->getOrCreateStateSet()->setAttributeAndModes(polymode,
                                                                         osg::StateAttribute::ON);
	
                boxGeometry->getOrCreateStateSet()->setAttributeAndModes(new osg::LineWidth(2.0),
                                                                         osg::StateAttribute::ON);
	
                osg::Geode * bodyGeode = new osg::Geode;
                //
                bodyGeode->addDrawable(boxGeometry);
	
                bodyGroup->addChild(bodyGeode);
            }
      
        } else {
      
            // no shape, just a point
            // ORSA_DEBUG("POINT...");
      
            // double radius = (*_b_it)->getInitialConditions().inertial->shape()->boundingRadius();
            double radius = (*_b_it)->getInitialConditions().inertial->localShape() ? (*_b_it)->getInitialConditions().inertial->localShape()->boundingRadius() : 0;
            if (radius == 0) {
                radius = FromUnits(100,Unit::METER); // hmm...
            }
      
            osg::ShapeDrawable* shape = new osg::ShapeDrawable(new osg::Sphere(osg::Vec3d(0,0,0),
                                                                               radius));
      
            shape->getOrCreateStateSet()->setMode(GL_LIGHTING,osg::StateAttribute::ON);
      
            shape->setColor(osg::Vec4d(0.6,0.6,0.6,1.0));
      
            osg::Geode * bodyGeode = new osg::Geode;
      
            bodyGeode->addDrawable(shape);
      
            bodyGroup->addChild(bodyGeode);
        }
    
        if (1) {
            // material
      
            ref_ptr<Material> m = new Material;
            //
            m->setColorMode(Material::DIFFUSE);
            m->setAmbient(Material::FRONT_AND_BACK, Vec4(0, 0, 0, 1));
            m->setSpecular(Material::FRONT_AND_BACK, Vec4(1, 1, 1, 1));
            m->setShininess(Material::FRONT_AND_BACK, 64.0f);
            //
            bodyGroup->getOrCreateStateSet()->setAttributeAndModes(m.get(), StateAttribute::ON);
      
        }
    
        {
      

            osg::PositionAttitudeTransform * attitudeTransform = new osg::PositionAttitudeTransform;
      
            attitudeTransform->setUpdateCallback(new orsaOSG::BodyAttitudeCallback(_bg.get(),
                                                                                   (*_b_it).get(),
                                                                                   _at.get()));
            attitudeTransform->addChild(bodyGroup);
      
            //
      
            _bodyAttitudeTransform[(*_b_it).get()] = attitudeTransform; 
      
            //
      
            osg::PositionAttitudeTransform * positionTransform = new osg::PositionAttitudeTransform;
      
            positionTransform->setUpdateCallback(new orsaOSG::BodyTranslationCallback(_bg.get(),
                                                                                      (*_b_it).get(),
                                                                                      _at.get()));
      
            positionTransform->addChild(attitudeTransform);
      
            //
      
            _bodyPositionTransform[(*_b_it).get()] = positionTransform; 
      
            //
      
            _bodySwitch[(*_b_it).get()]->addChild(positionTransform);
      
            //
      
            lastLightSource->addChild(_bodySwitch[(*_b_it).get()].get());
      
        }
    
        ++_b_it;
    }
  
    {
        // center of mass node
    
        osg::Group * centerOfMassGroup = new osg::Group;
    
        centerOfMassGroup->setName("centerOfMass");
    
        osg::PositionAttitudeTransform * positionTransform = new osg::PositionAttitudeTransform;
    
        positionTransform->setUpdateCallback(new orsaOSG::CenterOfMassPositionCallback(_bg.get(),
                                                                                       _at.get()));
        positionTransform->addChild(centerOfMassGroup);
    
        lastLightSource->addChild(positionTransform);
    }
  
    if (1) {
    
        // orbit
    
        BodyGroup::BodyList::const_iterator _b_it = _bg->getBodyList().begin();
        while (_b_it != _bg->getBodyList().end()) {
      
            osg::Geode * orbitGeode = new osg::Geode;
      
            orbitGeode->addDrawable(unitOffsetCircle);
            //
            /* for (unsigned int s=0; s<numSegments; ++s) {
               orbitGeode->addDrawable(unitOffsetCircleSegment[s]);
               }
            */
      
            orbitGeode->getOrCreateStateSet()->setAttributeAndModes(new osg::LineWidth(1),
                                                                    osg::StateAttribute::ON);
      
            // anti aliasing?

            /* 
               OLD CODE in ORSA 0.7
               glEnable(GL_LINE_SMOOTH);
               //
               glEnable(GL_BLEND);
               //
               glBlendFunc(GL_SRC_ALPHA,GL_ONE);
               //
               glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
               glDepthMask(GL_FALSE);
            */
            //
            /* orbitGeode->getOrCreateStateSet()->setMode(GL_LINE_SMOOTH, 
               osg::StateAttribute::ON);
               // Enable blending, select transparent bin.
               orbitGeode->getOrCreateStateSet()->setMode(GL_BLEND, 
               osg::StateAttribute::ON);
               orbitGeode->getOrCreateStateSet()->setRenderingHint(osg::StateSet::TRANSPARENT_BIN);
	 
               glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
            */
      
            // Enable depth test so that an opaque polygon will occlude a transparent one behind it.
            // orbitGeode->getOrCreateStateSet()->setMode( GL_DEPTH_TEST, osg::StateAttribute::ON );
      
            // Conversely, disable writing to depth buffer so that
            // a transparent polygon will allow polygons behind it to shine thru.
            // OSG renders transparent polygons after opaque ones.
            /* {
               osg::Depth * depth = new osg::Depth;
               depth->setWriteMask( false );
               orbitGeode->getOrCreateStateSet()->setAttributeAndModes( depth, osg::StateAttribute::ON );
               }
            */
      
            /* {
               osg::BlendFunc * bf = new
               osg::BlendFunc(osg::BlendFunc::SRC_ALPHA,
               osg::BlendFunc::ONE);
               orbitGeode->getOrCreateStateSet()->setAttributeAndModes(bf);
               }
            */
      
            // Disable conflicting modes.
            // orbitGeode->getOrCreateStateSet()->setMode( GL_LIGHTING, osg::StateAttribute::OFF );
      
      
            // orbitGeode->setComputeBoundingSphereCallback(new OrbitComputeBoundingSphereCallback);
      
            osg::MatrixTransform * mt = new osg::MatrixTransform;
      
            mt->setUpdateCallback(new OrbitNodeUpdateCallback(_bg.get(),
                                                              (*_b_it).get(),
                                                              _at.get()));
      
            mt->addChild(orbitGeode);
      
            _bodyPositionTransform[(*_b_it).get()]->addChild(mt);
      
            ++_b_it;
        }
    
    }
  
    // HUD
    lastLightSource->addChild(createHUD()); 
  
    return rootNode;    
}

osg::Node * Viz::createHUD() {
  
    osg::CameraNode * HUD = new orsaOSG::HUD(_at.get());
  
    return HUD;
}
