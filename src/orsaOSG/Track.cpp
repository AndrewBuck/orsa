#include <orsaOSG/Track.h>

// #include <orsa/attitude.h>
#include <orsa/util.h>

using namespace orsaOSG;

void TrackCallback::operator () (osg::Node * node, 
                                 osg::NodeVisitor * nv) {
  
    if (nv->getVisitorType()==osg::NodeVisitor::UPDATE_VISITOR && nv->getFrameStamp()) {
    
        osg::Geode * geode = dynamic_cast<osg::Geode * > (node);
        //
        if (geode) {
      
            const orsa::Time simulationTime = _at->getSimulationTime(nv->getFrameStamp()->getFrameNumber());
      
            const std::vector<TrackElement> & tv = _track->getTrackVector();
      
            _track->trackVectorMutex.lock();
            const unsigned int tvSize = tv.size();
            _track->trackVectorMutex.unlock();	  
      
            // ORSA_DEBUG("size: %i   old: %i",tvSize,oldTrackVectorSize);
      
            if ( (geode->getNumDrawables() == 0) || 
                 (tvSize > oldTrackVectorSize) ) {
	
                // ORSA_DEBUG("updating...");
	
                oldTrackVectorSize = tvSize;
	
                // remove previous geometries, if any
                if (_trackGeometry.get()) {
                    geode->removeDrawable(_trackGeometry.get());
                }
	
                _trackGeometry = new osg::Geometry;
	
                osg::Vec4Array * color = new osg::Vec4Array(1);
                (*color)[0] = osg::Vec4d(1.0,1.0,1.0,0.3);
                _trackGeometry->setColorArray(color);
                _trackGeometry->setColorBinding(osg::Geometry::BIND_OVERALL);
	
                _track->trackVectorMutex.lock();
                osg::Vec3Array * coords = new osg::Vec3Array(tv.size());
                //
                for (unsigned int j=0; j<tv.size(); ++j) {
                    (*coords)[j].set(tv[j].r.getVec3d());
                } 
                _track->trackVectorMutex.unlock();	  
                //
                _trackGeometry->setVertexArray(coords);
	
                _trackGeometry->getOrCreateStateSet()->setMode(GL_LIGHTING,osg::StateAttribute::OFF);
	
                geode->addDrawable(_trackGeometry.get());
	
            }
      
            _track->trackVectorMutex.lock();
            unsigned int count=0;
            std::vector<TrackElement>::const_iterator it = tv.begin();
            while (it != tv.end()) {
                if ((*it).t <= simulationTime) {
                    ++count;
                }
                ++it;
            }
            _track->trackVectorMutex.unlock();	 
      
            /* 
               {
               _track->trackVectorMutex.lock();
               ORSA_DEBUG("count: %i   size: %i   _trackGeometry->getNumPrimitiveSets(): %i",
               count,
               tv.size(),
               _trackGeometry->getNumPrimitiveSets());
               _track->trackVectorMutex.unlock();	
               }
            */
      
            if (_trackGeometry->getNumPrimitiveSets() == 1) {
                _trackGeometry->setPrimitiveSet(0,
                                                new osg::DrawArrays(osg::PrimitiveSet::LINE_STRIP,
                                                                    0,
                                                                    count));
            } else {
                _trackGeometry->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::LINE_STRIP,
                                                                    0,
                                                                    count));
            }
      
        } else {
            ORSA_DEBUG("problems");
        }
    
    }
  
    // must call any nested node callbacks and continue subgraph traversal.
    NodeCallback::traverse(node,nv);
}

////

TrackFillThread::TrackFillThread(orsa::BodyGroup  * bg,
                                 const orsa::Body * b,
                                 const orsa::Body * ref_b,
                                 const orsa::Time & dt,
                                 const bool         groundTrack) :
    QThread(),
    _bg(bg),
    _b(b),
    _ref_b(ref_b),
    _dt(dt),
    _groundTrack(groundTrack) {
  
}

void TrackFillThread::run() {
  
    dataMutex.lock();
    data.clear();
    dataMutex.unlock();
  
    doAbort = false;
  
    orsa::Time t_start, t_stop;
    _bg->getGlobalInterval(t_start,t_stop,false);
  
    orsa::Vector r_b, r_ref_b;
    orsa::Vector intersectionPoint;
    orsa::Vector intersectionNormal;
  
    TrackElement tmpTrackElement;
  
    orsa::Time t = t_start;
    while (t <= t_stop) {
    
        if (doAbort) {
            break;
        }
    
        if (_b.get() && _ref_b.get()) {
      
            if (_bg->getInterpolatedPosition(r_b,         _b.get(), t) &&
                _bg->getInterpolatedPosition(r_ref_b, _ref_b.get(), t) ) {
	
                const orsa::Vector dr = r_b - r_ref_b;
	
                if (_groundTrack) {
	  
                    if (_ref_b->getInitialConditions().inertial->localShape()) {
	    
                        /* 
                           const orsa::Vector dr_local = (_ref_b->getAttitude()) ? 
                           (_ref_b->getAttitude()->globalToLocal(t)*dr) : dr;
                        */
                        //
                        /* 
                           const orsa::Vector dr_local = (_bg->getInterpolatedAttitude(_ref_b.get(),t)) ? 
                           (_bg->getInterpolatedAttitude(_ref_b.get(),t)->globalToLocal(t)*dr) : dr;
                        */
                        //
                        // osg::ref_ptr<orsa::Attitude> attitude = new orsa::BodyAttitude(_ref_b.get(), _bg.get());
                        const orsa::Matrix g2l = orsa::globalToLocal(_ref_b.get(),_bg.get(),t);
                        const orsa::Vector dr_local = g2l * dr;
	    
                        if (_ref_b->getInitialConditions().inertial->localShape()->rayIntersection(intersectionPoint,
                                                                                                   intersectionNormal,
                                                                                                   dr_local,
                                                                                                   (-dr_local).normalized(),
                                                                                                   false)) {
	      
                            tmpTrackElement.t = t;
                            //
                            tmpTrackElement.r = 1.05*intersectionPoint;
                            // tmpTrackElement.r = orsa::FromUnits(300,orsa::Unit::KM)*(intersectionPoint.normalized());
                            //
                            // trackVector.push_back(tmpTrackElement);
                            //	      
                            dataMutex.lock();
                            data.push_back(tmpTrackElement);
                            dataMutex.unlock();
	      
                        } else {
                            ORSA_DEBUG("problems...");
                        }
	    
                    } else {
                        ORSA_DEBUG("groundTrack without a shape?");
                    }
	  
                } else {
	  
                    tmpTrackElement.t = t;
                    tmpTrackElement.r = dr;
                    //
                    // trackVector.push_back(tmpTrackElement);
                    //
                    dataMutex.lock();
                    data.push_back(tmpTrackElement);
                    dataMutex.unlock();
                }
	
            } else {
                ORSA_DEBUG("problems...");
            }
      
        } else if (_b.get()) {
      
            /* 
               if (_bg->getInterpolatedPosition(r_b, _b.get(), t)) {
	 
               tmpTrackElement.t = t;
               tmpTrackElement.r = r_b - _at->centralBodyPosition(t_start);
               //
               trackVector.push_back(tmpTrackElement);
	 
               ORSA_DEBUG("tmpTrackElement.r: %f",
               tmpTrackElement.r.length());
	 
               } else {
               ORSA_DEBUG("problems...");
               }
            */
            //
            ORSA_DEBUG("this case cannot be pre-computed, because of the moving AnimationTime::centralBodyPosition(...), and requires different code.");
      
        } else {
            ORSA_ERROR("problems");      
        }
    
        t += _dt;
    }
  
}

////

void Track::trackFill() {
  
    trackFillThread->dataMutex.lock();
  
    trackVectorMutex.lock();
  
    // this could be optimized a bit...
    trackVector = trackFillThread->data;
  
    trackVectorMutex.unlock();
  
    trackFillThread->dataMutex.unlock();
}
