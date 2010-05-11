#include <orsa/shape.h>
#include <orsa/unit.h>

#include <algorithm>

using namespace orsa;

// Shape

// nothing here...

// TriShape

/* 
   const Vector & TriShape::getNormal(const unsigned int face_index,
   const unsigned int vertex_index) const {
   switch(_normal_type) {
   case NORMAL_RADIAL:       return _getNormalRadial(vertex_index);      break;
   case NORMAL_FACE:         return _getNormalFace(face_index);          break;
   case NORMAL_FACE_AVERAGE: return _getNormalFaceAverage(vertex_index); break;
   };
   
   {
   ORSA_DEBUG("we should not be here...");
   _dummy_normal.set(0,0,1);
   return _dummy_normal;
   }
   }
*/

/* 
   bool TriShape::_updateNormal() const {
   switch(_normal_type) {
   case NORMAL_RADIAL:        break;
   case NORMAL_FACE:          break;
   case NORMAL_FACE_AVERAGE:  break;
   };
   return true;
   }
*/

const Vector & TriShape::_getVertexNormal(const unsigned int vertex_index) const {
    if (_vertex_normal.size() != _vertex.size()) {
        _vertex_normal.resize(_vertex.size());
        Vector _n;
        for (unsigned int _v=0; _v<_vertex.size(); ++_v) {
            _n.set(0,0,0);
            for (unsigned int _f=0; _f<_face.size(); ++_f) {
                if ( (_face[_f].i() == _v) ||
                     (_face[_f].j() == _v) ||
                     (_face[_f].k() == _v) ) {
                    _n += _getFaceNormal(_f);
                }
            }
            _vertex_normal[_v] = _n.normalized();
        }
    }
    return _vertex_normal[vertex_index];
}

const Vector & TriShape::_getFaceNormal(const unsigned int face_index) const {
    if (_face_normal.size() != _face.size()) {
        _face_normal.resize(_face.size());
        for (unsigned int _f=0; _f<_face.size(); ++_f) {
            const TriIndex & _t = _face[_f];
            _face_normal[_f] = externalProduct(_vertex[_t.k()]-_vertex[_t.j()],
                                               _vertex[_t.i()]-_vertex[_t.j()]).normalized();
        }
    }
    return _face_normal[face_index];
}

double TriShape::_getFaceArea(const unsigned int face_index) const {
    if (_face_area.size() != _face.size()) {
        _face_area.resize(_face.size());
        const double half = 0.5;
        for (unsigned int _f=0; _f<_face.size(); ++_f) {
            const TriIndex & _t = _face[_f];
            _face_area[_f] = half * externalProduct(_vertex[_t.k()]-_vertex[_t.j()],
                                                    _vertex[_t.i()]-_vertex[_t.j()]).length();
        }
    }
    return _face_area[face_index];
}

bool TriShape::_updateCache() const {
    if ((!_r_min.isSet()) || (!_r_max.isSet())) {
        VertexVector::const_iterator _it = _vertex.begin();
        double _d2, _d2_min, _d2_max;
        // init
        _d2_min = _d2_max = (*_it).lengthSquared();
        while (_it != _vertex.end()) {
            _d2 = (*_it).lengthSquared();
            if (_d2 < _d2_min) {
                _d2_min = _d2;
            }
            if (_d2 > _d2_max) {
                _d2_max = _d2;
            }      
            ++_it;
        }
        _r_min.set(sqrt(_d2_min));
        _r_max.set(sqrt(_d2_max));
        //
        // std::cerr << "r_min: " << FromUnits(_r_min.getRef(),Unit::KM,-1) << " KM" << std::endl;
    }
    //
    if (!_boundingBox.isSet()) {
        double xMin, xMax, yMin, yMax, zMin, zMax;
        xMin = yMin = zMin =  _r_max.get();
        xMax = yMax = zMax = -_r_max.get();
        VertexVector::const_iterator _it = _vertex.begin();
        while (_it != _vertex.end()) {
            const Vector & v = (*_it);
            if (v.getX() < xMin) xMin = v.getX();
            if (v.getX() > xMax) xMax = v.getX();
            if (v.getY() < yMin) yMin = v.getY();
            if (v.getY() > yMax) yMax = v.getY();
            if (v.getZ() < zMin) zMin = v.getZ();
            if (v.getZ() > zMax) zMax = v.getZ();
            ++_it;
        }
        _boundingBox.set(xMin, xMax, yMin, yMax, zMin, zMax);
    }
    //
    if (!_symmetricBoundingBox.isSet()) {
        const double xL = std::max(fabs(_boundingBox.getXMin()),
                                   fabs(_boundingBox.getXMax()));
        const double yL = std::max(fabs(_boundingBox.getYMin()),
                                   fabs(_boundingBox.getYMax()));
        const double zL = std::max(fabs(_boundingBox.getZMin()),
                                   fabs(_boundingBox.getZMax()));
        _symmetricBoundingBox.set(-xL,xL,-yL,yL,-zL,zL);
    }
    //
    if ((!_delta_min.isSet()) || (!_delta_max.isSet())) {
        FaceVector::const_iterator _it = _face.begin();
        double _d2_ij, _d2_ik, _d2_jk, _d2_tmp;
        // init
        double _d2_min = _r_max.getRef()*_r_max.getRef();
        double _d2_max = 0;
        while (_it != _face.end()) {
            _d2_ij = (_vertex[(*_it).i()]-_vertex[(*_it).j()]).lengthSquared();
            _d2_ik = (_vertex[(*_it).i()]-_vertex[(*_it).k()]).lengthSquared();
            _d2_jk = (_vertex[(*_it).j()]-_vertex[(*_it).k()]).lengthSquared();
            //
            _d2_tmp = std::min(_d2_ij,std::min(_d2_ik,_d2_jk));
            if (_d2_tmp < _d2_min) {
                _d2_min = _d2_tmp;
            }
            //
            _d2_tmp = std::max(_d2_ij,std::max(_d2_ik,_d2_jk));
            if (_d2_tmp > _d2_max) {
                _d2_max = _d2_tmp;
            }
            //
            ++_it;
        }
        _delta_min.set(sqrt(_d2_min));
        _delta_max.set(sqrt(_d2_max));
    } 
    //
    // ORSA_DEBUG("need to fill the _ref_point_inside_model vector somewhere...");
    //
    return true;
}

bool TriShape::isInside(const Vector & v) const {
    // choose one
    //
    // return _isInside_useLineMethod(v);
    // return _isInside_useNormalMethod(v);
    return _isInside_useFaceMethod(v);
}

bool TriShape::_isInside_useLineMethod(const Vector & v) const {
    _updateCache();
    if (v.lengthSquared() > (_r_max.getRef()*_r_max.getRef())) {
        return false;
    } else if (v.lengthSquared() < (_r_min.getRef()*_r_min.getRef())) {
        return true;
    }
    //
    if (!_boundingBox.isInside(v)) {
        return false;
    }
    //
    ORSA_DEBUG("check this method...");
    // search for closest ref. point
    if (_ref.size() == 0) {
        ORSA_ERROR("no reference points available");
        return false; 
    }
    double _min_d2  = (v-_ref[0]).lengthSquared();
    Vector _min_ref = _ref[0];
    double _tmp_dx, _tmp_d2;
    for (unsigned int _k=1; _k<_ref.size(); ++_k) {
        _tmp_dx = v.getX()-_ref[_k].getX();
        if ((_tmp_dx*_tmp_dx) > _min_d2) {
            // fast skip     
            continue;
        }
        //
        _tmp_d2 = (v-_ref[_k]).lengthSquared();
        if (_tmp_d2 < _min_d2) {
            _min_ref = _ref[_k];
            _min_d2 = _tmp_d2;
        }
    }
  
    // search for vertex point close to the line
    // passing through ref. point and v
    const double _max_d_from_line  = _delta_max.getRef();
    const double _max_d2_from_line = _max_d_from_line*_max_d_from_line;
    const Vector _unit_v = (v - _min_ref).normalized();
    // index of vertex close to the line
    std::vector<unsigned int> _index;
    for (unsigned int _j=0; _j<_vertex.size(); ++_j) {
        const Vector & _p = _vertex[_j];
        if ((_p-_min_ref)*_unit_v < 0) continue;
        {
            // fast skip test
            if (_unit_v.getX() > 0) {
                if ((_p.getX() - _min_ref.getX()) < -_max_d_from_line) continue;
                if ((_p.getX() - v.getX())        >  _max_d_from_line) continue;
            } else {
                if ((_p.getX() - _min_ref.getX()) >  _max_d_from_line) continue;
                if ((_p.getX() - v.getX())        < -_max_d_from_line) continue;
            }
      
            if (_unit_v.getY() > 0) {
                if ((_p.getY() - _min_ref.getY()) < -_max_d_from_line) continue;
                if ((_p.getY() - v.getY())        >  _max_d_from_line) continue;
            } else {
                if ((_p.getY() - _min_ref.getY()) >  _max_d_from_line) continue;
                if ((_p.getY() - v.getY())        < -_max_d_from_line) continue;
            }
      
            if (_unit_v.getZ() > 0) {
                if ((_p.getZ() - _min_ref.getZ()) < -_max_d_from_line) continue;
                if ((_p.getZ() - v.getZ())        >  _max_d_from_line) continue;
            } else {
                if ((_p.getZ() - _min_ref.getZ()) >  _max_d_from_line) continue;
                if ((_p.getZ() - v.getZ())        < -_max_d_from_line) continue;
            }
        }
    
        const Vector _w = _p - v;
        _tmp_d2 = (_w - (_w*_unit_v)*_unit_v).lengthSquared();
        if (_tmp_d2 < _max_d2_from_line) {
            _index.push_back(_j);
        }
    }
    // check if the line is passing too close to the edge of the shape
    if (_index.size() > 0) {
        const double _line_step = 0.5*_delta_min.getRef();
        Vector _p = _min_ref;
        while ((_p-_min_ref).lengthSquared() < (_min_d2+(_line_step*_line_step))) {
            if ((_p-_min_ref).lengthSquared() > _min_d2) {
                _p = v;
            }
      
            for (unsigned int _j=0; _j<_index.size(); ++_j) {
                const unsigned int _vertex_index = _index[_j];
                const Vector _vi = _vertex[_vertex_index];
                const double _ref_scalar_product = _getVertexNormal(_vertex_index)*(_min_ref - _vi);
                const double _test_point_scalar_product = _getVertexNormal(_vertex_index)*(_p - _vi);
	
                if (_ref_scalar_product*_test_point_scalar_product < 0) {
                    // different sign!
                    return false;
                }
            }
            _p += _unit_v*_line_step;
        }
        // all points checked, none close enough to the shape edge
        return true;
    } else {
        return true;
    }
    return false;
}

bool TriShape::_isInside_useNormalMethod(const Vector & v) const {
    _updateCache();
    if (v.lengthSquared() > (_r_max.getRef()*_r_max.getRef())) {
        /* 
           ORSA_DEBUG("fast out: v.length()=%f > _r_max.getRef()=%f",
           v.length(),
           _r_max.getRef());
        */
        return false;
    } else if (v.lengthSquared() < (_r_min.getRef()*_r_min.getRef())) {
        // ORSA_DEBUG("fast in");
        return true;
    }
  
    // search closest vertex
    // should check for _vertex.size() > 0...
    double _vertex_d2 = (v-_vertex[0]).lengthSquared();
    unsigned int _vertex_index = 0;
    double _tmp_dx, _tmp_d2;
    for (unsigned int _k=0; _k<_vertex.size(); ++_k) {
        _tmp_dx = v.getX()-_vertex[_k].getX();
        if ((_tmp_dx*_tmp_dx) > _vertex_d2) {
            // fast skip     
            continue;
        }
        //
        _tmp_d2 = (v-_vertex[_k]).lengthSquared();
        if (_tmp_d2 < _vertex_d2) {
            _vertex_index = _k;
            _vertex_d2 = _tmp_d2;
        }
    }
  
    // here's the test
    if (_getVertexNormal(_vertex_index)*(_vertex[_vertex_index]-v) > 0) {
        // ORSA_DEBUG("last in");
        return true;
    }
  
    // ORSA_DEBUG("last out");
    return false;
}

bool TriShape::_isInside_useFaceMethod(const Vector & v) const {
    _updateCache();
    if (v.lengthSquared() > (_r_max.getRef()*_r_max.getRef())) {
        /* 
           ORSA_DEBUG("fast out: v.length()=%f > _r_max.getRef()=%f",
           v.length(),
           _r_max.getRef());
        */
        return false;
    } else if (v.lengthSquared() < (_r_min.getRef()*_r_min.getRef())) {
        // ORSA_DEBUG("fast in");
        return true;
    }
  
    // for all faces containing the vertex closet to "v", check that
    // the relative distance between "v" and each face vertex lays inside the face
    unsigned int id = closestVertexIndex(v);
    for (unsigned int fi=0; fi<_face.size(); ++fi) {
        if ( (_face[fi].i() == id) ||
             (_face[fi].j() == id) ||
             (_face[fi].k() == id) ) {
            if (_getFaceNormal(fi)*(_vertex[_face[fi].i()]-v) > 0) { return true; }
            if (_getFaceNormal(fi)*(_vertex[_face[fi].j()]-v) > 0) { return true; }
            if (_getFaceNormal(fi)*(_vertex[_face[fi].k()]-v) > 0) { return true; }
        }
    }
  
    return false;
}



const Vector & TriShape::closestVertex(const Vector & v) const {
    return _vertex[closestVertexIndex(v)];
}

unsigned int TriShape::closestVertexIndex(const Vector & v) const {
    _updateCache();
    //
    /* unsigned int _vertex_index;
       if (_old_vertex_index < _vertex.size()) {
       _vertex_index = _old_vertex_index;
       } else {
       _vertex_index = 0;
       }
    */
    //
    unsigned int _vertex_index = _old_closest_vertex_index;
    double _vertex_d2 = (v-_vertex[_vertex_index]).lengthSquared();
    //
    // double _tmp_dx, _tmp_dy, _tmp_dz, _tmp_d2;
    double _tmp_dx, _tmp_d2;
    for (unsigned int _k=0; _k<_vertex.size(); ++_k) {
        const Vector _tmp_v = v - _vertex[_k];
        // _tmp_dx = v.getX()-_vertex[_k].getX();
        _tmp_dx = _tmp_v.getX();
        if ((_tmp_dx*_tmp_dx) > _vertex_d2) {
            // fast skip     
            continue;
        }
        //
        // _tmp_dy = v.getY()-_vertex[_k].getY();
        /* 
           _tmp_dy = _tmp_v.getY();
           if ((_tmp_dy*_tmp_dy) > _vertex_d2) {
           // fast skip     
           continue;
           }
        */
        //
        // _tmp_dz = v.getZ()-_vertex[_k].getZ();
        /* 
           _tmp_dz = _tmp_v.getZ();
           if ((_tmp_dz*_tmp_dz) > _vertex_d2) {
           // fast skip     
           continue;
           }
        */
        //
        // _tmp_d2 = (v-_vertex[_k]).lengthSquared();
        _tmp_d2 = _tmp_v.lengthSquared();
        if (_tmp_d2 < _vertex_d2) {
            _vertex_index = _k;
            _vertex_d2 = _tmp_d2;
        }
    }
    // update local static variable
    _old_closest_vertex_index = _vertex_index;
    return _vertex_index;
}

bool TriShape::rayIntersection(orsa::Vector & intersectionPoint,
                               orsa::Vector & normal,
                               const orsa::Vector & P,
                               const orsa::Vector & u,
                               const bool fullLine) const {
    for (unsigned int j=0; j<_face.size(); ++j) {
        // ORSA_DEBUG("considering face %i / %i",j,_face.size());
        const TriIndex & t = _face[j];
        if (rayIntersectsTriangle(intersectionPoint,
                                  P,
                                  u,
                                  _vertex[t.i()],
                                  _vertex[t.j()],
                                  _vertex[t.k()],
                                  fullLine)) {
      
            normal = _getFaceNormal(j);
            if ((normal*u) < 0) {
                return true;
            }
        }
    }
    return false;
}

// each vertex with delta > deltaMax is not computed
bool TriShape::vertexIlluminationAngles(const orsa::Vector & lightSource,
                                        const orsa::Vector & observerPosition,
                                        double       & phase,
                                        AngleVector        & i, 
                                        AngleVector        & e,
                                        AngleVector        & delta,
                                        const double & deltaMax,
                                        const bool includeShadows) const {
  
    phase = acos((lightSource.normalized())*(observerPosition.normalized()));
  
    const unsigned int size = _vertex.size();
    //
    i.resize(size);
    e.resize(size);
    delta.resize(size);
  
    for (unsigned int k=0; k<size; ++k) {

        const Vector & vertex            = _vertex[k];
        const Vector   observerDirection = (observerPosition - vertex).normalized();
        //
        delta[k] = acos(observerDirection*(observerPosition.normalized()));
        //
        //! REMEMBER: deltaMax is used to skip extra computations where it is already
        //! known that this vertex is not going to be used (out of instrument field).
        //
        if (delta[k] > deltaMax) {
            i[k] = e[k] = pi();
            continue;
        }
    
        Vector intersectionPoint;
    
        const Vector &   vertexNormal = _getVertexNormal(k);
        const Vector   lightDirection = (lightSource - vertex).normalized();
        //
        i[k] = acos(vertexNormal*lightDirection);
        e[k] = acos(vertexNormal*observerDirection);
    
        if (includeShadows) {
            if ( (i[k] < halfpi()) && 
                 (e[k] < halfpi()) ) {
                for (unsigned int j=0; j<_face.size(); ++j) {
                    const TriIndex & t = _face[j];
                    // skip all faces containing this vertex
                    if ( (t.i() == k) || 
                         (t.j() == k) || 
                         (t.k() == k) ) {
                        continue;
                    }
                    //
                    if (rayIntersectsTriangle(intersectionPoint,
                                              vertex,
                                              lightDirection,
                                              _vertex[t.i()],
                                              _vertex[t.j()],
                                              _vertex[t.k()],
                                              false)) {
                        i[k] = pi();
                        //
                        break;
                    }
                }
            }
        }
    }
  
    return true;
}

// each face with every vertex with delta > deltaMax is not computed
/* 
   bool TriShape::faceIlluminationAngles(const orsa::Vector & lightSource,
   const orsa::Vector & observerPosition,
   double       & phase,
   AngleVector        & i, 
   AngleVector        & e,
   AngleVector        & delta,
   const double & deltaMax,
   const bool includeShadows) const {
   
   phase = acos((lightSource.normalized())*(observerPosition.normalized()));
   
   const unsigned int size = _face.size();
   //
   i.resize(size);
   e.resize(size);
   delta.resize(size);
   
   for (unsigned int k=0; k<size; ++k) {
   
   delta[k] = pi();
   const TriIndex & t = _face[k];
   //
   {
   const Vector & vertex = _vertex[t.i()];
   const Vector   observerDirection = (observerPosition - vertex).normalized();
   
   delta[k] = std::min(acos(observerDirection*(observerPosition.normalized())),
   delta[k]);
   }
   //
   {
   const Vector & vertex = _vertex[t.j()];
   const Vector   observerDirection = (observerPosition - vertex).normalized();
   
   delta[k] = std::min(acos(observerDirection*(observerPosition.normalized())),
   delta[k]);
   }
   //
   {
   const Vector & vertex = _vertex[t.k()];
   const Vector   observerDirection = (observerPosition - vertex).normalized();
   
   delta[k] = std::min(acos(observerDirection*(observerPosition.normalized())),
   delta[k]);
   }
   
   //! REMEMBER: deltaMax is used to skip extra computations where it is already
   //! known that this vertex is not going to be used (out of instrument field).
   //
   if (delta[k] > deltaMax) {
   i[k] = e[k] = pi();
   continue;
   }
   
   Vector intersectionPoint;
   
   // const Vector &   vertexNormal = _getVertexNormal(k);
   const Vector & faceNormal = _getFaceNormal(k);
   //
   // averaged...
   const Vector observerDirection = ( (observerPosition - _vertex[t.i()]).normalized() +
   (observerPosition - _vertex[t.j()]).normalized() +
   (observerPosition - _vertex[t.k()]).normalized() ).normalized();
   //
   // const Vector   lightDirection = (lightSource - vertex).normalized();
   // averaged lightDirection, maybe something better is needed (weighted average?)
   const Vector lightDirection = ( (lightSource - _vertex[t.i()]).normalized() +
   (lightSource - _vertex[t.j()]).normalized() +
   (lightSource - _vertex[t.k()]).normalized() ).normalized();
   //
   i[k] = acos(faceNormal*lightDirection);
   e[k] = acos(faceNormal*observerDirection);
   
   if (includeShadows) {
   ORSA_ERROR("implementation needed...");
   }
   }
   
   return true;
   }
*/

bool orsa::rayIntersectsTriangle(orsa::Vector & intersectionPoint,
                                 const orsa::Vector & P,
                                 const orsa::Vector & u,
                                 const orsa::Vector & t1,
                                 const orsa::Vector & t2,
                                 const orsa::Vector & t3,
                                 const bool fullLine) {
    const Vector t21 = t2 - t1;
    const Vector t31 = t3 - t1;
    //
    const Vector planeNormal = externalProduct(t21,t31).normalized();
    //
    const double pDistance   = (t1-P)*planeNormal*((planeNormal*u > 0) ? 1 : -1);
    //
    // ORSA_DEBUG("pDistance: %f",pDistance);
    //
  
    if ((!fullLine) && (pDistance <= 0)) {
        // ORSA_DEBUG("negative distance without fullLine...");
        return false;
    }  
  
    // projection of P on the plane, along the u direction
    const Vector Q = P + u * pDistance / (fabs(u*planeNormal));
    //
    intersectionPoint = Q;
    //
    const Vector Qt1 = Q - t1;
  
    /* 
       {
       // debug
       ORSA_DEBUG("original normal: %f %f %f",
       planeNormal.getX(),
       planeNormal.getY(),
       planeNormal.getZ());
     
       const Vector newNormal = externalProduct(t2-Q,t3-Q).normalized();
     
       ORSA_DEBUG("new normal.....: %f %f %f",
       planeNormal.getX(),
       planeNormal.getY(),
       planeNormal.getZ());
       }
    */
  
    if (Qt1.lengthSquared() > std::max(t21.lengthSquared(),
                                       t31.lengthSquared())) {
        // Q is too far from t1
        // ORSA_DEBUG("Q is too far from t1");
    
        return false;
    
    } else {
    
        const Vector a = t21;
        const Vector b = t31;
        const Vector c = Qt1;
        //
        const double ab = a*b;
        const double ac = a*c;
        const double bc = b*c;
        //
        const double a2 = a*a;
        const double b2 = b*b;
        const double c2 = c*c;
        //
        // const double beta  = (ab*ac-bc*a2)/(ab*ab-a2*b2);
        //
        double beta_tmp;
        //
        {
            const double denom_ab = (ab*ab-a2*b2);
            const double denom_ac = (ac*ac-a2*c2);
            //
            if (denom_ab != 0) {
                beta_tmp = (ab*ac-bc*a2)/denom_ab;
            } else if (denom_ac != 0) {
                beta_tmp = (ab*ac-bc*a2)/denom_ac;
            } else {
                ORSA_ERROR("beta: singular case...");
                return false;
            }
        }
        //
        const double beta = beta_tmp;
        //
        // const double alpha = (bc-beta*b2)/ab;
        //
        double alpha_tmp;
        //
        if (ab != 0) {
            alpha_tmp = (bc-beta*b2)/ab;
        } else if (ac != 0) {
            alpha_tmp = (c2-beta*bc)/ac;
        } else {
            ORSA_ERROR("alpha: singular case...");
            return false;
        }
        //
        const double alpha = alpha_tmp;
    
        /* 
           {
           // debug
       
           ORSA_DEBUG("original Qt1: %f %f %f    l: %f",
           Qt1.getX(),
           Qt1.getY(),
           Qt1.getZ(),
           Qt1.length());
       
           orsa::Vector newQt1 = alpha * a + beta * b;
       
           ORSA_DEBUG("new Qt1.....: %f %f %f    l: %f",
           newQt1.getX(),
           newQt1.getY(),
           newQt1.getZ(),
           newQt1.length());
           }
        */
    
        if ( (alpha >= 0) &&
             (beta >= 0)  && 
             (alpha+beta <= 1) ) {
      
            /* 
               ORSA_DEBUG("alpha: %f   beta: %f   alpha+beta: %f",
               alpha,
               beta,
               alpha+beta);
            */
      
            return true;
        }
    }
  
    return false;
}

// LatLonShape

bool LatLonShape::_updateCache() const {
    if ((!_r_min.isSet()) || (!_r_max.isSet())) {  
        if (_rt.size() > 0) {
            if (_rt[0].size() > 0) {
                double _d, _d_min, _d_max;
                // init
                _d_min = _d_max = _rt[0][0];
                //
                for (unsigned int _j=0; _j<_rt.size(); ++_j) {
                    for (unsigned int _k=0; _k<_rt[_j].size(); ++_k) {
                        _d = _rt[_j][_k];
                        if (_d < _d_min) {
                            _d_min = _d;
                        }
                        if (_d > _d_max) {
                            _d_max = _d;
                        }      
                    }
                }
                _r_min.set(_d_min);
                _r_max.set(_d_max);
            } else {
                ORSA_ERROR("inconsistent data problem");
                return false;
            }
        } else {
            ORSA_ERROR("inconsistent data problem");
            return false;
        }
    }
    //
    if (!_boundingBox.isSet()) {
        ORSA_DEBUG("code needed here!!");
    }
    if (!_symmetricBoundingBox.isSet()) {
        ORSA_DEBUG("code needed here!!");
    }
    return true;
}

bool LatLonShape::isInside(const Vector &) const {
    ORSA_DEBUG("code needed here!");
    return false;
}

// EllipsoidShape

bool EllipsoidShape::isInside(const Vector & v) const {
    /* 
       return (((v.getX()*v.getX())/_a2 +
       (v.getY()*v.getY())/_b2 +
       (v.getZ()*v.getZ())/_c2 ) <= 1);
    */
    //
    return (((v.getX()*v.getX())*_am2 +
             (v.getY()*v.getY())*_bm2 +
             (v.getZ()*v.getZ())*_cm2 ) <= 1);
}

bool EllipsoidShape::_updateCache() const {
    if ((!_r_min.isSet()) || (!_r_max.isSet())) {
        _r_min.set(std::min(std::min(_a,_b),_c));
        _r_max.set(std::max(std::max(_a,_b),_c));
    }
    if (!_boundingBox.isSet()) {
        _boundingBox.set(-_a, _a, -_b, _b, -_c, _c);
    }
    if (!_symmetricBoundingBox.isSet()) {
        _symmetricBoundingBox = _boundingBox;
    }
    return true;
}

bool EllipsoidShape::rayIntersection(orsa::Vector & intersectionPoint,
                                     orsa::Vector & normal,
                                     const orsa::Vector & P,
                                     const orsa::Vector & u,
                                     const bool fullLine) const {
  
    const double polyA = 
        u.getX()*u.getX()*_am2 +
        u.getY()*u.getY()*_bm2 +
        u.getZ()*u.getZ()*_cm2 ;
  
    const double polyB = 2*( P.getX()*u.getX()*_am2 +
                             P.getY()*u.getY()*_bm2 +
                             P.getZ()*u.getZ()*_cm2 );
  
    const double polyC = 
        P.getX()*P.getX()*_am2 +
        P.getY()*P.getY()*_bm2 +
        P.getZ()*P.getZ()*_cm2 - 1;
  
    const double delta = polyB*polyB - 4*polyA*polyC;
  
    if (delta < 0) {
        return false;
    }
  
    // smallest gamma: Ellipsoid point closest to 'P' along ray 'u'
    //
    // const double gamma = (-polyB-sqrt(delta))/(2.0*polyA);
    //
    const double s1 = (-polyB+sqrt(delta))/(2*polyA);
    const double s2 = (-polyB-sqrt(delta))/(2*polyA);
  
    // intersectionPoint = P + gamma*u;
    //
    if (fullLine) {
        if (fabs(s1) < fabs(s2)) {
            intersectionPoint = P + s1*u;
        } else {
            intersectionPoint = P + s2*u;
        }
    } else {
        const double sMin = std::min(s1,s2);
        const double sMax = std::max(s1,s2);
        if (sMin >= 0) {
            intersectionPoint = P + sMin*u;
        } else if (sMax >= 0) {
            // inside the shape, going out... 
            intersectionPoint = P + sMax*u;
        } else {
            return false;
        }
    }
  
    normal = normalVector(intersectionPoint);
  
    return true;
}
