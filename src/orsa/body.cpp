#include <orsa/body.h>

#include <orsa/print.h>
#include <orsa/util.h>

using namespace orsa;

const orsa::Shape * InertialBodyProperty::localShape() const {
  
    if (originalShape() == 0) return 0;
  
    /* ORSA_DEBUG("%x",ls_cache.get());
       if (ls_cache.get()) {
       ORSA_DEBUG("%x",ls_cache->data.get());
       }
    */
  
    if (ls_cache.get() != 0) {
        if (ls_cache->data.get() != 0) {
            bool valid=true;
            //
            if (ls_cache->data->originalShape.get() != originalShape()) valid=false;
            if (ls_cache->data->cm                  != centerOfMass())  valid=false;
            if (ls_cache->data->s2l                 != shapeToLocal())  valid=false;
            //
            if (valid) {
                // ORSA_DEBUG("cached");
                return ls_cache->data->localShape.get();
            }
        } 
    }
  
    // OK, old local_shape not valid, need to compute a new one
  
    // ORSA_DEBUG("not-cached");
  
    // this should never been neded, defeats the purpose of having a cache
    // ls_cache = new LocalShapeCache;
  
    ls_cache->data = new LocalShapeData;
  
    ls_cache->data->originalShape = originalShape();
    ls_cache->data->cm            = centerOfMass();
    ls_cache->data->s2l           = shapeToLocal();
  
    const orsa::Vector cm  = centerOfMass();
    const orsa::Matrix s2l = shapeToLocal();
  
    // bad style for the moment...
    switch (originalShape()->getType()) {
        case orsa::Shape::SHAPE_TRI:
        {
            osg::ref_ptr<const orsa::TriShape> ts = (const orsa::TriShape *) originalShape();
      
            const TriShape::VertexVector ts_vertex = ts->getVertexVector();
            const TriShape::FaceVector   ts_face   = ts->getFaceVector();
      
            TriShape::VertexVector local_vertex;
      
            local_vertex.resize(ts_vertex.size());
      
            for (unsigned int j=0; j<ts_vertex.size(); ++j) {
                local_vertex[j] = s2l*(ts_vertex[j]-cm);
            }
      
            TriShape::FaceVector local_face = ts_face;
      
            ls_cache->data->localShape = new orsa::TriShape(local_vertex,local_face);
        }
        break;
        case orsa::Shape::SHAPE_ELLIPSOID:
        {
            ORSA_DEBUG("CODE NEEDED!!");
            ORSA_DEBUG("setting same ellipsoid for the moment (BAD BAD BAD)");
            ls_cache->data->localShape = originalShape();
        }
        break;
        default:
            ORSA_WARNING("switch case not handled yet...   CODE NEEDED!!!");
            ls_cache->data->localShape = originalShape();
            break;
    }
  
    return ls_cache->data->localShape.get();
}


// IBPS 

IBPS::IBPS() {
  
    // ORSA_DEBUG("creating new IBPS, address: %x",this);
  
    tmp = false;
}

IBPS::IBPS(const IBPS & ibps) {
  
    /* ORSA_DEBUG("creating IBPS, address: %x   copy of address: %x",this,&ibps);
       ORSA_DEBUG("%x.translational.get(): %x   %x.translational.get(): %x",
       &ibps,ibps.translational.get(),
       this,translational.get());
    */
  
    if (ibps.time.isSet()) time = ibps.time.getRef();
  
    if (ibps.inertial.get()) {
        if (ibps.inertial->dynamic()) {
            inertial = ibps.inertial->clone();
        } else {
            inertial = ibps.inertial;
        }
    } else {
        inertial = 0;
    }
  
    if (ibps.translational.get()) {
        if (ibps.translational->dynamic()) {
            translational = ibps.translational->clone();
        } else {
            translational = ibps.translational;
        }
    } else {
        translational = 0;
    }
  
    if (ibps.rotational.get()) {
        if (ibps.rotational->dynamic()) {
            rotational = ibps.rotational->clone();
        } else {
            rotational = ibps.rotational;
        }
    } else {
        rotational = 0;
    }
  
    tmp = ibps.tmp;
  
    /* ORSA_DEBUG(" --LEAVING-- %x.translational.get(): %x   %x.translational.get(): %x",
       &ibps,ibps.translational.get(),
       this,translational.get());
    */
}

IBPS::~IBPS() { 
    /* ORSA_DEBUG("destroying IBPS address: %x   %x.translational.get(): %x",
       this,
       this,translational.get());
    */
}

const IBPS & IBPS::operator = (const IBPS & ibps) {
  
    /* 
       ORSA_DEBUG("copying IBPS, from address: %x   to address %x",&ibps,this);
       ORSA_DEBUG("%x.translational.get(): %x   %x.translational.get(): %x",
       &ibps,ibps.translational.get(),
       this,translational.get());
    */
  
    if (ibps.time.isSet()) time = ibps.time.getRef();
    /* 
       if (ibps.dynamic()) {
       // time = ibps.time.getRef();
       // ORSA_DEBUG("copy operator, time: %.6f",time.getRef().get_d());
       } 
    */
    if (ibps.inertial.get()) {
        inertial = ibps.inertial->clone();
    } else {
        inertial = 0;
    }
    if (ibps.translational.get()) {
        translational = ibps.translational->clone();
    } else {
        translational = 0;
    }
    if (ibps.rotational.get()) {
        rotational = ibps.rotational->clone();
    } else {
        rotational = 0;
    }
    tmp = ibps.tmp;
  
    /* ORSA_DEBUG(" --LEAVING-- %x.translational.get(): %x   %x.translational.get(): %x",
       &ibps,ibps.translational.get(),
       this,translational.get());
    */
  
    return (*this);
}

//

orsa::Quaternion RotationalBodyProperty::qDot (const orsa::Quaternion & q,
                                               const orsa::Vector     & omega) {
    return (orsa::Quaternion(omega) * q);
}

orsa::Quaternion RotationalBodyProperty::qDotDot (const orsa::Quaternion & q,
                                                  const orsa::Vector     & omega,
                                                  const orsa::Vector     & omegaDot) {
    return (orsa::Quaternion(omegaDot) * q + orsa::Quaternion(omega) * RotationalBodyProperty::qDot(q,omega));
}

orsa::Vector RotationalBodyProperty::omega (const orsa::Quaternion & q,
                                            const orsa::Quaternion & qDot) {
    return (qDot * inverse(q)).getVector();
}

orsa::Vector RotationalBodyProperty::omegaDot (const orsa::Quaternion & q,
                                               const orsa::Quaternion & qDot,
                                               const orsa::Quaternion & qDotDot) {
    return ((qDotDot - orsa::Quaternion(RotationalBodyProperty::omega(q,qDot)) * qDot) * inverse(q)).getVector();
}

orsa::Quaternion RotationalBodyProperty::qFiniteRotation (const orsa::Quaternion & q,
                                                          const orsa::Vector     & omega,
                                                          const orsa::Time       & dt) {
    const double angle = omega.length() * dt.get_d() / 2;
    return unitQuaternion(Quaternion(cos(angle),sin(angle)*(omega.normalized()))*q);
}

orsa::Vector RotationalBodyProperty::newOmega (const orsa::Vector & omega,
                                               const orsa::Vector & omegaDot,
                                               const orsa::Time   & dt) {
  
    /* 
       ORSA_DEBUG("---prints---");
       print(omegaDot);
       print(dt);
    */
  
    if (omegaDot.length() < orsa::epsilon()) {
        // ORSA_DEBUG("omegaDot is really tiny...");
        return (omega + omegaDot * dt.get_d());
    }
  
    // test
    // return orsa::Vector(omega + omegaDot * dt.get_d());
  
    const double  omegaMagnitude = omega.length();
    const orsa::Vector uOmega          = omega.normalized();

    // const double  omegaDotMagnitude = omegaDot.length();
    const orsa::Vector uOmegaDot          = omegaDot.normalized();
  
    orsa::Vector newOmega;
  
    if (omegaMagnitude > orsa::epsilon()) {
    
        const double radialComponent = omegaDot * uOmega;
    
        // ORSA_DEBUG("radialComponent/omegaDotMagnitude: %g",radialComponent/omegaDotMagnitude);
    
        const orsa::Vector  tangentOmegaDot          = omegaDot - radialComponent * uOmega;
        const double  tangentOmegaDotMagnitude = tangentOmegaDot.length();
        const orsa::Vector uTangentOmegaDot          = tangentOmegaDot.normalized();
    
        const orsa::Vector uRotationAxis = (orsa::externalProduct(uOmega,uTangentOmegaDot)).normalized();
        const double rotationAngle  = tangentOmegaDotMagnitude * dt.get_d() / omegaMagnitude;
    
        // rotate old omega direction around uRotationAxis of angle rotationAngle
        const double     qAngle = rotationAngle / 2; // factor 2 due to the quaternion angle definition
        const orsa::Quaternion qRot   = unitQuaternion(Quaternion(cos(qAngle),sin(qAngle)*uRotationAxis));
    
        newOmega = (radialComponent * dt.get_d() + omegaMagnitude) * (qRot*uOmega*conjugate(qRot)).getVector().normalized();
        //
        /* 
           newOmega = 
           radialComponent * dt.get_d() * uOmega +
           omegaMagnitude * (qRot*uOmega*conjugate(qRot)).getVector().normalized();
        */
    
    } else {
    
        newOmega = omega + omegaDot * dt.get_d();
    }
  
    // ORSA_DEBUG("omega comparison...");
    // print(omega);
    // print(newOmega);
  
    /* 
       ORSA_DEBUG("scalar product: %20.16f",
       newOmega.normalized() * omega.normalized());
    */
  
    return newOmega;  
}

// Body

Body::BodyID Body::bodyID_counter = 0;

Body::Body() : osg::Referenced(true), _id(bodyID_counter++) {
    _init();
    // ORSA_DEBUG("created body id: %i",id());
}

void Body::_init() {
    // _validBodyType = false;
    /* 
       _validBirthTime = _validDeathTime = false;
    */
  
    // _ibps      = 0;
    // _attitude  = 0;
    // _shape     = 0;
    // _multipole = 0;
    // _bpvc      = 0;
    // _bac       = 0;
  
    // comet particle beta
    beta    = 0;
    betaSun = 0;
  
    // light source
    isLightSource = false;
  
    nonInteractingGroup = false;
}

Body::~Body() {
  
}

/* 
   bool Body::setBirthTime(const orsa::Time & t) {
   if (!_validBirthTime) {
   _birthTime = t;
   _validBirthTime = true;
   } else {
   // error...
   }
   return true;
   }
   
   bool Body::validBirthTime() const {
   return _validBirthTime;
   }
   
   const orsa::Time & Body::getBirthTime() const {
   if (!_validBirthTime) {
   // error...
   }
   return _birthTime;
   }
   
   bool Body::setDeathTime(const orsa::Time & t) {
   if (!_validDeathTime) {
   _deathTime = t;
   _validDeathTime = true;
   } else {
   // error...
   }
   return true;
   }
   
   bool Body::validDeathTime() const {
   return _validDeathTime;
   }
   
   const orsa::Time & Body::getDeathTime() const {
   if (!_validDeathTime) {
   // error...
   }
   return _deathTime;
   }
   
   bool Body::isActive(const orsa::Time &) {
   // code needed here!
   return true;
   }
*/

/* 
   bool Body::setMass(const double m) {
   _mass = m;
   return true;
   }
   
   double Body::getMass() const {
   return _mass;
   }
*/

bool Body::setName(const std::string & s) {
    _name = s;
    return true;
}

const std::string & Body::getName() const {
    return _name;
}
