#ifndef _ORSA_OSG_ORBIT_NODE_CALLBACK_H_
#define _ORSA_OSG_ORBIT_NODE_CALLBACK_H_

#include <osg/MatrixTransform>
#include <osg/NodeCallback>
#include <osg/NodeVisitor>

#include <orsaOSG/AnimationTime.h>

#include <orsa/bodygroup.h>
#include <orsa/orbit.h>
#include <orsa/unit.h>

namespace orsaOSG {
  
    class OrbitNodeUpdateCallback : public osg::NodeCallback {
    public:
        OrbitNodeUpdateCallback(orsa::BodyGroup  * bg,
                                const orsa::Body * b,
                                orsaOSG::AnimationTime * at) : 
            osg::NodeCallback(),
            _bg(bg),
            _b(b),
            _at(at) { 
            // _orbitProxy = new orsa::OrbitProxy(b,bg);
            _orbitProxy = new orsa::OrbitProxy(b,bg,0.001,orsa::Time(0,1,0,0,0));
        }
    protected:
        ~OrbitNodeUpdateCallback() { }
    
    public:
        void operator () (osg::Node * node, osg::NodeVisitor * nv) {
      
            /* 
               ORSA_DEBUG("body [%s]",
               _b->getName().c_str());
            */
      
            if (nv->getVisitorType()==osg::NodeVisitor::UPDATE_VISITOR && nv->getFrameStamp()) {
	
                osg::MatrixTransform * mt = dynamic_cast<osg::MatrixTransform *> (node);
                //
                if (mt) {
	  
                    bool needToClearOrbit=true;
	  
                    const orsa::Time simulationTime = _at->getSimulationTime(nv->getFrameStamp()->getFrameNumber());
	  
                    if (_b->alive(simulationTime)) {
	    
                        if ( (!(_lastSimulationTime.isSet())) ||
                             ((_lastSimulationTime.isSet()) && (simulationTime != _lastSimulationTime.getRef())) ) {
	      
                            _lastSimulationTime = simulationTime;
	      
                            orsa::Orbit orbit;
                            const bool goodOrbit = _orbitProxy->getOrbit(orbit,
                                                                         simulationTime);
	      
                            if (goodOrbit) {
		
                                if (orbit.e < 1) {
		  
                                    needToClearOrbit=false;
		  
                                    /* 
                                       ORSA_DEBUG("a: %f KM   e: %f",
                                       FromUnits(orbit.a,orsa::Unit::KM,-1)(),
                                       orbit.e());
                                    */
		  
                                    orsa::Matrix m = orsa::Matrix::identity();
                                    // rot
                                    m.rotZ(-orbit.omega_node);
                                    m.rotX(-orbit.i);
                                    m.rotZ(-orbit.omega_pericenter);
                                    // scale
                                    const orsa::Matrix m_scale(orbit.a, 0, 0,
                                                               0, orbit.a*sqrt(1-orbit.e*orbit.e),0,
                                                               0, 0, orbit.a);	
                                    //
                                    /* const orsa::Double orbit_a_AU = orsa::FromUnits(orbit.a,orsa::Unit::AU,-1);
                                       const orsa::Matrix m_scale(orbit_a_AU, 0, 0,
                                       0, orbit_a_AU*sqrt(orsa::one()-orbit.e*orbit.e),0,
                                       0, 0, orbit_a_AU);
                                    */
                                    m = m_scale * m;
                                    // ecc
                                    m.rotZ(-orbit.eccentricAnomaly());
		  
                                    mt->setMatrix(m.getMatrixd());
		  
                                } else {
		  
                                    needToClearOrbit=true;
                                } 		
		
                                /* 
                                   } else {
                                   needToClearOrbit=true;
                                   } 	      
                                   } else {
                                   needToClearOrbit=true;
                                   }
                                */
		
                            }
	      
                        }
                    }
	  
                    if (needToClearOrbit) {
                        orsa::Matrix m = 0.0*orsa::Matrix::identity();
                        mt->setMatrix(m.getMatrixd());
                    }

                    // ORSA_DEBUG("needToClearOrbit: %i",needToClearOrbit);
                }
            }
      
            // must call any nested node callbacks and continue subgraph traversal.
            NodeCallback::traverse(node,nv);
        }
    
    protected:    
        osg::ref_ptr<orsa::BodyGroup> _bg;
        osg::ref_ptr<const orsa::Body> _b;
    protected:
        osg::ref_ptr<AnimationTime> _at;
    protected:
        orsa::Cache<orsa::Time> _lastSimulationTime;
    protected:
        osg::ref_ptr<orsa::OrbitProxy> _orbitProxy;
    };
  
}; // namespace orsaOSG

#endif // _ORSA_OSG_ORBIT_NODE_CALLBACK_H_
