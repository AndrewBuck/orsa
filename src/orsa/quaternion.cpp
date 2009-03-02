#include <orsa/quaternion.h>

#include <orsa/print.h>

using namespace orsa;

Quaternion Quaternion::operator * (const Quaternion & rhs) const {
  /* 
     ORSA_DEBUG("...");
     print(*this);
     print(rhs);
     ORSA_DEBUG("...more...");
     print(_s*rhs._s);
     print(_v*rhs._v);
     print(_s*rhs._v);
     print(rhs._s*_v);
     print(orsa::externalProduct(_v,rhs._v));
  */
  
  return Quaternion(_s*rhs._s - _v*rhs._v,
		    _s*rhs._v + rhs._s*_v + orsa::externalProduct(_v,rhs._v));
}
