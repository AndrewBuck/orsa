#ifndef _ORSA_PRINT_
#define _ORSA_PRINT_

#include <orsa/datetime.h>
#include <orsa/matrix.h>
#include <orsa/quaternion.h>

namespace orsa {
  
  inline void print(const orsa::Double & d) {
    ORSA_DEBUG("d: %+20.12Fg",d.get_mpf_t());
  }
  
  inline void print(const orsa::Vector & v) {
    ORSA_DEBUG("v: length: %.12Fe [%+.12Fe,%+.12Fe,%+.12Fe]",
	       v.length().get_mpf_t(),
	       v.getX().get_mpf_t(),
	       v.getY().get_mpf_t(),
	       v.getZ().get_mpf_t());
  }
  
  inline void print(const orsa::Matrix & m) {
    ORSA_DEBUG("m: det(m) = %.12Fe\n"
	       "[%+.12Fe,%+.12Fe,%+.12Fe]\n"
	       "[%+.12Fe,%+.12Fe,%+.12Fe]\n"
	       "[%+.12Fe,%+.12Fe,%+.12Fe]",
	       m.determinant().get_mpf_t(),
	       m.getM11().get_mpf_t(),  m.getM12().get_mpf_t(),  m.getM13().get_mpf_t(), 
	       m.getM21().get_mpf_t(),  m.getM22().get_mpf_t(),  m.getM23().get_mpf_t(), 
	       m.getM31().get_mpf_t(),  m.getM32().get_mpf_t(),  m.getM33().get_mpf_t());  
  }
  
  inline void print(const orsa::Quaternion & q) {
    ORSA_DEBUG("q: s: %+9.6Fe   v: %9.6Fe x [%6.3Ff,%6.3Ff,%6.3Ff]  l: %9.6Fe",
	       q.getScalar().get_mpf_t(),
	       q.getVector().length().get_mpf_t(),
	       q.getVector().normalized().getX().get_mpf_t(),
	       q.getVector().normalized().getY().get_mpf_t(),
	       q.getVector().normalized().getZ().get_mpf_t(),
	       q.length().get_mpf_t());
  }
  
  inline void print(const orsa::Time & t) {
    ORSA_DEBUG("t: %20.14Ff [day]",orsa::FromUnits(t.asDouble(),orsa::Unit::DAY,-1).get_mpf_t());
  }
  
} // namespace orsa

#endif // _ORSA_PRINT_
