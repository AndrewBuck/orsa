#ifndef _VESTA_H_
#define _VESTA_H_

// #include <orsa/attitude.h>
#include <orsa/shape.h>
#include <orsa/unit.h>
#include <orsa/debug.h>

#include <string>

using namespace std;
using namespace orsa;

class VestaShape : public orsa::TriShape {
 public:
  VestaShape() : TriShape() {
    
  }
 public:
  bool read(const std::string & filename) {
    FILE * fp_model = fopen(filename.c_str(),"r");
    if (!fp_model) {
      ORSA_ERROR("can't open file %s",filename.c_str());
      return false;
    }
    //
    _vertex.clear();
    _face.clear();
    //
    char line[1024];
    unsigned int nv, nt;
    if (fgets(line,1024,fp_model) != 0) {
      sscanf(line,"%i %i",&nv,&nt);
    }
    double x,y,z;
    unsigned int i,j,k;
    // char type[16];
    string s_type;
    unsigned int count_v=0;
    while (fgets(line,1024,fp_model) != 0) {
      // sscanf(line,"%s",type);
      // s_type = type;
      // if (s_type == "v") {
      ++count_v;
      if (count_v <= nv) {
	if (3 == sscanf(line,"%lf %lf %lf",&x,&y,&z)) {
	  // x = FromUnits(x,Unit::KM);
	  // y = FromUnits(y,Unit::KM);
	  // z = FromUnits(z,Unit::KM);
	  _vertex.push_back(Vector(FromUnits(x,Unit::KM),
				   FromUnits(y,Unit::KM),
				   FromUnits(z,Unit::KM)));
	}
      } else {
	if (3 == sscanf(line,"%i %i %i",&i,&j,&k)) {
	  // this is needed sometimes...
	  // a flag fould be useful
	  // --i; --j; --k; // 1,2,3... to 0,1,2...
	  if ((i != j) && (j != k)) _face.push_back(TriIndex(i,j,k));
	  else {
	    ORSA_ERROR("repeated indexes...");
	    return false;
	  }	
	}
      } 
    }
    // post-read check: vertex numbering starting from 0?
    {
      bool start_from_zero = false;
      unsigned int q=0;
      while (q<_face.size()) {
	const TriIndex t = _face[q];
	if ((t.i() == 0) || (t.j() == 0) || (t.k() == 0)) {
	  start_from_zero = true;
	  break;
	}
	++q;
      }
      if (!start_from_zero) {
	for (unsigned int p=0;p<_face.size();++p) {
	  const TriIndex t = _face[p];
	  _face[p] = TriIndex(t.i()-1,t.j()-1,t.k()-1);
	}
      }
    }
    fclose(fp_model);
    
    return true;
  }	
};

#endif // _VESTA_H_
