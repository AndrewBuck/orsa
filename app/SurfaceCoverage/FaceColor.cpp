#include "FaceColor.h"

Qt::GlobalColor FaceColor (MainThread         * mainThread,
			   const unsigned int   face,
			   const orsa::Time   & t) {
  Qt::GlobalColor c;
  
  unsigned int goodAllFace=0;
  //
  for (unsigned int k=0; k<mainThread->allFaceAtOnce_out[face].size(); ++k) {
    if (mainThread->allFaceAtOnce_out[face][k] <= t) {
      ++goodAllFace;
    }
  }
  //
  if (goodAllFace > 0) {
    if (goodAllFace == 1) {
      c = Qt::green;
    } else {
      c = Qt::darkGreen;
    }
  } else { 
    unsigned int coverage=0;
    const orsa::TriShape::TriIndex tri = mainThread->faceVector_out[face];  
    
    if (coverage == 0) {
      for (unsigned int k=0; k<mainThread->vertexCoverage_out[tri.i()].size(); ++k) {
	if (mainThread->vertexCoverage_out[tri.i()][k] <= t) { 
	  ++coverage;
	  // must break, to count each vertex at most once
	  break;
	}
      }
    }
    
    if (coverage == 0) {
      for (unsigned int k=0; k<mainThread->vertexCoverage_out[tri.j()].size(); ++k) {
	if (mainThread->vertexCoverage_out[tri.j()][k] <= t) { 
	  ++coverage;
	  // must break, to count each vertex at most once
	  break;
	}
      }
    }
    
    if (coverage == 0) {
      for (unsigned int k=0; k<mainThread->vertexCoverage_out[tri.k()].size(); ++k) {
	if (mainThread->vertexCoverage_out[tri.k()][k] <= t) { 
	  ++coverage;
	  // must break, to count each vertex at most once
	  break;
	}
      }
    }
    
    if (coverage > 0) {
      c = Qt::yellow;
    } else {
      c = Qt::black;
    }
  }
  
  // ORSA_DEBUG("color: %i",c);
  
  return c;
}

bool CumulativeColorMatch (const Qt::GlobalColor faceColor,
			   const Qt::GlobalColor plotColor) {
  
  if (faceColor == plotColor) return true;
  
  switch (plotColor) {
  case Qt::yellow:
    if (faceColor == Qt::green)     return true;
    if (faceColor == Qt::darkGreen) return true;
    break;
  case Qt::green:
    if (faceColor == Qt::darkGreen) return true;
    break;
  case Qt::darkGreen:
    break;
  default:
    return false;
  }
  
  return false;
}
