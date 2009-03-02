#ifndef __SURFACE_COVERAGE_FACECOLOR_H__
#define __SURFACE_COVERAGE_FACECOLOR_H__

#include <Qt>

#include "mainThread.h"

Qt::GlobalColor FaceColor (MainThread         * mainThread,
			   const unsigned int   face,
			   const orsa::Time   & t);

bool CumulativeColorMatch (const Qt::GlobalColor faceColor,
			   const Qt::GlobalColor plotColor);

#endif // __SURFACE_COVERAGE_FACECOLOR_H__
