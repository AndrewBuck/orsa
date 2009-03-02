#include <orsa/reference_frame.h>

using namespace orsa;

// AbsoluteReferenceFrame

ReferenceFrame * AbsoluteReferenceFrame::_instance = 0;

AbsoluteReferenceFrame::AbsoluteReferenceFrame() : _identity(Matrix::identity()) { }
