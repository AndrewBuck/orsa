#include <orsa/matrix.h>

#include <osg/Matrixd>

using namespace orsa;

/* 
 * 11 12 13
 * 21 22 23
 * 31 32 33
 *
 */

Matrix::Matrix() {
  // m11 = m12 = m13 = m21 = m22 = m23 = m31 = m32 = m33 = 0.0;
}

Matrix::Matrix(const Matrix & m) {
  m11 = m.m11;
  m12 = m.m12;
  m13 = m.m13;
  m21 = m.m21;
  m22 = m.m22;
  m23 = m.m23;
  m31 = m.m31;
  m32 = m.m32;
  m33 = m.m33;
}

Matrix::Matrix(const double & _m11, 
	       const double & _m12,
	       const double & _m13,
	       const double & _m21,
	       const double & _m22,
	       const double & _m23,
	       const double & _m31,
	       const double & _m32,
	       const double & _m33) {
  m11 = _m11;
  m12 = _m12;
  m13 = _m13;
  m21 = _m21;
  m22 = _m22;
  m23 = _m23;
  m31 = _m31;
  m32 = _m32;
  m33 = _m33;
  //
  check();
}

void Matrix::set(const double & _m11, 
		 const double & _m12,
		 const double & _m13,
		 const double & _m21,
		 const double & _m22,
		 const double & _m23,
		 const double & _m31,
		 const double & _m32,
		 const double & _m33) {
  m11 = _m11;
  m12 = _m12;
  m13 = _m13;
  m21 = _m21;
  m22 = _m22;
  m23 = _m23;
  m31 = _m31;
  m32 = _m32;
  m33 = _m33;
  //
  check();
}

void Matrix::get(double & _m11, 
		 double & _m12,
		 double & _m13,
		 double & _m21,
		 double & _m22,
		 double & _m23,
		 double & _m31,
		 double & _m32,
		 double & _m33) const {
  _m11 = m11;
  _m12 = m12;
  _m13 = m13;
  _m21 = m21;
  _m22 = m22;
  _m23 = m23;
  _m31 = m31;
  _m32 = m32;
  _m33 = m33;
}

osg::Matrixd Matrix::getMatrixd() const {
  return osg::Matrixd(getM11(),getM12(),getM13(),0,
		      getM21(),getM22(),getM23(),0,
		      getM31(),getM32(),getM33(),0,
		      0,       0,       0,       1);
}

double Matrix::determinant() const {
  return (m11*(m22*m33-m23*m32)-m12*(m21*m33-m31*m23)+m13*(m21*m32-m31*m22));
}

Matrix Matrix::identity() {
  Matrix M;
  M.m11 = M.m22 = M.m33 = 1;
  M.m12 = M.m13 = M.m21 = M.m23 = M.m31 = M.m32 = 0;
  return M;
}

bool Matrix::invert(const Matrix & src, Matrix & inverse) {
  
  const unsigned int m_size = 3;
  
  double t;
  unsigned int i, j, k, swap;
  double tmp[m_size][m_size];
  double local_inverse[m_size*m_size];
  
  // inverse = identity
  for (i=0;i<m_size;++i) {
    for (j=0;j<m_size;++j) {
      if (i==j) {
	local_inverse[i*m_size+j] = 1;
      } else {
	local_inverse[i*m_size+j] = 0;
      }
    }
  }
  
  /* 
     for (i = 0; i < m_size; i++) {
     for (j = 0; j < m_size; j++) {
     tmp[i][j] = src[i*m_size+j];
     }
     }
  */
  tmp[0][0] = src.m11;
  tmp[0][1] = src.m12;
  tmp[0][2] = src.m13;
  tmp[1][0] = src.m21;
  tmp[1][1] = src.m22;
  tmp[1][2] = src.m23;
  tmp[2][0] = src.m31;
  tmp[2][1] = src.m32;
  tmp[2][2] = src.m33;
  
  for (i = 0; i < m_size; i++) {
    /* look for largest element in column. */
    swap = i;
    for (j = i + 1; j < m_size; j++) {
      if (fabs(tmp[j][i]) > fabs(tmp[i][i])) {
	swap = j;
      }
    }
      
    if (swap != i) {
      /* swap rows. */
      for (k = 0; k < m_size; k++) {
	t = tmp[i][k];
	tmp[i][k] = tmp[swap][k];
	tmp[swap][k] = t;
	  
	t = local_inverse[i*m_size+k];
	local_inverse[i*m_size+k] = local_inverse[swap*m_size+k];
	local_inverse[swap*m_size+k] = t;
      }
    }
      
    if (tmp[i][i] == 0) {
      /* no non-zero pivot.  the matrix is singular, which
	 shouldn't happen.  This means the user gave us a bad
	 matrix. */
      ORSA_ERROR("sorry, cannot invert this matrix");
      // src.print();
      return false;
    }
      
    t = tmp[i][i];
    for (k = 0; k < m_size; k++) {
      tmp[i][k] /= t;
      local_inverse[i*m_size+k] /= t;
    }
    for (j = 0; j < m_size; j++) {
      if (j != i) {
	t = tmp[j][i];
	for (k = 0; k < m_size; k++) {
	  tmp[j][k] -= tmp[i][k]*t;
	  local_inverse[j*m_size+k] -= local_inverse[i*m_size+k]*t;
	}
      }
    }
  }
    
  inverse.set(local_inverse[0],
	      local_inverse[1],
	      local_inverse[2],
	      local_inverse[3],
	      local_inverse[4],
	      local_inverse[5],
	      local_inverse[6],
	      local_inverse[7],
	      local_inverse[8]);
    
  return true;
}
  
void Matrix::transpose(const Matrix & src, Matrix & transposed) {
  transposed.m11 = src.m11;
  transposed.m22 = src.m22;
  transposed.m33 = src.m33;
  transposed.m12 = src.m21;
  transposed.m21 = src.m12;
  transposed.m13 = src.m31;
  transposed.m31 = src.m13;
  transposed.m23 = src.m32;
  transposed.m32 = src.m23;
}
  
void Matrix::OpenGLMatrix(const Matrix & src, double opengl_matrix[16]) {
  opengl_matrix[0]  = src.m11;
  opengl_matrix[1]  = src.m21;
  opengl_matrix[2]  = src.m31;
  opengl_matrix[3]  = 0;
  opengl_matrix[4]  = src.m12;
  opengl_matrix[5]  = src.m22;
  opengl_matrix[6]  = src.m32;
  opengl_matrix[7]  = 0;
  opengl_matrix[8]  = src.m13;
  opengl_matrix[9]  = src.m23;
  opengl_matrix[10] = src.m33;
  opengl_matrix[11] = 0;
  opengl_matrix[12] = 0;
  opengl_matrix[13] = 0;
  opengl_matrix[14] = 0;
  opengl_matrix[15] = 1;
  //
  /* { 
     ORSA_DEBUG("Matrix::OpenGLMatrix() src.determinant(): %g",
     src.determinant());
     src.print();
     }
  */
}

Matrix & Matrix::operator += (const Matrix & m) {
  m11 += m.m11;
  m12 += m.m12;
  m13 += m.m13;
  m21 += m.m21;
  m22 += m.m22;
  m23 += m.m23;
  m31 += m.m31;
  m32 += m.m32;
  m33 += m.m33;
  return (*this);
}

Matrix & Matrix::operator -= (const Matrix & m) {
  m11 -= m.m11;
  m12 -= m.m12;
  m13 -= m.m13;
  m21 -= m.m21;
  m22 -= m.m22;
  m23 -= m.m23;
  m31 -= m.m31;
  m32 -= m.m32;
  m33 -= m.m33;
  return (*this);
}

Matrix & Matrix::operator *= (const double & f) {
  m11 *= f;
  m12 *= f;
  m13 *= f;
  m21 *= f;
  m22 *= f;
  m23 *= f;
  m31 *= f;
  m32 *= f;
  m33 *= f;
  return (*this);
}

Matrix & Matrix::operator /= (const double & f) {
  m11 /= f;
  m12 /= f;
  m13 /= f;
  m21 /= f;
  m22 /= f;
  m23 /= f;
  m31 /= f;
  m32 /= f;
  m33 /= f;
  return (*this);
}

Matrix Matrix::operator + () const {
  return Matrix(*this);
}

Matrix Matrix::operator - () const {
  Matrix m;
  m.m11 = -m11;
  m.m12 = -m12;
  m.m13 = -m13;
  m.m21 = -m21;
  m.m22 = -m22;
  m.m23 = -m23;
  m.m31 = -m31;
  m.m32 = -m32;
  m.m33 = -m33;
  return m;
}

// from the OpenGL Red Book, third edition, Appendix F, page 672
Matrix Matrix::axisRotation(const Vector & v, const double & angle) {
  
  // NOTE: check if angle sign agrees with other rotations...
  
  const Vector u = v.normalized();
  //
  Matrix S;
  S.m11 = S.m22 = S.m33 = 0;
  S.m21 =  u.getZ();
  S.m12 = -u.getZ();
  S.m13 =  u.getY();
  S.m31 = -u.getY();
  S.m32 =  u.getX();
  S.m23 = -u.getX();
  //
  Matrix P;
  P.m11 = u.getX()*u.getX();
  P.m22 = u.getY()*u.getY();
  P.m33 = u.getZ()*u.getZ();
  P.m12 = P.m21 = u.getX()*u.getY();
  P.m13 = P.m31 = u.getX()*u.getZ();
  P.m23 = P.m32 = u.getY()*u.getZ();
  //
  double s,c;
  orsa::sincos(angle,&s,&c);
  //
  const Matrix M = P + c * (Matrix::identity() - P) + s * S;
  return M;
}

Matrix Matrix::rotX(const double & angle) {
  double s,c;
  orsa::sincos(angle,&s,&c);
  Matrix rot;
  rot.m12 = rot.m13 = rot.m21 = rot.m31 = 0;
  rot.m11 = 1;
  rot.m22 = rot.m33 = c;
  rot.m23 = -s;
  rot.m32 =  s;
  //
  /* Matrix u = (*this);
     Matrix p = rot*u;
     (*this) = p;
  */
  //
  (*this) = rot*(*this);
  //
  return (*this);
}

Matrix Matrix::rotY(const double & angle) {
  double s,c;
  orsa::sincos(angle,&s,&c);
  Matrix rot;
  rot.m12 = rot.m21 = rot.m23 = rot.m32 = 0;
  rot.m22 = 1;
  rot.m11 = rot.m33 = c;
  rot.m13 =  s;
  rot.m31 = -s;
  //
  /* Matrix u = (*this);
     Matrix p = rot*u;
     (*this) = p;
  */
  //
  (*this) = rot*(*this);
  //
  return (*this);
}

Matrix Matrix::rotZ(const double & angle) {
  double s,c;
  orsa::sincos(angle,&s,&c);
  Matrix rot;
  rot.m13 = rot.m23 = rot.m31 = rot.m32 = 0;
  rot.m33 = 1;
  rot.m11 = rot.m22 = c;
  rot.m12 = -s;
  rot.m21 =  s;
  //
  /* Matrix u = (*this);
     Matrix p = rot*u;
     (*this) = p;
  */
  //
  (*this) = rot*(*this);
  //
  return (*this);
}

/* 
   Matrix operator + (const Matrix & p, const Matrix & q) {
   return Matrix(p.m11+q.m11,
   p.m12+q.m12,
   p.m13+q.m13,
   p.m21+q.m21,
   p.m22+q.m22,
   p.m23+q.m23,
   p.m31+q.m31,
   p.m32+q.m32,
   p.m33+q.m33);
   }
   
   Matrix operator - (const Matrix & p, const Matrix & q) {
   return Matrix(p+(-q));
   }
   
   Matrix operator * (const double f, const Matrix & m) {
   Matrix q(m);
   q *= f;
   return q;
   }
   
   Matrix operator * (const Matrix & m, const double f) {
   Matrix q(m);
   q *= f;
   return q;
   }
   
   Matrix operator * (const Matrix & p, const Matrix & q) {
   return Matrix(p.m11*q.m11+p.m12*q.m21+p.m13*q.m31,
   p.m11*q.m12+p.m12*q.m22+p.m13*q.m32,
   p.m11*q.m13+p.m12*q.m23+p.m13*q.m33,
   p.m21*q.m11+p.m22*q.m21+p.m23*q.m31,
   p.m21*q.m12+p.m22*q.m22+p.m23*q.m32,
   p.m21*q.m13+p.m22*q.m23+p.m23*q.m33,
   p.m31*q.m11+p.m32*q.m21+p.m33*q.m31,
   p.m31*q.m12+p.m32*q.m22+p.m33*q.m32,
   p.m31*q.m13+p.m32*q.m23+p.m33*q.m33);
   }
   
   Vector operator * (const Matrix & m, const Vector & v) {
   return Vector (m.m11*v.getX()+m.m12*v.getY()+m.m13*v.getZ(),
   m.m21*v.getX()+m.m22*v.getY()+m.m23*v.getZ(),
   m.m31*v.getX()+m.m32*v.getY()+m.m33*v.getZ());
   }
   
   Vector operator * (const Vector & v, const Matrix & m) {
   return Vector (v.getX()*m.m11+v.getY()*m.m21+v.getZ()*m.m31,
   v.getX()*m.m12+v.getY()*m.m22+v.getZ()*m.m32,
   v.getX()*m.m13+v.getY()*m.m23+v.getZ()*m.m33);
   }
*/

/* 
   void Matrix::print() const {
   ORSA_DEBUG("Matrix::print():");
   ORSA_DEBUG("%g %g %g\n",m11,m12,m13);
   ORSA_DEBUG("%g %g %g\n",m21,m22,m23);
   ORSA_DEBUG("%g %g %g\n",m31,m32,m33);
   }
*/

/* 
   Matrix operator + (const Matrix & rhs) const {
   return Matrix(m11+rhs.m11,
   m12+rhs.m12,
   m13+rhs.m13,
   m21+rhs.m21,
   m22+rhs.m22,
   m23+rhs.m23,
   m31+rhs.m31,
   m32+rhs.m32,
   m33+rhs.m33);
   }
*/

Matrix Matrix::operator + (const Matrix & rhs) const {
  Matrix m = (*this);
  m += rhs;
  return m;
}

Matrix Matrix::operator - (const Matrix & rhs) const {
  Matrix m = (*this);
  m -= rhs;
  return m;
}

/* 
   Matrix operator * (const double & f, const Matrix & m) {
   Matrix q(m);
   q *= f;
   return q;
   }
   
   Matrix operator * (const Matrix & m, const double & f) {
   Matrix q(m);
   q *= f;
   return q;
   }
*/

Matrix Matrix::operator * (const Matrix & rhs) const {
  return Matrix(m11*rhs.m11+m12*rhs.m21+m13*rhs.m31,
		m11*rhs.m12+m12*rhs.m22+m13*rhs.m32,
		m11*rhs.m13+m12*rhs.m23+m13*rhs.m33,
		m21*rhs.m11+m22*rhs.m21+m23*rhs.m31,
		m21*rhs.m12+m22*rhs.m22+m23*rhs.m32,
		m21*rhs.m13+m22*rhs.m23+m23*rhs.m33,
		m31*rhs.m11+m32*rhs.m21+m33*rhs.m31,
		m31*rhs.m12+m32*rhs.m22+m33*rhs.m32,
		m31*rhs.m13+m32*rhs.m23+m33*rhs.m33);
}

void Matrix::check() const {
  orsa::check(m11);
  orsa::check(m12);
  orsa::check(m13);
  orsa::check(m21);
  orsa::check(m22);
  orsa::check(m23);
  orsa::check(m31);
  orsa::check(m32);
  orsa::check(m33); 
}

bool orsa::operator == (const orsa::Matrix & a, const orsa::Matrix & b) {
  if (a.getM11() != b.getM11()) return false;
  if (a.getM12() != b.getM12()) return false;
  if (a.getM13() != b.getM13()) return false;
  if (a.getM21() != b.getM21()) return false;
  if (a.getM22() != b.getM22()) return false;
  if (a.getM23() != b.getM23()) return false;
  if (a.getM31() != b.getM31()) return false;
  if (a.getM32() != b.getM32()) return false;
  if (a.getM33() != b.getM33()) return false;
  return true;
}
