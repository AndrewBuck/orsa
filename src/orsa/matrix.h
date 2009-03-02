#ifndef _ORSA_MATRIX_
#define _ORSA_MATRIX_

#include <orsa/cache.h>
#include <orsa/debug.h>
#include <orsa/double.h>
#include <orsa/vector.h>

#include <osg/Matrixd>

namespace orsa {
  
  class Matrix {
  public:
    Matrix();
    Matrix(const Matrix &);
    Matrix(const Double m11, 
	   const Double m12,
	   const Double m13,
	   const Double m21,
	   const Double m22,
	   const Double m23,
	   const Double m31,
	   const Double m32,
	   const Double m33);
    
  public:
    void set(const Double m11, 
	     const Double m12,
	     const Double m13,
	     const Double m21,
	     const Double m22,
	     const Double m23,
	     const Double m31,
	     const Double m32,
	     const Double m33); 

    void get(Double & m11, 
	     Double & m12,
	     Double & m13,
	     Double & m21,
	     Double & m22,
	     Double & m23,
	     Double & m31,
	     Double & m32,
	     Double & m33) const;

  public:
    const Double & getM11() const { return m11; }
    const Double & getM12() const { return m12; }
    const Double & getM13() const { return m13; }
    const Double & getM21() const { return m21; }
    const Double & getM22() const { return m22; }
    const Double & getM23() const { return m23; }
    const Double & getM31() const { return m31; }
    const Double & getM32() const { return m32; }
    const Double & getM33() const { return m33; }
    
  public:
    osg::Matrixd getMatrixd() const;
    
  public:
    Double determinant() const;
 
  public:
    static Matrix identity();
    
  public:
    static bool invert(const Matrix & src, Matrix & inverse);
    static void transpose(const Matrix & src, Matrix & transposed);
    
  public:
    static void OpenGLMatrix(const Matrix & src, Double opengl_matrix[16]);
    
    // operators
    Matrix & operator += (const Matrix &);
    Matrix & operator -= (const Matrix &);
    Matrix & operator *= (const Double &);
    Matrix & operator /= (const Double &);
    
    // sign
    Matrix operator + () const;
    Matrix operator - () const;
    
  public:
    Matrix rotX(const Double & alpha);
    Matrix rotY(const Double & alpha);
    Matrix rotZ(const Double & alpha);
    
  public:
    static Matrix axisRotation(const Vector & axis, const Double & alpha);
    
  public:
    Matrix operator + (const Matrix &) const;
    Matrix operator - (const Matrix &) const;
    Matrix operator * (const Matrix &) const;
    
    /* 
       friend Matrix operator + (const Matrix &, const Matrix &);
       friend Matrix operator - (const Matrix &, const Matrix &);
       friend Matrix operator * (const Matrix &, const Matrix &);
       friend Vector operator * (const Matrix &, const Vector &);
       friend Vector operator * (const Vector &, const Matrix &);
    */
    
    // friend Matrix operator * (const Double  , const Matrix &);
    // friend Matrix operator * (const Matrix &, const Double  );
    // friend Vector operator * (const Matrix &, const Vector &);
    // friend Vector operator * (const Vector &, const Matrix &);
    
  public:
    // util...
    // void print() const;
    
  protected:
    Double m11,m12,m13,m21,m22,m23,m31,m32,m33;
  };
  
  // Matrix operator * (const Double &, const Matrix &);
  // Matrix operator * (const Matrix &, const Double &);
  // Vector operator * (const Matrix &, const Vector &);
  // Vector operator * (const Vector &, const Matrix &);
  
  inline Matrix operator * (const Double & f, const Matrix & m) {
    Matrix q(m);
    q *= f;
    return q;
  }
  
  inline Matrix operator * (const Matrix & m, const Double & f) {
    Matrix q(m);
    q *= f;
    return q;
  }

  inline Vector operator * (const Matrix & m, const Vector & v) {
    return Vector (m.getM11()*v.getX()+m.getM12()*v.getY()+m.getM13()*v.getZ(),
		   m.getM21()*v.getX()+m.getM22()*v.getY()+m.getM23()*v.getZ(),
		   m.getM31()*v.getX()+m.getM32()*v.getY()+m.getM33()*v.getZ());
  }
  
  inline Vector operator * (const Vector & v, const Matrix & m) {
    return Vector (v.getX()*m.getM11()+v.getY()*m.getM21()+v.getZ()*m.getM31(),
		   v.getX()*m.getM12()+v.getY()*m.getM22()+v.getZ()*m.getM32(),
		   v.getX()*m.getM13()+v.getY()*m.getM23()+v.getZ()*m.getM33());
  }
  
  /* 
     Matrix operator + (const Matrix &, const Matrix &);
     Matrix operator - (const Matrix &, const Matrix &);
     Matrix operator * (const Double  , const Matrix &);
     Matrix operator * (const Matrix &, const Double  );
     Matrix operator * (const Matrix &, const Matrix &);
     Vector operator * (const Matrix &, const Vector &);
     Vector operator * (const Vector &, const Matrix &);
  */
  
} // namespace orsa

#endif // _ORSA_MATRIX_
