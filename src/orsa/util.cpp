#include <orsa/util.h>
// #include <orsa/print.h>
#include <orsa/statistic.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

using namespace orsa;

std::string & orsa::removeLeadingAndTrailingSpaces(std::string & s) {
  const int first = s.find_first_not_of(" ");
  s.erase(0,first);
  const int last  = s.find_last_not_of(" ");
  s.erase(last+1,s.size());
  return s;
}

std::string & orsa::removeAllSpaces(std::string & s) {
  s.erase(remove(s.begin(),s.end(),' '),s.end());
  return s;
}

std::string & orsa::removeLeadingPlusSign(std::string & s) {
  const int first = s.find_first_not_of("+");
  s.erase(0,first);
  return s;
}


/* This code was used in orsa/interaction.cpp at some point
   to determine on the fly a good value for dt. 
   The only unsolved problem is that of using non-tmp IBPS entries,
   as temporary entries are just wrong outside the integrator.
   
   if (ref_b->propulsion.get()) {
   const orsa::Vector propulsionThrust = ref_b->propulsion->getThrust(t);
   if (propulsionThrust.lengthSquared() > orsa::epsilon()*orsa::epsilon()) {
   const orsa::BodyGroup::BodyInterval * bi =
   bg->getBodyInterval(ref_b);
   const orsa::Time dt = 
   std::min(orsa::Time(0,0,0,5,0),
   std::min(bi->max().time.getRef()-t,
   t-bi->min().time.getRef()));
   // REMEMBER: now thrust is in FrenetSerret components
   orsa::Vector T, N, B;
   ORSA_DEBUG("bi->min().time.getRef().tmp: %i",bi->min().tmp);
   ORSA_DEBUG("bi->max().time.getRef().tmp: %i",bi->max().tmp);
   if (t > bi->min().time.getRef()) {
   FrenetSerret(ref_b, bg,
   t,
   -dt,
   T, N, B);
   } else if (t < bi->max().time.getRef()) {
   FrenetSerret(ref_b, bg,
   t,
   dt,
   T, N, B);
   } else {
   ORSA_DEBUG("interval smaller than dt");
   ORSA_DEBUG("--BODY-INTERVAL-TIME-RANGE--");
   print(bi->min().time.getRef());
   print(bi->max().time.getRef());
   ORSA_DEBUG("call time:");
   print(t);
   //
   T = orsa::Vector(1,0,0);
   N = orsa::Vector(0,1,0);
   B = orsa::Vector(0,0,1);
   }
   thrust = 
   propulsionThrust.getX()*T +	
   propulsionThrust.getY()*N +
   propulsionThrust.getZ()*B;
   }
   }
*/

bool orsa::FrenetSerret(const orsa::Body * b,
			orsa::BodyGroup  * bg,
			const orsa::Time & t,
			const orsa::Time & dt,
			orsa::Vector & T,
			orsa::Vector & N,
			orsa::Vector & B) {
  if (dt.getMuSec() == 0) return false;
  orsa::Vector r_t,   v_t;
  orsa::Vector r_tdt, v_tdt;
  if (bg->getInterpolatedPosVel(r_t,  v_t,  b,t) && 
      bg->getInterpolatedPosVel(r_tdt,v_tdt,b,t+dt)) {
    T = v_t.normalized();
    // const orsa::Vector dv = v_tdt-v_t;
    // should chech that velocity is not constant...
    N = ((v_tdt-v_t)/dt.get_d()).normalized();
    B = externalProduct(T,N).normalized();
    //
    /* 
       ORSA_DEBUG("b: [%s]",b->getName().c_str());
       print(t);
       print(dt);
       print(T);
       print(N);
       print(B);
    */
    //
    return true;
  } else {
    // ORSA_DEBUG("problems...");
    return false;
  }
}

bool orsa::eulerAnglesToMatrix(orsa::Matrix       & m,
			       const double & psi,
			       const double & theta,
			       const double & phi) {
  
  // ORSA_DEBUG("this code needs to be verified!!");
  // #warning "this code needs to be verified!!"
  
  m = orsa::Matrix::identity();
  
  m.rotZ(phi);
  
  m.rotX(theta);
  
  m.rotZ(psi);
  
  return true;
}	

bool orsa::matrixToEulerAngles(double       & psi,
			       double       & theta,
			       double       & phi,
			       const orsa::Matrix & m) {
  
  // const orsa::Matrix l2g = localToGlobal(t);
  
  // const double sinTheta = sqrt(1-int_pow(m.getM33(),2));
  const double cosTheta = m.getM33();
  
  // if (fabs(sinTheta) > epsilon()) {
  if ((1-cosTheta*cosTheta) > (epsilon()*epsilon())) {
    
    psi   = atan2(-m.getM13(),
		   m.getM23());
    
    phi   = atan2(-m.getM31(), 
		  -m.getM32());
    
    // should check this better...
    if (fabs(m.getM23()) > epsilon()) {
      theta = atan2(-m.getM23()/cos(psi),
		     m.getM33());
    } else if (fabs(m.getM31()) > epsilon()) {
      theta = atan2(m.getM31()/sin(phi),
		    m.getM33());
    } else {
      // we should not be here...
      ORSA_DEBUG("problems...");
      theta = 0;
    }
    
  } else {
    
    // ORSA_DEBUG("using singular code...");
    
    psi = theta = 0;
    
    phi = atan2(m.getM21(), 
		m.getM11());
    
  }
  
  if (0) {
    
    ORSA_DEBUG("psi..: %f",  psi);
    ORSA_DEBUG("theta: %f",theta);
    ORSA_DEBUG("phi..: %f",  phi);
    
    orsa::Matrix m = orsa::Matrix::identity();
    //
    m.rotZ(phi);
    m.rotX(theta);
    m.rotZ(psi);
    
    const orsa::Matrix dm = m-m;
    
    ORSA_DEBUG("m: %f %f %f %f %f %f %f %f %f",
	       m.getM11(),
	       m.getM12(),
	       m.getM13(),
	       m.getM21(),
	       m.getM22(),
	       m.getM23(),
	       m.getM31(),
	       m.getM32(),
	       m.getM33());
    
    ORSA_DEBUG("m: %f %f %f %f %f %f %f %f %f",
	       m.getM11(),
	       m.getM12(),
	       m.getM13(),
	       m.getM21(),
	       m.getM22(),
	       m.getM23(),
	       m.getM31(),
	       m.getM32(),
	       m.getM33());
    
    ORSA_DEBUG("dm: %f %f %f %f %f %f %f %f %f",
	       dm.getM11(),
	       dm.getM12(),
	       dm.getM13(),
	       dm.getM21(),
	       dm.getM22(),
	       dm.getM23(),
	       dm.getM31(),
	       dm.getM32(),
	       dm.getM33());
  }
  
  return true;
}	


orsa::Matrix orsa::QuaternionToMatrix (const orsa::Quaternion & q) {
  
  const double s = q.getScalar();
  const orsa::Vector v = q.getVector();
  
  const double q0 = s;
  const double q1 = v.getX();
  const double q2 = v.getY();
  const double q3 = v.getZ();
  
  // Eq. (7.7) book by J.B. Kuipers on Quaternions
  return orsa::Matrix(2*q0*q0-1+2*q1*q1,   2*q1*q2-2*q0*q3, 2*q1*q3+2*q0*q2,
		      2*q1*q2+2*q0*q3,   2*q0*q0-1+2*q2*q2, 2*q2*q3-2*q0*q1,
		      2*q1*q3-2*q0*q2,     2*q2*q3+2*q0*q1, 2*q0*q0-1+2*q3*q3);
}

orsa::Quaternion orsa::MatrixToQuaternion (const orsa::Matrix & m) {
  //! assumes that m is a rotation matrix, with det(m)=1;
  
  // From Sec. 7.9 of Quaternions book bu J.B. Kuipers
  
  const double q0 = sqrt(m.getM11()+m.getM22()+m.getM33()+1)/2;
  
  const double q1 = (m.getM23()-m.getM32())/(4*q0);
  const double q2 = (m.getM31()-m.getM13())/(4*q0);
  const double q3 = (m.getM12()-m.getM21())/(4*q0);
  
  return orsa::Quaternion(q0, orsa::Vector(q1,q2,q3));
}

orsa::Matrix orsa::localToGlobal(const orsa::Body       * b,
				 const orsa::BodyGroup  * bg,
				 const orsa::Time       & t) {
  
  if ((bg==0) || (b==0)) {
    ORSA_ERROR("invalid pointer...");
    return orsa::Matrix::identity();
  }
  
  orsa::IBPS ibps;
  if (bg->getInterpolatedIBPS(ibps, b, t)) { 
    if (ibps.rotational.get()) {
      const orsa::Matrix m = QuaternionToMatrix(ibps.rotational.get()->getQ());
      orsa::Matrix m_tr;
      orsa::Matrix::transpose(m, m_tr);
      return m_tr;
    } else {
      // ORSA_DEBUG("problems... body: [%s]",b->getName().c_str());
      return orsa::Matrix::identity();
    }
  } else {
    ORSA_DEBUG("problems...body: [%s]",b->getName().c_str());
    return orsa::Matrix::identity();
  }
}

orsa::Matrix orsa::globalToLocal(const orsa::Body       * b,
				 const orsa::BodyGroup  * bg,
				 const orsa::Time       & t) {
  
  if ((bg==0) || (b==0)) {
    ORSA_ERROR("invalid pointer...");
    return orsa::Matrix::identity();
  }
  
  orsa::IBPS ibps;
  if (bg->getInterpolatedIBPS(ibps, b, t)) { 
    if (ibps.rotational.get()) {
      return QuaternionToMatrix(ibps.rotational.get()->getQ());
    } else {
      // ORSA_DEBUG("problems... body: [%s]",b->getName().c_str());
      return orsa::Matrix::identity();
    }
  } else {
    ORSA_DEBUG("problems... body: [%s]",b->getName().c_str());
    return orsa::Matrix::identity();
  }
}

double orsa::asteroidDiameter(const double & p, 
			      const double & H) {
  return orsa::FromUnits(1329,orsa::Unit::KM)*pow(10,-0.2*H)/sqrt(p);
}

// ConstantZRotation

/* 
   bool ConstantZRotation::update(const orsa::Time & t) {
   const double _phi = _phi0 + _omega*(t-_t0).get_d();
   orsa::Matrix localMatrix = Matrix::identity();
   // same rotation as Attitude::localToGlobal(t)
   localMatrix.rotZ(_phi);
   _m = localMatrix;
   return true;
   }
*/


void orsa::principalAxis(orsa::Matrix & genericToPrincipal,
			 orsa::Matrix & principalInertiaMatrix,
			 const orsa::Matrix & inertiaMatrix) {
  
  const orsa::Matrix & I = inertiaMatrix;
  
  double data[] = { 
    I.getM11(), I.getM12(), I.getM13(), 
    I.getM21(), I.getM22(), I.getM23(), 
    I.getM31(), I.getM32(), I.getM33()};
  
  gsl_matrix_view m = gsl_matrix_view_array (data, 3, 3);
  
  gsl_vector * eval = gsl_vector_alloc (3);
  
  gsl_matrix * evec = gsl_matrix_alloc (3, 3);
  
  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (3);
  
  gsl_eigen_symmv (&m.matrix, eval, evec, w);
  
  gsl_eigen_symmv_free (w);
  
  gsl_eigen_symmv_sort (eval, evec, 
			GSL_EIGEN_SORT_ABS_ASC);
  
  principalInertiaMatrix.set(gsl_vector_get(eval,0),0,0,
			     0,gsl_vector_get(eval,1),0,
			     0,0,gsl_vector_get(eval,2));
  
  genericToPrincipal.set(gsl_matrix_get(evec,0,0),gsl_matrix_get(evec,0,1),gsl_matrix_get(evec,0,2),
			 gsl_matrix_get(evec,1,0),gsl_matrix_get(evec,1,1),gsl_matrix_get(evec,1,2),
			 gsl_matrix_get(evec,2,0),gsl_matrix_get(evec,2,1),gsl_matrix_get(evec,2,2));
  
  // force det=1
  // genericToPrincipal /= genericToPrincipal.determinant();
  
  gsl_vector_free (eval);
  gsl_matrix_free (evec);
}

/***/

RandomPointsInShape::RandomPointsInShape(const orsa::Shape * s,
					 const unsigned int N,
					 const int rs) : osg::Referenced(), shape(s), size(N), randomSeed(rs) {
  rng = new orsa::RNG(randomSeed);
  in.clear();
  in.reserve(N);
  const Box boundingBox = shape->boundingBox();
  counter = 0;
  while (counter < N) {
    in[counter] = shape->isInside(__randomVectorUtil(rng.get(),boundingBox));
    ++counter;
  }
  reset();
}

bool RandomPointsInShape::get(orsa::Vector & v) const {
  if (counter >= size) {
    // ORSA_DEBUG("forgot to call reset() ?");
    // reset();
    return false;
  }
  do { 
    v = __randomVectorUtil(rng,shape->boundingBox());
    if (in[counter]) {
      return true;
    }
    ++counter;
  } while (counter < size);
  ORSA_DEBUG("could not find any point inside the shape");
  v = orsa::Vector(0,0,0);
  return false;
}

orsa::Vector RandomPointsInShape::__randomVectorUtil(const orsa::RNG * rng,
						     const Box & boundingBox) {
  return Vector(boundingBox.getXMin()+(boundingBox.getXMax()-boundingBox.getXMin())*rng->gsl_rng_uniform(),
		boundingBox.getYMin()+(boundingBox.getYMax()-boundingBox.getYMin())*rng->gsl_rng_uniform(),
		boundingBox.getZMin()+(boundingBox.getZMax()-boundingBox.getZMin())*rng->gsl_rng_uniform());
}

/***/

double volume(const orsa::RandomPointsInShape * randomPointsInShape) {
  unsigned int in=0;
  orsa::Vector(v);
  randomPointsInShape->reset();
  while (randomPointsInShape->get(v)) {
    ++in;
  }
  return ((randomPointsInShape->shape->boundingBox().volume()*in)/randomPointsInShape->size);
}

orsa::Vector centerOfMass(const orsa::RandomPointsInShape * randomPointsInShape,
			  const orsa::MassDistribution * massDistribution) {
  
  osg::ref_ptr<orsa::Statistic<double> > stat_M   = new orsa::Statistic<double>;
  //
  osg::ref_ptr<orsa::Statistic<double> > stat_CMx = new orsa::Statistic<double>;
  osg::ref_ptr<orsa::Statistic<double> > stat_CMy = new orsa::Statistic<double>;
  osg::ref_ptr<orsa::Statistic<double> > stat_CMz = new orsa::Statistic<double>;
  
  orsa::Vector v;
  randomPointsInShape->reset();
  while (randomPointsInShape->get(v)) {
    const double density = massDistribution->density(v);
    if (density > 0) {
      stat_M->insert(density);
      //
      stat_CMx->insert(density*v.getX());
      stat_CMy->insert(density*v.getY());
      stat_CMz->insert(density*v.getZ());
    }
  }
  
  const orsa::Vector center_of_mass(stat_CMx->sum()/stat_M->sum(),
				    stat_CMy->sum()/stat_M->sum(),
				    stat_CMz->sum()/stat_M->sum());
  
  if (1) {
    // debug output
    const orsa::Vector center_of_mass_uncertainty(stat_CMx->averageError(),
						  stat_CMy->averageError(),
						  stat_CMz->averageError());
    ORSA_DEBUG("cm.x: %14.6e +/- %14.6e",
	       center_of_mass.getX(),
	       center_of_mass_uncertainty.getX());
    ORSA_DEBUG("cm.y: %14.6e +/- %14.6e",
	       center_of_mass.getY(),
	       center_of_mass_uncertainty.getY());
    ORSA_DEBUG("cm.z: %14.6e +/- %14.6e",
	       center_of_mass.getZ(),
	       center_of_mass_uncertainty.getZ());
  }
  
  return center_of_mass;
}

void diagonalizedInertiaMatrix(orsa::Matrix & shapeToLocal,
			       orsa::Matrix & localToShape,
			       orsa::Matrix & inertiaMatrix,
			       const orsa::Vector & centerOfMass,
			       const orsa::RandomPointsInShape * randomPointsInShape,
			       const orsa::MassDistribution * massDistribution) {
  
  shapeToLocal = orsa::Matrix::identity();
  localToShape = orsa::Matrix::identity();
  // two runs, one not rotated, one rotated
  bool rotatedRun = false;
  while (1) {
    
    osg::ref_ptr<orsa::Statistic<double> > stat_Ixx = new orsa::Statistic<double>;
    osg::ref_ptr<orsa::Statistic<double> > stat_Iyy = new orsa::Statistic<double>;
    osg::ref_ptr<orsa::Statistic<double> > stat_Izz = new orsa::Statistic<double>;
    osg::ref_ptr<orsa::Statistic<double> > stat_Ixy = new orsa::Statistic<double>;
    osg::ref_ptr<orsa::Statistic<double> > stat_Ixz = new orsa::Statistic<double>;
    osg::ref_ptr<orsa::Statistic<double> > stat_Iyz = new orsa::Statistic<double>;
    
    orsa::Vector v;
    randomPointsInShape->reset();
    while (randomPointsInShape->get(v)) {
      const double density = massDistribution->density(v);
      if (density > 0) {
	v -= centerOfMass;
	v = shapeToLocal*v;
	//
	stat_Ixx->insert(density*(v.getY()*v.getY()+v.getZ()*v.getZ()));
	stat_Iyy->insert(density*(v.getX()*v.getX()+v.getZ()*v.getZ()));
	stat_Izz->insert(density*(v.getX()*v.getX()+v.getY()*v.getY()));
	//
	stat_Ixy->insert(density*(-v.getX()*v.getY()));
	stat_Ixz->insert(density*(-v.getX()*v.getZ()));
	stat_Iyz->insert(density*(-v.getY()*v.getZ()));
      }
    }
    
    inertiaMatrix.set(stat_Ixx->average(),
		      stat_Ixy->average(),
		      stat_Ixz->average(),
		      stat_Ixy->average(),
		      stat_Iyy->average(),
		      stat_Iyz->average(),
		      stat_Ixz->average(),
		      stat_Iyz->average(),
		      stat_Izz->average());
    
    if (1) {
      // debug output    
      ORSA_DEBUG("rotated: %i",rotatedRun);
      //
      ORSA_DEBUG("Ixx:  %14.6e +/- %14.6e",
		 double(stat_Ixx->average()),
		 double(stat_Ixx->averageError()));
      ORSA_DEBUG("Iyy:  %14.6e +/- %14.6e",
		 double(stat_Iyy->average()),
		 double(stat_Iyy->averageError()));
      ORSA_DEBUG("Izz:  %14.6e +/- %14.6e",
		 double(stat_Izz->average()),
		 double(stat_Izz->averageError()));
      //
      ORSA_DEBUG("Ixy:  %14.6e +/- %14.6e",
		 double(stat_Ixy->average()),
		 double(stat_Ixy->averageError()));
      ORSA_DEBUG("Ixz:  %14.6e +/- %14.6e",
		 double(stat_Ixz->average()),
		 double(stat_Ixz->averageError()));
      ORSA_DEBUG("Iyz:  %14.6e +/- %14.6e",
		 double(stat_Iyz->average()),
		 double(stat_Iyz->averageError()));
    }
    
#warning invert??
    ORSA_DEBUG("invert??");
    orsa::principalAxis(localToShape,
		        inertiaMatrix,
			inertiaMatrix);
    orsa::Matrix::invert(localToShape,shapeToLocal);
    
#warning check largest along z, smallest along x
    ORSA_DEBUG("principal axis: check largest I along z, smallest along x");
    
    if (rotatedRun) break;
    rotatedRun = true;
  }
  
}



orsa::PaulMoment * orsa::computePaulMoment(const unsigned int order,
					   const orsa::Matrix & shapeToLocal,
					   const orsa::Matrix & localToShape,
					   const orsa::Vector & centerOfMass,
					   const orsa::RandomPointsInShape * randomPointsInShape,
					   const orsa::MassDistribution * massDistribution) {
  
  orsa::PaulMoment * pm = new orsa::PaulMoment(order);
  
  const unsigned int order_plus_one = order+1;
  //
  osg::ref_ptr<orsa::Statistic<double> > _stat_M = new orsa::Statistic<double>;
  //
  std::vector< std::vector< std::vector< osg::ref_ptr< orsa::Statistic<double> > > > > _stat_Mo(order_plus_one);
  //
  // std::vector< std::vector< std::vector< osg::ref_ptr< orsa::Statistic<double> > > > > _stat_Vo(order_plus_one);
  //
  _stat_Mo.resize(order_plus_one);
  // _stat_Vo.resize(order_plus_one);
  for (unsigned int i=0; i<order_plus_one; ++i) {
    _stat_Mo[i].resize(order_plus_one-i);
    // _stat_Vo[i].resize(order_plus_one-i);
    for (unsigned int j=0; j<order_plus_one-i; ++j) {
      _stat_Mo[i][j].resize(order_plus_one-i-j);
      // _stat_Vo[i][j].resize(order_plus_one-i-j);
      for (unsigned int k=0; k<order_plus_one-i-j; ++k) {
	_stat_Mo[i][j][k] = new orsa::Statistic<double>;
     	// _stat_Vo[i][j][k] = new orsa::Statistic<double>;
      }
    }
  }
  
  orsa::Vector v;
  randomPointsInShape->reset();
  while (randomPointsInShape->get(v)) {
    const double density = massDistribution->density(v);
    if (density > 0) {
      v -= centerOfMass;
      v = shapeToLocal*v;
      //
      _stat_M->insert(density);
      //
      for (unsigned int i=0; i<order_plus_one; ++i) {
	for (unsigned int j=0; j<order_plus_one-i; ++j) {
	  for (unsigned int k=0; k<order_plus_one-i-j; ++k) {
	    _stat_Mo[i][j][k]->insert(density*
				      int_pow(v.getX(),i)*
				      int_pow(v.getY(),j)*
				      int_pow(v.getZ(),k));
	  }
	}
      }
      
    }
  }
  
  /* Mo.resize(order_plus_one);
     Mo_uncertainty.resize(order_plus_one);
     for (unsigned int i=0; i<order_plus_one; ++i) {
     Mo[i].resize(order_plus_one-i);
     Mo_uncertainty[i].resize(order_plus_one-i);
     for (unsigned int j=0; j<order_plus_one-i; ++j) {
     Mo[i][j].resize(order_plus_one-i-j);
     Mo_uncertainty[i][j].resize(order_plus_one-i-j);
     }
     }
  */
  
  for (unsigned int i=0; i<order_plus_one; ++i) {
    for (unsigned int j=0; j<order_plus_one-i; ++j) {
      for (unsigned int k=0; k<order_plus_one-i-j; ++k) {
	// Mo[i][j][k] = _stat_Mo[i][j][k]->sum() / _stat_M->sum();
	// Mo_uncertainty[i][j][k] = _stat_Mo[i][j][k]->averageError();
	//
	// save
	// pm->setM(Mo[i][j][k],i,j,k);
	pm->setM(_stat_Mo[i][j][k]->sum()/_stat_M->sum(),
		 i,j,k);
	
      }
    }
  }
  
  return pm;
}
